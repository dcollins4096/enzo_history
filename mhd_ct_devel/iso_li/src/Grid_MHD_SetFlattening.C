#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "pout.h"
#include "MHD_Athena.h"

#ifdef ATHENA

int grid::MHD_SetFlattening(){

  //temporary.  References to pressure should be replaced 
  //by C^2 \rho for isothermal EOS.
  if( EquationOfState == 1 ){
    fprintf(stderr,"SetFlattening: Needs calls to MHD_Pressure replaced by C^2*Density for Isothermal\n");
    return FAIL;
  }
  //
  if( MyProcessorNumber != ProcessorNumber )
    return SUCCESS;

  if( MHD_Flattening == 0 )
    return SUCCESS;

  // things you need:
  int i,j,k,dim, field, indexF, index[5];
  // map of indicies:   |  3  |  1  | 0 |  2  |  4  |
  //                    | i-2 | i-1 | i | i+1 | i+2 |
  //                    Swap i with j or k as needed.
  //                    This swap is controlled by the Offset variable.
  //                    cell(i  ,j+1,k  ) = i + nx*(j+1) + nx*ny(k)  = i + nx*j + ny*k + nx
  //                                      = index[0] + nx
  //                    cell(i  ,j  ,k+2) = i + nx*j + nx*ny*(k+2) 
  //                                      = index[0] + 2*nx*ny 
  //                    etc.

  // In 1 or 2 dimensional sims, the boundary faces aren't necessary.
  int EndX=1, EndY=1,EndZ=1;

  if( GridRank < 3 ){
    EndZ = 0;
    if( GridRank < 2 ){
      EndY = 0;
    }
  }
  
  //variables I'll move later
    
  float dp1, dp2, Pp1,Pm1,Pp2,Pm2;      //pressure jumps, p(i+1), p(i-1), p(+-2)
  float divV;                 //divergence of V
  float wj;                   //w_j from CW84
  float EpsilonCW = 0.33;     //Epsilon from CW84, A.1
  float EpsilonB  = 6.923e-2; //Epsilon from Balsara94, B3
  float PressureRatio;        //Pressure rato used in \tilde{f}_j, CW84 A.2
  float Omega, OmegaTilde, Omega1 = 0.75, Omega2 = 10, Nu1 = 0.1; //more parameters from CW84: steepness of shock
  float Sigma, SigmaTilde, Sigma1 = 0.5, Sigma2 = 1.0; //Sigma is roughly the Wave Strength
  float Kappa, KappaTilde, Kappa1 = 2, Kappa2 = 0.01;  //Kappa is the wavelength of the noise (ish)
  float dE1, dE2, Eratio;
  float We, Wj; //Shock speed approximation.  Note that I have NOT updated this for MHD,
                //this shock approximation comes from the HYDRO Rankine Hugoniot relations, not MHD.
  float Cs;     //sound speed.
  int UpWind, DownWind; // j+2s, j-2s.
  //var

  for(k=AthStart[2]-2; k<= AthEnd[2]+2;k++)
    for(j=AthStart[1]-2; j<= AthEnd[1]+2;j++)
      for(i=AthStart[0]-2; i<= AthEnd[0]+2;i++){

	//set up the state variable.  This will be permuted in the dimension loop.
	for( field=0; field<3; field++)
	  State[  MapEtoS[0][ Ev[field]] ] = BaryonField[ Ev[field] ];

	divV = 
	  (  State[ Sv[0] ][ index0(i+1,j,k) ] - State[ Sv[0] ][ index0(i-1,j,k) ] );
	if( GridRank >= 2 )
	  divV += State[ Sv[1] ][ index0(i,j+1,k) ] - State[ Sv[1] ][ index0(i,j-1,k) ];
	if( GridRank >= 3 )	
	  divV += State[ Sv[2] ][ index0(i,j,k+1) ] - State[ Sv[2] ][ index0(i,j,k-1) ];

	for( dim = 0; dim < GridRank; dim++){
	  
	  if( dim == 0 )
	    if( (k == AthEnd[2]+1) || (j== AthEnd[1]+1) )
	      continue;
	  if( dim == 1 )
	    if( (i == AthEnd[0]+1) || (k== AthEnd[2]+1) )
	      continue;
	  if( dim == 2 )
	    if( (i == AthEnd[0]+1) || (j== AthEnd[1]+1) )
	      continue;
	  
	  index[0] = i + GridDimension[0]*(j + GridDimension[1]*k );
	  index[1] = index[0] - Offset[dim];
	  index[2] = index[0] + Offset[dim];
	  index[3] = index[0] - 2*Offset[dim];
	  index[4] = index[0] + 2*Offset[dim];
	  indexF = dim + GridRank*(i + GridDimension[0]*(j + GridDimension[1]*k));

	  
	  //Cast the State pointer to the right Vector field.  Hopefully
	  //this won't be too violent for the Cache performance.
	  //ick: Sv[0,1,2] = kinetic field (velocity or momentum) that the Solver sees.
	  //     Ev[0,1,2] = kinetic field that Enzo sees.
	  //     MapEtoS[dim][0,1,2] = the map (which is dim dependant) from one to the other.
	  //     Sb[0,1]   = transvers field components for the Solver
	  //     BNum[dim][0,1] = transverse fields from Enzo.  


	  for( field=0; field<3; field++)
	    State[  MapEtoS[dim][ Ev[field]] ] = BaryonField[ Ev[field] ];
	  for( field=0; field<2; field++)
	    State[ Sb[field] ] = CenteredB[ BNum[dim][field] ];
	  
	  if( AllCellCentered == 1 ){
	    State[ Sb[2] ] = CenteredB[dim];
	  }

	  Pp1 =  MHD_Pressure[ index[2] ];
	  Pm1 =  MHD_Pressure[ index[1] ];
	  Pp2 =  MHD_Pressure[ index[4] ];
	  Pm2 =  MHD_Pressure[ index[3] ];
	  dp1 = Pp1 - Pm1;
	  dp2 = Pp2 - Pm2;


	  wj = 0;
	  if( EquationOfState == 0 ){
	    if( divV < 0 && fabs(dp1)/min(Pp1,Pm1) > EpsilonCW )
	      wj = 1;
	  }else{
	    //see Balsara, ApJ 1994, B3
	    if( divV < 0 && (divV*divV)/(IsothermalSoundSpeed*IsothermalSoundSpeed) > EpsilonB )
	      wj = 1;
	  }
	  
	  if( dp2 / min(Pp2,Pm2) < EpsilonCW ){
	    PressureRatio = 1;
	  }else{
	    PressureRatio = fabs(dp1/dp2);	  
	  }

	  ShockDirection[indexF] = (( Pp1 < Pm1 ) ? 1 : -1);	  	  

	  switch( MHD_Flattening ){
	  case 1:

	    FlatteningField[indexF] = wj*max( 0, min( (PressureRatio - Omega1) * Omega2 , 1) );
	    break;
	  case 2:
	    fprintf(stderr," Due to lack of Lagrangian MHD, Flattening2 isn't an option.  Only 1 & 3.\n");
	    fprintf(stderr," This numbering is to be consistent with CW84 and PPM in Enzo.\n");
	    return FAIL;
	  case 3:
	    
	    // Wave Strength
	    SigmaTilde = wj * dp2 / min(Pp2,Pm2);
	    Sigma = max( 0, ( ( SigmaTilde - Sigma1 )/(SigmaTilde + Sigma2 ) ) );

	    //Steepness of shock
	    if( EquationOfState == 0 ){
	      dE1 = BaryonField[ Eeng ][ index[2] ] - BaryonField[ Eeng ][ index[1] ];
	      dE2 = BaryonField[ Eeng ][ index[4] ] - BaryonField[ Eeng ][ index[3] ];
	    }else{
	      dE1 = 0;
	      dE2 = 0;
	    }
	    if( dE2 < tiny_number )
	      Eratio = 0;
	    else
	      Eratio = dE1/dE2;

	    OmegaTilde = max( PressureRatio, Eratio );
	    Omega = max( 0, Omega2*(OmegaTilde - Omega1 ) ); //copying what's done in Enzo.
	                                                     //Suspected typo in CW84.

	    Wj = sqrt( (max(Pp2,Pm2) + ( (EquationOfState == 0 ) ? 0.5*(Gamma-1)*(Pp2+Pm2) : 0 ) )/
		       (max( 1/State[ Sden ][ index[3] ], 1/State[ Sden ][ index[4] ] )) );

	    //Upwind = j + 2*s_j
	    UpWind =   ( (ShockDirection[ indexF ] == 1 ) ? index[4] : index[3] );
	    DownWind = ( (ShockDirection[ indexF ] == 1 ) ? index[3] : index[4] );

	    We = ShockDirection[indexF]* Wj/( State[ Sden][ DownWind ] ) 
	      + State[ Sv[0] ][ DownWind ] + tiny_number;

	    //Note that here, MHD_Pressure is actually the Pressure.
	    Cs = ( (EquationOfState == 0 )? sqrt( Gamma*MHD_Pressure[ UpWind ]/State[Sden][ UpWind ] ) :
		   IsothermalSoundSpeed );
	    KappaTilde = fabs( (We - State[ Sv[0] ][ UpWind ] + ShockDirection[ indexF ] * Cs)/We );
	    Kappa = max( 0, (KappaTilde - Kappa1 )/(KappaTilde + Kappa2 ) );

	    FlatteningField[indexF] = min( wj * Omega, min(wj*Sigma, Kappa));
	    break;
	  }//switch
	  
	}//dim
      }//i,j,k

  

  return SUCCESS;
}

#endif //ATHENA
