


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


int grid::MHD_FluxDifference(float dT, float * Fluxes[], float RandomForcingNormalization){


  //loop variables
  int field, dim, i,j,k,index1, index2, index3;

  int kludge = 0; //0 for the real thing.  1 for saving the flux directly
  for(k=AthStart[2]; k<=AthEnd[2]; k++)
    for(j=AthStart[1]; j<=AthEnd[1]; j++)
      for(i=AthStart[0]; i<=AthEnd[0]; i++){
	
	index1 = i + GridDimension[0]*( j + GridDimension[1]*k);

	for(field=0;field<NumberOfBaryonFields;field++){
	  if( (field == Ev[0] || field == Ev[1] || field == Ev[2]) && MHD_ReconField[0] == 0 ){
	    BaryonField[field][index1] = OldBaryonField[field][index1] * OldBaryonField[Eden][index1];
	  }else{
	    BaryonField[field][index1] = OldBaryonField[field][index1];
	  }
	}

	if( AllCellCentered == 1 )
	  for(field=0;field<3;field++)
	    CenteredB[field][index1] = OldCenteredB[field][index1];
	
	for( dim = 0; dim < GridRank; dim++){
	//for(dim=0;dim< min( GridRank,1);dim++){
	  //fprintf(stderr,"kludge: diff only doing one dim\n");
	  index2 = index1*NumberOfMHDFluxes; // + Map[ field ];
	  index3 = (index1 + Offset[dim]) * NumberOfMHDFluxes;
	  

	  //<dbg>
	  //0 for the actual update
	  //1 for Just the right fluxes.

	  switch( kludge ) { 
	  case 0:

	    for(field=0;field<NumberOfBaryonFields;field++){
	      BaryonField[field][index1] -=
		dT*dxI[dim]*(Fluxes[dim][index3+MapEtoS[dim][field]] -Fluxes[dim][index2+MapEtoS[dim][field] ]);
	      //Some error checking.
	      //<dbg>
	      if( BaryonField[field][index1] != BaryonField[field][index1] ){
		fprintf(stderr," FluxDifference: field problem, field %d (i,j,k,dim) %d %d %d %d\n",
		       field,i,j,k,dim);
		SuggestFailure = TRUE;
		return SUCCESS;
	      }
		 //</dbg>
	    }//field
	    
	    if( AllCellCentered == 1 ){
	      CenteredB[ BNum[dim][0] ][index1 ] -= 
		dT*dxI[dim]*(Fluxes[dim][index3+Sb[0] ]
			     -Fluxes[dim][index2+Sb[0] ]);
	      CenteredB[ BNum[dim][1] ][index1 ] -= 
		dT*dxI[dim]*(Fluxes[dim][index3+Sb[1] ]
			     -Fluxes[dim][index2+Sb[1] ]);
	    }
	    
	    break;

	  case 1:
	    fprintf(stderr,"kludge: only fluxes, not differences.\n");
	    dT = 1.0;
	    dxI[dim] = 1.0;
	    
	    
	    for(field=0;field<NumberOfBaryonFields;field++){

	      BaryonField[field][index1] = Fluxes[dim][index3+MapEtoS[dim][field]];
	      //dT*dxI[dim]*(Fluxes[dim][index3+MapEtoS[dim][field]] /*-Fluxes[dim][index2+MapEtoS[dim][field] ] */);
	      
	      //We store the Velocity, but the method is based on Conservative variables.
	      /*
		if( (field == Ev[0]) || (field == Ev[1]) || (field == Ev[2]) && 0==1 )
		BaryonField[field][index1] /= BaryonField[Eden][index1];
	      */
	      
	    }
	    
	    if( AllCellCentered == 1 ){
	      CenteredB[ BNum[dim][0] ][index1 ] =  
		dT*dxI[dim]*(Fluxes[dim][index3+Sb[0] ]
			     /*-Fluxes[dim][index2+Sb[0] ]*/ );
	      CenteredB[ BNum[dim][1] ][index1 ] = 
		dT*dxI[dim]*(Fluxes[dim][index3+Sb[1] ]
			     /*-Fluxes[dim][index2+Sb[1] ]*/);
	    }
	  default:
	    break;
	  }//switch	  
	  
	  
	}//dim

	if (SelfGravity || UniformGravity || PointSourceGravity){
	  for( dim=0; dim<GridRank; dim++){
	    BaryonField[ Ev[dim] ][index1] += dT*AccelerationField[dim][index1]*
	      0.5*(BaryonField[ Eden ][index1 ] + OldBaryonField[ Eden ][index1 ]);

	    if( EquationOfState == 0 ){
	      BaryonField[Eeng][index1 ] += dT*AccelerationField[dim][index1]*
		0.5*(BaryonField[ Eden ][index1 ]*BaryonField[ Ev[dim] ][index1]
		     + OldBaryonField[ Eden ][index1 ]*OldBaryonField[ Ev[dim] ][index1]);
										  
	    }//eos = 0

	  }//dim
	}//gravity

	//We store the Velocity, but the method is based on Conservative variables.
	//This loop is written stupid.  It should get fixed.
	if( kludge ==  0 ){ //<dbg>	
	  for( field=0; field<NumberOfBaryonFields; field++){
	    if( ( (field == Ev[0]) || (field == Ev[1]) || (field == Ev[2]) ) && MHD_ReconField[0] == 0 ){
	      BaryonField[field][index1] /= BaryonField[Eden][index1];
	    }

	  }//field
	  }//<dbg>

	if( TRUE == RandomForcing ){
	  for(dim = 0; dim < GridRank; dim++){
	    if( 0 == EquationOfState )
	      BaryonField[ Seng ][index1] +=
		BaryonField[ Sden ][ index1 ]*
		(BaryonField[ Sv[dim] ][index1]*
		 RandomForcingField[dim][index1]*
		 RandomForcingNormalization+
		 0.5*RandomForcingField[dim][index1]*
		 RandomForcingField[dim][index1]* 
		 RandomForcingNormalization*RandomForcingNormalization);

	    BaryonField[ Sv[dim] ][index1] += 
	      RandomForcingField[dim][index1] * RandomForcingNormalization;
					       

									   
	  }//dim
	}//rf

	//Density Floor.  Kind of experimental.
	/*
	  if( BaryonField[ Sden ][ index1 ] <= tiny_density ){
	  fprintf( stderr, "tiny density: (%d, %d,%d)\n", i,j,k);
	  BaryonField[ Sden ][index1] = tiny_density ;
	  }
	*/
      }//i,j,k
  
  return SUCCESS;
}
#endif //ATHENA
