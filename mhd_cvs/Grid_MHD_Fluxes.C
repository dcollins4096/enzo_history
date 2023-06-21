

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
#include "MHDFluid.h"
#ifdef ATHENA


int MHD_1dFlux(float * Fluxes, float * Lhs);
//int MHD_Roe1(float * Fluxes, Fluid * L, Fluid * R );
int MHD_Roe1(float * Fluxes, float * Lhs, float * Rhs );
int MHD_ROE_NEW(float * Fluxes, Fluid * L, Fluid * R );
int MHD_RiemannFirst(float * Fluxes, Fluid * L, Fluid * R );
int MHD_HLLE(float * Fluxes, Fluid * L, Fluid * R );
int MHD_HLLC(float * Fluxes, Fluid * L, Fluid * R );
int MHD_HLLD(float * Fluxes, Fluid * L, Fluid * R );
int MHD_HLLD_Sergey(float * Fluxes, Fluid * L, Fluid * R );
int MHD_TestFlux(float * Fluxes, Fluid * L, Fluid * R );
//
// MHD_Recon and MHD_Riemann are the selected riemann solver and reconstruction
// method.  These may be altered by the solution sensitive switch, depending on 
// options selected.  
//

int grid::MHD_Fluxes(float * Fluxes[], int MHD_ReconInput, int MHD_RiemannInput, 
		     int MHD_DiffusionMethodInput, float dT ){

  //Left and right hand states, for the riemann solver.
  float Lhs[MAX_MHD_WAVES +1],Rhs[MAX_MHD_WAVES +1];

  //The new object:  I plan on rolling the entire solver into this object,
  //for now it's only used in the HLLE, HLLC solvers.
  Fluid L, R;

  //Used to check for failure of various solvers.
  int TheSolverError = SUCCESS;

  // In 1 or 2 dimensional sims, the boundary faces aren't necessary.
  int EndX=1, EndY=1,EndZ=1;

  if( GridRank < 3 ){
    EndZ = 0;
    if( GridRank < 2 ){
      EndY = 0;
    }
  }

  //Loop Variables
  int i,j,k,dim,field, index[5];
  // map of indicies:   |  3  |  1  | 0 |  2  |  4  |
  //                    | i-2 | i-1 | i | i+1 | i+2 |
  //                    Swap i with j or k as needed.
  //                    This swap is controlled by the Offset variable.
  //                    cell(i  ,j+1,k  ) = i + nx*(j+1) + nx*ny(k)  = i + nx*j + ny*k + nx
  //                                      = index[0] + nx
  //                    cell(i  ,j  ,k+2) = i + nx*j + nx*ny*(k+2) 
  //                                      = index[0] + 2*nx*ny 
  //                    etc.

  //
  // Loop over all cells.  
  // Compute Interface State
  // Solve Riemann Problem.
  //
  // Since this loop is computing fluxes, it's really best to think
  // of this loop as a loop over cell Interfaces. 
  // 

  for(k=AthStart[2]; k<= AthEnd[2]+EndZ;k++)
    for(j=AthStart[1]; j<= AthEnd[1]+EndY;j++)
      for(i=AthStart[0]; i<= AthEnd[0]+EndX;i++){

	//<dbg>
	//if( i == 12 && j == 21 && k == 40 )
	//fprintf( stderr," dude, I'm really going to tank.\n");
	//</dbg>

	for( dim = 0; dim < 3; dim++){
	//for(dim=0;dim< min( GridRank,1);dim++){


	  //fprintf(stderr,"kludge: lin only doing one dim.\n");
	  //Don't get faces outside of volume: 

	  if( dim == 0 )
	    if( (k == AthEnd[2]+1) || (j== AthEnd[1]+1) )
	      continue;
	  if( dim == 1 )
	    if( (i == AthEnd[0]+1) || (k== AthEnd[2]+1) )
	      continue;
	  if( dim == 2 )
	    if( (i == AthEnd[0]+1) || (j== AthEnd[1]+1) )
	      continue;


	  //Compute the various indicies.
	  index[0] = i + GridDimension[0]*(j + GridDimension[1]*k );
	  index[1] = index[0] - Offset[dim];
	  index[2] = index[0] + Offset[dim];
	  index[3] = index[0] - 2*Offset[dim];
	  index[4] = index[0] + 2*Offset[dim];
	  
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

	  //
	  // Solution Sensitive Switch.
	  //

	  MHD_SSSr(dim,index, MHD_ReconInput, MHD_RiemannInput);

	  //
	  // Reconstruction on Baryon Fields.
	  //

	  switch( ReconLocal ){
	    
	    //piecewise constant
	  case 0:

	    for( field=0; field< NumberReconstructed; field ++ ){
	      //Piecewise Constant.
	      Lhs[ field ] = State[field][ index[1] ];
	      Rhs[ field ] = State[field][ index[0] ];
	    }

	    break;
	  case 1:
	    //Primitive Reconstruction.
	    if( MHD_LinearSlope(Lhs, Rhs, index, dim) == FAIL ){
	      fprintf( stderr,"MHD_Fluxes: error in MHD_PLM\n");
	      return FAIL;
	    }

	    break;

	  case 2:

	    //Hancock, method 1.  (see Grid_MHD_Hancock.C for details)
	    if( MHD_Hancock1(Lhs, Rhs, index, dim, dT, 
			     MagneticField[dim][ indexba(i,j,k,dim) ]) == FAIL ){
	      fprintf( stderr,"MHD_Fluxes: error in MHD_PLM\n");
	      return FAIL;
	    }

	    break;

	  case 3:

	    //Hancock, method 2.  (see Grid_MHD_Hancock.C for details)
	    if( MHD_Hancock2(Lhs, Rhs, index, dim, dT,
			     MagneticField[dim][ indexba(i,j,k,dim) ]) == FAIL ){
	      fprintf( stderr,"MHD_Fluxes: error in MHD_PLM\n");
	      return FAIL;
	    }

	    break;
	  case -1:
	    //For Rank<3, one dimension will only have 1 set of zones, so only Rhs makes any sense.
	    for( field=0; field< NumberReconstructed; field ++ )
	      Rhs[ field ] = State[field][ index[0] ];

	    break;
	  default:
	    fprintf(stderr," need to define what you mean by MHD_Recon = %d\n", ReconLocal);
	    return FAIL;
	    break;
	  }//reconstruction switch
	  
	  //If we're using the cell centered field for the Logitudinal
	  //field, it needs to be reconstructed (done above). 
	  //Otherwise, we use the FaceCentered magnetic field (here.)
	  if( AllCellCentered == 0 ){
	    Lhs[ Sb[2] ] = MagneticField[dim][ indexba(i,j,k,dim) ];
	    Rhs[ Sb[2] ] = MagneticField[dim][ indexba(i,j,k,dim) ];
	  }


	  //
	  // Flattening.
	  //
	  if( ReconLocal != 0 && ReconLocal != 2){ //Not done for piecewise constant
	    if( MHD_Flattening != 0 ){
	      if( MHD_Flatten(Lhs, Rhs, index, dim) == FAIL )
		{fprintf(stderr, "error with MHD_Flatten\n"); return FAIL;}
	    }
	  }

	  //Once the other solvers have been Objectified, the fill commands should be
	  //moved higher in the chain, so the reconstruction is done on Fluid objects
	  if( ReconLocal != -1 ) L.Fill( Lhs );
	  R.Fill( Rhs );

	  //
	  // Riemann Solver.
	  //
	  
	  
	  if( RiemannLocal == 0 || RiemannLocal == 1 ){
	    fprintf(stderr, "MHD_Riemann = %d not defined for this HydroMethod\n", 
		    RiemannLocal);
	    return FAIL;
	  }
	  
	  if( RiemannLocal == 5 )
	    //if( MHD_Roe1(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), &L, &R) == FAIL )
	    TheSolverError = MHD_Roe1(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), Lhs, Rhs);

	  
	  if( RiemannLocal == 3 )
	    //For 2 + 1/2 d, or 1+1/2+1/2 d, only the 1d flux is needed, since the 
	    //higher order terms are mostly proportional to gradients.
	    TheSolverError = MHD_1dFlux(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), Rhs);
	  
	  if( RiemannLocal == 4 ){
	    TheSolverError = MHD_RiemannFirst(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), &L, &R);
	  }

	  if( RiemannLocal == 2 ||
	      RiemannLocal == 8 ||
	      RiemannLocal == 9 )
	    TheSolverError = MHD_ROE_NEW(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), &L, &R);

	  //<dbg>
	  if( CheckForSolverProblems == TRUE )
	    SolverProblems[dim][ index[0] ] = TheSolverError;
	  //</dbg>
	  
	  if( ( RiemannLocal == 6 ) ||
	      ( RiemannLocal == 8 &&
		(TheSolverError == 2 || TheSolverError == FAIL ) ) ){
	    //if( RiemannLocal == 8 ) fprintf(stderr,"parachute\n");
	    TheSolverError = MHD_HLLE(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), &L, &R);
	  }
	  if( ( RiemannLocal == 7 ) ||
	      ( RiemannLocal == 9 &&
		(TheSolverError == 2 || TheSolverError == FAIL ) ) ){
	    TheSolverError =  MHD_HLLC(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), &L, &R);
	    if( TheSolverError != FAIL ){
	      SuggestFailure = FALSE;
	    }else{
	      fprintf(stderr,"The backup solver failed, too.\n");
	    }
	  }

	  if( RiemannLocal == 10 ){
	    TheSolverError =  MHD_HLLD(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), &L, &R);
	  }

	  if( RiemannLocal == 11 ){
	    TheSolverError =  MHD_HLLD_Sergey(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), &L, &R);
	  }

	  //fprintf(stderr,"arse %.13g\n",Fluxes[dim][ index[0] * NumberOfMHDFluxes + Sv[2] ]);
	  if( RiemannLocal == 12 ){
	    TheSolverError =  MHD_TestFlux(&( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), &L, &R);

	  }

	  if( RiemannLocal > 12 ){
	    fprintf(stderr," MHD_Riemann = %d Not defined. Pick again\n", RiemannLocal);
	    return FAIL;
	  }

	  if( TheSolverError == FAIL ){
	    fprintf(stderr, "Solver error: error code %d \n",TheSolverError );
	    fprintf(stderr, "   (he he! my first cryptic error code!)\n");
	    return FAIL;
	  }

	  if( SuggestFailure == TRUE ) {
	    fprintf(stderr,"  Looks like the riemann solver fucked up.  That's too bad for you!\n");
	    fprintf(stderr," Just to be nice, here's some information.\n");
	    fprintf(stderr,"  LHS       RHS\n");
	    for(int Clown=0;Clown<NumberReconstructed; Clown++)
	      fprintf(stderr," %f    %f\n", Lhs[Clown], Rhs[Clown]);
	    fprintf(stderr, "(i,j,k,dim) = (%d,%d,%d,%d)\n",i,j,k,dim);

	    //<dbg>
	    if( CheckForSolverProblems == TRUE )
	      SolverProblems[dim][ index[0] ] = 0;
	    //</dbg>

	    return SUCCESS;
	  }

	  //
	  // Viscosity
	  //

	  if( MHD_DiffusionMethodInput > 0 && dim < GridRank )
	    if( MHD_Diffusion( &( Fluxes[dim][ index[0] * NumberOfMHDFluxes ] ), Lhs, Rhs, 
			       index, dim, MHD_DiffusionMethodInput, dT)
		== FAIL )
	      {fprintf(stderr,"Problem with MHD_Diffusion.\n"); return FAIL;}

	}//dim

      }//i,j,k
  
  return SUCCESS;
}

inline void grid::MHD_SSSr(int dim, int * index, int MHD_ReconInput, int MHD_RiemannInput){

  ReconLocal = MHD_ReconInput;
  RiemannLocal = MHD_RiemannInput;

  if( GridRank < 3  ){
    if( dim==2 ){
      ReconLocal = -1;   //piecewise constant, only right state.  For 2.5d
      RiemannLocal = 3; //Only magnetic flux, nothing fancy.  Calls MHD_Zflux2d
    }
    if( GridRank < 2 && dim == 1 ){
      ReconLocal = -1;   //piecewise constant, only right state.  For 2.5d
      RiemannLocal = 3; //Only magnetic flux, nothing fancy.  Calls MHD_Zflux2d
    }
  }
}
#endif //Athena
