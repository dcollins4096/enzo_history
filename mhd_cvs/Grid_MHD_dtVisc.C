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


// Limits the timestep based on the viscous criterion.
// From VonNeuman (sine wave) analysis, dt < 0.5 dX^/D,
// where D is the diffusion constant.
// In our case, D = q \rho \Delta V dX, so
// dt < 0.5 dX / (q \rho \Delta V)
//
// In general, we follow what was done in Colella and Woodward 1984 
// for the viscosity.  Nu = - min( \Div V, 0) * q
// Nu is used in place of \Delta V in the above.
//
// Currently this code used "diffusion method 3," which doesn't rely on Left and Right
// states.
int grid::MHD_dtVisc(float &dT){
  

  
  float MaxNu = 0.0;




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


  //Set up various variables needed by the Athena Unsplit machenery.
  MHD_AthenaSetup();

  State[ Sden ] = BaryonField[ Eden ];

  // In 1 or 2 dimensional sims, the boundary faces aren't necessary.
  int EndX=1, EndY=1,EndZ=1;
  if( GridRank < 3 ){
    EndZ = 0;
    if( GridRank < 2 ){
      EndY = 0;
    }
  }



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
	for( dim = 0; dim < GridRank; dim++){

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

	  //note that the last argument is "timestep," which goes in only as a failsafe.
	  
	  MaxNu = max(MHD_Diffusion(NULL,NULL,NULL,index, dim, 3, 0)
		      /CellWidth[dim][0],
		      MaxNu);

	  if( MaxNu == -1 ){
	    fprintf(stderr,"MHD_dtVisc: MHD_Diffusion failed.\n");
	    dT = -1;
	    return FAIL;
	  }
	  
	}//dim

	if( MaxNu > tiny_number )
	  dT = min( dT, 0.5/(State[ Sden ][ index[0] ] * MaxNu ) ); 
      }//i,j,k
  
  return SUCCESS;
}

#endif //ATHENA
