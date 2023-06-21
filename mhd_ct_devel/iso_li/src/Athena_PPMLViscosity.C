#ifdef PPML
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"

#include "AthenaObj.h"
#include "PPML.h"

#include "DaveTools.h"

//
// Wrapper to add the Athena Viscosity module to the PPML flux.
//

int AthenaObj::PPMLViscosity(float * FluxX,float * FluxY, float * FluxZ, int SweepDirection){
  /*
    grid::AddDiffusion( FluxX, FluxY, FluxX ){
    Define MapEtoS  steal from DaveThena
    Declare State
    float * FluxArray[3] = {fluxX ...}
    for( k,j,k)
    BaryonIndex = i + nx( j + ny k )
    FluxIndex = i + nx( j + ny ( k + nf field ) //wait, shit-- compute inside.
    for( dim )
     cast state
     MHD_Diffusion( Flux[dim], Null,Null, BaryonIndex, Dim,
                    DiffusionMethod, dT, state);
    } 
  */

  float * FluxArray[3] = {FluxX, FluxY, FluxZ};
  int i,j,k,dim,field, index[5];

  //For directionally unsplit runs, loop over all 3 components.
  //For directionally split runs, only do the 'relevant' one.
  int dimStart = ( (MHD_DirectionalSplitting == 1) ? SweepDirection : 0 );
  int dimEnd =  ( (MHD_DirectionalSplitting == 1) ? SweepDirection+1 : GridRank );
  // In 1 or 2 dimensional sims, the boundary faces aren't necessary.
  int EndX=1, EndY=1,EndZ=1;

  if( GridRank < 3 ){
    EndZ = 0;
    if( GridRank < 2 ){
      EndY = 0;
    }
  }

  for(k=GridStartIndex[2]; k<= GridEndIndex[2]+EndZ;k++)
  for(j=GridStartIndex[1]; j<= GridEndIndex[1]+EndY;j++)
  for(i=GridStartIndex[0]; i<= GridEndIndex[0]+EndX;i++){
    for(dim=dimStart;dim<dimEnd; dim++){    
      
      //Only the longitudinal face, for a given flux, needs to be fixed.
      if( dim == 0 && (k == AthEnd[2]+1) || (j== AthEnd[1]+1) )
	continue;
      if( dim == 1 && (i == AthEnd[0]+1) || (k== AthEnd[2]+1) )
	continue;
      if( dim == 2 && (i == AthEnd[0]+1) || (j== AthEnd[1]+1) )
	continue;

      //Rotate the vectors so that V[0] is along dim.

      RotateState(dim);

      //Loop Variables

      // map of indicies:   |  3  |  1  | 0 |  2  |  4  |
      //                    | i-2 | i-1 | i | i+1 | i+2 |
      //                    Swap i with j or k as needed.
      //                    This swap is controlled by the Offset variable.
      //                    cell(i  ,j+1,k  ) = i + nx*(j+1) + nx*ny(k)  = i + nx*j + ny*k + nx
      //                                      = index[0] + nx
      //                    cell(i  ,j  ,k+2) = i + nx*j + nx*ny*(k+2)
      //                                      = index[0] + 2*nx*ny
      //                    etc.


      //Compute the various indicies.
      index[0] = i + GridDimension[0]*(j + GridDimension[1]*k );
      index[1] = index[0] - Offset[dim];
      index[2] = index[0] + Offset[dim];
      index[3] = index[0] - 2*Offset[dim];
      index[4] = index[0] + 2*Offset[dim];

      if( MHD_Viscosity(FluxArray[dim], NULL, NULL, index, 
		    dim, MHD_DiffusionMethod[0], dtFixed) < 0 )
	return FAIL;;
      /*
float AthenaObj::MHD_Viscosity(float * Fluxes, float * Lhs, float * Rhs, int * index, 
			  int dim, int MHD_DiffusionMethodInput, float dT){
      */
    }//dim
  }//i,j,k

  return SUCCESS;

}

#endif //PPML
