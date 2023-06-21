
//
// routine to fill SGFE from solver fluxes.
//
// Basic flow:
// Loop over all grid points.
//  Compare point to left and right flux extents, for each subgrid.
//  If it's a match, copy from the Solver Flux to the Flux Estimate.
//  
//
// While it may be more direct to loop over only the subgrid extents,
// we felt looping over the grid points would ensure that the Solver Flux (a much
// larger, 3 dimensional object) was accessed with stride 1, to
// improve cache performance.  
//
// Variables:

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
int grid::FillFluxes( int NumberOfSubgrids, fluxes *SubgridFluxes[], float * Fluxes[], float dT ){

  //loop variables
  int i,j,k,dim,axis,nSg,field;

  //indexing variables

  //CellOffset  : Coordinate of first data cell (w/ ghost zones) WRT
  //              domain, in units of cells.
  //iSg, jSg,kSg: i,j,k WRT the current working subgrid;
  // SGindex    : index of subgrid cell in memory
  // PGindex    : index of parent grid cell in memory.
  int CellOffset[3], SGindex, PGindex;
  int iSg,jSg,kSg;

  int * FluxDims[3][3];

  for( dim=0;dim<3;dim++){
    //error check.  Don't think I really need it, but whatever.
    if( CellWidth[dim][0] != 0 ){
      CellOffset[dim] = nint( (GridLeftEdge[dim] - DomainLeftEdge[dim])/CellWidth[dim][0] )-GridStartIndex[dim];
    }else{
      CellOffset[dim] = 0;
    }
  }

  //Set up the flux sizes

  for(axis=0;axis<3;axis++)
    for(dim=0;dim<3; dim++){
      FluxDims[axis][dim] = new int[NumberOfSubgrids];
      for( nSg = 0;nSg < NumberOfSubgrids; nSg++){
	FluxDims[axis][dim][nSg] = 
	  SubgridFluxes[nSg]->LeftFluxEndGlobalIndex[axis][dim]-
	  SubgridFluxes[nSg]->LeftFluxStartGlobalIndex[axis][dim] + 1;
      }
    }
  //<dbg> inside the next loop
  //PGindex fill, steal from dif
  
  //loop over the grid, fill the relavent fluxes.


  for( k=GridStartIndex[2]; k<=GridEndIndex[2]; k++)
    for( j=GridStartIndex[1]; j<=GridEndIndex[1]; j++)
      for( i=GridStartIndex[0]; i<=GridEndIndex[0]; i++){


	// <dbg> fprintf(stderr,"crowd (%d,%d,%d)\n",i,j,k);
	for(axis=0;axis< GridRank; axis++){
	  for( nSg = 0;nSg < NumberOfSubgrids; nSg++){
	    
	    //Left Flux

	    iSg = i - ( SubgridFluxes[nSg]->LeftFluxStartGlobalIndex[axis][0] - CellOffset[0] );
	    jSg = j - ( SubgridFluxes[nSg]->LeftFluxStartGlobalIndex[axis][1] - CellOffset[1] );
	    kSg = k - ( SubgridFluxes[nSg]->LeftFluxStartGlobalIndex[axis][2] - CellOffset[2] );
	    
	    SGindex = iSg + FluxDims[axis][0][nSg]*(jSg + FluxDims[axis][1][nSg]*kSg);
	    PGindex = (i + GridDimension[0]*(j+GridDimension[1]*k) )*NumberOfMHDFluxes;
	    if( 
	       iSg >= 0 && iSg < FluxDims[axis][0][nSg] &&
	       jSg >= 0 && jSg < FluxDims[axis][1][nSg] &&
	       kSg >= 0 && kSg < FluxDims[axis][2][nSg] ){

	      for(field=0;field<NumberOfBaryonFields; field++)
		SubgridFluxes[nSg]->LeftFluxes[field][axis][SGindex] =
		  dT*dxI[axis]*Fluxes[axis][PGindex +MapEtoS[axis][field] ];

	    }//if

	    //Right Flux

	    iSg = i - ( SubgridFluxes[nSg]->RightFluxStartGlobalIndex[axis][0] - CellOffset[0] );
	    jSg = j - ( SubgridFluxes[nSg]->RightFluxStartGlobalIndex[axis][1] - CellOffset[1] );
	    kSg = k - ( SubgridFluxes[nSg]->RightFluxStartGlobalIndex[axis][2] - CellOffset[2] );
	    
	    SGindex = iSg + FluxDims[axis][0][nSg]*(jSg + FluxDims[axis][1][nSg]*kSg);
	    PGindex = (i + GridDimension[0]*(j+GridDimension[1]*k) + Offset[axis] )*NumberOfMHDFluxes;
	    
	    if( iSg >= 0 && iSg < FluxDims[axis][0][nSg] &&
	        jSg >= 0 && jSg < FluxDims[axis][1][nSg] &&
	        kSg >= 0 && kSg < FluxDims[axis][2][nSg] ){
	      
	      for(field=0;field<NumberOfBaryonFields; field++){
		SubgridFluxes[nSg]->RightFluxes[field][axis][SGindex] = 
		  dT*dxI[axis]*Fluxes[axis][PGindex + MapEtoS[axis][field] ];
	      }
	      
		
	      }//if

	  }//nSubgrids
	}//axis
	//</dbg> needs to be inside the loop
	
      }//i,j,k loop
  
  
  
  return SUCCESS;
}
#endif //ATHENA
