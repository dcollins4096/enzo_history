/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRID CLASS (CLEAR/ALLOCATE BOUNDARY FLUXES FOR THIS GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// Clear the boundary fluxes (allocate first if necesary).
//  Note that we assume the left and right flux boundaries are the same size.
//   
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void grid::ClearBoundaryFluxes()
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return;
    JBMEM_MESSAGE(MyProcessorNumber,"jb: BeforeClearBoundary");  
  int dim, i, field, size, Dim[3]={1,1,1}, EDim[3]={1,1,1};
  
  for (dim = 0; dim < GridRank; dim++) {
    
    /* compute size for this dim's face */
    
    size = 1;
    
    for (i = 0; i < GridRank; i++){
      Dim[i] = BoundaryFluxes.LeftFluxEndGlobalIndex[dim][i] -
	BoundaryFluxes.LeftFluxStartGlobalIndex[dim][i] + 1;
      size *= Dim[i];
    }
    /* allocate if necesary, and clear */

    for (field = 0; field < NumberOfBaryonFields; field++) {
      if (BoundaryFluxes.LeftFluxes[field][dim] == NULL) {
	BoundaryFluxes.LeftFluxes[field][dim] = new float[size];
	BoundaryFluxes.RightFluxes[field][dim] = new float[size];
      }
      for (i = 0; i < size; i++) {
	BoundaryFluxes.LeftFluxes[field][dim][i] = 0.0;
	BoundaryFluxes.RightFluxes[field][dim][i] = 0.0;
      }
    }
    /* set unused pointers to NULL */
  
    for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	 field++)
      for (i = 0; i < MAX_DIMENSION; i++) {
	BoundaryFluxes.LeftFluxes[field][i]  = NULL;
	BoundaryFluxes.RightFluxes[field][i] = NULL;
      }
  
    for(field=0;field<3;field++){
      
      if(MHD_Used == TRUE && field != dim){	
	size = 1;
	
	for(int j=0;j<3;j++){	  
	  EDim[j] = Dim[j] + ( ( ( j != dim ) && ( j != field ) ) ? 1:0 );
	  size *= EDim[j];
	}
	if( BoundaryFluxes.LeftElectric[field][dim] == NULL )
	  BoundaryFluxes.LeftElectric[field][dim]  = new float[size];
	if( BoundaryFluxes.RightElectric[field][dim] == NULL)
	  BoundaryFluxes.RightElectric[field][dim] = new float[size];
	BoundaryFluxes.ElectricSize[field][dim]=size;

	for(i=0;i<size; i++){
	  BoundaryFluxes.LeftElectric[field][dim][i] = 0.0;
	  BoundaryFluxes.RightElectric[field][dim][i] = 0.0;
	}
	  
      }
      
    }//E field
  } // end loop over dims

  /* set unused pointers to NULL */

  for (dim = GridRank; dim < MAX_DIMENSION; dim++)
    for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
      BoundaryFluxes.LeftFluxes[field][dim]   = NULL;
      BoundaryFluxes.RightFluxes[field][dim] = NULL;
    }
    JBMEM_MESSAGE(MyProcessorNumber,"jb: AfterClearBoundary");

}
