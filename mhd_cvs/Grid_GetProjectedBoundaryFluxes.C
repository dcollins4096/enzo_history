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
/  GRID CLASS: PROJECT (DOWNSAMPLE) THIS GRIDS FLUXES BY GIVEN REFINEMENT
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  NOTE: This routine assumes that the fluxes structure and current grid 
/        have the same baryon fields.
/
************************************************************************/

// Use the refinement factors in the arguement to down sample the fluxes
//   and put them into the fluxes structure in the argument.  Also,
//   fill out the rest of the fluxes structure.

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int CommunicationSendFluxes(fluxes *Fluxes, int ToProc, int NumberOfFields,
			    int Rank);
int CommunicationReceiveFluxes(fluxes *Fluxes, int FromProc, 
			       int NumberOfFields, int Rank);
void DeleteFluxes(fluxes *Fluxes);


int grid::GetProjectedBoundaryFluxes(grid *ParentGrid, fluxes &ProjectedFluxes)
{

  /* Return if this doesn't involve us. */

  if (MyProcessorNumber != ProcessorNumber && 
      MyProcessorNumber != ParentGrid->ProcessorNumber)
    return SUCCESS;


  if (CommunicationDirection == COMMUNICATION_SEND &&
      (MyProcessorNumber == ParentGrid->ProcessorNumber || 
       ProcessorNumber == ParentGrid->ProcessorNumber))
    return SUCCESS;

  if (CommunicationDirection == COMMUNICATION_RECEIVE &&
      MyProcessorNumber != ParentGrid->ProcessorNumber &&
      ProcessorNumber != ParentGrid->ProcessorNumber)
    return SUCCESS;

  /* Compute the subgrid's refinement factors (for step #19) */

  int RefinementFactors[MAX_DIMENSION];
  ParentGrid->ComputeRefinementFactors(this, RefinementFactors);

  int i, j, k, i1, j1, k1, dim, coord, Dims[3], field, size;
  int index1, index2;
  int ProjectedDims[MAX_DIMENSION];

  if (NumberOfBaryonFields > 0) {

    /* fill out Flux indicies */

    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < GridRank; i++) {
	ProjectedFluxes.LeftFluxStartGlobalIndex[dim][i] =
	BoundaryFluxes.LeftFluxStartGlobalIndex[dim][i]/ RefinementFactors[i];
	ProjectedFluxes.LeftFluxEndGlobalIndex[dim][i] = 
	BoundaryFluxes.LeftFluxEndGlobalIndex[dim][i] / RefinementFactors[i];
	ProjectedFluxes.RightFluxStartGlobalIndex[dim][i] = 
	BoundaryFluxes.RightFluxStartGlobalIndex[dim][i]/RefinementFactors[i];
	ProjectedFluxes.RightFluxEndGlobalIndex[dim][i] = 
	BoundaryFluxes.RightFluxEndGlobalIndex[dim][i]/RefinementFactors[i];
      }

    /* loop over all dimensions */

    for (dim = 0; dim < GridRank; dim++) {

      /* compute size of flux region */

      Dims[0] = Dims[1] = Dims[2] = 1;

      for (i = 0; i < GridRank; i++)
	Dims[i] = BoundaryFluxes.LeftFluxEndGlobalIndex[dim][i] - 
	          BoundaryFluxes.LeftFluxStartGlobalIndex[dim][i] + 1;

      /* compute ProjectedFlux dimensions (and TotalRefinement). */

      int TotalRefinement = 1;
      for (i = 0; i < 3; i++){
#ifdef DC_COSMOLOGY_FLUX
	  TotalRefinement *= RefinementFactors[i];
#endif
	if (i != dim) {
	  ProjectedDims[i] = Dims[i]/RefinementFactors[i];
#ifndef DC_COSMOLOGY_FLUX
	  TotalRefinement *= RefinementFactors[i];
#endif
	}
	else
	  ProjectedDims[i] = 1;
      }
      size = ProjectedDims[0]*ProjectedDims[1]*ProjectedDims[2];

      /* compute the fraction of each (current) grid flux cell
         that each subgrid's flux cell occupies.  Er, whatever. */

      float dArea = 1.0/float(TotalRefinement);

      /* loop over fields */
    
      for (field = 0; field < NumberOfBaryonFields; field++) {

      /* Allocate and clear Fluxes */

	ProjectedFluxes.LeftFluxes[field][dim] = new float[size];
	ProjectedFluxes.RightFluxes[field][dim] = new float[size];
	for (i = 0; i < size; i++) {
	  ProjectedFluxes.LeftFluxes[field][dim][i] = 0.0;
	  ProjectedFluxes.RightFluxes[field][dim][i] = 0.0;
	}

	/* if this dim is of length 0, then there is no Flux. */
	
	if (GridDimension[dim] > 1 && MyProcessorNumber == ProcessorNumber) {

	  /* project (downsample by RefinementFactors[i] Fluxes */

	  for (i = 0; i < Dims[0]; i++) {
	    i1 = i/RefinementFactors[0];
	    for (j = 0; j < Dims[1]; j++) {
	      j1 = j/RefinementFactors[1];
	      for (k = 0; k < Dims[2]; k++) {
		k1 = k/RefinementFactors[2];
		*(ProjectedFluxes.LeftFluxes[field][dim] +
		  i1 + j1*ProjectedDims[0] + 
		       k1*ProjectedDims[0]*ProjectedDims[1]) +=
		  (*(BoundaryFluxes.LeftFluxes[field][dim] +
			i + j*Dims[0] + k*Dims[0]*Dims[1])) * dArea;
		*(ProjectedFluxes.RightFluxes[field][dim] +
		  i1 + j1*ProjectedDims[0] + 
		       k1*ProjectedDims[0]*ProjectedDims[1]) +=
		  (*(BoundaryFluxes.RightFluxes[field][dim] +
			i + j*Dims[0] + k*Dims[0]*Dims[1])) * dArea;
	      }
	    }
	  }

	}  // end: if Dims[dim] > 1

      }  // next field

      /* set unused flux pointers to null (makes cleanup easier) */
    
      for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	   field++) {
	ProjectedFluxes.LeftFluxes[field][dim]  = NULL;
	ProjectedFluxes.RightFluxes[field][dim] = NULL;
      }

      if( MHD_Used == TRUE ){

	for(field=0;field<3;field++){
	  ProjectedFluxes.LeftElectric[field][dim]=NULL;
	  ProjectedFluxes.RightElectric[field][dim]=NULL;
	}
	for(field=0;field<3;field++)
	  if( field != dim ){
	    //Electric fields along the face aren't needed.
	    size=1;

	    for(coord=0;coord<3;coord++){
	      Dims[coord] = BoundaryFluxes.LeftFluxEndGlobalIndex[dim][coord] - 
	          BoundaryFluxes.LeftFluxStartGlobalIndex[dim][coord] + 1;

	      ProjectedDims[coord] = Dims[coord]/RefinementFactors[coord];
	      if(dim == coord) ProjectedDims[coord]=1;

	      if( coord != field && coord != dim ){
		Dims[coord]++;
		ProjectedDims[coord]++;
	      }
	      size*= ProjectedDims[coord];

	    }

	    ProjectedFluxes.LeftElectric[field][dim]=new float[size];
	    ProjectedFluxes.RightElectric[field][dim]=new float[size];

	    for(i=0;i<size;i++){
	      ProjectedFluxes.LeftElectric[field][dim][i]=0.0;
	      ProjectedFluxes.RightElectric[field][dim][i]=0.0;
	    }	
	  
	    if( ProcessorNumber == MyProcessorNumber ){
	      //Do the projection.
	      //Note that, since electric fields are co-located along the same edge as their
	      //children, some are skipped. 
	      //( if this ) ? <return this>: <otherwise, this>;
	      
	      for(k=0;k<Dims[2];k+= ((field==2)?1:RefinementFactors[2])){
		k1 = k / RefinementFactors[2];
		for(j=0;j<Dims[1];j+= ((field==1)?1:RefinementFactors[1])){
		  j1=j/RefinementFactors[1];
		  for(i=0;i<Dims[0]; i+= ((field==0)?1:RefinementFactors[0])){
		    i1=i/RefinementFactors[0];
		    
		    index1 = i1+ProjectedDims[0]*(j1+ProjectedDims[1]*k1);
		    index2 = i+Dims[0]*(j+Dims[1]*k);
		    
		    ProjectedFluxes.LeftElectric[field][dim][index1]+=
		      BoundaryFluxes.LeftElectric[field][dim][index2] // /1.0
		      /RefinementFactors[field];
		    
		    ProjectedFluxes.RightElectric[field][dim][index1]+=
		      BoundaryFluxes.RightElectric[field][dim][index2] // /1.0
		      /RefinementFactors[field];
		    
		    /*
		      fprintf(stderr, "gpbf: dim %d field %d ijk, %d %d %d left %f right %f\n",
		      dim, field, i, j, k,
		      BoundaryFluxes.LeftElectric[field][dim][index2],
		      BoundaryFluxes.RightElectric[field][dim][index2]);
		    */
		  }}}//i,j,k
	    }//proc
	  }//E field (field != dim)
      }//mhd correction


    }  // next dimension
  
    /* set unused pointers to NULL */

    for (dim = GridRank; dim < MAX_DIMENSION; dim++)
      for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
	ProjectedFluxes.LeftFluxes[field][dim]  = NULL;
	ProjectedFluxes.RightFluxes[field][dim] = NULL;
      }

    /* If appropriate, receive data and exit. */

    if (ProcessorNumber != MyProcessorNumber) {
      if (CommunicationReceiveFluxes(&ProjectedFluxes, ProcessorNumber,
				   NumberOfBaryonFields, GridRank) == FAIL) {
	fprintf(stderr, "Error in CommunicationReceiveFluxes.\n");
	return FAIL;
      }
      return SUCCESS;
    }

    /* Send fluxes and delete. */

    if (ParentGrid->ProcessorNumber != ProcessorNumber) {
      if (CommunicationSendFluxes(&ProjectedFluxes, 
				  ParentGrid->ProcessorNumber,
				  NumberOfBaryonFields, GridRank) == FAIL) {
	fprintf(stderr, "Error in CommunicationSendFluxes.\n");
	return FAIL;
      }
      DeleteFluxes(&ProjectedFluxes);
    }

  } // end: if (NumberOfBaryonFields > 0)

  return SUCCESS;

}
