/***********************************************************************
/
/  GRID CLASS (PROJECT SOLUTION IN CURRENT GRID TO PARENT GRID
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  NOTE: This routine assumes that the parent and current grids have the
/        same baryon fields.
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

/* function prototypes */

int FindField(int f, int farray[], int n);
extern "C" void FORTRAN_NAME(mult3d)(float *source, float *dest, 
				     int *sdim1, int *sdim2, int *sdim3, 
				     int *ddim1, int *ddim2, int *ddim3,
				     int *sstart1, int *sstart2, int *sstart3, 
				     int *dstart1, int *dstart2, int *dstart3);
extern "C" void FORTRAN_NAME(div3d)(float *source, float *dest, 
				    int *sdim1, int *sdim2, int *sdim3, 
				    int *ddim1, int *ddim2, int *ddim3,
				    int *sstart1, int *sstart2, int *sstart3, 
				    int *dstart1, int *dstart2, int *dstart3,
                                    int *rstart1, int *rstart2, int *rstart3,
                                    int *rend1, int *rend2, int *rend3);


int grid::ProjectSolutionToParentGrid(grid &ParentGrid)
{
  /* Return if this doesn't involve us. */

  if (MyProcessorNumber != ProcessorNumber && 
      MyProcessorNumber != ParentGrid.ProcessorNumber)
    return SUCCESS;

  if (CommunicationDirection == COMMUNICATION_SEND &&
      (MyProcessorNumber == ParentGrid.ProcessorNumber || 
       ProcessorNumber == ParentGrid.ProcessorNumber))
    return SUCCESS;

  if (CommunicationDirection == COMMUNICATION_RECEIVE &&
      MyProcessorNumber != ParentGrid.ProcessorNumber &&
      ProcessorNumber != ParentGrid.ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  this->DebugCheck("ProjectSolutionToParentGrid (before)");

  /* declarations */
    
  int i, j, k, dim, field, One = 1, Zero = 0, skipi, skipj, skipk;
  int ParentStartIndex[MAX_DIMENSION], ParentDim[MAX_DIMENSION],
      ParentEndIndex[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION], Dim[MAX_DIMENSION];

  /* compute size of current grid fields */
    
  int ParentSize = 1, Size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    Size *= GridDimension[dim];
    ParentSize *= ParentGrid.GridDimension[dim];
    ParentDim[dim] = ParentGrid.GridDimension[dim];
    Dim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
  }

  /* set values that are needed for triply-nested loops. */

  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    ParentDim[dim]        = 1;
    ParentStartIndex[dim] = 0;
    ParentEndIndex[dim]   = 0;
    Dim[dim]              = 1;
  }
    
  /* compute refinement factor */
    
  ParentGrid.ComputeRefinementFactors(this, Refinement);

  /* compute the offset (in parent grid index units) from the edge of the
     parent grid to the beginning of the active region in this grid. */

  for (dim = 0; dim < GridRank; dim++) {
    if (GridLeftEdge[dim] >= ParentGrid.GridRightEdge[dim] ||
	GridRightEdge[dim] <= ParentGrid.GridLeftEdge[dim])
      return SUCCESS;
    ParentStartIndex[dim] = nint((GridLeftEdge[dim] - 
				  ParentGrid.GridLeftEdge[dim])/
				 (ParentGrid.CellWidth[dim][0]))
                          + ParentGrid.GridStartIndex[dim];
    ParentEndIndex[dim] = ParentStartIndex[dim] + Dim[dim]/Refinement[dim] - 1;
  }
    
  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  /* Compute the ratio Volume[ThisGridCell]/Volume[ParentCell]. */

  float RelativeVolume = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    RelativeVolume /= float(Refinement[dim]);

  /* Multiply all fields by the density to get conserved quantities. */

  if (ProcessorNumber == MyProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++)
      if (FieldTypeIsDensity(FieldType[field]) == FALSE && (
	  (FieldType[field] < Velocity1 || FieldType[field] > Velocity3)
          || HydroMethod != Zeus_Hydro      )       )
      FORTRAN_NAME(mult3d)(BaryonField[DensNum], BaryonField[field],
			   &Size, &One, &One, &Size, &One, &One,
			   &Zero, &Zero, &Zero, &Zero, &Zero, &Zero);

  /* Allocate Parent temps if it is in another processor. */

  if (ParentGrid.ProcessorNumber != MyProcessorNumber) {
    int ParentSize = 1;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      ParentStartIndex[dim] = 0;
      ParentDim[dim]        = Dim[dim]/Refinement[dim];
      ParentEndIndex[dim]   = ParentDim[dim] - 1;
      ParentSize *= ParentDim[dim];
    }
    for (field = 0; field < NumberOfBaryonFields; field++) {
      delete ParentGrid.BaryonField[field];
      ParentGrid.BaryonField[field] = new float[ParentSize];
    }
  }

  /* For each field, zero the appropriate parental zones. */

  int i1, j1, k1, pindex, gindex;
  if (ProcessorNumber == MyProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++)
      for (k = ParentStartIndex[2]; k <= ParentEndIndex[2]; k++)
	for (j = ParentStartIndex[1]; j <= ParentEndIndex[1]; j++) {
	  pindex = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, pindex++)
	    ParentGrid.BaryonField[field][pindex] = 0.0;
	}

  /* For each field, accumulate it's conserved quantities in the parent 
     grid. */

  if (ProcessorNumber == MyProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++) {
      skipi = skipj = skipk = 1;
      float weight = RelativeVolume;
      if (HydroMethod == Zeus_Hydro) {
	if (FieldType[field] == Velocity1) skipi = Refinement[0];
	if (FieldType[field] == Velocity2) skipj = Refinement[1];
	if (FieldType[field] == Velocity3) skipk = Refinement[2];
      }
      if (skipi*skipj*skipk != 1)
	weight *= float(skipi*skipj*skipk);
      for (k = 0; k < Dim[2]; k += skipk) {
	k1 = k/Refinement[2];
	for (j = 0; j < Dim[1]; j += skipj) {
	  j1 = j/Refinement[1];

	  pindex = (0  + ParentStartIndex[0])                            + 
	           (j1 + ParentStartIndex[1])*ParentDim[0]               +
	           (k1 + ParentStartIndex[2])*ParentDim[0]*ParentDim[1];

	  gindex = 0+GridStartIndex[0]                                      + 
		  (j+GridStartIndex[1])*GridDimension[0]                    +
		  (k+GridStartIndex[2])*GridDimension[0]*GridDimension[1];

	  for (i = 0; i < Dim[0]; i += skipi) { 
	    i1 = i/Refinement[0];
	    ParentGrid.BaryonField[field][pindex+i1] += 
	      BaryonField[field][gindex+i]*weight;
	  }
	}
      }
    }

  /* If necessary, copy the projected field from the 'fake' ParentGrid to
     the real one. */

  if (ProcessorNumber != ParentGrid.ProcessorNumber) {
    int ParentRegionDim[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ParentRegionDim[dim] = ParentEndIndex[dim] - ParentStartIndex[dim] + 1;
    ParentGrid.CommunicationReceiveRegion(&ParentGrid, ProcessorNumber,
	      ALL_FIELDS, NEW_ONLY, ParentStartIndex, ParentRegionDim, TRUE);
  }

  /* Divide all fields by mass to return to original quantity. */
    
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (FieldTypeIsDensity(FieldType[field]) == FALSE && (
	  (FieldType[field] < Velocity1 || FieldType[field] > Velocity3)
          || HydroMethod != Zeus_Hydro      )       ) {
      if (ProcessorNumber == MyProcessorNumber)
	FORTRAN_NAME(div3d)(BaryonField[DensNum], BaryonField[field],
			    &Size, &One, &One, &Size, &One, &One,
			    &Zero, &Zero, &Zero, &Zero, &Zero, &Zero,
			    &Zero, &Zero, &Zero, &Size, &Zero, &Zero);
      if (ParentGrid.ProcessorNumber == MyProcessorNumber)
	FORTRAN_NAME(div3d)(ParentGrid.BaryonField[DensNum],
			    ParentGrid.BaryonField[field],
			    ParentGrid.GridDimension, 
			      ParentGrid.GridDimension+1,
			      ParentGrid.GridDimension+2,
			    ParentGrid.GridDimension, 
			      ParentGrid.GridDimension+1,
			      ParentGrid.GridDimension+2,
			    &Zero, &Zero, &Zero, &Zero, &Zero, &Zero,
			    ParentStartIndex, ParentStartIndex+1, 
                             ParentStartIndex+2,
			    ParentEndIndex, ParentEndIndex+1, ParentEndIndex+2);
			  
    }
    
  /* If appropriate, restore consistency between total and internal
     energy in projected regions. */

  if (ParentGrid.ProcessorNumber == MyProcessorNumber)
   if (DualEnergyFormalism)
    for (k = ParentStartIndex[2]; k <= ParentEndIndex[2]; k++)
      for (j = ParentStartIndex[1]; j <= ParentEndIndex[1]; j++) {

	i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];

	for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++)
	  ParentGrid.BaryonField[TENum][i1] = 
	    ParentGrid.BaryonField[GENum][i1] + 0.5*
	    ParentGrid.BaryonField[Vel1Num][i1] * 
	    ParentGrid.BaryonField[Vel1Num][i1];

	i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];

	if (GridRank > 1)
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++)
	    ParentGrid.BaryonField[TENum][i1] += 0.5*
	      ParentGrid.BaryonField[Vel2Num][i1] *
	      ParentGrid.BaryonField[Vel2Num][i1];

	i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];

	if (GridRank > 2)
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++)
	    ParentGrid.BaryonField[TENum][i1] += 0.5*
	      ParentGrid.BaryonField[Vel3Num][i1] *
	      ParentGrid.BaryonField[Vel3Num][i1];
		  
      } // end: loop over faces

  /* Clean up the fake ParentGrid. */

  if (ParentGrid.ProcessorNumber != MyProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++) {
      delete ParentGrid.BaryonField[field];
      ParentGrid.BaryonField[field] = NULL;
    }

  ParentGrid.DebugCheck("ProjectSolutionToParentGrid (Parent, after)");

  return SUCCESS;
}
