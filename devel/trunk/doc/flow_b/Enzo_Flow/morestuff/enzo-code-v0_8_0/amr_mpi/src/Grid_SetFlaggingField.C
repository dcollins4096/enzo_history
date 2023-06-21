/***********************************************************************
/
/  GRID CLASS (CLEAR THE FLAGGING FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* The following is defined in Grid_DepositParticlePositions.C. */

extern float DepositParticleMaximumParticleMass;


int grid::SetFlaggingField(int &NumberOfFlaggedCells, int level)
{

  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  NumberOfFlaggedCells = INT_UNDEFINED;
  int method;

  for (method = 0; method < MAX_FLAGGING_METHODS; method++) {

  /***********************************************************************/
  /* beginning of Cell flagging criterion routine                        */

    switch (CellFlaggingMethod[method]) {

    case 0:   /* no action */
      NumberOfFlaggedCells = (NumberOfFlaggedCells == INT_UNDEFINED ?
			      0 : NumberOfFlaggedCells);
      break;

      /* ==== METHOD 1: BY SLOPE ==== */
  
    case 1:
    
      /* flag all points needing extra resolution (FlagCellsToBeRefinedBySlop
	 returns the number of flagged cells). */
    
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedBySlope();
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedBySlope.\n");
	return FAIL;
      }
      break;

      /* ==== METHOD 2: BY BARYON MASS OR OVERDENSITY ==== */

    case 2:
    
      /* allocate and clear mass flagging field */
    
      this->ClearMassFlaggingField();
    
      /* baryons: add baryon density to mass flagging field (so the mass 
	 flagging field contains the mass in the cell (not the density) */
    
      if (this->AddFieldMassToMassFlaggingField() == FAIL) {
	fprintf(stderr, "Error in grid->AddFieldMassToMassFlaggingField.\n");
	return FAIL;
      }
    
      /* flag all points that need extra resolution (FlagCellsToBeRefinedByMass
	 return the number of flagged cells). */
    
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMass(level, method);
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByMass.\n");
	return FAIL;
      }
      break;

      /* ==== METHOD 3: BY SHOCKS ==== */

    case 3:

      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByShocks();
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByShocks.\n");
	return FAIL;
      }
      break;

      /* ==== METHOD 4: BY PARTICLE MASS ==== */

    case 4:
    
      /* Allocate and clear mass flagging field. */
    
      this->ClearMassFlaggingField();

      /* Set the maximum particle mass to be deposited (cleared below). */

      DepositParticleMaximumParticleMass = 
	0.99999*MinimumMassForRefinement[method]*POW(RefineBy, 
		    level*MinimumMassForRefinementLevelExponent[method]);

      /* Deposit particles in this grid to MassFlaggingField. */

      if (this->DepositParticlePositions(this, this->ReturnTime(),
		  		         MASS_FLAGGING_FIELD) == FAIL) {
	fprintf(stderr, "Error in grid->DepositParticlePositions.\n");
	return FAIL;
      }

      DepositParticleMaximumParticleMass = 0;
    
      /* Flag all points that need extra resolution (FlagCellsToBeRefinedByMass
	 return the number of flagged cells). */
    
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMass(level, method);
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByMass.\n");
	return FAIL;
      }
      break;

      /* ==== METHOD 6: BY JEANS LENGTH ==== */

    case 6:
    
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByJeansLength();
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByJeansLength.\n");
	return FAIL;
      }
      break;

      /* ==== METHOD 7: BY COOLING TIME < DX/SOUND SPEED ==== */

    case 7:
    
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByCoolingTime();
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByCoolingTime.\n");
	return FAIL;
      }
      break;

      /* ==== undefined ==== */

    case INT_UNDEFINED:
      break;

    default:
      fprintf(stderr, "CellFlaggingMethod[%d] = %d unknown\n", method,
	      CellFlaggingMethod[method]);
      return FAIL;

    }
  
    /* End of Cell flagging criterion routine                              */
    /***********************************************************************/

  } // end: loop over refinement methods

  if (NumberOfFlaggedCells == INT_UNDEFINED) {
    fprintf(stderr, "No valid CellFlaggingMethod specified.\n");
    return FAIL;
  }

#ifdef USE_MPI
#ifdef MPI_INSTRUMENTATION
  /* Zhiling Lan's instrumented part */
  counter[4]++;
  ZLAN_COUNT(4,NumberOfFlaggedCells);
#endif /* MPI_INSTRUMENTATION */
#endif

  if (debug)
    printf("SetFlaggingField: NumberOfFlaggedCells = %d.\n", 
	   NumberOfFlaggedCells);

  return SUCCESS;

}
