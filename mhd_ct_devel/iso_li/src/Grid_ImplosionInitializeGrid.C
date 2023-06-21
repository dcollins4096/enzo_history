/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR IMPLOSION PROBLEM)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, December 2004.
/
/  PURPOSE: Sets density and total energy in the lower left corner of 
/           the domain.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::ImplosionInitializeGrid(float ImplosionDiamondDensity,
				  float ImplosionDiamondTotalEnergy,
				  float ImplosionDensity,
				  float ImplosionTotalEnergy)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* set fields in a "diamond" region: y=0.15-x. */

#define OLDCOLD

#ifdef OLDCOLD

  printf ("DEBUG: using OLDCODE (delete when finished!)\n");
  int index, jndex, i;
  for (i = 0; i < size; i++) {
    index = i % GridDimension[0];
    jndex = (i-index)/GridDimension[0];
    if (*(CellLeftEdge[0] + index) + 0.5*(*(CellWidth[0] + index))
	     + *(CellLeftEdge[1] + jndex) + 0.5*(*(CellWidth[1] + jndex))
	< 0.1517) { // must be 0.15 but this creates an imperfect initial front
                   // 0.151 is good for 400^2 simulation to fix the 
      // imperfection. 0.1517 is good for L2x2 start-up to match the unigrid.
      BaryonField[0][i] = ImplosionDiamondDensity;
      BaryonField[1][i] = ImplosionDiamondTotalEnergy;
    }
  }
#else /* NEWCODE */

  printf ("DEBUG: using NEWCODE (delete OLDCODE when finished!)\n");

  int index, jndex, i;

  for (jndex = 0; jndex < GridDimension[1]; jndex++) {
    for (index = 0; index < GridDimension[0]; index++) {
      i = index + jndex*GridDimension[0];
      // ll = lower-left point; ur = upper-right point
      float ll = CellLeftEdge[0][index] + CellLeftEdge[1][jndex];
      float ur = ll + CellWidth[0][index] + CellWidth[1][jndex];
      if (ur < 0.15) {
        BaryonField[0][i] = ImplosionDiamondDensity;
        BaryonField[1][i] = ImplosionDiamondTotalEnergy;
      } else if (ll > 0.15) {
        BaryonField[0][i] = ImplosionDensity;
        BaryonField[1][i] = ImplosionTotalEnergy;
      } else {
        // (Rough) estimate of volume fraction overlap
        float infrac  = (0.15-ll)/(ur-ll);
        float outfrac = 1.0 - infrac;
        BaryonField[0][i] =  infrac*ImplosionDiamondDensity
          +                 outfrac*ImplosionDensity;
        BaryonField[1][i] =  infrac*ImplosionDiamondTotalEnergy
          +                 outfrac*ImplosionTotalEnergy;
      }
      printf ("DEBUG %d %d %g %g\n",index,jndex,BaryonField[0][i],BaryonField[1][i]);
    }
  }

#endif

  return SUCCESS;
}
