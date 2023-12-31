/***********************************************************************
/
/  GRID CLASS (SOLVE POISSON ON POTENTIAL FIELD)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int MultigridSolver(float *RHS, float *Solution, int Rank, int TopDims[],
		    float &norm, float &mean, int start_depth,
		    float tolerance, int max_iter);
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
extern "C" void FORTRAN_NAME(smooth2)(float *source, float *dest, int *ndim,
                                   int *sdim1, int *sdim2, int *sdim3);
 
#define TOLERANCE 5.0e-4
#define MAX_ITERATION 10
 
int grid::SolveForPotential(int level, FLOAT PotentialTime)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  JBPERF_START("grid_SolveForPotential");
 
  /* declarations */
 
  int dim, size = 1, i;
  float tol_dim = TOLERANCE * POW(0.1, 3-GridRank);
  //  if (GridRank == 3)
  //    tol_dim = 1.0e-5;
 
  /* Compute adot/a at time = t+1/2dt (time-centered). */
 
  if (PotentialTime < 0)
    PotentialTime = Time + 0.5*dtFixed;
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(PotentialTime, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      ENZO_FAIL("");
    }
 
  /* Compute right hand side. */
 
  float InverseVolumeElement = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GravitatingMassFieldDimension[dim];
    InverseVolumeElement *= (GravitatingMassFieldDimension[dim]-1);
  }
  tol_dim = max(sqrt(float(size))*1e-6, tol_dim);
 
  float *rhs = new float[size];
 
  float Constant = GravitationalConstant * InverseVolumeElement *
                   POW(GravitatingMassFieldCellSize, 2) / a;
 
#define NO_SMOOTH_SOURCE
#ifdef SMOOTH_SOURCE
 
  FORTRAN_NAME(smooth2)(GravitatingMassField, rhs, &GridRank,
			GravitatingMassFieldDimension,
			GravitatingMassFieldDimension+1,
			GravitatingMassFieldDimension+2);
#if 0
  FORTRAN_NAME(smooth2)(rhs, GravitatingMassField, &GridRank,
			GravitatingMassFieldDimension,
			GravitatingMassFieldDimension+1,
			GravitatingMassFieldDimension+2);
  FORTRAN_NAME(smooth2)(GravitatingMassField, rhs, &GridRank,
			GravitatingMassFieldDimension,
			GravitatingMassFieldDimension+1,
			GravitatingMassFieldDimension+2);
#endif
  for (i = 0; i < size; i++)
    rhs[i] *= Constant;
 
#else /* SMOOTH_SOURCE */
 
  for (i = 0; i < size; i++)
    rhs[i] = GravitatingMassField[i] * Constant;
 
#endif /* SMOOTH_SOURCE */
 
  /* Restrict fields to lower resolution if desired. */
 
  int GravitySmooth = max(level - MaximumGravityRefinementLevel, 0);
  GravitySmooth = 0;
 
  /* Iterate with multigrid. */
 
  float norm = huge_number, mean = norm;
#ifdef UNUSED
  int iteration = 0;
#endif /* UNUSED */
 
  if (MultigridSolver(rhs, PotentialField, GridRank,
		      GravitatingMassFieldDimension, norm, mean,
		      GravitySmooth, tol_dim, MAX_ITERATION) == FAIL) {
    fprintf(stderr, "Error in MultigridDriver.\n");
    ENZO_FAIL("");
  }
 
#ifdef UNUSED
  while (norm/mean > tol_dim) {
    if (MultigridSolver(rhs, PotentialField, GridRank,
			GravitatingMassFieldDimension, norm, mean,
			GravitySmooth) == FAIL) {
      fprintf(stderr, "Error in MultigridDriver.\n");
      ENZO_FAIL("");
    }
    printf("%"ISYM" %"GSYM"\n", iteration, norm/mean);
    if (iteration++ > MAX_ITERATION) {
      fprintf(stderr, "exceeding iteration count (%"ISYM")\n", iteration);
      ENZO_FAIL("");
    }
  }
#endif /* UNUSED */
 
  /* Clean up. */
 
  delete [] rhs;

#ifdef POTENTIALDEBUGOUTPUT
  for (int i=0;i<GridDimension[0]; i++) {
    int igrid = GRIDINDEX_NOGHOST(i,(GridEndIndex[0]+GridStartIndex[0])/2,(GridEndIndex[0]+GridStartIndex[0])/2);
    printf("i: %i \t SolvedSub %g\n", i, PotentialField[igrid]);
  }
  float maxPot=-1e30, minPot=1e30;    
  for (int i=0;i<size; i++) {
    maxPot = max(maxPot,PotentialField[i]);
    minPot = min(minPot,PotentialField[i]);
  }
  if (debug1) printf("SolvedPotential: Potential minimum: %g \t maximum: %g\n", minPot, maxPot);
#endif
 
  JBPERF_STOP("grid_SolveForPotential");
  return SUCCESS;
}
 
