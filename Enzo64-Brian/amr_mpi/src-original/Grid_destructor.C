/***********************************************************************
/
/  GRID CLASS (DESTRUCTOR)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
//
//  Grid destructor
//
#include <stdio.h>
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
void DeleteFluxes(fluxes *Fluxes);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
 
grid::~grid()
{
 
  int i;
 
  /* Error check. */
 
#ifdef UNUSED
  if (NumberOfParticles > 0) {
    fprintf(stderr, "warning: destroying live particles (%"ISYM").\n",
	    NumberOfParticles);
  /* exit(EXIT_FAILURE); */
  }
#endif /* UNUSED */
 
  for (i = 0; i < MAX_DIMENSION; i++) {
    delete CellLeftEdge[i];
    delete CellWidth[i];
    delete ParticlePosition[i];
    delete ParticleVelocity[i];
    delete ParticleAcceleration[i];
    delete AccelerationField[i];
    delete RandomForcingField[i];
  }
 
  delete ParticleAcceleration[MAX_DIMENSION];
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    delete BaryonField[i];
    delete OldBaryonField[i];
  }

  for (i = 0; i < MAX_DIMENSION; i++) {
    if(OldAccelerationField[i] != NULL ){
      delete OldAccelerationField[i];
      OldAccelerationField[i] = NULL;
    }
  }
 
  DeleteFluxes(BoundaryFluxes);
  delete BoundaryFluxes;
 
  delete ParticleMass;
  delete ParticleNumber;
  delete ParticleType;
  delete PotentialField;
  delete GravitatingMassField;
  delete GravitatingMassFieldParticles;
  delete FlaggingField;
  delete MassFlaggingField;
 
  for (i = 0; i < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; i++)
    delete ParticleAttribute[i];

/* 
  if (debug && GridRank > 0) {
    printf("grid->destructor: deleting grid with dims = ");
    WriteListOfInts(stdout, GridRank, GridDimension);
  }
*/
 
}
