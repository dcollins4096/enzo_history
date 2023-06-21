/***********************************************************************
/
/  GRID CLASS (CREATES TRACER PARTICLES IN GRID ACCORDING TO PARAMETERS)
/
/  written by: Greg Bryan
/  date:       April, 2004
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
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
 
 
/* function prototypes */
 
 
 
int grid::TracerParticleCreateParticles(FLOAT LeftEdge[], FLOAT RightEdge[],
					FLOAT Spacing, int &TotalParticleCount)
{
 
  /* Compute the number of new particles to create. */
 
  int i, j, k, dim, NumberOfTracerParticles = 1, Dims[] = {1,1,1};
  int count;
 
  for (dim = 0; dim < GridRank; dim++) {
    Dims[dim] = nint((RightEdge[dim] - LeftEdge[dim])/Spacing);
    NumberOfTracerParticles *= Dims[dim];
  }
 
/*
  if (ProcessorNumber == MyProcessorNumber) {
 
  fprintf(stderr, "PX %"ISYM" NOP %"ISYM" TPC %"ISYM" NOTP %"ISYM"\n", MyProcessorNumber,
    NumberOfParticles, TotalParticleCount, NumberOfTracerParticles);
  fprintf(stderr, "GL %"ISYM"  %12.3e  %12.3e  %12.3e\n", MyProcessorNumber,
    GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
  fprintf(stderr, "GR %"ISYM"  %12.3e  %12.3e  %12.3e\n", MyProcessorNumber,
    GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
 
  }
*/
 
  /* Exit if grid is not on this processor, but not before recording increase
     in the number of particles, both in this grid and in the metadata. */
 
  if (ParallelRootGridIO == FALSE) {
  if (MyProcessorNumber != ProcessorNumber) {
    NumberOfParticles += NumberOfTracerParticles;
    TotalParticleCount += NumberOfTracerParticles;
    return SUCCESS;
  }
  }
 
//  if (debug)
//    printf("Total number of new tracer particles = %"ISYM"\n", NumberOfTracerParticles);
 
  if (NumberOfTracerParticles > 1000000) {
    fprintf(stderr, "Too many new particles!\n");
    return FAIL;
  }
 
  if (ProcessorNumber == MyProcessorNumber) {
 
  /* Record pointers to current particles. */
 
  FLOAT *TempPos[MAX_DIMENSION];
  float *TempVel[MAX_DIMENSION], *TempMass,
        *TempAttribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  int *TempNumber, *TempType;
 
  TempMass = ParticleMass;
  TempNumber = ParticleNumber;
  TempType = ParticleType;
 
  for (dim = 0; dim < GridRank; dim++) {
    TempPos[dim] = ParticlePosition[dim];
    TempVel[dim] = ParticleVelocity[dim];
  }
  for (j = 0; j < NumberOfParticleAttributes; j++)
    TempAttribute[j] = ParticleAttribute[j];
 
  /* Allocate enough space for new particles. */
 
  this->AllocateNewParticles(NumberOfParticles + NumberOfTracerParticles);
 
  /* Copy and delete old particles. */
 
  for (i = 0; i < NumberOfParticles; i++) {
    ParticleMass[i]   = TempMass[i];
    ParticleNumber[i] = TempNumber[i];
    ParticleType[i]   = TempType[i];
  }
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < NumberOfParticles; i++) {
      ParticlePosition[dim][i] = TempPos[dim][i];
      ParticleVelocity[dim][i] = TempVel[dim][i];
    }
  for (j = 0; j < NumberOfParticleAttributes; j++)
    for (i = 0; i < NumberOfParticles; i++)
      ParticleAttribute[j][i] = TempAttribute[j][i];
	
  delete [] TempMass;
  delete [] TempNumber;
  delete [] TempType;
  for (dim = 0; dim < GridRank; dim++) {
    delete [] TempPos[dim];
    delete [] TempVel[dim];
  }
  for (j = 0; j < NumberOfParticleAttributes; j++)
    delete [] TempAttribute[j];
 
  /* Fill out tracer particle information. */
 
  int index = NumberOfParticles;
 
//  int count = 0;
  count = 0;
 
  FLOAT pos[MAX_DIMENSION];
 
  for (k = 0; k < Dims[2]; k++) {
 
    if (GridRank > 2)
      pos[2] = LeftEdge[2] + Spacing*(0.5+k);
 
    if ((GridLeftEdge[2] < pos[2]) && (pos[2] < GridRightEdge[2])) {
 
    for (j = 0; j < Dims[1]; j++) {
 
      if (GridRank > 1)
	pos[1] = LeftEdge[1] + Spacing*(0.5+j);
 
      if ((GridLeftEdge[1] < pos[1]) && (pos[1] < GridRightEdge[1])) {
 
      for (i = 0; i < Dims[0]; i++) {
 
	pos[0] = LeftEdge[0] + Spacing*(0.5+i);
 
        if ((GridLeftEdge[0] < pos[0]) && (pos[0] < GridRightEdge[0])) {
 
        count++;
 
	/* Set particle position. */
 
	for (dim = 0; dim < GridRank; dim++)
	  ParticlePosition[dim][index] = pos[dim];
 
	/* Set particle mass and type. */
 
	ParticleMass[index] = 0;
	ParticleType[index] = PARTICLE_TYPE_TRACER;
 
	/* Set particle index.  Note that there is currently a problem
	   with star particles, which may overlap this range of indexes
	   if star particles have already been created at this point. */
 
//	ParticleNumber[index] = TotalParticleCount +
//	                        (index-NumberOfParticles);
 
        ParticleNumber[index] = index;
 
	/* Set Particle attributes to FLOAT_UNDEFINED. */
 
	for (int atr = 0; atr < NumberOfParticleAttributes; atr++)
	  ParticleAttribute[atr][index] = FLOAT_UNDEFINED;
 
	/* update particle count index */
 
	index++;
 
        }
      } // end: loop over i
      }
    } // end: loop over j
    }
  } // end: loop over k
 
  fprintf(stderr, "PG %"ISYM" actual count %"ISYM"\n", MyProcessorNumber, count);
 
  /* Set particle velocities. */
 
//  this->TracerParticleSetVelocity();
 
  }
 
  /* Update number of particles. */
 
  if (ParallelRootGridIO == FALSE) {
  TotalParticleCount += NumberOfTracerParticles;
  NumberOfParticles += NumberOfTracerParticles;
  this->TracerParticleSetVelocity();
  }
 
  if (ParallelRootGridIO == TRUE) {
  if (MyProcessorNumber == ProcessorNumber) {
//  TotalParticleCount += count;
  NumberOfParticles += count;
  this->TracerParticleSetVelocity();
  }
  }
 
  return SUCCESS;
}
