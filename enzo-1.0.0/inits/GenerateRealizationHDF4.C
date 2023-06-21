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
/  GENERATES THE FIELD AND PARTICLE REALIZATIONS
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_HDF4

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "Parameters.h"

/* function prototypes */

extern "C" void FORTRAN_NAME(set_common)(float *lam0_in, float *omega0_in, 
					 float *zri_in, float *hub_in);
extern "C" float FORTRAN_NAME(calc_f)(float *aye);
extern "C" float FORTRAN_NAME(calc_ayed)(float *aye);
int GenerateField(int Rank, int Dims[3], int MaxDims[3], int WaveNumberCutoff,
		  float *Field, int FieldType, int NewCenter[3], 
		  int Refinement, int StartIndex[3], int Species);
int WriteFieldHDF4(int Rank, int Dims[3], float *Field, char *Name, int Type);
void RemoveSubGridParticles(float *From, float *To, int Dims[],
			   int Start[], int End[]);



int GenerateRealization(parmstruct *Parameters, parmstruct *SubGridParameters)
{

  float *ParticleField, *GridField, LeftEdge[3], RightEdge[3], ParticleOffset;
  float GrowthFunction, aye = 1.0, ayed, Temp;
  int i, j, k, dim, size, index, NumberOfParticles, MaxDimsThisLevel[3];

  /* Calculate some cosmological quantities for later use. */

  FORTRAN_NAME(set_common)(&OmegaLambdaNow, &OmegaMatterNow, &InitialRedshift,
			   &HubbleConstantNow);
  GrowthFunction = FORTRAN_NAME(calc_f)(&aye);  /* dlog(D)/dlog(a) */
  ayed           = FORTRAN_NAME(calc_ayed)(&aye);
  if (debug) printf("GrowthFunction: dlog(D+)/dlog(a) = %g\n", GrowthFunction);

  /* Compute the position of the left corner of the particle grid. */

  for (dim = 0; dim < Parameters->Rank; dim++) {
    if (Parameters->NewCenter[0] != INT_UNDEFINED && Parameters->MaxDims[dim]
	!= Parameters->ParticleDims[dim]*Parameters->ParticleRefinement)
      LeftEdge[dim] = 
	float((Parameters->StartIndex[dim] + Parameters->MaxDims[dim]/2-1 -
	       Parameters->NewCenter[dim] + Parameters->MaxDims[dim]      ) 
	      % Parameters->MaxDims[dim]) /
	float(Parameters->MaxDims[dim]);
    else
      LeftEdge[dim] = float(Parameters->StartIndex[dim]) / 
	              float(Parameters->MaxDims[dim]);
    RightEdge[dim] = LeftEdge[dim] + 
      float(Parameters->ParticleDims[dim]*Parameters->ParticleRefinement)/
      float(Parameters->MaxDims[dim]);
  }

  if (debug) {
    printf("ParticleSubgridLeftEdge  = %f %f %f\n", LeftEdge[0], LeftEdge[1], 
	   LeftEdge[2]);
    printf("ParticleSubgridRightEdge = %f %f %f\n", RightEdge[0], 
	   RightEdge[1], RightEdge[2]);
  }

  /* ------------------------------------------------------------------- */
  /* Set particles. */

  if (Parameters->InitializeParticles) {

    /* Compute required size and allocate field. */

    size = 1, NumberOfParticles = 1;
    for (dim = 0; dim < Parameters->Rank; dim++) {
      size *= (Parameters->ParticleDims[dim] + ((dim == 0) ? 2 : 0));
      NumberOfParticles *= Parameters->ParticleDims[dim];
    }
    if ((ParticleField = new float[size]) == NULL) {
      fprintf(stderr, "GenerateRealization: malloc failure (%d).\n", size);
      exit(EXIT_FAILURE);
    }

    /* Find out where to remove particles if SubGridParameters is set. */

    int SubStart[3], SubEnd[3], ParticlesRemoved = 1, Shift[3];
    if (SubGridParameters) {
      for (dim = 0; dim < Parameters->Rank; dim++) {
	if (Parameters->MaxDims[dim] != SubGridParameters->MaxDims[dim]) {
	  fprintf(stderr, "SubGridParameter MaxDims must match.\n");
	  exit(EXIT_FAILURE);
	}
	if (SubGridParameters->StartIndex[dim] % Parameters->ParticleRefinement
	    != 0) {
	  fprintf(stderr, "SubGrid StartIndex must be divisible by refinement.\n");
	  exit(EXIT_FAILURE);
	}

	/* If this is the top grid, the corner is shifted if recentering. */

	if (Parameters->MaxDims[dim] == 
	    Parameters->ParticleDims[dim]*Parameters->ParticleRefinement &&
	    Parameters->NewCenter[dim] != INT_UNDEFINED)
	  Shift[dim] = Parameters->NewCenter[dim] - 
	    (Parameters->MaxDims[dim]/2 - 1);
	else
	  Shift[dim] = 0;
	printf("Shift[%d] = %d\n", dim, Shift[dim]);

	/* Compute start and end indices subgrid region. */

	SubStart[dim] = ((SubGridParameters->StartIndex[dim]-
			         Parameters->StartIndex[dim] - Shift[dim]) 
			 % Parameters->MaxDims[dim])/
	  Parameters->ParticleRefinement;
	MaxDimsThisLevel[dim] = Parameters->MaxDims[dim]/
	                        Parameters->ParticleRefinement;
	SubStart[dim] = (SubStart[dim] + MaxDimsThisLevel[dim]) %
	                MaxDimsThisLevel[dim];

	SubEnd[dim] = ((SubGridParameters->StartIndex[dim] -
		               Parameters->StartIndex[dim] + 
		        SubGridParameters->ParticleDims[dim]*
		        SubGridParameters->ParticleRefinement - Shift[dim])
		       % Parameters->MaxDims[dim])/
	  Parameters->ParticleRefinement - 1;
	SubEnd[dim] = (SubEnd[dim] + MaxDimsThisLevel[dim]) %
	              MaxDimsThisLevel[dim];

	ParticlesRemoved *= SubEnd[dim] - SubStart[dim] + 1;
      } // end: loop over dims
      NumberOfParticles -= ParticlesRemoved;
      if (debug)
	printf("Removing Particle Region = %d %d %d -> %d %d %d\n",
	       SubStart[0], SubStart[1], SubStart[2],
	       SubEnd[0], SubEnd[1], SubEnd[2]);
    }
    if (debug) printf("NumberOfParticles = %d\n", NumberOfParticles);

      /* Loop over dimensions. */

    for (dim = 0; dim < Parameters->Rank; dim++) {
      if (debug) printf("GenerateRealization: particle dim %d.\n", dim);

      /* 1) velocities
	    -generate the displacement field (f delta_k -i vec(k)/k^2). 
	    -multiply by adot to get velocity */
	
      GenerateField(Parameters->Rank, Parameters->ParticleDims,
		    Parameters->MaxDims, Parameters->WaveNumberCutoff, 
		    ParticleField, 1+dim, Parameters->NewCenter, 
		    Parameters->ParticleRefinement, Parameters->StartIndex, 1);

      Temp = ayed * GrowthFunction;
      for (i = 0; i < size; i++)
	ParticleField[i] *= Temp;

      /* Remove subgrid particles if necessary (using temp field),
         and write out field. */

      if (SubGridParameters) {
	float *TempField = new float[size];
	RemoveSubGridParticles(ParticleField, TempField, 
			       Parameters->ParticleDims, SubStart, SubEnd);
	WriteFieldHDF4(1, &NumberOfParticles,
		   TempField, Parameters->ParticleVelocityName, dim);
	delete TempField;
      } 
      else 
	WriteFieldHDF4(1, &NumberOfParticles,
		   ParticleField, Parameters->ParticleVelocityName, dim);

	
      /* 2) Make the position field by converting velocity to displacement
	    and adding the initial position. */

      Temp = 1.0/ayed;
      for (i = 0; i < size; i++)
	ParticleField[i] *= Temp;

//      ParticleOffset = (Parameters->ParticleRefinement == 
//			Parameters->GridRefinement       ) ? 0.5 : 0.0;
        ParticleOffset = 0.5;
#ifdef SHIFT_FOR_LARS
        ParticleOffset = 1.0;
#endif /* SHIFT_FOR_LARS */
      if (debug) printf("ParticleOffset = %g\n", ParticleOffset);

      for (k = 0; k < Parameters->ParticleDims[2]; k++)
	for (j = 0; j < Parameters->ParticleDims[1]; j++) {
	  index = (k*Parameters->ParticleDims[1] + j)*
	          Parameters->ParticleDims[0];
	  Temp = float(Parameters->ParticleRefinement) / 
	         float(Parameters->MaxDims[dim]);
	  if (dim == 0)
	    for (i = 0; i < Parameters->ParticleDims[0]; i++)
	      ParticleField[index+i] += LeftEdge[dim] + 
		                 (float(i)+ParticleOffset)*Temp;
	  if (dim == 1)
	    for (i = 0; i < Parameters->ParticleDims[0]; i++)
	      ParticleField[index+i] += LeftEdge[dim] +
		                 (float(j)+ParticleOffset)*Temp;
	  if (dim == 2)
	    for (i = 0; i < Parameters->ParticleDims[0]; i++)
	      ParticleField[index+i] += LeftEdge[dim] +
		                 (float(k)+ParticleOffset)*Temp;
	  for (i = 0; i < Parameters->ParticleDims[0]; i++) {
	    if (ParticleField[index+i] <  0.0) ParticleField[index+i] += 1.0;
	    if (ParticleField[index+i] >= 1.0) ParticleField[index+i] -= 1.0;
	  }
	}

      /* Remove subgrid particles if necessary (using temp field),
         and write out field. */

      if (SubGridParameters) {
	float *TempField = new float[size];
	RemoveSubGridParticles(ParticleField, TempField, 
			       Parameters->ParticleDims, SubStart, SubEnd);
	WriteFieldHDF4(1, &NumberOfParticles,
		   TempField, Parameters->ParticlePositionName, dim);
	delete TempField;
      } 
      else 
	WriteFieldHDF4(1, &NumberOfParticles,
		   ParticleField, Parameters->ParticlePositionName, dim);

    } /* end: loop over dims */

    /* Generate and write out particle mass field. */

    if (Parameters->ParticleMassName != NULL) {
      float ParticleMass = (OmegaMatterNow-OmegaBaryonNow)/OmegaMatterNow;
      for (dim = 0; dim < Parameters->Rank; dim++)
	ParticleMass *= Parameters->ParticleRefinement/
	  Parameters->GridRefinement;
      for (i = 0; i < NumberOfParticles; i++)
	ParticleField[i] = ParticleMass;
      WriteFieldHDF4(1, &NumberOfParticles, ParticleField, 
		 Parameters->ParticleMassName, 1);
    }

    /* Free field. */

    delete ParticleField;

  } /* end: if (InitializeParticles) */

  /* ------------------------------------------------------------------- */
  /* Set grids. */

  if (Parameters->InitializeGrids) {

    /* Compute required size and allocate field. */

    size = 1;
    for (dim = 0; dim < Parameters->Rank; dim++)
      size *= (Parameters->GridDims[dim] + ((dim == 0) ? 2 : 0));
    if ((GridField = new float[size]) == NULL) {
      fprintf(stderr, "GenerateRealization: malloc failure (%d).\n", size);
      exit(EXIT_FAILURE);
    }

    /* 1) density (add one and multiply by mean density). */

    if (debug) printf("Generating grid densities.\n");
    GenerateField(Parameters->Rank, Parameters->GridDims,
		  Parameters->MaxDims, Parameters->WaveNumberCutoff, 
		  GridField, 0, Parameters->NewCenter, 
		  Parameters->GridRefinement, Parameters->StartIndex, 2);

    Temp = OmegaBaryonNow/OmegaMatterNow;
    for (i = 0; i < size; i++)
      GridField[i] = max(GridField[i] + 1.0, 0.1) * Temp;

    WriteFieldHDF4(Parameters->Rank, Parameters->GridDims,
	       GridField, Parameters->GridDensityName, 0);

    /* 2) velocities. */

    for (dim = 0; dim < Parameters->Rank; dim++) {
      if (debug) printf("GenerateRealization: grid velocity dim %d.\n", dim);

      /* 1) velocities
	    -generate the displacement field (f delta_k -i vec(k)/k^2). 
	    -multiply adot to get velocity */
	
      GenerateField(Parameters->Rank, Parameters->GridDims,
		    Parameters->MaxDims, Parameters->WaveNumberCutoff, 
		    GridField, 1+dim, Parameters->NewCenter, 
		    Parameters->GridRefinement, Parameters->StartIndex, 2);

      Temp = ayed * GrowthFunction;
      for (i = 0; i < size; i++)
	GridField[i] *= Temp;

      WriteFieldHDF4(Parameters->Rank, Parameters->GridDims,
		 GridField, Parameters->GridVelocityName, dim);
    }

  }

  return SUCCESS;
}


void RemoveSubGridParticles(float *From, float *To, int Dims[],
			   int Start[], int End[])
{
  int i, j, k, findex = 0, tindex = 0;

  for (k = 0; k < Dims[2]; k++)
    for (j = 0; j < Dims[1]; j++)
      for (i = 0; i < Dims[0]; i++, findex++)
	if (k < Start[2] || k > End[2] ||
	    j < Start[1] || j > End[1] || 
	    i < Start[0] || i > End[0])
	  To[tindex++] = From[findex];

}
			   
#endif /* USE_HDF4 */
