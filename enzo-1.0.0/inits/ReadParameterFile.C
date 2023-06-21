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
/  READ A PARAMETER FILE
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine reads the parameter file in the argument and sets parameters
//   based on it.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "global_data.h"
#include "Parameters.h"


int ReadParameterFile(FILE *fptr, parmstruct *Parameters)
{

  /* declarations */

  char line[MAX_LINE_LENGTH];
  int dim, ret;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  if (debug) printf("Reading parameter file\n");
  
  /* read until out of lines */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read Parameters. */

    ret += sscanf(line, "Rank = %d", &Parameters->Rank);
    ret += sscanf(line, "GridRefinement = %d", &Parameters->GridRefinement);
    ret += sscanf(line, "ParticleRefinement = %d", 
		  &Parameters->ParticleRefinement);
    ret += sscanf(line, "GridDims = %d %d %d", Parameters->GridDims,
		  Parameters->GridDims+1, Parameters->GridDims+2);
    ret += sscanf(line, "ParticleDims = %d %d %d", Parameters->ParticleDims,
		  Parameters->ParticleDims+1, Parameters->ParticleDims+2);
    ret += sscanf(line, "MaxDims = %d %d %d", Parameters->MaxDims,
		  Parameters->MaxDims+1, Parameters->MaxDims+2);
    ret += sscanf(line, "NewCenter = %d %d %d", Parameters->NewCenter,
		  Parameters->NewCenter+1, Parameters->NewCenter+2);
    ret += sscanf(line, "NewCenterFloat = %f %f %f", 
		  Parameters->NewCenterFloat,
		  Parameters->NewCenterFloat+1, Parameters->NewCenterFloat+2);
    ret += sscanf(line, "StartIndex = %d %d %d", Parameters->StartIndex,
		  Parameters->StartIndex+1, Parameters->StartIndex+2);
    ret += sscanf(line, "StartIndexInNewCenterTopGridSystem = %d %d %d", 
		  Parameters->StartIndexInNewCenterTopGridSystem,
		  Parameters->StartIndexInNewCenterTopGridSystem+1, 
		  Parameters->StartIndexInNewCenterTopGridSystem+2);
    ret += sscanf(line, "EndIndexInNewCenterTopGridSystem = %d %d %d", 
		  Parameters->EndIndexInNewCenterTopGridSystem,
		  Parameters->EndIndexInNewCenterTopGridSystem+1, 
		  Parameters->EndIndexInNewCenterTopGridSystem+2);
    ret += sscanf(line, "RootGridDims = %d %d %d\n",
		  Parameters->RootGridDims, Parameters->RootGridDims+1,
		  Parameters->RootGridDims+2);

    ret += sscanf(line, "WaveNumberCutoff = %d",&Parameters->WaveNumberCutoff);
    ret += sscanf(line, "InitializeParticles = %d", 
		  &Parameters->InitializeParticles);
    ret += sscanf(line, "InitializeGrids = %d", &Parameters->InitializeGrids);

    if (sscanf(line, "ParticlePositionName = %s", dummy) == 1)
      Parameters->ParticlePositionName = dummy;
    if (sscanf(line, "ParticleVelocityName = %s", dummy) == 1)
      Parameters->ParticleVelocityName = dummy;
    if (sscanf(line, "ParticleMassName = %s", dummy) == 1)
      Parameters->ParticleMassName = dummy;
    if (sscanf(line, "GridDensityName = %s", dummy) == 1)
      Parameters->GridDensityName = dummy;
    if (sscanf(line, "GridVelocityName = %s", dummy) == 1)
      Parameters->GridVelocityName = dummy;

    /* If the dummy char space was used, then make another. */

    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }

    /* check to see if the line belongs to one of the test problems */

    if (strstr(line, "Cosmology")           ) ret++;
    if (strstr(line, "PowerSpectrum")       ) ret++;

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s", line);

  }

  /* clean up */

  delete dummy;

  /* Set parameters that were left undefined. */

  for (dim = 0; dim < Parameters->Rank; dim++)
    if (Parameters->MaxDims[dim] == INT_UNDEFINED)
      Parameters->MaxDims[dim] = max(Parameters->ParticleDims[dim],
				     Parameters->GridDims[dim]);
  if (Parameters->WaveNumberCutoff == INT_UNDEFINED)
    Parameters->WaveNumberCutoff = Parameters->MaxDims[0]/2;

  for (dim = Parameters->Rank; dim < 3; dim++) {
    Parameters->MaxDims[dim] = 1;
    Parameters->ParticleDims[dim] = 1;
    Parameters->GridDims[dim] = 1;
  }

  /* If NewCenterFloat was set, then set NewCenter. */

  int SetNewCenter = FALSE;
  for (dim = 0; dim < Parameters->Rank; dim++)
    if (Parameters->NewCenterFloat[dim] != FLOAT_UNDEFINED) {
      if (Parameters->NewCenter[dim] != INT_UNDEFINED) {
	fprintf(stderr, "You have set both NewCenterFloat and NewCenter -- use only one.\n");
	return FAIL;
      }
      if (Parameters->NewCenterFloat[dim] < 0 || 
	  Parameters->NewCenterFloat[dim] > 1) {
	fprintf(stderr, "NewCenterFloat must be >= 0 and < 1.\n");
	return FAIL;
      }
      if (Parameters->RootGridDims[dim] == INT_UNDEFINED) {
	fprintf(stderr, "You must set RootGridDims.\n");
	return FAIL;
      }
      Parameters->NewCenter[dim] = 
	nint(Parameters->NewCenterFloat[dim]*Parameters->RootGridDims[dim])*
	Parameters->MaxDims[dim]/Parameters->RootGridDims[dim] - 1;
      Parameters->NewCenter[dim] = 
	(Parameters->NewCenter[dim] + Parameters->MaxDims[dim]) %
	Parameters->MaxDims[dim];
      SetNewCenter = TRUE;
    }

  if (SetNewCenter)
    printf("Setting NewCenter = %d %d %d\n",
	   Parameters->NewCenter[0], Parameters->NewCenter[1],
	   Parameters->NewCenter[2]);

  /* If Start/EndIndexInNewCenterTopGridSystem was set, then use it to
     set StartIndex and GridDims. */

  int SetWithNewSystem = FALSE;
  for (dim = 0; dim < Parameters->Rank; dim++)
    if (Parameters->StartIndexInNewCenterTopGridSystem[dim] != INT_UNDEFINED) {
      if (Parameters->RootGridDims[dim] == INT_UNDEFINED) {
	fprintf(stderr, "You must set RootGridDims.\n");
	return FAIL;
      }
      SetWithNewSystem = TRUE;
      int MoveBy = Parameters->NewCenter[dim] - 
	(Parameters->MaxDims[dim]/2 - 1);
      Parameters->StartIndex[dim] = 
	(Parameters->StartIndexInNewCenterTopGridSystem[dim]*
	 Parameters->MaxDims[dim])/Parameters->RootGridDims[dim] +
	MoveBy;
      Parameters->StartIndex[dim] = 
	(Parameters->StartIndex[dim] + Parameters->MaxDims[dim]) %
	Parameters->MaxDims[dim];
      if (Parameters->InitializeGrids)
	Parameters->GridDims[dim] = 
	  (Parameters->EndIndexInNewCenterTopGridSystem[dim] -
	   Parameters->StartIndexInNewCenterTopGridSystem[dim]+1)*
	  (Parameters->MaxDims[dim]/(Parameters->RootGridDims[dim]*
				     Parameters->GridRefinement));
      if (Parameters->InitializeParticles)
	Parameters->ParticleDims[dim] = 
	  (Parameters->EndIndexInNewCenterTopGridSystem[dim] -
	   Parameters->StartIndexInNewCenterTopGridSystem[dim]+1)*
	  (Parameters->MaxDims[dim]/(Parameters->RootGridDims[dim]*
				     Parameters->ParticleRefinement));
				     
    }

  if (SetWithNewSystem) {
    printf("Setting StartIndex = %d %d %d\n",
	   Parameters->StartIndex[0], Parameters->StartIndex[1],
	   Parameters->StartIndex[2]);
    printf("Setting GridDims = %d %d %d\n",
	   Parameters->GridDims[0], Parameters->GridDims[1],
	   Parameters->GridDims[2]);
    printf("Setting ParticleDims = %d %d %d\n",
	   Parameters->ParticleDims[0], Parameters->ParticleDims[1],
	   Parameters->ParticleDims[2]);
  }

  return SUCCESS;
}
