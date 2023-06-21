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
/  SET THE PARAMETER DEFAULTS
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

/* function prototypes */

int SetParameterDefaults(parmstruct *Parameters)
{

  int dim;

  char ppos_name[] = "ParticlePositions";
  char pvel_name[] = "ParticleVelocities";
  char gden_name[] = "GridDensity";
  char gvel_name[] = "GridVelocity";

  /* Set Defaults. */

  Parameters->Rank                = 3;
  Parameters->GridRefinement      = 1;
  Parameters->ParticleRefinement  = 1;
  Parameters->WaveNumberCutoff    = INT_UNDEFINED;
  Parameters->InitializeParticles = TRUE;
  Parameters->InitializeGrids     = TRUE;

  Parameters->ParticlePositionName = ppos_name;
  Parameters->ParticleVelocityName = pvel_name;
  Parameters->ParticleMassName     = NULL;
  Parameters->GridDensityName      = gden_name;
  Parameters->GridVelocityName     = gvel_name;

  for (dim = 0; dim < Parameters->Rank; dim++) {
    Parameters->GridDims[dim]     = INT_UNDEFINED;
    Parameters->ParticleDims[dim] = INT_UNDEFINED;
    Parameters->MaxDims[dim]      = INT_UNDEFINED;
    Parameters->NewCenter[dim]    = INT_UNDEFINED;
    Parameters->NewCenterFloat[dim] = FLOAT_UNDEFINED;
    Parameters->StartIndex[dim]   = 0;
    Parameters->StartIndexInNewCenterTopGridSystem[dim] = INT_UNDEFINED;
    Parameters->EndIndexInNewCenterTopGridSystem[dim] = INT_UNDEFINED;
    Parameters->RootGridDims[dim] = INT_UNDEFINED;
  }

  return SUCCESS;
}
