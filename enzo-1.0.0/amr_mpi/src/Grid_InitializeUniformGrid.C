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
/  GRID CLASS (INITIALIZE THE GRID TO A UNIFORM POOL OF GAS)
/
/  written by: Greg Bryan
/  date:       February, 1995
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

int grid::InitializeUniformGrid(float UniformDensity,
				float UniformTotalEnergy,
				float UniformInternalEnergy,
				float UniformVelocity[])
{
  /* declarations */

  int dim, i, size, field;

  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* compute size of fields */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* allocate fields */

  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];

  /* set density, total energy */

  for (i = 0; i < size; i++) {
    BaryonField[0][i] = UniformDensity;
    BaryonField[1][i] = UniformTotalEnergy;
  }

  /* set velocities */

  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      BaryonField[vel+dim][i] = UniformVelocity[dim];

  /* Set internal energy if necessary. */

  if (DualEnergyFormalism)
    for (i = 0; i < size; i++)
      BaryonField[2][i] = UniformInternalEnergy;

  return SUCCESS;
}
