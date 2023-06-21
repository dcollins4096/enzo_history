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
				float UniformVelocity[],
				float *UniformMagneticField) // default = NULL
{
  /* declarations */

  int dim, i, size, field;

#ifndef ATHENA
 int EquationOfState = 0;
#endif

  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  if( EquationOfState == 0 ){
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1 || MHD_Used == TRUE) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2 || MHD_Used == TRUE)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* compute size of fields */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* allocate fields */

  //Variable names.
  int Eeng, Eden, Ev[3], Egas; 
  if (this->IdentifyPhysicalQuantities(Eden, Egas, Ev[0], Ev[1], 
				       Ev[2], Eeng) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }  

  this->AllocateGrids();

  /* set density, total energy */

  for (i = 0; i < size; i++) {
    BaryonField[ Eden ][i] = UniformDensity;
    if( EquationOfState == 0 )
      BaryonField[ Eeng ][i] = UniformTotalEnergy;
  }

  /* set velocities */

  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      BaryonField[vel+dim][i] = UniformVelocity[dim];

  /* Set internal energy if necessary. */

  if (DualEnergyFormalism)
    for (i = 0; i < size; i++)
      BaryonField[2][i] = UniformInternalEnergy;


  if( MHD_Used ){
    for(field=0;field<3;field++){
      for(i=0;i<MagneticSize[field];i++)
	MagneticField[field][i] = 
	  ((UniformMagneticField != NULL)? UniformMagneticField[field]:1e-5);
      for(i=0;i<size;i++)
	CenteredB[field][i]=
	  ((UniformMagneticField != NULL)? UniformMagneticField[field]:1e-5);
    }      
  }
  
  fprintf(stderr,"Sedov: D %f E %f B %f %f %f\n",
	  BaryonField[0][0],BaryonField[1][0], MagneticField[0][0]
	  , MagneticField[1][0], MagneticField[2][0]);

  return SUCCESS;
}
