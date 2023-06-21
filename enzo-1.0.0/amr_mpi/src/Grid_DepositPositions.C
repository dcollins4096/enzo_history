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
/  GRID CLASS (DEPOSIT POSITIONS ONTO THE GRAVITATING MASS FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE: 
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

/* This variable is only set in ReadParameterFile. */

float DepositPositionsParticleSmoothRadius = 0;

/* function prototypes */

extern "C" void FORTRAN_NAME(cic_deposit)(FLOAT *posx, FLOAT *posy, 
			FLOAT *posz, int *ndim, int *npositions, 
                        float *densfield, float *field, FLOAT *leftedge, 
                        int *dim1, int *dim2, int *dim3, float *cellsize);

extern "C" void FORTRAN_NAME(smooth_deposit)(FLOAT *posx, FLOAT *posy, 
			FLOAT *posz, int *ndim, int *npositions, 
                        float *densfield, float *field, FLOAT *leftedge, 
                        int *dim1, int *dim2, int *dim3, float *cellsize,
			       float *rsmooth);


int grid::DepositPositions(FLOAT *Position[], float *Mass, int Number, 
			   int DepositField)
{
  if (Number == 0) return SUCCESS;

  /* DepositField specifies where the particles should go.  Set LeftEdge, 
     Dimension, CellSize, DepositFieldPointer according to it. */

  float *DepositFieldPointer, CellSize;
  FLOAT LeftEdge[MAX_DIMENSION];
  int   dim, Dimension[MAX_DIMENSION];

  /* 1) GravitatingMassField. */

  if (DepositField == GRAVITATING_MASS_FIELD) {
    DepositFieldPointer = GravitatingMassField;
    CellSize            = float(GravitatingMassFieldCellSize);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = GravitatingMassFieldLeftEdge[dim];
      Dimension[dim] = GravitatingMassFieldDimension[dim];
    }
  }

  /* 2) GravitatingMassFieldParticles. */

  else if (DepositField == GRAVITATING_MASS_FIELD_PARTICLES) {
    DepositFieldPointer = GravitatingMassFieldParticles;
    CellSize            = float(GravitatingMassFieldParticlesCellSize);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = GravitatingMassFieldParticlesLeftEdge[dim];
      Dimension[dim] = GravitatingMassFieldParticlesDimension[dim];
    }
  }

  /* 3) MassFlaggingField */

  else if (DepositField == MASS_FLAGGING_FIELD) {
    DepositFieldPointer = MassFlaggingField;
    CellSize            = float(CellWidth[0][0]);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = CellLeftEdge[dim][0];
      Dimension[dim] = GridDimension[dim];
    }
  }

  /* 4) error */

  else {
    fprintf(stderr, "DepositField = %d not recognized.\n", DepositField);
    return FAIL;
  }

  /* Error check. */

  if (DepositFieldPointer == NULL) {
    fprintf(stderr, "DepositFieldPointer (%d) is NULL, Number = %d.\n",
	    DepositField, Number);
    return FAIL;
  }

  if (GridRank != 3) {
    fprintf(stderr, "New gravity module currently supports only 3d.\n");
    return FAIL;
  }

  if (DepositPositionsParticleSmoothRadius < CellSize)

    /* Deposit to field using CIC. */

    FORTRAN_NAME(cic_deposit)(Position[0], Position[1], Position[2], &GridRank,
			      &Number, Mass, DepositFieldPointer, LeftEdge,
			      Dimension, Dimension+1, Dimension+2,
			      &CellSize);

  else

    /* Deposit to field using large-spherical CIC, with radius of
       DepositPositionsPartaicleSmoothRadius . */

    FORTRAN_NAME(smooth_deposit)(
			  Position[0], Position[1], Position[2], &GridRank,
			  &Number, Mass, DepositFieldPointer, LeftEdge,
			  Dimension, Dimension+1, Dimension+2,
			  &CellSize, &DepositPositionsParticleSmoothRadius);

  return SUCCESS;
}
