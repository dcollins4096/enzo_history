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
/  GRID CLASS (COPY PARTICLES INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"


int grid::CommunicationTransferParticles(grid* Grids[], int NumberOfGrids, 
		 int ToGrid[6], int NumberToMove[6], 
		 float_int *ParticleData[6], int CopyDirection)
{

  /* If there are no particles to move, we're done. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Declarations. */

  int i, j, k, dim, grid;

  /* ----------------------------------------------------------------- */
  /* Copy particle out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* Count particles to move (mark already counted by setting mass < 0). */

    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < NumberOfParticles; i++)
	if (ParticleMass[i] >= 0) {
	  if (ParticlePosition[dim][i] < GridLeftEdge[dim]) {
	    NumberToMove[dim*2+0]++;
	    ParticleMass[i] = -ParticleMass[i];
	  }
	  if (ParticlePosition[dim][i] > GridRightEdge[dim]) {
	    NumberToMove[dim*2+1]++;
	    ParticleMass[i] = -ParticleMass[i];
	  }
	}

    /* Allocate space. */

    int TotalToMove = 0, NumberOfParticleFields = 8+NumberOfParticleAttributes;
    for (i = 0; i < 6; i++) {
      TotalToMove += NumberToMove[i];
      if (NumberToMove[i] > 0) {
	ParticleData[i] = new float_int[NumberToMove[i]*
				        NumberOfParticleFields];
      } else {
	ParticleData[i] = NULL;
      }
    }

    /* Set ToGrid. */

    for (dim = 0; dim < GridRank; dim++) {
      int DimSize = nint((DomainRightEdge[dim] - 
		      DomainLeftEdge[dim])/CellWidth[dim][0]);

      /* Find Left grid */

      for (grid = 0; grid < NumberOfGrids; grid++) {
	int FoundIt = TRUE;
	for (i = 0; i < GridRank; i++) {
	  if (dim != i && nint(GridLeftEdge[i]/CellWidth[i][0]) !=
	      nint(Grids[grid]->GridLeftEdge[i]/CellWidth[i][0]))
	    FoundIt = FALSE;
	  if (dim == i && (nint(
               Grids[grid]->GridRightEdge[i]/CellWidth[i][0]) % DimSize)
	      != nint(GridLeftEdge[i]/CellWidth[i][0]))
	    FoundIt = FALSE;
	}
	if (FoundIt) {
	  ToGrid[dim*2+0] = grid;
	  break;
	}
      }

      /* Find right grid */

      for (grid = 0; grid < NumberOfGrids; grid++) {
	int FoundIt = TRUE;
	for (i = 0; i < GridRank; i++) {
	  if (dim != i && nint(GridLeftEdge[i]/CellWidth[i][0]) !=
	      nint(Grids[grid]->GridLeftEdge[i]/CellWidth[i][0]))
	    FoundIt = FALSE;
	  if (dim == i && (nint(
               GridRightEdge[i]/CellWidth[i][0]) % DimSize)
	      != nint(Grids[grid]->GridLeftEdge[i]/CellWidth[i][0]))
	    FoundIt = FALSE;
	}
	if (FoundIt) {
	  ToGrid[dim*2+1] = grid;
	  break;
	}
      }

    } // end loop over dims

    if (TotalToMove == 0)
      return SUCCESS;
    
    /* Move particles (mark those moved by setting mass = FLOAT_UNDEFINED). */

    for (dim = 0; dim < GridRank; dim++) {
      int n1 = 0, n2 = 0;
      for (i = 0; i < NumberOfParticles; i++) 
	if (ParticleMass[i] != FLOAT_UNDEFINED) {

	  /* shift left. */

	  if (ParticlePosition[dim][i] < GridLeftEdge[dim]) {
	    for (j = 0; j < GridRank; j++) {
	      ParticleData[dim*2][n1++].FVAL = ParticlePosition[j][i];
	      ParticleData[dim*2][n1++].fval = ParticleVelocity[j][i];
	    }
	    ParticleData[dim*2][n1++].fval = -ParticleMass[i];
	    ParticleData[dim*2][n1++].ival = ParticleNumber[i];
	    for (j = 0; j < NumberOfParticleAttributes; j++)
	      ParticleData[dim*2][n1++].fval = ParticleAttribute[j][i];
	    ParticleMass[i] = FLOAT_UNDEFINED;
	  }

	  /* shift right. */

	  if (ParticlePosition[dim][i] > GridRightEdge[dim]) {
	    for (j = 0; j < MAX_DIMENSION; j++) {
	      ParticleData[dim*2+1][n2++].FVAL = ParticlePosition[j][i];
	      ParticleData[dim*2+1][n2++].fval = ParticleVelocity[j][i];
	    }
	    ParticleData[dim*2+1][n2++].fval = -ParticleMass[i];
	    ParticleData[dim*2+1][n2++].ival = ParticleNumber[i];
	    for (j = 0; j < NumberOfParticleAttributes; j++)
	      ParticleData[dim*2+1][n2++].fval = ParticleAttribute[j][i];
	    ParticleMass[i] = FLOAT_UNDEFINED;
	  }

	} // end: if (ParticleMass[i] != FLOAT_UNDEFINED)
    } // end: loop over dims

  } // end: if (COPY_OUT)

  /* ----------------------------------------------------------------- */
  /* Copy particle back into grid. */

  else {

    /* Count up total number. */

    int TotalNumberOfParticles = NumberOfParticles;

    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleMass[i] == FLOAT_UNDEFINED)
	TotalNumberOfParticles--;

    for (i = 0; i < 6; i++)
      TotalNumberOfParticles += NumberToMove[i];

    if (TotalNumberOfParticles == 0 && NumberOfParticles == 0)
      return SUCCESS;

    /* Allocate space for the particles. */

    FLOAT *Position[MAX_DIMENSION];
    float *Velocity[MAX_DIMENSION], *Mass,
          *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    int *Number;

    Mass = new float[TotalNumberOfParticles];
    Number = new int[TotalNumberOfParticles];
    for (dim = 0; dim < GridRank; dim++) {
      Position[dim] = new FLOAT[TotalNumberOfParticles];
      Velocity[dim] = new float[TotalNumberOfParticles];
    }
    for (i = 0; i < NumberOfParticleAttributes; i++)
      Attribute[i] = new float[TotalNumberOfParticles];
    if (Velocity[GridRank-1] == NULL && TotalNumberOfParticles != 0) {
      fprintf(stderr, "malloc error (out of memory?)\n");
      return FAIL;
    }

    /* Copy this grid's particles to the new space (don't copy any erased 
       particles, those with Mass == FLOAT_UNDEFINED). */

    int n = 0;
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleMass[i] != FLOAT_UNDEFINED) {
	Mass[n]   = ParticleMass[i];
	Number[n] = ParticleNumber[i];

	for (dim = 0; dim < GridRank; dim++) {
	  Position[dim][n] = ParticlePosition[dim][i];
	  Velocity[dim][n] = ParticleVelocity[dim][i];
	}
	for (j = 0; j < NumberOfParticleAttributes; j++)
	  Attribute[j][n] = ParticleAttribute[j][i];
	n++;
      }

    /* Copy new particles (and wrap around edges). */

    for (j = 0; j < 6; j++)
      if (NumberToMove[j] > 0) {
	int n2 = 0;
	for (i = 0; i < NumberToMove[j]; i++) {

	  for (dim = 0; dim < GridRank; dim++) {
	    Position[dim][n] = ParticleData[j][n2++].FVAL;
	    if (Position[dim][n] > DomainRightEdge[dim])
	      Position[dim][n] -= DomainRightEdge[dim] - DomainLeftEdge[dim];
	    if (Position[dim][n] < DomainLeftEdge[dim])
	      Position[dim][n] += DomainRightEdge[dim] - DomainLeftEdge[dim];
	    Velocity[dim][n] = ParticleData[j][n2++].fval;
	  }

	  Mass[n]   = ParticleData[j][n2++].fval;
	  Number[n] = ParticleData[j][n2++].ival;
	  for (k = 0; k < NumberOfParticleAttributes; k++)
	    Attribute[k][n] = ParticleData[j][n2++].fval;

	  n++;
	}
      }

    /* Set new number of particles in this grid. */
    
    NumberOfParticles = TotalNumberOfParticles;

    /* Delete old pointers and copy new ones into place. */

    this->DeleteParticles();
    this->SetParticlePointers(Mass, Number, Position, Velocity, Attribute);

  } // end: if (COPY_IN)

  return SUCCESS;
}
