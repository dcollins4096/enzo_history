/***********************************************************************
/
/  GRID CLASS (CONSTRUCTOR)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
//
//  Grid constructor (Set all data to null/default state).
//

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

grid::grid()
{

  /* clear scalars */

  GridRank                              = 0;
  Time                                  = 0.0;
  OldTime                               = 0.0;
  NumberOfBaryonFields                  = 0;
  dtFixed                               = 0.0;
  NumberOfParticles                     = 0;
  GravitatingMassFieldCellSize          = FLOAT_UNDEFINED;
  GravitatingMassFieldParticlesCellSize = FLOAT_UNDEFINED;
  SubgridsAreStatic                     = FALSE;
  ProcessorNumber                       = ROOT_PROCESSOR;

  /* clear MAX_DIMENSION vectors */

  int i, j;
  for (i = 0; i < MAX_DIMENSION; i++) {
    GridDimension[i]                 = 1;
    GridStartIndex[i]                = 0;
    GridEndIndex[i]                  = 0;
    GridLeftEdge[i]                  = DomainLeftEdge[i];
    GridRightEdge[i]                 = DomainRightEdge[i];
    CellLeftEdge[i]                  = NULL;
    CellWidth[i]                     = NULL;
    ParticlePosition[i]              = NULL;
    ParticleVelocity[i]              = NULL;
    ParticleAcceleration[i]          = NULL;
    AccelerationField[i]             = NULL;
    GravitatingMassFieldDimension[i] = 0;
  }

  ParticleAcceleration[MAX_DIMENSION]      = NULL;

  /* clear MAX_NUMBER_OF_BARYON_FIELDS vectors & [][MAX_DIMENSION] matricies */

  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    BaryonField[i]          = NULL;
    OldBaryonField[i]       = NULL;
    FieldType[i]            = FieldUndefined;
    for (j = 0; j < MAX_DIMENSION; j++ ) {
      BoundaryFluxes.LeftFluxes[i][j]  = NULL;
      BoundaryFluxes.RightFluxes[i][j] = NULL;
    }
  }

  /* Clear miscelaneous pointers */

  ParticleMass                  = NULL;
  ParticleNumber                = NULL;
  PotentialField                = NULL;
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;
  GravityBoundaryType           = GravityUndefined;
  for (i = 0; i < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; i++)
    ParticleAttribute[i] = NULL;

  /* Clear flagging field pointers */

  MassFlaggingField             = NULL;
  FlaggingField                 = NULL;
  
}
