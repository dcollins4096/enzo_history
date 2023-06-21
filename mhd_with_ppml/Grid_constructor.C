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

#include <stdio.h>
#include <stdlib.h>
 
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
#ifndef PPML
    RandomForcingField[i]            = NULL;
#endif //! ppml
  }
#ifdef PPML
    for( i=0; i< MAX_NUMBER_RANDOM_FORCING ; i++)
      RandomForcingField[i]            = NULL;
#endif //PPML
    
#ifdef MHDF
  NFaces = 3;
  NEdges = 3;

  for( j=0; j<MAX_FACE_FIELDS; j++){
    MagneticField[j] = NULL;
    MagneticSize[j] = 1;
    for( i=0; i<MAX_DIMENSION; i++){
      MagneticDims[j][i]  = 1;
      MHDAdd[j][i] = ( j == i ) ? 1 : 0 ;
    }}
  for( j=0; j<MAX_EDGE_FIELDS; j++) {
    ElectricSize[j] =  1;
    ElectricField[j] = NULL;
    for( i=0; i<MAX_DIMENSION; i++){    
      ElectricDims[j][i] = 1;
    }
  }//j
#endif //MHDF

  ParticleAcceleration[MAX_DIMENSION]      = NULL;
 
  /* clear MAX_NUMBER_OF_BARYON_FIELDS vectors & [][MAX_DIMENSION] matricies */
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    BaryonField[i]          = NULL;
    OldBaryonField[i]       = NULL;
    FieldType[i]            = FieldUndefined;
  }

/*
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    for (j = 0; j < MAX_DIMENSION; j++ ) {
      BoundaryFluxes->LeftFluxes[i][j]  = NULL;
      BoundaryFluxes->RightFluxes[i][j] = NULL;
    }
  }
*/

  for (i = 0; i < MAX_DIMENSION; i++) {
    OldAccelerationField[i] = NULL;
  }

  AccelerationHack = FALSE;

  /* Clear miscelaneous pointers */
 
  ParticleMass                  = NULL;
  ParticleNumber                = NULL;
  ParticleType                  = NULL;
  PotentialField                = NULL;
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;
  GravityBoundaryType           = GravityUndefined;
  for (i = 0; i < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; i++)
    ParticleAttribute[i] = NULL;

  BoundaryFluxes                = NULL;
 
  /* Clear flagging field pointers */
 
  MassFlaggingField             = NULL;
  FlaggingField                 = NULL;

#ifdef ISO_GRAV
  /* Initialize top level parallelism information */
  for (i=0; i<MAX_DIMENSION; i++) {
    ProcLayout[i]       = 1;
    ProcLocation[i]     = 0;
    ProcNeighbors[i][0] = 0;
    ProcNeighbors[i][1] = 0;
  }
#endif

#ifdef PPML
  InterfaceStatesInitialized = FALSE;
  NumberOfFluidQuantities    = 0;
  PPML_NFaces                = 0;
#endif //PPML

}
