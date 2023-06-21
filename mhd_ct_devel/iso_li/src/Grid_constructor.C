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
    RandomForcingField[i]            = NULL;
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


  //Whether you use them or not, these should be nullified.
  for(i=0;i<3;i++){
    MagneticField[i] = NULL;
    CenteredB[i]     = NULL;
    ElectricField[i] = NULL;
    AvgElectricField[i] = NULL;
    OldMagneticField[i] = NULL;
    OldElectricField[i] = NULL;
    OldCenteredB[i] = NULL;
    OldAccelerationField[i] = NULL;
    Current[i] = NULL;
    MHDParentTemp[i] = NULL;

    for(j=0;j<3;j++){
      BoundaryFluxes.LeftElectric[i][j]=NULL;
      BoundaryFluxes.RightElectric[i][j]=NULL;
      BoundaryFluxes.ElectricSize[i][j]=0;
    }
  }


#ifdef ATHENA  //HX 06/28/2006
  MHD_Pressure    = NULL;
#endif //ATHENA
  AccelerationHack = FALSE;

  DivB = NULL;
  dtParent = -1;

  DyBx = NULL;
  DzBx = NULL;
  DyzBx = NULL;
  DBxFlag = NULL;

  DxBy = NULL;
  DzBy = NULL;
  DxzBy = NULL;
  DByFlag = NULL;

  DxBz = NULL;
  DyBz = NULL;
  DxyBz = NULL;
  DBzFlag = NULL;

  for(int field=0;field<3;field++){
    MagneticSize[field] = -100;
    ElectricSize[field] = -100;
    for(int dim=0;dim<3;dim++){
      MHDAdd[field][dim]=(field==dim) ? 1:0;
      MagneticDims[field][dim] = -100;
    }}
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
