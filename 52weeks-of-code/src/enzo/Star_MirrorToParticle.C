/***********************************************************************
/
/  REFLECTS CHANGES IN STAR PARTICLE IN NORMAL PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

void Star::MirrorToParticle(void)
{

  if (CurrentGrid == NULL)
    return;

  const double Msun = 1.989e33;
  int i, dim, place = -1;
  float MassConversion;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, CurrentGrid->Time);

  double dx = LengthUnits * CurrentGrid->CellWidth[0][0];
  MassConversion = (float) (dx*dx*dx * double(DensityUnits) / Msun);

  //  printf("CurrentGrid->NumberOfParticles = %d\n", CurrentGrid->NumberOfParticles); 

  // Find where this star particle is stored in main arrays
  for (i = 0; i < CurrentGrid->NumberOfParticles; i++) 
    if (CurrentGrid->ParticleNumber[i] == this->Identifier) {
      place = i;
      break;
    }

  if (place < 0) {
    ENZO_VFAIL("Star_MirrorToParticle: Cannot find matching particle.\n"
	   "StarID = %d, grid::NumberOfParticles = %d\n"
	   "Position        = %"FSYM" %"FSYM" %"FSYM"\n"
	   "Grid Left Edge  = %"FSYM" %"FSYM" %"FSYM"\n"
	   "Grid Right Edge = %"FSYM" %"FSYM" %"FSYM"\n",
	   this->Identifier, CurrentGrid->NumberOfParticles,
	   this->pos[0], this->pos[1], this->pos[2],
	   CurrentGrid->GridLeftEdge[0], CurrentGrid->GridLeftEdge[1],
	   CurrentGrid->GridLeftEdge[2],
	   CurrentGrid->GridRightEdge[0], CurrentGrid->GridRightEdge[1],
	   CurrentGrid->GridRightEdge[2])

  }

  // Change all particle data in favor of updated Star particle
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    CurrentGrid->ParticlePosition[dim][place] = this->pos[dim];
    CurrentGrid->ParticleVelocity[dim][place] = this->vel[dim];
  }
  CurrentGrid->ParticleMass[place] = (float)(this->Mass) / MassConversion;
  CurrentGrid->ParticleType[place] = this->type;
  CurrentGrid->ParticleAttribute[0][place] = this->BirthTime;
  CurrentGrid->ParticleAttribute[1][place] = this->LifeTime;
  
  return;
}
