/***********************************************************************
/
/  FROM THE ATTRIBUTES, DETERMINE THE STAR TYPE AND FEEDBACK MODE
/
/  written by: John Wise
/  date:       November, 2005
/  modified1: Ji-hoon Kim
/             July, 2009
/
/  NOTES:  When the star particle is created, it is assigned the 
/          ParticleType from the grid.  Change this to represent the 
/          proper star_type (typedefs.h)
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
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
#include "StarParticleData.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

void Star::SetFeedbackFlag(int flag)
{
  this->FeedbackFlag = flag;
  return;
}

#ifdef LARGE_INTS
void Star::SetFeedbackFlag(Eint32 flag)
{
  this->FeedbackFlag = flag;
  return;
}
#endif

int Star::SetFeedbackFlag(FLOAT Time)
{

  const float PISNLowerMass = 140, PISNUpperMass = 260;
  const float StarClusterSNeStart = 4.0;   // Myr after cluster is born
  const float StarClusterSNeEnd = 20.0; // Myr (lifetime of a 8 Msun star)
  const double G = 6.673e-8, k_b = 1.38e-16, m_h = 1.673e-24;
  const double Msun = 1.989e33;

  int abs_type;
  float AgeInMyr;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  abs_type = ABS(this->type);
  switch (abs_type) {

  case PopIII:
    if (this->type < 0) // birth
      this->FeedbackFlag = FORMATION;
    else if (Time > this->BirthTime + this->LifeTime) // endpoint
      if (this->Mass >= PISNLowerMass && this->Mass <= PISNUpperMass)
	this->FeedbackFlag = SUPERNOVA;
      else
	this->FeedbackFlag = NO_FEEDBACK; // BH formation
    else // main sequence
      this->FeedbackFlag = NO_FEEDBACK;
    break;
    
  case PopII:
    AgeInMyr = (Time - BirthTime) * TimeUnits / 3.15e13;
    if (this->type > 0)
      if (AgeInMyr > StarClusterSNeStart && AgeInMyr < StarClusterSNeEnd)
	this->FeedbackFlag = CONT_SUPERNOVA;
      else
	this->FeedbackFlag = NO_FEEDBACK;
    else
      this->FeedbackFlag = FORMATION;
    break;

  case BlackHole:
    this->FeedbackFlag = NO_FEEDBACK;
    break;

  /* For MBH particle. Even with the NO_FEEDBACK flag, 
     the particle still can act as a Radiation Source if RadiativeTransfer = 1. */  

  case MBH:
    AgeInMyr = (Time - BirthTime) * TimeUnits / 3.15e13;
    if (this->type > 0 && AgeInMyr > 0 && MBHFeedbackThermal)
      this->FeedbackFlag = MBH_THERMAL;
    else
      this->FeedbackFlag = NO_FEEDBACK; //It could still be a Radiation Source. 

#define NOT_SEDOV_TEST
#ifdef SEDOV_TEST
  //if (this->type > 0 && AgeInMyr > 0 && AgeInMyr < 0.001)  //for Sedov test (injecting for 1kyr)
    if (this->type > 0 && AgeInMyr > 0)                      //for Ostriker & McKee test (injecting continuously)
      this->FeedbackFlag = MBH_THERMAL;
    else
      this->FeedbackFlag = NO_FEEDBACK;      
#endif

    break;

  } // ENDSWITCH

  //this->type = abs_type;

  return SUCCESS;
}
