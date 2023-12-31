/***********************************************************************
/
/  FROM THE ATTRIBUTES, DETERMINE WHETHER THE STAR DIES
/
/  written by: John Wise
/  date:       April, 2009
/  modified1:
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

#define NO_DEATH 0
#define KILL_STAR 1
#define KILL_ALL 2

int Star::HitEndpoint(FLOAT Time)
{

  const float PISNLowerMass = 140, PISNUpperMass = 260;

  /* First check if the star's past its lifetime and then check other
     constrains based on its star type */

  int result = NO_DEATH;
  if ((Time > this->BirthTime + this->LifeTime) && this->type >=0)
    result = KILL_STAR;
  else
    return result;

  switch (this->type) {

  case PopIII:
    // If a Pop III star is going supernova, only kill it after it has
    // applied its feedback sphere
    if (this->Mass >= PISNLowerMass && this->Mass <= PISNUpperMass)
      if (this->FeedbackFlag == DEATH) {
	this->Mass = 1e-10;  // Needs to be non-zero

	// Set lifetime so the time of death is exactly now.
	this->LifeTime = Time - this->BirthTime;

	this->FeedbackFlag = NO_FEEDBACK;
	//result = KILL_STAR;
	result = NO_DEATH;
      } else
	result = NO_DEATH;

    // Check mass: Don't want to kill tracer SN particles formed
    // (above) in the previous timesteps.

    else if (this->Mass > 1e-9) {
      // Turn particle into a black hole (either radiative or tracer)
      if (PopIIIBlackHoles) {
	this->type = BlackHole;
	this->LifeTime = huge_number;
	this->FeedbackFlag = NO_FEEDBACK;
	result = NO_DEATH;
      } else {
	this->type = PARTICLE_TYPE_DARK_MATTER;
	result = KILL_STAR;
      }
    } else // SN tracers (must refine)
      result = NO_DEATH;
    break;
    
  case PopII:
    break;

  case BlackHole:
    break;
    
  case PopIII_CF:
    break;

  } // ENDSWITCH

  return result;
}
