/***********************************************************************
/
/  GRID CLASS (HANDLE CALLING AND SOLVING COOLING/CHEMISTRY)
/
/  written by: Matthew Turk
/  date:       June, 2009
/  modified1:
/
/  PURPOSE: Move logic for chemistry/cooling module selection here
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::MultiSpeciesHandler()
{
  if ((!MultiSpecies) && (!RadiativeCooling)) return SUCCESS; 
  if (GadgetEquilibriumCooling != 0) return SUCCESS;

  JBPERF_START("grid_MultiSpeciesHandler");

  if (MultiSpecies && RadiativeCooling) {
    this->SolveRateAndCoolEquations();
  } else {
    if (MultiSpecies)
      this->SolveRateEquations();
    if (RadiativeCooling)
      this->SolveRadiativeCooling();
  }

  if (ProblemType == 62)
    this->CoolingTestResetEnergies();

  JBPERF_STOP("grid_MultiSpeciesHandler");
  return SUCCESS;
}
