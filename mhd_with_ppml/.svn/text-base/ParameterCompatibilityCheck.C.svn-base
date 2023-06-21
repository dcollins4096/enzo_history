#include <math.h>
#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"


int ParameterCompatibilityCheck(TopGridData &MetaData){

  int Check = SUCCESS;
  if( HydroMethod != PPM_Local ){
    
    if( MHD_Used == TRUE ){
      fprintf(stderr,"Incompatable Parameters: MHD_Used = TRUE and HydroMethod %"ISYM"\n",
	      HydroMethod);
      Check = FAIL;
    }
    if( EquationOfState == 1 ){
      fprintf(stderr,"Incompatable Parameters: EquationOfState = 1 and HydroMethod %"ISYM"\n",
	      HydroMethod);
      Check = FAIL;
    }
    if( DEFAULT_GHOST_ZONES != 3 ){
      fprintf(stderr,"Incompatable Parameters: DEFAULT_GHOST_ZONES = %"ISYM" and HydroMethod %"ISYM"\n",
	      DEFAULT_GHOST_ZONES, HydroMethod);
      Check = FAIL;
    }
  }

  /* Unlock when Athena gets installed.
  if( EquationOfState == 1 && MHD_ReconField[1] > 0 ){
    fprintf(stderr, "Parameter Consistency Warning: MHD_ReconField[1] > 0 && EquaitonOfState == 1.\n");
    fprintf(stderr, "Setting MHD_ReconField[1] = -1.\n");
    MHD_ReconField[1] = -1;
  }
  */
  return Check;
}
