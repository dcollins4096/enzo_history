/***********************************************************************
/
/  WRITES UNITS TO AN OUTPUT FILE
/
/  written by: Elizabeth Tasker
/  date:       May, 2005
/  modified1:
/
/  PURPOSE:
/
/  NOTE: 
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "units.h"

int WriteUnits(FILE *fptr)
{

  /* write output to file */

  fprintf(fptr, "MassUnits = %g\n", GlobalMassUnits);
  fprintf(fptr, "DensityUnits    = %g\n", GlobalDensityUnits);
  fprintf(fptr, "TimeUnits    = %g\n", GlobalTimeUnits);
  fprintf(fptr, "LengthUnits  = %g\n", GlobalLengthUnits);
     
  return SUCCESS;
}
