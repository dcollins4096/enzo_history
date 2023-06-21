/***********************************************************************
/
/  OUTPUT THE FIELD TO A FILE (FlexIO VERSION)
/
/  written by: Greg Bryan
/  date:       August, 1999
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include "macros_and_parameters.h"
#include "IO.hh"
#include "IEEEIO.hh"

/* function prototypes */

int WriteField(int Rank, int Dims[3], float *Field, char *Name, int Type)
{

  /* Allocate IO writer (open file). */

  IObase *writer;
  if (Type == 0)
    writer = new IEEEIO(Name, IObase::Create);
  else
    writer = new IEEEIO(Name, IObase::Append);

  if (!writer->isValid()) {
    fprintf(stderr, "Error opening file %s\n", Name);
    return FAIL;
  }

  /* Output field. */

  writer->write(IObase::Float32, Rank, Dims, (void *) Field);

  /* Close file. */

  delete writer;

  return SUCCESS;
}
