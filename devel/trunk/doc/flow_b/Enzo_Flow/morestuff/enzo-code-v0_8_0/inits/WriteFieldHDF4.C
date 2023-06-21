/***********************************************************************
/
/  OUTPUT THE FIELD TO A FILE
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#ifdef USE_HDF4

#include <stdio.h>
#include <string.h>
#include <df.h>
#include "macros_and_parameters.h"

/* function prototypes */

int WriteFieldHDF4(int Rank, int Dims[3], float *Field, char *Name, int Type)
{

  int dim, ret;
  int32 OutDims[3];

  /* Reverse dim ordering since we are using fortran array ordering
     (actually, we don't have to do this here). */

  for (dim = 0; dim < Rank; dim++)
    OutDims[dim] = Dims[dim];
  /* OutDims[Rank-dim-1] = Dims[dim]; */

  /* Set dims. */

  if (DFSDsetdims(Rank, OutDims) == HDF_FAIL) {
    fprintf(stderr, "Error in DFSDsetdims.\n");
    exit(EXIT_FAILURE);
  }

  /* Output field. */

  if (Type == 0)
    ret = DFSDputdata(Name, Rank, OutDims, (VOIDP) Field);
  else
    ret = DFSDadddata(Name, Rank, OutDims, (VOIDP) Field);
  if (ret == HDF_FAIL) {
    fprintf(stderr, "Error in DFSDadddata.\n");
    exit(EXIT_FAILURE);
  }

  return SUCCESS;
}
#endif /* USE_HDF4 */
