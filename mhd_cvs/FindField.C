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
/  FIND FIELD FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: field index or -1 on failure
/
************************************************************************/

// Find field type field in array field_type, returning the index into the
//   field array or -1 if it is not there.

#include "macros_and_parameters.h"
#include "typedefs.h"


int FindField(int field, int farray[], int numfields)
{
  for (int i = 0; i < numfields; i++)
    if (field == farray[i]) 
      return i;

  /* not found */

  return -1;
}
