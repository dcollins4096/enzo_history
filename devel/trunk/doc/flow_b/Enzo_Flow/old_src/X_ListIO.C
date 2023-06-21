/***********************************************************************
/
/  READ/WRITE LIST OF INTS/FLOATS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"

int ReadListOfInts(FILE *fptr, int N, int nums[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%d", nums + i) != 1)
      return FAIL;

  fscanf(fptr, "\n");
  return SUCCESS;
}

void WriteListOfInts(FILE *fptr, int N, int nums[])
{
  for (int i = 0; i < N; i++)
    fprintf(fptr, "%d ", nums[i]);
  fprintf(fptr, "\n");
}

#if (defined(p4))

int ReadListOfFloats(FILE *fptr, int N, float floats[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%"FSYM, floats + i) != 1)
      return FAIL;

  fscanf(fptr, "\n");
  return SUCCESS;
}

void WriteListOfFloats(FILE *fptr, int N, float floats[])
{
  for (int i = 0; i < N; i++)
    fprintf(fptr, "%.7g ", floats[i]);
  fprintf(fptr, "\n");
}

#endif

/* If using double or quad precision, then define high-precision versions. */

#if (defined(p8) || defined(p16))

void WriteListOfFLOATS(FILE *fptr, int N, FLOAT floats[])
{
  for (int i = 0; i < N; i++)
    fprintf(fptr, "%"GOUTSYM" ", floats[i]);
  fprintf(fptr, "\n");
}

int ReadListOfFLOATS(FILE *fptr, int N, FLOAT floats[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%"PSYM, floats + i) != 1)
      return FAIL;

  fscanf(fptr, "\n");
  return SUCCESS;
}

#endif
