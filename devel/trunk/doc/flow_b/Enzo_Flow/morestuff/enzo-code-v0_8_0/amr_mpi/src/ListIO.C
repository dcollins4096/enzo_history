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

#ifdef FAIL
#undef FAIL
#endif
#define FAIL      0
#define SUCCESS   1

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

int ReadListOfFloats(FILE *fptr, int N, float floats[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%f", floats + i) != 1)
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


void WriteListOfFloats(FILE *fptr, int N, double floats[])
{
  for (int i = 0; i < N; i++)
    fprintf(fptr, "%.14g ", floats[i]);
  fprintf(fptr, "\n");
}

int ReadListOfFloats(FILE *fptr, int N, double floats[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%lf", floats + i) != 1)
      return FAIL;

  fscanf(fptr, "\n");
  return SUCCESS;
}

void WriteListOfFloats(FILE *fptr, int N, long double floats[])
{
  for (int i = 0; i < N; i++)
    fprintf(fptr, "%.21Lg ", floats[i]);
  fprintf(fptr, "\n");
}

int ReadListOfFloats(FILE *fptr, int N, long double floats[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%Lf", floats + i) != 1)
      return FAIL;

  fscanf(fptr, "\n");
  return SUCCESS;

}
