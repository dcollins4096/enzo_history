#include <stdio.h>
#include <stdlib.h>

#define FAIL -1
#define PSYM "f"
#define FSYM "f"

int ReadListOfInts(FILE *fptr, int N, int nums[]);
int ReadListOfFloats(FILE *fptr, int N, float floats[]);


int main()
{

FILE *fptr;

int filecounter;

int Grid;
int ThisGrid, NextGrid;
int GridRank;
int GridDimension[3];
int GridStartIndex[3];
int GridEndIndex[3];
float GridLeftEdge[3];
float GridRightEdge[3];
float Time;
int SubgridsAreStatic;
int NumberOfBaryonFields;
int FieldType[16];
char *BaryonFileName;
float CourantSafetyNumber;
int PPMFlatteningParameter;
int PPMDiffusionParameter;
int PPMSteepeningParameter;
int NumberOfParticles;
char *ParticleFileName;
int GravityBoundaryType;
int SelfGravity;

BaryonFileName = new char[101];
ParticleFileName = new char[101];

fptr = fopen("XXX","r");

for(filecounter = 0; filecounter < 4; filecounter++)
{

  fscanf(fptr, "\n");

  if (fscanf(fptr, "Grid = %d\n", &Grid) != 1) {
    fprintf(stderr, "Error reading Grid.\n");
    return FAIL;
  }

  printf("Grid = %d\n", Grid);

  if (fscanf(fptr, "GridRank = %d\n", &GridRank) != 1) {
    fprintf(stderr, "Error reading GridRank.\n");
    return FAIL;
  }

  printf("GridRank = %d\n", GridRank);

  if (fscanf(fptr, "GridDimension = ") != 0) {
    fprintf(stderr, "Error reading GridDimension(0).\n");
    return FAIL;
  }

  if (ReadListOfInts(fptr, GridRank, GridDimension) == FAIL) {
    fprintf(stderr, "Error reading GridDimension(1).\n");
    return FAIL;
  }

  fscanf(fptr, "GridStartIndex = ");

  if (ReadListOfInts(fptr, GridRank, GridStartIndex) == FAIL) {
    fprintf(stderr, "Error reading GridStartIndex.\n");
    return FAIL;
  }

  fscanf(fptr, "GridEndIndex = ");

  if (ReadListOfInts(fptr, GridRank, GridEndIndex) == FAIL) {
    fprintf(stderr, "Error reading GridEndIndex.\n");
    return FAIL;
  }

  fscanf(fptr, "GridLeftEdge = ");

  if (ReadListOfFloats(fptr, GridRank, GridLeftEdge) == FAIL) {
    fprintf(stderr, "Error reading GridLeftEdge.\n");
    return FAIL;
  }

  fscanf(fptr, "GridRightEdge = ");

  if (ReadListOfFloats(fptr, GridRank, GridRightEdge) == FAIL) {
    fprintf(stderr, "Error reading GridRightEdge.\n");
    return FAIL;
  }

  printf("GridLeftEdge[0] = %10.4f  GridRightEdge[0] = %10.4f\n", GridLeftEdge[0], GridRightEdge[0]);
  printf("GridLeftEdge[1] = %10.4f  GridRightEdge[1] = %10.4f\n", GridLeftEdge[1], GridRightEdge[1]);
  printf("GridLeftEdge[2] = %10.4f  GridRightEdge[2] = %10.4f\n", GridLeftEdge[2], GridRightEdge[2]);

  if (fscanf(fptr, "Time = %"PSYM"\n", &Time) != 1) {
    fprintf(stderr, "Error reading Time.\n");
    return FAIL;
  }

  printf("Time %10.4f\n", Time);

  if (fscanf(fptr, "SubgridsAreStatic = %d\n", &SubgridsAreStatic) != 1) {
    fprintf(stderr, "Error reading SubgridsAreStatic.\n");
    return FAIL;
  }

  if (fscanf(fptr, "NumberOfBaryonFields = %d\n", &NumberOfBaryonFields) != 1) {
    fprintf(stderr, "Error reading NumberOfBaryonFields.\n");
    return FAIL;
  }

  printf("NumberOfBaryonFields = %d\n", NumberOfBaryonFields);

  if (NumberOfBaryonFields > 0) {

    fscanf(fptr, "FieldType = ");

    if (ReadListOfInts(fptr, NumberOfBaryonFields, FieldType) == FAIL) {
      fprintf(stderr, "Error reading FieldType.\n");
      return FAIL;
    }

    if (fscanf(fptr, "BaryonFileName = %s\n", BaryonFileName) != 1) {
      fprintf(stderr, "Error reading BaryonFileName.\n");
      return FAIL;
    }

    printf("BaryonFileName = %s\n", BaryonFileName);

    fscanf(fptr, "CourantSafetyNumber    = %"FSYM"\n", &CourantSafetyNumber);
    fscanf(fptr, "PPMFlatteningParameter = %d\n", &PPMFlatteningParameter);
    fscanf(fptr, "PPMDiffusionParameter  = %d\n", &PPMDiffusionParameter);
    fscanf(fptr, "PPMSteepeningParameter = %d\n", &PPMSteepeningParameter);

  }

  if (fscanf(fptr, "NumberOfParticles = %d\n", &NumberOfParticles) != 1) {
    fprintf(stderr, "error reading NumberOfParticles.\n");
    return FAIL;
  }

  printf("NumberOfParticles = %d\n", NumberOfParticles);

  if (NumberOfParticles > 0) {

    if (fscanf(fptr, "ParticleFileName = %s\n", ParticleFileName) != 1) {
      fprintf(stderr, "Error reading ParticleFileName.\n");
      return FAIL;
    }

    printf("ParticleFileName = %s\n", ParticleFileName);

  }

//  if (SelfGravity)
    if (fscanf(fptr, "GravityBoundaryType = %d\n",&GravityBoundaryType) != 1) {
      fprintf(stderr, "Error reading GravityBoundaryType.\n");
      return FAIL;
    }

    fscanf(fptr, "Pointer: Grid[%d]->NextGridThisLevel = %d\n", &ThisGrid, &NextGrid);

}

fclose(fptr);

}



// ListIO.C


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
