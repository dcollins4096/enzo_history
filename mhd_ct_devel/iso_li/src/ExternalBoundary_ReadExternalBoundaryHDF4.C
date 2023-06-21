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
/  EXTERNAL BOUNDARY CLASS (READ THE EXTERNAL BOUNDARY VALUES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  James Bordner, June 2003   added USE_HDF4
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_HDF4

#include <string.h>
#include <stdio.h>
#include <df.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


// This routine reads the external boundary from the provided file pointer
//

/* function prototypes */

int ReadListOfInts(FILE *fptr, int N, int nums[]);
int ReadListOfFloats(FILE *fptr, int N, float floats[]);


int ExternalBoundary::ReadExternalBoundaryHDF4(FILE *fptr)
{

  /* declarations */

  int Dims[MAX_DIMENSION], index, size, i;
  int BoundaryValuePresent[2*MAX_DIMENSION];
  int field, TempInt, j;
  int32 TempIntArray[MAX_DIMENSION];
  float32 *buffer;
  char hdfname[MAX_LINE_LENGTH];

  /* read general class data */

  if (fscanf(fptr, "BoundaryRank = %d\n", &BoundaryRank) != 1) {
    fprintf(stderr, "Error reading BoundaryRank.\n");
    return FAIL;
  }

  fscanf(fptr, "BoundaryDimension =");
  if (ReadListOfInts(fptr, BoundaryRank, BoundaryDimension) == FAIL) {
    fprintf(stderr, "Error reading BoundaryDimension.\n");
    return FAIL;
  }

  /* read baryon field quantities */

  if (fscanf(fptr, "NumberOfBaryonFields = %d\n", 
	     &NumberOfBaryonFields) != 1) {
    fprintf(stderr, "Error reading NumberOfBaryonFields.\n");
    return FAIL;
  }

  /* Read particle boundary type. */

  if (fscanf(fptr, "ParticleBoundaryType = %d\n",&ParticleBoundaryType) != 1) {
    fprintf(stderr, "Error reading ParticleBoundaryType.\n");
    return FAIL;
  }

  if (NumberOfBaryonFields > 0) {

    /* read field types */

    fscanf(fptr, "BoundaryFieldType = ");
    if (ReadListOfInts(fptr, NumberOfBaryonFields, BoundaryFieldType) 
        == FAIL) {
      fprintf(stderr, "Error reading BoundaryFieldType.\n");
      return FAIL;
    }

    /* read hdf file name */

    if (fscanf(fptr, "BaryonFileName = %s\n", hdfname) != 1) {
      fprintf(stderr, "Error reading BaryonFileName.\n");
      return FAIL;
    }    

    /* read BoundaryValue present line */

    fscanf(fptr, "BoundaryValuePresent = ");
    if (ReadListOfInts(fptr, BoundaryRank*2, BoundaryValuePresent) == FAIL) {
      fprintf(stderr, "Error reading BoundaryValuePresent.\n");
      return FAIL;
    }

    /* loop over faces, reading each */

    for (int dim = 0; dim < BoundaryRank; dim++)
      if (BoundaryDimension[dim] > 1) {

	/* calculate size and dims of flux plane */
	
	index = 0;
	size  = 1;
	Dims[0] = 1;
	for (i = 0; i < BoundaryRank; i++)
	  if (i != dim) {
	    Dims[index++] = BoundaryDimension[i];
	    size *= BoundaryDimension[i];
	  }
	index = max(BoundaryRank-1, 1);   // make index at least 1

	/* Read HDF dims */

	if (DFSDgetdims(hdfname, &TempInt, TempIntArray, BoundaryRank) 
	    == HDF_FAIL) {
	  fprintf(stderr, "Error in DFSDgetdims.\n");
	  return FAIL;
	}

	/* Check rank and dimensions (dims are stored backwards for us). */

	if (TempInt != index) {
	  fprintf(stderr, "HDF file rank does not match BoundaryRank.\n");
	  return FAIL;
	}

	for (i = 0; i < index; i++)
	  if (TempIntArray[index-i-1] != Dims[i]) {
	    fprintf(stderr, "HDF file dims do not match BoundaryDims.\n");
	    fprintf(stderr, " Dims[%d] = %d   HDF Dims[%d] = %d\n", i, Dims[i],
		    index-i-1, TempIntArray[index-i-1]);
	    return FAIL;
	  }

	/* Allocate temporary space. */
	
	buffer = new float32[size];

	/* loop over fields, reading each */

	for (field = 0; field < NumberOfBaryonFields; field++)
	  for (i = 0; i < 2; i++) {

	    /* allocate room for BoundaryType */

	    BoundaryType[field][dim][i] = new boundary_type[size];

	    /* read BoundaryType (then convert to int) */

	    if (DFSDgetdata(hdfname, TempInt, TempIntArray, (VOIDP) buffer)
		== HDF_FAIL) {
	      fprintf(stderr, "Error in DFSDgetdata(0).\n");
	      return FAIL;
	    }

//RH
//          printf("REB read BoundaryType\n");
//          printf("REB getdata from %s\n",hdfname);
//          printf("REB getdata rank %d\n",TempInt);
//          for (i = 0; i < TempInt; i++)
//          {
//            printf("%d  %d\n",i,TempIntArray[i]);
//          }
//RH

	    for (j = 0; j < size; j++)
	      BoundaryType[field][dim][i][j] = 
		                      (boundary_type) nint(buffer[j]);

	    /* read BoundaryValue */

	    if (BoundaryValuePresent[2*dim+i]) {
	      BoundaryValue[field][dim][i] = new float[size];
	      if (DFSDgetdata(hdfname, TempInt, TempIntArray, (VOIDP)
			      buffer) == HDF_FAIL) {
		fprintf(stderr, "Error in DFSDgetdata(1).\n");
		fprintf(stderr, "dim = %d field = %d i = %d\n", dim, field, i);
		return FAIL;
	      }

//RH
//          printf("REB read BoundaryValue\n");
//          printf("REB getdata from %s\n",hdfname);
//          printf("REB getdata rank %d\n",TempInt);
//          for (i = 0; i < TempInt; i++)
//          {
//            printf("%d  %d\n",i,TempIntArray[i]);
//          }
//RH

	      for (j = 0; j < size; j++)
		BoundaryValue[field][dim][i][j] = float(buffer[j]);
	    }
	  }  // end of loop over fields

	delete buffer;
 	
      }   // end of loop over dims

  }

  return SUCCESS;

}

#else

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "ExternalBoundary.h"
#include "message.h"

// HDF4 is not used, so HDF4ReadExternalBoundary should not be called

int ExternalBoundary::ReadExternalBoundaryHDF4(FILE *fptr)
{ 
  ERROR_MESSAGE;
  return FAIL;
}
#endif
