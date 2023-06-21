/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (WRITE THE EXTERNAL BOUNDARY VALUES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

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

// This routine writes the external boundary to the provided file pointer
//

/* function prototypes */

void WriteListOfInts(FILE *fptr, int N, int nums[]);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);


int ExternalBoundary::WriteExternalBoundary(FILE *fptr, char *hdfname)
{

  /* declarations */

  int dim, field, i, j, index, ret, size;
  int BoundaryValuePresent[MAX_DIMENSION*2], Temp[MAX_DIMENSION];
  int32 OutDims[MAX_DIMENSION];
  float32 *buffer;

  /* Save general class data */

  fprintf(fptr, "BoundaryRank         = %d\n", BoundaryRank);

  fprintf(fptr, "BoundaryDimension    = ");
  WriteListOfInts(fptr, BoundaryRank, BoundaryDimension);

  /* save baryon field quantities */

  fprintf(fptr, "NumberOfBaryonFields = %d\n", NumberOfBaryonFields);

  /* Save particle boundary type. */

  fprintf(fptr, "ParticleBoundaryType = %d\n", ParticleBoundaryType);

  if (NumberOfBaryonFields > 0) {

    fprintf(fptr, "BoundaryFieldType    = ");
    WriteListOfInts(fptr, NumberOfBaryonFields, BoundaryFieldType);

    fprintf(fptr, "BaryonFileName       = %s\n", hdfname);

    /* write out information about the BoundaryValue fields. */

    for (dim = 0; dim < BoundaryRank; dim++)
      for (i = 0; i < 2; i++)  
	if (BoundaryValue[0][dim][i] == NULL)
	  BoundaryValuePresent[2*dim+i] = FALSE;
	else
	  BoundaryValuePresent[2*dim+i] = TRUE;
    fprintf(fptr, "BoundaryValuePresent = ");
    WriteListOfInts(fptr, BoundaryRank*2, BoundaryValuePresent);

    /* Write HDF files */

//RH
//  printf("WriteEB start\n");
//  printf("  NumberOfBaryonFields %d\n",NumberOfBaryonFields);
//  printf("  BoundaryRank %d\n", BoundaryRank);
//  for (dim = 0; dim < BoundaryRank; dim++)
//  {
//     printf("    BoundaryDimension[%d] %d\n",dim,BoundaryDimension[dim]);
//  }
//RH

    for (dim = 0; dim < BoundaryRank; dim++)
      if (BoundaryDimension[dim] > 1) {

	/* calculate size and dims of flux plane */
	
	index   = 0;
	size    = 1;
	Temp[0] = 1;
	for (i = 0; i < BoundaryRank; i++)
	  if (i != dim) {
	    Temp[index++] = BoundaryDimension[i];
	    size *= BoundaryDimension[i];
	  }
	index = max(BoundaryRank-1, 1);   // make index at least 1

	/* Reverse outdims (for HDF). */

	for (i = 0; i < index; i++)
	  OutDims[index-i-1] = Temp[i];

//RH
//      printf("    Index %d\n",index);
//      for (i = 0; i < index; i++)
//      {
//        printf("      OutDims[%d] %d\n",i,OutDims[i]);
//      }
//      printf("    Size %d\n",size);
//RH

	/* Allocate temporary space. */
	
	buffer = new float32[size];

	for (field = 0; field < NumberOfBaryonFields; field++)
	  for (i = 0; i < 2; i++) {

//RH
//        printf("        dim %d : field %d : i %d\n",dim,field,i);
//RH

            /* Set dimensions. */

	    if (DFSDsetdims(index, OutDims) == HDF_FAIL) {
	      fprintf(stderr, "Error in DFSDsetdims(EB).\n");
	      return FAIL;
	    }

	    /* write out BoundaryType (convert to float first) */

	    for (j = 0; j < size; j++)
	      buffer[j] = float32(BoundaryType[field][dim][i][j]);

	    if (field == 0 && i == 0 && dim == 0)
	      ret = DFSDputdata(hdfname, index, OutDims, (VOIDP) buffer);
	    else
	      ret = DFSDadddata(hdfname, index, OutDims, (VOIDP) buffer);
	    if (ret == HDF_FAIL) {
	      fprintf(stderr, "Error in DFSDputdata(1).\n");
	      return FAIL;
	    }

	    /* write out BoundaryValue */

            if (BoundaryValue[field][dim][i] != NULL) {
	      for (j = 0; j < size; j++)
		buffer[j] = float32(BoundaryValue[field][dim][i][j]);
	      if (DFSDadddata(hdfname, index, OutDims, (VOIDP) buffer)
		  == HDF_FAIL) {
		fprintf(stderr, "Error in DFSDputdata(2).\n");
		return FAIL;
	      }
	    }

	  }  // end of loop over fields

	delete buffer;

      }  // end of loop over dims

  }

  return SUCCESS;

}
