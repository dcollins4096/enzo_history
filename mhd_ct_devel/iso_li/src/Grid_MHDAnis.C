
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

// This looks for anisotropies in test problems.  Input is a string of characters
// that will be put to the screen when the check fails.
// Only usefull for problems that have both x and y translation symmetry.
// dcc sept. 04

int grid::MHDAnis(char * string){

  //fprintf(stderr,"kludge: no anis\n");
  return SUCCESS;
  if( !MHD_Used)
    return SUCCESS;
  if( ProcessorNumber != MyProcessorNumber ) 
    return SUCCESS;

  float tolerance = 1e-10;
  int anis = FALSE, i,j,k,field;

  for( field=0;field<NumberOfBaryonFields;field++)
    for( i=GridStartIndex[0];i<=GridEndIndex[0];i++)
      for( j=GridStartIndex[1];j<=GridEndIndex[1];j++)
	for( k=GridStartIndex[2];k<=GridEndIndex[2];k++){

	  if( fabs(BaryonField[field][index0(i,j,k)] 
	      -BaryonField[field][index0(GridStartIndex[0], GridStartIndex[1],k)])
	      > tolerance)  {
	    anis = TRUE;
	    goto Anis;
	  }//anis
	}//k

  if( anis == FALSE )
    for( field=0;field<3;field++)
      for( i=MHDStartIndex[field][0];i<=MHDEndIndex[field][0];i++)
	for( j=MHDStartIndex[field][1];j<=MHDEndIndex[field][1];j++)
	  for( k=MHDStartIndex[field][2];k<=MHDEndIndex[field][2];k++){

	    if( fabs(MagneticField[field][indexba(i,j,k,field)] 
		     -MagneticField[field][indexba(MHDStartIndex[field][0], 
						   MHDStartIndex[field][1],k,field)])
		> tolerance)  {

	      fprintf(stderr,"Magnetic %d %d %d %d",field,i,j,k );
	      fprintf(stderr," %f %f (The fuck?)\n",MagneticField[field][indexba(i,j,k,field)],
		      MagneticField[field][indexba(MHDStartIndex[field][0],
						   MHDStartIndex[field][1],k,field)]);
	      
	      anis = TRUE;
	      goto Anis;
	    }//anis
	  }//k

 Anis:

  if( anis == TRUE ){
    fprintf(stderr, "*************** AMR ANISOTROPY!!! ******\n");
    fprintf(stderr, "*************** %d %d %d %d ***************\n",field, i,j,k);
    fprintf(stderr, " %s \n", string);
    fprintf(stderr, " Baryons 3 3 k = %f , i j k = %f, diff = %f\n", 
	    BaryonField[field][index0(i,j,k)],
	    BaryonField[field][index0(GridStartIndex[0], GridStartIndex[1],k)],
	    BaryonField[field][index0(i,j,k)]-
	    BaryonField[field][index0(GridStartIndex[0], GridStartIndex[1],k)]);
    fprintf(stderr, " Mag 3 3 k = %f , i j k = %f, diff = %f\n", 
	    MagneticField[field][indexba(i,j,k,field)],
	    MagneticField[field][indexba(GridStartIndex[0], GridStartIndex[1],k,field)],
	    MagneticField[field][indexba(i,j,k,field)]-
	    MagneticField[field][indexba(GridStartIndex[0], GridStartIndex[1],k,field)]);
    fprintf(stderr, " Start %d %d %d \n", GridStartIndex[0], GridStartIndex[1], 
	    GridStartIndex[2]);
    fprintf(stderr, " End %d %d %d \n", GridEndIndex[0], GridEndIndex[1], 
	    GridEndIndex[2]);
    fprintf(stderr, " Dims %d %d %d \n", GridDimension[0], GridDimension[1], 
	    GridDimension[2]);
    
    // dump data
    if(1==0){
      char basename[30];
      sprintf(basename, "data0666.grid");
      FILE *failptr = fopen(basename,"a");
      
      WriteGrid(failptr, basename, 666);
    }
    return FAIL;
  }
  else
    return SUCCESS;
}
