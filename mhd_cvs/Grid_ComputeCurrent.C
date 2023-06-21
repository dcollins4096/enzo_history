

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

#ifdef PASSIVE

int grid::MHD_ComputeCurrent(void)
{

  if( MyProcessorNumber != ProcessorNumber )
    return SUCCESS;


   int Crash = FALSE;

  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  float TotalEnergy = 0;
  float Kinetic = 0;
  float MagneticCentered = 0;
  float MagneticFace = 0;
  float Mass = 0;
  float GasEnergy = 0;
  float AbsDivB = 0;
  float TotalDivB = 0;
  float dens, v1, v2, v3, eng, b1, b2, b3, divergence;
  float MaxDivB = -1;

  int i,j,k, index;

  int max[3]={-1,-1,-1}, min[3] = {GridDimension[0],GridDimension[1],GridDimension[2]};
  int Convex = 1;
  float MaxTolerance = 1e-8;

  if( DivB == NULL) {
    DivB = new float[size];
    for(i=0;i<size;i++) DivB[i] = 0.0;
    
  }

    for(i=0;i<3;i++)
      if( Current[i] == NULL ){
	Current[i] = new float[size];
	for(j=0;j<size;j++) Current[i][j] = 0.0;
     }

  

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
  
  
	
	
  //Compute current


 
  for(k=GridStartIndex[2]-2;k<=GridEndIndex[2]+2; k++)
    for(j=GridStartIndex[1]-2;j<=GridEndIndex[1]+2;j++)
      for(i=GridStartIndex[0]-2;i<=GridEndIndex[0]+2;i++) {
	int tmp1= indexb2(i,j,k-1);
	Current[0][index0(i,j,k)] = (1/CellWidth[1][0]* 
				     ( MagneticField[2][indexb3(i,j,k)]-MagneticField[2][indexb3(i,j-1,k)] )
				     -1/CellWidth[2][0]*
                                     ( MagneticField[1][indexb2(i,j,k)]-MagneticField[1][tmp1] ));
	//indexb2(i,j,k-1)] ));
	
	Current[1][index0(i,j,k)] = (1/CellWidth[2][0]* 
				     ( MagneticField[0][indexb1(i,j,k)]-MagneticField[0][indexb1(i,j,k-1)] )
				     -1/CellWidth[0][0]* 
				     ( MagneticField[2][indexb3(i,j,k)]-MagneticField[2][indexb3(i-1,j,k)] ));
	
	Current[2][index0(i,j,k)] = (1/CellWidth[0][0]* 
				     ( MagneticField[1][indexb2(i,j,k)]-MagneticField[1][indexb2(i-1,j,k)] )
				     -1/CellWidth[1][0]*
				     ( MagneticField[0][indexb1(i,j,k)]-MagneticField[0][indexb1(i,j-1,k)] ));
	
	//need to correct for the fact that Magneticfield isn't trustworthy in the ghost zones
	if(0==1){
	  if( k == GridStartIndex[2] || k == GridEndIndex[2] ){
	    Current[0][index0(i,j,k)] = 0.0;
	    Current[1][index0(i,j,k)] = 0.0;
	  }
	  
	  if( j == GridStartIndex[1] || j == GridEndIndex[1] ) {
	    Current[0][index0(i,j,k)] = 0.0;
	    Current[2][index0(i,j,k)] = 0.0;
	  }
	  
	  
	  if( i == GridStartIndex[0] || i == GridEndIndex[0] ) {
	    Current[1][index0(i,j,k)] = 0.0;
	    Current[2][index0(i,j,k)] = 0.0;
	  }
	}
	
    } // end current  
  

  return SUCCESS;

}

#endif //PASSIVE
