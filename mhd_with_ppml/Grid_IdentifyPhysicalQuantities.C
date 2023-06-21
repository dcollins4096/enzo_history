/***********************************************************************
/
/  GRID CLASS (IDENTIFY CERTAIN COMMONLY USED VARIABLES FROM THE LIST)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
#ifdef PPML
#include "PPML.h"
#endif //PPML
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
 
int grid::IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num,
				     int &Vel2Num, int &Vel3Num, int &TENum)
{
 
  DensNum = GENum = Vel1Num = Vel2Num = Vel3Num = TENum = 0;
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "GIPQ: Cannot find density.\n");
    return FAIL;
  }
 
  /* Find Total energy, if possible. */
 
#ifdef PPML
  if( EquationOfState == 0 )
    if ((TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find total energy.\n");
      return FAIL;
    }
#else
  if ((TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find total energy.\n");
    return FAIL;
  }
#endif //PPML 
  /* Find gas energy, if possible. */
 
  if (DualEnergyFormalism == TRUE)
    if ((GENum = FindField(InternalEnergy, FieldType,
			   NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find gas energy.\n");
      return FAIL;
    }
 
  /* Find Velocity1, if possible. */
 
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find Velocity1.\n");
    return FAIL;
  }
 
  /* Find Velocity2, if possible. */
 
  if (GridRank > 1)
    if ((Vel2Num = FindField(Velocity2, FieldType,
			     NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find Velocity2.\n");
      return FAIL;
    }
 
  /* Find Velocity3, if possible. */
 
  if (GridRank > 2)
    if ((Vel3Num = FindField(Velocity3, FieldType,
			     NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find Velocity3.\n");
      return FAIL;
    }
 
  return SUCCESS;
}
#ifdef PPML

//This is only better in that it's extendable to future physics packages.

int grid::IdentifyMHDQuantities(int &DensNum, int &GENum, int &Vel1Num,
				int &Vel2Num, int &Vel3Num, int &TENum,
				int & BX_Num, int &BY_Num, int & BZ_Num ){
  
  DensNum = GENum = Vel1Num = Vel2Num = Vel3Num = TENum = 0;
  
  /* Find Density, if possible. */
  
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Identify MHD QUantities: Cannot find density.\n");
    return FAIL;
  }
  
  /* Find Total energy, if possible. */
  if( EquationOfState == 0 ){
    if ((TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find total energy.\n");
      return FAIL;
    }
  }else{
    TENum = 0;
  }
  /* Find gas energy, if possible. */
  
  if (DualEnergyFormalism == TRUE)
    if ((GENum = FindField(InternalEnergy, FieldType,
			   NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find gas energy.\n");
      return FAIL;
    }
  
  /* Find Velocity1, if possible. */
  
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find Velocity1.\n");
    return FAIL;
  }
  
  /* Find Velocity2, if possible. */
  
  if (GridRank > 1
#ifdef PPML
      || MHD_Used
#endif //PPML
      )
    if ((Vel2Num = FindField(Velocity2, FieldType,
			     NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find Velocity2.\n");
      return FAIL;
    }
  
  /* Find Velocity3, if possible. */
  
  if (GridRank > 2
#ifdef PPML
      || MHD_Used
#endif //PPML
      )
    if ((Vel3Num = FindField(Velocity3, FieldType,
			     NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find Velocity3.\n");
      return FAIL;
    }

  if( MHD_Used == TRUE ){

    if ((BX_Num = FindField(Magnetic1, FieldType,
			    NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find BX\n");
      return FAIL;
    }

    if ((BY_Num = FindField(Magnetic2, FieldType,
			    NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find BY.\n");
      return FAIL;
    }

    if ((BZ_Num = FindField(Magnetic3, FieldType,
			    NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find BZ.\n");
      return FAIL;
    }
  }

  return SUCCESS;
}


int grid::IdentifyPhysicalQuantities_2(IndexPointerMap &ind){
  
  
  ind.D = ind.GE = ind.VX = ind.VY = ind.VZ = ind.TE = 0;
  ind.BX = ind.BY = ind.BZ;
  /* Find Density, if possible. */
  
  if ((ind.D = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "GIPQ_2: Cannot find density.\n");
    return FAIL;
  }
  
  if( EquationOfState == 0 ){
    /* Find Total energy, if possible. */
    
    if ((ind.TE = FindField(TotalEnergy, FieldType, NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find total energy.\n");
      return FAIL;
    }
  }
  /* Find gas energy, if possible. */
  
  if (DualEnergyFormalism == TRUE)
    if ((ind.GE = FindField(InternalEnergy, FieldType,
			    NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find gas energy.\n");
      return FAIL;
    }
  
  /* Find Velocity1, if possible. */
  
  if ((ind.VX= FindField(Velocity1, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find Velocity1.\n");
    return FAIL;
  }
  
  /* Find Velocity2, if possible. */
  
  if (GridRank > 1 || MHD_Used)
    if ((ind.VY = FindField(Velocity2, FieldType,
			    NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find Velocity2.\n");
      return FAIL;
    }
  
  /* Find Velocity3, if possible. */
  
  if (GridRank > 2 || MHD_Used)
    if ((ind.VZ = FindField(Velocity3, FieldType,
			    NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find Velocity3.\n");
      return FAIL;
    }
  
  
  if( MHD_Used == TRUE){
    if ((ind.BX = FindField(Magnetic1, FieldType,
			    NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find BX\n");
      return FAIL;
    }
    if ((ind.BY = FindField(Magnetic2, FieldType,
			    NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find BY.\n");
      return FAIL;
    }
    if ((ind.BZ = FindField(Magnetic3, FieldType,
			    NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find BZ.\n");
      return FAIL;
    }
  } 
    
  ind.V[0] = ind.VX;
  ind.V[1] = ind.VY;
  ind.V[2] = ind.VZ;
  ind.B[0] = ind.BX;
  ind.B[1] = ind.BY;
  ind.B[2] = ind.BZ;
#ifdef MHDF
  if( MHD_Used == TRUE )
    for(int field=0; field<3; field++)
      ind.CenteredB[field] = BaryonField[ ind.B[field] ];
#endif //MHDF
  return SUCCESS;
}
#endif //PPLM
