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
#include "IndexPointerMap.h"
/* function prototypes */

int FindField(int f, int farray[], int n);



int grid::IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, 
				     int &Vel2Num, int &Vel3Num, int &TENum)
{



  DensNum = GENum = Vel1Num = Vel2Num = Vel3Num = TENum = 0;
    
  /* Find Density, if possible. */

  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find density.\n");
    return FAIL;
  }

  /* Find Total energy, if possible. */

  if ((TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields)) < 0) {
#ifdef ATHENA
    if( EquationOfState != 1){
#endif //ATHENA       
    fprintf(stderr, "Cannot find total energy.\n");
    return FAIL;
#ifdef ATHENA
    }
#endif //ATHENA
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
  // dcc June 3 2006. MHD always needs all 3 velocities, even for 1d problems. 
  if (GridRank > 1 || MHD_Used == TRUE)
    if ((Vel2Num = FindField(Velocity2, FieldType, 
			     NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find Velocity2.\n");
      return FAIL;
    }

  /* Find Velocity3, if possible. */

  if (GridRank > 2 || MHD_Used == TRUE)
    if ((Vel3Num = FindField(Velocity3, FieldType, 
			     NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find Velocity3.\n");
      return FAIL;
    }

  return SUCCESS;
}


//New version, that uses a structure
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

  ind.V[0] = ind.VX;
  ind.V[1] = ind.VY;
  ind.V[2] = ind.VZ;
  //ind.B[0] = ind.BX;
  //ind.B[1] = ind.BY;
  //ind.B[2] = ind.BZ;

  
  /* for now, the magnetic field isn't stored in the Field Type.  
     Code stolen from the version that does (PPML)
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
  */

  return SUCCESS;
}
