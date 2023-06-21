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
/  GRID CLASS (ADD THE COMOVING EXPANSION TERMS TO VELOCITY AND ENERGY)
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
#include "CosmologyParameters.h"

#define USE_FORTRAN

#ifdef HAOXU
#undef USE_FORTRAN
#endif 

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
extern "C" void FORTRAN_NAME(expand_terms)(
   int *rank, int *isize, int *idual, float *coef, int *imethod, float *gamma,
   float *p, float *pdual, float *d, float *e, float *ge, 
      float *u, float *v, float *w,
   float *dold, float *eold, float *geold, float *uold, float *vold, 
      float *wold);


int grid::ComovingExpansionTerms()
{
    
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  this->DebugCheck("ComovingExpansionTerms");

  if (NumberOfBaryonFields > 0) {

    /* Compute adot/a at time = t-1/2dt (time-centered). */

    FLOAT a, dadt;
    if (CosmologyComputeExpansionFactor(0.5*(Time+OldTime), &a, &dadt) 
	== FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
    float Coefficient = dtFixed*dadt/a;

    /* Determine the size of the grids. */

    int dim, size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

    /* Find the density, gas energy, velocities & total energy
       (where appropriate). */

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      return FAIL;
    }

#ifdef HAOXU
    float *PressureDual,*OldPressureDual, *Pressure = new float[size],*OldPressure = new float[size];
    int i;
#else
    float *PressureDual, *Pressure = new float[size];
#endif /* HAOXU */

#ifdef USE_FORTRAN

    /* Compute the time-centered pressure for this grid. */

    if (DualEnergyFormalism) {
      PressureDual = new float[size];
      if (this->ComputePressureDualEnergyFormalism(0.5*(Time+OldTime),
						   PressureDual) == FAIL) {
//      if (this->ComputePressureDualEnergyFormalism(Time,
//						   PressureDual) == FAIL) {
	fprintf(stderr, "Error in ComputePressureDualEnergyFormalism.\n");
	return FAIL;
      }
    }

    if (this->ComputePressure(0.5*(Time+OldTime), Pressure) == FAIL) {
//    if (this->ComputePressure(Time, Pressure) == FAIL) {
      fprintf(stderr, "Error in ComputePressure.\n");
      return FAIL;
    }

    /* Call fortran routine to do the real work. */

    FORTRAN_NAME(expand_terms)(
              &GridRank, &size, &DualEnergyFormalism, &Coefficient, 
	          (int*) &HydroMethod, &Gamma,
	      Pressure, PressureDual,
                  BaryonField[DensNum], BaryonField[TENum], 
                  BaryonField[GENum], BaryonField[Vel1Num], 
                  BaryonField[Vel2Num], BaryonField[Vel3Num],
              OldBaryonField[DensNum], OldBaryonField[TENum], 
                  OldBaryonField[GENum], OldBaryonField[Vel1Num], 
                  OldBaryonField[Vel2Num], OldBaryonField[Vel3Num]);

    if (DualEnergyFormalism)
      delete PressureDual;

#else /* USE_FORTRAN */

#ifdef HAOXU
	if(MHD_Used){

#define pressure_method_1
           
	if(DualEnergyFormalism){
			 PressureDual = new float[size];
			 OldPressureDual = new float[size];
#ifdef pressure_method_1_no
  if(this->ComputePressureDualEnergyFormalism(0.5*(Time+OldTime),
               PressureDual) == FAIL) {
              fprintf(stderr, "Error in ComputePressure.\n");
             return FAIL;
             }
#else
           if(this->ComputePressureDualEnergyFormalism(Time,
               PressureDual) == FAIL) {
              fprintf(stderr, "Error in ComputePressure.\n");
             return FAIL;
             }
	if (this->ComputePressureDualEnergyFormalism(OldTime,
              OldPressureDual) == FAIL) {
            fprintf(stderr, "Error in ComputePressure.\n");
            return FAIL;
         }
#endif
       }

#ifdef pressure_method_1_no
   if (this->ComputePressure(0.5*(Time+OldTime), Pressure) == FAIL) {
            fprintf(stderr, "Error in ComputePressure.\n");
            return FAIL;
        }
#else
	if (this->ComputePressure(Time, Pressure) == FAIL) {
            fprintf(stderr, "Error in ComputePressure.\n");
            return FAIL;
        }
	if (this->ComputePressure(OldTime, OldPressure) == FAIL) {
              fprintf(stderr, "Error in ComputePressure.\n");
             return FAIL;
         }
#endif

	for(i = 0; i<size; i++) {
	  if(DualEnergyFormalism) {
#ifdef pressure_method_1
	     PressureDual[i] *= (1.0 -
               0.5*3.0*(Gamma-1.0)*Coefficient)/(1.0 +
               0.5*3.0*(Gamma-1.0)*Coefficient);
       	     PressureDual[i] = max(PressureDual[i],tiny_pressure);
  	     BaryonField[GENum][i] = PressureDual[i]/(BaryonField[DensNum][i]*(Gamma-1.0));
#endif //pressure_method_1

#ifdef pressure_method_2
     PressureDual[i] -= 2*Coefficient*PressureDual[i];
     BaryonField[GENum][i] = PressureDual[i]/(BaryonField[DensNum][i]*(Gamma-1.0));
#endif //pressure_method_2
           } //DualEnergyFormalism

#ifdef pressure_method_1
	Pressure[i] *= (1.0 -
                   0.5*3.0*(Gamma-1.0)*Coefficient)/(1.0  +
                   0.5*3.0*(Gamma-1.0)*Coefficient);
        Pressure[i] = max(Pressure[i],tiny_pressure);
#endif

#ifdef pressure_method_2
        Pressure[i] -= 2*Coefficient*Pressure[i];
#endif 

#define Velocity_Method_1

#ifdef Velocity_Method_1  //Self_implicit  
        BaryonField[Vel1Num][i] *= (1.0 -
            0.5*Coefficient)/(1.0+0.5*Coefficient);
	BaryonField[Vel2Num][i] *= (1.0 -
           0.5*Coefficient)/(1.0+0.5*Coefficient);
        BaryonField[Vel3Num][i] *= (1.0 -
           0.5*Coefficient)/(1.0+0.5*Coefficient);
#endif

#ifdef Velocity_Method_2 
         BaryonField[Vel1Num][i] -= Coefficient*0.5*(   BaryonField[Vel1Num][i]
           + OldBaryonField[Vel1Num][i]);
        BaryonField[Vel2Num][i] -= Coefficient*0.5*(   BaryonField[Vel2Num][i]
           + OldBaryonField[Vel2Num][i]);
        BaryonField[Vel3Num][i] -= Coefficient*0.5*(   BaryonField[Vel3Num][i]
            +       OldBaryonField[Vel3Num][i]);

#endif

// Expansion of Magnetic Fields
   if(MHD_Equation == 2) {
   CenteredB[0][i] *= (1.0 - 0.25*Coefficient)/(1.0+0.25*Coefficient);
   CenteredB[1][i] *= (1.0 - 0.25*Coefficient)/(1.0+0.25*Coefficient);
   CenteredB[2][i] *= (1.0 - 0.25*Coefficient)/(1.0+0.25*Coefficient);
   }

 if(DualEnergyFormalism) {
      BaryonField[TENum][i] = PressureDual[i]/(Gamma-1.0) +
0.5*BaryonField[DensNum][i]*
             (BaryonField[Vel1Num][i]*BaryonField[Vel1Num][i]
              +BaryonField[Vel2Num][i]*BaryonField[Vel2Num][i]
             +BaryonField[Vel3Num][i]*BaryonField[Vel3Num][i])
             + 0.5*(CenteredB[0][i]*CenteredB[0][i]+CenteredB[1][i]
             *CenteredB[1][i]+CenteredB[2][i]*CenteredB[2][i]);
 }else{
	BaryonField[TENum][i] = Pressure[i]/(Gamma-1.0) +
0.5*BaryonField[DensNum][i]*
             (BaryonField[Vel1Num][i]*BaryonField[Vel1Num][i]
              +BaryonField[Vel2Num][i]*BaryonField[Vel2Num][i]
	     +BaryonField[Vel3Num][i]*BaryonField[Vel3Num][i])
             + 0.5*(CenteredB[0][i]*CenteredB[0][i]+CenteredB[1][i]
             *CenteredB[1][i]+CenteredB[2][i]*CenteredB[2][i]);
 }//DualEnergyFormalism

/*
 // Compute Expansion of Energy directly
        BaryonField[TENum][i]-=0.5*(CenteredB[0][i]*CenteredB[0][i]+CenteredB[1][i]
             *CenteredB[1][i]+CenteredB[2][i]*CenteredB[2][i]);

        BaryonField[TENum][i]*=(1.0 -
               0.5*3.0*(Gamma-1.0)*Coefficient)/(1.0 +
             0.5*3.0*(Gamma-1.0)*Coefficient);

BaryonField[TENum][i]+=0.5*(CenteredB[0][i]*CenteredB[0][i]+CenteredB[1][i]
*CenteredB[1][i]+CenteredB[2][i]*CenteredB[2][i]);  
*/
	}
	if(DualEnergyFormalism) {
		delete PressureDual;
		delete OldPressureDual;
	}
  if ( MHD_Equation == 2){
     for(int field=0;field<3;field++)
     for(i=0;i<MagneticSize[field];i++)
     MagneticField[field][i]
      *=(1.0-0.25*Coefficient)/(1.0+0.25*Coefficient);
  }

	}else{  //MHD_used

    /* Compute the time-centered pressure for this grid. */
    if (this->ComputePressure(Time, Pressure) == FAIL) {
      fprintf(stderr, "Error in ComputePressure.\n");
      return FAIL;
    }

    /* Replace pressure with the combination 3*p/d. */
    for ( i = 0; i < size; i++)
      Pressure[i] = 3.0*Pressure[i] / BaryonField[DensNum][i];

    /* Apply the expansion terms using time-centered quantities
       (BaryonFields are at t, OldBaryonFields at t-dt).  We apply
       the expansion term over the entire grid because we can and
       its faster. */

    /* (B) Total energy (first add v^2 to sum in pressure field). */
    if (HydroMethod != Zeus_Hydro) {
      for (i = 0; i < size; i++)
	Pressure[i] += BaryonField[Vel1Num][i]*BaryonField[Vel1Num][i];
      
      if (Vel2Num != 0)
	for (i = 0; i < size; i++)
	  Pressure[i] += BaryonField[Vel2Num][i]*BaryonField[Vel2Num][i];
      
      if (Vel3Num != 0)
	for (i = 0; i < size; i++)
	  Pressure[i] += BaryonField[Vel3Num][i]*BaryonField[Vel3Num][i];
    }

    for (i = 0; i < size; i++)
      BaryonField[TENum][i] -= Coefficient*Pressure[i];


    /* (A) velocity terms. */

#define VELOCITY_METHOD2

#ifdef VELOCITY_METHOD1

    /*    i) semi-implicit way: */

    for (i = 0; i < size; i++)
      BaryonField[Vel1Num][i] *= 
	(1.0 - 0.5*Coefficient)/(1.0 + 0.5*Coefficient);

    if (Vel2Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel2Num][i] *= 
	  (1.0 - 0.5*Coefficient)/(1.0 + 0.5*Coefficient);

    if (Vel3Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel3Num][i] *= 
	  (1.0 - 0.5*Coefficient)/(1.0 + 0.5*Coefficient);

#endif /* VELOCITY_METHOD1 */

#ifdef VELOCITY_METHOD2

    /*    ii) old way: */

    for (i = 0; i < size; i++)
      BaryonField[Vel1Num][i] -= Coefficient*0.5*(   BaryonField[Vel1Num][i]
           + OldBaryonField[Vel1Num][i]);
    if (Vel2Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel2Num][i] -= Coefficient*0.5*(   BaryonField[Vel2Num][i] 
           + OldBaryonField[Vel2Num][i]);
    if (Vel3Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel3Num][i] -= Coefficient*0.5*(   BaryonField[Vel3Num][i]
            + 	    OldBaryonField[Vel3Num][i]);

#endif /* VELOCITY_METHOD2 */

#ifdef VELOCITY_METHOD3

    /*    iii) dumb way: */

    for (i = 0; i < size; i++)
      BaryonField[Vel1Num][i] -= Coefficient*OldBaryonField[Vel1Num][i];
    if (Vel2Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel2Num][i] -= 
	                           Coefficient*OldBaryonField[Vel2Num][i];
    if (Vel3Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel3Num][i] -= 
			           Coefficient*OldBaryonField[Vel3Num][i];

#endif /* VELOCITY_METHOD3 */

    /* (C) If we're using the dual energy formalism, we must recompute the
       pressure from the gas energy (for consistency). */

    if (DualEnergyFormalism) {

      if (this->ComputePressureDualEnergyFormalism(0.5*(Time+OldTime),
						   Pressure) == FAIL) {
	fprintf(stderr, "Error in ComputePressureDualEnergyFormalism.\n");
	return FAIL;
      }

      /* Replace pressure with the time-centered combination 3*p/d. */

      for (i = 0; i < size; i++)
	Pressure[i] = 6.0*Pressure[i] / 
              (BaryonField[DensNum][i] + OldBaryonField[DensNum][i]);

      /* Add the expansion terms to the gas energy. */

      for (i = 0; i < size; i++)
	BaryonField[GENum][i] -= Coefficient*Pressure[i];
    }
}
#else /* HAOXU */

    /* Compute the time-centered pressure for this grid. */

    if (this->ComputePressure(Time, Pressure) == FAIL) {
      fprintf(stderr, "Error in ComputePressure.\n");
      return FAIL;
    }

    /* Replace pressure with the combination 3*p/d. */

    for (int i = 0; i < size; i++)
      Pressure[i] = 3.0*Pressure[i] / BaryonField[DensNum][i];

    /* Apply the expansion terms using time-centered quantities
       (BaryonFields are at t, OldBaryonFields at t-dt).  We apply
       the expansion term over the entire grid because we can and
       its faster. */

    /* (B) Total energy (first add v^2 to sum in pressure field). */

    if (HydroMethod != Zeus_Hydro) {
      for (i = 0; i < size; i++)
	Pressure[i] += BaryonField[Vel1Num][i]*BaryonField[Vel1Num][i];
      
      if (Vel2Num != 0)
	for (i = 0; i < size; i++)
	  Pressure[i] += BaryonField[Vel2Num][i]*BaryonField[Vel2Num][i];
      
      if (Vel3Num != 0)
	for (i = 0; i < size; i++)
	  Pressure[i] += BaryonField[Vel3Num][i]*BaryonField[Vel3Num][i];
    }

    for (i = 0; i < size; i++)
      BaryonField[TENum][i] -= Coefficient*Pressure[i];

    /* (A) velocity terms. */

#define VELOCITY_METHOD2

#ifdef VELOCITY_METHOD1

    /*    i) semi-implicit way: */

    for (i = 0; i < size; i++)
      BaryonField[Vel1Num][i] *= 
	(1.0 - 0.5*Coefficient)/(1.0 + 0.5*Coefficient);

    if (Vel2Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel2Num][i] *= 
	  (1.0 - 0.5*Coefficient)/(1.0 + 0.5*Coefficient);

    if (Vel3Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel3Num][i] *= 
	  (1.0 - 0.5*Coefficient)/(1.0 + 0.5*Coefficient);

#endif /* VELOCITY_METHOD1 */

#ifdef VELOCITY_METHOD2

    /*    ii) old way: */

    for (i = 0; i < size; i++)
      BaryonField[Vel1Num][i] -= 
	                       Coefficient*0.5*(   BaryonField[Vel1Num][i] + 
						OldBaryonField[Vel1Num][i] );
    if (Vel2Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel2Num][i] -= 
			       Coefficient*0.5*(   BaryonField[Vel2Num][i] +
				                OldBaryonField[Vel2Num][i] );
    if (Vel3Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel3Num][i] -= 
			       Coefficient*0.5*(   BaryonField[Vel3Num][i] + 
					        OldBaryonField[Vel3Num][i] );
#endif /* VELOCITY_METHOD2 */

#ifdef VELOCITY_METHOD3

    /*    iii) dumb way: */

    for (i = 0; i < size; i++)
      BaryonField[Vel1Num][i] -= Coefficient*OldBaryonField[Vel1Num][i];
    if (Vel2Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel2Num][i] -= 
	                           Coefficient*OldBaryonField[Vel2Num][i];
    if (Vel3Num != 0)
      for (i = 0; i < size; i++)
	BaryonField[Vel3Num][i] -= 
			           Coefficient*OldBaryonField[Vel3Num][i];

#endif /* VELOCITY_METHOD3 */

    /* (C) If we're using the dual energy formalism, we must recompute the
       pressure from the gas energy (for consistency). */

    if (DualEnergyFormalism) {

      if (this->ComputePressureDualEnergyFormalism(0.5*(Time+OldTime),
						   Pressure) == FAIL) {
	fprintf(stderr, "Error in ComputePressureDualEnergyFormalism.\n");
	return FAIL;
      }

      /* Replace pressure with the time-centered combination 3*p/d. */

      for (i = 0; i < size; i++)
	Pressure[i] = 6.0*Pressure[i] / 
              (BaryonField[DensNum][i] + OldBaryonField[DensNum][i]);

      /* Add the expansion terms to the gas energy. */

      for (i = 0; i < size; i++)
	BaryonField[GENum][i] -= Coefficient*Pressure[i];

    }
#endif /* HAOXU */

#endif /* USE_FORTRAN */

    /* clean up */
#ifdef HAOXU
    delete Pressure;
    delete OldPressure;
#else
    delete Pressure;
#endif /* HAOXU */

  }

  return SUCCESS;
}
