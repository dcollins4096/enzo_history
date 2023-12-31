/***********************************************************************
/
/  GRID CLASS (COMPUTE TIME STEP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    dt   - timestep
/
************************************************************************/
 
// Compute the timestep from all the constrains for this grid.
//
// Somebody fix the error handling in this routine! please.
//
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "StarParticleData.h"
#include "RadiativeTransferParameters.h"
#include "hydro_rk/EOS.h"
#include "hydro_rk/tools.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionTimestep(FLOAT time, float *dtExpansion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
extern "C" void FORTRAN_NAME(calc_dt)(
                  int *rank, int *idim, int *jdim, int *kdim,
                  int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
			     hydro_method *ihydro, float *C2,
                  FLOAT *dx, FLOAT *dy, FLOAT *dz, float *vgx, float *vgy,
                             float *vgz, float *gamma, int *ipfree, float *aye,
                  float *d, float *p, float *u, float *v, float *w,
			     float *dt, float *dtviscous);
 
 
float grid::ComputeTimeStep()
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return huge_number;
 
  this->DebugCheck("ComputeTimeStep");
 
  /* initialize */
 
  float dt, dtTemp;
  float dtBaryons      = huge_number;
  float dtViscous      = huge_number;
  float dtParticles    = huge_number;
  float dtExpansion    = huge_number;
  float dtAcceleration = huge_number;
  float dtMHD          = huge_number;
  int dim, i, result;
 
  /* Compute the field size. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
 
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  float afloat = float(a);
 
  /* 1) Compute Courant condition for baryons. */
 
  if (NumberOfBaryonFields > 0 && (HydroMethod != HD_RK) && (HydroMethod != MHD_RK)) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
      exit(FAIL);
    }
 
    /* Compute the pressure. */
 
    float *pressure_field = new float[size];
    if (DualEnergyFormalism)
      result = this->ComputePressureDualEnergyFormalism(Time, pressure_field);
    else
      result = this->ComputePressure(Time, pressure_field);
 
    if (result == FAIL) {
      fprintf(stderr, "Error in grid->ComputePressure.\n");
      exit(EXIT_FAILURE);
    }
 
#ifdef UNUSED
    int Zero[3] = {0,0,0}, TempInt[3] = {0,0,0};
    for (dim = 0; dim < GridRank; dim++)
      TempInt[dim] = GridDimension[dim]-1;
#endif /* UNUSED */
 
    /* Call fortran routine to do calculation. */
 
    FORTRAN_NAME(calc_dt)(&GridRank, GridDimension, GridDimension+1,
                               GridDimension+2,
//                        Zero, TempInt, Zero+1, TempInt+1, Zero+2, TempInt+2,
                          GridStartIndex, GridEndIndex,
                               GridStartIndex+1, GridEndIndex+1,
                               GridStartIndex+2, GridEndIndex+2,
			       &HydroMethod, &ZEUSQuadraticArtificialViscosity,
                          CellWidth[0], CellWidth[1], CellWidth[2],
                               GridVelocity, GridVelocity+1, GridVelocity+2,
                               &Gamma, &PressureFree, &afloat,
                          BaryonField[DensNum], pressure_field,
                               BaryonField[Vel1Num], BaryonField[Vel2Num],
                               BaryonField[Vel3Num], &dtBaryons, &dtViscous);
 
    /* Clean up */
 
    delete [] pressure_field;
 
    /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
 
    dtBaryons *= CourantSafetyNumber;
 
  }


  if (NumberOfBaryonFields > 0 && 
      (HydroMethod == HD_RK)) {

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
      exit(FAIL);
    }

    FLOAT dxinv = 1.0 / CellWidth[0][0];
    FLOAT dyinv = (GridRank > 1) ? 1.0 / CellWidth[1][0] : 0.0;
    FLOAT dzinv = (GridRank > 2) ? 1.0 / CellWidth[2][0] : 0.0;
    float dt_temp = 1.e-20, dt_ltemp, dt_x, dt_y, dt_z;
    float rho, p, vx, vy, vz, v2, eint, etot, h, cs, dpdrho, dpde,
      v_signal_x, v_signal_y, v_signal_z;
    int n = 0;
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++, n++) {
	  rho = BaryonField[DensNum][n];
	  vx  = BaryonField[Vel1Num][n];
	  vy  = BaryonField[Vel2Num][n];
	  vz  = BaryonField[Vel3Num][n];

	  if (DualEnergyFormalism) {
	    eint = BaryonField[GENum][n];
	  }
	  else {
	    etot = BaryonField[TENum][n];
	    v2 = vx*vx + vy*vy + vz*vz;
	    eint = etot - 0.5*v2;
	  }

	  EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);

	  v_signal_x = (cs + fabs(vx));
	  v_signal_y = (cs + fabs(vy));
	  v_signal_z = (cs + fabs(vz));

	  dt_x = v_signal_x * dxinv;
	  dt_y = v_signal_y * dyinv;
	  dt_z = v_signal_z * dzinv;

	  dt_ltemp = my_MAX(dt_x, dt_y, dt_z);

	  if (dt_ltemp > dt_temp) {
	    dt_temp = dt_ltemp;
	  }
        }
      }
    }

    dtBaryons = CourantSafetyNumber / dt_temp;

  }


  // MHD
  if (NumberOfBaryonFields > 0 && HydroMethod == MHD_RK) {

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
      B1Num, B2Num, B3Num, PhiNum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      return FAIL;
    }

    FLOAT dxinv = 1.0 / CellWidth[0][0];
    FLOAT dyinv = (GridRank > 1) ? 1.0 / CellWidth[1][0] : 0.0;
    FLOAT dzinv = (GridRank > 2) ? 1.0 / CellWidth[2][0] : 0.0;
    float vxm, vym, vzm, Bm, rhom;
    float dt_temp = 1.e-20, dt_ltemp, dt_x, dt_y, dt_z;
    float rho, p, vx, vy, vz, v2, eint, etot, h, cs, cs2, dpdrho, dpde,
      v_signal_x, v_signal_y, v_signal_z, cf, cf2, temp1, Bx, By, Bz, B2, ca2;
    int n = 0;
    float rho_dt, B_dt, v_dt;
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++, n++) {
	  rho = BaryonField[DensNum][n];
	  vx  = BaryonField[Vel1Num][n];
	  vy  = BaryonField[Vel2Num][n];
	  vz  = BaryonField[Vel3Num][n];
	  Bx  = BaryonField[B1Num][n];
	  By  = BaryonField[B2Num][n];
	  Bz  = BaryonField[B3Num][n];

	  if (DualEnergyFormalism) {
	    eint = BaryonField[GENum][n];
	  }
	  else {
	    etot = BaryonField[TENum][n];
	    v2 = vx*vx + vy*vy + vz*vz;
	    B2 = Bx*Bx + By*By + Bz*Bz;
	    eint = etot - 0.5*v2 - 0.5*B2/rho;
	  }

	  EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
	  cs2 = cs*cs;
	  temp1 = cs2 + B2/rho;

	  ca2 = Bx*Bx/rho;
	  cf2 = 0.5 * (temp1 + sqrt(temp1*temp1 - 4.0*cs2*ca2));
	  cf = sqrt(cf2);
	  v_signal_x = (cf + fabs(vx));

	  ca2 = By*By/rho;
	  cf2 = 0.5 * (temp1 + sqrt(temp1*temp1 - 4.0*cs2*ca2));
	  cf = sqrt(cf2);
	  v_signal_y = (cf + fabs(vy));

	  ca2 = Bz*Bz/rho;
	  cf2 = 0.5 * (temp1 + sqrt(temp1*temp1 - 4.0*cs2*ca2));
	  cf = sqrt(cf2);
	  v_signal_z = (cf + fabs(vz));

	  dt_x = v_signal_x * dxinv;
	  dt_y = v_signal_y * dyinv;
	  dt_z = v_signal_z * dzinv;

	  dt_ltemp = my_MAX(dt_x, dt_y, dt_z);

	  if (dt_ltemp > dt_temp) {
	    dt_temp = dt_ltemp;
	    rho_dt = rho;
	    B_dt = sqrt(Bx*Bx+By*By+Bz*Bz);
	    v_dt = max(fabs(vx), fabs(vy));
	    v_dt = max(v_dt, fabs(vz));
	  }

        }
      }
    }
    dtMHD = CourantSafetyNumber / dt_temp;

    //    if (dtMHD*TimeUnits/yr < 5) {
    //float ca = B_dt/sqrt(rho_dt)*VelocityUnits;
    //printf("dt=%g, rho=%g, B=%g\n, v=%g, ca=%g, dt=%g", dtMHD*TimeUnits/yr, rho_dt*DensityUnits, B_dt*MagneticUnits, 
    //    v_dt*VelocityUnits/1e5, ca/1e5, LengthUnits/(dxinv*ca)/yr*CourantSafetyNumber);
    // }

  } // HydroMethod = MHD_RK

 
  /* 2) Calculate dt from particles. */
 
  if (NumberOfParticles > 0) {
 
    /* Compute dt constraint from particle velocities. */
 
    for (dim = 0; dim < GridRank; dim++) {
      float dCell = CellWidth[dim][0]*a;
      for (i = 0; i < NumberOfParticles; i++) {
        dtTemp = dCell/max(fabs(ParticleVelocity[dim][i]), tiny_number);
	dtParticles = min(dtParticles, dtTemp);
      }
    }
 
    /* Multiply resulting dt by ParticleCourantSafetyNumber. */
 
    dtParticles *= ParticleCourantSafetyNumber;
 
  }
 

  /* 3) Find dt from expansion. */
 
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionTimestep(Time, &dtExpansion) == FAIL) {
      fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
      exit(FAIL);
    }
 
  /* 4) Calculate minimum dt due to acceleration field (if present). */
 
  if (SelfGravity) {
    for (dim = 0; dim < GridRank; dim++)
      if (AccelerationField[dim])
	for (i = 0; i < size; i++) {
	  dtTemp = sqrt(CellWidth[dim][0]/
			fabs(AccelerationField[dim][i])+tiny_number);
	  dtAcceleration = min(dtAcceleration, dtTemp);
	}
    if (dtAcceleration != huge_number)
      dtAcceleration *= 0.5;
  }
 
  /* 5) calculate minimum timestep */
 
  dt = min(dtBaryons, dtParticles);
  dt = min(dt, dtMHD);
  dt = min(dt, dtViscous);
  dt = min(dt, dtAcceleration);
  dt = min(dt, dtExpansion);

#ifdef TRANSFER

  /* 6) If star formation, set a minimum timestep */

  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    ENZO_FAIL("");
  }

  float mindtNOstars;  // Myr
  const int NumberOfStepsInLifetime = 5;
  float dtStar = huge_number;

  if (STARFEED_METHOD(POP3_STAR))
    mindtNOstars = 3;  // Myr
  if (STARFEED_METHOD(STAR_CLUSTER))
    mindtNOstars = 10;  // Myr

  if (STARFEED_METHOD(POP3_STAR) || STARFEED_METHOD(STAR_CLUSTER))
    if (G_TotalNumberOfStars > 0 && minStarLifetime < 1e6)
      dtStar = minStarLifetime/NumberOfStepsInLifetime;
    else
      dtStar = 3.1557e13*mindtNOstars/TimeUnits;

  dt = min(dt, dtStar);



  /* 7) If using radiation pressure, calculate minimum dt */

  float dtRadPressure = huge_number;
  float absVel, absAccel;

  if (RadiationPressure && RadiativeTransfer) {

    int RPresNum1, RPresNum2, RPresNum3;
    if (IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3) 
	== FAIL) {
      fprintf(stdout, "Error in IdentifyRadiationPressureFields.\n");
      ENZO_FAIL("");
    }

    for (i = 0; i < size; i++)
      for (dim = 0; dim < GridRank; dim++) {
	dtTemp = sqrt(CellWidth[dim][0] / (fabs(BaryonField[RPresNum1+dim][i])+
					   tiny_number));
	dtRadPressure = min(dtRadPressure, dtTemp);
      }
    
    if (dtRadPressure < huge_number)
      dtRadPressure *= 0.5;

    dt = min(dt, dtRadPressure);

  } /* ENDIF RadiationPressure */

  /* 8) Safety Velocity to limit timesteps */

  float dtSafetyVelocity = huge_number;
  if (TimestepSafetyVelocity > 0)
    dtSafetyVelocity = CellWidth[0][0] / 
      (TimestepSafetyVelocity*1e5 / VelocityUnits);    // parameter in km/s

  dt = min(dt, dtSafetyVelocity);
  
#endif /* TRANSFER */
 
  /* Debugging info. */
 
  if (debug1) {
    printf("ComputeTimeStep = %"FSYM" (", dt);
    if (HydroMethod != MHD_RK && NumberOfBaryonFields > 0)
      printf("Bar = %"FSYM" ", dtBaryons);
    if (HydroMethod == MHD_RK)
      printf("dtMHD = %"FSYM" ", dtMHD);
    if (HydroMethod == Zeus_Hydro)
      printf("Vis = %"FSYM" ", dtViscous);
    if (ComovingCoordinates)
      printf("Exp = %"FSYM" ", dtExpansion);
    if (dtAcceleration != huge_number)
      printf("Acc = %"FSYM" ", dtAcceleration);
    if (NumberOfParticles)
      printf("Part = %"FSYM" ", dtParticles);
    printf(")\n");
  }
 
  return dt;
}
