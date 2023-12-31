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
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionTimestep(FLOAT time, float *dtExpansion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
extern "C" void FORTRAN_NAME(calc_dt)(
                  int *rank, int *idim, int *jdim, int *kdim,
                  int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
			     hydro_method *ihydro, float *C2,
                  FLOAT *dx, FLOAT *dy, FLOAT *dz, float *vgx, float *vgy,
                             float *vgz, float *gamma, int *ipfree, float *aye,
                  float *d, float *p, float *u, float *v, float *w,
			     float *dt, float *dtviscous);
 
#ifdef PPML
extern "C" void
FORTRAN_NAME(mhd_dt)(float *bxc, float *byc, float *bzc, 
		     float *vx, float *vy, float *vz,
		     float *d, float *p, float *gamma, float *dt, 
		     FLOAT *dx, FLOAT *dy, FLOAT *dz, 
		     int *idim, int *jdim, int *kdim, int * rank,
		     int *i1, int *i2,
		     int *j1, int *j2,
		     int *k1, int *k2);//, float* eng);


#endif //PPML
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
 
  if (NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
#ifdef PPML
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, BXNum, BYNum, BZNum;
    if (this->IdentifyMHDQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum,BXNum, BYNum, BZNum) == FAIL) {
      fprintf(stderr, "ComputeTimeStep: IdentifyMHDQuantities error.\n");
      exit(FAIL);
    }
#else //PPML
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
      exit(FAIL);
    }
#endif //PPM
 
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
#ifdef PPML
    if( MHD_Used == 0 ){
#endif //PPML
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

      
#ifdef PPML
    }else{


      FORTRAN_NAME(mhd_dt)(BaryonField[BXNum], BaryonField[BYNum], BaryonField[BZNum],
			   BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
			   BaryonField[DensNum], pressure_field, &Gamma, &dtBaryons, 
			   CellWidth[0], CellWidth[1], CellWidth[2],
			   GridDimension, GridDimension + 1, GridDimension +2,&GridRank,
			   GridStartIndex, GridEndIndex,
			   GridStartIndex+1, GridEndIndex+1,
			   GridStartIndex+2, GridEndIndex+2);//, BaryonField[TENum]);
      
    }
#endif //PPML

      
    /* Clean up */
 
    delete pressure_field;
 
    /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
 
    dtBaryons *= CourantSafetyNumber;
 
  }
 
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
  dt = min(dt, dtViscous);
  dt = min(dt, dtAcceleration);
  dt = min(dt, dtExpansion);
 
  /* Debugging info. */
 
//  if (debug || NumberOfProcessors > 1) {
  if (debug) {
    printf("ComputeTimeStep = %"FSYM" (", dt);
    if (NumberOfBaryonFields > 0)
      printf("Bar = %"FSYM" ", dtBaryons);
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
#ifdef PPML

  //Fixed Timestep.
  //Set to -1: ignored.
  //Set set to something > 0: Fixed Timestep is used (regardless of other dt's)
  //    Checks for FixedDT > Computed DT, issues warning only.
  //Set to -2: Looks for a file called TimeStep<i>.  Reads number from that file.
  //    Checks for ReadDT > Computed DT, issues warning only.
  if( FixedTimestep > -1e-6 ){
    if( FixedTimestep > dt ){
      fprintf(stderr, "WARNING: Fixed Timestep (%"GSYM") > Computed (%"GSYM").  Using (%"GSYM")\n",
	      FixedTimestep, dt, FixedTimestep);
    }
    dt = FixedTimestep;
  }//Fixed dt > 0

  if( fabs(FixedTimestep + 2) < 1e-6 ){  

    float TempInputStep;
    FILE * StepToRead;
    char FileName[20];
    sprintf(FileName,"step_sizes/TimeStep%d",dccCounter02);

    StepToRead = fopen(FileName, "r");

    if( StepToRead != NULL ){
      fscanf(StepToRead,"%"FSYM,&TempInputStep); 
      fclose(StepToRead);
      fprintf(stderr,"Read Timestep (%"FSYM") from file %s\n", TempInputStep,FileName);
      if( TempInputStep > dt ){
	fprintf(stderr, "WARNING: Read Timestep (%"GSYM") > Computed (%"GSYM").  Using (%"GSYM")\n",
		TempInputStep, dt, TempInputStep);
      }
      dt = TempInputStep;
    }else{
      fprintf(stderr,"Timestep File Open Error: %s\n", FileName);
    }
  }//FixedTimeStep == -2

#endif //PPML 
  return dt;
}
