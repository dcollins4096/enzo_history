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


extern "C" void
FORTRAN_NAME(mhd_dt)(float *bxc, float *byc, float *bzc, 
		     float *vx, float *vy, float *vz,
		     float *d, float *p, float *gamma, float *dt, 
		     FLOAT *dx, FLOAT *dy, FLOAT *dz, 
		     int *idim, int *jdim, int *kdim, int * rank,
		     int *i1, int *i2,
		     int *j1, int *j2,
		     int *k1, int *k2, float* eng);

void dump(float *A, int nx, int ny, int nz, int nb, char * filename);

float grid::ComputeTimeStep(int level)
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


  /* Isothermal timestepping is easier? */
  if( Gamma < 0 ){
    //dt = 0.003;
    //fprintf(stderr,"kludge; hacking dt = %f  Need isothermal sound business.\n", dt);
    //return 0.003;
  }

  /*
  dt = 1.0/510.0;
  fprintf(stderr,"kludge: hacking dt = %f\n", dt);
  return dt;
  */
  /* Compute the field size. */
  
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
  
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
	  if(CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL){
	  fprintf(stderr,"Error in CosmologyComputeExpansionFactor,\n");
	  exit(EXIT_FAILURE);
	  }
  float afloat = float(a);
  
  /* 1) Compute Courant condition for baryons. */
  
  if (NumberOfBaryonFields > 0) {
    
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
      fprintf(stderr, "Error in grid->ComputePressure, called in Grid_ComputeTimeStep.\n");
      exit(EXIT_FAILURE);
    }
    
    /* Call fortran routine to do calculation. */

    //MHD timesteps are stricly smaller than Hydro.  It's a waste of time to do this step
    //if MHD is in use.
    if( MHD_Used != TRUE )
      FORTRAN_NAME(calc_dt)(&GridRank, GridDimension, GridDimension+1,
			    GridDimension+2,
			    GridStartIndex, GridEndIndex,
			    GridStartIndex+1, GridEndIndex+1,
			    GridStartIndex+2, GridEndIndex+2,
			    &HydroMethod,
                            &ZEUSQuadraticArtificialViscosity,
			    CellWidth[0], CellWidth[1], CellWidth[2],
			    GridVelocity, GridVelocity+1, GridVelocity+2,
			    &Gamma, &PressureFree, &afloat,
			    BaryonField[DensNum], pressure_field,
			    BaryonField[Vel1Num], BaryonField[Vel2Num],
			    BaryonField[Vel3Num], &dtBaryons, &dtViscous);
    

#ifdef ATHENA
    /*
    //Viscous timestep criterion, as implemented by the ATHENA solver.
    if( MHD_DiffusionMethod[0] != 0 || MHD_DiffusionMethod[1] != 0 ){
      if( MHD_dtVisc( dtViscous ) == FAIL ){
	fprintf(stderr,"ComputeTimeStep: MHD_dtVisc failed\n");
      }
    }
    */
#endif //ATHENA

    /* Clean up */

    //for(dcci=0;dcci<size;dcci++)
    //  cerr << "PDEBUG2 "  << pressure_field[dcci] << endl;
    
    
    /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
    
    dtBaryons *= CourantSafetyNumber;
    
    if(MHD_Used){
      
      /* 1.5) Calculate minimum dt due to MHD: Maximum Fast MagnetoSonic Shock Speed */
      
      //Cosmos nees this, for some reason.
      if(GridRank < 3 ){
	if( CellWidth[2] == NULL ) CellWidth[2] = new float;
	CellWidth[2][0] = 1.0;
	if( GridRank < 2 ){
	  if( CellWidth[1] == NULL ) CellWidth[1] = new float;
	  CellWidth[1][0] = 1.0;
	}}

      FORTRAN_NAME(mhd_dt)(CenteredB[0], CenteredB[1], CenteredB[2],
			   BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
			   BaryonField[DensNum], pressure_field, &Gamma, &dtMHD, 
			   CellWidth[0], CellWidth[1], CellWidth[2],
			   GridDimension, GridDimension + 1, GridDimension +2,&GridRank,
			   GridStartIndex, GridEndIndex,
			   GridStartIndex+1, GridEndIndex+1,
			   GridStartIndex+2, GridEndIndex+2, BaryonField[TENum]);
      
      dtMHD *= CourantSafetyNumber;
#ifdef HAOXU
	  dtMHD *= afloat;  
#endif /* HAOXU */
    }//if MHD_Used
    
    delete pressure_field;
       
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
      if (AccelerationField[dim] != NULL){

	for (i = 0; i < size; i++) {

	  dtTemp = sqrt(CellWidth[dim][0]/
			(fabs(AccelerationField[dim][i])+tiny_number));
			
	  dtAcceleration = min(dtAcceleration, dtTemp);
	}
      }

    if (dtAcceleration != huge_number)
      dtAcceleration *= 0.5;
  }

  /* 6 ) calculate minimum timestep */
  
  dt = min(dtMHD, dtBaryons);
  dt = min(dt, dtParticles);
  dt = min(dt, dtViscous);
  dt = min(dt, dtAcceleration);
  dt = min(dt, dtExpansion);


  /* Debugging info. */

//  if (debug || NumberOfProcessors > 1) {

  if (debug) {
    printf("ComputeTimeStep = %e (", dt);
    if (NumberOfBaryonFields > 0)
      printf("Baryon = %e ", dtBaryons);
    if (HydroMethod == Zeus_Hydro 
#ifdef ATHENA
|| MHD_DiffusionMethod[0] != 0 || MHD_DiffusionMethod[1] != 0
#endif //ATHENA
)      printf("Vis = %e ", dtViscous);
    if (ComovingCoordinates)
      printf("Exp = %e ", dtExpansion);
    if (dtAcceleration != huge_number)
      printf("Acc = %e ", dtAcceleration);
    if (NumberOfParticles)
      printf("Part = %e ", dtParticles);
    if (MHD_Used) 
      printf("MHD = %e", dtMHD);
    //printf(" Ratio: dtMHD/dtBaryons = %f",dtMHD/dtBaryons);
    printf(")\n");
  }


  //Fixed Timestep.
  //Set to -1: ignored.
  //Set set to something > 0: Fixed Timestep is used (regardless of other dt's)
  //    Checks for FixedDT > Computed DT, issues warning only.
  //Set to -2: Looks for a file called TimeStep<i>.  Reads number from that file.
  //    Checks for ReadDT > Computed DT, issues warning only.
  if( MHD_FixedTimestep > -1e-6 ){
    if( MHD_FixedTimestep > dt ){
      fprintf(stderr, "WARNING: Fixed Timestep (%g) > Computed (%g).  Using (%g)\n",
	      MHD_FixedTimestep, dt, MHD_FixedTimestep);
    }
    dt = MHD_FixedTimestep;
  }//Fixed dt > 0

  if( fabs(MHD_FixedTimestep + 2) < 1e-6 ){  

    float TempInputStep;
    FILE * StepToRead;
    char FileName[20];
    sprintf(FileName,"step_sizes/TimeStep%d",dccCounter9);

    StepToRead = fopen(FileName, "r");

    if( StepToRead != NULL ){
      fscanf(StepToRead,"%"FSYM,&TempInputStep); 
      fclose(StepToRead);
      fprintf(stderr,"Read Timestep (%"FSYM") from file %s\n", TempInputStep,FileName);
      if( TempInputStep > dt ){
	fprintf(stderr, "WARNING: Read Timestep (%g) > Computed (%g).  Using (%g)\n",
		TempInputStep, dt, TempInputStep);
      }
      dt = TempInputStep;
    }else{
      fprintf(stderr,"Timestep File Open Error: %s\n", FileName);
    }
  }//FixedTimeStep == -2



  return dt;
}

