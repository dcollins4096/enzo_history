/***********************************************************************
/
/  EVOLVE ROUTINE FOR COUPLED RADIATION-HYDRODYNAMICS-CHEMICAL 
/  KINETICS SYSTEM
/
/  written by: Daniel Reynolds
/  date:       November, 2006
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
#ifdef RAD_HYDRO
 
#include <stdio.h>
#include <math.h>

#include "gFLDProblem_preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "InexactNewton.h"
#include "EnzoVector.h"
#include "gFLDProblem.h"



int EvolveRadHydro(HierarchyEntry *ThisGrid, gFLDProblem *FLDsolver)
{

  if (debug)
    fprintf(stdout,"Entering EvolveRadHydro routine\n");

  // Only continue if we own this grid
  if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber())
    return SUCCESS;


  // Get grid size
  int size=1, i;
  for (i=0; i<MAX_DIMENSION; i++) 
    size *= ThisGrid->GridData->GetGridDimension(i);

  // Get number of chemical species from FLDsolver
  int RadHydroNChem = FLDsolver->GetNumChemicalSpecies();

  // Set pointers to each variable (allocate fluid energy correction 
  // since that is not internal to Enzo)
  float *FluidEnergy = ThisGrid->GridData->AccessTotalEnergy();
  float *FluidEnergyCorrection = new float[size];
  for (i=0; i<size; i++)  FluidEnergyCorrection[i] = 0.0;
  float *RadiationEnergy = ThisGrid->GridData->AccessRadiationFrequency0();
  float **Ni = new float *[20];   // make large enough for future
  if (RadHydroNChem > 0) 
    Ni[0] = ThisGrid->GridData->AccessHIDensity();
  if (RadHydroNChem > 1) 
    Ni[1] = ThisGrid->GridData->AccessHeIDensity();
  if (RadHydroNChem > 2) 
    Ni[2] = ThisGrid->GridData->AccessHeIIDensity();

  // Get time-related information
  float told = ThisGrid->GridData->ReturnTime();
  float dt = ThisGrid->GridData->ReturnTimeStep();
 
  // Set up the gFLDProblem for this time step
  if (FLDsolver->Setup(*ThisGrid, dt, told, RadiationEnergy, 
		       FluidEnergyCorrection, Ni) == FAIL) {
    fprintf(stderr,"Error in gFLDProblem::Setup\n");
    return FAIL;
  }

  //   obtain example vector from FLDsolver
  EnzoVector *oldvec = FLDsolver->ExampleVector();

  //   initialize solution vector and copy previous state
  EnzoVector *newvec = oldvec->clone();
  newvec->copy(oldvec);

  //   initialize solver object
  InexactNewtonSolver *INSolve = new InexactNewtonSolver(newvec);

  //   set nonlinear solver parameters
  INSolve->SetMaxIters(FLDsolver->GetMaxNewtonIters());
  INSolve->SetInexactNewton(0, FLDsolver->GetInexactNewtonConstant(), 1.0);
  INSolve->SetNewtonTolerance(FLDsolver->GetNewtonTolerance());
  INSolve->SetNewtonNorm(FLDsolver->GetNewtonNormChoice());
  INSolve->SetMinLinesearch(FLDsolver->GetNewtonMinLinesearch());


  // Call nonlinear solver to compute updated time step
  if (INSolve->Solve(FLDsolver,newvec) == FAIL) {
    fprintf(stderr,"ERROR: INSolve failure\n");
    return FAIL;
  }

//   float onenorm = newvec->l1norm();
//   if (debug)
//     fprintf(stdout,"\nEvolveRadHydro:  ||newvec||_1 = %g\n",onenorm);

  // Update solution data with new values
  //   Radiation Energy, Chemical Species and Fluid Correction are in newvec
  oldvec->copy(newvec);
  //   add fluid correction to fluid energy field
  for (i=0; i<size; i++)
    FluidEnergy[i] += FluidEnergyCorrection[i];

//   onenorm = oldvec->l1norm();
//   if (debug)
//     fprintf(stdout,"\nEvolveRadHydro:  ||oldvec||_1 = %g\n",onenorm);
 
//   if (debug) {
//     float inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += FluidEnergyCorrection[i]*FluidEnergyCorrection[i];
//     fprintf(stdout,"\nEvolveRadHydro:  ||ec||_2 = %g\n",sqrt(inorm));

//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += FluidEnergy[i]*FluidEnergy[i];
//     fprintf(stdout,"EvolveRadHydro:  ||e||_2 = %g\n",sqrt(inorm));

//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += RadiationEnergy[i]*RadiationEnergy[i];
//     fprintf(stdout,"EvolveRadHydro:  ||Eg||_2 = %g\n",sqrt(inorm));
    
//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += Ni[0][i]*Ni[0][i];
//     fprintf(stdout,"EvolveRadHydro:  ||nHI||_2 = %g\n\n",sqrt(inorm));
//   }


  // Clean up
  delete INSolve;
  delete newvec;
  delete[] Ni;
  delete[] FluidEnergyCorrection;


  // Return
  return SUCCESS;
 
}

#endif   // RAD_HYDRO
