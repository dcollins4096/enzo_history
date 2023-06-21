/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Self-Gravity Problem Class setup routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Takes in relevant problem-defining parameters, as well as
/           Enzo data arrays.  This routine will be called repeatedly 
/           (once per time step), so it should NOT allocate any memory.  
/           Instead, it is designed to extract relevant information 
/           from other Enzo structures (cosmological expansion 
/           coefficient, gravitating mass field, etc.) that will be
/           used in defining the [non]linear residual routine at this 
/           time step.
/
/           Also use this stage as an opportunity to set up the linear 
/           self-gravity matrix and associated boundary conditions.
/
/           NOTE: this routine must be called *after* 
/           SelfGravityProblem::SetupBoundary has been called for 
/           all boundary conditions.
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#include "SelfGravityProblem.h"



// Function prototypes
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);



int SelfGravityProblem::Setup(HierarchyEntry &TopGrid) 
{

  // find the grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("Error: p%"ISYM" could not locate his grid\n",MyProcessorNumber);
    return FAIL;
  }

#ifdef USE_MPI
  //  check that MyProcessorNumber agrees with MPI process ID
  MPI_Arg MPI_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_id);
  if (MyProcessorNumber != MPI_id) {
    fprintf(stderr, "ERROR: Enzo PID %"ISYM" doesn't match MPI ID %"ISYM"\n", 
	    MyProcessorNumber, int(MPI_id));
    return FAIL;
  }
#endif

  // set up rhs vector
  //   Compute cosmology expansion term at time = t+1/2dt
  FLOAT a = 1, dadt, MidTime = ThisGrid->GridData->ReturnTime() + 
                               0.5*ThisGrid->GridData->ReturnTimeStep();
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(MidTime, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
  } // end: if (ComovingCoordinates)
  
  //    extract gravitating mass field from grid
  float *grav_mass_field = ThisGrid->GridData->GetGravitatingMassField();

  //    attach grav_mass_field to rhorhs
  rhorhs->SetData(0, grav_mass_field);


//   // TEMPORARY:  RE-SET RHORHS TO 1.0 (CURRENT ERROR IN PARALLELISM FROM ENZO)
//   rhorhs->constant(1.0);


  //    subtract rho0 from grav_mass_field, and scale by gravity constant
  rhorhs->addconst(-AverageDensity);
  rhorhs->scale(GravitationalConstant/a);

  // set up standard, 7pt Laplacian matrix in A
  //    destroy old matrix if necessary
  if (AInit == 1) {
    HYPRE_SStructMatrixDestroy(A);
    AInit = 0;
  }

  //    create the matrix, and set init flag
  HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);
  AInit = 1;

  //    set matrix storage type
  HYPRE_SStructMatrixSetObjectType(A, mattype);

  //    initialize matrix
  HYPRE_SStructMatrixInitialize(A);

  //    entries holds the stencil locations
  Eint32 entries[7] = {0, 1, 2, 3, 4, 5, 6};

  //    set matrix values over grid
  double dxi2 = 1.0/dx[0]/dx[0];
  double dyi2 = 1.0/dx[1]/dx[1];
  double dzi2 = 1.0/dx[2]/dx[2];
  double dval = -2.0*(dxi2 + dyi2 + dzi2); 
  // Nx, Ny, Nz include interior and boundaries, but no ghosts
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  double *tmpvec = new double[7*Nx*Ny*Nz];
  int idx=0, ix, iy, iz;
  for (iz=0; iz<Nz; iz++) {
    for (iy=0; iy<Ny; iy++) {
      for (ix=0; ix<Nx; ix++) {
	tmpvec[idx++] = dzi2;
	tmpvec[idx++] = dyi2;
	tmpvec[idx++] = dxi2;
	tmpvec[idx++] = dval;
	tmpvec[idx++] = dxi2;
	tmpvec[idx++] = dyi2;
	tmpvec[idx++] = dzi2;
      }
    }
  }
  //       x0L face adjustments
  if ((BdryType[0]==1) && OnBdry[0][0]) {
    ix=0;
    for (iz=0; iz<Nz; iz++) {
      for (iy=0; iy<Ny; iy++) {
	idx = 7*((iz*Ny + iy)*Nx + ix);
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = dval;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
      }
    }
  }
  //       x0R face adjustments
  if ((BdryType[0]==1) && OnBdry[0][1]) {
    ix=Nx-1;
    for (iz=0; iz<Nz; iz++) {
      for (iy=0; iy<Ny; iy++) {
	idx = 7*((iz*Ny + iy)*Nx + ix);
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = dval;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
      }
    }    
  }
  //       x1L face adjustments
  if ((BdryType[1]==1) && OnBdry[1][0]) {
    iy=0;
    for (iz=0; iz<Nz; iz++) {
      for (ix=0; ix<Nx; ix++) {
	idx = 7*((iz*Ny + iy)*Nx + ix);
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = dval;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
      }
    }    
  }
  //       x1R face adjustments
  if ((BdryType[1]==1) && OnBdry[1][1]) {
    iy=Ny-1;
    for (iz=0; iz<Nz; iz++) {
      for (ix=0; ix<Nx; ix++) {
	idx = 7*((iz*Ny + iy)*Nx + ix);
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = dval;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
      }
    }
  }
  //       x2L face adjustments
  if ((BdryType[2]==1) && OnBdry[2][0]) {
    iz=0;
    for (iy=0; iy<Ny; iy++) {
      for (ix=0; ix<Nx; ix++) {
	idx = 7*((iz*Ny + iy)*Nx + ix);
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = dval;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
      }
    }
  }
  //       x2R face adjustments
  if ((BdryType[2]==1) && OnBdry[2][1]) {
    iz=Nz-1;
    for (iy=0; iy<Ny; iy++) {
      for (ix=0; ix<Nx; ix++) {
	idx = 7*((iz*Ny + iy)*Nx + ix);
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = dval;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
	tmpvec[idx++] = 0.0;
      }
    }    
  }
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  Eint32 zed=0;
  HYPRE_SStructMatrixSetBoxValues(A, zed, ilower, iupper, zed, 
				  stSize, entries, tmpvec); 
  delete[] tmpvec;

  //    assemble matrix
  HYPRE_SStructMatrixAssemble(A);

//   // TEMPORARY: output matrix to file
//   if (debug)  printf("Writing out matrix to file selfgrav.mat\n");
//   HYPRE_SStructMatrixPrint("selfgrav.mat",A,0);

  // return success
  return SUCCESS;
}
#endif
