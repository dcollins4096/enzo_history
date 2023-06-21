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
/  Gray Flux-Limited Diffusion Implicit Problem Class setup routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Takes in relevant problem-defining parameters, as well as
/           Enzo data arrays.  These arrays may be converted from the 
/           finite-volume Enzo format to another format amenable to 
/           the implicit solve if required.  The computations to 
/           determine these transformations will have previously 
/           occured in the Problem constructor, but the actual data 
/           transformations/interpolations must occur here.  This 
/           routine will be called repeatedly (once per time step), so 
/           it should NOT allocate any memory.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"


// Function prototypes
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);



int gFLDProblem::Setup(HierarchyEntry &TopGrid, float Dt, float tval, 
		       float *Er, float *ec, float **ni) 
{
  
//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::Setup routine\n");

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

  // store input args into class
  dt = Dt;
  told = tval;
  tnew = tval+dt;

  // store cosmology expansion parameters
  aold = 1.0; adotold = 0.0;
  anew = 1.0; adotnew = 0.0;
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(told, &aold, &adotold) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
    if (CosmologyComputeExpansionFactor(tnew, &anew, &adotnew) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
  }

  // attach arrays to U0 vector
  U0->SetData(0, Er);
  U0->SetData(1, ec);
  int ns;
  for (ns=1; ns<=Nchem; ns++)  U0->SetData(ns+1, ni[ns-1]);
  U0->exchange();

  // access velocities, density and total fluid energy fields
  vx  = ThisGrid->GridData->AccessVelocity1();
  if (vx == NULL) {
    fprintf(stderr,"gFLDProblem Setup: could not obtain velocity1\n");
    return FAIL;
  }
  vy  = ThisGrid->GridData->AccessVelocity2();
  if (vy == NULL) {
    fprintf(stderr,"gFLDProblem Setup: could not obtain velocity2\n");
    return FAIL;
  }
  vz  = ThisGrid->GridData->AccessVelocity3();
  if (vz == NULL) {
    fprintf(stderr,"gFLDProblem Setup: could not obtain velocity3\n");
    return FAIL;
  }
  rho = ThisGrid->GridData->AccessDensity();
  if (rho == NULL) {
    fprintf(stderr,"gFLDProblem Setup: could not obtain density\n");
    return FAIL;
  }
  eh  = ThisGrid->GridData->AccessTotalEnergy();
  if (eh == NULL) {
    fprintf(stderr,"gFLDProblem Setup: could not obtain fluid energy\n");
    return FAIL;
  }
  if (Nchem > 0) {
    ne = ThisGrid->GridData->AccessElectronDensity();
    if (ne == NULL) {
      fprintf(stderr,"gFLDProblem Setup: could not obtain electron density\n");
      return FAIL;
    }
  }

//   if (debug) {
//     fprintf(stdout,"\n  gFLDProblem Setup: current internal quantities:\n");
//     fprintf(stdout,"        ec typical values = %g, %g, %g\n",ec[1],ec[10],ec[100]);
//     fprintf(stdout,"        Eg typical values = %g, %g, %g\n",Er[1],Er[10],Er[100]);
//     fprintf(stdout,"       nHI typical values = %g, %g, %g\n",ni[0][1],ni[0][10],ni[0][100]);
//     fprintf(stdout,"        vx typical values = %g, %g, %g\n",vx[1],vx[10],vx[100]);
//     fprintf(stdout,"        vy typical values = %g, %g, %g\n",vy[1],vy[10],vy[100]);
//     fprintf(stdout,"        vz typical values = %g, %g, %g\n",vz[1],vz[10],vz[100]);
//     fprintf(stdout,"       rho typical values = %g, %g, %g\n",rho[1],rho[10],rho[100]);
//     fprintf(stdout,"        eh typical values = %g, %g, %g\n",eh[1],eh[10],eh[100]);
//     fprintf(stdout,"        ne typical values = %g, %g, %g\n\n",ne[1],ne[10],ne[100]);
//   }

//   if (debug) {
//     float inorm = 0.0;
//     float size=ArrDims[0]*ArrDims[1]*ArrDims[2];
//     for (i=0; i<size; i++) 
//       inorm += ec[i]*ec[i];
//     fprintf(stdout,"\ngFLDProblem Setup:  ||ec||_2 = %g\n",sqrt(inorm));

//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += Er[i]*Er[i];
//     fprintf(stdout,"gFLDProblem Setup:  ||Eg||_2 = %g\n",sqrt(inorm));
    
//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += ni[0][i]*ni[0][i];
//     fprintf(stdout,"gFLDProblem Setup:  ||nHI||_2 = %g\n",sqrt(inorm));
    
//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += vx[i]*vx[i];
//     fprintf(stdout,"gFLDProblem Setup:  ||vx||_2 = %g\n",sqrt(inorm));
    
//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += vy[i]*vy[i];
//     fprintf(stdout,"gFLDProblem Setup:  ||vy||_2 = %g\n",sqrt(inorm));
    
//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += vz[i]*vz[i];
//     fprintf(stdout,"gFLDProblem Setup:  ||vz||_2 = %g\n",sqrt(inorm));

//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += rho[i]*rho[i];
//     fprintf(stdout,"gFLDProblem Setup:  ||rho||_2 = %g\n",sqrt(inorm));

//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += eh[i]*eh[i];
//     fprintf(stdout,"gFLDProblem Setup:  ||eh||_2 = %g\n",sqrt(inorm));

//     inorm = 0.0;
//     for (i=0; i<size; i++) 
//       inorm += ne[i]*ne[i];
//     fprintf(stdout,"gFLDProblem Setup:  ||ne||_2 = %g\n\n",sqrt(inorm));

//   }


  // set up rhs0 for the current state, at current time t0
  if (this->ComputeRHS(rhs0, told, U0) == FAIL) {
    fprintf(stderr,"gFLDProblem Setup: Error in ComputeRHS routine\n");
    return FAIL;
  }

//   if (debug)
//     fprintf(stdout,"gFLDProblem Setup complete!\n");

  // set prepared flag to true
  prepared = true;

  // return success
  return SUCCESS;
}
#endif
