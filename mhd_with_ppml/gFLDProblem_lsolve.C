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
/  Gray Flux-Limited Diffusion Implicit Problem linear Newton system 
/  solution function
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Solves the linear Newton system J(u)*s = b.  For the Gray 
/           FLD Problem (one radiation group) without advection, the 
/           problem may be solved via the following Schur-complement 
/           formulation.
/
/           The Jacobian may be written in the operator form: 
/                           [ L_ee  L_en      L_eE      ]   
/               J = I3 + dt*[ L_ne  L_nn      L_nE      ] = [ M U ]
/                           [ L_Ee  L_En  (L_EE + D_EE) ]   [ L D ]
/           where for the fluid energy e, chemical species n, and 
/           radiation group r,
/               L_ee = local Jacobian of e wrt e
/               L_en = local Jacobian of e wrt n_j, j=1:Nchem
/               L_eE = local Jacobian of e wrt E
/               L_ne = local Jacobian of n_i wrt e, i=1:Nchem
/               L_nn = local Jacobian of n_i wrt n_j, i,j=1:Nchem
/               L_nE = local Jacobian of n_i wrt E, i=1:Nchem
/               L_Ee = local Jacobian of E wrt e
/               L_En = local Jacobian of E wrt n_j, j=1:Nchem
/               L_EE = local Jacobian of E wrt E
/               D_EE = diffusive Jacobian of E wrt E
/                 M = I2 + dt*[ L_ee L_en ]
/                             [ L_ne L_nn ]
/                 U = dt*[ L_eE L_nE ]^T
/                 L = dt*[ L_Ee L_En ]
/                 D = I1 + dt*(L_EE + D_EE).
/            The Schur complement formulation provides that
/                 Ji = [ I -Mi*U ][ Mi 0  ][   I   0 ]
/                      [ 0    I  ][ 0  Pi ][ -L*Mi I ]
/            where we use 'i' to denote the inverse, e.g. Ai = A^{-1}, 
/            and where the Schur complement is formed as P = D-L*Mi*U.  
/            Therefore, the solve J*s = b, where 
/            s = (s_e s_n s_E)^T = (s_m s_E)^T  (other vectors similarly 
/            condense e and n into m) may be broken down into the stages:
/                 (1) Solve for c_m:  M*c_m = b_m             (local)
/                 (2) Solve for y_m:  M*y_m = U               (local)
/                 (3) Construct:  y_E = I + dt*L_EE - L*y_m   (local)
/                     note: as this matrix is diagonal, store in y_E 
/                 (4) Update: b_E = b_E - L*c_m               (local)
/                 (5) Construct: P = dt*D_EE + I*y_E          (local)
/                 (6) Solve for s_E:  P*s_E = b_E             (nonlocal)
/                 (7) Set s_m = c_m - y_m*s_E.                (local)
/            We note that all of the steps are completely local except 
/            for the step (6), which requires one scalar-valued diffusive 
/            solve, for which we use the HYPRE library's SysPFMG solver 
/            (extensible to vector-valued diffusive systems).
/
/            The matrix components L_* have already been computed in 
/            the lsetup routine, and these are assumed to contain 
/            the time step scaling terms dt.  Moreover, the component 
/            D_EE is also already assumed to contain the time step 
/            scaling dt.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"


int gFLDProblem::lsolve(EnzoVector *s, EnzoVector *b, 
			EnzoVector *u, float delta)
{
  
//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::lsolve routine\n");

  // check that the gFLDProblem has been set up
  if (!prepared) {
    fprintf(stderr,"lsolve error: gFLDProblem not yet prepared\n");
    return FAIL;
  }
  
//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- enforcing BCs on RHS\n");

//   // have b communicate neighbor information and enforce BCs
//   b->exchange();
//   this->EnforceBoundary(b,1);

  // check that b matches local vector size (including ghosts, etc)
  int vsz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  b->size(&vsz[0], &vsz[1], &vsz[2], &vsz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (vsz[0] != LocDims[0]) {
    fprintf(stderr,"lsolve error: x0 vector dims do not match\n");
    return FAIL;
  }
  if (vsz[1] != LocDims[1]) {
    fprintf(stderr,"lsolve error: x1 vector dims do not match\n");
    return FAIL;
  }
  if (vsz[2] != LocDims[2]) {
    fprintf(stderr,"lsolve error: x2 vector dims do not match\n");
    return FAIL;
  }
  if (vsz[3] != (2+Nchem)) {
    fprintf(stderr,"lsolve error: nspecies do not match\n");
    return FAIL;
  }
  if ((vsz[0]+ghXl+ghXr) != ArrDims[0]) {
    fprintf(stderr,"lsolve error: x0 vector sizes do not match\n");
    return FAIL;
  }
  if ((vsz[1]+ghYl+ghYr) != ArrDims[1]) {
    fprintf(stderr,"lsolve error: x1 vector sizes do not match\n");
    return FAIL;
  }
  if ((vsz[2]+ghZl+ghZr) != ArrDims[2]) {
    fprintf(stderr,"lsolve error: x2 vector sizes do not match\n");
    return FAIL;
  }

//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- creating Schur comp. temp vector\n");

  // create temporary vector yvec for Schur complement correction
  EnzoVector *yvec = u->clone();


//   fprintf(stdout,"Writing out state vector to Erad.vec\n");
//   u->write("Erad.vec",0);

//   fprintf(stdout,"  ||rhs||_rms = %g\n",b->rmsnorm());
//   fprintf(stdout,"  ||rhs||_inf = %g\n",b->infnorm());
//   fprintf(stdout,"Writing out residual to Eresid.vec\n");
//   b->write("Eresid.vec",0);

//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- performing steps 1-2\n");

  //////////////////////////////////////////////////////////////
  // steps (1) and (2): local solves 
  //          c_m = Mi*b_m   (store c_m in b_m) 
  //          y_m = Mi*U
  int ix, iy, iz, idx, irow, icol, size=Nchem+1, two=2;
  float M[(1+Nchem)*(1+Nchem)], bvec[2*(1+Nchem)], xvec[2*(1+Nchem)];
  float *Lblock, *tmp;
  for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
    for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
      for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;

	// set up local M matrix
	for (irow=0; irow<(1+Nchem); irow++) {
	  for (icol=0; icol<(1+Nchem); icol++) {

	    // set matrix element with correct Jacobian component
	    Lblock = L[irow+1]->GetData(icol+1);
	    M[icol*(1+Nchem)+irow] = Lblock[idx];
	  }

	  // add one to matrix diagonal for time-independent piece
	  M[irow*(1+Nchem)+irow] += 1.0;

	  // set up rhs vector bvec (1st column is b_m, 2nd is U)
	  //    1st column contains b_m
	  tmp = b->GetData(irow+1);
	  bvec[irow] = tmp[idx];
	  //    2nd column contains U
	  tmp = (L[irow+1])->GetData(0);
	  bvec[Nchem+1+irow] = tmp[idx];
	}

	// solve dense local systems
	if (this->BlockSolve(M, xvec, bvec, &size, &two) != SUCCESS) {
	  fprintf(stderr,"lsolve: Error in BlockSolve routine\n");
	  return FAIL;
	}

	// extract solution components to appropriate locations
	for (irow=0; irow<(1+Nchem); irow++) {
	  //    put c_m back in b_m
	  tmp = b->GetData(irow+1);
	  tmp[idx]= xvec[irow];

	  //    put Mi*U into y_m vector (all of yvec except species 0)
	  tmp = yvec->GetData(irow+1);
	  tmp[idx] = xvec[Nchem+1+irow];
	}
      }
    }
  }


//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- performing steps 3-4\n");

  //////////////////////////////////////////////////////////////
  // steps (3) and (4): Construct local update for P and adjust rhs b_E
  //     (3) construct:  y_E = I + L_EE - L*y_m
  //     (4) update:     b_E = b_E - L*c_m  (note: c_m stored in b_m)
  float *y_E = yvec->GetData(0);
  float *b_E = b->GetData(0);
  float *Ldiag = (L[0])->GetData(0);
  float *b_m, *y_m;
  for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
    for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
      for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	
	// set y_E with (I + L_EE) at first
	y_E[idx] = 1.0 + Ldiag[idx];
	
	// iterate over other fluid energy, chemistry 
	// (radiation in irow 0)
	for (irow=1; irow<(2+Nchem); irow++) {

	  // L_E* is contained in L[0]
	  Lblock = (L[0])->GetData(irow);

	  // update y_E component
	  //   Mi*U is contained in y_m (the rest of yvec)
	  y_m = yvec->GetData(irow);
	  y_E[idx] -= Lblock[idx]*y_m[idx];

	  // update rhs vector b_E
	  //   store c_m in b_m
	  b_m = b->GetData(irow);
	  b_E[idx] -= Lblock[idx]*b_m[idx];
	}
      }
    }
  }

//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- writing out schur_adjustment vector\n");
//   yvec->write("schur_adjustment",0);


//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- performing step 5\n");

  //////////////////////////////////////////////////////////////
  // step (5): Construct (locally) the Schur complement matrix 
  // P = D-L*Mi*U.  In the code, this corresponds to the matrix
  // P = D_EE + I*y_E
  //       create the matrix
  HYPRE_SStructMatrix P;
  HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &P);

  //       set matrix storage type
  HYPRE_SStructMatrixSetObjectType(P, mattype);

  //       initialize matrix
  HYPRE_SStructMatrixInitialize(P);

  //       entries holds the stencil locations
  Eint32 entries[7] = {0, 1, 2, 3, 4, 5, 6};

  //       set matrix values over grid
  float *u_E = u->GetData(0);
  float *u0_E = U0->GetData(0);
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  Eflt64 *tmpvec = new Eflt64[7*Nx*Ny*Nz];
  for (ix=0; ix<7*Nx*Ny*Nz; ix++)  tmpvec[ix]=0.0;
  //       NOTE: Temp, sigmaA and sigmaS are still valid from nlresid
  if (this->MatrixEntries(tmpvec, u_E, u0_E, Temp, 
			  sigmaA, sigmaS, y_E) != SUCCESS) {
    fprintf(stderr,"lsolve: Error in MatrixEntries routine\n");
    return FAIL;
  }
  Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
  Eint32 zed=0;
//   fprintf(stderr,"lsolve: calling HYPRE_SStructMatrixSetBoxValues\n");
  HYPRE_SStructMatrixSetBoxValues(P, zed, ilower, iupper, zed, 
				  stSize, entries, tmpvec); 
  delete[] tmpvec;

  //       assemble matrix
//   fprintf(stderr,"lsolve: calling HYPRE_SStructMatrixAssemble\n");
  HYPRE_SStructMatrixAssemble(P);

  // TEMPORARY: output matrix to file
  if (debug)  printf("Writing out matrix to file P.mat\n");
  HYPRE_SStructMatrixPrint("P.mat",P,0);

//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- performing step 6\n");

  //////////////////////////////////////////////////////////////
  // step (6): Solve the (nonlocal) Schur complement system,
  // i.e. solve for s_E:  P*s_E = b_E

  //       re-scale delta to relative residual and not actual
  delta /= b->rmsnorm();
  delta = min(delta, 1.0e-2);
//   if (debug)  fprintf(stdout,"   HYPRE delta = %g\n",delta);

  //       create the HYPRE vectors
  HYPRE_SStructVector rhsvec, solvec;
//   fprintf(stderr,"lsolve: calling HYPRE_SStructVectorCreate\n");
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &rhsvec);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &solvec);

  //       set vector storage type
//   fprintf(stderr,"lsolve: calling HYPRE_SStructVectorSetObjectType\n");
  HYPRE_SStructVectorSetObjectType(rhsvec, mattype);
  HYPRE_SStructVectorSetObjectType(solvec, mattype);

  //       initialize vectors
//   fprintf(stderr,"lsolve: calling HYPRE_SStructVectorInitialize\n");
  HYPRE_SStructVectorInitialize(rhsvec);
  HYPRE_SStructVectorInitialize(solvec);

  //       insert rhs, sol vectors into HYPRE vectors x and b
  ilower[0] = SolvIndices[0][0];
  iupper[0] = SolvIndices[0][1];
  int xBuff, yBuff, zBuff;
  xBuff = ghXl;
  yBuff = ghYl-SolvIndices[1][0];
  zBuff = ghZl-SolvIndices[2][0];
  float *s_E = s->GetData(0);
  int Zbl, Ybl;
//   fprintf(stderr,"lsolve: calling HYPRE_SStructVectorSetBoxValues\n");
  Eflt64 *buff = new Eflt64[SolvIndices[0][1]-SolvIndices[0][0]+1];
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
    ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];
      ilower[1] = iy;  iupper[1] = iy;
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	buff[ix] = b_E[Zbl+Ybl+xBuff+ix];
      HYPRE_SStructVectorSetBoxValues(rhsvec, zed, ilower, iupper, zed, buff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	buff[ix] = s_E[Zbl+Ybl+xBuff+ix];
      HYPRE_SStructVectorSetBoxValues(solvec, zed, ilower, iupper, zed, buff);
    }
  }

  //       assemble vectors
//   fprintf(stderr,"lsolve: calling HYPRE_SStructVectorAssemble\n");
  HYPRE_SStructVectorAssemble(solvec);
  HYPRE_SStructVectorAssemble(rhsvec);

  // TEMPORARY: output rhs to file
  if (debug)  printf("Writing out rhs to file rhs.vec\n");
  HYPRE_SStructVectorPrint("rhs.vec",rhsvec,zed);

  //       set up the solver [SysPFMG]
  //          create the solver
  HYPRE_SStructSolver solver;
  fprintf(stderr,"lsolve: calling HYPRE_SStructSysPFMGCreate\n");
  HYPRE_SStructSysPFMGCreate(MPI_COMM_WORLD, &solver);

  //          set solver options
//   fprintf(stderr,"lsolve: calling HYPRE_SStructSysPFMGSet*\n");
  HYPRE_SStructSysPFMGSetMaxIter(solver, sol_maxit);
  HYPRE_SStructSysPFMGSetRelChange(solver, sol_relch);
  HYPRE_SStructSysPFMGSetNumPreRelax(solver, sol_npre);
  HYPRE_SStructSysPFMGSetNumPostRelax(solver, sol_npost);
  HYPRE_SStructSysPFMGSetPrintLevel(solver, sol_printl);
  HYPRE_SStructSysPFMGSetLogging(solver, sol_log);
  if (delta != 0.0)   HYPRE_SStructSysPFMGSetTol(solver, Eflt64(delta));
  if (sol_zeroguess)  HYPRE_SStructSysPFMGSetZeroGuess(solver);
//   fprintf(stderr,"lsolve: calling HYPRE_SStructSysPFMGSetup\n");
  HYPRE_SStructSysPFMGSetup(solver, P, rhsvec, solvec);

  //       solve the linear system
//   fprintf(stderr,"lsolve: calling HYPRE_SStructSysPFMGSolve\n");
  HYPRE_SStructSysPFMGSolve(solver, P, rhsvec, solvec);

  //       extract solver statistics
  double finalresid;
  Eint32 its;
//   fprintf(stderr,"lsolve: calling HYPRE_SStructSysPFMGGet*\n");
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_SStructSysPFMGGetNumIterations(solver, &its);
  totIters += its;
  if (debug) 
    fprintf(stdout,"   HYPRE solve: iters = %i, relresid = %g\n",its,finalresid);

  //       gather the solution vector and extract values
//   fprintf(stderr,"lsolve: calling HYPRE_SStructVectorGather\n");
  HYPRE_SStructVectorGather(solvec);
//   fprintf(stderr,"lsolve: calling HYPRE_SStructVectorGetBoxValues\n");
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
    ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];
      ilower[1] = iy;  iupper [1] = iy;
      HYPRE_SStructVectorGetBoxValues(solvec, zed, ilower, iupper, zed, buff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	s_E[Zbl+Ybl+xBuff+ix] = buff[ix];
    }
  }
  delete[] buff;

  // TEMPORARY: output rhs to file
  if (debug)  printf("Writing out solution to file sol.vec\n");
  HYPRE_SStructVectorPrint("sol.vec",solvec,zed);

  //       destroy HYPRE vector and solver structures
//   fprintf(stderr,"lsolve: calling HYPRE_SStruct*Destroy\n");
  HYPRE_SStructSysPFMGDestroy(solver);
  HYPRE_SStructVectorDestroy(rhsvec);
  HYPRE_SStructVectorDestroy(solvec);
  HYPRE_SStructMatrixDestroy(P);


//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- performing step 7\n");

  //////////////////////////////////////////////////////////////
  // step (7) Set s_m = c_m - y_m*s_E  (Note: c_m stored in b_m)
  float *s_m;
  for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
    for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
      for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	for (irow=1; irow<(2+Nchem); irow++) {
	  b_m = b->GetData(irow);
	  y_m = yvec->GetData(irow);
	  s_m = s->GetData(irow);
	  s_m[idx] = b_m[idx] - y_m[idx]*s_E[idx];
	}
      }
    }
  }

//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- deleting Schur temp. vec.\n");

  // delete temporary vector yvec
  delete yvec;

//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- solution neighbor exchange\n");

  // have s communicate neighbor information
  s->exchange();

//   if (debug)
//     fprintf(stdout,"gFLDProblem::lsolve -- finished!\n");

  // return success
  return SUCCESS;
}
#endif
