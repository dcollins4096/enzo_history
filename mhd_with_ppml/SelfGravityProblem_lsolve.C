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
/  Self-Gravity Problem linear Newton system solution function
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Solves the linear Newton system J(u)*s = b.  If based on
/           an iterative solver, reduces the linear residual below 
/           tolerance delta.
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#include "SelfGravityProblem.h"



int SelfGravityProblem::lsolve(EnzoVector *s, EnzoVector *b, 
			       EnzoVector *u, float delta)
{
  
  // check that the SelfGravityProblem has been set up
  if (!prepared) {
    fprintf(stderr,"lsolve error: SelfGravityProblem not yet prepared\n");
    return FAIL;
  }

  // check that the matrix has been initialized
  if (AInit == 0) {
    fprintf(stderr,"lsolve error: self gravity matrix un-initialized\n");
    return FAIL;
  }

  // have b communicate neighbor information and enforce BCs
  b->exchange();
  this->EnforceBoundary(b,1);

  // re-scale delta so that it solves to actual residual and not relative
  delta /= b->rmsnorm();

  // check that b matches local vector size (including ghosts, etc)
  int vsize[4], vghosts[3][2];
  b->size(&vsize[0], &vsize[1], &vsize[2], &vsize[3], 
	  &vghosts[0][0], &vghosts[0][1], &vghosts[1][0], 
	  &vghosts[1][1], &vghosts[2][0], &vghosts[2][1]);
  if (vsize[0] != (EdgeIndices[0][1]-EdgeIndices[0][0]+1)) {
    fprintf(stderr,"lsolve error: x0 vector dims do not match\n");
    return FAIL;
  }
  if (vsize[1] != (EdgeIndices[1][1]-EdgeIndices[1][0]+1)) {
    fprintf(stderr,"lsolve error: x1 vector dims do not match\n");
    return FAIL;
  }
  if (vsize[2] != (EdgeIndices[2][1]-EdgeIndices[2][0]+1)) {
    fprintf(stderr,"lsolve error: x2 vector dims do not match\n");
    return FAIL;
  }
  if ((vsize[0]+vghosts[0][0]+vghosts[0][1]) != ArrDims[0]) {
    fprintf(stderr,"lsolve error: x0 vector sizes do not match\n");
    return FAIL;
  }
  if ((vsize[1]+vghosts[1][0]+vghosts[1][1]) != ArrDims[1]) {
    fprintf(stderr,"lsolve error: x1 vector sizes do not match\n");
    return FAIL;
  }
  if ((vsize[2]+vghosts[2][0]+vghosts[2][1]) != ArrDims[2]) {
    fprintf(stderr,"lsolve error: x2 vector sizes do not match\n");
    return FAIL;
  }

  // create the HYPRE vectors
  HYPRE_SStructVector bvec, xvec;
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &bvec);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &xvec);

  // set vector storage type
  HYPRE_SStructVectorSetObjectType(bvec, mattype);
  HYPRE_SStructVectorSetObjectType(xvec, mattype);

  // initialize vectors
  HYPRE_SStructVectorInitialize(bvec);
  HYPRE_SStructVectorInitialize(xvec);

  // insert rhs vector into HYPRE vectors x and b
  Eint32 ilower[3], iupper[3];
  ilower[0] = SolvIndices[0][0];
  iupper[0] = SolvIndices[0][1];
  int xBuff, yBuff, zBuff;
  xBuff = vghosts[0][0]+SolvIndices[0][0]-EdgeIndices[0][0];
  yBuff = vghosts[1][0]-EdgeIndices[1][0];
  zBuff = vghosts[2][0]-EdgeIndices[2][0];
  float *bdata = b->GetData(0);
  float *sdata = s->GetData(0);
  int ix, iy, iz, Zbl, Ybl;
  Eint32 zed=0;
  Eflt64 *buff = new Eflt64[SolvIndices[0][1]-SolvIndices[0][0]+1];
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
    ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];
      ilower[1] = iy;  iupper[1] = iy;
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	buff[ix] = bdata[Zbl+Ybl+xBuff+ix];
      HYPRE_SStructVectorSetBoxValues(bvec, zed, ilower, iupper, zed, buff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	buff[ix] = sdata[Zbl+Ybl+xBuff+ix];
      HYPRE_SStructVectorSetBoxValues(xvec, zed, ilower, iupper, zed, buff);
    }
  }

  // assemble vectors
  HYPRE_SStructVectorAssemble(xvec);
  HYPRE_SStructVectorAssemble(bvec);

  // set up the solver [SysPFMG]
  //    create the solver
  HYPRE_SStructSolver solver;
  HYPRE_SStructSysPFMGCreate(MPI_COMM_WORLD, &solver);

  //    set solver options
  HYPRE_SStructSysPFMGSetMaxIter(solver, sol_maxit);
  HYPRE_SStructSysPFMGSetRelChange(solver, sol_relch);
  HYPRE_SStructSysPFMGSetNumPreRelax(solver, sol_npre);
  HYPRE_SStructSysPFMGSetNumPostRelax(solver, sol_npost);
  HYPRE_SStructSysPFMGSetPrintLevel(solver, sol_printl);
  HYPRE_SStructSysPFMGSetLogging(solver, sol_log);
  if (delta != 0.0)  HYPRE_SStructSysPFMGSetTol(solver, delta);
  if (sol_zeroguess) HYPRE_SStructSysPFMGSetZeroGuess(solver);
  HYPRE_SStructSysPFMGSetup(solver, A, bvec, xvec);

  // solve the linear system
  HYPRE_SStructSysPFMGSolve(solver, A, bvec, xvec);

  // extract solver statistics
  Eflt64 finalresid;
  Eint32 its;
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_SStructSysPFMGGetNumIterations(solver, &its);
  totIters += its;
  if (debug) 
    fprintf(stdout,"   HYPRE solve: iters = %i, relresid = %g\n",its,finalresid);

  // gather the solution vector and extract values
  HYPRE_SStructVectorGather(xvec);
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
    ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];
      ilower[1] = iy;  iupper [1] = iy;
      HYPRE_SStructVectorGetBoxValues(xvec, zed, ilower, iupper, zed, buff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	sdata[Zbl+Ybl+xBuff+ix] = buff[ix];
    }
  }
  delete[] buff;

  // destroy HYPRE vector and solver structures
  HYPRE_SStructSysPFMGDestroy(solver);
  HYPRE_SStructVectorDestroy(bvec);
  HYPRE_SStructVectorDestroy(xvec);

  // have s communicate neighbor information
  s->exchange();

  // return success
  return SUCCESS;
}
#endif
