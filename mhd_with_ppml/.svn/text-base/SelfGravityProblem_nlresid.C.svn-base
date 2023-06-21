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
/  Self-Gravity Problem nonlinear residual function
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: Residual function that defines the implicit self-gravity 
/           problem.  As the problem is linear, just compute the 
/           linear residual using MGMPI with the already-defined matrix.
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#include "SelfGravityProblem.h"



int SelfGravityProblem::nlresid(EnzoVector *fu, EnzoVector *u)
{

  // check that the SelfGravityProblem has been prepared
  if (!prepared) {
    fprintf(stderr,"nlresid error: SelfGravityProblem not yet prepared\n");
    return FAIL;
  }
  
//   // write initial state to file
//   if (debug)
//     fprintf(stdout,"Writing initial state to u0_vec\n");
//   u->writeall("u0_vec",0);

  // have u communicate neighbor information
  u->exchange();

//   // write communicated state to file
//   if (debug)
//     fprintf(stdout,"Writing communicated state to u0comm_vec\n");
//   u->writeall("u0comm_vec",0);

  // enforce boundary conditions on state u
  this->EnforceBoundary(u,0);

//   // write BC-enforced state to file
//   if (debug)
//     fprintf(stdout,"Writing BC-enforced state to u0bc_vec\n");
//   u->writeall("u0bc_vec",0);

//   // write rhs data state to file
//   if (debug)
//     fprintf(stdout,"Writing right-hand-side to rhs_vec\n");
//   rhorhs->writeall("rhs_vec",0);

  // iterate over interior, computing poisson residual
  fu->constant(0.0);
  float *udata = u->GetData(0);
//   float *rdata = rhorhs->GetData(0);
  float *fdata = fu->GetData(0);
  int i, j, k, idx, idxl00, idxr00, idx0l0, idx0r0, idx00l, idx00r;
  int udims[4], ugh[3][2];
  u->size(&udims[0], &udims[1], &udims[2], &udims[3], 
	  &ugh[0][0], &ugh[0][1], &ugh[1][0], 
	  &ugh[1][1], &ugh[2][0], &ugh[2][1]);
  int x0len = udims[0] + ugh[0][0] + ugh[0][1];
  int x1len = udims[1] + ugh[1][0] + ugh[1][1];
//   int x2len = udims[2] + ugh[2][0] + ugh[2][1];
  for (k=ugh[2][0]; k<udims[2]+ugh[2][0]; k++) 
    for (j=ugh[1][0]; j<udims[1]+ugh[1][0]; j++)
      for (i=ugh[0][0]; i<udims[0]+ugh[0][0]; i++) {
	idx = (k*x1len + j)*x0len + i;
	idxl00 = (k*x1len + j)*x0len + i - 1;
	idxr00 = (k*x1len + j)*x0len + i + 1;
	idx0l0 = (k*x1len + j - 1)*x0len + i;
	idx0r0 = (k*x1len + j + 1)*x0len + i;
	idx00l = ((k-1)*x1len + j)*x0len + i;
	idx00r = ((k+1)*x1len + j)*x0len + i;
	fdata[idx] = 
	  (udata[idxl00] - 2.0*udata[idx] + udata[idxr00])/dx[0]/dx[0] + 
	  (udata[idx0l0] - 2.0*udata[idx] + udata[idx0r0])/dx[1]/dx[1] + 
	  (udata[idx00l] - 2.0*udata[idx] + udata[idx00r])/dx[2]/dx[2];
      }

//   // write Laplace information to file
//   if (debug)
//     fprintf(stdout,"Writing local Laplace to laplace_vec\n");
//   fu->writeall("laplace_vec",0);

  // subtract off the rhs to obtain residual
  fu->axpy(-1.0,rhorhs);

//   // write local residual to file
//   if (debug)
//     fprintf(stdout,"Writing local residual to resid_vec\n");
//   fu->writeall("resid_vec",0);

//   // have residual communicate neighbor information
//   fu->exchange();

//   // write communicated residual to file
//   if (debug)
//     fprintf(stdout,"Writing communicated residual to residcomm_vec\n");
//   fu->writeall("residcomm_vec",0);

  // return success
  return SUCCESS;

}
#endif
