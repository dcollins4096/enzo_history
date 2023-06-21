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
/  Self-Gravity Problem Boundary Condition Enforcement Routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Enforces boundary conditions on a Self-Gravity problem 
/           vector.  Depending on 'flag' this routine will perform
/           one of two BC-related actions:
/
/              flag=0: sets Dirichlet values into a given vector.    
/              Useful for enforcing conditions on a solution.
/
/              flag!=0: set zero-valued homogeneous dirichlet 
/              conditions on any external Dirichlet face.  Useful
/              for enforcing conditions into a Newton update, 
/              which should not interfere with Dirichlet values.
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#include "SelfGravityProblem.h"



int SelfGravityProblem::EnforceBoundary(EnzoVector *u, int flag)
{
  
  // check that the SelfGravityProblem has been prepared
  if (!prepared) {
    fprintf(stderr,"EnforceBoundary ERROR: SelfGravityProblem unprepared\n");
    return FAIL;
  }
  
  // get information about the vector u, and check against BC dims
  float *udata = u->GetData(0);
  int i, j, k, idx, idxbc;
  int udims[4], ugh[3][2];
  u->size(&udims[0], &udims[1], &udims[2], &udims[3], 
	  &ugh[0][0], &ugh[0][1], &ugh[1][0], 
	  &ugh[1][1], &ugh[2][0], &ugh[2][1]);
  if (udims[0] != LocDims[0]) {
    fprintf(stderr,"p%"ISYM" EnforceBC: mismatched x0 dims %"ISYM"!=%"ISYM"\n",
	    MyProcessorNumber,udims[0],LocDims[0]);
    return FAIL;
  }
  if (udims[1] != LocDims[1]) {
    fprintf(stderr,"p%"ISYM" EnforceBC: mismatched x1 dims %"ISYM"!=%"ISYM"\n",
	    MyProcessorNumber,udims[1],LocDims[1]);
    return FAIL;
  }
  if (udims[2] != LocDims[2]) {
    fprintf(stderr,"p%"ISYM" EnforceBC: mismatched x2 dims %"ISYM"!=%"ISYM"\n",
	    MyProcessorNumber,udims[2],LocDims[2]);
    return FAIL;
  }

  // set some shortcuts for the EnzoVector dimensions
  int x0len = udims[0] + ugh[0][0] + ugh[0][1];
  int x1len = udims[1] + ugh[1][0] + ugh[1][1];

  // flag == 0: set Dirichlet values into udata
  if (flag == 0) {

    // x0 left boundary
    if (OnBdry[0][0] && (BdryType[0]==1)) {
      i = -1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = BdryVals[0][0][idxbc];
	}
    }

    // x0 right boundary
    if (OnBdry[0][1] && (BdryType[0]==1)) {
      i = LocDims[0];
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = BdryVals[0][1][idxbc];
	}
    }

    // x1 left boundary
    if (OnBdry[1][0] && (BdryType[1]==1)) {
      j = -1;
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = BdryVals[1][0][idxbc];
	}
    }

    // x1 right boundary
    if (OnBdry[1][1] && (BdryType[1]==1)) {
      j = LocDims[1];
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = BdryVals[1][1][idxbc];
	}
    }

    // x2 left boundary
    if (OnBdry[2][0] && (BdryType[2]==1)) {
      k = -1;
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = BdryVals[2][0][idxbc];
	}
    }

    // x2 right boundary
    if (OnBdry[2][1] && (BdryType[2]==1)) {
      k = LocDims[2];
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = BdryVals[2][1][idxbc];
	}
    }
  }
  // flag != 0: enforce zero Dirichlet values
  else {
    
    // x0 left boundary
    if (OnBdry[0][0] && (BdryType[0]==1)) {
      i = ugh[0][0]-1;
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x0 right boundary
    if (OnBdry[0][1] && (BdryType[0]==1)) {
      i = LocDims[0]+ugh[0][0];
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x1 left boundary
    if (OnBdry[1][0] && (BdryType[1]==1)) {
      j = ugh[1][0]-1;
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x1 right boundary
    if (OnBdry[1][1] && (BdryType[1]==1)) {
      j = LocDims[1]+ugh[1][0];
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x2 left boundary
    if (OnBdry[2][0] && (BdryType[2]==1)) {
      k = ugh[2][0]-1;
      for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x2 right boundary
    if (OnBdry[2][1] && (BdryType[2]==1)) {
      k = LocDims[2]+ugh[2][0];
      for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }
  }

  // return success
  return SUCCESS;

}
#endif
