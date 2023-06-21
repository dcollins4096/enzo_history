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
/  Gray Flux-Limited Diffusion Implicit Problem Class 
/  EnforceBoundary routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Enforces boundary conditions on a FLD problem 
/           vector.  Depending on 'flag' this routine will perform
/           one of two BC-related actions:
/
/              flag=0: sets Dirichlet or Neumann values into a given 
/              vector.  Useful for enforcing conditions on a solution.
/
/              flag!=0: set zero-valued homogeneous Dirichlet or 
/              Neumann conditions on any external face.  Useful
/              for enforcing conditions into a Newton update, 
/              which should not interfere with BC values.
/
/           Note: Neumann values are enforced on the first layer of 
/           ghost zones using a second-order central difference 
/           approximation to the first (outward-normal) derivative.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"



int gFLDProblem::EnforceBoundary(EnzoVector *u, int flag)
{
  
//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::EnforceBoundary routine\n");

  // check that the gFLDProblem has been prepared
  if (!prepared) {
    fprintf(stderr,"EnforceBoundary ERROR: gFLDProblem unprepared\n");
    return FAIL;
  }
  
  // get information about the vector u, and check against BC dims
  int i, i2, j, j2, k, k2, idx, idx2, idxbc;
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
  if (udims[3] != (2+Nchem)) {
    fprintf(stderr,"p%"ISYM" EnforceBC: mismatched nspecies %"ISYM"!=3\n",
	    MyProcessorNumber,udims[3]);
    return FAIL;
  }

  // set some shortcuts for the EnzoVector dimensions
  int x0len = udims[0] + ugh[0][0] + ugh[0][1];
  int x1len = udims[1] + ugh[1][0] + ugh[1][1];

  // flag == 0: set BCs into udata
  float *udata;
  if (flag == 0) {

    //////////////////////////////////////////////
    // first do the radiation energy term
    udata = u->GetData(0);

    // x0 left boundary
    //   Dirichlet
    if (OnBdry[0][0] && (BdryType[0][0]==1)) {
//       i = -1;
      i = 0;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = EBdryVals[0][0][idxbc];
	}
    }
    //   Neumann
    if (OnBdry[0][0] && (BdryType[0][0]==2)) {
      i = -1;  i2 = i+1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dx[0]*EBdryVals[0][0][idxbc];
	}
    }

    // x0 right boundary
    //   Dirichlet
    if (OnBdry[0][1] && (BdryType[0][1]==1)) {
//       i = LocDims[0];
      i = LocDims[0]-1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = EBdryVals[0][1][idxbc];
	}
    }
    //   Neumann
    if (OnBdry[0][1] && (BdryType[0][1]==2)) {
      i = LocDims[0];  i2 = i-1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dx[0]*EBdryVals[0][1][idxbc];
	}
    }

    // x1 left boundary
    //   Dirichlet
    if (OnBdry[1][0] && (BdryType[1][0]==1)) {
//       j = -1;
      j = 0;
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = EBdryVals[1][0][idxbc];
	}
    }
    //   Neumann
    if (OnBdry[1][0] && (BdryType[1][0]==2)) {
      j = -1;  j2 = j+1;
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = udata[idx2] + dx[1]*EBdryVals[1][0][idxbc];
	}
    }

    // x1 right boundary
    //   Dirichlet
    if (OnBdry[1][1] && (BdryType[1][1]==1)) {
//       j = LocDims[1];
      j = LocDims[1]-1;
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = EBdryVals[1][1][idxbc];
	}
    }
    //   Neumann
    if (OnBdry[1][1] && (BdryType[1][1]==2)) {
      j = LocDims[1];  j2 = j-1;
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = udata[idx2] + dx[1]*EBdryVals[1][1][idxbc];
	}
    }

    // x2 left boundary
    //   Dirichlet
    if (OnBdry[2][0] && (BdryType[2][0]==1)) {
//       k = -1;
      k = 0;
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = EBdryVals[2][0][idxbc];
	}
    }
    //   Neumann
    if (OnBdry[2][0] && (BdryType[2][0]==2)) {
      k = -1;  k2 = k+1;
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = udata[idx2] + dx[2]*EBdryVals[2][0][idxbc];
	}
    }

    // x2 right boundary
    //   Dirichlet
    if (OnBdry[2][1] && (BdryType[2][1]==1)) {
//       k = LocDims[2];
      k = LocDims[2]-1;
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = EBdryVals[2][1][idxbc];
	}
    }
    //   Neumann
    if (OnBdry[2][1] && (BdryType[2][1]==2)) {
      k = LocDims[2];  k2 = k-1;
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = udata[idx2] + dx[2]*EBdryVals[2][1][idxbc];
	}
    }


    //////////////////////////////////////////////
    // second do the fluid energy correction term
    udata = u->GetData(1);

    /// GET THESE FROM ENZO
    // [I believe that these need not be enforced, since the implicit problem
    //  is purely local for this variable, so it is well-posed without any BCs]


    //////////////////////////////////////////////
    // last do the chemical species variables
    for (int ns=1; ns<=Nchem; ns++) {      
      udata = u->GetData(1+ns);

      /// GET THESE FROM ENZO
      // [I believe that these need not be enforced, since the implicit problem
      //  is purely local for this variable, so it is well-posed without any BCs]

    }
  }
  // flag != 0: enforce zero values at Dirichlet faces
  else {

    //////////////////////////////////////////////
    // first do the radiation energy term
    udata = u->GetData(0);

    // x0 left boundary
    if (OnBdry[0][0] && (BdryType[0][0]==1)) {
//       i = ugh[0][0]-1;
      i = ugh[0][0];
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x0 right boundary
    if (OnBdry[0][1] && (BdryType[0][1]==1)) {
//       i = LocDims[0]+ugh[0][0];
      i = LocDims[0]+ugh[0][0]-1;
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x1 left boundary
    if (OnBdry[1][0] && (BdryType[1][0]==1)) {
//       j = ugh[1][0]-1;
      j = ugh[1][0];
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x1 right boundary
    if (OnBdry[1][1] && (BdryType[1][1]==1)) {
//       j = LocDims[1]+ugh[1][0];
      j = LocDims[1]+ugh[1][0]-1;
      for (k=ugh[2][0]-1; k<LocDims[2]+ugh[2][0]+1; k++)
	for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x2 left boundary
    if (OnBdry[2][0] && (BdryType[2][0]==1)) {
//       k = ugh[2][0]-1;
      k = ugh[2][0];
      for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }

    // x2 right boundary
    if (OnBdry[2][1] && (BdryType[2][1]==1)) {
//       k = LocDims[2]+ugh[2][0];
      k = LocDims[2]+ugh[2][0]-1;
      for (j=ugh[1][0]-1; j<LocDims[1]+ugh[1][0]+1; j++)
	for (i=ugh[0][0]-1; i<LocDims[0]+ugh[0][0]+1; i++)
	  udata[(k*x1len + j)*x0len + i] = 0.0;
    }



    //////////////////////////////////////////////
    // second do the fluid energy correction term
    udata = u->GetData(1);

    /// GET THESE FROM ENZO
    // [I believe that these need not be enforced, since the implicit problem
    //  is purely local for this variable, so it is well-posed without any BCs]



    //////////////////////////////////////////////
    // last do the chemical species variables
    for (int ns=0; ns<Nchem; ns++) {      
      udata = u->GetData(2+ns);

      /// GET THESE FROM ENZO
      // [I believe that these need not be enforced, since the implicit problem
      //  is purely local for this variable, so it is well-posed without any BCs]

    }
  }

  // return success
  return SUCCESS;

}
#endif
