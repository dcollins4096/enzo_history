/***********************************************************************
/
/  INITIALIZE PROTOSTELLAR CORE COLLAPSE
/
/  written by: Daniel R. Reynolds
/  date:       October, 2005
/  modified1:  
/
/  PURPOSE:    This file sets up problem-specific gravity potential 
/              boundary conditions for a process-local subgrid.  As
/              this routine assumes the root grid has already been 
/              distributed among different processors, it must be 
/              called *after* CommunicationDistributeGrid.  This 
/              routine is designed to be called after the standard 
/              problem initialization routine, here 
/              ProtostellarCollapseInitialize.C
/
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "GravityPotentialBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"


/* Set up the Gravitational Potential Boundary: in order to use 
   non-triply periodic boundary conditions for the Gravity potential, 
   we need to set up the problem-specific boundary conditions here.  
   We use the boundary condition choices from MetaData, and 
   then must manually set the potential boundary values for the 
   boundaries owned by this subgrid.
   
   ** Currently, isolating values are all set to zero, but code is **
   ** in place to allow for location-specific isolating BC values, **
   ** with location computations already performed.                **
*/
#ifdef ISO_GRAV
int ProtostellarCollapseInitializeGravity(TopGridData &MetaData, 
					  GravityPotentialBoundary &PotBdry)
{

  if (debug)
    fprintf(stdout,"Entering ProtostellarCollapseInitializeGravity\n");

  /* get gravity boundary types from MetaData */
  int GravityBoundaries[3];
  GravityBoundaries[0] = MetaData.GravityBoundaryFaces[0];
  GravityBoundaries[1] = MetaData.GravityBoundaryFaces[1];
  GravityBoundaries[2] = MetaData.GravityBoundaryFaces[2];
 

  /* get gravity subdomain information from PotBdry */
  
  /*   x*leftedge gives the location of the left edge of the domain */
  FLOAT x0leftedge = PotBdry.GetBoundaryLeftEdge(0);
  FLOAT x1leftedge = PotBdry.GetBoundaryLeftEdge(1);
  FLOAT x2leftedge = PotBdry.GetBoundaryLeftEdge(2);
  
  /*   x*rightedge gives the location of the right edge of the domain */
  FLOAT x0rightedge = PotBdry.GetBoundaryRightEdge(0);
  FLOAT x1rightedge = PotBdry.GetBoundaryRightEdge(1);
  FLOAT x2rightedge = PotBdry.GetBoundaryRightEdge(2);
  
  /*   OnX*LBdry, OnX*RBdry gives whether proc owns part of bdry */
  bool OnX0LBdry = PotBdry.AmIOnBoundary(0,0);
  bool OnX0RBdry = PotBdry.AmIOnBoundary(0,1);
  bool OnX1LBdry = PotBdry.AmIOnBoundary(1,0);
  bool OnX1RBdry = PotBdry.AmIOnBoundary(1,1);
  bool OnX2LBdry = PotBdry.AmIOnBoundary(2,0);
  bool OnX2RBdry = PotBdry.AmIOnBoundary(2,1);
  
  /*   Nx* holds the dimensions of the local MGMPI gravity domain */
  int Nx0 = PotBdry.GetLocalDimension(0);
  int Nx1 = PotBdry.GetLocalDimension(1);
  int Nx2 = PotBdry.GetLocalDimension(2);
  
  /* dxyz gives grid cell size */
  float dxyz = PotBdry.GetCellSize();
  
  
  /* put together the external boundary conditions */
  int size = Nx0 * Nx1 * Nx2;
  int x0facesize = size/Nx0;
  int x1facesize = size/Nx1;
  int x2facesize = size/Nx2;
  float ZERO = 0.0;
  int i0, i1, i2;

  /* x*center gives the location of the centroid of the overall domain */
  float x0center = (DomainLeftEdge[0] + DomainRightEdge[0])/2.0;
  float x1center = (DomainLeftEdge[1] + DomainRightEdge[1])/2.0;
  float x2center = (DomainLeftEdge[2] + DomainRightEdge[2])/2.0;

  /* in each of the following loops, x*val gives the evaluation location 
     for the isolating BC */
  float x0val, x1val, x2val;
  float offset=0.5*dxyz;

  /* set x0L gravity boundaries, if non-periodic and on this processor */
  if ((GravityBoundaries[0] != 0) && (OnX0LBdry)) {
      
    /* set constant-valued gravity potential boundary values */
    // PotBdry.SetGravityBoundaryValues(0, 0, 1, &ZERO);

    /* set array-valued gravity potential boundary values */
    x0val = x0leftedge + offset;
    float *Lface12 = new float[Nx1*Nx2];
    for (i2 = 0; i2 < Nx2; i2++) {
      x2val = x2leftedge + offset + dxyz*i2;
      for (i1 = 0; i1 < Nx1; i1++) {
	x1val = x1leftedge + offset + dxyz*i1;
	Lface12[i2*Nx1 + i1] = 0.0;
      }
    }
    PotBdry.SetGravityBoundaryValues(0, 0, 0, Lface12);
    delete [] Lface12;
  }


  /* set x0R gravity boundaries, if non-periodic and on this processor */
  if ((GravityBoundaries[0] != 0) && (OnX0RBdry)) {
      
    /* set constant-valued gravity potential boundary values */
    // PotBdry.SetGravityBoundaryValues(0, 1, 1, &ZERO);

    /* set array-valued gravity potential boundary values */
    x0val = x0leftedge + offset + dxyz*Nx0;
    float *Rface12 = new float[Nx1*Nx2];
    for (i2 = 0; i2 < Nx2; i2++) {
      x2val = x2leftedge + offset + dxyz*i2;
      for (i1 = 0; i1 < Nx1; i1++) {
	x1val = x1leftedge + offset + dxyz*i1;
	Rface12[i2*Nx1 + i1] = 0.0;
      }
    }
    PotBdry.SetGravityBoundaryValues(0, 1, 0, Rface12);
    delete [] Rface12;
  }


  /* set x1L gravity boundaries, if non-periodic and on this processor */
  if ((GravityBoundaries[1] != 0) && (OnX1LBdry)) {
      
    /* set constant-valued gravity potential boundary values */
    // PotBdry.SetGravityBoundaryValues(1, 0, 1, &ZERO);

    /* set array-valued gravity potential boundary values */
    x1val = x1leftedge + offset;
    float *Lface20 = new float[Nx2*Nx0];
    for (i0 = 0; i0 < Nx0; i0++) {
      x0val = x0leftedge + offset + dxyz*i0;
      for (i2 = 0; i2 < Nx2; i2++) {
	x2val = x2leftedge + offset + dxyz*i2;
	Lface20[i0*Nx2 + i2] = 0.0;
      }
    }
    PotBdry.SetGravityBoundaryValues(1, 0, 0, Lface20);
    delete [] Lface20;
  }


  /* set x1R gravity boundaries, if non-periodic and on this processor */
  if ((GravityBoundaries[1] != 0) && (OnX1RBdry)) {
      
    /* set constant-valued gravity potential boundary values */
    // PotBdry.SetGravityBoundaryValues(1, 1, 1, &ZERO);

    /* set array-valued gravity potential boundary values */
    x1val = x1leftedge + offset + dxyz*Nx1;
    float *Rface20 = new float[Nx2*Nx0];
    for (i0 = 0; i0 < Nx0; i0++) {
      x0val = x0leftedge + offset + dxyz*i0;
      for (i2 = 0; i2 < Nx2; i2++) {
	x2val = x2leftedge + offset + dxyz*i2;
	Rface20[i0*Nx2 + i2] = 0.0;
      }
    }
    PotBdry.SetGravityBoundaryValues(1, 1, 0, Rface20);
    delete [] Rface20;
  }


  /* set x2L gravity boundaries, if non-periodic and on this processor */
  if ((GravityBoundaries[2] != 0) && (OnX2LBdry)) {
      
    /* set constant-valued gravity potential boundary values */
    // PotBdry.SetGravityBoundaryValues(2, 0, 1, &ZERO);

    /* set array-valued gravity potential boundary values */
    x2val = x2leftedge + offset;
    float *Lface01 = new float[Nx0*Nx1];
    for (i1 = 0; i1 < Nx1; i1++) {
      x1val = x1leftedge + offset + dxyz*i1;
      for (i0 = 0; i0 < Nx0; i0++) {
	x0val = x0leftedge + offset + dxyz*i0;
	Lface01[i1*Nx0 + i0] = 0.0;
      }
    }
    PotBdry.SetGravityBoundaryValues(2, 0, 0, Lface01);
    delete [] Lface01;
  }


  /* set x2R gravity boundaries, if non-periodic and on this processor */
  if ((GravityBoundaries[2] != 0) && (OnX2RBdry)) {
      
    /* set constant-valued gravity potential boundary values */
    // PotBdry.SetGravityBoundaryValues(2, 1, 1, &ZERO);

    /* set array-valued gravity potential boundary values */
    x2val = x2leftedge + offset + dxyz*Nx2;
    float *Rface01 = new float[Nx0*Nx1];
    for (i1 = 0; i1 < Nx1; i1++) {
      x1val = x1leftedge + offset + dxyz*i1;
      for (i0 = 0; i0 < Nx0; i0++) {
	x0val = x0leftedge + offset + dxyz*i0;
	Rface01[i1*Nx0 + i0] = 0.0;
      }
    }
    PotBdry.SetGravityBoundaryValues(2, 1, 0, Rface01);
    delete [] Rface01;
  }
 
  if (debug)
    fprintf(stdout,"Leaving ProtostellarCollapseInitializeGravity\n");

  return SUCCESS;

}
#endif
