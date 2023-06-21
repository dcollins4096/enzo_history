/*****************************************************************************
 *                                                                           *
 * Copyright 2005 Daniel R. Reynolds                                         *
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRAVITY POTENTIAL BOUNDARY CLASS
/
/  written by: Daniel Reynolds
/  date:       September, 2005
/  modified1:  
/
/  PURPOSE: This class is designed to handle setup and information 
/           about the boundary conditions for the self-gravity 
/           potential field solve.  It is only used for potential field
/           solves on the root grid, since subgrid solves are handled 
/           using a hard-coded multigrid solver with boundary conditions 
/           arising from the root grid solve.
/
/           In this class, we allow for specification of either 
/           periodic or isolating boundary conditions on any face 
/           pair, for example one may set the x1 and x2 boundaries to 
/           be periodic, with the x3 boundary set with isolating 
/           conditions.  Isolating conditions correspond to prescribing
/           values on the potential field at the outer-most layer of 
/           volume cells (Dirichlet BCs).  These values may be provided 
/           in either scalar (all values on the layer are the same) or 
/           array form (the values on the layer may change spatially).
/
/           It may be safely assumed that root-grid boundary conditions
/           set in this class will supercede those set directly in the 
/           Grid class and TopGridData object, although all efforts will
/           be made to adjust those if possible when boundary conditions 
/           in this class are changed.
/
************************************************************************/

#ifndef GRAVITY_POTENTIAL_BOUNDARY_DEFINED__
#define GRAVITY_POTENTIAL_BOUNDARY_DEFINED__

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

class GravityPotentialBoundary
{
 private:
  bool ReadyToGo;               // Flag denoting whether internal
                                // structures are ready to be used

  int  BoundaryRank;            // This is the rank of the top grid */

  int  LocalDimension[3];       // Dims of the local Gravity grid

  int  GlobalDimension[3];      // Dims of the overall Gravity grid

  int  BoundaryType[3];         // Type of BC on face:
                                //    periodic   = 0
                                //    isolating != 0 
                                // (to match current scheme)

  float *BoundaryValues[3][2];  // boundary values for isolating

  FLOAT GravityLeftEdge[3];     // Left  edge of this grid's gravity domain
  FLOAT GravityRightEdge[3];    // Right edge of this grid's gravity domain

  int LeftEdgeIndex[3];         // Global index in gravity domain of left edge
  int RightEdgeIndex[3];        // Global index in gravity domain of right edge

  bool  OwnBoundary[3][2];      // logic denoting if proc on boundary

  float GravityCellSize;        // uniform cell size for gravity domain

 public:
//
// Constructor (set rank and dims from TopGrid)
//
  GravityPotentialBoundary();

//
// Checks if the Gravity boundary is ready to be used
//
  int AmIReadyToGo() {return ReadyToGo;};

//
// Set one face of gravity boundary values to either a constant value, or 
// set values over the entire face
//   Returns SUCCESS or FAIL.
//
  int SetGravityBoundaryValues(int Dimension, int Face, 
			       int BdryConst, float *BdryValue);
//
// Return the BoundaryType in a given dimension
//
  int GetBoundaryType(int Dimension) {
    return BoundaryType[Dimension];};
//
// Return the BoundaryValue array pointer
//
  float* GetBoundaryValues(int Dimension, int Face) {
    return BoundaryValues[Dimension][Face];};
//
// Set up the gravity domain extents, etc., based on the Gravity 
// Boundary Types and local Enzo domain
//
  int SetupGravityDomain(HierarchyEntry &TopGrid, TopGridData &MetaData);
//
// Get the Gravity-local boundary size in a given dimension
//
  int GetLocalDimension(int Dimension) {
    return LocalDimension[Dimension];};
//
// Get the left and right Gravity-local grid edges
//
  FLOAT GetBoundaryLeftEdge(int Dimension) {
    return GravityLeftEdge[Dimension];};
  FLOAT GetBoundaryRightEdge(int Dimension) {
    return GravityRightEdge[Dimension];};
//
// Get whether my processor owns piece of selected boundary
//
  bool AmIOnBoundary(int Dimension, int Face) {
    return OwnBoundary[Dimension][Face];};
//
// Get gravity domain cell size
//
  float GetCellSize() {
    return GravityCellSize;};

//
// Write/Read a potential boundry
//
  int WritePotentialBoundary (FILE *fptr, char *hdfname);  
  int ReadPotentialBoundary (FILE *fptr);  

};




#endif
