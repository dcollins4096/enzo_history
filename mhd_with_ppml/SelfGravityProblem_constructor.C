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
/  Self Gravity Problem Class constructor routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Allocates all necessary internal memory for problem 
/           definition and associated linear solver.  This begins the
/           interface between Enzo and the solver module, so any
/           and all grid/index transformations must be performed and 
/           stored here.
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#include "SelfGravityProblem.h"


SelfGravityProblem::SelfGravityProblem(HierarchyEntry &TopGrid,
				       TopGridData &MetaData)
{

  // initialize prepared flag to false
  prepared = false;
  
  // find the grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("SelfGravity constructor ERROR: p%"ISYM" could not locate his grid\n",
	   MyProcessorNumber);
  }

  // set rank of self-gravity problem to 3
  rank = 3;
  
  // get processor layout from Grid
  for (dim=0; dim<rank; dim++) 
    layout[dim] = ThisGrid->GridData->GetProcessorLayout(dim);
  
  // get processor location in MPI grid
  for (dim=0; dim<rank; dim++) 
    location[dim] = ThisGrid->GridData->GetProcessorLocation(dim);

  // get neighbor information from grid
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      NBors[dim][face] = ThisGrid->GridData->GetProcessorNeighbors(dim,face);

  // store gravity domain boundary conditions
  for (dim=0; dim<rank; dim++) 
    BdryType[dim] = MetaData.GravityBoundaryFaces[dim];

  // check that these give appropriate values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim] < 0) || (BdryType[dim] > 1))  BdryType[dim] = 0;

  // ensure that new BdryVals array pointers are set to NULL
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      BdryVals[dim][face] = NULL;


  // set up subdomain information
  //   GravEdges gives the location of the left/right edge of the
  //      gravity domain (no bdry) -- start with Enzo grid size
  for (dim=0; dim<rank; dim++) {
    GravEdges[dim][0] = ThisGrid->GridData->GetGridLeftEdge(dim);
    GravEdges[dim][1] = ThisGrid->GridData->GetGridRightEdge(dim);
  }

  //   LocDims holds the dimensions of the local gravity 
  //        domain: we first obtain the Enzo grid size
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;

  // dx gives grid cell size
  for (dim=0; dim<rank; dim++)
    dx[dim] = (GravEdges[dim][1]-GravEdges[dim][0])/LocDims[dim];

  // set flags denoting if this processor is on the external boundary
  for (dim=0; dim<rank; dim++) {
    OnBdry[dim][0] = (location[dim] == 0);
    OnBdry[dim][1] = (location[dim] == layout[dim]-1);
  }

  // set gravity ghost zones
  GravGhosts = GRAVITY_BUFFER_SIZE + 1;

  //   for iso. domain, add ghosts at boundary, and unset neighbor info.
  for (dim=0; dim<rank; dim++) {
    if ((OnBdry[dim][0]) && (BdryType[dim] == 1)) {
      LocDims[dim] += GravGhosts-1;
      GravEdges[dim][0] = GravEdges[dim][0] - (GravGhosts-1)*dx[dim];
      NBors[dim][0] = MPI_PROC_NULL;
    }
    if ((OnBdry[dim][1]) && (BdryType[dim] == 1)) {
      LocDims[dim] += GravGhosts-1;
      GravEdges[dim][1] = GravEdges[dim][1] + (GravGhosts-1)*dx[dim];
      NBors[dim][1] = MPI_PROC_NULL;
    }
  }
  
  // compute global dimension information via hydro size + iso. BC adjustment
  for (dim=0; dim<rank; dim++) {
    GlobDims[dim] = MetaData.TopGridDims[dim];
    if (BdryType[dim] == 1)  GlobDims[dim] += 2*(GravGhosts-1);
  }

  // compute global index information for this gravity subdomain
  float fCellsLeft;
  for (dim=0; dim<rank; dim++) {

    // the global indexing is easy if we're at the left edge
    if (location[dim]==0)  EdgeIndices[dim][0]=0;

    // otherwise we compute the number of intervening cells to left edge
    else {

      // get floating point value for number of cells
      fCellsLeft = (GravEdges[dim][0] - DomainLeftEdge[dim])/dx[dim];

      // round floating point value to closest integer
      EdgeIndices[dim][0] = lround(fCellsLeft);

      // add padding if isolating (gravity domain larger than hydro)
      if (BdryType[dim]==1)  EdgeIndices[dim][0] += GravGhosts-1;
    }

    // add on local size to obtain right edge indices
    EdgeIndices[dim][1] = EdgeIndices[dim][0] + LocDims[dim] - 1;
  }  

  // determine this proc's left/right ghost zones
  int NGhosts[3][2];
  for (dim=0; dim<rank; dim++)
    for (face=0; face<2; face++) {
      if ((OnBdry[dim][face]) && (BdryType[dim] == 1))  NGhosts[dim][face] = 1;
      else NGhosts[dim][face] = GravGhosts;
    }

  // store local array sizes (active + ghost)
  for (dim=0; dim<rank; dim++)
    ArrDims[dim] = LocDims[dim] + NGhosts[dim][0] + NGhosts[dim][1];
 
  // set up vector container for gravitating mass field (NULL data pointer)
  float *ndata = NULL;
  rhorhs = new EnzoVector(LocDims[0], LocDims[1], LocDims[2], 
			  NGhosts[0][0], NGhosts[0][1], NGhosts[1][0], 
			  NGhosts[1][1], NGhosts[2][0], NGhosts[2][1], 1, 
			  NBors[0][0], NBors[0][1], NBors[1][0], 
			  NBors[1][1], NBors[2][0], NBors[2][1], &ndata);


  // initialize HYPRE stuff
  //    set the matrix type: HYPRE_SSTRUCT or HYPRE_PARCSR
  mattype = HYPRE_SSTRUCT;

  //    initialize the diagnostic information
  totIters = 0;

  //    set up the grid
  //       create the grid object
  HYPRE_SStructGridCreate(MPI_COMM_WORLD, 3, 1, &grid);

  //       set my grid extents as if we have one part with multiple boxes.
  //       Have each processor describe it's own global extents
  for (int dim=0; dim<3; dim++) {
    SolvIndices[dim][0] = EdgeIndices[dim][0];
    SolvIndices[dim][1] = EdgeIndices[dim][1];
    if (OnBdry[dim][0] && (BdryType[dim]==1))  SolvIndices[dim][0] -= 1;
    if (OnBdry[dim][1] && (BdryType[dim]==1))  SolvIndices[dim][1] += 1;
  }
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  HYPRE_SStructGridSetExtents(grid, 0, ilower, iupper);

  //       set grid variables for this part
  HYPRE_SStructVariable vartypes = HYPRE_SSTRUCT_VARIABLE_CELL;
  HYPRE_SStructGridSetVariables(grid, 0, 1, &vartypes);

  //       set grid periodicity
  Eint32 periodicity[3] = {0, 0, 0};
  if (BdryType[0] == 0)  periodicity[0] = GlobDims[0];
  if (BdryType[1] == 0)  periodicity[1] = GlobDims[1];
  if (BdryType[2] == 0)  periodicity[2] = GlobDims[2];
  HYPRE_SStructGridSetPeriodic(grid, 0, periodicity);
  
  //       assemble the grid
  HYPRE_SStructGridAssemble(grid);

  //   set up the stencil
  stSize = 7;
  HYPRE_SStructStencilCreate(3, stSize, &stencil);

  //      set stencil entries
  Eint32 offset[3];
  //         dependency to x2 left
  offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(stencil, 0, offset, 0);
  //         dependency to x1 left
  offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(stencil, 1, offset, 0);
  //         dependency to x0 left
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(stencil, 2, offset, 0);
  //         dependency to self
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(stencil, 3, offset, 0);
  //         dependency to x0 right
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(stencil, 4, offset, 0);
  //         dependency to x1 right
  offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(stencil, 5, offset, 0);
  //         dependency to x2 right
  offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(stencil, 6, offset, 0);

  //   set up the graph
  //      create the graph object
  HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);

  //      set graph type according to solver desired
  HYPRE_SStructGraphSetObjectType(graph, mattype);

  //      set stencils into graph
  HYPRE_SStructGraphSetStencil(graph, 0, 0, stencil);

  //      add any additional non-stencil entries into graph
  //      (none that I can think of)

  //      assemble the graph
  HYPRE_SStructGraphAssemble(graph);

  //   set matrix initialization flag to fail
  AInit = 0;

  //   set solver defaults
  sol_maxit     = 50;
  sol_relch     = 0;
  sol_rlxtype   = 1;
  sol_npre      = 1;
  sol_npost     = 1;
  sol_printl    = 1;
  sol_log       = 1;
  sol_zeroguess = 0;


  // store SelfGravityProblem preparedness status
  prepared = true;

}
#endif
