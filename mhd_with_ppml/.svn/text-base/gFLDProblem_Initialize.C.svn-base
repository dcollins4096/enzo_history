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
/  Gray Flux-Limited Diffusion Implicit Problem Class problem 
/  initialization routine
/
/  written by: Daniel Reynolds
/  date:       September, 2006
/  modified1:  
/
/  PURPOSE: Allocates all necessary internal memory for problem 
/           definition and associated linear solver.  This begins the
/           interface between Enzo and the FLD solver module, so any
/           and all grid/index transformations must be performed and 
/           stored here.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"


// character strings
EXTERN char outfilename[];


int gFLDProblem::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData)
{

//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::Initialize routine\n");

  // find the grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("FLD Initialize ERROR: p%"ISYM" could not locate his grid\n",
	   MyProcessorNumber);
    return FAIL;
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


  // set default module parameters
  Nchem  = 1;           // hydrogen only
  Model  = 0;           // no ec, Eg dependence on chemistry
  theta  = 1.0;         // backwards euler implicit time discret.
  LimImp = 0;           // lag implicit dependence of limiter in time
  BdryType[0][0] = 0;
  BdryType[0][1] = 0;   // set default radiation boundaries to 
  BdryType[1][0] = 0;   //   periodic in each direction
  BdryType[1][1] = 0;
  BdryType[2][0] = 0;
  BdryType[2][1] = 0;
  BoundaryFName = "radiation_bdry";

  // set default solver parameters
  newt_maxit         = 20;        // 20 Newton iterations
  newt_norm          = 0;         // standard RMS norm
  newt_INconst       = 1.0e-8;    // inexact Newton forcing constant
  newt_tol           = 1.0e-4;    // default nonlinear tolerance
  newt_MinLinesearch = 1.0e-10;   // minimum linesearch step length
  sol_relch          = 0;         // HYPRE relative change stopping crit.
  sol_printl         = 1;         // HYPRE print level
  sol_log            = 1;         // HYPRE logging level
  sol_zeroguess      = 0;         // HYPRE uses a zero initial guess
  sol_maxit          = 50;        // HYPRE max multigrid iters
  sol_rlxtype        = 1;         // HYPRE relaxation type
  sol_npre           = 1;         // HYPRE num pre-smoothing steps
  sol_npost          = 1;         // HYPRE num post-smoothing steps
  
  ////////////////////////////////
  // if input file present, over-write defaults with module inputs
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int ret;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  // check whether input file is non-null
  if (MetaData.RadHydroParameterFname != NULL) {
    if ((fptr = fopen(MetaData.RadHydroParameterFname, "r")) == NULL)
      fprintf(stderr,"Error opening RadHydro parameter file %s, using defaults\n",
	      MetaData.RadHydroParameterFname);
    else {

      // read until out of lines
      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	ret = 0;
	ret += sscanf(line, "RadHydroESpectrum = %"ISYM, &ESpectrum);
	ret += sscanf(line, "RadHydroChemistry = %"ISYM, &Nchem);
	ret += sscanf(line, "RadHydroModel = %"ISYM, &Model);
	ret += sscanf(line, "RadHydroErUnits = %"FSYM, &ErUnits);
	ret += sscanf(line, "RadHydroTheta = %"FSYM, &theta);
	ret += sscanf(line, "RadHydroImplicitLimiter = %"ISYM, &LimImp);
	ret += sscanf(line,"RadiationBoundaryX0Faces = %"ISYM" %"ISYM, 
		      BdryType[0], BdryType[0]+1);
	ret += sscanf(line,"RadiationBoundaryX1Faces = %"ISYM" %"ISYM,
		      BdryType[1], BdryType[1]+1);
	ret += sscanf(line,"RadiationBoundaryX2Faces = %"ISYM" %"ISYM,
		      BdryType[2], BdryType[2]+1);
	if (sscanf(line, "RadHydroBoundaryFName = %s", dummy) == 1)
	  BoundaryFName = dummy;
	ret += sscanf(line, "RadHydroNewtIters = %"ISYM, &newt_maxit);
	ret += sscanf(line, "RadHydroNewtNorm = %"ISYM, &newt_norm);
	ret += sscanf(line, "RadHydroINConst = %"FSYM, &newt_INconst);
	ret += sscanf(line, "RadHydroNewtTolerance = %"FSYM, &newt_tol);
	ret += sscanf(line, "RadHydroMinLinesearch = %"FSYM,
		      &newt_MinLinesearch);
	ret += sscanf(line, "RadHydroMaxMGIters = %i", &sol_maxit);
	ret += sscanf(line, "RadHydroMGRelaxType = %i", &sol_rlxtype);
	ret += sscanf(line, "RadHydroMGPreRelax = %i", &sol_npre);
	ret += sscanf(line, "RadHydroMGPostRelax = %i", &sol_npost);
	
      }  // end loop over file lines
    }  // end successfule file open
  }  // end if file name exists
 
  // clean up
  delete dummy;
  rewind(fptr);
  fclose(fptr);

  ////////////////////////////////

 
  // check that these give appropriate values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++)
      /// ADD NEW BOUNDARY CONDITION TYPES HERE!
      if ((BdryType[dim][face] < 0) || (BdryType[dim][face] > 2)) {
	fprintf(stderr,"gFLDProblem_Initialize Warning: re-setting BC to periodic, dim %"ISYM", face %"ISYM"\n",dim,face);
	BdryType[dim][face] = 0;
      }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim][0]*BdryType[dim][1] == 0) && 
	(BdryType[dim][0]+BdryType[dim][1] != 0)) {
      fprintf(stderr,"gFLDProblem_Initialize Warning: non-matching periodic BCs, dim %"ISYM"\n",dim);
      BdryType[dim][0] = 0;
      BdryType[dim][1] = 0;
    }

  // ensure that new EBdryVals array pointers are set to NULL
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      EBdryVals[dim][face] = NULL;


  // set up subdomain information
  //   EdgeVals gives the location of the left/right edge of the
  //      domain (no bdry) -- start with Enzo grid size
  for (dim=0; dim<rank; dim++) {
    EdgeVals[dim][0] = ThisGrid->GridData->GetGridLeftEdge(dim);
    EdgeVals[dim][1] = ThisGrid->GridData->GetGridRightEdge(dim);
  }

  //   LocDims holds the dimensions of the local domain: we 
  //      first obtain the Enzo grid size
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;

  // dx gives grid cell size
  for (dim=0; dim<rank; dim++)
    dx[dim] = (EdgeVals[dim][1]-EdgeVals[dim][0])/LocDims[dim];

  // dt gives the time step size (initialize to zero)
  dt = 0.0;

  // a, adot give cosmological expansion & rate at old/new time steps
  aold = 1.0;
  anew = 1.0;
  adotold = 0.0;
  adotnew = 0.0;

  // Nchem gives the number of chemical species
  if ((Nchem < 0) || (Nchem > 3)) {
    fprintf(stderr,"gFLDProblem Initialize: illegal Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"   re-setting Nchem to 1\n");
    Nchem = 1;  // default is Hydrogen only
  }

  // LimImp gives the implicitness of the radiation flux limiter (see header)
  if ((LimImp < 0) || (LimImp > 4)) {
    fprintf(stderr,"gFLDProblem Initialize: illegal LimImp = %"ISYM"\n",
	    LimImp);
    fprintf(stderr,"   re-setting LimImp to 0\n");
    LimImp = 0;  // default is time-lagged in implicit solve
  }

  // Theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"gFLDProblem Initialize: illegal theta = %g\n",
	    theta);
    fprintf(stderr,"   re-setting theta to 1.0 (Backwards Euler)\n");
    theta = 1.0;  // default is backwards Euler
  }

  // set flags denoting if this processor is on the external boundary
  for (dim=0; dim<rank; dim++) {
    if (layout[dim]==0) {
      OnBdry[dim][0] = OnBdry[dim][1] = true;
    }
    else {
      OnBdry[dim][0] = (location[dim] == 0);
      OnBdry[dim][1] = (location[dim] == layout[dim]-1);
    }
  }
  fprintf(stdout,"gFLDProblem::Initialize p%"ISYM": layout = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,layout[0],layout[1],layout[2]);
  fprintf(stdout,"gFLDProblem::Initialize p%"ISYM": location = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,location[0],location[1],location[2]);

  //   for non-periodic domain, unset neighbor info.
  for (dim=0; dim<rank; dim++) {
    if ((OnBdry[dim][0]) && (BdryType[dim][0] != 0))
      NBors[dim][0] = MPI_PROC_NULL;
    if ((OnBdry[dim][1]) && (BdryType[dim][1] != 0))
      NBors[dim][1] = MPI_PROC_NULL;
  }
  fprintf(stdout,"gFLDProblem::Initialize p%"ISYM": OnBdry = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,int(OnBdry[0][0]),int(OnBdry[0][1]),int(OnBdry[1][0]),int(OnBdry[1][1]),int(OnBdry[2][0]),int(OnBdry[2][1]));
  fprintf(stdout,"gFLDProblem::Initialize p%"ISYM": BdryType = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,BdryType[0][0],BdryType[0][1],BdryType[1][0],BdryType[1][1],BdryType[2][0],BdryType[2][1]);
  
  
  // compute global dimension information
  for (dim=0; dim<rank; dim++)
    GlobDims[dim] = MetaData.TopGridDims[dim];

  // compute global index information for this gravity subdomain
  float fCellsLeft;
  for (dim=0; dim<rank; dim++) {

    // the global indexing is easy if we're at the left edge
    if (location[dim]==0)  SolvIndices[dim][0]=0;

    // otherwise we compute the number of intervening cells to left edge
    else {

      // get floating point value for number of cells
      fCellsLeft = (EdgeVals[dim][0] - DomainLeftEdge[dim])/dx[dim];

      // round floating point value to closest integer
      SolvIndices[dim][0] = lround(fCellsLeft);

    }

    // add on local size to obtain right edge indices
    SolvIndices[dim][1] = SolvIndices[dim][0] + LocDims[dim] - 1;
  }  

  // store local array sizes (active + ghost)
  for (dim=0; dim<rank; dim++)
    ArrDims[dim] = LocDims[dim] + 2*DEFAULT_GHOST_ZONES;
 
  // set up vector container for previous time step (empty data)
  int ghosts = DEFAULT_GHOST_ZONES;
  int empty=1;
  U0 = new EnzoVector(LocDims[0], LocDims[1], LocDims[2], 
		      ghosts, ghosts, ghosts, ghosts, ghosts, ghosts, 2+Nchem, 
		      NBors[0][0], NBors[0][1], NBors[1][0], 
		      NBors[1][1], NBors[2][0], NBors[2][1], empty);
  GhDims[0][0] = ghosts;
  GhDims[0][1] = ghosts;
  GhDims[1][0] = ghosts;
  GhDims[1][1] = ghosts;
  GhDims[2][0] = ghosts;
  GhDims[2][1] = ghosts;

  // set up vectors for temporary storage and Jacobian components
  rhs  = U0->clone();
  rhs0 = U0->clone();
  Temp = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  sigmaA = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  sigmaS = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  L = (EnzoVector **) new EnzoVector*[2+Nchem];
  for (int i=0; i<(2+Nchem); i++)
    L[i] = U0->clone();


  // compute Radiation Energy spectrum integrals
  if (this->ComputeRadiationIntegrals() == FAIL) {
    fprintf(stderr,"gFLDProblem::Initialize Error in computing radiation spectrum integrals\n");
    return FAIL;
  }

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
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  HYPRE_SStructGridSetExtents(grid, 0, ilower, iupper);

  //       set grid variables for this part
  HYPRE_SStructVariable vartypes = HYPRE_SSTRUCT_VARIABLE_CELL;
  HYPRE_SStructGridSetVariables(grid, 0, 1, &vartypes);

  //       set grid periodicity
  Eint32 periodicity[3] = {0, 0, 0};
  if (BdryType[0][0] == 0)  periodicity[0] = GlobDims[0];
  if (BdryType[1][0] == 0)  periodicity[1] = GlobDims[1];
  if (BdryType[2][0] == 0)  periodicity[2] = GlobDims[2];
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


  //   check MG solver parameters
  if (sol_maxit < 0) {
    fprintf(stderr,"Illegal RadHydroMaxMGIters = %i. Setting to 20\n",
	    sol_maxit);
    sol_maxit = 20;
  }
  if ((sol_rlxtype<0) || (sol_rlxtype>3)) {
    fprintf(stderr,"Illegal RadHydroMGRelaxType = %i. Setting to 1\n",
	    sol_rlxtype);
    sol_rlxtype = 1;
  }
  if (sol_npre < 0) {
    fprintf(stderr,"Illegal RadHydroMGPreRelax = %i. Setting to 1\n",
	    sol_npre);
    sol_npre = 1;
  }
  if (sol_npost < 0) {
    fprintf(stderr,"Illegal RadHydroMGPostRelax = %i. Setting to 1\n",
	    sol_npost);
    sol_npost = 1;
  }

  //   check Newton solver parameters
  if (newt_maxit < 1) {
    fprintf(stderr,"Illegal RadHydroNewtIters = %"ISYM". Setting to 20\n",
	    newt_maxit);
    newt_maxit = 20;
  }
  if ((newt_norm < 0) || (newt_norm > 5)) {
    fprintf(stderr,"Illegal RadHydroNewtNorm = %"ISYM". Setting to 0\n",
	    newt_norm);
    newt_norm = 0;
  }
  if (newt_tol < 1.0e-15) {
    fprintf(stderr,"Illegal RadHydroNewtTolerance = %g. Setting to 1e-4\n",
	    newt_tol);
    newt_tol = 1.0e-4;
  }
  if ((newt_INconst < 1.0e-15) || (newt_INconst > 1.)) {
    fprintf(stderr,"Illegal RadHydroINConst = %g. Setting to 0.1\n",
	    newt_INconst);
    newt_INconst = 0.1;
  }
  if (newt_MinLinesearch < 1.0e-15) {
    fprintf(stderr,"Illegal RadHydroMinLinesearch = %g. Setting to 1e-12\n",
	    newt_MinLinesearch);
    newt_MinLinesearch = 1.0e-12;
  }


  // output RadHydro solver parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL) {
      fprintf(stderr,"Error opening parameter output file %s!!\n", 
	      outfilename);
      return FAIL;
    }
    else {
      fprintf(outfptr, "RadHydroESpectrum = %"ISYM"\n", ESpectrum);
      fprintf(outfptr, "RadHydroChemistry = %"ISYM"\n", Nchem);
      fprintf(outfptr, "RadHydroModel = %"ISYM"\n", Model);
      fprintf(outfptr, "RadHydroTheta = %"FSYM"\n", theta);
      fprintf(outfptr, "RadHydroImplicitLimiter = %"ISYM"\n", LimImp);
      fprintf(outfptr, "RadiationBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	      BdryType[0][0], BdryType[0][1]);
      fprintf(outfptr, "RadiationBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
	      BdryType[1][0], BdryType[1][1]);
      fprintf(outfptr, "RadiationBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
	      BdryType[2][0], BdryType[2][1]);
      fprintf(outfptr, "RadHydroBoundaryFName = %s\n", BoundaryFName);
      fprintf(outfptr, "RadHydroNewtIters = %"ISYM"\n", newt_maxit);    
      fprintf(outfptr, "RadHydroNewtNorm = %"ISYM"\n", newt_norm);    
      fprintf(outfptr, "RadHydroINConst = %"FSYM"\n", newt_INconst);    
      fprintf(outfptr, "RadHydroNewtTolerance = %"FSYM"\n", newt_tol);    
      fprintf(outfptr, "RadHydroMinLinesearch = %"FSYM"\n", 
	      newt_MinLinesearch);    
      fprintf(outfptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
      fprintf(outfptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
      fprintf(outfptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
      fprintf(outfptr, "RadHydroMGPostRelax = %i\n", sol_npost);    
      
      // close parameter file
      fclose(outfptr);
    }
  }

  return SUCCESS;
}
#endif
