/***********************************************************************
/
/  COMPUTE THE POTENTIAL FIELD
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:  Daniel R. Reynolds
/  date:       February, 2006
/  modified2:
/
/  PURPOSE:
/
************************************************************************/
/***********************************************************************
/
/  COMPUTE THE GRAVITATIONAL POTENTIAL FIELD
/
/  This file extends Greg Bryan's Original code that uses an FFT-based
/  algorithm to solve for the gravitational potential under periodic 
/  boundary conditions on the root (level 0) grid.
/
/  Additional functionality has been added to allow for isolating 
/  (Dirichlet) boundary conditions on the root grid.  This solve calls
/  the MGMPI library for solution of the poisson equation.
/
/  NOTE: both approaches compute and store relevant solver information 
/  during the first call to the routine.  Neither of these routines 
/  'clean up' after themselves upon exit of the program, i.e. memory 
/  is allocated but never freed, requiring that the compiler take care 
/  of the remaining data upon program completion.  A future version of 
/  this gravity solver module may include a C++ class for the solver, 
/  which is initialized at the same point as the Enzo grids, stores 
/  all necessary information internally and privately, and is cleared 
/  prior to exit of the overall Enzo program.
/
/  written by: Greg Bryan
/  date:       January, 1998
/
/  modified1:  Daniel R. Reynolds
/  date:       August, 2005
/
************************************************************************/
#ifdef ISO_GRAV
/* New includes for MGMPI */
#include "mgmpi.h"
#endif

/* Original includes for Enzo */
#include <stdio.h>
#include <mpi.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#ifdef ISO_GRAV
#include "GravityPotentialBoundary.h"
#endif
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
 

/* Function prototypes */
int CommunicationParallelFFT(region *InRegion, int NumberOfInRegions,
			     region **OutRegion, int *NumberOfOutRegions,
			     int DomainDim[], int Rank,
			     int direction, int TransposeOnCompletion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

#ifdef ISO_GRAV
int ComputePotentialFieldLevelZeroIso(TopGridData *MetaData,
				      HierarchyEntry *Grids[], int NumberOfGrids, 
				      GravityPotentialBoundary *PotBdry);
#endif

int ComputePotentialFieldLevelZeroPer(TopGridData *MetaData, 
				      HierarchyEntry *Grids[], int NumberOfGrids);

int MGMPIBdryExchange(float *mgmpivec, int *layout, int *location, 
		      int x0l, int x0r, int x1l, int x1r, int x2l, 
		      int x2r, int mdim0, int mdim1, int mdim2);

int GravityBdryExchange(float *mgmpivec, int *layout, int *location,
			int GravityGhosts, int edim0, int edim1, 
			int edim2, int x0l, int x0r, int x1l, int x1r, 
			int x2l, int x2r);

extern "C" void FORTRAN_NAME(mgmpi_to_enzo)(float *mgmpivec, float *enzovec, 
					    int *mdim1, int *mdim2, int *mdim3, 
					    int *x1Li, int *x1Ri, int *x2Li, 
					    int *x2Ri, int *x3Li, int *x3Ri,
					    int *edim1, int *edim2, int *edim3, 
					    int *es1, int *es2, int *es3, 
					    int *ef1, int *ef2, int *ef3, 
					    int *off1, int *off2, int *off3); 

extern "C" void FORTRAN_NAME(enzo_to_mgmpi)(float *enzovec, float *mgmpivec, 
					    int *edim1, int *edim2, int *edim3, 
					    int *mdim1, int *mdim2, int *mdim3, 
					    int *es1, int *es2, int *es3, 
					    int *ef1, int *ef2, int *ef3, 
					    int *off1, int *off2, int *off3);



/******************************************************************/
/* ComputePotentialFieldLevelZero observes teh root-grid boundary */
/* conditions and calls the corresponding potential field solver. */
/******************************************************************/
#ifdef ISO_GRAV
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids,
				   GravityPotentialBoundary *PotBdry)
{

  /* call the appropriate solver depending on the
     desired root grid boundary conditions */
  if (MetaData->GravityBoundary == TopGridPeriodic) {

    if (debug)
      fprintf(stdout, "ComputePotentialFieldLevelZero: TopGridPeriodic \n");
    /* Periodic root grid BC's */
    if (ComputePotentialFieldLevelZeroPer(MetaData, Grids, 
					  NumberOfGrids) == FAIL) {
      fprintf(stderr, "Error in ComputePotentialFieldLevelZeroPer.\n");
      return FAIL;
    }    
  }

  else if (MetaData->GravityBoundary == TopGridIsolated) {
    
    if (debug)
      fprintf(stdout, "ComputePotentialFieldLevelZero: TopGridIsolated \n");
    if (MetaData->TopGridRank == 3) {
      /* Isolating root grid BC's */
      if (ComputePotentialFieldLevelZeroIso(MetaData, Grids, NumberOfGrids,
					    PotBdry) == FAIL) {
	fprintf(stderr, "Error in ComputePotentialFieldLevelZeroIso.\n");
	return FAIL;
      }
    }
    else {
      fprintf(stderr, "Error in ComputePotentialFieldLevelZero: \n");
      fprintf(stderr, "  Isolating BC's are only allowd for 3D problems.\n");
    }
  }

  else {

    /* Unknown root grid BC's */
    fprintf(stderr, "Error in ComputePotentialFieldLevelZero: \n");
    fprintf(stderr, "  only Periodic and Isolating BC's are allowed.\n");
    return FAIL;
  }  // end: if (TopGridPeriodic)
  
  return SUCCESS;
}
#else
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], 
				   int NumberOfGrids)
{
  /* call the periodic solver */
  if (ComputePotentialFieldLevelZeroPer(MetaData, Grids, 
					NumberOfGrids) == FAIL) {
    fprintf(stderr, "Error in ComputePotentialFieldLevelZeroPer.\n");
    return FAIL;
  }    
  return SUCCESS;
}
#endif

  
#ifdef ISO_GRAV
/******************************************************************/
/*  ComputePotentialFieldLevelZeroIso performs a root-grid        */
/*  potential field solver using isolating (Dirichlet) boundary   */
/*  conditions, through calls to the external MGMPI library.      */
/******************************************************************/
int ComputePotentialFieldLevelZeroIso(TopGridData *MetaData,
				      HierarchyEntry *Grids[], 
				      int NumberOfGrids,
				      GravityPotentialBoundary *PotentialBdry)
{
  /* Static declarations (for MGMPI setup) */
  static int FirstCall = TRUE;
  static Eint32 MGparams;     /* Handle to parameters object */
  static Eint32 MGconcurr;    /* Handle to concurrent object */
  static Eint32 MGdistrib;    /* Handle to the data distribute object */
  static Eint32 MGgrid;       /* Handle to the grid object */
  static Eint32 MGbc;         /* Handle to the boundary conditions object */
  static Eint32 MGa;          /* Handle to the coefficient matrix object */
  static Eint32 MGsolver;     /* Handle to the linear solver object */
  static Eint32 MGstopping;   /* Handle to the stopping criteria object */
  static Eint32 MGx;          /* Handle to the solution vector object */
  static Eint32 MGb;          /* Handle to the RHS vector object */

  /* Declarations */
  int i, j, k;
  int mygrid = -1;

  /* Figure out which grid I own (as this MPI process) */
  for (i=0; i<NumberOfGrids; i++) {
    if (MyProcessorNumber == Grids[i]->GridData->ReturnProcessorNumber()) {
      mygrid = i;
      break;
    }
  }
  if (mygrid == -1) {
    fprintf(stderr, "ERROR: process %"ISYM" could not locate its grid.\n",
	    MyProcessorNumber);
    return FAIL;
  }


#ifdef USE_MPI
  //  check that MyProcessorNumber agrees with MPI process ID
  MPI_Arg MPI_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_id);
  if (MyProcessorNumber != MPI_id) {
    fprintf(stderr, "ERROR: Enzo processor ID %"ISYM" does not match MPI ID %"ISYM"\n", 
	   MyProcessorNumber, MPI_id);
    return FAIL;
  }
#endif


  /* Ensure that PotentialBdry has been set up already */
  if (PotentialBdry->AmIReadyToGo() == FALSE) {
    fprintf(stderr,"ERROR: GravityPotentialBoundary not yet set up\n");
    return FAIL;
  }
  
  /* get process parallelism information from grid */
  /*   layout gives Nprocs in each dim (1-based) */
  /*   location gives this proc's index in each dim (0-based) */
  /*   per* gives periodicity of each dimension (1=>periodic) */
  int layout[3], location[3];
  layout[0] = Grids[mygrid]->GridData->GetProcessorLayout(0);
  layout[1] = Grids[mygrid]->GridData->GetProcessorLayout(1);
  layout[2] = Grids[mygrid]->GridData->GetProcessorLayout(2);
  location[0] = Grids[mygrid]->GridData->GetProcessorLocation(0);
  location[1] = Grids[mygrid]->GridData->GetProcessorLocation(1);
  location[2] = Grids[mygrid]->GridData->GetProcessorLocation(2);
  int ZERO = 0;
  int per0 = (PotentialBdry->GetBoundaryType(0) == ZERO) ? 1 : 0;
  int per1 = (PotentialBdry->GetBoundaryType(1) == ZERO) ? 1 : 0;
  int per2 = (PotentialBdry->GetBoundaryType(2) == ZERO) ? 1 : 0;
  
  /*   store whether this proc is on a given external boundary */
  bool x0LBdry = PotentialBdry->AmIOnBoundary(0,0);
  bool x0RBdry = PotentialBdry->AmIOnBoundary(0,1);
  bool x1LBdry = PotentialBdry->AmIOnBoundary(1,0);
  bool x1RBdry = PotentialBdry->AmIOnBoundary(1,1);
  bool x2LBdry = PotentialBdry->AmIOnBoundary(2,0);
  bool x2RBdry = PotentialBdry->AmIOnBoundary(2,1);

  /* get gravity solve local grid dimensions */
  int x0loc = PotentialBdry->GetLocalDimension(0);
  int x1loc = PotentialBdry->GetLocalDimension(1);
  int x2loc = PotentialBdry->GetLocalDimension(2);

  /* the 'local' size for MGMPI doesn't include Dirichlet bdry nodes
     (those constitute padding).  Thus, if proc owns a piece of the 
     bdry for ISOLATING dims, reduce size by 1 and note buffer size. */
  int MGx0loc, MGx1loc, MGx2loc, MGx0buf, MGx1buf, MGx2buf, MGStartIndex;
  MGx0buf = ((per0 == 0) && x0LBdry) ? 1 : 0;
  MGx0loc = x0loc - MGx0buf;
  MGx1buf = ((per1 == 0) && x1LBdry) ? 1 : 0;
  MGx1loc = x1loc - MGx1buf;
  MGx2buf = ((per2 == 0) && x2LBdry) ? 1 : 0;
  MGx2loc = x2loc - MGx2buf;
  if ((per0 == 0) && x0RBdry)  MGx0loc-=1;
  if ((per1 == 0) && x1RBdry)  MGx1loc-=1;
  if ((per2 == 0) && x2RBdry)  MGx2loc-=1;
  MGStartIndex = (MGx2buf*(x1loc+2) + MGx1buf)*(x0loc+2) + MGx0buf;
  
  /* get enzo local grid dimensions (gravitating_mass_field, etc.) */
  int edim0 = Grids[mygrid]->GridData->GetGravitatingMassFieldDimension(0);
  int edim1 = Grids[mygrid]->GridData->GetGravitatingMassFieldDimension(1);
  int edim2 = Grids[mygrid]->GridData->GetGravitatingMassFieldDimension(2);
  
  /* set MGMPI grid dimensions */
  int mdim0, mdim1, mdim2;
  mdim0 = x0loc+2;  mdim1 = x1loc+2;  mdim2 = x2loc+2;

  /* get process parallelism neighbor information from grid */
  /* (these assume periodicity, only for the added info,    */
  /*  and return the grid index and not the process index)  */
  int x0Nb[2], x1Nb[2], x2Nb[2];
  x0Nb[0] = Grids[mygrid]->GridData->GetProcessorNeighbors(0,0);
  x0Nb[1] = Grids[mygrid]->GridData->GetProcessorNeighbors(0,1);
  x1Nb[0] = Grids[mygrid]->GridData->GetProcessorNeighbors(1,0);
  x1Nb[1] = Grids[mygrid]->GridData->GetProcessorNeighbors(1,1);
  x2Nb[0] = Grids[mygrid]->GridData->GetProcessorNeighbors(2,0);
  x2Nb[1] = Grids[mygrid]->GridData->GetProcessorNeighbors(2,1);

//   /*   translate neighbor grid IDs to MPI process IDs */
//   if (debug)
//     fprintf(stdout,"CompPotFieldL0Iso: getting neighbor proc IDs\n");
//   x0Nb[0] = Grids[x0Nb[0]]->GridData->ReturnProcessorNumber();
//   x0Nb[1] = Grids[x0Nb[1]]->GridData->ReturnProcessorNumber();
//   x1Nb[0] = Grids[x1Nb[0]]->GridData->ReturnProcessorNumber();
//   x1Nb[1] = Grids[x1Nb[1]]->GridData->ReturnProcessorNumber();
//   x2Nb[0] = Grids[x2Nb[0]]->GridData->ReturnProcessorNumber();
//   x2Nb[1] = Grids[x2Nb[1]]->GridData->ReturnProcessorNumber();
  
  /*   if proc on boundary for isolating dimension, set neighbor to null */
//   if (debug)
//     fprintf(stdout,"CompPotFieldL0Iso: resetting exterior neighbor IDs\n");
  if (x0LBdry && (per0==0))  x0Nb[0] = MPI_PROC_NULL;
  if (x1LBdry && (per1==0))  x1Nb[0] = MPI_PROC_NULL;
  if (x2LBdry && (per2==0))  x2Nb[0] = MPI_PROC_NULL;
  if (x0RBdry && (per0==0))  x0Nb[1] = MPI_PROC_NULL;
  if (x1RBdry && (per1==0))  x1Nb[1] = MPI_PROC_NULL;
  if (x2RBdry && (per2==0))  x2Nb[1] = MPI_PROC_NULL;


  /* ------------------------------------------------------------------- */
  /* If this is the first time this routine has been called, then set up */
  /* the required MGMPI data structures describing the problem.          */
  if (FirstCall) {

    /* allocate and initialize the PotentialField in mygrid */
    /* (this will be re-used throughout the simulation)     */
    float *PotentialField = new float[edim0*edim1*edim2];
    Grids[mygrid]->GridData->SetPotentialField(PotentialField);

    /* allocate and initialize the MGMPI rhs and solution vectors 
       in mygrid (these will be reused throughout the simulation) */
    float *MGMPIrhsvec = new float[mdim0*mdim1*mdim2];
    PotentialBdry->SetMGMPIrhs(MGMPIrhsvec);
    /* Note: solvec currently requires additional padding for access
       within MGMPI (hence the additional mdim0*(mdim1+1)+1 in length) */
    float *MGMPIsolvec = new float[mdim0*mdim1*mdim2+mdim0*(mdim1+1)+1];
    PotentialBdry->SetMGMPIsol(MGMPIsolvec);
    

    /* Initialize MGMPI */
    MGMPI_Init();

    /* create and set up MGMPI parameters object */
    if (debug)
      fprintf(stdout,"CompPotFieldL0Iso: loading MGparams\n");
    MGMPI_ParametersCreate(&MGparams);
    MGMPI_ParametersReadFile(MGparams, "MGMPI_GravParams.in");

    /* create MGMPI concurrency object */
    if (debug)
      fprintf(stdout,"CompPotFieldL0Iso: creating MGconcurr object\n");
    MGMPI_ConcurrentCreateCartBypass(&MGconcurr, MGparams, MPI_COMM_WORLD, 
				     Eint32(layout[0]), Eint32(layout[1]), 
				     Eint32(layout[2]), Eint32(location[0]), 
				     Eint32(location[1]), Eint32(location[2]),
				     Eint32(x0Nb[0]), Eint32(x0Nb[1]), 
				     Eint32(x1Nb[0]), Eint32(x1Nb[1]), 
				     Eint32(x2Nb[0]), Eint32(x2Nb[1]), 
				     Eint32(per0), Eint32(per1), Eint32(per2));

    /* create and set up MGMPI distribute object */
    if (debug)
      fprintf(stdout, "CompPotFieldL0Iso: creating MGdistrib object\n");
    MGMPI_DistributeCreate(&MGdistrib, MGparams, MGconcurr);
    MGMPI_DistributeSetLocalSize(MGdistrib, Eint32(MGx0loc), 
				 Eint32(MGx1loc), Eint32(MGx2loc));

    /* create and set up MGMPI boundary conditions object, */
    /*   (save setting of actual BC values for each solve) */
    if (debug)
      fprintf(stdout, "CompPotFieldL0Iso: creating MGbc object\n");
    MGMPI_BcCreate(&MGbc, MGparams, MGdistrib);
    if (per0 == 0) {
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_LOWER_X, MGMPI_BC_DIRICHLET);
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_UPPER_X, MGMPI_BC_DIRICHLET);
      if (debug)
	fprintf(stdout, "  Self-Gravity BCs in x0 set to isolating.\n");
    }
    else {
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_LOWER_X, MGMPI_BC_PERIODIC);
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_UPPER_X, MGMPI_BC_PERIODIC);
      if (debug)
	fprintf(stdout, "  Self-Gravity BCs in x0 set to periodic.\n");
    }
    if (per1 == 0) {
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_LOWER_Y, MGMPI_BC_DIRICHLET);
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_UPPER_Y, MGMPI_BC_DIRICHLET);
      if (debug)
	fprintf(stdout, "  Self-Gravity BCs in x1 set to isolating.\n");
    }
    else {
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_LOWER_Y, MGMPI_BC_PERIODIC);
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_UPPER_Y, MGMPI_BC_PERIODIC);
      if (debug)
	fprintf(stdout, "  Self-Gravity BCs in x1 set to periodic.\n");
    }
    if (per2 == 0) {
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_LOWER_Z, MGMPI_BC_DIRICHLET);
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_UPPER_Z, MGMPI_BC_DIRICHLET);
      if (debug)
	fprintf(stdout, "  Self-Gravity BCs in x2 set to isolating.\n");
    }
    else {
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_LOWER_Z, MGMPI_BC_PERIODIC);
      MGMPI_BcSetFaceType(MGbc, MGMPI_FACE_UPPER_Z, MGMPI_BC_PERIODIC);
      if (debug)
	fprintf(stdout, "  Self-Gravity BCs in x2 set to periodic.\n");
    }


    /* create and set up MGMPI grid object */
    if (debug)
      fprintf(stdout, "CompPotFieldL0Iso: creating MGgrid object\n");
    MGMPI_GridCreate(&MGgrid, MGparams, MGdistrib);
    
    /* create matrix, rhs and solution vectors on the grid */
    if (debug)
      fprintf(stdout, "CompPotFieldL0Iso: creating matrix/vector objects\n");
    MGMPI_MatrixCreate(&MGa, MGparams, MGgrid, MGgrid);
    MGMPI_VectorCreate(&MGx, MGparams, MGgrid);
    MGMPI_VectorCreate(&MGb, MGparams, MGgrid);

    /* set up standard, 7pt Laplacian matrix in MGa   */
    /* (since matrix is constant-coefficient, we only */
    /*  need to set one element per stencil entry)    */
    float dx = Grids[0]->GridData->GetGravitatingMassFieldCellSize();
    float dy = dx;
    float dz = dx;
    if (debug)
      fprintf(stdout, "CompPotFieldL0Iso: setting matrix entries, dx=%g\n",dx);
    float dxi2 = 1.0/dx/dx;
    float dyi2 = 1.0/dy/dy;
    float dzi2 = 1.0/dz/dz;
    float diagval = -2.0*(dxi2 + dyi2 + dzi2); 
    Eint32 zero32=0; 
    Eint32 one32=1; 
    Eint32 mone32=-1;
    MGMPI_MatrixSetUnknown(MGa, zero32, zero32, zero32, 
			   zero32, zero32, zero32, diagval);
    MGMPI_MatrixSetUnknown(MGa, zero32, zero32, zero32, 
			   mone32, zero32, zero32, dxi2);
    MGMPI_MatrixSetUnknown(MGa, zero32, zero32, zero32,  
			   one32, zero32, zero32, dxi2);
    MGMPI_MatrixSetUnknown(MGa, zero32, zero32, zero32, 
			   zero32, mone32, zero32, dyi2);
    MGMPI_MatrixSetUnknown(MGa, zero32, zero32, zero32, 
			   zero32,  one32, zero32, dyi2);
    MGMPI_MatrixSetUnknown(MGa, zero32, zero32, zero32, 
			   zero32, zero32, mone32, dzi2);
    MGMPI_MatrixSetUnknown(MGa, zero32, zero32, zero32, 
			   zero32, zero32,  one32, dzi2);

    /* create and set up MGMPI stopping criteria */
    if (debug)
      fprintf(stdout, "CompPotFieldL0Iso: creating MGstopping object\n");
    MGMPI_StoppingCreate(&MGstopping, MGparams);

    /* create and set up MGMPI solver object */
    if (debug)
      fprintf(stdout, "CompPotFieldL0Iso: creating MGsolver object\n");
    MGMPI_SolverCreate(&MGsolver, MGparams, MGstopping);

    /* unset FirstCall flag */
    FirstCall = FALSE;

  } // end: if (FirstCall)


  /* set up rhs vector */
  /*   Compute cosmology expansion term at time = t+1/2dt */
  FLOAT a = 1, dadt, MidTime = Grids[mygrid]->GridData->ReturnTime() + 
                               0.5*Grids[mygrid]->GridData->ReturnTimeStep();
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(MidTime, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
  } // end: if (ComovingCoordinates)


  /* insert boundary conditions into matrix (done here to allow for
     changing BCs in time) */
//   if (debug)
//     fprintf(stdout, "CompPotFieldL0Iso: inserting local BC values\n");
  float *BCptr;
  /*    grid owns part of x0L boundary */
  if ((per0 == 0) && x0LBdry) {
    BCptr = PotentialBdry->GetBoundaryValues(0,0);
    if (BCptr != NULL)
      MGMPI_BcSetFaceValues(MGbc, MGMPI_FACE_LOWER_X, BCptr, 
			    Eint32(x1loc), Eint32(x2loc));
    else {
      fprintf(stderr, "Error: grid %"ISYM", x0L boundary pointer NULL", mygrid);
      return FAIL;
    }
  }
  /*    grid owns part of x0R boundary */
  if ((per0 == 0) && x0RBdry) {
    BCptr = PotentialBdry->GetBoundaryValues(0,1);
    if (BCptr != NULL)
      MGMPI_BcSetFaceValues(MGbc, MGMPI_FACE_UPPER_X, BCptr, 
			    Eint32(x1loc), Eint32(x2loc));
    else {
      fprintf(stderr, "Error: grid %"ISYM", x0R boundary pointer NULL", mygrid);
      return FAIL;
    }
  }
  /*    grid owns part of x1L boundary */
  if ((per1 == 0) && x1LBdry) {
    BCptr = PotentialBdry->GetBoundaryValues(1,0);
    if (BCptr != NULL)
      MGMPI_BcSetFaceValues(MGbc, MGMPI_FACE_LOWER_Y, BCptr, 
			    Eint32(x2loc), Eint32(x0loc));
    else {
      fprintf(stderr, "Error: grid %"ISYM", x1L boundary pointer NULL", mygrid);
      return FAIL;
    }
  }
  /*    grid owns part of x1R boundary */
  if ((per1 == 0) && x1RBdry) {
    BCptr = PotentialBdry->GetBoundaryValues(1,1);
    if (BCptr != NULL)
      MGMPI_BcSetFaceValues(MGbc, MGMPI_FACE_UPPER_Y, BCptr, 
			    Eint32(x2loc), Eint32(x0loc));
    else {
      fprintf(stderr, "Error: grid %"ISYM", x1R boundary pointer NULL", mygrid);
      return FAIL;
    }
  }
  /*    grid owns part of x2L boundary */
  if ((per2 == 0) && x2LBdry) {
    BCptr = PotentialBdry->GetBoundaryValues(2,0);
    if (BCptr != NULL)
      MGMPI_BcSetFaceValues(MGbc, MGMPI_FACE_LOWER_Z, BCptr, 
			    Eint32(x0loc), Eint32(x1loc));
    else {
      fprintf(stderr, "Error: grid %"ISYM", x2L boundary pointer NULL", mygrid);
      return FAIL;
    }
  }
  /*    grid owns part of x2R boundary */
  if ((per2 == 0) && x2RBdry) {
    BCptr = PotentialBdry->GetBoundaryValues(2,1);
    if (BCptr != NULL)
      MGMPI_BcSetFaceValues(MGbc, MGMPI_FACE_UPPER_Z, BCptr, 
			    Eint32(x0loc), Eint32(x1loc));
    else {
      fprintf(stderr, "Error: grid %"ISYM", x2R boundary pointer NULL", mygrid);
      return FAIL;
    }
  }
  MGMPI_MatrixSetBc(MGa, MGbc, MGbc);


  /*    Compute rhs coefficient = 4*pi*G/a  */
  float coef = GravitationalConstant/a;


  /*    Scale density by coef, and interpolate to FD grid for system rhs  */
  /*       x*Li, x*Ri denote if the L/R face is internal to the domain    */
  /*       es*, ef* denote where the enzo active extents start and finish */
  /*       off* denotes the offset between the enzo and mgmpi grids       */
//   if (debug)
//     fprintf(stdout, "CompPotFieldL0Iso: setting up system rhs\n");
  float *grav_mass_field = Grids[mygrid]->GridData->GetGravitatingMassField();
  int x0Li, x0Ri, x1Li, x1Ri, x2Li, x2Ri;
  int es0, es1, es2, ef0, ef1, ef2, off0, off1, off2;
  x0Li = 0;  x0Ri = 0;  x1Li = 0;  x1Ri = 0;  x2Li = 0;  x2Ri = 0;  
  float *rhsvec = PotentialBdry->GetMGMPIrhs();
  int GravityGhosts = Grids[mygrid]->GridData->GetGravityGhostZones();
  es0 = es1 = es2 = 1 + GravityGhosts;
  off0 = off1 = off2 = -GravityGhosts;
  ef0 = edim0 - GravityGhosts;  
  ef1 = edim1 - GravityGhosts;  
  ef2 = edim2 - GravityGhosts;
  if (x0LBdry)
    if (per0 == 1)  off0 = 1-GravityGhosts;
    else  {es0 = 1;  off0 = 1;}
  else  x0Li = 1;

  if (x1LBdry)
    if (per1 == 1)  off1 = 1-GravityGhosts;
    else  {es1 = 1;  off1 = 1;}
  else  x1Li = 1;

  if (x2LBdry)
    if (per2 == 1)  off2 = 1-GravityGhosts;
    else  {es2 = 1;  off2 = 1;}
  else  x2Li = 1;
  
  if (x0RBdry)
    { if (per0 == 0)  ef0 = edim0; }
  else  x0Ri = 1;

  if (x1RBdry) 
    { if (per1 == 0)  ef1 = edim1; }
  else  x1Ri = 1;

  if (x2RBdry) 
    { if (per2 == 0)  ef2 = edim2; }
  else  x2Ri = 1;

  /* get background density rho0 */
  float rho0 = PotentialBdry->GetOmegaBaryonNow();
//   if (debug)
//     fprintf(stdout, "CompPotFieldL0Iso: OmegaBaryonNow = %g\n",rho0);

  /* subtract rho0 from grav_mass_field, and scale by gravity constant */
  int idx;
  for (idx=0; idx<(edim0*edim1*edim2); idx++)
    grav_mass_field[idx] = (grav_mass_field[idx] - rho0)*coef;

  /* in case of non-conforming hydro/gravity BCs, zero out 
     exterior enzo ghost cells for isolating solve (in case 
     of conforming iso/iso BCs, these will already be zero) */
  /*   x0, left external boundary */
  if (x0LBdry  &&  (per0 == 0)) {
    for (k=0; k<edim2; k++)
      for (j=0; j<edim1; j++)
	for (i=0; i<GravityGhosts; i++)
	  grav_mass_field[(k*edim1 + j)*edim0 + i] = 0.0;
  }
  /*   x0, right external boundary */
  if (x0RBdry  &&  (per0 == 0)) {
    for (k=0; k<edim2; k++)
      for (j=0; j<edim1; j++)
	for (i=edim0-GravityGhosts; i<edim0; i++)
	  grav_mass_field[(k*edim1 + j)*edim0 + i] = 0.0;
  }
  /*   x1, left external boundary */
  if (x1LBdry  &&  (per1 == 0)) {
    for (k=0; k<edim2; k++)
      for (j=0; j<GravityGhosts; j++)
	for (i=0; i<edim0; i++)
	  grav_mass_field[(k*edim1 + j)*edim0 + i] = 0.0;
  }
  /*   x1, right external boundary */
  if (x1RBdry  &&  (per1 == 0)) {
    for (k=0; k<edim2; k++)
      for (j=edim1-GravityGhosts; j<edim1; j++)
	for (i=0; i<edim0; i++)
	  grav_mass_field[(k*edim1 + j)*edim0 + i] = 0.0;
  }
  /*   x2, left external boundary */
  if (x2LBdry  &&  (per2 == 0)) {
    for (k=0; k<GravityGhosts; k++)
      for (j=0; j<edim1; j++)
	for (i=0; i<edim0; i++)
	  grav_mass_field[(k*edim1 + j)*edim0 + i] = 0.0;
  }
  /*   x2, right external boundary */
  if (x2RBdry && (per2 == 0)) {
    for (k=edim2-GravityGhosts; k<edim2; k++)
      for (j=0; j<edim1; j++)
	for (i=0; i<edim0; i++)
	  grav_mass_field[(k*edim1 + j)*edim0 + i] = 0.0;
  }


  /* extend lower proc internal extents by one for interpolation */
   if (!x0RBdry)  ef0 += 1;
   if (!x1RBdry)  ef1 += 1;
   if (!x2RBdry)  ef2 += 1;

   /* interpolate grav_mass_field (Enzo grid) to rhsvec (MGMPI grid) */
//   if (debug)
//     fprintf(stdout, "CompPotFieldL0Iso: interpolating rhs to MGMPI grid\n");
  FORTRAN_NAME(enzo_to_mgmpi)(grav_mass_field, rhsvec, &edim0, &edim1, 
			      &edim2, &mdim0, &mdim1, &mdim2, &es0, &es1, 
			      &es2, &ef0, &ef1, &ef2, &off0, &off1, &off2);

#ifdef USE_MPI
  /* communicate boundary data in MGMPI grid */
  if (MGMPIBdryExchange(rhsvec, layout, location, x0Nb[0], 
			x0Nb[1], x1Nb[0], x1Nb[1], x2Nb[0], x2Nb[1], 
			mdim0, mdim1, mdim2) != SUCCESS) {
    fprintf(stderr,"Error in MGMPIBdryExchange, Grid %"ISYM"\n",mygrid);
    return FAIL;
  }
#endif

  /* reset internal extents */
   if (!x0RBdry)  ef0 -= 1;
   if (!x1RBdry)  ef1 -= 1;
   if (!x2RBdry)  ef2 -= 1;


  /* alias rhs vector to MGMPI vector MGb (send in first 'active' index) */
//   if (debug)
//     fprintf(stdout, "CompPotFieldL0Iso: aliasing rhsvec to MGb\n");
  MGMPI_VectorAliasVertices(MGb, &(rhsvec[MGStartIndex]), Eint32(mdim0), 
			    Eint32(mdim1), Eint32(mdim2));

  /* create solution vector and alias to MGMPI vector MGx */
//   if (debug)
//     fprintf(stdout, "CompPotFieldL0Iso: aliasing solvec to MGx\n");
  float *solvec = PotentialBdry->GetMGMPIsol();
  for (i=0; i<mdim0*mdim1*mdim2; i++)  solvec[i] = 0.0;  // initialize to 0.0
  MGMPI_VectorAliasVertices(MGx, &(solvec[MGStartIndex]), Eint32(mdim0), 
			    Eint32(mdim1), Eint32(mdim2));


  /* perform gravity solve for potential */
//   if (debug)
//     fprintf(stdout, "CompPotFieldL0Iso: applying MGMPI solver\n");
  MGMPI_SolverApply(MGsolver, MGa, MGx, MGb);



//   // output solution
//   FILE *fptr;
//   if (debug)
//     fprintf(stdout,"CompPotFieldL0Iso: outputting solution\n");
//   switch (mygrid) {
//   case 0:
//     printf("  Grid 0, outputting solvec to solinit_0\n");
//     fptr = fopen("solinit_0", "w");
//     for (k=0, idx=0; k<mdim2; k++)
//       for (j=0; j<mdim1; j++)
//         for (i=0; i<mdim0; i++, idx++)
//           fprintf(fptr, "  %"ISYM"  %"ISYM"  %"ISYM"  %g\n",
// 		  i+1,j+1,k+1,solvec[idx]);
//     fclose(fptr);
//     break;
    
//   case 1:
//     printf("  Grid 1, outputting solvec to solinit_1\n");
//     fptr = fopen("solinit_1", "w");
//     for (k=0, idx=0; k<mdim2; k++)
//       for (j=0; j<mdim1; j++)
//         for (i=0; i<mdim0; i++, idx++)
//           fprintf(fptr, "  %"ISYM"  %"ISYM"  %"ISYM"  %g\n",
// 		  i+1,j+1,k+1,solvec[idx]);
//     fclose(fptr);
//     break;
    
//   case 2:
//     printf("  Grid 2, outputting solvec to solinit_2\n");
//     fptr = fopen("solinit_2", "w");
//     for (k=0, idx=0; k<mdim2; k++)
//       for (j=0; j<mdim1; j++)
//         for (i=0; i<mdim0; i++, idx++)
//           fprintf(fptr, "  %"ISYM"  %"ISYM"  %"ISYM"  %g\n",
// 		  i+1,j+1,k+1,solvec[idx]);
//     fclose(fptr);
//     break;
    
//   case 3:
//     printf("  Grid 3, outputting solvec to solinit_3\n");
//     fptr = fopen("solinit_3", "w");
//     for (k=0, idx=0; k<mdim2; k++)
//       for (j=0; j<mdim1; j++)
//         for (i=0; i<mdim0; i++, idx++)
//           fprintf(fptr, "  %"ISYM"  %"ISYM"  %"ISYM"  %g\n",
// 		  i+1,j+1,k+1,solvec[idx]);
//     fclose(fptr);
//     break;
    
//   default:
//     printf("  Grid %"ISYM", not outputting solvec\n",mygrid);
//   }


  /* interpolate MGMPI solution vector to ENZO gravity potential */
//   if (debug)
//     fprintf(stdout, "CompPotFieldL0Iso: interpolating solvec to Enzo grid\n");
  float *potential_field = Grids[mygrid]->GridData->GetPotentialField();
  FORTRAN_NAME(mgmpi_to_enzo)(solvec, potential_field, &mdim0, &mdim1, 
			      &mdim2, &x0Li, &x0Ri, &x1Li, &x1Ri, &x2Li, 
			      &x2Ri, &edim0, &edim1, &edim2, &es0, &es1, 
			      &es2, &ef0, &ef1, &ef2, &off0, &off1, &off2);

#ifdef USE_MPI
  /* update gravitational potential on neighbors */
  if (GravityBdryExchange(potential_field, layout, location, GravityGhosts,
			  edim0, edim1, edim2, x0Nb[0], x0Nb[1], x1Nb[0], 
			  x1Nb[1], x2Nb[0], x2Nb[1]) != SUCCESS) {
    fprintf(stderr,"Error in GravityBdryExchange, Grid %"ISYM"\n",mygrid);
    return FAIL;
  }
#endif


//   /* Clean up. */ 
//   if (debug)
//     fprintf(stdout, "CompPotFieldL0Iso: FINISHED!\n");

  return SUCCESS;
}
#endif




/******************************************************************/
/*  ComputePotentialFieldLevelZeroPer performs a root-grid        */
/*  potential field solver using periodic boundary conditions,    */
/*  via an FFT-based solution strategy.  This solver just calls   */
/*  the pre-existing code that Greg Bryan wrote.                  */
/******************************************************************/
int ComputePotentialFieldLevelZeroPer(TopGridData *MetaData,
				      HierarchyEntry *Grids[], int NumberOfGrids)
{
  /* Static declarations (for Green's function). */
 
  static int FirstCall = TRUE, NumberOfGreensRegions;
  static region *GreensRegion;
 
  /* Declarations. */
 
  region *OutRegion = NULL;
  int NumberOfOutRegions, DomainDim[MAX_DIMENSION];
  int i, j, n, grid, grid2;
 
  /* Allocate space for grid info. */
 
  int NumberOfRegions = NumberOfGrids;
  region *InitialRegion = new region[NumberOfRegions];
 
  /* Compute adot/a at time = t+1/2dt (time-centered). */
 
  FLOAT a = 1, dadt, MidTime = Grids[0]->GridData->ReturnTime() +
                           0.5*Grids[0]->GridData->ReturnTimeStep();
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(MidTime, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
 
  /* ------------------------------------------------------------------- */
  /* If this is the first time this routine has been called, then generate
     the Green's function. */
 
  if (FirstCall) {
 
    if (MetaData->GravityBoundary == TopGridPeriodic) {
 
      /* Periodic -- Prepare in k-space. */
 
      NumberOfGreensRegions = NumberOfGrids;
      GreensRegion = new region[NumberOfGreensRegions];
      for (grid = 0; grid < NumberOfGrids; grid++)
	if (Grids[grid]->GridData->PreparePeriodicGreensFunction(
					     &(GreensRegion[grid])) == FAIL) {
	  fprintf(stderr, "Error in grid->PreparePeriodicGreensFunction.\n");
	  return FAIL;
	}
 
    } else {
 
      fprintf(stderr, "Isolated BC's not yet implemented.\n");
      return FAIL;
 
#ifdef UNUSED
 
      for (grid = 0; grid < NumberOfGrids; grid++) {
 
	/* Generate Green's function in real space (doesn't work!). */
	
	if (Grids[grid]->GridData->PrepareGreensFunction() == FAIL) {
	  fprintf(stderr, "Error in grid->PrepareGreensFunction.\n");
	  return FAIL;
	}
 
	/* Turn it into regions. */
	
	if (Grids[grid]->GridData->PrepareFFT(&InitialRegion[grid],
					      POTENTIAL_FIELD, DomainDim)
	    == FAIL) {
	  fprintf(stderr, "Error in grid->PrepareFFT.\n");
	  return FAIL;
	}
 
      } // end loop over grids
 
      /* Forward FFT Green's function. */
 
      if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
				   &GreensRegion, &NumberOfGreensRegions,
				   DomainDim, MetaData->TopGridRank,
				   FFT_FORWARD, FALSE) == FAIL) {
	fprintf(stderr, "Error in CommunicationParallelFFT.\n");
	return FAIL;
      }
 
#endif /* UNUSED */
 
    } // end: if (Periodic)
 
    FirstCall = FALSE;
 
  } // end: if (FirstCall)
 
  /* ------------------------------------------------------------------- */
  /* Generate FFT regions for density field. */
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->PrepareFFT(&InitialRegion[grid],
					  GRAVITATING_MASS_FIELD, DomainDim)
	== FAIL) {
      fprintf(stderr, "Error in grid->PrepareFFT.\n");
      return FAIL;
    }
 
  /* Forward FFT density field. */
 
  if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
			       &OutRegion, &NumberOfOutRegions,
			       DomainDim, MetaData->TopGridRank,
			       FFT_FORWARD, TRUE) == FAIL) {
    fprintf(stderr, "Error in CommunicationParallelFFT.\n");
    return FAIL;
  }
 
  /* Quick error check. */
 
  if (NumberOfOutRegions != NumberOfGreensRegions) {
    fprintf(stderr, "OutRegion(%"ISYM") != GreensRegion(%"ISYM")\n", NumberOfOutRegions,
	    NumberOfGreensRegions);
    return FAIL;
  }
 
  /* Compute coefficient for Greens function. */
 
  float coef = GravitationalConstant/a;
  //  for (int dim = 0; dim < MetaData->TopGridRank; dim++)
  //    coef *= (DomainRightEdge[dim] - DomainLeftEdge[dim])/float(DomainDim[dim]);
			
  /* Multiply density by Green's function to get potential. */
 
  for (i = 0; i < NumberOfGreensRegions; i++)
    if (OutRegion[i].Data != NULL) {
      int size = OutRegion[i].RegionDim[0]*OutRegion[i].RegionDim[1]*
	         OutRegion[i].RegionDim[2];
      for (n = 0, j = 0; j < size; j += 2, n++) {
	OutRegion[i].Data[j  ] *= coef*GreensRegion[i].Data[n];
	OutRegion[i].Data[j+1] *= coef*GreensRegion[i].Data[n];
      }
    }
 
  /* Inverse FFT potential field. */
 
  if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
			       &OutRegion, &NumberOfOutRegions,
			       DomainDim, MetaData->TopGridRank,
			       FFT_INVERSE, TRUE) == FAIL) {
    fprintf(stderr, "Error in CommunicationParallelFFT.\n");
    return FAIL;
  }
 
  /* Copy Potential in active region into whole grid. */
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->FinishFFT(&InitialRegion[grid], POTENTIAL_FIELD,
			       DomainDim) == FAIL) {
      fprintf(stderr, "Error in grid->FinishFFT.\n");
      return FAIL;
    }
 
  /* Update boundary regions of potential
     (first set BCTempL/R which are fluid BC's because that's the format
      that CheckForOverlap takes). */
 
  boundary_type BCTempLeft[MAX_DIMENSION], BCTempRight[MAX_DIMENSION];
  if (Grids[0]->GridData->ReturnGravityBoundaryType() == TopGridPeriodic) {
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      BCTempLeft[dim] = BCTempRight[dim] = periodic;
  } else {
    fprintf(stderr, "Only periodic gravity BC's allowed.\n");
    return FAIL;
  }
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				      BCTempLeft, BCTempRight,
     	                              &grid::CopyPotentialField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	return FAIL;
      }
 
  /* Clean up. */
 
  delete [] InitialRegion;
  if (OutRegion != InitialRegion)
    delete [] OutRegion;
 
  if (CopyGravPotential)
    for (grid = 0; grid < NumberOfGrids; grid++)
    {
      fprintf(stderr, "Call CP from ComputePotentialFieldLevelZero\n");
      Grids[grid]->GridData->CopyPotentialToBaryonField();
    }
 
  return SUCCESS;
}

/******************************************************************/
