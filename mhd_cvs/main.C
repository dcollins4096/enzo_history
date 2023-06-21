
/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  AMR MAIN CODE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/    This is main() for the amr code.  It interprets the arguments and
/      then calls the appropriate routines depending on whether we are
/      doing a new run, a restart, an extraction, or a projection.
/
************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <exception>
#include "performance.h"
#include <string.h>
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */
#include "macros_and_parameters.h"
#include "typedefs.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "flowdefs.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"


#ifdef USE_FLEXIO
#include "WriteEvolution.h"
#endif

#include "pout.h"
#undef DEFINE_STORAGE
#ifdef USE_PABLO
extern "C"
{
#include "ProcTrace.h"
}
#endif

/* function prototypes */

int InitializeNew(  char *filename, HierarchyEntry &TopGrid, TopGridData &tgd, 
		    ExternalBoundary &Exterior, float *Initialdt);
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
int EvolveHierarchy(                HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, LevelHierarchyEntry *Array[],
		    float Initialdt);
void ExtractSection(HierarchyEntry &TopGrid, TopGridData &tgd, 
		    LevelHierarchyEntry *Array[], ExternalBoundary *Exterior,
		    int ExtractStart[], int ExtractEnd[], 
		    FLOAT ProjectStartCoordinates[],
		    FLOAT ProjectEndCoordinates[], int ExtractLevel);
int OutputLevelInformation(FILE *fptr, TopGridData &tgd, 
			   LevelHierarchyEntry *Array[]);
int ProjectToPlane(TopGridData &MetaData, LevelHierarchyEntry *LevelArray[],
		   int ProjectStart[], int ProjectEnd[], 
		   FLOAT ProjectStartCoordinates[], 
		   FLOAT ProjectEndCoordinates[], int ProjectLevel,
		   int ProjectionDimension, char *ProjectionFileName,
		   int ProjectionSmooth, ExternalBoundary *Exterior);
int OutputAsParticleData(TopGridData &MetaData, 
			 LevelHierarchyEntry *LevelArray[],
			 int RegionStart[], int RegionEnd[], 
			 FLOAT RegionStartCoordinates[], 
			 FLOAT RegionEndCoordinates[], int RegionLevel,
			 char *OutputFileName);
int InterpretCommandLine(int argc, char *argv[], char *myname,
			 int &restart, int &debug, int &extract,
			 int &ForceAnalysisOnly, char *AnalysisParameterFileName[],
			 int &InformationOutput,
			 int &OutputAsParticleDataFlag,
			 int &project, int &ProjectionDimension, 
			 int &ProjectionSmooth,
			 char *ParameterFile[],
			 int RegionStart[], int RegionEnd[], 
			 FLOAT RegionStartCoordinates[],
			 FLOAT RegionEndCoordinates[], 
			 int &Level, int MyProcessorNumber);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
int CommunicationInitialize(int *argc, char **argv[]);
int CommunicationFinalize();
int CommunicationPartitionGrid(HierarchyEntry *Grid);
void DoAllAnalysis(HierarchyEntry *TopGrid, TopGridData &MetaData);		 
void my_exit(int status);

#include "message.h"
#include "error.h"

int main(int argc, char *argv[])
{

  // SetupNormal
  // MHD_Roe1
  // MHD_HLLE
  // MHD_Fluxes
  // ComputeAccelerations
  // SetAccelerationBoundary
  // MHD_Athena
  // FlagCellsToBeRefinedByMHD
  // UpdateFromFinerGrids
  // interp3d
  // RebuildHierarchy
  // InterpolateFieldValues
  // InterpolateBoundaryFromParent  
  // mhd_interpolate
  // MHDShockInitializeGrid
  // MHDOrszagTangInitGrid
  // MHDBlastInitialize
  // MHDBlastInitializeGrid
  // SolveMHDEquations()
  // EvolveLevel
  // EvolveHierarchy
  // ComputeTimeStep
  // InitializeNew Prepare
  // mhd_harten
  // mhdtvd
  // CopyZonesFromGrid
  // MHD_Diagnose
  // PrepareGridDerivedQuantities
  // PrepareGrid
  // CorrectForRefinedFluxes  
  /* To Do List
  fprintf(stderr,"there's a bunch of ifs in mhd_interpolate.  TestX.  Get rid of them.\n");
  fprintf(stderr,"test for CVS\n");
  fprintf(stderr,"----- Notes to dcc: -----\n");
  fprintf(stderr,"dcc: irradicate the OldCenterdB crap.\n");
  fprintf(stderr,"dcc: fix where the E project happens: you should just clean it from the baryon\n");
  fprintf(stderr,"dcc: clean up non-processor fields everywhere.");
  fprintf(stderr,"dcc: pout->Pig for Send Old Fields, Interpolate crap.\n");
  fprintf(stderr,"dcc: pout->Prong for Project Face crap\n");
  fprintf(stderr,"dcc: clean up calls to mhd_interpolate.\n");
  fprintf(stderr,"DCC: there's a redundant call to SetBoundaryConditions in EvolveLevel?\n");
  fprintf(stderr,"DCC: fix double counting for flux correction.\n");
  fprintf(stderr,"DCC: turned off the WARNING! difusion crap in ReadParameter\n");
  //fprintf(stderr,"kludge: Maximum Subgrid Size to 1000 (Proto Subgrid_Acceptable)\n");
  fprintf(stderr,"dcc: Kludge in ProtoSubgrid_AcceptableSubgrid: No Processor number Check\n");
  fprintf(stderr,"kludge: ExternalBoundary_SetMagneticBoundary exits\n");
  fprintf(stderr,"kludge: mucked up InitStyle = 8, Magnetic Field\n");
  fprintf(stderr,"kludge: not exiting when DivB != 0\n");
  fprintf(stderr,"kludge: movie dump for martin\n");
  //fprintf(stderr,"dcc: No WriteInThis output in InterpBoundary\n");
  fprintf(stderr,"dcc: Strang?\n");
  fprintf(stderr,"dcc: Remove MHDcUnits\n");
  fprintf(stderr,"dcc: Be careful about where Bc is centered.\n");
  fprintf(stderr,"kludge: dccCounter12 = -12 initially.  DON'T USE\n");
  fprintf(stderr,"kludge: mhd_interpolate, cycle .eq. -12 crap: hard wired for r=2\n");
  fprintf(stderr,"kludge: no divergence warnings in mhd_interpolate. \n");
  fprintf(stderr,"kludge: THERE IS DIVERGENCE IN ALL BOUNDARIES.  \n");
  fprintf(stderr,"kludge: the fix for non-periodic boundaries in CZFG removed.\n");
  */

  try {

  /* Initialize Communications package. */

    CommunicationInitialize(&argc, &argv);
    wall_time_start();
    wall_time("Start Main");    
#ifdef OOC_BOUNDARY
    ExternalBoundaryIO = TRUE;
    ExternalBoundaryTypeIO = FALSE;
    ExternalBoundaryValueIO = FALSE;
#else
    //ExternalBoundaryIO = FALSE;
    //ExternalBoundaryTypeIO = FALSE;
    //ExternalBoundaryValueIO = FALSE;
#endif


    /* Print out splash screen, just in case people forgot who wrote the thing... */

    if (MyProcessorNumber == 0) {
      printf("\n\n\n");
      printf("*********************************************************************\n");
      printf("*                                                                   *\n");
      printf("*  Enzo is an Eulerian adaptive mesh refinement cosmology code      *\n");
      printf("*  developed and maintained by the Laboratory for Computational     *\n");
      printf("*  Astrophysics at the University of California at San Diego under  *\n");
      printf("*  the direction of Michael L. Norman.                              *\n");
      printf("*                                                                   *\n");
      printf("*  Enzo was originally written by Greg Bryan at the National Center *\n");
      printf("*  for Supercomputing Applications at the University of Illinois in *\n");
      printf("*  Urbana-Champaign.  For more inforamtion, see the Enzo website at *\n");
      printf("*  http://cosmos.ucsd.edu/enzo/.                                    *\n");
      printf("*                                                                   *\n");
      printf("*********************************************************************\n");
      printf("\n\n\n");
      fflush(stdout);
    }

#ifdef USE_PABLO
   HDFinitTrace( "myTrace",ID_ALLHDF,MPI_RUNTIME_TRACE);
#endif
#ifdef USE_JBPERF
    jb::perf.init ();
    jb::perf.user ("mpi-send");
    jb::perf.user ("mpi-recv");
    jb::perf.user ("mpi-gather");
    jb::perf.user ("mpi-reduce");
    jb::perf.user ("mpi-barrier");
    jb::perf.user ("hdf5-read");
    jb::perf.user ("hdf5-write");
#endif

#ifdef MPI_INSTRUMENTATION
  Start_Wall_Time = MPI_Wtime(); // Zhiling Lan's instrumented
#endif

  /* Main declarations */

  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

  /* Initialize */

  int restart                  = FALSE,
      extract                  = FALSE,
      OutputAsParticleDataFlag = FALSE,
      InformationOutput        = FALSE,
      project                  = FALSE, 
      ProjectionDimension      = INT_UNDEFINED,
      ProjectionSmooth         = FALSE;
  debug                        = FALSE;
  trace                        = FALSE;
  ForceAnalysisOnly            = FALSE;
  AnalysisParameterFileName    = NULL;
  char *ParameterFile          = NULL;
  char *myname                 = argv[0];
  int RegionStart[MAX_DIMENSION], RegionEnd[MAX_DIMENSION], RegionLevel;
  FLOAT RegionStartCoordinates[MAX_DIMENSION],
        RegionEndCoordinates[MAX_DIMENSION];
  float Initialdt              = 0;

#ifdef TRACE
  char pid[5];
  sprintf(pid, "%4.4d", MyProcessorNumber);
  char *tracelog = new char[strlen(pid)+6+1];
  strcpy(tracelog, "Trace.");
  strcat(tracelog, pid);
  trace_fptr = fopen( tracelog, "w" );
  trace_level=0;
  trace = TRUE;
#endif

  ZLAN_INIT;

  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;

  /* Interpret command-line arguments. */

  // Metadata.commandline: options for the WriteEvolution module.
  // MetaData.commandline.parseAndEat(argc,argv); // parse cmdln options
  // and remove them so as not to confuse InterpretCommandLine
  //  MetaData.commandline.print(stderr); // print results of cmdln parse
  
  if (InterpretCommandLine(argc, argv, myname, restart, debug, extract,
			   ForceAnalysisOnly, &AnalysisParameterFileName,
			   InformationOutput, OutputAsParticleDataFlag,
			   project, ProjectionDimension, ProjectionSmooth,
			   &ParameterFile,
			   RegionStart, RegionEnd, 
			   RegionStartCoordinates, RegionEndCoordinates,
			   RegionLevel, MyProcessorNumber) == FAIL)
    my_exit(EXIT_FAILURE);

  // =========================================================
  // Initialize jbPerf
  // =========================================================

  if( ForceAnalysisOnly){
    fprintf(stderr,"Forcing parameter file %s\n",AnalysisParameterFileName);
  }
#ifdef USE_JBPERF
    jb::perf.begin ();
    {
      char ipstr[10] = "0";
#   ifdef USE_MPI
      int ip,ierr;
      CHECK_MPI_ERROR(MPI_Comm_rank (MPI_COMM_WORLD,&ip));
      sprintf (ipstr,"%d",ip);
#   endif
      jb::perf.category("ip",ipstr);
    }
    jb::perf.start ("main");
#endif
  // =========================================================

  /* If we need to read the parameter file as a restart file, do it now. */

  if (project || OutputAsParticleDataFlag || extract || 
      InformationOutput || restart || ForceAnalysisOnly) {
    SetDefaultGlobalValues(MetaData);

    if (debug) printf("reading parameter file %s\n", ParameterFile);

    if (ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Error in ParameterFile %s.\n", ParameterFile);
      my_exit(EXIT_FAILURE);
    }

    if (!ParallelRootGridIO && restart && TopGrid.NextGridThisLevel == NULL)
      CommunicationPartitionGrid(&TopGrid); // partition top grid if necessary

    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Successfully read ParameterFile %s.\n", ParameterFile);

    AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels
  }

  /* Information dump. */

  if (InformationOutput) {
    OutputLevelInformation(stdout, MetaData, LevelArray);
    my_exit(EXIT_SUCCESS);
  }

  /* Other analysis.  Might be redundant with extract or whatnot. Whatever.*/
  if( ForceAnalysisOnly){
    DoAllAnalysis(&TopGrid, MetaData);
    my_exit(EXIT_SUCCESS);
  }
  /* Extract a grid subsection (doesn't return). */

  if (extract)
    ExtractSection(TopGrid, MetaData, LevelArray, &Exterior,
		   RegionStart, RegionEnd, 
		   RegionStartCoordinates, RegionEndCoordinates, RegionLevel);

  /* Output fields as particle data. */

  if (OutputAsParticleDataFlag)
    if (OutputAsParticleData(MetaData, LevelArray, RegionStart, RegionEnd, 
			     RegionStartCoordinates, RegionEndCoordinates,
			     RegionLevel, "amr.particles") == FAIL)
      my_exit(EXIT_FAILURE);
    else
      my_exit(EXIT_SUCCESS);

  /* Project 3D field to 2D plane. */

  if (project) 
    if (ProjectToPlane(MetaData, LevelArray, RegionStart, RegionEnd, 
		     RegionStartCoordinates, RegionEndCoordinates,
		     RegionLevel, ProjectionDimension, "amr.project",
		     ProjectionSmooth, &Exterior) == FAIL)
      my_exit(EXIT_FAILURE);
    else
      my_exit(EXIT_SUCCESS);

  /* Open and read parameter file (if not restarting). */
  //wall_time("Start InitializeNew");
  //wall_time("Start InitializeNew");
  if (!restart)
    if (InitializeNew(ParameterFile, TopGrid, MetaData, Exterior, &Initialdt) 
	== FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Error in Parameter File %s.\n", ParameterFile);
      my_exit(EXIT_FAILURE);
    }
    else {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Successfully read in parameter file %s.\n", 
		ParameterFile);

      AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels
    }
  //wall_time("End InitializeNew");
  //wall_time("End InitializeNew");
  /* Call the main evolution routine. */



  // OK, we are going to evolve, so we must initialize the 
  // WriteEvolution module
  //  WriteEvolution wEvolve(MetaData); // init on construction
  // Can we now pass this into EvolveHierarchy safely?

  if (EvolveHierarchy(TopGrid, MetaData, &Exterior, LevelArray, Initialdt) 
      == FAIL) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Error in EvolveHierarchy.\n");
    my_exit(EXIT_FAILURE);
  }
  else
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Successful run, exiting.\n");

  /* done */
  wall_time("End Main");
  wall_time_stop();
#ifdef MPI_INSTRUMENTATION
  /* Zhiling Lan's instrumented part */
  End_Wall_Time = MPI_Wtime();
  WallTime = End_Wall_Time - Start_Wall_Time;

  /* Zhiling Lan's instrumented part */

  char *message[] = {
    /* 0 */ "total interval time used between RebuildHierarchy"
  };

  fprintf(filePtr, "elapsed wall time is %ef .\n", WallTime);
  fprintf(filePtr, "communication time is %ef .\n", CommunicationTime);
  fprintf(filePtr, "comm. time for transferring region is %ef (%d times).\n", timer[5], counter[5]);
  fprintf(filePtr, "comm. time for sending particles is %ef (%d times).\n", timer[7], counter[7]);
  fprintf(filePtr, "comm. time for transferring particles is %ef (%d times).\n", timer[9], counter[9]);
  fprintf(filePtr, "comm. time for UpdateStarParticleCount is %ef (%d times).\n", timer[11], counter[11]);
  fprintf(filePtr, "comm. time for transferring Fluxes is %ef (%d times).\n", timer[12], counter[12]);
  fprintf(filePtr, "comm. time for ShareGrids is %ef (%d times).\n", timer[13], counter[13]);
  fprintf(filePtr, "comm. time for Transpose is %ef (%d times).\n", timer[14], counter[14]);
  fprintf(filePtr, "comm. time for BroadcastValue is %ef (%d times).\n", timer[15], counter[15]);
  fprintf(filePtr, "comm. time for MinValue is %ef (%d times).\n", timer[16], counter[16]);
  fprintf(filePtr, "total time used for RebuildHierarchy is %ef (%d times).\n", timer[1], counter[1]);
  fprintf(filePtr, message[0]);
  fprintf(filePtr, " is %ef (%d times).\n", timer[0], counter[0]);
  fprintf(filePtr, "total time used for loadbalancing is %ef (%d times).\n", timer[2], counter[2]);
  fprintf(filePtr, "total number of successful loadbalancing is %d/%d (LOAD_BALANCE_RATIO=%f).\n",counter[3], counter[2], timer[3]);
  fprintf(filePtr, "total number of flagging cells is %ef (%d times).\n",timer[4], counter[4]);
  fprintf(filePtr, "total size of region transferred is %ef (%d times). \n", timer[6], counter[5]);
  fprintf(filePtr, "total size of particles sended is %ef (%d times). \n", timer[8],counter[7]);
  fprintf(filePtr, "total size of particles transferred is %ef (%d times). \n", timer[10],counter[9]);
  fprintf(filePtr, "global communication time is %ef .\n", GlobalCommunication);
  fprintf(filePtr, "MPI_Recv communication time is %ef .\n", RecvComm);
  if ( flagging_count != 0 ) 
    fprintf(filePtr, "average percentage of flagging cells is %ef(= %ef/%d).\n",flagging_pct/flagging_count, flagging_pct, flagging_count);
  else 
    fprintf(filePtr, "average percentage of flagging cells is 0.\n");
  if ( moving_count != 0 )
    fprintf(filePtr, "average percentage of moving cells is %ef(= %ef/%d).\n",moving_pct/moving_count, moving_pct,moving_count);
  else
    fprintf(filePtr, "average percentage of moving cells is 0.\n");

  fprintf(filePtr,"***********************************\n");

  fprintf(filePtr, "square time for transferring region is %ef (%d times).\n", timer[25], counter[5]);
  fprintf(filePtr, "square time for sending particles is %ef (%d times).\n", timer[27], counter[7]);
  fprintf(filePtr, "square time for transferring particles is %ef (%d times).\n", timer[29], counter[9]);
  fprintf(filePtr, "square time for UpdateStarParticleCount is %ef (%d times).\n", timer[31], counter[11]);
  fprintf(filePtr, "square time for transferring Fluxes is %ef (%d times).\n", timer[32], counter[12]);
  fprintf(filePtr, "square time for ShareGrids is %ef (%d times).\n", timer[33], counter[13]);
  fprintf(filePtr, "square time for Transpose is %ef (%d times).\n", timer[34], counter[14]);
  fprintf(filePtr, "square time for BroadcastValue is %ef (%d times).\n", timer[35], counter[15]);
  fprintf(filePtr, "square time for MinValue is %ef (%d times).\n", timer[36], counter[16]);

  fprintf(filePtr, "square time used for RebuildHierarchy is %ef (%d times).\n", timer[21], counter[1]);
  fprintf(filePtr, "square interval time used between RebuildHierarchy is %ef (%d times).\n", timer[20],counter[0]);
  fprintf(filePtr, "square time used for loadbalancing is %ef (%d times).\n", timer[22], counter[2]);
  fprintf(filePtr, "square number of flagging cells is %ef (%d times).\n",timer[24], counter[4]);
  fprintf(filePtr, "square size of region transferred is %ef (%d times). \n", timer[26], counter[5]);
  fprintf(filePtr, "square size of particles sended is %ef (%d times). \n", timer[28],counter[7]);
  fprintf(filePtr, "square size of particles transferred is %ef (%d times). \n", timer[30],counter[9]);

  fclose(filePtr); 

#endif /* MPI_INSTRUMENTATION */

#ifdef TRACE
    fclose(trace_fptr);
#endif
  
  // =========================================================
  // Finalize jbPerf
  // =========================================================
#ifdef USE_JBPERF
  jb::perf.stop("main");
  jb::perf.write ();
  jb::perf.end();
#endif
  // =========================================================

#ifdef USE_PABLO
  HDFendTrace();
#endif
  my_exit(EXIT_SUCCESS);
  }

  catch (std::exception& e) {
    fprintf (stderr,"ENZO ERROR: main() caught STL exception %s!\n",e.what());
    ERROR_MESSAGE;
  }
  catch (...) {
    fprintf (stderr,"ENZO ERROR: main() caught an unknown exception!\n");
    ERROR_MESSAGE;
  }

}


void my_exit(int exit_status)
{
  CommunicationFinalize();
  if (exit_status==EXIT_FAILURE) {
    ERROR_MESSAGE;  // Call MPI_Abort!
  }
  exit(exit_status);
}

