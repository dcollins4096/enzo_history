/***********************************************************************
/
/  AMR MAIN CODE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       January, 2006
/
/  PURPOSE:
/    This is main() for the amr code.  It interprets the arguments and
/    then calls the appropriate routines depending on whether we are
/    doing a new run, a restart, an extraction, or a projection.
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#endif
#include "svn_version.def"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "flowdefs.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#ifdef ISO_GRAV
#include "GravityPotentialBoundary.h"
#endif
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#ifdef RAD_HYDRO
#include "gFLDProblem.h"
#endif
#undef DEFINE_STORAGE
#ifdef PPML
#include "DaveTools.h"
#endif //PPML 
// Function prototypes
 
int InitializeNew(  char *filename, HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary &Exterior, float *Initialdt);
#ifdef ISO_GRAV
#ifdef RAD_HYDRO
int InitializeLocal(int restart, HierarchyEntry &TopGrid, 
		    TopGridData &MetaData, 
		    GravityPotentialBoundary &PotentialBdry,
		    gFLDProblem &FLDsolver);
#else
int InitializeLocal(int restart, HierarchyEntry &TopGrid, 
		    TopGridData &MetaData, 
		    GravityPotentialBoundary &PotentialBdry);
#endif
#else
#ifdef RAD_HYDRO
int InitializeLocal(int restart, HierarchyEntry &TopGrid, 
		    TopGridData &MetaData, gFLDProblem &FLDsolver);
#else
int InitializeLocal(int restart, HierarchyEntry &TopGrid, 
		    TopGridData &MetaData);
#endif
#endif

int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
#ifdef ISO_GRAV
#ifdef RAD_HYDRO
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, 
		    GravityPotentialBoundary *PotentialBdry, 
		    gFLDProblem &FLDsolver,
		    LevelHierarchyEntry *Array[], float Initialdt);
#else
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, 
		    GravityPotentialBoundary *PotentialBdry, 
		    LevelHierarchyEntry *Array[], float Initialdt);
#endif
#else
#ifdef RAD_HYDRO
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, 
		    gFLDProblem &FLDsolver,
		    LevelHierarchyEntry *Array[], float Initialdt);
#else
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, 
		    LevelHierarchyEntry *Array[], float Initialdt);
#endif
#endif
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
int CommunicationInitialize(Eint32 *argc, char **argv[]);
int CommunicationFinalize();
int CommunicationPartitionGrid(HierarchyEntry *Grid);
int ENZO_OptionsinEffect(void);
void my_exit(int status);
#if defined(ISO_GRAV) || defined(RAD_HYDRO)
int DetermineParallelism(HierarchyEntry *TopGrid, TopGridData &MetaData);
#endif

#ifdef MEM_TRACE
Eint64 mused(void);
#endif 
 
 
Eint32 main(Eint32 argc, char *argv[])
{
#ifdef DCC_SCALING
  wall_time ("Start Main");
#endif //DCC_SCALING

#ifdef MEM_TRACE
    Eint64 MemInUse;
#endif
 
  // Initialize Communications
    
  CommunicationInitialize(&argc, &argv);
 
  int int_argc;
  int_argc = argc;
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("=========================\n");
    printf("Enzo SVN Branch   %s\n",ENZO_SVN_BRANCH);
    printf("Enzo SVN Revision %d\n",ENZO_SVN_REVISION);
    printf("=========================\n");
    fflush(stdout);
  }

#ifdef MPI_INSTRUMENTATION
  Start_Wall_Time = MPI_Wtime();
#endif

#ifdef OOC_BOUNDARY
  ExternalBoundaryIO = TRUE;
  ExternalBoundaryTypeIO = FALSE;
  ExternalBoundaryValueIO = FALSE;
#else
  ExternalBoundaryIO = FALSE;
  ExternalBoundaryTypeIO = FALSE;
  ExternalBoundaryValueIO = FALSE;
#endif

  ENZO_OptionsinEffect();
 
  // Main declarations
 
  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
#ifdef ISO_GRAV
  GravityPotentialBoundary PotentialBdry;
#endif
#ifdef RAD_HYDRO
  gFLDProblem FLDsolver;
#endif
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
 
  // Initialize
 
  int restart                  = FALSE,
      OutputAsParticleDataFlag = FALSE,
      InformationOutput        = FALSE,
      project                  = FALSE,
      ProjectionDimension      = INT_UNDEFINED,
      ProjectionSmooth         = FALSE;
  debug                        = FALSE;
  extract                      = FALSE;
  flow_trace_on                = FALSE;
  char *ParameterFile          = NULL;
  char *myname                 = argv[0];
  int RegionStart[MAX_DIMENSION],
      RegionEnd[MAX_DIMENSION],
      RegionLevel;
  FLOAT RegionStartCoordinates[MAX_DIMENSION],
        RegionEndCoordinates[MAX_DIMENSION];
  float Initialdt              = 0;
 
#ifdef FLOW_TRACE
  char pid[5];
  char flow_trace_log[16];
  sprintf(pid, "%4.4"ISYM, MyProcessorNumber);
  strcpy(flow_trace_log, "FlowTrace.");
  strcat(flow_trace_log, pid);
  flow_trace_fptr = fopen( flow_trace_log, "w" );
  flow_trace_level = 0;
  flow_trace_on = TRUE;
#endif
 
  if (flow_trace_on) flow_trace1("ENZO");
 
#ifdef MPI_INSTRUMENTATION
  char perfname[20];
 
  flagging_count = 0;
  moving_count = 0;
  in_count = 0;
  out_count = 0;
  flagging_pct = 0.0;
  moving_pct = 0.0;
 
  GlobalCommunication = 0.0;
  RecvComm = 0.0;
  WaitComm = 0.0;
 
  for (int i=0; i < 40; i++) {
    timer[i] = 0.0;
    counter[i] = 0;
  }
 
  sprintf(perfname, "perfdata_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  perfname[strlen(perfname)] = '\0';
  filePtr = fopen(perfname, "w");
 
  sprintf(tracename, "trace_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  tracename[strlen(tracename)] = '\0';

#ifdef MPI_TRACE
  tracePtr = fopen(tracename, "w");
  traceMPI = TRUE;
#else
  traceMPI = FALSE;
#endif /* MPI_TRACE */

#endif /* MPI_INSTRUMENTATION */

#ifdef MEM_TRACE
  sprintf(memtracename, "mem_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  memtracename[strlen(memtracename)] = '\0';
  memtracePtr = fopen(memtracename, "w");
  traceMEM = TRUE;
#else
  traceMEM = FALSE;
#endif /* MEM_TRACE */

 
  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
 
  // Interpret command-line arguments
 
  // Metadata.commandline: options for the WriteEvolution module.
  // MetaData.commandline.parseAndEat(argc,argv); // parse cmdln options
  // and remove them so as not to confuse InterpretCommandLine
  // MetaData.commandline.print(stderr); // print results of cmdln parse
 
  if (InterpretCommandLine(int_argc, argv, myname, restart, debug, extract,
			   InformationOutput, OutputAsParticleDataFlag,
			   project, ProjectionDimension, ProjectionSmooth,
			   &ParameterFile,
			   RegionStart, RegionEnd,
			   RegionStartCoordinates, RegionEndCoordinates,
			   RegionLevel, MyProcessorNumber) == FAIL) {
    WARNING_MESSAGE;
    fflush(stdout);
    my_exit(EXIT_FAILURE);
  }
 

  // If we need to read the parameter file as a restart file, do it now
 
  if (project || OutputAsParticleDataFlag || extract ||
      InformationOutput || restart) {
 
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


  // Information dump
 
  if (InformationOutput) {
    OutputLevelInformation(stdout, MetaData, LevelArray);
    my_exit(EXIT_SUCCESS);
  }
 
 
  // Extract a grid subsection (doesn't return)
 
  if (extract)
    ExtractSection(TopGrid, MetaData, LevelArray, &Exterior,
		   RegionStart, RegionEnd,
		   RegionStartCoordinates, RegionEndCoordinates, RegionLevel);
 
 
  // Output fields as particle data
 
  if (OutputAsParticleDataFlag)
    if (OutputAsParticleData(MetaData, LevelArray, RegionStart, RegionEnd,
			     RegionStartCoordinates, RegionEndCoordinates,
			     RegionLevel, "amr.particles") == FAIL)
      my_exit(EXIT_FAILURE);
    else
      my_exit(EXIT_SUCCESS);
 
 
  // Project 3D field to 2D plane
 
  if (project)
    if (ProjectToPlane(MetaData, LevelArray, RegionStart, RegionEnd,
		     RegionStartCoordinates, RegionEndCoordinates,
		     RegionLevel, ProjectionDimension, "amr.project",
		     ProjectionSmooth, &Exterior) == FAIL)
      my_exit(EXIT_FAILURE);
    else
      my_exit(EXIT_SUCCESS);
 
 
  // Normal start: Open and read parameter file

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Call initialize %16"ISYM" \n", MemInUse);
#endif
 

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
 
#ifdef USE_JBPERF

    // Initialize jbPerf performance collecting

    jbPerfInitialize(MaximumRefinementLevel);

#endif


#if defined(ISO_GRAV) || defined(RAD_HYDRO)
  // Determine Top-Grid Parallelism Information (store in grids)

  if (DetermineParallelism(&TopGrid, MetaData) == FAIL) {
    fprintf(stderr, "Error in DetermineParallelism.\n");
    my_exit(EXIT_FAILURE);
  }
#endif

  // Perform local initialization (even for restart, may just return, depends on problem)
#ifdef ISO_GRAV
#ifdef RAD_HYDRO
  if (InitializeLocal(restart, TopGrid, MetaData, 
		      PotentialBdry, FLDsolver) == FAIL) {
#else
  if (InitializeLocal(restart, TopGrid, MetaData, PotentialBdry) == FAIL) {
#endif
#else
#ifdef RAD_HYDRO
  if (InitializeLocal(restart, TopGrid, MetaData, FLDsolver) == FAIL) {
#else
  if (InitializeLocal(restart, TopGrid, MetaData) == FAIL) {
#endif
#endif
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Error in Local Initialization.\n");
    my_exit(EXIT_FAILURE);
  }
#ifdef MEM_TRACE
  MemInUse = mused();
  fprintf(memtracePtr, "Call evolve hierarchy %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif
 
  // Call the main evolution routine
 
#ifdef ISO_GRAV
#ifdef RAD_HYDRO
  if (EvolveHierarchy(TopGrid, MetaData, &Exterior, &PotentialBdry, 
		      FLDsolver, LevelArray, Initialdt) == FAIL) {
#else
  if (EvolveHierarchy(TopGrid, MetaData, &Exterior, &PotentialBdry, 
		      LevelArray, Initialdt) == FAIL) {
#endif
#else
#ifdef RAD_HYDRO
  if (EvolveHierarchy(TopGrid, MetaData, &Exterior, FLDsolver, 
		      LevelArray, Initialdt) == FAIL) {
#else
  if (EvolveHierarchy(TopGrid, MetaData, &Exterior, 
		      LevelArray, Initialdt) == FAIL) {
#endif
#endif
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Error in EvolveHierarchy.\n");
    my_exit(EXIT_FAILURE);
  }
  else
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Successful run, exiting.\n");
  
  
  if (flow_trace_on) flow_trace2("ENZO");
#ifdef DCC_SCALING
  wall_time ("End Main");
#endif //DCC_SCALING
  
#ifdef MPI_INSTRUMENTATION
  End_Wall_Time = MPI_Wtime();
  WallTime = End_Wall_Time - Start_Wall_Time;
 
  fprintf(filePtr, "Elapsed wall time:                   %12.6e\n", WallTime);
  fprintf(filePtr, "Communication time:                  %12.6e\n", CommunicationTime);
  fprintf(filePtr, "Global communication time:           %12.6e\n", GlobalCommunication);
  fprintf(filePtr, "Receive communication time:          %12.6e\n", RecvComm);
  fprintf(filePtr, "Waiting communication time:          %12.6e\n", WaitComm);
 
  fprintf(filePtr, "\n\n");
  fprintf(filePtr, "Transferring region       (%8"ISYM" times) %12.6e\n", counter[5],  timer[5]);
  fprintf(filePtr, "Sending particles         (%8"ISYM" times) %12.6e\n", counter[7],  timer[7]);
  fprintf(filePtr, "Transferring particles    (%8"ISYM" times) %12.6e\n", counter[9],  timer[9]);
  fprintf(filePtr, "Transferring Fluxes       (%8"ISYM" times) %12.6e\n", counter[12], timer[12]);
  fprintf(filePtr, "ShareGrids                (%8"ISYM" times) %12.6e\n", counter[13], timer[13]);
  fprintf(filePtr, "Transpose                 (%8"ISYM" times) %12.6e\n", counter[14], timer[14]);
  fprintf(filePtr, "BroadcastValue            (%8"ISYM" times) %12.6e\n", counter[15], timer[15]);
  fprintf(filePtr, "MinValue                  (%8"ISYM" times) %12.6e\n", counter[16], timer[16]);
  fprintf(filePtr, "UpdateStarParticleCount   (%8"ISYM" times) %12.6e\n", counter[11], timer[11]);
 
  fprintf(filePtr, "\n\n");
  fprintf(filePtr, "RebuildHierarchy          (%8"ISYM" times) %12.6e\n", counter[1],  timer[1]);
  fprintf(filePtr, "RebuildHierarchy interval (%8"ISYM" times) %12.6e\n", counter[0],  timer[0]);
  fprintf(filePtr, "Load balancing            (%8"ISYM" times) %12.6e\n", counter[2],  timer[2]);
  fprintf(filePtr, "Region transfer size      (%8"ISYM" times) %12.6e\n", counter[5],  timer[6]);
  fprintf(filePtr, "Particles sent            (%8"ISYM" times) %12.6e\n", counter[7],  timer[8]);
  fprintf(filePtr, "Particle transfer size    (%8"ISYM" times) %12.6e\n", counter[9],  timer[10]);
 
  fprintf(filePtr, "\n\n");
  fprintf(filePtr, "Number of load balancing calls %"ISYM"/%"ISYM" (LOAD_BALANCE_RATIO=%"FSYM")\n",counter[3], counter[2], timer[3]);
  fprintf(filePtr, "Number of flagging cells  (%8"ISYM" times) %12.6e\n", counter[4], timer[4]);
 
  fprintf(filePtr, "\n\n");
  if ( flagging_count != 0 )
    fprintf(filePtr, "Average percentage of flagging cells %12.6e(= %12.6e/%"ISYM")\n",
            flagging_pct/flagging_count, flagging_pct, flagging_count);
  else
    fprintf(filePtr, "Average percentage of flagging cells 0\n");
 
  if ( moving_count != 0 )
    fprintf(filePtr, "Average percentage of moving cells   %12.6e(= %12.6e/%"ISYM")\n",
            moving_pct/moving_count, moving_pct,moving_count);
  else
    fprintf(filePtr, "Average percentage of moving cells 0\n");
 
  fclose(filePtr);
 
#ifdef MPI_TRACE
  fclose(tracePtr);
#endif /* MPI_TRACE */

#endif /* MPI_INSTRUMENTATION */

#ifdef MEM_TRACE
  fclose(memtracePtr);
#endif /* MEM_TRACE */
 
#ifdef FLOW_TRACE
    fclose(flow_trace_fptr);
#endif
 
  my_exit(EXIT_SUCCESS);
 
}
 
void my_exit(int status)
{
  if (status == EXIT_FAILURE) {
    ERROR_MESSAGE;
  } else {
    CommunicationFinalize();
    exit(0);
  }
}
