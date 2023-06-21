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
#ifndef __TopGridData_h_
#define __TopGridData_h_
/***********************************************************************
/
/  TOP GRID DATA STRUCTURE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#ifdef USE_FLEXIO
#include "Commandline.h"
#endif
struct TopGridData
{

  /* Counters for the TopGrid. */

  int   CycleNumber;         // Number of top grid timestep performed
  FLOAT Time;                // Current problem time
  float CPUTime;             // Current CPU time used

  /* Stopping criteria for TopGrid. */

  FLOAT StopTime;            // time to stop at
  int   StopCycle;           // timestep number to stop at
  float StopCPUTime;         // Maximum CPU time to be used


  /* Parameters governing when analysis is done. */

  FLOAT TimeLastAnalysis;   // CPU time of the last restart dump (seconds)
  FLOAT dtAnalysis;         // CPU time between restart dumps (0 = never)

  /* Parameters governing when output is done. */

  FLOAT TimeLastRestartDump;  // CPU time of the last restart dump (seconds)
  FLOAT dtRestartDump;        // CPU time between restart dumps (0 = never)

  FLOAT TimeLastDataDump;     // Problem time of the last data dump
  FLOAT dtDataDump;           // Problem time between data dumps (0 = never)

  FLOAT TimeLastHistoryDump; // Problem time of the last history (small) dump
  FLOAT dtHistoryDump;       // Problem time between history dumps (0 = never)

  FLOAT TimeLastMovieDump;   // Problem time of the last movie dump
  FLOAT dtMovieDump;          // Problem time between movie dumps (0 = never)

  FLOAT MovieRegionLeftEdge[MAX_DIMENSION];  // region for movie output
  FLOAT MovieRegionRightEdge[MAX_DIMENSION];

  int CycleLastRestartDump;  // Cycle of the last restart dump (seconds)
  int CycleSkipRestartDump;  // Cycles between restart dumps (0 = never)

  int CycleLastDataDump;     // Cycle of the last data dump
  int CycleSkipDataDump;     // Number of cycles between data dumps (0 = never)

  int CycleLastHistoryDump;  // Cycle of the last history (small) dump
  int CycleSkipHistoryDump;  // Cycles between history dumps (0 = never)
  int CycleSkipGlobalDataDump;//AK Cycles between global data dumps (0 = never)

  int OutputFirstTimeAtLevel; // Outputs when a new level is generated
  int StopFirstTimeAtLevel;   // Stops when this level is first reached

  /* Parameters governing output names. */

  int RestartDumpNumber;    // number appended to end of restart dump name
  int DataDumpNumber;       // number appended to end of data dump name
  int HistoryDumpNumber;    // number appended to end of history dump
  int MovieDumpNumber;      // number appended to end of movie dump name

  char *RestartDumpName;    // restart dump base name
  char *DataDumpName;       // data dump base name
  char *HistoryDumpName;    // history dump base name
  char *MovieDumpName;      // movie dump base name
  char *RedshiftDumpName;   // redshift dump base name

  char *AnalysisName;        // analysis dump base name
  /* TopGrid Parameters governing hierarchy */

  int StaticHierarchy;     // TRUE for static mesh refinement

  /* Some grid defining data
     These are here out of convenience, the real ones are in the grids. */

  int TopGridRank;
  int TopGridDims[MAX_DIMENSION];
  boundary_type  LeftFaceBoundaryCondition[MAX_DIMENSION],
                RightFaceBoundaryCondition[MAX_DIMENSION];
  char *BoundaryConditionName;

  /* Gravity data.
     These are here out of convenience, the real ones are in the grids. */

  gravity_boundary_type GravityBoundary;

  /* Particle and Particle boundary data. (real one in ExternalBoundary). */

  boundary_type ParticleBoundaryType;
  int           NumberOfParticles;

  /* Hydro Parameters.  
     These are here out of convenience, the real ones are in the grids. */

  float  CourantSafetyNumber;                       // Hydro parameter
  int    PPMFlatteningParameter;                    // PPM parameter
  int    PPMDiffusionParameter;                     // PPM parameter
  int    PPMSteepeningParameter;                    // PPM parameter
#ifdef USE_FLEXIO
  Commandline commandline; // Commandline options for WriteEvolution module.
#endif

  ProjectionParameters *ProjectionParams;
  PDFGeneratorParameters *PDFGeneratorParams;

};

#endif
