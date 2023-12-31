/***********************************************************************
/
/  WRITE OUT ALL THE DATA (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       February, 2004
/              October, 2004
/              Direct or indirect SRB driver
/              Local or Global file system, or cwd
/              July , 2006
/              Assemble group file in-core
/  modified2:  Robert Harkness
/  date:       April 2008
/
/  PURPOSE:
/
************************************************************************/
#define SYSCALL
 
// This function writes out the data hierarchy (TopGrid), the External
//   Boundary (Exterior), the TopGridData, and the global_data.
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <hdf5.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "CommunicationUtilities.h"

void my_exit(int status);
 
// HDF5 function prototypes
 

 
// function prototypes
 
int SysMkdir(char *startdir, char *directory);
 
int WriteDataCubes(HierarchyEntry *TopGrid, int TDdims[], char *gridbasename, int &GridID, FLOAT WriteTime);
int Group_WriteDataHierarchy(FILE *fptr, TopGridData &MetaData, HierarchyEntry *TopGrid,
		       char *gridbasename, int &GridID, FLOAT WriteTime, hid_t file_id);
int WriteMemoryMap(FILE *fptr, HierarchyEntry *TopGrid,
		   char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteConfigure(FILE *optr);
int WriteTaskMap(FILE *fptr, HierarchyEntry *TopGrid,
		 char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteParameterFile(FILE *fptr, TopGridData &MetaData);
int WriteStarParticleData(FILE *fptr);
int WriteRadiationData(FILE *fptr);
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int CommunicationCombineGrids(HierarchyEntry *OldHierarchy,
			      HierarchyEntry **NewHierarchyPointer,
			      FLOAT WriteTime);
void DeleteGridHierarchy(HierarchyEntry *GridEntry);
void ContinueExecution(void);
int CreateSmoothedDarkMatterFields(TopGridData &MetaData, HierarchyEntry *TopGrid);
 
 
extern char BCSuffix[];
extern char GridSuffix[];
extern char HierarchySuffix[];
extern char hdfsuffix[];
extern char RadiationSuffix[];
extern char TaskMapSuffix[];
extern char MemoryMapSuffix[];
extern char ConfigureSuffix[];

char CPUSuffix[]       = ".cpu";
 
extern char LastFileNameWritten[MAX_LINE_LENGTH];
 
 
 
 
int Group_WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, FLOAT WriteTime = -1)
{
 
  char id[MAX_CYCLE_TAG_SIZE], *cptr, name[MAX_LINE_LENGTH];
  char dumpdirname[MAX_LINE_LENGTH];
  char dumpdirroot[MAX_LINE_LENGTH];
  char unixcommand[MAX_LINE_LENGTH];
  char gridbasename[MAX_LINE_LENGTH];
  char hierarchyname[MAX_LINE_LENGTH];
  char radiationname[MAX_LINE_LENGTH];
  char taskmapname[MAX_LINE_LENGTH];
  char memorymapname[MAX_LINE_LENGTH];
  char configurename[MAX_LINE_LENGTH];
  char groupfilename[MAX_LINE_LENGTH];
 
  int unixresult;
  int status;
  int local, global;
  int file_status;
  int ii, pe, nn;
 
#ifdef USE_MPI
  double io_start, io_stop;
  double dc_start, dc_stop;
  double ttenter, ttexit;
  double iot1a, iot1b, iot2a, iot2b, iot3a, iot3b, iot4a, iot4b;
  char io_logfile[MAX_NAME_LENGTH];
  FILE *xptr;
#endif /* USE_MPI */
  char pid[MAX_TASK_TAG_SIZE];
 
  FILE *fptr;
  FILE *sptr;
  FILE *gptr;
  FILE *tptr;
  FILE *mptr;
  FILE *optr;
 
  hid_t       file_id;
  hid_t       file_acc_template;
  size_t      memory_increment; // in bytes
  hbool_t     dump_flag;
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int GridID = 1;
  int GridJD = 1;
  int GridKD = 1;
  int GridLD = 1;
 
#ifdef USE_MPI
  ttenter = MPI_Wtime();
#endif /* USE_MPI */

  /* If this is an interpolated time step, then temporary replace  the time
     in MetaData.  Note:  Modified 6 Feb 2006 to fix interpolated  data outputs. */

  FLOAT SavedTime = MetaData.Time;
  MetaData.Time = (WriteTime < 0) ? MetaData.Time : WriteTime;

  /* If we're writing interpolated dark matter fields, create them now. */

  CreateSmoothedDarkMatterFields(MetaData, TopGrid);

  // Global or local filesystem?
 
  local = 0;
  global = 0;
 
  if (MetaData.LocalDir != NULL)
  {
     local = 1;
     strcpy(dumpdirroot, MetaData.LocalDir);
     // fprintf(stdout, "XXXX local dir: %s\n", MetaData.LocalDir);
     // Create on node - locking?
  }
 
  if (MetaData.GlobalDir != NULL)
  {
     global = 1;
     strcpy(dumpdirroot, MetaData.GlobalDir);
     // fprintf(stdout, "XXXX global dir: %s\n", MetaData.GlobalDir);
     // Create on task 0 only
  }
  
  if (( local == 1) && (global == 1))
    fprintf(stdout, "Local AND Global !!\n");
 
  // Create main name
 
  if (ComovingCoordinates && (cptr = strstr(name, "RRRR"))) {
    FLOAT a, dadt;
    CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
    sprintf(cptr, "%"CYCLE_TAG_FORMAT""ISYM, nint(100*((1 + InitialRedshift)/a - 1)));
  } else {
 
    sprintf(id, "%"CYCLE_TAG_FORMAT""ISYM, filenumber);
    sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);

    /******************** CYCLE / DT BASED OUTPUTS ********************/
 
    if ( (cptr = strstr(basename, MetaData.DataDumpName)) ) {
 
      if (MetaData.DataDumpDir != NULL)
      {
        if (MetaData.LocalDir != NULL) {
          // Local fs
          strcpy(dumpdirname, MetaData.LocalDir);
          strcat(dumpdirname, "/");
          strcat(dumpdirname, MetaData.DataDumpDir);
          strcat(dumpdirname, id);
  
          // Create once per node...
#ifdef USE_NODE_LOCAL
            strcat(dumpdirname, "/mpi");
            strcat(dumpdirname, pid);
#endif /* USE_NODE_LOCAL */
 
          strcpy(name, dumpdirname);
          strcat(name, "/");
          strcat(name, basename);
        } // if LocalDir
 
        else
 
        {
          if (MetaData.GlobalDir != NULL) {
            // Global fs
            strcpy(dumpdirname, MetaData.GlobalDir);
            strcat(dumpdirname, "/");
            strcat(dumpdirname, MetaData.DataDumpDir);
            strcat(dumpdirname, id);
            // Do mkdir on cpu 0 only
            strcpy(name, dumpdirname);
            strcat(name, "/");
            strcat(name, basename);
          } // if GlobalDir
 
          else
 
          {
            // No local or global specified
            strcpy(name, basename);
          } // else GlobalDir
 
        } // else LocalDir
 
      } // if DataDumpDir
 
      else
 
      {
        strcpy(name, basename);
      } // else DataDumpDir
 
      if (debug) fprintf(stdout, "DATA dump: %s\n", name);
 
    } // if DataDumpName

    /******************** RESTART BASED OUTPUTS ********************/
 
    if ( (cptr = strstr(basename, MetaData.RedshiftDumpName)) ) {
 
      if (MetaData.RedshiftDumpDir != NULL)
      {
        if (MetaData.LocalDir != NULL) {
          // Local fs
          strcpy(dumpdirname, MetaData.LocalDir);
          strcat(dumpdirname, "/");
          strcat(dumpdirname, MetaData.RedshiftDumpDir);
          strcat(dumpdirname, id);
  
          // Create once per node...
#ifdef USE_NODE_LOCAL
          strcat(dumpdirname, "/mpi");
          strcat(dumpdirname, pid);
#endif /* USE_NODE_LOCAL */
 
          strcpy(name, dumpdirname);
          strcat(name, "/");
          strcat(name, basename);
        } // if LocalDir
 
        else
 
        {
          if (MetaData.GlobalDir != NULL) {
            // Global fs
            strcpy(dumpdirname, MetaData.GlobalDir);
            strcat(dumpdirname, "/");
            strcat(dumpdirname, MetaData.RedshiftDumpDir);
            strcat(dumpdirname, id);
            // Do mkdir on cpu 0 only
            strcpy(name, dumpdirname);
            strcat(name, "/");
            strcat(name, basename);
          } // if GlobalDir
 
          else
 
          {
            // No local or global specified
            strcpy(name, basename);
          } // else GlobalDir
 
        } // else LocalDir
 
      } // if RedshiftDumpDir
 
      else
 
      {
        strcpy(name, basename);
      } // else RedshiftDumpDir

      if (debug)
	fprintf(stdout, "REDSHIFT dump: %s\n", name);
 
    } // if RedshiftDumpName

    /******************** RESTART BASED OUTPUTS ********************/

    if ( (cptr = strstr(basename, MetaData.RestartDumpName)) ) {
 
      if (MetaData.RestartDumpDir != NULL)
      {
        if (MetaData.LocalDir != NULL) {
          // Local fs
          strcpy(dumpdirname, MetaData.LocalDir);
          strcat(dumpdirname, "/");
          strcat(dumpdirname, MetaData.RestartDumpDir);
          strcat(dumpdirname, id);
  
          // Create once per node...
#ifdef USE_NODE_LOCAL
          strcat(dumpdirname, "/mpi");
          strcat(dumpdirname, pid);
#endif /* USE_NODE_LOCAL */
 
          strcpy(name, dumpdirname);
          strcat(name, "/");
          strcat(name, basename);
        } // if LocalDir
 
        else
 
        {
          if (MetaData.GlobalDir != NULL) {
            // Global fs
            strcpy(dumpdirname, MetaData.GlobalDir);
            strcat(dumpdirname, "/");
            strcat(dumpdirname, MetaData.RestartDumpDir);
            strcat(dumpdirname, id);
            // Do mkdir on cpu 0 only
            strcpy(name, dumpdirname);
            strcat(name, "/");
            strcat(name, basename);
          } // if GlobalDir
 
          else
 
          {
            // No local or global specified
            strcpy(name, basename);
          } // else GlobalDir
 
        } // else LocalDir
 
      } // if RedshiftDumpDir
 
      else
 
      {
        strcpy(name, basename);
      } // else RedshiftDumpDir

      if (debug)
	fprintf(stdout, "RESTART dump: %s\n", name);
 
    } // if RestartDumpName
 
    if (filenumber >= 0)
      strcat(name, id);
  }
 
  strcpy(LastFileNameWritten, name);
 
  strcpy(groupfilename, name);
  strcat(groupfilename, CPUSuffix);
  strcat(groupfilename, pid);
 
  if (debug)
    fprintf(stdout, "WriteAllData: writing group file %s\n", groupfilename);
 
//  Synchronization point for directory creation
 
#ifdef USE_MPI
  iot1a = MPI_Wtime();
  CommunicationBarrier();
  dc_start = MPI_Wtime();
  iot1b = MPI_Wtime();
#endif /* USE_MPI */
 
//  Get cwd
//  Generate command
//  Execute system call
 
    if ( local )
    {

      MPI_Arg mpi_size;
      MPI_Arg mpi_rank;

#ifdef USE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else
      mpi_rank = 0;
      mpi_size = 1;
#endif

      pe = mpi_rank;
      nn = mpi_size;

 
      for ( ii = 0; ii < nn; ii++ )
      {

        CommunicationBarrier();
        if( pe == ii )
        {
 
          if ( (cptr = strstr(basename, MetaData.DataDumpName)) ) {
            if (MetaData.DataDumpDir != NULL) {
#ifdef SYSCALL
              unixresult = SysMkdir("", dumpdirname);
              if (debug) fprintf(stdout, "DATA dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
              strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
              unixresult = system(unixcommand);
              if (debug) fprintf(stdout, "DATA dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
            }
          } // ENDIF datadump
 
          if ( (cptr = strstr(basename, MetaData.RedshiftDumpName)) ) {
            if (MetaData.RedshiftDumpDir != NULL) {
#ifdef SYSCALL
              unixresult = SysMkdir("", dumpdirname);
              fprintf(stdout, "REDSHIFT dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
              strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
              unixresult = system(unixcommand);
              fprintf(stdout, "REDSHIFT dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
            }
          } // ENDIF redshift

          if ( (cptr = strstr(basename, MetaData.RestartDumpName)) ) {
            if (MetaData.RestartDumpDir != NULL) {
#ifdef SYSCALL
              unixresult = SysMkdir("", dumpdirname);
              fprintf(stdout, "RESTART dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
              strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
              unixresult = system(unixcommand);
              fprintf(stdout, "RESTART dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
            }
          } // ENDIF restart

        } // ENDIF pe == ii
      } // ENDFOR ii
    } // ENDIF local
 
    if ( global )
    {
      if ( MyProcessorNumber == ROOT_PROCESSOR )
      {
 
        if ( (cptr = strstr(basename, MetaData.DataDumpName)) ) {
          if (MetaData.DataDumpDir != NULL) {
#ifdef SYSCALL
            unixresult = SysMkdir("", dumpdirname);
            if (debug) fprintf(stdout, "DATA dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
            strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
            unixresult = system(unixcommand);
            if (debug) fprintf(stdout, "DATA dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
          }
        } // ENDIF datadump
 
        if ( (cptr = strstr(basename, MetaData.RedshiftDumpName)) ) {
          if (MetaData.RedshiftDumpDir != NULL) {
#ifdef SYSCALL
            unixresult = SysMkdir("", dumpdirname);
            fprintf(stdout, "REDSHIFT dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
            strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
            unixresult = system(unixcommand);
            fprintf(stdout, "REDSHIFT dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
          }
        } // ENDIF redshift

        if ( (cptr = strstr(basename, MetaData.RestartDumpName)) ) {
          if (MetaData.RestartDumpDir != NULL) {
#ifdef SYSCALL
            unixresult = SysMkdir("", dumpdirname);
            fprintf(stdout, "RESTART dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
            strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
            unixresult = system(unixcommand);
            fprintf(stdout, "RESTART dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
          }
        } // ENDIF restart
 
      }
    }


  
//  fprintf(stdout, "Sync point ok\n");
 
#ifdef USE_MPI
  iot2a = MPI_Wtime();
  CommunicationBarrier();
  dc_stop = MPI_Wtime();
  iot2b = MPI_Wtime();
#endif /* USE_MPI */
 
//  Start I/O timing
 
#ifdef USE_MPI
  io_start = MPI_Wtime();
#endif /* USE_MPI */

#ifdef USE_HDF5_OUTPUT_BUFFERING

  memory_increment = 1024*1024;
  dump_flag = 1;

  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
    if( file_acc_template == h5_error ){my_exit(EXIT_FAILURE);}

  h5_status = H5Pset_fapl_core(file_acc_template, memory_increment, dump_flag);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

  file_id = H5Fcreate(groupfilename, H5F_ACC_TRUNC, H5P_DEFAULT, file_acc_template);
    if( file_id == h5_error ){my_exit(EXIT_FAILURE);}

#else

  file_id = H5Fcreate(groupfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //  h5_status = H5Fclose(file_id);
  //    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#endif
 
  // Set MetaData.BoundaryConditionName
 
  if (MetaData.BoundaryConditionName != NULL)
    delete [] MetaData.BoundaryConditionName;
  MetaData.BoundaryConditionName = new char[MAX_LINE_LENGTH];
  strcpy(MetaData.BoundaryConditionName, name);
  strcat(MetaData.BoundaryConditionName, BCSuffix);
 
  // Output TopGrid data
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((fptr = fopen(name, "w")) == NULL) {
      fprintf(stdout, "Error opening output file %s\n", name);
      ENZO_FAIL("");
    }
    if (WriteTime >= 0)
      fprintf(fptr, "# WARNING! Interpolated output: level = %"ISYM"\n",
	      MetaData.OutputFirstTimeAtLevel-1);
    if (WriteParameterFile(fptr, MetaData) == FAIL)
      ENZO_FAIL("Error in WriteParameterFile");
    fclose(fptr);
  
  }

 
  // Output Boundary condition info
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((fptr = fopen(MetaData.BoundaryConditionName, "w")) == NULL) {
      fprintf(stdout, "Error opening boundary condition file: %s\n",
	      MetaData.BoundaryConditionName);
      ENZO_FAIL("");
    }
    strcat(MetaData.BoundaryConditionName, hdfsuffix);
    if (Exterior->WriteExternalBoundary(fptr, MetaData.BoundaryConditionName)
	== FAIL)
      ENZO_FAIL("Error in WriteExternalBoundary");
    fclose(fptr);
 
  }
 
  // Create hierarchy name and grid base name
 
  strcpy(hierarchyname, name);
  strcat(hierarchyname, HierarchySuffix);
 
  strcpy(gridbasename, name);
//  strcat(gridbasename, GridSuffix);

  strcpy(taskmapname, name);
  strcat(taskmapname, TaskMapSuffix);
  strcat(taskmapname, pid);

  strcpy(memorymapname, name);
  strcat(memorymapname, MemoryMapSuffix);

  strcpy(configurename, name);
  strcat(configurename, ConfigureSuffix);

 
  /* Combine the top level grids into a single grid for output
     (TempTopGrid is the top of an entirely new hierarchy). */
 
  HierarchyEntry *TempTopGrid;
  CommunicationCombineGrids(TopGrid, &TempTopGrid, WriteTime);
 
  // Output Data Hierarchy
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    if ((fptr = fopen(hierarchyname, "w")) == NULL) {
      fprintf(stdout, "Error opening hierarchy file %s\n", hierarchyname);
      ENZO_FAIL("");
    }
 
  if (Group_WriteDataHierarchy(fptr, MetaData, TempTopGrid, gridbasename, GridID, WriteTime, file_id) == FAIL)
    ENZO_FAIL("Error in Group_WriteDataHierarchy");

  // At this point all the grid data has been written

  h5_status = H5Fclose(file_id);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#ifdef USE_HDF5_OUTPUT_BUFFERING

  h5_status = H5Pclose(file_acc_template);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#endif


  if (MyProcessorNumber == ROOT_PROCESSOR)
    if ((mptr = fopen(memorymapname, "w")) == NULL) {
      fprintf(stdout, "Error opening memory map file %s\n", memorymapname);
      ENZO_FAIL("");
    }

  if (WriteMemoryMap(mptr, TempTopGrid, gridbasename, GridKD, WriteTime) == FAIL)
    ENZO_FAIL("Error in WriteMemoryMap");

  // Output configure

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((optr = fopen(configurename, "w")) == NULL) {
      fprintf(stdout, "Error opening configure file %s\n", configurename);
      fprintf(stdout, "Not crucial but worrysome. Will continue.\n" );
      //      ENZO_FAIL("");
    }

    WriteConfigure(optr);

    fclose(optr);
  }

  // Output task map

  if ((tptr = fopen(taskmapname, "w")) == NULL) {
    fprintf(stdout, "Error opening task map file %s\n", taskmapname);
    ENZO_FAIL("");
  }

  if (WriteTaskMap(tptr, TempTopGrid, gridbasename, GridLD, WriteTime) == FAIL)
    ENZO_FAIL("Error in WriteTaskMap");
 
  int TGdims[3];
 
  TGdims[0] = MetaData.TopGridDims[0];
  TGdims[1] = MetaData.TopGridDims[1];
  TGdims[2] = MetaData.TopGridDims[2];
 
  //  fprintf(stdout, "TGdims  %"ISYM"  %"ISYM"  %"ISYM"\n", TGdims[0], TGdims[1], TGdims[2]);
 
  if (CubeDumpEnabled == 1)
    if (WriteDataCubes(TempTopGrid, TGdims, name, GridJD, WriteTime) == FAIL)
      ENZO_FAIL("Error in WriteDataCubes");
 
  // Clean up combined top level grid, and first two levels of hierarchy
 
  if (TempTopGrid != TopGrid) {
    if (TempTopGrid->NextGridNextLevel != NULL)
      DeleteGridHierarchy(TempTopGrid->NextGridNextLevel);
    delete TempTopGrid->GridData;
    delete TempTopGrid;
  }
 
  // Output StarParticle data (actually just number of stars)
 
  if (WriteStarParticleData(fptr) == FAIL)
    ENZO_FAIL("Error in WriteStarParticleData");
 
  // Create radiation name and write radiation data
 
  if (RadiationFieldType >= 10 && RadiationFieldType <= 11 &&
      MyProcessorNumber == ROOT_PROCESSOR) {
 
    FILE *Radfptr;
 
    strcpy(radiationname, name);
    strcat(radiationname, RadiationSuffix);
 
    if ((Radfptr = fopen(radiationname, "w")) == NULL) {
      fprintf(stdout, "Error opening radiation file %s\n", radiationname);
      ENZO_FAIL("");
    }
    if (WriteRadiationData(Radfptr) == FAIL)
      ENZO_FAIL("Error in WriteRadiationData");
 
    fclose(Radfptr);
 
  }
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fclose(fptr);
    fclose(mptr);
  }

  fclose(tptr);
 
  // Replace the time in metadata with the saved value (above)
 
  MetaData.Time = SavedTime;

  CommunicationBarrier();
// if (debug)
    //  fprintf(stdout, "WriteAllData: finished writing data\n");
 
//  Stop I/O timing
 
#ifdef USE_MPI
  io_stop = MPI_Wtime();
#endif /* USE_MPI */
 
//  Synchronization point for SRB
 
#ifdef USE_MPI
  iot3a = MPI_Wtime();
  CommunicationBarrier();
  iot3b = MPI_Wtime();
#endif /* USE_MPI */
 
  ContinueExecution();
 
#ifdef USE_MPI
  iot4a = MPI_Wtime();
  CommunicationBarrier();
  iot4b = MPI_Wtime();
#endif /* USE_MPI */

  if ( MyProcessorNumber == ROOT_PROCESSOR ){
    sptr = fopen("OutputLog", "a");
    fprintf(sptr, "DATASET WRITTEN %s \n", name);
    fclose(sptr);
  }
 
#ifdef USE_MPI
  ttexit = MPI_Wtime();
#endif /* USE_MPI */
 
#ifdef USE_MPI
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
  strcpy(io_logfile, "IO_perf.");
  strcat(io_logfile, pid);
  xptr = fopen(io_logfile, "a");
  fprintf(xptr, "IO %12.4e  %s\n", (io_stop-io_start), name);
  fprintf(xptr, "DC %12.4e  %s\n", (dc_stop-dc_start), name);
  fprintf(xptr, "XX %12.4e  %12.4e  %12.4e  %12.4e\n",
                (iot1b-iot1a), (iot2b-iot2a), (iot3b-iot3a), (iot4b-iot4a));
  fprintf(xptr, "TT %12.4e\n", (ttexit-ttenter));
  fclose(xptr);
#endif /* USE_MPI */
 
  //  fprintf(stdout,"Safe exit from WriteAllData\n");
 
  return SUCCESS;
}

/* 
void DeleteGridHierarchy(HierarchyEntry *GridEntry)
{
  if (GridEntry->NextGridThisLevel != NULL)
     DeleteGridHierarchy(GridEntry->NextGridThisLevel);
 
  delete GridEntry;
 
  return;
}
*/
