/* ---------------------------------------------------------------------- *

                 CosmicCollage - Enzo projection stacker
		            Prototype version (serial)

   Written by:  bwo, 28 March 2006 (and thereabouts)

   Purpose:  Converts Enzo projections into a stacked image.  Various details
     are necessary to get this actually work right, of course, and are discussed
     below.  This particular version is actually intended to be a prototype for
     a future OpenMP-parallel version of the code, pending our happiness with 
     it and complaints regarding its overall speed.  Direct all 
     complaints/suggestions/verbal abuse to bwoshea@lanl.gov

   Compile:  Figure out where the HDF5 library is, then type "make" or perhaps
     "gmake."  This will produce a binary called "collage".

   Usage:  The user has to set a bunch of stuff in the file UserDefines.C. 
     Compile, then go to the directory where your enzo projections are and run
     the collage binary.  Warning:  This may require a LOT of memory!

   modified1: BWO, 20 september 2006.  Now this code also allows the user to
     read in files full of cluster information and makes SZ y-decrement maps
     which contain _only_ cluster information and ones where the cluster
     information was masked out.

 * ---------------------------------------------------------------------- */

#include "standard_includes.h"
#define DEFINE_STORAGE
#include "global_variables.h"
#undef DEFINE_STORAGE
#include "subroutines.h"


int main(int argc, char *argv[]){

  FILE *thelogfile;

  /* time variables:  we will keep track of the overall amount of time
     we require to do the various components of the collage-making.
     In principle this will allow us to assess how well it parallelizes,
     and also figure out what needs to be tweaked.  Print out the time
     information incrementally, and also print out a final report at the end
     detailing the main usages of time.  Note that the time outputs only record
     time in one-second increments, so there will be lots of zeros for routines that
     take little time.  */
  time_t start_time, stop_time, total_read_time, total_rebin_time, 
    total_correction_time, total_add_time, total_write_time, 
    total_overall_time, total_leftover_time, total_rebin_writeout_time,
    time_1, time_2, total_cluster_time;

  total_read_time = total_rebin_time = total_correction_time = 
    total_add_time = total_write_time = total_leftover_time = 
    total_rebin_writeout_time = total_cluster_time = time_t(0.0);

  FILE *timefile;

  timefile = fopen("timefile.txt","w");

  // load user defines -- this absolutely, positively must be done first!
  UserDefines();

#ifdef USE_CLUSTER_SPLIT
  // notify the user that the USE_CLUSTER_SPLIT junk is on
  printf("******  USE_CLUSTER_SPLIT is ON!  ******\n");
  fflush(stdout);

  // USE_CLUSTER_SPLIT stuff should _only_ be used with projection type
  // 1 (SZ Y effect) so make sure that we don't ever use it with other 
  // types of light cones.
  if(ProjectionType != 1){
    fprintf(stderr,"USE_CLUSTER_SPLIT is on, but ProjectionType is != 1!  This is wrong!\n");
    exit(-999);
  }
#endif 

  // read in random seed on user flag in UserDefines.C - no error checking, don't be stupid!
  if(read_in_random_seed==1){

    random_seed = atoi(argv[1]);
    printf("Just read in random seed.  Value is %d\n",random_seed);

  }

  // read in number of projections on user flag in UserDefines.C - no error checking, don't be stupid!
  if(read_in_num_projections==1){

    NumberOfProjections = atoi(argv[2]);
    printf("Just read in number of projections.  Value is %d\n",NumberOfProjections);

  }

  printf("opening file %s  %s\n",LogFileName,InputDataFileName);
  fflush(stdout);

  // the thelogfile keeps track of the axis of projection, offsets, etc. for future reference
  thelogfile = fopen(LogFileName,"w");

  start_time = time(NULL);

  // seed random number generators (one for ints, one for doubles)
  // The seed is provided by the user
  srand( random_seed );  // for integars
  srand48( random_seed );  // for doubles

  // output everything to the log file so that we can reconstruct
  // what happened later.
  fprintf(thelogfile,"Misc. information:\n");
  fprintf(thelogfile,"Random seed:                  %d\n", random_seed);
  fprintf(thelogfile,"Input dataset precision:      %d bits\n", InputDatasetPrecision);
  fprintf(thelogfile,"Total number of projections:  %d\n", NumberOfProjections);
  fprintf(thelogfile,"Output projection grid size:  %d\n", OutputProjectionGridSize);
  fprintf(thelogfile,"Projection type:              %d\n", ProjectionType);
  fprintf(thelogfile,"Angular projection size:      %.12lf (degrees)\n", AngularProjSize);
  fprintf(thelogfile,"debug flags:                  %d  %d\n", debug, verbosedebug);
  fprintf(thelogfile,"Input file list:              %s\n",InputProjectionFileNameList);
  fprintf(thelogfile,"Input dataset name:           %s\n",InputProjectionFileDatasetName);
  fprintf(thelogfile,"Input file info:              %s\n",InputDataFileName);
  fprintf(thelogfile,"Output projection name:       %s\n",OutputProjectionFileName);
  fprintf(thelogfile,"Output dataset name:          %s\n",OutputProjectionFileDatasetName);
  fprintf(thelogfile,"\n\n");

  // read in user data from text files
  ReadInUserData();

  // Initialize master projection array and zero it out!
  InitializeProjectionArray();

  // Loop over total number of projections, doing the following things:
  //  1) create and zero out loop arrays (array for read-in and rebinning)
  //  2) choose random axis and x and y offsets
  //  3) read in projection arrays
  //  3a) mask out clusters in projection arrays (optional, for SZ y effect)
  //  4) rebin projection arrays to temporary output grid
  //  5) apply corrections to rebinned array
  //  6) add rebinned and corrected array to master output array
  //  7) delete and set to null loop arrays
  for(int thisproj=0; thisproj < NumberOfProjections; thisproj++){
    int axis;
    double xoffset, yoffset, *xoffsetptr, *yoffsetptr;

    // create and zero projection field array and rebinned field array
    CreateAndZeroLoopArrays(thisproj);
    
    // choose random axis (x,y,z) for projection
    axis = chooseaxis(thisproj);

    // choose random x and y offsets for bottom left corner of input grid (between 0-1)
    xoffsetptr = &xoffset;  // set pointers
    yoffsetptr = &yoffset;
    chooseoffset(xoffsetptr,yoffsetptr, thisproj);  // actually call routine
    if(debug){ printf("x, y offsets:  %lf %lf\n",xoffset,yoffset); fflush(stdout); }

    // print out all relevant information to logfile so that we can reconstruct the
    // projection later, if necessary
    fprintf(thelogfile,"Projection step %d of %d:\n",thisproj+1, NumberOfProjections);
    fprintf(thelogfile,"    Axis:                  %d\n", axis);
    fprintf(thelogfile,"    x,y offsets:           %.12lf  %.12lf\n", xoffset, yoffset);
    fprintf(thelogfile,"    Projection file name:  %s\n", InputProjectionFileNames[axis][thisproj]);
    fprintf(thelogfile,"    Redshift:              %.12lf\n", projredshift[thisproj]);
    fprintf(thelogfile,"    Input Grid ratio:      %.12lf\n", gridratio[thisproj] );
    fprintf(thelogfile,"    Input Grid size:       %d\n", InputProjectionGridSize[thisproj] );
    fprintf(thelogfile,"\n\n" );

    // read in desired field of projection for appropriate axis 
    // (open and close file in here
    time_1 = time(NULL);

    ReadInProjectionArray(axis, thisproj);

    time_2 = time(NULL);
    fprintf(timefile,"spent %e seconds reading in projection, loop %d of %d\n",
	    difftime(time_2,time_1),thisproj+1,NumberOfProjections);
    total_read_time += time_t(difftime(time_2,time_1));

#ifdef USE_CLUSTER_SPLIT
    /* this chunk of the code reads in the cluster file (which contains halo positions,
       SZ y dec., etc) and then creates arrays of "cluster only" and "no cluster" y
       decrement maps. */
    time_1 = time(NULL);

    NumberOfClusters = 0;  // keep track of number of clusters for each redshift bin

    // read cluster file to get halo positions, virial radii, redshifts, y decrement, etc.
    ReadClusterFile(axis, thisproj);

    // modify the arrays which contain y dec to subtract out cluster information in the
    // "no cluster" dataset and add it to the "clusters only" dataset.
    ModifyArraysWithClusterInfo(axis, thisproj);

    // clean up global arrays which are created in ReadClusterFile and set them to null.
    DeleteClusterArrays();

    time_2 = time(NULL);
    fprintf(timefile,"spent %e seconds zeroing out clusters, loop %d of %d\n",
	    difftime(time_2,time_1),thisproj+1,NumberOfProjections);
    total_cluster_time += time_t(difftime(time_2,time_1));

#endif // CLUST_SP

    // rebin read-in projection field in order to add it to the "master" (summed) 
    // projection, taking into account resizing and offsets
    time_1 = time(NULL);

    // Rebin input projection array so that it's the same angular resolution as the output projection array. 
    RebinProjectionArray(thisproj,xoffset,yoffset);

    time_2 = time(NULL);
    fprintf(timefile,"spent %e seconds rebinning projection, loop %d of %d\n",
	    difftime(time_2,time_1),thisproj+1,NumberOfProjections);
    total_rebin_time += time_t(difftime(time_2,time_1));

    // apply various corrections to the rebinned field:  correct for X-ray luminosity, etc.
    // this has the potential to get complicated.
    time_1 = time(NULL);

    // apply corrections to the rebinned array
    ApplyCorrectionsToRebinnedArray(thisproj);

    time_2 = time(NULL);
    fprintf(timefile,"spent %e seconds applying corrections to projection, loop %d of %d\n",
	    difftime(time_2,time_1),thisproj+1,NumberOfProjections);
    total_correction_time += time_t(difftime(time_2,time_1));

    // Write out the rebinned, corrected array
    time_1 = time(NULL);
    WriteOutRebinnedProjectionArray(thisproj);
    time_2 = time(NULL);

    fprintf(timefile,"spent %e seconds writing out rebinned, corrected projections, loop %d of %d\n",
	    difftime(time_2,time_1),thisproj+1,NumberOfProjections);
    total_rebin_writeout_time += time_t(difftime(time_2,time_1));

    // actually add it to the projection!
    time_1 = time(NULL);

    // adds rebinned projections to master projections
    AddRebinnedProjectionToMaster();

    time_2 = time(NULL);
    fprintf(timefile,"spent %e seconds adding projections, loop %d of %d\n\n",
	    difftime(time_2,time_1),thisproj+1,NumberOfProjections);
    total_add_time += time_t(difftime(time_2,time_1));

    // delete and set to null projection field array and rebinned field array
    DeleteLoopArrays();

  }

  time_1 = time(NULL);

  // Write out final projection
  WriteOutProjectionArray();

  time_2 = time(NULL);
  fprintf(timefile,"spent %e seconds writing out final projection\n", difftime(time_2,time_1));
  total_write_time = time_t(difftime(time_2,time_1));

  // Erase projection array
  EraseProjectionArray();

  stop_time = time(NULL);

  // print out final time report
  fprintf(timefile,"------------------- FINAL REPORT ---------------- \n");
  fprintf(timefile,"spent %e total seconds reading files\n", double(total_read_time));
  fprintf(timefile,"spent %e total seconds rebinning projections\n", double(total_rebin_time));
  fprintf(timefile,"spent %e total seconds correcting projections\n", double(total_correction_time));
  fprintf(timefile,"spent %e total seconds writing out corrected projections.\n", double(total_rebin_writeout_time));
  fprintf(timefile,"spent %e total seconds adding together projections\n", double(total_add_time));
  fprintf(timefile,"spent %e total seconds writing out new projection\n", double(total_write_time));

#ifdef USE_CLUSTER_SPLIT
  fprintf(timefile,"spent %e total seconds dealing with masking out clusters\n", double(total_cluster_time));
#endif

  total_overall_time = time_t(difftime(stop_time, start_time));  // this is the overall time the code uses
  fprintf(timefile,"spent %e seconds overall in this program\n", double(total_overall_time));

  // calculate misc. time (that doesn't fall into the bins I just wrote out)
  total_leftover_time = total_overall_time - total_read_time - total_rebin_time
    - total_correction_time - total_add_time - total_write_time;

#ifdef USE_CLUSTER_SPLIT
  total_leftover_time -= total_cluster_time;
#endif

  fprintf(timefile,"spent %e seconds doing misc. things that don't fall into the above categories\n",
	  double(total_leftover_time));

  fclose(timefile);

  fclose(thelogfile);

  // all done!
  return SUCCESS;
}


