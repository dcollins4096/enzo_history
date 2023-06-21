#include "standard_includes.h"
#include "global_variables.h"
#include "subroutines.h"


/* ---------------------------------------------------------------------- *
   The user sets all of their parameters in this routine.  See comments 
   below for more information.
 * ---------------------------------------------------------------------- */
void UserDefines(void){

  /* -----------  USER SETS ALL OF THE PARAMETERS BELOW  ---------------- */
  random_seed = 42273645;   // must be a positive int, less than nine digits
                            // this is used for both position shift and axis

  read_in_random_seed = 0;  // set to 1 or 0.  If set to 1, read in the random seed from the command line.

  NumberOfProjections = 26;  // number of individual redshifts at which we
                             // have dumped out data (int value)

  read_in_num_projections = 0;  // set to 0 or 1.  If set to 1, read in the random seed from the command line.

  read_in_random_shifts = 0;  // set to 0 or 1.  If set to 1, read in the random shifts and axis from a file
                              // called "random_shifts.txt"

  OutputProjectionGridSize = 2048;  // pixels per edge of the output projection
                                    // (int value)

  ProjectionType = 1; // projection type: 
                      //   1 is SZ Y effect
                      //   2 is SZ Kinetic effect
                      //   3 is Baryon density
                      //   4 is Dark Matter density
                      //   5 is X-ray luminosity
                      // (int value)

  // size of projection in degrees on the sky 
  // (degrees on a side, not square degrees)
  // (floating point value)
  AngularProjSize = 10.0;  

  // in Comoving Mpc/h - at some point this is going to have to be changed
  // to a line in the read-in file
  ComovingBoxSize = 512.0;  

  // in units of 100 km/s/Mpc - this will also have to be changed to be
  // read in at some point.
  Hubble0 = 0.7;  

  // File containing list of input projection names - this is explained 
  // more fully in the comments above the function ReadInUserData
  // in ProjectionRoutines.C
  // (character string)
  InputProjectionFileNameList = "filenames.txt";

  // precision of input projections (64 or 32 bits).  Outputs are
  // at present always 32 bits.
  InputDatasetPrecision = 64; // can be 32 or 64 (bits)

  // Projection dataset name (in input file)
  // (character string)

  InputProjectionFileDatasetName = "/SZ_Y_Effect";
  //InputProjectionFileDatasetName = "/X_Ray_Luminosity";

  // file containing list of input projection information.  This is
  // explained more fully in comments above the function ReadInUserData
  // in ProjectionRoutines.C
  // (character string)
  InputDataFileName = "inputdata.txt";

  // This file contains a list of the file names of the output,
  // rebinned projection arrays (which will be summed up to
  // make the output light cone)
  RebinFileNameList = "rebinned_projection_names.txt";

#ifdef USE_CLUSTER_SPLIT

  // file containing the names of rebinned tiles for the "no cluster" output.
  RebinNoClustersFileNameList = "rebinned_no_clusters_projection_names.txt";

  // file containing the names of rebinned tiles for the "cluster only" output
  RebinJustClustersFileNameList = "rebinned_just_clusters_projection_names.txt";

  // file containing the names of the "cluster files" - those containing all of
  // the halo information such as position and y-decrement
  ClusterPositionFileNameList = "cluster_position_files.txt";

  // projection name for the "no clusters" light cone
  NoClustersProjectionFileName = "SZY.NoClusters.hdf5";

  // projection name for the "clusters only" light cone.
  JustClustersProjectionFileName = "SZY.JustClusters.hdf5";

#endif // USE_CLUSTER_SPLIT

  // Name of output projection file (HDF5 file)
  // (character string)

  OutputProjectionFileName = "SZY.lightcone.hdf5";
  //OutputProjectionFileName = "Xray.lightcone.hdf5";

  // Name of HDF5 dataset in output projection file
  // (character string)
  OutputProjectionFileDatasetName = "/Projection_Dataset";

  // name of log file that axes, offsets, etc. are written to
  LogFileName = "lightcone_log.txt";

  // debug flag:  1 = on, 0 = off
  // (int value)
  debug = 1;    

  // verbose debug flag:  1 = on, 0 = off
  // this is VERY noisy and will hugely slow down the calculation!
  // (int value)
  verbosedebug = 0; 

  return;
}

