#include "standard_includes.h"
#include "global_variables.h"
#include "subroutines.h"

/* -------------------------------------------------------------------- *
   This routine initializes the master projection array and zeros it
   out.  The initializations are broken up as much as possible to simplify
   later finding of seg faults when somebody inevitably tries to run this
   on a system with far too little memory.
   Note that the master projection array is the one that will be written
   out to disk at the end!

   If the USE_CLUSTER_SPLIT stuff is enabled, it initializes and zeros
   the "no clusters" and "just clusters" arrays as well.
 * -------------------------------------------------------------------- */
void InitializeProjectionArray(void){

  if(debug){
    printf("InitializeProjectionArray: about to declare master projection array: total size %d bytes\n",
	   OutputProjectionGridSize*OutputProjectionGridSize*sizeof(double) );
    fflush(stdout);
  }

  // actually declare the array
  MasterProjectionArray = new float[OutputProjectionGridSize*OutputProjectionGridSize];

#ifdef USE_CLUSTER_SPLIT
  MasterProjectionArrayNoClusters = new float[OutputProjectionGridSize*OutputProjectionGridSize];
  MasterProjectionArrayJustClusters = new float[OutputProjectionGridSize*OutputProjectionGridSize];
#endif // USE_CLUSTER_SPLIT

  if(debug){ printf("InitializeProjectionArray: Array initialized, zeroing.\n"); fflush(stdout); }

  // Loop over newly-initialized array and zero it out.  In principle, if we don't
  // have enough memory yet try to declare the array anyway, the code should crash. 
  // However, this tends to happen on some systems only when the declared (and too-large)
  // array is written to.  So if things crash at this line, probably the user is attempting
  // to use more memory than they have available to them.  The "unlimit memoryuse" command
  // is your friend.
  for(int i=0; i<OutputProjectionGridSize*OutputProjectionGridSize; i++)
#ifdef USE_CLUSTER_SPLIT
    MasterProjectionArrayNoClusters[i] = MasterProjectionArrayJustClusters[i] = 
#endif // USE_CLUSTER_SPLIT
    MasterProjectionArray[i] = 0.0;

  if(debug){ printf("InitializeProjectionArray: Done zeroing array, returning to main.\n");  fflush(stdout); }

  return;
}



/* -------------------------------------------------------------------- *
   Frees up memory and sets to null the master projection array pointer.  
   Probably unnecessary, but good coding form.  Does the same for the 
   "no cluster" and "just cluster" arrays if necessary.
 * -------------------------------------------------------------------- */
void EraseProjectionArray(void){
  if(debug){ printf("EraseProjectionArray: entering.\n"); fflush(stdout); }

  // free up memory usage
  delete [] MasterProjectionArray;

#ifdef USE_CLUSTER_SPLIT
  delete [] MasterProjectionArrayNoClusters;
  delete [] MasterProjectionArrayJustClusters;
#endif // USE_CLUSTER_SPLIT

  // set it to null
  MasterProjectionArray = NULL;

#ifdef USE_CLUSTER_SPLIT
  MasterProjectionArrayNoClusters = NULL;
  MasterProjectionArrayJustClusters = NULL;
#endif // USE_CLUSTER_SPLIT

  if(debug){ printf("EraseProjectionArray: leaving.\n");  fflush(stdout); }
  return;
}



/* -------------------------------------------------------------------- *
   Declares the arrays which contain the rebinned data and input projection
   data.  See ranting in InitializeProjectionArray for potentially useful
   information if you're seeing seg faults.  Does the same thing for 
   the "no cluster" and "cluster only" routines if that option is enabled.
 * -------------------------------------------------------------------- */
void CreateAndZeroLoopArrays(int thisprojection){

  if(debug){ 
    printf("CreateAndZeroLoopArrays: about to declare rebin array: total size %d bytes\n",
	   OutputProjectionGridSize*OutputProjectionGridSize*sizeof(float) );
    fflush(stdout);
  }

  // declare rebinned array - this is the size of the output array!
  RebinArray = new float[OutputProjectionGridSize*OutputProjectionGridSize];

#ifdef USE_CLUSTER_SPLIT
  RebinArrayNoClusters = new float[OutputProjectionGridSize*OutputProjectionGridSize];
  RebinArrayJustClusters = new float[OutputProjectionGridSize*OutputProjectionGridSize];
#endif // USE_CLUSTER_SPLIT

  if(debug){ 
    printf("CreateAndZeroLoopArrays: Rebin array initialized, zeroing.\n"); 
    fflush(stdout); 
  }

  // zero the array out.  If you see a seg fault here, you probably don't
  // have enough memory allocated, for reasons which are discussed in
  // the InitializeProjectionArray routine, above
  for(int i=0; i<OutputProjectionGridSize*OutputProjectionGridSize; i++)
#ifdef USE_CLUSTER_SPLIT
    RebinArrayNoClusters[i] = RebinArrayJustClusters[i] =
#endif
    RebinArray[i] = 0.0;

  if(debug){ 
    printf("CreateAndZeroLoopArrays: Done zeroing array, moving to next one.\n");

    if(InputDatasetPrecision==32){
      printf("CreateAndZeroLoopArrays: about to declare input projection array: total size %d bytes\n",
	     OutputProjectionGridSize*OutputProjectionGridSize*sizeof(float) );
    } else {
      printf("CreateAndZeroLoopArrays: about to declare input projection array: total size %d bytes\n",
	     OutputProjectionGridSize*OutputProjectionGridSize*sizeof(double) );
    }

    fflush(stdout);
  }

  // declare input projection array - this does not have to be the same size as the
  // output array, and can be either float or double precision
  if(InputDatasetPrecision==32){  // 32 bit (float) inputs

    InputProjectionArrayFloat = new float[InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection]];

#ifdef USE_CLUSTER_SPLIT
    InputProjectionArrayNoClustersFloat = new float[InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection]];
    InputProjectionArrayJustClustersFloat = new float[InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection]];
#endif // USE_CLUSTER_SPLIT

  } else if(InputDatasetPrecision==64){  // 64 bit (double) inputs

    InputProjectionArrayDouble = new double[InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection]];

#ifdef USE_CLUSTER_SPLIT
    InputProjectionArrayNoClustersDouble = new double[InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection]];
    InputProjectionArrayJustClustersDouble = new double[InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection]];
#endif // USE_CLUSTER_SPLIT

  } else {  // if it's not 32 or 64, it's wrong!

    fprintf(stderr,"CreateAndZeroLoopArrays:  Precision is %d, expect 32 or 64.  Exiting.\n",InputDatasetPrecision);
    exit(-1);
  }
  if(debug){ 

    printf("CreateAndZeroLoopArrays: Array initialized, zeroing.\n");
    fflush(stdout);
  }

  // zero array out.  If you see a seg fault here, you probably don't
  // have enough memory allocated, for reasons which are discussed above.
  if(InputDatasetPrecision==32){

    for(int i=0; i<InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection]; i++)
#ifdef USE_CLUSTER_SPLIT
      InputProjectionArrayNoClustersFloat[i] = InputProjectionArrayJustClustersFloat[i] = 
#endif // USE_CLUSTER_SPLIT
      InputProjectionArrayFloat[i] = 0.0;

  } else if(InputDatasetPrecision==64){

    for(int i=0; i<InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection]; i++)
#ifdef USE_CLUSTER_SPLIT
      InputProjectionArrayNoClustersDouble[i] = InputProjectionArrayJustClustersDouble[i] = 
#endif // USE_CLUSTER_SPLIT
      InputProjectionArrayDouble[i] = 0.0;

  } else {

    fprintf(stderr,"problems with input dataset precision, can only be 32 or 64:  %d\n",InputDatasetPrecision);
    exit(-111);

  }

  if(debug){ 
    printf("CreateAndZeroLoopArrays: Done zeroing second array, returning to main.\n");
    fflush(stdout);
  }

  return;
}



/* -------------------------------------------------------------------- *
   Frees up memory from arrays which were declared in CreateAndZeroLoopArrays.
   These are the arrays which contain the rebinned data and the input 
   projection data.  Clean up "no cluster" and "cluster only" arrays if
   that option is enabled.
 * -------------------------------------------------------------------- */
void DeleteLoopArrays(void){

  if(debug){ printf("DeleteLoopArrays: entering\n"); fflush(stdout); }

  delete [] RebinArray;
  RebinArray = NULL;

#ifdef USE_CLUSTER_SPLIT
  delete [] RebinArrayNoClusters;
  delete [] RebinArrayJustClusters;

  RebinArrayNoClusters = NULL;
  RebinArrayJustClusters = NULL;
#endif // USE_CLUSTER_SPLIT

  if(InputDatasetPrecision==32){

    delete [] InputProjectionArrayFloat;
    InputProjectionArrayFloat = NULL;

#ifdef USE_CLUSTER_SPLIT
    delete [] InputProjectionArrayNoClustersFloat;
    delete [] InputProjectionArrayJustClustersFloat;

    InputProjectionArrayNoClustersFloat = NULL;
    InputProjectionArrayJustClustersFloat = NULL;
#endif // USE_CLUSTER_SPLIT

  } else {

    delete [] InputProjectionArrayDouble;
    InputProjectionArrayDouble = NULL;

#ifdef USE_CLUSTER_SPLIT
    delete [] InputProjectionArrayNoClustersDouble;
    delete [] InputProjectionArrayJustClustersDouble;

    InputProjectionArrayNoClustersDouble = NULL;
    InputProjectionArrayJustClustersDouble = NULL;
#endif // USE_CLUSTER_SPLIT

  }

  if(debug){ printf("DeleteLoopArrays: exiting\n"); fflush(stdout); }

  return;
}



/* -------------------------------------------------------------------- *
   This simply adds the rebinned, corrected projection array to the
   master projection array.  This routine is amenable to parallelization.
 * -------------------------------------------------------------------- */
void AddRebinnedProjectionToMaster(void){

  if(debug){ printf("AddRebinnedProjectionToMaster: entering\n"); fflush(stdout); }  

  for(int i=0; i< OutputProjectionGridSize*OutputProjectionGridSize; i++){
    MasterProjectionArray[i] += RebinArray[i];
#ifdef USE_CLUSTER_SPLIT
    MasterProjectionArrayNoClusters[i] += RebinArrayNoClusters[i];
    MasterProjectionArrayJustClusters[i] += RebinArrayJustClusters[i];
#endif // USE_CLUSTER_SPLIT
  }

  if(debug){ printf("AddRebinnedProjectionToMaster: exiting\n"); fflush(stdout); }  

  return;
}



/* -------------------------------------------------------------------- *
   Writes out the projection array in hdf5 format.  This can be easily
   modified to other (less pedantic) formats if necessary or desired.
   Does the same thing for the "no clusters" and "cluster only" arrays
   if that option is enabled.
 * -------------------------------------------------------------------- */
void WriteOutProjectionArray(void){

  if(debug){ printf("WriteOutProjectionArray:  entering.\n"); fflush(stdout); }

  // hdf 5 declarations - only one set necessary
  hid_t outputfile_id;
  hsize_t outputdims[2];
  herr_t h5_status;
  herr_t h5_error = -1;

  hid_t output_dset_id, output_dsp_id;

  /* ------------------ write out master projection ----------------- */

  // declare dimensions of output array
  outputdims[0] = outputdims[1] = OutputProjectionGridSize;

  // create output file
  outputfile_id = H5Fcreate(OutputProjectionFileName, H5F_ACC_TRUNC, 
			    H5P_DEFAULT, H5P_DEFAULT);
  if(debug){ 
    printf("Output file is named %s\n",OutputProjectionFileName);
    printf("status for outputfile_id is %d\n",(int) outputfile_id);
    fflush(stdout);
  }

  assert( outputfile_id != h5_error );

  // -------------- create dataspace
  output_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(debug){ 
    printf("status for output_dsp_id is %d\n",(int) output_dsp_id);
    fflush(stdout);
  }

  assert( output_dsp_id != h5_error);

  // -------------- create dataset
  output_dset_id=H5Dcreate(outputfile_id, OutputProjectionFileDatasetName, 
			   HDF5_FILE_I4, output_dsp_id, H5P_DEFAULT);
  if(debug){ 
    printf("read status for output_dset_id is: %d\n",(int) output_dset_id);
    fflush(stdout);
  }

  assert(output_dset_id != h5_error);

   // -------------- write dataset
  h5_status = H5Dwrite(output_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT, MasterProjectionArray);
  if(debug){
    printf("status for H5Dwrite is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);


  // -------------- close dataset 
  h5_status = H5Dclose(output_dset_id);
  if(debug){ 
    printf("status for H5Dclose is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);

  // -------------- close file
  h5_status = H5Fclose(outputfile_id);
  if(debug){
    printf("status for H5Fclose is %d\n",(int) h5_status);
    fflush(stdout);
  }
  assert( h5_status != h5_error );

#ifdef USE_CLUSTER_SPLIT

  /* ------------------ write out "no cluster" projection ---------------- */

  // create output file
  outputfile_id = H5Fcreate(NoClustersProjectionFileName, H5F_ACC_TRUNC, 
			    H5P_DEFAULT, H5P_DEFAULT);
  if(debug){ 
    printf("Output file is named %s\n",NoClustersProjectionFileName);
    printf("status for outputfile_id is %d\n",(int) outputfile_id);
    fflush(stdout);
  }

  assert( outputfile_id != h5_error );

  // -------------- create dataspace
  output_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(debug){ 
    printf("status for output_dsp_id is %d\n",(int) output_dsp_id);
    fflush(stdout);
  }

  assert( output_dsp_id != h5_error);

  // -------------- create dataset
  output_dset_id=H5Dcreate(outputfile_id, OutputProjectionFileDatasetName, 
			   HDF5_FILE_I4, output_dsp_id, H5P_DEFAULT);
  if(debug){ 
    printf("read status for output_dset_id is: %d\n",(int) output_dset_id);
    fflush(stdout);
  }

  assert(output_dset_id != h5_error);

   // -------------- write dataset
  h5_status = H5Dwrite(output_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		       MasterProjectionArrayNoClusters);
  if(debug){
    printf("status for H5Dwrite is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);


  // -------------- close dataset 
  h5_status = H5Dclose(output_dset_id);
  if(debug){ 
    printf("status for H5Dclose is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);

  // -------------- close file
  h5_status = H5Fclose(outputfile_id);
  if(debug){
    printf("status for H5Fclose is %d\n",(int) h5_status);
    fflush(stdout);
  }
  assert( h5_status != h5_error );



  /* ------------------ write out "just clusters" projection ---------------- */

  // create output file
  outputfile_id = H5Fcreate(JustClustersProjectionFileName, H5F_ACC_TRUNC, 
			    H5P_DEFAULT, H5P_DEFAULT);
  if(debug){ 
    printf("Output file is named %s\n",JustClustersProjectionFileName);
    printf("status for outputfile_id is %d\n",(int) outputfile_id);
    fflush(stdout);
  }

  assert( outputfile_id != h5_error );

  // -------------- create dataspace
  output_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(debug){ 
    printf("status for output_dsp_id is %d\n",(int) output_dsp_id);
    fflush(stdout);
  }

  assert( output_dsp_id != h5_error);

  // -------------- create dataset
  output_dset_id=H5Dcreate(outputfile_id, OutputProjectionFileDatasetName, 
			   HDF5_FILE_I4, output_dsp_id, H5P_DEFAULT);
  if(debug){ 
    printf("read status for output_dset_id is: %d\n",(int) output_dset_id);
    fflush(stdout);
  }

  assert(output_dset_id != h5_error);

   // -------------- write dataset
  h5_status = H5Dwrite(output_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		       MasterProjectionArrayJustClusters);
  if(debug){
    printf("status for H5Dwrite is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);


  // -------------- close dataset 
  h5_status = H5Dclose(output_dset_id);
  if(debug){ 
    printf("status for H5Dclose is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);

  // -------------- close file
  h5_status = H5Fclose(outputfile_id);
  if(debug){
    printf("status for H5Fclose is %d\n",(int) h5_status);
    fflush(stdout);
  }
  assert( h5_status != h5_error );


#endif

  if(debug){ printf("WriteOutProjectionArray:  leaving.\n"); fflush(stdout); }

  return;
}



/* -------------------------------------------------------------------- *
   Reads in the projection array, assuming it's a standard HDF5 dataset.
   This assumes that the filenames were given by the user in a text file
   of some sort, a while back.
 * -------------------------------------------------------------------- */
void ReadInProjectionArray(int axis, int boxnumber){
  if(debug){
    printf("ReadInProjectionArray: entering.  %d  %d\n", axis, boxnumber);
    fflush(stdout);
  }

  // we only need one of all of the following
  hid_t file_id;
  hid_t       mem_type_id;
  hid_t       mem_dsp_id, file_dsp_id;
  hid_t       dsp_id;
  hid_t       typ_id;

  hsize_t     size;
  hsize_t     dims[4];
  hsize_t     xdims[4];
  hsize_t     maxdims[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  // need one of these per dataset!
  hid_t input_dset_id;

  int i, ndims;

  // open file - only once
  if(debug){
    printf("ReadInProjectionArray:  about to open file %s\n",InputProjectionFileNames[axis][boxnumber]);    
    fflush(stdout);
  }
  file_id = H5Fopen(InputProjectionFileNames[axis][boxnumber], H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );

  // open dataset dataset
  input_dset_id = H5Dopen(file_id, InputProjectionFileDatasetName);
  assert( input_dset_id != h5_error );

  // open dataspace (to get dimensions) 
  dsp_id = H5Dget_space(input_dset_id);
  assert( dsp_id != h5_error );

  // get data type
  typ_id = H5Dget_type(input_dset_id);
  assert( typ_id != h5_error );

  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

 // from the dimensional information, calculate the size of the buffer.
  size = 1;
  if(debug){ printf("Ndims %d\n",ndims); fflush(stdout); }

  // calculate size of dataset
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(debug){ printf(" Dim %d\n", (int) xdims[i]);  fflush(stdout); }
    }
  if(debug){ printf("Size %d\n", (int) size); fflush(stdout); }

  // make sure everything is the expected size!
  assert( (InputProjectionGridSize[boxnumber]*InputProjectionGridSize[boxnumber]) == ((int) size) );

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
  assert( file_dsp_id != h5_error );

  mem_dsp_id = H5Screate_simple(1, &size, NULL);
  assert( mem_dsp_id != h5_error );

  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) ){

    fprintf(stderr,"opening 32 bit!\n");

    mem_type_id = H5T_NATIVE_FLOAT;

    // read input field into an array
    h5_status = H5Dread(input_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, InputProjectionArrayFloat);
    if(debug){ printf("float read status %d for input field\n", 
		      (int) h5_status); fflush(stdout); }
    assert( h5_status != h5_error ); 

  } else if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) ){
      
    fprintf(stderr,"opening 64 bit!\n");

    mem_type_id = H5T_NATIVE_DOUBLE;

    // read input field into an array
    h5_status = H5Dread(input_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, InputProjectionArrayDouble);
    if(debug){ printf("float read status %d for input field\n", 
		      (int) h5_status); fflush(stdout); }
    assert( h5_status != h5_error ); 
  } else {
    fprintf(stderr,"problems reading dataset, aborting.\n");
    exit(-123);
  }


  // ---------- close hdf5 file, doing appropriate error checking
  h5_status = H5Dclose(input_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  if(debug){
    printf("ReadInProjectionArray: exiting.\n");
    fflush(stdout);
  }

#ifdef USE_CLUSTER_SPLIT

  // copy over all of the read-in data to the "NoClusters" dataset.  We will then
  // subtract the cluster from that datset, and add them to the "Clusters" dataset.
  if(InputDatasetPrecision==32){
    for(i=0; i<int(size); i++)
      InputProjectionArrayNoClustersFloat[i] = InputProjectionArrayFloat[i];
  } else if (InputDatasetPrecision==64){
    for(i=0; i<int(size); i++)
      InputProjectionArrayNoClustersDouble[i] = InputProjectionArrayDouble[i];
  } else {
    fprintf(stderr,"aaaah!  wrong precision!\n");
    exit(-398);
  }

#endif // USE_CLUSTER_SPLIT

  return;
}



/* -------------------------------------------------------------------- *
   This routine reads in the user data from three separate files (not 
   counting the files used for the cluster positioning stuff - see below).  

   The first file is the list of projections, and assumes that the reader 
   has put the file names in this order:

   1.  Box 1 x-axis
   2.  Box 1 y-axis
   3.  Box 1 z-axis
   4.  Box 2 x-axis
   5.  Box 2 y-axis
   6.  Box 2 z-axis

   etc., up to Box N z-axis.  Therefore, there are 3*NumberOfProjections
   lines in this file (where NumberOfProjections is set in UserDefines.C).

   The second file contains columns of ascii data with the following
   information:

   1)  Projection redshift (redshift at which corrections should be 
       applied)
   2)  grid ratio (defined as the ratio of the angular size that this box
       subtends on the sky to AngularProjSize, which is set in 
       UserDefines.C.  Values < 1 indicate a angular size less than
       AngularProjSize, values > 1 indicate an angular size greater than
       that.
   3)  Projection size in terms of number of cells along one edge.  We
       always assume that our projections are squares.
 
   The third file contains a list of the file names for the rebinned,
   corrected tiles which will get added up to become the light cone.  We
   want the individual tiles for debugging purposes and as yet to be 
   determined analysis.

   Note that when we want to cut out clusters from our dataset, we have
   to read in three more files, with the names of the individual output
   data dumps for the "No clusters" and "clusters only" datasets, and
   with the list of names for the input clusters.  These are in the
   ifdef statements at the bottom.

 * -------------------------------------------------------------------- */
void ReadInUserData(void){

  FILE *fptr;

  char *line = new char[MAX_LINE_LENGTH];

  // open up list of projection files
  fptr = fopen(InputProjectionFileNameList,"r");

  // error check
  assert( NumberOfProjections <= MAX_PROJECTION_FILES );

  // loop through total number of projection files, assuming user put them in
  // the right order (no real way to check this, but we could in principle do 
  // some input tests here)
  for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++)
    for(int thisaxis = 0; thisaxis < 3; thisaxis++){

      InputProjectionFileNames[thisaxis][thisprojfile] = new char[MAX_LINE_LENGTH];      

      fgets(line, MAX_LINE_LENGTH, fptr);

      sscanf(line,"%s",InputProjectionFileNames[thisaxis][thisprojfile]);
    }


  // if user wants debug information, print out the projection file, etc.
  if(verbosedebug){
    for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++)
      for(int thisaxis = 0; thisaxis < 3; thisaxis++){

	printf("%d %d -- %s\n",thisprojfile, thisaxis, InputProjectionFileNames[thisaxis][thisprojfile]);

      }  
    fflush(stdout);
  }

  fclose(fptr);

  
  // open up file with columns of user-defined data
  fptr = fopen(InputDataFileName,"r");

  // loop through total number of projection outputs and read in redshift, grid ratio, grid size
  for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++){

      fgets(line, MAX_LINE_LENGTH, fptr);

      // actually read in data
      sscanf(line,"%lf %lf %d", &projredshift[thisprojfile],
	     &gridratio[thisprojfile], &InputProjectionGridSize[thisprojfile]);

      // very basic error checking on gridratio
      if(gridratio[thisprojfile] <= 0.0){
	fprintf(stderr,"error in ReadInUserData:  gridratio[%d] = %lf is less than or equal to zero!\n",
		thisprojfile,gridratio[thisprojfile]);
	exit(-1);
      }

      // very basic error checking on redshifts
      if(projredshift[thisprojfile] < 0.0){
	fprintf(stderr,"error in ReadInUserData:  projredshift[%d] = %lf is less than zero!\n",
		thisprojfile,projredshift[thisprojfile]);
	exit(-1);
      }

      // very basic error checking on input grid size
      if(InputProjectionGridSize[thisprojfile] < 0){
	fprintf(stderr,"error in ReadInUserData:  InputProjectionGridSize[%d] = %lf is less than zero!\n",
		thisprojfile,InputProjectionGridSize[thisprojfile]);
	exit(-1);
      }

      // print out data for debug
      if(debug){
	printf("redshift: %lf   grid ratio: %lf   grid size: %d    projection: %d\n",
	       projredshift[thisprojfile], gridratio[thisprojfile],
	       InputProjectionGridSize[thisprojfile],thisprojfile);
	fflush(stdout);
      }
  }

  fclose(fptr);


  // open up file containing list of rebinned, corrected tiles
  fptr = fopen(RebinFileNameList,"r");

  // loop through total number of rebinned, corrected tiles and read them in
  for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++){
    RebinnedProjectionFileNames[thisprojfile] = new char[MAX_LINE_LENGTH];      

    fgets(line, MAX_LINE_LENGTH, fptr);

    sscanf(line,"%s",RebinnedProjectionFileNames[thisprojfile]);
  }

  // if user wants debug information, print out the projection file, etc.
  if(verbosedebug){
    for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++)
      printf("%d -- %s\n",thisprojfile, RebinnedProjectionFileNames[thisprojfile]);

    fflush(stdout);
  }

  fclose(fptr);


#ifdef USE_CLUSTER_SPLIT

  // open up file containing list of rebinned, corrected "No cluster" tiles
  fptr = fopen(RebinNoClustersFileNameList,"r");

  // loop through total number of rebinned, corrected tiles and read them in
  for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++){
    RebinnedProjectionNoClustersFileNames[thisprojfile] = new char[MAX_LINE_LENGTH];      

    fgets(line, MAX_LINE_LENGTH, fptr);

    sscanf(line,"%s",RebinnedProjectionNoClustersFileNames[thisprojfile]);
  }

  // if user wants debug information, print out the projection file, etc.
  if(verbosedebug){
    for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++)
      printf("NC: %d -- %s\n",thisprojfile, RebinnedProjectionNoClustersFileNames[thisprojfile]);

    fflush(stdout);
  }

  fclose(fptr);


  // open up file containing list of rebinned, corrected "Just Clusters" tiles
  fptr = fopen(RebinJustClustersFileNameList,"r");

  // loop through total number of rebinned, corrected tiles and read them in
  for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++){
    RebinnedProjectionJustClustersFileNames[thisprojfile] = new char[MAX_LINE_LENGTH];      

    fgets(line, MAX_LINE_LENGTH, fptr);

    sscanf(line,"%s",RebinnedProjectionJustClustersFileNames[thisprojfile]);
  }

  // if user wants debug information, print out the projection file, etc.
  if(verbosedebug){
    for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++)
      printf("JC: %d -- %s\n",thisprojfile, RebinnedProjectionJustClustersFileNames[thisprojfile]);

    fflush(stdout);
  }

  fclose(fptr);

  fptr = fopen(ClusterPositionFileNameList,"r");

  // get list of names for halo cluster file!
  for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++)
    for(int thisaxis = 0; thisaxis < 3; thisaxis++){

      ClusterPositionFileNames[thisaxis][thisprojfile] = new char[MAX_LINE_LENGTH];      

      fgets(line, MAX_LINE_LENGTH, fptr);

      sscanf(line,"%s",ClusterPositionFileNames[thisaxis][thisprojfile]);
    }

  fclose(fptr);

  if(verbosedebug){
    for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++)
      for(int thisaxis = 0; thisaxis < 3; thisaxis++){

	printf("CPFN:  %d %d -- %s\n",thisprojfile, thisaxis, ClusterPositionFileNames[thisaxis][thisprojfile]);

      }  
    fflush(stdout);
  }

#endif // USE_CLUSTER_SPLIT


  if(read_in_random_shifts == 1){

    random_axes = new int[NumberOfProjections];
    xshifts = new double[NumberOfProjections];
    yshifts = new double[NumberOfProjections];

    fptr = fopen("random_shifts.txt","r");

    for(int thisprojfile = 0; thisprojfile < NumberOfProjections; thisprojfile++){

      fgets(line, MAX_LINE_LENGTH, fptr);

      // actually read in data
      sscanf(line,"%lf %lf %d", &xshifts[thisprojfile],
	     &yshifts[thisprojfile], &random_axes[thisprojfile]);

      if(verbosedebug){ 
	printf("shifts, axis:  %lf %lf %d\n",xshifts[thisprojfile],yshifts[thisprojfile],random_axes[thisprojfile]);
	fflush(stdout);
      }

    }
 
   fclose(fptr);

  }

  return;
}



/* -------------------------------------------------------------------- *
   Goes through read-in projection array and rebins it to the correct
   size.  This is more of a pain than one might think, since the field of
   view we ask for is a given angular size and we have to possibly tile 
   the input projection array as well as arbitrarily shifting it along
   the x and y coordinates, and do bilinear interpolation to get the 
   values of whatever quantity we want.  Edge effects are important as
   well.

   We assume that the ratio of the spatial extent of the projection array 
   to that of the master array are given at the redshift in question.  This 
   ratio can be greater or less than one - if less than one, we need to
   tile, and if greater than one, we choose some subregion.

   Conventions:
   The double-precision variables "xoffset" and "yoffset" are the offsets
   of the zero-coordinate of the bottom left hand (0,0) coordinate in 
   our read-in projection.  

   The double-precision array gridratio[NumberOfProjections] contains 
   information about the relative size of each projection compared to 
   the angular scale demanded by the user.
 * -------------------------------------------------------------------- */
void RebinProjectionArray(int thisprojection, double xoffset, double yoffset){
  int i,j, rebin_index, input_index,xind_input,yind_input,xind_input_orig, yind_input_orig, index11,
    index12,index21,index22;
  double dtheta_master, dtheta_readin, theta_xout, theta_yout, x1,y1,x2,y2,x,y;

  if(debug){
    printf("RebinProjectionArray: entered.\n");
    fflush(stdout);
  }

  // angular resolution of each grid cell in the master ("rebinned") array.
  dtheta_master = AngularProjSize / double(OutputProjectionGridSize);

  // angular resolution of each grid cell in the input array
  dtheta_readin = AngularProjSize*gridratio[thisprojection]/double(InputProjectionGridSize[thisprojection]);

  if(debug){
    printf("RebinProjectionArray: dtheta master, readin:  %e %e %e.\n",
	   dtheta_master,dtheta_readin, gridratio[thisprojection]);
    fflush(stdout);
  }

  for(i=0; i<OutputProjectionGridSize; i++)
    for(j=0; j<OutputProjectionGridSize; j++){

      // cell index in output array
      rebin_index = j*OutputProjectionGridSize+i;
    
      // calculate position of center of output array cell in our funky angular spatial coordinate
      theta_xout = (double(i)+0.5)*dtheta_master;
      theta_yout = (double(j)+0.5)*dtheta_master;

      // what index number is this in the input array?
      // note that this does not have to be in the range 0-InputProjectionGridSize-1,
      // it can be wildly different since we're also correcting for the x and y offsets
      // and for the grid ratios.
      xind_input_orig = xind_input = int((theta_xout + xoffset*AngularProjSize)/dtheta_readin);
      yind_input_orig = yind_input = int((theta_yout + yoffset*AngularProjSize)/dtheta_readin);

      // loop until the x input index is ok
      while(xind_input < 0 || xind_input >= InputProjectionGridSize[thisprojection]){
	if(xind_input < 0) xind_input += InputProjectionGridSize[thisprojection];
	if(xind_input >= InputProjectionGridSize[thisprojection]) 
	  xind_input -= InputProjectionGridSize[thisprojection];
      }

      // loop until the y input index is ok
      while(yind_input < 0 || yind_input >= InputProjectionGridSize[thisprojection]){
	if(yind_input < 0) yind_input += InputProjectionGridSize[thisprojection];
	if(yind_input >= InputProjectionGridSize[thisprojection]) 
	  yind_input -= InputProjectionGridSize[thisprojection];
      }

      // now interpolate onto rebinned grid by bilinear interpolation.  There are multiple 
      //     cases:
      // 1)  input indices != 0 or != InputProjectionGridSize-1: simple center values,
      //     go to town
      // 2)  x or y index = 0 or InputProjectionGridSize-1:  edge case, if it's zero we're
      //     ok and if it's IPGS-1, need to be periodic.

      // along right edge?  If so, need to wrap around.
      if(xind_input==InputProjectionGridSize[thisprojection]-1 || 
	 yind_input==InputProjectionGridSize[thisprojection]-1){

	// along right side (x) but not at the top
	if((xind_input==InputProjectionGridSize[thisprojection]-1) && 
	   (yind_input < InputProjectionGridSize[thisprojection] - 1)){

	  // calculate indices!
	  // x,y - standard
	  index11 = yind_input*InputProjectionGridSize[thisprojection] + xind_input;

	  // x, y+dy - standard
	  index12 = (yind_input+1)*InputProjectionGridSize[thisprojection] + xind_input;

	  // x+dx, y
	  index21 = yind_input*InputProjectionGridSize[thisprojection] + 0;

	  // x+dx, y+dy
	  index22 = (yind_input+1)*InputProjectionGridSize[thisprojection] + 0;

	  // along the top but not at the right side
	} else if((xind_input < InputProjectionGridSize[thisprojection]-1) && 
		  (yind_input==InputProjectionGridSize[thisprojection]-1)){

	  // calculate indices!
	  // x,y - standard
	  index11 = yind_input*InputProjectionGridSize[thisprojection] + xind_input;

	  // x, y+dy
	  index12 = 0*InputProjectionGridSize[thisprojection] + xind_input;

	  // x+dx, y - standard
	  index21 = yind_input*InputProjectionGridSize[thisprojection] + (xind_input+1);

	  // x+dx, y+dy
	  index22 = 0*InputProjectionGridSize[thisprojection] + (xind_input+1);

	  // top right hand corner
	} else {

	  // calculate indices!
	  // x,y
	  index11 = yind_input*InputProjectionGridSize[thisprojection] + xind_input;

	  // x, y+dy
	  index12 = 0*InputProjectionGridSize[thisprojection] + xind_input;

	  // x+dx, y
	  index21 = yind_input*InputProjectionGridSize[thisprojection] + 0;

	  // x+dx, y+dy
	  index22 = 0*InputProjectionGridSize[thisprojection] + 0;

	}

	// otherwise we're in the center or along the left edge, and 
	// we're cool.
      } else {
	  
	//fprintf(stderr,"we're cool!\n");
	// calculate indices!

	// x,y
	index11 = yind_input*InputProjectionGridSize[thisprojection] + xind_input;

	// x, y+dy
	index12 = (yind_input+1)*InputProjectionGridSize[thisprojection] + xind_input;

	// x+dx, y
	index21 = yind_input*InputProjectionGridSize[thisprojection] + (xind_input+1);

	// x+dx, y+dy
	index22 = (yind_input+1)*InputProjectionGridSize[thisprojection] + (xind_input+1);

      }

      // no matter what, our assumed spatial locations are always the same
      x1 = (double(xind_input_orig)+0.0)*dtheta_readin;
      x2 = (double(xind_input_orig)+1.0)*dtheta_readin;

      y1 = (double(yind_input_orig)+0.0)*dtheta_readin;
      y2 = (double(yind_input_orig)+1.0)*dtheta_readin;

      // reset these to x,y for convenience
      x = theta_xout+xoffset*AngularProjSize;
      y = theta_yout+yoffset*AngularProjSize;


      // error check to make sure our values of x aren't insane
      if( x  < x1 || x > x2 || y < y1 || y > y2 ){ 
	printf("ah shit!  %e %e %e %e %e %e\n",x,x1,x2,y,y1,y2);
	printf("%e %e %e %e\n",x-x1,x2-x,y-y1,y2-y);
	printf("%d %d %d %d (%d, %d) -- %lf %lf %lf %lf -- %lf %lf --\n",index11,index21,index12,index22,
	       InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection],
	       OutputProjectionGridSize*OutputProjectionGridSize,
	       x1,x2,y1,y2,x,y
	       );

	fflush(stdout);
	exit(-1);
      }

      // verbose output:  this will really slow the code down!
      if(verbosedebug){
	printf("%d %d %d %d (%d, %d) -- %lf %lf %lf %lf -- %lf %lf --\n",index11,index21,index12,index22,
	       InputProjectionGridSize[thisprojection]*InputProjectionGridSize[thisprojection],
	       OutputProjectionGridSize*OutputProjectionGridSize,
	       x1,x2,y1,y2,x,y
	       );
	fflush(stdout);
      }
      
      // do bilinear interpolation into rebinned array, using cell-centered
      // position of the input array and vertex edges of the output array.  This
      // implicitly shifts everything by half a cell, but that's ok because
      // we're randomly shifting everything anyway, and it's always done exactly
      // the same way.
      // Note that RebinArray is always float, the positions are always double, and
      // InputProjectionArray can be float or double (as named).  So we need to 
      // explicitly case all of the input floats to doubles and then recast it all
      // to float in the end.  Hopefully this won't be too much of a performance hit.
      
      if(InputDatasetPrecision == 32){ // 32 bit - explicitly case input projeciton to double

	RebinArray[rebin_index] =  float((
					  (x2-x)*(y2-y)*double(InputProjectionArrayFloat[index11]) + 
					  (x-x1)*(y2-y)*double(InputProjectionArrayFloat[index21]) + 
					  (x2-x)*(y-y1)*double(InputProjectionArrayFloat[index12]) + 
					  (x-x1)*(y-y1)*double(InputProjectionArrayFloat[index22])
					  ) / 
					 ( (y2-y1)*(x2-x1) ) );

#ifdef USE_CLUSTER_SPLIT
	// rebin "no cluster" and "cluster only" data as well if this option is enabled
	RebinArrayNoClusters[rebin_index] =  float((
						    (x2-x)*(y2-y)*double(InputProjectionArrayNoClustersFloat[index11]) + 
						    (x-x1)*(y2-y)*double(InputProjectionArrayNoClustersFloat[index21]) + 
						    (x2-x)*(y-y1)*double(InputProjectionArrayNoClustersFloat[index12]) + 
						    (x-x1)*(y-y1)*double(InputProjectionArrayNoClustersFloat[index22])
						    ) / 
						   ( (y2-y1)*(x2-x1) ) );

	RebinArrayJustClusters[rebin_index] =  float((
						      (x2-x)*(y2-y)*double(InputProjectionArrayJustClustersFloat[index11]) + 
						      (x-x1)*(y2-y)*double(InputProjectionArrayJustClustersFloat[index21]) + 
						      (x2-x)*(y-y1)*double(InputProjectionArrayJustClustersFloat[index12]) + 
						      (x-x1)*(y-y1)*double(InputProjectionArrayJustClustersFloat[index22])
						      ) / 
						     ( (y2-y1)*(x2-x1) ) );
#endif // USE_CLUSTER_SPLIT

      } else if (InputDatasetPrecision == 64){  // 64 bit inputs - only recast at the end.

	RebinArray[rebin_index] =  float((
					  (x2-x)*(y2-y)*InputProjectionArrayDouble[index11] + 
					  (x-x1)*(y2-y)*InputProjectionArrayDouble[index21] + 
					  (x2-x)*(y-y1)*InputProjectionArrayDouble[index12] + 
					  (x-x1)*(y-y1)*InputProjectionArrayDouble[index22]
					  ) / 
					 ( (y2-y1)*(x2-x1) ) );

#ifdef USE_CLUSTER_SPLIT
	// rebin "no cluster" and "cluster only" data as well if this option is enabled
	RebinArrayNoClusters[rebin_index] =  float((
						    (x2-x)*(y2-y)*InputProjectionArrayNoClustersDouble[index11] + 
						    (x-x1)*(y2-y)*InputProjectionArrayNoClustersDouble[index21] + 
						    (x2-x)*(y-y1)*InputProjectionArrayNoClustersDouble[index12] + 
						    (x-x1)*(y-y1)*InputProjectionArrayNoClustersDouble[index22]
						    ) / 
						   ( (y2-y1)*(x2-x1) ) );

	RebinArrayJustClusters[rebin_index] =  float((
						      (x2-x)*(y2-y)*InputProjectionArrayJustClustersDouble[index11] + 
						      (x-x1)*(y2-y)*InputProjectionArrayJustClustersDouble[index21] + 
						      (x2-x)*(y-y1)*InputProjectionArrayJustClustersDouble[index12] + 
						      (x-x1)*(y-y1)*InputProjectionArrayJustClustersDouble[index22]
						      ) / 
						     ( (y2-y1)*(x2-x1) ) );
#endif // USE_CLUSTER_SPLIT

      } else {
	
	fprintf(stderr,"Problem with InputDatasetPrecision:  %d expected 32 or 64\n",InputDatasetPrecision);
	exit(-123);

      }

    } // for(j=0; ...

  if(debug){
    printf("RebinProjectionArray: exiting.\n");
    fflush(stdout);
  }

  return;
}


/* -------------------------------------------------------------------- *
   All this routine does is apply a simple correction to the rebinned
   array to take into account redshift/distance effects.  For the most 
   part we're worried about things that are redshift independent (like 
   the SZ effect) or that are very simple and are really just tests of
   the algorithm (like baryon and dark matter density).  The only one
   that really matters is the X-ray luminosity, which should have a 
   correction of 4*pi*(1+z)^4.
 * -------------------------------------------------------------------- */
void ApplyCorrectionsToRebinnedArray(int thisprojection){

  if(debug){
    printf("ApplyCorrectionToRebinnedArray: entering.\n");
    fflush(stdout);
  }

  float correction_factor;

  switch(ProjectionType){ 

  case 1:  // 1: SZ y effect:  do nothing!

    if(debug){
      printf("ApplyCorrectionToRebinnedArray: SZ Y effect.\n");
      printf("       Correction:  None whatsoever!\n");
      fflush(stdout);
    }

    break;

  case 2:  // 2:  SZ kinetic effect: do nothing!

    if(debug){
      printf("ApplyCorrectionToRebinnedArray: SZ kinetic effect.\n");
      printf("       Correction:  None whatsoever!\n");
      fflush(stdout);
    }

    break;

  case 3:  // 3: baryon density: do nothing!

    if(debug){
      printf("ApplyCorrectionToRebinnedArray: baryon density.\n");
      printf("       Correction:  None whatsoever!\n");
      fflush(stdout);
    }

    break;

  case 4:  // 4: dark matter density:  do nothing!

    if(debug){
      printf("ApplyCorrectionToRebinnedArray: dark matter density.\n");
      printf("       Correction:  None whatsoever!\n");
      fflush(stdout);
    }

    break;

  case 5:  // 5: X-ray luminosity:  divide by 4*pi*(1+z)^4, where this z
    //    is the redshift at which the slice is taken (user input)

    correction_factor = (1.0 / 4.0 / 3.1415926 / pow( (1.0+float(projredshift[thisprojection])), 4.0 ) );

    if(debug){
      printf("ApplyCorrectionToRebinnedArray: X-ray luminosity.\n");
      printf("       Correction:  divide by 4*pi*(1+z)^4!\n");
      printf("       correction:  %e   redshift:  %lf\n",correction_factor, projredshift[thisprojection]);
      fflush(stdout);
    }


    // actually apply correction!
    for(int i=0; i<OutputProjectionGridSize; i++)
      for(int j=0; j<OutputProjectionGridSize; j++){
	RebinArray[j*OutputProjectionGridSize+i] *= correction_factor;

#ifdef USE_CLUSTER_SPLIT
	RebinArrayNoClusters[j*OutputProjectionGridSize+i] *= correction_factor;
	RebinArrayJustClusters[j*OutputProjectionGridSize+i] *= correction_factor;
#endif // USE_CLUSTER_SPLIT

      }  // end of double for loop

    break;
  }

  if(debug){
    printf("ApplyCorrectionToRebinnedArray: exiting.\n");
    fflush(stdout);
  }

  return;
}

/* -------------------------------------------------------------------- *
   This routine writes out the rebinned, corrected projection arrays, 
   with file names given by the user and a fixed dataset name.  Write
   out "no cluster" and "clusters only" stuff at the end as well.
 * -------------------------------------------------------------------- */
void WriteOutRebinnedProjectionArray(int thisprojection){

  if(debug){ printf("WriteOutRebinnedProjectionArray:  entering.\n"); fflush(stdout); }

  // hdf 5 declarations - only one set necessary
  hid_t outputfile_id;
  hsize_t outputdims[2];
  herr_t h5_status;
  herr_t h5_error = -1;

  hid_t output_dset_id, output_dsp_id;

  // declare dimensions of output array
  outputdims[0] = outputdims[1] = OutputProjectionGridSize;

  /* ------------------- master output file (per projection output) -------------- */
  // create output file
  outputfile_id = H5Fcreate(RebinnedProjectionFileNames[thisprojection], H5F_ACC_TRUNC, 
			    H5P_DEFAULT, H5P_DEFAULT);
  if(debug){ 
    printf("Output file is named %s (tile %d)\n", 
	   RebinnedProjectionFileNames[thisprojection], thisprojection);
    printf("status for outputfile_id is %d\n",(int) outputfile_id);
    fflush(stdout);
  }

  assert( outputfile_id != h5_error );

  // -------------- create dataspace
  output_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(debug){ 
    printf("status for output_dsp_id is %d\n",(int) output_dsp_id);
    fflush(stdout);
  }

  assert( output_dsp_id != h5_error);

  // -------------- create dataset
  output_dset_id=H5Dcreate(outputfile_id, "/RebinnedProjection",
			   HDF5_FILE_I4, output_dsp_id, H5P_DEFAULT);
  if(debug){ 
    printf("read status for output_dset_id is: %d\n",(int) output_dset_id);
    fflush(stdout);
  }

  assert(output_dset_id != h5_error);

   // -------------- write dataset
  h5_status = H5Dwrite(output_dset_id, HDF5_FILE_I4, H5S_ALL, H5S_ALL, H5P_DEFAULT, RebinArray);
  if(debug){
    printf("status for H5Dwrite is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);

  // -------------- close dataset 
  h5_status = H5Dclose(output_dset_id);
  if(debug){ 
    printf("status for H5Dclose is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);

  // -------------- close file
  h5_status = H5Fclose(outputfile_id);
  if(debug){
    printf("status for H5Fclose is %d\n",(int) h5_status);
    fflush(stdout);
  }
  assert( h5_status != h5_error );


#ifdef USE_CLUSTER_SPLIT
  if(debug){ printf("writing out no clusters and clusters-only stuff!\n"); fflush(stdout); }
  
  /* ------------------- "no clusters" output file (per projection output) -------------- */
  // create output file
  outputfile_id = H5Fcreate(RebinnedProjectionNoClustersFileNames[thisprojection], H5F_ACC_TRUNC, 
			    H5P_DEFAULT, H5P_DEFAULT);
  if(debug){ 
    printf("Output file is named %s (tile %d)\n", 
	   RebinnedProjectionNoClustersFileNames[thisprojection], thisprojection);
    printf("status for outputfile_id is %d\n",(int) outputfile_id);
    fflush(stdout);
  }

  assert( outputfile_id != h5_error );

  // -------------- create dataspace
  output_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(debug){ 
    printf("status for output_dsp_id is %d\n",(int) output_dsp_id);
    fflush(stdout);
  }

  assert( output_dsp_id != h5_error);

  // -------------- create dataset
  output_dset_id=H5Dcreate(outputfile_id, "/RebinnedProjection",
			   HDF5_FILE_I4, output_dsp_id, H5P_DEFAULT);
  if(debug){ 
    printf("read status for output_dset_id is: %d\n",(int) output_dset_id);
    fflush(stdout);
  }

  assert(output_dset_id != h5_error);

   // -------------- write dataset
  h5_status = H5Dwrite(output_dset_id, HDF5_FILE_I4, H5S_ALL, H5S_ALL, H5P_DEFAULT, RebinArrayNoClusters);
  if(debug){
    printf("status for H5Dwrite is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);

  // -------------- close dataset 
  h5_status = H5Dclose(output_dset_id);
  if(debug){ 
    printf("status for H5Dclose is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);

  // -------------- close file
  h5_status = H5Fclose(outputfile_id);
  if(debug){
    printf("status for H5Fclose is %d\n",(int) h5_status);
    fflush(stdout);
  }
  assert( h5_status != h5_error );



  /* ------------------- "just clusters" output file (per projection output) -------------- */
  // create output file
  outputfile_id = H5Fcreate(RebinnedProjectionJustClustersFileNames[thisprojection], H5F_ACC_TRUNC, 
			    H5P_DEFAULT, H5P_DEFAULT);
  if(debug){ 
    printf("Output file is named %s (tile %d)\n", 
	   RebinnedProjectionJustClustersFileNames[thisprojection], thisprojection);
    printf("status for outputfile_id is %d\n",(int) outputfile_id);
    fflush(stdout);
  }

  assert( outputfile_id != h5_error );

  // -------------- create dataspace
  output_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(debug){ 
    printf("status for output_dsp_id is %d\n",(int) output_dsp_id);
    fflush(stdout);
  }

  assert( output_dsp_id != h5_error);

  // -------------- create dataset
  output_dset_id=H5Dcreate(outputfile_id, "/RebinnedProjection",
			   HDF5_FILE_I4, output_dsp_id, H5P_DEFAULT);
  if(debug){ 
    printf("read status for output_dset_id is: %d\n",(int) output_dset_id);
    fflush(stdout);
  }

  assert(output_dset_id != h5_error);

   // -------------- write dataset
  h5_status = H5Dwrite(output_dset_id, HDF5_FILE_I4, H5S_ALL, H5S_ALL, H5P_DEFAULT, RebinArrayJustClusters);
  if(debug){
    printf("status for H5Dwrite is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);

  // -------------- close dataset 
  h5_status = H5Dclose(output_dset_id);
  if(debug){ 
    printf("status for H5Dclose is %d\n",(int) h5_status);
    fflush(stdout);
  }

  assert( h5_status != h5_error);

  // -------------- close file
  h5_status = H5Fclose(outputfile_id);
  if(debug){
    printf("status for H5Fclose is %d\n",(int) h5_status);
    fflush(stdout);
  }
  assert( h5_status != h5_error );

#endif // USE_CLUSTER_SPLIT

  if(debug){ printf("WriteOutRebinnedProjectionArray:  leaving.\n"); fflush(stdout); }

  return;
}


/* -------------------------------------------------------------------- *
   Reads the 20-column cluster information file which contains cluster
   information such as radius, mass, position, redshift, and y-decrement
   for various times.  This is only called if the USE_CLUSTER_SPLIT stuff
   is actually on.
 * -------------------------------------------------------------------- */
void ReadClusterFile(int axis, int boxnumber){

  if(debug){ 
    printf("Entering ReadClusterFile: %d halos, reading %s\n",
	   NumberOfClusters, ClusterPositionFileNames[axis][boxnumber]); 
    fflush(stdout); 
  }

  double fjunk;
  FILE *fptr;
  char *line = new char[MAX_LINE_LENGTH];

  fptr = fopen(ClusterPositionFileNames[axis][boxnumber],"r");

  fgets(line, MAX_LINE_LENGTH, fptr);  // read junk line

  // read once to get total number of clusters
  while(fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk,
	       &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk, &fjunk) == 20)
    NumberOfClusters++;

  fclose(fptr);

  if(debug){ printf("Read in file:  There are %d clusters\n",NumberOfClusters); fflush(stdout); };

  // declare arrays (all double)
  haloxpos = new double[NumberOfClusters];
  haloypos = new double[NumberOfClusters];
  halozpos = new double[NumberOfClusters];
  haloRvir = new double[NumberOfClusters];
  haloredshift = new double[NumberOfClusters];

  if(HALO_YDEC_BINS != 10){
    fprintf(stderr,"we're assuming HALO_YDEC_BINS is always 10!\n");
    exit(-661);
  }

  for(int i=0; i < HALO_YDEC_BINS; i++)
    haloYdec[i] = new double[NumberOfClusters];

  int numclust = 0;

  // read a second time to get all values
  fptr = fopen(ClusterPositionFileNames[axis][boxnumber],"r");

  fgets(line, MAX_LINE_LENGTH, fptr);  // read junk line

  /* columns are (20 cols total, dear god):
     1) current redshift
     2) xpos
     3) ypos
     4) zpos 
     5) rvir
     6) mvir
     7) mvir_gas
     8) mvir_dm
     9) mvir_star
     10) tvir_predicted
     11-20) SZ y-decrement
  */

  while(fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &haloredshift[numclust], &haloxpos[numclust], &haloypos[numclust], &halozpos[numclust],
	       &haloRvir[numclust], &fjunk, &fjunk, &fjunk, &fjunk, &fjunk,
	       &haloYdec[0][numclust],&haloYdec[1][numclust],&haloYdec[2][numclust],&haloYdec[3][numclust],
	       &haloYdec[4][numclust],&haloYdec[5][numclust],&haloYdec[6][numclust],&haloYdec[7][numclust],
	       &haloYdec[8][numclust],&haloYdec[9][numclust])==20)
    numclust++;
  fclose(fptr);

  // renormalize radii to be in units of box size
  double BoxSizeProper;

  // normalize virial radius so that it's in units of the box length 
  // (0.0-1.0 inclusive)
  for(int i=0; i<NumberOfClusters; i++){
    BoxSizeProper = ComovingBoxSize / Hubble0 / (1.0 + haloredshift[i]);
    haloRvir[i] /= BoxSizeProper;
  }

  if(verbosedebug){
    printf("there are %d halos!\n",NumberOfClusters);
    printf("#redshift xpos ypos zpos rvir ydec_0 ydec_1 ydec_8 ydec_9\n");

    for(int i=0; i<NumberOfClusters; i++)
      printf("%lf %lf %lf %lf %lf %e %e %e %e\n",haloredshift[i],haloxpos[i],haloypos[i],halozpos[i],haloRvir[i],
	     haloYdec[0][i],haloYdec[1][i],haloYdec[8][i],haloYdec[9][i]);
    fflush(stdout);
  }


  if(debug){ printf("Exiting ReadClusterArrays\n"); fflush(stdout); }
  return;
}


/* -------------------------------------------------------------------- *
   Delete arrays which are declared in ReadClusterFile and set them
   to null.
 * -------------------------------------------------------------------- */
void DeleteClusterArrays(void){
  int i;

  if(debug){ printf("entering DeleteClusterArrays\n"); fflush(stdout); }

  delete [] haloxpos;
  delete [] haloypos;
  delete [] halozpos;
  delete [] haloRvir;
  delete [] haloredshift;

  haloxpos = NULL;
  haloypos = NULL;
  halozpos = NULL;
  haloRvir = NULL;
  haloredshift = NULL;
  
  for(i=0; i<HALO_YDEC_BINS; i++){
    delete [] haloYdec[i];
    haloYdec[i] = NULL;
  }

  if(debug){ printf("Exiting DeleteClusterArrays\n"); fflush(stdout); }

  return;
}

/* -------------------------------------------------------------------- *
   Given the list of cluster positions and various other quantities 
   (most usefully the y-decrement), run through list of halos and subtract
   the y-dec from the "No clusters" array, then add it to the "clusters 
   only" array.  This is only called if USE_CLUSTER_SPLIT is enabled.
 * -------------------------------------------------------------------- */
void ModifyArraysWithClusterInfo(int axis, int boxnumber){
  if(debug){ printf("Entering ModifyArraysWithClusterInfo %d %d %d\n",axis,boxnumber,NumberOfClusters); fflush(stdout); }

  int xpix, ypix, Rvirpix, thishalo, i, j, index;
  double xposonmap, yposonmap, cellxpos, cellypos, radiusinrvir, cellsizeinrvir, ydec;

  double getydec(double radinrvir, double cellinrvir, int halonum);

  // loop over each of the halos and do all sorts of clever stuff using it
  for(thishalo=0; thishalo < NumberOfClusters; thishalo++){

    // get position in (x,y) coords of the final image (which are between 0 and 1, inclusive)
    if(axis==0){
      xposonmap = haloypos[thishalo];
      yposonmap = halozpos[thishalo];
    }

    if(axis==1){
      xposonmap = haloxpos[thishalo];
      yposonmap = halozpos[thishalo];
    }

    if(axis==2){
      xposonmap = haloxpos[thishalo];
      yposonmap = haloypos[thishalo];
    }

    // calculate the x and y pixels (on the pixel map) of the center of the halo in question
    xpix = int( xposonmap * double(InputProjectionGridSize[boxnumber]) );
    if(xpix < 0) xpix = 0;
    if(xpix >= InputProjectionGridSize[boxnumber]) xpix = InputProjectionGridSize[boxnumber]-1;

    ypix = int( yposonmap * double(InputProjectionGridSize[boxnumber]) );
    if(ypix < 0) ypix = 0;
    if(ypix >= InputProjectionGridSize[boxnumber]) ypix = InputProjectionGridSize[boxnumber]-1;

    // calculate virial radius in pixels
    Rvirpix = int(haloRvir[thishalo]*double(InputProjectionGridSize[boxnumber]) );

    // cell size in units of the virial radius of the halo
    cellsizeinrvir = 1.0 / double(InputProjectionGridSize[boxnumber]) / haloRvir[thishalo];

    // this thing must be at least one cell on a side!
    if(Rvirpix < 1) Rvirpix = 1;

    int xlow, xhigh, ylow, yhigh;

    // set high and low bounds for the pixel range we're iterating over.
    xlow = xpix-Rvirpix-1;
    xhigh = xpix+Rvirpix+1;
    ylow = ypix-Rvirpix-1;
    yhigh = ypix+Rvirpix+1;

    // error check on bounds - can't be out of range
    if(xlow < 0) xlow = 0;
    if(xhigh >= InputProjectionGridSize[boxnumber]) xhigh = InputProjectionGridSize[boxnumber]-1;

    if(ylow < 0) ylow = 0;
    if(yhigh >= InputProjectionGridSize[boxnumber]) yhigh = InputProjectionGridSize[boxnumber]-1;

    // loop over all possible pixels which may be included in the y decrement from a given halo
    for(i=xlow; i<= xhigh; i++)
      for(j=ylow; j<=yhigh; j++){

	// get x,y position of center of cell
	cellxpos = (double(i) + 0.5) / double(InputProjectionGridSize[boxnumber]);
	cellypos = (double(j) + 0.5) / double(InputProjectionGridSize[boxnumber]);

	// distance between cell center and halo center in units of virial radius
	radiusinrvir = sqrt( pow( (cellxpos-xposonmap), 2.0) + pow( (cellypos-yposonmap), 2.0) ) / haloRvir[thishalo];

	// only consider a cell if it is within the virial radius!
	// Then get y decrement for that cell (averaged in a vaguely clever way)
	// and remove that from the no cluster image/add to just cluster image.
	if(radiusinrvir < 1.0){

	  // call snazzy function to get y decrement
	  ydec = getydec(radiusinrvir, cellsizeinrvir, thishalo);

	  index = j*InputProjectionGridSize[boxnumber] + i;

	  // subtract y decrement from "no clusters" array, and then
	  // add that to the "just clusters" array.
	  if(InputDatasetPrecision==32){

	    InputProjectionArrayNoClustersFloat[index] -= float(ydec);
	    InputProjectionArrayJustClustersFloat[index] += float(ydec);

	    if(InputProjectionArrayNoClustersFloat[index]<0.0) InputProjectionArrayNoClustersFloat[index]=0.0;
	    
	  } else { // must be 64 bit

	    InputProjectionArrayNoClustersDouble[index] -= ydec;
	    InputProjectionArrayJustClustersDouble[index] += ydec;

	    if(InputProjectionArrayNoClustersDouble[index]<0.0) InputProjectionArrayNoClustersDouble[index]=0.0;
	    
	  }

	} // if(radiusinrvir < 1.0)

	//fprintf(stderr,"(4c) %d %d\n",i,j);

      } // end of double for loop

    //fprintf(stderr,"(5)\n");

  }  //   for(thishalo=0; thishalo < NumberOfClusters; thishalo++)

  if(debug){ printf("Exiting ModifyArraysWithClusterInfo\n"); fflush(stdout); }
  return;
}


/* -------------------------------------------------------------------- *
   This function calculates the surface-area averaged y decrement of 
   a given cell.  Takes in the cell size, the distance from halo center,
   and the halo number, and uses the global array for the halo y decrement
   info to calculate average.

   inputs:
   radinrvir = radius in units of the virial radius
   cellinrvir = cell size in units of the virial radius
   halonum = this halo number
 * -------------------------------------------------------------------- */
double getydec(double radinrvir, double cellinrvir, int halonum){
  double thisydec=0.0, volweight=0.0, thisbinweight, Rvirfrac;

  // loop over bins (fractions of Rvir)
  for(int i=0; i < HALO_YDEC_BINS; i++){

    Rvirfrac = double(i) / double(HALO_YDEC_BINS);

    // Only average over bins which conceivably could overlap this cell!  
    // note that this is slightly in error - if the cell is not on the same axis
    // (i.e. offset in both x and y) there's more overlap, but not a ton.
    if( (radinrvir-cellinrvir/2.0 <= Rvirfrac) && (Rvirfrac < radinrvir+cellinrvir/2.0) ) {

      if(i==0){ // center bin - weight is pi r^2 , where r = 0.05

	thisydec += haloYdec[i][halonum]*3.14159*0.05*0.05;
	volweight += 3.14159*0.05*0.05;

      } else {  // not center bin - weight is 2 pi R dr (dr=0.1)

	thisydec += haloYdec[i][halonum]*2.0*3.14159*(Rvirfrac)*0.1;
	volweight += 2.0*3.14159*(Rvirfrac)*0.1;

      }

    } // if(radinrvir <= Rvirfrac){

  }  // for(int i=0; i < HALO_YDEC_BINS; i++)

  // don't divide by zero so as not to get hilarious (and irritating
  // to debug) results.
  if(volweight <= 1.0e-20){
    thisydec = 0.0;
  } else {
    thisydec /= volweight;
  }

  // return surface area-weighted y-decrement
  return thisydec;
}
