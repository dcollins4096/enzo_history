#include "EnzoExtract.h"

/* --------------------- CreateExtractionArrays ---------------------
 *  Creates Extraction arrays.
 *  This is the routine that must be modified to ADD more projection
 *  arrays.
 * ------------------------------------------------------------------ */
void CreateExtractionArrays(void){
  if(DEBUG) fprintf(stderr,"CreateExtractionArrays:  entering\n");
  int i,j;
  // calculate the number of cells the extraction arrays take up in the x direction
  xextractnumcells = (int) 
    ((xend - xstart) * ((double) rootgridsize)*((double) pow((double) refineby,(double) outputmaxlevel)));//AK
  if(xextractnumcells < 1) xextractnumcells = 1;
  
  // calculate the number of cells the extraction arrays take up in the y direction
  yextractnumcells = (int) 
    ((yend - ystart) * ((double) rootgridsize)*((double) pow((double) refineby,(double) outputmaxlevel)));//AK
  if(yextractnumcells < 1) yextractnumcells = 1;
  
  // calculate the number of cells the extraction arrays take up in the z direction
  zextractnumcells = (int) 
    ((zend - zstart) * ((double) rootgridsize)*((double) pow((double) refineby,(double) outputmaxlevel)));//AK
  if(zextractnumcells < 1) zextractnumcells = 1;
  
  // extraction array is this size
  extractionarraysize = xextractnumcells * yextractnumcells * zextractnumcells ;
  
  if(DEBUG) fprintf(stderr,"CreateExtractionArrays: %d %d %d %d\n",
		    xextractnumcells, yextractnumcells, zextractnumcells,
		    extractionarraysize);

  // create arrays!

  numberofextractionfields = 0;

  ext_fieldnames = new char*[MAX_EXTRACTION_FIELDS];

#ifdef FIELD_VALUES_DOUBLE
  ext_fieldvalues = new double*[MAX_EXTRACTION_FIELDS];
#else
  ext_fieldvalues = new float*[MAX_EXTRACTION_FIELDS];
#endif



  for(i=0; i<MAX_EXTRACTION_FIELDS; i++)
    ext_fieldnames[i] = new char[MAX_NAME_LENGTH];

  /* ---------------- THIS IS WHERE NEW FIELDS SHOULD BE ADDED --------
     Make follow along in the same format.  Make sure that the field 
     name (ext_fieldnames[numberofextractionfields]) is the same as the
     field name in the grid files.

     Logic can also be added to figure out what fields are in the grid
     file - perhaps in ReadParameterFile().
     ------------------------------------------------------------------ */

  // density
  ext_fieldnames[numberofextractionfields] = "Density";
#ifdef FIELD_VALUES_DOUBLE
  ext_fieldvalues[numberofextractionfields] = new double[extractionarraysize];
#else
  ext_fieldvalues[numberofextractionfields] = new float[extractionarraysize];
#endif
  if(DEBUG) fprintf(stderr,"CreateExtractionArrays: field number %d is named %s\n",
		    numberofextractionfields,ext_fieldnames[numberofextractionfields]);
  numberofextractionfields++;



  
  // temperature
  ext_fieldnames[numberofextractionfields] = "Temperature";
#ifdef FIELD_VALUES_DOUBLE
  ext_fieldvalues[numberofextractionfields] = new double[extractionarraysize];
#else
  ext_fieldvalues[numberofextractionfields] = new float[extractionarraysize];
#endif
  if(DEBUG) fprintf(stderr,"CreateExtractionArrays: field number %d is named %s\n",
		    numberofextractionfields,ext_fieldnames[numberofextractionfields]);
  numberofextractionfields++;

  // x-velocity
  ext_fieldnames[numberofextractionfields] = "x-velocity";
#ifdef FIELD_VALUES_DOUBLE
  ext_fieldvalues[numberofextractionfields] = new double[extractionarraysize];
#else
  ext_fieldvalues[numberofextractionfields] = new float[extractionarraysize];
#endif
  if(DEBUG) fprintf(stderr,"CreateExtractionArrays: field number %d is named %s\n",
		    numberofextractionfields,ext_fieldnames[numberofextractionfields]);
  numberofextractionfields++;

  // y-velocity
  ext_fieldnames[numberofextractionfields] = "y-velocity";
#ifdef FIELD_VALUES_DOUBLE
  ext_fieldvalues[numberofextractionfields] = new double[extractionarraysize];
#else
  ext_fieldvalues[numberofextractionfields] = new float[extractionarraysize];
#endif
  if(DEBUG) fprintf(stderr,"CreateExtractionArrays: field number %d is named %s\n",
		    numberofextractionfields,ext_fieldnames[numberofextractionfields]);
  numberofextractionfields++;

  // z-velocity
  ext_fieldnames[numberofextractionfields] = "z-velocity";
#ifdef FIELD_VALUES_DOUBLE
  ext_fieldvalues[numberofextractionfields] = new double[extractionarraysize];
#else
  ext_fieldvalues[numberofextractionfields] = new float[extractionarraysize];
#endif
  if(DEBUG) fprintf(stderr,"CreateExtractionArrays: field number %d is named %s\n",
		    numberofextractionfields,ext_fieldnames[numberofextractionfields]);
  numberofextractionfields++;
 

  for(i=0; i<numberofextractionfields; i++)
    for(j=0; j < extractionarraysize; j++)
      ext_fieldvalues[i][j] = -1.0;

  if(DEBUG) fprintf(stderr,"CreateExtractionArrays:  exiting\n");

}



/* ----------------------- DeleteExtractionArrays ----------------------
 *  Deletes Extraction Arrays (woo)
 * ------------------------------------------------------------------ */
void DeleteExtractionArrays(void){
  if(DEBUG) fprintf(stderr,"DeleteExtractionArrays:  entering\n");
  int i;

  if(DEBUG) fprintf(stderr,"DeleteExtractionArrays:  deleting extraction field values\n");
  for(i=0; i < numberofextractionfields; i++)
    delete [] ext_fieldvalues[i];

  /*  TEMPORARILY COMMENTED OUT - PROBABLY NOT A BIG DEAL ANYWAY
    if(DEBUG) fprintf(stderr,"DeleteExtractionArrays:  deleting extraction field names\n");
    for(i=0; i < MAX_EXTRACTION_FIELDS; i++)
    delete [] ext_fieldnames[i];
  */

  if(DEBUG) fprintf(stderr,"DeleteExtractionArrays:  deleting arrays of pointers\n");
  delete [] ext_fieldvalues;
  //delete [] ext_fieldnames;


  if(DEBUG) fprintf(stderr,"DeleteExtractionArrays:  exiting\n");
}



/* ---------------------- SwitchArraysToRowOrder -----------------------
 *  Switches output arrays from column-major order to row-major order.
 *  C(i,j,k) for a 3D array of dimensions Nx,Ny,Nz is:
 *
 *  Row Major (i varies slowest):
 *      i*Ny*Nz + j*Nz +k
 *
 *  Column Major (k varies slowest):
 *      k*Nx*Ny + j*Nx + i
 *
 * ------------------------------------------------------------------ */
void SwitchArraysToRowOrder(void){
  if(DEBUG) fprintf(stderr,"SwitchArraysToRowOrder:  entering\n");
  int i,j,k,cell,rowindex,columnindex,field;

#ifdef FIELD_VALUES_DOUBLE
  double *tempbuffer = new double[extractionarraysize];
#else 
  float *tempbuffer = new float[extractionarraysize];
#endif

  // loop over all fields!
  for(field=0; field < numberofextractionfields; field++){

    if(DEBUG) fprintf(stderr,"SwitchArraysToRowOrder:  Switching %s field\n",
		      ext_fieldnames[field]);
    
    // make buffer values really, really negative so we can later
    // make sure that we've allocated cells correctly
    for(i=0; i< extractionarraysize; i++) tempbuffer[i] = -HUGENUMBER;
    
    // loop over entire projection volume and switch cells
    for(i=0; i < xextractnumcells; i++)
      for(j=0; j < xextractnumcells; j++)
	for(k=0; k < xextractnumcells; k++){
	  
	// calculate row-major index
	  rowindex = i * yextractnumcells * zextractnumcells 
	    + j * zextractnumcells + k;
	  
	  // calculate column-major index
	  columnindex = k * xextractnumcells * yextractnumcells
	    + j * xextractnumcells + i; 
	  
	  // copy values from ext_fieldvalues[field] to temporary buffer
	  tempbuffer[rowindex] = ext_fieldvalues[field][columnindex];

	}
    
    // now copy back to the field value
    for(i=0; i< extractionarraysize; i++){
      ext_fieldvalues[field][i]=tempbuffer[i];
      if(ext_fieldvalues[field][i]==-HUGENUMBER){
	fprintf(stderr,"WTF?  SwitchArraysToRowOrder -- field not allocated corrrectly. %d %d\n",field,i);
      } 
    }
    
  } // for(field=0; field < numberofextractionfields; field++){
  
  if(DEBUG) fprintf(stderr,"SwitchArraysToRowOrder:  exiting\n");
}



/* ---------------------- WriteExtractionArrays ------------------------
 *  Write out extraction files - this was the ugliest routine (arguably)
 *  in EnzoProj, but I figured out how to neaten it up quite a bit for
 *  EnzoExtract.  Loop over all numberofextractionfields fields, outputting
 *  each one into its very own HDF 5 file which has the same name as
 *  the dataset.
 * ------------------------------------------------------------------ */
void WriteExtractionArrays(void){
  if(DEBUG) fprintf(stderr,"WriteExtractionArrays:  entering\n");

  // hdf 5 declarations 
  hid_t extract_file_id, extract_dset_id, extract_dsp_id;
  hsize_t outputdims[3];
  herr_t h5_status, h5_error=-1;

  int field;

  // set output dimensions
  if(outputrowmajor){ // for row major format
    outputdims[0] = (hsize_t) zextractnumcells;
    outputdims[1] = (hsize_t) yextractnumcells;
    outputdims[2] = (hsize_t) xextractnumcells;
  } else {  // otherwise it's column major
    outputdims[0] = (hsize_t) xextractnumcells;
    outputdims[1] = (hsize_t) yextractnumcells;
    outputdims[2] = (hsize_t) zextractnumcells;
  }

  // loop over fields, creating a file and adding values to each one.
  for(field=0; field < numberofextractionfields; field++){
    if(DEBUG) fprintf(stderr,"WriteExtractionArrays:  writing field %s\n",ext_fieldnames[field]);

    // create file (which is the same as the field name)
    extract_file_id = H5Fcreate(ext_fieldnames[field], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert( extract_file_id != h5_error);

    // -------------- create dataspace for field
    extract_dsp_id=H5Screate_simple(3, outputdims, NULL);
    assert( extract_dsp_id != h5_error);

    // -------------- create dataset for field
#ifdef FIELD_VALUES_DOUBLE
    extract_dset_id=H5Dcreate(extract_file_id,ext_fieldnames[field], HDF5_FILE_I8, extract_dsp_id, H5P_DEFAULT);
#else
    extract_dset_id=H5Dcreate(extract_file_id,ext_fieldnames[field], HDF5_FILE_I4, extract_dsp_id, H5P_DEFAULT);
#endif
    assert(extract_dset_id != h5_error);
    
    // -------------- write dataset 
#ifdef FIELD_VALUES_DOUBLE
    h5_status = H5Dwrite(extract_dset_id, HDF5_FILE_I8,H5S_ALL, H5S_ALL, H5P_DEFAULT,ext_fieldvalues[field]);
#else
    h5_status = H5Dwrite(extract_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,ext_fieldvalues[field]);

#endif
    assert( h5_status != h5_error);
    
    // -------------- close dataset 
    h5_status = H5Dclose(extract_dset_id);
    assert( h5_status != h5_error);
    
    // -------------- close dataspace
    h5_status = H5Sclose(extract_dsp_id);
    assert( h5_status != h5_error);

    // -------------- close file 
    h5_status = H5Fclose(extract_file_id);
    assert( h5_status != h5_error );
  }

  if(DEBUG) fprintf(stderr,"WriteExtractionArrays:  exiting\n");
}
