/* --------------------------------------------------------------------------   

     enzo_dataset_reader.C

     Purpose:
        This is an example enzo dataset reader


     Written by:  Brian O'Shea  (bwoshea@cosmos.ucsd.edu)

     Date:  Sept. 23, 2004

     Notes:
       This is an example enzo dataset reader.  The redshift to be read
       in is specified at compile time for convenience.  All datasets are
       read in and converted to proper (not comoving!) CGS coordinates.
       A section at the end of the code is marked - this is where the user
       would put in their analysis and IO routines.

     Compile notes:
       This code requires the HDF 5 library, available from 
       http://hdf.ncsa.uiuc.edu/HDF5.  I suggest you use the 1.4.5-post2
       version of the code, as we have had compatibility issues with the 
       1.6.x versions of the code.  You will have to modify the makefile 
       to point to the correct library and also modify the defined redshift
       parameter and file names in the source (below) so that the correct
       files are accessed.

     Usage notes:
       Run this code from the directory in which the data cubes reside.
       As of now, all of the instructions on which datasets to use, etc.
       are added during compilation and controlled by the use-set 
       parameters below

    ---------------------------------------------------------------------- */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#include "extern_hdf5.h"  // hdf 5 prototypes

#define DEBUG 1  // this flag turns on and off the verbose debug options:
                 // 1 is on, 0 is off

//-------------------- USER SETS THESE ----------------------------

// DATASET FILENAMES
char  *densityfile="RedshiftOutput0007.l6.Density",
  *temperaturefile="RedshiftOutput0007.l6.Temperature",
  *xvelfile="RedshiftOutput0007.l6.x-velocity",
  *yvelfile="RedshiftOutput0007.l6.y-velocity",
  *zvelfile="RedshiftOutput0007.l6.z-velocity";

char *outputfilename ="z_1_00_L6.dat";

double redshift = 0.06;  // REDSHIFT

//----------------------------------------------------------------

int main(int argc, char *argv[]){

  // hdf5 declarations
  hid_t file_id;
  hid_t       mem_type_id;
  hid_t       mem_dsp_id, file_dsp_id;
  hid_t       dsp_id;
  hid_t       typ_id;

  hsize_t     size;
  hsize_t     dims[3];
  hsize_t     xdims[3];
  hsize_t     maxdims[3];

  herr_t      h5_status;
  herr_t      h5_error = -1;
  
  hid_t    input_dset_id, output_dset_id, dataset_id, dataspace_id;
  int ndims, gridsize, i,j,k;

  // dataset names (user does not modify)
  char  *densdataset, *tempdataset, *xveldataset, *yveldataset, *zveldataset;

  // stuff used for conversion factors, etc.
  double densityconversion, lengthconversion, timeconversion, velocityconversion,
    omegamatter,comovingboxsize,hubbleconstantnow,initialredshift,rho0;

  // pointers for datasets
  float *density=NULL, 
    *temperature=NULL, 
    *xvelocity=NULL, 
    *yvelocity=NULL, 
    *zvelocity=NULL;

  // leave these alone - dataset names and simulation parameters
  densdataset = "Density";
  tempdataset = "Temperature";
  xveldataset = "x-velocity";
  yveldataset = "y-velocity";
  zveldataset = "z-velocity";

  FILE *outfile;
  
  omegamatter = 0.3;
  comovingboxsize = 256.0;  // in Mpc/h
  hubbleconstantnow = 0.7; // in units of 100 km/s/mpc
  initialredshift = 30.0;  // this is the redshift that the simulation was initialized
  gridsize = 256;          // number of cells in top grid


  // ---------------------- get data out of density file -------------------

  if(DEBUG) fprintf(stderr,"reading Density dataset\n");

  // open file
  file_id = H5Fopen(densityfile, H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );

  // open dataset
  input_dset_id = H5Dopen(file_id, densdataset);
  assert( input_dset_id != h5_error);

  // open sataspace (to get dimensions) 
  dsp_id = H5Dget_space(input_dset_id);
  assert( dsp_id != h5_error );

  // get data type 
  typ_id = H5Dget_type(input_dset_id);
  assert( typ_id != h5_error );

  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  // from the dimensional information, calculate the size of the buffer.
  size = 1;

  if(DEBUG) fprintf(stderr,"Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(DEBUG) fprintf(stderr," Dim %d\n", (int) xdims[i]);
    }

   if(DEBUG) fprintf(stderr,"Size %d\n", (int) size);

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
  assert( file_dsp_id != h5_error );

  mem_dsp_id = H5Screate_simple(1, &size, NULL);
  assert( mem_dsp_id != h5_error );

  // allocate buffer
  density = new float[(int) size];

  mem_type_id = H5T_IEEE_F32BE;
  
  // read field into an array
  h5_status = H5Dread(input_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, density);
   if(DEBUG) fprintf(stderr,"float read status %d for input field\n", 
		    (int) h5_status);
  assert( h5_status != h5_error ); 

  // close dataset
  h5_status = H5Dclose(input_dset_id);
  assert( h5_status != h5_error );

  // close hdf5 file, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );

  // close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );


  if(DEBUG) fprintf(stderr,"%f %f %f\n",density[0],density[1],density[(int) size - 1]);




  // ---------------------- get data out of temperature file -------------------

  if(DEBUG) fprintf(stderr,"reading Temperature dataset\n");

  // open file
  file_id = H5Fopen(temperaturefile, H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );

  // open dataset
  input_dset_id = H5Dopen(file_id, tempdataset);
  assert( input_dset_id != h5_error);

  // open sataspace (to get dimensions) 
  dsp_id = H5Dget_space(input_dset_id);
  assert( dsp_id != h5_error );

  // get data type 
  typ_id = H5Dget_type(input_dset_id);
  assert( typ_id != h5_error );

  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  // from the dimensional information, calculate the size of the buffer.
  size = 1;

  if(DEBUG) fprintf(stderr,"Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(DEBUG) fprintf(stderr," Dim %d\n", (int) xdims[i]);
    }

   if(DEBUG) fprintf(stderr,"Size %d\n", (int) size);

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
  assert( file_dsp_id != h5_error );

  mem_dsp_id = H5Screate_simple(1, &size, NULL);
  assert( mem_dsp_id != h5_error );

  // allocate buffer
  temperature = new float[(int) size];

  mem_type_id = H5T_IEEE_F32BE;
  
  // read field into an array
  h5_status = H5Dread(input_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, temperature);
   if(DEBUG) fprintf(stderr,"float read status %d for input field\n", 
		    (int) h5_status);
  assert( h5_status != h5_error ); 

  // close dataset
  h5_status = H5Dclose(input_dset_id);
  assert( h5_status != h5_error );

  // close hdf5 file, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );

  // close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );


  if(DEBUG) fprintf(stderr,"%f %f %f\n",temperature[0],temperature[1],temperature[(int) size - 1]);



  // ---------------------- get data out of x-velocity file -------------------

  if(DEBUG) fprintf(stderr,"reading x-velocity dataset\n");

  // open file
  file_id = H5Fopen(xvelfile, H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );

  // open dataset
  input_dset_id = H5Dopen(file_id, xveldataset);
  assert( input_dset_id != h5_error);

  // open sataspace (to get dimensions) 
  dsp_id = H5Dget_space(input_dset_id);
  assert( dsp_id != h5_error );

  // get data type 
  typ_id = H5Dget_type(input_dset_id);
  assert( typ_id != h5_error );

  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  // from the dimensional information, calculate the size of the buffer.
  size = 1;

  if(DEBUG) fprintf(stderr,"Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(DEBUG) fprintf(stderr," Dim %d\n", (int) xdims[i]);
    }

   if(DEBUG) fprintf(stderr,"Size %d\n", (int) size);

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
  assert( file_dsp_id != h5_error );

  mem_dsp_id = H5Screate_simple(1, &size, NULL);
  assert( mem_dsp_id != h5_error );

  // allocate buffer
  xvelocity = new float[(int) size];

  mem_type_id = H5T_IEEE_F32BE;
  
  // read field into an array
  h5_status = H5Dread(input_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, xvelocity);
   if(DEBUG) fprintf(stderr,"float read status %d for input field\n", 
		    (int) h5_status);
  assert( h5_status != h5_error ); 

  // close dataset
  h5_status = H5Dclose(input_dset_id);
  assert( h5_status != h5_error );

  // close hdf5 file, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );

  // close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );


  if(DEBUG) fprintf(stderr,"%f %f %f\n",xvelocity[0],xvelocity[1],xvelocity[(int) size - 1]);



  // ---------------------- get data out of y-velocity file -------------------

  if(DEBUG) fprintf(stderr,"reading y-velocity dataset\n");

  // open file
  file_id = H5Fopen(yvelfile, H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );

  // open dataset
  input_dset_id = H5Dopen(file_id, yveldataset);
  assert( input_dset_id != h5_error);

  // open sataspace (to get dimensions) 
  dsp_id = H5Dget_space(input_dset_id);
  assert( dsp_id != h5_error );

  // get data type 
  typ_id = H5Dget_type(input_dset_id);
  assert( typ_id != h5_error );

  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  // from the dimensional information, calculate the size of the buffer.
  size = 1;

  if(DEBUG) fprintf(stderr,"Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(DEBUG) fprintf(stderr," Dim %d\n", (int) xdims[i]);
    }

   if(DEBUG) fprintf(stderr,"Size %d\n", (int) size);

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
  assert( file_dsp_id != h5_error );

  mem_dsp_id = H5Screate_simple(1, &size, NULL);
  assert( mem_dsp_id != h5_error );

  // allocate buffer
  yvelocity = new float[(int) size];

  mem_type_id = H5T_IEEE_F32BE;
  
  // read field into an array
  h5_status = H5Dread(input_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, yvelocity);
   if(DEBUG) fprintf(stderr,"float read status %d for input field\n", 
		    (int) h5_status);
  assert( h5_status != h5_error ); 

  // close dataset
  h5_status = H5Dclose(input_dset_id);
  assert( h5_status != h5_error );

  // close hdf5 file, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );

  // close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );


  if(DEBUG) fprintf(stderr,"%f %f %f\n",yvelocity[0],yvelocity[1],yvelocity[(int) size - 1]);



  // ---------------------- get data out of z-velocity file -------------------

  if(DEBUG) fprintf(stderr,"reading z-velocity dataset\n");

  // open file
  file_id = H5Fopen(zvelfile, H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );

  // open dataset
  input_dset_id = H5Dopen(file_id, zveldataset);
  assert( input_dset_id != h5_error);

  // open sataspace (to get dimensions) 
  dsp_id = H5Dget_space(input_dset_id);
  assert( dsp_id != h5_error );

  // get data type 
  typ_id = H5Dget_type(input_dset_id);
  assert( typ_id != h5_error );

  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  // from the dimensional information, calculate the size of the buffer.
  size = 1;

  if(DEBUG) fprintf(stderr,"Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(DEBUG) fprintf(stderr," Dim %d\n", (int) xdims[i]);
    }

   if(DEBUG) fprintf(stderr,"Size %d\n", (int) size);

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
  assert( file_dsp_id != h5_error );

  mem_dsp_id = H5Screate_simple(1, &size, NULL);
  assert( mem_dsp_id != h5_error );

  // allocate buffer
  zvelocity = new float[(int) size];

  mem_type_id = H5T_IEEE_F32BE;
  
  // read field into an array
  h5_status = H5Dread(input_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, zvelocity);
   if(DEBUG) fprintf(stderr,"float read status %d for input field\n", 
		    (int) h5_status);
  assert( h5_status != h5_error ); 

  // close dataset
  h5_status = H5Dclose(input_dset_id);
  assert( h5_status != h5_error );

  // close hdf5 file, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );

  // close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  if(DEBUG) fprintf(stderr,"%f %f %f\n",zvelocity[0],zvelocity[1],zvelocity[(int) size - 1]);

  if(DEBUG) fprintf(stderr,"done reading in datasets, now calculating conversion factors.\n");

  //------------------ calculate conversion factors here!
  rho0 = omegamatter * 1.8788e-29 * hubbleconstantnow * hubbleconstantnow;

  // enzo conversions:  see http://cosmos.ucsd.edu/enzo/amr_guide/output.html for full
  // descriptions of how these are set up!
  lengthconversion = comovingboxsize/hubbleconstantnow*(3.0857e+24)/(1.0+redshift);
  densityconversion = rho0 * pow( (1.0+redshift)  , 3.0);
  timeconversion = 1.0 / sqrt( 4.0 * 3.14159 * 6.67e-8 * rho0 
			       * pow( (1.0+initialredshift), 3.0) );
  velocityconversion = lengthconversion/timeconversion*(1.0+redshift)/(1.0+initialredshift);

  // print out for user's convenience
  if(DEBUG){
    fprintf(stderr,"\n\n");
    fprintf(stderr,"comoving box size:         %f  (Mpc/h)\n",comovingboxsize);
    fprintf(stderr,"omega matter:              %f\n",omegamatter);
    fprintf(stderr,"hubble constant (z=0):     %f  (units of 100 km/s/Mpc)\n",hubbleconstantnow);
    fprintf(stderr,"current redshift:          %f\n",redshift);
    fprintf(stderr,"sim. initial redshift:     %f\n",initialredshift);
    fprintf(stderr,"\n\n");
    fprintf(stderr,"conversion factors (to convert enzo to proper CGS units):\n");
    fprintf(stderr,"density:   %e\n",densityconversion);
    fprintf(stderr,"length:    %e\n",lengthconversion);
    fprintf(stderr,"time:      %e\n",timeconversion);
    fprintf(stderr,"velocity:  %e\n",velocityconversion);
    fprintf(stderr,"\n\n");
   }

  if(DEBUG) fprintf(stderr,"first, last 5 cell values:\n\n");

  // actually modify density, velocity datasets (temperature is unmodified)
  for(i=0; i < (gridsize*gridsize*gridsize); i++){
    density[i] = density[i]* ((float) densityconversion);
    xvelocity[i] = xvelocity[i] * ((float) velocityconversion);
    yvelocity[i] = yvelocity[i] * ((float) velocityconversion);
    zvelocity[i] = zvelocity[i] * ((float) velocityconversion);
 
    // if debug is turned on, print out first and last 5 values
    // (just so we can check out what's going on)
    if(DEBUG)
      if( i <= 4 || i >= gridsize*gridsize*gridsize-5)
	fprintf(stderr,"%d:   %e %e %e %e %e\n",i, density[i],temperature[i],
		xvelocity[i],yvelocity[i],zvelocity[i]);
      
  } // end of for(i=0; i < (gridsize*gridsize*gridsize); i++){

  if(DEBUG){
    fprintf(stderr,"\n\n");
    fprintf(stderr,"done calculating conversion factors and converting variables.\n");
    fprintf(stderr,"beginning analysis.\n");
    fprintf(stderr,"\n\n");
  }

  // ------------------------- USER DOES ANALYSIS HERE --------------------- 

  /* at this point there are 5 buffers:  density[], temperature[], xvelocity[],
     yvelocity[] and zvelocity[], which are all in proper (NOT comoving!) CGS units.  
     So, density is in g/cm^3, temperature is in Kelvin, and velocities are in 
     cm/s.

     All of the datasets provided are cubes of size 256^3, and the size of the 
     grid is held in the int variable gridsize.  The "l6" datasets cover a 
     spatial extent of 8 Mpc/h (comoving) on a side and the "l7" datasets cover 
     a spatial extent of 4 Mpc/h on a side, so the cells have comoving resolutions 
     of 31.25 Kpc/h and 15.625 Kpc/h on a side, respectively.  For reference, the 
     virial radius of the cluster is 2.8 Mpc at z=0.

     The datasets are in column-major ("fortran ordered") format.  The index of
     cell (i,j,k) in the large 1D array is therefore:

            cellindex = k*gridsize*gridsize + j*gridsize + i;

     That should be all of the information that is necessary to get started 
     with data analysis.  If you have problems or questions, please email me
     (Brian O'Shea) at bwoshea@cosmos.ucsd.edu.
  */


  fprintf(stderr,"about to open file %s\n",outputfilename);

  outfile = fopen(outputfilename,"w");

  fprintf(stderr,"writing to file %s, don't hold your breath\n",outputfilename);

  // actually modify density, velocity datasets (temperature is unmodified)
  for(i=0; i < (gridsize*gridsize*gridsize); i++)
	fprintf(outfile,"%d %e %e %e %e %e\n",i, density[i],temperature[i],
		xvelocity[i],yvelocity[i],zvelocity[i]);

  fclose(outfile);

  fprintf(stderr,"done writing to file %s\n",outputfilename);







  // --------------- DONE DOING ANALYSIS AND OUTPUT ------------------------ 

  // clean up arrays
  delete [] density;
  delete [] temperature;
  delete [] xvelocity;
  delete [] yvelocity;
  delete [] zvelocity;

  // exit successfully
  return 1;

}
