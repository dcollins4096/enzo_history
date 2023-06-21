/*-----------------------------------------------------------

Finds grids with highest baryon and dark matter density in 
the simulation and returns the (x,y,z) position in double
precision.  This is calculated for both dark matter AND
baryons.

Usage:

denpeak <amr file name>

BWO, June 19 2004

Note:  June 19 is not the original date of writing, but it's
the day I added 64-bit capability to the position and also
put dark matter and baryon in the same file.
-----------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#include "extern_hdf5.h"  // hdf 5 prototypes

#define MAX_LINE_LENGTH 256
#define DEBUG 0

int main(int argc, char *argv[]){

  FILE *hierfile, *headerfile;

  char *hierarchyfilename=NULL,*datafilename=NULL;

  char *line = new char[MAX_LINE_LENGTH];
  char *gfname = new char[MAX_LINE_LENGTH];

  int griddimx,griddimy,griddimz,  // grid xyz dims
    gridsex,gridsey,gridsez,  // grid start indices
    grideex,grideey,grideez,  // grid end indices
    gdx,gdy,gdz,  // grid dims - no buffers
    cellindex;

  double glex,gley,glez,  // grid left edge
    grex,grey,grez,  // grid right edge
    ccx,ccy,ccz;    // cell center

float
    fjunk,
    *densbuff,*dmdensbuff;

  float maxdens_baryon, maxdens_dm;

  double max_bar_xpos, max_bar_ypos, max_bar_zpos,
    max_dm_xpos, max_dm_ypos, max_dm_zpos;

  int i,j,k,ndims,m,numpart;
  
  // hdf 5 stuff
  hid_t       file_id;
  hid_t      dens_dset_id,dmdens_dset_id;

  hid_t       mem_dsp_id, file_dsp_id;
  hid_t       mem_type_id;
  hid_t       dsp_id;
  hid_t       typ_id;

  hsize_t     size;
  hsize_t     dims[4];

  hsize_t     xdims[4];
  hsize_t     maxdims[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  // read in command-line arguments
  datafilename=argv[1];

  maxdens_baryon = maxdens_dm = -1.0;

  headerfile = fopen(datafilename,"r");

  // get hierarchy file name
  hierarchyfilename = datafilename;
  hierarchyfilename=strcat(hierarchyfilename,".hierarchy");
  fprintf(stderr,"hierarchy file:  %s\n",hierarchyfilename);
  fprintf(stderr,"data file:  %s\n",datafilename);

  // Open up hierarchy file 
  hierfile = fopen(hierarchyfilename,"r");

  /* read through the hierarchy file. and do stuff */

  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){

    /* when the "Grid = " line is encountered, this is 
       the beginning of a grid entry.  Proceed to read 
       in the pertinent info.    */
    if(strncmp(line,"Grid = ",7)==0){


      fgets(line, MAX_LINE_LENGTH, hierfile);  // junk line
      fgets(line, MAX_LINE_LENGTH, hierfile);  // grid dim
      sscanf(line,"GridDimension     = %d %d %d",
	     &griddimx,&griddimy,&griddimz);
      fgets(line, MAX_LINE_LENGTH, hierfile);  // start index
      sscanf(line,"GridStartIndex    = %d %d %d",
	     &gridsex,&gridsey,&gridsez);
      fgets(line, MAX_LINE_LENGTH, hierfile);  // end index
      sscanf(line,"GridEndIndex      = %d %d %d",
	     &grideex,&grideey,&grideez);
      fgets(line, MAX_LINE_LENGTH, hierfile);  // left edge
      sscanf(line,"GridLeftEdge      = %lf %lf %lf",&glex,&gley,&glez);
      fgets(line, MAX_LINE_LENGTH, hierfile);  // right edge
      sscanf(line,"GridRightEdge     = %lf %lf %lf",&grex,&grey,&grez);

      // "junk" lines
      for(i=0;i<5;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
      
      // get name of grid file
      sscanf(line,"BaryonFileName = %s",gfname);

      // "junk" lines
      for(i=0;i<5;i++)
	fgets(line, MAX_LINE_LENGTH, hierfile);

      sscanf(line,"NumberOfParticles   = %d",&numpart);

      if(DEBUG) fprintf(stderr,"%d %d %d %d %d %d %d %d %d\n",
		       griddimx,griddimy,griddimz,
		       gridsex,gridsey,gridsez,
		       grideex,grideey,grideez);
      if(DEBUG) fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n",
		       glex,gley,glez,
		       grex,grey,grez);
      if(DEBUG) fprintf(stderr,"%s\n",gfname);

      gdx = 1+grideex-gridsex;
      gdy = 1+grideey-gridsey;
      gdz = 1+grideez-gridsez;

      if(DEBUG) fprintf(stderr,"Grid is %d %d %d without buffers\n",gdx,gdy,gdz);
      if(DEBUG) fprintf(stderr,"that comes out to %d total cells\n",gdx*gdy*gdz);
      
      // open grid file
      file_id = H5Fopen(gfname, H5F_ACC_RDWR, H5P_DEFAULT);
      assert( file_id != h5_error );
      
      // open density dataset
      dens_dset_id = H5Dopen(file_id,"Density");
      assert( dens_dset_id != h5_error );

      if(numpart>0){
	// open density dataset
	dmdens_dset_id = H5Dopen(file_id,"Dark_Matter_Density");
	assert( dmdens_dset_id != h5_error );
      }

      // open dataspace (to get dimensions)
      dsp_id = H5Dget_space(dens_dset_id);
      assert( dsp_id != h5_error );
	
      // get data type
      typ_id = H5Dget_type(dens_dset_id);
      assert( typ_id != h5_error );
	
      // get dimensional information from dataspace
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

      if ( H5Tequal( typ_id, H5T_NATIVE_FLOAT ) ){
	
	densbuff = new float[(int) size];

	if(numpart>0)
	  dmdensbuff = new float[(int) size];

	// read density field
	h5_status = H5Dread(dens_dset_id,  H5T_IEEE_F32BE, mem_dsp_id, 
			    file_dsp_id, H5P_DEFAULT, densbuff);
	assert( h5_status != h5_error );


	if(numpart>0){
	  // read dark matter density field
	  h5_status = H5Dread(dmdens_dset_id,  H5T_IEEE_F32BE, mem_dsp_id, 
			      file_dsp_id, H5P_DEFAULT, dmdensbuff);
	  assert( h5_status != h5_error );
	}

      }
	
      for(k=0; k < gdz; k++)
	for(j=0; j < gdy; j++)
	  for(i=0; i < gdx; i++){

	    cellindex = k*gdx*gdy + j*gdx + i;	      
	      
	    // calculate center of this grid cell (i,j,k) in double precision
	    ccx = glex + ( (grex-glex) / ((double) gdx) ) * ( ((double) i) + 0.5 );
	    ccy = gley + ( (grey-gley) / ((double) gdy) ) * ( ((double) j) + 0.5 );
	    ccz = glez + ( (grez-glez) / ((double) gdz) ) * ( ((double) k) + 0.5 );

	    if( densbuff[cellindex] > maxdens_baryon){
	      maxdens_baryon = densbuff[cellindex];
	      max_bar_xpos = ccx;
	      max_bar_ypos = ccy;
	      max_bar_zpos = ccz;
	    }

	    if(numpart > 0)
	      if( dmdensbuff[cellindex] > maxdens_dm){
		maxdens_dm = dmdensbuff[cellindex];
		max_dm_xpos = ccx;
		max_dm_ypos = ccy;
		max_dm_zpos = ccz;
	      }
	    
	  }
	      
      delete [] densbuff;

      if(numpart>0)
	delete [] dmdensbuff;

      // close everything, doing appropriate error checking
      h5_status = H5Sclose(dsp_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Tclose(typ_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Sclose(mem_dsp_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Sclose(file_dsp_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Dclose(dens_dset_id);
      assert( h5_status != h5_error );

      if(numpart>0){
	h5_status = H5Dclose(dmdens_dset_id);
	assert( h5_status != h5_error );
      }
      
      h5_status = H5Fclose(file_id);
      assert( h5_status != h5_error );
    }  
  }  // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){
  
  fclose(hierfile);  // close hierarchy file

  printf("\n\nthe max baryon density is %e and is at\n\n",maxdens_baryon);
  printf("      x,y,z =   %.12lf   %.12lf   %.12lf\n",
	 max_bar_xpos, max_bar_ypos, max_bar_zpos);

  printf("\n\nthe max dark matter density is %e and is at\n\n",maxdens_dm);
  printf("      x,y,z =   %.12lf   %.12lf   %.12lf\n\n",
	 max_dm_xpos, max_dm_ypos, max_dm_zpos);

  exit(0);
}

