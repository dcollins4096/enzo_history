#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#include "extern_hdf5.h"  // hdf 5 prototypes

#define MAX_LINE_LENGTH 256
#define HDF5_I4 H5T_NATIVE_INT
//#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I4 H5T_IEEE_F32BE  

#define DEBUG 1  // 1 on, 0 off
#define VERBOSEDEBUG 1 // 1 on, 0 off

// return calls
#define SUCCESS 1
#define FAILURE 0


int total_number_particles = 0, total_number_grids;

int NumberOfGrids(char *hierfilename);
int GetGridInfo(int numberofgrids,char *hierfilename);
int AddParticlesToArrays(int gridnum);

float *xpositions, *ypositions, *zpositions, *particlemasses,
  *xvelocities, *yvelocities, *zvelocities,
  *particlecrtime,*particledyntime,*particlemetalfrac;

int *griddx;

float *gridlex,*gridrex;

int *gridnump, array_position=0, particle_count=0;
char **gridfilenames;

int main(int argc, char *argv[]){

  char *inputfilename=argv[1];  // name of Enzo hierarchy
  int i;
  FILE *outputfile;



  printf("%s\n",inputfilename);
  fflush(stdout);

  // get hierarchy file name
  int inlen = (int) strlen(inputfilename);
  char *hierfile = new char[inlen+12];
  strcpy(hierfile,inputfilename);
  hierfile=strcat(hierfile,".hierarchy");

  int total_number_grids = NumberOfGrids(hierfile);

  printf("there are %i total particles\n",total_number_particles);
  fflush(stdout);

  fprintf(stderr,"**starting**\n");

  xpositions = new float[total_number_particles];
  ypositions = new float[total_number_particles];
  zpositions = new float[total_number_particles];
  xvelocities = new float[total_number_particles];
  yvelocities = new float[total_number_particles];
  zvelocities = new float[total_number_particles];
  particlemasses = new float[total_number_particles];
  particlecrtime = new float[total_number_particles];
  particledyntime = new float[total_number_particles];
  particlemetalfrac = new float[total_number_particles];

  fprintf(stderr,"**ending**\n");


  for(i=0; i<total_number_particles; i++){

    //fprintf(stderr,"%i\n",i);

    xpositions[i] = ypositions[i] = zpositions[i] = 
    xvelocities[i] = yvelocities[i] = zvelocities[i] = 
      particlemasses[i] = particlecrtime[i] = 
      particledyntime[i] = particlemetalfrac[i] = -1.0;
  }

  fprintf(stderr,"**done setting values**\n");


  // grid file name info
  gridfilenames = new char*[total_number_grids];

  // number of particles info
  gridnump = new int[total_number_grids];

  // grid info
  griddx = new int[total_number_grids];
  gridlex = new float[total_number_grids];
  gridrex = new float[total_number_grids];


  fprintf(stderr,"**(1)**\n");

  for(i=0; i<total_number_grids; i++){
    gridnump[i]=griddx[i]=-1;
    gridlex[i]=gridrex[i]=-1.0;
    gridfilenames[i] = new char[MAX_LINE_LENGTH];
   }

  GetGridInfo(total_number_grids,hierfile);

  for(i=0; i<total_number_grids; i++)
    if(gridnump[i]>0)
      AddParticlesToArrays(i);


  for(i=0; i<total_number_particles; i++)
    if(xpositions[i] < 0.0 || xpositions[i] > 1.0 ||
       ypositions[i] < 0.0 || ypositions[i] > 1.0 ||
       zpositions[i] < 0.0 || zpositions[i] > 1.0)
      fprintf(stderr,"uh-oh %f %f %f %i %i\n",xpositions[i],ypositions[i],zpositions[i],i,total_number_particles);

   
  outputfile = fopen("starinfo.dat","w");

  fprintf(stderr,"writing file starinfo.dat\n") ;

  
    for(i=0; i<total_number_particles; i++)
      if(particlecrtime[i] > 0.0)
	fprintf(outputfile,"%f %f %f %e %e %e %e %f %f %e\n",
		xpositions[i],ypositions[i],zpositions[i],
		xvelocities[i],yvelocities[i],zvelocities[i],
		particlemasses[i], particlecrtime[i], 
		particledyntime[i], particlemetalfrac[i]);
    
  fclose(outputfile);

  delete [] xpositions;
  delete [] ypositions;
  delete [] zpositions;

  delete [] xvelocities;
  delete [] yvelocities;
  delete [] zvelocities;

  delete [] particlemasses;
  delete [] particlecrtime;
  delete [] particledyntime;
  delete [] particlemetalfrac;

}


/*------------------------- NumberOfGrids(char) --------------- 
 *
 *   Reads the hierarchy file and counts the number of grids.
 *   This will be used shortly thereafter to get all of the 
 *   grid info.
 *
 *-------------------------------------------------------------*/
int NumberOfGrids(char *hierfilename){
  if(DEBUG) fprintf(stderr,"in NumberOfGrids\n");

  FILE *hierfile;
  char *line = new char[MAX_LINE_LENGTH];
  int numgrids=0,npartthisgrid;

  if(DEBUG) fprintf(stderr,"hierarchy file: %s\n",hierfilename);
    
  hierfile = fopen(hierfilename,"r");

  // read through hierarchy file, counting # of grids
  // lines that start with "Grid = " are the beginning of a new
  // piece of grid info
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){
    if(strncmp(line,"Grid = ",7)==0) numgrids++;

    if(strncmp(line,"NumberOfParticles",17)==0){
      sscanf(line,"NumberOfParticles   = %d",
	     &npartthisgrid);
      total_number_particles += npartthisgrid;
    }
  }
    
  fclose(hierfile);

  if(DEBUG) fprintf(stderr,"NumberOfGrids:  there are %i grids\n",numgrids);
  
  // clean up dynamically allocated stuff
  delete [] line;

  if(DEBUG) fprintf(stderr,"exiting NumberOfGrids\n");

  // return # of grids
  return numgrids;  
}


/*----------------- GetGridInfo(int,char) ------------------------
 *
 *   Reads through grid file and retrieves all important information -
 *   grid dimensions, bounds, level, and file name, and puts them all 
 *   into arrays which were allocated in main().
 *
 *---------------------------------------------------------------- */
int GetGridInfo(int numberofgrids,char *hierfilename){
  if(DEBUG) fprintf(stderr,"in GetGridInfo\n");
  FILE *hierfile;
  char *line = new char[MAX_LINE_LENGTH];
  char *gfname = new char[MAX_LINE_LENGTH];
  int grid=0,i;

  int griddimx,griddimy,griddimz,gridsex,gridsey,gridsez,
    grideex,grideey,grideez,nbfields, particleonly, ijunk;

  float fjunk;

  if(DEBUG) fprintf(stderr,"hierarchy file: %s\n",hierfilename);

  // open hierarchy file
  hierfile = fopen(hierfilename,"r");

  // read through hierarchy file, get grid info
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){

    // lines that start with "Grid =" indicate a new set of grid info
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
      sscanf(line,"GridLeftEdge      = %f %f %f",
	     &gridlex[grid],&fjunk,&fjunk);

      fgets(line, MAX_LINE_LENGTH, hierfile);  // right edge
      sscanf(line,"GridRightEdge     = %f %f %f",
	     &gridrex[grid],&fjunk,&fjunk);
      
      for(i=0;i<3;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
      sscanf(line,"NumberOfBaryonFields = %d",&nbfields);

      if(nbfields==0){
	particleonly=1;

	fgets(line, MAX_LINE_LENGTH, hierfile);
	sscanf(line,"NumberOfParticles   = %d",
	       &gridnump[grid]);  // number of particles on this grid

	fgets(line, MAX_LINE_LENGTH, hierfile);
	sscanf(line,"ParticleFileName = %s",gridfilenames[grid]);

	fprintf(stderr,"GetGridInfo:  No baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

      } else {
	// "junk" lines
	for(i=0;i<2;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
	
	// grid file name
	sscanf(line,"BaryonFileName = %s",gridfilenames[grid]);

	// "junk" lines
	for(i=0;i<5;i++) fgets(line, MAX_LINE_LENGTH, hierfile);

	sscanf(line,"NumberOfParticles   = %d",
	       &gridnump[grid]);  // number of particles on this grid

	fprintf(stderr,"GetGridInfo: WITH baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

      }

      // grid dims - buffers stripped
      griddx[grid] = 1+grideex-gridsex;

      grid++;  // increment grid number!
    }
  }

  fclose(hierfile);

  // clean up dynamically allocated stuff
  delete [] line;
  delete [] gfname;

  if(DEBUG) fprintf(stderr,"exiting GetGridInfo\n");

  // return code - number of grids!
  return grid;
}

int AddParticlesToArrays(int gridnum){
  if(DEBUG) fprintf(stderr,"in AddParticlesToArrays %d\n",gridnum);


  // are there no particles on this grid?  if not, exit!
  if( gridnump[gridnum] == 0) return SUCCESS;

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

  hid_t       pposx_dset_id,pposy_dset_id, pposz_dset_id,
    pxvel_dset_id, pyvel_dset_id, pzvel_dset_id,
    pmass_dset_id,pcrtime_dset_id,
    pdyntime_dset_id,pmfrac_dset_id;

  int i, ndims;
  float *particle_mass_buffer, *particle_crtime_buffer,
    *particle_xvel_buffer, *particle_yvel_buffer,
    *particle_zvel_buffer,
    *particle_dyntime_buffer, *particle_mfrac_buffer;

  double *particle_posx_buffer, *particle_posy_buffer, 
    *particle_posz_buffer;
  int numberofparticles,totalnumberofstars=0;

  float cellvolume, deltax;

  deltax = (gridrex[gridnum]-gridlex[gridnum])/((float)griddx[gridnum]);

  cellvolume = deltax*deltax*deltax;

  if(DEBUG) fprintf(stderr,"GridFoo:  %f %f %d %e %e\n",
		    gridlex[gridnum],gridrex[gridnum],griddx[gridnum],
		    deltax,cellvolume);

  // open file - only once
  fprintf(stderr,"AddParticlesToArrays %s %d\n",gridfilenames[gridnum],
	  gridnump[gridnum]);

  if(DEBUG) fprintf(stderr,"about to open %s\n",gridfilenames[gridnum]);

  file_id = H5Fopen(gridfilenames[gridnum], H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );

  if(DEBUG) fprintf(stderr,"file opened %s\n",gridfilenames[gridnum]);

  /*------------------------------ NOTE -------------------------------
    The particle position stuff is 64 bit, whereas the creation time and
    particle masses are 32 bit.  Get particle positions and then go back and
    get the creation time and particle mass info.  */

  // get particle positions - 64 bit!
  pposx_dset_id = H5Dopen(file_id, "particle_position_x");
  assert( pposx_dset_id != h5_error);

  pposy_dset_id = H5Dopen(file_id, "particle_position_y");  
  assert( pposy_dset_id != h5_error);

  pposz_dset_id = H5Dopen(file_id, "particle_position_z");  
  assert( pposz_dset_id != h5_error);

  // open particle pos x dataspace (to get dimensions) 
  dsp_id = H5Dget_space(pposx_dset_id);
  assert( dsp_id != h5_error );

  // get data type (only once!)
  typ_id = H5Dget_type(pposx_dset_id);
  assert( typ_id != h5_error );

  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  // from the dimensional information, calculate the size of the buffer.
  // only once!
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

  if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) )
    {

      if(DEBUG) fprintf(stderr,"(3)\n");
      particle_posx_buffer = new double[(int) size];
      particle_posy_buffer = new double[(int) size];
      particle_posz_buffer = new double[(int) size];

      mem_type_id = H5T_IEEE_F64BE;

      // read particle position x field into an array
      h5_status = H5Dread(pposx_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_posx_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle pos x field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      // read particle position y field into an array
      h5_status = H5Dread(pposy_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_posy_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle pos y field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      // read particle position z field into an array
      h5_status = H5Dread(pposz_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_posz_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle pos z field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

    }  //   if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) )...

  // close hdf5 file, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );

  // close all position datasets
  h5_status = H5Dclose(pposx_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pposy_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pposz_dset_id);
  assert( h5_status != h5_error );


  // ------------------ get particle mass, dyn time, creat time, metal fraction info 

  pmass_dset_id = H5Dopen(file_id, "particle_mass");  
  assert( pmass_dset_id != h5_error);

  pxvel_dset_id = H5Dopen(file_id, "particle_velocity_x");  
  assert( pxvel_dset_id != h5_error);

  pyvel_dset_id = H5Dopen(file_id, "particle_velocity_y");  
  assert( pyvel_dset_id != h5_error);

  pzvel_dset_id = H5Dopen(file_id, "particle_velocity_z");  
  assert( pzvel_dset_id != h5_error);

  pcrtime_dset_id = H5Dopen(file_id, "creation_time");  
  assert( pmass_dset_id != h5_error);

  pdyntime_dset_id = H5Dopen(file_id, "dynamical_time");  
  assert( pmass_dset_id != h5_error);

  pmfrac_dset_id = H5Dopen(file_id, "metallicity_fraction");  
  assert( pmass_dset_id != h5_error);


  // open particle mass dataspace (to get dimensions) 
  // only once!
  dsp_id = H5Dget_space(pmass_dset_id);
  assert( dsp_id != h5_error );

  // get data type (only once!)
  typ_id = H5Dget_type(pmass_dset_id);
  assert( typ_id != h5_error );

  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  // from the dimensional information, calculate the size of the buffer.
  // only once!
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

  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
    {

      particle_xvel_buffer = new float[(int) size];
      particle_yvel_buffer = new float[(int) size];
      particle_zvel_buffer = new float[(int) size];
      particle_mass_buffer = new float[(int) size];
      particle_crtime_buffer = new float[(int) size];
      particle_dyntime_buffer = new float[(int) size];
      particle_mfrac_buffer = new float[(int) size];

      mem_type_id = H5T_IEEE_F32BE;

      // read particle mass field into an array
      h5_status = H5Dread(pmass_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_mass_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle mass field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      // read particle xvel field into an array
      h5_status = H5Dread(pxvel_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_xvel_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle xvel field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      // read particle yvel field into an array
      h5_status = H5Dread(pyvel_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_yvel_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle yvel field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      // read particle zvel field into an array
      h5_status = H5Dread(pzvel_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_zvel_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle zvel field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      // read particle crtime field into an array
      h5_status = H5Dread(pcrtime_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_crtime_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle crtime field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      // read particle dyntime field into an array
      h5_status = H5Dread(pdyntime_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_dyntime_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle dyntime field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      // read particle mfrac field into an array
      h5_status = H5Dread(pmfrac_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_mfrac_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle mfrac field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


    }  //  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )

  // close hdf5 data, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Dclose(pmass_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pxvel_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pyvel_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pzvel_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pcrtime_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pdyntime_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pmfrac_dset_id);
  assert( h5_status != h5_error );

  // close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  numberofparticles = ((int) size);

  assert( numberofparticles == gridnump[gridnum] );

  for(i=0; i<numberofparticles; i++)
    if(array_position < total_number_particles){
      //  for(i=array_position; i<array_position+numberofparticles; i++){

      xpositions[array_position] = (float) particle_posx_buffer[i];
      ypositions[array_position] = (float) particle_posy_buffer[i];
      zpositions[array_position] = (float) particle_posz_buffer[i];

      xvelocities[array_position] = particle_xvel_buffer[i];
      yvelocities[array_position] = particle_yvel_buffer[i];
      zvelocities[array_position] = particle_zvel_buffer[i];

      particlemasses[array_position] = particle_mass_buffer[i]*cellvolume;
      particlecrtime[array_position] = particle_crtime_buffer[i];
      particledyntime[array_position] = particle_dyntime_buffer[i];
      particlemetalfrac[array_position] = particle_mfrac_buffer[i];
      
      array_position++;
      particle_count++;

      if(array_position >= total_number_particles)
	fprintf(stderr,"HEY KIDS!  TOTAL NUM PARTICLES REACHED!\n\n");
    }

  // clean up arrays!
  delete [] particle_posx_buffer;
  delete [] particle_posy_buffer;
  delete [] particle_posz_buffer;

  delete [] particle_xvel_buffer;
  delete [] particle_yvel_buffer;
  delete [] particle_zvel_buffer;

  delete [] particle_mass_buffer;
  delete [] particle_crtime_buffer;
  delete [] particle_dyntime_buffer;
  delete [] particle_mfrac_buffer;


  if(DEBUG) fprintf(stderr,"exiting AddParticlesToArrays\n");
  return SUCCESS;


}
