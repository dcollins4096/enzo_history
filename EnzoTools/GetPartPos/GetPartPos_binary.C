/* GetPartPos_binary.C
   written by BWO

   this reads in Enzo HDF5 datasets (packed or not packed) and writes out a bunch
   of dark matter particle information in a manner that is consistent with Katrin's
   "cosmo" format.  This has to be further modified by an IDL script (which should be
   in the same directory as this program) called convdm.pro

   Note that if the data is in PACKED AMR FORMAT you may have to modify the line that says
   "#define PACK_AMR" below.

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>

#define hssize_t hsize_t

#define PACK_AMR  // set to PACK_AMR or NO_PACK_AMR  -------------

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

float *xpositions, *ypositions, *zpositions, *xvelocities, *yvelocities, *zvelocities;
int *gridnump, array_position=0, particle_count=0,totalparticlesfromuser=0;
char **gridfilenames;

int main(int argc, char *argv[]){

  char *inputfilename=argv[1];  // name of Enzo hierarchy
  int numpart;
  int i;
  FILE *outputfile;

  numpart=atoi( argv[2] );   // number of particles (total)

  if(argc != 3){
    fprintf(stderr,"I need two arguments!  getpartpos_bin <paramfile> <numpart>\n");
    exit(-123);
  }

  /* this is done because sometimes the "real" number of particles is NOT
     N^3 - there are spurious particles.  So, we're going to just throw out
     the last few particles on the list. */
  totalparticlesfromuser=numpart;  

  printf("%s %i\n",inputfilename,totalparticlesfromuser);
  fflush(stdout);

  //ReadParameterFile(inputfilename);

  // get hierarchy file name
  int inlen = (int) strlen(inputfilename);
  char *hierfile = new char[inlen+13];
  strcpy(hierfile,inputfilename);
  hierfile=strcat(hierfile,".hierarchy");

  int total_number_grids = NumberOfGrids(hierfile);

  printf("there are %i total particles\n",total_number_particles);

  xpositions = new float[totalparticlesfromuser];
  ypositions = new float[totalparticlesfromuser];
  zpositions = new float[totalparticlesfromuser];


  xvelocities = new float[totalparticlesfromuser];
  yvelocities = new float[totalparticlesfromuser];
  zvelocities = new float[totalparticlesfromuser];

  for(i=0; i<totalparticlesfromuser; i++)
    xpositions[i] = ypositions[i] = zpositions[i] = 
      xvelocities[i] = yvelocities[i] = zvelocities[i] = -1.0;



  // grid file name info
  gridfilenames = new char*[total_number_grids];

  // number of particles info
  gridnump = new int[total_number_grids];

  for(i=0; i<total_number_grids; i++){
    gridnump[i]=-1;
    gridfilenames[i] = new char[MAX_LINE_LENGTH];
   }

  GetGridInfo(total_number_grids,hierfile);

  for(i=0; i<total_number_grids; i++)
    if(gridnump[i]>0)
      AddParticlesToArrays(i);


  for(i=0; i<totalparticlesfromuser; i++)
    if(xpositions[i] < 0.0 || xpositions[i] > 1.0 ||
       ypositions[i] < 0.0 || ypositions[i] > 1.0 ||
       zpositions[i] < 0.0 || zpositions[i] > 1.0)
      fprintf(stderr,"uh-oh %f %f %f\n",xpositions[i],ypositions[i],zpositions[i]);

   
  outputfile = fopen("dmposinfo.dat","w");

  fprintf(stderr,"writing file dmposinfo.dat\n") ;

  
  for(i=0; i<totalparticlesfromuser; i++){
    fwrite(&xpositions[i],sizeof(float),1,outputfile);
    fwrite(&xvelocities[i],sizeof(float),1,outputfile);
    fwrite(&ypositions[i],sizeof(float),1,outputfile);
    fwrite(&yvelocities[i],sizeof(float),1,outputfile);
    fwrite(&zpositions[i],sizeof(float),1,outputfile);
    fwrite(&zvelocities[i],sizeof(float),1,outputfile);
  }
   
  fclose(outputfile);

  delete [] xpositions;
  delete [] ypositions;
  delete [] zpositions;
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

  int nbfields,particleonly;
  
  if(DEBUG) fprintf(stderr,"hierarchy file: %s\n",hierfilename);

  // open hierarchy file
  hierfile = fopen(hierfilename,"r");

  // read through hierarchy file, get grid info
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){

    // lines that start with "Grid =" indicate a new set of grid info
    if(strncmp(line,"Grid = ",7)==0){
      
      for(i=0;i<9;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
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
  hid_t       file_id, group_id;
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

  hid_t       pposx_dset_id,pposy_dset_id,
    pposz_dset_id,pmass_dset_id,pcrtime_dset_id,  pvelx_dset_id,pvely_dset_id,
    pvelz_dset_id;

  int i, ndims;
  float *particle_mass_buffer, *particle_crtime_buffer;
  double *particle_posx_buffer, *particle_posy_buffer, 
    *particle_posz_buffer;

  float  *particle_velx_buffer, *particle_vely_buffer, 
    *particle_velz_buffer;

  int numberofparticles,totalnumberofstars=0;


  // open file - only once
  fprintf(stderr,"AddParticlesToArrays %s %d\n",gridfilenames[gridnum],
	  gridnump[gridnum]);

  file_id = H5Fopen(gridfilenames[gridnum], H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );


#ifdef PACK_AMR
  //AK open group
  char id[10];
  sprintf(id, "%8.8d", gridnum+1);
  char name[MAX_LINE_LENGTH];
  strcpy(name, "/Grid");
  strcat(name, id);
  printf("H5Gopen with Name %s\n", name);
  group_id = H5Gopen(file_id, name);
  assert( group_id != h5_error );
#endif

  /*------------------------------ NOTE -------------------------------
    The particle position stuff is 64 bit, whereas the creation time and
    particle masses are 32 bit.  Get particle positions and then go back and
    get the creation time and particle mass info.  */

  // get particle positions - 64 bit!
#ifdef PACK_AMR

  pposx_dset_id = H5Dopen(group_id, "particle_position_x");
  assert( pposx_dset_id != h5_error);

  pposy_dset_id = H5Dopen(group_id, "particle_position_y");  
  assert( pposy_dset_id != h5_error);

  pposz_dset_id = H5Dopen(group_id, "particle_position_z");  
  assert( pposz_dset_id != h5_error);

#else /* PACK_AMR */

  pposx_dset_id = H5Dopen(file_id, "particle_position_x");
  assert( pposx_dset_id != h5_error);

  pposy_dset_id = H5Dopen(file_id, "particle_position_y");  
  assert( pposy_dset_id != h5_error);

  pposz_dset_id = H5Dopen(file_id, "particle_position_z");  
  assert( pposz_dset_id != h5_error);

#endif /* PACK_AMR */

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


  // ------------------ get particle mass and creation time info 


#ifdef PACK_AMR

  pmass_dset_id = H5Dopen(group_id, "particle_mass");  
  assert( pmass_dset_id != h5_error);


  pvelx_dset_id = H5Dopen(group_id, "particle_velocity_x");  
  assert( pvelx_dset_id != h5_error);

  pvely_dset_id = H5Dopen(group_id, "particle_velocity_y");  
  assert( pvely_dset_id != h5_error);

  pvelz_dset_id = H5Dopen(group_id, "particle_velocity_z");  
  assert( pvelz_dset_id != h5_error);

#else /* PACK_AMR */

  pmass_dset_id = H5Dopen(file_id, "particle_mass");  
  assert( pmass_dset_id != h5_error);


  pvelx_dset_id = H5Dopen(file_id, "particle_velocity_x");  
  assert( pvelx_dset_id != h5_error);

  pvely_dset_id = H5Dopen(file_id, "particle_velocity_y");  
  assert( pvely_dset_id != h5_error);

  pvelz_dset_id = H5Dopen(file_id, "particle_velocity_z");  
  assert( pvelz_dset_id != h5_error);

#endif /* PACK_AMR */

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

      particle_mass_buffer = new float[(int) size];
      particle_velx_buffer = new float[(int) size];
      particle_vely_buffer = new float[(int) size];
      particle_velz_buffer = new float[(int) size];

      mem_type_id = H5T_IEEE_F32BE;

      // read particle mass field into an array
      h5_status = H5Dread(pmass_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_mass_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle mass field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      h5_status = H5Dread(pvelx_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_velx_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle x-velocity \n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      h5_status = H5Dread(pvely_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_vely_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle y-velocity \n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      h5_status = H5Dread(pvelz_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_velz_buffer);
      if(DEBUG) fprintf(stderr,"float read status %d for particle z-velocity \n", 
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

  // close group
#ifdef PACK_AMR
  h5_status = H5Gclose(group_id);
  assert( h5_status != h5_error );
#endif

  // close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  numberofparticles = ((int) size);

  assert( numberofparticles == gridnump[gridnum] );

  for(i=0; i<numberofparticles; i++)
    if(array_position < totalparticlesfromuser){
      //  for(i=array_position; i<array_position+numberofparticles; i++){

      xpositions[array_position] = (float) particle_posx_buffer[i];
      ypositions[array_position] = (float) particle_posy_buffer[i];
      zpositions[array_position] = (float) particle_posz_buffer[i];

      xvelocities[array_position] = (float) particle_velx_buffer[i];
      yvelocities[array_position] = (float) particle_vely_buffer[i];
      zvelocities[array_position] = (float) particle_velz_buffer[i];
      
      array_position++;
      particle_count++;

      if(array_position >= totalparticlesfromuser)
	fprintf(stderr,"HEY KIDS!  TOTAL NUM PARTICLES REACHED!\n\n");
    }

  // clean up arrays!
  delete [] particle_posx_buffer;
  delete [] particle_posy_buffer;
  delete [] particle_posz_buffer;

  delete [] particle_velx_buffer;
  delete [] particle_vely_buffer;
  delete [] particle_velz_buffer;

  delete [] particle_mass_buffer;

  if(DEBUG) fprintf(stderr,"exiting AddParticlesToArrays\n");
  return SUCCESS;


}
