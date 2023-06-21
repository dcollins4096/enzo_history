/***************************************************************************
 *                                                                         *
 * Copyright 2006 Brian O'Shea                                             *
 * Copyright 2006 Laboratory for Computational Astrophysics                *
 * Copyright 2006 Regents of the University of California                  *
 *                                                                         *
 * This software is released under the terms of the "University of         *
 * California/BSD License" in the accompanying LICENSE file.               *
 *                                                                         *
 ***************************************************************************/
#include "EnzoProj.h"

// external variables declared in other files
extern int projaxis,projlevel,*gridlevel,*griddx,*griddy,*griddz, *gridnump,
  rootgridsize,xprojnumcells,yprojnumcells,zprojnumcells,starformation;

extern double *gridlex,*gridley,*gridlez, *gridrex,*gridrey,*gridrez, 
  xstart,ystart,zstart, xend,yend,zend;

extern float boxsize,redshift,initialredshift,
  hubble,omegamatter,*projbuff_star,*projbuff_dmdensity;

extern char **gridfilenames;

// external variables declared in this file
float *particle_mass_buffer, *particle_crtime_buffer;
double *particle_posx_buffer, *particle_posy_buffer, 
  *particle_posz_buffer;
int numberofparticles,totalnumberofstars=0;

// functions in this file
int ProjectParticlesToPlane_X(int gridnum);
int ProjectParticlesToPlane_Y(int gridnum);
int ProjectParticlesToPlane_Z(int gridnum);

/*------------------ AddParticlesToProjection -----------------------
 *
 * we want ALL grids, not just grids with higher levels!  Use dark matter
 * and star particles, assuming that both exist.
 *
 *--------------------------------------------------------------*/
int AddParticlesToProjection(int gridnum){
  if(DEBUG) fprintf(stderr,"in AddParticlesToProjection %d\n",gridnum);

  // is any part of grid within bounds of projection?  if not, exit!
  if( !((gridrex[gridnum] > xstart) && (gridlex[gridnum] < xend) && 
	(gridrey[gridnum] > ystart) && (gridley[gridnum] < yend) && 
	(gridrez[gridnum] > zstart) && (gridlez[gridnum] < zend))
      ) return SUCCESS;

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
    pposz_dset_id,pmass_dset_id,pcrtime_dset_id;

  int i, ndims, m;

  // open file - only once
  fprintf(stderr,"AddParticlesToProjection %s %d\n",gridfilenames[gridnum],
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
#else
  pposx_dset_id = H5Dopen(file_id, "particle_position_x");
#endif
  assert( pposx_dset_id != h5_error);

#ifdef PACK_AMR
  pposy_dset_id = H5Dopen(group_id, "particle_position_y");  
#else
  pposy_dset_id = H5Dopen(file_id, "particle_position_y");  
#endif
  assert( pposy_dset_id != h5_error);

#ifdef PACK_AMR
  pposz_dset_id = H5Dopen(group_id, "particle_position_z");  
#else
  pposz_dset_id = H5Dopen(file_id, "particle_position_z");  
#endif
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

  //  if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) )
  //  {

  if(DEBUG) fprintf(stderr,"(3)\n");
  particle_posx_buffer = new double[(int) size];
  particle_posy_buffer = new double[(int) size];
  particle_posz_buffer = new double[(int) size];

  mem_type_id = H5T_NATIVE_DOUBLE;
  double *d_read_buff = new double[(int) size];

  // read particle position x field into an array
  h5_status = H5Dread(pposx_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, d_read_buff);
  if(DEBUG) fprintf(stderr,"float read status %d for particle pos x field\n", 
		    (int) h5_status);
  assert( h5_status != h5_error ); 
  for(m=0; m < size; m++) { particle_posx_buffer[m] = d_read_buff[m]; }

  fprintf(stderr,"aaah:  %.12lf %.12lf %.12lf %.12lf\n",particle_posx_buffer[0],particle_posx_buffer[1],
	  particle_posx_buffer[size-2],particle_posx_buffer[size-1]);

  // read particle position y field into an array
  h5_status = H5Dread(pposy_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, d_read_buff);
  if(DEBUG) fprintf(stderr,"float read status %d for particle pos y field\n", 
		    (int) h5_status);
  assert( h5_status != h5_error ); 
  for(m=0; m < size; m++) { particle_posy_buffer[m] = d_read_buff[m]; }

  // read particle position z field into an array
  h5_status = H5Dread(pposz_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, d_read_buff);
  if(DEBUG) fprintf(stderr,"float read status %d for particle pos z field\n", 
		    (int) h5_status);
  assert( h5_status != h5_error ); 
  for(m=0; m < size; m++) { particle_posz_buffer[m] = d_read_buff[m]; }
  
  //    }  //   if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) )...

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
#else
  pmass_dset_id = H5Dopen(file_id, "particle_mass");  
#endif
  assert( pmass_dset_id != h5_error);

  if(starformation==1){
#ifdef PACK_AMR
    pcrtime_dset_id = H5Dopen(group_id, "creation_time");  
#else
    pcrtime_dset_id = H5Dopen(file_id, "creation_time");  
#endif
    assert( pcrtime_dset_id != h5_error);
  }

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

#ifdef R8
  double *read_buff = new double[(int) size];
  mem_type_id = H5T_NATIVE_DOUBLE;
  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )  {
    fprintf(stderr, "TYPE MISMATCH!\n");
    fprintf(stderr, "Input Type is F32BE\n");
    exit(1);
  }
#endif

#ifdef R4
  float *read_buff = new float[(int) size];
  mem_type_id = H5T_NATIVE_FLOAT;
  if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) )  {
    fprintf(stderr, "TYPE MISMATCH!\n");
    fprintf(stderr, "Input Type is F64BE\n");
    exit(1);
  }
#endif

  //  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
  //{

  particle_mass_buffer = new float[(int) size];
  if(starformation==1)
    particle_crtime_buffer = new float[(int) size];
  
  //      mem_type_id = H5T_IEEE_F32BE;

  if(starformation==1){
    // read particle creation time field into an array
    h5_status = H5Dread(pcrtime_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(DEBUG) fprintf(stderr,"float read status %d for particle crtime field\n", 
		      (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { particle_crtime_buffer[m] = (float) read_buff[m]; }
  }

  // read particle mass field into an array
  h5_status = H5Dread(pmass_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, read_buff);
  if(DEBUG) fprintf(stderr,"float read status %d for particle mass field\n", 
		    (int) h5_status);
  assert( h5_status != h5_error ); 
  for(m=0; m < size; m++) { particle_mass_buffer[m] = (float) read_buff[m]; }
  
  //    }  //  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )

  // close hdf5 data, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );

  // close datasets
  if(starformation==1){
    h5_status = H5Dclose(pcrtime_dset_id);
    assert( h5_status != h5_error );
  }
  
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

  // project to appropriate axis!
  if(projaxis==0) ProjectParticlesToPlane_X(gridnum);
  if(projaxis==1) ProjectParticlesToPlane_Y(gridnum);
  if(projaxis==2) ProjectParticlesToPlane_Z(gridnum);

  // clean up arrays!
  delete [] particle_posx_buffer;
  delete [] particle_posy_buffer;
  delete [] particle_posz_buffer;
  delete [] particle_mass_buffer;

  if(starformation==1)
    delete [] particle_crtime_buffer;

  delete [] d_read_buff;
  delete [] read_buff;

  if(DEBUG) fprintf(stderr,"exiting AddParticlesToProjection\n");
  return SUCCESS;
}


/*-------------------- ProjectParticlesToPlane_X ------------------------
 *
 * This routine adds information to projections along the x axis.
 * 
 *--------------------------------------------------------------*/
int ProjectParticlesToPlane_X(int gridnum){
  if(DEBUG) fprintf(stderr,"in ProjectParticlesToPlane_X\n");
  int i, starprojind_y,starprojind_z,starsthisgrid=0;
  double proj_delta,grid_delta;

  double DensityUnits, LengthUnits, DensityConversion, MassConversion;

  LengthUnits = (((double) boxsize) * MPC_CM / ((double) hubble) / (1.0 + ((double) redshift)));
  DensityUnits = ((double) omegamatter) * RHOCRIT * pow( (1.0 + ((double) redshift) ), 3.0)
    *((double) hubble) * ( (double) hubble);

  // grid spacing of the projections (code units)
  proj_delta = (yend - ystart) / ((double) yprojnumcells);

  // spacing of the grid (code units)
  grid_delta = (gridrey[gridnum] - gridley[gridnum]) / ( (double) griddy[gridnum]);

  // density conversion
  DensityConversion = proj_delta * DensityUnits * LengthUnits / 
    MSOLAR_G * MPC_CM * MPC_CM;

  /* mass conversion - particle masses are stored as densities - therefore, since
     we're going from cells of different sizes, we need to multiply by the ratio
     of their volumes to renormalize the "masses".  This makes no assumption about
     the size of the grid cell compared to the size of the projection cell */

  MassConversion = pow( (grid_delta/proj_delta),3.0);

  if(DEBUG) fprintf(stderr,"ProjectParticlesToPlane_X: %lf %lf %e %e %i\n", 
		    proj_delta,grid_delta,DensityConversion,MassConversion,numberofparticles);

  /* loop over all particles.  If the particle is within the projection
     volume, calculate the grid cell that it will be in, and then
     depending on whether or not it's a star particle or dark matter
     particle, put it in the correct projection array.  */

  for(i=0; i < numberofparticles; i++) 
    if(( particle_posx_buffer[i] >= xstart ) && 
       ( particle_posx_buffer[i] <  xend ) && 
       ( particle_posy_buffer[i] >= ystart ) && 
       ( particle_posy_buffer[i] <  yend ) && 
       ( particle_posz_buffer[i] >= zstart ) && 
       ( particle_posz_buffer[i] <  zend )){

      // calculate indices
      starprojind_y = ((int) (( particle_posy_buffer[i] - ystart) / (yend - ystart) 
			      * ((double) yprojnumcells) ));
      starprojind_z = ((int) (( particle_posz_buffer[i] - zstart) / (zend - zstart) 
		       * ((double) zprojnumcells) ));

      // Is star formation turned on?  If so, use creation time as the
      // way to determine whether a particle is DM or star.
      if(starformation==1){ 
	if(particle_crtime_buffer[i] >= 0.0){  // then it's a star particle

	  // weight the appropriate grid cell, playing appropriate stupid 
	  // tricks with precision
	  projbuff_star[ starprojind_z * yprojnumcells + starprojind_y ] +=
	    ((float)
	     (((double)particle_mass_buffer[i])* 
	      MassConversion*DensityConversion));
	  
	  // keep track of the star particles formed!
	  starsthisgrid++;
	  totalnumberofstars++;
	} else { // Otherwise it's a Dm particle
	  	  
	  projbuff_dmdensity[starprojind_z * yprojnumcells + starprojind_y] +=
	    ((float)
	     (((double)particle_mass_buffer[i])* 
	      MassConversion*DensityConversion));

	}
      } else {  // if starformation == 0, then there are only dark matter ptcles
	
	projbuff_dmdensity[starprojind_z * yprojnumcells + starprojind_y] +=
	  ((float)
	   (((double)particle_mass_buffer[i])* 
	    MassConversion*DensityConversion));
      }
    } // for(i=0; i < numberofparticles; i++)
  
  if(DEBUG) fprintf(stderr,"ProjectParticlesToPlane_X:  stars this grid:  %i  total stars:  %i\n",
		    starsthisgrid,totalnumberofstars);
  if(DEBUG) fprintf(stderr,"leaving ProjectParticlesToPlane_X\n");
  return SUCCESS;
}


/*-------------------- ProjectsStarsToPlane_Y ------------------------
 *
 * This routine adds information to projections along the y axis.
 * 
 *--------------------------------------------------------------*/
int ProjectParticlesToPlane_Y(int gridnum){
  if(DEBUG) fprintf(stderr,"in ProjectParticlesToPlane_Y\n");

  int i, starprojind_x,starprojind_z,starsthisgrid=0;
  double proj_delta,grid_delta;

  double DensityUnits, LengthUnits, DensityConversion, MassConversion;

  LengthUnits = (((double) boxsize) * MPC_CM / ((double) hubble) / (1.0 + ((double) redshift)));
  DensityUnits = ((double) omegamatter) * RHOCRIT * pow( (1.0 + ((double) redshift) ), 3.0)
    *((double) hubble) * ( (double) hubble);

  // grid spacing of the projections (code units)
  proj_delta = (xend - xstart) / ((double) xprojnumcells);

  // spacing of the grid (code units)
  grid_delta = (gridrex[gridnum] - gridlex[gridnum]) / ( (double) griddx[gridnum]);

  // density conversion
  DensityConversion =  proj_delta * DensityUnits * LengthUnits / 
    MSOLAR_G * MPC_CM * MPC_CM;

  /* mass conversion - particle masses are stored as densities - therefore, since
     we're going from cells of different sizes, we need to multiply by the ratio
     of their volumes to renormalize the "masses".  This makes no assumption about
     the size of the grid cell compared to the size of the projection cell */

  MassConversion = pow( ((double) (grid_delta/proj_delta)),3.0);

  if(DEBUG) fprintf(stderr,"ProjectParticlesToPlane_X: %lf %lf %e %e %i\n", 
		    proj_delta,grid_delta,DensityConversion,MassConversion,numberofparticles);

  /* loop over all particles.  If the particle is within the projection
     volume, calculate the grid cell that it will be in, and then
     depending on whether or not it's a star particle or dark matter
     particle, put it in the correct projection array.  */

  for(i=0; i < numberofparticles; i++) 
    if(( particle_posx_buffer[i] >= xstart ) && 
       ( particle_posx_buffer[i] <  xend ) && 
       ( particle_posy_buffer[i] >= ystart ) && 
       ( particle_posy_buffer[i] <  yend ) && 
       ( particle_posz_buffer[i] >= zstart ) && 
       ( particle_posz_buffer[i] <  zend )){

      // calculate indices
      starprojind_x = ((int) (( particle_posx_buffer[i] - xstart) / (xend - xstart) 
			      * ((double) xprojnumcells) ));
      starprojind_z = ((int) (( particle_posz_buffer[i] - zstart) / (zend - zstart) 
		       * ((double) zprojnumcells) ));

      // Is star formation turned on?  If so, use creation time as the
      // way to determine whether a particle is DM or star.
      if(starformation==1){ 
	if(particle_crtime_buffer[i] >= 0.0){  // then it's a star particle

	  // weight the appropriate grid cell, playing appropriate stupid 
	  // tricks with precision
	  projbuff_star[ starprojind_z * xprojnumcells + starprojind_x ] +=
	    ((float)
	     (((double)particle_mass_buffer[i])* 
	      MassConversion*DensityConversion));
	  
	  // keep track of the star particles formed!
	  starsthisgrid++;
	  totalnumberofstars++;
	} else { // Otherwise it's a Dm particle
	  	  
	  projbuff_dmdensity[starprojind_z * xprojnumcells + starprojind_x] +=
	    ((float)
	     (((double)particle_mass_buffer[i])* 
	      MassConversion*DensityConversion));

	}	
      } else {  // if starformation == 0, then there are only dark matter ptcles
	
	projbuff_dmdensity[starprojind_z * xprojnumcells + starprojind_x] +=
	  ((float)
	   (((double)particle_mass_buffer[i])* 
	    MassConversion*DensityConversion));
      }
    }  //   for(i=0; i < numberofparticles; i++) ...

  if(DEBUG) fprintf(stderr,"leaving ProjectParticlesToPlane_Y\n");
  return SUCCESS;
}


/*-------------------- ProjectParticlesToPlane_Z ------------------------
 *
 * This routine adds information to projections along the z axis.
 * 
 *--------------------------------------------------------------*/
int ProjectParticlesToPlane_Z(int gridnum){
  if(DEBUG) fprintf(stderr,"in ProjectParticlesToPlane_Z\n");

  int i, starprojind_y,starprojind_x,starsthisgrid=0;
  double proj_delta,grid_delta;

  double DensityUnits, LengthUnits, DensityConversion, MassConversion;

  LengthUnits = (((double) boxsize) * MPC_CM / ((double) hubble) / (1.0 + ((double) redshift)));
  DensityUnits = ((double) omegamatter) * RHOCRIT * pow( (1.0 + ((double) redshift) ), 3.0)
    *((double) hubble) * ( (double) hubble);

  // grid spacing of the projections (code units)
  proj_delta = (yend - ystart) / ((double) yprojnumcells);

  // spacing of the grid (code units)
  grid_delta = (gridrey[gridnum] - gridley[gridnum]) / ( (double) griddy[gridnum]);

  // density conversion
  DensityConversion =  proj_delta * DensityUnits * LengthUnits / 
    MSOLAR_G * MPC_CM * MPC_CM;

  /* mass conversion - particle masses are stored as densities - therefore, since
     we're going from cells of different sizes, we need to multiply by the ratio
     of their volumes to renormalize the "masses".  This makes no assumption about
     the size of the grid cell compared to the size of the projection cell */

  MassConversion = pow( ((double) (grid_delta/proj_delta)),3.0);

  if(DEBUG) fprintf(stderr,"ProjectParticlesToPlane_X: %lf %lf %e %e %i\n", 
		    proj_delta,grid_delta,DensityConversion,MassConversion,numberofparticles);

  /* loop over all particles.  If the particle is within the projection
     volume, calculate the grid cell that it will be in, and then
     depending on whether or not it's a star particle or dark matter
     particle, put it in the correct projection array.  */

  for(i=0; i < numberofparticles; i++) 
    if((particle_posx_buffer[i] >= xstart ) && 
       (particle_posx_buffer[i] <  xend ) && 
       (particle_posy_buffer[i] >= ystart ) && 
       (particle_posy_buffer[i] <  yend ) && 
       (particle_posz_buffer[i] >= zstart ) && 
       (particle_posz_buffer[i] <  zend )){

      // calculate indices
      starprojind_y = ((int) (( particle_posy_buffer[i] - ystart) / (yend - ystart) 
			      * ((double) yprojnumcells) ));
      starprojind_x = ((int) (( particle_posx_buffer[i] - xstart) / (xend - xstart) 
			      * ((double) xprojnumcells) ));

      // Is star formation turned on?  If so, use creation time as the
      // way to determine whether a particle is DM or star.
      if(starformation==1){ 
	if(particle_crtime_buffer[i] >= 0.0){  // then it's a star particle

	  // weight the appropriate grid cell, playing appropriate stupid 
	  // tricks with precision
	  projbuff_star[ starprojind_y * xprojnumcells + starprojind_x ] +=
	    ((float)
	     (((double)particle_mass_buffer[i])* 
	      MassConversion*DensityConversion));
	  
	  // keep track of the star particles formed!
	  starsthisgrid++;
	  totalnumberofstars++;
	} else { // Otherwise it's a Dm particle
	  	  
	  projbuff_dmdensity[starprojind_y * xprojnumcells + starprojind_x] +=
	    ((float)
	     (((double)particle_mass_buffer[i])* 
	      MassConversion*DensityConversion));

	}	
      } else {  // if starformation == 0, then there are only dark matter ptcles
	
	projbuff_dmdensity[starprojind_y * xprojnumcells + starprojind_x] +=
	  ((float)
	   (((double)particle_mass_buffer[i])* 
	    MassConversion*DensityConversion));
      }
    }  //   for(i=0; i < numberofparticles; i++) ...

  if(DEBUG) fprintf(stderr,"leaving ProjectParticlesToPlane_Z\n");
  return SUCCESS;
}
