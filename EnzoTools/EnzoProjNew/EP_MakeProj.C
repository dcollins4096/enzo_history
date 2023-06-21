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
extern int projaxis,projlevel,*gridlevel,*griddx,*griddy,*griddz,
  rootgridsize,xprojnumcells,yprojnumcells,zprojnumcells,
  metal_flag,z1_flag,z2_flag,star_flag,multispecies;

extern double *gridlex,*gridley,*gridlez, *gridrex,*gridrey,*gridrez, 
  xstart,ystart,zstart, xend,yend,zend;

extern float boxsize,redshift,initialredshift,
  hubble,omegamatter,*projbuff_density,*projbuff_dmdensity,*projbuff_xraylum,
  *projbuff_tempxray, *projbuff_level, *projbuff_szy, *projbuff_szkin,
  *projbuff_metal, *projbuff_z1, *projbuff_z2, *projbuff_star,
  *projbuff_HIdensity,*projbuff_HIIdensity,*projbuff_HeIdensity,
  *projbuff_HeIIdensity,*projbuff_HeIIIdensity,*projbuff_electrondensity,
  *projbuff_HMdensity,*projbuff_H2Idensity,*projbuff_H2IIdensity, 
  *projbuff_xrayem, *projbuff_spectemp, *projbuff_spectemp_weight,
  *projbuff_masstemp, *projbuff_masstemp_weight;

extern char **gridfilenames;

// externs declared in THIS file
float *densitybuff,*dmdensbuff,*tempbuff,*metalbuff,*z1buff,*z2buff,
  *velocitybuff,*HIbuff,*HIIbuff,*HeIbuff,*HeIIbuff,*HeIIIbuff,*elecbuff,
  *HMbuff,*H2Ibuff,*H2IIbuff, *xrayem_buff;
int *flagbuffer,dmonthisgrid;

// functions that are only for use in this file!
int FlagGridCells(int gridnum,int total_number_grids);
int ProjectToPlane_X(int gridnum);
int ProjectToPlane_Y(int gridnum);
int ProjectToPlane_Z(int gridnum);


/*------------------ AddGridToProjection -----------------------
 *
 *  This function takes in a grid (assumed to be at level <= 
 *  the maximum projection level) and adds information from that
 *  grid to the 2D projection, after first checking to make sure
 *  that the grid actually lies within the bounds of the
 *  projection area.
 *
 *--------------------------------------------------------------*/
int AddGridToProjection(int gridnum,int total_number_grids){

  if(DEBUG) fprintf(stderr,"in AddGridToProjection %d\n",gridnum);

  if(gridlevel[gridnum] > projlevel) return SUCCESS;

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

  // need one of these per dataset!
  hid_t dens_dset_id,dmdens_dset_id,temp_dset_id,
    metal_dset_id,z1_dset_id,z2_dset_id,
    velocity_dset_id,HI_dset_id,HII_dset_id,HeI_dset_id,
    HeII_dset_id,HeIII_dset_id,e_dset_id,HM_dset_id,
    H2I_dset_id,H2II_dset_id;

  int i, ndims,m;

  // is any part of grid within bounds of projection?  if not, exit!
  if( !((gridrex[gridnum] > xstart) && (gridlex[gridnum] < xend) && 
	(gridrey[gridnum] > ystart) && (gridley[gridnum] < yend) && 
	(gridrez[gridnum] > zstart) && (gridlez[gridnum] < zend))
      ) return SUCCESS;

  // open grid file and extract various quantities of interest

  // open file - only once
  fprintf(stderr,"AddGridToProjection %s\n",gridfilenames[gridnum]);
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

  // open density dataset
#ifdef PACK_AMR
  dens_dset_id = H5Dopen(group_id, "Density");
#else
  dens_dset_id = H5Dopen(file_id, "Density");
#endif
  assert( dens_dset_id != h5_error );

  // open temperature dataset
#ifdef PACK_AMR
  temp_dset_id = H5Dopen(group_id, "Temperature");
#else
  temp_dset_id = H5Dopen(file_id, "Temperature");
#endif
  assert( dens_dset_id != h5_error );

  // open velocity dataset based on projection axis! -- for SZ calculations
  if(projaxis==0){
#ifdef PACK_AMR
    velocity_dset_id = H5Dopen(group_id, "x-velocity");
#else
    velocity_dset_id = H5Dopen(file_id, "x-velocity");
#endif
    assert( velocity_dset_id != h5_error );
  }

  if(projaxis==1){
#ifdef PACK_AMR
    velocity_dset_id = H5Dopen(group_id, "y-velocity");
#else
    velocity_dset_id = H5Dopen(file_id, "y-velocity");
#endif
    assert( velocity_dset_id != h5_error );
  }

  if(projaxis==2){
#ifdef PACK_AMR
    velocity_dset_id = H5Dopen(group_id, "z-velocity");
#else
    velocity_dset_id = H5Dopen(file_id, "z-velocity");
#endif
    assert( velocity_dset_id != h5_error );
  }

  if(metal_flag==1){
    // open metal density dataset
#ifdef PACK_AMR
    metal_dset_id = H5Dopen(group_id, "Metal_Density");
#else
    metal_dset_id = H5Dopen(file_id, "Metal_Density");
#endif
    assert( metal_dset_id != h5_error );
  }

  if(z1_flag==1){
    // open Z1 metal density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening Z_Field1\n");
#ifdef PACK_AMR
    z1_dset_id = H5Dopen(group_id, "Z_Field1");
#else
    z1_dset_id = H5Dopen(file_id, "Z_Field1");
#endif
    assert( z1_dset_id != h5_error );
  }

  if(z2_flag==1){
    // open Z2 metal density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening Z_Field2\n");
#ifdef PACK_AMR
    z2_dset_id = H5Dopen(group_id, "Z_Field2");
#else
    z2_dset_id = H5Dopen(file_id, "Z_Field2");
#endif
    assert( z2_dset_id != h5_error );
  }

  if((multispecies==1) || (multispecies==2)){
    // open neutral hydrogen density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening HI_Density\n");
#ifdef PACK_AMR
    HI_dset_id = H5Dopen(group_id, "HI_Density");
#else
    HI_dset_id = H5Dopen(file_id, "HI_Density");
#endif
    assert( HI_dset_id != h5_error );

    // open ionized hydrogen density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening HII_Density\n");
#ifdef PACK_AMR
    HII_dset_id = H5Dopen(group_id, "HII_Density");
#else
    HII_dset_id = H5Dopen(file_id, "HII_Density");
#endif
    assert( HII_dset_id != h5_error );

    // open neutral helium density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening HeI_Density\n");
#ifdef PACK_AMR
    HeI_dset_id = H5Dopen(group_id, "HeI_Density");
#else
    HeI_dset_id = H5Dopen(file_id, "HeI_Density");
#endif
    assert( HeI_dset_id != h5_error );

    // open singly ionized helium density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening HeII_Density\n");
#ifdef PACK_AMR
    HeII_dset_id = H5Dopen(group_id, "HeII_Density");
#else
    HeII_dset_id = H5Dopen(file_id, "HeII_Density");
#endif
    assert( HeII_dset_id != h5_error );

    // open doubly ionized helium density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening HeIII_Density\n");
#ifdef PACK_AMR
    HeIII_dset_id = H5Dopen(group_id, "HeIII_Density");
#else
    HeIII_dset_id = H5Dopen(file_id, "HeIII_Density");
#endif
    assert( HeIII_dset_id != h5_error );

    // open electron density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening Electron_Density\n");
#ifdef PACK_AMR
    e_dset_id = H5Dopen(group_id, "Electron_Density");
#else
    e_dset_id = H5Dopen(file_id, "Electron_Density");
#endif
    assert( e_dset_id != h5_error );
  }

  if(multispecies==2){
    // open negatively charged helium dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening HM_Density\n");
#ifdef PACK_AMR
    HM_dset_id = H5Dopen(group_id, "HM_Density");
#else
    HM_dset_id = H5Dopen(file_id, "HM_Density");
#endif
    assert( HM_dset_id != h5_error );

    // open molecular hydrogen density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening H2I_Density\n");
#ifdef PACK_AMR
    H2I_dset_id = H5Dopen(group_id, "H2I_Density");
#else
    H2I_dset_id = H5Dopen(file_id, "H2I_Density");
#endif
    assert( H2I_dset_id != h5_error );

    // open singly ionized molecular hydrogen density dataset
    if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection:  Opening H2II_Density\n");
#ifdef PACK_AMR
    H2II_dset_id = H5Dopen(group_id, "H2II_Density");
#else
    H2II_dset_id = H5Dopen(file_id, "H2II_Density");
#endif
    assert( H2II_dset_id != h5_error );
  }

  // open density dataspace (to get dimensions) 
  // only once!
  dsp_id = H5Dget_space(dens_dset_id);
  assert( dsp_id != h5_error );

  // get data type (only once!)
  typ_id = H5Dget_type(dens_dset_id);
  assert( typ_id != h5_error );

  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )  {
    fprintf(stderr, "Input Type is F32BE\n");
  }

  if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) )  {
    fprintf(stderr, "Input Type is F64BE\n");
  }
  
  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  // from the dimensional information, calculate the size of the buffer.
  // only once!
  size = 1;
  if(VERBOSEDEBUG) fprintf(stderr,"Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(VERBOSEDEBUG) fprintf(stderr," Dim %d\n", (int) xdims[i]);
    }
  if(VERBOSEDEBUG) fprintf(stderr,"Size %d\n", (int) size);

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
//    {
 
  // allocate buffers - one per projection buffer!
  densitybuff = new float[(int) size];
  tempbuff = new float[(int) size];
  velocitybuff  = new float[(int) size];
  if(metal_flag==1) metalbuff = new float[(int) size];
  if(z1_flag==1) z1buff = new float[(int) size];
  if(z2_flag==1) z2buff = new float[(int) size];
  if((multispecies==1) || (multispecies==2) ){
    HIbuff = new float[(int) size];
    HIIbuff = new float[(int) size];
    HeIbuff = new float[(int) size];
    HeIIbuff = new float[(int) size];
    HeIIIbuff = new float[(int) size];
    elecbuff = new float[(int) size];
  }

  if(multispecies==2){
    HMbuff = new float[(int) size];
    H2Ibuff = new float[(int) size];
    H2IIbuff = new float[(int) size];
  }

#ifdef MEKAL
  xrayembuff = new float[(int) size]
#endif

  //      mem_type_id = H5T_IEEE_F32BE;

  // read density field into an array
  h5_status = H5Dread(dens_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, read_buff);
  if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for density field\n", 
			   (int) h5_status);
  assert( h5_status != h5_error ); 
  for(m=0; m < size; m++) { densitybuff[m] = (float) read_buff[m]; }
  
  // read temperature field into an array
  h5_status = H5Dread(temp_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, read_buff);
  if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for temperature field\n", 
			   (int) h5_status);
  assert( h5_status != h5_error );      
  for(m=0; m < size; m++) { tempbuff[m] = (float) read_buff[m]; }
  
  // read velocity field into an array
  h5_status = H5Dread(velocity_dset_id, mem_type_id, 
		      mem_dsp_id, file_dsp_id, 
		      H5P_DEFAULT, read_buff);
  if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for velocity field\n", 
			   (int) h5_status);
  assert( h5_status != h5_error );      
  for(m=0; m < size; m++) { velocitybuff[m] = (float) read_buff[m]; }
  
  if(metal_flag==1){
    // read metal field into an array
    h5_status = H5Dread(metal_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for metal density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error );      
    for(m=0; m < size; m++) { metalbuff[m] = (float) read_buff[m]; }
  }
  
  if(z1_flag==1){
    // read z1 field into an array
    h5_status = H5Dread(z1_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for ZField_1 field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { z1buff[m] = (float) read_buff[m]; }
  }

  if(z2_flag==1){
    // read z2 field into an array
    h5_status = H5Dread(z2_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for ZField_2 field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error );      
    for(m=0; m < size; m++) { z2buff[m] = (float) read_buff[m]; }
  }
  

  if((multispecies==1) || (multispecies==2)){
    // read HI field into an array
    h5_status = H5Dread(HI_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for HI_Density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { HIbuff[m] = (float) read_buff[m]; }
    
    
    // read HII field into an array
    h5_status = H5Dread(HII_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for HII_Density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { HIIbuff[m] = (float) read_buff[m]; }
    
    // read HeI field into an array
    h5_status = H5Dread(HeI_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for HeI_Density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { HeIbuff[m] = (float) read_buff[m]; }
    
    // read HeII field into an array
    h5_status = H5Dread(HeII_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for HeII_Density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { HeIIbuff[m] = (float) read_buff[m]; }
    
    // read HeIII field into an array
    h5_status = H5Dread(HeIII_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for HeIII_Density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { HeIIIbuff[m] = (float) read_buff[m]; }
    
    
    // read e field into an array
    h5_status = H5Dread(e_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for Electron_Density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { elecbuff[m] = (float) read_buff[m]; }
    
  }


  if(multispecies==2){
    // read HM field into an array
    h5_status = H5Dread(HM_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for HM_Density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { HMbuff[m] = (float) read_buff[m]; }
    
    // read H2I field into an array
    h5_status = H5Dread(H2I_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for H2I_Density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { H2Ibuff[m] = (float) read_buff[m]; }
    
    // read H2II field into an array
    h5_status = H5Dread(H2II_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, read_buff);
    if(VERBOSEDEBUG) fprintf(stderr,"float read status %d for H2II_Density field\n", 
			     (int) h5_status);
    assert( h5_status != h5_error ); 
    for(m=0; m < size; m++) { H2IIbuff[m] = (float) read_buff[m]; }
    
  }
  
  //    } // end of type declaration
  
  // ---------- close hdf5 file, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );


  // ---------- must close each dataset - one per projection buffer!
  h5_status = H5Dclose(dens_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(temp_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(velocity_dset_id);
  assert( h5_status != h5_error );

  if(metal_flag==1){
    h5_status = H5Dclose(metal_dset_id);
    assert( h5_status != h5_error );
  }

  if(z1_flag==1){
    h5_status = H5Dclose(z1_dset_id);
    assert( h5_status != h5_error );
  }

  if(z2_flag==1){
    h5_status = H5Dclose(z2_dset_id);
    assert( h5_status != h5_error );
  }

  if((multispecies==1) || (multispecies==2)){
    h5_status = H5Dclose(HI_dset_id);
    assert( h5_status != h5_error );

    h5_status = H5Dclose(HII_dset_id);
    assert( h5_status != h5_error );

    h5_status = H5Dclose(HeI_dset_id);
    assert( h5_status != h5_error );

    h5_status = H5Dclose(HeII_dset_id);
    assert( h5_status != h5_error );

    h5_status = H5Dclose(HeIII_dset_id);
    assert( h5_status != h5_error );

    h5_status = H5Dclose(e_dset_id);
    assert( h5_status != h5_error );
  }

  if(multispecies==2){
    h5_status = H5Dclose(HM_dset_id);
    assert( h5_status != h5_error );

    h5_status = H5Dclose(H2I_dset_id);
    assert( h5_status != h5_error );

    h5_status = H5Dclose(H2II_dset_id);
    assert( h5_status != h5_error );
  }

  // ---------- close group
#ifdef PACK_AMR
  h5_status = H5Gclose(group_id);
  assert( h5_status != h5_error );
#endif

  // ---------- close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  // ---------- create flag field and zero it out
  flagbuffer = new int[(int) size];
  for(i=0; i < ((int) size) ; i++) flagbuffer[i] = 0;

  // ---------- check density buffer - debug operation
  for(i=0; i < ((int) size) ; i++){
    if(densitybuff[i]<=0.0) fprintf(stderr,"AddGridToProjection: negative or zero density value!\n");
  }

  // if this grid is NOT on maxprojectionlevel, flag cells that
  // have a grid below them
  if(gridlevel[gridnum]<projlevel) FlagGridCells(gridnum,total_number_grids);

  // add cells to hierarchy (based on projection axis)
  if(projaxis == 0) ProjectToPlane_X(gridnum);
  if(projaxis == 1) ProjectToPlane_Y(gridnum);
  if(projaxis == 2) ProjectToPlane_Z(gridnum);

  // clean up dynamically allocated arrays - must erase everything!
  if(VERBOSEDEBUG) fprintf(stderr,"AddGridToProjection: cleaning up\n");
 
  delete [] flagbuffer;
  delete [] densitybuff;
  delete [] tempbuff;
  delete [] velocitybuff;
  if(metal_flag == 1) delete [] metalbuff;
  if(z1_flag == 1) delete [] z1buff;
  if(z2_flag == 1) delete [] z2buff;

  if((multispecies==1) || (multispecies==2)){
    delete [] HIbuff;
    delete [] HIIbuff;
    delete [] HeIbuff;
    delete [] HeIIbuff;
    delete [] HeIIIbuff;
    delete [] elecbuff;
  }

  if(multispecies==2){
    delete [] HMbuff;
    delete [] H2Ibuff;
    delete [] H2IIbuff;
  }

#ifdef MEKAL
  delete [] xrayem_buff;
#endif

  delete [] read_buff;

  if(DEBUG) fprintf(stderr,"exiting AddGridToProjection\n");
  return SUCCESS;
}


/*--------------------- FlagGridCells --------------------------
 *
 *  This function is handed the number of a grid, and there is a
 *  buffer for flags.  We go through each of the OTHER grids, 
 *  only looking for grids that are at level L+1 (assuming this
 *  grid is on level L) and that is overlapping the grid in
 *  question in some way.  If this is true, we go over all of the
 *  cells in the grid and flag each one individually if that cell
 *  happens to be either a) covering a more refined region or
 *  b) if that cell is completely outside the bounds of the
 *  desired projection volume.
 *
 *  The algorithm chosen is probably not the most effective one,
 *  though I can't think of a better one offhand.  If you've
 *  actually taken the time to read the source and can
 *  significantly improve on this, please email me your idea at
 *  bwoshea@cosmos.ucsd.edu, and I'll be very, very grateful.
 *
 *--------------------------------------------------------------*/
int FlagGridCells(int gridnum, int total_number_grids){

  if(VERBOSEDEBUG) fprintf(stderr,"in FlagGridCells: %i %i\n",gridnum,total_number_grids);

  int counter, flagged_cells_this_grid,i,j,k,cellindex;

  double clex, cley, clez, crex, crey, crez;  // cell left and right edges (x,y,z)

  flagged_cells_this_grid=0;  // reset counter to zero

  /* loop over all grids and flag cells that are covered by a grid of
     higher resolution.  We can reject a lot of grids offhand because
     they are:
     1) not of level L+1 (where the grid of interest is level L) or
     2) outside of the bounds of the grid of interest
   */
  for(counter=0; counter < total_number_grids; counter++){
    
    // if the grid of interest is on level l, we only want to
    // look at grids on level l+1 to avoid working too hard
    // (and to make our lives simpler)
    // also takes care of getting rid of grid[gridnum]
    // so we don't flag ALL of the cells!
    if(gridlevel[counter] != (gridlevel[gridnum]+1)) continue;

    // does grid i overlap our grid of interest?  If not, skip it.
    if( !((gridrex[counter] > gridlex[gridnum]) && (gridlex[counter] < gridrex[gridnum]) &&
	  (gridrey[counter] > gridley[gridnum]) && (gridley[counter] < gridrey[gridnum]) &&
	  (gridrez[counter] > gridlez[gridnum]) && (gridlez[counter] < gridrez[gridnum])) 
	) continue;

    // if we've passed all of those tests, grid[counter] is
    //   1)  on level L+1 (where grid[gridnum] is on L)
    //   2)  guaranteed to overlap grid[gridnum] in some way
    // so we loop over every cell in grid[gridnum] and flag all cells that
    // grid[counter] overlaps, and also if there are any cells that are
    // outside the bounds.

    for(i=0; i < griddx[gridnum]; i++){
      for(j=0; j < griddy[gridnum]; j++){
	for(k=0; k < griddz[gridnum]; k++){

	  // calculate cell index
	  cellindex = k*griddx[gridnum]*griddy[gridnum] + j*griddx[gridnum] + i;

	  clex = gridlex[gridnum] + (( (double) i ) * (gridrex[gridnum] - gridlex[gridnum])/
				     ((double) griddx[gridnum]));

	  crex = gridlex[gridnum] + (( (double) (i+1) ) * (gridrex[gridnum] - gridlex[gridnum])/
				     ((double) griddx[gridnum]));

	  cley = gridley[gridnum] + (( (double) j ) * (gridrey[gridnum] - gridley[gridnum])/
				     ((double) griddy[gridnum]));

	  crey = gridley[gridnum] + (( (double) (j+1) ) * (gridrey[gridnum] - gridley[gridnum])/
				     ((double) griddy[gridnum]));

	  clez = gridlez[gridnum] + (( (double) k ) * (gridrez[gridnum] - gridlez[gridnum])/
				     ((double) griddz[gridnum]));

	  crez = gridlez[gridnum] + (( (double) (k+1) ) * (gridrez[gridnum] - gridlez[gridnum])/
				     ((double) griddz[gridnum]));

	  /*  grid[counter] (ie, grid on level L+1) can either only partially fill the 
	     cell in question on grid[gridnum] (grid on level L) or can totally encompass
	     it (much less likely - it would have to be a 2x2x2 grid!) - therefore, see if
	     grid[counter] overlaps grid[gridnum] in ANY way, and flag the cell if it hasn't been
	     flagged before! */

	  if( (gridlex[counter] < crex) && (gridrex[counter] > clex) &&
	      (gridley[counter] < crey) && (gridrey[counter] > cley) &&
	      (gridlez[counter] < crez) && (gridrez[counter] > clez)
	      )
	    // flag cell if it hasn't yet been flagged, and increment counter
	    if(flagbuffer[cellindex] != -1){
	      flagbuffer[cellindex]=-1;
	      flagged_cells_this_grid++;
	    }
	  
	  /* we're going with the approach where, if ANY part of the grid cell overlaps
	     the user-defined projection boundaries, we want to include it - that way
	     there aren't any funky edge effects.  
	     but, only flag if the cell hasn't been flagged already (so as to not screw
	     up the flagged_cells counter)
	  */
	  if( (crex <= xstart) || (crey <= ystart) || (crez <= zstart) ||
	      (clex >= xend) || (cley >= yend) || (clez >= zend) ) 
	    if(flagbuffer[cellindex] != -1){
	      flagbuffer[cellindex]=-1;
	      flagged_cells_this_grid++;
	    }
	  
	}  // for(k=...
      } // for(j=0...
    }  // for(i=0...

    if(VERBOSEDEBUG) fprintf(stderr,"FlagGridCells:  In Grid: %i (level %i) Grid Looked At: %i (level %i) cells flagged this grid (ttl): %i\n",
		      gridnum,gridlevel[gridnum],counter,gridlevel[counter],flagged_cells_this_grid);
    if(VERBOSEDEBUG) fprintf(stderr,"FlagGridCells: %lf %lf %lf   %lf %lf %lf\n",gridlex[counter],gridley[counter],gridlez[counter],
		      gridrex[counter],gridrey[counter],gridrez[counter]);

  }  // for(counter=0...

  if(VERBOSEDEBUG) fprintf(stderr,"exiting FlagGridCells...\n");
  return SUCCESS;
}


/*-------------------- ProjectToPlane_X ------------------------
 *
 * This routine adds information to projections along the x axis.
 * 
 *--------------------------------------------------------------*/
int ProjectToPlane_X(int gridnum){
  if(DEBUG) fprintf(stderr,"in ProjectToPlane_X\n");

  int xbeg,ybeg,zbeg,xfin,yfin,zfin,i,j,k,xprojeffnumcells,
    gridindex_y,gridindex_z,gridindex_x,gridcellindex,projcellindex;

  double proj_delta,y_phys,z_phys,x_phys;

  double sigma_thompson = 6.65e-25, mh = 1.67e-24, me = 9.11e-28,
    kboltz = 1.38e-16, clight = 3.00e10, csquared = 8.99e20;

  double DensityUnits,
    LengthUnits,
    VelocityUnits,
    DensityConversion,
    XrayConversion,
    TempXrayConversion,
    SZYConversion,SZKinConversion;

  int numcellswritten;

  LengthUnits = (((double) boxsize) * MPC_CM / ((double) hubble) / (1.0 + ((double) redshift)));
  DensityUnits = ((double) omegamatter) * RHOCRIT * pow( (1.0 + ((double) redshift) ), 3.0)
    *((double) hubble) * ( (double) hubble);
  VelocityUnits = ((double) boxsize) * MPC_CM * 
    pow( (4.0 * 3.14 * 6.67 * pow(10.0,-8.0) * ((double) omegamatter) * 
	  RHOCRIT * (1.0 + ((double) initialredshift))),
	 0.5);

  // effective number of cells in x-direction
  xprojeffnumcells = (int) (
			    (xend - xstart) * 
			    ((double) rootgridsize)*
			     pow(2.0,((double) projlevel))
			    );
  if(xprojeffnumcells < 1) xprojeffnumcells = 1;

  /* calculate the beginning and ending cells in all three dimensions,
     for the projection arrays.  We know that:
     (a) grid resolution is <= projection resolution, and the difference
         is always a factor of two to some integar power.
     (b) at the maximum grid resolution allowed by the projection, the 
         edges of the grid cells match up to the projection cells.

     We use the "effective number of cells" (calculated below)
     for the x-axis because we want to weight it appropriately.
  */
  xbeg=(int) ( (( max(xstart,gridlex[gridnum]) - xstart ) / (xend-xstart) ) * 
	       (((double) xprojeffnumcells)) );
  ybeg=(int) ( (( max(ystart,gridley[gridnum]) - ystart ) / (yend-ystart) ) * 
	       (((double) yprojnumcells)) );
  zbeg=(int) ( (( max(zstart,gridlez[gridnum]) - zstart ) / (zend-zstart) ) * 
	       (((double) zprojnumcells)) );

  xfin=(int) ( (( min(xend,gridrex[gridnum]) - xstart) / (xend-xstart) ) * 
	       (((double) xprojeffnumcells)-1.0) );
  yfin=(int) ( (( min(yend,gridrey[gridnum]) - ystart) / (yend-ystart) ) * 
	       (((double) yprojnumcells)-1.0) );
  zfin=(int) ( (( min(zend,gridrez[gridnum]) - zstart) / (zend-zstart) ) * 
	       (((double) zprojnumcells)-1.0) );

  // grid spacing of the projections
  proj_delta = (yend - ystart) / ((double) yprojnumcells);

  // density conversion
  DensityConversion = proj_delta * DensityUnits * LengthUnits / 
    MSOLAR_G * MPC_CM * MPC_CM;

  // xray luminosity conversion
  XrayConversion = proj_delta *  pow(10.0,-20.0) * DensityUnits * 
    DensityUnits / MSOLAR_G / MSOLAR_G * pow( MPC_CM,6.0);

  // temperature conversion
  TempXrayConversion = XrayConversion;

  // sz y effect conversion

  SZYConversion = DensityUnits *0.88/mh
                       *kboltz/(me*csquared)*sigma_thompson
		       *LengthUnits*proj_delta;

  SZKinConversion =  VelocityUnits*DensityUnits*0.88
                       *sigma_thompson/mh/clight
		       *LengthUnits*proj_delta;

  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_X: %lf %i %i %i %i %i %i %i\n",
	  proj_delta,xbeg,ybeg,zbeg,xfin,yfin,zfin,gridnum);
  if(VERBOSEDEBUG) fprintf(stderr,"ProjectTOPlane_X: %lf %lf %lf %lf %lf %lf %i %i %i\n",
	  xstart,ystart,zstart,
	  gridlex[gridnum],gridley[gridnum],gridlez[gridnum],
	  xprojeffnumcells,yprojnumcells,zprojnumcells);

  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_X: conversions: %e %e %e %e %e\n",
	  ((float) DensityConversion),((float) XrayConversion),((float) TempXrayConversion),
	  ((float) SZYConversion),((float) SZKinConversion));
  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_X: units:  %e %e %e %e %e\n",((float) LengthUnits),
	  ((float) DensityUnits),((float) VelocityUnits),MPC_CM,MSOLAR_G);

  if( (xbeg > xfin) || (ybeg > yfin) || (zbeg > zfin) ){
    fprintf(stderr,"error in ProjectToPlane_X: %i %i %i %i %i %i.  Exiting...\n",
		   xbeg,ybeg,zbeg,xfin,yfin,zfin);
    exit(-3);
  }

  numcellswritten=0;

  for(j=ybeg; j <= yfin; j++){
    for(k=zbeg; k <= zfin; k++){

      // calculate 'physical' center of projection array cell in y,z space
      y_phys = ( ((double) j) + 0.5 ) * proj_delta + ystart;
      z_phys = ( ((double) k) + 0.5 ) * proj_delta + zstart;

      // from that, calculate the y,z indices of the grid!
      gridindex_y= (int) (( (y_phys - gridley[gridnum])/(gridrey[gridnum]-gridley[gridnum]) ) * 
			  ((double) griddy[gridnum]) ); // removed -1
      gridindex_z= (int) (( (z_phys - gridlez[gridnum])/(gridrez[gridnum]-gridlez[gridnum]) ) * 
			  ((double) griddz[gridnum]) );  // removed -1

      if(gridindex_z==griddz[gridnum]){
	fprintf(stderr,"ProjectToPlane_X:  array bounds issue %i %i\n",
		gridindex_z,griddz[gridnum]);
	fprintf(stderr,"ProjectToPlane_X:  %i %f %f %f %i %f\n",
		k,proj_delta,zstart,z_phys,zfin,( ((double) k) + 0.5 ));
      }

      projcellindex = k * yprojnumcells + j;
      
      for(i=xbeg; i <= xfin; i++){

	// calculate 'physical' center of x projection array cell
	x_phys = ( ((double) i) + 0.5 ) * proj_delta + xstart;
	
	// calculate which grid cell that is in along the x axis
	gridindex_x= (int) (( (x_phys - gridlex[gridnum])/
			       (gridrex[gridnum]-gridlex[gridnum]) ) 
			    * ((double) griddx[gridnum]) );  // removed -1
	
	gridcellindex = gridindex_z*griddx[gridnum]*griddy[gridnum] + 
	  gridindex_y*griddx[gridnum] + gridindex_x;

	// increment projection buffer values
	if(flagbuffer[gridcellindex] != -1){

	  // baryon density
	  projbuff_density[projcellindex] += densitybuff[gridcellindex] 
	    * ((float) DensityConversion);

	  // dm density
	  //if(dmonthisgrid==1)
	  //projbuff_dmdensity[projcellindex] += dmdensbuff[gridcellindex]
	  //  * ((float) DensityConversion);
	 
	  // xray luminosity
	  projbuff_xraylum[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), 0.5))
	    * ((float) XrayConversion);

	  // xray weighted temperature
	  projbuff_tempxray[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), 1.5))
	    * ((float) TempXrayConversion);

	  // spectral emissivity weighted temperature
	  // this is taken from Rasia et al. 2004, ApJL, 618, 1-4.
	  // Note that we don't use the volume weighting because all cells have the same
	  // size.  Woot.
	  projbuff_spectemp[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), 0.25));
	  projbuff_spectemp_weight[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), -0.75));

	  // mass-weighted temperature.  Note that we don't multiply by dV because
	  // all cells have the same size.
	  projbuff_masstemp[projcellindex] += densitybuff[gridcellindex]*tempbuff[gridcellindex];
	  projbuff_masstemp_weight[projcellindex] += densitybuff[gridcellindex];

	  // x-ray emissivity calculated with the MEKAL tables.
#ifdef USE_MEKAL
	  projbuff_xrayem[projcellindex] += 
	    ComputeSpectralEmissivity(densitybuff[gridcellindex],tempbuff[gridcellindex]);	  
#endif // USE_MEKAL

	  // sz y effect
	  projbuff_szy[projcellindex] += densitybuff[gridcellindex] * tempbuff[gridcellindex]
	    * ((float) SZYConversion);

	  // sz kinematic effect
	  projbuff_szkin[projcellindex] += densitybuff[gridcellindex] * velocitybuff[gridcellindex]
	    *((float) SZKinConversion);

	  // grid level
	  if(projbuff_level[projcellindex] < gridlevel[gridnum]) 
	    projbuff_level[projcellindex] = gridlevel[gridnum];

	  // metal density
	  if(metal_flag==1) projbuff_metal[projcellindex] += metalbuff[gridcellindex]
			      * ((float) DensityConversion);

	  // z1 field density
	  if(z1_flag==1) projbuff_z1[projcellindex] += z1buff[gridcellindex]
			      * ((float) DensityConversion);

	  // z2 field density
	  if(z2_flag==1) projbuff_z2[projcellindex] += z2buff[gridcellindex]
			      * ((float) DensityConversion);
	  
	  if( (multispecies==1) || (multispecies==2) ){
	    // HI density
	    projbuff_HIdensity[projcellindex] += HIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HII density
	    projbuff_HIIdensity[projcellindex] += HIIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HeI density
	    projbuff_HeIdensity[projcellindex] += HeIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HeII density
	    projbuff_HeIIdensity[projcellindex] += HeIIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HeIII density
	    projbuff_HeIIIdensity[projcellindex] += HeIIIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // electron density
	    projbuff_electrondensity[projcellindex] += elecbuff[gridcellindex] 
	      * ((float) DensityConversion);
	  }

	  if(multispecies==2){
	    // HM density
	    projbuff_HMdensity[projcellindex] += HMbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // H2I density
	    projbuff_H2Idensity[projcellindex] += H2Ibuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // H2II density
	    projbuff_H2IIdensity[projcellindex] += H2IIbuff[gridcellindex] 
	      * ((float) DensityConversion);
	  }

	  numcellswritten++;

	}

      } // for(i=xbeg...
            
    }  // for(k=zbeg...
  } // for(j=ybeg...

  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_X: %i cells written\n",numcellswritten);
  if(DEBUG) fprintf(stderr,"leaving ProjectToPlane_X\n");
  return SUCCESS;
}


/*-------------------- ProjectToPlane_Y ------------------------
 *
 * This routine adds information to projections along the y axis.
 * 
 *--------------------------------------------------------------*/
int ProjectToPlane_Y(int gridnum){
  if(DEBUG) fprintf(stderr,"in ProjectToPlane_Y\n");

  int xbeg,ybeg,zbeg,xfin,yfin,zfin,i,j,k,yprojeffnumcells,
    gridindex_y,gridindex_z,gridindex_x,gridcellindex,projcellindex;

  double proj_delta,y_phys,z_phys,x_phys;

  double sigma_thompson = 6.65e-25, mh = 1.67e-24, me = 9.11e-28,
    kboltz = 1.38e-16, clight = 3.00e10, csquared = 8.99e20;

  double DensityUnits,
    LengthUnits,
    VelocityUnits,
    DensityConversion,
    XrayConversion,
    TempXrayConversion,
    SZYConversion,SZKinConversion;

  int numcellswritten;

  LengthUnits = (((double) boxsize) * MPC_CM / ((double) hubble) / (1.0 + ((double) redshift)));
  DensityUnits = ((double) omegamatter) * RHOCRIT * pow( (1.0 + ((double) redshift) ), 3.0)
    *((double) hubble) * ( (double) hubble);
  VelocityUnits = ((double) boxsize) * MPC_CM * 
    pow( (4.0 * 3.14 * 6.67 * pow(10.0,-8.0) * ((double) omegamatter) * 
	  RHOCRIT * (1.0 + ((double) initialredshift))),
	 0.5);

  // effective number of cells in y-direction
  yprojeffnumcells = (int) (
			    (yend - ystart) * 
			    ((double) rootgridsize)*
			     pow(2.0,((double) projlevel))
			    );
  if(yprojeffnumcells < 1) yprojeffnumcells = 1;

  /* calculate the beginning and ending cells in all three dimensions,
     for the projection arrays.  We know that:
     (a) grid resolution is <= projection resolution, and the difference
         is always a factor of two to some integar power.
     (b) at the maximum grid resolution allowed by the projection, the 
         edges of the grid cells match up to the projection cells.

     We use the "effective number of cells" (calculated below)
     for the y-axis because we want to weight it appropriately.
  */
  xbeg=(int) ( (( max(xstart,gridlex[gridnum]) - xstart ) / (xend-xstart) ) * 
	       (((double) xprojnumcells)) );
  ybeg=(int) ( (( max(ystart,gridley[gridnum]) - ystart ) / (yend-ystart) ) * 
	       (((double) yprojeffnumcells)) );
  zbeg=(int) ( (( max(zstart,gridlez[gridnum]) - zstart ) / (zend-zstart) ) * 
	       (((double) zprojnumcells)) );

  xfin=(int) ( (( min(xend,gridrex[gridnum]) - xstart) / (xend-xstart) ) * 
	       (((double) xprojnumcells)-1.0) );
  yfin=(int) ( (( min(yend,gridrey[gridnum]) - ystart) / (yend-ystart) ) * 
	       (((double) yprojeffnumcells)-1.0) );
  zfin=(int) ( (( min(zend,gridrez[gridnum]) - zstart) / (zend-zstart) ) * 
	       (((double) zprojnumcells)-1.0) );

  // grid spacing of the projections
  proj_delta = (xend - xstart) / ((double) xprojnumcells);

  // density conversion
  DensityConversion = proj_delta * DensityUnits * LengthUnits / 
    MSOLAR_G * MPC_CM * MPC_CM;

  // xray luminosity conversion
  XrayConversion = proj_delta *  pow(10.0,-20.0) * DensityUnits * 
    DensityUnits / MSOLAR_G / MSOLAR_G * pow( MPC_CM,6.0);

  // temperature conversion
  TempXrayConversion = XrayConversion;

  // sz y effect conversion

  SZYConversion = DensityUnits *0.88/mh
                       *kboltz/(me*csquared)*sigma_thompson
		       *LengthUnits* proj_delta;

  SZKinConversion =  VelocityUnits*DensityUnits*0.88
                       *sigma_thompson/mh/clight
		       *LengthUnits*proj_delta;

  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Y: %lf %i %i %i %i %i %i %i\n",
	  proj_delta,xbeg,ybeg,zbeg,xfin,yfin,zfin,gridnum);
  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Y: %lf %lf %lf %lf %lf %lf %i %i %i\n",
	  xstart,ystart,zstart,
	  gridlex[gridnum],gridley[gridnum],gridlez[gridnum],
	  xprojnumcells,yprojeffnumcells,zprojnumcells);

  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Y: conversions: %e %e %e %e %e\n",
	  ((float) DensityConversion),((float) XrayConversion),((float) TempXrayConversion),
	  ((float) SZYConversion),((float) SZKinConversion));
  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Y: units:  %e %e %e %e %e\n",((float) LengthUnits),
	  ((float) DensityUnits),((float) VelocityUnits),MPC_CM,MSOLAR_G);

  if( (xbeg > xfin) || (ybeg > yfin) || (zbeg > zfin) ){
    fprintf(stderr,"error in ProjectToPlane_Y: %i %i %i %i %i %i.  Exiting...\n",
		   xbeg,ybeg,zbeg,xfin,yfin,zfin);
    exit(-3);
  }

  numcellswritten=0;

  for(i=xbeg; i <= xfin; i++){
    for(k=zbeg; k <= zfin; k++){

      // calculate 'physical' center of projection array cell in y,z space
      x_phys = ( ((double) i) + 0.5 ) * proj_delta + xstart;
      z_phys = ( ((double) k) + 0.5 ) * proj_delta + zstart;

      // from that, calculate the y,z indices of the grid!
      gridindex_x= (int) (( (x_phys - gridlex[gridnum])/(gridrex[gridnum]-gridlex[gridnum]) ) * 
			  ((double) griddx[gridnum]) ); // removed -1
      gridindex_z= (int) (( (z_phys - gridlez[gridnum])/(gridrez[gridnum]-gridlez[gridnum]) ) * 
			  ((double) griddz[gridnum]) );  // removed -1

      if(gridindex_z==griddz[gridnum]){
	fprintf(stderr,"ProjectToPlane_Y:  array bounds issue %i %i\n",
		gridindex_z,griddz[gridnum]);
	fprintf(stderr,"ProjectToPlane_Y:  %i %f %f %f %i %f\n",
		k,proj_delta,zstart,z_phys,zfin,( ((double) k) + 0.5 ));
      }

      projcellindex = k * xprojnumcells + i;
      
      for(j=ybeg; j <= yfin; j++){

	// calculate 'physical' center of x projection array cell
	y_phys = ( ((double) j) + 0.5 ) * proj_delta + ystart;
	
	// calculate which grid cell that is in along the x axis
	gridindex_y= (int) (( (y_phys - gridley[gridnum])/
			       (gridrey[gridnum]-gridley[gridnum]) ) 
			    * ((double) griddy[gridnum]) );  // removed -1
	
	gridcellindex = gridindex_z*griddx[gridnum]*griddy[gridnum] + 
	  gridindex_y*griddx[gridnum] + gridindex_x;

	// increment projection buffer values
	if(flagbuffer[gridcellindex] != -1){

	  // baryon density
	  projbuff_density[projcellindex] += densitybuff[gridcellindex] 
	    * ((float) DensityConversion);

	  // dm density
	  //if(dmonthisgrid==1)
	  //projbuff_dmdensity[projcellindex] += dmdensbuff[gridcellindex]
	  //  * ((float) DensityConversion);
	 
	  // xray luminosity
	  projbuff_xraylum[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), 0.5))
	    * ((float) XrayConversion);

	  // xray weighted temperature
	  projbuff_tempxray[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), 1.5))
	    * ((float) TempXrayConversion);

	  // spectral emissivity weighted temperature
	  // this is taken from Rasia et al. 2004, ApJL, 618, 1-4.
	  // Note that we don't use the volume weighting because all cells have the same
	  // size.  Woot.
	  projbuff_spectemp[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), 0.25));
	  projbuff_spectemp_weight[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), -0.75));

	  // mass-weighted temperature.  Note that we don't multiply by dV because
	  // all cells have the same size.
	  projbuff_masstemp[projcellindex] += densitybuff[gridcellindex]*tempbuff[gridcellindex];
	  projbuff_masstemp_weight[projcellindex] += densitybuff[gridcellindex];

	  // x-ray emissivity calculated with the MEKAL tables.
#ifdef USE_MEKAL
	  projbuff_xrayem[projcellindex] += 
	    ComputeSpectralEmissivity(densitybuff[gridcellindex],tempbuff[gridcellindex]);	  
#endif // USE_MEKAL

	  // sz y effect
	  projbuff_szy[projcellindex] += densitybuff[gridcellindex] * tempbuff[gridcellindex]
	    * ((float) SZYConversion);

	  // sz kinematic effect
	  projbuff_szkin[projcellindex] += densitybuff[gridcellindex] * velocitybuff[gridcellindex]
	    *((float) SZKinConversion);

	  // grid level
	  if(projbuff_level[projcellindex] < gridlevel[gridnum]) 
	    projbuff_level[projcellindex] = gridlevel[gridnum];

	  // metal density
	  if(metal_flag==1) projbuff_metal[projcellindex] += metalbuff[gridcellindex]
			      * ((float) DensityConversion);

	  // z1 field density
	  if(z1_flag==1) projbuff_z1[projcellindex] += z1buff[gridcellindex]
			      * ((float) DensityConversion);

	  // z2 field density
	  if(z2_flag==1) projbuff_z2[projcellindex] += z2buff[gridcellindex]
			      * ((float) DensityConversion);
	  

	  if( (multispecies==1) || (multispecies==2) ){
	    // HI density
	    projbuff_HIdensity[projcellindex] += HIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HII density
	    projbuff_HIIdensity[projcellindex] += HIIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HeI density
	    projbuff_HeIdensity[projcellindex] += HeIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HeII density
	    projbuff_HeIIdensity[projcellindex] += HeIIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HeIII density
	    projbuff_HeIIIdensity[projcellindex] += HeIIIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // electron density
	    projbuff_electrondensity[projcellindex] += elecbuff[gridcellindex] 
	      * ((float) DensityConversion);
	  }

	  if(multispecies==2){
	    // HM density
	    projbuff_HMdensity[projcellindex] += HMbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // H2I density
	    projbuff_H2Idensity[projcellindex] += H2Ibuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // H2II density
	    projbuff_H2IIdensity[projcellindex] += H2IIbuff[gridcellindex] 
	      * ((float) DensityConversion);
	  }


	  numcellswritten++;

	}

      } // for(j=ybeg...
            
    }  // for(k=zbeg...
  } // for(i=xbeg...

  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Y: %i cells written\n",numcellswritten);
  if(DEBUG) fprintf(stderr,"leaving ProjectToPlane_Y\n");
  return SUCCESS;
}


/*-------------------- ProjectToPlane_Z ------------------------
 *
 * This routine adds information to projections along the z axis.
 * 
 *--------------------------------------------------------------*/
int ProjectToPlane_Z(int gridnum){
  if(DEBUG) fprintf(stderr,"in ProjectToPlane_Z\n");

  int xbeg,ybeg,zbeg,xfin,yfin,zfin,i,j,k,zprojeffnumcells,
    gridindex_y,gridindex_z,gridindex_x,gridcellindex,projcellindex;

  double proj_delta,y_phys,z_phys,x_phys;

  double sigma_thompson = 6.65e-25, mh = 1.67e-24, me = 9.11e-28,
    kboltz = 1.38e-16, clight = 3.00e10, csquared = 8.99e20;

  double DensityUnits,
    LengthUnits,
    VelocityUnits,
    DensityConversion,
    XrayConversion,
    TempXrayConversion,
    SZYConversion,SZKinConversion;

  int numcellswritten;

  LengthUnits = (((double) boxsize) * MPC_CM / ((double) hubble) / (1.0 + ((double) redshift)));
  DensityUnits = ((double) omegamatter) * RHOCRIT * pow( (1.0 + ((double) redshift) ), 3.0)
    *((double) hubble) * ( (double) hubble);
  VelocityUnits = ((double) boxsize) * MPC_CM * 
    pow( (4.0 * 3.14 * 6.67 * pow(10.0,-8.0) * ((double) omegamatter) * 
	  RHOCRIT * (1.0 + ((double) initialredshift))),
	 0.5);

  // effective number of cells in z-direction
  zprojeffnumcells = (int) (
			    (zend - zstart) * 
			    ((double) rootgridsize)*
			     pow(2.0,((double) projlevel))
			    );
  if(zprojeffnumcells < 1) zprojeffnumcells = 1;

  /* calculate the beginning and ending cells in all three dimensions,
     for the projection arrays.  We know that:
     (a) grid resolution is <= projection resolution, and the difference
         is always a factor of two to some integar power.
     (b) at the maximum grid resolution allowed by the projection, the 
         edges of the grid cells match up to the projection cells.

     We use the "effective number of cells" (calculated below)
     for the x-axis because we want to weight it appropriately.
  */
  xbeg=(int) ( (( max(xstart,gridlex[gridnum]) - xstart ) / (xend-xstart) ) * 
	       (((double) xprojnumcells)) );
  ybeg=(int) ( (( max(ystart,gridley[gridnum]) - ystart ) / (yend-ystart) ) * 
	       (((double) yprojnumcells)) );
  zbeg=(int) ( (( max(zstart,gridlez[gridnum]) - zstart ) / (zend-zstart) ) * 
	       (((double) zprojeffnumcells)) );

  xfin=(int) ( (( min(xend,gridrex[gridnum]) - xstart) / (xend-xstart) ) * 
	       (((double) xprojnumcells)-1.0) );
  yfin=(int) ( (( min(yend,gridrey[gridnum]) - ystart) / (yend-ystart) ) * 
	       (((double) yprojnumcells)-1.0) );
  zfin=(int) ( (( min(zend,gridrez[gridnum]) - zstart) / (zend-zstart) ) * 
	       (((double) zprojeffnumcells)-1.0) );

  // grid spacing of the projections
  proj_delta = (yend - ystart) / ((double) yprojnumcells);

  // density conversion
  DensityConversion = proj_delta * DensityUnits * LengthUnits / 
    MSOLAR_G * MPC_CM * MPC_CM;

  // xray luminosity conversion
  XrayConversion = proj_delta *  pow(10.0,-20.0) * DensityUnits * 
    DensityUnits / MSOLAR_G / MSOLAR_G * pow( MPC_CM,6.0);

  // temperature conversion
  TempXrayConversion = XrayConversion;

  // sz y effect conversion

  SZYConversion = DensityUnits *0.88/mh
                       *kboltz/(me*csquared)*sigma_thompson
		       *LengthUnits*proj_delta;

  SZKinConversion =  VelocityUnits*DensityUnits*0.88
                       *sigma_thompson/mh/clight
		       *LengthUnits* proj_delta;

  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Z: %lf %i %i %i %i %i %i %i\n",
	  proj_delta,xbeg,ybeg,zbeg,xfin,yfin,zfin,gridnum);
  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Z: %lf %lf %lf %lf %lf %lf %i %i %i\n",
	  xstart,ystart,zstart,
	  gridlex[gridnum],gridley[gridnum],gridlez[gridnum],
	  xprojnumcells,yprojnumcells,zprojeffnumcells);

  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Z: conversions: %e %e %e %e %e\n",
	  ((float) DensityConversion),((float) XrayConversion),((float) TempXrayConversion),
	  ((float) SZYConversion),((float) SZKinConversion));
  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Z: units:  %e %e %e %e %e\n",((float) LengthUnits),
	  ((float) DensityUnits),((float) VelocityUnits),MPC_CM,MSOLAR_G);

  if( (xbeg > xfin) || (ybeg > yfin) || (zbeg > zfin) ){
    fprintf(stderr,"error in ProjectToPlane_Z: %i %i %i %i %i %i.  Exiting...\n",
		   xbeg,ybeg,zbeg,xfin,yfin,zfin);
    exit(-3);
  }

  numcellswritten=0;

  for(j=ybeg; j <= yfin; j++){
    for(i=xbeg; i <= xfin; i++){

      // calculate 'physical' center of projection array cell in y,z space
      y_phys = ( ((double) j) + 0.5 ) * proj_delta + ystart;
      x_phys = ( ((double) i) + 0.5 ) * proj_delta + xstart;

      // from that, calculate the y,z indices of the grid!
      gridindex_y= (int) (( (y_phys - gridley[gridnum])/(gridrey[gridnum]-gridley[gridnum]) ) * 
			  ((double) griddy[gridnum]) ); // removed -1
      gridindex_x= (int) (( (x_phys - gridlex[gridnum])/(gridrex[gridnum]-gridlex[gridnum]) ) * 
			  ((double) griddx[gridnum]) );  // removed -1

      if(gridindex_x==griddx[gridnum]){
	fprintf(stderr,"ProjectToPlane_Z:  array bounds issue %i %i\n",
		gridindex_x,griddx[gridnum]);
	fprintf(stderr,"ProjectToPlane_Z:  %i %f %f %f %i %f\n",
		i,proj_delta,xstart,x_phys,xfin,( ((double) i) + 0.5 ));
      }

      projcellindex = j * xprojnumcells + i;
      
      for(k=zbeg; k <= zfin; k++){

	// calculate 'physical' center of x projection array cell
	z_phys = ( ((double) k) + 0.5 ) * proj_delta + zstart;
	
	// calculate which grid cell that is in along the z axis
	gridindex_z= (int) (( (z_phys - gridlez[gridnum])/
			       (gridrez[gridnum]-gridlez[gridnum]) ) 
			    * ((double) griddz[gridnum]) );  // removed -1
	
	gridcellindex = gridindex_z*griddx[gridnum]*griddy[gridnum] + 
	  gridindex_y*griddx[gridnum] + gridindex_x;

	// increment projection buffer values
	if(flagbuffer[gridcellindex] != -1){


	  // baryon density
	  projbuff_density[projcellindex] += densitybuff[gridcellindex] 
	    * ((float) DensityConversion);

	  // dm density
	  //if(dmonthisgrid==1)
	  //projbuff_dmdensity[projcellindex] += dmdensbuff[gridcellindex]
	  //  * ((float) DensityConversion);
	 
	  // xray luminosity
	  projbuff_xraylum[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), 0.5))
	    * ((float) XrayConversion);

	  // xray weighted temperature
	  projbuff_tempxray[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), 1.5))
	    * ((float) TempXrayConversion);

	  // spectral emissivity weighted temperature
	  // this is taken from Rasia et al. 2004, ApJL, 618, 1-4.
	  // Note that we don't use the volume weighting because all cells have the same
	  // size.  Woot.
	  projbuff_spectemp[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), 0.25));
	  projbuff_spectemp_weight[projcellindex] += densitybuff[gridcellindex]*
	    densitybuff[gridcellindex] * ((float) pow( ((double) tempbuff[gridcellindex]), -0.75));

	  // mass-weighted temperature.  Note that we don't multiply by dV because
	  // all cells have the same size.
	  projbuff_masstemp[projcellindex] += densitybuff[gridcellindex]*tempbuff[gridcellindex];
	  projbuff_masstemp_weight[projcellindex] += densitybuff[gridcellindex];

	  // x-ray emissivity calculated with the MEKAL tables.
#ifdef USE_MEKAL
	  projbuff_xrayem[projcellindex] += 
	    ComputeSpectralEmissivity(densitybuff[gridcellindex],tempbuff[gridcellindex]);	  
#endif // USE_MEKAL

	  // sz y effect
	  projbuff_szy[projcellindex] += densitybuff[gridcellindex] * tempbuff[gridcellindex]
	    * ((float) SZYConversion);

	  // sz kinematic effect
	  projbuff_szkin[projcellindex] += densitybuff[gridcellindex] * velocitybuff[gridcellindex]
	    *((float) SZKinConversion);

	  // grid level
	  if(projbuff_level[projcellindex] < gridlevel[gridnum]) 
	    projbuff_level[projcellindex] = gridlevel[gridnum];

	  // metal density
	  if(metal_flag==1) projbuff_metal[projcellindex] += metalbuff[gridcellindex]
			      * ((float) DensityConversion);

	  // z1 field density
	  if(z1_flag==1) projbuff_z1[projcellindex] += z1buff[gridcellindex]
			      * ((float) DensityConversion);

	  // z2 field density
	  if(z2_flag==1) projbuff_z2[projcellindex] += z1buff[gridcellindex]
			      * ((float) DensityConversion);
	  
	  if( (multispecies==1) || (multispecies==2) ){
	    // HI density
	    projbuff_HIdensity[projcellindex] += HIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HII density
	    projbuff_HIIdensity[projcellindex] += HIIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HeI density
	    projbuff_HeIdensity[projcellindex] += HeIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HeII density
	    projbuff_HeIIdensity[projcellindex] += HeIIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // HeIII density
	    projbuff_HeIIIdensity[projcellindex] += HeIIIbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // electron density
	    projbuff_electrondensity[projcellindex] += elecbuff[gridcellindex] 
	      * ((float) DensityConversion);
	  }

	  if(multispecies==2){
	    // HM density
	    projbuff_HMdensity[projcellindex] += HMbuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // H2I density
	    projbuff_H2Idensity[projcellindex] += H2Ibuff[gridcellindex] 
	      * ((float) DensityConversion);

	    // H2II density
	    projbuff_H2IIdensity[projcellindex] += H2IIbuff[gridcellindex] 
	      * ((float) DensityConversion);
	  }

	  numcellswritten++;
	  
	}

      } // for(i=xbeg...
            
    }  // for(k=zbeg...
  } // for(j=ybeg...

  if(VERBOSEDEBUG) fprintf(stderr,"ProjectToPlane_Z: %i cells written\n",numcellswritten);
  if(DEBUG) fprintf(stderr,"leaving ProjectToPlane_Z\n");
  return SUCCESS;
}
