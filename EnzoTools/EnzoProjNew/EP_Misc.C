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

// global variables originally declared in another file
extern int projaxis,projlevel,  maxlevel,rootgridsize,
  outputfile_newname;
extern double xstart,ystart,zstart,xend,yend,zend;
extern char **gridfilenames, *outputfilename;

// global vars declared here for the first time
int xprojnumcells,yprojnumcells,zprojnumcells,metal_flag=0,z1_flag=0,z2_flag=0;
float *projbuff_density,*projbuff_dmdensity,*projbuff_xraylum,*projbuff_tempxray,
  *projbuff_level,*projbuff_szy,*projbuff_szkin,*projbuff_metal,*projbuff_z1,
  *projbuff_z2,*projbuff_star,*projbuff_HIdensity,*projbuff_HIIdensity,*projbuff_HeIdensity,
  *projbuff_HeIIdensity,*projbuff_HeIIIdensity,*projbuff_electrondensity,
  *projbuff_HMdensity,*projbuff_H2Idensity,*projbuff_H2IIdensity,
  *projbuff_xrayem, *projbuff_spectemp, *projbuff_spectemp_weight,
  *projbuff_masstemp, *projbuff_masstemp_weight;

/*----------------- CreateProjectionArrays   ---------------------
 *
 * Figure out the size of the projection arrays given their bounds
 * and the maximum level, then allocate memory for the arrays and
 * zero them out.
 *
 * Note that I had to play a bit of a trick with the pow() functions.
 * Their outputs are float, but the arguments are double.  This is 
 * due to a strange implementation of pow() in the IBMS where there
 * are problems with 32 bit/64 bit floats.  Really, this should 
 * probably be some sort of macro.
 *
 *---------------------------------------------------------------- */
int CreateProjectionArrays(void){
  int i;

  if(DEBUG) fprintf(stderr,"in CreateProjectionArrays\n");

  if(VERBOSEDEBUG) fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n",
		    xstart,ystart,zstart,xend,yend,zend);
  if(VERBOSEDEBUG) fprintf(stderr,"%i %i\n",maxlevel,projlevel); 
  
  // given projection axis, calculate size of arrays:
  // number of root grid cells * 2^maxlevel per dimension
  // this works even if numrootgridcells is < 1
  if(projaxis==0){  // proj. along x axis

    // number of cells - x axis
    xprojnumcells = 1;

    // number of cells - y axis
    yprojnumcells = (int) 
      ((yend - ystart) * ((double) rootgridsize)* pow(2.0,(double) projlevel));
    if(yprojnumcells < 1) yprojnumcells = 1;

    // number of cells - z axis
    zprojnumcells = (int)
      ((zend - zstart) * ((double) rootgridsize)* pow(2.0,(double) projlevel));
    if(zprojnumcells < 1) zprojnumcells = 1;

    if(DEBUG) fprintf(stderr,"num. cells: %i %i %i total: %i\n",xprojnumcells,
		      yprojnumcells,zprojnumcells,
		      xprojnumcells*yprojnumcells*zprojnumcells);
  } 

  if(projaxis==1){  // proj. along y axis

    // number of cells - x axis
    xprojnumcells = (int) 
      ((xend - xstart) * ((double) rootgridsize)* pow(2.0,(double) projlevel));
    if(xprojnumcells < 1) xprojnumcells = 1;

    // number of cells - y axis
    yprojnumcells = 1;

    // number of cells - z axis
    zprojnumcells = (int)
      ((zend - zstart) * ((double) rootgridsize)* pow(2.0,(double) projlevel));
    if(zprojnumcells < 1) zprojnumcells = 1;

    if(DEBUG) fprintf(stderr,"num. cells: %i %i %i total: %i\n",xprojnumcells,
		      yprojnumcells,zprojnumcells,
		      xprojnumcells*yprojnumcells*zprojnumcells);
  }

  if(projaxis==2){  // proj. along z axis

    // number of cells - x axis
    xprojnumcells = (int) 
      ((xend - xstart) * ((double) rootgridsize)* pow(2.0,(double) projlevel));
    if(xprojnumcells < 1) xprojnumcells = 1;

    // number of cells - y axis
    yprojnumcells = (int) 
      ((yend - ystart) * ((double) rootgridsize)* pow(2.0,(double) projlevel));
    if(yprojnumcells < 1) yprojnumcells = 1;

    // number of cells - z axis
    zprojnumcells = 1;

    if(DEBUG) fprintf(stderr,"num. cells: %i %i %i total: %i\n",xprojnumcells,
		      yprojnumcells,zprojnumcells,
		      xprojnumcells*yprojnumcells*zprojnumcells);
  }

  fprintf(stderr,"CreateProjectionArrays: num. cells: %i %i %i total: %i\n",xprojnumcells,
		    yprojnumcells,zprojnumcells,
		    xprojnumcells*yprojnumcells*zprojnumcells);
  
  // allocate projection buffers as 1D arrays
  projbuff_density = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_dmdensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_xraylum = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_tempxray = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_level = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_szy = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_szkin = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_metal = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_z1 = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_z2 = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_star = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_HIdensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_HIIdensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_HeIdensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_HeIIdensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_HeIIIdensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_electrondensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_HMdensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_H2Idensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_H2IIdensity = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_xrayem = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_spectemp = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_spectemp_weight = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_masstemp = new float[xprojnumcells*yprojnumcells*zprojnumcells];
  projbuff_masstemp_weight = new float[xprojnumcells*yprojnumcells*zprojnumcells]; 

  // zero out projection arrays
  for(i=0; i < xprojnumcells*yprojnumcells*zprojnumcells ; i++){
    projbuff_density[i]=projbuff_dmdensity[i]=projbuff_xraylum[i]=
      projbuff_tempxray[i]=projbuff_level[i]=projbuff_szy[i]=
      projbuff_szkin[i]=projbuff_metal[i]=projbuff_z1[i]=
      projbuff_z2[i]=projbuff_star[i]=projbuff_HIdensity[i]=
      projbuff_HIIdensity[i]=projbuff_HeIdensity[i]=projbuff_HeIIdensity[i]=
      projbuff_HeIIIdensity[i]=projbuff_electrondensity[i]=projbuff_HMdensity[i]=
      projbuff_H2Idensity[i]=projbuff_H2IIdensity[i]= projbuff_xrayem[i] =
      projbuff_spectemp[i] = projbuff_spectemp_weight[i] = 
      projbuff_masstemp[i] = projbuff_masstemp_weight[i] = 
      0.0;
  }

  if(DEBUG) fprintf(stderr,"CreateProjectionArrays: leaving successfully\n");
  return SUCCESS;
}


/*----------------- NormalizeProjectionArrays   ------------------
 *
 * Normalize a bunch of the projection arrays by the values that
 * they are weighted by!
 * (currently when I say 'a bunch' I mean 'one')
 *
 *---------------------------------------------------------------- */
int NormalizeProjectionArrays(void){
  int i;

  for(i=0; i<xprojnumcells*yprojnumcells*zprojnumcells; i++){
    if(projbuff_xraylum[i]>0.0) projbuff_tempxray[i] /= projbuff_xraylum[i];
    if(projbuff_spectemp_weight[i] > 0.0) projbuff_spectemp[i] / projbuff_spectemp_weight[i];
    if(projbuff_masstemp_weight[i] > 0.0) projbuff_masstemp[i] / projbuff_masstemp_weight[i];    
  }

  return SUCCESS;
}


/*----------------- DeleteProjectionArrays   ---------------------
 *
 * Deletes all of the projection arrays created in 
 * CreateProjectionArrays.  There should be one of these for each 
 * projection buffer.
 * 
 *---------------------------------------------------------------- */
int DeleteProjectionArrays(void){

  delete [] projbuff_density;
  delete [] projbuff_dmdensity;
  delete [] projbuff_xraylum;
  delete [] projbuff_tempxray;
  delete [] projbuff_level;
  delete [] projbuff_szy;
  delete [] projbuff_szkin;
  delete [] projbuff_metal;
  delete [] projbuff_z1;
  delete [] projbuff_z2;
  delete [] projbuff_star;
  delete [] projbuff_HIdensity;
  delete [] projbuff_HIIdensity;
  delete [] projbuff_HeIdensity;
  delete [] projbuff_HeIIdensity;
  delete [] projbuff_HeIIIdensity;
  delete [] projbuff_electrondensity;
  delete [] projbuff_HMdensity;
  delete [] projbuff_H2Idensity;
  delete [] projbuff_H2IIdensity;
  delete [] projbuff_spectemp;
  delete [] projbuff_spectemp_weight;
  delete [] projbuff_masstemp;
  delete [] projbuff_masstemp_weight;
  delete [] projbuff_xrayem;

  return SUCCESS;
}


/*----------------- WriteProjectionArrays   ---------------------
 *
 * Writes the projection information into an hdf5 file, named
 * enzo.project .  This file can then be analyzed using your 2D
 * data visualization tool of choice.
 * 
 *---------------------------------------------------------------- */
int WriteProjectionArrays(void){

  char *outputfile = "enzo.project";

  if(DEBUG) fprintf(stderr,"in WriteProjectionArrays\n");

  // hdf 5 declarations - only one set necessary
  hid_t outputfile_id;
  hsize_t outputdims[2];
  herr_t h5_status;
  herr_t h5_error = -1;

  // hdf 5 declarations - need one per dataset
  hid_t gas_dset_id,gas_dsp_id,  // baryon gas density
    dm_dset_id,dm_dsp_id,  // dark matter density
    xlum_dset_id,xlum_dsp_id,  // xray luminosity
    txray_dset_id,txray_dsp_id,  // xray lum weighted temperature
    spectemp_dset_id, spectemp_dsp_id, // spectral-weighted temperature
    masstemp_dset_id, masstemp_dsp_id, // mass-weighted temperature
    xrayem_dset_id, xrayem_dsp_id, // MEKAL-weighted x-ray emissivity
    level_dset_id,level_dsp_id,  // level
    szy_dset_id,szy_dsp_id, // sz y effect
    szkin_dset_id,szkin_dsp_id,  // kinetic sz effect
    metal_dset_id,metal_dsp_id,  // metal density
    z1_dset_id,z1_dsp_id,  // z1 metal density
    z2_dset_id,z2_dsp_id,  // z2 metal density
    star_dset_id,star_dsp_id,  // star particle metal density
    HI_dset_id,HI_dsp_id,  // neutral hydrogen density
    HII_dset_id,HII_dsp_id,  // ionized hydrogen density
    HeI_dset_id,HeI_dsp_id,  // neutral helium
    HeII_dset_id,HeII_dsp_id,  //  singly ionized helium
    HeIII_dset_id,HeIII_dsp_id,  // doubly ionized helium
    e_dset_id,e_dsp_id,  //  electron density
    HM_dset_id,HM_dsp_id,  //  negatively charged hydrogen
    H2I_dset_id,H2I_dsp_id, //  neutral molecular hydrogen
    H2II_dset_id,H2II_dsp_id; //  singly ionized molecular hydrogen

  // create data space for dataset (depending on axes in question)
  if(projaxis==0){  // x projection
    outputdims[0]= (hsize_t) yprojnumcells;
    outputdims[1]= (hsize_t) zprojnumcells;
  }

  if(projaxis==1){  // y projection
    outputdims[0]= (hsize_t) xprojnumcells;
    outputdims[1]= (hsize_t) zprojnumcells;
  }

  if(projaxis==2){  // z projection
    outputdims[0]= (hsize_t) xprojnumcells;
    outputdims[1]= (hsize_t) yprojnumcells;
  }

  if(DEBUG) fprintf(stderr,"output dimensions are:  %d %d\n",
		    (int) outputdims[0], (int) outputdims[1]);

  // -------------- open output file

  if(outputfile_newname==0){
  outputfile_id = H5Fcreate(outputfile, H5F_ACC_TRUNC, 
			    H5P_DEFAULT, H5P_DEFAULT);
  printf("Output file is named %s\n",outputfile);

  } else {
  outputfile_id = H5Fcreate(outputfilename, H5F_ACC_TRUNC, 
			    H5P_DEFAULT, H5P_DEFAULT);
  printf("Output file is named %s\n",outputfilename);

  }

  if(DEBUG) fprintf(stderr,"status for outputfile_id is %d\n",(int) outputfile_id);
  assert( outputfile_id != h5_error);


  // -------------- create dataspaces (1 per dataset)
  gas_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for gas_dsp_id is %d\n",(int) gas_dsp_id);
  assert( gas_dsp_id != h5_error);

  dm_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for dm_dsp_id is %d\n",(int) dm_dsp_id);
  assert( dm_dsp_id != h5_error);

  xlum_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for xlum_dsp_id is %d\n",(int) xlum_dsp_id);
  assert( xlum_dsp_id != h5_error);

  txray_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for txray_dsp_id is %d\n",(int) txray_dsp_id);
  assert( txray_dsp_id != h5_error);

  spectemp_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for spectemp_dsp_id is %d\n",(int) spectemp_dsp_id);
  assert( spectemp_dsp_id != h5_error);

  masstemp_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for masstemp_dsp_id is %d\n",(int) masstemp_dsp_id);
  assert( masstemp_dsp_id != h5_error);
 
  xrayem_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for xrayem_dsp_id is %d\n",(int) xrayem_dsp_id);
  assert( xrayem_dsp_id != h5_error);

  level_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for level_dsp_id is %d\n",(int) level_dsp_id);
  assert( level_dsp_id != h5_error);

  szy_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for szy_dsp_id is %d\n",(int) szy_dsp_id);
  assert( szy_dsp_id != h5_error);

  szkin_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for szkin_dsp_id is %d\n",(int) szkin_dsp_id);
  assert( szkin_dsp_id != h5_error);

  metal_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for metal_dsp_id is %d\n",(int) metal_dsp_id);
  assert( metal_dsp_id != h5_error);

  z1_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for z1_dsp_id is %d\n",(int) z1_dsp_id);
  assert( z1_dsp_id != h5_error);

  z2_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for z2_dsp_id is %d\n",(int) z2_dsp_id);
  assert( z2_dsp_id != h5_error);

  star_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for star_dsp_id is %d\n",(int) star_dsp_id);
  assert( star_dsp_id != h5_error);

  HI_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for HI_dsp_id is %d\n",(int) HI_dsp_id);
  assert( HI_dsp_id != h5_error);

  HII_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for HII_dsp_id is %d\n",(int) HII_dsp_id);
  assert( HII_dsp_id != h5_error);

  HeI_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for HeI_dsp_id is %d\n",(int) HeI_dsp_id);
  assert( HeI_dsp_id != h5_error);

  HeII_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for HeII_dsp_id is %d\n",(int) HeII_dsp_id);
  assert( HeII_dsp_id != h5_error);

  HeIII_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for HeIII_dsp_id is %d\n",(int) HeIII_dsp_id);
  assert( HeIII_dsp_id != h5_error);

  e_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for e_dsp_id is %d\n",(int) e_dsp_id);
  assert( e_dsp_id != h5_error);

  HM_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for HM_dsp_id is %d\n",(int) HM_dsp_id);
  assert( HM_dsp_id != h5_error);

  H2I_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for H2I_dsp_id is %d\n",(int) H2I_dsp_id);
  assert( H2I_dsp_id != h5_error);

  H2II_dsp_id=H5Screate_simple(2, outputdims, NULL);
  if(DEBUG) fprintf(stderr,"status for H2II_dsp_id is %d\n",(int) H2II_dsp_id);
  assert( H2II_dsp_id != h5_error);


  // -------------- create datasets (1 per dataset)
  gas_dset_id=H5Dcreate(outputfile_id, "/Density", HDF5_FILE_I4, gas_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for gas_dset_id is: %d\n",(int) gas_dset_id);
  assert(gas_dset_id != h5_error);

  dm_dset_id=H5Dcreate(outputfile_id, "/Dark_Matter_Density", HDF5_FILE_I4, dm_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for dm_dset_id is: %d\n",(int) dm_dset_id);
  assert(dm_dset_id != h5_error);
  
  xlum_dset_id=H5Dcreate(outputfile_id, "/X_Ray_Luminosity", HDF5_FILE_I4, xlum_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for xlum_dset_id is: %d\n",(int) xlum_dset_id);
  assert(xlum_dset_id != h5_error);
  
  txray_dset_id=H5Dcreate(outputfile_id, "/Temp_X_Ray_Weighted", HDF5_FILE_I4, txray_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for txray_dset_id is: %d\n",(int) txray_dset_id);
  assert(txray_dset_id != h5_error);

  spectemp_dset_id=H5Dcreate(outputfile_id, "/Temp_Spectral_Weighted", HDF5_FILE_I4, spectemp_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for spectemp_dset_id is: %d\n",(int) spectemp_dset_id);
  assert(spectemp_dset_id != h5_error);

  masstemp_dset_id=H5Dcreate(outputfile_id, "/Temp_Mass_Weighted", HDF5_FILE_I4, masstemp_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for masstemp_dset_id is: %d\n",(int) masstemp_dset_id);
  assert(masstemp_dset_id != h5_error);

  xrayem_dset_id=H5Dcreate(outputfile_id, "/Xray_Emissivity_MEKAL", HDF5_FILE_I4, xrayem_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for xrayem_dset_id is: %d\n",(int) xrayem_dset_id);
  assert(xrayem_dset_id != h5_error);
  
  level_dset_id=H5Dcreate(outputfile_id, "/Level", HDF5_FILE_I4, level_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for level_dset_id is: %d\n",(int) level_dset_id);
  assert(level_dset_id != h5_error);
  
  szy_dset_id=H5Dcreate(outputfile_id, "/SZ_Y_Effect", HDF5_FILE_I4, szy_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for szy_dset_id is: %d\n",(int) szy_dset_id);
  assert(szy_dset_id != h5_error);
  
  szkin_dset_id=H5Dcreate(outputfile_id, "/SZ_Kinetic", HDF5_FILE_I4, szkin_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for szkin_dset_id is: %d\n",(int) szkin_dset_id);
  assert(szkin_dset_id != h5_error);
  
  metal_dset_id=H5Dcreate(outputfile_id, "/Metal_Density", HDF5_FILE_I4, metal_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for metal_dset_id is: %d\n",(int) metal_dset_id);
  assert(metal_dset_id != h5_error);
  
  z1_dset_id=H5Dcreate(outputfile_id, "/Z1_Metal_Density", HDF5_FILE_I4, z1_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for z1_dset_id is: %d\n",(int) z1_dset_id);
  assert(z1_dset_id != h5_error);
  
  z2_dset_id=H5Dcreate(outputfile_id, "/Z2_Metal_Density", HDF5_FILE_I4, z2_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for z2_dset_id is: %d\n",(int) z2_dset_id);
  assert(z2_dset_id != h5_error);
  
  star_dset_id=H5Dcreate(outputfile_id, "/Star_Density", HDF5_FILE_I4, star_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for star_dset_id is: %d\n",(int) star_dset_id);
  assert(star_dset_id != h5_error);

  HI_dset_id=H5Dcreate(outputfile_id, "/HI_Density", HDF5_FILE_I4, HI_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for HI_dset_id is: %d\n",(int) HI_dset_id);
  assert(HI_dset_id != h5_error);

  HII_dset_id=H5Dcreate(outputfile_id, "/HII_Density", HDF5_FILE_I4, HII_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for HII_dset_id is: %d\n",(int) HII_dset_id);
  assert(HII_dset_id != h5_error);

  HeI_dset_id=H5Dcreate(outputfile_id, "/HeI_Density", HDF5_FILE_I4, HeI_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for HeI_dset_id is: %d\n",(int) HeI_dset_id);
  assert(HeI_dset_id != h5_error);

  HeII_dset_id=H5Dcreate(outputfile_id, "/HeII_Density", HDF5_FILE_I4, HeII_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for HeII_dset_id is: %d\n",(int) HeII_dset_id);
  assert(HeII_dset_id != h5_error);

  HeIII_dset_id=H5Dcreate(outputfile_id, "/HeIII_Density", HDF5_FILE_I4, HeIII_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for HeIII_dset_id is: %d\n",(int) HeIII_dset_id);
  assert(HeIII_dset_id != h5_error);

  e_dset_id=H5Dcreate(outputfile_id, "/Electron_Density", HDF5_FILE_I4, e_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for e_dset_id is: %d\n",(int) e_dset_id);
  assert(e_dset_id != h5_error);

  HM_dset_id=H5Dcreate(outputfile_id, "/HM_Density", HDF5_FILE_I4, HM_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for HM_dset_id is: %d\n",(int) HM_dset_id);
  assert(HM_dset_id != h5_error);

  H2I_dset_id=H5Dcreate(outputfile_id, "/H2I_Density", HDF5_FILE_I4, H2I_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for H2I_dset_id is: %d\n",(int) H2I_dset_id);
  assert(H2I_dset_id != h5_error);

  H2II_dset_id=H5Dcreate(outputfile_id, "/H2II_Density", HDF5_FILE_I4, H2II_dsp_id, H5P_DEFAULT);
  if(DEBUG) fprintf(stderr,"read status for H2II_dset_id is: %d\n",(int) H2II_dset_id);
  assert(H2II_dset_id != h5_error);

  
  // -------------- write dataset (1 per dataset)
  SwapArray(outputdims[0],outputdims[1],projbuff_density);  // swap arrays first!
  h5_status = H5Dwrite(gas_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_density);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_dmdensity);  // swap arrays first!
  h5_status = H5Dwrite(dm_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_dmdensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_xraylum);  // swap arrays first!
  h5_status = H5Dwrite(xlum_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_xraylum);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_tempxray);  // swap arrays first!
  h5_status = H5Dwrite(txray_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_tempxray);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_spectemp);  // swap arrays first!
  h5_status = H5Dwrite(spectemp_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_spectemp);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_masstemp);  // swap arrays first!
  h5_status = H5Dwrite(masstemp_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_masstemp);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);


  SwapArray(outputdims[0],outputdims[1],projbuff_xrayem);  // swap arrays first!
  h5_status = H5Dwrite(xrayem_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_xrayem);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_level);  // swap arrays first!
  h5_status = H5Dwrite(level_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_level);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_szy);  // swap arrays first!
  h5_status = H5Dwrite(szy_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_szy);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_szkin);  // swap arrays first!
  h5_status = H5Dwrite(szkin_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_szkin);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_metal);  // swap arrays first!
  h5_status = H5Dwrite(metal_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_metal);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_z1);  // swap arrays first!
  h5_status = H5Dwrite(z1_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_z1);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_z2);  // swap arrays first!
  h5_status = H5Dwrite(z1_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_z2);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_star);  // swap arrays first!
  h5_status = H5Dwrite(star_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_star);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_HIdensity);  // swap arrays first!
  h5_status = H5Dwrite(HI_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_HIdensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_HIIdensity);  // swap arrays first!
  h5_status = H5Dwrite(HII_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_HIIdensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_HeIdensity);  // swap arrays first!
  h5_status = H5Dwrite(HeI_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_HeIdensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_HeIIdensity);  // swap arrays first!
  h5_status = H5Dwrite(HeII_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_HeIIdensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_HeIIIdensity);  // swap arrays first!
  h5_status = H5Dwrite(HeIII_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_HeIIIdensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_electrondensity);  // swap arrays first!
  h5_status = H5Dwrite(e_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_electrondensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_HMdensity);  // swap arrays first!
  h5_status = H5Dwrite(HM_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_HMdensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_H2Idensity);  // swap arrays first!
  h5_status = H5Dwrite(H2I_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_H2Idensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  SwapArray(outputdims[0],outputdims[1],projbuff_H2IIdensity);  // swap arrays first!
  h5_status = H5Dwrite(H2II_dset_id, HDF5_FILE_I4,H5S_ALL, H5S_ALL, H5P_DEFAULT,projbuff_H2IIdensity);
  if(DEBUG) fprintf(stderr,"status for H5Dwrite is %d\n",(int) h5_status);
  assert( h5_status != h5_error);


  // -------------- close dataset (1 per dataset)
  h5_status = H5Dclose(gas_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(dm_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(xlum_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(txray_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(spectemp_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(masstemp_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(xrayem_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(level_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(szy_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(szkin_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(metal_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(z1_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(z2_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(star_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(HI_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(HII_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(HeI_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(HeII_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(HeIII_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(e_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(HM_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(H2I_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Dclose(H2II_dset_id);
  if(DEBUG) fprintf(stderr,"status for H5Dclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);


  // -------------- close dataspace (1 per dataset)
  h5_status = H5Sclose(gas_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(dm_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(xlum_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(txray_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(spectemp_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(masstemp_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(xrayem_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(level_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(szy_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(szkin_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(metal_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(z1_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(z2_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(star_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(HI_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(HII_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(HeI_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(HeII_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(HeIII_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(e_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(HM_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(H2I_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);

  h5_status = H5Sclose(H2II_dsp_id);
  if(DEBUG) fprintf(stderr,"status for H5Sclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error);


  // -------------- close file - only one!
  h5_status = H5Fclose(outputfile_id);
  if(DEBUG) fprintf(stderr,"status for H5Fclose is %d\n",(int) h5_status);
  assert( h5_status != h5_error );

  if(DEBUG) fprintf(stderr,"WriteProjectionArrays: exiting\n");
  return SUCCESS;
}


/*-------------------- CheckBoundaryValues -----------------------
 *
 * If the user sets the bounds (ie, if the bounds are not default
 * of [0,0,0] - [1,1,1]) we want to check to make sure that the
 * bounds line up with the edge of the most highly refined level of
 * prjection, because there are assumptions in ProjectToPlane that
 * are dependent upon this!
 *
 *---------------------------------------------------------------- */
int CheckBoundaryValues(void){
  if(DEBUG) fprintf(stderr,"in CheckBoundaryValues %i %i\n",projlevel,rootgridsize);
  double delta_grid,totalcells,xstartorig,ystartorig,zstartorig,
    xendorig,yendorig,zendorig;
  long int tempcellnum;

  if(DEBUG) fprintf(stderr,"%.16lf %.16lf %.16lf\n",xend-xstart,yend-ystart,zend-zstart);

  // if boundaries are default, don't mess around with them.
  if( (xstart==0.0) && (ystart==0.0) && (zstart==0.0) &&
      (xend==1.0) && (yend==1.0) && (zend==1.0) ) return SUCCESS;
  
  /* The bounds were checked in ParseArgs(), so this is superfluous.
     However, one day somebody will probably erase one or the other,
     so it's good to check. */

  // check the (modified) boundaries to make sure that the 
  if( (xstart < 0.0) || (ystart < 0.0) || (zstart < 0.0) ||
      (xend > 1.0) || (yend > 1.0) || (zend > 1.0) ||
      (xstart >= xend) || (ystart >= yend) || (zstart >= zend) ){
    fprintf(stderr,"CheckBoundaryValues:  Error in boundary values. %lf %lf %lf %lf %lf %lf\n",
	    xstart,ystart,zstart,xend,yend,zend);
    return FAILURE;    
  }

  // calculate the grid spacing of the highest level of refinement
  // in the projection
  delta_grid = 1.0 / ( (double) rootgridsize ) / 
     pow( 2.0, ((double) projlevel)  );

  fprintf(stderr,"delta_grid = %e\n",delta_grid);

  // calculate the total number of grid cells in the cube
  // it's float instead of int in order to make the rest of our
  // calculations straightforward
  totalcells = ((double) rootgridsize) *  pow( 2.0, ((double) projlevel) );

  // save original boundary values;
  xstartorig = xstart;
  ystartorig = ystart;
  zstartorig = zstart;
  xendorig = xend;
  yendorig = yend;
  zendorig = zend;
  
  // x axis values
  if(xstart != 0.0){

    // calculate the value of the cell that xstart should be in
    // force it to be an integar!
    tempcellnum = ((long int) (xstart * totalcells));

    // calculate xstart.  since tempcellnum was forced to be an
    // integar, we now can be assured that xstart is exactly at
    // at the boundary of the most highly refined level of
    // projection
    xstart = ((double) tempcellnum ) / totalcells;
  }

  if(xend != 1.0){
    tempcellnum = ((long int) (xend * totalcells));
    xend = ((double) tempcellnum ) / totalcells;
    // to ensure that we have at least one grid cell
    // thickness
    if(xend==xstart) xend += delta_grid;
  }

  // y axis values - see x stuff for explanation
  if(ystart != 0.0){
    tempcellnum = ((long int) (ystart * totalcells));
    ystart = ((double) tempcellnum ) / totalcells;
  }

  if(yend != 1.0){
    tempcellnum = ((long int) (yend * totalcells));
    yend = ((double) tempcellnum ) / totalcells;
    if(yend==ystart) yend += delta_grid;
  }

  // z axis values - see x stuff for explanation
  if(zstart != 0.0){
    tempcellnum = ((long int) (zstart * totalcells));
    zstart = ((double) tempcellnum ) / totalcells;
  }

  if(zend != 1.0){
    tempcellnum = ((long int) (zend * totalcells));
    zend = ((double) tempcellnum ) / totalcells;
    if(zend==zstart) zend += delta_grid;
  }

  // if the value WAS changed, print out so that the user knows!
  if((xstartorig != xstart))
    fprintf(stderr,"CheckBoundaryValues: xstart was %.16lf and is now %.16lf (%e diff)\n",
	    xstartorig,xstart,xstart-xstartorig);

  if((xendorig != xend))
    fprintf(stderr,"CheckBoundaryValues: xend was %.16lf and is now %.16lf (%e diff)\n",
	    xendorig,xend, xend-xstartorig);

  if((ystartorig != ystart))
    fprintf(stderr,"CheckBoundaryValues: ystart was %.16lf and is now %.16lf (%e diff)\n",
	    ystartorig,ystart,ystart-ystartorig);

  if((yendorig != yend))
    fprintf(stderr,"CheckBoundaryValues: yend was %.16lf and is now %.16lf (%e diff)\n",
	    yendorig,yend, yend-ystartorig);

  if((zstartorig != zstart))
    fprintf(stderr,"CheckBoundaryValues: zstart was %.16lf and is now %.16lf (%e diff)\n",
	    zstartorig,zstart,zstart-zstartorig);

  if((zendorig != zend))
    fprintf(stderr,"CheckBoundaryValues: zend was %.16lf and is now %.16lf (%e diff)\n",
	    zendorig,zend, zend-zstartorig);

  if(DEBUG) fprintf(stderr,"numcells(x):  %lf\n",(xend-xstart)/delta_grid);
  if(DEBUG) fprintf(stderr,"numcells(y):  %lf\n",(yend-ystart)/delta_grid);
  if(DEBUG) fprintf(stderr,"numcells(z):  %lf\n",(zend-zstart)/delta_grid);

  // leave
  if(DEBUG) fprintf(stderr,"leaving CheckBoundaryValues\n");
  return SUCCESS;
}


/*-------------------- CheckForMetals ----------------------------
 *
 * This routine checks to see if the various metal fields exist
 * in the simulation in question.  If so, toggle flags on, if not,
 * toggle off.  There's probably a prettier way to do this, but I
 * don't know what it is.  As it is, there are lots of nasty error
 * messages that appear if the datasets don't exist.  Blah.
 *
 * This could be made significantly less kluged together if we start 
 * printing out flags in the parameter file to tell which metal fields 
 * exist - perhaps when we add more metal fields.
 *
 *---------------------------------------------------------------- */
int CheckForMetals(void){
  if(DEBUG) fprintf(stderr,"in CheckForMetals: %s\n",gridfilenames[0]);

  hid_t     file_id;
  herr_t    h5_status;
  herr_t    h5_error = -1;

  hid_t     test_dset_id;

  // open file
  file_id = H5Fopen(gridfilenames[0], H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );

  // attempt to open each dataset - if it DOESN'T fail,
  // then the dataset exists.  Set flag to 1.
  test_dset_id = H5Dopen(file_id, "Metal_Density");
  if(test_dset_id != h5_error){
    metal_flag = 1;
    h5_status = H5Dclose(test_dset_id);
    assert( h5_status != h5_error);
  }

  test_dset_id = H5Dopen(file_id, "Z_Field1");
  if(test_dset_id != h5_error){
    z1_flag = 1;
    h5_status = H5Dclose(test_dset_id);
    assert( h5_status != h5_error);
  }

  test_dset_id = H5Dopen(file_id, "Z_Field2");
  if(test_dset_id != h5_error){
    z2_flag = 1;
    h5_status = H5Dclose(test_dset_id);
    assert( h5_status != h5_error);
  }

  // print out flags
  fprintf(stderr,"CheckForMetals:  Z: %i  Z1: %i  Z2: %i\n",
	  metal_flag,z1_flag,z2_flag);

  // close file!
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  if(DEBUG) fprintf(stderr,"leaving CheckForMetals\n");

  return SUCCESS;
}


int SwapArray(int xcells, int ycells, float *arraytoswap){

#ifdef SWAP_ARRAYS  // if the user wants us to swap arrays (defined in EnzoProj.h)

  int i,j,rowindex,colindex;
  float *swaparray;

  swaparray = new float[xcells*ycells];

  for(i=0; i<xcells; i++)
    for(j=0; j<ycells; j++){

      rowindex = i*ycells + j;
      colindex = j*xcells + i;

      swaparray[rowindex] = arraytoswap[colindex];

    }

  for(i=0; i<xcells*ycells; i++)
    arraytoswap[i] = swaparray[i];


  delete [] swaparray;

#endif // SWAP_ARRAYS

  return SUCCESS;
}

