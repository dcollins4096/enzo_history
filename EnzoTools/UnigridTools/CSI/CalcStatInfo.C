/*-----------------------------------------------------------

calcsi (CalcStatInfo.C)


This program goes through all of the grid files in a given
data dump (from the .hierarchy file) and calculates the
following information:

1.  1D mass- and volume-weighted distribution functions for
    metallicity (complete, Z1, Z2), temperature, entropy,
    and density

2.  2D distribution functions for density vs. temperature,
    entropy vs. temperature, and all metal and pop III metal
    values vs. density and temperature
    
3.  various mean values

All of this information is put into .dat files.

Usage:

    calcsi <data dump name>

Example:

    calcsi RedshiftOutput0011

Brian O'Shea
bwoshea@cosmos.ucsd.edu
3 Feb 2002

-----------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#include "extern_hdf5.h"  // hdf 5 prototypes

#define MAX_LINE_LENGTH 256
#define HDF5_I4 H5T_NATIVE_INT
#define HDF5_FILE_I4 H5T_STD_I32BE

#define DEBUG 1  // 1 on, 0 off
#define VDEBUG 0 // verbose debug - 1 on, 0 off

// these define the bounds of our various arrays
//   (all in log10 space)
#define ZMIN -12.0  // metallicity
#define ZMAX 2.0    

// These numbers are the minimum metallicity used
// for calculating the volume filling factor (VFF)
// for the various metal quantities.
#define Z_VFF_MIN  0.000001
#define Z1_VFF_MIN 0.000001
#define Z2_VFF_MIN 0.000001

#define NUMBEROFBINS 100  // number of bins in the arrays

// these are for the 2D distribution functions
#define SMIN 17.0  // entropy 
#define SMAX 30.0
#define DMIN -4.0
#define DMAX 8.0
#define TMIN -1.5
#define TMAX 9.0

#define PIXELS 300
#define RHOCRIT_CGS  1.8788e-29 // * hubble * hubble
#define OMEGA_BARYON 0.04  // probably should fix THIS at some point

int main(int argc, char *argv[]){

  FILE *hierfile, *outputfile,  *headerfile, //*evr_vol,*evr_mass,*tvr_vol,*tvr_mass,
    *projname;

  // file names and such
  char *hierarchyfilename=NULL,*datafilename=NULL;
  char *line = new char[MAX_LINE_LENGTH];
  char *gfname = new char[MAX_LINE_LENGTH];

  // field names to extract
  char *densfield = "Density";
  char *metal = "Metal_Density";
  char *tracer1 = "Z_Field1";
  char *tracer2 = "Z_Field2";
  char *tempfield = "Temperature";

  float
    *densbuff,*metalbuff,*z1buff,*z2buff,*tempbuff,
    glex,gley,glez,  // grid left edge
    grex,grey,grez,  // grid right edge
    hubble,omegamatter,boxsize, redshift,
    densconstant,
    massconv,
    celllength,
    cellvol,

    metalval[NUMBEROFBINS],    // arrays for metal field
    metalfrac_mass[NUMBEROFBINS],
    metalfrac_vol[NUMBEROFBINS],
    metalfrac_mass_cum[NUMBEROFBINS],
    metalfrac_vol_cum[NUMBEROFBINS],
    metalfrac_mass_total,
    metalfrac_vol_total,

    z1frac_mass[NUMBEROFBINS],    // arrays for tracer field 1
    z1frac_vol[NUMBEROFBINS],
    z1frac_mass_cum[NUMBEROFBINS],
    z1frac_vol_cum[NUMBEROFBINS],
    z1frac_mass_total,
    z1frac_vol_total,

    z2frac_mass[NUMBEROFBINS],     // arrays for tracer field 2
    z2frac_vol[NUMBEROFBINS],
    z2frac_mass_cum[NUMBEROFBINS],
    z2frac_vol_cum[NUMBEROFBINS],
    z2frac_mass_total,
    z2frac_vol_total,

    mfracbin,
    total_volume,
    total_mass,
    mean_metal_mass, mean_metal_vol,
    mean_z1_mass,mean_z1_vol,
    mean_z2_mass,mean_z2_vol,
    metallicity,
    z1_metallicity,
    z2_metallicity,
    totalmetalmass,
    totalgasmass,
    totalz1mass,
    totalz2mass,
    z_min,z_max,
    z1_min,z1_max,
    z2_min,z2_max,

    TVR_MassweightedTable[PIXELS][PIXELS],  // mass vs. density 2d dist'n fctns
    TVR_VolumeweightedTable[PIXELS][PIXELS],
    EVR_MassweightedTable[PIXELS][PIXELS],  // entropy vs. density
    EVR_VolumeweightedTable[PIXELS][PIXELS],

    MVT_MassweightedTable[PIXELS][PIXELS],  // metals vs. temperature
    MVT_VolumeweightedTable[PIXELS][PIXELS],

    MVR_MassweightedTable[PIXELS][PIXELS],  // metals vs. density
    MVR_VolumeweightedTable[PIXELS][PIXELS],

    P3VT_MassweightedTable[PIXELS][PIXELS],  // pop III metals vs. temperature
    P3VT_VolumeweightedTable[PIXELS][PIXELS],

    P3VR_MassweightedTable[PIXELS][PIXELS],  // pop III metals vs. density
    P3VR_VolumeweightedTable[PIXELS][PIXELS],

    maxdens,maxtemp,entropy,
    
    
    mean_entropy_mass,mean_entropy_volume,mean_temp_mass,mean_temp_volume,
    sum_density,sum_density_squared,clumping_factor,
    ent_mass[NUMBEROFBINS],ent_vol[NUMBEROFBINS],
    dens_mass[NUMBEROFBINS],dens_vol[NUMBEROFBINS],
    temp_mass[NUMBEROFBINS],temp_vol[NUMBEROFBINS],
    ent_mass_cum[NUMBEROFBINS],ent_vol_cum[NUMBEROFBINS],
    dens_mass_cum[NUMBEROFBINS],dens_vol_cum[NUMBEROFBINS],
    temp_mass_cum[NUMBEROFBINS],temp_vol_cum[NUMBEROFBINS],
    dens_mass_total,dens_vol_total,temp_mass_total,temp_vol_total,ent_mass_total,ent_vol_total,
    ent[NUMBEROFBINS],temperature[NUMBEROFBINS],density[NUMBEROFBINS],
    tbin1D,dbin1D,sbin1D,tbin2D,dbin2D,sbin2D,mbin2D
    ;


  int
    griddimx,griddimy,griddimz,  // grid xyz dims
    gridsex,gridsey,gridsez,  // grid start indices
    grideex,grideey,grideez,  // grid end indices
    gdx,gdy,gdz,  // grid dims - no buffers
    i,j, ndims,cellcounter,
    metalbin,
    z1bin,
    z2bin,
    vffcells,
    z1_vffcells,
    z2_vffcells,
    dummy,
    tempbin,entbin,densbin
    ;


  // hdf 5 stuff
  hid_t       file_id;
  hid_t       dens_dset_id, metal_dset_id,
    z1_dset_id,z2_dset_id,temp_dset_id;
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
  
  // output file names
  //outfilename="metal_distfctn.dat";
  //outfile1name="meaninfo.dat";
  //outfile2name="baryon_distfctn.dat";

  printf("data dump name:          %s\n",datafilename);
  //printf("dist fctn. file name:    %s\n",outfilename);
  //printf("mean data file name:     %s\n",outfile1name);

  // ------------ set up distribution function arrays 

  // 1D bin size (in log10 space)
  mfracbin=(ZMAX-ZMIN)/NUMBEROFBINS;
  sbin1D = (SMAX-SMIN)/NUMBEROFBINS;
  dbin1D = (DMAX-DMIN)/NUMBEROFBINS;
  tbin1D = (TMAX-TMIN)/NUMBEROFBINS;

  // calculate values in the center of the bin so that
  // the plots look pretty 'n stuff
  for(i=0;i<NUMBEROFBINS;i++){
    metalfrac_mass[i]=metalfrac_vol[i]=0.0;
    z1frac_mass[i]=z1frac_vol[i]=z2frac_mass[i]=z2frac_vol[i]=
      ent_mass[i] = ent_vol[i] = dens_mass[i] = dens_vol[i]
      = temp_mass[i] = temp_vol[i] = 0.0;

    if(i==0){
      metalval[i]=mfracbin/2.0+ZMIN;
      ent[i] = sbin1D/2.+SMIN;
      temperature[i] = tbin1D/2.+TMIN;
      density[i] = dbin1D/2.+DMIN;
    } else{
      metalval[i]=((float) i + 0.5)*mfracbin + ZMIN;
      ent[i] = ((float) i + .5)*sbin1D + SMIN;
      temperature[i] = ((float) i + .5)*tbin1D + TMIN;
      density[i] = ((float) i + .5)*dbin1D + DMIN;

    }
  }

  // zero mean quantities
  total_volume=total_mass=mean_metal_mass=mean_metal_vol=
    mean_z1_mass=mean_z1_vol=mean_z2_mass=mean_z2_vol=
    totalz1mass=totalz2mass=
    totalmetalmass=totalgasmass=
    mean_entropy_mass
    = mean_entropy_volume = mean_temp_mass = mean_temp_volume
    = sum_density = sum_density_squared = 0.0;


  // count cells to make sure we have the right number!
  cellcounter=0;
  
  // keep track of volume filling factor 
  vffcells=z1_vffcells=z2_vffcells=0;

  // calculate bin sizes for 2D distribution functions
  // set # of bins
  sbin2D = (SMAX-SMIN)/PIXELS;
  dbin2D = (DMAX-DMIN)/PIXELS;
  tbin2D = (TMAX-TMIN)/PIXELS;
  mbin2D = (ZMAX-ZMIN)/PIXELS;
  
  // zero out 2D dist fctn tables
  for(i=0; i<PIXELS;i++){
    for(j=0; j<PIXELS;j++){
      EVR_MassweightedTable[i][j]=0.0;
      EVR_VolumeweightedTable[i][j]=0.0;

      TVR_MassweightedTable[i][j]=0.0;
      TVR_VolumeweightedTable[i][j]=0.0;

      MVT_MassweightedTable[i][j]=0.0;
      MVT_VolumeweightedTable[i][j]=0.0;

      MVR_MassweightedTable[i][j]=0.0;
      MVR_VolumeweightedTable[i][j]=0.0;

      P3VT_MassweightedTable[i][j]=0.0;
      P3VT_VolumeweightedTable[i][j]=0.0;

      P3VR_MassweightedTable[i][j]=0.0;
      P3VR_VolumeweightedTable[i][j]=0.0;


    }
  }


  // set initial (out-of-bounds) values for z, z1, z2 min/max
  z_min=z1_min=z2_min=100.0;
  z_max=z1_max=z2_max=-1.0;

  // set initial (out-of-bounds) values for max temp,dens
  maxdens=maxtemp=-1.0;


  // open up data output parameter file, get redshift, omega_matter
  // and hubble constant and box size in Mpc/h
  headerfile = fopen(datafilename,"r");

  while( fgets(line, MAX_LINE_LENGTH, headerfile) != NULL){

    // current redshift
    if(strncmp(line,"CosmologyCurrentRedshift",24)==0)
      sscanf(line,"CosmologyCurrentRedshift = %f",
	     &redshift);

    // comoving box size (in Mpc/h)
    if(strncmp(line,"CosmologyComovingBoxSize",24)==0)
      sscanf(line,"CosmologyComovingBoxSize   = %f",
	     &boxsize);
    // hubble constant (in units of 100 km/s/Mpc)
    if(strncmp(line,"CosmologyHubbleConstantNow",26)==0)
      sscanf(line,"CosmologyHubbleConstantNow = %f",
	     &hubble);

    // comoving box size
    if(strncmp(line,"CosmologyOmegaMatterNow",23)==0)
      sscanf(line,"CosmologyOmegaMatterNow = %f",
	     &omegamatter);

  } // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)

  fclose(headerfile);

  if(DEBUG) printf("\tredshift:  %f\n\thubble:  %f\n\tbox size:\n\t:  %f\n\tomega_matter:  %f\n",
		   redshift,hubble,boxsize,omegamatter);

  // calculate density constant - omit factors of redshift since they fall
  // out when we calculate masses anyway
  //densconstant = omegamatter * 2.78*pow(10.0,11.0) * hubble * hubble;
  densconstant = omegamatter * 2.78 * 100000000000.0 * hubble * hubble;

  if(VDEBUG) printf("density constant is %e\n",densconstant);

  // get hierarchy file name
  hierarchyfilename = datafilename;
  hierarchyfilename=strcat(hierarchyfilename,".hierarchy");
  printf("hierarchy file:          %s\n",hierarchyfilename);
  //printf("data file:  %s\n",datafilename);

  // Open up hierarchy file 
  hierfile = fopen(hierarchyfilename,"r");

  /* read through the hierarchy file.  At every grid, extract all of the
     particle quantities mentioned above.  Then go through the list of creation
     times.  If T_create for a given particle is greater than zero, it's a 
     star particle, so add it to the output file. */
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
      sscanf(line,"GridLeftEdge      = %f %f %f",&glex,&gley,&glez);
      fgets(line, MAX_LINE_LENGTH, hierfile);  // right edge
      sscanf(line,"GridRightEdge     = %f %f %f",&grex,&grey,&grez);
      
      // "junk" lines
      for(i=0;i<5;i++) fgets(line, MAX_LINE_LENGTH, hierfile);

      // get name of grid file
      sscanf(line,"BaryonFileName = %s",gfname);

      if(VDEBUG) printf("%d %d %d %d %d %d %d %d %d\n",
		       griddimx,griddimy,griddimz,
		       gridsex,gridsey,gridsez,
		       grideex,grideey,grideez);
      if(VDEBUG) printf("%f %f %f %f %f %f\n",
		       glex,gley,glez,
		       grex,grey,grez);
      if(VDEBUG) printf("%s\n",gfname);

      // calculate grid size w/out buffers
      gdx = 1+grideex-gridsex;
      gdy = 1+grideey-gridsey;
      gdz = 1+grideez-gridsez;

      // calculate mass conversion factor
      celllength = ((grex-glex) / ((float) gdx) );

      //cellvol = pow( celllength, 3.0);
      cellvol = celllength*celllength*celllength;

      //massconv = densconstant * pow(  (celllength*boxsize/hubble),3.0);
      massconv = densconstant * cellvol * boxsize * boxsize * boxsize 
	/ hubble / hubble / hubble;

      if(VDEBUG) printf("Grid is %d %d %d without buffers\n",gdx,gdy,gdz);
      if(VDEBUG) printf("that comes out to %d total cells\n",gdx*gdy*gdz);
      if(VDEBUG) printf("Mass conversion factor: %f\n",massconv);

      /* NOW OPERATE ON THIS GRID */

      // open grid file
      file_id = H5Fopen(gfname, H5F_ACC_RDWR, H5P_DEFAULT);
      assert( file_id != h5_error );

      dens_dset_id = H5Dopen(file_id,densfield);
      assert( dens_dset_id != h5_error );

      metal_dset_id = H5Dopen(file_id,metal);
      assert( metal_dset_id != h5_error );

      z1_dset_id = H5Dopen(file_id,tracer1);
      assert( z1_dset_id != h5_error);

      z2_dset_id = H5Dopen(file_id,tracer2);
      assert( z2_dset_id != h5_error);

      temp_dset_id = H5Dopen(file_id,tempfield);
      assert( temp_dset_id != h5_error);

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
      if(VDEBUG) printf("Ndims %d\n",ndims);
      for ( i = 0; i < ndims; i++)
	{
	  dims[i] = xdims[i];
	  size = size * dims[i];
	  if(VDEBUG) printf(" Dim %d\n", (int) xdims[i]);
	}
      if(VDEBUG) printf("Size %d\n", (int) size);

      file_dsp_id = H5Screate_simple(ndims, dims, NULL);
      assert( file_dsp_id != h5_error );
      mem_dsp_id = H5Screate_simple(1, &size, NULL);
      assert( mem_dsp_id != h5_error );

      // now read arrays into memory
      if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
	{
	  densbuff = new float[(int) size];
	  metalbuff = new float[(int) size];
	  z1buff = new float[(int) size];
	  z2buff = new float[(int) size];
	  tempbuff = new float[(int) size];

	  mem_type_id = H5T_NATIVE_FLOAT;

	  // read dens field into an array
	  h5_status = H5Dread(dens_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, densbuff);
	  if(VDEBUG) printf("float read status %d for dens field\n", (int) h5_status);
	  assert( h5_status != h5_error );

	  // read metal field into an array
	  h5_status = H5Dread(metal_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, metalbuff);
	  if(VDEBUG) printf("float read status %d for metal field\n", (int) h5_status);
	  assert( h5_status != h5_error );

	  // read z1 field into an array
	  h5_status = H5Dread(z1_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, z1buff);
	  if(VDEBUG) printf("float read status %d for z1 field\n", (int) h5_status);
	  assert( h5_status != h5_error );

	  // read z2 field into an array
	  h5_status = H5Dread(z2_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, z2buff);
	  if(VDEBUG) printf("float read status %d for z2 field\n", (int) h5_status);
	  assert( h5_status != h5_error );

	  // read temperature field into array
	  h5_status = H5Dread(temp_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, tempbuff);
	  if(VDEBUG) printf("float read status %d for temp field\n", (int) h5_status);
	  assert( h5_status != h5_error );

	}

      // calculate distribution functions and weights

      // go through buffers
      for(i=0; i<size; i++){
	// get metallicity in terms of solar metallicity
	metallicity = metalbuff[i]/densbuff[i]/0.02;
	z1_metallicity = z1buff[i]/densbuff[i]/0.02;
	z2_metallicity = z2buff[i]/densbuff[i]/0.02;
	entropy = tempbuff[i] * 
	  (
	   (float) pow( ((double) (densbuff[i]*omegamatter*RHOCRIT_CGS*hubble*hubble)), 
				  ((double) (-2./3.)))
	   / pow( ((double) (1.0+redshift)),2.0)
	   
	   );

	if(VDEBUG) fprintf(stderr,"entropy: %e %e %e %i\n",entropy,tempbuff[i],
			   densbuff[i]*omegamatter*RHOCRIT_CGS*hubble*hubble,i);

	if(i==0 && DEBUG) printf("%e %e\n",z1_metallicity,z2_metallicity);

	// check to make sure that there isn't more metal than there is gas!
	if(metallicity > (1.0/0.02)) fprintf(stderr,"wtf?  metal value > one:  %f %f\n",
				      metalbuff[i],densbuff[i]);

	// keep track of maximum and minimum values of metallicity
	// (for all metal fields)
	if(metallicity > z_max) z_max = metallicity;
	if(metallicity < z_min) z_min = metallicity;

	if(z1_metallicity > z1_max) z1_max = z1_metallicity;
	if(z1_metallicity < z1_min) z1_min = z1_metallicity;

	if(z2_metallicity > z2_max) z2_max = z2_metallicity;
	if(z2_metallicity < z2_min) z2_min = z2_metallicity;

	if(densbuff[i]>maxdens) maxdens=densbuff[i];
	if(tempbuff[i]>maxtemp) maxtemp=tempbuff[i];

	// increment mean quantities
	mean_metal_mass += metallicity * densbuff[i] * cellvol;
	mean_metal_vol += metallicity * cellvol;

	mean_z1_mass += z1_metallicity * densbuff[i] * cellvol;
	mean_z1_vol += z1_metallicity * cellvol;

	mean_z2_mass += z2_metallicity * densbuff[i] * cellvol;
	mean_z2_vol += z2_metallicity * cellvol;

	mean_entropy_mass += entropy*densbuff[i] * cellvol;
	mean_entropy_volume += entropy*cellvol;

	mean_temp_mass += tempbuff[i]*densbuff[i] * cellvol;
	mean_temp_volume += tempbuff[i]*cellvol;

	sum_density += densbuff[i]*cellvol;
	sum_density_squared += densbuff[i]*densbuff[i]*cellvol;

	// increment total "mass" (ie, density) and volume
	total_mass += densbuff[i] * cellvol;
	total_volume += cellvol;

	// calculate which metallicity bin this cell belongs to
	// (all three metal variables)
	metalbin=(int) ( (log10(metallicity) - ZMIN)/mfracbin);
	z1bin=(int) ( (log10(z1_metallicity) - ZMIN)/mfracbin);
	z2bin=(int) ( (log10(z2_metallicity) - ZMIN)/mfracbin);

	// calculate 1d bins for temp,dens,entropy
	tempbin = (int) ( (log10(tempbuff[i]) - TMIN)/tbin1D);
	densbin = (int) ( (log10(densbuff[i]*omegamatter/OMEGA_BARYON) - DMIN)/dbin1D);
	entbin = (int) ( (log10(entropy) - SMIN)/sbin1D);

	// error check - arrays in bounds?
	// it's just a check - if they DON'T fit, we just put it in the lowest
	// bin.
	if(metalbin < 0 || metalbin >= NUMBEROFBINS || 
	   z1bin >= NUMBEROFBINS || z2bin >= NUMBEROFBINS){
	  fprintf(stderr,"wtf? %d %d %d %d %e %e %e\n",metalbin,z1bin,z2bin,NUMBEROFBINS,
		  metallicity,z1_metallicity,z2_metallicity);
	  //exit(10);
	}

	// error check on bins
	if(tempbin<0 || densbin<0 || entbin < 0 || 
	   tempbin >= NUMBEROFBINS || densbin >= NUMBEROFBINS
	   || entbin >= NUMBEROFBINS){
	  fprintf(stderr,"what the heck? %d %d %d %d %s\n",tempbin,densbin,
		  entbin,NUMBEROFBINS,gfname);
	  exit(10);
	}

	// increment volume filling factor # cells
	if(metallicity >= Z_VFF_MIN) vffcells++;
	if(z1_metallicity >= Z1_VFF_MIN) z1_vffcells++;
	if(z2_metallicity >= Z2_VFF_MIN) z2_vffcells++;

	// check to make sure that the metallicity isn't zero anywhere - 
	// though it really doesn't make much of a difference
	if(metallicity==0.0 || z1_metallicity == 0.0 || z2_metallicity == 0.0)
	  fprintf(stderr,"wtf?  metallicity is zero! %f %f %f %d %s\n",
		  metallicity, z1_metallicity, z2_metallicity,
		  i,gfname);

	// add weights to bins
	if(metalbin >=0){
	  metalfrac_mass[metalbin] += densbuff[i]*cellvol;
	  metalfrac_vol[metalbin] += cellvol;
	} else{
	  metalfrac_mass[0] += densbuff[i]*cellvol;
	  metalfrac_vol[0] += cellvol;
	}

	if(z1bin >= 0){
	  z1frac_mass[z1bin] += densbuff[i]*cellvol;
	  z1frac_vol[z1bin] += cellvol;
	} else{
	  z1frac_mass[0] += densbuff[i]*cellvol;
	  z1frac_vol[0] += cellvol;
	}

	if(z2bin >= 0){
	  z2frac_mass[z2bin] += densbuff[i]*cellvol;
	  z2frac_vol[z2bin] += cellvol;
	} else{
	  z2frac_mass[0] += densbuff[i]*cellvol;
	  z2frac_vol[0] += cellvol;
	}

	// increment 1D dist'n fctns
	ent_mass[entbin] += densbuff[i]*cellvol;
	ent_vol[entbin]+= cellvol;
	dens_mass[densbin] += densbuff[i]*cellvol;
	dens_vol[densbin] += cellvol;
	temp_mass[tempbin] += densbuff[i]*cellvol;
	temp_vol[tempbin] += cellvol;
	
	//----------------- 2D distribution function stuff	

	// calculate bins for 2D distribution functions
	entbin = (int) ( (log10(entropy) - SMIN)/sbin2D);
	densbin = (int) ( (log10(densbuff[i]*omegamatter/OMEGA_BARYON) - DMIN)/dbin2D);
	tempbin = (int) ( (log10(tempbuff[i]) - TMIN)/tbin2D);
	metalbin = (int) ( (log10(metallicity) - ZMIN)/mbin2D);
	z1bin = (int) ( (log10(z1_metallicity) - ZMIN)/mbin2D);


	if(z1bin < 0) z1bin=0;
	if(metalbin < 0) metalbin=0;

	// error check on bins
	if(entbin<0 || densbin<0 || entbin >= PIXELS || densbin >= PIXELS || tempbin<0 || tempbin >= PIXELS
	   || metalbin < 0 || z1bin < 0 || metalbin >= PIXELS || z1bin >= PIXELS
	   ){
	  fprintf(stderr,"what the heck? %d %d %d %d %d %i %s\n",entbin,tempbin,
		  densbin,metalbin,z1bin, PIXELS,gfname);
	  exit(10);
	}

	// entropy vs. density, mass and volume weighted
	EVR_MassweightedTable[entbin][densbin]+= densbuff[i]*cellvol;
	EVR_VolumeweightedTable[entbin][densbin]+= cellvol;

	// temp vs density, mass and volume weighted
	TVR_MassweightedTable[tempbin][densbin]+= densbuff[i]*cellvol;
	TVR_VolumeweightedTable[tempbin][densbin]+= cellvol;

	// metallicity vs. density, mass and volume weighted
	MVR_MassweightedTable[metalbin][densbin]+= densbuff[i]*cellvol;
	MVR_VolumeweightedTable[metalbin][densbin]+= cellvol;

	// metallicity vs. temperature, mass and volume weighted
	MVT_MassweightedTable[metalbin][tempbin]+= densbuff[i]*cellvol;
	MVT_VolumeweightedTable[metalbin][tempbin]+= cellvol;

	// pop III metallicity vs. density, mass and volume weighted
	P3VR_MassweightedTable[z1bin][densbin]+= densbuff[i]*cellvol;
	P3VR_VolumeweightedTable[z1bin][densbin]+= cellvol;

	// pop III metallicity vs. temperature, mass and volume weighted
	P3VT_MassweightedTable[z1bin][tempbin]+= densbuff[i]*cellvol;
	P3VT_VolumeweightedTable[z1bin][tempbin]+= cellvol;

	// keep track of totals
	cellcounter++;
	totalmetalmass += metalbuff[i];
	totalgasmass += densbuff[i];


      } // end of for(i=0; i<size; i++){

      // clean up
      delete [] densbuff;
      delete [] metalbuff;
      delete [] z1buff;
      delete [] z2buff;
      delete [] tempbuff;

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

      h5_status = H5Dclose(metal_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(z1_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(z2_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(temp_dset_id);
      assert( h5_status != h5_error);

    }  // end of if(strncmp(line,"Grid = ",7)==0){
  }  // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){

  fclose(hierfile);  // close hierarchy file

  // divide through by weighting quantities to get mean vals
  mean_metal_mass /= total_mass;
  mean_metal_vol /= total_volume;

  mean_z1_mass /= total_mass;
  mean_z1_vol /= total_volume;

  mean_z2_mass /= total_mass;
  mean_z2_vol /= total_volume;

  mean_entropy_mass/=total_mass;
  mean_entropy_volume/=total_volume;
  mean_temp_mass/=total_mass;
  mean_temp_volume/=total_volume;
  sum_density/=total_volume;
  sum_density_squared/=total_volume;

  // calculate clumping factor:  <rho^2> / <rho>^2
  clumping_factor = sum_density_squared / sum_density / sum_density;

  // output mean info and various statistics
  outputfile=fopen("meaninfo.dat","w");

  // metallicity info
  fprintf(outputfile,"Z: mass weighted mean:      %e\n",mean_metal_mass);
  fprintf(outputfile,"Z: volume weighted mean:    %e\n",mean_metal_vol);
  fprintf(outputfile,"Z: max metallicity:         %e\n",z_max);
  fprintf(outputfile,"Z: min metallicity:         %e\n",z_min);
  fprintf(outputfile,"Z: volume filling factor    %e\n\n", 
	 (((float) vffcells)/((float) cellcounter)));

  // Z1 info
  fprintf(outputfile,"Z1: mass weighted mean:      %e\n",mean_z1_mass);
  fprintf(outputfile,"Z1: volume weighted mean:    %e\n",mean_z1_vol);
  fprintf(outputfile,"Z1: max metallicity:         %e\n",z1_max);
  fprintf(outputfile,"Z1: min metallicity:         %e\n",z1_min);
  fprintf(outputfile,"Z1: volume filling factor:   %e\n\n", 
	 (((float) z1_vffcells)/((float) cellcounter)));

  // Z2 info
  fprintf(outputfile,"Z2: mass weighted mean:      %e\n",mean_z2_mass);
  fprintf(outputfile,"Z2: volume weighted mean:    %e\n",mean_z2_vol);
  fprintf(outputfile,"Z2: max metallicity:         %e\n",z2_max);
  fprintf(outputfile,"Z2: min metallicity:         %e\n",z2_min);
  fprintf(outputfile,"Z2: volume filling factor:   %e\n\n",
	 (((float) z2_vffcells)/((float) cellcounter)));

  fprintf(outputfile,"total volume: %f\ntotal mass:  %f\n",total_volume,total_mass);
  fprintf(outputfile,"number of cells:  %d\n",cellcounter);
  fprintf(outputfile,"total gas mass:  %f\ntotal metal mass: %f\nM_metal/M_gas:  %e\n",
	 totalgasmass,totalmetalmass,(totalmetalmass/totalgasmass));
  fprintf(outputfile,"total z1 mass:  %e  total z2 mass:  %e\n",totalz1mass,totalz2mass);
  fprintf(outputfile,"M_metal/M_gas/.02:  %e\n",
	  (totalmetalmass/totalgasmass/0.02));

  fprintf(outputfile,"MeanEntropyMassWeighted = %e\n",mean_entropy_mass);
  fprintf(outputfile,"MeanEntropyVolumeWeighted = %e\n",mean_entropy_volume);
  fprintf(outputfile,"MeanTempMassWeighted = %e\n",mean_temp_mass);
  fprintf(outputfile,"MeanTempVolumeWeighted = %e\n",mean_temp_volume);
  fprintf(outputfile,"ClumpingFactor = %f\n",clumping_factor);
  fprintf(outputfile,"MeanDensity = %f\n",total_mass/total_volume);
  fprintf(outputfile,"max temp:  %f\n",maxtemp);
  fprintf(outputfile,"max dens:  %f\n",maxdens);
  fprintf(outputfile,"BinSizeLogDens = %e\n",dbin1D);
  fprintf(outputfile,"BinSizeLogTemp = %e\n",tbin1D);
  fprintf(outputfile,"BinSizeLogEntropy = %e\n",sbin1D);

  // close mean data output file
  fclose(outputfile);

  // zero arrays and counters for cumulative dist'n fctns
  for(i=0; i < NUMBEROFBINS; i++) metalfrac_mass_cum[i]=metalfrac_vol_cum[i]=
				    z1frac_mass_cum[i]=z1frac_vol_cum[i]=
				    z2frac_mass_cum[i]=z2frac_vol_cum[i]=
				    ent_vol_cum[i]=ent_mass_cum[i]=
				    temp_vol_cum[i]=temp_mass_cum[i]=
				    dens_vol_cum[i]=dens_mass_cum[i]=
				    0.0;

  // keep track of total so we normalize to 1.0 - zero first
  metalfrac_mass_total=metalfrac_vol_total=
    z1frac_mass_total=z1frac_vol_total=
    z2frac_mass_total=z2frac_vol_total=
    dens_mass_total=dens_vol_total=
    temp_mass_total=temp_vol_total=
    ent_mass_total=ent_vol_total=
    0.0;

  // calculate cumulative dist'n fctns, counting backwards from the highest bin
  for(i=NUMBEROFBINS-1; i>= 0; i--){
    // get total (from bin i to NUMBEROFBINS -1
    metalfrac_mass_total+=metalfrac_mass[i];
    metalfrac_vol_total+=metalfrac_vol[i];

    z1frac_mass_total+=z1frac_mass[i];
    z1frac_vol_total+=z1frac_vol[i];

    z2frac_mass_total+=z2frac_mass[i];
    z2frac_vol_total+=z2frac_vol[i];

    dens_mass_total+=dens_mass[i];
    dens_vol_total+=dens_vol[i];

    temp_mass_total+=temp_mass[i];
    temp_vol_total+=temp_vol[i];

    ent_mass_total+=ent_mass[i];
    ent_vol_total+=ent_vol[i];

    // store in cumulative bin
    metalfrac_mass_cum[i]=metalfrac_mass_total;
    metalfrac_vol_cum[i]=metalfrac_vol_total;

    z1frac_mass_cum[i]=z1frac_mass_total;
    z1frac_vol_cum[i]=z1frac_vol_total;

    z2frac_mass_cum[i]=z2frac_mass_total;
    z2frac_vol_cum[i]=z2frac_vol_total;

    dens_mass_cum[i]=dens_mass_total;
    dens_vol_cum[i]=dens_vol_total;

    ent_mass_cum[i]=ent_mass_total;
    ent_vol_cum[i]=ent_vol_total;

    temp_mass_cum[i]=temp_mass_total;
    temp_vol_cum[i]=temp_vol_total;

  }

  // output metal distribution functions
  outputfile = fopen("metal_distfctn.dat","w");
  for(i=0; i<NUMBEROFBINS;i++)
    fprintf(outputfile,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n",
	    metalval[i],
	    (metalfrac_mass[i]/total_mass/mfracbin),
	    (metalfrac_vol[i]/total_volume/mfracbin),
	    (metalfrac_mass_cum[i]/metalfrac_mass_total),
	    (metalfrac_vol_cum[i]/metalfrac_vol_total), 

	    (z1frac_mass[i]/total_mass/mfracbin),
	    (z1frac_vol[i]/total_volume/mfracbin),
	    (z1frac_mass_cum[i]/z1frac_mass_total),
	    (z1frac_vol_cum[i]/z1frac_vol_total), 

	    (z2frac_mass[i]/total_mass/mfracbin),
	    (z2frac_vol[i]/total_volume/mfracbin),
	    (z2frac_mass_cum[i]/z2frac_mass_total),
	    (z2frac_vol_cum[i]/z2frac_vol_total) 
	    );

  fclose(outputfile);  // close output file

  // output density 1D distribution functions
  outputfile = fopen("dens_distfctn.dat","w");
  for(i=0; i<NUMBEROFBINS;i++)
    fprintf(outputfile,"%e %e %e %e %e\n",
	    density[i],
	    (dens_mass[i]/total_mass/dbin1D),
	    (dens_vol[i]/total_volume/dbin1D),
	    (dens_mass_cum[i]/dens_mass_total),
	    (dens_vol_cum[i]/dens_vol_total)
	    );

  fclose(outputfile);

  // output entropy 1D distribution fctns
  outputfile = fopen("ent_distfctn.dat","w");
  for(i=0; i<NUMBEROFBINS;i++)
    fprintf(outputfile,"%e %e %e %e %e\n",
	    ent[i],
	    (ent_mass[i]/total_mass/dbin1D),
	    (ent_vol[i]/total_volume/dbin1D),
	    (ent_mass_cum[i]/ent_mass_total),
	    (ent_vol_cum[i]/ent_vol_total)
	    );

  fclose(outputfile);

  // output temperature 1D distribution fctns
  outputfile = fopen("temp_distfctn.dat","w");
  for(i=0; i<NUMBEROFBINS;i++)
    fprintf(outputfile,"%e %e %e %e %e\n",
	    temperature[i],
	    (temp_mass[i]/total_mass/dbin1D),
	    (temp_vol[i]/total_volume/dbin1D),
	    (temp_mass_cum[i]/temp_mass_total),
	    (temp_vol_cum[i]/temp_vol_total)
	    );

  fclose(outputfile);

  dummy=PIXELS;

  // output all of the 2D distribution functions
  projname = fopen("entvsrho_vol.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(EVR_VolumeweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("entvsrho_mass.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(EVR_MassweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("tempvsrho_vol.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(TVR_VolumeweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("tempvsrho_mass.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(TVR_MassweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("metalvsrho_vol.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(MVR_VolumeweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("metalvsrho_mass.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(MVR_MassweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("metalvstemp_vol.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(MVT_VolumeweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("metalvstemp_mass.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(MVT_MassweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("pop3metrho_vol.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(P3VR_VolumeweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("pop3metrho_mass.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(P3VR_MassweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("pop3mettemp_vol.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(P3VT_VolumeweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  projname = fopen("pop3mettemp_mass.dat","w");
  fwrite(&dummy,sizeof(int),1,projname);
  fwrite(P3VT_MassweightedTable,sizeof(float),PIXELS*PIXELS,projname);
  fclose(projname);

  // everything's happy, exit!
  exit(0);
}
