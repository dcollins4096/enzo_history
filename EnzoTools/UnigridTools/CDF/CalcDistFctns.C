/*-----------------------------------------------------------

calcdf (CalcDistFctns.C)

This program goes through all of the grid files in a given
data dump (from the .hierarchy file) and calculates the
following information:

1.  The 1d mass- and volume-weighted metallicity distribution
    functions (differential and cumulative)

2.  volume filling factor as a function of metallicity

3.  various mean values

This program extracts this info and puts it into
a text file named distfctninfo.dat

Usage:

    calcdf <data dump name>

Example:

    calcdf RedshiftOutput0011

Brian O'Shea
bwoshea@cosmos.ucsd.edu
1/13/2002

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

#define DEBUG 0  // 1 on, 0 off
#define VDEBUG 0 // verbose debug - 1 on, 0 off

// these define the bounds of our various arrays
//   (all in log10 space)
#define ZMIN -10.0  // metallicity
#define ZMAX 2.0    

#define METAL_MULTIPLY


// These numbers are the minimum metallicity used
// for calculating the volume filling factor (VFF)
// for the various metal quantities.
//#define Z_VFF_MIN  0.000001
//#define Z1_VFF_MIN 0.000001
//#define Z2_VFF_MIN 0.000001

#define NUMBEROFBINS 100  // number of bins in the arrays

int main(int argc, char *argv[]){

  FILE *hierfile, *outputfile, *outputfile1, *headerfile;

  // file names and such
  char *outfilename=NULL,*outfile1name=NULL,*hierarchyfilename=NULL,*datafilename=NULL;
  char *line = new char[MAX_LINE_LENGTH];
  char *gfname = new char[MAX_LINE_LENGTH];

  // field names to extract
  char *density = "Density";
  char *metal = "Metal_Density";
  char *tracer1 = "Z_Field1";
  char *tracer2 = "Z_Field2";

  float
    *densbuff,*metalbuff,*z1buff,*z2buff,
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
    z2_min,z2_max,vff_min[7];


  int
    griddimx,griddimy,griddimz,  // grid xyz dims
    gridsex,gridsey,gridsez,  // grid start indices
    grideex,grideey,grideez,  // grid end indices
    gdx,gdy,gdz,  // grid dims - no buffers
    i, ndims,cellcounter,m,
    metalbin,
    z1bin,
    z2bin,
    vffcells[7],
    z1_vffcells[7],
    z2_vffcells[7];

  float
    mffcells[7],
    z1_mffcells[7],
    z2_mffcells[7];


  // hdf 5 stuff
  hid_t       file_id;
  hid_t       dens_dset_id, metal_dset_id,
    z1_dset_id,z2_dset_id;
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
  outfilename="distfctninfo.dat";
  outfile1name="meaninfo.dat";

  printf("data dump name:          %s\n",datafilename);
  printf("dist fctn. file name:    %s\n",outfilename);
  printf("mean data file name:     %s\n",outfile1name);

  // ------------ set up distribution function arrays 

  // bin size (in log10 space)
  mfracbin=(ZMAX-ZMIN)/NUMBEROFBINS;

  // calculate values in the center of the bin so that
  // the plots look pretty 'n stuff
  for(i=0;i<NUMBEROFBINS;i++){
    metalfrac_mass[i]=metalfrac_vol[i]=0.0;
    z1frac_mass[i]=z1frac_vol[i]=0.0;
    z2frac_mass[i]=z2frac_vol[i]=0.0;

    if(i==0){
      metalval[i]=mfracbin/2.0+ZMIN;
    } else{
      metalval[i]=((float) i + 0.5)*mfracbin + ZMIN;
    }
  }

  // zero mean quantities
  total_volume=total_mass=mean_metal_mass=mean_metal_vol=
    mean_z1_mass=mean_z1_vol=mean_z2_mass=mean_z2_vol=
    totalz1mass=totalz2mass=totalmetalmass=totalgasmass=0.0;

  // count cells to make sure we have the right number!
  cellcounter=0;
 

  //------ REMOVE ------
  if(DEBUG) fprintf(stderr,"(1)\n");
  //------ REMOVE ------

  // keep track of volume filling factor  
  for(i=0; i<7; i++){
    vff_min[i] = (float) pow( 10.0, ((double) (i-7)) );
    vffcells[i]=z1_vffcells[i]=z2_vffcells[i]=0;
    mffcells[i]=z1_mffcells[i]=z2_mffcells[i]=0.0;
  }


  // set initial (out-of-bounds) values for z, z1, z2 min/max
  z_min=z1_min=z2_min=100.0;
  z_max=z1_max=z2_max=-1.0;

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

      dens_dset_id = H5Dopen(file_id,density);
      assert( dens_dset_id != h5_error );

      metal_dset_id = H5Dopen(file_id,metal);
      assert( metal_dset_id != h5_error );

      z1_dset_id = H5Dopen(file_id,tracer1);
      assert( z1_dset_id != h5_error);

      z2_dset_id = H5Dopen(file_id,tracer2);
      assert( z2_dset_id != h5_error);


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
	}

      // calculate distribution functions and weights

      // go through buffers
      for(i=0; i<size; i++){
	// get metallicity in terms of solar metallicity

	metallicity = metalbuff[i]/densbuff[i]/0.02;
	z1_metallicity = z1buff[i]/densbuff[i]/0.02;
	z2_metallicity = z2buff[i]/densbuff[i]/0.02;

#ifdef METAL_MULTIPLY
	metallicity *= 130.0 / 50.0;
	z1_metallicity *= 130.0 / 50.0;
	z2_metallicity *= 130.0 / 50.0;
#endif

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


	// increment mean quantities
	mean_metal_mass += metallicity * densbuff[i] * cellvol;
	mean_metal_vol += metallicity * cellvol;

	mean_z1_mass += z1_metallicity * densbuff[i] * cellvol;
	mean_z1_vol += z1_metallicity * cellvol;

	mean_z2_mass += z2_metallicity * densbuff[i] * cellvol;
	mean_z2_vol += z2_metallicity * cellvol;

	// increment total "mass" (ie, density) and volume
	total_mass += densbuff[i] * cellvol;
	total_volume += cellvol;

	// calculate which metallicity bin this cell belongs to
	// (all three metal variables)
	metalbin=(int) ( (log10(metallicity) - ZMIN)/mfracbin);
	z1bin=(int) ( (log10(z1_metallicity) - ZMIN)/mfracbin);
	z2bin=(int) ( (log10(z2_metallicity) - ZMIN)/mfracbin);
	// error check - arrays in bounds?
	// it's just a check - if they DON'T fit, we just put it in the lowest
	// bin.
	if(metalbin < 0 || metalbin >= NUMBEROFBINS || 
	   z1bin >= NUMBEROFBINS || z2bin >= NUMBEROFBINS){
	  fprintf(stderr,"wtf? %i %i %i %i %e %e %e\n",metalbin,z1bin,z2bin,NUMBEROFBINS,
		  metallicity,z1_metallicity,z2_metallicity);
	  //exit(10);
	}

	// increment volume filling factor # cells
	for(m=0; m < 7; m++){
	  if(metallicity >= vff_min[m]){
	    vffcells[m]++;
	    mffcells[m]+=densbuff[i];
	  }
	  if(z1_metallicity >= vff_min[m]){
	    z1_vffcells[m]++;
	    z1_mffcells[m]+=densbuff[i];
	  }
	  if(z2_metallicity >= vff_min[m]){
	    z2_vffcells[m]++;
	    z2_mffcells[m]+=densbuff[i];
	  }
	}


	
	//if(metallicity >= Z_VFF_MIN) vffcells++;
	//if(z1_metallicity >= Z1_VFF_MIN) z1_vffcells++;
	//if(z2_metallicity >= Z2_VFF_MIN) z2_vffcells++;
	
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

  // open up mean data output file
  outputfile1 = fopen(outfile1name,"w");

  //fprintf(outputfile1,"mean metallicity (in units of Z/Z_solar):\n");
  fprintf(outputfile1,"Z: mass weighted mean:      %e\n",mean_metal_mass);
  fprintf(outputfile1,"Z: volume weighted mean:    %e\n",mean_metal_vol);
  fprintf(outputfile1,"Z: max metallicity:         %e\n",z_max);
  fprintf(outputfile1,"Z: min metallicity:         %e\n",z_min);
  for(m=0; m <7; m++)
    fprintf(outputfile1,"Z: volume filling factor    %e  above %e  Z/Zsolar\n", 
	    (((float) vffcells[m])/((float) cellcounter)), vff_min[m] );
  for(m=0; m <7; m++)
    fprintf(outputfile1,"Z: mass filling factor    %e  above %e  Z/Zsolar\n", 
	    (((float) mffcells[m])/totalgasmass), vff_min[m] );


  fprintf(outputfile1,"\n");


  //fprintf(outputfile1,"mean z1 metallicity (in units of Z/Z_solar):\n");
  fprintf(outputfile1,"Z1: mass weighted mean:      %e\n",mean_z1_mass);
  fprintf(outputfile1,"Z1: volume weighted mean:    %e\n",mean_z1_vol);
  fprintf(outputfile1,"Z1: max metallicity:         %e\n",z1_max);
  fprintf(outputfile1,"Z1: min metallicity:         %e\n",z1_min);
  for(m=0; m <7; m++)
    fprintf(outputfile1,"Z1: volume filling factor:   %e  above %e  Z/Zsolar\n", 
	    (((float) z1_vffcells[m])/((float) cellcounter)), vff_min[m] );
  for(m=0; m <7; m++)
    fprintf(outputfile1,"Z1: mass filling factor    %e  above %e  Z/Zsolar\n", 
	    (((float) z1_mffcells[m])/totalgasmass), vff_min[m] );

  fprintf(outputfile1,"\n");


  //fprintf(outputfile1,"mean z2 metallicity (in units of Z/Z_solar):\n");
  fprintf(outputfile1,"Z2: mass weighted mean:      %e\n",mean_z2_mass);
  fprintf(outputfile1,"Z2: volume weighted mean:    %e\n",mean_z2_vol);
  fprintf(outputfile1,"Z2: max metallicity:         %e\n",z2_max);
  fprintf(outputfile1,"Z2: min metallicity:         %e\n",z2_min);
  for(m=0; m <7; m++)
    fprintf(outputfile1,"Z2: volume filling factor:   %e  above %e  Z/Zsolar\n",
	    (((float) z2_vffcells[m])/((float) cellcounter)), vff_min[m] );
  for(m=0; m <7; m++)
    fprintf(outputfile1,"Z2: mass filling factor    %e  above %e  Z/Zsolar\n", 
	    (((float) z2_mffcells[m])/totalgasmass), vff_min[m] );

  fprintf(outputfile1,"\n");
  fprintf(outputfile1,"\n");
  fprintf(outputfile1,"total volume: %f\ntotal mass:  %f\n",total_volume,total_mass);
  fprintf(outputfile1,"number of cells:  %i\n",cellcounter);
  fprintf(outputfile1,"total gas mass:  %f\ntotal metal mass: %f\nM_metal/M_gas:  %e\n",
	 totalgasmass,totalmetalmass,(totalmetalmass/totalgasmass));
  fprintf(outputfile1,"M_metal/M_gas/.02:  %e\n",
	  (totalmetalmass/totalgasmass/0.02));

  // close mean data output file
  fclose(outputfile1);

  // zero arrays and counters for cumulative dist'n fctns
  for(i=0; i < NUMBEROFBINS; i++) metalfrac_mass_cum[i]=metalfrac_vol_cum[i]=
				    z1frac_mass_cum[i]=z1frac_vol_cum[i]=
				    z2frac_mass_cum[i]=z2frac_vol_cum[i]=0.0;

  // keep track of total so we normalize to 1.0
  metalfrac_mass_total=metalfrac_vol_total=
    z1frac_mass_total=z1frac_vol_total=
    z2frac_mass_total=z2frac_vol_total=0.0;

  // calculate cumulative dist'n fctns
  for(i=NUMBEROFBINS-1; i>= 0; i--){
    metalfrac_mass_total+=metalfrac_mass[i];
    metalfrac_vol_total+=metalfrac_vol[i];

    z1frac_mass_total+=z1frac_mass[i];
    z1frac_vol_total+=z1frac_vol[i];

    z2frac_mass_total+=z2frac_mass[i];
    z2frac_vol_total+=z2frac_vol[i];

    metalfrac_mass_cum[i]=metalfrac_mass_total;
    metalfrac_vol_cum[i]=metalfrac_vol_total;

    z1frac_mass_cum[i]=z1frac_mass_total;
    z1frac_vol_cum[i]=z1frac_vol_total;

    z2frac_mass_cum[i]=z2frac_mass_total;
    z2frac_vol_cum[i]=z2frac_vol_total;
  }

  // open up distribution function output file
  outputfile = fopen(outfilename,"w");

  //------------- write out all of the stuff to file
  for(i=0; i<NUMBEROFBINS;i++)
    fprintf(outputfile,"%f %f %f %f %f %f %f %f %f %f %f %f %f\n",
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
  exit(0);
}
