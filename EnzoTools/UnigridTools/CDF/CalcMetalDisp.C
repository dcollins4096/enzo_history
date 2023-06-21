/*-----------------------------------------------------------


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
#define ODMIN -5.0
#define ODMAX 7.0

#define NUMBEROFBINS 1220  // number of bins in the arrays

// function declarations
void ReadGridFiles(int OperationFlag, char *hierarchyfilename);
void GetDistFctns(int ArraySize, float CellVolume);
void GetDispersion(int Arraysize, float CellVolume);
void PrintResults(char *outfilename);
void CreateArrays(void);
void DivideArrays(void);
void FinalizeArrays(void);

// global variable declarations
float  *densbuff,*metalbuff,*z1buff,*z2buff,
  logdrho,
  meanbaryondensity = 0.13333333,
  overdensity[NUMBEROFBINS],

  metal_mass[NUMBEROFBINS],
  metal_vol[NUMBEROFBINS],
  metal_mass_weights[NUMBEROFBINS],
  metal_vol_weights[NUMBEROFBINS],
  metal_mass_disp[NUMBEROFBINS],
  metal_vol_disp[NUMBEROFBINS],

  z1_mass[NUMBEROFBINS],
  z1_vol[NUMBEROFBINS],
  z1_mass_weights[NUMBEROFBINS],
  z1_vol_weights[NUMBEROFBINS],
  z1_mass_disp[NUMBEROFBINS],
  z1_vol_disp[NUMBEROFBINS],

  z2_mass[NUMBEROFBINS],
  z2_vol[NUMBEROFBINS],
  z2_mass_weights[NUMBEROFBINS],
  z2_vol_weights[NUMBEROFBINS],
  z2_mass_disp[NUMBEROFBINS],
  z2_vol_disp[NUMBEROFBINS]
  ;






int main(int argc, char *argv[]){
  // declare variables


  // file names and such
  char *outfilename=NULL,*hierarchyfilename=NULL,*datafilename=NULL;
  char *line = new char[MAX_LINE_LENGTH];
  char *gfname = new char[MAX_LINE_LENGTH];

  // read in command line arguments
  datafilename=argv[1];
  outfilename="dispersion.dat";

  printf("data dump name:          %s\n",datafilename);
  printf("dist fctn. file name:    %s\n",outfilename);


  // generate and zero arrays
  CreateArrays();

  // get hierarchy file name
  hierarchyfilename = datafilename;
  hierarchyfilename=strcat(hierarchyfilename,".hierarchy");
  printf("hierarchy file:          %s\n",hierarchyfilename);

  // read in grids one time, to get Z-Rho distribution functions
  ReadGridFiles(1,hierarchyfilename);

  DivideArrays();

  // read in grids second time, to get Z dispersion.
  ReadGridFiles(2,hierarchyfilename);

  FinalizeArrays();


  // print out all results!
  PrintResults(outfilename);

  return 1;
}

/*------------------------------------------------------------------*/
void ReadGridFiles(int OperationFlag,char *hierarchyfilename){

  if(DEBUG) fprintf(stderr,"ReadGridFiles: %d %s\n",OperationFlag,hierarchyfilename);

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

  int
    griddimx,griddimy,griddimz,  // grid xyz dims
    gridsex,gridsey,gridsez,  // grid start indices
    grideex,grideey,grideez,  // grid end indices
    gdx,gdy,gdz,  // grid dims - no buffers
    i,ndims;

  FILE *hierfile;
  
  char *line = new char[MAX_LINE_LENGTH];
  char *gfname = new char[MAX_LINE_LENGTH];
  
  // field names to extract
  char *density = "Density";
  char *metal = "Metal_Density";
  char *tracer1 = "Z_Field1";
  char *tracer2 = "Z_Field2";


 float celllength,cellvol,
    glex,gley,glez,  // grid left edge
    grex,grey,grez;  // grid right edge



  // Open up hierarchy file 
  hierfile = fopen(hierarchyfilename,"r");


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

      if(DEBUG) fprintf(stderr,"ReadGridFiles: opening %s\n",gfname);

      // calculate grid size w/out buffers
      gdx = 1+grideex-gridsex;
      gdy = 1+grideey-gridsey;
      gdz = 1+grideez-gridsez;

      // calculate mass conversion factor
      celllength = ((grex-glex) / ((float) gdx) );

 
      cellvol = celllength*celllength*celllength;

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

      if(OperationFlag==1) GetDistFctns( (int) size,cellvol );
      if(OperationFlag==2) GetDispersion( (int) size,cellvol );

      // clean up
      delete [] densbuff;
      delete [] metalbuff;
      delete [] z1buff;
      delete [] z2buff;

    }  // end of if(strncmp(line,"Grid = ",7)==0){
  }  // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){

  fclose(hierfile);  // close hierarchy file

}


/*------------------------------------------------------------------*/
void GetDistFctns(int ArraySize, float CellVolume){
  int i, densitybin;
  float overdensity,logod;

  for(i=0; i<ArraySize; i++){
    logod = log10(densbuff[i]/meanbaryondensity);

    densitybin = (int) ((logod - ODMIN) / logdrho);

    // error check!
    if( (densitybin < 0) || (densitybin>NUMBEROFBINS) ){
      fprintf(stderr, "GetDistFctns:  WTF? %d %d %f %f %f %f\n",densitybin,
	      NUMBEROFBINS, (densbuff[i]/meanbaryondensity), 
	      log10(densbuff[i]/meanbaryondensity),ODMAX,ODMIN);
      assert(1==0);  // force crash
    }

    // mass weighted distribution is just weighted by mass, which is
    // density times volume, but this is unigrid so the volumes are 
    // always the same
    metal_mass[densitybin] += metalbuff[i];
    metal_mass_weights[densitybin] += densbuff[i];

    z1_mass[densitybin] += z1buff[i];
    z1_mass_weights[densitybin] += densbuff[i];

    z2_mass[densitybin] += z2buff[i];
    z2_mass_weights[densitybin] += densbuff[i];

    // weight metal fraction by volume!
    metal_vol[densitybin] += metalbuff[i]/densbuff[i] * CellVolume;
    metal_vol_weights[densitybin] += CellVolume;

    z1_vol[densitybin] += z1buff[i]/densbuff[i] * CellVolume;
    z1_vol_weights[densitybin] += CellVolume;
    z2_vol[densitybin] += z2buff[i]/densbuff[i] * CellVolume;
    z2_vol_weights[densitybin] += CellVolume;

  }

}

/*------------------------------------------------------------------*/
void GetDispersion(int Arraysize, float CellVolume){

  int i, densitybin;
  float overdensity,logod,metallicity, z1_metallicity, z2_metallicity;

  for(i=0; i<Arraysize; i++){
    logod = log10(densbuff[i]/meanbaryondensity);

    densitybin = (int) ((logod - ODMIN) / logdrho);
    
    // error check!
    if( (densitybin < 0) || (densitybin>NUMBEROFBINS) ){
      fprintf(stderr, "GetDispersion:  WTF? %d %d %f %f %f %f\n",densitybin,
	      NUMBEROFBINS, (densbuff[i]/meanbaryondensity), 
	      log10(densbuff[i]/meanbaryondensity),ODMAX,ODMIN);
      assert(1==0);  // force crash
    }

    metallicity = metalbuff[i]/densbuff[i]/0.02;  // in units of Z/Zsolar
    z1_metallicity = z1buff[i]/densbuff[i]/0.02;  // in units of Z/Zsolar
    z2_metallicity = z2buff[i]/densbuff[i]/0.02;  // in units of Z/Zsolar


    // get dispersions - 'dispersion' arrays are summed 
    // ( val - <val> )^2 * weight for each cell.  The 'weight' values are
    // already calculated from GetDistFctns so we don't calculate them here.

    metal_mass_disp[densitybin] += 
      (metallicity-metal_mass[densitybin])*(metallicity-metal_mass[densitybin])*
      densbuff[i];
    
    z1_mass_disp[densitybin] += 
      (z1_metallicity-z1_mass[densitybin])*(z1_metallicity-z1_mass[densitybin])*
      densbuff[i];

    z2_mass_disp[densitybin] += 
      (z2_metallicity-z2_mass[densitybin])*(z2_metallicity-z2_mass[densitybin])*
      densbuff[i];

    metal_vol_disp[densitybin] += 
      (metallicity-metal_vol[densitybin])*(metallicity-metal_vol[densitybin])*
      CellVolume;
    
    z1_vol_disp[densitybin] += 
      (z1_metallicity-z1_vol[densitybin])*(z1_metallicity-z1_vol[densitybin])*
      CellVolume;

    z2_vol_disp[densitybin] += 
      (z2_metallicity-z2_vol[densitybin])*(z2_metallicity-z2_vol[densitybin])*
      CellVolume;

  }

}

void DivideArrays(void){
  int i;

  // divide through by weights to get mean metallicity as a 
  // function of overdensity!
  // also divide by 0.02 to get everything in terms of solar metallicity!
  for(i = 0; i < NUMBEROFBINS; i++){
    if(metal_mass_weights[i] > 0.0) metal_mass[i] /= (metal_mass_weights[i]*0.02);
    if(z1_mass_weights[i] > 0.0) z1_mass[i] /= (z1_mass_weights[i]*0.02);
    if(z2_mass_weights[i] > 0.0) z2_mass[i] /= (z2_mass_weights[i]*0.02);

    if(metal_vol_weights[i] > 0.0) metal_vol[i] /= (metal_vol_weights[i]*0.02);
    if(z1_vol_weights[i] > 0.0) z1_vol[i] /= (z1_vol_weights[i]*0.02);
    if(z2_vol_weights[i] > 0.0) z2_vol[i] /= (z2_vol_weights[i]*0.02);
  }

}


/*------------------------------------------------------------------*/
void FinalizeArrays(void){
  int i;

  for(i = 0; i < NUMBEROFBINS; i++){
    if(metal_mass_weights[i] > 0.0) 
      metal_mass_disp[i] = sqrt(  metal_mass_disp[i] / metal_mass_weights[i] );

    if(z1_mass_weights[i] > 0.0) 
      z1_mass_disp[i] = sqrt(  z1_mass_disp[i] / z1_mass_weights[i] );

    if(z2_mass_weights[i] > 0.0) 
      z2_mass_disp[i] = sqrt(  z2_mass_disp[i] / z2_mass_weights[i] );

    if(metal_vol_weights[i] > 0.0) 
      metal_vol_disp[i] = sqrt(  metal_vol_disp[i] / metal_vol_weights[i] );

    if(z1_vol_weights[i] > 0.0) 
      z1_vol_disp[i] = sqrt(  z1_vol_disp[i] / z1_vol_weights[i] );

    if(z2_vol_weights[i] > 0.0) 
      z2_vol_disp[i] = sqrt(  z2_vol_disp[i] / z2_vol_weights[i] );

  }



}

/*------------------------------------------------------------------*/
void PrintResults(char *outfilename){
  FILE *outfile;
  int i;
  outfile=fopen(outfilename,"w");

  for(i=0; i<NUMBEROFBINS; i++){

    fprintf(outfile,"%f %e %e %e %e %e %e %e %e %e %e %e %e\n",
	    overdensity[i],metal_mass[i],metal_mass_disp[i],
	    metal_vol[i],metal_vol_disp[i],z1_mass[i],z1_mass_disp[i],
	    z1_vol[i],z1_vol_disp[i],z2_mass[i],z2_mass_disp[i],
	    z2_vol[i],z2_vol_disp[i]);

    if(DEBUG)  printf("%f %e %e %e %e %e %e %e %e %e %e %e %e\n",
	    overdensity[i],metal_mass[i],metal_mass_disp[i],
	    metal_vol[i],metal_vol_disp[i],z1_mass[i],z1_mass_disp[i],
	    z1_vol[i],z1_vol_disp[i],z2_mass[i],z2_mass_disp[i],
	    z2_vol[i],z2_vol_disp[i]);

  }
  fclose(outfile);
}


/*------------------------------------------------------------------*/
void CreateArrays(void){
  int i;

  logdrho = (ODMAX - ODMIN)/NUMBEROFBINS; // get bin size

  for(i=0; i < NUMBEROFBINS; i++){
    metal_mass[i]=metal_vol[i]=metal_mass_weights[i]=
      metal_vol_weights[i]=metal_mass_disp[i]=metal_vol_disp[i]=0.0;
   
    if(i==0){
      overdensity[i] = logdrho/2.0 + ODMIN;
    } else{
      overdensity[i] = ((float) i + 0.5)*logdrho + ODMIN;
    }
  }
}
