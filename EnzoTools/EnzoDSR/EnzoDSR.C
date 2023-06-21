/* ---------------------------------------------------------------------- *

EnzoDSR.C

Sample enzo dataset reader

bwoshea, 22 June 2005

This is a sample Enzo dataset reader/writer for Gabe & Chris.  It's written
to be a serial code, so I imagine that there will have to be some fiddling 
done to parallelize it.

The way it works:
1) open up enzo parameter file and hierarchy file, get all of the necessary
   grid information, etc.
2) calculate conversion factors to go from Enzo internal code units to 
   proper (NOT comoving) CGS units
3) serially go through grids.  Open each one, extract all relevant datasets
   (basically all baryon quantities), modify these, write them back into files.
4) profit!

There are some user defined parameters in the header file, in lieu of having
command line arguments.

 * ---------------------------------------------------------------------- */

#include "EnzoDSR.h"

int main(){

  int i;

  // read Enzo restart dump parameter file
  ReadParameterFile(amrparfilename);

  numberofgrids=0;  // this is the total number of grids, set
		    // in GetGridInfo

  // read in enzo hierarchy file
  GetGridInfo(amrhierfilename);

  // set conversion factors
  SetConversionFactors();

  // loop through grids
  for(i=0; i<numberofgrids; i++)
    SetGridValues(i);

  return SUCCESS;
}


/* --------------------------- ReadParameterFile --------------------
 *  Parse parameter file looking for the top grid size and the maximum
 *  level.  For the purposes of an extraction we don't actually care
 *  about any cosmological parameters, so we don't get those.
 * ------------------------------------------------------------------ */
void ReadParameterFile(char *filename){
  fprintf(stderr,"ReadParameterFile:  Reading %s\n",filename);
  FILE *headerfile;
  //char *line = new char[MAX_LINE_LENGTH];
  char line[MAX_LINE_LENGTH];
  int ijunk;

  // open up header file and go through line by line
  headerfile=fopen(filename,"r");

  while( fgets(line, MAX_LINE_LENGTH, headerfile) != NULL){
    
    // top grid dimension
    if(strncmp(line,"TopGridDimensions",17)==0)
      sscanf(line,"TopGridDimensions   = %d %d %d",
	     &rootgridsize,&ijunk,&ijunk);

    // initial redshift
    if(strncmp(line,"CosmologyInitialRedshift",24)==0)
      sscanf(line,"CosmologyInitialRedshift = %lf",
             &initial_redshift);

    // omega matter
    if(strncmp(line,"CosmologyOmegaMatterNow",23)==0)
      sscanf(line,"CosmologyOmegaMatterNow    = %lf",
             &omegamatter);

    // comoving box size (in Mpc/h)
    if(strncmp(line,"CosmologyComovingBoxSize",24)==0)
      sscanf(line,"CosmologyComovingBoxSize   = %lf",
             &boxsize);

    // current redshift
    if(strncmp(line,"CosmologyCurrentRedshift",24)==0)
      sscanf(line,"CosmologyCurrentRedshift = %lf",
	     &redshift);

    // hubble constant (in units of 100 km/s/Mpc)
    if(strncmp(line,"CosmologyHubbleConstantNow",26)==0)
      sscanf(line,"CosmologyHubbleConstantNow = %lf",
	     &hubble);

  } // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)

  fclose(headerfile);


  //delete [] line;

  fprintf(stderr,"current redshift:  %g  hubble:  %g\n",omegamatter,hubble);
  fprintf(stderr,"Initial redshift:  %g  Box size:  %g  Omega Matter %g\n",
	 initial_redshift, boxsize, omegamatter);
  fprintf(stderr,"ReadParameterFile:  exiting\n");
}



/* --------------------------- GetGridInfo -----------------------------
 *  This routine goes through the hierarchy file twice.  The first time,
 *  it simple counts up the number of grid files so we can allocate memory
 *  for all of the arrays.  Then we create the arrays and read through the
 *  hierarchy file again to get info such as grid sizes (in number of
 *  cells), the start and end indices, and the level of the grid.
 * ------------------------------------------------------------------ */
void GetGridInfo(char *hierfilename){

  if(DEBUG) fprintf(stderr,"GetGridInfo:  entering.  %s. %i\n",hierfilename,numberofgrids);

  FILE *hierfile;
  int i,grid=0,griddimx,griddimy,griddimz,gridsex,gridsey,gridsez,
    grideex,grideey,grideez;
  char *line = new char[MAX_LINE_LENGTH];

  if(DEBUG) fprintf(stderr,"GetGridInfo: about to open hierarchy file.\n");

  /* open the hierarchy file once, and count the number of grids.
     Then close the hierarchy file */

  if(DEBUG) fprintf(stderr,"hierfilename:  %s\n",hierfilename);

  hierfile=fopen(hierfilename,"r");

  if(DEBUG) fprintf(stderr,"GetGridInfo: file opened, counting grids\n");
  // count grids
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){
    if(strncmp(line,"Grid = ",7)==0) numberofgrids++;
  }

  if(DEBUG) fprintf(stderr,"GetGridInfo: done counting, closing file.\n");

  fclose(hierfile);

  if(DEBUG) fprintf(stderr,"GetGridInfo:  There are %d grids\n",numberofgrids);

  
  /* initialize arrays which grid information is stored in
     (information such as level, grid bounds and number of
     cells, etc.) */
  gridlevel = new int[numberofgrids];  // level info

  // grid number of cells info
  griddx = new int[numberofgrids];  
  griddy = new int[numberofgrids];
  griddz = new int[numberofgrids];

  // grid bounds info
  gridlex = new double[numberofgrids];
  gridley = new double[numberofgrids];
  gridlez = new double[numberofgrids];
  gridrex = new double[numberofgrids];
  gridrey = new double[numberofgrids];
  gridrez = new double[numberofgrids];

  // grid file name info
  gridfilenames = new char*[numberofgrids];

  /* initialize all of the arrays to negative values 
     (in the case of the int and float arrays) or
     to a string in the case of the filename array. */
  for(i=0; i<numberofgrids; i++){
    gridlevel[i]=griddx[i]=griddy[i]=griddz[i]=-1;
    
    gridlex[i]=gridley[i]=gridlez[i]=
      gridrex[i]=gridrey[i]=gridrez[i]=-1.0;

    gridfilenames[i]= new char[MAX_LINE_LENGTH];
  }
  
  
  /* open the hierarchy file again.  This time, we read in all of 
     the various values that we want for the grids!
     
  */
  hierfile=fopen(hierfilename,"r");
  
  grid=0;
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
      sscanf(line,"GridLeftEdge      = %lf %lf %lf",
	     &gridlex[grid],&gridley[grid],&gridlez[grid]);

      fgets(line, MAX_LINE_LENGTH, hierfile);  // right edge
      sscanf(line,"GridRightEdge     = %lf %lf %lf",
	     &gridrex[grid],&gridrey[grid],&gridrez[grid]);

      for(i=0;i<3;i++) fgets(line, MAX_LINE_LENGTH, hierfile);

      // "junk" lines
      for(i=0;i<2;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
      

      // grid file name
      sscanf(line,"BaryonFileName = %s",gridfilenames[grid]);

      // "junk" lines
      for(i=0;i<5;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
      
      // grid dims - buffers stripped
      griddx[grid] = 1+grideex-gridsex;
      griddy[grid] = 1+grideey-gridsey;
      griddz[grid] = 1+grideez-gridsez;
      
      // calculate level from grid bounds, etc.
      gridlevel[grid] = (int) (log10(
				     (1.0 / ((double) rootgridsize)) / 
				     ( (gridrex[grid]-gridlex[grid]) / 
				       ( (double) griddx[grid] ) )
				     ) / log10(2.0)); 
      
      if(VERBOSEDEBUG) fprintf(stderr,"%lf %lf %lf    %lf %lf %lf    %d %d %d    %d    %s\n",
			       gridlex[grid],gridley[grid],gridlez[grid],
			       gridrex[grid],gridrey[grid],gridrez[grid],
			       griddx[grid],griddy[grid],griddz[grid],
			       gridlevel[grid],gridfilenames[grid]);

      grid++;  // increment grid number!

    }
    
  }

  fclose(hierfile);
   
  if(DEBUG) fprintf(stderr,"GetGridInfo:  exiting\n");
}


/* ------------------------------- diff -----------------------------
  calculates distance between (x1,y1,z1) and (x2,y2,z2)
 * ------------------------------------------------------------------ */
double diff(double x1, double y1, double z1, double x2, double y2, double z2){

  return sqrt( ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) );

}


/* ----------------------- SetConversionFactors ---------------------
   Set factors for conversion from "enzo units", which are somewhat
   insane (unless you do lots of acid, or use enzo far too much for
   the good of your mental health), into proper (NOT COMOVING) CGS
   units

   Conversion factors are defined at:
     http://cosmos.ucsd.edu/enzo/amr_guide/output.html
 * ------------------------------------------------------------------ */
void SetConversionFactors(void){

  density_conversion = omegamatter * RHOCRIT * hubble * hubble * 
    pow( (1.0+redshift), 3.0);

  length_conversion = boxsize * MPC_TO_CM / hubble / (1.0 + redshift) ;

  time_conversion = 4.0 * PI * GRAVC * omegamatter * RHOCRIT *
    pow( (1.0+initial_redshift), 3.0 );

  time_conversion = pow( time_conversion, -0.5);

  velocity_conversion = length_conversion/time_conversion * (1.0+redshift) / 
    (1.0 + initial_redshift);

  if(DEBUG){
    fprintf(stderr,"SetConversionFactors:\n\n");
    fprintf(stderr,"       density conversion:     %e\n", density_conversion);
    fprintf(stderr,"       length conversion:      %e\n", length_conversion);
    fprintf(stderr,"       time conversion:        %e\n", time_conversion);
    fprintf(stderr,"       velocity conversion:    %e\n", velocity_conversion);
  }

}



/* --------------------------- SetGridValues --------------------------
   Go through and put the actual supernova values in over the original
   grid values.  Extract arrays of grid information, modify them, 
   and then put the modified buffers back into the grids.
 * ------------------------------------------------------------------ */
void SetGridValues(int gridnumber){

  int i,j,k,ndims,m,buffer_index,cellindex,sn_int_radius,
    hlx,hly,hlz, 
    hrx,hry,hrz,sn_nx,sn_ny,sn_nz;
  
  double cellsize, ccx, ccy, ccz, supernova_radius, radial_distance;

  // hdf 5 stuff
  hid_t       file_id;
  hid_t       velx_dset_id, vely_dset_id, velz_dset_id,
    dens_dset_id, HI_dset_id, HII_dset_id, te_dset_id,
    ge_dset_id, temp_dset_id, elec_dset_id,HeI_dset_id,
    HeII_dset_id, HeIII_dset_id, H2I_dset_id, H2II_dset_id, 
    HM_dset_id, metal_dset_id;

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

  // buffers for grid info
  float *velxbuff,*velybuff,*velzbuff,
    *densbuff,*HIbuff,*HIIbuff,*tebuff,*gebuff,
    *tempbuff, *electronbuff,*HMbuff,
    *HeIbuff, *HeIIbuff,*HeIIIbuff, *H2Ibuff,*H2IIbuff,*metalbuff;

  // field names that will be extracted
  char *velfield_x = "x-velocity";
  char *velfield_y = "y-velocity";
  char *velfield_z = "z-velocity";
  char *densfield = "Density";
  char *HIfield = "HI_Density";
  char *HIIfield = "HII_Density";
  char *tefield = "Total_Energy";
  char *gefield = "Gas_Energy";
  char *tempfield = "Temperature";
  char *electronfield = "Electron_Density";
  char *HeIfield = "HeI_Density";
  char *HeIIfield = "HeII_Density";
  char *HeIIIfield = "HeIII_Density";
  char *H2Ifield = "H2I_Density";
  char *H2IIfield = "H2II_Density";
  char *HMfield = "HM_Density";
  char *metalfield = "Metal_Density";

  // GABE:  these are placeholder 32-bit floats for your data values.  I'm assuming
  // that all of your values are in CGS
  float SN_velx_cgs, SN_vely_cgs, SN_velz_cgs,   // x,y,z supernova velocities
	SN_density_cgs,		// total baryon density 
	SN_HIdens_cgs,		// neutral hydrogen density
	SN_HIIdens_cgs,		// ionized hydrogen density
	SN_totenergy_cgs,	// specific total energy (assumed to be in ergs/gram)
	SN_gasenergy_cgs,	// specific internal energy (assumed to be in ergs/gram)
	SN_temperature_cgs,	// temperature (Kelvin)
	SN_elecdens_cgs,	// electronn density - see note below about units
	SN_HeIdens_cgs,		// neutral helium density
	SN_HeIIdens_cgs,	// singly-ionized helium density
	SN_HeIIIdens_cgs,	// doubly ionized helium density
	SN_H2Idens_cgs,		// neutral molecular hydrogen
	SN_H2IIdens_cgs,	// ionized molecular hydrogen
	SN_HMdens_cgs,		// negatively charged atomic hydrogen
	SN_metaldens_cgs;	// metal field density

  supernova_radius = 1.5*SN_RADIUS/length_conversion;

  /*  
    Does the grid overlap the supernova?  If so, figure out which cells 
    actually overlap, and figure out their locations.  THEN we want to modify
    the density/energy/velocity/metal values in those cells, and then store 
    the values back in the grids
  */

    // does possible extent of sn overlap with this grid?
    if( ((SN_XPOS+supernova_radius) > gridlex[gridnumber] ) && ((SN_XPOS-supernova_radius) < gridrex[gridnumber] ) &&
	((SN_YPOS+supernova_radius) > gridley[gridnumber] ) && ((SN_YPOS-supernova_radius) < gridrey[gridnumber] ) &&
	((SN_ZPOS+supernova_radius) > gridlez[gridnumber] ) && ((SN_ZPOS-supernova_radius) < gridrez[gridnumber] )){
      
      if(DEBUG) fprintf(stderr,"Supernova overlaps grid %s\n",gridfilenames[gridnumber]);

      cellsize = (gridrex[gridnumber]-gridlex[gridnumber])/((float) griddx[gridnumber]);  
      
      // open grid file
      file_id = H5Fopen(gridfilenames[gridnumber], H5F_ACC_RDWR, H5P_DEFAULT);
      assert( file_id != h5_error );
      
      // open x velocity dataset
      velx_dset_id= H5Dopen(file_id,velfield_x);
      assert( velx_dset_id != h5_error );
      
      // open y velocity dataset
      vely_dset_id= H5Dopen(file_id,velfield_y);
      assert( vely_dset_id != h5_error );
      
      // open z velocity dataset
      velz_dset_id= H5Dopen(file_id,velfield_z);
      assert( velz_dset_id != h5_error );
	
      // open density dataset
      dens_dset_id = H5Dopen(file_id, densfield);
      assert( dens_dset_id != h5_error );
      
      // open HI density dataset
      HI_dset_id = H5Dopen(file_id, HIfield);
      assert( HI_dset_id != h5_error );
      
      // open HII density dataset
      HII_dset_id = H5Dopen(file_id, HIIfield);
      assert( HII_dset_id != h5_error );
      
      // open total energy dataset
      te_dset_id = H5Dopen(file_id, tefield);
      assert( te_dset_id != h5_error );
      
      // open gas energy dataset
      ge_dset_id = H5Dopen(file_id, gefield);
      assert( ge_dset_id != h5_error );
      
      // open temperature dataset
      temp_dset_id = H5Dopen(file_id, tempfield);
      assert( temp_dset_id != h5_error );
      
      // open electron density dataset
      elec_dset_id = H5Dopen(file_id, electronfield);
      assert( elec_dset_id != h5_error );
      
      // open neutral helium density dataset
      HeI_dset_id = H5Dopen(file_id, HeIfield);
      assert( HeI_dset_id != h5_error );
      
      // open singly ionized helium density dataset
      HeII_dset_id = H5Dopen(file_id, HeIIfield);
      assert( HeII_dset_id != h5_error );
      
      // open doubly ionized helium density dataset
      HeIII_dset_id = H5Dopen(file_id, HeIIIfield);
      assert( HeIII_dset_id != h5_error );
      
      // open neutral molecular hydrogen density dataset
      H2I_dset_id = H5Dopen(file_id, H2Ifield);
      assert( H2I_dset_id != h5_error );
      
      // open singly ionized molecular hydrogen density dataset
      H2II_dset_id = H5Dopen(file_id, H2IIfield);
      assert( H2II_dset_id != h5_error );
      
      // open negatively charged hydrogen density dataset
      HM_dset_id = H5Dopen(file_id, HMfield);
      assert( HM_dset_id != h5_error );

      // open metal density dataset
      metal_dset_id = H5Dopen(file_id, metalfield);
      assert(metal_dset_id != h5_error );
      
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
      
      if ( H5Tequal( typ_id, H5T_NATIVE_FLOAT ) ){
	
	// get memory for all of the buffers
	velxbuff = new float[(int) size];
	velybuff = new float[(int) size];
	velzbuff = new float[(int) size];
	densbuff = new float[(int) size];
	HIbuff = new float[(int) size];
	HIIbuff = new float[(int) size];
	tebuff = new float[(int) size];
	gebuff = new float[(int) size];
	tempbuff = new float[(int) size];
	electronbuff = new float[(int) size];
	HeIbuff = new float[(int) size];
	HeIIbuff = new float[(int) size];
	HeIIIbuff = new float[(int) size];
	H2Ibuff = new float[(int) size];
	H2IIbuff = new float[(int) size];
	HMbuff = new float[(int) size];
	metalbuff = new float[(int) size];
	
	mem_type_id = H5T_NATIVE_FLOAT;
	
	// read x velocity field
	h5_status = H5Dread(velx_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, velxbuff);
	assert( h5_status != h5_error ); 

	// read y velocity field 
	h5_status = H5Dread(vely_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, velybuff);
	assert( h5_status != h5_error );
	  
	// read z velocity field
	h5_status = H5Dread(velz_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, velzbuff);
	assert( h5_status != h5_error );
	
	// read density field
	h5_status = H5Dread(dens_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, densbuff);
	assert( h5_status != h5_error );
	
	// read HI density field
	h5_status = H5Dread(HI_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HIbuff);
	assert( h5_status != h5_error );
	  
	// read HII density field
	h5_status = H5Dread(HII_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HIIbuff);
	assert( h5_status != h5_error );
	
	// read total energy field
	h5_status = H5Dread(te_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, tebuff);
	assert( h5_status != h5_error );
	
	// read gas energy field
	h5_status = H5Dread(ge_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, gebuff);
	assert( h5_status != h5_error );

	// read temperature field
	h5_status = H5Dread(temp_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, tempbuff);
	assert( h5_status != h5_error );
	
	// read electron density field
	h5_status = H5Dread(elec_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, electronbuff);
	assert( h5_status != h5_error );
	
	// read HeI density field
	h5_status = H5Dread(HeI_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HeIbuff);
	assert( h5_status != h5_error );
	
	// read HeII density field
	h5_status = H5Dread(HeII_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HeIIbuff);
	assert( h5_status != h5_error );
	
	// read HeIII density field
	h5_status = H5Dread(HeIII_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HeIIIbuff);
	assert( h5_status != h5_error );
	
	// read neutral H2 density field
	h5_status = H5Dread(H2I_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, H2Ibuff);
	assert( h5_status != h5_error );
	  
	// read  H2II density field
	h5_status = H5Dread(H2II_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, H2IIbuff);
	assert( h5_status != h5_error );
	
	// read  HM density field
	h5_status = H5Dread(HM_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HMbuff);
	assert( h5_status != h5_error );

	// read  metal density field
	h5_status = H5Dread(metal_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, metalbuff);
	assert( h5_status != h5_error );
      }

      // calculate which cell the sn is in (ie, integer cell indices)
      sn_nx = (int) ((SN_XPOS - gridlex[gridnumber])/cellsize);
      sn_ny = (int) ((SN_YPOS - gridley[gridnumber])/cellsize);
      sn_nz = (int) ((SN_ZPOS - gridlez[gridnumber])/cellsize);
	
      // calculate sn radius in units of cell size
      // (plus some padding in case radius < cellsize)
      sn_int_radius = (int) (supernova_radius / cellsize + 1.0);
	
      /*  set left and right bounds that we're going to add metals in.  Double
	  sn radius to ensure that we get everything - it's no big deal b/c 
	  sn_int_radius is not that big.  Check to make sure that hl* and hr* 
          are within the bounds of the grid cell.
      */
	
      hlx = ( (sn_nx - 2*sn_int_radius) > 0 ) ? (sn_nx - 2*sn_int_radius) : 0;
      hly = ( (sn_ny - 2*sn_int_radius) > 0 ) ? (sn_ny - 2*sn_int_radius) : 0;
      hlz = ( (sn_nz - 2*sn_int_radius) > 0 ) ? (sn_nz - 2*sn_int_radius) : 0;
	
      hrx = ( (sn_nx + 2*sn_int_radius) < (griddx[gridnumber]-1) ) ? (sn_nx + 2*sn_int_radius) : (griddx[gridnumber]-1);
      hry = ( (sn_ny + 2*sn_int_radius) < (griddy[gridnumber]-1) ) ? (sn_ny + 2*sn_int_radius) : (griddy[gridnumber]-1);
      hrz = ( (sn_nz + 2*sn_int_radius) < (griddz[gridnumber]-1) ) ? (sn_nz + 2*sn_int_radius) : (griddz[gridnumber]-1);
	
      if(VERBOSEDEBUG) fprintf(stderr,"overlap: %d %d %d    %d %d %d\n",hlx,hly,hlz,hrx,hry,hrz);
      if(VERBOSEDEBUG) fprintf(stderr,"overlap grid dims: %d %d %d\n",
				griddx[gridnumber],griddy[gridnumber],griddz[gridnumber]);
      if(VERBOSEDEBUG) fprintf(stderr,"overlap sn pos: %d %d %d\n",sn_nx,sn_ny,sn_nz);
      if(VERBOSEDEBUG) fprintf(stderr,"overlap sn pos (double): %lf %lf %lf\n",
			(SN_XPOS - gridlex[gridnumber])/cellsize,
			(SN_YPOS - gridley[gridnumber])/cellsize,
			(SN_ZPOS - gridlez[gridnumber])/cellsize);

      
      // loop over all cells that are covered by the supernova
      for(k=hlz; k <= hrz; k++){	
	for(j=hly; j <= hry; j++){
	  for(i=hlx; i <= hrx; i++){
	      	    
	    // cell index (3d -> 1d conversion - USES FORTRAN ORDERING)
	    cellindex = k*griddx[gridnumber]*griddy[gridnumber] + j*griddx[gridnumber] + i;

	    // calculate center of this grid cell (i,j,k) in double precision
	    ccx = gridlex[gridnumber] + ( (gridrex[gridnumber]-gridlex[gridnumber]) / 
					((double) griddx[gridnumber]) ) * ( ((double) i) + 0.5 );
	    ccy = gridley[gridnumber] + ( (gridrey[gridnumber]-gridley[gridnumber]) / 
					((double) griddy[gridnumber]) ) * ( ((double) j) + 0.5 );
	    ccz = gridlez[gridnumber] + ( (gridrez[gridnumber]-gridlez[gridnumber]) / 
					((double) griddz[gridnumber]) ) * ( ((double) k) + 0.5 );

	    // calculate distance between center of cell and the sn
	    radial_distance = diff(SN_XPOS,SN_YPOS,SN_ZPOS,ccx,ccy,ccz);

            if(VERBOSEDEBUG) fprintf(stderr, "overlap:  cell at %lf %lf %lf is within sn radius %lf of sn center %lf %lf %lf\n",
				ccx,ccy,ccz,supernova_radius,SN_XPOS,SN_YPOS,SN_ZPOS);


	    // Now that we know where everything is, we have to input all of the values.
	    // GABE:  You'd probably need to put your stuff here!



	    /* -------- GABE:  This is where values need to be set! ------- */

	    // gas velocity
	    velxbuff[cellindex] = SN_velx_cgs/velocity_conversion_float;
	    velybuff[cellindex] = SN_vely_cgs/velocity_conversion_float;
	    velzbuff[cellindex] = SN_velz_cgs/velocity_conversion_float;
		  
	    // atomic hydrogen densities
	    HIbuff[cellindex] = SN_HIdens_cgs / density_conversion_float;
	    HIIbuff[cellindex] = SN_HIIdens_cgs / density_conversion_float;
	    HMbuff[cellindex] = SN_HMdens_cgs / density_conversion_float;

	    // molecular hydrogen densities
	    H2Ibuff[cellindex] = SN_H2Idens_cgs / density_conversion_float;
	    H2IIbuff[cellindex] = SN_H2IIdens_cgs / density_conversion_float;

            // helium densities
	    HeIbuff[cellindex] = SN_HeIdens_cgs / density_conversion_float;
	    HeIIbuff[cellindex] = SN_HeIIdens_cgs / density_conversion_float;
	    HeIIIbuff[cellindex] = SN_HeIIIdens_cgs / density_conversion_float;

	    // gas total and internal specific energy  (energy/particle)
	    tebuff[cellindex] = SN_totenergy_cgs / (velocity_conversion_float*velocity_conversion_float);
	    gebuff[cellindex] = SN_gasenergy_cgs / (velocity_conversion_float*velocity_conversion_float);

	    // gas temperature
	    tempbuff[cellindex] = SN_temperature_cgs;

	    // metal densities
	    metalbuff[cellindex] = SN_metaldens_cgs / density_conversion_float;


	    /* --------- GABE: LEAVE THESE ALONE  ---------*/

	    // electrons -- note that electron 'code units'
	    // are different than code units for everything else!
	    // don't use Zeus electron dens. values in order to maintain
	    // charge conservation - compute it from other stuff.
	    electronbuff[cellindex] = HIIbuff[cellindex] + 
	      (HeIIbuff[cellindex] / 4.0) +
	      (HeIIIbuff[cellindex] / 2.0) +
	      (H2IIbuff[cellindex] / 2.0) -
	       HMbuff[cellindex];

	    // total density -- computed from partial
	    // densities for accuracy!
	    densbuff[cellindex] = (HIbuff[cellindex] + HIIbuff[cellindex] +
	      HMbuff[cellindex] + H2Ibuff[cellindex] + H2IIbuff[cellindex] +
	      HeIbuff[cellindex] + HeIIbuff[cellindex] + HeIIIbuff[cellindex]);

	  } // for i
	}  // for j
      }  // for k


      /* ------------- DONE DOING STUFF TO THIS GRID - NOW WRITE DATA ----------- */

      // write modified x velocity buffer values to grid file
      h5_status = H5Dwrite(velx_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, velxbuff);
      assert( h5_status != h5_error );
	
      // write modified y velocity buffer values to grid file
      h5_status = H5Dwrite(vely_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, velybuff);
      assert( h5_status != h5_error );
      
      // write modified z velocity buffer values to grid file
      h5_status = H5Dwrite(velz_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, velzbuff);
      assert( h5_status != h5_error );
	
      // write modified density buffer values to grid file
      h5_status = H5Dwrite(dens_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, densbuff);
      assert( h5_status != h5_error );
	
      // write modified HI density buffer values to grid file
      h5_status = H5Dwrite(HI_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HIbuff);
      assert( h5_status != h5_error );
	
      // write modified HII density buffer values to grid file
      h5_status = H5Dwrite(HII_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HIIbuff);
      assert( h5_status != h5_error );
	
      // write modified total energy buffer values to grid file
      h5_status = H5Dwrite(te_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, tebuff);
      assert( h5_status != h5_error );
	
      // write modified gas energy buffer values to grid file
      h5_status = H5Dwrite(ge_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, gebuff);
      assert( h5_status != h5_error );

      // write modified temperature buffer values to grid file
      h5_status = H5Dwrite(temp_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, tempbuff);
      assert( h5_status != h5_error );

      // write modified electron density buffer values to grid file
      h5_status = H5Dwrite(elec_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, electronbuff);
      assert( h5_status != h5_error );

      // write modified HeI density buffer values to grid file
      h5_status = H5Dwrite(HeI_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HeIbuff);
      assert( h5_status != h5_error );

      // write modified HeII density buffer values to grid file
      h5_status = H5Dwrite(HeII_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HeIIbuff);
      assert( h5_status != h5_error );

      // write modified HeIII density buffer values to grid file
      h5_status = H5Dwrite(HeIII_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HeIIIbuff);
      assert( h5_status != h5_error );

      // write modified neutral H2 density buffer values to grid file
      h5_status = H5Dwrite(H2I_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, H2Ibuff);
      assert( h5_status != h5_error );

      // write modified negative H2 density buffer values to grid file
      h5_status = H5Dwrite(H2II_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, H2IIbuff);
      assert( h5_status != h5_error );

      // write modified negative H density buffer values to grid file
      h5_status = H5Dwrite(HM_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, HMbuff);
      assert( h5_status != h5_error );

      // write modified negative H density buffer values to grid file
      h5_status = H5Dwrite(metal_dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, metalbuff);
      assert( h5_status != h5_error );

      // cleanup
      delete [] velxbuff;
      delete [] velybuff;
      delete [] velzbuff;
      delete [] densbuff;
      delete [] HIbuff;
      delete [] HIIbuff;
      delete [] tebuff;
      delete [] gebuff;
      delete [] tempbuff;
      delete [] electronbuff;
      delete [] HeIbuff;
      delete [] HeIIbuff;
      delete [] HeIIIbuff;
      delete [] H2Ibuff;
      delete [] H2IIbuff;
      delete [] HMbuff;
      delete [] metalbuff;
      
      // close everything, doing appropriate error checking
      h5_status = H5Sclose(dsp_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Tclose(typ_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Sclose(mem_dsp_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Sclose(file_dsp_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Dclose(velx_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(vely_dset_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Dclose(velz_dset_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Dclose(dens_dset_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Dclose(HI_dset_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Dclose(HII_dset_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Dclose(te_dset_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Dclose(ge_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(temp_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(elec_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(HeI_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(HeII_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(HeIII_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(H2I_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(H2II_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(HM_dset_id);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(metal_dset_id);
      assert( h5_status != h5_error );
	
      h5_status = H5Fclose(file_id);
      assert( h5_status != h5_error );
      
    } // if( ((SN_XPOS ..... loop

}
