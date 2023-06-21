#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#include "extern_hdf5.h"  // hdf 5 prototypes

#define MAX_LINE_LENGTH 256
#define OMEGA_BARYON 0.04
#define OMEGA_MATTER 0.30
#define DEBUG 0
#define VERBOSEDEBUG 0
#define SUCCESS 1
#define FAILURE 0

// -------------- USER MUST SET THIS ---------------
// in Msolar, not Msolar/h

float dmpmass = 166.17; // 256^3 DM particles in a 0.3 Mpc box, 
                        // Omega_CDM = 0.26 and h=0.7


// global hierarchy file values
int *gridlevel,*griddx,*griddy,*griddz,*gridnump;
float *gridlex,*gridley,*gridlez,
  *gridrex,*gridrey,*gridrez;
char **gridfilenames;
float *densitybuff,*dmdensitybuff,*H2buff,*HIIbuff,boxsize;
int *flagbuffer;
int rootgridsize, particleonly;
float mean_density,omegamatter,omegalambda,massconv,hubble,redshift;


float *haloxpos,*haloypos,*halozpos,*halomaxrad,*halogasmass,*halodmmass,*haloH2mass, 
  *haloHIImass, *halotmass;
int numberofhalos,dark_matter_exists;


int ReadParameterFile(char *filename);
int GetGridInfo(int numberofgrids,char *hierfilename);
int NumberOfGrids(char *hierfilename);
int GetCellInformationFromGrid(int gridnum,int total_number_grids);
int GetParticleInformationFromGrid(int gridnum);
int FlagGridCells(int gridnum,int total_number_grids);
int GetDensityInfo(int gridnum);
float diff(float x1, float y1, float z1, float x2, float y2, float z2);


void ReadHaloFile(char *halofilename);
void CalculateVirialRadii(void);

int main(int argc, char *argv[]){
  char *inputfilename=NULL;
  char *halofilename=NULL;
  int junk,i,return_value;
  int total_number_grids,total_uo_cells;
  inputfilename=argv[1];
  halofilename=argv[2];

  FILE *output;


  boxsize = omegamatter = omegalambda = hubble = -1.0;

  junk =  ReadParameterFile(inputfilename);

  massconv = 2.78 * ((float) pow(10.0,11.0) ) * boxsize * boxsize * boxsize * omegamatter / hubble;

  fprintf(stderr,"%e %f %f %f %f %f\n", 2.78 * ((float) pow(10.0,11.0) ),boxsize,omegamatter,omegalambda,redshift,hubble);
  fprintf(stderr,"Mass conversion is:  %e\n",massconv);

  ReadHaloFile(halofilename);

  CalculateVirialRadii();

  // get hierarchy file name
  int inlen = (int) strlen(inputfilename);
  char *hierfile = new char[inlen+1];
  strcpy(hierfile,inputfilename);
  hierfile=strcat(hierfile,".hierarchy");

  // get total number of grids - to generate grid information
  total_number_grids = NumberOfGrids(hierfile);

  if(DEBUG) fprintf(stderr,"there are %i grids\n",total_number_grids);

  /* initialize arrays which grid information is stored in
     (information such as level, grid bounds and number of
     cells, etc.) */
  gridlevel = new int[total_number_grids];  // level info

  // grid number of cells info
  griddx = new int[total_number_grids];  
  griddy = new int[total_number_grids];
  griddz = new int[total_number_grids];

  // grid bounds info
  gridlex = new float[total_number_grids];
  gridley = new float[total_number_grids];
  gridlez = new float[total_number_grids];
  gridrex = new float[total_number_grids];
  gridrey = new float[total_number_grids];
  gridrez = new float[total_number_grids];

  // grid file name info
  gridfilenames = new char*[total_number_grids];

  // number of particles info
  gridnump = new int[total_number_grids];

  /* initialize all of the arrays to negative values 
     (in the case of the int and float arrays) or
     to a string in the case of the filename array. */
  for(i=0; i<total_number_grids; i++){
    gridlevel[i]=griddx[i]=griddy[i]=griddz[i]=gridnump[i]=-1;

    gridlex[i]=gridley[i]=gridlez[i]=
      gridrex[i]=gridrey[i]=gridrez[i]=-1.0;

    gridfilenames[i]= new char[MAX_LINE_LENGTH];
  }

  // get grid information from hierarchy file
  return_value = GetGridInfo(total_number_grids,hierfile);
  assert( return_value == total_number_grids );


  for(i=0; i<total_number_grids; i++){
    return_value=GetCellInformationFromGrid(i,total_number_grids);
    assert( return_value != FAILURE );    
  }

  /*
  for(i=0; i<total_number_grids; i++){
    return_value=GetParticleInformationFromGrid(i);
    assert( return_value != FAILURE );   
  }
  */


  //PrintOutValues();

  // clean up dynamically allocated stuff
  for(i=0; i<total_number_grids; i++)
    delete [] gridfilenames[i];
  
  delete [] gridfilenames;
  delete [] gridlevel;
  delete [] griddx;
  delete [] griddy;
  delete [] griddz;
  delete [] gridlex;
  delete [] gridley;
  delete [] gridlez;
  delete [] gridrex;  
  delete [] gridrey;  
  delete [] gridrez;  
  delete [] hierfile;
  delete [] gridnump;

  output = fopen("massfrac.dat","w");

  for(i=0; i<numberofhalos;i++)
    fprintf(output,"%e %e %e %e %e %f %e %e\n",
	    halogasmass[i],
	    halodmmass[i],
	    haloH2mass[i],
	    haloHIImass[i],
	    halotmass[i],
	    (halogasmass[i]/(halodmmass[i] + halogasmass[i])/(0.04/0.3)),
	    haloH2mass[i]/halogasmass[i],
	    haloHIImass[i]/halogasmass[i]
	    );

  fclose(output);

  if(DEBUG) fprintf(stderr,"exiting...\n");
  exit(0);
}

int GetParticleInformationFromGrid(int gridnum){
  if(DEBUG) fprintf(stderr,"in GetParticleInformationFromGrid %d\n",gridnum);

  // in case there aren't any particles!
  if(gridnump[gridnum]==0) return SUCCESS;

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

  // one per dataset
  hid_t ppx_dset_id, ppy_dset_id, ppz_dset_id;

  double *particlex,*particley,*particlez;

  int i, j, numpart,  ndims;

  // open file - only once
  fprintf(stderr,"GetParticleInformationFromGrid %s\n",gridfilenames[gridnum]);
  file_id = H5Fopen(gridfilenames[gridnum], H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );
  
  // open part. pos. x
  ppx_dset_id = H5Dopen(file_id, "particle_position_x");
  assert( ppx_dset_id != h5_error );

  // open part. pos. y
  ppy_dset_id = H5Dopen(file_id, "particle_position_y");
  assert( ppx_dset_id != h5_error );

  // open part. pos. z
  ppz_dset_id = H5Dopen(file_id, "particle_position_z");
  assert( ppx_dset_id != h5_error );


  // open ppx dataspace (to get dimensions) 
  // only once!
  dsp_id = H5Dget_space(ppx_dset_id);
  assert( dsp_id != h5_error );

  // get data type (only once!)
  typ_id = H5Dget_type(ppx_dset_id);
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

      particlex = new double[(int) size];
      particley = new double[(int) size];
      particlez = new double[(int) size];


      mem_type_id = H5T_IEEE_F64BE;

      // read part. pos. x field into an array
      h5_status = H5Dread(ppx_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, particlex);
      if(DEBUG) fprintf(stderr,"float read status %d for ppx field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 
     
      // read part. pos. y field into an array
      h5_status = H5Dread(ppy_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, particley);
      if(DEBUG) fprintf(stderr,"float read status %d for ppy field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

       // read part. pos. z field into an array
      h5_status = H5Dread(ppz_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, particlez);
      if(DEBUG) fprintf(stderr,"float read status %d for ppz field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

    }



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
  h5_status = H5Dclose(ppx_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(ppy_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(ppz_dset_id);
  assert( h5_status != h5_error );

  // ---------- close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );


  numpart = (int) size;  // # of particles


  // loop over all halos and particles on this grid.  If the halo and particle are
  // within halomaxrad of each other, add a dark matter particle mass to the total.
  // booyah.
  for(i=0; i < numberofhalos; i++)
    for(j=0; j < numpart; j++)
      if( diff( ((float) particlex[j]), ((float) particley[j]), ((float) particlez[j]), 
		haloxpos[i], haloypos[i], halozpos[i]) < halomaxrad[i] ) 
	halodmmass[i] += dmpmass;

  delete [] particlex;
  delete [] particley;
  delete [] particlez;


  if(DEBUG) fprintf(stderr,"exiting GetParticleInformationFromGrid\n");
  return SUCCESS;

}



/*------------------ GetCellInformationFromGrid -----------------------
 *
 *
 *--------------------------------------------------------------*/
int GetCellInformationFromGrid(int gridnum,int total_number_grids){

  if(DEBUG) fprintf(stderr,"in GetCellInformationFromGrid %d\n",gridnum);



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

  // need one of these per dataset!
  hid_t dens_dset_id,dmdens_dset_id,temp_dset_id,
     H2_dset_id,
    metal_dset_id,z1_dset_id,z2_dset_id,
    velocity_dset_id,HI_dset_id,HII_dset_id,HeI_dset_id,
    HeII_dset_id,HeIII_dset_id,e_dset_id,HM_dset_id,
    H2I_dset_id,H2II_dset_id;

  int i, ndims;

  dark_matter_exists = 0;

  // open grid file and extract various quantities of interest

  // open file - only once
  fprintf(stderr,"GetCellInformationFromGrid %s\n",gridfilenames[gridnum]);
  file_id = H5Fopen(gridfilenames[gridnum], H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );
  
  // open density dataset
  dens_dset_id = H5Dopen(file_id, "Density");
  assert( dens_dset_id != h5_error );

  // dark matter density dataset
  dmdens_dset_id = H5Dopen(file_id,"Dark_Matter_Density");
  //assert( dmdens_dset_id != h5_error );

  // open HII density dataset
  HII_dset_id = H5Dopen(file_id, "HII_Density");
  assert( HII_dset_id != h5_error );

  // open H2 density dataset
  H2_dset_id = H5Dopen(file_id, "H2I_Density");
  assert( H2_dset_id != h5_error );


  if( dmdens_dset_id != h5_error) dark_matter_exists = 1;
  

  // open density dataspace (to get dimensions) 
  // only once!
  dsp_id = H5Dget_space(dens_dset_id);
  assert( dsp_id != h5_error );

  // get data type (only once!)
  typ_id = H5Dget_type(dens_dset_id);
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
      // allocate buffers - one per projection buffer!
      densitybuff = new float[(int) size];
      dmdensitybuff = new float[(int) size];
      H2buff = new float[(int) size];
      HIIbuff = new float[(int) size];

      mem_type_id = H5T_IEEE_F32BE;

      // read density field into an array
      h5_status = H5Dread(dens_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, densitybuff);
      if(DEBUG) fprintf(stderr,"float read status %d for density field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      // read HII density field into an array
      h5_status = H5Dread(HII_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, HIIbuff);
      if(DEBUG) fprintf(stderr,"float read status %d for HII density field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      // read H2 density field into an array
      h5_status = H5Dread(H2_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, H2buff);
      if(DEBUG) fprintf(stderr,"float read status %d for H2 density field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      if(dark_matter_exists){
      // read dark matter density field into an array
      h5_status = H5Dread(dmdens_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, dmdensitybuff);
      if(DEBUG) fprintf(stderr,"float read status %d for dark matter density field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 
      }


      
    }

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

  h5_status = H5Dclose(HII_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(H2_dset_id);
  assert( h5_status != h5_error );

  // ---------- must close each dataset - one per projection buffer!
  h5_status = H5Dclose(dmdens_dset_id);
  //assert( h5_status != h5_error );

  // ---------- close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  // ---------- create flag field and zero it out
  flagbuffer = new int[(int) size];
  for(i=0; i < ((int) size) ; i++) flagbuffer[i] = 0;

  // ---------- check density buffer - debug operation
  for(i=0; i < ((int) size) ; i++){
    if(densitybuff[i]<=0.0) 
      fprintf(stderr,"GetCellInformationFromGrid: negative or zero density value!\n");
  }

  FlagGridCells(gridnum,total_number_grids);

  // ----------------------- GET INFORMATION! ----------
  GetDensityInfo(gridnum);


  delete [] flagbuffer;
  delete [] densitybuff;
  delete [] HIIbuff;
  delete [] H2buff;
  delete [] dmdensitybuff;

  if(DEBUG) fprintf(stderr,"exiting GetCellInformationFromGrid\n");
  return SUCCESS;
}



/*-------------------- GetDensityInfo ------------------------
 *
 * This routine adds information to projections along the x axis.
 * 
 *--------------------------------------------------------------*/
int GetDensityInfo(int gridnum){
  if(DEBUG) fprintf(stderr,"in GetDensityInfo\n");
  float cellsize,gridmass,massdensconv,halopos_x,halopos_y,halopos_z;
  int i,j,k,cellindex,totalcells,cellcounter,odbin,massbin,m;
  int halo_nx, halo_ny, halo_nz,hlx,hly,hlz,hrx,hry,hrz,haloradius;
  float ccx,ccy,ccz;
  if(DEBUG) fprintf(stderr,"1\n");

  cellcounter = 0;

  totalcells = griddx[gridnum]*griddy[gridnum]*griddz[gridnum];

  // size of grid cell, in code units
  cellsize = (gridrex[gridnum] - gridlex[gridnum]) / ((float) griddx[gridnum] );

  // mass conversion - all you need to do is multiply by density (in code units) 
  // and you've got a mass.
  massdensconv = cellsize * cellsize * cellsize * massconv;


  for(m=0; m<numberofhalos; m++)
    for(halopos_x = (haloxpos[m]-1.0); halopos_x <= (haloxpos[m]+1.001); halopos_x+=1.0 )
      for(halopos_y = (haloypos[m]-1.0); halopos_y <= (haloypos[m]+1.001); halopos_y+=1.0 )
	for(halopos_z = (halozpos[m]-1.0); halopos_z <= (halozpos[m]+1.001); halopos_z+=1.0 )
	  // does halo overlap with this grid?
	  if( ((halopos_x+halomaxrad[m]) > gridlex[gridnum] ) && ((halopos_x-halomaxrad[m]) < gridrex[gridnum] ) &&
	      ((halopos_y+halomaxrad[m]) > gridley[gridnum] ) && ((halopos_y-halomaxrad[m]) < gridrey[gridnum] ) &&
	      ((halopos_z+halomaxrad[m]) > gridlez[gridnum] ) && ((halopos_z-halomaxrad[m]) < gridrez[gridnum] )
	      ){  // if...
	    
	    
	    
	    // calculate which cell the halo is in
	    halo_nx = (int) ((halopos_x - gridlex[gridnum])/cellsize);
	    halo_ny = (int) ((halopos_y - gridley[gridnum])/cellsize);
	    halo_nz = (int) ((halopos_z - gridlez[gridnum])/cellsize);
	    
	    // calculate halo radius in units of cell size
	    // (plus some padding in case radius < cellsize)
	    haloradius = (int) (halomaxrad[m] / cellsize + 1.0);
	    
	    /*  set left and right bounds that we're going to add metals in.  Double
		halo radius to ensure that we get everything - it's no big deal b/c 
		haloradius is on the order of a few cells.  Check to make sure that
		hl* and hr* are within the bounds of the grid cell.
	    */
	    
	    hlx = ( (halo_nx - 2*haloradius) > 0 ) ? (halo_nx - 2*haloradius) : 0;
	    hly = ( (halo_ny - 2*haloradius) > 0 ) ? (halo_ny - 2*haloradius) : 0;
	    hlz = ( (halo_nz - 2*haloradius) > 0 ) ? (halo_nz - 2*haloradius) : 0;
	    
	    hrx = ( (halo_nx + 2*haloradius) < (griddx[gridnum]-1) ) ? (halo_nx + 2*haloradius) : (griddx[gridnum]-1);
	    hry = ( (halo_ny + 2*haloradius) < (griddy[gridnum]-1) ) ? (halo_ny + 2*haloradius) : (griddy[gridnum]-1);
	    hrz = ( (halo_nz + 2*haloradius) < (griddz[gridnum]-1) ) ? (halo_nz + 2*haloradius) : (griddz[gridnum]-1);
	    

	    if(VERBOSEDEBUG){
	      printf("overlap: %d %d %d    %d %d %d\n",hlx,hly,hlz,hrx,hry,hrz);
	      printf("overlap grid dims: %d %d %d\n",griddx[gridnum],griddy[gridnum],griddz[gridnum]);
	      printf("overlap halo pos: %d %d %d\n",halo_nx,halo_ny,halo_nz);
	      printf("overlap halo pos (float): %f %f %f\n",(halopos_x - gridlex[gridnum])/cellsize,
		     (halopos_y - gridley[gridnum])/cellsize,
		     (halopos_z - gridlez[gridnum])/cellsize);
	    }

	    // loop over all cells that are covered by the halo
	    for(i=hlx; i <= hrx; i++){
	      for(j=hly; j <= hry; j++){
		for(k=hlz; k <= hrz; k++){
		  
		  // cell index (3d -> 1d conversion - USES FORTRAN ORDERING)
		  cellindex = k*griddx[gridnum]*griddy[gridnum] + j*griddx[gridnum] + i;
		  
		  // calculate center of cell (i,j,k)
		  ccx = gridlex[gridnum] + ( (gridrex[gridnum]-gridlex[gridnum]) / ((float) griddx[gridnum]) ) * ( ((float) i) + 0.5 );
		  ccy = gridley[gridnum] + ( (gridrey[gridnum]-gridley[gridnum]) / ((float) griddy[gridnum]) ) * ( ((float) j) + 0.5 );
		  ccz = gridlez[gridnum] + ( (gridrez[gridnum]-gridlez[gridnum]) / ((float) griddz[gridnum]) ) * ( ((float) k) + 0.5 );
		  
		  // is the center of this cell within the halo radius?
		  // if so, increment metal density and tracer density by rho_z_halo
		  //if( diff(halopos_x,ccx,halopos_y,ccy,halopos_z,ccz) < halomaxrad[m]){
		  if( diff(halopos_x,halopos_y,halopos_z,ccx,ccy,ccz) < halomaxrad[m]){

		    if(VERBOSEDEBUG){
		      printf("overlap:  cell at %f %f %f is within halo radius %f of halo center %f %f %f\n",
			     ccx,ccy,ccz,halomaxrad[m],halopos_x,halopos_y,halopos_z);
		      printf("Grid coords:  %d %d %d %d %d %d %d %d\n",i,j,k,
			     griddx[gridnum],griddy[gridnum],griddz[gridnum],cellindex,totalcells);
		    }

		    if(flagbuffer[cellindex] != -1){ // read just non-overlapped grid cells
		      halogasmass[m] += densitybuff[cellindex]*massdensconv;
		      haloH2mass[m] += H2buff[cellindex]*massdensconv;
		      haloHIImass[m] += HIIbuff[cellindex]*massdensconv;
		      if(dark_matter_exists) halodmmass[m] += dmdensitybuff[cellindex]*massdensconv;

		    }
		  } // if( diff
		  
		}  // for(k=hlz; k <= hrz; k++){
	      }  // for(j=hly; j <= hry; j++){
	    }  // for(i=hlx; i <= hrx; i++){
	  } // if( ((halopos_x+halomaxrad[m]) > gridlex[gridnum] )
	    
	    
  if(DEBUG) fprintf(stderr,"GetDensityInfo:  Cells counted: %d  Total Cells:  %d\n",cellcounter,totalcells);
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

  if(DEBUG) fprintf(stderr,"in FlagGridCells: %i %i\n",gridnum,total_number_grids);

  int counter, flagged_cells_this_grid,i,j,k,cellindex;

  float clex, cley, clez, crex, crey, crez;  // cell left and right edges (x,y,z)

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

	  clex = gridlex[gridnum] + (( (float) i ) * (gridrex[gridnum] - gridlex[gridnum])/
				     ((float) griddx[gridnum]));

	  crex = gridlex[gridnum] + (( (float) (i+1) ) * (gridrex[gridnum] - gridlex[gridnum])/
				     ((float) griddx[gridnum]));

	  cley = gridley[gridnum] + (( (float) j ) * (gridrey[gridnum] - gridley[gridnum])/
				     ((float) griddy[gridnum]));

	  crey = gridley[gridnum] + (( (float) (j+1) ) * (gridrey[gridnum] - gridley[gridnum])/
				     ((float) griddy[gridnum]));

	  clez = gridlez[gridnum] + (( (float) k ) * (gridrez[gridnum] - gridlez[gridnum])/
				     ((float) griddz[gridnum]));

	  crez = gridlez[gridnum] + (( (float) (k+1) ) * (gridrez[gridnum] - gridlez[gridnum])/
				     ((float) griddz[gridnum]));

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

	  //if(flagbuffer[cellindex] != -1){
	  //  flagbuffer[cellindex]=-1;
	  //  flagged_cells_this_grid++;
	  //}
	  
	}  // for(k=...
      } // for(j=0...
    }  // for(i=0...

    if(DEBUG) fprintf(stderr,"FlagGridCells:  In Grid: %i (level %i) Grid Looked At: %i (level %i) cells flagged this grid (ttl): %i\n",
		      gridnum,gridlevel[gridnum],counter,gridlevel[counter],flagged_cells_this_grid);
    if(DEBUG) fprintf(stderr,"FlagGridCells: %f %f %f   %f %f %f\n",gridlex[counter],gridley[counter],gridlez[counter],
		      gridrex[counter],gridrey[counter],gridrez[counter]);

  }  // for(counter=0...

  if(DEBUG) fprintf(stderr,"exiting FlagGridCells...\n");
  return SUCCESS;
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
    grideex,grideey,grideez,nbfields;
  
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
	     &gridlex[grid],&gridley[grid],&gridlez[grid]);

      fgets(line, MAX_LINE_LENGTH, hierfile);  // right edge
      sscanf(line,"GridRightEdge     = %f %f %f",
	     &gridrex[grid],&gridrey[grid],&gridrez[grid]);
      
      for(i=0;i<3;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
      sscanf(line,"NumberOfBaryonFields = %d",&nbfields);

      if(nbfields==0){
	particleonly=1;

	fgets(line, MAX_LINE_LENGTH, hierfile);
	sscanf(line,"NumberOfParticles   = %d",
	       &gridnump[grid]);  // number of particles on this grid

	fgets(line, MAX_LINE_LENGTH, hierfile);
	sscanf(line,"ParticleFileName = %s",gridfilenames[grid]);

	//fprintf(stderr,"GetGridInfo:  No baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

      } else {
	// "junk" lines
	for(i=0;i<2;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
	
	// grid file name
	sscanf(line,"BaryonFileName = %s",gridfilenames[grid]);

	// "junk" lines
	for(i=0;i<5;i++) fgets(line, MAX_LINE_LENGTH, hierfile);

	sscanf(line,"NumberOfParticles   = %d",
	       &gridnump[grid]);  // number of particles on this grid

	//fprintf(stderr,"GetGridInfo: WITH baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

      }

      // grid dims - buffers stripped
      griddx[grid] = 1+grideex-gridsex;
      griddy[grid] = 1+grideey-gridsey;
      griddz[grid] = 1+grideez-gridsez;

      // calculate level from grid bounds, etc.
      gridlevel[grid] = (int) (log10(
			      (1.0 / ((float) rootgridsize)) / 
			      ( (gridrex[grid]-gridlex[grid]) / 
				( (float) griddx[grid] ) )
			      ) / log10(2.0)); 

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
  int numgrids=0;

  if(DEBUG) fprintf(stderr,"hierarchy file: %s\n",hierfilename);
    
  hierfile = fopen(hierfilename,"r");

  // read through hierarchy file, counting # of grids
  // lines that start with "Grid = " are the beginning of a new
  // piece of grid info
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)
    if(strncmp(line,"Grid = ",7)==0) numgrids++;
    
  fclose(hierfile);

  if(DEBUG) fprintf(stderr,"NumberOfGrids:  there are %i grids\n",numgrids);
  
  // clean up dynamically allocated stuff
  delete [] line;

  if(DEBUG) fprintf(stderr,"exiting NumberOfGrids\n");

  // return # of grids
  return numgrids;  
}



/* ---------------------------- ReadParameterFile -----------------------------
 * 
 *
 *
 * ---------------------------------------------------------------------------- */

int ReadParameterFile(char *filename){
  FILE *headerfile;
  char *line = new char[MAX_LINE_LENGTH];
  int ijunk;

  if(DEBUG) fprintf(stderr,"ReadParameterFile:  reading %s\n",filename);

  // open up header file and go through line by line
  headerfile=fopen(filename,"r");

  while( fgets(line, MAX_LINE_LENGTH, headerfile) != NULL){

    // comoving box size (in Mpc/h)
    if(strncmp(line,"CosmologyComovingBoxSize",24)==0){
      sscanf(line,"CosmologyComovingBoxSize   = %g",
	     &boxsize);
      fprintf(stderr,"found boxsize:  %s\n",line);
      fprintf(stderr,"box size is: %g\n",boxsize);
    }
    
    // top grid dimension
    if(strncmp(line,"TopGridDimensions",17)==0)
      sscanf(line,"TopGridDimensions   = %d %d %d",
	     &rootgridsize,&ijunk,&ijunk);

    // omega lambda
    if(strncmp(line,"CosmologyOmegaLambdaNow",23)==0)
      sscanf(line,"CosmologyOmegaLambdaNow    = %f",
	     &omegalambda);

    // omega matter
    if(strncmp(line,"CosmologyOmegaMatterNow",23)==0)
      sscanf(line,"CosmologyOmegaMatterNow    = %f",
	     &omegamatter);

    // redshift
    if(strncmp(line,"CosmologyCurrentRedshift",24)==0)
      sscanf(line,"CosmologyCurrentRedshift   = %f",
	     &redshift);

    // hubble constant (in units of 100 km/s/Mpc)
    if(strncmp(line,"CosmologyHubbleConstantNow",26)==0)
      sscanf(line,"CosmologyHubbleConstantNow = %f",
	     &hubble);


  } // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)

  fclose(headerfile);

  if(DEBUG) fprintf(stderr,"ReadParameterFile:  root grid is %d\n",rootgridsize);

  return 1;

}


void ReadHaloFile(char *halofilename){

  FILE *halofile;
  char *line = new char[MAX_LINE_LENGTH];
  int ijunk,halo,i,*halonumparticles;
  float fjunk;

  numberofhalos=0;

  if(DEBUG) fprintf(stderr,"ReadHaloFile:  reading %s\n",halofilename);

  // open up header file and go through line by line
  halofile=fopen(halofilename,"r");

  fgets(line, MAX_LINE_LENGTH, halofile); // read junk line at the top!

  // count number of halos
  while(fscanf(halofile,"%d %f %d %e %f %f %f",&ijunk,&fjunk,&ijunk,&fjunk,&fjunk,&fjunk,&fjunk)==7) numberofhalos++;
  fclose(halofile);

  haloxpos = new float[numberofhalos];
  haloypos = new float[numberofhalos];
  halozpos = new float[numberofhalos];
  halomaxrad = new float[numberofhalos];
  halogasmass = new float[numberofhalos];
  halodmmass = new float[numberofhalos];
  haloH2mass = new float[numberofhalos];
  haloHIImass = new float[numberofhalos];
  halotmass= new float[numberofhalos];
  halonumparticles = new int[numberofhalos];

  // open up header file, again, and go through line by line reading in values
  halofile=fopen(halofilename,"r");

  fgets(line, MAX_LINE_LENGTH, halofile); // read junk line at the top!

  // read in halos - column order is  intjunk, #numpart, intjunk, floatjunk, xpos,ypos,zpos
  // (we're talking about Hop files here)
  halo = 0;
  while(fscanf(halofile,"%d %d %d %f %f %f %f",
	       &ijunk,&halonumparticles[halo],&ijunk,&fjunk,
	       &haloxpos[halo],&haloypos[halo],&halozpos[halo])==7) halo++;

  // &ijunk,&halomaxrad[halo],&ijunk,&halotmass[halo],
  // &haloxpos[halo],&haloypos[halo],&halozpos[halo])==7) halo++;

  fclose(halofile);

  if(halo != numberofhalos) fprintf(stderr,"ReadHaloFile:  WTF?  %d %d\n",halo, numberofhalos);

  for(i=0; i<numberofhalos; i++){
    halomaxrad[i]=0.1;
    halotmass[i] = dmpmass * ((float) halonumparticles[i]);

    // zero these out.
    halogasmass[i] = halodmmass[i] = haloH2mass[i] = haloHIImass[i] = 0.0;

    fprintf(stderr,"%f %f %f %f %e %f %f %f %f\n",
	    haloxpos[i], haloypos[i], halozpos[i],
	    halomaxrad[i], halotmass[i], halodmmass[i],
	    halogasmass[i], haloH2mass[i], haloHIImass[i]);

  }



  fprintf(stderr,"ReadHaloFile:  there are %d halos\n",numberofhalos);


  //fclose(halofile);

  return;
}


float diff(float x1, float y1, float z1, float x2, float y2, float z2){

  return sqrt( ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) );

}

void CalculateVirialRadii(void){

  int i;
  float omega_mz, d, delta_c, pi;
  float old_radius, virial_radius;
  float properbox;

  pi = 3.1415926;

  omega_mz = omegamatter * pow( (1.0+redshift),3.0) / 
    ( omegamatter * pow( (1.0+redshift),3.0) + omegalambda );

  d = omega_mz - 1.0;

  delta_c = 18.0*pi*pi + 82.0*d - 39.0*d*d;

  // get box size in proper kpc/h

  properbox = boxsize * 1000.0 / (1.0+redshift);

  for(i=0; i < numberofhalos; i++){
    // old_radius stores "max radius" in units of box size
    old_radius = halomaxrad[i];

    // calculate virial radius in kpc/h as per
    // Barkana & Loeb 2000, eqtn. 24
    virial_radius = 0.784 * pow( ( halotmass[i] *hubble / 1.0e+8), 0.3333 )
      * pow( ( omegamatter*delta_c / omega_mz / 18.0 / pi / pi), -0.3333 )
      * pow( ( (1.0+redshift)/10.0 ), -1.0);
    
    // get virial radius in code units
    virial_radius /= properbox;

    // multiple of virial radius

    virial_radius *= 4.0;

    //if(DEBUG) 
    fprintf(stderr,"old, new halo radii: %f %f\n",old_radius,virial_radius);

    // new radius in box units!
    halomaxrad[i] = virial_radius;
  
  }

}
