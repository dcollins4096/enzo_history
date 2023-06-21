/* ------------------------------- DensVar.C ------------------------- *
   Calculates the density variance within a given region.
 * ------------------------------------------------------------------- */

#include "DensVar.h"
/*------------------------------------------------------------------
 *MAKE SURE TO ADD NEW 1D and 2D Distributions for multispecies = 1!
 *----------------------------------------------------------------*/



/*--------------------------- main -----------------------------
 *
 *
 *--------------------------------------------------------------*/
int main(int argc, char *argv[]){
  char *inputfilename=NULL;

  int junk,i,return_value;
  inputfilename=argv[1];


  boxsize = omegamatter = omegalambda = hubble = -1.0;
  rootgridsize = multispecies = -1;

  // read parameter file, get some values
  return_value =  ReadParameterFile(inputfilename);

  if(DEBUG) fprintf(stderr,"%lf %lf %lf %lf %lf %d %d %lf\n", boxsize,omegamatter,
		    omegalambda,redshift,hubble,rootgridsize,multispecies, currenttime);

  if(DEBUG) fprintf(stderr,"Main: about to get hierarchy info\n");

  // get hierarchy file name
  int inlen = (int) strlen(inputfilename);
  char *hierfile = new char[inlen+6];
  strcpy(hierfile,inputfilename);
  hierfile=strcat(hierfile,".hierarchy");

  if(DEBUG) fprintf(stderr,"Main: hierarchy file name is:  %s\n",hierfile);

  // get total number of grids - to generate grid information
  total_number_grids = NumberOfGrids(hierfile);

  if(DEBUG) fprintf(stderr,"Main: there are %i grids\n",total_number_grids);

  // declares grid arrays - basically just uses new and sets things to zero
  DeclareGridArrays();

  if(DEBUG) fprintf(stderr,"Main: entering GetGridInfo\n");

  // get grid information from hierarchy file
  return_value = GetGridInfo(total_number_grids,hierfile);
  assert( return_value == total_number_grids );

  delete [] hierfile;

  meanoverdensity = 0.0;
  meanoverdensity_weight = 0.0;
  totalmass = totalvolume = 0.0;

  // Loop through grids and get information for each grid, making sure
  // to zero under subgrids
  if(DEBUG) fprintf(stderr,"Main: entering GetCellInformation\n");
  for(i=0; i<total_number_grids; i++){
    return_value=GetCellInformationFromGrid(i,total_number_grids,DENSITY);
    assert( return_value != FAILURE );    
  }

  meanoverdensity /= meanoverdensity_weight;

  totalvariance = totalvariance_weight = 0.0;

  // Loop through grids and get information for each grid, making sure
  // to zero under subgrids
  if(DEBUG) fprintf(stderr,"Main: entering GetCellInformation (Variance)\n");
  for(i=0; i<total_number_grids; i++){
    return_value=GetCellInformationFromGrid(i,total_number_grids,VARIANCE);
    assert( return_value != FAILURE );    
  }

  totalvariance = pow( totalvariance / totalvariance_weight, 0.5);

  // Actually write everything out
  if(DEBUG) fprintf(stderr,"Main: entering OutputAllInformation\n");
  OutputAllInformation(argv[1]);

  // Out of here.
  if(DEBUG) fprintf(stderr,"Main: exiting...\n");
  exit(0);
}


/*------------------ OutputAllInformation -----------------------
 *
 *  Pretty straightforward, really.  Writes all of the information
 *  out to files.
 *
 *--------------------------------------------------------------*/
void OutputAllInformation(char *infilename){
  if(DEBUG) fprintf(stderr,"In OutputAllInformation  \n");
  int i, numpixels;
  FILE *output;

  // misc. information, mean information
  output = fopen("varianceinfo.dat","w");
  fprintf(output,"# mean information file\n");
  fprintf(output,"InputFileName       = %s\n",infilename);
  fprintf(output,"Redshift            = %lf\n",redshift);
  fprintf(output,"CurrentTime         = %lf\n",currenttime);
  fprintf(output,"BoxSize             = %lf\n",boxsize);
  fprintf(output,"RootGridSize        = %d\n",rootgridsize);
  fprintf(output,"TotalNumberOfGrids  = %d\n",total_number_grids);
  fprintf(output,"TotalMassInSim      = %e\n",totalmass);
  fprintf(output,"TotalVolumeInSim    = %e\n",totalvolume);
  fprintf(output,"\n");
  fprintf(output,"mean overdensity     = %e\n",meanoverdensity);
  fprintf(output,"overdensity variance = %e\n",totalvariance);
 
  if(DEBUG) fprintf(stderr,"Leaving OutputAllInformation  \n");
  return;
}



/*--------------------- DeclareGridArrays ----------------------
 *
 *  This simple function just uses new to declare the arrays for
 *  grid information like position, level, etc.  Also set various
 *  values to absurd values so we can error-check later.
 *
 *--------------------------------------------------------------*/
void DeclareGridArrays(void){
  if(DEBUG) fprintf(stderr,"In DeclareGridArrays\n");

  int i;

  /* initialize arrays which grid information is stored in
     (information such as level, grid bounds and number of
     cells, etc.) */
  gridlevel = new int[total_number_grids];  // level info

  // grid number of cells info
  griddx = new int[total_number_grids];  
  griddy = new int[total_number_grids];
  griddz = new int[total_number_grids];

  // grid bounds info
  gridlex = new double[total_number_grids];
  gridley = new double[total_number_grids];
  gridlez = new double[total_number_grids];
  gridrex = new double[total_number_grids];
  gridrey = new double[total_number_grids];
  gridrez = new double[total_number_grids];

  // grid file name info
  gridfilenames = new char*[total_number_grids];

  // number of particles info
  gridnump = new int[total_number_grids];

  /* initialize all of the arrays to negative values 
     (in the case of the int and double arrays) or
     to a string in the case of the filename array. */
  for(i=0; i<total_number_grids; i++){
    gridlevel[i]=griddx[i]=griddy[i]=griddz[i]=gridnump[i]=-1;

    gridlex[i]=gridley[i]=gridlez[i]=
      gridrex[i]=gridrey[i]=gridrez[i]=-1.0;

    gridfilenames[i]= new char[MAX_LINE_LENGTH];
  }

  if(DEBUG) fprintf(stderr,"Leaving DeclareGridArrays  \n");
  return;
}

/*---------------------- CleanGridArrays -----------------------
 *
 *  Cleans up the memory consumed by the grid arrays created in 
 *  DeclareGridArrays.  Nothing exciting here.
 *  
 *--------------------------------------------------------------*/
void CleanGridArrays(void){
  if(DEBUG) fprintf(stderr,"In CleanGridArrays  \n");

  int i;

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
  delete [] gridnump;

  if(DEBUG) fprintf(stderr,"Leaving CleanGridArrays\n");
  return;
}


/*------------------ GetCellInformationFromGrid ---------------------------
 *
 *  This function opens up the grid, extracts all of the baryon quantities
 *  we care about.  It then calls FlagGridCells and GetDensityInfo, which is
 *  the routine that actually calculates all of the quantities that we care
 *  about.
 *
 *-------------------------------------------------------------------------*/
int GetCellInformationFromGrid(int gridnum,int total_number_grids, int whichroutine){

  if(DEBUG) fprintf(stderr,"in GetCellInformationFromGrid %d\n",gridnum);

  // is any part of grid within bounds of projection?  if not, exit!
  if( !((gridrex[gridnum] > BEGIN_POS) && (gridlex[gridnum] < END_POS) && 
        (gridrey[gridnum] > BEGIN_POS) && (gridley[gridnum] < END_POS) && 
        (gridrez[gridnum] > BEGIN_POS) && (gridlez[gridnum] < END_POS))
      ) return SUCCESS;


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
  hid_t dens_dset_id;

  int i, ndims;

  // open grid file and extract various quantities of interest

  // open file - only once
  if(VERBOSEDEBUG) fprintf(stderr,"GetCellInformationFromGrid %s\n",gridfilenames[gridnum]);
  file_id = H5Fopen(gridfilenames[gridnum], H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );
  
  // open density dataset
  dens_dset_id = H5Dopen(file_id, "Density");
  assert( dens_dset_id != h5_error );

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

  // get arrays
  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
    {
      // allocate buffers - one per projection buffer!
      densitybuff = new float[(int) size];

      mem_type_id = H5T_NATIVE_FLOAT;

      // read density field into an array
      h5_status = H5Dread(dens_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, densitybuff);
      if(VERBOSEDEBUG) fprintf(stderr,"double read status %d for density field\n", 
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


  // ---------- must close each dataset - one per buffer extracted!
  h5_status = H5Dclose(dens_dset_id);
  assert( h5_status != h5_error );

  // ---------- close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  // ---------- create flag field and zero it out
  flagbuffer = new int[(int) size];
  for(i=0; i < ((int) size) ; i++) flagbuffer[i] = 0;

  // ---------- check density buffer as a debug operation
  for(i=0; i < ((int) size) ; i++){
    if(densitybuff[i]<=0.0) 
      fprintf(stderr,"GetCellInformationFromGrid: negative or zero density value!\n");
  }

  // flag grid cells here!
  FlagGridCells(gridnum,total_number_grids);

  // CALL THE ROUTINE THAT ACTUALLY GETS THE INFORMATION! 
  if(whichroutine==DENSITY)
    GetDensityInfo(gridnum);

  if(whichroutine==VARIANCE)
    GetVarianceInfo(gridnum);

  // ------------------- erase buffers - once per buffer!
  delete [] flagbuffer;
  delete [] densitybuff;

  if(VERBOSEDEBUG) fprintf(stderr,"exiting GetCellInformationFromGrid\n");
  return SUCCESS;
}



/*------------------------- GetDensityInfo ----------------------------
 *
 *  This is the routine that actually calculates all of the quantites
 *  we care about.  Loop through grid cells and for all cells that are
 *  the most highly-resolved, we calculate all of the quantites shown
 *  at the top of the file.
 * 
 *---------------------------------------------------------------------*/
int GetDensityInfo(int gridnum){
  if(VERBOSEDEBUG) fprintf(stderr,"in GetDensityInfo\n");

  double cellsize,cellvol,cellmass, ccx,ccy,ccz, H2frac, efrac, HMfrac,
    overdensity,entropy,jeansmass,HIIfrac;

  int i,j,k,cellindex,totalcells, logodbin, logtempbin, logjeansbin, 
    logefracbin, logH2fracbin, logHminusfracbin,logentbin, logHIIfracbin;


  totalcells = griddx[gridnum]*griddy[gridnum]*griddz[gridnum];

  // size of grid cell, in code units
  cellsize = (gridrex[gridnum] - gridlex[gridnum]) / ((double) griddx[gridnum] );

  // cell volume in code units
  cellvol = cellsize*cellsize*cellsize;

  // density -> mass conversion, cell volume included
  massconv = densconstant * cellvol * boxsize * boxsize * boxsize 
    / hubble / hubble / hubble;
  
  // loop over all cells
  for(k=0; k < griddz[gridnum]; k++)
    for(j=0; j < griddy[gridnum]; j++)
      for(i=0; i < griddx[gridnum]; i++){

	// cell index (3d -> 1d conversion - USES FORTRAN ORDERING)
	cellindex = k*griddx[gridnum]*griddy[gridnum] + j*griddx[gridnum] + i;

	if(flagbuffer[cellindex] != -1){ // read just non-overlapped grid cells

	  // calculate center of cell (i,j,k)
	  ccx = gridlex[gridnum] + ( (gridrex[gridnum]-gridlex[gridnum]) / 
				     ((double) griddx[gridnum]) ) 
	    * ( ((double) i) + 0.5 );
	  ccy = gridley[gridnum] + ( (gridrey[gridnum]-gridley[gridnum]) / 
				     ((double) griddy[gridnum]) ) 
	    * ( ((double) j) + 0.5 );
	  ccz = gridlez[gridnum] + ( (gridrez[gridnum]-gridlez[gridnum]) / 
				     ((double) griddz[gridnum]) ) 
	    * ( ((double) k) + 0.5 );

	  // calculate mass of gas in cells
	  cellmass = ((double) densitybuff[cellindex])*cellvol;  // cell mass in code units

	  // calculate overdensity - this is what we're really interested
	  // in, typically.
	  overdensity = ((double) densitybuff[cellindex]) / (OMEGA_BARYON/OMEGA_MATTER);

	  // increment total mass, volume: these are the weights for
	  // most of the quantities below.
	  totalmass += cellmass;
	  totalvolume += cellvol;

	  meanoverdensity += overdensity*cellvol;
	  meanoverdensity_weight += cellvol;
	  


	} // if(flagbuffer[cellindex] != -1){

      } // end of triply-nested for loop

  if(VERBOSEDEBUG) fprintf(stderr,"leaving GetDensityInfo\n");
  return SUCCESS;
}



/*------------------------- GetVarianceInfo ----------------------------
 *
 * 
 *---------------------------------------------------------------------*/
int GetVarianceInfo(int gridnum){
  if(VERBOSEDEBUG) fprintf(stderr,"in GetVarianceInfo\n");

  double cellsize,cellvol,cellmass, ccx,ccy,ccz, H2frac, efrac, HMfrac,
    overdensity,entropy,jeansmass,HIIfrac, thisvariance;

  int i,j,k,cellindex,totalcells, logodbin, logtempbin, logjeansbin, 
    logefracbin, logH2fracbin, logHminusfracbin,logentbin, logHIIfracbin;


  totalcells = griddx[gridnum]*griddy[gridnum]*griddz[gridnum];

  // size of grid cell, in code units
  cellsize = (gridrex[gridnum] - gridlex[gridnum]) / ((double) griddx[gridnum] );

  // cell volume in code units
  cellvol = cellsize*cellsize*cellsize;

  // density -> mass conversion, cell volume included
  massconv = densconstant * cellvol * boxsize * boxsize * boxsize 
    / hubble / hubble / hubble;
  
  // loop over all cells
  for(k=0; k < griddz[gridnum]; k++)
    for(j=0; j < griddy[gridnum]; j++)
      for(i=0; i < griddx[gridnum]; i++){

	// cell index (3d -> 1d conversion - USES FORTRAN ORDERING)
	cellindex = k*griddx[gridnum]*griddy[gridnum] + j*griddx[gridnum] + i;

	if(flagbuffer[cellindex] != -1){ // read just non-overlapped grid cells

	  // calculate center of cell (i,j,k)
	  ccx = gridlex[gridnum] + ( (gridrex[gridnum]-gridlex[gridnum]) / 
				     ((double) griddx[gridnum]) ) 
	    * ( ((double) i) + 0.5 );
	  ccy = gridley[gridnum] + ( (gridrey[gridnum]-gridley[gridnum]) / 
				     ((double) griddy[gridnum]) ) 
	    * ( ((double) j) + 0.5 );
	  ccz = gridlez[gridnum] + ( (gridrez[gridnum]-gridlez[gridnum]) / 
				     ((double) griddz[gridnum]) ) 
	    * ( ((double) k) + 0.5 );

	  // calculate mass of gas in cells
	  cellmass = ((double) densitybuff[cellindex])*cellvol;  // cell mass in code units

	  // calculate overdensity - this is what we're really interested
	  // in, typically.
	  overdensity = ((double) densitybuff[cellindex]) / (OMEGA_BARYON/OMEGA_MATTER);

	  // increment total mass, volume: these are the weights for
	  // most of the quantities below.
	  totalmass += cellmass;
	  totalvolume += cellvol;

	  thisvariance = meanoverdensity - overdensity;

	  totalvariance += thisvariance*thisvariance*cellvol;

	  totalvariance_weight += cellvol;


	} // if(flagbuffer[cellindex] != -1){

      } // end of triply-nested for loop

  if(VERBOSEDEBUG) fprintf(stderr,"leaving GetVarianceInfo\n");
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

	  //if(flagbuffer[cellindex] != -1){
	  //  flagbuffer[cellindex]=-1;
	  //  flagged_cells_this_grid++;
	  //}
	  
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
      sscanf(line,"GridLeftEdge      = %lf %lf %lf",
	     &gridlex[grid],&gridley[grid],&gridlez[grid]);

      fgets(line, MAX_LINE_LENGTH, hierfile);  // right edge
      sscanf(line,"GridRightEdge     = %lf %lf %lf",
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
			      (1.0 / ((double) rootgridsize)) / 
			      ( (gridrex[grid]-gridlex[grid]) / 
				( (double) griddx[grid] ) )
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
 *  Reads the parameter file, extracts stuff needed in the rest of the program.
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
      sscanf(line,"CosmologyComovingBoxSize   = %lf",
	     &boxsize);
      if(DEBUG){
	fprintf(stderr,"found boxsize:  %s\n",line);
	fprintf(stderr,"box size is: %g\n",boxsize);
      }
    }
    
    // top grid dimension
    if(strncmp(line,"TopGridDimensions",17)==0)
      sscanf(line,"TopGridDimensions   = %d %d %d",
	     &rootgridsize,&ijunk,&ijunk);

    // multispecies
    if(strncmp(line,"MultiSpecies",12)==0)
      sscanf(line,"MultiSpecies                   = %d",
	     &multispecies);

    // omega lambda
    if(strncmp(line,"CosmologyOmegaLambdaNow",23)==0)
      sscanf(line,"CosmologyOmegaLambdaNow    = %lf",
	     &omegalambda);

    // omega matter
    if(strncmp(line,"CosmologyOmegaMatterNow",23)==0)
      sscanf(line,"CosmologyOmegaMatterNow    = %lf",
	     &omegamatter);

    // current redshift
    if(strncmp(line,"CosmologyCurrentRedshift",24)==0)
      sscanf(line,"CosmologyCurrentRedshift   = %lf",
	     &redshift);

    // current time
    if(strncmp(line,"InitialTime",11)==0)
      sscanf(line,"InitialTime         = %lf",
	     &currenttime);


    // hubble constant (in units of 100 km/s/Mpc)
    if(strncmp(line,"CosmologyHubbleConstantNow",26)==0)
      sscanf(line,"CosmologyHubbleConstantNow = %lf",
	     &hubble);


  } // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)

  fclose(headerfile);

  if(DEBUG) fprintf(stderr,"ReadParameterFile:  root grid is %d\n",rootgridsize);

  return 1;

}

// end of the file, hooray
