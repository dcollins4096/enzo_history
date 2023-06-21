#include "EnzoExtract.h"

// function declarations 
int FlagGridCells(int gridnum);
int AddGridCellsToExtractionVolume(int gridnum);
int AddSubVolumeToExtraction(int gridnum);

/*-------------------------- Extract ---------------------------
 *
 *  This function is basically the wrapper for the function that
 *  flags grid cells which have overlapping cells underneath them
 *  and the function that actually adds stuff to the extraction 
 *  volume.  Hopefully this promotes reusability...  Or something.
 *
 *--------------------------------------------------------------*/
int Extract(int gridnum){
  if(VERBOSEDEBUG) fprintf(stderr,"Extract:  entering\n");
  int i=0;

  // don't even look at this grid if it's too high of a level
  if(gridlevel[gridnum] > outputmaxlevel) return SUCCESS;

  // is any part of grid within bounds of projection?  if not, exit!
  if( !((gridrex[gridnum] > xstart) && (gridlex[gridnum] < xend) && 
	(gridrey[gridnum] > ystart) && (gridley[gridnum] < yend) && 
	(gridrez[gridnum] > zstart) && (gridlez[gridnum] < zend))
      ) return SUCCESS;

  // allocate buffer to flag cells
  numberofgridcells = griddx[gridnum]*griddy[gridnum]*griddz[gridnum];
  flagbuffer = new int[numberofgridcells];

  for(i=0; i<numberofgridcells; i++)
    flagbuffer[i]=0;

  // flag grid cells that have more highly refined cells under them,
  //  assuming the more highly refined cells are also above
  //  extractionmaxlevel.  Only flag cells which are less refined
  //  than the maximum extraction level, though!
  if(gridlevel[gridnum] < outputmaxlevel) FlagGridCells(gridnum);

  AddGridCellsToExtractionVolume(gridnum);

  // delete flagging buffer
  delete [] flagbuffer;

  if(VERBOSEDEBUG) fprintf(stderr,"Extract:  exiting\n");
  return SUCCESS;
}



/*------------- AddGridCellsToExtractionVolume -----------------
 *
 *  This function goes through the grid and opens each dataset 
 *  that we want to extract, adds it to a set of arrays, and then
 *  finally adds all of the stuff to the extraction arrays.
 *
 *--------------------------------------------------------------*/
int AddGridCellsToExtractionVolume(int gridnum){
  if(VERBOSEDEBUG) fprintf(stderr,"AddGridCellsToExtractionVolume %s\n",
			   gridfilenames[gridnum]);

  int i,ndims,havedimensions=0,field;

#ifdef FIELD_VALUES_DOUBLE
  gridbuffers = new double*[numberofextractionfields];
#else 
  gridbuffers = new float*[numberofextractionfields];
#endif

  // hdf 5 stuff
  hid_t   file_id, mem_type_id, mem_dsp_id, file_dsp_id,
    dsp_id,typ_id, grid_dset_id, group_id;
  hsize_t     size, dims[4],xdims[4],maxdims[4];
  herr_t      h5_status, h5_error = -1;


  // ---------------------------- SETUP -------------------------

  // open file 
  if(VERBOSEDEBUG) fprintf(stderr,"AddGridCells:  about to open file %s\n",gridfilenames[gridnum]);
  file_id = H5Fopen(gridfilenames[gridnum], H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );

#ifdef PACK_AMR
  //AK open group
  char id[10];
  sprintf(id, "%8.8d", gridnum+1);
  char name[MAX_LINE_LENGTH];
  strcpy(name, "/Grid");
  strcat(name, id);
  if(VERBOSEDEBUG) printf("H5Gopen with Name %s\n", name);
  group_id = H5Gopen(file_id, name);
  assert( group_id != h5_error );
#endif

  /* open up each dataset that we want extractions of and read that buffer
     into an array of arrays (gridbuffers[]) */
  for(field=0; field<numberofextractionfields; field++){

    // open dataset 
#ifdef PACK_AMR
    grid_dset_id = H5Dopen(group_id, ext_fieldnames[field]);//AK
#else
    grid_dset_id = H5Dopen(file_id, ext_fieldnames[field]);
#endif
    assert( grid_dset_id != h5_error );

    // we need to get information about the dimensionality,
    // but we only do it once.
    if(!havedimensions){
      havedimensions=1;

      // open dataspace (to get dimensions) 
      dsp_id = H5Dget_space(grid_dset_id);
      assert( dsp_id != h5_error );

      // get data type (only once!)
      typ_id = H5Dget_type(grid_dset_id);
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

    } // if(!havedimensions){

    // allocate field buffer
#ifdef FIELD_VALUES_DOUBLE
    gridbuffers[field] = new double[(int) size];
#else 
    gridbuffers[field] = new float[(int) size];
#endif


    if(DEBUG) fprintf(stderr,"****** SIZE %d  NUMGRIDCELLS %d ******\n",(int)size,numberofgridcells);
    assert( ((int) size) == numberofgridcells);

#ifdef FIELD_VALUES_DOUBLE
    mem_type_id = H5T_IEEE_F64BE;
#else 
    mem_type_id = H5T_IEEE_F32BE;
#endif

    // read buffer from file
    h5_status = H5Dread(grid_dset_id, mem_type_id, 
			mem_dsp_id, file_dsp_id, 
			H5P_DEFAULT, gridbuffers[field]);
    assert( h5_status != h5_error ); 


    // close stuff - this only gets closed on the last field
    if(field==(numberofextractionfields-1) ){
      h5_status = H5Sclose(dsp_id);
      assert( h5_status != h5_error );
      
      h5_status = H5Tclose(typ_id);
      assert( h5_status != h5_error );
      
      h5_status = H5Sclose(mem_dsp_id);
      assert( h5_status != h5_error );
      
      h5_status = H5Sclose(file_dsp_id);
      assert( h5_status != h5_error );
    }

    // close dataset
    h5_status = H5Dclose(grid_dset_id);
    assert( h5_status != h5_error );

  } // for(field=0; field<numberofextractionfields; field++){

#ifdef PACK_AMR
  h5_status = H5Gclose(group_id);  //AK
  assert( h5_status != h5_error );
#endif

  // close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  // -------------------- ACTUALLY DO STUFF -----------------

  AddSubVolumeToExtraction(gridnum);

  // --------------------- CLEANUP --------------------------

  for(field=0; field < numberofextractionfields; field++)
    delete [] gridbuffers[field];

  delete [] gridbuffers;

  return SUCCESS;
}



/*---------------- AddSubVolumeToExtraction --------------------
 *
 *
 *
 *
 *
 *------------------------------------------------------------*/
int AddSubVolumeToExtraction(int gridnum){
  int xbeg,ybeg,zbeg,xfin,yfin,zfin,i,j,k,
    gridindex_y,gridindex_z,gridindex_x,
    gridcellindex,extractcellindex,
    numcellswritten,field;

  double proj_delta,y_phys,z_phys,x_phys;

  /* calculate the beginning and ending cells in all three dimensions,
     for the extraction arrays.  We know that:
     (a) grid resolution is <= extraction resolution, and the difference
         is always a factor of two to some integar power.
     (b) at the maximum grid resolution allowed by the extraction, the 
         edges of the grid cells match up to the extraction cells.
  */
  xbeg=(int) ( (( max(xstart,gridlex[gridnum]) - xstart ) / (xend-xstart) ) * 
	       (((double) xextractnumcells)) );
  ybeg=(int) ( (( max(ystart,gridley[gridnum]) - ystart ) / (yend-ystart) ) * 
	       (((double) yextractnumcells)) );
  zbeg=(int) ( (( max(zstart,gridlez[gridnum]) - zstart ) / (zend-zstart) ) * 
	       (((double) zextractnumcells)) );

  xfin=(int) ( (( min(xend,gridrex[gridnum]) - xstart) / (xend-xstart) ) * 
	       (((double) xextractnumcells)-1.0) );
  yfin=(int) ( (( min(yend,gridrey[gridnum]) - ystart) / (yend-ystart) ) * 
	       (((double) yextractnumcells)-1.0) );
  zfin=(int) ( (( min(zend,gridrez[gridnum]) - zstart) / (zend-zstart) ) * 
	       (((double) zextractnumcells)-1.0) );

  // grid spacing of the extractions
  proj_delta = (xend - xstart) / ((double) xextractnumcells);

  if(VERBOSEDEBUG){
    fprintf(stderr,"AddSubVolumeToExtraction: %lf %i %i %i %i %i %i %i\n",
	    proj_delta,xbeg,ybeg,zbeg,xfin,yfin,zfin,gridnum);
    fprintf(stderr,"AddSubVolumeToExtraction: %lf %lf %lf %lf %lf %lf %i %i %i\n",
	    xstart,ystart,zstart,
	    gridlex[gridnum],gridley[gridnum],gridlez[gridnum],
	    xextractnumcells,yextractnumcells,zextractnumcells);
  }

  if( (xbeg > xfin) || (ybeg > yfin) || (zbeg > zfin) ){
    fprintf(stderr,"error in AddSubVolumeToExtraction: %i %i %i %i %i %i.  Exiting...\n",
		   xbeg,ybeg,zbeg,xfin,yfin,zfin);
    exit(FAILURE);
  }

  numcellswritten=0;

  for(i=xbeg; i <= xfin; i++)
    for(j=ybeg; j <= yfin; j++)
      for(k=zbeg; k <= zfin; k++){
	
	// calculate 'physical' center of projection array cell 
	x_phys = ( ((double) i) + 0.5 ) * proj_delta + xstart;
	y_phys = ( ((double) j) + 0.5 ) * proj_delta + ystart;
	z_phys = ( ((double) k) + 0.5 ) * proj_delta + zstart;
	
	// from that, calculate the x,y,z indices of the grid!
	gridindex_x= (int) (( (x_phys - gridlex[gridnum])/(gridrex[gridnum]-gridlex[gridnum]) ) * 
			    ((double) griddx[gridnum]) ); 
	gridindex_y= (int) (( (y_phys - gridley[gridnum])/(gridrey[gridnum]-gridley[gridnum]) ) * 
			    ((double) griddy[gridnum]) ); 
	gridindex_z= (int) (( (z_phys - gridlez[gridnum])/(gridrez[gridnum]-gridlez[gridnum]) ) * 
			    ((double) griddz[gridnum]) ); 

	extractcellindex =  k * yextractnumcells * xextractnumcells 
	  + j * xextractnumcells +i;
	
	gridcellindex = gridindex_z*griddx[gridnum]*griddy[gridnum] + 
	  gridindex_y*griddx[gridnum] + gridindex_x;
	
	if(flagbuffer[gridcellindex] != -1)
	  for(field=0; field<numberofextractionfields; field++)
	    ext_fieldvalues[field][extractcellindex] = gridbuffers[field][gridcellindex];

      }

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
 *  I'll even buy you a beer.
 *
 *--------------------------------------------------------------*/
int FlagGridCells(int gridnum){

  if(VERBOSEDEBUG) fprintf(stderr,"in FlagGridCells: %i %i\n",gridnum,numberofgrids);

  int counter, flagged_cells_this_grid,i,j,k,cellindex;

  double clex, cley, clez, crex, crey, crez;  // cell left and right edges (x,y,z)

  flagged_cells_this_grid=0;  // reset counter to zero

  /* loop over all grids and flag cells that are covered by a grid of
     higher resolution.  We can reject a lot of grids offhand because
     they are:
     1) not of level L+1 (where the grid of interest is level L) or
     2) outside of the bounds of the grid of interest
   */
  for(counter=0; counter < numberofgrids; counter++){
    
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


    // debug info
    if(VERBOSEDEBUG){
      fprintf(stderr,"FlagGridCells:  In Grid: %i (level %i)",
	      gridnum,gridlevel[gridnum]);
      fprintf(stderr,"Grid Looked At: %i (level %i)",counter,gridlevel[counter]);
      fprintf(stderr,"cells flagged this grid (ttl): %i\n",flagged_cells_this_grid);
      fprintf(stderr,"FlagGridCells: %lf %lf %lf   %lf %lf %lf\n",
	      gridlex[counter],gridley[counter],gridlez[counter],
	      gridrex[counter],gridrey[counter],gridrez[counter]);
    }

  }  // for(counter=0...

  if(VERBOSEDEBUG) fprintf(stderr,"exiting FlagGridCells...\n");
  return SUCCESS;
}
