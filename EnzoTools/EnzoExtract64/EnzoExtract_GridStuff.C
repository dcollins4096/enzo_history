#include "EnzoExtract.h"

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
  int i,grid=0;
  char *line = new char[MAX_LINE_LENGTH];

  int griddimx,griddimy,griddimz,gridsex,gridsey,gridsez,
    grideex,grideey,grideez;

  if(DEBUG) fprintf(stderr,"GetGridInfo: about to open hierarchy file.\n");

  /* open the hierarchy file once, and count the number of grids.
     Then close the hierarchy file */

  
  if(DEBUG) fprintf(stderr,"hierfilename:  %s\n",hierfilename);

  hierfile=fopen(hierfilename,"r");

  if(DEBUG) fprintf(stderr,"GetGridInfo: file opened, counting grids\n");

  // count grids
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)
    if(strncmp(line,"Grid = ",7)==0) numberofgrids++;
  
  if(DEBUG) fprintf(stderr,"GetGridInfo: done counting, closing file.\n");

  fclose(hierfile);

  if(DEBUG) fprintf(stderr,"GetGridInfo: done counting grids, leaving. \n");

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


      
      // grid file name
#ifdef HIER_RICK
      for(i=0;i<3;i++) fgets(line, MAX_LINE_LENGTH, hierfile); // junk line

      sscanf(line,"FileName = %s",gridfilenames[grid]);
#else
      for(i=0;i<3;i++) fgets(line, MAX_LINE_LENGTH, hierfile);

      // "junk" lines
      for(i=0;i<2;i++) fgets(line, MAX_LINE_LENGTH, hierfile);

      sscanf(line,"BaryonFileName = %s",gridfilenames[grid]);
#endif
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
				     ) / log10((double) refineby)); //AK 
      
      if(VERBOSEDEBUG) fprintf(stderr,"%lf %lf %lf %lf %lf %lf %d %d %d %d FN: %s\n",
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


