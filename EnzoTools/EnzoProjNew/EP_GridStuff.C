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

// declared in other files
extern int *gridlevel,*griddx,*griddy,*griddz,*gridnump,rootgridsize,particleonly;
extern double *gridlex,*gridley,*gridlez,
  *gridrex,*gridrey,*gridrez;
extern char **gridfilenames;

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

	if(VERBOSEDEBUG) fprintf(stderr,"GetGridInfo:  No baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

      } else {
	// "junk" lines
	for(i=0;i<2;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
	
	// grid file name
	sscanf(line,"BaryonFileName = %s",gridfilenames[grid]);

	// "junk" lines
	for(i=0;i<5;i++) fgets(line, MAX_LINE_LENGTH, hierfile);

	sscanf(line,"NumberOfParticles   = %d",
	       &gridnump[grid]);  // number of particles on this grid

	if(VERBOSEDEBUG) fprintf(stderr,"GetGridInfo: WITH baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

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
  if(DEBUG) fprintf(stderr,"about to open hierarchy file!\n");  
  hierfile = fopen(hierfilename,"r");
  if(DEBUG) fprintf(stderr,"just opened hierarchy file!\n");  
  // read through hierarchy file, counting # of grids
  // lines that start with "Grid = " are the beginning of a new
  // piece of grid info
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){
    if(VERBOSEDEBUG) fprintf(stderr,". ");
    if(strncmp(line,"Grid = ",7)==0) numgrids++;
    if(VERBOSEDEBUG) fprintf(stderr," - ");
  }

  if(DEBUG) fprintf(stderr,"finished reading file!\n");
    
  fclose(hierfile);

  if(DEBUG) fprintf(stderr,"just closed file!\n");

  if(DEBUG) fprintf(stderr,"NumberOfGrids:  there are %i grids\n",numgrids);
  
  // clean up dynamically allocated stuff
  delete [] line;

  if(DEBUG) fprintf(stderr,"exiting NumberOfGrids\n");

  // return # of grids
  return numgrids;  
}
