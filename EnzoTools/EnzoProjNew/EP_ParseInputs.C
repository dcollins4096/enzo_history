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

// external variables declared in other files!
extern int projaxis,projlevel,maxlevel,rootgridsize,
  statichierarchy,starformation,multispecies,outputfile_newname;
extern double xstart,ystart,zstart,xend,yend,zend;
extern float omegamatter, boxsize,hubble,redshift,initialredshift;
extern char *inputfilename, *outputfilename, *MEKAL_Filename;


/*---------------------- ParseArgs ----------------------------
 *
 *  This function (called from main()) parses the command-line
 *  arguments handed in by the user.  If the user makes a
 *  mistake then various return codes go back - otherwise,
 *  return SUCCESS.
 *
 *-------------------------------------------------------------*/
int ParseArgs(int numarg,char *arguments[]){
  int i;
  if(DEBUG) fprintf(stderr,"ParseArgs:  entering\n");
  // if there aren't enough arguments display help file and exit
  if(numarg < 2) return -1;

  /* if there is only one argument, it's either a request for the
     help file or the hierarchy name.  Display help file and exit,
     or store hierarchy name and return to main */
  if(numarg==2){
    if(strcmp(arguments[1],"-h")==0 
       || strcmp(arguments[1],"-help")==0){ 
      return -1;
    } else{
      inputfilename=arguments[1];
      return SUCCESS;
    }
  }

  /* if there is more than one argument, loop through arguments,
     assuming that the last one is the file name.  Error check
     as appropriate. */
  for(i=1; i<numarg-1;i++){
    if(strcmp(arguments[i],"-l")==0){  // level
      i++;
      projlevel=atoi(arguments[i]);
      if(projlevel < 0){
	fprintf(stderr,"max level of projection must be >= 0!\n");
	fprintf(stderr,"value is:  %i\nExiting...\n\n",projlevel);
	return -2;
      }
    } else if(strcmp(arguments[i],"-p")==0){  // projection axis
      i++;
      if(strcmp(arguments[i],"x")==0){ 
	projaxis=0;
      } else if(strcmp(arguments[i],"y")==0){
	projaxis = 1;
      } else if(strcmp(arguments[i],"z")==0){
	projaxis = 2;
      } else{
	fprintf(stderr,"I don't know what this is (after -p): %s\n",arguments[i]);
	fprintf(stderr,"Type enzoproj -h for help file.\nExiting...\n\n");
	return -2;
      }
    } else if(strcmp(arguments[i],"-b")==0){  // beginning corner (bottom left front)
 
#ifdef NEW_INPUT_STYLE     
      xstart= strtod(arguments[++i], NULL);
      ystart= strtod(arguments[++i], NULL);
      zstart= strtod(arguments[++i], NULL);

      if(DEBUG) fprintf(stderr," -b %.16lf %.16lf %.16lf\n", 
			xstart, ystart, zstart);
#else
      xstart= (double) atof(arguments[++i]); 
      ystart= (double) atof(arguments[++i]);
      zstart= (double) atof(arguments[++i]);
#endif

      if(xstart < 0.0 || xstart >= 1.0 || ystart < 0.0 || ystart >= 1.0
	 || zstart < 0.0 || zstart >= 1.0){
	fprintf(stderr,"arguments following -b are out of bounds:  %lf %lf %lf\n",
		xstart,ystart,zstart);
	fprintf(stderr,"they must be between 0.0 and 1.0, inclusive.\nExiting...\n\n");
	return -2;
      }
      
    } else if(strcmp(arguments[i],"-f")==0){  // final corner (top left back)

#ifdef NEW_INPUT_STYLE
      xend= strtod(arguments[++i], NULL);
      yend= strtod(arguments[++i], NULL);
      zend= strtod(arguments[++i], NULL);

      if(DEBUG) fprintf(stderr," -f %.16lf %.16lf %.16lf\n", 
			xend, yend, zend);

#else
      xend= (double) atof(arguments[++i]); 
      yend= (double) atof(arguments[++i]);
      zend= (double) atof(arguments[++i]);
#endif

      if(xend <= 0.0 || xend > 1.0 || yend <= 0.0 || yend > 1.0
	 || zend <= 0.0 || zend > 1.0){
	fprintf(stderr,"arguments following -f are out of bounds:  %lf %lf %lf\n",
		xend,yend,zend);
	fprintf(stderr,"they must be between 0.0 and 1.0, inclusive.\nExiting...\n\n");
	return -2;
       }
    } else if(strcmp(arguments[i],"-h")==0 || strcmp(arguments[i],"-help")==0){
      fprintf(stderr,"-h or -help %s\n",arguments[i]);
      return -1;
    } else if(strcmp(arguments[i],"-o")==0){ // output file name
      outputfilename = arguments[++i];
      outputfile_newname = 1;
      
    } else if(strcmp(arguments[i],"-m")==0){ // MEKAL file name
      MEKAL_Filename = arguments[++i];
    } else{
      fprintf(stderr,"I don't recognize this:  %s\n",arguments[i]);
      fprintf(stderr,"type enzoproj -h for the help file.\nExiting...\n\n");
      return -2;
    }

  } // for(i=1; i<numarg-2;i++){

  // hierarchy name is the last argument - if it isn't the program
  // will crash somewhere else in the code.
  inputfilename=arguments[numarg-1];

  // error check boundaries
  if(xstart >= xend || ystart >= yend || zstart >= zend){
    fprintf(stderr,"Check your boundaries!  Start corner cannot\n");
    fprintf(stderr,"be greater than end corner!\n");
    fprintf(stderr,"start:  %lf %lf %lf\n",xstart,ystart,zstart);
    fprintf(stderr,"end:    %lf %lf %lf\n",xend,yend,zend);
    fprintf(stderr,"exiting...\n\n");
    return -2;
  }

  if(DEBUG) fprintf(stderr,"ParseArgs:  leaving successfully\n");
  return SUCCESS;  
}


/*--------------------- ReadParameterFile() ---------------------
 *
 *  This function reads in the enzo parameter file and extracts
 *  a bunch of information about the simulation: namely, root
 *  grid information, level information, and various cosmological
 *  parameters that are needed in other parts of the simulation.
 *
 *---------------------------------------------------------------*/
int ReadParameterFile(char *filename){
  FILE *headerfile;
  char *line = new char[MAX_LINE_LENGTH];
  int ijunk;

  if(DEBUG) printf("ReadParameterFile:  reading %s\n",filename);

  // open up header file and go through line by line
  headerfile=fopen(filename,"r");

  while( fgets(line, MAX_LINE_LENGTH, headerfile) != NULL){
    
    // top grid dimension
    if(strncmp(line,"TopGridDimensions",17)==0)
      sscanf(line,"TopGridDimensions   = %d %d %d",
	     &rootgridsize,&ijunk,&ijunk);

    // current redshift
    if(strncmp(line,"CosmologyCurrentRedshift",24)==0)
      sscanf(line,"CosmologyCurrentRedshift = %f",
	     &redshift);

    // current redshift
    if(strncmp(line,"CosmologyInitialRedshift",24)==0)
      sscanf(line,"CosmologyInitialRedshift = %f",
	     &initialredshift);

    // comoving box size (in Mpc/h)
    if(strncmp(line,"CosmologyComovingBoxSize",24)==0)
      sscanf(line,"CosmologyComovingBoxSize   = %f",
	     &boxsize);

    // hubble constant (in units of 100 km/s/Mpc)
    if(strncmp(line,"CosmologyHubbleConstantNow",26)==0)
      sscanf(line,"CosmologyHubbleConstantNow = %f",
	     &hubble);

    // omega matter
    if(strncmp(line,"CosmologyOmegaMatterNow",23)==0)
      sscanf(line,"CosmologyOmegaMatterNow    = %f",
	     &omegamatter);

    // maximum refinement level
    if(strncmp(line,"MaximumRefinementLevel",22)==0)
      sscanf(line,"MaximumRefinementLevel = %d",
	     &maxlevel);

    // maximum refinement level
    if(strncmp(line,"StaticHierarchy",15)==0)
      sscanf(line,"StaticHierarchy     = %d",
	     &statichierarchy);

    // multispecies
    if(strncmp(line,"MultiSpecies",12)==0)
      sscanf(line,"MultiSpecies                   = %d",
	     &multispecies);


    if(strncmp(line,"StarParticleCreation",20)==0){
      sscanf(line,"StarParticleCreation                  = %d",&ijunk);
      if(ijunk != 0){ 
	starformation = 1;
	if(DEBUG) fprintf(stderr,"ReadParameterFile:  Star formation is ON\n");
      } else starformation = 0;
    }


  } // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)

  fclose(headerfile);

  // make sure that all of the values are sensical
  if(maxlevel >= 0 && redshift >= 0.0
     && rootgridsize > 1 && omegamatter >= 0.0
     && boxsize > 0.0 && hubble > 0.0
     && statichierarchy >= 0){
    if(DEBUG) fprintf(stderr,"ReadParameterFile:  exiting successfully\n");
    return SUCCESS;
  } else{
    return -1;
  }
}


/*----------------------- HelpMe()  ----------------------------
 *
 *  This function prints out all of the help information, along 
 *  with an example call to the program.  Any call to HelpMe 
 *  automatically exits the program.
 *
 *--------------------------------------------------------------*/
void HelpMe(void){
  fprintf(stderr,"enzoproj - enzo piecewise projection tool\n\n");
  fprintf(stderr,"Usage:  enzoproj [-p <axis>] [-b <x1> <y1> <z1>] [-f <x2> <y2> <z2>] [-l <level>] <hierarchy name>\n\n");
  fprintf(stderr,"Options:\n\n");
  fprintf(stderr,"\t-p    projection axis\n");
  fprintf(stderr,"\t\t\tvalid arguments:  x,y,z\n");
  fprintf(stderr,"\t\t\tdefault argument: x\n");
  fprintf(stderr,"\t-l    max projection level\n");
  fprintf(stderr,"\t\t\tvalid arguments:  0 - Maxlevel\n");
  fprintf(stderr,"\t\t\tdefault argument: 0\n");
  fprintf(stderr,"\t-b    begin flag\n");
  fprintf(stderr,"\t\t\tvalid args:  0.0 - 1.0 per axis\n");
  fprintf(stderr,"\t\t\tdefault args: 0.0 0.0 0.0\n");
  fprintf(stderr,"\t-f    finish flag\n");
  fprintf(stderr,"\t\t\tvalid args:  0.0 - 1.0 per axis\n");
  fprintf(stderr,"\t\t\tdefault argument: 1.0 1.0 1.0\n");
  fprintf(stderr,"\t-o <filename> outputfile name\n");
  fprintf(stderr,"\t-m <MEKALfile> MEKAL table file name (x-ray emissivities)\n");
  fprintf(stderr,"\t\t\tdefault file name: enzo.project\n");
  fprintf(stderr,"\t-h or -help    print help message\n\n");
  fprintf(stderr,"Example:  enzoproj -p y  -l 2  -b 0.3 0.3 0.4  -f 0.5 0.6 0.5 -m datafile.mekal  -o test.project RedshiftOutput0002\n");
  fprintf(stderr,"\n");

  exit(0);
}
