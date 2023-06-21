#include "EnzoExtract.h"

/* --------------------------- ParseInputs --------------------------
 *  This routine reads in the command line inputs and makes stores
 *  values in appropriate variables.  Some error checking too.
 * ------------------------------------------------------------------ */
void ParseInputs(int numarg,char *arguments[]){
  if(DEBUG) fprintf(stderr,"ParseInputs:  entering\n");
  int i;
  // set default values
  xstart=ystart=zstart=0.0;
  xend=yend=zend=1.0;
  outputrowmajor=0;
  outputmaxlevel=0;
  refineby = 2;
  parameterfile=NULL;


  // no arguments?  print command list and exit.
  if(numarg < 2){
    HelpMe();
    exit(FAILURE);
  } 
  // one argument?  if it's help request, print command list.
  // otherwise it's the input filename, and get that.
  else if(numarg==2){
    if(strcmp(arguments[1],"-h")==0 
       || strcmp(arguments[1],"-help")==0){ 
      HelpMe();
      exit(SUCCESS);
    } else{
      parameterfile=arguments[1];
    }
  } 
  // more than one argument?  parse the arguments!
  else {
    for(i=1; i<numarg-1;i++){

    if(strcmp(arguments[i],"-l")==0){  // level
      i++;
      outputmaxlevel=atoi(arguments[i]);
      if(outputmaxlevel < 0){
	fprintf(stderr,"max level of extraction must be >= 0!\n");
	fprintf(stderr,"value is:  %i\nExiting...\n\n",outputmaxlevel);
	exit(FAILURE);
      }
    } else if(strcmp(arguments[i],"-b")==0){  // beginning corner (bottom left front)
      
      xstart=(double) atof(arguments[++i]);
      ystart=(double) atof(arguments[++i]);
      zstart=(double) atof(arguments[++i]);
      
    }  else if(strcmp(arguments[i],"-f")==0){  // final corner (top left back)

      xend=(double) atof(arguments[++i]);
      yend=(double) atof(arguments[++i]);
      zend=(double) atof(arguments[++i]);
    } else if(strcmp(arguments[i],"-row")==0){ // row major outputs?
      outputrowmajor=1;
    } else if(strcmp(arguments[i],"-h")==0 || strcmp(arguments[i],"-help")==0){
      fprintf(stderr,"-h or -help %s\n",arguments[i]);
      exit(FAILURE);
    } else{
      fprintf(stderr,"I don't recognize this:  %s\n",arguments[i]);
      fprintf(stderr,"type enzoextract -h for the help file.\nExiting...\n\n");
      exit(FAILURE);
    }

    } // for(i=1; i<numarg-1;i++){
  } 
  
  // hierarchy name is the last argument - if it isn't the program
  // will crash somewhere else in the code.
  parameterfile=arguments[numarg-1];

  if(DEBUG){
    fprintf(stderr,"ParseInputs:  Input values are:\n");
    fprintf(stderr,"\txstart:  %lf   ystart:  %lf   zstart:  %lf\n",xstart,ystart,zstart);
    fprintf(stderr,"\txend:  %lf   yend:  %lf   zend:  %lf\n",xend,yend,zend);
    fprintf(stderr,"\tmaxlevel:   %d\n",outputmaxlevel);
    fprintf(stderr,"\toutputrowmajor:  %d\n",outputrowmajor);
    fprintf(stderr,"\tfile name:  %s\n",parameterfile);
  }

  //fprintf(stderr,"parsinputs:  %s\n",parameterfile);
  if(DEBUG) fprintf(stderr,"ParseInputs:  exiting\n");

}



/* --------------------------- ReadParameterFile --------------------
 *  Parse parameter file looking for the top grid size and the maximum
 *  level.  For the purposes of an extraction we don't actually care
 *  about any cosmological parameters, so we don't get those.
 * ------------------------------------------------------------------ */
void ReadParameterFile(char *filename){
  if(DEBUG) fprintf(stderr,"ReadParameterFile:  Reading %s\n",filename);
  FILE *headerfile;
  //char *line = new char[MAX_LINE_LENGTH];
  char line[MAX_LINE_LENGTH];
  int ijunk;

  rootgridsize = maxlevel = -1;

  // open up header file and go through line by line
  headerfile=fopen(filename,"r");

  while( fgets(line, MAX_LINE_LENGTH, headerfile) != NULL){
    
    // top grid dimension
    if(strncmp(line,"TopGridDimensions",17)==0)
      sscanf(line,"TopGridDimensions   = %d %d %d",
	     &rootgridsize,&ijunk,&ijunk);

    // maximum refinement level
    if(strncmp(line,"MaximumRefinementLevel",22)==0)
      sscanf(line,"MaximumRefinementLevel = %d",
	     &maxlevel);

    // refinement factor //AK
    if(strncmp(line,"RefineBy",8)==0)
      sscanf(line,"RefineBy = %d",
	     &refineby);

  } // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)

  fclose(headerfile);

  if( (rootgridsize==-1) || (maxlevel == -1) ){
    fprintf(stderr,"ReadParameterFile:  Something wrong!  %d %d\n",
	    rootgridsize,maxlevel);
    exit(FAILURE);
  }

  //delete [] line;

  if(DEBUG) fprintf(stderr,"ReadParameterFile:  root grid size: %d\n",rootgridsize);
  if(DEBUG) fprintf(stderr,"ReadParameterFile:  maximum level: %d\n",maxlevel);
  if(DEBUG) fprintf(stderr,"ReadParameterFile:  exiting\n");
}



/* --------------------------- CheckBoundaryValues -------------------
 *  Checks boundary information to make sure that:
 *    a)  starting and ending indices are within bounds
 *    b)  starting indices < corresponding end indices
 *    c)  starting and ending indices line up with a grid edge (on the
 *        highest level of extraction).
 *
 *  If (a) or (b) are not true, exit.  If (c) is not true, we just correct
 *  so that the start and end indices line up with a grid edge and print
 *  out a message informing the user.
 * ------------------------------------------------------------------ */
void CheckBoundaryValues(void){
  if(DEBUG) fprintf(stderr,"CheckBoundaryValues:  entering\n");

  double delta_grid,totalcells,xstartorig,ystartorig,zstartorig,
    xendorig,yendorig,zendorig;
  int tempcellnum;

  // are start or end indices out of bounds?  If so, print an error message and exit.
  if( (xstart < 0.0) || (ystart < 0.0) || (zstart < 0.0) || 
      (xend > 1.0) || (yend > 1.0) || (zend > 1.0) ||
      (xstart >= xend) || (ystart > yend) || (zstart > zend) ){
    fprintf(stderr,"CheckBoundaryValues: ERROR. start and end indices must be\n");
    fprintf(stderr,"between 0.0 and 1.0, and end indices must be larger than\n");
    fprintf(stderr,"start indices.  Your values are:\n");
    fprintf(stderr,"\txstart: %lf\n",xstart);
    fprintf(stderr,"\tystart: %lf\n",ystart);
    fprintf(stderr,"\tzstart: %lf\n",zstart);
    fprintf(stderr,"\txend: %lf\n",xend);
    fprintf(stderr,"\tyend: %lf\n",yend);
    fprintf(stderr,"\tzend: %lf\n",zend);
    fprintf(stderr,"\nexiting...\n");
    exit(SUCCESS);
  }

  // need to check boundary values here to make sure that they line up with edges!

  // calculate the grid spacing of the highest level of refinement
  // in the projection
  delta_grid = 1.0 / ( (double) rootgridsize ) / 
    ( (double) pow( (double) refineby, ((double) outputmaxlevel) ) );

  // calculate the total number of grid cells in the cube
  // it's float instead of int in order to make the rest of our
  // calculations straightforward
  totalcells = (double) rootgridsize * ( (double) pow( (double) refineby, ((double) outputmaxlevel) ) );//AK

  // save original boundary values;
  xstartorig = xstart;
  ystartorig = ystart;
  zstartorig = zstart;
  xendorig = xend;
  yendorig = yend;
  zendorig = zend;
  
  // x axis values
  if(xstart != 0.0){

    // calculate the value of the cell that xstart should be in
    // force it to be an integar!
    tempcellnum = ((int) (xstart * totalcells));

    // calculate xstart.  since tempcellnum was forced to be an
    // integar, we now can be assured that xstart is exactly at
    // at the boundary of the most highly refined level of
    // projection
    xstart = ((double) tempcellnum ) / totalcells;
  }

  if(xend != 1.0){
    tempcellnum = ((int) (xend * totalcells));
    xend = ((double) tempcellnum ) / totalcells;
    // to ensure that we have at least one grid cell
    // thickness
    if(xend==xstart) xend += delta_grid;
  }

  // y axis values - see x stuff for explanation
  if(ystart != 0.0){
    tempcellnum = ((int) (ystart * totalcells));
    ystart = ((double) tempcellnum ) / totalcells;
  }

  if(yend != 1.0){
    tempcellnum = ((int) (yend * totalcells));
    yend = ((double) tempcellnum ) / totalcells;
    if(yend==ystart) yend += delta_grid;
  }

  // z axis values - see x stuff for explanation
  if(zstart != 0.0){
    tempcellnum = ((int) (zstart * totalcells));
    zstart = ((double) tempcellnum ) / totalcells;
  }

  if(zend != 1.0){
    tempcellnum = ((int) (zend * totalcells));
    zend = ((double) tempcellnum ) / totalcells;
    if(zend==zstart) zend += delta_grid;
  }

  // if the value WAS changed, print out so that the user knows!
  if((xstartorig != xstart))
    fprintf(stderr,"CheckBoundaryValues: xstart was %lf and is now %lf\n",
	    xstartorig,xstart);

  if((xendorig != xend))
    fprintf(stderr,"CheckBoundaryValues: xend was %lf and is now %lf\n",
	    xendorig,xend);

  if((ystartorig != ystart))
    fprintf(stderr,"CheckBoundaryValues: ystart was %lf and is now %lf\n",
	    ystartorig,ystart);

  if((yendorig != yend))
    fprintf(stderr,"CheckBoundaryValues: yend was %lf and is now %lf\n",
	    yendorig,yend);

  if((zstartorig != zstart))
    fprintf(stderr,"CheckBoundaryValues: zstart was %lf and is now %lf\n",
	    zstartorig,zstart);

  if((zendorig != zend))
    fprintf(stderr,"CheckBoundaryValues: zend was %lf and is now %lf\n",
	    zendorig,zend);

  if(DEBUG) fprintf(stderr,"CheckBoundaryValues:  exiting\n");
} 




/* --------------------------- HelpMe -------------------------------
 *  Prints out helpful usage information to the user and exits.
 * ------------------------------------------------------------------ */
void HelpMe(void){

  fprintf(stderr,"enzoextract -- enzo extraction utility\n\n");
  fprintf(stderr,"Usage:  enzoextract [-b <x> <y> <z>] [-f <x> <y> <z>] [-row] [-l <level>] <hierarchy name>\n");
  fprintf(stderr,"Options:\n\n");
  fprintf(stderr,"\t-b    begin flag\n");
  fprintf(stderr,"\t\t\tvalid args:  0.0 - 1.0 per axis\n");
  fprintf(stderr,"\t\t\tdefault args: 0.0 0.0 0.0\n");
  fprintf(stderr,"\t-f    finish flag\n");
  fprintf(stderr,"\t\t\tvalid args:  0.0 - 1.0 per axis\n");
  fprintf(stderr,"\t\t\tdefault argument: 1.0 1.0 1.0\n");
  fprintf(stderr,"\t-l    max projection level\n");
  fprintf(stderr,"\t\t\tvalid arguments:  0 - Maxlevel\n");
  fprintf(stderr,"\t\t\tdefault argument: 0\n");
  fprintf(stderr,"\t-row  output arrays in row-major format\n");
  fprintf(stderr,"\t-h or -help    print help message\n\n");
  fprintf(stderr,"Example:  enzoextract -b 0.3 0.3 0.4 -f 0.5 0.6 0.6 -l 2  RedshiftOutput0002\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"****** COMPILATION NOTES ******\n");
#ifdef PACK_AMR
  fprintf(stderr,"PACK_AMR is           ON.\n");
#else
  fprintf(stderr,"PACK_AMR is           OFF.\n");
#endif

#ifdef FIELD_VALUES_DOUBLE
  fprintf(stderr,"64-bit float IO is    ON\n");
#else 
  fprintf(stderr,"32-bit float IO is    ON\n");
#endif

#ifdef DEBUG
  fprintf(stderr,"DEBUG flag is         ON\n");
#else
  fprintf(stderr,"DEBUG flag is         OFF\n");
#endif

#ifdef VERBOSEDEBUG
  fprintf(stderr,"VERBOSEDEBUG flag is  ON\n");
#else
  fprintf(stderr,"VERBOSEDEBUG flag is  OFF\n");
#endif


#ifdef HIER_RICK
  fprintf(stderr,"Hierarchy file is     RICK\n");
#else
  fprintf(stderr,"Hierarchy file is     STANDARD\n");
#endif

  fprintf(stderr,"\nYou will have to recompile if you want these changed.\n");
  fprintf(stderr,"****************\n\n");

  exit(SUCCESS);

}
