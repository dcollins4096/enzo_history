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
/*-----------------------------------------------------------
 *
 *                            EnzoProj
 *
 *  EnzoProj is an exciting and fun-filled routine that makes
 *  projections of enzo5 (the hdf5 version of enzo) data.  It
 *  is intended to fully replace the projection feature built
 *  into the code, ie, any routine that goes as "enzo5 -p ...".
 *
 *  This version is much more suited to large simulations, 
 *  distributed-memory architectures where you can't load the
 *  entire hierarchy into the memory on a single node, or really
 *  any situation where you're dealing with non-trivial datasets
 *  and would like to make a whole boatload of images in a short
 *  time.
 *
 *  One thing that is important to note is that all of the fields
 *  that enzoproj supports are always written out, regardless of
 *  whether the simulation in question actually HAS those fields.
 *  The 'empty' datasets are zeroed.  This is so analysis tools
 *  for IDL and other programs can be written simply.
 *
 *  Compilation Instructions:
 *    compile with gnu make (gmake on some systems):  
 *    "make -f MakeEnzoProj".
 *    you may have to tweak the makefile in order to find the hdf
 *    libraries.  Also, you have to change the top line of the 
 *    makefile to reflect the site you're at - this will hopefully
 *    change soon.  You will have to change the second parameter of the
 *    makefile to reflect the precision of the dataset you are 
 *    reading from (32 or 64 bit internal precision).  The third
 *    parameter of the makefile tells the code whether you want to use
 *    "packed" AMR or not.
 *
 *  Usage instructions:
 *    After compilation, type "enzoproj -h" or "enzoproj -help"
 *    for instructions, or see the enzo user's manual and cookbook.
 *    You must have the entire enzo hierarchy (all grid files, the
 *    hierarchy file and the parameter file) in the directory that
 *    you run enzoproj in (though of course the executable doesn't
 *    have to be in that directory).  The output file, enzo.project,
 *    will also appear in that directory.
 *
 *  Written by Brian W. O'Shea (bwoshea@cosmos.ucsd.edu)
 *    February 2003
 *
 *  Revision history:
 *
 *    July 2003  -  Finally got off my butt and added all sorts of
 *      slick new stuff.  Now all multispecies fields that are 
 *      currently functioning in enzo (e-, HI, HII, HeI, HeII, HeIII
 *      H-, H2I, H2II densities) are included.  Particle fields are
 *      now included - rho_dm and rho_star are both calculated from
 *      particle data.  Originally rho_dm was calculated from field
 *      data, and there were artifacts, which are now gone.  The code
 *      has also been cleaned up and made a bit more readable. -- BWO
 *
 *    June 2004 - Fixing EnzoProj so that it can make images of very
 *      high resolution simulations - essentially just modifying absolute
 *      positioning information (grid and projection pos. info) to be 64
 *      bit instead of 32 bit, which is how it should have been done from
 *      the get-go.  To-do: add some sort of cloud-in-cell smoothing
 *      to the dark matter density profile -- BWO
 *
 *    June 2006 - Rolled many changes into EnzoProj.  Output arrays are
 *      now row-major instead of column major, which makes outputs look
 *      right in IDL.  Can specify any output name we want.  Some
 *      debug flags fixed.  A few other things improved as well,
 *      but nothing major except for now particle and grid positions
 *      are all represented as 64-bit floats.
 *
 *    June 18, 2006 - Made a couple more changes in EnzoProj, primarily
 *      so that extremely high resolution calculations will now correctly
 *      generate projections.
 *
 *    August 2, 2006 - added Robert Harkness' changes for packed AMR
 *      and for the ability to read 32 vs. 64-bit precision in each
 *      dataset.  These are set in the makefile.
 *
 *    September 21, 2006 - added ability to read in MEKAL file tables,
 *      also added ability to compute spectral-weighted temperature.
 *      The usage of a MEKAL file table is set in the Makefile.
 *
 *--------------------------------------------------------------*/

#include "EnzoProj.h"

// global values from input arguments
int projaxis=0,projlevel=0,particleonly=0;
double xstart=0.0,ystart=0.0,zstart=0.0,
  xend=1.0,yend=1.0,zend=1.0;
char *inputfilename=NULL;
char *outputfilename=NULL;
char *MEKAL_Filename=NULL;
int outputfile_newname = 0;

// global values from parameter file
int maxlevel=-1,rootgridsize=-1,statichierarchy=-1,starformation=-1,
  multispecies=-1;
float omegamatter=-1.0,boxsize=-1.0,hubble=-1.0,redshift=-1.0,
  initialredshift=-1.0;

// global hierarchy file values
int *gridlevel,*griddx,*griddy,*griddz,*gridnump;
double *gridlex,*gridley,*gridlez,
  *gridrex,*gridrey,*gridrez;
char **gridfilenames;


/*------------------------ main() ------------------------------
 *
 *  main function - controls everything.  Not too much to say
 *  about that.  All of the really exciting (and sneaky) stuff
 *  takes place in the various subroutines (declared in 
 *  EnzoProj.h).
 *
 * -------------------------------------------------------------*/

int main(int argc, char *argv[]){
  int parse_return, param_file_return,total_number_grids,i,return_value;

#ifdef PACK_AMR
  printf("PACK_AMR is ON\n");
  fflush(stdout);
#else
  printf("PACK_AMR is OFF\n");
  fflush(stdout);
#endif

#ifdef R4
  printf("File precision is ___32 BIT___\n");
  fflush(stdout);
#endif

#ifdef R8
  printf("File precision is ___64 BIT___\n");
  fflush(stdout);
#endif

#ifdef USE_MEKAL
  printf("Using the MEKAL tables to calculate X-ray emissivity\n");
  fflush(stdout);
#endif // USE_MEKAL

  // look at command line arguments and figure out what the user wants
  parse_return = ParseArgs(argc,argv);

  // check return values of ParseArgs
  if(parse_return ==SUCCESS){  // normal, print out values for user (if debug flag on)
    if(DEBUG){
      if(projaxis==0) printf("Projection axis:        x\n");
      if(projaxis==1) printf("Projection axis:        y\n");
      if(projaxis==2) printf("Projection axis:        z\n");
      printf("Projection max level:   %i\n",projlevel);
      printf("Projection begin edge:  %lf  %lf  %lf\n",xstart,ystart,zstart);
      printf("Projection end edge:    %lf  %lf  %lf\n",xend,yend,zend);
      printf("Input file hierarchy:   %s\n",inputfilename);
      if(outputfile_newname==1)
	printf("Output file name: %s\n",outputfilename);
    }
  } else if(parse_return == -1){  // print help file, which exits.
    HelpMe();
  } else if(parse_return == -2){  // just exit
    exit(-1);
  } else{  // what the hell?
    fprintf(stderr,"function parse_return returned an unexpected value of %i.  Exiting.\n\n",parse_return);
    exit(-1);
  }

  // read the enzo parameter file to get important info about simulation
  param_file_return = ReadParameterFile(inputfilename);

  // check return values
  if(param_file_return==SUCCESS){  // everything's good
    if(DEBUG){
      printf("Maximum Refinement Level:        %i\n",maxlevel);
      printf("Static Hierarchy?:               %i\n",statichierarchy);
      printf("Root Grid # of cells:            %i\n",rootgridsize);
      printf("Box Size (Mpc/h):                %f\n",boxsize);
      printf("Omega Matter:                    %f\n",omegamatter);
      printf("Hubble constant (100 km/s/mpc):  %f\n",hubble);
      printf("Redshift:                        %f\n",redshift);
      printf("MultiSpecies:                    %i\n",multispecies);
    }
  } else if(param_file_return == -1){  // bad value in parameter file, dump and exit.
    fprintf(stderr,"read_parameter_file:  bad value in %s\n",inputfilename);
    fprintf(stderr,"Maximum Refinement Level:        %i\n",maxlevel);
    fprintf(stderr,"Static Hierarchy?:               %i\n",statichierarchy);
    fprintf(stderr,"Root Grid # of cells:            %i\n",rootgridsize);
    fprintf(stderr,"Box Size (Mpc/h):                %f\n",boxsize);
    fprintf(stderr,"Omega Matter:                    %f\n",omegamatter);
    fprintf(stderr,"Hubble constant (100 km/s/mpc):  %f\n",hubble);
    fprintf(stderr,"Redshift:                        %f\n",redshift);
    fprintf(stderr,"Exiting...\n\n");
    exit(-1);
  } else{  // something's seriously screwed - shouldn't ever see this!
    fprintf(stderr,"read_parameter_file returned a bad value of %i.\n",param_file_return);
    fprintf(stderr,"Exiting...\n\n");
    exit(-1);
  }

  // check user-defined boundary values, just in case
  return_value = CheckBoundaryValues();
  assert( return_value != FAILURE );

#ifdef USE_MEKAL
  return_value = ReadMEKALTable(MEKAL_Filename);
  assert( return_value != FAILURE );
#endif // USE_MEKAL

  // get hierarchy file name
  int inlen = (int) strlen(inputfilename);
  fprintf(stderr,"inlen is: %i\n",inlen);
  char *hierfile = new char[inlen+6];  // no fucking clue why this works
  strcpy(hierfile,inputfilename);
  hierfile=strcat(hierfile,".hierarchy");

  if(DEBUG) fprintf(stderr,"hierarchy name is %s\n",hierfile);

  // get total number of grids - to generate grid information
  total_number_grids = NumberOfGrids(hierfile);
  assert( total_number_grids != FAILURE );

  if(DEBUG) fprintf(stderr,"there are %i grids\n",total_number_grids);

  /* initialize arrays which grid information is stored in
     (information such as level, grid bounds and number of
     cells, etc.) */
  gridlevel = new int[total_number_grids];  // level info

  // grid number of cells info
  griddx = new int[total_number_grids];  
  griddy = new int[total_number_grids];
  griddz = new int[total_number_grids];

  // grid bounds info -- make double precision for high-rez calcs (bwo, june '04)
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

  // create arrays for projection based on begin/end bounds, proj axis and # of levels
  return_value = CreateProjectionArrays();
  assert( return_value != FAILURE );


  // do this stuff only if there are baryon files!
  if(!particleonly){


#ifdef PACK_AMR
    fprintf(stderr,"PACK_AMR is on - not looking for metals right now\n");
#else
    fprintf(stderr,"PACK_AMR is OFF - looking for metals right now\n");
    // figure out which metal fields are included (if any) and flag for them.
    return_value = CheckForMetals();
    assert( return_value != FAILURE );
#endif


    // iterate through grids, skipping grids that are > maxlevel
    //   flag cells that have grids below them that are <= maxlevel!
    //   read in all data along appropriate axis - calculate quantities as nec.
    for(i=0; i<total_number_grids; i++){
      return_value=AddGridToProjection(i,total_number_grids);
      assert( return_value != FAILURE );

    }

  }

  // Add particles to projection arrays - this includes both dark matter
  // and star particles (if star formation is on).  AddParticlesToProjection
  // goes through ALL grids, since particles are kept on the highest
  // resolution grid for a given piece of space.
  for(i=0; i<total_number_grids; i++){
    return_value=AddParticlesToProjection(i);
    assert( return_value != FAILURE );
  }

  // normalize some of the projection arrays
  return_value=NormalizeProjectionArrays();
  assert( return_value != FAILURE );

  // write out projections to file
  return_value=WriteProjectionArrays();
  assert( return_value != FAILURE );

  // smooth dark matter density projection here - TO DO!


  // delete projection arrays (dynamic memory alloc. cleanup)
  return_value = DeleteProjectionArrays();
  assert( return_value != FAILURE );

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

  if(DEBUG) fprintf(stderr,"exiting enzoproj...\n");

  exit(0);
}

