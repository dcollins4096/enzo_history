#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#include "extern_hdf5.h"  // hdf 5 prototypes


/*-------------------------- USER DEFINES ----------------------*/

#define SN_XPOS 0.4941383507  // x,y,z position (in code units) of the supernova center
#define SN_YPOS 0.4969428387  // This is the point of max baryon density in the enzo
#define SN_ZPOS 0.4943933878  // calculations

#define SN_RADIUS 1.0         // in proper parsecs - GABE: set this!
  
char *amrparfilename = "RS0014";
char *amrhierfilename = "RS0014.hierarchy";

/*-------------------------------------------------------------*/

#define MAX_LINE_LENGTH 256

#define DEBUG 1         // 1 on, 0 off
#define VERBOSEDEBUG 1  // 1 on, 0 off

#define MPC_TO_CM 3.0857e+24  // megaparsecs in cgs
#define RHOCRIT 1.8788e-29    // critical density in cgs (times h^2 to get actual value)
#define GRAVC  6.67e-8        // gravitational constant in cgs
#define PI  3.1415926   
#define MPROTON 1.66e-24      // proton mass, cgs
#define MELECTRON 9.11e-28    // electron mass, cgs
#define MSOLAR 1.989e+33      // solar mass in cgs

// return calls
#define SUCCESS 1
#define FAILURE 0


int numberofgrids;

// pointers to arrays of info about the individual grids - this is all set
// in the routine GetGridInfo and then never changed
int   *gridlevel,
      *griddx,*griddy,*griddz;  // grid size (number of cells)

double *gridlex,*gridley,*gridlez,  // grid bottom left corner (code units)
       *gridrex,*gridrey,*gridrez;  // grid top right corner (code units)

char **gridfilenames;  // grid file names


// important simulation parameters - these are set in the routine
// ReadParameterFile and never changed
int rootgridsize;
double redshift, hubble, boxsize, initial_redshift,omegamatter;


// conversion factors from enzo code units to proper (NOT COMOVING)
// cgs units, double and float
double density_conversion, length_conversion, 
       velocity_conversion, time_conversion;

float density_conversion_float, length_conversion_float, 
      velocity_conversion_float, time_conversion_float;


// function declarations
void GetGridInfo(char *hierfilename);
void ReadParameterFile(char *filename);
void SetGridValues(int gridnumber);
void SetConversionFactors(void);
double diff(double x1, double y1, double z1, double x2, double y2, double z2);
