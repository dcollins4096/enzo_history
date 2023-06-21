#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>

#define hssize_t hsize_t

#include "extern_hdf5.h"  // hdf 5 prototypes

#define MAX_LINE_LENGTH 256
#define SUCCESS 1
#define FAILURE 0



// ---------- USER DEFINES STUFF IN THIS SECTION ---------
#define DEBUG 1              // debug flag - 1 is on, 0 is off
#define VERBOSEDEBUG 0       // very verbose debug flat

#define OMEGA_BARYON 0.04
#define OMEGA_MATTER 0.3

#define BEGIN_POS  0.375
#define END_POS    0.625

#define DENSITY 1
#define VARIANCE 2

double meanoverdensity, meanoverdensity_weight,
  totalvariance, totalvariance_weight,
  totalmass, totalvolume;

#define RHOCRIT_CGS  1.8788e-29 // * hubble * hubble

#define HUGE_NUMBER 1.0e20 
#define HUGE_NEGATIVE_NUMBER 1.0e20

int particleonly=0;


// global hierarchy file values
int *gridlevel,*griddx,*griddy,*griddz,*gridnump;
double *gridlex,*gridley,*gridlez,
  *gridrex,*gridrey,*gridrez;
char **gridfilenames;
int *flagbuffer;


float *densitybuff;


// simulation parameters
double boxsize;
int rootgridsize, multispecies;
double mean_density,omegamatter,omegalambda,massconv,hubble,redshift,currenttime;

double densconstant, jeansmass_prefactor;


int numberofhalos,dark_matter_exists, total_number_grids;

int ReadParameterFile(char *filename);
int GetGridInfo(int numberofgrids,char *hierfilename);
int NumberOfGrids(char *hierfilename);
int GetCellInformationFromGrid(int gridnum,int total_number_grids, int whichroutine);
int FlagGridCells(int gridnum,int total_number_grids);
int GetDensityInfo(int gridnum);
int GetVarianceInfo(int gridnum);

void DeclareGridArrays(void);
void CleanGridArrays(void);

void OutputAllInformation(char *infilename);


