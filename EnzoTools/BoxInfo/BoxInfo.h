#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>

//#define hssize_t hsize_t

#include "extern_hdf5.h"  // hdf 5 prototypes

#define MAX_LINE_LENGTH 256
#define SUCCESS 1
#define FAILURE 0



// ---------- USER DEFINES STUFF IN THIS SECTION ---------
#define DEBUG 0              // debug flag - 1 is on, 0 is off
#define VERBOSEDEBUG 0       // very verbose debug flat

#define USE_METALS   // set to USE_METALS or NO_USE_METALS

#define NUMBINS 300          // total number of bins in 1D, 2D arrays

#define OMEGA_BARYON 0.0418
#define OMEGA_MATTER 0.2383

// these are for the 1D and 2D distribution functions: log bins for
// all of the various quantities we're interested in (log base 10)
#define SMIN 11.0  // entropy 
#define SMAX 21.0
#define DMIN -2.0  // overdensity
#define DMAX 7.0
#define TMIN 0.0   // temperature
#define TMAX 7.0
#define H2MIN -6.0  // molecular hydrogen fraction
#define H2MAX 0.0
#define HIIMIN -6.0  // Hplus fraction
#define HIIMAX 0.0
#define EMIN -11.0  // electron fraction
#define EMAX -1.0
#define HMMIN -20.0 // Hminus fraction
#define HMMAX -9.0
#define JMIN 0.0    // jeans mass
#define JMAX 8.0

// ------------------------------------------------------


double maxbaryondens,
  maxbar_temp, maxbar_H2frac,
  maxbar_HIIfrac,
  maxbar_xpos, maxbar_ypos, maxbar_zpos,
  rhobar_mean_mw, rhobarsquared_mean_mw,
  rhobar_mean_vw, rhobarsquared_mean_vw,
  clumping_factor_mw, clumping_factor_vw,
  meantemp_mw, meantemp_vw,
  meanent_mw, meanent_vw,
  meanH2frac_mw, meanH2frac_vw,
  meanHIIfrac_mw, meanHIIfrac_vw,
  meanelecfrac_mw, meanelecfrac_vw,
  meanHMfrac_mw, meanHMfrac_vw,
  maxtemp, mintemp,
  maxdens, mindens,
  maxent, minent,
  maxH2frac, minH2frac,
  maxHIIfrac, minHIIfrac,
  maxelecfrac, minelecfrac,
  maxHMfrac, minHMfrac,
  maxjeans, minjeans,
  totalmass, totalvolume,
    metalmff, metalvff, metalmass, metalvol,
  temp_od_2ddist_mw[NUMBINS][NUMBINS],
  temp_od_2ddist_vw[NUMBINS][NUMBINS],

  ent_od_2ddist_mw[NUMBINS][NUMBINS],
  ent_od_2ddist_vw[NUMBINS][NUMBINS],

  temp_ent_2ddist_mw[NUMBINS][NUMBINS],
  temp_ent_2ddist_vw[NUMBINS][NUMBINS],

  jeans_od_2ddist_mw[NUMBINS][NUMBINS],
  jeans_od_2ddist_vw[NUMBINS][NUMBINS],

  H2frac_od_2ddist_mw[NUMBINS][NUMBINS],
  H2frac_od_2ddist_vw[NUMBINS][NUMBINS],
  
  HIIfrac_od_2ddist_mw[NUMBINS][NUMBINS],
  HIIfrac_od_2ddist_vw[NUMBINS][NUMBINS],

  HMfrac_od_2ddist_mw[NUMBINS][NUMBINS],
  HMfrac_od_2ddist_vw[NUMBINS][NUMBINS],

  efrac_od_2ddist_mw[NUMBINS][NUMBINS],
  efrac_od_2ddist_vw[NUMBINS][NUMBINS],

  efrac_H2frac_2ddist_mw[NUMBINS][NUMBINS],
  efrac_H2frac_2ddist_vw[NUMBINS][NUMBINS],

  HMfrac_H2frac_2ddist_mw[NUMBINS][NUMBINS],
  HMfrac_H2frac_2ddist_vw[NUMBINS][NUMBINS],

  efrac_HMfrac_2ddist_mw[NUMBINS][NUMBINS],
  efrac_HMfrac_2ddist_vw[NUMBINS][NUMBINS],

  temp_1ddist_value[NUMBINS],
  temp_1ddist_mw[NUMBINS],
  temp_1ddist_vw[NUMBINS],
  temp_1ddist_mw_cum[NUMBINS],
  temp_1ddist_vw_cum[NUMBINS],

  od_1ddist_value[NUMBINS],
  od_1ddist_mw[NUMBINS],
  od_1ddist_vw[NUMBINS],
  od_1ddist_mw_cum[NUMBINS],
  od_1ddist_vw_cum[NUMBINS],

  ent_1ddist_value[NUMBINS],
  ent_1ddist_mw[NUMBINS],
  ent_1ddist_vw[NUMBINS],
  ent_1ddist_mw_cum[NUMBINS],
  ent_1ddist_vw_cum[NUMBINS],

  H2frac_1ddist_value[NUMBINS],
  H2frac_1ddist_mw[NUMBINS],
  H2frac_1ddist_vw[NUMBINS],
  H2frac_1ddist_mw_cum[NUMBINS],
  H2frac_1ddist_vw_cum[NUMBINS],

  HIIfrac_1ddist_value[NUMBINS],
  HIIfrac_1ddist_mw[NUMBINS],
  HIIfrac_1ddist_vw[NUMBINS],
  HIIfrac_1ddist_mw_cum[NUMBINS],
  HIIfrac_1ddist_vw_cum[NUMBINS],

  efrac_1ddist_value[NUMBINS],
  efrac_1ddist_mw[NUMBINS],
  efrac_1ddist_vw[NUMBINS],
  efrac_1ddist_mw_cum[NUMBINS],
  efrac_1ddist_vw_cum[NUMBINS],

  HMfrac_1ddist_value[NUMBINS],
  HMfrac_1ddist_mw[NUMBINS],
  HMfrac_1ddist_vw[NUMBINS],
  HMfrac_1ddist_mw_cum[NUMBINS],
  HMfrac_1ddist_vw_cum[NUMBINS],

  jeans_1ddist_value[NUMBINS],
  jeans_1ddist_mw[NUMBINS],
  jeans_1ddist_vw[NUMBINS],
  jeans_1ddist_mw_cum[NUMBINS],
  jeans_1ddist_vw_cum[NUMBINS],

  dbin, tbin, sbin, H2bin, ebin, HMbin,jbin, HIIbin;

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


float *densitybuff,*temperaturebuff,*H2Ibuff,*Hminusbuff,*electronbuff,*HIIbuff, *metalbuff;


// simulation parameters
double boxsize;
int rootgridsize, multispecies;
double mean_density,omegamatter,omegalambda,massconv,hubble,redshift,currenttime;

double densconstant, jeansmass_prefactor;


int numberofhalos,dark_matter_exists, total_number_grids;

int ReadParameterFile(char *filename);
int GetGridInfo(int numberofgrids,char *hierfilename);
int NumberOfGrids(char *hierfilename);
int GetCellInformationFromGrid(int gridnum,int total_number_grids);
int FlagGridCells(int gridnum,int total_number_grids);
int GetDensityInfo(int gridnum);

void DeclareGridArrays(void);
void CleanGridArrays(void);

void OutputAllInformation(char *infilename);

void SetAllArraysAndValues(void);

void CleanUpArraysAndValues(void);

