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

// global variables originally declared in another file
extern int projaxis,projlevel,  maxlevel,rootgridsize,
  outputfile_newname;
extern double xstart,ystart,zstart,xend,yend,zend;
extern char **gridfilenames, *outputfilename;
extern int xprojnumcells,yprojnumcells,zprojnumcells,metal_flag,z1_flag,z2_flag;
extern float *projbuff_xrayem, *tempbuff, *densitybuff;

// global vars declared here for the first time




/*------------------------ ReadMEKALTable() --------------------
 *
 * reads in a MEKAL table which is used in ComputeSpectralEmissivity
 * to calculate spectral emissivities (go figure).
 *
 * -------------------------------------------------------------*/

int ReadMEKALTable(char *mekal_filename){
  if(DEBUG) fprintf(stderr,"in ReadMEKALTable:  %s\n",mekal_filename);
  int returnval=SUCCESS;
  FILE *fptr;

  fptr = fopen(mekal_filename,"r");





  fclose(fptr);

  if(DEBUG) fprintf(stderr,"leaving ReadMEKALTable\n");
  return returnval;
}


/*---------------- ComputeSpectralEmissivity() -----------------
 *
 * Calculates spectral emissivities using the tables obtained in
 * ReadMEKALTable()
 *
 * -------------------------------------------------------------*/
float ComputeSpectralEmissivity(float density, float temperature){
  if(VERBOSEDEBUG) fprintf(stderr,"in ComputeSpectramEmissivity, %e %e\n",density, temperature);
  int returnval=SUCCESS;
  float emissivity=0.0;




  if(VERBOSEDEBUG) fprintf(stderr,"leaving ComputeSpectramEmissivity\n");
  return emissivity;
}
