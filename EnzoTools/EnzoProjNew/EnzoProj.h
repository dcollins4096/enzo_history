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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>

#define hssize_t hsize_t // hdf5 1.6.4 and greater

#include "extern_hdf5.h"  // hdf 5 prototypes

#define NEW_INPUT_STYLE  // this is either NEW_INPUT_STYLE or OLD_INPUT_STYLE
                         // for the -b and -f input pieces.

#define MAX_LINE_LENGTH 256
#define HDF5_I4 H5T_NATIVE_INT
//#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I4 H5T_IEEE_F32BE  

#define DEBUG 1  // 1 on, 0 off
#define VERBOSEDEBUG 0 // 1 on, 0 off

// return calls
#define SUCCESS 1
#define FAILURE 0

#define SWAP_ARRAYS 1 //  // if 1, swaps column-major output arrays to row major.


// physical constants
#define MPC_CM 3.0824e+24  // megaparsec to centimeter conversion
#define MSOLAR_G 1.989e+33 // solar mass to grams conversion
#define RHOCRIT 1.8788e-29 // critical density in g/cm^3
 
#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))

// in file EP_ParseInputs.C
int ParseArgs(int argc,char *argv[]);
int ReadParameterFile(char *filename);
int CheckBoundaryValues(void);
void HelpMe(void);

// in file EP_Misc.C
int CreateProjectionArrays(void);
int DeleteProjectionArrays(void);
int NormalizeProjectionArrays(void);
int WriteProjectionArrays(void);
int CheckForMetals(void);
int SwapArray(int xcells, int ycells, float *arraytoswap);

// in file EP_GridStuff.C
int GetGridInfo(int numberofgrids,char *hierfilename);
int NumberOfGrids(char *hierfilename);

// in file EP_MakeProj.C
int AddGridToProjection(int gridnum,int total_number_grids);

// in file EP_PartProj.C
int AddParticlesToProjection(int gridnum);

// in file EP_MEKAL.C
int ReadMEKALTable(char *mekalfilename);
float ComputeSpectralEmissivity(float density, float temperature);

