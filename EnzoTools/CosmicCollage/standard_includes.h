#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#include <time.h>
//#include "extern_hdf5.h"  // hdf 5 prototypes

#define MAX_LINE_LENGTH 256
#define HDF5_I4 H5T_NATIVE_INT
//#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I4 H5T_IEEE_F32BE  

// return calls
#define SUCCESS 0
#define FAILURE 1

// physical constants
#define MPC_CM 3.0824e+24  // megaparsec to centimeter conversion
#define MSOLAR_G 1.989e+33 // solar mass to grams conversion
#define RHOCRIT_C0S 2.7754e11  // comoving critical density in Msolar/Mpc^3 - needs h^2
#define RHOCRIT_CGS 1.8788e-29 // comoving critical density in g/cm^3 - needs h^2
#define MPROTON 1.6726e-24 // proton mass in cgs
#define MELECTRON 9.1094e-28 // electron mass in cgs
#define G_CGS 6.6726e-8  // gravitational constant in cgs
#define C_CGS 2.9979e10  // speed of light in cgs

#define PI 3.1415926

#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))

#define MAX_PROJECTION_FILES 100

#define HALO_YDEC_BINS 10
