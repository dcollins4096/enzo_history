#ifndef __typedefs_h_
#define __typedefs_h_
/***********************************************************************
/
/  MISCELANEOUS TYPEDEFS AND ENUMERATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include "CoolData.h"
#include "RateData.h"
#include "RadiationFieldData.h"

/* These are the different types of baryon fields. */

#ifdef SMALL_INTS
typedef int field_type;
typedef int boundary_type;
typedef int gravity_boundary_type;
typedef int interpolation_type;
typedef int hydro_method;
#ifdef MHDF
typedef int MHD_divbmethod;
typedef int MHD_Centering;
#endif
#endif

#ifdef LARGE_INTS
typedef long_int field_type;
typedef long_int boundary_type;
typedef long_int gravity_boundary_type;
typedef long_int interpolation_type;
typedef long_int hydro_method;
#ifdef MHDF
typedef long_int MHD_divbmethod;
typedef long_int MHD_Centering;
#endif
#endif

const field_type 
  Density         = 0,
  TotalEnergy     = 1,
  InternalEnergy  = 2,
  Pressure        = 3,
  Velocity1       = 4,
  Velocity2       = 5,
  Velocity3       = 6,
  ElectronDensity = 7,
  HIDensity       = 8,
  HIIDensity      = 9,
  HeIDensity      = 10,
  HeIIDensity     = 11,
  HeIIIDensity    = 12,
  HMDensity       = 13,
  H2IDensity      = 14,
  H2IIDensity     = 15,
  DIDensity       = 16,
  DIIDensity      = 17,
  HDIDensity      = 18,
  Metallicity     = 19,
  ExtraType0      = 20,
  ExtraType1      = 21,
  GravPotential   = 22,
  Acceleration0   = 23,
  Acceleration1   = 24,
  Acceleration2   = 25,
  RadiationFreq0  = 26,
  RadiationFreq1  = 27,
  RadiationFreq2  = 28,
  RadiationFreq3  = 29,
  RadiationFreq4  = 30,
  RadiationFreq5  = 31,
#ifdef PPML
  Face_X_L_D       = 32,
  Face_X_L_TE      = 33,
  Face_X_L_VX      = 34,
  Face_X_L_VY      = 35,
  Face_X_L_VZ      = 36,
  Face_X_L_BX      = 37,
  Face_X_L_BY      = 38,
  Face_X_L_BZ      = 39,
  
  Face_X_R_D       = 40,
  Face_X_R_TE      = 41,
  Face_X_R_VX      = 42,
  Face_X_R_VY      = 43,
  Face_X_R_VZ      = 44,
  Face_X_R_BX      = 45,
  Face_X_R_BY      = 46,
  Face_X_R_BZ      = 47,
  
  Face_Y_L_D       = 48,
  Face_Y_L_TE      = 49,
  Face_Y_L_VX      = 50,
  Face_Y_L_VY      = 51,
  Face_Y_L_VZ      = 52,
  Face_Y_L_BX      = 53,
  Face_Y_L_BY      = 54,
  Face_Y_L_BZ      = 55,
  
  Face_Y_R_D       = 56,
  Face_Y_R_TE      = 57,
  Face_Y_R_VX      = 58,
  Face_Y_R_VY      = 59,
  Face_Y_R_VZ      = 60,
  Face_Y_R_BX      = 61,
  Face_Y_R_BY      = 62,
  Face_Y_R_BZ      = 63,
  
  Face_Z_L_D       = 64,
  Face_Z_L_TE      = 65,
  Face_Z_L_VX      = 66,
  Face_Z_L_VY      = 67,
  Face_Z_L_VZ      = 68,
  Face_Z_L_BX      = 69,
  Face_Z_L_BY      = 70,
  Face_Z_L_BZ      = 71,
  
  Face_Z_R_D       = 72,
  Face_Z_R_TE      = 73,
  Face_Z_R_VX      = 74,
  Face_Z_R_VY      = 75,
  Face_Z_R_VZ      = 76,
  Face_Z_R_BX      = 77,
  Face_Z_R_BY      = 78,
  Face_Z_R_BZ      = 79,

  Magnetic1        = 80,
  Magnetic2        = 81,
  Magnetic3        = 82,
  
  FieldUndefined  = 83;
#else
FieldUndefined  = 32;
#endif //PPML


#ifdef MHDF
const MHD_divbmethod 
  NoDivB = 0,
  BalsaraSpicer = 1,
  Athena_LF = 2,
  Athena_Switch = 3,
  BalsaraToth = 4;

const MHD_Centering
  MHD_CenteringUndefined=0, 
  MHD_SplitWithEng = 1, 
  MHD_Split = 2, 
  MHD_Volumetric = 3, 
  MHD_dtOnly = 4, 
  MHD_none = 5, 
  MHD_Non3d = 6;

#endif //MHDF

/*
enum field_type {Density, TotalEnergy, InternalEnergy, Pressure,
		 Velocity1, Velocity2, Velocity3, 
		 ElectronDensity, HIDensity, HIIDensity,  HeIDensity, 
		 HeIIDensity, HeIIIDensity, HMDensity, H2IDensity, 
		 H2IIDensity, DIDensity, DIIDensity, HDIDensity,
                 Metallicity, ExtraType0, ExtraType1, GravPotential,
		 Acceleration0, Acceleration1,Acceleration2,
		 RadiationFreq0, RadiationFreq1, RadiationFreq2, 
		 RadiationFreq3, RadiationFreq4, RadiationFreq5,
		 FieldUndefined};
*/

#define FieldTypeIsDensity(A) (((A) >= TotalEnergy && (A) <= Velocity3) ? FALSE : TRUE)

/* These are the different types of fluid boundary conditions. */

const boundary_type
  reflecting        = 0,
  outflow           = 1,
  inflow            = 2,
  periodic          = 3,
  BoundaryUndefined = 4;

// enum boundary_type {reflecting, outflow, inflow, periodic, BoundaryUndefined};

/* These are the different types of gravity boundary conditions. */

const gravity_boundary_type
  TopGridPeriodic  = 0,
  TopGridIsolated  = 1,
  SubGridIsolated  = 2,
  GravityUndefined = 3;

// enum gravity_boundary_type {TopGridPeriodic, TopGridIsolated, 
// 				    SubGridIsolated, GravityUndefined};

/* Interpolation types. */

const interpolation_type
  ThirdOrderA            = 0,
  SecondOrderA           = 1,
  SecondOrderB           = 2,
  SecondOrderC           = 3,
  FirstOrderA            = 4,
  InterpolationUndefined = 5;


// enum interpolation_type {ThirdOrderA, SecondOrderA, SecondOrderB, SecondOrderC,
// 			 FirstOrderA, InterpolationUndefined};

/* Hydrodynamics methods. */

const hydro_method
  PPM_DirectEuler      = 0,
  PPM_LagrangeRemap    = 1,
  Zeus_Hydro           = 2,
  PPM_Local            = 3,
  MHD_Test             = 4,
#ifdef MHDF
  MHD_Athena           = 5,
  HydroMethodUndefined = 6;
#else
  HydroMethodUndefined = 5;
#endif
// enum hydro_method {PPM_DirectEuler, PPM_LagrangeRemap, Zeus_Hydro, PPM_Local};

/* Define a float/int union. */

union float_int {
  long_int ival;
  float fval;
  FLOAT FVAL;
};

#endif
