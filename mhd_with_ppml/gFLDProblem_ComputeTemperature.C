/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Implicit Problem Class temperature 
/  calculation routine (nearly identical to 
/  Grid_ComputeTemperatureField, but allows input arguments based on u 
/  fields, as opposed to just extracting all information out of the grid.
/
/  written by: Daniel Reynolds
/  date:       October, 2006
/  modified1:  
/
/  PURPOSE: Takes in EnzoVector and returns temperature field at 
/           desired time.  This routine will be called repeatedly, so 
/           it should NOT result in any net allocation of memory.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"



/* default constants */
#define DEFAULT_MU 0.6       // mean molecular mass
#define MIN_TEMP 1.0     // minimum temperature [K]


/* function prototypes */
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);




int gFLDProblem::ComputeTemperature(float *TempArr, float time, 
				    FLOAT a, EnzoVector *u) 
{

//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::ComputeTemperature routine\n");

  // get local mesh description
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) {
    fprintf(stderr,"Temperature error: x0 vector dims do not match\n");
    return FAIL;
  }
  if (usz[1] != LocDims[1]) {
    fprintf(stderr,"Temperature error: x1 vector dims do not match\n");
    return FAIL;
  }
  if (usz[2] != LocDims[2]) {
    fprintf(stderr,"Temperature error: x2 vector dims do not match\n");
    return FAIL;
  }
  if (usz[3] != (2+Nchem)) {
    fprintf(stderr,"Temperature error: nspecies dims do not match\n");
    return FAIL;
  }
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) {
    fprintf(stderr,"Temperature error: x0 vector sizes do not match\n");
    return FAIL;
  }
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) {
    fprintf(stderr,"Temperature error: x1 vector sizes do not match\n");
    return FAIL;
  }
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) {
    fprintf(stderr,"Temperature error: x2 vector sizes do not match\n");
    return FAIL;
  }

  // extract fluid energy, radiation energy and chemistry arrays
  float *ec = u->GetData(1);
  float *ni[Nchem];
  for (int i=0; i<Nchem; i++)  ni[i] = u->GetData(2+i);


  // Compute the size of the fields
  int size=1;
  int i, dim;
  for (dim = 0; dim<3; dim++)  size *= ArrDims[dim];


  // Find the temperature units if we are using comoving coordinates
  float TempUnits, RhoUnits, LengthUnits, VelUnits, TimeUnits;
  TempUnits = RhoUnits = LengthUnits = VelUnits = TimeUnits = 1.0;
  if (ComovingCoordinates)
    if (CosmologyGetUnits(&RhoUnits, &LengthUnits, &TempUnits,
			  &TimeUnits, &VelUnits, time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }
 
  ////////////////////////////
  // Compute the pressure first, storing in Temperature array
  if (DualEnergyFormalism) {
    // compute pressure using ideal gas law
    // (convert density from comoving to proper, scale to CGS)
    for (i=0; i<size; i++) {
      TempArr[i] = (Gamma - 1.0) * rho[i]/a/a/a * (eh[i] + ec[i]);
      // correct for near-zero temperature
      TempArr[i] = max(TempArr[i], tiny_number);
    }

    // Grid_ComputePressure here corrects the pressure due to
    // H2 fields, which we do not consider here, so move along
  }
  else {
    // compute pressure using ideal gas law
    // (convert density from comoving to proper)
    for (i=0; i<size; i++) {
      TempArr[i] = (Gamma - 1.0) * rho[i]/a/a/a * 
	((eh[i]+ec[i]) - 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]));
      // correct for near-zero temperature 
      TempArr[i] = max(TempArr[i], tiny_number);
    }

    // Grid_ComputePressure here corrects the pressure due to
    // H2 fields, which we do not consider here, so move along
  }
 

  //////////////////////
  // NOTE: For some reason the computations of chemical quantities 
  //       (nHII, nHI, etc.) give terribly incorrect values for the
  //       temperature.  For now, stick with the standard approach.
  //////////////////////

//   float num_density, nHI, nHII, nHeI, nHeII, nHeIII;
//   float mp = 1.67262171e-24;  // mass of a proton [g]
  // for no chemistry, compute temperature as pressure/density
  // (convert density from comoving to proper)
//   if (Nchem == 0) {
    for (i=0; i<size; i++)
      TempArr[i] = max((TempUnits*TempArr[i]*DEFAULT_MU
			/max(rho[i]/a/a/a,tiny_number)),MIN_TEMP);
//   }

//   // otherwise compute temperature with mu calculated directly
//   else {

//     // Hydrogen only
//     if (Nchem == 1) {
//       for (i=0; i<size; i++) {

// 	// compute the overall number density
// 	// (convert density from comoving to proper)
// 	nHI = ni[0][i];
// 	nHII = rho[i]/mp - nHI;
// 	num_density = (nHI + nHII + ne[i])/a/a/a;

// 	// compute Temperature
// 	TempArr[i] *= TempUnits/max(num_density, tiny_number);
// 	TempArr[i] = max(TempArr[i], MIN_TEMP);

//       }  // end: loop over i
//     }  // end: Nchem == 1

//     // Hydrogen and Helium
//     else {
//       for (i=0; i<size; i++) {

// 	// compute the overall number density
// 	// (convert density from comoving to proper)
// 	nHI = ni[0][i];
// 	nHII = 0.75*rho[i]/mp - nHI;
// 	nHeI = ni[1][i];
// 	nHeII = ni[2][i];
// 	nHeIII = 0.25*rho[i]/mp - nHeI - nHeII;
// 	num_density = (0.25*(nHeI+nHeII+nHeIII) + nHI + nHII + ne[i])/a/a/a;

// 	// compute Temperature
// 	TempArr[i] *= TempUnits/max(num_density, tiny_number);
// 	TempArr[i] = max(TempArr[i], MIN_TEMP);

//       }  // end: loop over i
//     }  // end: Hydrogen and Helium
//   }  // end: Nchem != 0


  // return success
  return SUCCESS;
}
#endif
