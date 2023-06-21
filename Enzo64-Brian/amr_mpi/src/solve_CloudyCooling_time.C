/********************************************************************************
 *                        solve_CloudyCooling_time.C                            *
 *                       Britton Smith - April, 2006                            *
 *                                                                              *
 * Compute cooling time using cooling method prescribed by                      *
 * solve_CloudyCooling.C.                                                       *
 *                                                                              *
 *        March, 2007: changed calculation of H number density for cloudy       *
 *                     cooling.                                                 *
 *                                                                              *
 *      October, 2006: added feature to prevent gas from cooling below CMB      *
 *                     temperature by calculating cooling as C(T)-C(T_CMB).     *
 * Updated June, 2006: added metallicity field to handle higher dimension of    *
 *                     interpolation.                                           *
 *******************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#define MH 1.67e-24
#define DEFAULT_MU 1.22
#define H_MASS_FRACTION 0.76

float coolingGridInterpolate0D_float(float temperature,float *dataField);
float coolingGridInterpolate1D_float(float parameter1,float temperature,float *dataField);
float coolingGridInterpolate2D_float(float parameter1,float parameter2,float temperature,
				       float *dataField);
float coolingGridInterpolate1D_float(float parameter1,float temperature,float *dataField);

int solve_CloudyCooling_time(float *density,float *totalenergy,float *gasenergy,
			     float *velocity1,float *velocity2,float *velocity3,
			     float *metallicity,
			     float *cooling_time,
			     int *GridDimension,int GridRank,float dtFixed,
			     float afloat,float TemperatureUnits,float aUnits,
			     float DensityUnits)
{

  /* Compute size (in floats) of the current grid. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  // Loop over all cells, solve cooling, update energies

  int i,iter;
  float mu,mu_old,mu_tol,energy;
  float properDensity,numberDensity,cellTemperature;
  double cooling,heating,edot,coolingAtTCMB;

  // CMB Temperature
  float compton2;

  double CurrentRedshift = 1.0/(afloat*aUnits) - 1.0;

  // Set compton cooling

  if (ComovingCoordinates) {
    compton2 = 2.73 * (1 + CurrentRedshift);
  }

  mu_tol = .001;
  mu = 1.0;

  for (i = 0;i < size;i++) {

    // convert density from comoving to proper

    properDensity = H_MASS_FRACTION * density[i] / pow(afloat,3);

    // calculate cell temperature

    if (DualEnergyFormalism) {
      energy = gasenergy[i];
    }
    else {
      energy = totalenergy[i] - 0.5*velocity1[i]*velocity1[i];
      if (GridRank > 1)
	energy -= 0.5*velocity2[i]*velocity2[i];
      if (GridRank > 2)
	energy -= 0.5*velocity3[i]*velocity3[i];
    }

    cellTemperature = (Gamma - 1.0) * energy * TemperatureUnits * DEFAULT_MU;

    /* In case we ever want to solve for mmw accurately, uncomment this block. *
    // converge on mean molecular weight and number density

    mu_old = 0;
    iter = 0;
    while ((fabs(mu-mu_old)/mu > mu_tol) && (iter < 10)) {
      mu_old = mu;
      numberDensity = log10(properDensity*DensityUnits/mu/MH);
      mu = coolingGridInterpolate1D_float(numberDensity,cellTemperature,
    				    CloudyCoolingData.coolingGridMeanMolecularWeight);
      iter++;
    }
    numberDensity += log10(mu_old/mu);
    ****************************************************************************/

    // just take mu = 1.22 for now
    mu = DEFAULT_MU;
    numberDensity = log10(properDensity*DensityUnits/mu/MH);

    if (CloudyCoolingData.CloudyCoolingGridRank == 0) {

      // Get Cloudy cooling

      cooling = coolingGridInterpolate0D_float(cellTemperature,
						CloudyCoolingData.coolingGridCooling);

      // include heating only if requested

      if (CloudyCoolingData.IncludeCloudyHeating) {
	heating = coolingGridInterpolate0D_float(cellTemperature,
						  CloudyCoolingData.coolingGridHeating);

	edot = -(pow(10,cooling) - pow(10,heating)) * properDensity;
      }

      // otherwise, just cooling

      else {
	edot = -(pow(10,cooling)) * properDensity;
      }

      // If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

      if (CMBTemperatureFloor) {
	coolingAtTCMB = coolingGridInterpolate0D_float(compton2,
						       CloudyCoolingData.coolingGridCooling);
	edot += (pow(10,coolingAtTCMB)) * properDensity;
      }

    }

    // Interpolate over density and temperature.

    else if (CloudyCoolingData.CloudyCoolingGridRank == 1) {

      // Get Cloudy cooling

      cooling = coolingGridInterpolate1D_float(numberDensity,cellTemperature,
						CloudyCoolingData.coolingGridCooling);

      // include heating only if requested

      if (CloudyCoolingData.IncludeCloudyHeating) {
	heating = coolingGridInterpolate1D_float(numberDensity,cellTemperature,
						  CloudyCoolingData.coolingGridHeating);

	edot = -(pow(10,cooling) - pow(10,heating)) * properDensity;
      }

      // otherwise, just cooling

      else {
	edot = -(pow(10,cooling)) * properDensity;
      }

      // If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

      if (CMBTemperatureFloor) {
	coolingAtTCMB = coolingGridInterpolate1D_float(numberDensity,compton2,
						       CloudyCoolingData.coolingGridCooling);
	edot += (pow(10,coolingAtTCMB)) * properDensity;
      }

    }

    // Interpolate over density, metallicity, and temperature.

    else if (CloudyCoolingData.CloudyCoolingGridRank == 2) {

      // Get Cloudy cooling

      cooling = coolingGridInterpolate2D_float(numberDensity,metallicity[i],cellTemperature,
						CloudyCoolingData.coolingGridCooling);

      // include heating only if requested

      if (CloudyCoolingData.IncludeCloudyHeating) {
	heating = coolingGridInterpolate2D_float(numberDensity,metallicity[i],cellTemperature,
						  CloudyCoolingData.coolingGridHeating);

	edot = -(pow(10,cooling) - pow(10,heating)) * properDensity;
      }

      // otherwise, just cooling

      else {
	edot = -(pow(10,cooling)) * properDensity;
      }

      // If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

      if (CMBTemperatureFloor) {
	coolingAtTCMB = coolingGridInterpolate2D_float(numberDensity,metallicity[i],compton2,
						       CloudyCoolingData.coolingGridCooling);
	edot += (pow(10,coolingAtTCMB)) * properDensity;
      }

    }

    else {

      fprintf(stderr,"CloudyCoolingData.CloudyCoolingGridRank must be 0, 1, or 2.\n");
      return FAIL;

    }

    // compute cooling time

    cooling_time[i] = -energy/edot;

  } // for (i = 0;i < size;i++)

  return SUCCESS;

}
