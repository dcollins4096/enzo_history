/********************************************************************************
 *                       multi_CloudyCooling_time.C                             *
 *                      Britton Smith - April, 2006                             *
 *                                                                              *
 * Compute cooling time using method prescribed by mulit_CloudyCooling.C.       *
 *                                                                              *
 *        March, 2007: changed calculation of H number density for cloudy       *
 *                     cooling.                                                 *
 *                                                                              *
 *      October, 2006: added feature to prevent gas from cooling below CMB      *
 *                     temperature by calculating cooling as C(T)-C(T_CMB).     *
 *                                                                              *
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

int multi_CloudyCooling_time(float *density,float *totalenergy,float *gasenergy,
			     float *velocity1,float *velocity2,float *velocity3,
			     float *De,float *HI,float *HII,
			     float *HeI,float *HeII,float *HeIII,
			     float *HM,float *H2I,float *H2II,
			     float *metallicity,
			     float *cooling_time,
			     int *GridDimension,int GridRank,float dtFixed,
			     float afloat,float TemperatureUnits,float LengthUnits,
			     float aUnits,float DensityUnits,float TimeUnits)
{

  /* Compute size (in floats) of the current grid. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  // Declare some variables.

  int i;
  float mu,energy;
  float gamma2,nH2,nOther,x;
  float logNumberDensity,cellTemperature;
  double cooling,heating,edot,coolingAtTCMB,edotMetals;

  // Variables for interpolation over enzo cooling data

  float logtem0 = log10(CoolData.TemperatureStart);
  float logtem9 = log10(CoolData.TemperatureEnd);
  float dlogtem = (logtem9 - logtem0)/(CoolData.NumberOfTemperatureBins - 1);
  float logCellTemperature,logTIndex;
  int enzoCoolIndex;

  // H/He atomic cooling variables

  double ceHI,ceHeI,ceHeII,ciHI,ciHeI,ciHeIS,ciHeII,
    reHII,reHeII1,reHeII2,reHeIII,brem;

  // H2 cooling variables

  float kH_H2,kH2_H2,vibrationH,vibrationL,rotationH,rotationL,Qn;

  // H2 cooling variables for Galli & Palla 1998

  float GPLow,GPHigh,GPHigh1;

  // Compton cooling variables

  float compton1,compton2;

  /* Get conversion units (ripped from InitializeEquilibriumCoolData.C) */

  double dom    = DensityUnits*pow(afloat,3)/MH;
  double tbase1 = TimeUnits;
  double xbase1 = LengthUnits/(afloat*aUnits);
  double dbase1 = DensityUnits * POW(afloat*aUnits, 3);
  double CoolUnit = (POW(aUnits,5) * POW(xbase1,2) * POW(MH,2)) /
                    (POW(tbase1,3) * dbase1);
  double CurrentRedshift = 1.0/(afloat*aUnits) - 1.0;

  // Set compton cooling

  if (ComovingCoordinates) {
    compton1 = CoolData.comp * pow((1 + CurrentRedshift),4);
    compton2 = 2.73 * (1 + CurrentRedshift);
  }
  else {
    compton1 = tiny_number;
    compton2 = tiny_number;
  }

  // Loop over all grid cells.

  for (i = 0;i < size;i++) {

    edot = 0;

    // convert density from comoving to proper

    density[i] /= pow(afloat,3);
    De[i]      /= pow(afloat,3);
    HI[i]      /= pow(afloat,3);
    HII[i]     /= pow(afloat,3);
    HeI[i]     /= pow(afloat,3);
    HeII[i]    /= pow(afloat,3);
    HeIII[i]   /= pow(afloat,3);
    HM[i]      /= pow(afloat,3);
    H2I[i]     /= pow(afloat,3);
    H2II[i]    /= pow(afloat,3);

    // calculate mu
    mu = density[i] / (De[i] + HI[i] + HII[i] + 
		       (HeI[i] + HeII[i] + HeIII[i])/4.0 +
		       HM[i] + (H2I[i] + H2II[i])/2.0);

    logNumberDensity = log10(density[i]*H_MASS_FRACTION*DensityUnits/MH);

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

    cellTemperature = (Gamma - 1.0) * energy * TemperatureUnits * mu;

    // Correct temperature for gamma from H2 (ripped from cool1d_multi.src)

    nH2 = (H2I[i] + H2II[i])/2.0;
    nOther = De[i] + HI[i] + HII[i] + (HeI[i] + HeII[i] + HeIII[i])/4.0; // why no H-?

    if (nH2/nOther > 1e-3) {
      x = 6100 / cellTemperature;
      if (x > 10) {
	gamma2 = 2.5;
      }
      else {
	gamma2 = 0.5*(5.0 + 2.0*pow(x,2)*exp(x)/pow(exp(x)-1,2));
      }
    }
    else {
      gamma2 = 2.5;
    }
    gamma2 = 1.0 + (nH2 + nOther)/(nH2*gamma2 + nOther/(Gamma-1.0));
    cellTemperature *= (gamma2 - 1.0)/(Gamma - 1.0);
    
    // Get Index for enzo interpolation

    logCellTemperature = log10(cellTemperature);

    enzoCoolIndex = int((logCellTemperature - logtem0)/dlogtem);
    enzoCoolIndex = max(enzoCoolIndex,0);
    enzoCoolIndex = min(enzoCoolIndex,CoolData.NumberOfTemperatureBins - 2);
    logTIndex = enzoCoolIndex*dlogtem + logtem0;

    // If cloudy is supplying only metals, get all H/He cooling

    if (CloudyCooling == 2) {
      ceHI = CoolData.ceHI[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.ceHI[enzoCoolIndex+1] - CoolData.ceHI[enzoCoolIndex])/dlogtem;
      ceHeI = CoolData.ceHeI[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.ceHeI[enzoCoolIndex+1] - CoolData.ceHeI[enzoCoolIndex])/dlogtem;
      ceHeII = CoolData.ceHeII[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.ceHeII[enzoCoolIndex+1] - CoolData.ceHeII[enzoCoolIndex])/dlogtem;
      ciHI = CoolData.ciHI[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.ciHI[enzoCoolIndex+1] - CoolData.ciHI[enzoCoolIndex])/dlogtem;
      ciHeI = CoolData.ciHeI[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.ciHeI[enzoCoolIndex+1] - CoolData.ciHeI[enzoCoolIndex])/dlogtem;
      ciHeIS = CoolData.ciHeIS[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.ciHeIS[enzoCoolIndex+1] - CoolData.ciHeIS[enzoCoolIndex])/dlogtem;
      ciHeII = CoolData.ciHeII[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.ciHeII[enzoCoolIndex+1] - CoolData.ciHeII[enzoCoolIndex])/dlogtem;
      reHII = CoolData.reHII[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.reHII[enzoCoolIndex+1] - CoolData.reHII[enzoCoolIndex])/dlogtem;
      reHeII1= CoolData.reHeII1[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.reHeII1[enzoCoolIndex+1]- CoolData.reHeII1[enzoCoolIndex])/dlogtem;
      reHeII2= CoolData.reHeII2[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.reHeII2[enzoCoolIndex+1]- CoolData.reHeII2[enzoCoolIndex])/dlogtem;
      reHeIII= CoolData.reHeIII[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.reHeIII[enzoCoolIndex+1]- CoolData.reHeIII[enzoCoolIndex])/dlogtem;
      brem = CoolData.brem[enzoCoolIndex] + (logCellTemperature - logTIndex)
	*(CoolData.brem[enzoCoolIndex+1] - CoolData.brem[enzoCoolIndex])/dlogtem;

      edot += (

	       // Collisional excitations

	       - ceHI  *HI  [i]*De[i]             // ce of HI
	       - ceHeI *HeII[i]*De[i]*De[i]*dom/4.0  // ce of HeI
	       - ceHeII*HeII[i]*De[i]/4.0         // ce of HeII

	       // Collisional ionizations

	       - ciHI  *HI  [i]*De[i]             // ci of HI
	       - ciHeI *HeI [i]*De[i]/4.0         // ci of HeI
	       - ciHeII*HeII[i]*De[i]/4.0         // ci of HeII
	       - ciHeIS*HeII[i]*De[i]*De[i]*dom/4.0  // ci of HeIS

	       // Recombinations

	       - reHII  *HII  [i]*De[i]          // re of HII
	       - reHeII1*HeII [i]*De[i]/4.0      // re of HeII
	       - reHeII2*HeII [i]*De[i]/4.0      // re of HeII
	       - reHeIII*HeIII[i]*De[i]/4.0      // re of HeIII

	       // X-ray compton heating

	       - CoolData.comp_xray * (cellTemperature-CoolData.temp_xray)*De[i]/dom

	       // Bremsstrahlung

	       - brem*(HII[i]+HeII[i]/4.0+HeIII[i])*De[i]

	       )/density[i];
    }

#define USE_GALLI_PALLA1999
#ifdef USE_GALLI_PALLA1999

	// Get H2 cooling from Galli & Palla 1998

	GPLow = (logCellTemperature - logTIndex)*
	  (CoolData.GP99LowDensityLimit[enzoCoolIndex+1]-
	   CoolData.GP99LowDensityLimit[enzoCoolIndex])/dlogtem
	  + CoolData.GP99LowDensityLimit[enzoCoolIndex];

	GPHigh = (logCellTemperature - logTIndex)*
	  (CoolData.GP99HighDensityLimit[enzoCoolIndex+1]-
	   CoolData.GP99HighDensityLimit[enzoCoolIndex])/dlogtem
	  + CoolData.GP99HighDensityLimit[enzoCoolIndex];

	GPHigh1 = GPHigh/(HI[i]*dom);

	edot -= (H2I[i]/2.0/dom/density[i]) * GPHigh/(1 + GPHigh1/GPLow);

#else /* USE_GALLI_PALLA1999 */

    // Get H2 cooling from Lepp & Shull 1983

    kH_H2 = (logCellTemperature - logTIndex)*
      (CoolData.hyd01k[enzoCoolIndex+1]-CoolData.hyd01k[enzoCoolIndex])/dlogtem +
      CoolData.hyd01k[enzoCoolIndex];
    kH2_H2 = (logCellTemperature - logTIndex)*
      (CoolData.h2k01[enzoCoolIndex+1]-CoolData.h2k01[enzoCoolIndex])/dlogtem +
      CoolData.h2k01[enzoCoolIndex];
    vibrationH = (logCellTemperature - logTIndex)*
      (CoolData.vibh[enzoCoolIndex+1]-CoolData.vibh[enzoCoolIndex])/dlogtem +
      CoolData.vibh[enzoCoolIndex];
    rotationH = (logCellTemperature - logTIndex)*
      (CoolData.roth[enzoCoolIndex+1]-CoolData.roth[enzoCoolIndex])/dlogtem +
      CoolData.roth[enzoCoolIndex];
    rotationL = (logCellTemperature - logTIndex)*
      (CoolData.rotl[enzoCoolIndex+1]-CoolData.rotl[enzoCoolIndex])/dlogtem +
      CoolData.rotl[enzoCoolIndex];
    
    Qn = 1.2*pow((HI[i]*dom),0.77) + pow((H2I[i]*dom/2.0),0.77);
    rotationL *= Qn;
    vibrationL = (HI[i]*kH_H2 + H2I[i]/2.0*kH2_H2)*dom*8.18e-13;

    edot -= (H2I[i]/2.0/dom/density[i]) * (rotationH/(1 + rotationH/rotationL) + 
					   vibrationH/(1 + vibrationH/vibrationL));

#endif /* USE_GALLI_PALLA1999 */

    // Just interpolate over temperature.

    if (CloudyCoolingData.CloudyCoolingGridRank == 0) {

      // Get Cloudy cooling

      cooling = coolingGridInterpolate0D_float(cellTemperature,
					       CloudyCoolingData.coolingGridCooling);

      edotMetals = -(pow(10,cooling)) * density[i] * H_MASS_FRACTION;

      // If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

      if (CMBTemperatureFloor) {
	coolingAtTCMB = coolingGridInterpolate0D_float(compton2,
						       CloudyCoolingData.coolingGridCooling);
	if (cooling > coolingAtTCMB) {
	  edotMetals += (pow(10,coolingAtTCMB)) * density[i] * H_MASS_FRACTION;
	}
	else {
	  edotMetals = 0.0;
	}
      }

      // include heating only if requested

      if (CloudyCoolingData.IncludeCloudyHeating) {
	heating = coolingGridInterpolate0D_float(cellTemperature,
						 CloudyCoolingData.coolingGridHeating);

	edotMetals += ((pow(10,heating)) * density[i]) * H_MASS_FRACTION;
      }

      edot += edotMetals;

    }

    // Interpolate over density and temperature.

    else if (CloudyCoolingData.CloudyCoolingGridRank == 1) {

      // Get Cloudy cooling

      cooling = coolingGridInterpolate1D_float(logNumberDensity,cellTemperature,
					       CloudyCoolingData.coolingGridCooling);

      edotMetals = -(pow(10,cooling)) * density[i] * H_MASS_FRACTION;

      // If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

      if (CMBTemperatureFloor) {
	coolingAtTCMB = coolingGridInterpolate1D_float(logNumberDensity,compton2,
						       CloudyCoolingData.coolingGridCooling);
	if (cooling > coolingAtTCMB) {
	  edotMetals += (pow(10,coolingAtTCMB)) * density[i] * H_MASS_FRACTION;
	}
	else {
	  edotMetals = 0.0;
	}
      }

      // include heating only if requested

      if (CloudyCoolingData.IncludeCloudyHeating) {
	heating = coolingGridInterpolate1D_float(logNumberDensity,cellTemperature,
						 CloudyCoolingData.coolingGridHeating);

	edotMetals += ((pow(10,heating)) * density[i]) * H_MASS_FRACTION;
      }

      edot += edotMetals;

    }

    // Interpolate over density, metallicity, and temperature.

    else if (CloudyCoolingData.CloudyCoolingGridRank == 2) {

      // Get Cloudy cooling

      cooling = coolingGridInterpolate2D_float(logNumberDensity,metallicity[i],cellTemperature,
					       CloudyCoolingData.coolingGridCooling);

      edotMetals = -(pow(10,cooling)) * density[i] * H_MASS_FRACTION;
 
      // If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

      if (CMBTemperatureFloor) {
	coolingAtTCMB = coolingGridInterpolate2D_float(logNumberDensity,metallicity[i],compton2,
						       CloudyCoolingData.coolingGridCooling);
	if (cooling > coolingAtTCMB) {
	  edotMetals += (pow(10,coolingAtTCMB)) * density[i] * H_MASS_FRACTION;
	}
	else {
	  edotMetals = 0.0;
	}
      }

      // include heating only if requested

      if (CloudyCoolingData.IncludeCloudyHeating) {
	heating = coolingGridInterpolate2D_float(logNumberDensity,metallicity[i],cellTemperature,
						 CloudyCoolingData.coolingGridHeating);

	edotMetals += ((pow(10,heating)) * density[i]) * H_MASS_FRACTION;
      }

      edot += edotMetals;

    }

    else {

      fprintf(stderr,"CloudyCoolingData.CloudyCoolingGridRank must be 0, 1, or 2.\n");
      return FAIL;

    }
    
    // Add Compton cooling
    
    edot -= compton1*(cellTemperature-compton2)*De[i]/dom/density[i];
      
    // compute cooling time

    cooling_time[i] = -energy/edot;

    // convert density back to comoving

    density[i] *= pow(afloat,3);
    De[i]      *= pow(afloat,3);
    HI[i]      *= pow(afloat,3);
    HII[i]     *= pow(afloat,3);
    HeI[i]     *= pow(afloat,3);
    HeII[i]    *= pow(afloat,3);
    HeIII[i]   *= pow(afloat,3);
    HM[i]      *= pow(afloat,3);
    H2I[i]     *= pow(afloat,3);
    H2II[i]    *= pow(afloat,3);

  } // for (i = 0;i < size;i++)

  return SUCCESS;

}
