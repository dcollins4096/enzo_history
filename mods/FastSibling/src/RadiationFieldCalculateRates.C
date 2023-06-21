/***********************************************************************
/
/  SETS THE MULTI-SPECIES RATES BASED ON THE EXTERNAL RADIATION FIELD
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);


int RadiationFieldCalculateRates(FLOAT Time)
{

  /* Return if there is no radiation (rates should be all zero). */

  if (RadiationFieldType == 0)
    return SUCCESS;

  /* Set units. */

  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;    

  if (!ComovingCoordinates) {
    fprintf(stderr, "RadiationField only defined for cosmology.\n");
    return FAIL;
  }

  CosmologyComputeExpansionFactor(Time, &a, &dadt);
  CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		    &TimeUnits, &VelocityUnits, Time);

  aUnits = 1.0/(1.0 + InitialRedshift);
  float Redshift = 1.0/(a*aUnits) - 1;

  double xbase1 = LengthUnits/(a*aUnits);
  double dbase1 = DensityUnits*pow(a*aUnits, 3);
  double mh     = 1.67e-24;
  double CoolingUnits = (pow(aUnits, 5) * xbase1*xbase1 * mh*mh) /
                       (pow(double(TimeUnits), 3) * dbase1);

  /* ------------------------------------------------------------------ */
  /* First, calculate the ramp value, a number between 0 and 1 which
     is used as an external control to the radiation. 
     (Only used if RadiationFieldType = 1 to 4). */

  float Ramp = 0;

  if (Redshift < CoolData.RadiationRedshiftOn && 
      Redshift > CoolData.RadiationRedshiftOff) {

    if (Redshift > CoolData.RadiationRedshiftFullOn)
      Ramp = 0.5 - 0.5*tanh(15.0*(Redshift - 0.5*
	    (CoolData.RadiationRedshiftOn+CoolData.RadiationRedshiftFullOn)));
    else if (Redshift < CoolData.RadiationRedshiftDropOff)
      Ramp = (Redshift - CoolData.RadiationRedshiftOff + CoolData.f0to3*
	                 (CoolData.RadiationRedshiftDropOff - Redshift)) /
             (CoolData.RadiationRedshiftDropOff - 
	      CoolData.RadiationRedshiftOff);
    else
      Ramp = 1.0;

  }

  float exp_arg = -1.0 * pow(Redshift-2.3, 2);

  /* ------------------------------------------------------------------ */
  /* 1) For the Haardt and Madau (1996) quasar spectrum (alpha_q = 1.5) */

  if (RadiationFieldType == 1) {
    RateData.k24 = 6.7e-13 * pow(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 * TimeUnits * Ramp;
    RateData.k25 = 6.3e-15 * pow(1.0+Redshift, 0.51) * exp(exp_arg/2.35)
                 * TimeUnits * Ramp;
    RateData.k26 = 3.2e-13 * pow(1.0+Redshift, 0.50) * exp(exp_arg/2.00)
                 * TimeUnits * Ramp;
    CoolData.piHI   = 4.7e-24 * pow(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 / CoolingUnits * Ramp;
    CoolData.piHeI  = 8.2e-24 * pow(1.0+Redshift, 0.50) * exp(exp_arg/2.00)
                 / CoolingUnits * Ramp;
    CoolData.piHeII = 1.6e-25 * pow(1.0+Redshift, 0.51) * exp(exp_arg/2.35)
                 / CoolingUnits * Ramp;
  }

  /* ------------------------------------------------------------------ */
  /* 2) For the Haardt and Madau (1996) quasar spectrum (alpha_q = 1.8) */

  if (RadiationFieldType == 2) {
    RateData.k24 = 5.6e-13 * pow(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 * TimeUnits * Ramp;
    RateData.k25 = 3.2e-15 * pow(1.0+Redshift, 0.30) * exp(exp_arg/2.60)
                 * TimeUnits * Ramp;
    RateData.k26 = 4.8e-13 * pow(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 * TimeUnits * Ramp;
    CoolData.piHI   = 3.9e-24 * pow(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 / CoolingUnits * Ramp;
    CoolData.piHeI  = 6.4e-24 * pow(1.0+Redshift, 0.43) * exp(exp_arg/2.10)
                 / CoolingUnits * Ramp;
    CoolData.piHeII = 8.7e-26 * pow(1.0+Redshift, 0.30) * exp(exp_arg/2.70)
                 / CoolingUnits * Ramp;
  }

  /* ------------------------------------------------------------------ */
  /* 3) This is a modified version of (1) but with the HeII heating rate
     multiplied by 1.8. */

  if (RadiationFieldType == 3) {
    RateData.k24 = 6.7e-13 * pow(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 * TimeUnits * Ramp;
    RateData.k25 = 6.3e-15 * pow(1.0+Redshift, 0.51) * exp(exp_arg/2.35)
                 * TimeUnits * Ramp;
    RateData.k26 = 3.2e-13 * pow(1.0+Redshift, 0.50) * exp(exp_arg/2.00)
                 * TimeUnits * Ramp;
    CoolData.piHI   = 4.7e-24 * pow(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 / CoolingUnits * Ramp;
    CoolData.piHeI  = 8.2e-24 * pow(1.0+Redshift, 0.50) * exp(exp_arg/2.00)
                 / CoolingUnits * Ramp;
    CoolData.piHeII = 1.6e-25 * pow(1.0+Redshift, 0.51) * exp(exp_arg/2.35)
                 / CoolingUnits * Ramp;
    CoolData.piHeII *= 1.8;
  }

  /* ------------------------------------------------------------------ */
  /* 4) For the Haardt and Madau (1996) quasar spectrum (alpha_q = 1.5) 
        with X-ray Compton heating from Madau & Efstathiou (9902080). */

  if (RadiationFieldType == 4) {
    RateData.k24 = 6.7e-13 * pow(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 * TimeUnits * Ramp;
    RateData.k25 = 6.3e-15 * pow(1.0+Redshift, 0.51) * exp(exp_arg/2.35)
                 * TimeUnits * Ramp;
    RateData.k26 = 3.2e-13 * pow(1.0+Redshift, 0.50) * exp(exp_arg/2.00)
                 * TimeUnits * Ramp;
    CoolData.piHI   = 4.7e-24 * pow(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 / CoolingUnits * Ramp;
    CoolData.piHeI  = 8.2e-24 * pow(1.0+Redshift, 0.50) * exp(exp_arg/2.00)
                 / CoolingUnits * Ramp;
    CoolData.piHeII = 1.6e-25 * pow(1.0+Redshift, 0.51) * exp(exp_arg/2.35)
                 / CoolingUnits * Ramp;

    /* This is sigma_thompson * c * (effective <h \mu>/<m_e c^2>) *
       U_xray * 1eV.  U_xray is the energy density of XRB in , <h \mu> is the
       average photon energy in keV, corrected for relativistic effects. */

    float RedshiftXrayCutoff = 5.0;
    /*
    CoolData.comp_xray = 6.65e-25 * 3.0e10 * 
                        (31.8*pow(1.0+Redshift, 0.3333)/511.0) * 
                        (6.3e-5 * 1.6e-12) * 
                        pow(1.0 + Redshift, 4) * 
                        exp(-pow(Redshift/RedshiftXrayCutoff, 2)) / 
                        CoolingUnits; */
    CoolData.comp_xray = 6.65e-25 * 3.0e10 * 
                        (1.0/511.0e3) * 
                        (4.0 * 1.38e-16/1.6e-12) *
                        (6.3e-5 * 1.6e-12) * 
                        pow(1.0 + Redshift, 4) * 
                        exp(-pow(Redshift/RedshiftXrayCutoff, 2)) / 
                        CoolingUnits;

    /* The effective temperature (in K). */

    CoolData.temp_xray = 31.8e3*pow(1.0+Redshift, 0.3333)*1.6e-12/
                         (4.0*1.38e-16);
  }


  /* ------------------------------------------------------------------ */
  /* 5) Featureless power-law spectrum (see calc_rates.src). */

  if (RadiationFieldType == 5)
    ;

  /* ------------------------------------------------------------------ */
  /* 8) An absorbed (hard) quasar-like spectrum plus molecular H constant
     photo-dissociation  (the rates are calculated in calc_rates.src). */

  if (RadiationFieldType == 8) {

    /* Insert redshift-dependent term here */

    /* molecular hydrogen constant photo-dissociation 
       rate is 1.13e-8 * F_LW  (flux in Lyman-Werner bands) 
       Note: this is hard-coded to F_LW = 1e-21; flux normalization controls
       the hard-radiation component (i.e. > 13.6 eV)*/

    RateData.k31 = 1.13e8 * 1.0e-21 * TimeUnits;
  }

  /* ------------------------------------------------------------------ */
  /* 9) molecular hydrogen constant photo-dissociation only! 
     rate is 1.13e-8 * F_LW  (flux in Lyman-Werner bands) */

  if (RadiationFieldType == 9)
    RateData.k31 = 1.13e8 * CoolData.f3 * TimeUnits;

  /* ------------------------------------------------------------------ */
  /* 10 & 11) - internally-computed radiation field.  Most of the rates
     are calculated in RadiationFieldUpdate, but the comp_xray rate is
     calculated here, from the X-ray energy density and temperature. */

  if (RadiationFieldType >= 10 && RadiationFieldType <= 11) {
    CoolData.comp_xray = 8.0/3.0*0.665e6/(9.1e0*3.0e0) *
      RadiationData.ComptonXrayEnergyDensity*1.38e2*
      (double(1.0e-30)/CoolingUnits);
    CoolData.temp_xray = RadiationData.ComptonXrayTemperature;
  }

  if (RadiationFieldType < 0 || RadiationFieldType > 11) {
    fprintf(stderr, "RadiationFieldType %d not recognized.\n", 
	    RadiationFieldType);
    return FAIL;
  }

  return SUCCESS;
}


