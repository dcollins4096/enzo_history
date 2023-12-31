/***********************************************************************
/
/  EVALUATE THE POWER SPECTRUM AT THE SPECIFIED K VALUE (IN 1/MPC)
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "PowerSpectrumParameters.h"

/* function prototypes */

extern "C" float FORTRAN_NAME(psfunc)
             (float *k, int *ispectrum, float *omega0,
	         float *hub, float *omega_nu,
	      float *omega_lam, float *psindex, float *omegab0, float *z, 
	         float *gamma,
	      float *psnorm, float *growth_factor, float *kcutoff);
extern "C" void FORTRAN_NAME(spline)(float *x, float *y, int *n, float *yp1,
				     float *ypn, float *y2);
extern "C" void FORTRAN_NAME(splint)(float *xa, float *ya, float *y2a, 
				     int *n, float *x, float *y);

static int CMBFastAlreadyInput = {0};
static int CMBFastNumberOfPoints[2];
static float *CMBFastWaveNumber[2];
static float *CMBFastTransferFunction[2][MAX_SPECIES];
static float *CMBFastTransferSpline[2][MAX_SPECIES];

/* Prototypes for Eisenstein & Hu fitting functions. */

extern int TFmdm_set_cosm(float omega_matter, float omega_baryon,
			      float omega_hdm, int degen_hdm, 
			      float omega_lambda, float hubble, 
			      float redshift); 
extern float TFmdm_onek_mpc(float kk);


float EvaluatePowerSpectrum(float k, int Species)
{

  /* Species indicates which p.s. to use:
     0 - total (mean) p.s.
     1 - CDM
     2 - baryon */

  int i, j, Table, n1;
  float psval, dummy, x0, x1, y0, y1;
  FILE *fptr;

  /* ------------------------------------------------------------------------
     1-10) The first 10 types are reserved for use in the generalized routine
     psfunc. */

  if (PowerSpectrumType < 10)
      psval = FORTRAN_NAME(psfunc)
	                  (&k, &PowerSpectrumType, &OmegaMatterNow,
			      &HubbleConstantNow, &OmegaHDMNow, 
			   &OmegaLambdaNow, &PrimordialIndex, &OmegaBaryonNow,
			      &Redshift, &Gamma, 
			   &Normalization, &GrowthFactor, &kcutoff);

  /* ------------------------------------------------------------------------
     11) This is the Eisenstein & Hu fitting function with 1 neutrino
     12)   (with two neutrinos). */
			   
  else if (PowerSpectrumType == 11 || PowerSpectrumType == 12) {

    /* Error check. */

    if (kcutoff != 0) {
      fprintf(stderr, "This power spectrum does not support kcutoff.\n");
      exit(EXIT_FAILURE);
    }

    /* Initialize cosmology. */

    int NumberOfNeutrinoSpecies = PowerSpectrumType - 10;
    if (TFmdm_set_cosm(OmegaMatterNow, OmegaBaryonNow, OmegaHDMNow,
		       NumberOfNeutrinoSpecies, OmegaLambdaNow, 
		       HubbleConstantNow, Redshift) != 0) {
      fprintf(stderr, "Error in TFmdm_set_cosm.\n");
      exit(EXIT_FAILURE);
    }

    /* Compute transfer function. */

    psval = TFmdm_onek_mpc(k);

    /* Multiply transfer function by power spectrum */

    psval = Normalization * POW(k, PrimordialIndex) * 
                            POW(GrowthFactor*psval, 2);

  }

  /* -----------------------------------------------------------------
     Warm Dark Matter - use modified Eisenstein & Hu with 
     one neutrino species!
     ----------------------------------------------------------------- */
  else if (PowerSpectrumType == 13){

    /* Error check. */
    
    if (kcutoff != 0) {
      fprintf(stderr, "This power spectrum does not support kcutoff.\n");
      exit(EXIT_FAILURE);
    }

    /* Initialize cosmology. */

    int NumberOfNeutrinoSpecies = 1;
    if (TFmdm_set_cosm(OmegaMatterNow, OmegaBaryonNow, OmegaHDMNow,
		       NumberOfNeutrinoSpecies, OmegaLambdaNow, 
		       HubbleConstantNow, Redshift) != 0) {
      fprintf(stderr, "Error in TFmdm_set_cosm.\n");
      exit(EXIT_FAILURE);
    }

    /* Compute transfer function. */

    psval = TFmdm_onek_mpc(k);

    /* compute fitting parameter for WDM power spectrum */
    float WDMAlpha = 0.05 *
      POW( (OmegaWDMNow/0.4), 0.15 ) *
      POW( (HubbleConstantNow/0.65), 1.3 ) *
      POW( WDMPartMass, -1.15 ) *
      POW( (1.5/WDMg_x), 0.29 );

    /* Multiply transfer function by power spectrum */
    psval = Normalization * POW(k, PrimordialIndex) * 
                            POW(GrowthFactor*psval, 2);

    psval *= POW( (1.0 + (WDMAlpha*k)*(WDMAlpha*k)), -10.0);

  }

  /* -----------------------------------------------------------------
     Warm Dark Matter - Equivalent to PowerSpectrumType = 1, but 
     using the WDM power spectrum modifier.
     ----------------------------------------------------------------- */

  else if (PowerSpectrumType == 14){

    int WDMPSType = 1;  // corresponds to power spec. type #1

      psval = FORTRAN_NAME(psfunc)
	                  (&k, &WDMPSType, &OmegaMatterNow,
			      &HubbleConstantNow, &OmegaHDMNow, 
			   &OmegaLambdaNow, &PrimordialIndex, &OmegaBaryonNow,
			      &Redshift, &Gamma, 
			   &Normalization, &GrowthFactor, &kcutoff);

    /* compute fitting parameter for WDM power spectrum */
    float WDMAlpha = 0.05 *
      POW( (OmegaWDMNow/0.4), 0.15 ) *
      POW( (HubbleConstantNow/0.65), 1.3 ) *
      POW( WDMPartMass, -1.15 ) *
      POW( (1.5/WDMg_x), 0.29 );

    psval *= POW( (1.0 + (WDMAlpha*k)*(WDMAlpha*k)), -10.0);

  }

  /* ------------------------------------------------------------------------
     20) Read in transfer function from CMBFAST. */

  else if (PowerSpectrumType == 20) {

    if (!CMBFastAlreadyInput) {

      char char_dummy[MAX_LINE_LENGTH];

      /* Loop over the required two tables, 0 is for z=0, 1 is for z=zinit. */

      for (i = 0; i < 2; i++) {

	/* Open CMBFast transfer function. */

	if ((fptr = fopen(PSFileName[i], "r")) == NULL) {
	  fprintf(stderr, "Error opening CMBFast[%d]: %s\n", i, PSFileName[i]);
	  exit(EXIT_FAILURE);
	}

	/* Count lines and allocate space. */
	
	CMBFastNumberOfPoints[i] = 0;
	while (fscanf(fptr, "%f", &dummy) > 0) {
	  CMBFastNumberOfPoints[i]++;
	  fgets(char_dummy, MAX_LINE_LENGTH, fptr);
	}
	if (debug) printf("CMBFast[%d] NumberOfPoints = %d\n", i, 
			  CMBFastNumberOfPoints[i]);
	CMBFastWaveNumber[i]          = new float[CMBFastNumberOfPoints[i]];
	for (j = 0; j < 3; j++) {
	  CMBFastTransferFunction[i][j] = new float[CMBFastNumberOfPoints[i]];
	  CMBFastTransferSpline[i][j] = new float[CMBFastNumberOfPoints[i]];
	}
	rewind(fptr);

	/* Read data: k/h, T(CDM), T(baryon), and compute T(mean) =
	   f_cdm*T(CDM) + f_baryon*T(baryon). */

	j = 0;
	while (fscanf(fptr, "%f %f %f", CMBFastWaveNumber[i]+j,
		      CMBFastTransferFunction[i][1]+j,
		      CMBFastTransferFunction[i][2]+j) == 3) {
	  CMBFastTransferFunction[i][0][j] = 
	    (OmegaMatterNow-OmegaBaryonNow)/OmegaMatterNow*
	        CMBFastTransferFunction[i][1][j] +
	    (               OmegaBaryonNow)/OmegaMatterNow*
		CMBFastTransferFunction[i][2][j];
	  CMBFastWaveNumber[i][j++] *= HubbleConstantNow;
	  fgets(char_dummy, MAX_LINE_LENGTH, fptr);  // get rid of rest of line
	}
	if (j != CMBFastNumberOfPoints[i]) {
	  fprintf(stderr, "Counting error (%d/%d) in %s\n", j, 
		  CMBFastNumberOfPoints[i], PSFileName[i]);
	  exit(EXIT_FAILURE);
	}

	/* generate spline table. */

	float BigNum = 1e32;
	for (j = 0; j < 3; j++)
	  FORTRAN_NAME(spline)(CMBFastWaveNumber[i], 
		  CMBFastTransferFunction[i][j], &CMBFastNumberOfPoints[i],
		  &BigNum, &BigNum, CMBFastTransferSpline[i][j]);

	/* Close file. */

	fclose(fptr);

      } // end loop over i

      CMBFastAlreadyInput = TRUE;

    } // end: if (!CMBFastAlreadyInput)

    /* If the redshift = 0, use table 0, otherwise use table 1. */

    Table = (Redshift == 0) ? 0 : 1;
    n1 = CMBFastNumberOfPoints[Table];

    /* Do a spline look-up to get transfer function. 
       (if k is too low, just use 1, if too high to a log extrapolation). */

    if (k < CMBFastWaveNumber[Table][0])
      psval = 1;
    else if (k > CMBFastWaveNumber[Table][n1-1]) {
      x0 = log(CMBFastWaveNumber[Table][n1-2]);
      x1 = log(CMBFastWaveNumber[Table][n1-1]);
      y0 = log(CMBFastTransferFunction[Table][Species][n1-2]);
      y1 = log(CMBFastTransferFunction[Table][Species][n1-1]);
      psval = exp(y0 + (log(k)-x0)/(x1-x0)*(y1-y0));
    }
    else
      FORTRAN_NAME(splint)(CMBFastWaveNumber[Table], 
			   CMBFastTransferFunction[Table][Species], 
			   CMBFastTransferSpline[Table][Species],
			   &CMBFastNumberOfPoints[Table], &k, &psval);

    /* Multiply transfer function by power spectrum. */

//    printf("%g %g\n", k, psval);
    psval = Normalization * POW(k, PrimordialIndex) * 
      (GrowthFactor*psval)*(GrowthFactor*psval); /*(before 7/5/98)*/
    //    psval = Normalization * POW(k, PrimordialIndex) * psval * psval;

  } // end: if (PowerSpectrumType == 20)

  /* ------------------------------------------------------------------------
     Unidentified. */

  else {
    fprintf(stderr, "PowerSpectrumType = %d unknown.\n", PowerSpectrumType);
    exit(EXIT_FAILURE);
  }

  return psval;
}
