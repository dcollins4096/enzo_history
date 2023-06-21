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
/  Gray Flux-Limited Diffusion Implicit Problem Class 
/  Species cross-section evaluation routine 
/
/  written by: Daniel Reynolds
/  date:       March, 2007
/  modified1:  
/
/  PURPOSE: Takes in a frequency and returns the appropriate cross
/  section at that frequency.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"

 
 

float gFLDProblem::CrossSections(float nu, int species)
{

  // set necessary constants
  float h = 6.6260693e-27;   // Planck's constant [ergs*s]
  float nu0_HI   = 13.6/h;   // ionization threshold of HI
  float nu0_HeI  = 20.4/h;   // ionization threshold of HeI
  float nu0_HeII = 54.4/h;   // ionization threshold of HeII
  float nuscaled;            // normalized frequency
  float eps;                 // constant in cross section definition
  float pi = 4.0*atan(1.0);
  float sigma;

  // evaluate the radiation spectrum based on the internal ESpectrum parameter
  switch (species) {


  case 0:  // sigma_HI

    // check that frequency is within the allowed range
    if (nu < nu0_HI) {
      fprintf(stderr,"gFLDProblem::CrossSections Error: (nu = %g) < (nu0 = %g)\n",
	      nu,nu0_HI);
      sigma = -1.0;
      break;
    }
    nuscaled = nu/nu0_HI;
    eps = sqrt(nuscaled - 1.0);
    sigma = (6.30e-18)*POW(nuscaled,-4.0)
      *exp(4.0-4.0*atan(eps)/eps)/(1.0-exp(-2.0*pi/eps));
    break;


  case 1:  // sigma_HeI

    // check that frequency is within the allowed range
    if (nu < nu0_HeI) {
      fprintf(stderr,"gFLDProblem::CrossSections Error: (nu = %g) < (nu0 = %g)\n",
	      nu,nu0_HeI);
      sigma = -1.0;
      break;
    }
    nuscaled = nu/nu0_HeI;
    sigma = (7.42e-18)*(1.66*POW(nuscaled,-2.05) - 0.66*POW(nuscaled,-3.05));
    break;


  case 2:  // sigma_HeII

    // check that frequency is within the allowed range
    if (nu < nu0_HeII) {
      fprintf(stderr,"gFLDProblem::CrossSections Error: (nu = %g) < (nu0 = %g)\n",
	      nu,nu0_HeII);
      sigma = -1.0;
      break;
    }
    nuscaled = nu/nu0_HeII;
    eps = sqrt(nuscaled - 1.0);
    sigma = (1.575e-18)*POW(nuscaled,-4.0)
      *exp(4.0-4.0*atan(eps)/eps)/(1.0-exp(-2.0*pi/eps));
    break;


  default:   // illegal input
    fprintf(stderr,"gFLDProblem::CrossSections Error: species %"ISYM" undefined\n",
	    species);
    sigma = -1.0;
  }

  return sigma;
}

#endif
