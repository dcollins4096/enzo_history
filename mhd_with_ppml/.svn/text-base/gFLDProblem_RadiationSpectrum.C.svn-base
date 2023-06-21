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
/  Radiation Spectrum Evaluation routine 
/
/  written by: Daniel Reynolds
/  date:       March, 2007
/  modified1:  
/
/  PURPOSE: Takes in a frequency and returns the assumed radiation 
/  energy spectrum at that frequency.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"

 
 

float gFLDProblem::RadiationSpectrum(float nu)
{

  // set necessary constants
  float h = 6.6260693e-27;   // Planck's constant [ergs*s]
  float nu0 = 13.6/h;        // ionization threshold of Hydrogen
  float sigma;

  // check that frequency is within the allowed range
  if (nu < nu0) {
    fprintf(stderr,"gFLDProblem::RadiationSpectrum Error: (nu = %g) < (nu0 = %g)\n",
	    nu,nu0);
    return -1.0;
  }

  // evaluate the radiation spectrum based on the internal ESpectrum parameter
  switch (ESpectrum) {

  // Add new spectrum choices here
  // case 0:
  //   ...
  //   break;

  default:
    // simple power law spectrum with power -1.5
    sigma = POW((nu/nu0),-1.5);
  }

  return sigma;
}

#endif
