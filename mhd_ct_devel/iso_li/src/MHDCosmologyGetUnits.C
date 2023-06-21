/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  COMPUTE AND RETURN THE MHD COSMOLOGY UNITS
/
/
/  Modified from CosmologyGetUnits.C
/  written by: Greg Bryan
/  date:       April, 1995
/  modified by: Hao Xu
/  date: March,2006
/
/  PURPOSE:  Returns the cosmology units:
/
/         time:        utim = 1 / sqrt(4 * \pi * G * \rho_0 * (1+zri)^3)
/         density:     urho = \rho_0 * (1+z)^3
/         length:      uxyz = (1 Mpc) * box / h / (1+z)
/         velocity:    uvel = uaye * uxyz / utim  (since u = a * dx/dt)
/    (*)  temperature: utem = m_H * \mu / k * uvel**2
/         a(t):        uaye = 1 / (1 + zri)
/
/         Magnetic Field: ub = sqrt(\rho_0) * uvel * (1+z)^2 
/
/           where:
/             box     - size of simulation box in Mpc/h
/             zri     - initial redshift (start of simulation)
/             \rho_0  = 3*\Omega_0*H_0^2/(8*\pi*G)
/             Omega_0 - the fraction of non-relativistic matter at z=0
/
/           Note that two definitions are dependent on redshift (urho
/             and uxyz) so make sure to call this routine immediately
/             before writing.
/
/           * - the utem given below assumes that \mu = 1, so you must
/               multiply the resulting temperature field by \mu.
/
/
/
/  NOTE: 
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "CosmologyParameters.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int MHDCosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time,
                      float *BFieldUnits)
{

  /* From the time, compute the current redshift. */

  FLOAT a, dadt;
  if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
    fprintf(stderr, "Error in ComputeExpansionFactor.\n");
    return FAIL;
  }

  /* Compute the current redshift (remember a(init) = 1). */

  FLOAT CurrentRedshift = (1 + InitialRedshift)/a - 1;

  /* Determine the units. */

  *DensityUnits     = 1.88e-29*OmegaMatterNow*POW(HubbleConstantNow,2)*
                      POW(1 + CurrentRedshift,3);

  *LengthUnits      = 3.086e24*ComovingBoxSize/HubbleConstantNow/
                      (1 + CurrentRedshift);

  *TemperatureUnits = 1.88e6*POW(ComovingBoxSize,2)*OmegaMatterNow*
                      (1 + InitialRedshift);

  *TimeUnits        = 2.52e17/sqrt(OmegaMatterNow)/HubbleConstantNow/
                      POW(1 + InitialRedshift,FLOAT(1.5));

  *VelocityUnits    = 1.225e7*ComovingBoxSize*sqrt(OmegaMatterNow)*
                      sqrt(1 + InitialRedshift);

//  *BFieldUnits      = sqrt((*DensityUnits)/POW(1 +CurrentRedshift,3))
//                      *(*VelocityUnits)*POW(1 + CurrentRedshift,2);

    *BFieldUnits      = sqrt(4*3.14159*(*DensityUnits))*(*VelocityUnits)/sqrt(a);
  if(MHD_Equation == 2) 
   *BFieldUnits = (*BFieldUnits)*sqrt(a);  

  return SUCCESS;
}
