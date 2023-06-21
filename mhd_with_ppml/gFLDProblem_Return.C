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
/  Gray Flux-Limited Diffusion Implicit Problem Class return routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Returns relevant data arrays relevant to the external 
/           Enzo solver.  These arrays may need to be converted 
/           from a solver-specified format back to the finite-volume 
/           Enzo format if required.  This routine will be called 
/           repeatedly (once per time step), so it should NOT delete 
/           any problem memory.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"


int gFLDProblem::Return(EnzoVector *sol, float *Er, float *ec, float **ni) 
{

//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::Return routine\n");

  // the Enzo arrays are currently attached to the U0 vector, so 
  // unless the updated solution must be copied into U0, there 
  // is nothing to do here.  However, if this is called, just 
  // extract solution and place it in corresponding arrays.
  float **data = (float **) new float*[2+Nchem];  // create data array
  data[0] = Er;                      // attach radiation energy
  data[1] = ec;                      // attach energy correction
  for (int ns=1; ns<=Nchem; ns++)
    data[1+ns] = ni[ns-1];           // attach chemistry
  EnzoVector *tmpvec = sol->clone(data);  // create temporary vector
  tmpvec->copy(sol);                 // copy sol to data arrays
  delete tmpvec;                     // delete temporary vector;
  delete[] data;

  // return success
  return SUCCESS;
}
#endif
