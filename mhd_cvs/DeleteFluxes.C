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
/  DELETE FLUXES FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// Delete all the data associated with the fluxes structure in the argument

#include "macros_and_parameters.h"
#include "Fluxes.h"
#include "typedefs.h"
#include "global_data.h"

void DeleteFluxes(fluxes *Fluxes)
{

  int field, dim;
  for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++){
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      if (Fluxes->LeftFluxes[field][dim] != NULL) {
	delete Fluxes->LeftFluxes[field][dim];
      }
      if (Fluxes->RightFluxes[field][dim] != NULL) {
	delete Fluxes->RightFluxes[field][dim];
      }
      Fluxes->LeftFluxes[field][dim]  = NULL;
      Fluxes->RightFluxes[field][dim] = NULL;

    }
  }

if (MHD_Used)
  for(field = 0;field<3;field++){
    for(dim =0;dim<3;dim++){
      if(Fluxes->LeftElectric[field][dim] != NULL )
	delete Fluxes->LeftElectric[field][dim];
      Fluxes->LeftElectric[field][dim] = NULL;	
      
      if(Fluxes->RightElectric[field][dim] != NULL )
	delete Fluxes->RightElectric[field][dim];
      Fluxes->RightElectric[field][dim] = NULL;	
    }
  }
  
}
