/***********************************************************************
/
/  CHECKS PARAMETER FILE FOR TRACER PARTICLE CREATION PARAMETERS AND
/    IF PRESENT, CREATES TRACER PARTICLES.
/
/  written by: Greg Bryan
/  date:       April, 2004
/  modified1:
/
/  PURPOSE:  This routine scans through the parameter file pointed to
/     by the file pointer passed in and looks for parameters specifying
/     that tracer particles should be created.  If present, it creates
/     the tracer particles.  Note that unlike most parameters, these
/     parameters are not recorded because tracer particle creation is
/     a one-time event.
/
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

/* function prototypes */



int TracerParticleCreation(FILE *fptr, HierarchyEntry &TopGrid,
			   TopGridData &MetaData)
{
  /* declarations */

  char line[MAX_LINE_LENGTH];
  int dim;

  /* Set default values for parameters. */

  FLOAT TracerParticleCreationLeftEdge[MAX_DIMENSION];
  FLOAT TracerParticleCreationRightEdge[MAX_DIMENSION];
  FLOAT TracerParticleCreationSpacing = FLOAT_UNDEFINED;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    TracerParticleCreationLeftEdge[dim] = FLOAT_UNDEFINED;
    TracerParticleCreationRightEdge[dim] = FLOAT_UNDEFINED;
  }
  
  /* read until out of lines */

  rewind(fptr);
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    /* read tracer particle parameters */
    
    sscanf(line, "TracerParticleCreation = %d", &MetaData.CycleNumber);
    sscanf(line, "TracerParticleCreationSpacing = %"FSYM, 
	   &TracerParticleCreationSpacing);
    sscanf(line, "TracerParticleCreationLeftEdge = %"FSYM" %"FSYM" %"FSYM, 
		  TracerParticleCreationLeftEdge,
		  TracerParticleCreationLeftEdge+1,
		  TracerParticleCreationLeftEdge+2);
    sscanf(line, "TracerParticleCreationRightEdge = %"FSYM" %"FSYM" %"FSYM, 
		  TracerParticleCreationRightEdge,
		  TracerParticleCreationRightEdge+1,
		  TracerParticleCreationRightEdge+2);

  }

  /* If spacing is non-zero, then create particles. */

  if (TracerParticleCreationSpacing > 0)
    if (TopGrid.GridData->TracerParticleCreateParticles(
                                     TracerParticleCreationLeftEdge,
				     TracerParticleCreationRightEdge,
				     TracerParticleCreationSpacing,
				     MetaData.NumberOfParticles) == FAIL) {
      fprintf(stderr, "Error in grid->TracerParticleCreateParticles.\n");
      return FAIL;
    }

  /* clean up */

  rewind(fptr);

  return SUCCESS;
}
