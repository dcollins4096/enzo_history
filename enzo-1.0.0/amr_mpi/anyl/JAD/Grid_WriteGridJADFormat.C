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
/  GRID CLASS (WRITE OUT GRID IN JOHN'S AMR DATA (JAD) FORMAT)
/
/  written by: Greg Bryan
/  date:       April, 1996
/  modified1:
/
/  PURPOSE: see http://bach.ncsa.uiuc.edu:80/IO
/
************************************************************************/

//  Write grid to file pointer fptr
//     (we assume that the grid is at an appropriate stopping point,
//      where the Old values aren't required)
//

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "df.h"
#include "IEEEIO.hh"
#include "AMRwriter.hh"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


int grid::WriteGridJADFormat(IObase *filehandle, float TopGridCellWidth, 
			     int level, int GridID, int Resolution, 
			     int FinestLevel)
{

  /* Declarations */

  int i, j, k, dim, field, ActiveDim[MAX_DIMENSION];
  int32 OutDims[MAX_DIMENSION];

  /* Compute the active dimensions. */

  int OutSize = 1, GridSize = 1;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    ActiveDim[dim]        = GridEndIndex[dim] - GridStartIndex[dim] +1;
    OutDims[dim]          = int32(ActiveDim[dim]);
    OutSize *= ActiveDim[dim];
    GridSize *= GridDimension[dim];
  }

  /* Set top grid info on AMRwriter. */

  double origin[MAX_DIMENSION], delta[MAX_DIMENSION], timestep;
  int leftedge[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    origin[dim] = double(DomainLeftEdge[dim]);
    delta[dim] = double(TopGridCellWidth);
    leftedge[dim] = nint((GridLeftEdge[dim]-DomainLeftEdge[dim])/
			 CellWidth[dim][0]);
  }
  timestep = double(1.0);

  AMRwriter *amrfile = new AMRwriter(*filehandle);
  amrfile->setRank(GridRank);
  amrfile->setToplevelParameters(origin, delta, timestep);
  amrfile->setLevelParameters(level, Resolution, Resolution, NULL);
  
  if (NumberOfBaryonFields > 0) {

    /* Loop over fields, writing each one. */

    float32 *temp = new float32[OutSize];

    for (field = 0; field < NumberOfBaryonFields; field++) {

      /* copy active part of field into grid */

      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           + 
	       (j-GridStartIndex[1])*ActiveDim[0]              + 
	       (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		       float32(
	     BaryonField[field][i + j*GridDimension[0] +
		                     k*GridDimension[0]*GridDimension[1]]
                              );

      /* Write data to ieeeio file. */

      amrfile->setDims(OutDims);
      amrfile->setOrigin(leftedge);
      amrfile->write((VOIDP) temp);

    }   // end of loop over fields

    /* If this is cosmology, compute the temperature field as well since
       its such a pain to compute after the fact. */

    if (ComovingCoordinates) {

      /* Allocate field and compute temperature. */

      float *temperature = new float[GridSize];
      if (this->ComputeTemperatureField(temperature) == FAIL) {
	fprintf(stderr, "Error in grid->ComputeTemperatureField.\n");
	return FAIL;
      }

      /* Copy active part of field into grid */

      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           + 
	         (j-GridStartIndex[1])*ActiveDim[0]              + 
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     float32(
		   temperature[(k*GridDimension[1] + j)*GridDimension[0] + i]
			     );

      /* Write out temperature. */

      amrfile->setDims(OutDims);
      amrfile->setOrigin(leftedge);
      amrfile->write((VOIDP) temp);

      delete temperature;

    } // end: if (NumberOfBaryonFields > 0)

    /* Make sure that there is a copy of dark matter field to save
       (and at the right resolution). */

#ifdef UNUSED
    if (SelfGravity && NumberOfParticles > 0) {
      float SaveGravityResolution = GravityResolution;
      GravityResolution = 1;
      this->InitializeGravitatingMassFieldParticles();
      this->ClearGravitatingMassFieldParticles();
      this->DepositParticlePositions(this, Time, 
				     GRAVITATING_MASS_FIELD_PARTICLES);
      GravityResolution = SaveGravityResolution;
    }
#endif /* UNUSED */

    /* If present, write out the GravitatingMassFieldParticles. */

    if (GravitatingMassFieldParticles != NULL) {

      /* Set dimensions. */

      int StartIndex[MAX_DIMENSION], EndIndex[MAX_DIMENSION];
      for (dim = 0; dim < 3; dim++)
	if (GravityBoundaryType == SubGridIsolated) {
	  StartIndex[dim] = GridStartIndex[dim];
	  EndIndex[dim]   = GridEndIndex[dim];
	}
	else {
	  StartIndex[dim] = 0;
	  EndIndex[dim]   = ActiveDim[dim] - 1;
	}

      /* Copy active part of field into grid */

      for (k = StartIndex[2]; k <= EndIndex[2]; k++)
	for (j = StartIndex[1]; j <= EndIndex[1]; j++)
	  for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	    temp[(i-StartIndex[0])                           + 
	         (j-StartIndex[1])*ActiveDim[0]              + 
	         (k-StartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     float32(
			     GravitatingMassFieldParticles[ i + 
			       j*GravitatingMassFieldParticlesDimension[0] +
			       k*GravitatingMassFieldParticlesDimension[0]*
			         GravitatingMassFieldParticlesDimension[1]]
			     );

      /* Write out DM field. */
      /* Write data to amr file. */

      amrfile->setDims(OutDims);
      amrfile->setOrigin(leftedge);
      amrfile->write((VOIDP) temp);

      /* Clean up if we modified the resolution. */

      if (SelfGravity && GravityResolution != 1)
	this->DeleteGravitatingMassFieldParticles();

    } // end of (if GravitatingMassFieldParticles != NULL)

    delete temp;

  }

#ifdef UNUSED

  /* 3) Save particle quantities. */
  
  fprintf(fptr, "NumberOfParticles   = %d\n", NumberOfParticles);

  if (NumberOfParticles > 0) {

    fprintf(fptr, "ParticleFileName = %s\n", name);

    /* Sort particles according to their identifier. */

    float *DragList[2*MAX_DIMENSION+1];
    DragList[0] = ParticleMass;
    for (dim = 0; dim < GridRank; dim++) {
      DragList[2*dim+1] = ParticlePosition[dim];
      DragList[2*dim+2] = ParticleVelocity[dim];
    }
    QuickSortAndDrag(ParticleNumber, 0, NumberOfParticles-1, 2*GridRank+1,
		     DragList);
    
    /* Clear any HDF parameters set previously. */

    DFSDclear();

    /* Create a temporary buffer (32 bit). */

    float32 *temp = new float32[NumberOfParticles];
    int32 TempIntArray[1];
    TempIntArray[0] = int32(NumberOfParticles);

    /* Copy particle positions to temp and write them. */

    for (dim = 0; dim < GridRank; dim++) {

      for (i = 0; i < NumberOfParticles; i++)
	*(temp + i) = float32(*(ParticlePosition[dim] + i));

      if (DFSDadddata(name, 1, TempIntArray, (VOIDP) temp) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDadddata (particle position).\n");
	return FAIL;
      }

    }

    /* Copy particle velocities to temp and write them. */

    for (dim = 0; dim < GridRank; dim++) {

      for (i = 0; i < NumberOfParticles; i++)
	*(temp + i) = float32(*(ParticleVelocity[dim] + i));

      if (DFSDadddata(name, 1, TempIntArray, (VOIDP) temp) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDadddata (particle position).\n");
	return FAIL;
      }

    }

    /* Copy mass to temp and write it. */

    for (i = 0; i < NumberOfParticles; i++)
      *(temp + i) = float32(*(ParticleMass + i));

    if (DFSDadddata(name, 1, TempIntArray, (VOIDP) temp) == HDF_FAIL) {
      fprintf(stderr, "Error in DFSDadddata (particles).\n");
      return FAIL;
    }

    /* Copy number (ID) to temp and write it. */

    int32 *tempint = new int32[NumberOfParticles];
    for (i = 0; i < NumberOfParticles; i++)
      *(tempint + i) = int32(*(ParticleNumber + i));

    if (DFSDadddata(name, 1, TempIntArray, (VOIDP) tempint) == HDF_FAIL) {
      fprintf(stderr, "Error in DFSDadddata (particles).\n");
      return FAIL;
    }

    /* clean up */

    delete temp;
    delete tempint;

  }

  /* 4) Save Gravity info. */

  if (SelfGravity)
    fprintf(fptr, "GravityBoundaryType = %d\n", GravityBoundaryType);

#endif /* UNUSED */

  delete amrfile;
  
  return SUCCESS;

}
