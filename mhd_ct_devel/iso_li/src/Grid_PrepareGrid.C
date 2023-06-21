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
/  GRID CLASS (PREPARES GRID FOR INITIALIZATION)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:  Sets some basic grid quantities based on the arguments.
/
/  INPUTS:
/    GridDim   - int array of dimensions (includes ghost zones unless dim = 1)
/    LeftEdge  - float array of left edge starting positions in terms of
/                problem space (i.e. starting 'vector').  Doesn't include
/                ghost zones.
/    RightEdge - stop (right) 'vector'.
/
************************************************************************/

//
//  Assign basic values to a grid (allocate fields)
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);

void grid::PrepareGrid(int Rank, int GridDim[], 
		       FLOAT LeftEdge[], FLOAT RightEdge[], int NumParticles)
{

  int field, dim;

  /* debugging code */

  /* Set Particle quantities. */

  NumberOfParticles = NumParticles;

  /* Set global grid quantities.
     (Start/End index are zero based). */

  GridRank = Rank;
  
  for (dim = 0; dim < GridRank; dim++) {
    GridDimension[dim]  = GridDim[dim];
    GridStartIndex[dim] = min(DEFAULT_GHOST_ZONES, GridDim[dim]-1);
    GridEndIndex[dim]   = min(abs(GridDim[dim]-DEFAULT_GHOST_ZONES-1), 
			      GridDim[dim]-1);
    GridLeftEdge[dim]   = LeftEdge[dim];
    GridRightEdge[dim]  = RightEdge[dim];
  }

  if(MHD_Used)
    {
      
      for(field=0; field<3; field++){
	MagneticSize[field] = 1;
	ElectricSize[field] = 1;

	for(dim=0; dim<3; dim++){
          MagneticDims[field][dim] = GridDimension[dim];
          ElectricDims[field][dim] = GridDimension[dim] +1;
	  

	  MHDStartIndex[field][dim] = GridStartIndex[dim];
	  MHDEndIndex[field][dim] = GridEndIndex[dim];

	  MHDeStartIndex[field][dim] = GridStartIndex[dim];
	  MHDeEndIndex[field][dim] = GridEndIndex[dim]+1;

          if( field == dim )
            {MagneticDims[field][dim]++;
            ElectricDims[field][dim]--;
	    MHDEndIndex[field][dim]++;
	    MHDeEndIndex[field][dim]--;}
	  

	  MagneticSize[field] *= MagneticDims[field][dim];
	  ElectricSize[field] *= ElectricDims[field][dim];
	}
      }

      //fprintf(stderr, "\n\n");
      for( field=0;field<3;field++)
	for(dim=0;dim<3;dim++){
	  MHDAdd[field][dim]=MagneticDims[field][dim]-GridDimension[dim];
	  //fprintf(stderr, "MHDAdd[%d][%d] = %d ", field, dim, MHDAdd[field][dim]);
	}
      //fprintf(stderr, "\n\n");

      
      

      
    }



    /* compute derived quantites */

  this->PrepareGridDerivedQuantities();

}

