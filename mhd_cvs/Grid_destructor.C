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
/  GRID CLASS (DESTRUCTOR)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//
//  Grid destructor
//
#include <stdio.h>
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

void DeleteFluxes(fluxes *Fluxes);
void WriteListOfInts(FILE *fptr, int N, int nums[]);


grid::~grid()
{

  int i;

  /* Error check. */

#ifdef UNUSED
  if (NumberOfParticles > 0) {
    fprintf(stderr, "warning: destroying live particles (%d).\n", 
	    NumberOfParticles);
  /* exit(EXIT_FAILURE); */
  }
#endif /* UNUSED */

  for (i = 0; i < MAX_DIMENSION; i++) {
    delete CellLeftEdge[i];
    delete CellWidth[i];
    delete ParticlePosition[i];
    delete ParticleVelocity[i];
    delete ParticleAcceleration[i];
    delete AccelerationField[i];
    delete RandomForcingField[i];
  }

  delete ParticleAcceleration[MAX_DIMENSION];

  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    if( BaryonField[i] != NULL ){
      delete BaryonField[i];
      BaryonField[i] = NULL;
    }
    if( OldBaryonField[i] != NULL ){
      delete OldBaryonField[i];
      OldBaryonField[i] = NULL;
    }
  }


  for(i=0; i<3; i++)
    if(OldAccelerationField[i] != NULL ){
      delete OldAccelerationField[i];
      OldAccelerationField[i] = NULL;
    }

  //MHD stuff 

  if( MHD_Used){
    for(i=0;i<3;i++){

      if(MagneticField[i] != NULL ){
	delete MagneticField[i];
	MagneticField[i] = NULL;
      }
      if(OldMagneticField[i] != NULL ){
	delete OldMagneticField[i];
	OldMagneticField[i] = NULL;
      }

      if( CenteredB[i] != NULL ){
	delete CenteredB[i];
	CenteredB[i] = NULL;
      }

      if( DivB != NULL ){
        delete DivB;
        DivB = NULL;
      }

      if(ElectricField[i] != NULL){
	delete ElectricField[i];
	ElectricField[i] = NULL;
      }

      if(AvgElectricField[i] != NULL){
	delete AvgElectricField[i];
	AvgElectricField[i] = NULL;
      }

      if( DivB != NULL ){
	delete DivB;
	DivB = NULL;
      }
      if( Current[i] != NULL ){
	delete Current[i];
	Current[i] = NULL;
      }
    }

    MHDCleanUpTemp();

  }

  DeleteFluxes(&BoundaryFluxes);
  delete ParticleMass;
  delete ParticleNumber;
  delete PotentialField;
  delete GravitatingMassField;
  delete GravitatingMassFieldParticles;
  delete FlaggingField;
  delete MassFlaggingField;


  for (i = 0; i < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; i++)
    delete ParticleAttribute[i];

  if (debug && GridRank > 0) {
    printf("grid->destructor: deleting grid with dims = ");
    WriteListOfInts(stdout, GridRank, GridDimension);
  }

}
