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
/  EVOLVE LEVEL ROUTINES (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modified1:
/
/  PURPOSE:  This is a collection of routines called by EvolveLevel.
/            These have been optimized for enhanced message passing 
/            performance by performing two passes -- one which generates
/            sends and the second which receives them.
/
************************************************************************/


#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "pout.h"
/* function prototypes */

int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
int PrepareGravitatingMassField1(HierarchyEntry *Grid);
int PrepareGravitatingMassField2(HierarchyEntry *Grid, int grid1,
				 SiblingGridList SiblingList[],
				 TopGridData *MetaData, int level);
#else
int PrepareGravitatingMassField(HierarchyEntry *Grid, TopGridData *MetaData,
                                LevelHierarchyEntry *LevelArray[], int level);
#endif
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);



extern int CopyPotentialFieldAverage;



/* ======================================================================= */
/* This routine sets all the boundary conditions for Grids by either
   interpolating from their parents or copying from sibling grids. */

int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
                          SiblingGridList SiblingList[],
#endif
			  int level, TopGridData *MetaData, 
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level)
{
  int grid, grid2;
#ifndef OLD_CENTER    
  for (grid = 0; grid < NumberOfGrids; grid++) {

    //dcc 01/24/06.  MagneticField needs to be centered before the boundary condition 
    //               are set for data consistancy.  

    if( MHD_Used == TRUE ){
      if( Grids[grid]->GridData->CenterMagneticField() == FAIL ){
	fprintf( stderr, "error in Center MagneticField\n");
	return FAIL;
      }
    }
  }//grids
#endif //OLD_CENTER
  /* -------------- FIRST PASS ----------------- */
  CommunicationDirection = COMMUNICATION_SEND;

  for (grid = 0; grid < NumberOfGrids; grid++) {

    /* a) Interpolate boundaries from the parent grid or set external
       boundary conditions. */
    
    if (level == 0) {
      
      if (Grids[grid]->GridData->SetExternalBoundaryValues(Exterior) 
	  == FAIL) {
	fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
	return FAIL;
      }

#ifdef EMF_BOUNDARY_off
    //update magnetic fileds after aplying emf boundary conditions
     if(EMF_Boundary==1)
         if (Grids[grid]->GridData->MHD_UpdateMagneticField(0, 0) == FAIL) {
      fprintf(stderr, "Error in grid->MHD_UpdateMagneticField.\n");
      return FAIL;
    }
#endif //EMF_BOUNDARY

      
    }
    else {

      if ((Grids[grid]->GridData->InterpolateBoundaryFromParent
	   (Grids[grid]->ParentGrid->GridData, Level)) == FAIL) {
	fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");
	return FAIL;
      }
    }
    
    /* b) Copy any overlapping zones for sibling grids.  */
    
      
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH

    for (grid2 = 0; grid2 < SiblingList[grid].NumberOfSiblings; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(SiblingList[grid].GridList[grid2],
						 MetaData->LeftFaceBoundaryCondition,
						 MetaData->RightFaceBoundaryCondition,
						 &grid::CopyZonesFromGrid)
	  == FAIL) {
	fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
      }

#else

    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
						 MetaData->LeftFaceBoundaryCondition,
						 MetaData->RightFaceBoundaryCondition,
						 &grid::CopyZonesFromGrid)
	  == FAIL) {
	fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
      }
#endif

  
    
    
  } // end loop over grids
  //fprintf(stderr, "XXX: SBC: EndOfFirstPasss\n");
  /* -------------- SECOND PASS ----------------- */
  
  CommunicationDirection = COMMUNICATION_RECEIVE;
  

    for (grid = 0; grid < NumberOfGrids; grid++) {

	/* a) Interpolate boundaries from the parent grid or set external
	   boundary conditions. */
	
	if (level > 0){
	  
	  if ((Grids[grid]->GridData->InterpolateBoundaryFromParent
	       (Grids[grid]->ParentGrid->GridData,Level)) == FAIL) {
	    fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");
	    return FAIL;
	  }
	}
	
      /* b) Copy any overlapping zones for sibling grids.  */
      
    	
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
	for (grid2 = 0; grid2 < SiblingList[grid].NumberOfSiblings; grid2++)
	  if (Grids[grid]->GridData->CheckForOverlap(SiblingList[grid].GridList[grid2],
						     MetaData->LeftFaceBoundaryCondition,
						     MetaData->RightFaceBoundaryCondition,
						     &grid::CopyZonesFromGrid)
	      == FAIL) {
	    fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
	  }
	
#else
	for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	  if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
						     MetaData->LeftFaceBoundaryCondition,
						     MetaData->RightFaceBoundaryCondition,
						     &grid::CopyZonesFromGrid)
	      == FAIL) {
	    fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
	  }
	
#endif
    
    } // end loop over grids
    
    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
    
    // some IO, for debugging purposes.
    if(WriteInThisF(60)==TRUE ){
      
      char basename[30];
      fprintf(stderr," SBC: NumberOfGrids = %d \n", NumberOfGrids);
      for(grid=0;grid<NumberOfGrids;grid++){
	
	//data60<iteration><level><grid>.grid<ProcessorNumber>
	sprintf(basename,"data60%d%d%d.grid", dccCounter3, level,grid);
	
	FILE * dccptr = fopen(basename,"a");
	fprintf(stderr,"%s counter %d proc %d\n",basename, dccCounter3, MyProcessorNumber);
	Grids[grid]->GridData->WriteGrid(dccptr,basename,(MyProcessorNumber+1));
	fclose(dccptr);
      }    
      dccCounter3++;
    }
    
    return SUCCESS;
}

int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
			    SiblingGridList SiblingList[],
#endif
			    int level, TopGridData *MetaData, 
			    ExternalBoundary *Exterior, LevelHierarchyEntry * Level,
			    int CycleNumber){
	    

  int DontIndentThisFloat;	  
  //The previous float is here to keep the automatic indentin in emacs working.  
  //The ifdefs for the FastSiblingLocator make it to go haywire.

  //Set the boundary on the Acceleration field.  Reuse SetBoundaryConditions.  
  //Juggle pointers around.

  int grid, ConservativeTruth, MHD_Truth;
  char basename[30];  
  
  //We don't want conservative interpolation actually being done for the acceleration field.
  //We also don't want to bother re-interpolating the MHD fields.
  ConservativeTruth = ConservativeInterpolation;
  ConservativeInterpolation = FALSE;
  MHD_Truth = MHD_Used;
  MHD_Used = FALSE;

  for (grid = 0; grid < NumberOfGrids; grid++) {
    if( Grids[grid]->GridData->AttachAcceleration() == FAIL ) {
      fprintf(stderr,"Error in AttachAcceleration \n");
      return FAIL;
    }

    if( level > 0 ){
      if( Grids[grid]->ParentGrid->GridData->AttachAcceleration() ==FAIL ){
	fprintf(stderr,"Error in AttachAcceleration, Parent \n");
	return FAIL;
      }
    }


    if( level > 0 ){
      if( WriteInThisF(61) == TRUE) {
	sprintf(basename, "data61%d%d.grid",level+1,CycleNumber);
	
	FILE *dummy = fopen(basename, "a");    
	if( Grids[grid]->ParentGrid->GridData->WriteGrid(dummy, basename, grid+1) == FAIL ){
	  fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
	  return FAIL;
	}
	fclose(dummy);
	
      }  
    }//level > 0

    if( WriteInThisF(62) == TRUE) {
      sprintf(basename, "data62%d%d.grid",level+1,CycleNumber);
      
      FILE *dummy = fopen(basename, "a");    
      if( Grids[grid]->GridData->WriteGrid(dummy, basename, grid+1) == FAIL ){
	fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
	return FAIL;
      }
      fclose(dummy);
      
    }  
  }

  if (SetBoundaryConditions(Grids, NumberOfGrids, 
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
			    SiblingList,
#endif
			    level, MetaData, 
			    Exterior, NULL) == FAIL)
    return FAIL;


  for (grid = 0; grid < NumberOfGrids; grid++) {

    if( level > 0 ){
      if( WriteInThisF(80) == TRUE) {
	sprintf(basename, "data80%d%d.grid",level+1,CycleNumber);
	FILE *dummy = fopen(basename, "a");    
	if( Grids[grid]->ParentGrid->GridData->WriteGrid(dummy, basename, grid+1) == FAIL ){
	  fprintf(stderr, "Shit.  Problem with Write Grid in SAB.\n");
	  return FAIL;
	}
	fclose(dummy);
      }  
    }

    if( WriteInThisF(81) == TRUE) {
      sprintf(basename, "data81%d%d.grid",level+1,CycleNumber);
      FILE *dummy = fopen(basename, "a");    
      if( Grids[grid]->GridData->WriteGrid(dummy, basename, grid+1) == FAIL ){
	fprintf(stderr, "Shit.  Problem with Write Grid in SAB.\n");
	return FAIL;
      }
      fclose(dummy);
    }  

    if( Grids[grid]->GridData->DetachAcceleration() == FAIL ) {
      fprintf(stderr,"Error in DetachAcceleration\n");
      return FAIL;
    }

    if( level > 0 ){
      if( Grids[grid]->ParentGrid->GridData->DetachAcceleration() == FAIL ) {
	fprintf(stderr,"Error in DetachAcceleration, parent\n");
	return FAIL;
      }
    }



  }


  ConservativeInterpolation = ConservativeTruth;
  MHD_Used = MHD_Truth;

  return SUCCESS;

}

/* ======================================================================= */
/* This routine prepares the density field for all the grids on this level,
   both particle and baryonic densities.  It also calculates the potential
   field if this is level 0 (since this involves communication). */


int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
                          SiblingGridList SiblingList[],
#endif
                        int level, TopGridData *MetaData)
{


  //=================================================================

  int grid, grid2;

  /* Set the time for evaluation of the fields, etc. */

  FLOAT EvaluateTime = LevelArray[level]->GridData->ReturnTime() +
                   0.5*LevelArray[level]->GridData->ReturnTimeStep();

  /* If level is above MaximumGravityRefinementLevel, then just update the
     gravity at the MaximumGravityRefinementLevel. */

  int reallevel = level;
  level = min(level, MaximumGravityRefinementLevel);

  /* Create an array (Grids) of all the grids. */

  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

  /* Grids: Deposit particles in their GravitatingMassFieldParticles. */

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  CommunicationDirection = COMMUNICATION_SEND;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (DepositParticleMassField(Grids[grid], EvaluateTime) == FAIL) {
      fprintf(stderr, "Error in DepositParticleMassField.\n");
      return FAIL;
    }

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  CommunicationDirection = COMMUNICATION_RECEIVE;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (DepositParticleMassField(Grids[grid], EvaluateTime) == FAIL) {
      fprintf(stderr, "Error in DepositParticleMassField.\n");
      return FAIL;
    }

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  /* Grids: compute the GravitatingMassField (baryons & particles). */

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
  const int passes = 2;
#else
  const int passes = 1;
#endif

  for (int pass = 0; pass < passes; pass++) {
    
    CommunicationDirection = COMMUNICATION_SEND;
    for (grid = 0; grid < NumberOfGrids; grid++) {
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
      if ((pass==0 && PrepareGravitatingMassField1(Grids[grid]) == FAIL) ||
          (pass==1 && PrepareGravitatingMassField2(Grids[grid], grid, SiblingList,MetaData, level) == FAIL))
#else
      if (PrepareGravitatingMassField(Grids[grid], MetaData, LevelArray,
                                        level) == FAIL) 
#endif
        {
          fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
          return FAIL;
        }
      }

    CommunicationDirection = COMMUNICATION_RECEIVE;
    for (grid = 0; grid < NumberOfGrids; grid++) {
      if 
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
        ((pass==0 && PrepareGravitatingMassField1(Grids[grid]) == FAIL) ||
         (pass==1 && PrepareGravitatingMassField2(Grids[grid], grid, SiblingList,MetaData, level) == FAIL))
#else
        (PrepareGravitatingMassField(Grids[grid], MetaData, LevelArray,
                                    level) == FAIL)
#endif
          {
          fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
          return FAIL;
        }
    }
  }

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif


  /* Copy overlapping mass fields to ensure consistency and B.C.'s. */

  CommunicationDirection = COMMUNICATION_SEND;
  for (grid = 0; grid < NumberOfGrids; grid++)
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
    for (grid2 = 0; grid2 < SiblingList[grid].NumberOfSiblings; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(SiblingList[grid].GridList[grid2],
#else

    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
#endif
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyOverlappingMassField) == FAIL) {
        fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
        return FAIL;
      }

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  CommunicationDirection = COMMUNICATION_RECEIVE;
  for (grid = 0; grid < NumberOfGrids; grid++)
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
    for (grid2 = 0; grid2 < SiblingList[grid].NumberOfSiblings; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(SiblingList[grid].GridList[grid2],
#else

    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
#endif
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyOverlappingMassField) == FAIL) {
        fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
        return FAIL;
      }

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  /* Compute the potential for the top grid. */

  if (level == 0)
    if (ComputePotentialFieldLevelZero(MetaData, Grids, NumberOfGrids)
        == FAIL) {
      fprintf(stderr, "Error in ComputePotentialFieldLevelZero.\n");
      return FAIL;
    }

  /* Compute a first iteration of the potential and share BV's. */

#define ITERATE_POTENTIAL
#ifdef ITERATE_POTENTIAL
      if (level > 0) {
        CopyPotentialFieldAverage = 1;
        for (int iterate = 0; iterate < 1; iterate++) {

          if (iterate > 0)
            CopyPotentialFieldAverage = 2;

          int Dummy, grid2;
          for (grid = 0; grid < NumberOfGrids; grid++)
            if (Grids[grid]->GridData->SolveForPotential(Dummy, level,
                                                         EvaluateTime) 
                == FAIL) {
              fprintf(stderr, "Error in grid->SolveForPotential.\n");
              return FAIL;
            }

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
          CommunicationDirection = COMMUNICATION_SEND;
          for (grid = 0; grid < NumberOfGrids; grid++)
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
    for (grid2 = 0; grid2 < SiblingList[grid].NumberOfSiblings; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(SiblingList[grid].GridList[grid2],
#else

    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
#endif
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyPotentialField) == FAIL) {
               fprintf(stderr, "Error in grid->CopyPotentialField.\n");
               return FAIL;
             }

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
          CommunicationDirection = COMMUNICATION_RECEIVE;
          for (grid = 0; grid < NumberOfGrids; grid++)
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
    for (grid2 = 0; grid2 < SiblingList[grid].NumberOfSiblings; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(SiblingList[grid].GridList[grid2],
#else

    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
#endif
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyPotentialField) == FAIL) {
               fprintf(stderr, "Error in grid->CopyPotentialField.\n");
               return FAIL;
             }
#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
          CommunicationDirection = COMMUNICATION_SEND_RECEIVE;


        }
        CopyPotentialFieldAverage = 0;
      }
#endif /* ITERATE_POTENTIAL */

  /* if level > MaximumGravityRefinementLevel, then do final potential
     solve (and acceleration interpolation) here rather than in the main
     EvolveLevel since it involves communications. */

  if (reallevel > MaximumGravityRefinementLevel) {

    /* compute potential and acceleration on coarser level [LOCAL]
       (but only if there is at least a subgrid -- it should be only
        if there is a subgrrid on reallevel, but this is ok). */

    for (grid = 0; grid < NumberOfGrids; grid++) 
      if (Grids[grid]->NextGridNextLevel != NULL) {
        Grids[grid]->GridData->SolveForPotential(level,
                                               MaximumGravityRefinementLevel);
        Grids[grid]->GridData->ComputeAccelerationField(
                           (HydroMethod == Zeus_Hydro) ? ZEUS_GRIDS : GRIDS,
                                               MaximumGravityRefinementLevel);
      }

    /* Interpolate potential for reallevel grids from coarser grids. */

    int Dummy;
    LevelHierarchyEntry *Temp = LevelArray[reallevel];

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    CommunicationDirection = COMMUNICATION_SEND;
    while (Temp != NULL) {
      HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
      for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
        Temp3 = Temp3->ParentGrid;
      if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) == FAIL) {
        fprintf(stderr, "Error in grid->InterpolateAccelerations.\n");
        return FAIL;
      }
      Temp = Temp->NextGridThisLevel;
    }

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    CommunicationDirection = COMMUNICATION_RECEIVE;
    Temp = LevelArray[reallevel];
    while (Temp != NULL) {
      HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
      for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
        Temp3 = Temp3->ParentGrid;
      if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) == FAIL) {
        fprintf(stderr, "Error in grid->InterpolateAccelerations.\n");
        return FAIL;
      }
      Temp = Temp->NextGridThisLevel;
    }

#ifdef JB_OPT_BARRIER
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  } // end: if (reallevel > MaximumGravityRefinementLevel)

  return SUCCESS;
}


/* ======================================================================= */
/* This routines does the flux correction and project for all grids on this
   level from the list of subgrids. */

int UpdateFromFinerGrids(HierarchyEntry *Grids[], LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData,
			 int NumberOfGrids,
			 int NumberOfSubgrids[], 
			 fluxes **SubgridFluxesEstimate[]
#ifdef JB_OPT_FLUXES_FIX
			 , LevelHierarchyEntry* SUBlingList[]
#endif
			 ,int Cycle //dcc this is for debugging only.
)

{
  int grid, subgrid;
  HierarchyEntry *NextGrid;
  LevelHierarchyEntry *NextSubgrid;


#ifdef JB_OPT_FLUXES_FIX
  int SUBlingGrid;
  LevelHierarchyEntry *NextSUBling;
#endif


  /* Define a temporary flux holder for the refined fluxes. */

  fluxes SubgridFluxesRefined;

  /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */

  /* -------------- FIRST PASS ----------------- */

  CommunicationDirection = COMMUNICATION_SEND;

  for (grid = 0; grid < NumberOfGrids; grid++) {

    /* Loop over subgrids for this grid. */

    NextGrid = Grids[grid]->NextGridNextLevel;
    subgrid = 0;


    while (NextGrid != NULL && FluxCorrection) {

      /* Project subgrid's refined fluxes to the level of this grid. */

      Grids[0]->GridData->TotalMass("0 GPBF");
      if (NextGrid->GridData->GetProjectedBoundaryFluxes(
		      Grids[grid]->GridData, SubgridFluxesRefined) == FAIL) {
	fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	return FAIL;
      }
      Grids[0]->GridData->TotalMass("1 GPBF");

      NextGrid = NextGrid->NextGridThisLevel;
      subgrid++;
    }
#ifdef JB_OPT_FLUXES_FIX /* Loop over sublings */
    NextSUBling = SUBlingList[grid];
    
    while (NextSUBling != NULL && FluxCorrection) {
      /* make sure this isn't a "proper" subgrid */
      // dcc; note that this doesn't actually do anything, as ParentGrid= NULL 
      // in the Fast SUBling search.
      if( !(NextSUBling->GridHierarchyEntry->ParentGrid == Grids[grid]) ){
	/* Project subgrid's refined fluxes to the level of this grid. */

	Grids[0]->GridData->TotalMass("0 GPBF SUB");
	if (NextSUBling->GridData->
	    GetProjectedBoundaryFluxes( Grids[grid]->GridData, 
					SubgridFluxesRefined ) == FAIL) {
	  fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	  return FAIL;
	}
	Grids[0]->GridData->TotalMass("1 GPBF SUB");
      }
      NextSUBling = NextSUBling->NextGridThisLevel;
    }
#endif /* JB_OPT_FLUXES_FIX */

    
    /* Loop over subgrids for this grid: replace solution. */
    
    NextGrid=Grids[grid]->NextGridNextLevel;
    while( NextGrid != NULL ){

      Grids[0]->GridData->TotalMass("0 ProjSol");
      if (NextGrid->GridData->ProjectSolutionToParentGrid
	  (*Grids[grid]->GridData) == FAIL) {
	fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	return FAIL;
      }
      Grids[0]->GridData->TotalMass("1 ProjSol");      
      NextGrid = NextGrid->NextGridThisLevel;
      
    }

    if(MHD_Used ){
      //dcc 10/15/05 make this SUBlingList[0]
#ifdef  JB_OPT_FLUXES_FIX
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
      NextSubgrid = SUBlingList[grid];
#else  /* Neighbor */
      NextSubgrid=LevelArray[level+1];
#endif /* Neighbor */
#else /* FluxFix */
      NextSubgrid = NULL;
#endif /*FluxFix */
      while( NextSubgrid != NULL ){
	
	if (NextSubgrid->GridData->MHD_ProjectFace
	    (*Grids[grid]->GridData, MetaData->LeftFaceBoundaryCondition,
	     MetaData->RightFaceBoundaryCondition  ) == FAIL) {
	  fprintf(stderr, "Error in grid->MHD_ProjectFace, Send Pass.\n");
	  return FAIL;
	}
	
	NextSubgrid = NextSubgrid->NextGridThisLevel;
      }
    }
   
    Grids[grid]->GridData->CheckForNans("End Of Send Loop");    
  } // end of loop over grids

  /* this output thing.


  */
  /* -------------- SECOND PASS ----------------- */
  CommunicationDirection = COMMUNICATION_RECEIVE;

  dccCounter8++;
  dccCounter5 = -1; // CFR output control, grid number  Initialized to -1 to make incrementing easier.
  for (grid = 0; grid < NumberOfGrids; grid++) {
    dccCounter5++;
    dccCounter6 = -1; // CFR output control, subgrid (subling) number
    //Generate the field that avoids double counting.

    if( MHD_Used == TRUE && MHD_FluxCorrection == TRUE )
      Grids[grid]->GridData->NewElectricFlag();

    NextGrid = Grids[grid]->NextGridNextLevel;
    subgrid = 0;
    
    while (NextGrid != NULL && FluxCorrection) {
      dccCounter6++;
      /* Project subgrid's refined fluxes to the level of this grid. */

      if (NextGrid->GridData->GetProjectedBoundaryFluxes(
			 Grids[grid]->GridData, SubgridFluxesRefined) == FAIL) {
	fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	return FAIL;
      }

      /* Correct this grid for the refined fluxes (step #19)
	 (this also deletes the fields in SubgridFluxesRefined). */
      if(0==0){//this flag should be gone if you're not debugging the flux correction.
	if (Grids[grid]->GridData->CorrectForRefinedFluxes
	    (SubgridFluxesEstimate[grid][subgrid], &SubgridFluxesRefined, 
	     SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1]     
#ifdef JB_OPT_FLUXES_FIX /* Pass extra information */
	     , FALSE, MetaData
#endif
	     ) == FAIL) {
	  fprintf(stderr, "Error in grid->CorrectForRefinedFluxes.\n");
	  return FAIL;
	}
      }else{
	fprintf(stderr,"kludge: no CFR for Subgrids while debugging SUBling\n");
      }
      NextGrid = NextGrid->NextGridThisLevel;
      subgrid++;

    }

#ifdef JB_OPT_FLUXES_FIX /* Do flux correction from Sublings */
    
    NextSUBling = SUBlingList[grid];
    dccCounter6=-1;
    while (NextSUBling != NULL && FluxCorrection) {

      dccCounter6++;

      /* make sure this isn't a "proper" subgrid */
      if( NextSUBling->GridHierarchyEntry->ParentGrid != Grids[grid] ){
	
	/* Project subgrid's refined fluxes to the level of this grid. */
	if (NextSUBling->GridData->
	    GetProjectedBoundaryFluxes( Grids[grid]->GridData, 
					SubgridFluxesRefined ) == FAIL) {
	  fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	  return FAIL;
	}

	/* Correct this grid for the refined fluxes (step #19)
	   (this also deletes the fields in SubgridFluxesRefined). */

	if (Grids[grid]->GridData->CorrectForRefinedFluxes
	    (SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1], 
	     &SubgridFluxesRefined, 
	     SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1],
	     TRUE, MetaData) == FAIL) {
	  fprintf(stderr, "Error in grid->CorrectForRefinedFluxes.\n");
	  return FAIL;
	}


      }

      NextSUBling = NextSUBling->NextGridThisLevel;
    }
    
#endif /* JB_OPT_FLUXES_FIX */

  
    NextGrid = Grids[grid]->NextGridNextLevel;
    dccCounter6 = -1; // CFR output control, subgrid (subling) number

    while (NextGrid != NULL) {
      dccCounter6++;
      /* Project the subgrid solution into this grid. */

      if (NextGrid->GridData->ProjectSolutionToParentGrid
	  (*Grids[grid]->GridData) == FAIL) {
        fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
        return FAIL;
      }

      NextGrid = NextGrid->NextGridThisLevel;
    }
  
    if(MHD_Used){
      //dcc 10/15/5 make this SUBlingList.

#ifdef  JB_OPT_FLUXES_FIX
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
      NextSubgrid = SUBlingList[grid];
#else  /* Neighbor */
      NextSubgrid=LevelArray[level+1];
#endif /* Neighbor */
#else /* FluxFix */
      NextSubgrid = NULL;
#endif /*FluxFix */


      while( NextSubgrid != NULL ){
	
	if (NextSubgrid->GridData->MHD_ProjectFace
	    (*Grids[grid]->GridData, MetaData->LeftFaceBoundaryCondition,
	     MetaData->RightFaceBoundaryCondition  ) == FAIL) {
	  fprintf(stderr, "Error in grid->MHD_ProjectFace, Receive Pass.\n");
	  return FAIL;
	}
	
	NextSubgrid = NextSubgrid->NextGridThisLevel;
      }
    }//mhd_used
    
    if( MHD_Used == TRUE && MHD_FluxCorrection == TRUE ){
      Grids[grid]->GridData->DeleteElectricFlag();

      //not sure why this is here any more.  dcc 11/1/5
      Grids[grid]->GridData->MHDCleanUpTemp();
    }

  } // end of loop over grids

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  return SUCCESS;
}



/* ======================================================================= */
/* This routine simply converts a linked list of grids into an array of
   pointers. */

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[])
{

  /* Count the number of grids on this level. */

  int NumberOfGrids = 0, counter = 0;
  LevelHierarchyEntry *Temp = LevelArray[level];
  while (Temp != NULL) {
    NumberOfGrids++;
    Temp             = Temp->NextGridThisLevel;
  }

  /* Create a list of pointers and number of subgrids (and fill it out). */

  *Grids = new HierarchyEntry *[NumberOfGrids];
  Temp = LevelArray[level];
  while (Temp != NULL) {
    (*Grids)[counter++] = Temp->GridHierarchyEntry;
    Temp              = Temp->NextGridThisLevel;
  }

  return NumberOfGrids;
}
  
