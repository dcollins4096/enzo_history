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
/  INITIALIZE A NEW SIMULATION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include <string.h>
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
#include "StarParticleData.h"

/* function prototypes */

int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt);
int WriteParameterFile(FILE *fptr, TopGridData &MetaData);
void ConvertTotalEnergyToGasEnergy(HierarchyEntry *Grid);
int SetDefaultGlobalValues(TopGridData &MetaData);
int CommunicationPartitionGrid(HierarchyEntry *Grid);

/* Initialization function prototypes. */

int ShockTubeInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid);
int WavePoolInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData);
int ShockPoolInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			TopGridData &MetaData);
int ImplosionInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			TopGridData &MetaData);
int SedovBlastInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry
			 &TopGrid, TopGridData &MetaData);
int DoubleMachInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			TopGridData &MetaData, ExternalBoundary &Exterior);
int ShockInABoxInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData, ExternalBoundary &Exterior);
int ZeldovichPancakeInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid);
int PressurelessCollapseInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid, TopGridData &MetaData);
int AdiabaticExpansionInitialize(FILE *fptr, FILE *Outfptr, 
				 HierarchyEntry &TopGrid);
int TestGravityInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData);
int TestGravitySphereInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid, TopGridData &MetaData);
int SphericalInfallInitialize(FILE *fptr, FILE *Outfptr, 
			      HierarchyEntry &TopGrid, TopGridData &MetaData);
int GravityEquilibriumTestInitialize(FILE *fptr, FILE *Outfptr, 
			      HierarchyEntry &TopGrid, TopGridData &MetaData);
int CollapseTestInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData);
int TestGravityMotion(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData);
int CosmologySimulationInitialize(FILE *fptr, FILE *Outfptr, 
				  HierarchyEntry &TopGrid,
				  TopGridData &MetaData);
int CosmologySimulationReInitialize(HierarchyEntry *TopGrid, 
				    TopGridData &MetaData);
int TurbulenceSimulationInitialize(FILE *fptr, FILE *Outfptr, 
                                  HierarchyEntry &TopGrid,
                                  TopGridData &MetaData);
int ProtostellarCollapseInitialize(FILE *fptr, FILE *Outfptr,
                                  HierarchyEntry &TopGrid,
                                  TopGridData &MetaData);
int TurbulenceSimulationReInitialize(HierarchyEntry *TopGrid, 
                                    TopGridData &MetaData);
int SupernovaRestartInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid, TopGridData &MetaData,
			       ExternalBoundary &Exterior);
//dcc
int MHDBlastInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData, ExternalBoundary &Exterior);
int MHDTestInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData, ExternalBoundary &Exterior);
int MHDShockInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData, ExternalBoundary &Exterior);
int MHDOrszagTangInit(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		      TopGridData &MetaData, ExternalBoundary &Exterior);
int MHDLoopInit(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		TopGridData &MetaData, ExternalBoundary &Exterior);

int MHDCollapseTestInitialize(FILE *fptr, FILE *Outfptr,
                          HierarchyEntry &TopGrid, TopGridData &MetaData);

int DiskInitialize(FILE *fptr, FILE *Outfptr,
		   HierarchyEntry &TopGrid, TopGridData &MetaData);

#ifdef HAOXU
int MHDZeldovichPancakeInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid);
int MHDSphericalInfallInitialize(FILE *fptr, FILE *Outfptr, 
		      HierarchyEntry &TopGrid, TopGridData &MetaData);
int MHDCausticInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid);
int MHDAdiabaticExpansionInitialize(FILE *fptr, FILE *Outfptr,
                                 HierarchyEntry &TopGrid);
int MHDCosmologySimulationInitialize(FILE *fptr, FILE *Outfptr,
                               HierarchyEntry &TopGrid, TopGridData &MetaData);

int MHDCosmologySimulationReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);

int BiermannBatteryInitialize(FILE *fptr, FILE *Outfptr,
                               HierarchyEntry &TopGrid, TopGridData &MetaData);

int BiermannBatteryReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);

#endif /* HAOXU */

/* a test for my code browser
class test{
int fucko();
}
*/


/* character strings */

char outfilename[] = "amr.out";


int InitializeNew(char *filename, HierarchyEntry &TopGrid, 
		  TopGridData &MetaData, ExternalBoundary &Exterior,
		  float *Initialdt)
{

  /* declarations */

  FILE *fptr, *BCfptr, *Outfptr;
  float Dummy[MAX_DIMENSION];
  int dim, i;


  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Dummy[dim] = 0.0;

  /* Open parameter file */

  if ((fptr = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "Error opening parameter file.\n");
    return FAIL;
  }

  /* Open output file */

  if (MyProcessorNumber == ROOT_PROCESSOR)
    if ((Outfptr = fopen(outfilename, "w")) == NULL) {
      fprintf(stderr, "Error opening parameter output file %s\n", outfilename);
      return FAIL;
    }

  /* set the default MetaData values. */

  SetDefaultGlobalValues(MetaData);

  /* Read the MetaData/global values from the Parameter file. */

  if (ReadParameterFile(fptr, MetaData, Initialdt) == FAIL) {
    fprintf(stderr, "Error in ReadParameterFile.\n");
    return FAIL;
  }
  fprintf(stderr,"Left: %f %f %f\n", DomainLeftEdge[0], DomainLeftEdge[1], DomainLeftEdge[2]);
  fprintf(stderr,"Right: %f %f %f\n", DomainRightEdge[0], DomainRightEdge[1], DomainRightEdge[2]);
  /* Set the number of particle attributes, if left unset. */

  if (NumberOfParticleAttributes == INT_UNDEFINED)
    if (StarParticleCreation || StarParticleFeedback)
      NumberOfParticleAttributes = 3;
    else
      NumberOfParticleAttributes = 0;

  /* Give unset parameters their default values. */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    if (RefineRegionLeftEdge[dim] == FLOAT_UNDEFINED)
      RefineRegionLeftEdge[dim]   = DomainLeftEdge[dim];
    if (RefineRegionRightEdge[dim] == FLOAT_UNDEFINED)
      RefineRegionRightEdge[dim]  = DomainRightEdge[dim];
    if (MetaData.MovieRegionLeftEdge[dim] == FLOAT_UNDEFINED)
      MetaData.MovieRegionLeftEdge[dim]   = DomainLeftEdge[dim];
    if (MetaData.MovieRegionRightEdge[dim] == FLOAT_UNDEFINED)
      MetaData.MovieRegionRightEdge[dim]  = DomainRightEdge[dim];
  }

#ifdef ATHENA
  //Some error handling.
  if( (EquationOfState == 1) ){
    if( (ProblemType != 100) &&
	(ProblemType != 103) &&
	(ProblemType != 104) &&
	(ProblemType != 60)  &&
	(ProblemType != 123) ){
      fprintf(stderr,"\n Catestrophic Error!!!!!!\n");
      fprintf(stderr," You've selected Isothermal Equation of State and Problem type %d,\n", ProblemType );
      fprintf(stderr," which hasn't been retro-fitted for the Isothermal Solver. \n");
      fprintf(stderr," See the MHD Documetation on notes for the retrofit.\n\n");
      return FAIL;
    }
    if(  ( HydroMethod != Athena ) &&
	 ( HydroMethod != MHD_Li ) ){
      fprintf(stderr, " Catestrophic Error!  Abandon all hope for the future!\n");
      fprintf(stderr, " Only the Athena HydroMethod (4) and Li (6) have been implemented with EquationOfState == 1\n");
      return FAIL;
    }
    

  }
#endif //ATHENA


  /* If the problem reads in a restart dump, then skip over the following. */

  if (ProblemType != 40) {

  /* Error check the rank. */

  if (MetaData.TopGridRank < 0 || MetaData.TopGridRank > 3) {
    fprintf(stderr, "TopGridRank = %d ill defined.\n", MetaData.TopGridRank);
    return FAIL;
  }

  /* Error check the dimensions and at the same time add ghost zones. */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    if (MetaData.TopGridDims[dim] < 1 ) {
      fprintf(stderr, "Fatal Error: TopGridDims[%d] = %d <1\n", dim, 
	      MetaData.TopGridDims[dim]);
      return FAIL;
    }
    if (MetaData.TopGridDims[dim] > 2048) {
      fprintf(stderr, "Fatal Error: TopGridDims[%d] = %d > 2048.  Come back when you get it right..\n", dim, 
	      MetaData.TopGridDims[dim]);
      return FAIL;
    }
    MetaData.TopGridDims[dim] = (MetaData.TopGridDims[dim] > 1) ?
                     MetaData.TopGridDims[dim] + 2*DEFAULT_GHOST_ZONES : 1;
  }

  /* Create the top grid, prepare it, set the time and parameters. */

  TopGrid.GridData = new grid;
  TopGrid.GridData->PrepareGrid(MetaData.TopGridRank, MetaData.TopGridDims, 
				DomainLeftEdge, DomainRightEdge,
				MetaData.NumberOfParticles);
  TopGrid.GridData->SetTime(MetaData.Time);
  TopGrid.GridData->SetHydroParameters(MetaData.CourantSafetyNumber,
				       MetaData.PPMFlatteningParameter,
				       MetaData.PPMDiffusionParameter,
				       MetaData.PPMSteepeningParameter);
  TopGrid.GridData->SetGravityParameters(MetaData.GravityBoundary);

  /* Repair TopGridDims (subtract ghost zones added earlier). */

  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    MetaData.TopGridDims[dim] = max(MetaData.TopGridDims[dim] - 
				    2*DEFAULT_GHOST_ZONES, 1);

  /* Set TopGrid Hierarchy Entry */

  TopGrid.NextGridThisLevel = NULL;  // always true
  TopGrid.ParentGrid        = NULL;  // always true
  TopGrid.NextGridNextLevel = NULL;  // can be reset by initializer

  } // end: if (ProblemType != 40)

  /* Call problem initializer */

  if (ProblemType == 0) {
    fprintf(stderr, "No problem specified.\n");
    return FAIL;
  }
  int ret = INT_UNDEFINED;

  if (debug)
  printf("InitializeNew: Starting problem initialization.\n");

  ////////////////////////////////////////////////////////
  /* 1) Shocktube problem */

  if (ProblemType == 1)
    ret = ShockTubeInitialize(fptr, Outfptr, TopGrid);

  /* 2) Wave pool */

  if (ProblemType == 2)
    ret = WavePoolInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 3) Shock pool */

  if (ProblemType == 3)
    ret = ShockPoolInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 4) Double Mach reflection */

  if (ProblemType == 4)
    ret = DoubleMachInitialize(fptr, Outfptr, TopGrid, MetaData, Exterior);

  /* 5) ShockInABox */

  if (ProblemType == 5)
    ret = ShockInABoxInitialize(fptr, Outfptr, TopGrid, MetaData, Exterior);

  /* 6) Implosion */

  if (ProblemType == 6)
    ret = ImplosionInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 7) SedovBlast */

  if (ProblemType == 7)
    ret = SedovBlastInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 20) Zeldovich Pancake */

  if (ProblemType == 20)
    ret = ZeldovichPancakeInitialize(fptr, Outfptr, TopGrid);

  /* 21) 1D Pressureless collapse. */

  if (ProblemType == 21)
    ret = PressurelessCollapseInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 22) Adiabatic expansion. */

  if (ProblemType == 22)
    ret = AdiabaticExpansionInitialize(fptr, Outfptr, TopGrid);

  /* 23) GravityTest. */

  if (ProblemType == 23)
    ret = TestGravityInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 24) Spherical Infall. */

  if (ProblemType == 24)
    ret = SphericalInfallInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 25) TestGravitySphere. */

  if (ProblemType == 25)
    ret = TestGravitySphereInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 26) GravityEquilibriumTest. */

  if (ProblemType == 26)
    ret = GravityEquilibriumTestInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 27) CollapseTest. */

  if (ProblemType == 27)
    ret = CollapseTestInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 28) TestGravityMotion. */

  if (ProblemType == 28)
    ret = TestGravityMotion(fptr, Outfptr, TopGrid, MetaData);

  /* 30) Cosmology Simulation. */

  if (ProblemType == 30)
    ret = CosmologySimulationInitialize(fptr, Outfptr, TopGrid, MetaData);

  /* 40) Supernova Explosion from restart. */
  
  if (ProblemType == 40)
    ret = SupernovaRestartInitialize(fptr, Outfptr, TopGrid, MetaData, 
				     Exterior);
  /* 60) Turbulence Simulation. */

  if (ProblemType == 60)
    ret = TurbulenceSimulationInitialize(fptr, Outfptr, TopGrid, MetaData);

  if (ProblemType == 61)
    ret = ProtostellarCollapseInitialize(fptr, Outfptr, TopGrid, MetaData);  

  /* 100) MHD Blast Wave.  Not working */
  if (ProblemType == 100)
    ret = MHDBlastInitialize(fptr, Outfptr, TopGrid, MetaData, Exterior);

  /* 101) MHD Test routines. No physics. */
  if (ProblemType == 101)
    ret = MHDTestInitialize(fptr, Outfptr, TopGrid, MetaData, Exterior);

  if (ProblemType == 102) //MHDShock
    ret = MHDShockInitialize(fptr, Outfptr, TopGrid, MetaData, Exterior);

  if (ProblemType == 103) //This doesn't actually need all those arguments
    ret = MHDOrszagTangInit(fptr, Outfptr, TopGrid, MetaData, Exterior);


  if (ProblemType ==104)
    ret = MHDLoopInit(fptr, Outfptr, TopGrid, MetaData, Exterior);
  
  // Insert New problem intializer here.
  ////////////////////////////////////////////////////////

  if(ProblemType == 123) //MHD Collapse Test
    ret = MHDCollapseTestInitialize(fptr, Outfptr, TopGrid, MetaData);

#ifdef HAOXU
  if(ProblemType == 120) //MHD Zeldovich Pancake 
    ret =  MHDZeldovichPancakeInitialize(fptr, Outfptr, TopGrid);

 if(ProblemType == 121) //MHD Spherical Infall 
    ret =  MHDSphericalInfallInitialize(fptr, Outfptr, TopGrid, MetaData);

  if(ProblemType == 122) //MHD Caustic
   ret = MHDCausticInitialize(fptr, Outfptr, TopGrid);


  if(ProblemType == 124) //MHD Adiabatic Expansion
   ret = MHDAdiabaticExpansionInitialize(fptr,Outfptr,TopGrid);
 
  if(ProblemType == 130)  //MHD Cosmology Simulation
   ret = MHDCosmologySimulationInitialize(fptr,Outfptr,TopGrid, MetaData);

  if(ProblemType == 131)  //Biermann Battery
   ret = BiermannBatteryInitialize(fptr, Outfptr,
                              TopGrid, MetaData);

#endif /* HAOXU */

  if(ProblemType == 132)  //Disk problems.
   ret = DiskInitialize(fptr, Outfptr,
			TopGrid, MetaData);
  
  if (ret == INT_UNDEFINED) {
    fprintf(stderr, "Shit. Problem Type %d undefined.\n", ProblemType);
    return FAIL;
  }

  if (ret == FAIL) {
    fprintf(stderr, "Error in problem initialization.\n");
    return FAIL;
  }

  if (debug)
    printf("InitializeNew: Finished problem initialization.\n");

  /* Do some error checking. */

  if (MetaData.StopTime == FLOAT_UNDEFINED) {
    fprintf(stderr, "StopTime never set.\n");
    return FAIL;
  }

  /* Initialize the exterior (unless it was set in the problem initializer) */

  if (Exterior.AmIPrepared() == FALSE) {

    Exterior.Prepare(TopGrid.GridData);   // set rank and dims

    if (MetaData.BoundaryConditionName != NULL) {

      if ((BCfptr = fopen(MetaData.BoundaryConditionName, "r")) == NULL) {
	fprintf(stderr, "Error opening BC file: %s\n", 
		MetaData.BoundaryConditionName);
	return FAIL;
      }
      if (Exterior.ReadExternalBoundary(BCfptr) == FAIL) {
	fprintf(stderr, "Error in ReadExternalBoundary.\n");
	return FAIL;
      }
      fclose(BCfptr);
    }

    else {

#ifdef OOC_BOUNDARY
      SimpleConstantBoundary = TRUE;

      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
        if (MetaData.LeftFaceBoundaryCondition[dim] != periodic ||
            MetaData.RightFaceBoundaryCondition[dim] != periodic) {
          SimpleConstantBoundary = FALSE;
        }
      }

      if (SimpleConstantBoundary) {
        fprintf(stderr, "SimpleConstantBoundary TRUE\n");
      } else {
        fprintf(stderr, "SimpleConstantBoundary FALSE\n");
      }
#endif

      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	if (Exterior.InitializeExternalBoundaryFace(dim, 
				    MetaData.LeftFaceBoundaryCondition[dim],
				    MetaData.RightFaceBoundaryCondition[dim],
				    Dummy, Dummy)
	    == FAIL) {
	  fprintf(stderr, "Error in InitializeExternalBoundaryFace.\n");
	  return FAIL;
	}

      /* Initialize particle boundary conditions. */

      Exterior.InitializeExternalBoundaryParticles(
					  MetaData.ParticleBoundaryType);

    }  // end: if (MetaData.BoundaryConditionName != NULL)

  }  // end of set Exterior

  /* Set values that were left undefined (above). */

  if (MetaData.TimeLastDataDump == FLOAT_UNDEFINED)
    MetaData.TimeLastDataDump = MetaData.Time - MetaData.dtDataDump*1.00001;
  if (MetaData.TimeLastHistoryDump == FLOAT_UNDEFINED)
    MetaData.TimeLastHistoryDump = MetaData.Time - MetaData.dtHistoryDump;
  if (MetaData.TimeLastMovieDump == FLOAT_UNDEFINED)
    MetaData.TimeLastMovieDump = MetaData.Time - MetaData.dtMovieDump;

  if (MetaData.CycleLastDataDump == INT_UNDEFINED)
    MetaData.CycleLastDataDump = MetaData.CycleNumber - 
                                 MetaData.CycleSkipDataDump;
  if (MetaData.CycleLastHistoryDump == INT_UNDEFINED)
    MetaData.CycleLastHistoryDump = MetaData.CycleNumber - 
                                    MetaData.CycleSkipHistoryDump;

  if (MetaData.TimeLastAnalysis == FLOAT_UNDEFINED)
    MetaData.TimeLastAnalysis = MetaData.Time - MetaData.dtAnalysis*1.00001;

  /* Make changes required for Zeus solver, and turn the TotalEnergy
     variable (should be renamed just Energy) into GasEnergy. */

  if (HydroMethod == Zeus_Hydro && ProblemType != 20 && ProblemType != 30 &&
      ProblemType != 27)
    ConvertTotalEnergyToGasEnergy(&TopGrid);

  /* If using StarParticles, set the number to zero. */

  if (StarParticleCreation || StarParticleFeedback)
    NumberOfStarParticles = 0;

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  for (i = 0; i < MAX_FLAGGING_METHODS; i++)
    if (MinimumMassForRefinement[i] == FLOAT_UNDEFINED) {
      MinimumMassForRefinement[i] = MinimumOverDensityForRefinement[i];
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	MinimumMassForRefinement[i] *=
	  (DomainRightEdge[dim]-DomainLeftEdge[dim])/
	  float(MetaData.TopGridDims[dim]);
    }

  /* Write the MetaData/global values to the Parameter file. */

  if (MyProcessorNumber == ROOT_PROCESSOR)
    if (WriteParameterFile(Outfptr, MetaData) == FAIL) {
      fprintf(stderr, "Error in WriteParameterFile.\n");
      return FAIL;
    }

  /* Close parameter files */

  fclose(fptr);
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fclose(Outfptr);

  /* Partition the top grid. */
  
  if( CommunicationPartitionGrid(&TopGrid) == FAIL ){
    fprintf(stderr," Error in CommunicationPartitioGrid.\n");
    return FAIL;
  }
  /* For problem 30, using ParallelGridIO, read in data only after
     partitioning grid. */

  if (ParallelRootGridIO == TRUE && ProblemType == 30)
    if (CosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL) {
      fprintf(stderr, "Error in CosmologySimulationReInitialize.\n");
      return FAIL;
    }

   /* For problem 60, using ParallelGridIO, read in data only after
     partitioning grid. */
     
  if (ParallelRootGridIO == TRUE && ProblemType == 60)
    if (TurbulenceSimulationReInitialize(&TopGrid, MetaData) == FAIL) {
      fprintf(stderr, "Error in TurbulenceSimulationReInitialize.\n");
      return FAIL;
    }

#ifdef HAOXU
   /* For problem 130 and 131, using ParallelGridIO, read in data only after
     partitioning grid. */

  if(ParallelRootGridIO == TRUE && ProblemType == 130)
     if (MHDCosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL) {
      fprintf(stderr, "Error in MHDCosmologySimulationReInitialize.\n");
      return FAIL;
    }             

  if(ParallelRootGridIO == TRUE && ProblemType == 131)
     if (BiermannBatteryReInitialize(&TopGrid, MetaData) == FAIL) {
      fprintf(stderr, "Error in BiermanBatteryReInitialize.\n");
      return FAIL;
    }
#endif  

  fprintf(stderr, "End of InitializeNew.C\n");

  return SUCCESS;
}


void ConvertTotalEnergyToGasEnergy(HierarchyEntry *Grid)
{
  if (Grid != NULL) {
    Grid->GridData->ConvertTotalEnergyToGasEnergy();
    ConvertTotalEnergyToGasEnergy(Grid->NextGridThisLevel);
    ConvertTotalEnergyToGasEnergy(Grid->NextGridNextLevel);
  }
}
