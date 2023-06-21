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
/  INITIALIZE A COSMOLOGY SIMULATION
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:  Initialize for cosmology simulations.  Reads in a number
/      of initial grids.  If more than one, then they are all numbered.
/      We currently assume that all are subgrids except the first.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include <string.h>
#include <stdio.h>
#include <math.h>
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "fortran.def"
#include "error.h"

/* function prototypes */

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
int CommunicationAllSumIntegerValues(int *Values, int Number);

/* Cosmology Parameters (that need to be shared) */

static float CosmologySimulationOmegaBaryonNow       = 1.0;  // standard
static float CosmologySimulationOmegaCDMNow          = 0.0;  // no dark matter
static float CosmologySimulationInitialTemperature   = FLOAT_UNDEFINED;

static char *CosmologySimulationDensityName          = NULL;
static char *CosmologySimulationTotalEnergyName      = NULL;
static char *CosmologySimulationGasEnergyName        = NULL;
static char *CosmologySimulationParticlePositionName = NULL;
static char *CosmologySimulationParticleVelocityName = NULL;
static char *CosmologySimulationParticleMassName     = NULL;
static char *CosmologySimulationVelocityNames[MAX_DIMENSION];

static int   CosmologySimulationSubgridsAreStatic    = TRUE;
static int   CosmologySimulationNumberOfInitialGrids = 1;

static float CosmologySimulationInitialFractionHII   = 1.2e-5;
static float CosmologySimulationInitialFractionHeII  = 1.0e-14;
static float CosmologySimulationInitialFractionHeIII = 1.0e-17;
static float CosmologySimulationInitialFractionHM    = 2.0e-9;
static float CosmologySimulationInitialFractionH2I   = 2.0e-20;
static float CosmologySimulationInitialFractionH2II  = 3.0e-14;
static int   CosmologySimulationUseMetallicityField  = FALSE;

#define MAX_INITIAL_GRIDS 10

int CosmologySimulationInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  //
  // Initialize variables.
  //

  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";
  char *MetalName = "Metal_Density";

  char *ExtraNames[2] = {"Z_Field1", "Z_Field2"};

  /* declarations */

  char line[MAX_LINE_LENGTH];
  int i, j, dim, gridnum, ret, SubgridsAreStatic, region;
  HierarchyEntry *Subgrid;

  char *DensityName = NULL, *TotalEnergyName = NULL, *GasEnergyName = NULL,
       *ParticlePositionName = NULL, *ParticleVelocityName = NULL, 
       *ParticleMassName = NULL, *VelocityNames[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    VelocityNames[dim] = NULL;

  /* Set default parameters: parameters, names and subgrid info */

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    CosmologySimulationVelocityNames[dim]       = NULL;

  int   CosmologySimulationGridDimension[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  int   CosmologySimulationGridLevel[MAX_INITIAL_GRIDS];
  FLOAT CosmologySimulationGridLeftEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  FLOAT CosmologySimulationGridRightEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  for (i = 0; i < MAX_INITIAL_GRIDS; i++)
    CosmologySimulationGridLevel[i] = 1;
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    CosmologySimulationGridLeftEdge[0][dim] = DomainLeftEdge[dim];
    CosmologySimulationGridRightEdge[0][dim] = DomainRightEdge[dim];
    CosmologySimulationGridDimension[0][dim] = MetaData.TopGridDims[dim];
  }
  CosmologySimulationGridLevel[0] = 0;

  /* Error check. */

  if (!ComovingCoordinates) {
    fprintf(stderr, "ComovingCoordinates must be TRUE!\n");
    return FAIL;
  }

  if (DualEnergyFormalism == FALSE && HydroMethod != Zeus_Hydro)
    fprintf(stderr, "CosmologySimulation: DualEnergyFormalism is off!\n");
  if (!SelfGravity)
    fprintf(stderr, "CosmologySimulation: gravity is off!?!\n");


  //
  // Read parameter file. 
  // 


  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* Read parameters */

    ret += sscanf(line, "CosmologySimulationOmegaBaryonNow = %"FSYM, 
		  &CosmologySimulationOmegaBaryonNow);
    ret += sscanf(line, "CosmologySimulationOmegaCDMNow = %"FSYM, 
		  &CosmologySimulationOmegaCDMNow);
    ret += sscanf(line, "CosmologySimulationInitialTemperature = %"FSYM, 
		  &CosmologySimulationInitialTemperature);
    
    if (sscanf(line, "CosmologySimulationDensityName = %s", dummy) == 1)
      CosmologySimulationDensityName = dummy;
    if (sscanf(line, "CosmologySimulationTotalEnergyName = %s", dummy) == 1)
      CosmologySimulationTotalEnergyName = dummy;
    if (sscanf(line, "CosmologySimulationGasEnergyName = %s", dummy) == 1)
      CosmologySimulationGasEnergyName = dummy;
    if (sscanf(line, "CosmologySimulationVelocity1Name = %s", dummy) == 1)
      CosmologySimulationVelocityNames[0] = dummy;
    if (sscanf(line, "CosmologySimulationVelocity2Name = %s", dummy) == 1)
      CosmologySimulationVelocityNames[1] = dummy;
    if (sscanf(line, "CosmologySimulationVelocity3Name = %s", dummy) == 1)
      CosmologySimulationVelocityNames[2] = dummy;
    if (sscanf(line, "CosmologySimulationParticlePositionName = %s", dummy) 
	== 1) CosmologySimulationParticlePositionName = dummy;
    if (sscanf(line, "CosmologySimulationParticleVelocityName = %s", dummy) 
	== 1) CosmologySimulationParticleVelocityName = dummy;
    if (sscanf(line, "CosmologySimulationParticleMassName = %s", dummy) == 1)
      CosmologySimulationParticleMassName = dummy;

    ret += sscanf(line, "CosmologySimulationNumberOfInitialGrids = %d",
		  &CosmologySimulationNumberOfInitialGrids);
    ret += sscanf(line, "CosmologySimulationSubgridsAreStatic = %d",
		  &CosmologySimulationSubgridsAreStatic);

    if (sscanf(line, "CosmologySimulationGridLeftEdge[%d]", &gridnum) > 0)
      ret += sscanf(line, "CosmologySimulationGridLeftEdge[%d] = %"PSYM" %"PSYM" %"PSYM,
		    &gridnum, &CosmologySimulationGridLeftEdge[gridnum][0],
		    &CosmologySimulationGridLeftEdge[gridnum][1],
		    &CosmologySimulationGridLeftEdge[gridnum][2]);
    if (sscanf(line, "CosmologySimulationGridRightEdge[%d]", &gridnum) > 0)
      ret += sscanf(line, "CosmologySimulationGridRightEdge[%d] = %"PSYM" %"PSYM" %"PSYM,
		    &gridnum, &CosmologySimulationGridRightEdge[gridnum][0],
		    &CosmologySimulationGridRightEdge[gridnum][1],
		    &CosmologySimulationGridRightEdge[gridnum][2]);
    if (sscanf(line, "CosmologySimulationGridDimension[%d]", &gridnum) > 0)
      ret += sscanf(line, "CosmologySimulationGridDimension[%d] = %d %d %d",
		    &gridnum, &CosmologySimulationGridDimension[gridnum][0],
		    &CosmologySimulationGridDimension[gridnum][1],
		    &CosmologySimulationGridDimension[gridnum][2]);
    if (sscanf(line, "CosmologySimulationGridLevel[%d]", &gridnum) > 0)
      ret += sscanf(line, "CosmologySimulationGridLevel[%d] = %d",
		    &gridnum, &CosmologySimulationGridLevel[gridnum]);

    ret += sscanf(line, "CosmologySimulationInitialFractionHII = %"FSYM,
		  &CosmologySimulationInitialFractionHII);
    ret += sscanf(line, "CosmologySimulationInitialFractionHeII = %"FSYM,
		  &CosmologySimulationInitialFractionHeII);
    ret += sscanf(line, "CosmologySimulationInitialFractionHeIII = %"FSYM,
		  &CosmologySimulationInitialFractionHeIII);
    ret += sscanf(line, "CosmologySimulationInitialFractionHM = %"FSYM,
		  &CosmologySimulationInitialFractionHM);
    ret += sscanf(line, "CosmologySimulationInitialFractionH2I = %"FSYM,
		  &CosmologySimulationInitialFractionH2I);
    ret += sscanf(line, "CosmologySimulationInitialFractionH2II = %"FSYM,
		  &CosmologySimulationInitialFractionH2II);
    ret += sscanf(line, "CosmologySimulationUseMetallicityField = %d",
		  &CosmologySimulationUseMetallicityField);

    /* If the dummy char space was used, then make another. */

    if (dummy[0] != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "CosmologySimulation") &&
	line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }

  //
  // More error checking. 
  //

  if (CosmologySimulationDensityName == NULL && 
      CosmologySimulationParticlePositionName == NULL) {
    fprintf(stderr, "Missing initial data.\n");
    return FAIL;
  }

  if (CosmologySimulationDensityName != NULL && CellFlaggingMethod[0] != 2)
      fprintf(stderr, "CosmologySimulation: check CellFlaggingMethod.\n");

  if (CosmologySimulationDensityName == NULL && CellFlaggingMethod[0] != 4)
      fprintf(stderr, "CosmologySimulation: check CellFlaggingMethod.\n");

  if (CosmologySimulationNumberOfInitialGrids > MAX_INITIAL_GRIDS) {
    fprintf(stderr, "Too many InitialGrids! increase MAX_INITIAL_GRIDS\n");
    return FAIL;
  }

  if (CosmologySimulationDensityName == NULL && 
      MultiSpecies+RadiativeCooling > 0) {
    fprintf(stderr, "warning: no density field; setting MultiSpecies/RadiativeCooling = 0\n");
    MultiSpecies = RadiativeCooling = 0;
  }


  // If temperature is left unset, set it assuming that T=550 K at z=200.

  if (CosmologySimulationInitialTemperature == FLOAT_UNDEFINED)
    CosmologySimulationInitialTemperature = 550.0 *
      POW((1.0 + InitialRedshift)/(1.0 + 200), 2);

  // ---------------------------------------------
  // Generate the grids and set-up the hierarchy. 
  // ---------------------------------------------

  HierarchyEntry *GridsList[MAX_INITIAL_GRIDS];
  GridsList[0] = &TopGrid;
  for (gridnum = 1; gridnum < CosmologySimulationNumberOfInitialGrids; 
       gridnum++) {

    // Create a spot in the hierarchy. 

    Subgrid    = new HierarchyEntry;

    //
    // Find Parent Grid
    //

    int ParentGrid = INT_UNDEFINED;
    for (i = 0; i < gridnum; i++)
      if (CosmologySimulationGridLevel[i] == 
	  CosmologySimulationGridLevel[gridnum]-1)
	for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	  if (CosmologySimulationGridLeftEdge[gridnum][dim] <
	      CosmologySimulationGridLeftEdge[i][dim]       ||
	      CosmologySimulationGridRightEdge[gridnum][dim] >
	      CosmologySimulationGridRightEdge[i][dim]       )
	    break;
	  ParentGrid = i;
	}

    if (ParentGrid == INT_UNDEFINED) {
      fprintf(stderr, "Grid %d has no valid parent.\n", gridnum);
      return FAIL;
    }

    //
    // Add to linked list
    //

    GridsList[gridnum] = Subgrid;
    Subgrid->NextGridNextLevel = NULL;
    Subgrid->NextGridThisLevel = GridsList[ParentGrid]->NextGridNextLevel;
    Subgrid->ParentGrid        = GridsList[ParentGrid];
    GridsList[ParentGrid]->NextGridNextLevel = Subgrid;

    /* Error check for consistency and add ghost zones to dimension. */

    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      FLOAT SubgridCellSize = (DomainRightEdge[dim] - DomainLeftEdge[dim])/
	FLOAT(MetaData.TopGridDims[dim]*
	      POW(FLOAT(RefineBy), CosmologySimulationGridLevel[gridnum]));
      printf("%g\n", SubgridCellSize);
      printf("%g\n", POW(FLOAT(RefineBy), CosmologySimulationGridLevel[gridnum]));
      printf("%d %d\n", MetaData.TopGridDims[dim], dim);
      printf("%g %g\n", CosmologySimulationGridRightEdge[gridnum][dim],
	     CosmologySimulationGridLeftEdge[gridnum][dim]);
      printf("%d\n", nint((CosmologySimulationGridRightEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[gridnum][dim]   )
			  /SubgridCellSize));

      /* Check if declared size matches left/right edges. */

      if (nint((CosmologySimulationGridRightEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[gridnum][dim]   )
	       /SubgridCellSize) !=
	  CosmologySimulationGridDimension[gridnum][dim]) {
	fprintf(stderr, "Subgrid inconsistency: grid %d, dim %d\n",
		gridnum, dim);
	fprintf(stderr, " subgrid: %"GOUTSYM" -> %"GOUTSYM", CellSize = %"GOUTSYM"\n",
	      CosmologySimulationGridLeftEdge[gridnum][dim],
	      CosmologySimulationGridRightEdge[gridnum][dim], SubgridCellSize);
	return FAIL;
      }

      /* Check if left/right edge fall on Parent cell boundary. */

      if (nint((CosmologySimulationGridLeftEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[ParentGrid][dim])/
	       SubgridCellSize) % RefineBy != 0 ||
	  nint((CosmologySimulationGridRightEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[ParentGrid][dim])/
	       SubgridCellSize) % RefineBy != 0 ) {
	fprintf(stderr, "Subgrid inconsistency: grid %d, dim %d\n",
		gridnum, dim);
	fprintf(stderr, "left or right edges are not on parent cell edge.\n");
	return FAIL;
      }

      /* Add ghost zones. */

      CosmologySimulationGridDimension[gridnum][dim] += 2*DEFAULT_GHOST_ZONES;
    }

    /* Create a new subgrid and initialize it */

    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(Subgrid->ParentGrid->GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, 
				   CosmologySimulationGridDimension[gridnum],
				   CosmologySimulationGridLeftEdge[gridnum],
				   CosmologySimulationGridRightEdge[gridnum],
				   0);

    /* If subgrids are static, convert to static regions. */

    if (CosmologySimulationSubgridsAreStatic == TRUE) {
      for (region = 0; region < MAX_STATIC_REGIONS; region++)
	if (StaticRefineRegionLevel[region] == INT_UNDEFINED) {
	  StaticRefineRegionLevel[region] = 
	    CosmologySimulationGridLevel[gridnum] - 1;
	  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	    StaticRefineRegionLeftEdge[region][dim] = 
	      CosmologySimulationGridLeftEdge[gridnum][dim];
	    StaticRefineRegionRightEdge[region][dim] = 
	      CosmologySimulationGridRightEdge[gridnum][dim];
	  }
	  for (dim = MetaData.TopGridRank; dim < MAX_DIMENSION; dim++) {
	    StaticRefineRegionLeftEdge[region][dim] = DomainLeftEdge[dim];
	    StaticRefineRegionRightEdge[region][dim] = DomainRightEdge[dim];
	  }
	  break;
	}
      if (region == MAX_STATIC_REGIONS) {
	fprintf(stderr, "Out of static refine regions\n");
	return FAIL;
      }
    }

    //
    //Remove ghost zones from dim, in order to write them out.
    //

    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      CosmologySimulationGridDimension[gridnum][dim] -= 
	2*DEFAULT_GHOST_ZONES;

  } // end: loop over gridnums

  /* Initialize the previously-generated grids. */

  for (gridnum = 0; gridnum < CosmologySimulationNumberOfInitialGrids; 
       gridnum++) {

    /* If there is more than one grid, add the grid number to the name. */

    if (CosmologySimulationNumberOfInitialGrids > 1) {

      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("CosmologySimulation: Initializing grid %d\n", gridnum);

      if (CosmologySimulationDensityName) 
	sprintf(DensityName = new char[MAX_LINE_LENGTH], "%s.%1d", 
		CosmologySimulationDensityName, gridnum);
      if (CosmologySimulationTotalEnergyName)
	sprintf(TotalEnergyName = new char[MAX_LINE_LENGTH], "%s.%1d", 
		CosmologySimulationTotalEnergyName, gridnum);
      if (CosmologySimulationGasEnergyName)
	sprintf(GasEnergyName = new char[MAX_LINE_LENGTH], "%s.%1d", 
		CosmologySimulationGasEnergyName, gridnum);
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	if (CosmologySimulationVelocityNames[dim])
	  sprintf(VelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1d", 
		  CosmologySimulationVelocityNames[dim], gridnum);
      if (CosmologySimulationParticlePositionName)
	sprintf(ParticlePositionName = new char[MAX_LINE_LENGTH], "%s.%1d", 
		CosmologySimulationParticlePositionName, gridnum);
      if (CosmologySimulationParticleVelocityName)
	sprintf(ParticleVelocityName = new char[MAX_LINE_LENGTH], "%s.%1d", 
		CosmologySimulationParticleVelocityName, gridnum);
      if (CosmologySimulationParticleMassName)
	sprintf(ParticleMassName = new char[MAX_LINE_LENGTH], "%s.%1d", 
		CosmologySimulationParticleMassName, gridnum);

    } else {
      DensityName            = CosmologySimulationDensityName;
      TotalEnergyName        = CosmologySimulationTotalEnergyName;
      GasEnergyName          = CosmologySimulationGasEnergyName;
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	VelocityNames[dim]   = CosmologySimulationVelocityNames[dim];
      ParticlePositionName   = CosmologySimulationParticlePositionName;
      ParticleVelocityName   = CosmologySimulationParticleVelocityName;
      ParticleMassName       = CosmologySimulationParticleMassName;
    }

    /* If there is a subgrid, use CosmologySimulationSubgridsAreStatic,
       otherwise just set to false. */

    SubgridsAreStatic = (GridsList[gridnum]->NextGridNextLevel == NULL) ?
      FALSE : CosmologySimulationSubgridsAreStatic;

    /* Initialize the grid by reading in data. */

    int TotalRefinement = nint(POW(FLOAT(RefineBy),
				   CosmologySimulationGridLevel[gridnum]));
    if (GridsList[gridnum]->GridData->CosmologySimulationInitializeGrid(
			     CosmologySimulationOmegaBaryonNow,
			       CosmologySimulationOmegaCDMNow,
			       CosmologySimulationInitialTemperature,
			     DensityName, TotalEnergyName,
			       GasEnergyName, VelocityNames,
			       ParticlePositionName, ParticleVelocityName,
			       ParticleMassName,
			     SubgridsAreStatic, TotalRefinement,
			     CosmologySimulationInitialFractionHII,
			     CosmologySimulationInitialFractionHeII,
			     CosmologySimulationInitialFractionHeIII,
			     CosmologySimulationInitialFractionHM,
			     CosmologySimulationInitialFractionH2I,
			     CosmologySimulationInitialFractionH2II,
			     CosmologySimulationUseMetallicityField,
			     MetaData.NumberOfParticles
						       ) == FAIL) {
      fprintf(stderr, "Error in grid->CosmologySimulationInitializeGrid.\n");
      return FAIL;
    }

    /* Set boundary conditions if necessary. */
  } // end loop over initial grids

  /* -------------------------------------------------------------------- */
  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set).
     Note: multiply MinimumMassForRefinement by the OmegaBaryonNow since the
     routine that uses this parameter only counts baryonic mass. */

  for (i = 0; i < MAX_FLAGGING_METHODS; i++) {
    if (MinimumMassForRefinement[i] == FLOAT_UNDEFINED) {

      MinimumMassForRefinement[i] = CosmologySimulationOmegaBaryonNow/
	                            OmegaMatterNow;
      if (CellFlaggingMethod[i] == 4)
	MinimumMassForRefinement[i] = CosmologySimulationOmegaCDMNow/
	                              OmegaMatterNow;

      MinimumMassForRefinement[i] *= MinimumOverDensityForRefinement[i];
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	MinimumMassForRefinement[i] *= 
	  (DomainRightEdge[dim]-DomainLeftEdge[dim])/
	  float(MetaData.TopGridDims[dim]);
    }
  }

  /* set up field names and units */

  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;
  if (MultiSpecies) {
    DataLabel[i++] = ElectronName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = HIIName;
    DataLabel[i++] = HeIName;
    DataLabel[i++] = HeIIName;
    DataLabel[i++] = HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[i++] = HMName;
      DataLabel[i++] = H2IName;
      DataLabel[i++] = H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
    }
  }

  if (CosmologySimulationUseMetallicityField) {
    DataLabel[i++] = MetalName;
    DataLabel[i++] = ExtraNames[0];
    DataLabel[i++] = ExtraNames[1];
  }

  for (j = 0; j < i; j++)
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */
  
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "CosmologySimulationOmegaBaryonNow       = %f\n", 
	    CosmologySimulationOmegaBaryonNow);
    fprintf(Outfptr, "CosmologySimulationOmegaCDMNow          = %f\n", 
	    CosmologySimulationOmegaCDMNow);
    fprintf(Outfptr, "CosmologySimulationInitialTemperature   = %f\n\n", 
	    CosmologySimulationInitialTemperature);

    fprintf(Outfptr, "CosmologySimulationDensityName          = %s\n",
	    CosmologySimulationDensityName);
    if (CosmologySimulationTotalEnergyName)
    fprintf(Outfptr, "CosmologySimulationTotalEnergyName      = %s\n",
	    CosmologySimulationTotalEnergyName);
    if (CosmologySimulationGasEnergyName)
    fprintf(Outfptr, "CosmologySimulationGasEnergyName        = %s\n",
	    CosmologySimulationGasEnergyName);
    fprintf(Outfptr, "CosmologySimulationVelocity1Name        = %s\n",
	    CosmologySimulationVelocityNames[0]);
    fprintf(Outfptr, "CosmologySimulationVelocity2Name        = %s\n",
	    CosmologySimulationVelocityNames[1]);
    fprintf(Outfptr, "CosmologySimulationVelocity3Name        = %s\n",
	    CosmologySimulationVelocityNames[2]);
    if (CosmologySimulationParticlePositionName)
    fprintf(Outfptr, "CosmologySimulationParticlePositionName = %s\n",
	    CosmologySimulationParticlePositionName);
    if (CosmologySimulationParticleVelocityName)
    fprintf(Outfptr, "CosmologySimulationParticleVelocityName = %s\n",
	    CosmologySimulationParticleVelocityName);
    if (CosmologySimulationParticleMassName)
    fprintf(Outfptr, "CosmologySimulationParticleMassName     = %s\n\n",
	    CosmologySimulationParticleMassName);

    fprintf(Outfptr, "CosmologySimulationNumberOfInitialGrids = %d\n",
	    CosmologySimulationNumberOfInitialGrids);
    fprintf(Outfptr, "CosmologySimulationSubgridsAreStatic    = %d\n",
	    CosmologySimulationSubgridsAreStatic);
    for (gridnum = 1; gridnum < CosmologySimulationNumberOfInitialGrids; 
	 gridnum++) {
      fprintf(Outfptr, "CosmologySimulationGridLeftEdge[%d]     = ", gridnum);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CosmologySimulationGridLeftEdge[gridnum]);
      fprintf(Outfptr, "CosmologySimulationGridRightEdge[%d]    = ", gridnum);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CosmologySimulationGridRightEdge[gridnum]);
      fprintf(Outfptr, "CosmologySimulationGridDimension[%d]    = ", gridnum);
      WriteListOfInts(Outfptr, MetaData.TopGridRank,
		      CosmologySimulationGridDimension[gridnum]);
    }
    fprintf(Outfptr, "\n");

    fprintf(Outfptr, "CosmologySimulationInitialFractionHII   = %g\n",
	    CosmologySimulationInitialFractionHII);
    fprintf(Outfptr, "CosmologySimulationInitialFractionHeII  = %g\n",
	    CosmologySimulationInitialFractionHeII);
    fprintf(Outfptr, "CosmologySimulationInitialFractionHeIII = %g\n",
	    CosmologySimulationInitialFractionHeIII);
    fprintf(Outfptr, "CosmologySimulationInitialFractionHM    = %g\n",
	    CosmologySimulationInitialFractionHM);
    fprintf(Outfptr, "CosmologySimulationInitialFractionH2I   = %g\n",
	    CosmologySimulationInitialFractionH2I);
    fprintf(Outfptr, "CosmologySimulationInitialFractionH2II  = %g\n",
	    CosmologySimulationInitialFractionH2II);
    fprintf(Outfptr, "CosmologySimulationUseMetallicityField  = %d\n\n",
	    CosmologySimulationUseMetallicityField);
  }

  /* Clean up. */

  delete dummy;

  return SUCCESS;
}


void RecursivelySetParticleCount(HierarchyEntry *GridPoint, int *Count);


/* -------------------------------------------------------------------- 
   Re-call the initializer on level zero grids.  Used in case of 
   ParallelRootGridIO. */

int CosmologySimulationReInitialize(HierarchyEntry *TopGrid, 
				    TopGridData &MetaData)
{

  /* Declarations. */

  int dim, gridnum = 0;
  char *DensityName = NULL, *TotalEnergyName = NULL, *GasEnergyName = NULL,
       *ParticlePositionName = NULL, *ParticleVelocityName = NULL, 
       *ParticleMassName = NULL, *VelocityNames[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    VelocityNames[dim] = NULL;

  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("CosmologySimulation: ReInitializing grid %d\n", gridnum);

  /* If there is more than one grid, add the grid number to the name. */

  if (CosmologySimulationNumberOfInitialGrids > 1) {

    if (CosmologySimulationDensityName) 
      sprintf(DensityName = new char[MAX_LINE_LENGTH], "%s.%1d", 
	      CosmologySimulationDensityName, gridnum);
    if (CosmologySimulationTotalEnergyName)
      sprintf(TotalEnergyName = new char[MAX_LINE_LENGTH], "%s.%1d", 
	      CosmologySimulationTotalEnergyName, gridnum);
    if (CosmologySimulationGasEnergyName)
      sprintf(GasEnergyName = new char[MAX_LINE_LENGTH], "%s.%1d", 
	      CosmologySimulationGasEnergyName, gridnum);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      if (CosmologySimulationVelocityNames[dim])
	sprintf(VelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1d", 
		CosmologySimulationVelocityNames[dim], gridnum);
    if (CosmologySimulationParticlePositionName)
      sprintf(ParticlePositionName = new char[MAX_LINE_LENGTH], "%s.%1d", 
	      CosmologySimulationParticlePositionName, gridnum);
    if (CosmologySimulationParticleVelocityName)
      sprintf(ParticleVelocityName = new char[MAX_LINE_LENGTH], "%s.%1d", 
	      CosmologySimulationParticleVelocityName, gridnum);
    if (CosmologySimulationParticleMassName)
      sprintf(ParticleMassName = new char[MAX_LINE_LENGTH], "%s.%1d", 
	      CosmologySimulationParticleMassName, gridnum);

  } else {
    DensityName            = CosmologySimulationDensityName;
    TotalEnergyName        = CosmologySimulationTotalEnergyName;
    GasEnergyName          = CosmologySimulationGasEnergyName;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      VelocityNames[dim]   = CosmologySimulationVelocityNames[dim];
    ParticlePositionName   = CosmologySimulationParticlePositionName;
    ParticleVelocityName   = CosmologySimulationParticleVelocityName;
    ParticleMassName       = CosmologySimulationParticleMassName;
  }

  /* If there is a subgrid, use CosmologySimulationSubgridsAreStatic,
     otherwise just set to false. */

  int SubgridsAreStatic = (TopGrid->NextGridNextLevel == NULL) ?
    FALSE : CosmologySimulationSubgridsAreStatic;

  /* Call grid initializer.  Use TotalRefinement = -1 to flag real read. */

  int TotalRefinement = -1;

  /* Loop over level zero grid. */

  HierarchyEntry *Temp = TopGrid;
  while (Temp != NULL) {

    if (Temp->GridData->CosmologySimulationInitializeGrid(
			     CosmologySimulationOmegaBaryonNow,
			       CosmologySimulationOmegaCDMNow,
			       CosmologySimulationInitialTemperature,
			     DensityName, TotalEnergyName,
			       GasEnergyName, VelocityNames,
			       ParticlePositionName, ParticleVelocityName,
			       ParticleMassName,
			     SubgridsAreStatic, TotalRefinement,
			     CosmologySimulationInitialFractionHII,
			     CosmologySimulationInitialFractionHeII,
			     CosmologySimulationInitialFractionHeIII,
			     CosmologySimulationInitialFractionHM,
			     CosmologySimulationInitialFractionH2I,
			     CosmologySimulationInitialFractionH2II,
			     CosmologySimulationUseMetallicityField,
			     MetaData.NumberOfParticles
						       ) == FAIL) {
      fprintf(stderr, "Error in grid->CosmologySimulationInitializeGrid.\n");
      return FAIL;
    }

    Temp = Temp->NextGridThisLevel;
  }

  int LocalNumberOfParticles;
//  int GlobalSum;

  Temp = TopGrid;
  while (Temp != NULL) {

    LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();
//    printf("OldLocalParticleCount: %d\n", LocalNumberOfParticles );

//    CHECK_MPI_ERROR(MPI_Allreduce(&LocalNumberOfParticles, &GlobalSum, 1, 
//                                  MPI_INT, MPI_SUM, MPI_COMM_WORLD));
//    printf("GlobalSum: %d\n", GlobalSum);
//    Temp->GridData->SetNumberOfParticles( GlobalSum );

    CommunicationAllSumIntegerValues(&LocalNumberOfParticles, 1);
    Temp->GridData->SetNumberOfParticles(LocalNumberOfParticles);

    LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();
//    printf("NewLocalParticleCount: %d\n", LocalNumberOfParticles );
    
    Temp = Temp->NextGridThisLevel;
  }



  /* Loop over grids and set particle ID number. */

  Temp = TopGrid;
  int ParticleCount = 0;
  RecursivelySetParticleCount(Temp, &ParticleCount);
  if (debug)
    printf("FinalParticleCount = %d\n", ParticleCount);
  MetaData.NumberOfParticles = 0;

  return SUCCESS;
}

void RecursivelySetParticleCount(HierarchyEntry *GridPoint, int *Count)
{
  /* Add Count to the particle id's on this grid (which start from zero
     since we are doing a parallel root grid i/o). */

  GridPoint->GridData->AddToParticleNumber(Count);

  /* Recursively apply this to siblings and childern. */

  if (GridPoint->NextGridThisLevel != NULL)
    RecursivelySetParticleCount(GridPoint->NextGridThisLevel, Count);

  if (GridPoint->NextGridNextLevel != NULL)
    RecursivelySetParticleCount(GridPoint->NextGridNextLevel, Count);

  return;
}
