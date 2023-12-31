/***********************************************************************
/
/  INITIALIZE A TURBULENCE SIMULATION
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/  modified1:
/
/  PURPOSE:  Initialize a turbulence simulation.  Reads in initial data
/            for the root grid.
/
/  DEFAULT:  Quasi-isothermal forced turbulence.
/            Requires RAREFACTION1 (or 0), defined in euler.src.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
 
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "fortran.def"
#include "error.h"
#include "message.h"
 
/* function prototypes */
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
int CommunicationAllSumIntegerValues(int *Values, int Number);
 
/* Turbulence Parameters (that need to be shared) */
 
static float TurbulenceSimulationInitialDensity       = FLOAT_UNDEFINED;
static float TurbulenceSimulationInitialTemperature   = FLOAT_UNDEFINED;
 
static char *TurbulenceSimulationDensityName          = NULL;
static char *TurbulenceSimulationTotalEnergyName      = NULL;
static char *TurbulenceSimulationGasEnergyName        = NULL;
static char *TurbulenceSimulationVelocityNames[MAX_DIMENSION];
static char *TurbulenceSimulationRandomForcingNames[MAX_DIMENSION];
 
static int   TurbulenceSimulationSubgridsAreStatic    = TRUE;
static int   TurbulenceSimulationNumberOfInitialGrids = 1;
 
 
#define MAX_INITIAL_GRIDS 10
 
int TurbulenceSimulationInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int i, j, dim, gridnum, ret, SubgridsAreStatic, region;
  HierarchyEntry *Subgrid;
 
  char *DensityName = NULL, *TotalEnergyName = NULL, *GasEnergyName = NULL,
       *VelocityNames[MAX_DIMENSION],
       *RandomForcingNames[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    VelocityNames[dim] = NULL;
    RandomForcingNames[dim] = NULL;
  }
 
  /* Set default parameters and names */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    TurbulenceSimulationVelocityNames[dim]       = NULL;
 
  int   TurbulenceSimulationGridDimension[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  int   TurbulenceSimulationGridLevel[MAX_INITIAL_GRIDS];
  FLOAT TurbulenceSimulationGridLeftEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  FLOAT TurbulenceSimulationGridRightEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  for (i = 0; i < MAX_INITIAL_GRIDS; i++)
    TurbulenceSimulationGridLevel[i] = 1;
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    TurbulenceSimulationGridLeftEdge[0][dim] = DomainLeftEdge[dim];
    TurbulenceSimulationGridRightEdge[0][dim] = DomainRightEdge[dim];
    TurbulenceSimulationGridDimension[0][dim] = MetaData.TopGridDims[dim];
  }
  TurbulenceSimulationGridLevel[0] = 0;
 
 
  /* Error check. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if (DualEnergyFormalism == TRUE && HydroMethod != Zeus_Hydro)
      fprintf(stderr, "TurbulenceSimulation: DualEnergyFormalism is ON.\n");
    if (!SelfGravity)
      fprintf(stderr, "TurbulenceSimulation: SelfGravity is OFF\n");
    if (RandomForcing && HydroMethod == Zeus_Hydro) {
      fprintf(stderr, "RandomForcing for Zeus hydro is not tested!\n");
      // ERROR_MESSAGE;
    }
  }
 
  /* Read input from file. */
 
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* Read parameters */
 
    if (sscanf(line, "TurbulenceSimulationDensityName = %s", dummy) == 1)
      TurbulenceSimulationDensityName = dummy;
    if (sscanf(line, "TurbulenceSimulationTotalEnergyName = %s", dummy) == 1)
      TurbulenceSimulationTotalEnergyName = dummy;
    if (sscanf(line, "TurbulenceSimulationGasEnergyName = %s", dummy) == 1)
      TurbulenceSimulationGasEnergyName = dummy;
    if (sscanf(line, "TurbulenceSimulationVelocity1Name = %s", dummy) == 1)
      TurbulenceSimulationVelocityNames[0] = dummy;
    if (sscanf(line, "TurbulenceSimulationVelocity2Name = %s", dummy) == 1)
      TurbulenceSimulationVelocityNames[1] = dummy;
    if (sscanf(line, "TurbulenceSimulationVelocity3Name = %s", dummy) == 1)
      TurbulenceSimulationVelocityNames[2] = dummy;
    if (sscanf(line, "TurbulenceSimulationRandomForcing1Name = %s", dummy) ==1)
      TurbulenceSimulationRandomForcingNames[0] = dummy;
    if (sscanf(line, "TurbulenceSimulationRandomForcing2Name = %s", dummy) ==1)
      TurbulenceSimulationRandomForcingNames[1] = dummy;
    if (sscanf(line, "TurbulenceSimulationRandomForcing3Name = %s", dummy) ==1)
      TurbulenceSimulationRandomForcingNames[2] = dummy;
 
    ret += sscanf(line, "TurbulenceSimulationInitialTemperature = %"FSYM,
                  &TurbulenceSimulationInitialTemperature);
    ret += sscanf(line, "TurbulenceSimulationInitialDensity = %"FSYM,
                  &TurbulenceSimulationInitialDensity);
    ret += sscanf(line, "TurbulenceSimulationNumberOfInitialGrids = %"ISYM,
                  &TurbulenceSimulationNumberOfInitialGrids);
    ret += sscanf(line, "TurbulenceSimulationSubgridsAreStatic = %"ISYM,
                  &TurbulenceSimulationSubgridsAreStatic);
 
    if (sscanf(line, "TurbulenceSimulationGridLeftEdge[%"ISYM"]", &gridnum) > 0)
      ret += sscanf(line, "TurbulenceSimulationGridLeftEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
                    &gridnum, &TurbulenceSimulationGridLeftEdge[gridnum][0],
                    &TurbulenceSimulationGridLeftEdge[gridnum][1],
                    &TurbulenceSimulationGridLeftEdge[gridnum][2]);
    if (sscanf(line, "TurbulenceSimulationGridRightEdge[%"ISYM"]", &gridnum) > 0)
      ret += sscanf(line, "TurbulenceSimulationGridRightEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
                    &gridnum,
		    &TurbulenceSimulationGridRightEdge[gridnum][0],
                    &TurbulenceSimulationGridRightEdge[gridnum][1],
                    &TurbulenceSimulationGridRightEdge[gridnum][2]);
    if (sscanf(line, "TurbulenceSimulationGridDimension[%"ISYM"]", &gridnum) > 0)
      ret += sscanf(line, "TurbulenceSimulationGridDimension[%"ISYM"] = %"ISYM" %"ISYM" %"ISYM,
                    &gridnum,
		    &TurbulenceSimulationGridDimension[gridnum][0],
                    &TurbulenceSimulationGridDimension[gridnum][1],
                    &TurbulenceSimulationGridDimension[gridnum][2]);
    if (sscanf(line, "TurbulenceSimulationGridLevel[%"ISYM"]", &gridnum) > 0)
      ret += sscanf(line, "TurbulenceSimulationGridLevel[%"ISYM"] = %"ISYM,
                    &gridnum, &TurbulenceSimulationGridLevel[gridnum]);
 
 
    /* If the dummy char space was used, then make another. */
 
    if (dummy[0] != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "TurbulenceSimulation")
	&& line[0] != '#')
      fprintf(stderr,
   "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* More error checking. */
 
  if (TurbulenceSimulationVelocityNames[0] == NULL) {
    fprintf(stderr, "Missing initial data.\n");
    ENZO_FAIL("");
  }
 
  if (CellFlaggingMethod[0] != 3)
      fprintf(stderr, "TurbulenceSimulation: check CellFlaggingMethod.\n");
 
  /* If density(temperature, actually, c^2) is left unset, set it = 1. */
 
  if (TurbulenceSimulationInitialDensity == FLOAT_UNDEFINED)
    TurbulenceSimulationInitialDensity = 1.0;
  if (TurbulenceSimulationInitialTemperature == FLOAT_UNDEFINED)
    TurbulenceSimulationInitialTemperature = 1.0;
 
  /* Check/define RandomForcing parameters [Mac Low 1999, ApJ 524, 169]. */
 
  if (RandomForcing)
    if (RandomForcingMachNumber <= 0.0) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Warning: RandomForcing is OFF\n");
      RandomForcing = 0;
    }
 
  /* If RandomForcingEdot is not set in the parameter file, get it
     from [MacLow1999] formula. Note: his formula is calibrated for
     general random forcing fields; coefficient 0.81 can potentially
     be inappropriate for a purely solenoidal forcing; also our
     Gamma is not quite 1.0. */
 
  if (RandomForcing && RandomForcingEdot < 0.0) {
    float BoxSize = TurbulenceSimulationGridRightEdge[0][0] -
                    TurbulenceSimulationGridLeftEdge[0][0];
    float BoxMass = (BoxSize*BoxSize*BoxSize)*
                    TurbulenceSimulationInitialDensity;
    float Vrms    = RandomForcingMachNumber/
      sqrt(TurbulenceSimulationInitialTemperature);
    RandomForcingEdot = 0.81/BoxSize*BoxMass*Vrms*Vrms*Vrms;
 
  /* Approximate correction to the MacLow's factor (see eqs (7) - (8))
     for **this PPM implementation**. Seems to be OK for 64^3, 128^3 and 256^3
     Mach=3,6,10 simulations of **solenoidally** driven turbulence. */
 
    RandomForcingEdot *= 0.8;
 
  }
 
  if (RandomForcing)
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("RandomForcingEdot: %"GSYM"\n", RandomForcingEdot);
 
  /* -------------------------------------------------------------------- */
  /* Generate the root grid and set-up the hierarchy. */
 
  HierarchyEntry *GridsList;
  GridsList = &TopGrid;
 
  /* Initialize the root grid. */
 
  DensityName            = TurbulenceSimulationDensityName;
  TotalEnergyName        = TurbulenceSimulationTotalEnergyName;
  GasEnergyName          = TurbulenceSimulationGasEnergyName;
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    VelocityNames[dim]   = TurbulenceSimulationVelocityNames[dim];
    RandomForcingNames[dim] = TurbulenceSimulationRandomForcingNames[dim];
  }
 
  /* Initialize the root grid by reading in data. */
 
  int TotalRefinement = nint(POW(FLOAT(RefineBy),
                                   TurbulenceSimulationGridLevel[gridnum]));
  if (GridsList->GridData->TurbulenceSimulationInitializeGrid(
			   TurbulenceSimulationInitialDensity,
			   TurbulenceSimulationInitialTemperature,
			   DensityName, TotalEnergyName,
			   GasEnergyName, VelocityNames, RandomForcingNames,
			   TurbulenceSimulationSubgridsAreStatic,
			   TotalRefinement) == FAIL) {
      fprintf(stderr, "Error in grid->TurbulenceSimulationInitializeGrid.\n");
      ENZO_FAIL("");
  }
 
  /* Set boundary conditions if necessary. */
 
  // this will be done in EvolveHierarchy()
 
  /* set up field names and units */
 
  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;
 
  for (j = 0; j < i; j++)
    DataUnits[j] = NULL;
 
  /* Write parameters to parameter output file.
     ATTENTION: printf %s on sun fails if passed a NULL pointer. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "TurbulenceSimulationInitialDensity   = %"FSYM"\n\n",
	    TurbulenceSimulationInitialTemperature);
 
    fprintf(Outfptr, "TurbulenceSimulationInitialTemperature   = %"FSYM"\n\n",
	    TurbulenceSimulationInitialTemperature);
 
    if (TurbulenceSimulationDensityName)
    fprintf(Outfptr, "TurbulenceSimulationDensityName          = %s\n",
	    TurbulenceSimulationDensityName);
    if (TurbulenceSimulationTotalEnergyName)
    fprintf(Outfptr, "TurbulenceSimulationTotalEnergyName      = %s\n",
	    TurbulenceSimulationTotalEnergyName);
    if (TurbulenceSimulationGasEnergyName)
    fprintf(Outfptr, "TurbulenceSimulationGasEnergyName        = %s\n",
	    TurbulenceSimulationGasEnergyName);
    fprintf(Outfptr, "TurbulenceSimulationVelocity1Name        = %s\n",
	    TurbulenceSimulationVelocityNames[0]);
    fprintf(Outfptr, "TurbulenceSimulationVelocity2Name        = %s\n",
	    TurbulenceSimulationVelocityNames[1]);
    fprintf(Outfptr, "TurbulenceSimulationVelocity3Name        = %s\n",
	    TurbulenceSimulationVelocityNames[2]);
    if (TurbulenceSimulationRandomForcingNames[0])
      fprintf(Outfptr, "TurbulenceSimulationRandomForcing1Name = %s\n",
	    TurbulenceSimulationRandomForcingNames[0]);
    if (TurbulenceSimulationRandomForcingNames[1])
      fprintf(Outfptr, "TurbulenceSimulationRandomForcing2Name = %s\n",
	    TurbulenceSimulationRandomForcingNames[1]);
    if (TurbulenceSimulationRandomForcingNames[2])
      fprintf(Outfptr, "TurbulenceSimulationRandomForcing3Name = %s\n",
	    TurbulenceSimulationRandomForcingNames[2]);
  }
 
  /* Clean up. */
 
  delete dummy;
 
  return SUCCESS;
}
 
 
/* --------------------------------------------------------------------
   Re-call the initializer on level zero grids.  Used in case of
   ParallelRootGridIO. */
 
int TurbulenceSimulationReInitialize(HierarchyEntry *TopGrid,
				    TopGridData &MetaData)
{
 
  /* Declarations. */
 
  int dim, gridnum = 0;
  char *DensityName = NULL, *TotalEnergyName = NULL, *GasEnergyName = NULL,
       *ParticlePositionName = NULL, *ParticleVelocityName = NULL,
       *ParticleMassName = NULL, *VelocityNames[MAX_DIMENSION],
       *RandomForcingNames[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    VelocityNames[dim] = NULL;
    RandomForcingNames[dim] = NULL;
  }
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("TurbulenceSimulation: ReInitializing grid %"ISYM"\n", gridnum);
 
  /* If there is more than one grid, add the grid number to the name. */
 
  if (TurbulenceSimulationNumberOfInitialGrids > 1) {
 
    if (TurbulenceSimulationDensityName)
      sprintf(DensityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      TurbulenceSimulationDensityName, gridnum);
    if (TurbulenceSimulationTotalEnergyName)
      sprintf(TotalEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      TurbulenceSimulationTotalEnergyName, gridnum);
    if (TurbulenceSimulationGasEnergyName)
      sprintf(GasEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      TurbulenceSimulationGasEnergyName, gridnum);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      if (TurbulenceSimulationVelocityNames[dim])
	sprintf(VelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		TurbulenceSimulationVelocityNames[dim], gridnum);
      if (TurbulenceSimulationRandomForcingNames[dim])
	sprintf(RandomForcingNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		TurbulenceSimulationRandomForcingNames[dim], gridnum);
    }
 
  } else {
 
    DensityName            = TurbulenceSimulationDensityName;
    TotalEnergyName        = TurbulenceSimulationTotalEnergyName;
    GasEnergyName          = TurbulenceSimulationGasEnergyName;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      VelocityNames[dim]   = TurbulenceSimulationVelocityNames[dim];
      RandomForcingNames[dim] = TurbulenceSimulationRandomForcingNames[dim];
    }
 
  }
 
  /* Call grid initializer.  Use TotalRefinement = -1 to flag real read. */
 
  int TotalRefinement = -1;
 
  /* Loop over level zero grid. */
 
  HierarchyEntry *Temp = TopGrid;
  while (Temp != NULL) {
 
    if (Temp->GridData->TurbulenceSimulationInitializeGrid(
		        TurbulenceSimulationInitialDensity,
		        TurbulenceSimulationInitialTemperature,
		        DensityName, TotalEnergyName,
		        GasEnergyName, VelocityNames, RandomForcingNames,
			TurbulenceSimulationSubgridsAreStatic,
			TotalRefinement) == FAIL) {
      fprintf(stderr, "Error in grid->TurbulenceSimulationInitializeGrid.\n");
      ENZO_FAIL("");
    }
 
    Temp = Temp->NextGridThisLevel;
  }
 
  return SUCCESS;
}
 
