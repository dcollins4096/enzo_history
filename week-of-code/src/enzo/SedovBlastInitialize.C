/***********************************************************************
/
/  INITIALIZE SEDOV BLAST WAVE
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:  Alexei Kritsuk, January 2005. 
/
/  PURPOSE:
/
/   REFERENCE: Self-similar solution: L.I. Sedov (1946); 
/              see also: Sedov (1959), Similarity and Dimensional Methods
/              in Mechanics, pp. 210, 219, 228;
/              see also: Landau & Lifshitz, Fluid Dynamics, Sect. 99 
/              "The Propagation of Strong Shock Waves" (1959).
/              Experiments, terrestrial/numerical: Taylor (1941, 1949).
/
/   Two dimensional parameters: explosion energy E and ambient density rho_1
/   Two independent variables: radius r, time t
/   One dimensionless combination: r*(rho_1/E/t^2)^(1/5)
/
/
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.

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
#define DEFINE_STORAGE
#include "SedovBlastGlobalData.h"
#undef DEFINE_STORAGE

int SedovBlastInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";

  /* parameter declarations */

  FLOAT SedovBlastSubgridLeft, SedovBlastSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  
  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
                          SubgridDims[MAX_DIMENSION];
 
  /* make sure this is 2D or 3D */

  if (MetaData.TopGridRank < 2 || MetaData.TopGridRank > 3) {
    printf("Cannot model SedovBlast in %"ISYM" dimension(s)\n", MetaData.TopGridRank);
    ENZO_FAIL("");
  }    

  /* There are four parameters:
  
     1) geometry (2D-cylindrical or 3D-spherical), 
     2) gamma, 
     3) E, 
     4) rho_1.

     Set their default values here. */

  float Pi                      = 3.14159;
  SedovBlastFullBox             = 0;                 // full box or one quadrant
  float SedovBlastVelocity[3]   = {0.0, 0.0, 0.0};   // gas initally (t=0) at rest
  float SedovBlastBField[3]   = {0.0, 0.0, 0.0};   // gas initally (t=0) at rest
  FLOAT SedovBlastInitialTime   = 0.0;
  float SedovBlastPressure      = 1e-5;              // can be arbitrarily small
  SedovBlastDensity             = 1.0;
  float SedovBlastEnergy        = 1.0;
  float dx = (DomainRightEdge[0] - DomainLeftEdge[0])/
                                                   MetaData.TopGridDims[0];

  /* Use 3.5 zones on the finest level to resolve the initial explosion at t=0. */

  float dr = 3.5*dx*max(POW(RefineBy,-MaximumRefinementLevel), 0.25);

  /* set no subgrids by default. */

  SedovBlastSubgridLeft         = 0.0;    // start of subgrid(s)
  SedovBlastSubgridRight        = 0.0;    // end of subgrid(s)

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "SedovBlastFullBox  = %"ISYM, &SedovBlastFullBox);
    ret += sscanf(line, "SedovBlastInitialTime  = %"FSYM, &SedovBlastInitialTime);
    ret += sscanf(line, "SedovBlastDensity  = %"FSYM, &SedovBlastDensity);
    ret += sscanf(line, "SedovBlastPressure = %"FSYM, &SedovBlastPressure);
    ret += sscanf(line, "SedovBlastEnergy   = %"FSYM, &SedovBlastEnergy);
    ret += sscanf(line, "SedovBlastSubgridLeft = %"FSYM, 
		        &SedovBlastSubgridLeft);
    ret += sscanf(line, "SedovBlastSubgridRight = %"FSYM, 
		        &SedovBlastSubgridRight);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "SedovBlast") && 
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, 
	 "warning: the following parameter line was not interpreted:\n%s\n", 
	      line);

  } // end input from parameter file


  /* Set up current problem time, ambient total energy. */

  MetaData.Time         = SedovBlastInitialTime;
  SedovBlastTotalEnergy = SedovBlastPressure/((Gamma - 1.0)*SedovBlastDensity);

  /* compute p_2 as a function of explosion energy E, initial explosion 
     radius dr, and gamma.

     2D:  p_2 = (gamma-1)*E/(pi*r^2)*rho_in
     3D:  p_2 = (gamma-1)*E/(4/3*pi*r^3)*rho_in
          rho_2 = 1 = rho_1
  */

  float SedovBlastInnerPressure = 1.0;
  if (SedovBlastInitialTime == 0.0) {

    SedovBlastInnerPressure = 3.0*(Gamma-1.0)*SedovBlastEnergy/
      (MetaData.TopGridRank + 1.0)/
      POW(dr,MetaData.TopGridRank)/Pi;

    /* Check the self-similarity condition: p2/p1 >> (gamma+1)/(gamma-1). */

    float pjump = SedovBlastInnerPressure/SedovBlastPressure;
    if ( pjump < 10.0*(Gamma+1)/(Gamma-1) )
      printf("SBI: WARNING! No self-similarity. Pressure jump %"GSYM".\n", pjump);

    /* Compute total energy in the explosion region. */

    SedovBlastInnerTotalEnergy = SedovBlastInnerPressure/((Gamma - 1.0)*
							  SedovBlastDensity);
  }

  /* set periodic boundaries (if FullBox=1),
     otherwise keep reflecting (the default) */

  if (SedovBlastFullBox)
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
      MetaData.RightFaceBoundaryCondition[dim] = periodic;
    }

  /* set up uniform grid as of before explosion */

  if (TopGrid.GridData->InitializeUniformGrid(SedovBlastDensity, 
					      SedovBlastTotalEnergy,
					      SedovBlastTotalEnergy,
					      SedovBlastVelocity,
					      SedovBlastBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }

  /* Create as many subgrids as there are refinement levels 
     needed to resolve the initial explosion region upon the start-up. */

  HierarchyEntry ** Subgrid;
  if (MaximumRefinementLevel > 0) 
    Subgrid   = new HierarchyEntry*[MaximumRefinementLevel];

  /* Create new HierarchyEntries. */

  int lev;
  for (lev = 0; lev < MaximumRefinementLevel; lev++) 
    Subgrid[lev] = new HierarchyEntry;

  for (lev = 0; lev < MaximumRefinementLevel; lev++) {

    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      NumberOfSubgridZones[dim] =
	nint((SedovBlastSubgridRight - SedovBlastSubgridLeft)/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *POW(RefineBy, lev + 1);

    if (debug)
      printf("SedovBlast:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1, 
	     NumberOfSubgridZones[0]);

    if (NumberOfSubgridZones[0] > 0) {

      /* fill them out */

      if (lev == 0)
	TopGrid.NextGridNextLevel  = Subgrid[0];
      Subgrid[lev]->NextGridThisLevel = NULL;
      if (lev == MaximumRefinementLevel-1)
	Subgrid[lev]->NextGridNextLevel = NULL;
      else
	Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];
      if (lev == 0)
	Subgrid[lev]->ParentGrid        = &TopGrid;
      else
	Subgrid[lev]->ParentGrid        = Subgrid[lev-1];

      /* compute the dimensions and left/right edges for the subgrid */

      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*DEFAULT_GHOST_ZONES;
	LeftEdge[dim]    = SedovBlastSubgridLeft;
	RightEdge[dim]   = SedovBlastSubgridRight;
      }

      /* create a new subgrid and initialize it */

      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				     LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(SedovBlastDensity,
						   SedovBlastTotalEnergy,
						   SedovBlastTotalEnergy,
							SedovBlastVelocity,
							SedovBlastBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }

      /* set up the initial explosion area on the finest resolution subgrid */

      if (lev == MaximumRefinementLevel - 1)
	if (SedovBlastInitialTime == 0.0) {
	  if (Subgrid[lev]->GridData->SedovBlastInitializeGrid(dr,
							       SedovBlastInnerTotalEnergy) 
	      == FAIL) {
	    	    ENZO_FAIL("Error in SedovBlastInitialize[Sub]Grid.");
	  }
	}
	else
	  if (Subgrid[lev]->GridData->SedovBlastInitializeGrid3D("sedov.in") 
	      == FAIL) {
	    	    ENZO_FAIL("Error in SedovBlastInitialize3D[Sub]Grid.");
	  }
    }
    else
      printf("SedovBlast: single grid start-up.\n");
  }

  /* set up subgrids from level 1 to max refinement level -1 */

  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
				       *(Subgrid[lev-1]->GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
  
  /* set up the root grid */

  if (SedovBlastInitialTime != 0.0) 
    printf("SBI: Setting up the top grid at time = %"GSYM" based on data from %s.\n",
	   SedovBlastInitialTime, "sedov.in");

  if (MaximumRefinementLevel > 0) {
    if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
  }
  else
    if (SedovBlastInitialTime == 0.0) {
      if (TopGrid.GridData->SedovBlastInitializeGrid(dr,
						     SedovBlastInnerTotalEnergy) == FAIL) {
		ENZO_FAIL("Error in SedovBlastInitializeGrid.");
      }
    }
    else
      if (TopGrid.GridData->SedovBlastInitializeGrid3D("sedov.in") == FAIL) {
		ENZO_FAIL("Error in SedovBlastInitializeGrid3D.");
      }
  

  /* set up field names and units */
  int i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;

  for (int j=0; j< i; j++) 
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "SedovBlastFullBox         = %"ISYM"\n"  , SedovBlastFullBox);
    fprintf(Outfptr, "SedovBlastDensity         = %"FSYM"\n"  , SedovBlastDensity);
    fprintf(Outfptr, "SedovBlastPressure        = %"FSYM"\n"  , SedovBlastPressure);
    fprintf(Outfptr, "SedovBlastEnergy          = %"FSYM"\n"  , SedovBlastEnergy);
    fprintf(Outfptr, "SedovBlastInnerPressure   = %"FSYM"\n"  , 
	    SedovBlastInnerPressure);
  }

  return SUCCESS;

}
