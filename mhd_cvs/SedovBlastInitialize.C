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
#define DEFINE_STORAGE
#include "SedovBlastGlobalData.h"
#undef DEFINE_STORAGE

int SedovBlastInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData)
{

  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";

  /* set up field names and units */
#ifdef ATHENA
  if( EquationOfState == 0 ) {
#endif
    DataLabel[0] = DensName;
    DataLabel[1] = TEName;
    DataLabel[2] = Vel1Name;
    DataLabel[3] = Vel2Name;
    DataLabel[4] = Vel3Name;
#ifdef ATHENA
  }else if ( EquationOfState == 1 ){
    DataLabel[0] = DensName;
    DataLabel[1] = Vel1Name;
    DataLabel[2] = Vel2Name;
    DataLabel[3] = Vel3Name;
  }
#endif

  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  
  if(MHD_Used == TRUE){
    
    MHDcLabel[0] = "MagneticField_C_1";
    MHDcLabel[1] = "MagneticField_C_2";
    MHDcLabel[2] = "MagneticField_C_3";
    
    MHDLabel[0] = "MagneticField_F_1";
    MHDLabel[1] = "MagneticField_F_2";
    MHDLabel[2] = "MagneticField_F_3";
    
    MHDeLabel[0] = "ElectricField_1";
    MHDeLabel[1] = "ElectricField_2";
    MHDeLabel[2] = "ElectricField_3";
    
    MHDcUnits[0] = "FourPiGauss";
    MHDcUnits[1] = "FourPiGauss";
    MHDcUnits[2] = "FourPiGauss";
    
    MHDUnits[0] = "FourPiGauss";
    MHDUnits[1] = "FourPiGauss";
    MHDUnits[2] = "FourPiGauss";
    
    MHDeUnits[0] = "FourPiGauss";
    MHDeUnits[1] = "FourPiGauss";
    MHDeUnits[2] = "FourPiGauss";
    
    CurrentLabel[0] = "Current_1";
    CurrentLabel[1] = "Current_2";
    CurrentLabel[2] = "Current_3";
  }

  /* parameter declarations */

  FLOAT SedovBlastSubgridLeft, SedovBlastSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  
  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
                          SubgridDims[MAX_DIMENSION];
 
  /* make sure it is 2D or 3D */

  if (MetaData.TopGridRank < 2 || MetaData.TopGridRank > 3) {
    printf("Cannot do SedovBlast in %d dimension(s)\n", MetaData.TopGridRank);
    return FAIL;
  }    

  /* There are four parameters:  geometry (cylindrical or spherical symmetry), 
                                 gamma, 
				 E, 
				 rho_1.
     Set their default values */

  float Pi                      = 3.14159;
  float SedovBlastVelocity[3]   = {0.0, 0.0, 0.0};   // gas initally at rest

  SedovBlastDensity             = 1.0;

  float SedovBlastPressure      = 1e-5;

  //If SedovBlastInnerPressure isn't read from the parameter file,
  //calculate it from the Energy.
  int   PressureGiven           = FALSE;
  float SedovBlastInnerPressure;  
  float SedovBlastEnergy        = 1.0;

  float SedovBlastMHDField[3]    = {1e-5, 1e-5,1e-5};// Magnetic Field.  Beta = 2/3

  float SedovBlastRadius        = 3.5;
  float dx = (DomainRightEdge[0] - DomainLeftEdge[0])/
                                                   MetaData.TopGridDims[0];

  /* 3.5 zones on the finest level to resolve the initial explosion */

  float dr = SedovBlastRadius*dx*max(POW(RefineBy,-MaximumRefinementLevel), 0.25);
  /* set no subgrids by default. */

  SedovBlastSubgridLeft         = 0.0;    // start of subgrid(s)
  SedovBlastSubgridRight        = 0.0;    // end of subgrid(s)

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "SedovBlastDensity  = %"FSYM, &SedovBlastDensity);
    ret += sscanf(line, "SedovBlastRadius   = %"FSYM, &SedovBlastRadius);
    ret += sscanf(line, "SedovBlastPressure = %"FSYM, &SedovBlastPressure);
    ret += sscanf(line, "SedovBlastEnergy   = %"FSYM, &SedovBlastEnergy);
    ret += sscanf(line, "SedovBlastMHDField  = %"FSYM" %"FSYM" %"FSYM,
		  SedovBlastMHDField,SedovBlastMHDField+1,SedovBlastMHDField+2);
    ret += sscanf(line, "SedovBlastSubgridLeft = %"FSYM, 
		        &SedovBlastSubgridLeft);
    ret += sscanf(line, "SedovBlastSubgridRight = %"FSYM, 
		        &SedovBlastSubgridRight);
    if(sscanf(line, "SedovBlastInnerPressure = %"FSYM, &SedovBlastInnerPressure) == 1 ){
      PressureGiven = TRUE;
      ret++;
    }
    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "SedovBlast") && 
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, 
	 "warning: the following parameter line was not interpreted:\n%s\n", 
	      line);

  } // end input from parameter file

  //Re scale the subgrid edges to line up with the parent grid.

  int nCells= nint(( SedovBlastSubgridLeft - DomainLeftEdge[0] )/
		   (DomainRightEdge[0]-DomainLeftEdge[0])*MetaData.TopGridDims[0]) - 1;
  SedovBlastSubgridLeft=nCells*(DomainRightEdge[0]-DomainLeftEdge[0])/MetaData.TopGridDims[0];

  nCells = nint(( SedovBlastSubgridRight - DomainLeftEdge[0] )/
		(DomainRightEdge[0]-DomainLeftEdge[0])*MetaData.TopGridDims[0])+1;

  SedovBlastSubgridRight = nCells*(DomainRightEdge[0]-DomainLeftEdge[0])/MetaData.TopGridDims[0];

  
  /* compute p_2 as a function of explosion energy E, initial explosion 
     radius dr, and gamma.

     2D:  p_2 = (gamma-1)*E/(pi*r^2)*rho_in
     3D:  p_2 = (gamma-1)*E/(4/3*pi*r^3)*rho_in
          rho_2 = 1 = rho_1

     If p_2 (inner) is too high, the code crashes in euler.src (dnu<0) after 
     a few time steps. 
     In 2D it is stabe for SedovBlastEnergy <= 0.025 (tested with uniform grid 200^2). 
  */


  if( PressureGiven == FALSE ){
    SedovBlastInnerPressure = 3.0*(Gamma-1.0)*SedovBlastEnergy/
      (MetaData.TopGridRank + 1.0)/
      POW(dr,MetaData.TopGridRank)/Pi;
  }
  /* Check the self-similarity condition: p2/p1 >> (gamma+1)/(gamma-1). */

  float pjump = SedovBlastInnerPressure/SedovBlastPressure;
  if ( pjump < 10.*(Gamma+1)/(Gamma-1) )
    printf("SBI: WARNING! No self-similarity. Pressure jump %g.\n", pjump);

  /* Compute total energies */

  SedovBlastTotalEnergy = SedovBlastPressure/((Gamma - 1.0)*SedovBlastDensity);
  SedovBlastInnerTotalEnergy = SedovBlastInnerPressure/((Gamma - 1.0)*
						   SedovBlastDensity);

  /* Note: As of this writing, the MHD solver still uses Energy Density, not
     specific energy.  Thus, the Density in the following expressoion is only ok 
     because it's 1. Soon I will revert the MHD solver to take specific energy,
     conforming to Enzo.  dcc 2/28/05
  */
  if( MHD_Used ){
    float MagEng = 0.0;
    for(int i=0;i<3;i++)
      MagEng+= SedovBlastMHDField[i]*SedovBlastMHDField[i]/SedovBlastDensity*0.5;

    SedovBlastTotalEnergy+=MagEng;
    SedovBlastInnerTotalEnergy+=MagEng;
  }

  /* set the periodic boundaries */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }

  /* set up uniform grid as of before explosion */

  fprintf(stderr,"SedovInits: Density %f TotalEnergy %f, Velocity %f %f %f B %f %f %f\n",
	  SedovBlastDensity,SedovBlastTotalEnergy,
	  SedovBlastVelocity[0],SedovBlastVelocity[1],SedovBlastVelocity[2],
	  SedovBlastMHDField[0],SedovBlastMHDField[1],SedovBlastMHDField[2]);

  if (TopGrid.GridData->InitializeUniformGrid(SedovBlastDensity, 
					      SedovBlastTotalEnergy,
					      SedovBlastTotalEnergy,
					      SedovBlastVelocity,
					      SedovBlastMHDField) == FAIL) {
    fprintf(stderr, "Error in InitializeUniformGrid.\n");
    return FAIL;
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
      printf("SedovBlast:: Level[%d]: NumberOfSubgridZones[0] = %d\n", lev+1, 
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
							SedovBlastMHDField) == FAIL) {
	fprintf(stderr, "Error in InitializeUniformGrid (subgrid).\n");
	return FAIL;
      }

      /* set up the initial explosion area on the finest resolution subgrid */

      if (lev == MaximumRefinementLevel - 1)
	if (Subgrid[lev]->GridData->SedovBlastInitializeGrid(dr,
				    SedovBlastInnerTotalEnergy, TopGrid.GridData) 
	    == FAIL) {
	  fprintf(stderr, "Error in SedovBlastInitialize[Sub]Grid.\n");
	  return FAIL;
	}
    }
    else
      printf("SedovBlast: single grid start-up.\n");
  }

  /* set up subgrids from level 1 to max refinement level -1 */

  
  int MHD_ProjectEtmp = MHD_ProjectE;
  int MHD_ProjectBtmp = MHD_ProjectB;
  MHD_ProjectE=FALSE;
  MHD_ProjectB=TRUE;

  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
				       *(Subgrid[lev-1]->GridData))
	== FAIL) {
      fprintf(stderr, "Error in ProjectSolutionToParentGrid.\n");
      return FAIL;
    }

  /* set up the root grid */

  if (MaximumRefinementLevel > 0) {
    if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	== FAIL) {
      fprintf(stderr, "Error in ProjectSolutionToParentGrid.\n");
      return FAIL;
    }
  }
  else
    if (TopGrid.GridData->SedovBlastInitializeGrid(dr,
			  SedovBlastInnerTotalEnergy,TopGrid.GridData) == FAIL) {
      fprintf(stderr, "Error in SedovBlastInitializeGrid.\n");
      return FAIL;
    }

  MHD_ProjectE = MHD_ProjectEtmp;
  MHD_ProjectB = MHD_ProjectBtmp;


  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "SedovBlastDensity         = %f\n"  , SedovBlastDensity);
    fprintf(Outfptr, "SedovBlastPressure        = %f\n"  , SedovBlastPressure);
    fprintf(Outfptr, "SedovBlastEnergy          = %f\n"  , SedovBlastEnergy);
    fprintf(Outfptr, "SedovBlastInnerPressure   = %f\n"  , 
	    SedovBlastInnerPressure);

    if(MHD_Used==TRUE){
      fprintf(Outfptr,"SedovBlastMagField        = %f %f %f\n"  , 
	      SedovBlastMHDField[0],SedovBlastMHDField[1],SedovBlastMHDField[2]);
    }
    
  }

  return SUCCESS;

}
