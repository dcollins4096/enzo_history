#include <math.h>
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

//
// Rotating Disk initializer.
// Sets up a uniform density, pressure, and magnetic field.
// Velocity field has rotation (DiskRotation), DiskCompression, DiskTranslation.
// Rotated by an angle DiskTheta about the Y axis 
// With Gaussian falloffs of DiskRotationExponent and DiskCompressionExponent
// 


int DiskInitialize(FILE *fptr, FILE *Outfptr,
		   HierarchyEntry &TopGrid, TopGridData &MetaData)
{

#ifndef ATHENA   
   int  EquationOfState = 0 ;
   float IsothermalSoundSpeed = 0;
#endif

  //Physics Parameters.
  float 
    DiskTheta = 0,
    DiskWaveNumber = 1,
    DiskRotation = 0,
    DiskCompression = 0.01,
    DiskMagneticField[3] = {2.0, 0.0, 2.0};
  float DiskTranslation[3] = {2.0, 0.0, 2.0};
  float DiskCenter[3] = {0.5,0.5,0.5};
  float DiskRotationExp = 3,
    DiskCompressionExp = 6,
    DiskDensity  = 30,
    DiskPressure = 30;
  
  /*
  // Setup parameters
  int DiskRefineOnStartup = 0;
  float DiskSubgridLeft[3], DiskSubgridRight[3]
  */
  //
  //
  // Labels and Units.  (For IO.)
  // 
  
  char *DensName = "Density";
  
  char *TEName = "Total_Energy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  
#ifdef HAOXU
  if(DualEnergyFormalism ){
    char *GEName = "GasEnergy";
    DataLabel[5] = GEName;
    DataUnits[5] = NULL;   
  }
#endif
  
  if( EquationOfState == 0 ){
    DataLabel[0] = DensName;
    DataLabel[1] = TEName;
    DataLabel[2] = Vel1Name;
    DataLabel[3] = Vel2Name;
    DataLabel[4] = Vel3Name;
  }else if (EquationOfState == 1){
    DataLabel[0] = DensName;
    DataLabel[1] = Vel1Name;
    DataLabel[2] = Vel2Name;
    DataLabel[3] = Vel3Name;
  }
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  
  
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
  
  //Control for IO.

  char line[MAX_LINE_LENGTH];
  int ret=0;
  
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;
    ret += sscanf(line, "DiskPressure = %"PSYM, &DiskPressure);
    ret += sscanf(line, "DiskDensity = %"PSYM, &DiskDensity);
    ret += sscanf(line, "DiskTheta = %"PSYM, &DiskTheta);
    ret += sscanf(line, "DiskWaveNumber = %"PSYM, &DiskWaveNumber);
    ret += sscanf(line, "DiskRotation = %"PSYM, &DiskRotation);
    ret += sscanf(line, "DiskCompression = %"PSYM, &DiskCompression);
    ret += sscanf(line, "DiskRotationExp = %"PSYM, &DiskRotationExp);
    ret += sscanf(line, "DiskCompressionExp = %"PSYM, &DiskCompressionExp);
    ret += sscanf(line, "DiskTranslation = %"PSYM" %"PSYM" %"PSYM, 
		  DiskTranslation, DiskTranslation+1, DiskTranslation+2);
    ret += sscanf(line, "DiskMagneticField = %"PSYM" %"PSYM" %"PSYM, 
		  DiskMagneticField, DiskMagneticField+1, DiskMagneticField+2);
    ret += sscanf(line, "DiskCenter = %"PSYM" %"PSYM" %"PSYM, 
		  DiskCenter, DiskCenter+1, DiskCenter+2);

    if (strstr(line, "Disc") ){
      fprintf(stderr,"Warning: DiskInitialize spells things with a k:\n");
      fprintf(stderr,"         %s\n", line);
    }
  }//read loop



  if( TopGrid.GridData->DiskInitializeGrid(DiskDensity,
					   DiskPressure,
					   DiskTheta,
					   DiskWaveNumber,
					   DiskRotation,
					   DiskCompression,
					   DiskTranslation,
					   DiskMagneticField,
					   DiskCenter,
					   DiskRotationExp,
					   DiskCompressionExp) == FAIL ){
    fprintf(stderr,"Error: DiskInitializeGrid Failed.\n");
    return FAIL;
  }

  /*
  //
  // Generate Hierarchy.  Clone only, will need to be finished.
  //
  if( RefineOnStartup == 1 ){
  //Create as many subgrids as there are refinement levels 
  //needed to resolve the initial explosion region upon the start-up. 
  
  HierarchyEntry ** Subgrid;
  if (MaximumRefinementLevel > 0) 
  Subgrid   = new HierarchyEntry*[MaximumRefinementLevel];
  
  //
  //Create new HierarchyEntries.  Note that 'lev' loops are only for the SUBGRIDS.
  //
  
  int lev;
  int NumberOfSubgridZones[3], SubgridDims[3];
  
  for (lev = 0; lev < MaximumRefinementLevel; lev++) 
  Subgrid[lev] = new HierarchyEntry;
  
  for (lev = 0; lev < MaximumRefinementLevel; lev++) {
  
  //Calculate number of cells on this level.
  
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
  NumberOfSubgridZones[dim] = nCells[dim]*POW(RefineBy, lev + 1);
  
  fprintf(stderr,"DiskInit:: Level[%d]: NumberOfSubgridZones[0] = %d\n", lev+1, 
  NumberOfSubgridZones[0]);
  
  if (NumberOfSubgridZones[0] > 0) {
  
  // fill them out 
  
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
  
  //  compute the dimensions and left/right edges for the subgrid 
  
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
  SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*DEFAULT_GHOST_ZONES;
  }
  
  // create a new subgrid and initialize it 
  
  Subgrid[lev]->GridData = new grid;
  Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
  Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
  MHDBlastSubgridLeft,MHDBlastSubgridRight, 0);
  
  if( TopGrid.GridData->DiskInitializeGrid(DiskDensity,
  DiskPressure,
  DiskTheta,
  DiskWaveNumber,
  DiskRotation,
  DiskCompression,
  DiskTranslation,
  DiskMagneticField,
  DiskCenter,
  DiskRotationExp,
  DiskCompressionExp) == FAIL ){
  fprintf(stderr,"Error: DiskInitializeGrid Failed: SubGrid.\n");
  return FAIL;
  }	
  
  
  }//NumberOfSubgridZones > 0
  else{
  printf("DiskInit: single grid start-up.\n");
  }
  
  }//level
  
  // Make sure each grid has the best data with respect to the finer grids.
  // This projection juggle is to ensure that, regardless of how the hierarchy is evolved, the field gets projected
  // properly here.
  
  
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
  
  // set up the root grid 
  
  if (MaximumRefinementLevel > 0) {
  if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
  == FAIL) {
  fprintf(stderr, "Error in ProjectSolutionToParentGrid.\n");
  return FAIL;
  }
  }
  
  // Put the projection options back to the inital.
  MHD_ProjectE = MHD_ProjectEtmp;
  MHD_ProjectB = MHD_ProjectBtmp;
  
  }//RefineOnStartup
  */
  return SUCCESS;

}
