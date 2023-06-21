/***********************************************************************
/
/  OUTPUTS GRID DATA AS PARTICLE DATA
/
/  written by: Jm
/  date:       August, 1996
/  modified1:
/
/  PURPOSE: Calls OutputAsParticleDataHDF4 or OutputAsParticleDataHDF5
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "ExternalBoundary.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"

int OutputAsParticleDataHDF4(TopGridData &MetaData, 
			     LevelHierarchyEntry *LevelArray[],
			     int RegionStart[], int RegionEnd[], 
			     FLOAT RegionStartCoordinates[], 
			     FLOAT RegionEndCoordinates[], int RegionLevel,
			     char *OutputFileName);

int OutputAsParticleDataHDF5(TopGridData &MetaData, 
			     LevelHierarchyEntry *LevelArray[],
			     int RegionStart[], int RegionEnd[], 
			     FLOAT RegionStartCoordinates[], 
			     FLOAT RegionEndCoordinates[], int RegionLevel,
			     char *OutputFileName);

int OutputAsParticleData(TopGridData &MetaData, 
				LevelHierarchyEntry *LevelArray[],
				int RegionStart[], int RegionEnd[], 
				FLOAT RegionStartCoordinates[], 
				FLOAT RegionEndCoordinates[], int RegionLevel,
				char *OutputFileName)
{
#if defined (USE_HDF4)
  return OutputAsParticleDataHDF4 (MetaData, 
				   LevelArray,
				   RegionStart, RegionEnd,
				   RegionStartCoordinates,
				   RegionEndCoordinates, RegionLevel,
				   OutputFileName);
				   
#elif defined (USE_HDF5)
  return OutputAsParticleDataHDF5 (MetaData, 
				   LevelArray,
				   RegionStart, RegionEnd,
				   RegionStartCoordinates,
				   RegionEndCoordinates, RegionLevel,
				   OutputFileName);
				   
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if WRITEGRID_* not set 
#endif

}
