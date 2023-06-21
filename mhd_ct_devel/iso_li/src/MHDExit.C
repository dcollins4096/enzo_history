/* Not all features of Enzo are supported by the MHD code.
   This routine exits enzo if any of them are used.*/

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

int MHDExit(grid *TopGrid, TopGridData * MetaData)
{
  //test comment
  //another test comment
  if( ! MHD_Used ) return SUCCESS;

#if defined (USE_HDF4)
  fprintf(stderr, "\n\n++++++++++++++++++ MHD doesn't support HDF4 (MHDExit.C)+++++++++++++\n\n");
  return FAIL;
#endif
#if defined (WRITEGRID_FORTRAN)
  fprintf(stderr, "\n\n+++++++++++++++++ Output must be in HDF5, or write your own. (MHDExit.C)+++++++++++++\n\n");
  return FAIL;
#endif

#ifdef ATHENA
  if( TopGrid->GridRank != 3 && HydroMethod != Athena && HydroMethod != MHD_None )
    {
      fprintf(stderr," Currenly, only the Athena MHD method works with rank != 3.\n");

       return FAIL;
    }
#else
  if( TopGrid->GridRank != 3 )
    {
      fprintf(stderr," Currenly, only the Athena MHD method works with rank != 3.\n");

       return FAIL;
    }

#endif
  if( MHD_DivB == MHD_DivB_Poisson )
    {
      fprintf(stderr, "\n\n ++++++++++++++++ You need to write the divergence cleaner +++++++++++\n");
      return FAIL;
    }

 if( PressureFree == TRUE )
    {
      fprintf(stderr, "\n\n +++++++++++++++ Pressure Free Not Supported by this release of MHD \n");
      return FAIL;
    }

#ifdef HAOXU
 
#else /* HAOXU */
if(DualEnergyFormalism == TRUE )
    {
      fprintf(stderr, "\n\n +++++++++++++++ Dual Energy Formalism isn't supported in this MHD release \n");
      return FAIL;
    }

 
  if( ComovingCoordinates == TRUE )
    {
      fprintf(stderr, "\n\n +++++++++++++++ Comoving Coordinates not supported by this MHD release \n");
      return FAIL;
    }
#endif /* HAOXU */

  if( RefineBy != 2 ){
    fprintf(stderr, "No refinement other than 2. You have %d\n", RefineBy);
    return FAIL;
  }

  for( int dim=0;dim<TopGrid->GridRank;dim++)
    if( TopGrid->GridDimension[dim] < DEFAULT_GHOST_ZONES 
	&& MetaData-> LeftFaceBoundaryCondition[dim] == 3 && MetaData-> RightFaceBoundaryCondition[dim] == 3  ){
      fprintf(stderr," SEVERE WARNING!!!\n");
      fprintf(stderr," Your top grid is smaller than the number of ghost zones, and your boundary is periodic. \n");
  fprintf(stderr," This cause significant, wierd, random problems because the boundary wants more data than is available");
      fprintf(stderr," Dim Grid[%d] = %d, GhostZones = %d\n",dim, TopGrid->GridDimension[dim], DEFAULT_GHOST_ZONES);

      return FAIL;
    }


  if( MaximumRefinementLevel != 0 ){
    if (MHD_ProjectE != 1 ){
      fprintf(stderr," Error:  You have AMR on and ProjectE not on.  Either this is a mistake, or you've finished\n");
      fprintf(stderr," the flux correction routine.  In that case, remove this warning.\n");
      return FAIL;
    }
    if (MHD_ProjectB == 1){
      fprintf(stderr," Error:  You have AMR on and ProjectB  on.  Either this is a mistake, or you've finished\n");
      fprintf(stderr," the flux correction routine.  In that case, remove this warning.\n");
      return FAIL;
    }
  }//amr?
  return SUCCESS;



}
