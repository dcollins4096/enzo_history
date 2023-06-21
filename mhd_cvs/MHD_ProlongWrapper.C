#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "pout.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

// This is a wrapper for the prolongation routine.  Since the prolongation algorithm must
// have the most resolved information for the grid, which may not always be contained in one
// other grid, we must first calculate the derivatives that will be used in mhd_interpolate,
// then do the actual interpolation.
// For historic reasons, this is done entirely in one routine, mhd_interpolate (called from
// grid::ProlongFineGrid)

int grid::MHD_ProlongWrapper(LevelHierarchyEntry * Level){

  LevelHierarchyEntry * OtherLHE;
  int Step;

  // Step = 1 is a flag to mhd_interpolate to only create derivatives.

  MHD_ProlongAllocate(GridDimension);

  //The first loop over the old grids only calculates derivatives.
  //The second uses those derivatives to re-interpolate problem cells.

  Step = 1;

  //fprintf(stderr,"kludge: only doing one grid pair!  This effects MHD_ProlongWrapper as well as REbuild Hierarchy\n");

  OtherLHE = Level;
  while( OtherLHE != NULL){

    if( MHD_ProlongFineGrid( OtherLHE->GridData, Step, "turd" ) == FAIL ){
      fprintf(stderr, "MHDProlongWrapper Derivative Loop.\n");
      return FAIL;
    }

    //if( dccCounter ==0 ) 
    //  OtherLHE = OtherLHE->NextGridThisLevel;
    //else
    //  OtherLHE = NULL;
    OtherLHE = OtherLHE->NextGridThisLevel;
      
  }


  Step = 2;

  OtherLHE = Level;
  while( OtherLHE != NULL){

    if( MHD_ProlongFineGrid( OtherLHE->GridData, Step, "turd" ) == FAIL ){
      fprintf(stderr, "MHDProlongWrapper Fill Loop.\n");
      return FAIL;
    }

    //OtherLHE = OtherLHE->NextGridThisLevel;
    //if( dccCounter ==0 ) 
    //  OtherLHE = OtherLHE->NextGridThisLevel;
    //else
    //  OtherLHE = NULL;
    OtherLHE = OtherLHE->NextGridThisLevel;
      
  }


  //Delete all the temporary memory for the prolongation.
  MHDCleanUpTemp();

  return SUCCESS;
}


