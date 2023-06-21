//
// Allocate space for the Flux Correction structure SubgridFluxes.
// This structure is set up in EvolveLevel, allocated in
// SolveHydroEquations (which ever variant is relevant), used in
// UpdateFromFinerGrids (in several places) and deleted in
// DeleteFluxes at the end of EvolveLevel, after UpdateFromFinerGrids.
// Currently, only grid::NewSMHD uses this routine, all other routines
// do this allocation manually.

#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "pout.h"


#ifdef ATHENA
int grid:: AllocateFluxes( int NumberOfSubgrids, fluxes *SubgridFluxes[] ){
  //fprintf(stderr," Allocate Fluxes\n");
  //fprintf(stderr," Check before using.\n");

  //
  //
  //Set up and allocate Subgrid Flux arrays.  These will be filled and
  //dealt with in the solver pass itslef.
  //
  //

  //Subgrid Setup variables
  int subgrid, Axis, Face, field, cell, size=1;


  for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++) {
    for (Axis = 0; Axis < GridRank; Axis++)  {
      
      size = 1;
      for (Face = 0; Face < GridRank; Face++)
	size *= SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[Axis][Face] -
	  SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[Axis][Face] + 1;
      
      
      //Set unused dims to zero.
      for (Face = GridRank; Face < 3; Face++) {
	SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[Axis][Face] = 0;
	SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[Axis][Face] = 0;
	SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[Axis][Face] = 0;
	SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[Axis][Face] = 0;
      }
      
      //Allocate Space and zero out data.
      //JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDBeforeFluxAllocate loop");  
      for (field = 0; field < NumberOfBaryonFields; field++) {
	if (SubgridFluxes[subgrid]->LeftFluxes[field][Axis] == NULL)
	  SubgridFluxes[subgrid]->LeftFluxes[field][Axis]  = new float[size];
	if (SubgridFluxes[subgrid]->RightFluxes[field][Axis] == NULL)
	  SubgridFluxes[subgrid]->RightFluxes[field][Axis] = new float[size];
	//dbg 
	int fucker = 1;
	//
	for (cell = 0; cell < size; cell++) {
	  //<dbg>
	  //SubgridFluxes[subgrid]->LeftFluxes[field][Axis][cell] = fucker;
	  //SubgridFluxes[subgrid]->RightFluxes[field][Axis][cell] =fucker++;
	  //</dbg>

	  SubgridFluxes[subgrid]->LeftFluxes[field][Axis][cell] = 0.0;
	  SubgridFluxes[subgrid]->RightFluxes[field][Axis][cell] = 0.0;

	}

      }


      //JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDBeforeFluxAllocate loop");        
      //Double check that the unused arrays are NULL.
      for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	   field++) {
	SubgridFluxes[subgrid]->LeftFluxes[field][Axis] = NULL;
	SubgridFluxes[subgrid]->RightFluxes[field][Axis] = NULL;
      }

    }//Axis loop

    //For all remaining dimensions (if not a 3d problem)
    for (Axis = GridRank; Axis < 3; Axis++)
      for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
	SubgridFluxes[subgrid]->LeftFluxes[field][Axis] = NULL;
	SubgridFluxes[subgrid]->RightFluxes[field][Axis] = NULL;
      }

    //Now do the ElectricField.
    //(Even though this isn't actually used, since the flux correction
    // for MHD isn't finished.)
    for( Axis=0;Axis < GridRank; Axis++){
      for(field=0;field<3;field++){
	
	size = 1;
	for (Face = 0; Face < GridRank; Face++){
	  size *= (SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[Axis][Face] -
		   SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[Axis][Face] + 1
		   +( ( Axis != field ) ? ( (Face != field)? 1:0) : 0 ) ) ;
	  
	}
	
	SubgridFluxes[subgrid]->LeftElectric[field][Axis] = new float[size];
	SubgridFluxes[subgrid]->RightElectric[field][Axis] = new float[size];	
	
	for( cell=0;cell<size; cell++){
	  SubgridFluxes[subgrid]->LeftElectric[field][Axis][cell] = cell;
	  SubgridFluxes[subgrid]->RightElectric[field][Axis][cell] = cell;
	}//cell
      }//Electric field
    }//Electric Axis
      
      
  } // end of loop over subgrids



  return SUCCESS;
}

#endif //ATHENA
