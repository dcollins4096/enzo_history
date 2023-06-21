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
/  GRID CLASS (COPY BARYON FIELD TO OLD BARYON FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

// Copy the current baryon fields to the old baryon fields
//   (allocate old baryon fields if they don't exist).


#include "performance.h"

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
int grid::CopyBaryonFieldToOldBaryonField()
{
  JBMEM_MESSAGE(MyProcessorNumber,"jb: CB1");
  int i, field;

  /* update the old baryon field time */

  OldTime = Time;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* compute the field size */

  /* dcc june 3 2006.  Yeah, this is just wrong.
  int size = 1;
  for (int dim = 0; dim < GridDimension[dim]; dim++)
    size *= GridDimension[dim];
  */
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* copy fields */
    
  for (field = 0; field < NumberOfBaryonFields; field++) {

    /* Check to make sure BaryonField exists. */

    if (BaryonField[field] == NULL) {
      fprintf(stderr, "BaryonField missing.\n");
      return FAIL;
    }

    /* Create OldBaryonField if necessary. */

    if ((OldBaryonField[field] == NULL))
      OldBaryonField[field] = new float[size];

    /* Copy. */

    for (i = 0; i < size; i++)
      OldBaryonField[field][i] = BaryonField[field][i];

  } // end loop over fields
  JBMEM_MESSAGE(MyProcessorNumber,"jb: CB3");

  if(MHD_Used){   
    for(field=0;field<3;field++){

      if(MagneticField[field] == NULL ){
	fprintf(stderr, "MagneticField mising in CopyBaryonFieldToOld \n");
	return FAIL;
      }
      
      if(OldMagneticField[field] == NULL) {
	OldMagneticField[field] = new float[MagneticSize[field]];
      }
      
      for(i=0;i<MagneticSize[field];i++){
	OldMagneticField[field][i] = MagneticField[field][i];
      }

      if(CenteredB[field] == NULL) {
	fprintf(stderr, "CenteredB missing in CopyBaryonFieldToOld... \n");
	return FAIL;
      }
      
      if(OldCenteredB[field] == NULL) {
	OldCenteredB[field] = new float[size];
      }
      
      for(i=0;i<size;i++){
	OldCenteredB[field][i] = CenteredB[field][i];
      }
/*      
	if(ElectricField[field] == NULL ){
	fprintf(stderr, "ElectricField Missing in CopyBaryonFieldToOld... \n");
	return FAIL;
	}
	
	if(OldElectricField[field] == NULL) {
	OldElectricField[field] = new float[ElectricSize[field]];
	}
	
	for(i=0;i<ElectricSize[field];i++){
	OldElectricField[field][i] = ElectricField[field][i];
	}
*/

    }//for(field < 3;)
  }//end if(MHD_Used)

  //for acceleration hacks:

  if( (SelfGravity || UniformGravity || PointSourceGravity) )
    for(field=0;field<GridRank;field++){
      if(AccelerationField[field] != NULL){
	if( OldAccelerationField[field] == NULL ){
	  OldAccelerationField[field] = new float[size];

	}
	for(i=0;i<size;i++)
	  OldAccelerationField[field][i] = AccelerationField[field][i];

      }else{
	fprintf(stderr,"Shit-- in CopyBF to Old, no AccelerationField.\n");
//	return FAIL;
      }
    }//field, GravityFlags.



  return SUCCESS;

}
