

//  grid::AttachAcceleration  and grid::DetachAcceleration().
// Pointer juggling for the boundary set of the acceleration field.

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

// Begin the pointer juggle to set the boundary on the acceleration field.
// Save all BaryonField pointers in temporary array, and set them to be Acceleration Field
// pointers.  This lets the SetBoundary condition machenery operate without heft code rewrite.


int grid::AttachAcceleration(){


  //This redundancy check is for the parent grid.  Multiple subgrids will have the same 
  //parent grid.

  if( AccelerationHack == TRUE )
    return SUCCESS;
  else
    AccelerationHack = TRUE;

  //fprintf(stderr,"SAB: Attaching\n");


  ActualNumberOfBaryonFields = NumberOfBaryonFields;
  NumberOfBaryonFields = GridRank; 

  for(int field=0;field<ActualNumberOfBaryonFields; field++){
    ActualBaryonField[field] = BaryonField[field];
    ActualOldBaryonField[field] = OldBaryonField[field];
    ActualFieldType[field] = FieldType[field];


    if( field<GridRank ){

      BaryonField[field] = AccelerationField[field];
      OldBaryonField[field] = OldAccelerationField[field];

    }else{
      BaryonField[field] = NULL;
      OldBaryonField[field] = NULL;
    }

    FieldType[field]=FieldUndefined;

  }

  FieldType[0] = ((GridRank >= 1 ) ? Acceleration0 : FieldUndefined );
  FieldType[1] = ((GridRank >= 2 ) ? Acceleration1 : FieldUndefined );
  FieldType[2] = ((GridRank >= 3 ) ? Acceleration2 : FieldUndefined );

  
  return SUCCESS;
}


// end pointer juggle for Boundary Set of AccelerationField.
// Return saved BaryonField pointers to their rightful position.

int grid::DetachAcceleration(){


  if( AccelerationHack == FALSE )
    return SUCCESS;  // the detachment has already been done.
  else
    AccelerationHack = FALSE;
    
  //fprintf(stderr,"SAB: Detach\n");

  NumberOfBaryonFields=ActualNumberOfBaryonFields;

  for( int field=0;field<NumberOfBaryonFields; field++){
    
    BaryonField[field]=ActualBaryonField[field];
    OldBaryonField[field] = ActualOldBaryonField[field];
    FieldType[field] = ActualFieldType[field];
  }


  return SUCCESS;
}
