//
// PPML Interface Initialization Routines.
//
// Contains:
//   int grid::PPML_InitInterfaceTypesAndLabels
//          Fills the Grid.FieldType array with the indicies for the Left and Right State pointers.
//          Also fills DataUnits, DataLabel for interface states.
//          Also sets NumberOfFluidQuantites
//          Also sets PPML_NFaces.

//   int PPML_InitInterfaceDataLoop( HierarchyEntry *TopGrid )
//          Loops over top grid and initializes interface data on each grid.
//   int grid::PPML_InitInterfaceDataGrid()
//          Initializes the Left and Right state data.
//          Data must be allocated by Grid_ProblemInitGrid
//          3 options, defined by PPML_InitInterfaceMethod
//          0: No init.  Assumes it's done by the problem initializer (if not done, states are uninitialized)
//          1: Piecewise Constant ( Density_X_Left = Density, etc.)
//          2: Piecewise Linear ( X_Left[i,j,k] = 0.5*(Field[ i,j,k ] + Field[i-1,j,k])


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

#include "PPML.h"

#ifdef SIB2
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#else
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          int level, TopGridData *MetaData,
                          ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#endif
int FindField(int field, int farray[], int numfields);

int PPML_InitInterfaceDataLoop(HierarchyEntry **Grids, int NumberOfGrids, int level,
			       TopGridData *MetaData,ExternalBoundary *Exterior,LevelHierarchyEntry *Level
#ifdef SIB2
			     ,SiblingGridList SiblingList[]
#endif //SIB2
			     ){

  int grid1;

  if( HydroMethod != PPM_Local )
    return SUCCESS;

#ifdef SIB2
    if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			      level, MetaData, Exterior, Level ) == FAIL)
      return FAIL;
#else
    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, Level ) == FAIL)
      return FAIL;
#endif
  

  for( grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (Grids[grid1]->GridData->PPML_InitInterfaceDataGrid() == FAIL ){
      fprintf(stderr," PPML_InitInterfaceDataGrid failed.\n");
      return FAIL;
    }

#ifdef SIB2
    if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			      level, MetaData, Exterior, Level ) == FAIL)
      return FAIL;
#else
    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, Level ) == FAIL)
      return FAIL;
#endif

return SUCCESS;
  
}

 
int grid::PPML_InitInterfaceDataGrid(){

  if( HydroMethod != PPM_Local )
    return SUCCESS;
  if( MyProcessorNumber != ProcessorNumber )
    return SUCCESS;

  if( InterfaceStatesInitialized == TRUE )
    return SUCCESS;

  InterfaceStatesInitialized = TRUE;

  //Set up pointer set.
  PPML_InterfacePointerBundle Face( this );
#ifdef trash
  if( this->ReturnInterfacePointers( Face ) == FAIL ){
    fprintf(stderr," ReturnInterfacePointers failed.\n"); return FAIL;}
#endif //trash
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, BC1Num, BC2Num, BC3Num;
  if (this->IdentifyMHDQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, BC1Num, BC2Num, BC3Num) == FAIL) {
    fprintf(stderr, "PPML_InitializeRoutines (InitInterfaceTypes) Error in IdentifyMHDQuantities.\n");
    return FAIL; }

  int face,field,i,j,k, index1, index2, dim;
  int size=1;
  for( dim=0; dim<GridRank; dim++)
    size *= GridDimension[dim];
  
  //Check for allocation
  for( face=0; face< PPML_NFaces; face++)
    for( field=0; field< NumberOfFluidQuantities; field++)
      if( Face.All[ face ][field] == NULL ){
	fprintf(stderr," PPML_InitInterfaceMethod Problem: Fields must be allocated before calling.\n");
	return FAIL;
      }

  switch( PPML_InitInterfaceMethod ){

    //
    // Disable auto set
    //

  case 0:
    for( face=0; face< PPML_NFaces; face++)
      for( field=0; field< NumberOfFluidQuantities; field++)
	if( Face.All[ face ][field] == NULL ){
	  fprintf(stderr," PPML_InitInterfaceMethod Problem: Auto Init deactivated, but no interface defined.\n");
	  return FAIL;
	}
    break;

    //
    // Piecewise Constant
    //

  case 1:
    for( field=0;field<NumberOfFluidQuantities;field++){
      // X Left
      for( index1 =0; index1<size; index1++)
	Face.X_L[ field ][ index1 ]= BaryonField[ field ][ index1 ];

      // X Right
      for( index1 =0; index1<size; index1++)
	Face.X_R[ field ][ index1 ]= BaryonField[ field ][ index1 ];
      
      // Y Left
      for( index1 =0; index1<size; index1++)
	Face.Y_L[ field ][ index1 ]= BaryonField[ field ][ index1 ];
      

      // Y Right
      for( index1 =0; index1<size; index1++)
	Face.Y_R[ field ][ index1 ]= BaryonField[ field ][ index1 ];

      // Z Left
      for( index1 =0; index1<size; index1++)
	Face.Z_L[ field ][ index1 ]= BaryonField[ field ][ index1 ];
      
      // Z Right
      for( index1 =0; index1<size; index1++)
	Face.Z_R[ field ][ index1 ]= BaryonField[ field ][ index1 ];
      
    }//field
    break;
    //
    // Direct Average.
    //
  case 2:
    for( field=0;field<NumberOfFluidQuantities;field++){
      // X Left
      for( k=0; k<GridDimension[2]; k++)
	for(j=0; j<GridDimension[1]; j++)
	  for(i=1; i<GridDimension[0]; i++){
	    index1 = i     + GridDimension[0]*(j + GridDimension[1]*k);
	    index2 = (i-1) + GridDimension[0]*(j + GridDimension[1]*k);
	    Face.X_L[ field ][ index1 ]= 
	      0.5*(BaryonField[ field ][ index1 ] + BaryonField[ field ][ index2 ]);
	  }
      // X Right
      for( k=0; k<GridDimension[2]; k++)
	for(j=0; j<GridDimension[1]; j++)
	  for(i=0; i<GridDimension[0]-1; i++){
	    index1 = i     + GridDimension[0]*(j + GridDimension[1]*k);
	    index2 = (i+1) + GridDimension[0]*(j + GridDimension[1]*k);
	    Face.X_R[ field ][ index1 ]= 
	      0.5*(BaryonField[ field ][ index1 ] + BaryonField[ field ][ index2 ]);
	  }
      // Y Left
      for( k=0; k<GridDimension[2]; k++)
	for(j=1; j<GridDimension[1]; j++)
	  for(i=0; i<GridDimension[0]; i++){
	    index1 = i + GridDimension[0]*( j    + GridDimension[1]*k);
	    index2 = i + GridDimension[0]*((j-1) + GridDimension[1]*k);
	    Face.Y_L[ field ][ index1 ]= 
	      0.5*(BaryonField[ field ][ index1 ] + BaryonField[ field ][ index2 ]);
	  }
      // Y Right
      for( k=0; k<GridDimension[2]; k++)
	for(j=0; j<GridDimension[1]-1; j++)
	  for(i=0; i<GridDimension[0]; i++){
	    index1 = i + GridDimension[0]*( j    + GridDimension[1]*k);
	    index2 = i + GridDimension[0]*((j+1) + GridDimension[1]*k);
	    Face.Y_R[ field ][ index1 ]= 
	      0.5*(BaryonField[ field ][ index1 ] + BaryonField[ field ][ index2 ]);
	  }
      // Z Left
      for( k=1; k<GridDimension[2]; k++)
	for(j=0; j<GridDimension[1]; j++)
	  for(i=0; i<GridDimension[0]; i++){
	    index1 = i + GridDimension[0]*(j + GridDimension[1]*k);
	    index2 = i + GridDimension[0]*(j + GridDimension[1]*(k-1));
	    Face.Z_L[ field ][ index1 ]= 
	      0.5*(BaryonField[ field ][ index1 ] + BaryonField[ field ][ index2 ]);
	  }
      // Z Right
      for( k=0; k<GridDimension[2]-1; k++)
	for(j=0; j<GridDimension[1]; j++)
	  for(i=0; i<GridDimension[0]; i++){
	    index1 = i + GridDimension[0]*(j + GridDimension[1]*k);
	    index2 = i + GridDimension[0]*(j + GridDimension[1]*(k+1));
	    Face.Z_R[ field ][ index1 ]= 
	      0.5*(BaryonField[ field ][ index1 ] + BaryonField[ field ][ index2 ]);
	  }

    }//field
    break;
  default:
    fprintf(stderr," PPML_InitInterfaceMethod not defined: %d\n", PPML_InitInterfaceMethod);
    return FAIL;
  }//switch
  
  if( RandomForcing == TRUE && PPML_InitInterfaceMethod != 0){
    fprintf(stderr," Re-copying BaryonField to RandomForcingField. \n");
    int counter = 0;
    int vel = Vel1Num; //Kind of lazy.
    for(dim=0; dim < GridRank; dim++){
      for (i = 0; i < size; i++)
	RandomForcingField[counter][i] =  BaryonField[vel+dim][i];
      counter++;
    }
    for(dim=0; dim < GridRank; dim++){
      for (i = 0; i < size; i++)
	RandomForcingField[counter][i] = Face.X_L[vel+dim][i];
      counter++;
    }
    for(dim=0; dim < GridRank; dim++){
      for (i = 0; i < size; i++)
	RandomForcingField[counter][i] = Face.X_R[vel+dim][i];
      counter++;
    }
    for(dim=0; dim < GridRank; dim++){
      for (i = 0; i < size; i++)
	RandomForcingField[counter][i] = Face.Y_L[vel+dim][i];
      counter++;
    }   
    for(dim=0; dim < GridRank; dim++){
      for (i = 0; i < size; i++)
	RandomForcingField[counter][i] = Face.Y_R[vel+dim][i];
      counter++;
    }
    for(dim=0; dim < GridRank; dim++){
      for (i = 0; i < size; i++)
	RandomForcingField[counter][i] = Face.Z_L[vel+dim][i];
      counter++;
    }
    for(dim=0; dim < GridRank; dim++){
      for (i = 0; i < size; i++)
	RandomForcingField[counter][i] = Face.Z_R[vel+dim][i];
      counter++;
    }
    
  }//Random Forcing.
    
  return SUCCESS;

}


int grid::PPML_InitInterfaceTypesAndLabels(){

// Sets up DataLabel, FieldType, and DataUnits.
// Increments NumberOfBaryonFields,
// and sets NumberOfFluidQuantities.
// And number of faces
// Must be called AFTER FieldType is filled
// Must be called on ALL processors, on ALL grids (just like setting NumberOfBaryonFields or FieldType[])

  if( HydroMethod != PPM_Local )
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, BC1Num, BC2Num, BC3Num;
  if (this->IdentifyMHDQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, BC1Num, BC2Num, BC3Num) == FAIL) {
    fprintf(stderr, "PPML_InitializeRoutines (InitInterfaceTypes) Error in IdentifyMHDQuantities.\n");
    return FAIL; }

  NumberOfFluidQuantities = 2; //Density and Vx.
  if( EquationOfState == 0 )  NumberOfFluidQuantities++;
  if( GridRank >= 2 || MHD_Used == 1 )   NumberOfFluidQuantities++;
  if( GridRank >= 3 || MHD_Used == 1 )   NumberOfFluidQuantities++;
  if( MHD_Used == TRUE )   NumberOfFluidQuantities += 3;

  PPML_NFaces = 6;
  
  //Yup.  Horrible.  
  DataLabel[NumberOfBaryonFields]   = "Face_X_L_D";
  FieldType[NumberOfBaryonFields++] =  Face_X_L_D;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ DensNum ];
  if( EquationOfState == 0 ){
    DataLabel[NumberOfBaryonFields]   = "Face_X_L_TE";
    FieldType[NumberOfBaryonFields++] =  Face_X_L_TE;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ TENum ];}
  DataLabel[NumberOfBaryonFields]   = "Face_X_L_VX";
  FieldType[NumberOfBaryonFields++] =  Face_X_L_VX;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel1Num ];
  if( GridRank >= 2 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_X_L_VY";
    FieldType[NumberOfBaryonFields++] =  Face_X_L_VY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel2Num ];}
  if( GridRank >= 3 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_X_L_VZ";
    FieldType[NumberOfBaryonFields++] =  Face_X_L_VZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel3Num ];}
  if( MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_X_L_BX";
    FieldType[NumberOfBaryonFields++] =  Face_X_L_BX;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC1Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_X_L_BY";
    FieldType[NumberOfBaryonFields++] =  Face_X_L_BY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC2Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_X_L_BZ";
    FieldType[NumberOfBaryonFields++] =  Face_X_L_BZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC3Num ];}
  DataLabel[NumberOfBaryonFields]   = "Face_X_R_D";
  FieldType[NumberOfBaryonFields++] =  Face_X_R_D;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ DensNum ];
  if( EquationOfState == 0 ){
    DataLabel[NumberOfBaryonFields]   = "Face_X_R_TE";
    FieldType[NumberOfBaryonFields++] =  Face_X_R_TE;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ TENum ];}
  DataLabel[NumberOfBaryonFields]   = "Face_X_R_VX";
  FieldType[NumberOfBaryonFields++] =  Face_X_R_VX;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel1Num ];
  if( GridRank >= 2 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_X_R_VY";
    FieldType[NumberOfBaryonFields++] =  Face_X_R_VY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel2Num ];}
  if( GridRank >= 3 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_X_R_VZ";
    FieldType[NumberOfBaryonFields++] =  Face_X_R_VZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel3Num ];}
  if( MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_X_R_BX";
    FieldType[NumberOfBaryonFields++] =  Face_X_R_BX;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC1Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_X_R_BY";
    FieldType[NumberOfBaryonFields++] =  Face_X_R_BY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC2Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_X_R_BZ";
    FieldType[NumberOfBaryonFields++] =  Face_X_R_BZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC3Num ];}
  DataLabel[NumberOfBaryonFields]   = "Face_Y_L_D";
  FieldType[NumberOfBaryonFields++] =  Face_Y_L_D;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ DensNum ];
  if( EquationOfState == 0 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Y_L_TE";
    FieldType[NumberOfBaryonFields++] =  Face_Y_L_TE;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ TENum ];}
  DataLabel[NumberOfBaryonFields]   = "Face_Y_L_VX";
  FieldType[NumberOfBaryonFields++] =  Face_Y_L_VX;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel1Num ];
  if( GridRank >= 2 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Y_L_VY";
    FieldType[NumberOfBaryonFields++] =  Face_Y_L_VY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel2Num ];}
  if( GridRank >= 3 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Y_L_VZ";
    FieldType[NumberOfBaryonFields++] =  Face_Y_L_VZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel3Num ];}
  if( MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Y_L_BX";
    FieldType[NumberOfBaryonFields++] =  Face_Y_L_BX;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC1Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_Y_L_BY";
    FieldType[NumberOfBaryonFields++] =  Face_Y_L_BY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC2Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_Y_L_BZ";
    FieldType[NumberOfBaryonFields++] =  Face_Y_L_BZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC3Num ];}
  DataLabel[NumberOfBaryonFields]   = "Face_Y_R_D";
  FieldType[NumberOfBaryonFields++] =  Face_Y_R_D;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ DensNum ];
  if( EquationOfState == 0 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Y_R_TE";
    FieldType[NumberOfBaryonFields++] =  Face_Y_R_TE;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ TENum ];}
  DataLabel[NumberOfBaryonFields]   = "Face_Y_R_VX";
  FieldType[NumberOfBaryonFields++] =  Face_Y_R_VX;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel1Num ];
  if( GridRank >= 2 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Y_R_VY";
    FieldType[NumberOfBaryonFields++] =  Face_Y_R_VY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel2Num ];}
  if( GridRank >= 3 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Y_R_VZ";
    FieldType[NumberOfBaryonFields++] =  Face_Y_R_VZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel3Num ];}
  if( MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Y_R_BX";
    FieldType[NumberOfBaryonFields++] =  Face_Y_R_BX;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC1Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_Y_R_BY";
    FieldType[NumberOfBaryonFields++] =  Face_Y_R_BY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC2Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_Y_R_BZ";
    FieldType[NumberOfBaryonFields++] =  Face_Y_R_BZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC3Num ];}
  DataLabel[NumberOfBaryonFields]   = "Face_Z_L_D";
  FieldType[NumberOfBaryonFields++] =  Face_Z_L_D;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ DensNum ];
  if( EquationOfState == 0 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Z_L_TE";
    FieldType[NumberOfBaryonFields++] =  Face_Z_L_TE;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ TENum ];}
  DataLabel[NumberOfBaryonFields]   = "Face_Z_L_VX";
  FieldType[NumberOfBaryonFields++] =  Face_Z_L_VX;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel1Num ];
  if( GridRank >= 2 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Z_L_VY";
    FieldType[NumberOfBaryonFields++] =  Face_Z_L_VY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel2Num ];}
  if( GridRank >= 3 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Z_L_VZ";
    FieldType[NumberOfBaryonFields++] =  Face_Z_L_VZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel3Num ];}
  if( MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Z_L_BX";
    FieldType[NumberOfBaryonFields++] =  Face_Z_L_BX;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC1Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_Z_L_BY";
    FieldType[NumberOfBaryonFields++] =  Face_Z_L_BY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC2Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_Z_L_BZ";
    FieldType[NumberOfBaryonFields++] =  Face_Z_L_BZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC3Num ];}
  DataLabel[NumberOfBaryonFields]   = "Face_Z_R_D";
  FieldType[NumberOfBaryonFields++] =  Face_Z_R_D;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ DensNum ];
  if( EquationOfState == 0 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Z_R_TE";
    FieldType[NumberOfBaryonFields++] =  Face_Z_R_TE;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ TENum ];}
  DataLabel[NumberOfBaryonFields]   = "Face_Z_R_VX";
  FieldType[NumberOfBaryonFields++] =  Face_Z_R_VX;
  DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel1Num ];
  if( GridRank >= 2 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Z_R_VY";
    FieldType[NumberOfBaryonFields++] =  Face_Z_R_VY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel2Num ];}
  if( GridRank >= 3 || MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Z_R_VZ";
    FieldType[NumberOfBaryonFields++] =  Face_Z_R_VZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ Vel3Num ];}
  if( MHD_Used == 1 ){
    DataLabel[NumberOfBaryonFields]   = "Face_Z_R_BX";
    FieldType[NumberOfBaryonFields++] =  Face_Z_R_BX;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC1Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_Z_R_BY";
    FieldType[NumberOfBaryonFields++] =  Face_Z_R_BY;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC2Num ];
    DataLabel[NumberOfBaryonFields]   = "Face_Z_R_BZ";
    FieldType[NumberOfBaryonFields++] =  Face_Z_R_BZ;
    DataUnits[NumberOfBaryonFields]   =  DataUnits[ BC3Num ];}
  
  return SUCCESS;
}


/*
Code for generating the above atrocity:
 #! /bin/tcsh

foreach i ( `grep Face typedefs.h | awk '{print $1}'` )

echo DataLabel\[NumberOfBaryonFields\]'   '= '"'$i '"'\;
echo FieldType\[NumberOfBaryonFields++\]' '= ' '$i\;
echo DataUnits\[NumberOfBaryonFields\]'   '= ' 'DataUnits\[ \]\;

end
*/ 
