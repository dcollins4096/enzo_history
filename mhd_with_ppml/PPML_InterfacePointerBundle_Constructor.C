
//
//  Two routines that create and cleanup the PPML_InterfactBundle.
//  

//   struct PPML_InterfacePointerBundle
//          Set of pointers for PPML Interface States.  
//          X_L[ field ] = left x field (where field is density, energy, vx, vy, vz, bx, by,bz)
//          X_R[ field ] = right x field.
//          Y_L, Y_R, Z_L, Z_R: [y,z] [left,right] face 
//          All = {X_L, X_R, Y_L, Y_R, Z_L, Z_R} (for looping over all fields)
//          In the rest of enzo: density is gotten by BaryonField[ DensNum ][ ]
//          X_L[ DensNum ] gives Density on the Left.  (Same indexing strategy)

//   PPML_InterfacePointerBundle::ReturnInterfacePointers( grid * Grid )
//          Identifies Physical quantities (cell centered)
//          Finds indecies for each variable on each face.
//          Assigns X_L[ Density ] = BaryonField[ X_Left_Density ]
//          Makes dealing with these pointers easier.



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

int FindField(int field, int farray[], int numfields);

//Default constructor.  
PPML_InterfacePointerBundle::PPML_InterfacePointerBundle()
{
  PPML_NFaces = 0;
  NumberOfFluidQuantities = 0;
  Error = FALSE;
  All = NULL;
  X_L = NULL;
  X_R = NULL;
  Y_L = NULL;
  Y_R = NULL;
  Z_L = NULL;
  Z_R = NULL;
}

//Here is something dangerous:  A delete extremely close to the acutal data.  
//Like sleeping with a handgun under your pillow.

PPML_InterfacePointerBundle::~PPML_InterfacePointerBundle(){

  int face, field;

  if( this->All != NULL ){
    for( face = 0; face<this->PPML_NFaces; face++){
      if( this->All[face] != NULL ){
	delete this->All[face] ;
	this->All[face] = NULL;
      }
    }
    
    delete this->All;
    this->All = NULL;
  }

}

PPML_InterfacePointerBundle::PPML_InterfacePointerBundle( grid * Grid ){

  if( MyProcessorNumber != Grid->ProcessorNumber ){
    PPML_NFaces = 0;
    NumberOfFluidQuantities = 0;
    Error = TRUE;  //This should never be set without a proper grid.
    All = NULL;
    X_L = NULL;
    X_R = NULL;
    Y_L = NULL;
    Y_R = NULL;
    Z_L = NULL;
    Z_R = NULL;
  }else{
    //Working directly with the BaryonFields for interfaces is asking for trouble.
    //This returns pointers to the Interface States, indexed in the same way as the BaryonFields. 
    //It's a god awful amount of data, so it' going to be god awful code.  Sorry.
    // The order is ALWAYS X_L, X_R, Y_L, Y_R, Z_L, Z_R.
    
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, BX_Num, BY_Num, BZ_Num;
    int face, field;
    
    if( Grid->IdentifyMHDQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				    Vel3Num, TENum,BX_Num,BY_Num,BZ_Num) == FAIL ){
      fprintf(stderr, "Error in Identify Magnetic Quantities.\n");Error=TRUE;  }
    
    
    //<dbg>
    /*
      fprintf(stderr,"monkey d %d  e %d x %d y %d z %d x %d y %d z %d\n", 
      DensNum, TENum, Vel1Num, Vel2Num, Vel3Num, BX_Num, BY_Num, BZ_Num);
      fprintf(stderr,"NFQ %"ISYM"\n", NumberOfFluidQuantities);
    */
    
    //I'm going to loop over the fields, so I need an pointer to pointers to pointers.  Yay!
    this->PPML_NFaces = Grid->PPML_NFaces;
    this->NumberOfFluidQuantities = Grid->NumberOfFluidQuantities;
    this->Error = FALSE;
    this->All = new float ** [PPML_NFaces];
    
    for( face=0; face<PPML_NFaces; face++){
      this->All[face] = new  float *[ NumberOfFluidQuantities];
      for( field=0; field<NumberOfFluidQuantities; field++)
	this->All[face][field] = NULL;
    }
    
    //Find all the indicies for all the various types
    //Not all types are always used, so for Non-MHD runs the BX[] array will be filled w/ -1.
    
    int Densities[6] = { FindField( Face_X_L_D, Grid->FieldType, Grid->NumberOfBaryonFields ),
			 FindField( Face_X_R_D, Grid->FieldType, Grid->NumberOfBaryonFields ),
			 FindField( Face_Y_L_D, Grid->FieldType, Grid->NumberOfBaryonFields ),
			 FindField( Face_Y_R_D, Grid->FieldType, Grid->NumberOfBaryonFields ),
			 FindField( Face_Z_L_D, Grid->FieldType, Grid->NumberOfBaryonFields ),
			 FindField( Face_Z_R_D, Grid->FieldType, Grid->NumberOfBaryonFields ) };
    
    int TE[6] = { FindField( Face_X_L_TE, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_X_R_TE, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_L_TE, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_R_TE, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_L_TE, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_R_TE, Grid->FieldType, Grid->NumberOfBaryonFields ) };
    
    int VX[6] = { FindField( Face_X_L_VX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_X_R_VX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_L_VX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_R_VX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_L_VX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_R_VX, Grid->FieldType, Grid->NumberOfBaryonFields ) };
    
    int VY[6] = { FindField( Face_X_L_VY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_X_R_VY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_L_VY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_R_VY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_L_VY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_R_VY, Grid->FieldType, Grid->NumberOfBaryonFields ) };
    
    int VZ[6] = { FindField( Face_X_L_VZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_X_R_VZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_L_VZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_R_VZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_L_VZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_R_VZ, Grid->FieldType, Grid->NumberOfBaryonFields ) };  
    
    int BX[6] = { FindField( Face_X_L_BX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_X_R_BX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_L_BX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_R_BX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_L_BX, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_R_BX, Grid->FieldType, Grid->NumberOfBaryonFields ) };  
    
    int BY[6] = { FindField( Face_X_L_BY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_X_R_BY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_L_BY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_R_BY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_L_BY, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_R_BY, Grid->FieldType, Grid->NumberOfBaryonFields ) };  
    
    int BZ[6] = { FindField( Face_X_L_BZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_X_R_BZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_L_BZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Y_R_BZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_L_BZ, Grid->FieldType, Grid->NumberOfBaryonFields ),
		  FindField( Face_Z_R_BZ, Grid->FieldType, Grid->NumberOfBaryonFields ) };  
    
    
    //Error check. 
    for( face=0;face<PPML_NFaces; face++){
      if( Densities[face] == -1 ) {
	fprintf(stderr,"PPML_ReturnPointers: Density Failure, face %d\n",face);Error=TRUE;  }
      if( EquationOfState == 0 ){
	if( TE[face] == -1 ) {
	  fprintf(stderr,"PPML_ReturnPointers: TotalEnergy Failure, face %d\n",face);Error=TRUE;  }}
      if( VX[face] == -1 ) {
	fprintf(stderr,"PPML_ReturnPointers: VX Failure, face %d\n",face);Error=TRUE;  }
      if( Grid->GridRank >= 2 || MHD_Used == 1 ){
	if( VY[face] == -1 ) {
	  fprintf(stderr,"PPML_ReturnPointers: VY Failure, face %d\n",face);Error=TRUE;   }}
      if( Grid->GridRank >= 3 || MHD_Used == 1 ){
	if( VX[face] == -1 ) {
	  fprintf(stderr,"PPML_ReturnPointers:  VZ Failure, face %d\n",face);Error=TRUE;   }}
      if( MHD_Used == 1 ){
	if( BX[face] == -1 ) {
	  fprintf(stderr,"PPML_ReturnPointers: BX Failure, face %d\n",face);Error=TRUE;   }
	if( BY[face] == -1 ) {
	  fprintf(stderr,"PPML_ReturnPointers: BY Failure, face %d\n",face);Error=TRUE;   }
	if( BZ[face] == -1 ) {
	  fprintf(stderr,"PPML_ReturnPointers: BZ Failure, face %d\n",face);Error=TRUE;   }}
    }//face loop
    
    for( face = 0; face < PPML_NFaces; face++ ){
      
      this->All[ face ][ DensNum ] = Grid->BaryonField[ Densities[face ] ];
      if( EquationOfState == 0 )
	this->All[ face ][ TENum ] = Grid->BaryonField[ TE[face ] ];
      this->All[ face ][ Vel1Num ] = Grid->BaryonField[ VX[face ] ];
      if( Grid->GridRank >= 2 || MHD_Used == 1 )
	this->All[ face ][ Vel2Num ] = Grid->BaryonField[ VY[face ] ];
      if( Grid->GridRank >= 3 || MHD_Used == 1 )
	this->All[ face ][ Vel3Num ] = Grid->BaryonField[ VZ[face ] ];
      if( MHD_Used == 1 ){
	this->All[ face ][ BX_Num ] = Grid->BaryonField[ BX[face ] ];
	this->All[ face ][ BY_Num ] = Grid->BaryonField[ BY[face ] ];
	this->All[ face ][ BZ_Num ] = Grid->BaryonField[ BZ[face ] ];}
      
      // For debugging.  These should be sequential.
      //fprintf(stderr,"%d Joker d  %d\n", face, Densities[face]);
      //fprintf(stderr,"%d Joker te %d\n", face, TE[face]);
      //fprintf(stderr,"%d Joker vx %d\n", face, VX[face]);
      //fprintf(stderr,"%d Joker vy %d\n", face, VY[face]);
      //fprintf(stderr,"%d Joker vz %d\n", face, VZ[face]);
      
      //For debugging.  These should match output from AllocateGrids.
      /*
	fprintf(stderr,"Joker d  this->All[ %d ][ %d ] %p \n", face, DensNum, this->All[ face ][ DensNum ] );
	if( EquationOfState == 0 )
	fprintf(stderr,"Joker te this->All[ %d ][ %d ] %p \n", face, TENum,   this->All[ face ][ TENum ]);
	fprintf(stderr,"Joker vx this->All[ %d ][ %d ] %p \n", face, Vel1Num, this->All[ face ][ Vel1Num]);
	fprintf(stderr,"Joker vy this->All[ %d ][ %d ] %p \n", face, Vel2Num, this->All[ face ][ Vel2Num]);
	fprintf(stderr,"Joker vz this->All[ %d ][ %d ] %p \n", face, Vel3Num, this->All[ face ][ Vel3Num]);
	if( MHD_Used == TRUE ){
	fprintf(stderr,"Joker bx this->All[ %d ][ %d ] %p \n", face, BX_Num, this->All[ face ][ BX_Num ]);
	fprintf(stderr,"Joker by this->All[ %d ][ %d ] %p \n", face, BY_Num, this->All[ face ][ BY_Num ]);
	fprintf(stderr,"Joker bz this->All[ %d ][ %d ] %p \n", face, BZ_Num, this->All[ face ][ BZ_Num ]);
	}
      */
    }
    
    this->X_L = this->All[0];
    this->X_R = this->All[1];
    this->Y_L = this->All[2];
    this->Y_R = this->All[3];
    this->Z_L = this->All[4];
    this->Z_R = this->All[5];
  }//processor check
}

