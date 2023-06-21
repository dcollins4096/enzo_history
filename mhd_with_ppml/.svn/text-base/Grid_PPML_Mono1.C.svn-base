#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "PPML.h"

extern "C" void FORTRAN_NAME(monotx)(int * mk,int * Nx,float * Qpr,float * Qpl,float * Qpm,
				     float * a2_sound,float * c0);

int grid::PPML_Mono1(){

  if( ProcessorNumber != MyProcessorNumber )
    return SUCCESS;

  //fprintf(stderr," PPM Mono1 OFF\n");
  //return SUCCESS;
  int i,j,k,field, index, *MapToPPML;

  int dim;
  float * Qpr, * Qpl,* Qpm; //Right, Left and Center stripes.
  float a2 = IsothermalSoundSpeed*IsothermalSoundSpeed;
  //
  // Layer between Physical quantities and Data structures.
  // Not in the slightest bit pretty.
  //

  PPML_InterfacePointerBundle Face( this );
  if( Face.Error == TRUE ){
    fprintf(stderr, "Error initializing interface pointer bundle.\n");
    return FAIL;
  }
#ifdef trash
  if( this->ReturnInterfacePointers( Face ) == FAIL ){
    fprintf(stderr," ReturnInterfacePointers failed.\n"); return FAIL;}
#endif //trash
  IndexPointerMap ind;
  if( this->IdentifyPhysicalQuantities_2(ind) == FAIL ){
    fprintf(stderr," IdentifyPhysicalQuantities_2 failed\n");return FAIL;}

  //
  // To ease our coding, we also make a map between PPML indexing and Enzo indexing.
  // It takes the PPML index, and returns the Enzo index.

  /*
    fprintf(stderr,"monkey d %d  e %d x %d y %d z %d x %d y %d z %d\n", 
    ind.D, ind.TE, ind.VX, ind.VY, ind.VZ, ind.BX, ind.BY, ind.BZ);
  */

  MapToPPML = new int[ NumberOfFluidQuantities] ;
  i = 0;
  MapToPPML[ i++ ] = ind.D;
  //This will have to move.
  if( EquationOfState == 0 )
    {MapToPPML[ i++ ] = ind.TE;
    fprintf(stderr," Mono1: PPML doesn't yet support TotalEnergy, so you'd better not have it in the map.\n");}
  MapToPPML[ i++ ] = ind.VX;
  MapToPPML[ i++ ] = ind.VY;
  MapToPPML[ i++ ] = ind.VZ;
  if( MHD_Used == TRUE ){
    MapToPPML[ i++ ] = ind.BX;
    MapToPPML[ i++ ] = ind.BY;
    MapToPPML[ i++ ] = ind.BZ;
  }

  /*
    for(j=0;j<7;j++)
    fprintf(stderr,"MapToPPML[%d] = %d\n", j, MapToPPML[j]);
  */

  dim = 1;
  Qpr = new float[ GridDimension[0]*NumberOfFluidQuantities ];
  Qpl = new float[ GridDimension[0]*NumberOfFluidQuantities ];
  Qpm = new float[ GridDimension[0]*NumberOfFluidQuantities ];

  for(k=0; k<GridDimension[2]; k++){
    for(j=0; j<GridDimension[1]; j++){

      for(field=0;field<NumberOfFluidQuantities; field++)
	for(i=0; i<GridDimension[0]; i++){
	  index = i + GridDimension[0]*(j + GridDimension[1]*k);
	  Qpr[i + GridDimension[0]*field] = Face.X_R[ MapToPPML[ field ] ][ index ];
	  Qpl[i + GridDimension[0]*field] = Face.X_L[ MapToPPML[ field ] ][ index ];
	  Qpm[i + GridDimension[0]*field] = BaryonField[ MapToPPML[ field ] ][ index ];
	}//field, i

      FORTRAN_NAME(monotx)(&dim, GridDimension, Qpr, Qpl, Qpm, &a2, &CourantSafetyNumber);

      for(field=0;field<NumberOfFluidQuantities; field++)
	for(i=0; i<GridDimension[0]; i++){
	  index = i + GridDimension[0]*(j + GridDimension[1]*k);
	  Face.X_R[ MapToPPML[ field ] ][ index ] = Qpr[i + GridDimension[0]*field];
	  Face.X_L[ MapToPPML[ field ] ][ index ] = Qpl[i + GridDimension[0]*field];
	}//field, i

	
    }//j
  }//k

  delete Qpr;
  delete Qpl;
  delete Qpm;


  //
  // Y face
  //

  if( GridRank > 1 ){
    dim = 2; //this is used by Monotx
    Qpr = new float[ GridDimension[1]*NumberOfFluidQuantities ];
    Qpl = new float[ GridDimension[1]*NumberOfFluidQuantities ];
    Qpm = new float[ GridDimension[1]*NumberOfFluidQuantities ];
    
    for(k=0; k<GridDimension[2]; k++){
      for(i=0; i<GridDimension[0]; i++){
	
	for(field=0;field<NumberOfFluidQuantities; field++)
	  for(j=0; j<GridDimension[1]; j++){
	    index = i + GridDimension[0]*(j + GridDimension[1]*k);
	    Qpr[j + GridDimension[1]*field] = Face.Y_R[ MapToPPML[ field ] ][ index ];
	    Qpl[j + GridDimension[1]*field] = Face.Y_L[ MapToPPML[ field ] ][ index ];
	    Qpm[j + GridDimension[1]*field] = BaryonField[ MapToPPML[ field ] ][ index ];
	  }//field, i
	
	FORTRAN_NAME(monotx)(&dim, GridDimension +1, Qpr, Qpl, Qpm, &a2, &CourantSafetyNumber);
	
	for(field=0;field<NumberOfFluidQuantities; field++)
	  for(j=0; j<GridDimension[1]; j++){
	    index = i + GridDimension[0]*(j + GridDimension[1]*k);
	    Face.Y_R[ MapToPPML[ field ] ][ index ] = Qpr[j + GridDimension[1]*field] ;
	    Face.Y_L[ MapToPPML[ field ] ][ index ] = Qpl[j + GridDimension[1]*field];
	  }//field, j
	
      }//j
    }//k
    
    delete Qpr;
    delete Qpl;
    delete Qpm;
  }//rank


  //
  // Z sweep
  //
  if( GridRank > 2 ){
    dim = 3; //this is used by Monotx
    Qpr = new float[ GridDimension[2]*NumberOfFluidQuantities ];
    Qpl = new float[ GridDimension[2]*NumberOfFluidQuantities ];
    Qpm = new float[ GridDimension[2]*NumberOfFluidQuantities ];
    
    
    for(j=0; j<GridDimension[1]; j++){
      for(i=0; i<GridDimension[0]; i++){
	
	for(field=0;field<NumberOfFluidQuantities; field++)
	  for(k=0; k<GridDimension[2]; k++){
	    index = i + GridDimension[0]*(j + GridDimension[1]*k);
	    Qpr[k + GridDimension[2]*field] = Face.Z_R[ MapToPPML[ field ] ][ index ];
	    Qpl[k + GridDimension[2]*field] = Face.Z_L[ MapToPPML[ field ] ][ index ];
	    Qpm[k + GridDimension[2]*field] = BaryonField[ MapToPPML[ field ] ][ index ];
	  }//field, i
	
	FORTRAN_NAME(monotx)(&dim, GridDimension+2, Qpr, Qpl, Qpm, &a2, &CourantSafetyNumber);
	for(field=0;field<NumberOfFluidQuantities; field++)
	  for(k=0; k<GridDimension[2]; k++){
	    index = i + GridDimension[0]*(j + GridDimension[1]*k);
	    Face.Z_R[ MapToPPML[ field ] ][ index ] =Qpr[k + GridDimension[2]*field];
	    Face.Z_L[ MapToPPML[ field ] ][ index ] =Qpl[k + GridDimension[2]*field];
	  }//field, i
	
	
      }//j
    }//k
    
    delete Qpr;
    delete Qpl;
    delete Qpm;
  }
  delete MapToPPML;

#ifdef trash
  if( this->DeleteInterfacePointers( Face ) == FAIL ){
    fprintf(stderr,"Grid_PPML_Mono1: error with Delete InterfacePointers");
    return FAIL;
  }
#endif //trash
  return SUCCESS;
}
