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
extern "C" void FORTRAN_NAME(ppml_barth)(int * nx,int * ny,int * nz,
					 float * Qplx, float * Qprx,
					 float * Qply, float * Qpry,
					 float * Qplz, float * Qprz,
					 float * Qpm);
int grid::PPML_Mono2(){

  if( ProcessorNumber != MyProcessorNumber )
    return SUCCESS;

  if( GridRank != 3 ){
    fprintf(stderr,"PPML_Mono2: Rank != 3.  Exiting until proper fix to the limiter.\n");
    return SUCCESS;
  }
  
  int field;
  //fprintf(stderr," PPM Mono2 OFF\n");
  //return SUCCESS;
  PPML_InterfacePointerBundle Face( this );
  if( Face.Error == TRUE )
    ERROR_MESSAGE;
#ifdef trash
  if( this->ReturnInterfacePointers( Face ) == FAIL ){
    fprintf(stderr," ReturnInterfacePointers failed.\n"); return FAIL;}
#endif //trash  
  for(field=0; field<NumberOfFluidQuantities; field++)
    FORTRAN_NAME(ppml_barth)(GridDimension,GridDimension+1,GridDimension+2,
			     Face.X_L[field], Face.X_R[field],
			     Face.Y_L[field], Face.Y_R[field],
			     Face.Z_L[field], Face.Z_R[field],
			     BaryonField[field]);

#ifdef trash
  if( this->DeleteInterfacePointers( Face ) == FAIL ){
    fprintf(stderr,"Grid_PPML_Mono1: error with Delete InterfacePointers");
    return FAIL;
  }
#endif //trash
  return SUCCESS;
}
