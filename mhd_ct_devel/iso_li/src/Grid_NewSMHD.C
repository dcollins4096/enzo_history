//
// New version of SolveMHDEquations. This version is for codes that
// only call C solvers.  It may or may not get wrapped into the
// regular SMHD, but it's easier for me to develop without having to
// deal with switching out the code I don't want.
//
// Flow:
// 1.) Set up Subgrid Flux arrays
// 2.) Set up Electric Field
// 3.) Call Solver
// 4.) Clean up leftovers.
//


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
int grid::NewSMHD(int CycleNumber, int NumberOfSubgrids, 
		  fluxes *SubgridFluxes[], ExternalBoundary *Exterior, int level, int grid,
		  float * RandomForcingNormalization, float TopGridTimeStep){


  //This routine only gets called by processors that have data.
  if( ProcessorNumber != MyProcessorNumber )
    return SUCCESS;
  
  //
  //Set up variables we'll need.
  //

  int size, dim;

  //Label for WriteInThis
  char basename[20];
  int MyGridNumber = 1; //THe grid number for output.  Use an integer, level, processor, bananna, whatever.

  //Some IO, for flow tracing.
  fprintf(stderr,"===  NewSMHD n = %d L = %d g = %d proc = %d dt = %15.12e ===\n",
	  CycleNumber, level, grid, MyProcessorNumber, dtFixed);
  Pout(" +++++ Solve MHD Equations (cycle, level, grid ) +++++ ", CycleNumber, level, grid);
  wall_time("Start NewSMHD");

  //
  //Allocate subgrid fluxes. (NOT the solver fluxes.)
  //

  AllocateFluxes(NumberOfSubgrids, SubgridFluxes);

  //
  // Write In This is good for debugging.
  //


  if( WriteInThisF(30) == TRUE) {
    sprintf(basename, "data30%d%d%d.grid",CycleNumber, level, grid);
    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, MyGridNumber) == FAIL ){
      fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
      return FAIL;
    }
    fclose(dummy);
  }  

  //
  // Other debugging things
  //


  //
  // Call the solver itself.
  //

  switch( HydroMethod ){

  case Athena:

    if( MHD_Athena(CycleNumber, level, grid,NumberOfSubgrids, SubgridFluxes,
		   RandomForcingNormalization,TopGridTimeStep) == FAIL ) 
      { fprintf(stderr,"shit.  Athena failed.\n"); return FAIL;}
    break;


  case MHD_None:
    fprintf(stderr,"=============== NASTY KLUDGE!!! NO SOLVER!!! (NewSMHD) ==================\n");

    break;

  default:
    if(MyProcessorNumber == ROOT_PROCESSOR) 
      fprintf(stderr, "SolveMHDEquations:  Hydro Method is not a defined MHD Method\n");
    return FAIL;
    
  }


  //
  // Write In This is good for debugging.
  //

  if( WriteInThisF(39) == TRUE){
    int writetmp = MHD_WriteElectric;
    MHD_WriteElectric=TRUE;
    sprintf(basename, "data39%d%d%d.grid",CycleNumber, level,grid);
    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, MyGridNumber) == FAIL ){
      fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
      return FAIL;
    }
    fclose(dummy);
    MHD_WriteElectric=writetmp;
  }  


  wall_time("End NewSMHD");

  return SUCCESS;
}
#endif //ATHENA
