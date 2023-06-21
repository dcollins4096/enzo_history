#ifdef PPML
//
// Computes global statistical information.
// Currently called every timestep.  Expensive-- more elaborate machinery will be needed in the future.
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "GlobalStats.h"

int ComputeGlobalStats(LevelHierarchyEntry *LevelArray[],
		       int level, TopGridData *MetaData){

  /* Work flow:
     1.) Flag for top grid: not defined for subgrids.
     2.) Allocate static GlobalStats Stats.  Pray that things work properly.
     3.) Loop with level array, call Grid_InlineStatistics( GlobalStats Stats)
     4.) Sum over summed variables, max over maxed variables.
     5.) Sqrts for those that need it.
     6.) Open file, append.
     7.) write data to file
     8.) close file.

  */

  //1.) Unclear how to do this in AMR.  Probably use Ricks tools.
  if( level != 0 ) {
    return SUCCESS;
  }

  //2.) Declared static in order to waste memory.  
  //    I mean, in order to keep track of when the first call is made, so header info
  //    is only printed once.
  static GlobalStats Stats;
  Stats.Clear();
  //2.5) Register all quantities.
  //     Adding a new quantity is easy:  Declare it's name (used here
  //and in Grid_ComputeGlobalStatsGrid for reference, and it labels
  //the output) , outputformat (just like in printf.  Don't forget a
  //space afterwards!) , and type of communcation here (options are
  //min,max,rms,avg,sum, which uses MPI to perform those operations,
  //and root_singleton for non-communicated values).
  //     Then compute the variable in Grid_ComputeGlobalStatsGrid (see
  //that routine for details.)  That's it!  Fun and easy!
  //Variables get output in the order they're registered.


  //Please put a space after the format. It's for your own good.
  Stats.RegisterNewStat("Cycle","%4.0f ", root_singleton);

  //Direct set of these variables without checking the existance is dangerous-- don't do it.
  //I've checked this code, so it's ok here.
  Stats.elem("Cycle")->value = MetaData->CycleNumber;
  Stats.RegisterNewStat("Time", "%9.6e ", root_singleton);
  Stats.elem("Time")->value = MetaData->Time;  
  Stats.RegisterNewStat("Timestep", "%9.6e ", root_singleton);
  Stats.elem("Timestep")->value = LevelArray[level]->GridData->ReturnTimeStep();

  Stats.RegisterNewStat("V-rms","%9.6e ", rms); 
  Stats.RegisterNewStat("B-rms","%9.6e ", rms);

  Stats.RegisterNewStat("AvgBx","%9.6e ", avg);
  Stats.RegisterNewStat("AvgBy","%9.6e ", avg);
  Stats.RegisterNewStat("AvgBz","%9.6e ", avg);
  Stats.RegisterNewStat("MaxB^2","%9.6e ", max);
  Stats.RegisterNewStat("RMS-mach","%9.6e ", rms);
  Stats.RegisterNewStat("AvgKE","%9.6e ", avg);
  Stats.RegisterNewStat("MaxV^2","%9.6e ", max);
  Stats.RegisterNewStat("AvgVx","%9.6e ", avg);
  if( MetaData->TopGridRank>1 || MHD_Used )
    Stats.RegisterNewStat("AvgVy","%9.6e ", avg);
  if( MetaData->TopGridRank>2 || MHD_Used )
    Stats.RegisterNewStat("AvgVz","%9.6e ", avg);  


  Stats.RegisterNewStat("DensityVariance","%9.6e ", rms);
  Stats.RegisterNewStat("MaxD","%9.6e ", max);
  Stats.RegisterNewStat("MinD","%9.6e ", min);

  Stats.RegisterNewStat("RMS-AlvMach","%9.6e ", rms);
  Stats.RegisterNewStat("AvgBeta","%9.6e ", avg);

  //3.) loop over grids.
  LevelHierarchyEntry * Temp = LevelArray[level];
  while( Temp != NULL ){
    if( Temp->GridData->ComputeGlobalStatsGrid(&Stats) == NULL )
      {fprintf(stderr,"Error in ComputeGlobalStatsGrid\n"); return FAIL;}
    Temp = Temp->NextGridThisLevel;
  }

  //4.) Trade the answers (Sum the sums, min the Mins, max the Maxs.)
  Stats.Communicate();

  //5.) Apply all devides, sqrts, etc. that you need.
  Stats.FinalAnalysis( MetaData );

  /* This will become a loop over summed, rms quantities

  int numberOfGridZones = 1;
  
  for (int dim = 0; dim <  MetaData->TopGridRank; dim++)
  numberOfGridZones *= MetaData->TopGridDims[dim];
  
  float OneOverN = 1./numberOfGridZones;
  
  Stats.RMSVelocity = sqrt( Stats.RMSVelocity * OneOverN );
  Stats.RMSMagneticField = sqrt( Stats.RMSMagneticField * OneOverN );
  for( int dim=0; dim < 3; dim++ ){
  Stats.AvgMagneticField[dim] = Stats.AvgMagneticField[dim] * OneOverN;
  }
  Stats.RMSMachNumber = sqrt( Stats.RMSMachNumber * OneOverN );
  Stats.RMSAlvMach  = sqrt( Stats.RMSAlvMach * OneOverN );
  Stats.AvgBeta     = Stats.AvgBeta * OneOverN;
  //Stats.MassWeightedMach = sqrt( Stats.MassWeightedMach * OneOverN ); coming soon!
  Stats.DensityVariance = sqrt( Stats.DensityVariance * OneOverN );
  Stats.AvgKinetic = Stats.AvgKinetic * OneOverN;
  
  Stats.MaxSpeed = sqrt(Stats.MaxSpeed);
  Stats.MaxMagneticField = sqrt(  Stats.MaxMagneticField );
  */

  /*
    if (MetaData->CycleSkipGlobalDataDump != 0)
    if (MyProcessorNumber == ROOT_PROCESSOR &&
    MetaData->CycleNumber % MetaData->CycleSkipGlobalDataDump == 0.0) {
  */
  int WriteDataThisTime = TRUE;
  if( TRUE == WriteDataThisTime  ){
    
    //6.) Open file, append.
    FILE * Fptr;
    if ((Fptr = fopen("GlobalStats.out", "a")) == NULL)
      ERROR_MESSAGE;
    
    //7.) write data to file
    Stats.print( Fptr );
    
    //8.) close file.
    fclose( Fptr );
    
  }

  return SUCCESS;
}
#endif //PPML
