
//
// GlobalStats constructor and destructor and clear function.
// Clear is separate from the constructor; since the object is declared as static, the constructor
// is only called once.
//

#include <stdio.h>
#include <math.h>
#include <string.h>
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

int   CommunicationAllSumValues(float *Values, int Number);
float CommunicationMinValue(float Value);
float CommunicationMaxValue(float Value);

GlobalStats::~GlobalStats(){
  //actually, there's nothing to destruct.
}
GlobalStats::GlobalStats(){
    PrintedAtFirstStep = FALSE;  //Flag to determine if the header line should be printed.
    N_Avg= N_RMS= N_Max= N_Min=0, N_Sum =0;
    for(int i=0;i<MAX_STATS;i++){
      Stat[i] =
	Min_stats[i]=
	Max_stats[i]=
	Avg_stats[i]=
	RMS_stats[i]=
	Sum_stats[i]=NULL;
    }
}

Stat_Obj::Stat_Obj(){
  value = 0.0;
  sprintf(Name,"NoName");
  sprintf(OutputMethod,"%f");
  Comm = undef;
}
Stat_Obj::Stat_Obj(char * input_name, char * input_output_method, GlobalStats_Comm input_comm){
  strcpy(Name, input_name);
  strcpy(OutputMethod, input_output_method);
  switch( input_comm ){
  case min:
    value = huge_number;
    break;
  case max:
    value = -huge_number;
    break;
  default:
    value = 0.0;
    break;
  }
  Comm = input_comm;
}

void GlobalStats::Clear(){
  for(int i=0; i<MAX_STATS; i++){
    if( Stat[i] == NULL ) continue;
    switch( Stat[i]->Comm ){
    case min:
      Stat[i]->value = huge_number;
      break;
    case max:
      Stat[i]->value = -huge_number;
      break;
    default:
      Stat[i]->value = 0.0;
      break;
    }//switch
  }//loop
}

int GlobalStats::RegisterNewStat(stat_id nameID,char * name, char * output, GlobalStats_Comm comm){

  if( this->elem(nameID) == NULL){

    if( nameID > MAX_STATS ){
      fprintf(stderr, "GlobalStats: stat_id larger than MAX_STATS.  Please Fix\n");
      fprintf(stderr, "GlobalStats: stat_id larger than MAX_STATS.  Please Fix\n");
      fprintf(stderr, "GlobalStats: stat_id larger than MAX_STATS.  Please Fix\n");
      fprintf(stderr, "GlobalStats: stat_id larger than MAX_STATS.  Please Fix\n");
      fprintf(stderr, "GlobalStats: stat_id larger than MAX_STATS.  Please Fix\n");
      fprintf(stderr, "GlobalStats: stat_id larger than MAX_STATS.  Please Fix\n");
      fprintf(stderr, "GlobalStats: stat_id larger than MAX_STATS.  Please Fix\n");
      fprintf(stderr, "GlobalStats: stat_id larger than MAX_STATS.  Please Fix\n");
      return SUCCESS;
    }
    
    Stat_Obj * ThisStat  = new Stat_Obj(name, output, comm);  
    
    Stat[nameID] = ThisStat;

    switch( comm ){
    case rms:
      RMS_stats[N_RMS++] = ThisStat;
      break;
    case avg:
      Avg_stats[N_Avg++] = ThisStat;
      break;
    case min:
      Min_stats[N_Min++] = ThisStat;
      break;
    case max:
      Max_stats[N_Max++] = ThisStat;
      break;
    case sum:
      Sum_stats[N_Sum++] = ThisStat;
      break;
    case root_singleton:
      break;
    case undef:
      break;
    default:
      fprintf(stderr, "GlobalStats::RegisterNewStat: ERROR!!!  You didn't define a communicator for %s\n",name);
      
    }
  }
  return SUCCESS;
}
int GlobalStats::Communicate(){
  //Copy all "Summed" quantities (RMS, Avg, Sum) to an array, communicate, copy back.
  //Then do a global Min and Max for each element in the Min_stats, Max_stats array.
  int N, TotalSummed = 0;
  float * AllSummed = new float[ N_RMS + N_Sum + N_Avg ];
  
  for( N = 0; N<N_RMS; N++)
    AllSummed[TotalSummed++] = RMS_stats[N]->value;
  for( N = 0; N<N_Sum; N++)
    AllSummed[TotalSummed++] = Sum_stats[N]->value;
  for( N = 0; N<N_Avg; N++)
    AllSummed[TotalSummed++] = Avg_stats[N]->value;

  CommunicationAllSumValues(AllSummed, TotalSummed);
  TotalSummed = 0;
  for( N = 0; N<N_RMS; N++)
    RMS_stats[N]->value=AllSummed[TotalSummed++];
  for( N = 0; N<N_Sum; N++)
    Sum_stats[N]->value=AllSummed[TotalSummed++];
  for( N = 0; N<N_Avg; N++)
    Avg_stats[N]->value=AllSummed[TotalSummed++];
  
  for( int N=0; N<N_Max; N++)
    Max_stats[N]->value = CommunicationMaxValue( Max_stats[N]->value );
  for( int N=0; N<N_Min; N++)
    Min_stats[N]->value = CommunicationMinValue( Min_stats[N]->value );

  return SUCCESS;
}
int GlobalStats::FinalAnalysis(TopGridData *MetaData){
  int N;
  int numberOfGridZones = 1;

  for (int dim = 0; dim <  MetaData->TopGridRank; dim++)
    numberOfGridZones *= MetaData->TopGridDims[dim];

  float OneOverN = 1./numberOfGridZones;

  for( N = 0; N<N_RMS; N++)
    RMS_stats[N]->value = sqrt(RMS_stats[N]->value*OneOverN);
  for( N = 0; N<N_Avg; N++)
    Avg_stats[N]->value = Avg_stats[N]->value*OneOverN;
  return SUCCESS;
}
Stat_Obj * GlobalStats::elem(stat_id nameID){

  return Stat[nameID];

}
void GlobalStats::print(FILE * FPTR){
  int i;
  //The first time this is called, output the name of all registered stats.
  if( MyProcessorNumber == ROOT_PROCESSOR ){
    if( FALSE == PrintedAtFirstStep ){
      PrintedAtFirstStep = TRUE;
      //The # is for gnuplot.
      fprintf(FPTR,"#GLOBAL_STATS\n#");
      for( i=0; i<MAX_STATS; i++)
	if( Stat[i] )
	  fprintf(FPTR, "#%d.) %s \n",i+1,Stat[i]->Name);

    }
    for( i=0; i<MAX_STATS; i++){
      if( Stat[i] )
	fprintf(FPTR, Stat[i]->OutputMethod,Stat[i]->value);
    }
    fprintf(FPTR,"\n");    

    }//root proc

}
