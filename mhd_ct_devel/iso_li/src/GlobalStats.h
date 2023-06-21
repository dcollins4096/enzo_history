
//
// class GlobalStats 
//
// Object to deal with global (and simple) statists for runs.  
// Only houses information for ease of manipulation.
#define MAX_STATS 45
#define NAME_SIZE 30
enum GlobalStats_Comm{undef,sum,rms,avg,min,max,root_singleton};

enum stat_id {Cycle, Time, Timestep, V_rms, B_rms, AvgBx, AvgBy, AvgBz,
	      MaxB2, RMS_mach, AvgKE, MaxV2, AvgVx, AvgVy, AvgVz,
	      DensityVariance, MaxD, MinD, RMSAlvMach, AvgBeta, C_rho,
	      C_ke, C_mach, C_px, C_py, C_pz, C_be, C_bx, C_by, C_bz, C_beta,MaxDivB,
	      AvgTE,AvgGE,AvgAlfv,AvgBE,px,py,pz,bulkX, bulkY,bulkZ};


struct TopGridData;
class Stat_Obj{
 public:
  Stat_Obj();
  Stat_Obj(char * name, char * output, GlobalStats_Comm comm);
  //~Stat_Obj();

  char Name[NAME_SIZE];
  float value;
  GlobalStats_Comm Comm;
  char OutputMethod[NAME_SIZE];

};

class GlobalStats{

 public:
  GlobalStats();
  ~GlobalStats();
  void Clear();   //Zeros data.
  int Communicate();
  int FinalAnalysis(TopGridData *MetaData);

  void print( FILE * fptr);


  //The array of objects
  int RegisterNewStat(stat_id nameID, char * name, char * output, GlobalStats_Comm comm);
  Stat_Obj *Stat[MAX_STATS];   
  
  int PrintedAtFirstStep;  //Flag to determine if the header line should be printed.

  //Indicies for each communicator.  Indicies will relate to the Stat array.
  int N_Avg, N_RMS, N_Sum, N_Max, N_Min; //For communication: number of each type of MPI call to make.
  Stat_Obj *Min_stats[MAX_STATS], *Max_stats[MAX_STATS],
    *Avg_stats[MAX_STATS], *RMS_stats[MAX_STATS], *Sum_stats[MAX_STATS];

  Stat_Obj * elem(stat_id nameID);

};

