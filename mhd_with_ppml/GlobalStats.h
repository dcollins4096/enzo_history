
//
// class GlobalStats 
//
// Object to deal with global (and simple) statists for runs.  
// Only houses information for ease of manipulation.
#define MAX_STATS 30
#define NAME_SIZE 30
enum GlobalStats_Comm{undef,sum,rms,avg,min,max,root_singleton};
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
  int RegisterNewStat(char * name, char * output, GlobalStats_Comm comm);
  Stat_Obj *Stat[MAX_STATS];   
  
  int PrintedAtFirstStep;  //Flag to determine if the header line should be printed.

  //Indicies for each communicator.  Indicies will relate to the Stat array.
  int N_Registered,N_Avg, N_RMS, N_Sum, N_Max, N_Min; //For communication: number of each type of MPI call to make.
  Stat_Obj *Min_stats[MAX_STATS], *Max_stats[MAX_STATS],
    *Avg_stats[MAX_STATS], *RMS_stats[MAX_STATS], *Sum_stats[MAX_STATS];

  Stat_Obj * elem(char * name);
  /*
  //step and timestep 
  int CycleNumber;
  float Timestep, Time;
  
  //Summed quantities.  The pointers above will be used to loop over these quantities,
  //so Summed quantities must be contiguous in the class definition.
  //If adding fields, change the count in the constructor!
  float RMSVelocity;         // sqrt( <v^2> )
  float RMSMagneticField;    // sqrt( <B^2> ) Small scale field.
  float RMSAlvMach;          // sqrt( <v^2/(B^2/rho) > ) Avg Alfven Mach Number
  float AvgBeta;             // sqrt( <2*P/B^2 > ) Avg Plasma Beta
  float AvgKinetic;          // Or do I want avg kinetic energy?
  float AvgMagneticField[3]; // Average vector quantity
  float RMSMachNumber;       // sqrt( < v^2/c^2 > )
  float DensityVariance;     // ?

  //Max quantities.  Quantites over which the max will be taken.  Also must be contiguous.
  float MaxDensity;
  float MaxSpeed;
  float MaxMagneticField;

  //Min quantities.  Quantites over which the min will be taken.  Also must be contiguous.
  float MinDensity;
  */

};


