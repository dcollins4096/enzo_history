#ifndef __COMMANDLINE_HH_
#define __COMMANDLINE_HH_

#include <stdio.h>
#include <FlexArrayTmpl.H>

struct Commandline {
  char *programname;
  char *outfileprefix,*filename;
  char *path;
  typedef char String32[32];
  // JAD-related Options-----------------------------------------------
  struct JadOpt {
    int enabled;
    int dumplevel; // dumprate is proportional to this level timestep
    
    FlexArray<char*> fieldnames; // fields to dump
    // FlexArray<String32> particlenames; // particlefields to dump or *all*
    int startstep,endstep; // in terms of simulationsteps
    float starttime,endtime; // in terms of simulationtime
    char fileext[32];
    char fileprefix[128];
    char *filename; // construct full filename based on prefix etc..??
    enum FileType {HDF=1,HDF4=1,HDF5=2,Raw=3,IEEE=3,Default=3};
    FileType filetype;
    // path?
    JadOpt(){ setDefaults(); }
    ~JadOpt(){ 
      for(int i=0;i<fieldnames.getSize();i++) 
	if(fieldnames[i]) delete fieldnames[i];
      if(filename) delete filename;
    }
    void setDefaults();
    void print(FILE *fp=stderr);
    int parseOption(int argc,char *argv[],int &i);
  };
  JadOpt jad;

  struct StarOpt {
    int enabled;
    int dumplevel;
    // always dump all fields for stars
    int startstep,endstep;
    float starttime,endtime;
    char fileext[32];
    char fileprefix[128];
    char *filename;
    // char *filename; // construct full filename based on prefix etc..??
    enum FileType {HDF=1,HDF4=1,HDF5=2,Raw=3,IEEE=3,Default=3};
    FileType filetype;
    StarOpt(){setDefaults();}
    ~StarOpt(){if(filename) delete filename;}
    //~StarOpt(){ for(int i=0;i<fieldnames.getSize();i++) 
    //  if(fieldnames[i]) delete fieldnames[i];}
    // currently we only dump *all* fields for stars.
    void setDefaults();
    void print(FILE *fp=stderr);
    int parseOption(int argc,char *argv[],int &i);
  };
  StarOpt stars;
  
  struct DMOpt {
    int enabled;
    int dumplevel;
    // always dump all fields for stars
    int startstep,endstep;
    float starttime,endtime;
    char fileext[32];
    char fileprefix[128];
    char *filename;
    // char *filename; // construct full filename based on prefix etc..??
    enum FileType {HDF=1,HDF4=1,HDF5=2,Raw=3,IEEE=3,Default=3};
    FileType filetype;
    DMOpt(){ setDefaults(); }
    ~DMOpt(){if(filename) delete filename;}
    void setDefaults();
    void print(FILE *fp=stderr);
    int parseOption(int argc,char *argv[],int &i);
  };
  DMOpt dm;
  
  // more options for later--------------------------------------------
  int checkpoint;  // generate checkpoint files
  int checkpoint_sim_step_fre; // stepfreq for checkpoints
  float checkpoint_sim_time_freq; // checkpoint based on simtime
  float checkpoint_elapsed_time; // checkpoint after every elapsed time
  // unlike regular outputs, the checkpoint will output *all* fields
  
  
  // File stuff is going to require a global datastructure
  void setDefaults(){
    programname="";
    filename=0;
    outfileprefix=0;
    path = 0;
    stars.setDefaults();
    jad.setDefaults();
    dm.setDefaults();
  } 
  Commandline(){
    setDefaults();
  }
  Commandline(int argc,char *argv[]){
    setDefaults();
    parse(argc,argv);
  }
  ~Commandline(){
    if(path) delete path;
  }
  void print(FILE *fp=stderr){
    fprintf(fp,"Program %s: JAD/Star/DM options are\n",programname);
    jad.print(fp);
    stars.print(fp);
    dm.print(fp);
  }
  void parse(int argc,char *argv[]);
  int parseOption(int argc,char *argv[],int &i);
  void parseAndEat(int &argc,char *argv[]);
};

#endif
