// AmrFileReader
#ifndef __AMRSTREAMINGFILEREADER_HH_
#define __AMRSTREAMINGFILEREADER_HH_
#include <stdio.h>
#include <IO.hh>
#include "AmrGridReader.hh"
#include "FlexArrayTmpl.H"

class AmrStreamingFileReader {
protected:
  int gridloading;
  AmrGridReader gridreader;
  FlexArray<int> activeGrids;
  FlexArray<AmrGrid> grids; // perhaps use a hashtable to find grids?
  // Change to AmrGridFile structures
  FlexArray<int> levelmask;
  IObase::DataType datatype;
  int maxlevel,maxtimeres,mintime,maxtime;
  int current_time;
  // Internal Utility methods
  void buildInfoTable();
  void loadGrids();
  void reclaimGrids();
  void purgeGrids();
  void printGridInfo(AmrGrid &g);
public:
  int debug;
  void printGridInfo();
  void printActiveGrids();
  AmrStreamingFileReader(IObase &f);
  int getNumLevels(){ return maxlevel+1; }
  void getTimeRange(int &min,int &max){
    min=mintime;
    max=maxtime;
  }
  void setTime(int timestep);
  // starts out with all selected
  void showLevel(int level=-1){ // default is all (-1)
    if(level>=levelmask.getSize() || level<0){
      printf("AmrConvert::showLevel(%u) :  Level out of range 0:%u\n",
	     level,levelmask.getSize()-1);
    }
    else
      levelmask[level]=1;
  }
  void showAllLevels();
  void hideLevel(int level=-1){  // default is all (-1)
    if(level>=levelmask.getSize() || level<0){
      printf("AmrConvert::showLevel(%u) :  Level out of range 0:%u\n",
	     level,levelmask.getSize()-1);
    }
    else
      levelmask[level]=0;
  }
  int nLevels(){ return maxlevel+1; }
  IObase::DataType getDataType(){return datatype;}
  // For C interface
  int getNumGrids(){ // number of active grids
    return activeGrids.getSize();
  }
  int getGrids(AmrGrid *g);
  // For C++ interface
  int getGrids(FlexArray<AmrGrid> &g);
  void setDataLoadingOff(){ gridloading=0; purgeGrids();}
  void setDataLoadingOn(){ gridloading=1; loadGrids();}
};

#endif
