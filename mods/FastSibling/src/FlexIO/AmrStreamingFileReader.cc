#include <stdio.h>
#include <stdlib.h>
#include "AmrStreamingFileReader.hh"
#include "FlexArrayTmpl.H"
void AmrStreamingFileReader::printGridInfo(AmrGrid &g){
    printf("Grid level=%u step=%u maxtime=%u\n",g.level,g.timestep,g.maxtime);
    printf("\trank=%u dims[",g.rank);
    for(int i=0;i<g.rank-1;i++) printf("%u,",g.dims[i]);
    printf("%u]\n",g.dims[g.rank-1]);
    printf("\tTimerefine=%u dx[0]=%lf origin[",g.timerefinement,g.delta[0]);
    for(i=0;i<g.rank-1;i++) printf("%lf,",g.origin[i]);
    printf("%lf]\n",g.origin[g.rank-1]);
    printf("\tData Pointer is %u\n",(unsigned long)(g.data));
    if(g.data){
      float *fdata = (float *)(g.data);
      for(int i=0;i<3;i++)
	printf("\t\tData[%u]=%f\n",i,fdata[i]);
    }
  }
 void AmrStreamingFileReader::printGridInfo(){
    // print it all out...
    printf("MaxLevel=%u\n",maxlevel);
    for(int i=0;i<grids.getSize();i++){
      printf("Grid[%u]--------------------\n",i);
      printGridInfo(grids[i]);
    }
  }
  void AmrStreamingFileReader::printActiveGrids(){
    // print it all out...
    printf("Num Active Grids=%u\n",activeGrids.getSize());
    for(int i = 0; i < activeGrids.getSize(); i++){
      printf("  Grid[%u] mapped from %u --------------------\n",
	     i,activeGrids[i]);
      printGridInfo(grids[activeGrids[i]]);
    }
  }

void AmrStreamingFileReader::buildInfoTable(){
  // Load all grids Info
  int index=0;
  AmrGrid g;
  while(gridreader.getGridInfo(g,index++)){
    if(this->debug) printf("buildInfoTable: getGrid index=%u\n",index);
    if(this->debug) printGridInfo(g);
    int i=grids.getSize();
    g.data=0; // zero out the data
    grids.append(g);
    if(g.level>maxlevel)
      maxlevel=g.level;
    if(!i){
      mintime = g.timestep;
      maxtime = g.timestep;
      maxtimeres = g.timerefinement;
      maxlevel= g.level;
    }
    if(g.timestep<mintime) 
      mintime=g.timestep;
    if(g.timestep>maxtime) 
      maxtime=g.timestep;
    if(g.timerefinement>maxtimeres) 
      maxtimeres=g.timerefinement;
  }
}

void AmrStreamingFileReader::loadGrids(){
  if(!gridloading) return;
  for(int i=0;i<activeGrids.getSize();i++){
    if(this->debug) printf("buildInfoTable: getGrid index=%u activegridindex %u\n",i,activeGrids[i]);
    if(this->debug) printGridInfo(grids[activeGrids[i]]);
    gridreader.getGridData(grids[activeGrids[i]],activeGrids[i]);
  }
}

void AmrStreamingFileReader::reclaimGrids(){
  if(!gridloading) return;
  for(int i=0;i<grids.getSize();i++){
    int f=0;
    for(int j=0;j<activeGrids.getSize();j++){
      if(activeGrids[j]==i){
	f=1;
	break;
      }
    }
    if(!f){
      free((grids[i]).data);
      (grids[i]).data=0;
    }
  }
}

void AmrStreamingFileReader::purgeGrids(){
  for(int i=0;i<grids.getSize();i++){
    if((grids[i]).data)
      free((grids[i]).data);
    (grids[i]).data=0;
  }
}

AmrStreamingFileReader::AmrStreamingFileReader(IObase &f):gridreader(f),debug(0),gridloading(1){
  // we need to build a table from the file
  // then select grids based on the timestep
  buildInfoTable(); // initialize the convertor
  levelmask.setSize(maxlevel+1);
  for(int i=0;i<=maxlevel;i++)
    levelmask[i]=1; // all levels visible is default
  showAllLevels();
  setTime(mintime);
}

void AmrStreamingFileReader::setTime(int timestep){
  // Make Grid Selections
  if(timestep<mintime || timestep>maxtime){
    printf("timestep %u is out of range %u:%u\n",
	   timestep,mintime,maxtime);
    return;
  }
  activeGrids.purge();
  current_time=timestep;
  if(this->debug) printf("setTime(%u): mintime=%u maxtime=%u\n",current_time,mintime,maxtime);
  for(int i=0;i<grids.getSize();i++){
    if(this->debug) printf("\tgrids[%u].timestep=%u maxtime=%u\n",i,grids[i].timestep,
	   grids[i].maxtime);
    if(current_time>=grids[i].timestep && 
       current_time<grids[i].maxtime && levelmask[grids[i].level]){
      activeGrids.append(i);
      if(this->debug) printf("\t\tAppendGrid number %u\n",i);
    }
  }
  if(this->debug) puts("load grids");
  loadGrids();
  if(this->debug) puts("reclaim grids");
  reclaimGrids();
}
	
void AmrStreamingFileReader::showAllLevels(){
  for(int i=0;i<levelmask.getSize();i++) levelmask[i]=1;
}
	
// For C interface

int AmrStreamingFileReader::getGrids(AmrGrid *g){ // assumes number of 
  for(int i=0;i<grids.getSize();i++)
    g[i]=grids[activeGrids[i]];
  return activeGrids.getSize();
}

// For C++ interface
int AmrStreamingFileReader::getGrids(FlexArray<AmrGrid> &g){
  g.setSize(activeGrids.getSize());
  for(int i=0;i<g.getSize();i++)
    g[i]=grids[activeGrids[i]];
  return activeGrids.getSize();
}
