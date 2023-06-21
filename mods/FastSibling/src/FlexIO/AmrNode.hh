#ifndef __AMRNODE_HH_
#define __AMRNODE_HH_

//#include "config.h"
#include "Bounds.hh"
#include <IO.hh>

struct sPoint { double x,y,z; };
union pPoint {
  double array[3];
  sPoint cartesian;
};

struct AmrNode {
  pPoint location;
  struct AmrNode *parent,*child;
  int index;
  int gridID;
  int level;
  float *data;
  /*
   IObase::DataType datatype;
   union {
     void *vdata;
     float *fdata;
     double *ddata;
     int *idata;
     char *cdata;
     short *sdata;
  }; */
  // DATATYPE data; 
  // my local ref to data... 
  // The grid is going to need the smarts to assign this...
  // So data segment length is important.  
  // Or we need a struct with a copy operator
  void setIndex(int i){
    index=i;
    if(parent) parent->setIndex(i);
  }
};

#endif
