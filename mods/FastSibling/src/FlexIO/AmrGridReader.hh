#ifndef __AMRGRIDREADER_HH_
#define __AMRGRIDREADER_HH_

#include <IO.hh>
#include "AmrGrid.h"

class AmrGridReader {
  IObase &file;
public:
  AmrGridReader(IObase &f):file(f){}
  // get specific grid
  AmrGrid *getGrid(AmrGrid &g,int index){
    if(file.seek(index)<index)
      return 0; // don't load past end
    getGridInfo(g,index);
    getGridData(g,index);
    return &g;
  }
  AmrGrid *getGrid(int index){
    AmrGrid *g=new AmrGrid;
    return this->getGrid(*g,index);
  }
  // Other stuff
  AmrGrid *getGridInfo(AmrGrid &g,int index);
  AmrGrid *getGridData(AmrGrid &g,int index);
};


#endif // __AMRGRIDREADER_HH_
