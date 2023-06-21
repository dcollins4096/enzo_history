// AmrGridReader
#include "AmrGridReader.hh"
#include <stdio.h>
#include <stdlib.h>

// Other stuff
AmrGrid *AmrGridReader::getGridInfo(AmrGrid &g,int index){
  g.dataveclen=1;
  if(file.seek(index)<index)
    return 0; // fail if index past end
  IObase::DataType dt;
  file.readInfo(dt,g.rank,g.dims);
  g.datatype = dt;
  g.nbytes = IObase::nBytes(dt,g.rank,g.dims);
  // find the deepest level (finest time resolution)
  // Attrib Names?
  IObase::DataType atype;
  int length;
  int attrnum=file.readAttributeInfo("level",atype,length);
  if(attrnum>=0){
    int lev; // should be Int level
    file.readAttribute(attrnum,&lev);
    if(lev>g.maxlevel) g.maxlevel=lev;
    g.level=lev;
  }
  attrnum=file.readAttributeInfo("time_refinement",atype,length);
  if(attrnum>=0){
    file.readAttribute(attrnum,&(g.timerefinement));
  }
  attrnum=file.readAttributeInfo("timestep",atype,length);
  if(attrnum>=0){
    file.readAttribute(attrnum,&(g.timestep));
  }
  attrnum=file.readAttributeInfo("origin",atype,length);
  if(attrnum>=0)
    file.readAttribute(attrnum,(g.origin));
  attrnum=file.readAttributeInfo("delta",atype,length);
  if(attrnum>=0)
    file.readAttribute(attrnum,(g.delta));
  attrnum=file.readAttributeInfo("min_ext",atype,length);
  if(attrnum>=0)
    file.readAttribute(attrnum,(g.minext));
  attrnum=file.readAttributeInfo("max_ext",atype,length);
  if(attrnum>=0)
    file.readAttribute(attrnum,(g.maxext));
  attrnum=file.readAttributeInfo("persistence",atype,length);
  if(attrnum>=0){
    file.readAttribute(attrnum,&(g.persistence));
    g.maxtime = g.timestep + g.persistence;
  }
  g.data=0;
  return &g;
} // done 
  
AmrGrid *AmrGridReader::getGridData(AmrGrid &g,int index){
  IObase::DataType atype;
  // if(data) free(data); data=0; // make certain it is empty first
  g.data = malloc(g.nbytes);
  file.seek(index);
  file.readInfo(atype,g.rank,g.dims);
  g.datatype=atype;
  file.read(g.data);
  return &g;
}
