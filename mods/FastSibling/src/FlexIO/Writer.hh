#ifndef __UNIWRITER_HH_
#define __UNIWRITER_HH_

#include "IO.hh"

class Writer {
protected:
  IObase &file;
  int drank;
  //int *ddims;
  int ddims[5];
  IObase::DataType dtypeID;
  // double *ddelta,*dorigin,*dext;
  double ddelta[5],dorigin[5],dext[5];
  int extsize;
  virtual void writeBounds();
public:
  Writer(IObase &outfile);
  virtual ~Writer();
  virtual void setRank(int rank); // more efficient to set rank once for all
  virtual void setDims(int *dims);//only if rank has been set(should have flag)
  virtual void setDims(int rank,int *dims);
  virtual void setType(IObase::DataType typeID);
  virtual void setOrigin(double *origin);
  virtual void setDelta(double *delta);
  virtual void setParams(int rank,int *dims,IObase::DataType typeID,
		 double *origin,double *delta);
  virtual void write(void *data);
  virtual void reserveChunk();
  virtual void writeChunk(int *chunkdims,int *chunkorigin,void *data);
};

extern "C" {
#include "Writer.h"
}

#endif
