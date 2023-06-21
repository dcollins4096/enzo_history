#include <stdio.h>
#include "Writer.hh"

void Writer::writeBounds(){
  if(!ddelta || !dorigin) return;
  /*
  if(extsize!=drank){
    if(dext) delete dext;
    dext=new double[drank];
  }*/
  for(int i=0;i<drank;i++) // always recompute extents
    dext[i]=dorigin[i]+ddelta[i]*(double)(ddims[i]);
  file.writeAttribute("origin",IObase::Float64,drank,dorigin);
  file.writeAttribute("delta",IObase::Float64,drank,ddelta);
  file.writeAttribute("min_ext",IObase::Float64,drank,dorigin);
  file.writeAttribute("max_ext",IObase::Float64,drank,dext);
}

Writer::Writer(IObase &outfile):
  file(outfile),dtypeID(IObase::Float32),drank(3),extsize(0) {
    // zero the values here...  Set to sensible initial values
    for(int i=0;i<5;i++){
      dorigin[i]=0;
      ddims[i]=0;
      dext[i]=1.0;
      ddelta[i]=1.0;
    }
}

Writer::~Writer() {
/*
  if(ddims) delete ddims;
  if(ddelta) delete ddelta;
  if(dorigin) delete dorigin;
  if(dext) delete dext;
*/
}

void Writer::setRank(int rank){
  drank=rank;
  /*
  if(ddims) delete ddims;
  ddims=new int[rank];
  */
}

void Writer::setDims(int *dims){
  for(int i=0;i<drank;i++) ddims[i]=dims[i];
}

void Writer::setDims(int rank,int *dims){
  setRank(rank);
  setDims(dims);
}

void Writer::setType(IObase::DataType typeID){
  dtypeID=typeID;
}

void Writer::setOrigin(double *origin){
/*
  if(dorigin) delete dorigin;
  dorigin = new double[drank]; */
  for(int i=0;i<drank;i++) dorigin[i]=origin[i];
}

void Writer::setDelta(double *delta){
/*
  if(ddelta) delete ddelta;
  ddelta = new double[drank]; */
  for(int i=0;i<drank;i++) ddelta[i]=delta[i];
}

void Writer::setParams(int rank,int *dims,IObase::DataType typeID,
			      double *origin,double *delta){
  setDims(rank,dims);
  setType(typeID);
  setOrigin(origin);
  setDelta(delta);
}

void Writer::write(void *data){
  file.write(dtypeID,drank,ddims,data);
  //puts("writebounds after data");
  writeBounds();
}

void Writer::reserveChunk(){
  file.reserveChunk(dtypeID,drank,ddims);
  writeBounds(); // write the boundary information
}

void Writer::writeChunk(int *dims,int *origin,void *data){
  // make certain the chunk is reserved? 
  file.writeChunk(dims,origin,data);
}
//===========C Interface======================================
// should have an RTTI interface and inherit everything from
// IOobject which contains the RTTI isOfType() information.
// isOfType() should propagate recursively to determin type info
// match.  must grab typeID from floating point pool.  And then
// we need a static initializer for everything.

// How does performer/inventor do this?
// if it hits an undefined interface (-1 default static value) then
// nobody has constructed an interface of that type yet.  So that would
// determine that the class is unrelated...

WRFile WRbeginFile(IOFile descriptor){
  IObase *io = (IObase*)descriptor;
  return (WRFile)(new Writer(*io));
}

void WRendFile(WRFile afile){
  Writer *w = (Writer*)afile;
  delete w;
}

void WRsetRank(WRFile afile,int rank){
  Writer *w = (Writer*)afile;
  w->setRank(rank);
}

void WRsetType(WRFile afile,int numbertype){
  Writer *w = (Writer*)afile;
  w->setType(IObase::Int2DataType(numbertype));
}

void WRsetParams(WRFile afile,
		 int rank,int *dims,int type,
		 double *origin,double *delta){	     
  Writer *w = (Writer*)afile;
  w->setParams(rank,dims,IObase::Int2DataType(type),origin,delta);
}

void WRsetDims(WRFile afile,int *dims){
  Writer *w = (Writer*)afile;
  w->setDims(dims);
}

void WRsetRankDims(WRFile afile,int rank, int *dims){
  Writer *w = (Writer*)afile;
  w->setDims(rank,dims);
}

void WRsetOrigin(WRFile afile,double *origin){
  Writer *w = (Writer*)afile;
  w->setOrigin(origin);
}

void WRsetDelta(WRFile afile,double *delta){
  Writer *w = (Writer*)afile;
  w->setOrigin(delta);
}

void WRwrite(WRFile afile,void *data){
  Writer *w = (Writer*)afile;
  w->write(data);
}

void WRwriteChunk(WRFile afile){
  Writer *w = (Writer*)afile;
  w->reserveChunk();
}

void WRwriteChunk(WRFile afile,int *dims,int *origin,void *data){
  Writer *w = (Writer*)afile;
  w->writeChunk(dims,origin,data);
}
