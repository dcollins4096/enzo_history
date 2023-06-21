#ifndef __HDFIO_HH_
#define __HDFIO_HH_
#include "Arch.h"
#include "IO.hh"
#include <hdf.h>

class HDFIO : public IObase {
  int32 fid,sid,current_rank;
  int32 DataType2HDF(IObase::DataType nt);
  IObase::DataType HDF2DataType(int32 nt);
  int chunkdims[5];
  int hasread;
  void create(int rank,CONST int *dims,DataType nt);
  void select(int i);
  void endaccess();
public:
  HDFIO(CONST char *fname,AccessMode access);
  virtual ~HDFIO();
  virtual int isValid();
  virtual int write(IObase::DataType typeID,int rank,CONST int *dims,void *data);
  virtual int readInfo(IObase::DataType &typeID,int &rank,int *dims,int maxdims=3);
  virtual int readInfo(char *name,IObase::DataType &typeID,int &rank,int *dims,int maxdims=3);
  virtual int read(void *data);
  virtual int seek(int i);
  virtual int nDatasets();
  virtual int writeAnnotation(CONST char *annotation);
  virtual int readAnnotationInfo(int number,int &length); // NEW
  virtual int readAnnotation(int index,char *annotation,int maxlen); // changed
  virtual int nAnnotations(); // New

  virtual int writeAttribute(CONST char *name,IObase::DataType typeID,Long length,void *data); // New
  virtual int readAttributeInfo(int number,char *name,IObase::DataType &typeID,Long &nelem,int maxnamelen=128); // New
  virtual int readAttributeInfo(CONST char *name,IObase::DataType &typeID,Long &nelem); // New
  virtual int readAttribute(int number,void *data); // New
  // virtual Long readAttribute(CONST char *name,void *data);
  virtual int nAttributes(); // New
  
  //-----------------Chunking Utilities..................
  virtual int reserveChunk(IObase::DataType typeID,int rank,CONST int *dims);
  virtual int writeChunk(CONST int *chunkdims,CONST int *chunkorigin,void *data);
  virtual int readChunk(CONST int *chunkdims,CONST int *chunkorigin,void *data);
  //-----------------Streaming interface (for PANDA, Sockets & MPIO).........
  virtual int reserveStream(IObase::DataType typeID,int rank,CONST int *dims) {return -1; }
  virtual int writeStream(void *data,int length){ return -1; } // not implemented yet
  virtual int readStream(void *data,int length) { return -1; } // not implemented yet
  //----------------Special HDF-specific methods.........
  int isCoord(); // Tells whether the dataset is a coordinate or a dataset
  int readDimInfo(int dimnumber,char *name, IObase::DataType &datatype, int &length);
  int readDim(int dimnumber, void *dim);
  int writeDim(int dimnumber,IObase::DataType datatype,int length,void *dim);
  int writeDimName(int dimnumber,CONST char *name);
};

#define f_hdf_open F77NAME(hdf_open_,hdf_open,HDF_OPEN)
#define f_hdf_openr F77NAME(hdf_openr_,hdf_openr,HDF_OPENR)
#define f_hdf_openw F77NAME(hdf_openw_,hdf_openw,HDF_OPENW)
#define f_hdf_opena F77NAME(hdf_opena_,hdf_opena,HDF_OPENA)

extern "C" {
#ifdef CRAY // Note: This isn't really implemented yet...
#include <fortran.h>
  Long8 f_hdf_open(_fcd fcfilename,_fcd accessname);
  Long8 f_hdf_openr(_fcd fcfilename);
  Long8 f_hdf_openw(_fcd fcfilename);
  Long8 f_hdf_opena(_fcd fcfilename);
#else
  Long8 f_hdf_open(char *file,char *access,int flen,int alen);
  Long8 f_hdf_openr(char *file,int flen);
  Long8 f_hdf_openw(char *file,int flen);
  Long8 f_hdf_opena(char *file,int flen);
#endif
#include "HDFIO.h"
}
/*
  long io_open_ieee_(constCONST char *file,constCONST char *access,int flen,int alen);
*/

#endif

