#ifndef __IO_HH_
#define __IO_HH_
#include "Arch.h"

class IObase {
public:
  //-----------------------Public Enums....................
  enum AccessMode {Read=0,Write=1,Create=1,Append=2};
  enum DataType {Byte=0,Int8=0,Int16=1,Int32=2,Int64=3,
		 Float32=4,Float64=5,
		 uInt8=6,uChar=6,uInt16=7,uInt32=8,uInt64=9,
		 Char=10,Char8=10,String=10,Unicode=11,Char16=11, // special string types
		 Error=-1};
protected:
  //-----------------------State variables.................
  int index,nitems;
  AccessMode accessmode;
  char filename[128];
public:
  //------------------------core stuff.....................
  IObase(CONST char *fname,AccessMode access);
  virtual ~IObase() {}
  virtual int isValid() { return 0; }
  // could use overloading to differentiate type here... (but I'm going simple)
  virtual int write(DataType typeID,int rank,CONST int *dims,void *data)=0;
  virtual int readInfo(DataType &typeID,int &rank,int *dims,int maxdims=3)=0;
  virtual int read(void *data)=0;
  virtual int seek(int dataset_index)=0;
  virtual int nDatasets()=0;
  virtual int writeAnnotation(CONST char *annotation)=0;
  virtual int readAnnotationInfo(int number,int &length)=0; // returns length (-1 if none left)
  virtual int readAnnotation(int number,char *annotation,int maxsize=128)=0;
  virtual int nAnnotations()=0;
  
  virtual int writeAttribute(CONST char *name,IObase::DataType typeID,Long length,void *data)=0;
  // returns number
  virtual int readAttributeInfo(int number,char *name,IObase::DataType &typeID,Long &nelem,int maxnamelen=128)=0;
  virtual int readAttributeInfo(CONST char *name,IObase::DataType &typeID,Long &nelem)=0; // returns number
  virtual int readAttribute(int number,void *data)=0;
  // virtual Long readAttribute(char *name,void *data);
  virtual int nAttributes()=0;

  //-----------------Chunking Features (for MPI/HPF)................
  virtual int reserveChunk(IObase::DataType typeID,int rank,CONST int *dims)=0;
  
  virtual int writeChunk(CONST int *chunkdims,CONST int *chunkorigin,void *data)=0;
  // virtual int writeStridedChunk()=0;
  // virtual int writeHPF(int processorID,int *proclayout,IObase::HPFlayout *arraylayout,void *data)=0;
  virtual int readChunk(CONST int *chunkdims,CONST int *chunkorigin,void *data)=0;
  // virtual int readStrided()=0;
  // virtual int readHPF(int processorID,int *proclayout,IObase::HPFlayout *arraylayout,void *data)=0;
  //-----------------Streaming interface (for PANDA, Sockets & MPIO).........
  virtual int reserveStream(IObase::DataType typeID,int rank,CONST int *dims)=0;
  virtual int writeStream(void *data,int length)=0;
  virtual int readStream(void *data,int length)=0;
  //-----------------------Utilities........................
  // unfortunately you can cast enums to int but not the reverse
  static DataType Int2DataType(int dt);
  static int sizeOf(DataType dt);
  static int nElements(int rank,CONST int *dims); // returns number of elements in dataset
  static int nBytes(DataType dt,int rank,CONST int *dims); // returns number of bytes (good for malloc)
  CONST char *name() { return filename; }
  virtual int pause()   { return 0; }
  virtual int resume() { return 0; }
};

//===================Fortran77 Interface 
#define f_io_close F77NAME(io_close_,io_close,IO_CLOSE)
#define f_io_isvalid F77NAME(io_isvalid_,io_isvalid,IO_ISVALID)
#define f_io_sizeof F77NAME(io_sizeof_,io_sizeof,IO_SIZEOF)
#define f_io_nbytes F77NAME(io_nbytes_,io_nbytes,IO_NBYTES)
#define f_io_nelements F77NAME(io_nelements_,io_nelements,IO_NELEMENTS)
#define f_io_write F77NAME(io_write_,io_write,IO_WRITE)
#define f_io_readinfo F77NAME(io_readinfo_,io_readinfo,IO_READINFO)
#define f_io_read F77NAME(io_read_,io_read,IO_READ)
#define f_io_seek F77NAME(io_seek_,io_seek,IO_SEEK)
#define f_io_numdata F77NAME(io_numdata_,io_numdata,IO_NUMDATA)
#define f_io_writenote F77NAME(io_writenote_,io_writenote,IO_WRITENOTE)
#define f_io_readnote F77NAME(io_readnote_,io_readnote,IO_READNOTE)
#define f_io_numnote F77NAME(io_numnote_,io_numnote,IO_NUMNOTE)
#define f_io_writeatt F77NAME(io_writeatt_,io_writeatt,IO_WRITEATT)
#define f_io_attinfo F77NAME(io_attinfo_,io_attinfo,IO_ATTINFO)
#define f_io_iattinfo F77NAME(io_iattinfo_,io_iattinfo,IO_IATTINFO)
#define f_io_readatt F77NAME(io_readatt_,io_readatt,IO_READATT)
#define f_io_numatt F77NAME(io_numatt_,io_numatt,IO_NUMATT)
#define f_io_reserveck F77NAME(io_reserveck_,io_reserveck,IO_RESERVECK)
#define f_io_writeck F77NAME(io_writeck_,io_writeck,IO_WRITECK)
#define f_io_readck F77NAME(io_readck_,io_readck,IO_READCK)
#define f_io_writestrm F77NAME(io_writestrm_,io_writestrm,IO_WRITESTRM)
#define f_io_readstrm F77NAME(io_readstrm_,io_readstrm,IO_READSTRM)
#define f_io_pause F77NAME(io_pause_,io_pause,IO_PAUSE)
#define f_io_resume F77NAME(io_resume_,io_resume,IO_RESUME)
extern "C"{ 
//================Ansi C interface
#include "IO.h"
//==== f77 interface
#ifdef CRAY
  // Cray fortran uses an _fcd data structure to pass strings.  This is in preparation
  // for using this on Crays, but since nobody is using them yet, I've left the code for
  // these versions of the routines blank.  It'll be filled in later if the need arises
#include <fortran.h>
  int f_io_close (Long8 *deviceID);
  int f_io_isvalid (Long8 *deviceID);
  int f_io_sizeof (int *datatype);
  int f_io_nelements (int *rank,int *dims);
  int f_io_nbytes (int *datatype,int *rank,int *dims);
  int f_io_write (Long8 *deviceID,int *typeID,int *rank,int *dims,void *data);
  int f_io_readinfo (Long8 *deviceID,int *typeID,int *rank,int *dims,int *maxdims);
  int f_io_read (Long8 *deviceID,void *data);
  int f_io_seek (Long8 *deviceID,int *dataset_index);
  int f_io_numdata (Long8 *deviceID);
  int f_io_writenote (Long8 *deviceID,_fcd fcannotation);
  int f_io_noteinfo (Long8 *deviceID,int *number,int *length);
  int f_io_readnote (Long8 *deviceID,int *number,_fcd fcannotation);
  int f_io_numnote (Long8 *deviceID);

  int f_io_writeatt (Long8 *deviceID,_fcd fcname,
		   int *datatype,Long *nelements,void *data);
  int f_io_attinfo (Long8 *deviceID,_fcd fcname,
		  int *datatype,Long *nelements);
  int f_io_iattinfo (Long8 *deviceID,int *index,_fcd fcname,
		  int *datatype,Long *nelements);
  int f_io_readatt (Long8 *deviceID,int *number,void *data);
  int f_io_numatt (Long8 *deviceID);
  int f_io_reserveck(Long8 *deviceID,int *typeID,int *rank,int *dims);
  int f_io_writeck(Long8 *deviceID,int *chunkdims,int *chunkorigin,void *data);
  int f_io_readck(Long8 *deviceID,int *chunkdims,int *chunkorigin,void *data);
  int f_io_pause(Long8 *deviceID);
  int f_io_resume(Long8 *deviceID);
#else
  int f_io_close (Long8 *deviceID);
  int f_io_isvalid (Long8 *deviceID);
  int f_io_sizeof (int *datatype);
  int f_io_nelements (int *rank,int *dims);
  int f_io_nbytes (int *datatype,int *rank,int *dims);
  int f_io_write (Long8 *deviceID,int *typeID,int *rank,int *dims,void *data);
  int f_io_readinfo (Long8 *deviceID,int *typeID,int *rank,int *dims,int *maxdims);
  int f_io_read (Long8 *deviceID,void *data);
  int f_io_seek (Long8 *deviceID,int *dataset_index);
  int f_io_numdata (Long8 *deviceID);

  int f_io_writenote (Long8 *deviceID,char *annotation,int size);
  int f_io_noteinfo (Long8 *deviceID,int *number,int length);
  int f_io_readnote (Long8 *deviceID,int *number,char *annotation,int maxsize);
  int f_io_numnote (Long8 *deviceID);

  int f_io_writeatt (Long8 *deviceID,char *name,
		   int *datatype,Long *nelements,void *data,int namesize);
  int f_io_attinfo (Long8 *deviceID,char *name,
		  int *datatype,Long *nelements,int namesize);
  int f_io_iattinfo (Long8 *deviceID,int *index,char *name,
		  int *datatype,Long *nelements,int namesize);
  int f_io_readatt (Long8 *deviceID,int *number,void *data);
  int f_io_numatt (Long8 *deviceID);
  int f_io_reserveck(Long8 *deviceID,int *typeID,int *rank,int *dims);
  int f_io_writeck(Long8 *deviceID,int *chunkdims,int *chunkorigin,void *data);
  int f_io_readck(Long8 *deviceID,int *chunkdims,int *chunkorigin,void *data);
  int f_io_writestrm(Long8 *deviceID,void *data,int *length);
  int f_io_readstrm(Long8 *deviceID,void *data,int *length);
  int f_io_pause(Long8 *deviceID);
  int f_io_resume(Long8 *deviceID);
#endif
}

#endif



