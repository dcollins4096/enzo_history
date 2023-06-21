#ifndef __IEEEIO_HH_
#define __IEEEIO_HH_
#ifndef WIN32
#include <stdio.h>
#include <unistd.h>
#include <sys/file.h>
#else
#include <io.h>
#endif

#include "IO.hh"
#include "FlexArrayTmpl.H"
#include "Arch.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#if defined(WIN32) || defined(SOLARIS)
// I'm a WIN32 box and I ignore standards
// or I'm a Sun box and I hate BSD.
#define bcopy(src,dest,len) memcpy(dest,src,len)
#else
#include <strings.h> // available non-ANSI/Posix compliant 
#endif

#ifdef	FFIO
#include <ffio.h>
#endif

// Win32/Solaris lseek crap...
#if defined(WIN32) || defined(SOLARIS)
#define L_INCR SEEK_CUR
#ifndef L_SET // for that wanky Solaris compiler
#define L_SET  SEEK_SET
#endif
#define L_XTND SEEK_END
#endif
// T3E lseek crap... 
#ifdef T3E // I'm a T3E and I ignore ANSI/Posix standards
#define L_INCR SEEK_CUR
#define L_SET  SEEK_SET
#define L_XTND SEEK_END
#endif

#define IEEE_NATIVE_MAGIC 0x01020304
#define IEEE_REVERSE_MAGIC 0x04030201
#define IEEE_ALTERNATING_MAGIC 0x02010403

#define IEEE_BIG_ENDIAN IEEE_NATIVE_MAGIC
#define IEEE_LITTLE_ENDIAN IEEE_REVERSE_MAGIC
#define IEEE_VAX_ENDIAN IEEE_ALTERNATING_MAGIC

#define IEEE_MAJORVERSION 1
#define IEEE_MINORVERSION 0

class IEEEIO : public IObase {
public:
  enum RecordType {DataRecord=1,AnnotationRecord,AttributeRecord}; 
  union btorder {
    char c[4]; // write as Int, read as Int
    Int  i;
  };
#define FileHdrSize (4+4+4+2+2) // T3E Kludge
  struct FileHdr {
    // do magic twice so we can figure out which byte-order
    // the current machine is in addition (without hardcoding) in
    // addition to which byte order the file is in
    btorder magic;
    btorder byteorder;
    Long ndatasets;
    // T3E will be stupid for this, but thats OK
    short majorversion;
    short minorversion;
    // ** Inline these because they are 
    // ** actually kludges for the T3E
    static int write(IEEEIO *ieeeio,FileHdr &rec){
      int r=0;
#ifdef T3E // data structures are not packed contiguously in memory
      r+=ieeeio->write((rec.magic.c),4);
      r+=ieeeio->write((rec.byteorder.c),4);
      r+=ieeeio->write(&(rec.ndatasets),4);
      // this may be incorrect, but it will work for now
      // It assumes the low order bytes will be first to
      // be written.
      r+=ieeeio->write(&(rec.majorversion),2);
      r+=ieeeio->write(&(rec.minorversion),2);
#else
      // problem if this includes pointers to static methods
      // will require an extra degree of inheritance. (inherit
      // from a base structure.
      r=ieeeio->write(&rec,FileHdrSize);
#endif
      return r;
    }
    static int read(IEEEIO *ieeeio,FileHdr &rec){
      int r=0;
#ifdef T3E // data structures are not packed contiguously in memory
      char buffer[FileHdrSize];
      r=ieeeio->read(buffer,FileHdrSize); // read full buffer in and then copy to struct
      ::bcopy(buffer,(rec.magic.c),4);
      ::bcopy(buffer+4,(rec.byteorder.c),4);
      ::bcopy(buffer+8,&(rec.ndatasets),4);
      // May be invalid if low order bytes are not the
      // proper ones to load. (eg. may fail on big-endian systems)
      ::bcopy(buffer+12,&(rec.majorversion),2);
      ::bcopy(buffer+14,&(rec.minorversion),2);
#else
      r=ieeeio->read(&rec,FileHdrSize);
#endif
      return r;
    }
    static void byteswap(FileHdr &hdr){
      IEEEIO::byteswapBuffer(&(hdr.ndatasets),1,sizeof(hdr.ndatasets));
      // because the t3e doesn't have a 2-byte datatype, must do this manually
      IEEEIO::byteswapBuffer(&(hdr.majorversion),1,2);
      IEEEIO::byteswapBuffer(&(hdr.minorversion),1,2);
    }
    // May be a problem to have a nonstatic class member here because
    // this will create a vtable entry.
    void byteswap(){
      IEEEIO::byteswapBuffer(&ndatasets,1,sizeof(ndatasets));
      // because the t3e doesn't have a 2-byte datatype, must do this manually
      IEEEIO::byteswapBuffer(&majorversion,1,2);
      IEEEIO::byteswapBuffer(&minorversion,1,2);
    }
  };
#define RecordHdrSize (4+4+4) // T3E Kludge
  struct RecordHdr {
    Int recordtype;
    Long recordsize;
    Int sequenceID;
    // ** Inline these because they are 
    // ** actually kludges for the T3E
    static int write(IEEEIO *ieeeio,RecordHdr &rec){
      int r=0;
      // copy into buffer
      if(ieeeio->swapbytes){
	char buffer[RecordHdrSize];
	IEEEIO::byteswapBuffer(&(rec.recordtype),buffer,1,4);
	IEEEIO::byteswapBuffer(&(rec.recordsize),buffer+4,1,4);
	IEEEIO::byteswapBuffer(&(rec.sequenceID),buffer+8,1,4);
	r=ieeeio->write(buffer,RecordHdrSize);
      }
      else {
#ifdef T3E // data structures are not packed contiguously in memory
	r+=ieeeio->write(&(rec.recordtype),4);
	r+=ieeeio->write(&(rec.recordsize),4);
	r+=ieeeio->write(&(rec.sequenceID),4);
#else
	r=ieeeio->write(&rec,RecordHdrSize);
	//printf("Wrote rechdr %u,%u,%u\n",rec.recordtype,rec.recordsize,rec.sequenceID);
#endif
      }
      return r;
    }
    static int read(IEEEIO *ieeeio,RecordHdr &rec){
      char buffer[RecordHdrSize];
      int r=ieeeio->read(buffer,RecordHdrSize);
      if(ieeeio->swapbytes)
	IEEEIO::byteswapBuffer(buffer,3,4);//3 elements of 4 bytes 
#ifdef T3E // data structures are not packed contiguously in memory
      // do it the stupid way first 
      // would be smarter to pack buffers really
      ::bcopy(buffer,&(rec.recordtype),4);
      ::bcopy(buffer+4,&(rec.recordsize),4);
      ::bcopy(buffer+8,&(rec.sequenceID),4);
#else
      ::bcopy(buffer,(char*)&rec,RecordHdrSize);
#endif
      return r;
    }
  };
  
#define DataRecordHdrSize (4+4+4) // T3E Kludge   
  struct DataRecordHdr {
    Long datasize;
    Int numbertype;
    Int rank;
    // ** Inline these because they are 
    // ** actually kludges for the T3E
    static int write(IEEEIO *ieeeio,DataRecordHdr &rec){
      int r=0;
      // copy into buffer
      if(ieeeio->swapbytes){
	char buffer[DataRecordHdrSize];
	IEEEIO::byteswapBuffer(&(rec.datasize),buffer,1,4);
	IEEEIO::byteswapBuffer(&(rec.numbertype),buffer+4,1,4);
	IEEEIO::byteswapBuffer(&(rec.rank),buffer+8,1,4);
	r=ieeeio->write(buffer,DataRecordHdrSize);
      }
      else {
#ifdef T3E
	r+=ieeeio->write(&(rec.datasize),4);
	r+=ieeeio->write(&(rec.numbertype),4);
	r+=ieeeio->write(&(rec.rank),4);
#else
	r=ieeeio->write(&rec,DataRecordHdrSize);
#endif
      }
      return r;
    }
    static int read(IEEEIO *ieeeio,DataRecordHdr &rec){	
      char buffer[DataRecordHdrSize];
      int r=ieeeio->read(buffer,DataRecordHdrSize);
      if(ieeeio->swapbytes)
	IEEEIO::byteswapBuffer(buffer,3,4);//3 elements of 4 bytes 
#ifdef T3E
      // do it the stupid way first 
      // would be smarter to pack buffers really	
      ::bcopy(buffer,&(rec.datasize),4);
      ::bcopy(buffer+4,&(rec.numbertype),4);
      ::bcopy(buffer+8,&(rec.rank),4);
#else
      ::bcopy(buffer,(char*)&rec,DataRecordHdrSize);
#endif
      return r;
    }
  };
#define AttributeRecordHdrSize (4+4+4) // T3E Kludge  
  struct AttributeRecordHdr {
    Long datasize;
    Int numbertype;
    Int namesize; // attributes are named
    // ** Inline these because they are 
    // ** actually kludges for the T3E
    static int write(IEEEIO *ieeeio,AttributeRecordHdr &rec){	
      int r=0;
      // copy into buffer
      if(ieeeio->swapbytes){ // pack a buffer
	char buffer[AttributeRecordHdrSize];
	IEEEIO::byteswapBuffer(&(rec.datasize),buffer,1,4);
	IEEEIO::byteswapBuffer(&(rec.numbertype),buffer+4,1,4);
	IEEEIO::byteswapBuffer(&(rec.namesize),buffer+8,1,4);
	r=ieeeio->write(buffer,RecordHdrSize);
      }
      else {
#ifdef T3E
	r+=ieeeio->write(&(rec.datasize),4);
	r+=ieeeio->write(&(rec.numbertype),4);
	r+=ieeeio->write(&(rec.namesize),4);
#else
	//printf("AttribRec::writing ds,nt,ns %d,%d,%d\n",rec.datasize,rec.numbertype,rec.namesize);
	r = ieeeio->write(&rec,AttributeRecordHdrSize);
	//printf("Wrote attribhdr %d,%d,%d\n",rec.datasize,rec.numbertype,rec.namesize);
	//printf("\tr=%u\n",r);
#endif
      }
      return r;
    }
    static int read(IEEEIO *ieeeio,AttributeRecordHdr &rec){
      int r=0;
      char buffer[AttributeRecordHdrSize];
      r=ieeeio->read(buffer,AttributeRecordHdrSize);
      if(ieeeio->swapbytes)
	IEEEIO::byteswapBuffer(buffer,3,4);
#ifdef T3E
      ::bcopy(buffer,&(rec.datasize),4);
      ::bcopy(buffer+4,&(rec.numbertype),4);
      ::bcopy(buffer+8,&(rec.namesize),4);
#else
      ::bcopy(buffer,(char*)&rec,AttributeRecordHdrSize);
      //printf("AttribRec::reading ds,nt,ns %d,%d,%d\n",rec.datasize,rec.numbertype,rec.namesize);
#endif
      return r;
    }
  };

  struct RecRef {
    RecordHdr rec;
    long offset;
    //virtual ~RecRef();
    // ** Inline these because they are 
    // ** actually kludges for the T3E
  };
  
  struct AttribRef : public RecRef {
    //FlexArray<char> name;
    char name[64];
    AttribRef();
    AttribRef(AttribRef &src);
    //virtual ~AttribRef();
    AttribRef &operator=(AttribRef &src);
    // ** Inline these because they are 
    // ** actually kludges for the T3E
  };
  
  struct DataRef : public RecRef {
    FlexArray<RecRef> annotations;
    FlexArray<AttribRef> attributes;
    long end_offset;
    DataRef();
    //virtual ~DataRef();
    DataRef(DataRef &src);
    DataRef &operator=(DataRef &src);
  };
private:
  friend class AttribRef;
  friend class DataRef;
  friend class FileHdr;
  friend class RecordHdr;
  friend class DataRecordHdr;
  friend class AttributeRecordHdr;
  // convenience utilities
#ifdef T3E
  int ffioIsTrue();
#endif
  void initPosition();
  inline long getPosition(){return this->lseek(0,L_INCR);}
  void restart();
  int nextRecord();
  int nextDataRecord();
  inline long getLength(){ return file_length; }
  static void byteswapBuffer(void *buf,long nelements,int elementsize);
  static void byteswapBuffer(void *source,void *dest,long nelements,int elementsize);
  int writeFileHeader();
  int readFileHeader();
  void rebuildFileHeader();
  void appendRecordTable();
  void buildRecordTable();
  void clearChunk(int nbytes);
  void openFile(CONST char *fname,IObase::AccessMode access,int swbytes);
  // annotations don't need a header
  // for c, size includes null terminator.
  // must decrement size for f77
  //inline int write(int fid,void *data,int size);
protected:
  FlexArray<DataRef> rec_table;
  FlexArray<Int> chunkdims;
  FileHdr file_header;
  RecordHdr current_rec,current_data_rec;
  DataRecordHdr current_dat;
  long current_rec_offset; // offset past the record header
  long stream_offset; // kludge to allow index cache to be turned off
  long current_dat_offset; // to datarec
  int fid,datasetnumber,ndatasets,swapbytes,current_reserved_chunk;
  int hasread,streaming,cur_type_size;
  IEEEIO *masterfile;
  int writebuffersize,writebuffercount;
  char *writebuffer;
  long savedposition;
  long virtual_position,actual_position;
  long file_length;
protected:
  inline int write(const void *data,size_t size);
  inline void flush();
  inline int lseek(long len, int direction);
  inline int read(void *data,size_t size);
public:
  //------------------------core stuff.....................
  enum AccessMode {Read=0,Write=1,Create=1,Append=2,SharedRead=4};
  IObase::AccessMode mode(IEEEIO::AccessMode amode){
    switch(amode){
    case Read:
      return IObase::Read;
    case Write:
      return IObase::Write;
    case Append:
      return IObase::Append;
    default:
      return IObase::Read;
    }
  }
  // the IEEEIO:: commented out to satisfy Microsoft Compiler
  IEEEIO(CONST char *fname, /*IEEEIO::*/AccessMode access,int swbytes=0);
  IEEEIO(CONST char *fname,IObase::AccessMode access,int swbytes=0);
  IEEEIO(IEEEIO *file); // read-only dup an existing open file
  virtual ~IEEEIO();
  virtual int isValid();
  // could use overloading to differentiate type here... (but I'm going simple)
  virtual int write(IObase::DataType typeID,int rank,CONST int *dims,void *data);
  virtual int readInfo(IObase::DataType &typeID,int &rank,int *dims,int maxdims=3);
  virtual int read(void *data);
  //virtual int readChunk(int dims,int origin,int stride,void *data)=0;
  virtual int seek(int dataset_index);
  virtual int nDatasets();
  
  virtual int writeAnnotation(CONST char *annotation);
  virtual int readAnnotationInfo(int number,int &length); // returns length (-1 if none left)
  virtual int readAnnotation(int number,char *annotation,int maxsize=128);
  virtual int nAnnotations();

  virtual int writeAttribute(CONST char *name,IObase::DataType typeID,Long length,void *data);
  // returns number
  virtual int readAttributeInfo(int number,char *name,IObase::DataType &typeID,Long &nelem,int maxnamelen=128);
  virtual int readAttributeInfo(CONST char *name,IObase::DataType &typeID,Long &nelem); // returns number
  virtual int readAttribute(int number,void *data);
  // virtual Long readAttribute(CONST char *name,void *data);
  virtual int nAttributes();

  //-----------------Chunking Utilities..................
  virtual int reserveChunk(IObase::DataType typeID,int rank,CONST int *dims);
  virtual int writeChunk(CONST int *chunkdims,CONST int *chunkorigin,void *data);
  virtual int readChunk(CONST int *chunkdims,CONST int *chunkorigin,void *data);
  // Streaming interface is for support of PANDA, Sockets etc..
  virtual int reserveStream(IObase::DataType typeID,int rank,CONST int *dims);
  virtual int writeStream(void *data,int length);
  virtual int readStream(void *data,int length);
  virtual int pause();
  virtual int resume();
  void bufferOn(long size=-1);
  void bufferOff();
};
/*
int bcopy(char *src,char *dst,int len){
  // copying bytes on the DEC Alpha processor can be costly due
  // to the penalty of non-aligned memory accesses.  This tries
  // to align the copys.

  You can use the program copyperf to see how bcopy() improves performance
  You must build copyperf separately using
     gmake copyperf
  
  The following will test copy a 64meg buffer 1 time.
     copyperf 64m
  The following will test copy a 32k buffer 1000 times
     copyperf 32k 1000
  The performance differences are 4-to-1 or greater in favor of bcopy().
}
*/
int IEEEIO::write(const void *data,size_t size){
  //puts("***********write()");
#ifdef T3E
  static int ffio_inline_true=0;
  if(ffio_inline_true >= 0){
    // do this test only once on first write
    // this just gives someone a clue that they've done the
    // wrong thing when they compiled
#ifdef FFIO
    if(!ffioIsTrue()){
      fprintf(stderr,"Error!!! Inconsistent FFIO flags on a T3E!\n");
      fprintf(stderr,"Your libraries have been compiled without -DFFIO, but your executable has been compled and linked with them with -DFFIO.\n");
      fprintf(stderr,"\t*****Please either rebuild the IEEEIO library with -DFFIO or *remove* -DFFIO from your own build.******\n");
    }
#else
    if(ffioIsTrue()){
      fprintf(stderr,"Error!!! Inconsistent FFIO flags on a T3E!\n");
      fprintf(stderr,"Your libraries have been compiled with -DFFIO, but your executable has been compled and linked with them without -DFFIO\n");
      fprintf(stderr,"\t*****Please recompile with -DFFIO******\n");
    }
#endif
    ffio_inline_true = -1;
  }
#endif
  if(!writebuffer){
    register long r=0;
#ifdef FFIO
    r = ::ffwrite(fid,data,size);
#else
    r = ::write(fid,data,size);
#endif
    if(r>0) actual_position += r;
    virtual_position = actual_position;
    // if(r>0) virtual_position = actual_position = (actual_position+r);
    if(virtual_position>file_length) file_length = virtual_position;
    return r;
  }
  // otherwise, do a buffered write
  int retval=0;
  char *datap = (char *)data;
  if(writebuffercount>0){
    // copy as much as you can to the writebuffer
    // should use bcopy... its a whole lot faster!!!
    long copysize = (size>(writebuffersize-writebuffercount))?(writebuffersize-writebuffercount):size;
    bcopy(datap,writebuffer+writebuffercount,copysize);
    writebuffercount += copysize;
    virtual_position += copysize;
    retval += copysize;
    if(writebuffercount>=writebuffersize)
      this->flush(); // dump immediately (updates both virtual and actual pos)
    size-=copysize; datap+=copysize;
    if(size<=0) {
      // virtual_position = actual_position + writebuffercount;
      if(virtual_position>file_length) file_length=virtual_position;
      return retval;
    }
  }
  // at this point size>0 and writebuffercount=0
  if(size<=0 || writebuffercount!=0){
    printf("***Write Failed: size=%u writebuffercount=%u\n",
	   (unsigned int)size,writebuffercount);
  }
  if(size<writebuffersize){// store to buffer (size>0 is implicit in loop)
    bcopy(datap,writebuffer+writebuffercount,size);
    writebuffercount += size;
    virtual_position += size;
    retval += size;
  }
  else {
    register long r=0;
    // its larger than buffer, so bypass it and
    // write it directly to the disk.
#ifdef FFIO
    r = ::ffwrite(fid,datap,size);
#else
    r = ::write(fid,datap,size);
#endif
    actual_position += r;
    virtual_position = actual_position; // essentially equiv to flush
    if(r>0) retval+=r;
  }
  // double-check virtual position
  // virtual_position = actual_position + writebuffercount;
  if(virtual_position>file_length) file_length = virtual_position;
  return retval;
}
void IEEEIO::flush(){
  register long r=0;
  if(writebuffercount){
#ifdef  FFIO
    r = ::ffwrite(fid,writebuffer,writebuffercount);
#else
    r = ::write(fid,writebuffer,writebuffercount);
#endif
  }
  writebuffercount=0;
  if((r+actual_position)!=virtual_position){
    fprintf(stderr,"IEEEIO::flush() : inconsistent file positions! r=%ld v=%ld a=%ld\n",
	 r,virtual_position,actual_position);
  }
  actual_position = virtual_position;
  if(virtual_position>file_length) file_length=virtual_position;
}
int IEEEIO::lseek(long len, int direction){
  // compute new pos
  register long npos=0;
  // can seek beyond end of file
  switch(direction){
  case L_XTND: /* should do special case for 0 len */
    npos = file_length + len;
    break;
  case L_SET:
    npos = len;
    break;
  case L_INCR:;
    npos = virtual_position + len;
    break;
  }
  if(npos>file_length)
    file_length = npos;
  if(virtual_position != npos) {
    if(writebuffer)
      this->flush(); // must flush current buffer contents
    // printf("Actually seeking from %ld to %ld\n",virtual_position,npos);
    virtual_position = actual_position = npos;
#ifdef  FFIO
    // if(!writebuffer) 
    return ::ffseek(fid,npos,L_SET);
#else
    // if(!writebuffer) 
    return ::lseek(fid,npos,L_SET);
#endif    
  }
  return npos; // otherwise no change required
}
/* Will eventually want to buffer the reads as well... */
int IEEEIO::read(void *data,size_t size){
  register int r=0;
  if(writebuffercount>0) 
    flush();
#ifdef  FFIO
  r = ::ffread(fid,data,size);
#else
  r = ::read(fid,data,size);
#endif
  if(r>0) virtual_position += r;
  actual_position=virtual_position;
  // can't read past end of file of course so no need to update file_length;
  return r;
}


//==========F77 Interface=============
#define f_ieee_open F77NAME(ieee_open_,ieee_open,IEEE_OPEN)
#define f_ieee_openr F77NAME(ieee_openr_,ieee_openr,IEEE_OPENR)
#define f_ieee_openw F77NAME(ieee_openw_,ieee_openw,IEEE_OPENW)
#define f_ieee_opena F77NAME(ieee_opena_,ieee_opena,IEEE_OPENA)
#define f_ieee_bufon F77NAME(ieee_bufon_,ieee_bufon,IEEE_BUFON)
#define f_ieee_bufoff F77NAME(ieee_bufoff_,ieee_bufoff,IEEE_BUFOFF)
extern "C"{
#ifdef CRAY // Note: This isn't really implemented yet...
#include <fortran.h>
  Long8 f_ieee_open (_fcd fcfilename,_fcd fcaccessmode);
  Long8 f_ieee_openr (_fcd fcfilename);
  Long8 f_ieee_openw (_fcd fcfilename);
  Long8 f_ieee_opena (_fcd fcfilename);
#else
  Long8 f_ieee_open (char *filename,char *accessmode,int namelen,int accesslength);
  Long8 f_ieee_openr (char *filename,int namelen);
  Long8 f_ieee_openw (char *filename,int namelen); // actually IObase::Create
  Long8 f_ieee_opena (char *filename,int namelen);
#endif
  void f_ieee_bufon (Long8 *fileID,int *bufsize);
  void f_ieee_bufoff (Long8 *fileID);
  //=====================ANSI C interface
#include "IEEEIO.h" // ANSI C interface
	  }

#endif
