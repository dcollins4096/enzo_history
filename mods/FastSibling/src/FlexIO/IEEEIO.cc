#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "IEEEIO.hh"

#ifdef WIN32
// Are we Microsoft VC++ 5.0 or 6.0?
#if defined(_MSC_VER) && ( (_MSC_VER == 1100) || (_MSC_VER == 1200) ) // yes we are
#include <fcntl.h>
#define O_RDONLY _O_RDONLY|_O_BINARY
#define O_RDWR _O_RDWR|_O_BINARY
#define O_WRONLY _O_WRONLY|_O_BINARY
#define O_CREAT _O_CREAT|_O_BINARY
#define O_TRUNC _O_TRUNC|_O_BINARY
#define O_APPEND _O_APPEND|_O_BINARY
#else // not an MSC compiler (use the old defines)
#define O_RDONLY _O_RDONLY
#define O_RDWR _O_RDWR
#define O_WRONLY _O_WRONLY
#define O_CREAT _O_CREAT
#define O_TRUNC _O_TRUNC
#define O_APPEND _O_APPEND
#endif // not an MSC compiler
#endif // WIN32

#ifdef FFIO
#include <ffio.h>
#define open(x,y,z) ffopens(x,y,z, 0, 0, NULL, "bufa.bufsize=256.num_buffers=4")
#define close(x) ffclose(x)
#endif

#ifdef T3E
int IEEEIO::ffioIsTrue(){
#ifdef FFIO
  return 1;
#else
  return 0;
#endif
}
#endif


// int IEEEIO::nextRecord() (stores state in current_record)
// also garauntees integrity of datasetnumber and 
// the current_rec and current_dat

IEEEIO::AttribRef::AttribRef(IEEEIO::AttribRef &src){ 
  //puts("AttribRef::copy constructor"); 
  rec=src.rec;
  offset=src.offset;
  strcpy(name,src.name);
}

IEEEIO::AttribRef::AttribRef(){
  //puts("empty attribref constructor");
}
/*
IEEEIO::AttribRef::~AttribRef(){
  //puts("delete an attrib ref");
}
IEEEIO::RecRef::~RecRef(){
  //puts("delete a recref");
}
*/

IEEEIO::AttribRef &IEEEIO::AttribRef::operator=(IEEEIO::AttribRef &src){
  if(this != &src) {
    //puts("\tAttribRef::copy operator");
    rec=src.rec;
    offset=src.offset;
    //name=src.name; // uses the flexarray copy constructor
    strcpy(name,src.name);
  }
  return *this;
}

IEEEIO::DataRef::DataRef():end_offset(0){
  //puts("Empty DataRef Constructor!!");
}


IEEEIO::DataRef::DataRef(IEEEIO::DataRef &src){
  //puts("DataRef:: copy constructor");
  rec=src.rec;
  offset=src.offset;
  end_offset=src.end_offset;
  annotations = src.annotations;
  attributes = src.attributes; 
  //puts("DataRef:: copy constructor DONE");
}
/*
IEEEIO::DataRef::~DataRef(){
  // puts("dataref destructor");
}*/

IEEEIO::DataRef &IEEEIO::DataRef::operator=(IEEEIO::DataRef &src){
  if(this != &src){
    rec=src.rec;
    offset=src.offset;
    end_offset=src.end_offset;
    annotations = src.annotations;
    attributes = src.attributes;
  }
  return *this;
}

void IEEEIO::initPosition(){
  // init
#ifdef FFIO
  actual_position = ::ffseek(fid,0,L_INCR);
#else 
  actual_position = ::lseek(fid,0,L_INCR);
#endif
  virtual_position = actual_position + writebuffercount;
  // get actual filesize
  
#ifdef FFIO
  file_length = ::ffseek(fid,0,L_XTND);
  ::ffseek(fid,actual_position,L_SET); // seek back to original pos
#else
  file_length = ::lseek(fid,0,L_XTND);
  ::lseek(fid,actual_position,L_SET); // seek back to original pos
#endif
  if(virtual_position > file_length) 
    file_length = virtual_position;
  // printf("Init position: a=%ld v=%ld f=%ld\n",
  //	 actual_position,virtual_position,file_length);
}


//long IEEEIO::getPosition(int immediate) {
//#ifdef FFIO
//  if(!writebuffer || immediate)
//    return ::ffseek(fid,0,L_INCR);
//  else 
//    return ::ffseek(fid,0,L_INCR)+writebuffercount;
//#else
//  if(!writebuffer || immediate)
//    return ::lseek(fid,0,L_INCR);
//  else 
//    return ::lseek(fid,0,L_INCR)+writebuffercount;
//#endif
//}

void IEEEIO::restart(){
  // this->initPosition();
  //puts("restart seek");
  IEEEIO::lseek(FileHdrSize,L_SET); // T3E Kludge
  //puts("restart");
  current_rec_offset=current_dat_offset=FileHdrSize; // T3E Kludge
  current_rec.recordsize=0;
  current_dat.datasize=0;
  current_dat.rank=0;
  datasetnumber=0; // initial
  hasread=0; // and it hasn't been read yet...
}

int IEEEIO::nextRecord(){
  RecordHdr rec;
  // should have it check the datasetnumber before seeking
  // this will require consistancy checks for the datasetnumber
  // elsewhere in the code
  long src=current_rec_offset;
  long dst=src+current_rec.recordsize;
  current_rec_offset=IEEEIO::lseek(dst,L_SET); // seek from start
  // or use rec.read(rec,fid)
  // or IEEEIO::RecordHdr::read(rec,fid)
  if(RecordHdr::read(this,rec) <= 0)
    return 0;
  current_rec_offset=dst+RecordHdrSize;
  current_rec=rec;
  return 1;
}

int IEEEIO::nextDataRecord(){
  hasread=0; // reset the hasread field
  if(current_rec.recordtype==DataRecord && 
     datasetnumber!=current_rec.sequenceID){
    IEEEIO::lseek(current_rec_offset,L_SET); // seek to record start
  }
  else {
    do { // scan through records for the next DataRecord
      if(!nextRecord()) 
	return 0; // reached end of file
    } while(current_rec.recordtype!=DataRecord);
  }
  // dat_offset is same as current_rec_offset for the data record
  // so nextRecord() can be used for scanning for annotations
  current_dat_offset=current_rec_offset;
  DataRecordHdr::read(this,current_dat);
  datasetnumber=current_rec.sequenceID;
  return 1;
}

//long IEEEIO::getLength(int immediate){
//#ifdef FFIO
//  long currentpos=::ffseek(fid,0,L_INCR); // find currentposition
//  long auspos = ::ffseek(fid,0,L_XTND); // seek to end to find length
//#else
//  long currentpos=::lseek(fid,0,L_INCR); // find currentposition
//  long auspos = ::lseek(fid,0,L_XTND); // seek to end to find length
//#endif
//  if(writebuffer && writebuffercount && !immediate){
//    // must compute real end of file
//    register long pos1=currentpos+writebuffercount;
//    if(pos1>auspos) auspos=pos1;
//  }
//#ifdef FFIO
//  ::ffseek(fid,currentpos,L_SET); // seek back to original position
//#else
//  ::lseek(fid,currentpos,L_SET); // seek back to original position
//#endif
//  return auspos; // return end position
//}

void IEEEIO::byteswapBuffer(void *buf,long nelements,int elementsize){
  char *buffer=(char *)buf; // treat as a character buffer
  if(elementsize<=1) return;
  for(long i=0;i<nelements;i++,buffer+=elementsize){
    register int s,d;
    register char c;
    // do the swap thang on each element
    for(s=0,d=elementsize-1;s<d;s++,d--){
      c=buffer[s];
      buffer[s]=buffer[d];
      buffer[d]=c;
    }
  }
}

void IEEEIO::byteswapBuffer(void *source,void *dest,long nelements,int elementsize){
  if(elementsize<=1) return;
  // Lets optimize for integers
  switch(elementsize){
  case 4:
    {
      Int *src=(Int*)source;
      Int *dst=(Int*)dest;
      for(long i=0;i<nelements;i++){
	union IntChar {
	  int i;
	  char a[4];
	};
	register IntChar ts,td;
	// do the swap thang on each element
	td.i=ts.i=src[i];
	td.a[0]=ts.a[3];
	td.a[1]=ts.a[2];
	td.a[2]=ts.a[1];
	td.a[3]=ts.a[0];
	dst[i]=td.i;
      }
    }
    break;
#ifndef WIN32 // Win32 can't have 8byte integers
  case 8:
    {
      Long *src=(Long*)source;
      Long *dst=(Long*)dest;
      for(long i=0;i<nelements;i++){
	union LongChar {
	  Long l;
	  char a[8];
	};
	register LongChar ts,td;
	// do the swap thang on each element
	td.l=ts.l=src[i];
	td.a[0]=ts.a[7];
	td.a[1]=ts.a[6];
	td.a[2]=ts.a[5];
	td.a[3]=ts.a[4];
	td.a[4]=ts.a[3];
	td.a[5]=ts.a[2];
	td.a[6]=ts.a[1];
	td.a[7]=ts.a[0];
	dst[i]=td.l;
      }
    }
    break;
#endif
  default:
    {
      char *src=(char *)source; // treat as a character buffer
      char *dst=(char *)dest;
      for(long i=0;i<nelements;i++,src+=elementsize,dst+=elementsize){
	register int s,d;
	// do the swap thang on each element
	for(s=0,d=elementsize-1;s<d;s++,d--){
	  dst[d]=src[s];
	  dst[s]=src[d];
	}
      }
    }
    break;
  }
}

int IEEEIO::writeFileHeader(){
  // long pos=getPosition(); // use virtual_position
  // this->flush();
  //puts("WriteHeader");
  for(char c=0;c<4;c++) // make sure magic is correct
    file_header.magic.c[c]=c+1;
  if(swapbytes)
    file_header.byteorder.i=IEEE_REVERSE_MAGIC;
  else 
    file_header.byteorder.i=IEEE_NATIVE_MAGIC;
  file_header.majorversion = IEEE_MAJORVERSION;
  file_header.minorversion = IEEE_MINORVERSION;
  IEEEIO::lseek(0,L_SET); // will force a flush (if needed)
  //puts("done writeheader seek");
  if(swapbytes) file_header.byteswap();
  // T3E Kludge
  int sz=FileHdr::write(this,file_header);
  if(swapbytes) file_header.byteswap();
  // if(pos>=FileHdrSize)  // T3E Kludge
  //  IEEEIO::lseek(pos,L_SET);
  return sz;
}

// reads the start of the file (including the magic number)
// and determines if this is 
//       1) A valid IEEEIO file
//       2) what the byte order is
//       3) How many datasets are contained
int IEEEIO::readFileHeader(){
  long pos=getPosition(); // use virtual_position
  IEEEIO::lseek(0,L_SET);
  // T3E Kludge
  int sz=FileHdr::read(this,file_header);
  IEEEIO::lseek(pos,L_SET); // return to original position
  for(char c=0;c<4;c++){
    if((c+1)!=file_header.magic.c[c]){
      file_header.ndatasets=0;
      fprintf(stderr,"IEEEIO: File %s has the wrong magic number\n",filename);
      return -1; // bad magic
    }
  }
  // now compare byte order
  if(file_header.byteorder.i==IEEE_NATIVE_MAGIC)
    swapbytes=0;
  else{
    swapbytes=1;
    //puts("byteswapping is On !!!!!!!!!!!!!!!");
  }
  if(swapbytes) file_header.byteswap();
  ndatasets=file_header.ndatasets;
  return 1;
}

void IEEEIO::rebuildFileHeader(){
  int lndatasets;
  long pos=getPosition(); // save file position (use virtual position
  restart();
  // steal seek loop from the nextDataRecord()
  for(lndatasets=0;
      nextDataRecord()>0;
      lndatasets++){
    // now check the length of the record
    // for each record.
  }
  file_header.ndatasets=lndatasets;
  writeFileHeader();
  IEEEIO::lseek(pos,L_SET);
}

void IEEEIO::appendRecordTable(){
  //printf("appending records==================\n");
  while(nextRecord()>0){
    //printf("appending record\n");
    switch(current_rec.recordtype){
    case DataRecord:
      //fprintf(stderr,"DataRecord");
      (rec_table[current_rec.sequenceID]).rec=current_rec;
      // use virtual_position
      (rec_table[current_rec.sequenceID]).offset=getPosition()-RecordHdrSize; 
      break;
    case AnnotationRecord:{
      //fprintf(stderr,"\tAnnotationRecord\n");
      RecRef ref;
      ref.rec=current_rec;
      ref.offset=getPosition()-	RecordHdrSize; // T3E Kludge
      (rec_table[current_rec.sequenceID]).annotations.append(ref);
    }
    break;
    case AttributeRecord:{
      //fprintf(stderr,"\tAttributeRecord\n");
      AttribRef ref,*refp;
      AttributeRecordHdr attribhdr;
      ref.rec=current_rec;
      ref.offset=getPosition()- RecordHdrSize; // T3E Kludge
      (rec_table[current_rec.sequenceID]).attributes.append(ref);
      int idx=(rec_table[current_rec.sequenceID]).attributes.getSize()-1;
      refp = &((rec_table[current_rec.sequenceID]).attributes[idx]);
      AttributeRecordHdr::read(this,attribhdr);
      IEEEIO::read(refp->name,attribhdr.namesize);
    }
    break;
    default:
      fprintf(stderr,"\tIEEEIO::Error UNKNOWN RECORD TYPE [%d]... recovering.",
	      (int)(current_rec.recordtype));
      return;
    }
    (rec_table[current_rec.sequenceID]).end_offset = 
      current_rec_offset + current_rec.recordsize;
  }
}

void IEEEIO::buildRecordTable(){
  if(file_header.ndatasets<0)
    return; // failure
  restart();
  rec_table.setSize(file_header.ndatasets);
  appendRecordTable();
}

void IEEEIO::openFile(CONST char *fname,IObase::AccessMode access,int swbytes){
  switch(access){
  case IObase::Read:
    fid=open(fname,O_RDONLY,0644);
    if(fid<0){
      fprintf(stderr,"IEEEIO: Failed to open %s for reading\n",fname);
      return;
    }
    initPosition();
    if(readFileHeader()<=0){
      close(fid);
      fprintf(stderr,"IEEEIO: File %s is empty (opened for reading)\n",fname);
      fid=-1; // invalid file
      return;
    }
    if(file_header.ndatasets<0){
      // must recover from crash
      close(fid);
      fid=open(fname,O_RDWR,0644);
      rebuildFileHeader();
      close(fid);
      fid=open(fname,O_RDONLY,0644);
      readFileHeader(); // reread..
    }
    buildRecordTable();
    restart(); // go to start of file
    ndatasets = file_header.ndatasets;
    nextRecord(); // prime the recordnumber
    break;
  case IObase::Write: // truncates
    fid=open(fname,O_WRONLY|O_CREAT|O_TRUNC,0644);
    if(fid<0){
      fprintf(stderr,"IEEEIO: Failed to create %s for writing\n",fname);
      return;
    }
    initPosition();
    file_header.ndatasets=-1;
    ndatasets=0;
    writeFileHeader();
    restart();
    break;
  case IObase::Append:
  default:
    fid=open(fname,O_RDWR,0644); // PW: We *don't* want O_APPEND here.
				// O_Append moves the file ptr to end
				// of file before a write; this means
				// when we do the writeHeader() below
				// this is stuck at the end of this file
				// which makes it puke.
    if(fid<0){
      fprintf(stderr,"IEEEIO: Failed to open %s for append\n",fname);
      return;
    }
    initPosition();
    if(readFileHeader()<=0){
      close(fid);
      fprintf(stderr,"IEEEIO: File %s is empty (opened for append)\n",fname);
      fid=-1; // invalid file
      return;
    }
    if(file_header.ndatasets<0){
      // must recover from crash
      close(fid);
      fid=open(fname,O_RDWR,0644);
      initPosition();
      rebuildFileHeader();
      close(fid);
      fid=open(fname,O_RDWR,0644);	// PW: We want to open RDWR again
      initPosition();
      readFileHeader(); // reread..
    }
    buildRecordTable();
    restart(); // go to start of file
    ndatasets = file_header.ndatasets;
    // now we sync up in the way we would for a write
    file_header.ndatasets=-1;
    writeFileHeader(); // write header with -1 to indicate we are writing
    // could refcount by increasing the - of numbers in header 
    // (needs extra flag for synced)
    restart();
    nextRecord(); // prime the recordnumber
    //seek(ndatasets);
    break;
  }
}

IEEEIO::IEEEIO(IEEEIO *file): // read only dup
  IObase(file->filename,IObase::Read),fid(dup(file->fid)),
  swapbytes(file->swapbytes),datasetnumber(0),
  ndatasets(file->ndatasets),hasread(0),writebuffer(0),
  writebuffersize(0),writebuffercount(0),savedposition(-1),
  file_length(-1),actual_position(-1),virtual_position(-1)
{
  if(file->masterfile) masterfile=file->masterfile;
  else masterfile=file;
  initPosition(); // initial positional pointers
  ndatasets = file_header.ndatasets = file->ndatasets;
  buildRecordTable();
  restart(); // go to start of file
  nextRecord(); // prime the recordnumber 
}

// IEEEIO::AccesMode has scope commented out to satisfy Microsoft compiler
IEEEIO::IEEEIO(CONST char *fname,/*IEEEIO::*/AccessMode access,int swbytes):
  IObase(fname,mode(access)),fid(-1),swapbytes(swbytes),datasetnumber(0),
  ndatasets(0),hasread(0),masterfile(0),writebuffer(0),
  writebuffersize(0),writebuffercount(0),savedposition(-1),
  file_length(-1),actual_position(-1),virtual_position(-1)
{
  if(access==IEEEIO::SharedRead){
    masterfile=this; // we are multi-reading
    fid=open(fname,O_RDONLY,0644);
    if(fid<0){
      fprintf(stderr,"IEEEIO: Failed to open %s for reading\n",fname);
      return;
    }
    this->initPosition();
    if(readFileHeader()<=0){
      close(fid);
      fprintf(stderr,"IEEEIO: File %s is empty (opened for reading)\n",fname);
      fid=-1; // invalid file
      return;
    }
    buildRecordTable();
    restart(); // go to start of file
    ndatasets = file_header.ndatasets;
    nextRecord(); // prime the recordnumber
  }
  else
    openFile(fname,mode(access),swbytes);
}

IEEEIO::IEEEIO(CONST char *fname,IObase::AccessMode access,int swbytes):
  IObase(fname,access),fid(-1),swapbytes(swbytes),datasetnumber(0),
  ndatasets(0),hasread(0),masterfile(0),writebuffer(0),
  writebuffersize(0),writebuffercount(0),savedposition(-1)
{
  //  long fpos;
  openFile(fname,access,swbytes);
}

IEEEIO::~IEEEIO(){
  if(fid<0 && savedposition>=0) resume();
  // resume IO only if it is paused so that
  // we can do the final writes to the file before
  // closing
  //puts("do bufferoff");
  if(writebuffer) bufferOff(); // automatically flushes buffer
  //puts("now rewrite header");
  if(fid>=0){
    if(accessmode!=IObase::Read){
      file_header.ndatasets=ndatasets;
      writeFileHeader();
    }
    close(fid);
  }
  fid=-1;
}

int IEEEIO::isValid(){
  if(fid>=0) return 1;
  else return 0;
}

int IEEEIO::write(IObase::DataType typeID,int rank,CONST int *dims,void *data){
  int i;
  RecordHdr rec;
  DataRecordHdr hdr;
  // make sure its the EOF
  if(accessmode==IObase::Read) return 0;
  hasread=0; // reset hasread; (JMS: changed from local Thu Mar 12 13:53:13 CST 1998)
  if(rec_table.getSize()>0){
    long endpos = rec_table[rec_table.getSize()-1].end_offset;
    IEEEIO::lseek(endpos,L_SET);  // don't know if this is costly
  }
  else 
    IEEEIO::lseek(0,L_XTND); // seek to end of file
  // if 0-length seeks are costly, then alternative
  // logic can be constructed to ensure writes to end of file
  for(i=0,hdr.datasize=1;i<rank;i++) 
    hdr.datasize*=dims[i];

  if(chunkdims.getSize()<rank)
    chunkdims.setSize(rank);
  for(i=0;i<rank;i++) chunkdims[i]=dims[i];
  
  hdr.datasize*=sizeOf(typeID);
  hdr.numbertype=typeID;
  hdr.rank=rank;
  // if last annotation slot is filled, it is a pointer
  // to another block of 8 annotation pointers.
  rec.recordtype = DataRecord;
  rec.recordsize = hdr.datasize + 
    DataRecordHdrSize +
    sizeof(Int) * rank;
  rec.sequenceID = datasetnumber = ndatasets++;
  // need to byteswap the header and records if swapbytes==True
  current_dat=hdr; // first copy native info to current
  current_rec=rec;
  RecordHdr::write(this,rec);
  current_dat_offset=current_rec_offset=getPosition();
  DataRecordHdr::write(this,hdr);
  // write the dims.... (byteswap if necessary)
  if(swapbytes) byteswapBuffer(chunkdims.getData(),rank,sizeof(Int));
  IEEEIO::write(chunkdims.getData(),sizeof(Int)*rank);
  if(swapbytes) byteswapBuffer(chunkdims.getData(),rank,sizeof(Int)); // swap back
  if(swapbytes) byteswapBuffer(data,hdr.datasize/sizeOf(typeID),sizeOf(typeID));
  int sz = IEEEIO::write(data,hdr.datasize);
  // yeah! double bytswap seem stupid, but consider the alternatives
  // write one byte at a time (swapping as we go) == systemcall overhead
  // copy to a temporary buffer and swap bytes == malloc, copy, and free!
  // strangely its faster to swap twice
  if(swapbytes) byteswapBuffer(data,hdr.datasize/sizeOf(typeID),sizeOf(typeID));
  DataRef dref;
  dref.rec=current_rec;
  // Must find exact size!!!
  dref.offset = current_rec_offset - RecordHdrSize;
  dref.end_offset = current_rec_offset + current_rec.recordsize;
  if(accessmode!=IObase::Write) // do not store if in write-only mode
    rec_table.append(dref);
  return sz;
}

//setCurrentRec(offset);
// storing in 8*1024 byte blocks (always aligned to chunk boundaries)
// will result in more efficient seeking behavior due to 
// disk block size... thus file is always aligned to unix
// filesystem blocks. (for now, we'll use die method maximo-stupido...
int IEEEIO::readInfo(IObase::DataType &typeID,int &rank,int *dims,int maxdims){
  if(accessmode!=IObase::Read && accessmode!=IObase::Append)
    return 0;
  //int sz;
  if(hasread) datasetnumber++; // increment to the next one if this has been read
  if(datasetnumber>=rec_table.getSize()) return 0; // end of file
  // read the record + header
  hasread=1; // we've read this, so next time we hit it, we'll increment again
  IEEEIO::lseek(rec_table[datasetnumber].offset,L_SET);
  // T3E Kludge
  RecordHdr::read(this,current_rec);
  DataRecordHdr::read(this,current_dat);
  rank=current_dat.rank;
  if(chunkdims.getSize()<rank) chunkdims.setSize(rank);
  typeID=Int2DataType(current_dat.numbertype);
  IEEEIO::read(chunkdims.getData(),sizeof(Int)*rank);
  if(swapbytes) byteswapBuffer(chunkdims.getData(),(rank>maxdims)?maxdims:rank,sizeof(Int));
  
  for(int i=0;i<maxdims && i<rank;i++) dims[i]=chunkdims[i];
  streaming=0;
  return 1;
}

/*
  It does not appear that the bytswapping is really configured here...
 */
int IEEEIO::read(void *data){
  if(accessmode!=IObase::Read && accessmode!=IObase::Append)
    return 0; // can't readif write only
  if(datasetnumber>=rec_table.getSize()) return 0; // seek past end.
  long datapos=rec_table[datasetnumber].offset+RecordHdrSize+DataRecordHdrSize+
    sizeof(Int)*current_dat.rank;
  // should make certain current file position is correct
  // if(getPosition() != datapos) (redundant call to lseek())
  IEEEIO::lseek(datapos,L_SET); // seek to position
  int sz=IEEEIO::read(data,current_dat.datasize); // read nbytes
  int typelen = sizeOf(Int2DataType(current_dat.numbertype));// compute elem sz
  if(swapbytes) byteswapBuffer(data,current_dat.datasize/typelen,typelen);
  return sz;
}

int IEEEIO::seek(int idx){ 
  if(accessmode!=IObase::Read &&
     accessmode!=IObase::Append) // can't seek unless readable
    return -1; // failed to seek
  if (idx >= ndatasets || idx < 0) {
    return -1;	// do a bounds check (probably should just clip)
  }
  // bound the index
  if((1+idx)>rec_table.getSize()) idx=rec_table.getSize();
  if(idx<0) idx=0;
  index = idx;
  IEEEIO::lseek(rec_table[index].offset,L_SET);
  // T3E Kludge
  RecordHdr::read(this,current_rec);
  current_rec_offset=getPosition();
  datasetnumber=current_rec.sequenceID; // Changed:  had -1
  hasread=0; // new file position
  return current_rec.sequenceID;
}

int IEEEIO::nDatasets() {
  // can work across threads or across processes since
  // it gathers information directly from the file.
  // it does not need to share info with masterfile
  // in fact you can delete the masterfile
  if(masterfile){ // we are passive read-only. update nDatasets
    // seek-scan to find datasets?
    int oldlength = (rec_table.getSize())?(rec_table[rec_table.getSize()-1].end_offset):0;
    //printf("oldlength=%u newlength=%u\n",oldlength,(unsigned int)(getLength()));
    if(getLength()>oldlength){
      // lets scan to find the end
      int lndatasets;
      puts("We should clearly not be here!!!");
      restart(); // go to beginning
      // steal seek loop from the nextDataRecord()
      for(lndatasets=0;
	  nextDataRecord()>0;
	  lndatasets++){
	//printf("scan dataset[%u]\n",lndatasets);
      }
      file_header.ndatasets = lndatasets;
      rec_table.setSize(lndatasets);
      // printf("counted %u datasets. old was %u datasets\n",lndatasets,ndatasets);
      seek(ndatasets); // seek to current last
      ndatasets = lndatasets; // then reset the last
      appendRecordTable(); // append new stuff to record table
    }
  }
  return ndatasets; 
}

int IEEEIO::writeAnnotation(CONST char *annotation){
  int stringlen;
  RecRef aref;
  RecordHdr rec;
  if(accessmode==IObase::Read || !annotation) return 0;
  stringlen=strlen(annotation)+1;

  if(rec_table.getSize()>0){
    long endpos = rec_table[rec_table.getSize()-1].end_offset;
    IEEEIO::lseek(endpos,L_SET);  // don't know if this is costly
  }
  else 
    IEEEIO::lseek(0,L_XTND); // seek to end of file

  rec.recordtype=AnnotationRecord;
  rec.recordsize=stringlen;
  if(datasetnumber>=0) rec.sequenceID=datasetnumber;
  else rec.sequenceID=current_rec.sequenceID; // a kludge for error immunity
  current_rec = rec;
  // T3E Kludge
  current_rec_offset = getPosition() + RecordHdrSize;
  // T3E Kludge
  RecordHdr::write(this,rec);
  // if(swapbytes) byteswap(rec);
  int sz=IEEEIO::write(annotation,stringlen);

  aref.rec=current_rec;
  // T3E Kludge
  aref.offset=current_rec_offset-RecordHdrSize;
  if(accessmode!=IObase::Write){ // don't write to index cache if in write-only mode
    rec_table[datasetnumber].annotations.append(aref); 
    rec_table[datasetnumber].end_offset=getPosition();
  }
  return sz;
}
/*
	If objects could have all attributes implicitly availible for querying
	when passed to subroutines.  Basic object methods for data movement like
	linearize() or object.linearized[index] object.linearized.size.
	Then for objects that a receiver can't deal with either
		1: have compiletime "gatekeepers" that restrict allowable types
			that it can receive
		2: Do a treesearch through object space to find a convertor sequence
			that can convert the current type into the one that the object
			accepts.  If no such object can be found, indicate runtime
			object adaptor failure.
*/
int IEEEIO::readAnnotationInfo(int number,int &length){
  if(datasetnumber<0 ||
     number>=rec_table[datasetnumber].annotations.getSize()){ 
    length=0;
    return -1;
  }
  length=rec_table[datasetnumber].annotations[number].rec.recordsize;
  return length;
}

int IEEEIO::readAnnotation(int number,char *annotation,int maxsize){
  // returns actual size of annotation or -1 if error
  if(datasetnumber<0 || accessmode==IObase::Write ||
     number>=rec_table[datasetnumber].annotations.getSize()){ 
    return -1;
  }
  IEEEIO::lseek(rec_table[datasetnumber].annotations[number].offset,L_SET);
  // T3E Kludge
  RecordHdr::read(this,current_rec);
  IEEEIO::read(annotation,(maxsize<rec_table[datasetnumber].annotations[number].rec.recordsize)?
       maxsize:(rec_table[datasetnumber].annotations[number].rec.recordsize));
  annotation[maxsize-1]='\0';
  return rec_table[datasetnumber].annotations[number].rec.recordsize;
}

int IEEEIO::nAnnotations(){
  if(datasetnumber<0 || accessmode==IObase::Write) 
    return -1;
  return rec_table[datasetnumber].annotations.getSize();
}

// for attributes
int IEEEIO::writeAttribute(CONST char *name,IObase::DataType typeID,Long length,void *data){
  int stringlen;
  AttribRef aref;
  RecordHdr rec;
  AttributeRecordHdr attrib;

  if(datasetnumber<0){
    fprintf(stderr,"IEEEIO::writeAttribute():  Error, cannot write attribute before any datasets have been written!\n");
    return 0;
  }
  if(accessmode==IObase::Read) return 0;
  stringlen=strlen(name)+1;

  if(rec_table.getSize()>0){
    long endpos = rec_table[rec_table.getSize()-1].end_offset;
    IEEEIO::lseek(endpos,L_SET);  // don't know if this is costly
  }
  else 
    IEEEIO::lseek(0,L_XTND); // seek to end of file
  attrib.datasize=length*sizeOf(typeID);
  attrib.namesize=stringlen;
  attrib.numbertype=typeID;
  rec.recordtype=AttributeRecord;
  rec.recordsize=attrib.datasize+attrib.namesize+AttributeRecordHdrSize;
  if(datasetnumber>=0) rec.sequenceID=datasetnumber;
  else rec.sequenceID=current_rec.sequenceID; // a kludge for error immunity
  current_rec=rec;
  current_rec_offset=getPosition() + RecordHdrSize;
  aref.rec=current_rec;
  aref.offset=current_rec_offset-RecordHdrSize;
  // no copy constructor for aref, so must do manually
  if(accessmode != IObase::Write){
    int lastindex;
    rec_table[datasetnumber].attributes.append(aref);
    lastindex=rec_table[datasetnumber].attributes.getSize()-1;
    strcpy(rec_table[datasetnumber].attributes[lastindex].name,name);
  }
  // T3E Kludge
  RecordHdr::write(this,rec);
  AttributeRecordHdr::write(this,attrib);
  //printf("\tnow write the string\n");
  IEEEIO::write(name,stringlen);
  //printf("\tdone IEEEIO::stringlen=%u typeid=%u\n",attrib.namesize,attrib.numbertype);
  // data is not byte-reversed...
  int sz=0;
  if(swapbytes) byteswapBuffer(data,length,sizeOf(typeID));
  IEEEIO::write(data,length*sizeOf(typeID));
  if(swapbytes) byteswapBuffer(data,length,sizeOf(typeID)); // swap back
  if(accessmode != IObase::Write)
    rec_table[datasetnumber].end_offset=getPosition();
  // doesn't appear to store attributes properly
  return sz;
}

int IEEEIO::readAttributeInfo(int number,char *name,IObase::DataType &typeID,
			      Long &nelem,int maxnamelen){
  if(accessmode==IObase::Write) return -1;
  FlexArray<AttribRef> *attribs= &(rec_table[datasetnumber].attributes);
  if(number>=(*attribs).getSize()) return -1; // > number of attributes
  AttribRef *attrib=&((*attribs)[number]);
  if(strlen(attrib->name)>maxnamelen){
    strncpy(name,attrib->name,maxnamelen); // don't we want to copy the other way?
    name[maxnamelen-1]='\0';
  }
  else strcpy(name,attrib->name);
  name[maxnamelen-1]='\0'; // make certain it is null capped
  AttributeRecordHdr attribhdr;
  IEEEIO::lseek((*attribs)[number].offset,L_SET);
  RecordHdr::read(this,current_rec);
  AttributeRecordHdr::read(this,attribhdr);
//printf("IEEEIO:attribrechdr = ds,nt,ns %d,%d,%d\n",
//	 attribhdr.datasize,
//	 attribhdr.numbertype,
//	 attribhdr.namesize);
  typeID = Int2DataType(attribhdr.numbertype);
  nelem = attribhdr.datasize/sizeOf(typeID);
  return -1;
}

int IEEEIO::readAttributeInfo(CONST char *name,IObase::DataType &typeID,Long &nelem){
  // by name
  if(accessmode==IObase::Write) return -1;
  FlexArray<AttribRef> *attribs= &(rec_table[datasetnumber].attributes);
  for(int i=0;i<attribs->getSize();i++){
    if(!strcmp(name,(*attribs)[i].name)){
      // must read the record to get this info
      AttributeRecordHdr attribhdr;
      IEEEIO::lseek((*attribs)[i].offset,L_SET);
      RecordHdr::read(this,current_rec);
      AttributeRecordHdr::read(this,attribhdr);
      typeID = Int2DataType(attribhdr.numbertype);
      nelem = attribhdr.datasize/sizeOf(typeID);
      return i;
    }
  }
  return -1; // Attribute not found
}

int IEEEIO::readAttribute(int number,void *data){
  if(accessmode==IObase::Write) return -1;
  FlexArray<AttribRef> *attribs= &(rec_table[datasetnumber].attributes);
  if(number>=(*attribs).getSize()) return -1; // > number of attributes
  AttributeRecordHdr attribhdr;
  IEEEIO::lseek((*attribs)[number].offset,L_SET);
  RecordHdr::read(this,current_rec);
  AttributeRecordHdr::read(this,attribhdr);
  IObase::DataType typeID = Int2DataType(attribhdr.numbertype);
  if(typeID==IObase::Error) return -1; // read failed due to bad datatype
  long nelem = attribhdr.datasize/sizeOf(typeID);
  IEEEIO::lseek(attribhdr.namesize,L_INCR);
  int sz=IEEEIO::read(data,attribhdr.datasize)/sizeOf(typeID);
  if(swapbytes) byteswapBuffer(data,sz,sizeOf(typeID));
  if(typeID==IObase::String){
    char *cdata=(char *)data;
    cdata[sz]='\0'; /* Null Terminate String data */
  }
  return sz; // returns number of elements read
}

int IEEEIO::nAttributes(){
  if(datasetnumber<0 || accessmode == IObase::Write) return -1;
  return rec_table[datasetnumber].attributes.getSize();
}

//================Chunking Interface-----------------------
void  IEEEIO::clearChunk(int nbytes){
#ifdef T3E
#define DISKBLOCKSIZE 128*1024
#else
#define DISKBLOCKSIZE 8192
#endif
// This is a wasteful feature that will be eliminated in later releases
// This is currently done because I needed the disk space to be zero'ed 
// prior to writing for debugging purposes...  Could have a writeStream 
// class.
  int nwritten;
  char dummy[DISKBLOCKSIZE];
  for(int i=0;i<DISKBLOCKSIZE;i++) 
    dummy[i]=0; // bzero it (reserve space with zero'ed data)
  while(nbytes>DISKBLOCKSIZE){
    nwritten=IEEEIO::write(dummy,DISKBLOCKSIZE);
    if(nwritten<=0) return; // this aint working
    nbytes-=nwritten;
  }
  if(nbytes>0) IEEEIO::write(dummy,nbytes); // write the remaining bytes
#undef DISKBLOCKSIZE
}

int IEEEIO::reserveChunk(IObase::DataType typeID,int rank,CONST int *dims){
  reserveStream(typeID,rank,dims); // sets up the chunking state machine
  // clearing does not appear to be needed... just need header
  //clearChunk(IObase::nBytes(typeID,rank,dims)); // this pre-clears the reserved area to 0. (inefficient)
  return 1;
}

// Chunking should probably include striding....
// But thats a pain since you'll need to read + write to
// Accomplish that (yuck!).  Lots of seeking.
// Lets leave striding out for now.

// should also check for contiguous data (eg. slicing) and optimize for that layout
// for now, the streaming interface best serves contiguous data.
int IEEEIO::writeChunk(CONST int *dims,CONST int *origin,void *data){
  int i,sz;
  // make sure its the EOF
  if(accessmode==IObase::Read) {
    fprintf(stderr,"IEEEIO::writeChunk(): Error!  Access is ReadOnly\n");
    return 0;
  }
  // should have a "chunk reserved" flag set.
  if(current_reserved_chunk!=datasetnumber){
    fprintf(stderr,"IEEEIO::writeChunk(): Error! You forgot to reserve space for the chunk using IO::reserveChunk()\n");
    return 0;
  }
  
  // now we need to seek to the position of the data
  int rank = current_dat.rank;
  IObase::DataType typeID = IObase::Int2DataType(current_dat.numbertype);
  int typesize=sizeOf(typeID);
  // long basefileoffset = rec_table[datasetnumber].offset+
  long basefileoffset = stream_offset +
    RecordHdrSize+DataRecordHdrSize+sizeof(Int)*current_dat.rank; // T3E Kludge
  long chunkcolumnsize=dims[0]*typesize; // stride between columns in chunk
  long filecolumnsize=chunkdims[0]*typesize; // stride between columns on disk (full array size)
  long originoffset;
  long accumdims;
  for(originoffset=0,accumdims=typesize,i=0;i<rank;i++){
    originoffset += (origin[i]*accumdims);
    accumdims *= chunkdims[i];
  }
  for(i=0;i<rank;i++){
    if((origin[i]+dims[i])>chunkdims[i]){
      fprintf(stderr,"IEEEIO::writeChunk(): ERROR!!  specified dims and origin exceed reserved block size\n");
      fprintf(stderr,"\tfor dimension %u origin=%u and dims=%u\n",i,origin[i],
	      dims[i]);
      fprintf(stderr,"\torigin+dims=%u whereas the maximum must be less than %u\n",origin[i]+dims[i],chunkdims[i]);
    }
  }
  originoffset+=basefileoffset; // for absolute seeking
  long maxindex,minindex;
  minindex=basefileoffset;
  maxindex=basefileoffset+current_dat.datasize;
  int ncolumns; // computed in loop below
  for(ncolumns=1,i=1;i<rank;i++) ncolumns*=dims[i];
  if(swapbytes) 
    byteswapBuffer(data,IObase::nElements(rank,dims),typesize);
  for(i=0;i<ncolumns;i++){ // read the columns
    if(originoffset<minindex){
      fprintf(stderr,"WriteChunk() inconsistency. Writing less than min index\n");
      fprintf(stderr,"\tCol[%u]: Requested %u, but minidex= %u\n",
	      i,(unsigned int)originoffset,(unsigned int)minindex);
    }
    if(originoffset>maxindex){
      fprintf(stderr,"WriteChunk() inconsistency. Writing greater than maximum index\n");
      fprintf(stderr,"\tCol[%u]: Requested %u, but maxindex= %u\n",
	      i,(unsigned int)originoffset,(unsigned int)maxindex);
    }
    if((originoffset+chunkcolumnsize)>maxindex){
      fprintf(stderr,"WriteChunk() inconsistency. This write will overrun the reserved data\n");
      fprintf(stderr,"\tCol[%u]: Requested %u, maxindex %u, and the write of %u will run to %u\n",
	      i,(unsigned int)originoffset,(unsigned int)maxindex,(unsigned int)chunkcolumnsize,(unsigned int)(originoffset+chunkcolumnsize));
    }
    IEEEIO::lseek(originoffset,L_SET);
    sz=IEEEIO::write(((char*)data)+i*chunkcolumnsize,(int)chunkcolumnsize);
    originoffset+=filecolumnsize;
    for(long j=1,idx=dims[1],planesize=filecolumnsize;
	j<(rank-1) && !((i+1)%idx);
	idx*=dims[++j]){
      long extraoffset=planesize*(chunkdims.getData()[j]-dims[j]);
      originoffset+=extraoffset;
      planesize*=dims[j];
    }
  }
  // now swap the data back to native order...
  if(swapbytes) 
    byteswapBuffer(data,IObase::nElements(rank,dims),typesize);

  return sz;
}
int IEEEIO::readChunk(CONST int *dims,CONST int *origin,void *data){ 
  int sz,i;
  if(accessmode!=IObase::Read && accessmode!=IObase::Append)
    return 0;
  // gonna have to stride through this sucker... (yuck)
  // use the same chunkdims for this.  Must set during readinfo...
  // now we need to seek to the position of the data
  int rank = current_dat.rank;
  IObase::DataType typeID = IObase::Int2DataType(current_dat.numbertype);
  int typesize=sizeOf(typeID);
  long basefileoffset = rec_table[datasetnumber].offset+
    // sizeof(RecordHdr)+sizeof(DataRecordHdr)+sizeof(Int)*current_dat.rank;
    RecordHdrSize+DataRecordHdrSize+sizeof(Int)*current_dat.rank; // T3E Kludge
  long chunkcolumnsize=dims[0]*typesize; // stride between columns in chunk
  long filecolumnsize=chunkdims[0]*typesize; // stride between columns on disk (full array size)
  long originoffset;
  long accumdims;
  // compute the offset into the data on disk required by the chunk origin (initial offset)
  for(originoffset=0,accumdims=typesize,i=0;i<rank;i++){
    originoffset += (origin[i]*accumdims);
    accumdims *= chunkdims[i];
  }
  originoffset+=basefileoffset; // for absolute seeking
  long ncolumns; // computed in loop below
  for(ncolumns=1,i=1;i<rank;i++) ncolumns*=dims[i];
  for(i=0;i<ncolumns;i++){ // read the columns
    IEEEIO::lseek(originoffset,L_SET);
    sz=IEEEIO::read(((char*)data)+i*chunkcolumnsize,(int)chunkcolumnsize); 
    originoffset+=filecolumnsize;
    for(long j=1,idx=dims[1],planesize=filecolumnsize;
	j<(rank-1) && !((i+1)%idx);
	idx*=dims[++j]){
      long extraoffset=planesize*(chunkdims.getData()[j]-dims[j]);
      originoffset+=extraoffset;
      planesize*=dims[j];
    }
  }
  // now swap the data back to native order...
  if(swapbytes)
    byteswapBuffer(data,IObase::nElements(rank,dims),typesize);
  return sz;
}
// Nearly identical to reserveChunk except that we don't need to pre-clear
// the space to 0 for streamed data.
int IEEEIO::reserveStream(IObase::DataType typeID,int rank,CONST int *dims){
  int i;
  RecordHdr rec;
  DataRecordHdr hdr;
  // make sure its the EOF
  if(accessmode==IObase::Read) return 0;
  
  if(rec_table.getSize()>0){
    long endpos = rec_table[rec_table.getSize()-1].end_offset;
    IEEEIO::lseek(endpos,L_SET);  // don't know if this is costly
  }
  else {
    IEEEIO::lseek(0,L_XTND); // seek to end of file
  }
  // if 0-length seeks are costly, then alternative
  // logic can be constructed to ensure writes to end of file
  if(chunkdims.getSize() != rank)
    chunkdims.setSize(rank); // KLUDGE!
  for(i=0,hdr.datasize=1;i<rank;i++) {
    hdr.datasize*=dims[i];
    chunkdims[i]=dims[i];
  }
  cur_type_size=sizeOf(typeID); // KLUDGE: for writeStream()
  hdr.datasize*=sizeOf(typeID);
  hdr.numbertype=typeID;
  hdr.rank=rank;
  // if last annotation slot is filled, it is a pointer
  // to another block of 8 annotation pointers.
  rec.recordtype = DataRecord;
  rec.recordsize = hdr.datasize + 
    DataRecordHdrSize +  // T3E Kludge
    sizeof(Int) * rank;
  rec.sequenceID = datasetnumber = ndatasets++;
  // need to byteswap the header and records if swapbytes==True
  current_dat=hdr; // first copy native info to current
  current_rec=rec;
  RecordHdr::write(this,rec);
  current_dat_offset=current_rec_offset=getPosition();
  DataRecordHdr::write(this,hdr);
  if(swapbytes) byteswapBuffer(chunkdims.getData(),rank,sizeof(Int));
  IEEEIO::write(chunkdims.getData(),sizeof(Int)*rank); // for better for worse...
  if(swapbytes) byteswapBuffer(chunkdims.getData(),rank,sizeof(Int)); // swap back
  // OK, now writing 8k at a time, reserve the space with 0's
  current_reserved_chunk=datasetnumber;
  DataRef dref;
  dref.rec=current_rec;
  dref.offset=current_rec_offset-RecordHdrSize;
  dref.end_offset=current_rec_offset+current_rec.recordsize;
  if(accessmode!=IObase::Write)
    rec_table.append(dref);
  // else (JMS debugging change 5/12/98)
  stream_offset = current_rec_offset-RecordHdrSize;
  streaming=0;
  return 1;
}

int IEEEIO::writeStream(void *data,int length){
  int len;
  if(!streaming){// seek to starting offset
    long basefileoffset = stream_offset +
      RecordHdrSize+DataRecordHdrSize+sizeof(Int)*current_dat.rank;
    IEEEIO::lseek(basefileoffset,L_SET);
    streaming=1;
  }
  // need to do a bounds check on the stream... for now it'll play dumb...
  long nbytes=cur_type_size * length;
  len = IEEEIO::write((char *)data,nbytes); // record file position??
  if(len<0)
    return -1; // write failure
  return len;
}

int IEEEIO::readStream(void *data,int length){
  int len;
  if(!streaming){// seek to starting offset
    long basefileoffset = rec_table[datasetnumber].offset+
      RecordHdrSize+DataRecordHdrSize+sizeof(Int)*current_dat.rank;
    IEEEIO::lseek(basefileoffset,L_SET);
    streaming=1;
  }
  // need to do a bounds check on the stream... for now it'll play dumb...
  int typesize= sizeOf(Int2DataType(current_dat.numbertype));
  long nbytes = typesize * length;
  len = IEEEIO::read((char *)data,nbytes);
  if(len<0)
    return -1; // read failure
  else
    if(swapbytes) byteswapBuffer(data,length,typesize);
  return len;
}

void IEEEIO::bufferOn(long size){
  if(writebuffer) {
    if(size==0)
      bufferOff();
    return;
  }
  else {
#if !defined(WIN32) && !defined(T3E)// if not Windows, use fstat() to choose optimal blocksiz
    if(size<0){
      struct stat mystat;
      fstat(fid,&mystat);
      size=mystat.st_blksize;
    }
#endif 
    // if windows then take a wild guess
#ifdef T3E
    size=(size<=64)?8*1024*1024:size; // 8Meg is a good guess for T3E
#else
    size=(size<=64)?8*1024:size; //(8k is a good guess for workstations/WinTel)
#endif
    writebuffer = new char[size];
  }
  writebuffersize=size;
  writebuffercount=0;
}

void IEEEIO::bufferOff(){
  if(!writebuffer) return;
  if(fid<0 && writebuffercount>0){
    fprintf(stderr,"IEEEIO::bufferOff() ERROR!!!!");
    fprintf(stderr,"\tFile %s is either paused or invalid, so I can't deallocate the write buffer safely\n",filename);
    return;
  }
  //printf("bufferoff on a=%ld v=%ld f=%ld b=%d\n",
  //	 actual_position,virtual_position,file_length,writebuffercount);
  this->flush();
  delete writebuffer;
  writebuffer=0;
}

int IEEEIO::pause(){
  if(fid<0) return 0; // fail
  this->flush(); // flush the buffer to pause
  savedposition=getPosition();
  close(fid);
  fid=-1;
  return 1; // success
}

int IEEEIO::resume(){
  if(fid>=0 || savedposition<0) return 0; // fail
  switch(IObase::accessmode){
  case IObase::Read:
    fid=open(filename,O_RDONLY,0644);
    break;
  case IObase::Write:
    fid=open(filename,O_WRONLY,0644);
    break;
  case IObase::Append:
    fid=open(filename,O_RDWR,0644);
    break;
  }
  if(fid<0){
    printf("IEEEIO::resume failed for file %s\n",filename);
    return -1;
  }
  ::lseek(fid,savedposition,L_SET);  // this is OK..
  savedposition=-1;
  return 1; // success
}

//*********************Fortran Interface**********************
Long8 f_ieee_open (char *file,char *accessname,int flen,int alen){
  // would have used tolower(), but it doesn't exist everywhere.... :(
  IObase::AccessMode mode;
  if(*accessname=='R' || *accessname=='r')
    mode=IObase::Read;
  else if(*accessname=='W' || *accessname=='w' ||
	  *accessname=='C' || *accessname=='c')
    mode=IObase::Write;
  else if(*accessname=='A' || *accessname=='a')
    mode=IObase::Append;
  else {
    fprintf(stderr,"IEEEopen(): Error unknown option [%s] to open file %s\n",
	    accessname,file);
    return 0;
  }
  IObase *fid=new IEEEIO(file,mode);
  if(fid->isValid()) 
    return (Long8)fid;
  else
    delete fid; // file open failed
  return 0;
}

Long8 f_ieee_openr(char *filename,int namelen){
  filename[namelen]='\0';
  return (Long8)(new IEEEIO(filename,IObase::Read));
}

Long8 f_ieee_openw(char *filename,int namelen){
  filename[namelen]='\0';
  return (Long8)(new IEEEIO(filename,IObase::Write));
}

Long8 f_ieee_opena(char *filename,int namelen){
  filename[namelen]='\0';
  return (Long8)(new IEEEIO(filename,IObase::Append));
}

void ieee_bufon(Long8 *fileID,int bufsize){
  IEEEIO *f=(IEEEIO*)(*fileID);
  if(bufsize<0) f->bufferOn(); // use default size
  else f->bufferOn(bufsize); // or specify
}

void ieee_bufoff(Long8 *fileID){
  IEEEIO *f=(IEEEIO*)(*fileID);
  f->bufferOff();
}

IOFile IEEEopen(char *filename,char *accessname){
  // Parse all of the ansi stdio access option strings
  IObase::AccessMode mode;
  if(!strcmp(accessname,"read") ||
     !strcmp(accessname,"r") ||
     !strcmp(accessname,"rb"))
    mode=IObase::Read;
  else if(*accessname=='a')
    mode=IObase::Append;
  else if(!strcmp(accessname,"write") || 
	  !strcmp(accessname,"create") ||
	  !strcmp(accessname,"w") ||
	  !strcmp(accessname,"wb") ||
	  !strcmp(accessname,"w+") ||
	  !strcmp(accessname,"w+b") ||
	  !strcmp(accessname,"wb+"))
    mode=IObase::Write;
  else{
    fprintf(stderr,"IEEEopen(): Error unknown option [%s] to open file %s\n",
	    accessname,filename);
    return 0;
  }
  IObase *fid=new IEEEIO(filename,mode);
  if(fid->isValid()) 
    return (IOFile)fid;
  else
    delete fid; // file open failed
  return 0; // unknown option
}

IOFile IEEEopenRead(char *filename){
  return (IOFile)(new IEEEIO(filename,IObase::Read));
}

IOFile IEEEopenWrite(char *filename){
  return (IOFile)(new IEEEIO(filename,IObase::Create));
}

IOFile IEEEopenAppend(char *filename){
  return (IOFile)(new IEEEIO(filename,IObase::Append));
}

void IEEEbufferOn(IOFile fileID,int bufsize){
  IEEEIO *f=(IEEEIO*)fileID;
  if(bufsize<0) f->bufferOn(); // use default size
  else f->bufferOn(bufsize); // or specify
}

void IEEEbufferOff(IOFile fileID){
  IEEEIO *f=(IEEEIO*)fileID;
  f->bufferOff();
}

