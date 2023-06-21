
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "IO.hh"

IObase::IObase(CONST char *fname,AccessMode access):index(0),nitems(0),accessmode(access){
  strcpy(filename,fname);
}

IObase::DataType IObase::Int2DataType(int dt){
  if(dt==Byte)
    return Byte;
  else if(dt==Int16)
    return Int16;
  else if(dt==Int32)
    return Int32;
  else if(dt==Int64)
    return Int64;
  else if(dt==Float32)
    return Float32;
  else if(dt==Float64)
    return Float64;
  else if(dt==uInt8)
    return uInt8;
  else if(dt==uInt16)
    return uInt16;
  else if(dt==uInt32)
    return uInt32;
  else if(dt==uInt64)
    return uInt64;
  else if(dt==Char8)
    return Char8;
  else if(dt==Char16)
    return Char16;
  else return Error;
}

int IObase::sizeOf(IObase::DataType dt){
  switch(dt){
  case Int8:
  case uInt8:
  case Char8:
    return 1;
  case Int16:
  case uInt16:
  case Char16:
    return 2;
  case Int32:
  case uInt32:
  case Float32:
    return 4;
  case Int64:
  case uInt64:
  case Float64:
    return 8;
  default:
    return 0;
  }
}

int IObase::nElements(int rank,CONST int *dims){
  int nelem,i;
  for(i=0,nelem=1;i<rank;i++) nelem*=dims[i];
  return nelem;
}

int IObase::nBytes(IObase::DataType dt, int rank,CONST int *dims){
  return nElements(rank,dims)*sizeOf(dt);
}

int f_io_close(Long8 *deviceID){
  IObase *dev=(IObase*)(*deviceID);
  delete dev;
  *deviceID=0;
  return 1;
}

int f_io_isvalid(Long8 *deviceID){
  if(!deviceID || !*deviceID) return 0;
  IObase *dev=(IObase*)(*deviceID);
  return dev->isValid();
}

int f_io_sizeof(int *type){
  return IObase::sizeOf(IObase::Int2DataType(*type));
}

int f_io_nelements(int *rank,int *dims){
  return IObase::nElements(*rank,dims);
}

int f_io_nbytes(int *type,int *rank,int *dims){
  return IObase::nBytes(IObase::Int2DataType(*type),*rank,dims);
}

int f_io_write(Long8 *deviceID,int *typeID,int *rank,int *dims,void *data){
  IObase *dev=(IObase*)(*deviceID);
  dev->write(dev->Int2DataType(*typeID),*rank,dims,data);
  return 1;
}

int f_io_readinfo(Long8 *deviceID,int *typeID,int *rank,int *dims,int *maxdims){
  IObase *dev=(IObase*)(*deviceID);
  IObase::DataType tid;
  if(*maxdims<0) *maxdims=0;
  int r=dev->readInfo(tid,*rank,dims,*maxdims);
  //int r=dev->readInfo(tid,*rank,dims);
  *typeID=tid;
  return r;
}

int f_io_numdata(Long8 *deviceID){ 
  IObase *dev=(IObase*)(*deviceID);
  return dev->nDatasets();
}
int f_io_read(Long8 *deviceID,void *data){
  IObase *dev=(IObase*)(*deviceID);
  return dev->read(data);
}

int f_io_seek(Long8 *deviceID,int *dataset_index){
  IObase *dev=(IObase*)(*deviceID);
  return dev->seek(*dataset_index);
}
//--------------Annotations
int f_io_writenote(Long8 *deviceID,char *annotation,int size){
  IObase *dev=(IObase*)(*deviceID);
  annotation[size]='\0'; // Yeah! its unsafe I know..!! (but it works!)
  return dev->writeAnnotation(annotation);
}

int f_io_readnote(Long8 *deviceID,int *index,char *annotation,int maxsize){ 
  IObase *dev=(IObase*)(*deviceID);
  return dev->readAnnotation(*index,annotation,maxsize);
}

int f_io_noteinfo(Long8 *deviceID,int *index,int *length){ 
  IObase *dev=(IObase*)(*deviceID);
  return dev->readAnnotationInfo(*index,*length);
}

int f_io_numnote(Long8 *deviceID){ 
  IObase *dev=(IObase*)(*deviceID);
  return dev->nAnnotations();
}

//---------------Attributes
int f_io_writeatt(Long8 *deviceID,char *name,
		 int *datatype,Long *nelements,void *data,int namesize){
  IObase *dev=(IObase*)(*deviceID);
  name[namesize]='\0'; // cap name (to be certain).. unsafe but it works
  // should copy into a flexarray which will be destroyed on return
  IObase::DataType typeID = IObase::Int2DataType(*datatype);
  return dev->writeAttribute(name,typeID,*nelements,data);
}

int f_io_attinfo(Long8 *deviceID,char *name,
		int *datatype,Long *nelements,int namesize){
  IObase *dev=(IObase*)(*deviceID);
  IObase::DataType typeID;
  if(namesize>0)
    name[namesize]='\0';
  int i=dev->readAttributeInfo(name,typeID,*nelements);
  *datatype=typeID;
  return i;
}

int f_io_iattinfo(Long8 *deviceID,int *index,char *name,
		 int *datatype,Long *nelements,int namesize){
  int i;
  IObase *dev=(IObase*)(*deviceID);
  IObase::DataType typeID;
  for(i=0;i<namesize;i++) name[i]='\0';
  //printf("io_iattinfo(): Namesize=%u Name=[%s]\n",namesize,name);
  i=dev->readAttributeInfo(*index,name,typeID,*nelements,namesize);
  // need to zero the array
  //printf("io_iattinfo(): Newname=[%s]\n",name); 
  *datatype=typeID;
  //printf("io_iattinfo(): attribs are index=%u type=%u nelements=%u namesize=%u\n",
  //	 *index,*datatype,*nelements,namesize);
  return i;
}

int f_io_readatt(Long8 *deviceID,int *number,void *data){
  IObase *dev=(IObase*)(*deviceID);
  return dev->readAttribute(*number,data);
}

int f_io_numatt(Long8 *deviceID){
  IObase *dev=(IObase*)(*deviceID);
  return dev->nAttributes();
}

//==========F77 Chunking interface--------------------
int f_io_reserveck(Long8 *deviceID,int *typeID,int *rank,int *dims){
  IObase *dev=(IObase*)(*deviceID);
  return dev->reserveChunk(dev->Int2DataType(*typeID),*rank,dims);
}

int f_io_writeck(Long8 *deviceID,int *chunkdims,int *chunkorigin,void *data){
  IObase *dev=(IObase*)(*deviceID);
  return dev->writeChunk(chunkdims,chunkorigin,data);
}

int f_io_readck(Long8 *deviceID,int *chunkdims,int *chunkorigin,void *data){
  IObase *dev=(IObase*)(*deviceID);
  return dev->readChunk(chunkdims,chunkorigin,data);
}
int f_io_writestrm(IOFile deviceID,void *data,int *length){
  IObase *dev=(IObase*)(deviceID);
  return dev->writeStream(data,*length);
}

int f_io_readstrm(IOFile deviceID,void *data,int *length){
  IObase *dev=(IObase*)(deviceID);
  return dev->readStream(data,*length);
}

int f_io_pause(Long8 *deviceID){
  IObase *dev=(IObase*)(deviceID);
  return dev->pause();
}

int f_io_resume(Long8 *deviceID){
  IObase *dev=(IObase*)(deviceID);
  return dev->resume();
}

//====================C Interface========================
int IOclose(IOFile deviceID){
  IObase *dev=(IObase*)(deviceID);
  delete dev;
  deviceID=0;
  return 1;
}

int IOisValid(IOFile deviceID){
  if(!deviceID) return 0;
  IObase *dev=(IObase*)(deviceID);
  return dev->isValid();
}

int IOsizeOf(int type){
  return IObase::sizeOf(IObase::Int2DataType(type));
}

int IOnElements(int rank,int *dims){
  return IObase::nElements(rank,dims);
}

int IOnBytes(int type,int rank,int *dims){
  return IObase::nBytes(IObase::Int2DataType(type),rank,dims);
}

int IOwrite(IOFile deviceID,int typeID,int rank,int *dims,void *data){
  IObase *dev=(IObase*)(deviceID);
  return dev->write(dev->Int2DataType(typeID),rank,dims,data);
}

int IOreadInfo(IOFile deviceID,int *typeID,int *rank,int *dims,int maxdims){
  IObase *dev=(IObase*)(deviceID);
  IObase::DataType tid;
  int r=dev->readInfo(tid,*rank,dims,maxdims);
  *typeID=tid;
  return r;
}

int IOread(IOFile deviceID,void *data){
  IObase *dev=(IObase*)(deviceID);
  return dev->read(data);
}

int IOseek(IOFile deviceID,int dataset_index){
  IObase *dev=(IObase*)(deviceID);
  return dev->seek(dataset_index);
}

int IOnDatasets(IOFile deviceID){
  IObase *dev=(IObase*)(deviceID);
  return dev->nDatasets();
}

int IOwriteAnnotation(IOFile deviceID,char *annotation){
  IObase *dev=(IObase*)(deviceID);
  //annotation[size]='\0'; // Yeah! its unsafe I know..!! (but it works!)
  return dev->writeAnnotation(annotation);
}
 
int IOreadAnnotation(IOFile deviceID,int index,char *annotation,int maxsize){
  IObase *dev=(IObase*)(deviceID);
  return dev->readAnnotation(index,annotation,maxsize);
}

int IOreadAnnotationInfo(IOFile deviceID,int index,int *size){
  IObase *dev=(IObase*)(deviceID);
  return dev->readAnnotationInfo(index,*size);
}

int IOnAnnotations(IOFile deviceID){
  IObase *dev=(IObase*)(deviceID);
  return dev->nAnnotations();
}

int IOwriteAttribute(IOFile deviceID,char *name,int type,Long length,void *data){
  IObase *dev=(IObase*)(deviceID);
  IObase::DataType typeID = IObase::Int2DataType(type);
  return dev->writeAttribute(name,typeID,length,data);
}

int IOreadIndexedAttributeInfo(IOFile deviceID,int number,char *name,
			int *type,Long *nelem,int maxnamelen){
  IObase *dev=(IObase*)(deviceID);
  IObase::DataType typeID;
  int i=dev->readAttributeInfo(number,name,typeID,*nelem,maxnamelen);
  *type=typeID; // convert from enum
  return i;
}

int IOreadAttributeInfo(IOFile deviceID,char *name,int *type,Long *nelem){
  IObase *dev=(IObase*)(deviceID);
  IObase::DataType typeID;
  int i=dev->readAttributeInfo(name,typeID,*nelem);
  *type=typeID; // convert from enum
  return i;
}  

int IOreadAttribute(IOFile deviceID,int number,void *data){
  IObase *dev=(IObase*)(deviceID);
  return dev->readAttribute(number,data);
}  

int IOnAttributes(IOFile deviceID){
  IObase *dev=(IObase*)(deviceID);
  return dev->nAttributes();
}

//=============C Chunking interface---------------
int IOreserveChunk(IOFile deviceID,int typeID,int rank,int *dims){
  IObase *dev=(IObase*)(deviceID);
  return dev->reserveChunk(dev->Int2DataType(typeID),rank,dims);
}

int IOwriteChunk(IOFile deviceID,int *chunkdims,int *chunkorigin,void *data){
  IObase *dev=(IObase*)(deviceID);
  return dev->writeChunk(chunkdims,chunkorigin,data);
}

int IOreadChunk(IOFile deviceID,int *chunkdims,int *chunkorigin,void *data){
  IObase *dev=(IObase*)(deviceID);
  return dev->readChunk(chunkdims,chunkorigin,data);
}

int IOwriteStream(IOFile deviceID,void *data,int length){
  IObase *dev=(IObase*)(deviceID);
  return dev->writeStream(data,length);
}

int IOreadStream(IOFile deviceID,void *data,int length){
  IObase *dev=(IObase*)(deviceID);
  return dev->readStream(data,length);
}
int IOpause(IOFile deviceID){
  IObase *dev=(IObase*)(deviceID);
  return dev->pause();
}

int IOresume(IOFile deviceID){
  IObase *dev=(IObase*)(deviceID);
  return dev->resume();
}
