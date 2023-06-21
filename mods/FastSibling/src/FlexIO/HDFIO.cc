#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <hdf.h>
#include "HDFIO.hh"
#include <mfhdf.h>

int32 HDFIO::DataType2HDF(IObase::DataType nt){
  switch(nt){
  case Int8:
    return DFNT_INT8; // means data
  case Char8:  // distinct from INT8..
  	return DFNT_CHAR8; // means string
  case Float32:
    return DFNT_FLOAT32;
  case Float64:
    return DFNT_FLOAT64;
  case Int32:
    return DFNT_INT32;
  case Int64:
    return DFNT_INT64;
  case Int16:
    return DFNT_INT16;
  case uInt8:
    return DFNT_UINT8;
  case uInt16:
    return DFNT_UINT16;
  case uInt32:
    return DFNT_UINT32;
  case uInt64:
    return DFNT_UINT64;
  case Char16: // unicode character type
  	return DFNT_CHAR16;
  }
  printf("HDFIO::HDF2DataType(): Don't recognize type %d\n",(int)nt);
  return -1;
}

IObase::DataType HDFIO::HDF2DataType(int32 nt){ 
  switch(nt){
  case DFNT_INT8:
  	return Int8;
  case DFNT_CHAR8:
    return Char8;
  case DFNT_FLOAT32:
    return Float32;
  case DFNT_FLOAT64:
    return Float64;
  case DFNT_INT32:
    return Int32;
  case DFNT_INT64:
    return Int64;
  case DFNT_INT16:
    return Int16;
  case DFNT_UINT8:
  case DFNT_UCHAR8:
    return uInt8;
  case DFNT_UINT16:
    return uInt16;
  case DFNT_UINT32:
    return uInt32;
  case DFNT_UINT64:
    return uInt64;
  case DFNT_CHAR16:
  	return Char16;
  }
  printf("HDFIO::HDF2DataType(): Don't recognize type %d\n",(int)nt);
  return Error;
}

void HDFIO::create(int rank,CONST int *dims,IObase::DataType nt){
  int32 hdims[5];
  if(sid>=0) SDendaccess(sid);
  nitems++;
  index=nitems-1;
  for(int i=0;i<rank;i++) hdims[i]=dims[i];
  sid=SDcreate(fid,"",DataType2HDF(nt),rank,hdims);
}

void HDFIO::select(int i){
  if(index==i && sid>=0) 
    return;
  if(index>=nitems) index=nitems-1;
  if(index<0) index=0;
  if(sid>=0) SDendaccess(sid);
  sid=SDselect(fid,i);
  index=i;
}
void HDFIO::endaccess(){
  if(sid<0) return;
  SDendaccess(sid);
  sid=-1;
}

HDFIO::HDFIO(CONST char *fname,AccessMode access):IObase(fname,access),sid(-1),fid(-1),hasread(0){
  switch(accessmode){
  case Read:
    fid=SDstart(filename,DFACC_RDONLY);
    break;
  case Write:
    fid=SDstart(filename,DFACC_CREATE);
    break;
  case Append:
    fid=SDstart(filename,DFACC_RDWR);
    break;
  default:
	puts("HDFIO: constructor... invalid accessmode");
	break;
  }
}

HDFIO::~HDFIO(){
  // printf("Destroying HDF file fid=%u, sid=%u\n",fid,sid);
  endaccess();
  if(fid>=0) SDend(fid);
  // SDend(fid);
}

int HDFIO::isValid() { if(fid>=0) return 1; else return 0; }

int HDFIO::write(IObase::DataType typeID,int rank,CONST int *dims,void *data){
  int32 origin[5]={0,0,0,0,0};
  int32 stride[5]={1,1,1,1,1}; // kludge... we'll fix later
  int32 hdims[5]={0,0,0,0,0};
  hasread=0;
  for(int i=0;i<rank;i++) hdims[i]=dims[i];
  create(rank,dims,typeID);
  //printf("write: sdsid=%u index=%u\n",sid,index);
  current_rank=rank;
  return (int)SDwritedata(sid,origin,stride,hdims,data);
}

int HDFIO::readInfo(char *name,IObase::DataType &typeID,int &rank,int *dims,int maxdims){
  int32 nt,nattr,hrank;
  int32 hdims[5];
  if(hasread)
    select(index+1);
  else
    select((int)index);
  SDgetinfo(sid,name,&hrank,hdims,&nt,&nattr);
  rank=(int)hrank;
  for(int i=0;i<rank && i<maxdims;i++) dims[i]=(int)(hdims[i]);
  typeID=HDF2DataType(nt);
  hasread=1;
  return 1;
}

int HDFIO::readInfo(IObase::DataType &typeID,int &rank,int *dims,int maxdims){
  int32 nt,nattr,hrank;
  int32 hdims[5];
  char name[128];
  if(hasread)
    select(index+1);
  else
    select((int)index);
  SDgetinfo(sid,name,&hrank,hdims,&nt,&nattr);
  rank=(int)hrank;
  for(int i=0;i<rank && i<maxdims;i++) dims[i]=(int)(hdims[i]);
  typeID=HDF2DataType(nt);
  // index++?
  hasread=1;
  return 1;
}

int HDFIO::read(void *data){
  int32 origin[5]={0,0,0,0,0};
  int32 stride[5]={1,1,1,1,1}; // kludge... we'll fix later
  int32 rank,dims[5],nt,natt;
  char name[128];
  select(index); // make certain its selected
  SDgetinfo(sid,name,&rank,dims,&nt,&natt);
  hasread=1;
  return (int)SDreaddata(sid,origin,stride,dims,data);
}

int HDFIO::seek(int i) { select((int)i); hasread=0; return index; }

int HDFIO::nDatasets(){ // not completely correct due to coordvar's
  // must scan for coordvar's and eliminate them from the count.
  int32 ndatasets,nattribs;
  SDfileinfo(fid,&ndatasets,&nattribs);
  return (int)ndatasets; //?
}

int HDFIO::writeAnnotation(CONST char *annotation){
  select((int)index); // select if not already selected
  int32 ref=SDidtoref(sid);
  return (int)DFANputlabel(filename,DFTAG_NDG,ref,(char*)annotation);
}

int HDFIO::readAnnotationInfo(int number,int &length){
  select(index);
  int32 ref=SDidtoref(sid);
  length=(int)DFANgetlablen(filename,DFTAG_NDG,ref);
  return length;
}

int HDFIO::readAnnotation(int number,char *annotation,int maxlen){
  // number=0; // How do I get the number of annotations availible?
  // use get lablist to get list of tags
  number=0; // number is ALWAYS 0 for hdf files
  select(index);
  int32 ref=SDidtoref(sid);
  return (int)DFANgetlabel(filename,DFTAG_NDG,ref,annotation,maxlen);
}


int HDFIO::nAnnotations(){
  select(index);
  int32 ref=SDidtoref(sid);
  if(DFANgetlablen(filename,DFTAG_NDG,ref)<=0) return 0; // no labels found
  return 1; // always 1 annotation per object limit for HDF is appears 
}

int HDFIO::writeAttribute(CONST char *name,IObase::DataType typeID,Long length,void *data){
  select(index); // select if not already selected
  //printf("write attrib: sdsid=%u index=%u\n",sid,index);
  return (int)SDsetattr(sid,(char*)name,DataType2HDF(typeID),(int32)length,data);
}

int HDFIO::readAttributeInfo(int number,char *name,IObase::DataType &typeID,Long &nelem,int maxnamelen){
  select(index);
  int32 numbertype;
  int32 nelements;
  int sz= SDattrinfo(sid,number,name,&numbertype,&nelements);
  typeID=HDF2DataType(numbertype);
  nelem=(int)nelements;
  return sz;
}

int HDFIO::readAttributeInfo(CONST char *name,IObase::DataType &typeID,Long &nelem){
  char fakename[128];
  select(index);
  int number=(int)SDfindattr(sid,(char*)name);
  if(number>0)
    return readAttributeInfo(number,fakename,typeID,nelem);
  else return -1;
}

int HDFIO::readAttribute(int number,void *data){
  select(index);
  return (int)SDreadattr(sid,number,data);
}

int HDFIO::nAttributes(){
  int32 nt,nattr,rank,dims[5];
  char name[128];
  select(index);
  SDgetinfo(sid,name,&rank,dims,&nt,&nattr);
  return (int)nattr;
}
//================Chunking Interface-----------------------
int HDFIO::reserveChunk(IObase::DataType typeID,int rank,CONST int *dims){
  //int32 origin[5]={0,0,0,0,0};
  //int32 stride[5]={1,1,1,1,1}; // kludge... we'll fix later
  //int32 hdims[5]={0,0,0,0,0};
  hasread=0;
  for(int i=0;i<rank;i++) chunkdims[i]=dims[i];
  create(rank,dims,typeID);
  current_rank=rank;
  return 1;
}

int HDFIO::writeChunk(CONST int *dims,CONST int *origin,void *data){
  int32 horigin[5]={0,0,0,0,0};
  int32 stride[5]={1,1,1,1,1}; // kludge... we'll fix later
  int32 hdims[5]={0,0,0,0,0};
  int32 rank = current_rank;
  for(int i=0;i<rank;i++) { hdims[i]=dims[i]; horigin[i]=origin[i]; }
  return (int)SDwritedata(sid,horigin,stride,hdims,data);
}


int HDFIO::readChunk(CONST int *dims,CONST int *origin,void *data){
  int32 horigin[5]={0,0,0,0,0};
  int32 stride[5]={1,1,1,1,1}; // kludge... we'll fix later
  int32 rank,nt,natt,hdims[5]={0,0,0,0,0};
  char name[128];
  select(index);
  SDgetinfo(sid,name,&rank,hdims,&nt,&natt);
  for(int i=0;i<rank;i++) {hdims[i]=dims[i]; horigin[i]=origin[i];}
  return (int)SDreaddata(sid,horigin,stride,hdims,data);
}
int HDFIO::isCoord() { // for HDFIO only!!
  select(index); // make sure it has been selected
  return SDiscoordvar(sid);
}

int HDFIO::readDimInfo(int dimnumber, char *name, IObase::DataType &datatype, int &length){
  int32 len,nt,nattribs;
  length=-1; // initial value
  // printf("HDFIO:readDimInfo\n");
  if(sid<0) return -1;
  int32 dim_id = SDgetdimid(sid,dimnumber);
  if(dim_id<0) return -1;
  SDdiminfo(dim_id,name,&len,&nt,&nattribs);
  //printf("HDFIO: dim_id=%u name=%s len=%u nt=%u nattribs=%u\n",
  //	  dim_id,name,len,nt,nattribs);
  length=len; 
  if(!nt) {
    datatype=IObase::Float32;
  }
  else {// compute custom datatype
    datatype = HDF2DataType(nt);
  }
  return length;
}

int HDFIO::readDim(int dimnumber, void *dim){
  if(sid<0) return -1;
  int32 dim_id = SDgetdimid(sid,dimnumber);
  if(dim_id<0) return -1;
  SDgetdimscale(dim_id,dim);
  return 1;
}

int HDFIO::writeDim(int dimnumber,IObase::DataType datatype,int length,void *dim){
  if(sid<0) return -1;
  int32 dim_id = SDgetdimid(sid,dimnumber);
  if(dim_id<0) return -1;
  SDsetdimscale(dim_id,length,DataType2HDF(datatype),dim);
  return 1;
}

int HDFIO::writeDimName(int dimnumber,CONST char *name){
  if(sid<0) return -1;
  int32 dim_id = SDgetdimid(sid,dimnumber);
  if(dim_id<0) return -1;
  SDsetdimname(dim_id,(char *)name);
  return 1;
}

//===============F77 Interface
Long8 f_hdf_open (char *file,char *accessname,int flen,int alen){
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
  IObase *fid=new HDFIO(file,mode);
  if(fid->isValid()) 
    return (Long8)fid;
  else
    delete fid; // file open failed
  return 0;
}

Long8 f_hdf_openr (char *file,int flen){
  file[flen]='\0'; // null terminate
  return (Long8)(new HDFIO(file,IObase::Read));
}

Long8 f_hdf_openw (char *file,int flen){
  file[flen]='\0'; // null terminate
  return (Long8)(new HDFIO(file,IObase::Create));
}

Long8 f_hdf_opena (char *file,int flen){
  file[flen]='\0'; // null terminate
  return (Long8)(new HDFIO(file,IObase::Append));
}

IOFile HDFIOopen (char *file,char *accessname){
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
	  !strcmp(accessname,"wb"))
    mode = IObase::Write;
  else if(!strcmp(accessname,"w+") ||
	  !strcmp(accessname,"w+b") ||
	  !strcmp(accessname,"wb+"))
    mode=IObase::Append;
  else{
    fprintf(stderr,"IEEEopen(): Error unknown option [%s] to open file %s\n",
	    accessname,file);
    return 0;
  }
  IObase *fid=new HDFIO(file,mode);
  if(fid->isValid())
    return (IOFile)fid;
  else
    delete fid; // file open failed
  return 0; // unknown option
}

IOFile HDFIOopenRead (char *file){
  return (IOFile)(new HDFIO(file,IObase::Read));
}

IOFile HDFIOopenWrite (char *file){
  return (IOFile)(new HDFIO(file,IObase::Write));
}

IOFile HDFIOopenAppend (char *file){
  return (IOFile)(new HDFIO(file,IObase::Append));
}

