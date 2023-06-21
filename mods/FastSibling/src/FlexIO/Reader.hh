#ifndef __READER_HH_
#define __READER_HH_

#include "IO.hh"
#include "FlexArrayTmpl.H"

struct AnnotationInfo{
  int length;
};

class IOdataset;

struct AttributeInfo{
  Long nelements;
  IObase::DataType datatype;
  char name[128];
  IOdataset *dataset;
  int index;
  int read(void *data);
};

class AttribInfoGroup {
  FlexArray<AttributeInfo> attribs;
  IOdataset *dataset;
public:
  AttributeInfo nilattrib;
  AttribInfoGroup(IOdataset *ds);
  ~AttribInfoGroup(){/*puts("Delete AttribInfoGroup");*/ }
  inline void setSize(int sz) {
    attribs.setSize(sz);
  }
  inline int getSize() { 
    return attribs.getSize(); 
  }
  inline AttributeInfo &operator[](int index){
    attribs[index].index=index;
    attribs[index].dataset = dataset;
    return attribs[index];
  }
  AttributeInfo &operator[](char *name);
  AttribInfoGroup(AttribInfoGroup &src);
  AttribInfoGroup &operator=(AttribInfoGroup &src);
};

class IOdataset {
  IObase &file;
  void update();
public:
  int index;
  int rank; // should put in constructor for const
  int dims[5];
  long nbytes;
  long nelements;
  IObase::DataType datatype;
  int typesize;
  

  FlexArray<AnnotationInfo> annotation;
  int nannotations;
  //FlexArray<AttributeInfo> attribute;
  AttribInfoGroup attribute;
  int nattributes;

  IOdataset(IObase &infile,int idx);
  IOdataset(IOdataset &src);
  ~IOdataset(){/*puts("delete IOdataset");*/}
  IOdataset &operator=(IOdataset &src);
  int readAttribute(int index,void *data);
  int attributeIndex(char *name);
  int read(void *data);
}; 

class Reader {
  IObase &file;
public:
  int ndatasets;
  Reader(IObase &infile):file(infile){
    ndatasets=file.nDatasets();
  }
  ~Reader(){/*puts("Delete Reader");*/}
  inline IOdataset operator[](int index){
    return IOdataset(file,index);
  }
};

#endif
