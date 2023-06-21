#include <stdio.h>
#include <string.h>
#include "Reader.hh"

AttribInfoGroup::AttribInfoGroup(IOdataset *ds):dataset(ds){
  nilattrib.name[0]='\0';
  nilattrib.datatype=IObase::Error;
  nilattrib.nelements=0;
}

AttributeInfo &AttribInfoGroup::operator[](char *name){
  for(int i=0;i<attribs.getSize();i++){
    if(!strcmp((*this)[i].name,name)){
      return attribs[i];
    }
  }
  return nilattrib;
}

AttribInfoGroup::AttribInfoGroup(AttribInfoGroup &src){
  //printf("AttribInfoGroup: copy constructor nattribs=%u\n",src.attribs.getSize());
  nilattrib=src.nilattrib;
  attribs=src.attribs;
  dataset=src.dataset;
}

AttribInfoGroup &AttribInfoGroup::operator=(AttribInfoGroup &src){
  //puts("AttribInfo = operator");
  if(&src==this) return *this;
  //printf("\tAttribInfoGroup: nattribs=%u\n",src.attribs.getSize());
  nilattrib=src.nilattrib;
  attribs=src.attribs;
  dataset=src.dataset;
  return *this;
}

int AttributeInfo::read(void *data){
  return dataset->readAttribute(index,data);
}

void IOdataset::update(){
  int i;
  // read values
  file.seek(index);
  file.readInfo(datatype,rank,dims);
  typesize=IObase::sizeOf(datatype);
  for(i=0,nelements=1;i<rank;i++) nelements*=dims[i];
  nbytes=nelements*typesize;
  annotation.setSize(nannotations=file.nAnnotations());
  attribute.setSize(nattributes=file.nAttributes());
  for(i=0;i<file.nAnnotations();i++){
    file.readAnnotationInfo(i,annotation[i].length);
  }
  for(i=0;i<file.nAttributes();i++){
    file.readAttributeInfo(i,
			   attribute[i].name,
			   attribute[i].datatype,
			   attribute[i].nelements,
			   127);
  }
}

IOdataset::IOdataset(IObase &infile,int idx):file(infile),index(idx),attribute(this){
  //puts("IOdataset constructor");
  update();
  // printf("\tIOdataset nattribs=%u:%u\n",attribute.getSize(),nattributes);
}

IOdataset::IOdataset(IOdataset &src):file(src.file),attribute(this){ 
  //printf("\tIOdataset: copy constructor start nattribs=%u:%u\n",src.attribute.getSize(),src.nattributes);
  index=src.index;
  datatype=src.datatype;
  rank=src.rank;
  for(int i=0;i<rank;i++) dims[i]=src.dims[i];
  nannotations=src.nannotations;
  annotation=src.annotation;
  nattributes=src.nattributes;
  attribute=src.attribute;
  nbytes=src.nbytes;
  nelements=src.nelements;
  typesize=src.typesize;
  // printf("\tIOdataset: copy constructor end nattribs=%u:%u\n",attribute.getSize(),nattributes);
}

IOdataset &IOdataset::operator=(IOdataset &src){
  //puts("IOdataset = operator");
  if(&src==this) return *this;
  //printf("\tIOdataset: = operator start nattribs=%u:%u\n",src.attribute.getSize(),src.nattributes);
  file=src.file;
  index=src.index;
  datatype=src.datatype;
  rank=src.rank;
  for(int i=0;i<rank;i++) dims[i]=src.dims[i];
  nannotations=src.nannotations;
  annotation=src.annotation;
  nattributes=src.nattributes;
  attribute=src.attribute;
  nbytes=src.nbytes;
  nelements=src.nelements;
  typesize=src.typesize;
  //printf("\tIOdataset: = operator end nattribs=%u:%u\n",attribute.getSize(),nattributes);
  return *this;
}

int IOdataset::readAttribute(int i,void *data){
  file.seek(index);
  file.readInfo(datatype,rank,dims);
  return file.readAttribute(i,data);
}

int IOdataset::attributeIndex(char *name){
  for(int i=0;i<attribute.getSize();i++){
    if(!strcmp(attribute[i].name,name)){
      return i;
    }
  }
  return -1; // failed
}

int IOdataset::read(void *data){
  file.seek(index);
  file.readInfo(datatype,rank,dims);
  return file.read(data);
}


