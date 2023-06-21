#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "sds.h"

 /*@@
   @file      sds.c
   @date      Wed Feb 26 20:05:00 1996
   @author   John Shalf 
   @desc 
   The SDS routines are a wrapper for the HDF SD interface
   that allows you to refer to files by name instead of
   explicitly opening them and refering to them by fileID.
   This file contain all of the source code that implements
   these wrappers.
   @enddesc 
   @includes sds.h
   @seeheader sds.h
   @version 0.9a
 @@*/

static int32 DefaultNumberType=DFNT_FLOAT32;
static int32 nfiles=0,maxfiles=54;

/* Need to set the default NumberType!!!! */

typedef struct fileID {
  char componentname[128];
  char filename[256];
  char sds_name[MAX_NC_NAME];  /* to keep track of SDS name internally */
  int32 id,sds_id,sdsnumber,nsds,nattrib,numbertype;
  int32 accessflags;
  clock_t lastaccess; /* in clock ticks.  Uses clock().  CLOCKS_PER_SEC */
  struct fileID *next;
}fileID;  

static fileID *files=NULL;

 /*@@
   @routine    deactivateOldFiles
   @date       Wed Feb 28 15:37:43 1996
   @desc 
	Finds the filesID's with longest time since last reference
	and deactivates the worst 1/8 of them to free up some
	file descriptors.
   @enddesc 
   @seeroutine SDSdeactivate SDSreactivate getFileID
   @calledby   NewFileID

   @history 
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf 
   @hdesc Initial Version   
   @endhistory 
 @@*/
static
void deactivateOldFiles() 
{
  fileID *tmpid;
  int i,ntimes;
  clock_t *worst_times;
  ntimes=maxfiles>>3; /* Ok, divide by 8 */
  ntimes=ntimes?ntimes:1; /* make sure its >0 */
  /* It was the best of times... it was the worst of times... */
  worst_times=(clock_t *)calloc(ntimes,sizeof(clock_t));
  /* first search for the top-1/8 longest times since last access */
  for(tmpid=files;tmpid;tmpid=tmpid->next){
    clock_t t=tmpid->lastaccess;
    if(tmpid->id<0) continue;
    for(i=0;i<ntimes;i++){	
      if(t<=worst_times[i]){ /* swap */
	clock_t tmp=worst_times[i];
	worst_times[i]=t;
	t=tmp; /* now looking for second-worst time */
      }
    }
  }
  /* next deactivate those top-10% */
  for(tmpid=files;tmpid;tmpid=tmpid->next)
    if(tmpid->id && tmpid->lastaccess<=worst_times[ntimes-1]) 
      SDSdeactivate(tmpid->filename); /* deactivate it */
  free(worst_times);
}

 /*@@
   @routine    NewFileID
   @date       Wed Feb 28 15:37:43 1996
   @desc 
	Creates a new fileID structure, opens the file and
	initializes all of the structure elements. 
   @enddesc 
   @calledby   SDSopen getFileID
   @calls deactivateOldFiles

   @par     name
   @pdesc   Name of the HDF file to open or create
   @ptype   char*
   @endpar 

   @par     access
   @pdesc   Access mode for the file (uses ANSI C access modes)
   @ptype   char*
   @pvalues "r","w","rw","w+","ra"
   @endpar 

   @history 
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf 
   @hdesc Initial Version   
   @endhistory 
 @@*/
fileID *NewFileID(char *name,char *access) {
  fileID *id=(fileID *)calloc(sizeof(fileID),1);
  strcpy(id->componentname,name);
  strcpy(id->filename,name);
  /* make certain that there are enough file descriptors around */
  if(nfiles==maxfiles) /* gotta deactivate a few */
    deactivateOldFiles();
  
  /* open SDS */
  if(*access=='r'){
    if(access[1]=='w'){ 
      id->id=SDstart(id->filename,id->accessflags=DFACC_RDWR);
      if(id->id<0){
	/* kludge around HDF bug */
	id->id=SDstart(id->filename,id->accessflags=DFACC_CREATE);
      }
      SDfileinfo(id->id,&(id->nsds),&(id->nattrib));
      id->sdsnumber=0; /* end of file */
      if(id->accessflags!=DFACC_CREATE)
	id->sds_id=SDselect(id->id,0); /* goto first */
      else
	id->sds_id=-1;
    }
    else if (access[1]=='a' || access[1]=='+') {
      id->id=SDstart(id->filename,id->accessflags=DFACC_RDWR);
      nfiles++; 
      if(id->id<0){
	/* kludge around HDF bug */
	id->id=SDstart(id->filename,id->accessflags=DFACC_CREATE);
      }
      SDfileinfo(id->id,&(id->nsds),&(id->nattrib));
      id->sdsnumber=id->nsds; /* end of file */
      if(id->accessflags!=DFACC_CREATE)
	id->sds_id=SDselect(id->id,id->sdsnumber); /* goto last */
      else
	id->sds_id=-1;
    }
    else {
      id->id=SDstart(id->filename,id->accessflags=DFACC_RDONLY);
      SDfileinfo(id->id,&(id->nsds),&(id->nattrib));
      id->sdsnumber=0; /* end of file */
      id->sds_id=SDselect(id->id,id->sdsnumber); /* means, need to create */
    }
  }
  else if(*access=='w'){
    if(access[1]=='\0'){
      id->id=SDstart(id->filename,id->accessflags=DFACC_CREATE);
      id->sdsnumber=id->nsds=0; /* end of file */
      id->sds_id=-1; /* means, need to create */
    }
    else{
      id->id=SDstart(id->filename,id->accessflags=DFACC_ALL);
      SDfileinfo(id->id,&(id->nsds),&(id->nattrib));
      id->sdsnumber=id->nsds;
      id->sds_id=SDselect(id->id,id->sdsnumber); /* goto last */
    }
  }
  if(id->id >= 0)
    nfiles++; /* it has been opened successfully */
  id->lastaccess=clock();
  id->numbertype=DefaultNumberType;
  id->next=files;
  files=id;
  return id;
}

 /*@@
   @routine    findIDbyComponent
   @date       Wed Feb 28 20:37:43 1996
   @desc
	Finds a File ID given the name of the component
	that the file contains.  This is to support a
	user interface style in which the componentnames
	can be used as alias for the actual filename.
   @enddesc

   @par     name
   @pdesc   Name of the component.
   @ptype   char*
   @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/

fileID *findIDbyComponent(name) 
char *name;
{
  fileID *id;
  for(id=files;id!=NULL;id=id->next)
    if(!strcmp(name,id->componentname))
      break;
  return id;
}

 /*@@
   @routine    findIDbyFilename
   @date       Wed Feb 28 20:37:43 1996
   @desc
        Finds a File ID given the name of the component
        that the file contains.  This is to support a
        user interface style in which the componentnames
        can be used as alias for the actual filename.
   @enddesc
   @calledby   SDSopen SDSclose SDSdeactivate SDSreactivate getFileID

   @par     name
   @pdesc   Name of the HDF file.
   @ptype   char*
   @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
fileID *findIDbyFilename(name) 
char *name;
{
  fileID *id;
  for(id=files;id!=NULL;id=id->next)
    if(!strcmp(name,id->filename))
      break;
  return id;
}

 /*@@
   @routine    removeFileID
   @date       Wed Feb 28 20:37:43 1996
   @desc
	Remove a file ID reference and free associated memory.
	It does not close the HDF file if is left open though.
   @enddesc
   @calledby   SDSclose getFileID

   @par     name
   @pdesc   Name of the file.
   @ptype   char*
   @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
void removeFileID(filename)
char *filename;
{
  fileID *id,*lastid;

  for(lastid=0,id=files;id;lastid=id,id=id->next) {
    if(!strcmp(id->filename,filename)){
      if(!lastid){
	files=id->next;
      }
      else {
	lastid->next=id->next;
      }
      free(id);
    }
  }
}

/*@@
 @routine    SDScomputeSize
   @date       Wed Feb 28 15:37:43 1996
   @desc
        Computes the Size in bytes required to store a dataset
	based on its dimensions and datatype.
   @enddesc
   @seeroutine SDSgetDims SDSreadData

   @par rank
   @pdesc Number of dimensions in the dataset
   @ptype int32
   @endpar

   @par dims
   @pdesc The size of each dimension in the dataset (in elements).
   @ptype int32*
   @endpar

   @par datatype
   @pdesc The HDF datatype for each element in the dataset.
   @ptype int32*
   @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
@@*/
long SDScomputeSize(rank,dims,datatype)
int32 rank;
int32 *dims;
int32 datatype;
{
  register int i,n;
  register long nbytes;
  for(nbytes=1,n=rank,i=0;i<n;i++){
    nbytes*=dims[i];
  }
  nbytes*=DFKNTsize(datatype);
  return nbytes;
}

 /*@@
   @routine    SDScloseAll
   @date       Wed Feb 28 15:37:43 1996
   @desc
        Closes and all SDS files and deallocates all
	in-memory datastructures used by the SDS routines.
   @enddesc
   @seeroutine SDSclose

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
void SDScloseAll()
{
  fileID *id,*nextid;
  for(id=files;id;id=nextid) {
    nextid=id->next;
    if(id->sds_id>=0)
      SDendaccess(id->sds_id); 
    id->sds_id=-1;
    if(id->id>=0){
      nfiles--;
      SDend(id->id);
    }
    id->id=-1;
    free(id);
  }
  files=NULL;
}

 /*@@
   @routine    SDSdeactivate
   @date       Wed Feb 28 15:37:43 1996
   @desc
        Closes and HDF file but preserves its access state so
	that it can be reactivated later.  Used to free up 
	file descriptors temporarily when more files are 
	active than there are file descriptors availible.
	Used for "file descriptor juggling".
   @enddesc
   @seeroutine SDSreactivate
   @calledby   getFileID
   @calls      findIDbyFilename

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/

/* make file handle inactive (use if shortage of descriptors) */
void SDSdeactivate(filename) 
char *filename;
{
  fileID *id=findIDbyFilename(filename);
  if(id){
    if(id->sds_id)
      SDendaccess(id->sds_id);
    SDend(id->id);
    id->id=-1;
    nfiles--;
    if(id->accessflags==DFACC_CREATE) 
      id->accessflags=DFACC_RDWR; /* make sure it doesn't purge the file
				     on reactivation */
  }
}


/*@@
   @routine    SDSreactivate
   @date       Wed Feb 28 15:37:43 1996
   @desc
        Reopens an HDF file closed by SDSdeactivate and
	returns it to its original state.
   @enddesc 
   @seeroutine SDSdeactivate
   @calledby   getFileID
   @calls      findIDbyFilename

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
void SDSreactivate(filename)
char *filename;
{
  fileID *id=findIDbyFilename(filename);
  if(id && id->id < 0){
    id->id=SDstart(filename,id->accessflags);
    if(id->id)
      nfiles++;
    if(id->sds_id>=0)
      id->sds_id=SDselect(id->id,id->sdsnumber);
  }
}

/*@@
   @routine    SDSclose
   @date       Wed Feb 28 15:37:43 1996
   @desc
        Closes an HDF file and removes its state description
	from memory.
   @enddesc 
   @seeroutine SDSopen SDSdeactivate SDSreactivate
   @calls      findIDbyFilename

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
void SDSclose(filename)
char *filename;
{
  fileID *id=findIDbyFilename(filename);
  if(id) {
    SDend(id->id);
    nfiles--;
    removeFileID(filename);
  }
}

 /*@@
   @routine    getFileID
   @date       Wed Feb 28 20:37:43 1996
   @desc
	Checks to see if the specified file is open.
	If it isn't already availbile, it opens it
	in read-write mode.  Otherwise return the existing
	fileID structure.
   @enddesc
   @calls findIDbyFilename
   @calledby   SDSopen SDSflush SDSpurge SDSseek SDSseekName SDSgetIndex
               SDSisCoord SDSgetFileInfo SDSgetDataName SDSgetNumDatasets
	       SDSgetDims SDSgetNT SDSsetNT SDSreadData SDSreadChunk
	       SDSallocateDataset SDSwriteData SDSwriteChunk
	       SDSaddAnnotation SDSgetAnnotationSize SDSgetAnnotation
	       SDSgetNumAttribs SDSwriteAttrib SDSreadAttrib 
	       SDSfindAttribInfo SDSgetAttribInfo SDSgetDataStrs
	       

   @par     name
   @pdesc   Name of the HDF file.
   @ptype   char*
   @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
fileID *getFileID(filename)
char *filename;
{
  fileID *id=findIDbyFilename(filename);
  if(!id)
    id=NewFileID(filename,"rw");
  else
    id->lastaccess=clock();
  if(id->id<0)
    SDSreactivate(filename); /* could be more efficient */
  return id;
}

/*@@
   @routine    SDSopen
   @date       Wed Feb 28 15:37:43 1996
   @desc
         Opens an HDF file in specified access mode.
   @enddesc 
   @seeroutine SDSclose SDSdeactivate SDSreactivate
   @calls      findIDbyFilename NewFileID

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar

   @par     access
   @pdesc   Access mode for the file (uses ANSI C access modes)
   @ptype   char*
   @pvalues "r","w","rw","w+","ra"
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
void SDSopen(filename,access)
char *filename,*access;
{
  fileID *id=findIDbyFilename(filename);
  if(!id)
    id=NewFileID(filename,access);
}

/*@@
   @routine    SDSflush
   @date       Wed Feb 28 15:37:43 1996
   @desc
         Flushes the HDF file so that if a crash occurs, 
	 the contents of the file will be preserved. (currently disabled)
   @enddesc 
   @seeroutine SDSclose SDSdeactivate SDSreactivate
   @calls      findIDbyFilename SDSdeactivate SDSreactivate

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @hdate Tues March 8 
   @hdesc Commented out deactivate/reactivate due to file access mode problems (must fix later)
   @endhistory
 @@*/
void SDSflush(filename)
char *filename;
{
  fileID *id=findIDbyFilename(filename);
  if(id && id->sds_id>=0)
    SDendaccess(id->sds_id);
  /*
  SDSdeactivate(filename);
  SDSreactivate(filename);
  */
}

/*@@
   @routine    SDSpurge
   @date       Wed Feb 28 15:37:43 1996
   @desc
         Purges an HDF file, removing all of its contents.
	 Can Use to garauntee that you are writing into a
	 fresh HDF file since the default way SDS opens 
	 a file is Read/Write/Append mode.
         Potentially dangerous routine...
   @enddesc 
   @calls      SDSclose SDSopen

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
void SDSpurge(filename)
char *filename;
{
  SDSclose(filename);
  SDSopen(filename,"w");
}

/*@@
   @routine    SDSseek
   @date       Wed Feb 28 15:37:43 1996
   @desc
         Allows random acces in HDF files.
	 Moves to the sds in the file specified by setnum.
   @enddesc 
   @calls      getFileID

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     setnum
   @pdesc   Index of the dataset in the HDf file
   @ptype   int
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
void SDSseek(filename,setnum) 
char *filename;
int setnum;
{
  fileID *id=getFileID(filename);
  id->sdsnumber=setnum;
  id->sds_id=SDselect(id->id,id->sdsnumber);
}

/*@@
   @routine    SDSseekName
   @date       Wed Feb 28 15:37:43 1996
   @desc
         Allows random acces in HDF files.
	 Moves to the sds in the file with the name specified by "dataname".
   @enddesc 
   @calls      getFileID

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     dataname
   @pdesc   Name of the dataset in the HDf file
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSseekName(filename,dataname)
char *filename,*dataname;
{
  fileID *id=getFileID(filename);
  id->sdsnumber=SDnametoindex(id->id,dataname);
  id->sds_id=SDselect(id->id,id->sdsnumber);
  return id->sdsnumber;
}

/*@@
   @routine    SDSgetNT
   @date       Wed Feb 28 15:37:43 1996
   @desc
           Gets the numbertype of the current SDS.
   @enddesc 
   @calls      getFileID

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
/* analogue to DFSDgetDims()
   it gets the AMR attributes of the next SDS */
int32 SDSgetNT(filename)
char *filename;
{
  int32 numbertype,attribs,rank,dims[MAX_VAR_DIMS];
  fileID *id=getFileID(filename);
  if(id->sds_id<0)
    id->sds_id=SDselect(id->id,id->sdsnumber);
  if(id->sds_id<0)
    return -1;  
  SDgetinfo(id->sds_id,NULL,&rank,dims,&numbertype,&attribs);
  id->numbertype=numbertype;
  return numbertype;
}

/*@@
   @routine    SDSisCoord
   @date       Wed Feb 28 15:37:43 1996
   @desc
           Returns 1 if the current dataset is a coordinate array
	   and false if it is not.  This is a direct passthru of
	   the function SDisCoordVar().  This function has become
	   necessary because in order to make HDF compatible with
	   NetCDF, the files the designers had to disable its ability
	   to tell the difference between datasets and stored coordinates
	   for the datasets.
   @enddesc 
   @calls      getFileID

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSisCoord(filename)
char *filename;
{
  fileID *id=getFileID(filename);
  if(id->sds_id<0)
    id->sds_id=SDselect(id->id,id->sdsnumber);
  if(id->sds_id<0)
    return -1;  
  return SDiscoordvar(id->sds_id);
}

/*@@
   @routine    SDSsetNT
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Sets the default numbertype for a particular
	  HDF file.  All further write operations will assume
	  this numbertype until SDSsetNT is used again to
	  select a new numbertype.<p>

	  The default numbertype is DFNT_FLOAT32 (32 bit float)	  
   @enddesc 
   @calls      getFileID

   @par     name
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     numbertype
   @pdesc   Numbertype for subsequent data writes
   @ptype   char*
   @pvalues DFNT_FLOAT32 (32-bit float), DFNT_FLOAT64 (64-bit float/double), 
            DFNT_INT8 (char), DFNT_INT16 (short int), DFNT_INT32,(int),
	    DFNT_UINT8 (unsigned char), DFNT_UNINT16 (unsigned short), 
	    DFNT_UINT32 (unsigned int),
	    DFNT_NATIVE (if bit-OR'ed with the above, uses machine native format)
	    
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
void SDSsetNT(filename,numbertype)
char *filename;
int32 numbertype;
{
   fileID *id=getFileID(filename);
   id->numbertype=numbertype;
}


/*@@
   @routine    SDSgetDims
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Reads data from an HDF file.
	  This must be preceded by an SDSgetDims() which
	  returns the dimensions of the data so that the user
	  can allocate the space to read in the data.
   @enddesc
   @seeroutine SDSreadData
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     name
   @pdesc   Name of the dataset (read from the file)
   @ptype   char*
   @endpar 

   @par     rank
   @pdesc   Number of dimension in the dataset
   @ptype   int32
   @endpar 

   @par     dims
   @pdesc   Dimensions of the dataset
   @ptype   int32*
   @endpar 

   @par     data
   @pdesc   Pointer to the actual dataset (preallocated by the user)
   @ptype   VOIDP
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
long SDSgetDims (filename,dataname,rank,dims)
char *filename;
char *dataname;
int32 *rank;
int32 *dims;
{
  int32 numbertype,attribs;
  long status;
  fileID *id=getFileID(filename);
  if(id->sds_id>=0)
    SDendaccess(id->sds_id);
  id->sds_id=SDselect(id->id,(id->sdsnumber)++);
  if(id->sds_id<0)
    return -1;  
  status=SDgetinfo(id->sds_id,dataname,rank,dims,&numbertype,&attribs);
  if(status>=0){
    id->numbertype=numbertype;
    SDScomputeSize(*rank,dims,numbertype);
  }
  else
    return -1;
}
/*@@
   @routine    SDSreadData
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Reads data from an HDF file.
	  This must be preceded by an SDSgetDims() which
	  returns the dimensions of the data so that the user
	  can allocate the space to read in the data.
   @enddesc
   @seeroutine SDSgetDims
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     dataname
   @pdesc   Name of the dataset (read from the file)
   @ptype   char*
   @endpar 

   @par     rank
   @pdesc   Number of dimension in the dataset
   @ptype   int32
   @endpar 

   @par     dims
   @pdesc   Dimensions of the dataset
   @ptype   int32*
   @endpar 

   @par     data
   @pdesc   Pointer to the actual dataset (preallocated by the user)
   @ptype   VOIDP
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int SDSreadData(filename,rank,dims,data)
char *filename;
int32 rank,*dims;
VOIDP data;
{
  static int32 originS[5]={0,0,0,0,0},*origin;
  int32 natt,nt;
  int status;
  fileID *id=getFileID(filename);
  if(id->sds_id<0)
    id->sds_id=SDselect(id->id,id->sdsnumber);

  if(rank > 5) {
    origin=(int32 *)calloc(sizeof(int),rank); /* origin is zeroed */
  }
  else 
    origin=originS;
  if(id->sds_id<0)
    id->sds_id=SDselect(id->id,id->sdsnumber);
  if(id->sds_id<0)
    return -1;
  /*SDgetinfo(id->sds_id,NULL,&rank,dims,&(id->numbertype),&(id->nattrib));*/
  status = SDreaddata(id->sds_id,origin,NULL,dims,data);
  if(origin!=originS) 
    free(origin);
  if(status>=0)
    return SDScomputeSize(rank,dims,id->numbertype);
  else
    return -1;
}



/*@@
   @routine    SDSreadChunk
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Reads a subsection (chunk) of data from a selected dataset.
	  You can specify any origin, stride, and region-
	  of-interest within the total size of the dataset.

	  A dataset must have been selected by an SDSgetDims() which
	  returns the total dimensions of the data.  Once selected
	  though, it can be followed by any number of reads.
   @enddesc
   @seeroutine SDSgetDims
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     rank
   @pdesc   Number of dimension in the dataset
   @ptype   int32
   @endpar 

   @par     dims
   @pdesc   Dimensions of the dataset to read
   @ptype   int32*
   @endpar 

   @par     origin
   @pdesc   Origin of the read within the dataset
   @ptype   int32*
   @endpar 

   @par     stride
   @pdesc   Stride of the read in each dimension
   @ptype   int32*
   @endpar 

   @par     data
   @pdesc   Pointer to the actual dataset (preallocated by the user)
   @ptype   VOIDP
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   need to eliminate chunk reading inefficiency.
   Basically, it needs to re-read the chunkname on every cycle
   which results in a lot of unneccessary seeking.  To
   optimize, will need to store the sds name internally and
   only get it again if its a different sds.  Either that
   or just separate readNamedChunk from this one.
   @hdate Tues Mar 8 @hauthor
   @hdesc Eliminated "name" from parameter list.
   Depends on SDSgetDims to read the name now.
   @endhistory
 @@*/
int SDSreadChunk(filename,rank,dims,
		      origin,stride,data)
char *filename;
int32 rank,*dims;
int32 *origin,*stride;
VOIDP data;
{
  static int32 originS[5]={0,0,0,0,0};
  int32 natt,nt;
  int status;
  char name[MAX_NC_NAME];
  fileID *id=getFileID(filename);
  if(!origin){
    if(rank > 5) {
      origin=(int32 *)calloc(sizeof(int),rank); /* origin is zeroed */
    }
    else 
      origin=originS;
  }
  if(id->sds_id<0){
    id->sds_id=SDselect(id->id,id->sdsnumber);
    if(id->sds_id<0)
      return -1;
    SDgetinfo(id->sds_id,name,&rank,dims,&(id->numbertype),&(id->nattrib));
  }
  status = SDreaddata(id->sds_id,origin,stride,dims,data);
  if(origin!=originS) 
    free(origin);
  return status;
}


/*@@
   @routine    SDSgetFileInfo
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Gets information relevant to the entire file.
	  Namely the number datasets contained in the file as
	  well as the number of attributes associated with the file.
   @enddesc
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     nsds
   @pdesc   Number of datasets (sds's) in the hdf file
   @ptype   int32*
   @endpar 

   @par     nattrib
   @pdesc   Number of netcdf attributes for the file
   @ptype   int32*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
void SDSgetFileInfo (filename,nsds,nattrib)
char *filename;
int32 *nsds;
int32 *nattrib;
{
  fileID *id=getFileID(filename);
  if(nsds)
    *nsds=id->nsds;
  if(nattrib)
    *nattrib=id->nattrib;
}

/*@@
   @routine    SDSgetDataName
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Gets the name of the current dataset
   @enddesc
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     dataname
   @pdesc   Name of the current dataset
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSgetDataName (filename,name)
char *filename;
char *name;
{
  int32 status;
  int32 rank,dims[MAX_NC_DIMS];
  fileID *id=getFileID(filename);
  if(id->sds_id<0)
    id->sds_id=SDselect(id->id,id->sdsnumber);
  if(id->sds_id<0)
    return -1;
  status=SDgetinfo(id->sds_id,name,rank,dims,
		   &(id->numbertype),
		   &(id->nattrib));
  return status;
}

/*@@
   @routine    SDSgetIndex
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Returns the index of the current dataset
   @enddesc
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32  SDSgetIndex (filename)
char *filename;
{
  fileID *id=getFileID(filename);
  return id->sdsnumber;
}

/*@@
   @routine    SDSgetDataStrs
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Retrieves the data strings associated with the current dataset.
   @enddesc
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     label
   @pdesc   Label of the current dataset
   @ptype   char*
   @endpar 

   @par     units
   @pdesc   Units of the current dataset (ie. meters, millibars, foot-pounds)
   @ptype   char*
   @endpar   

   @par     format
   @pdesc   Format of the current dataset
   @ptype   char*
   @endpar 

   @par     coordsys
   @pdesc   Coordinate system for the current dataset.
   @ptype   char*
   @endpar 

   @par     maxlen
   @pdesc   Maximum length of a string the variables can receive.
   @ptype   int32
   @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/

int32  SDSgetDataStrs (filename,label,units,format,coordsys,maxlen)
char *filename,*label,*units,*format,*coordsys;
int maxlen;
{
  fileID *id=getFileID(filename);
  if(id->sds_id<0)
    id->sds_id=SDselect(id->id,id->sdsnumber);
  if(id->sds_id<0)
    return -1;
  return SDgetdatastrs(id->id,label,units,format,coordsys,maxlen);
}

/*@@
   @routine    SDSgetNumDatasets
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Returns the number of SDS's in the current dataset
   @enddesc
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32  SDSgetNumDatasets (filename)
char *filename;
{
  fileID *id=getFileID(filename);
  return id->nsds;
}

/*@@
   @routine    SDSallocateDataset
   @date       Wed Feb 28 15:37:43 1996
   @desc
            Creates and reserves a dataset to receive data in chunks.
	    Use in conjunction with SDSwriteChunk.
   @enddesc
   @seeroutine SDSwriteChunk
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     dataname
   @pdesc   Name of the SDS
   @ptype   char*
   @endpar 

   @par     rank
   @pdesc   Number of dimensions for the dataset
   @ptype   char*
   @endpar 

   @par     dims
   @pdesc   The maximum dimension sizes for the dataset 
            (can be SD_UNLIMITED in any dimension to allow the
	     space availible to automatically grow as you write chunks)
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSallocateDataset(filename,dataname,rank,dims)
char *filename,*dataname;
int32 rank,*dims;
{
  static int32 originS[5]={0,0,0,0,0};
  int status;
  fileID *id=getFileID(filename);
  if(id->accessflags==DFACC_RDONLY)
    return -1; /* cant write to this file */
  if(id->sds_id>=0)
   SDendaccess(id->sds_id); /* may no longer be a valid operation */

  id->sds_id=SDcreate(id->id,dataname,id->numbertype,rank,dims);
  if(id->sds_id>=0)
    (id->nsds)++;

  return id->sdsnumber=SDreftoindex(id->id,SDidtoref(id->sds_id));
}

/*@@
   @routine    SDSwriteChunk
   @date       Wed Feb 28 15:37:43 1996
   @desc
            Writes a chunk of data to an sds with any
	    stride, origin, and dimension so long as the
	    chunk falls within the bounds of the dataset
	    reserved with SDSallocateDataset().  If any
	    of the dimensions were declared SD_UNLIMITED
	    for size, then the size of the chunk is only
	    constrained by system and file size limits.
   @enddesc
   @seeroutine SDSallocateDataset
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     rank
   @pdesc   Number of dimensions for the dataset
   @ptype   char*
   @endpar 

   @par     dims
   @pdesc   The maximum dimension sizes for the dataset 
            (can be SD_UNLIMITED in any dimension to allow the
	     space availible to automatically grow as you write chunks)
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int SDSwriteChunk (filename,rank,dims,
			origin,stride,data)
char *filename;
int32 rank,*dims;
int32 *origin,*stride;
VOIDP data;
{
  static int32 originS[5]={0,0,0,0,0};
  int status;
  fileID *id=getFileID(filename);
  if(id->accessflags==DFACC_RDONLY)
    return -1; /* cant write to this file */
  
  if(!origin){
    if(rank > 5) {
      origin=(int32 *)calloc(sizeof(int32),rank); /* origin is zeroed */
    }
    else 
      origin=originS;
  }
  if(id->sds_id<0)
	return -1; /* have to reserve the chunk first!! */
/*
  if(id->sds_id<0)  only create if it doesn't already exist 
    id->sds_id=SDcreate(id->id,id->dataname,id->numbertype,rank,dims);
*/
  if(id->sds_id>=0)
    (id->nsds)++;
  id->sdsnumber=SDreftoindex(id->id,SDidtoref(id->sds_id));
  status=SDwritedata(id->sds_id,origin,stride,dims,data);
  /* SDendaccess(id->sds_id); (might append annotation) */
  if(origin!=originS) 
    free(origin);
  return status;
}

/*@@
   @routine    SDSwriteData
   @date       Wed Feb 28 15:37:43 1996
   @desc
            Writes a chunk of data to an sds with any
	    stride, origin, and dimension so long as the
	    chunk falls within the bounds of the dataset
	    reserved with SDSallocateDataset().  If any
	    of the dimensions were declared SD_UNLIMITED
	    for size, then the size of the chunk is only
	    constrained by system and file size limits.
   @enddesc
   @seeroutine SDSreadData SDSsetNT()
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     dataname
   @pdesc   Name for the dataset
   @ptype   char*
   @endpar 

   @par     rank
   @pdesc   Number of dimensions for the dataset
   @ptype   char*
   @endpar 

   @par     dims
   @pdesc   The maximum dimension sizes for the dataset 
   @ptype   char*
   @endpar 

   @par     data
   @pdesc   Pointer to the dataset to write.  
             Default type is DFNT_FLOAT32 (32-bit float).
	     This default can be changed using SDSsetNT().
   @ptype   VOIDP
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int SDSwriteData (filename,dataname,rank,dims,data)
char *filename,*dataname;
int32 rank,*dims;
VOIDP data;
{
  static int32 originS[5]={0,0,0,0,0},*origin;
  int status;
  fileID *id=getFileID(filename);
  if(id->accessflags==DFACC_RDONLY)
    return -1; /* cant write to this file */
  if(id->sds_id>=0)
   SDendaccess(id->sds_id);
    
  if(rank > 5) {
    origin=(int32 *)calloc(sizeof(int32),rank); /* origin is zeroed */
  }
  else 
    origin=originS;

  id->sds_id=SDcreate(id->id,dataname,id->numbertype,rank,dims);
  if(id->sds_id>=0)
    (id->nsds)++;
  id->sdsnumber=SDreftoindex(id->id,SDidtoref(id->sds_id));
  status=SDwritedata(id->sds_id,origin,NULL,dims,data);
  /* SDendaccess(id->sds_id); (might append annotation) */
  if(origin!=originS) 
    free(origin);
  return status;
}

/*@@
   @routine    SDSgetAnnotationSize
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Returns the size of the current annotation attached to
	  the current dataset
   @enddesc
   @seeroutine SDSgetAnnotation
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSgetAnnotationSize(filename)
char *filename;
{
  int32 ref;
  fileID *id=getFileID(filename);
  if(id->accessflags==DFACC_RDONLY)
    return -1; /* cant write to this file */
  if(id->sds_id<0)
   id->sds_id=SDselect(id->id,id->sdsnumber);
  ref=SDidtoref(id->sds_id);
  return DFANgetlablen(id->filename,DFTAG_NDG,ref);
}

/*@@
   @routine    SDSgetAnnotation
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Returns the size of the current annotation attached to
	  the current dataset.
   @enddesc
   @seeroutine SDSgetAnnotationSize SDSaddAnnotation
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     annotation
   @pdesc   buffer to store the annotation in
   @ptype   char*
   @endpar 

   @par     maxlen
   @pdesc   Maximum length of the annotation (ie. the buffer size)
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSgetAnnotation(filename,annotation,maxlen)
char *filename;
char *annotation;
int32 maxlen;
{
  int32 ref;
  fileID *id=getFileID(filename);
  if(id->accessflags==DFACC_RDONLY)
    return -1; /* cant write to this file */
  if(id->sds_id<0)
   id->sds_id=SDselect(id->id,id->sdsnumber);
  ref=SDidtoref(id->sds_id);
  return DFANgetlabel(id->filename,DFTAG_NDG,ref,annotation,maxlen);
}

/*@@
   @routine    SDSaddAnnotation
   @date       Wed Feb 28 15:37:43 1996
   @desc
         Adds a text string annotation to the currently selected SDS
   @enddesc
   @seeroutine SDSgetAnnotation SDSgetAnnotationSize
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     annotation
   @pdesc   buffer to store the annotation in
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSaddAnnotation(filename,annotation)
char *filename;
char *annotation;
{
  int32 ref;
  fileID *id=getFileID(filename);
  if(id->accessflags==DFACC_RDONLY)
    return -1; /* cant write to this file */
  if(id->sds_id<0)
   id->sds_id=SDselect(id->id,id->sdsnumber);
  ref=SDidtoref(id->sds_id);
  return DFANputlabel(id->filename,DFTAG_NDG,ref,annotation);
}

/*@@
   @routine    SDSgetNumAttribs
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Returns the number of attributes associated with the
	  current dataset.  The attributes can be retreived by
	  name (if you already know which attributes you are
	  looking for) or by index (if you want to step through them
	  in order).  Returns an error if no SDS is currently selected.
   @enddesc
   @seeroutine SDSreadAttrib SDSfindAttribInfo SDSgetAttribInfo
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSgetNumAttribs(filename)
char *filename;
{
  int nattribs;
  fileID *id=getFileID(filename);
  if(id->sds_id>=0)
    return id->nattrib;
  else
    return -1; /* or should it return 0???? */
}

/*@@
   @routine    SDSwriteAttrib
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Writes an attribute associated with the current SDS
	  to the HDF file.  The attribute can be referred to
	  by its name or by its index when read back in.
   @enddesc
   @seeroutine SDSreadAttrib SDSfindAttribInfo SDSgetAttribInfo SDSgetNumAttribs
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     attribname
   @pdesc   Name of the attribute.
   @ptype   char*
   @endpar 

   @par     numbertype
   @pdesc   Numbertype of the attribute
   @pvalues DFNT_FLOAT32 (32-bit float), DFNT_FLOAT64 (64-bit float/double), 
            DFNT_INT8 (char), DFNT_INT16 (short int), DFNT_INT32,(int),
	    DFNT_UINT8 (unsigned char), DFNT_UNINT16 (unsigned short), 
	    DFNT_UINT32 (unsigned int),
	    DFNT_NATIVE (if bit-OR'ed with the above, uses machine native format)
   @ptype   int32
   @endpar 

   @par     nelements
   @pdesc   Number of elements (of the type specified by the numbertype param)
            to store in the attribute.
   @ptype   int32
   @endpar 

   @par     buffer
   @pdesc   Pointer to the attribute data.
   @ptype   VOIDP
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSwriteAttrib(filename,attribname,numbertype,nelements,buffer)
char *filename;
char *attribname;
int32 numbertype;
int32 nelements;
VOIDP buffer;
{
  int32 ref;
  fileID *id=getFileID(filename);
  if(id->accessflags==DFACC_RDONLY)
    return -1; /* cant write to this file */
  if(id->sds_id<0)
    return -1; /* needs to be existing writeable SDS ? */
  /*
   id->sds_id=SDselect(id->id,id->sdsnumber);
  if(id->sds_id<0)
    return -1;
  */
  return SDsetattr(id->sds_id,attribname,numbertype,nelements,buffer);
}
/*@@
   @routine    SDSreadAttrib
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Reads the attribute with the specified "attribname" associated 
	  with the current SDS to the HDF file.  If there is no attribute
	  by that name, it returns 0.  If an attribute is found by that
	  name, it returns the number of elements read.
	  The attribute can be referred to by its name or by its 
	  index when read back in.<p>

	  If you want to pre-allocate the buffer to match the size of
	  the attribute being read in, then you should use 
	  SDSfindAttribInfo() or SDSgetAttribInfo() to get the 
	  exact size of the data contained in the attribute before
	  reading.  If the size of the buffer is smaller than the
	  amount of data that needs to be read, the read will be
	  truncated to fit into the buffer.<p>
	  
	  If the stored data is a character string and the returned
	  size==the size of your buffer, you will need to null
	  terminate the buffer yourself after it is returned.
   @enddesc
   @seeroutine SDSwriteAttrib SDSfindAttribInfo SDSgetAttribInfo SDSgetNumAttribs
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     attribname
   @pdesc   Name of the attribute.
   @ptype   char*
   @endpar 

   @par     nelements
   @pdesc   The maximimum number of elements that can be read by the
            from the attribute (ie. the size of the buffer).
   @ptype   int32
   @endpar 

   @par     buffer
   @pdesc   Pointer to the attribute data.
   @ptype   VOIDP
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSreadAttrib(filename,attribname,nelements,buffer)
char *filename;
char *attribname;
int32 nelements;
VOIDP buffer;
{
  int status;
  int32 ref;
  int index;
  int32 nt,count,rtval;
  fileID *id=getFileID(filename);
  if(id->sds_id<0)
   id->sds_id=SDselect(id->id,id->sdsnumber);
  if(id->sds_id<0)
    return -1;
  index=SDfindattr(id->sds_id,attribname);
  /* Check to make sure the number of elements in the attribute
     is less than the number of elements in the buffer */     
  status=SDattrinfo(id->sds_id,index,attribname,&nt,&count);
  if(status<0)
    return status;
  if(count>nelements){
    char *tmpdata=(char *)malloc(count*DFKNTsize(nt));
    rtval=SDreadattr(id->sds_id,index,tmpdata);
    bcopy(tmpdata,buffer,DFKNTsize(nt)*nelements);
    free(tmpdata);
    /* if(nt==DFNT_INT8){
       might be a string, so null terminate 
       buffer[nelements-1]='\0';
       }*/
  }  
  else
    rtval=SDreadattr(id->sds_id,index,buffer);
  return (rtval>=0)?count:rtval;
}

/*@@
   @routine    SDSfindAttribInfo
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Finds the named attribute associated with the current dataset
	  by name and retrieves information about it.  It will return
	  -1 if no attribute with that name is found.
   @enddesc
   @seeroutine SDSreadAttrib SDSgetAttribInfo SDSgetAttribInfo SDSgetNumAttribs
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     attribname
   @pdesc   Name of the attribute.
   @ptype   char*
   @endpar 

   @par     numbertype
   @pdesc   Storage in which to return the numbertype of the attribute
   @pvalues DFNT_FLOAT32 (32-bit float), DFNT_FLOAT64 (64-bit float/double), 
            DFNT_INT8 (char), DFNT_INT16 (short int), DFNT_INT32,(int),
	    DFNT_UINT8 (unsigned char), DFNT_UNINT16 (unsigned short), 
	    DFNT_UINT32 (unsigned int),
	    DFNT_NATIVE (if bit-OR'ed with the above, uses machine native format)
   @ptype   int32*
   @endpar 

   @par     nelements
   @pdesc   Storage in which to return the number of elements in the attribute
   of the type defined by the numbertype param.
            
   @ptype   int32*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSfindAttribInfo(filename,attribname,numbertype,nelements)
char *filename;
char *attribname;
int32 *numbertype;
int32 *nelements;
{
  int32 ref;
  int index;
  fileID *id=getFileID(filename);
  if(id->sds_id<0)
   id->sds_id=SDselect(id->id,id->sdsnumber);
  if(id->sds_id<0)
    return -1;
  index=SDfindattr(id->sds_id,attribname);
  return SDattrinfo(id->sds_id,index,attribname,numbertype,nelements);
}
/*@@
   @routine    SDSgetAttribInfo
   @date       Wed Feb 28 15:37:43 1996
   @desc
           Allows you to get the attribute information by index
	   instead of by name.
   @enddesc
   @seeroutine SDSreadAttrib SDSwriteAttrib SDSfindAttribInfo SDSgetNumAttribs
   @calls      getFileID

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     index
   @pdesc   Index of the attribute attached to the SDS
   @ptype   int32
   @endpar

   @par     attribname
   @pdesc   Storage area for name of the attribute (so that it can be
   returned to the user)
   @ptype   char*
   @endpar 

   @par     numbertype
   @pdesc   Storage in which to return the numbertype of the attribute
   @pvalues DFNT_FLOAT32 (32-bit float), DFNT_FLOAT64 (64-bit float/double), 
            DFNT_INT8 (char), DFNT_INT16 (short int), DFNT_INT32,(int),
	    DFNT_UINT8 (unsigned char), DFNT_UNINT16 (unsigned short), 
	    DFNT_UINT32 (unsigned int),
	    DFNT_NATIVE (if bit-OR'ed with the above, uses machine native format)
   @ptype   int32*
   @endpar 

   @par     nelements
   @pdesc   Storage in which to return the number of elements in the attribute
   of the type defined by the numbertype param.
            
   @ptype   int32*
   @endpar 

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int32 SDSgetAttribInfo(filename,index,attribname,numbertype,nelements)
char *filename;
int32 index;
char *attribname;
int32 *numbertype;
int32 *nelements;
{
  int32 ref;
  fileID *id=getFileID(filename);
  if(id->sds_id<0)
   id->sds_id=SDselect(id->id,id->sdsnumber);
  if(id->sds_id<0)
    return -1;
  return SDattrinfo(id->sds_id,index,attribname,numbertype,nelements);
}
