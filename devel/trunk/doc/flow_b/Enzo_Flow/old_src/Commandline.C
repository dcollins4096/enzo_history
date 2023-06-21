#ifdef USE_FLEXIO
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Commandline.h"
void Commandline::JadOpt::setDefaults(){
  enabled=0;
  dumplevel=-1; // -1 means dump to deepest level
  startstep=endstep=-1;
  starttime=endtime=-1;
  filetype = Default;
  strcpy(fileext,"raw");
  strcpy(fileprefix,"enzodata");
  filename = new char[sizeof("enzodata.amr.raw")+1];
  strcpy(filename,"enzodata.amr.raw");
}
void Commandline::JadOpt::print(FILE *fp){
  fprintf(fp,"Jad Output Options------------\n");
  if(enabled){
    fprintf(fp,"\tEnabled\n");
    fprintf(fp,"\tFilename prefix=[%s] ext=[%s]\n",fileprefix,fileext);
    fprintf(fp,"\t  for a filename of %s.amr.%s\n",fileprefix,fileext);
    if(dumplevel>=0) fprintf(fp,"\tdumplevel=%u\n",dumplevel);
    else fprintf(fp,"\tdumplevel=MaxDepthOfHierarchy\n");
    if(startstep>=0) fprintf(fp,"\tstartstep=%u : ",startstep);
    else fprintf(fp,"\tstartstep=Undefined : ",startstep);
    if(endstep>=0) fprintf(fp,"\tendstep=%u\n",endstep);
    else fprintf(fp,"\tendstep=Undefined\n",endstep);
    if(starttime>=0) fprintf(fp,"\tstarttime=%f : ",starttime);
    else fprintf(fp,"\tstarttime=Undefined : ",starttime);
    if(endtime>=0) fprintf(fp,"\tendtime=%f\n",endtime);
    else fprintf(fp,"\tendtime=Undefined\n",endtime);
  }
  else {
    printf("\tDisabled\n");
  }
}
int Commandline::JadOpt::parseOption(int argc,char *argv[],int &i){ 
      char *str=argv[i];
      if(*str=='-'){
	//puts("here1");
	str++;
	if(!strncmp(str,"jad",3)){ // verify this is JAD option
	  //puts("here2");
	  // create a copy of this arg
	  char *argp,*arg = new char[strlen(str)+1];
	  argp=arg;
	  char *nextarg;
	  strcpy(arg,str);
	  // lets get our suboptions
	  enabled=1; // enable jad first
	  // chop up using strtok on ":"
	  if(arg=strchr(arg,':')){
	    arg++;
	    // while((nextarg=strchr(arg,':')))
	    do{
	      char *val;
	      //printf("Arg=[%s]\n",arg);
	      nextarg=strchr(arg,':');
	      if(nextarg){
		*nextarg='\0';
		nextarg++;
	      }
	      val = strchr(arg,'=');
	      //puts("did strchr of arg to make val");
	      *val='\0';
	      val++;
	      //puts("after");
	      if(!strncmp(arg,"filep",5)){ // fileprefix
		// filename.amr.filetype
		strcpy(fileprefix,val);
		// now construct filename
		if(filename) delete filename;
		filename=new char[strlen(fileprefix)+5+strlen(fileext)+2];
		sprintf(filename,"%s.amr.%s",fileprefix,fileext);
	      }
	      else if(!strncmp(arg,"filet",5)){ // filetype
		// Recognizes hdf,raw,hdf4,hdf5,ieee,raw
		// save as file-extension as well as 
		// an enumerated filetype
		if(!strcmp(val,"HDF") || 
		   !strcmp(val,"hdf") ||
		   !strcmp(val,"hdf4") ||
		   !strcmp(val,"HDF4")){
		  filetype=HDF4;
		  strcpy(fileext,val);
		}
		else if(!strcmp(val,"HDF5") ||
			!strcmp(val,"hdf5")){
		  fprintf(stderr,"JadParser: Suboption filetype=%s.  %s is unsupported right now\n");
		  /*
		    filetype=HDF5;
		    strcpy(fileext,val);
		  */
		}
		else if(!strcmp(val,"IEEE") ||
			!strcmp(val,"ieee") ||
			!strcmp(val,"raw")){
		  filetype=IEEE;
		  strcpy(fileext,val);
		}
		else {
		  fprintf(stderr,"JadParser: Suboption filetype=%s failed\n",
			  val);
		  fprintf(stderr,"\tAvailable options are hdf,hdf4,hdf5,raw,ieee\n");
		}
		// now construct filename
		if(filename) delete filename;
		filename=new char[strlen(fileprefix)+5+strlen(fileext)+2];
		sprintf(filename,"%s.amr.%s",fileprefix,fileext);
	      }
	      else if(!strcmp(arg,"field") || !strcmp(arg,"fieldname")){
		// expect , separated list of fields
		char *s=strtok(val,",");
		while(s){
		  char *ns=new char[strlen(s)+1];
		  fieldnames.append(ns);
		  s=strtok(0,","); // next token
		}
	      }
	      else if(!strcmp(arg,"dumplevel") || !strcmp(arg,"level")){ 
		// what level to dump
		dumplevel=atoi(val);
	      }
	      else if(!strcmp(arg,"starttime")){
		starttime=atof(val);
	      }
	      else if(!strcmp(arg,"endtime")){
		starttime=atof(val);
	      }
	      else if(!strncmp(arg,"startiter",9)||!strcmp(arg,"startstep")){
		startstep=atoi(val);
	      }
	      else if(!strncmp(arg,"enditer",7)||!strcmp(arg,"endstep")){
		endstep=atoi(val);
	      }
	      else {
		fprintf(stderr,"JAD cmdln parser: Unrecognized token \"%s\"\n",arg);
	      }
	      arg = nextarg;
	    } while(arg);
	  }
	  delete argp;
	  return 1;
	}
	else return 0;
      }
      return 0;
 }
void Commandline::StarOpt::setDefaults(){
  enabled=0;
  dumplevel=-1; // -1 means dump to deepest level
  startstep=endstep=-1;
  starttime=endtime=-1;
  filetype = Default;
  strcpy(fileext,"raw");
  strcpy(fileprefix,"enzodata");
  filename = new char[sizeof("enzodata.stars.raw")+1];
  strcpy(filename,"enzodata.stars.raw");
}
void Commandline::StarOpt::print(FILE *fp){
  fprintf(fp,"Star Output Options------------\n");
  if(enabled){
    fprintf(fp,"\tEnabled\n");
    fprintf(fp,"\tFilename prefix=[%s] ext=[%s]\n",fileprefix,fileext);
    fprintf(fp,"\t  for a filename of %s.star.%s\n",fileprefix,fileext);
    if(dumplevel>=0) fprintf(fp,"\tdumplevel=%u\n",dumplevel);
    else fprintf(fp,"\tdumplevel=MaxDepthOfHierarchy\n");
    if(startstep>=0) fprintf(fp,"\tstartstep=%u : ",startstep);
    else fprintf(fp,"\tstartstep=Undefined : ",startstep);
    if(endstep>=0) fprintf(fp,"\tendstep=%u\n",endstep);
    else fprintf(fp,"\tendstep=Undefined\n",endstep);
    if(starttime>=0) fprintf(fp,"\tstarttime=%f : ",starttime);
    else fprintf(fp,"\tstarttime=Undefined : ",starttime);
    if(endtime>=0) fprintf(fp,"\tendtime=%f\n",endtime);
    else fprintf(fp,"\tendtime=Undefined\n",endtime);
  }
  else {
    printf("\tDisabled\n");
  }
}

int Commandline::StarOpt::parseOption(int argc,char *argv[],int &i){
      char *str=argv[i];
      if(*str=='-'){
	char *val=argv[i+1];
	str++;
	if(!strncmp(str,"star",3)){ // verify this is a star option
	  // create a copy of this arg
	  char *argp,*arg = new char[strlen(str)+1];
	  argp=arg;
	  char *nextarg;
	  strcpy(arg,str);
	  // lets get our suboptions
	  enabled=1; // enable jad first
	  // chop up using strtok on ":"
	  if(arg=strchr(arg,':')){
	    arg++;
	    // while((nextarg=strchr(arg,':'))){
	    do{
	      nextarg=strchr(arg,':');
	      if(nextarg){
		*nextarg='\0';
		nextarg++;
	      }
	      val = strchr(arg,'=');
	      *val='\0';
	      val++;
	      if(!strncmp(arg,"filep",5)){ // fileprefix
		// filename.amr.filetype
		strcpy(fileprefix,val);
		// now construct filename
		if(filename) delete filename;
		filename=new char[strlen(fileprefix)+8+strlen(fileext)+2];
		sprintf(filename,"%s.stars.%s",fileprefix,fileext);
	      }
	      else if(!strncmp(arg,"filet",5)){ // filetype
		// Recognizes hdf,raw,hdf4,hdf5,ieee,raw
		// save as file-extension as well as 
		// an enumerated filetype
		if(!strcmp(val,"HDF") || 
		   !strcmp(val,"hdf") ||
		   !strcmp(val,"hdf4") ||
		   !strcmp(val,"HDF4")){
		  filetype=HDF4;
		  strcpy(fileext,val);
		}
		else if(!strcmp(val,"HDF5") ||
			!strcmp(val,"hdf5")){
		  fprintf(stderr,"StarParser: Suboption filetype=%s.  %s is unsupported right now\n");
		  /*
		    filetype=HDF5;
		    strcpy(fileext,val);
		  */
		}
		else if(!strcmp(val,"IEEE") ||
			!strcmp(val,"ieee") ||
			!strcmp(val,"raw")){
		  filetype=IEEE;
		  strcpy(fileext,val);
		}
		else {
		  fprintf(stderr,"StarParser: Suboption filetype=%s failed\n",
			  val);
		  fprintf(stderr,"\tAvailable options are hdf,hdf4,hdf5,raw,ieee\n");
		}
		// now construct filename
		if(filename) delete filename;
		filename=new char[strlen(fileprefix)+8+strlen(fileext)+2];
		sprintf(filename,"%s.stars.%s",fileprefix,fileext);
	      }
	      else if(!strcmp(arg,"dumplevel") || !strcmp(arg,"level")){ 
		// what level to dump
		dumplevel=atoi(val);
	      }
	      else if(!strcmp(arg,"starttime")){
		starttime=atof(val);
	      }
	      else if(!strcmp(arg,"endtime")){
		starttime=atof(val);
	      }
	      else if(!strncmp(arg,"startiter",9)||!strcmp(arg,"startstep")){
		startstep=atoi(val);
	      }
	      else if(!strncmp(arg,"enditer",7)||!strcmp(arg,"endstep")){
		endstep=atoi(val);
	      }
	      else {
		fprintf(stderr,"Star cmdln parser: Unrecognized token %s\n",arg);
	      }
	      arg = nextarg;
	    }while(arg);
	  }
	  delete argp;
	  return 1;
	}
	else return 0;
      }
      return 0;
}
void Commandline::DMOpt::setDefaults(){
  enabled=0;
  dumplevel=-1; // -1 means dump to deepest level
  startstep=endstep=-1;
  starttime=endtime=-1;
  filetype = Default;
  strcpy(fileext,"raw");
  strcpy(fileprefix,"enzodata");
  filename = new char[sizeof("enzodata.stars.raw")+1];
  strcpy(filename,"enzodata.dm.raw");
}
void Commandline::DMOpt::print(FILE *fp){
  fprintf(fp,"DarkMatter Particle Output Options------------\n");
  if(enabled){
    fprintf(fp,"\tEnabled\n");
    fprintf(fp,"\tFilename prefix=[%s] ext=[%s]\n",fileprefix,fileext);
    fprintf(fp,"\t  for a filename of %s.dm.%s\n",fileprefix,fileext);
    if(dumplevel>=0) fprintf(fp,"\tdumplevel=%u\n",dumplevel);
    else fprintf(fp,"\tdumplevel=MaxDepthOfHierarchy\n");
    if(startstep>=0) fprintf(fp,"\tstartstep=%u : ",startstep);
    else fprintf(fp,"\tstartstep=Undefined : ",startstep);
    if(endstep>=0) fprintf(fp,"\tendstep=%u\n",endstep);
    else fprintf(fp,"\tendstep=Undefined\n",endstep);
    if(starttime>=0) fprintf(fp,"\tstarttime=%f : ",starttime);
    else fprintf(fp,"\tstarttime=Undefined : ",starttime);
    if(endtime>=0) fprintf(fp,"\tendtime=%f\n",endtime);
    else fprintf(fp,"\tendtime=Undefined\n",endtime);
  }
  else {
    printf("\tDisabled\n");
  }
}

int Commandline::DMOpt::parseOption(int argc,char *argv[],int &i){
      char *str=argv[i];
      if(*str=='-'){
	char *val=argv[i+1];
	str++;
	if(!strncmp(str,"dm",3)){ // must check the debug option for
	  // create a copy of this arg
	  char *argp;
	  char *arg = new char[strlen(str)+1];
	  argp=arg; // placeholder for arg pointer
	  char *nextarg;
	  strcpy(arg,str);
	  // lets get our suboptions
	  enabled=1; // enable jad first
	  // chop up using strtok on ":"
	  if(arg=strchr(arg,':')){
	    arg++;
	    // while((nextarg=strchr(arg,':'))){
	    do{
	      nextarg=strchr(arg,':');
	      if(nextarg){
		*nextarg='\0';
		nextarg++;
	      }
	      val = strchr(arg,'=');
	      *val='\0';
	      val++;
	      if(!strncmp(arg,"filep",5)){ // fileprefix
		// filename.amr.filetype
		strcpy(fileprefix,val);
		// now construct filename
		if(filename) delete filename;
		filename=new char[strlen(fileprefix)+5+strlen(fileext)+2];
		sprintf(filename,"%s.dm.%s",fileprefix,fileext);
	      }
	      else if(!strncmp(arg,"filet",5)){ // filetype
		// Recognizes hdf,raw,hdf4,hdf5,ieee,raw
		// save as file-extension as well as 
		// an enumerated filetype
		if(!strcmp(val,"HDF") || 
		   !strcmp(val,"hdf") ||
		   !strcmp(val,"hdf4") ||
		   !strcmp(val,"HDF4")){
		  filetype=HDF4;
		  strcpy(fileext,val);
		}
		else if(!strcmp(val,"HDF5") ||
			!strcmp(val,"hdf5")){
		  fprintf(stderr,"DMParser: Suboption filetype=%s.  %s is unsupported right now\n");
		  /*
		    filetype=HDF5;
		    strcpy(fileext,val);
		  */
		}
		else if(!strcmp(val,"IEEE") ||
			!strcmp(val,"ieee") ||
			!strcmp(val,"raw")){
		  filetype=IEEE;
		  strcpy(fileext,val);
		}
		else {
		  fprintf(stderr,"DMParser: Suboption filetype=%s failed\n",
			  val);
		  fprintf(stderr,"\tAvailable options are hdf,hdf4,hdf5,raw,ieee\n");
		}
		// now construct filename
		if(filename) delete filename;
		filename=new char[strlen(fileprefix)+5+strlen(fileext)+2];
		sprintf(filename,"%s.dm.%s",fileprefix,fileext);
	      }
	      else if(!strcmp(arg,"dumplevel") || !strcmp(arg,"level")){ 
		// what level to dump
		dumplevel=atoi(val);
	      }
	      else if(!strcmp(arg,"starttime")){
		starttime=atof(val);
	      }
	      else if(!strcmp(arg,"endtime")){
		starttime=atof(val);
	      }
	      else if(!strncmp(arg,"startiter",9)||!strcmp(arg,"startstep")){
		startstep=atoi(val);
	      }
	      else if(!strncmp(arg,"enditer",7)||!strcmp(arg,"endstep")){
		endstep=atoi(val);
	      }
	      else {
		fprintf(stderr,"DM cmdln parser: Unrecognized token %s\n",arg);
	      }
	      arg = nextarg;
	    } while(arg);
	  }
	  delete argp;
	  return 1;
	}
	else return 0;
      }
      return 0;
}

int Commandline::parseOption(int argc,char *argv[],int &i){
  char *str=argv[i];
  if(*str=='-'){
    char *val=argv[i+1];
    str++;
    if(stars.parseOption(argc,argv,i)){ // delegate to star subparser
      return 1;
    }
    else if(jad.parseOption(argc,argv,i)){ // delegate to jad subparser
      return 1;
    }
    else return 0;
  }
  return 0; // no option parsed
#if 0
  else{
    // its the infilename which should be a .hierarchy file
    filename = str; // just point directly (why copy?)
    // got to get more sophisticated here
    if(strchr(filename,'/')) {
      path = new char[strlen(filename)+1];
      strcpy(path,filename);
      char *s = strrchr(path,'/');
      s[1] = '\0';
      printf("path=[%s] filename=[%s]\n",path,filename);
    }
    else {
      path = new char[1];
      *path = '\0';
      printf("path=[%s] filename=[%s]\n",path,filename);
    }
  }
#endif
}
void Commandline::parse(int argc,char *argv[]){
  programname=argv[0];
  // printf("parse file\n");
  for(int i=1;i<argc;i++){
    parseOption(argc,argv,i);
  }
}
void Commandline::parseAndEat(int &argc,char *argv[]){
  programname=argv[0];
  int nargc=argc;
  // printf("parse file\n");
  for(int o=1,i=1;i<argc;i++){
    if(parseOption(argc,argv,i)) // eat that option
      nargc--;
    else 
      o++;
    // copy down
    if((i+1) < argc)
      argv[o] = argv[i+1];
  }
  argc=nargc;
}
#endif /* USE_FLEXIO */
