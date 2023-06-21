#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/times.h>
#include <limits.h>
#include <time.h>
#include "Arch.h"
#include "IO.hh"
#include "IEEEIO.hh"
#include "HDFIO.hh"
#include "FlexArrayTmpl.H"

#define f_openf77 F77NAME(openf77_,openf77,OPENF77)
#define f_writef77 F77NAME(writef77_,writef77,WRITEF77)
#define f_closef77 F77NAME(closef77_,closef77,CLOSEF77)
extern "C" {
void f_openf77();
void f_writef77(double *data);
void f_closef77();
}

//typedef FlexArray<float> FloatVector;
//
//class fVectorStatistics : public FloatVector {
//  fVectorStatistics(int sz=0):FloatVector(sz){}
//};

struct ResultRecord {
  float usertime,systemtime,realtime,combined,nbytes,megs_per_second;
  // ResultRecord():usertime(0),systemtime(0),realtime(0),
  //  combined(0),nbytes(0),megs_per_second(0){}
};

#define ClearLine fprintf(stderr,"\r                                                                                         \r")

struct ResultRecordList {
  FlexArray<ResultRecord> results;
  ResultRecord average_results;
  ResultRecordList(){}
  void append(ResultRecord &rec);
  ResultRecord &operator[](int index){
    return results[index];
  }
  
};

void ResultRecordList::append(ResultRecord &rec){
  int idx=results.getSize();
  float nrecs=(float)idx;
  results.append(rec);
  average_results.usertime=(average_results.usertime*nrecs+rec.usertime)/(nrecs+1);
  average_results.systemtime=(average_results.systemtime*nrecs+rec.systemtime)/(nrecs+1);
  average_results.realtime=(average_results.realtime*nrecs+rec.realtime)/(nrecs+1);
  average_results.combined=(average_results.combined*nrecs+rec.combined)/(nrecs+1);
  average_results.nbytes=(average_results.nbytes*nrecs+rec.nbytes)/(nrecs+1);
  rec.megs_per_second = (double)(rec.nbytes)/(1024.0L*1024.0L*(double)rec.realtime);
  average_results.megs_per_second = (average_results.megs_per_second*nrecs+rec.megs_per_second)/(nrecs+1);
}

void main(int argc,char *argv[]){
  double data[64*64*64];
  int dims[3]={64,64,64};
  int i,nds;
  int ntests;
  struct tms stms,etms;
  struct timeval stmv,etmv;
  double srt,ert;
  ResultRecordList ieeeio_results,f77_results,hdf_results;

  if(argc>1){
    ntests=atoi(argv[1]);
    printf("Ntests=%u\n",ntests);
  }
  else {
    ntests=1;
  }

  for(int tst=0;tst<ntests;tst++){
    ResultRecord results;
    
    printf("\nTrial %u --------------\n",tst); // next line
    puts("--------------IEEE--------------");
    for(i=0,nds=5;i<8;i++,nds+=5){
      fprintf(stderr,"IEEEIO <open speed.raw> ");
      times(&stms);
      gettimeofday(&stmv);
      IObase *file = new IEEEIO("speed.raw",IObase::Write);
      fprintf(stderr,"Write %2u datasets:",nds);
      for(int n=0;n<nds;n++){
	fprintf(stderr,"*");
	file->write(IObase::Float64,3,dims,data);
      }
      delete file;
      gettimeofday(&etmv);
      times(&etms); // times after close to account for buffer flushing
      puts("");
      ert=(double)etmv.tv_sec + (double)etmv.tv_usec/1000000.0L;
      srt=(double)stmv.tv_sec + (double)stmv.tv_usec/1000000.0L;
      results.realtime = ert-srt;
      results.usertime = (float)(etms.tms_utime-stms.tms_utime)/CLK_TCK;
      results.systemtime = (float)(etms.tms_stime-stms.tms_stime)/CLK_TCK;
      results.combined = (float)(etms.tms_utime-stms.tms_utime+
				 etms.tms_stime-stms.tms_stime)/CLK_TCK;
      results.nbytes = (float)(nds*IObase::nBytes(IObase::Float64,3,dims));
      ieeeio_results.append(results);
    }
    puts("--------------HDF---------------");
    for(i=0,nds=5;i<8;i++,nds+=5){
      fprintf(stderr,"HDF <open speed.hdf> ");
      times(&stms);
      gettimeofday(&stmv);
      IObase *file = new HDFIO("speed.hdf",IObase::Write);
      fprintf(stderr,"Write %2u datasets:",nds);
      for(int n=0;n<nds;n++){
	fprintf(stderr,"*");
	file->write(IObase::Float64,3,dims,data);
      }
      delete file;
      gettimeofday(&etmv);
      times(&etms); // times after close to account for buffer flushing
      puts("");
      ert=(double)etmv.tv_sec + (double)etmv.tv_usec/1000000.0L;
      srt=(double)stmv.tv_sec + (double)stmv.tv_usec/1000000.0L;
      results.realtime = ert-srt;
      results.usertime = (float)(etms.tms_utime-stms.tms_utime)/CLK_TCK;
      results.systemtime = (float)(etms.tms_stime-stms.tms_stime)/CLK_TCK;
      results.combined = (float)(etms.tms_utime-stms.tms_utime+
				      etms.tms_stime-stms.tms_stime)/CLK_TCK;
      results.nbytes = (float)(nds*IObase::nBytes(IObase::Float64,3,dims)); 
      hdf_results.append(results);
    }  
    puts("-------------F77 Unformatted---------------");
    for(i=0,nds=5;i<8;i++,nds+=5){
      fprintf(stderr,"F77 Unf <open f77speed.unf> ");
      times(&stms);
      gettimeofday(&stmv);
      f_openf77();
      fprintf(stderr,"Write %2u datasets:",nds);
      for(int n=0;n<nds;n++){
	fprintf(stderr,"*");
	f_writef77(data);
      }
      f_closef77();
      gettimeofday(&etmv);
      times(&etms); // times after close to account for buffer flushing
      puts("");
      ert=(double)etmv.tv_sec + (double)etmv.tv_usec/1000000.0L;
      srt=(double)stmv.tv_sec + (double)stmv.tv_usec/1000000.0L;
      results.realtime = ert-srt;
      results.usertime = (float)(etms.tms_utime-stms.tms_utime)/CLK_TCK;
      results.systemtime = (float)(etms.tms_stime-stms.tms_stime)/CLK_TCK;
      results.combined = (float)(etms.tms_utime-stms.tms_utime+
				      etms.tms_stime-stms.tms_stime)/CLK_TCK;
      results.nbytes = (float)(nds*IObase::nBytes(IObase::Float64,3,dims)); 
      f77_results.append(results);
    }
  }

  // Now for all of the results
  printf("-------IEEE: Average of %u trials------------\n",ntests);
  printf("\tRealtime=%f, UserTime=%f, SystemTime=%f, Combined=%f Megs/sec=%f\n",
	 ieeeio_results.average_results.realtime,
	 ieeeio_results.average_results.usertime,
	 ieeeio_results.average_results.systemtime,
	 ieeeio_results.average_results.combined,
	 ieeeio_results.average_results.megs_per_second);
  printf("-------HDF: Average of %u trials------------\n",ntests);
  printf("\tRealtime=%f, UserTime=%f, SystemTime=%f, Combined=%f Megs/sec=%f\n",
	 hdf_results.average_results.realtime,
	 hdf_results.average_results.usertime,
	 hdf_results.average_results.systemtime,
	 hdf_results.average_results.combined,
	 hdf_results.average_results.megs_per_second);
  printf("-------F77 UNF: Average of %u trials------------\n",ntests);
  printf("\tRealtime=%f, UserTime=%f, SystemTime=%f, Combined=%f Megs/sec=%f\n",
	 f77_results.average_results.realtime,
	 f77_results.average_results.usertime,
	 f77_results.average_results.systemtime,
	 f77_results.average_results.combined,
	 f77_results.average_results.megs_per_second);
}


