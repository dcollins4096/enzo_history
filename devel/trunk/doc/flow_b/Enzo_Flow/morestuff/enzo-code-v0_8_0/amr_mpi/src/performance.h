#ifndef PERFORMANCE_H
#define PERFORMANCE_H

//======================================================================
//
// File:        performance.h
//
// Description: Interface layer between Enzo and jbPerf
//
// Requires:    jbPerf if USE_JBPERF is defined
// Requires:     jbMem if USE_JBMEM is defined
// Requires:      PAPI if USE_PAPI is defined
//
//----------------------------------------------------------------------
//
// James Bordner (jbordner@cosmos.ucsd.edu)
// 2003-06-20
//
//======================================================================

#define PERF_EVOLVE /* Compile EvolveLevel calls */
#define PERF_MPI    /* Compile MPI_* jbPerf calls */
#define PERF_COMM   /* Compile Communication* jbPerf calls */

//----------------------------------------------------------------------
/* Zhiling Lan's instrumented part */
//----------------------------------------------------------------------

#ifdef USE_MPI
#define MPI_INSTRUMENTATION
#endif /* USE_MPI */

#ifdef MPI_INSTRUMENTATION

#  define ZLAN_MAX_INDEX 20

//----------------------------------------------------------------------
#  define ZLAN_INIT \
   char filename[ZLAN_MAX_INDEX]; \
   flagging_count=0; \
   moving_count=0; \
   in_count=0; \
   out_count=0; \
   flagging_pct=0.0; \
   moving_pct=0.0; \
   GlobalCommunication = 0.0; \
   RecvComm = 0.0; \
   for (int i=0; i<2*ZLAN_MAX_INDEX; i++) { \
     timer[i] = 0.0; \
     counter[i] = 0; \
   } \
   sprintf(filename,"perfdata_%d.%d",NumberOfProcessors,MyProcessorNumber); \
   filename[strlen(filename)]='\0'; \
   filePtr=fopen(filename,"w");

//----------------------------------------------------------------------
#  define ZLAN_START starttime = MPI_Wtime();

//----------------------------------------------------------------------
#  define ZLAN_STOP(index) \
   endtime = MPI_Wtime(); \
   timer[index]+= endtime-starttime; \
   counter[index] ++; \
   timer[index+ZLAN_MAX_INDEX] += (endtime-starttime)*(endtime-starttime)

//----------------------------------------------------------------------
#  define ZLAN_STOP_GLOBAL(index) \
   ZLAN_STOP(index); \
   GlobalCommunication += ReturnCPUTime() - time1

//----------------------------------------------------------------------
#  define ZLAN_STOP_RECV(index) \
   ZLAN_STOP(index); \
   RecvComm += ReturnCPUTime() - time1

//----------------------------------------------------------------------
#  define ZLAN_COUNT(index,value) \
   timer[index] += double(value); \
   timer[30] += double(value*value);


#else

//----------------------------------------------------------------------
#  define ZLAN_INIT ;
#  define ZLAN_START ;
#  define ZLAN_STOP(index) ;
#  define ZLAN_STOP_GLOBAL(index) ;
#  define ZLAN_STOP_RECV(index) ;
#  define ZLAN_COUNT(index,value) ;

#endif


#ifdef USE_JBMEM
#  include "jbMem.h"

#  define JBMEM_MESSAGE(ip) \
      { \
	long long bcg,bmg,bcl,bml; \
	bcg = jb::mem.bytes(); \
	bmg = jb::mem.bytes_max(); \
	bcl = jb::mem.bytes_local(); \
	bml = jb::mem.bytes_max_local(); \
	if (ip==0) printf ("jbMem: %lld %lld %lld %lld\n", \
			   bcg,bmg,bcl,bml); \
      }


/* Definitions for using jbMem NEW and DELETE functions */

// #   define new NEW
// #   define delete DELETE
// extern void * NEW (size_t bytes) throw (std::bad_alloc);
// extern void DELETE (void *p) throw ();
#else
#  define JBMEM_MESSAGE(ip) ;
#endif

#ifdef USE_JBPERF

#   include "jbPerf.h"

#   define JBPERF_INIT \
     const int jbLength = 20; \
     char jbRegion[jbLength];

#   define JBPERF_START(REGION) \
     jb::perf.start(REGION);

#   define JBPERF_STOP(REGION) \
     jb::perf.stop(REGION);

#   define JBPERF_STOP_BYTES(REGION,COUNT,TYPE) \
     { \
       int bytes; \
       MPI_Type_size(TYPE,&bytes); \
       jb::perf.stop(REGION,COUNT*bytes); \
     }

#   define JBPERF_WRITE \
       jb::perf.write ();

#   define JBPERF_ADVANCE \
       jb::perf.advance ();

#   ifdef PERF_COMM
#      define JBPERF_COMM_START(REGION) \
        jb::perf.start(REGION);

#      define JBPERF_COMM_STOP(REGION) \
        jb::perf.stop(REGION);
#   else
#      define JBPERF_COMM_START(REGION) ;
#      define JBPERF_COMM_STOP(REGION) ;
#   endif

#   ifdef PERF_MPI
#      define JBPERF_MPI_START(REGION) \
        jb::perf.start(REGION);
#      define JBPERF_MPI_STOP(REGION) \
        jb::perf.stop(REGION);
#      define JBPERF_MPI_STOP_BYTES(REGION,COUNT,SIZE) \
        { \
          int bytes; \
          MPI_Type_size(TYPE,&bytes); \
          jb::perf.stop(REGION,COUNT*bytes); \
        }
#   else
#      define JBPERF_MPI_START(REGION) ;
#      define JBPERF_MPI_STOP(REGION) ;
#      define JBPERF_MPI_STOP_BYTES(REGION) ;
#   endif

#   ifdef PERF_EVOLVE
#      define JBPERF_EVOLVE_START(REGION) \
        jb::perf.start(REGION);

#      define JBPERF_EVOLVE_STOP(REGION) \
        jb::perf.stop(REGION);
#   else
#      define JBPERF_EVOLVE_START(REGION) ;
#      define JBPERF_EVOLVE_STOP(REGION) ;
#   endif

#   define JBPERF_NEXT(REGION) \
     jb::perf.next(REGION);

#   define JBPERF_START_LEVEL(REGION) \
     if (snprintf (jbRegion,jbLength,REGION".%d",level) >= jbLength) { \
       WARNING_MESSAGE; \
     } \
     jb::perf.start(jbRegion);

#   define JBPERF_STOP_LEVEL(REGION) \
     if (snprintf (jbRegion,jbLength,REGION".%d",level) >= jbLength) { \
       WARNING_MESSAGE; \
     } \
     jb::perf.stop(jbRegion);
#   define JBPERF_NEXT_LEVEL(REGION) \
     if (snprintf (jbRegion,jbLength,REGION".%d",level) >= jbLength) { \
       WARNING_MESSAGE; \
     } \
     jb::perf.next(jbRegion);

#else

#   define JBPERF_INIT ;
#   define JBPERF_START(X) ;
#   define JBPERF_STOP(X) ;
#   define JBPERF_STOP_BYTES(X,COUNT,TYPE) ;
#   define JBPERF_WRITE ;
#   define JBPERF_ADVANCE ;
#   define JBPERF_COMM_START(X) ;
#   define JBPERF_COMM_STOP(X) ;
#   define JBPERF_MPI_START(X) ;
#   define JBPERF_MPI_STOP(X) ;
#   define JBPERF_EVOLVE_START(X) ;
#   define JBPERF_EVOLVE_STOP(X) ;
#   define JBPERF_NEXT(X) ;
#   define JBPERF_START_LEVEL(X) ;
#   define JBPERF_STOP_LEVEL(X) ;
#   define JBPERF_NEXT_LEVEL(X) ;

#endif

#endif
