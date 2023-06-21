/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
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

// Local defines

#undef USE_JBPERF_MPI
#undef USE_JBPERF_LEVEL

//----------------------------------------------------------------------
// Zhiling Lan's instrumented part
//----------------------------------------------------------------------

#ifdef USE_MPI
#define MPI_INSTRUMENTATION
#endif /* USE_MPI */

#ifndef MPI_INSTRUMENTATION

#  define ZLAN_INIT ;
#  define ZLAN_START ;
#  define ZLAN_STOP(index) ;
#  define ZLAN_STOP_GLOBAL(index) ;
#  define ZLAN_STOP_RECV(index) ;
#  define ZLAN_COUNT(index,value) ;

# else

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
#endif

//----------------------------------------------------------------------
// jbMem
//----------------------------------------------------------------------
#define JBMEM_LOCAL
#ifdef JBMEM_LOCAL

#ifndef USE_JBMEM
#  define JBMEM_MESSAGE(ip,MESSAGE) ;
#else
#  include "jbMem.h"
#  define JBMEM_MESSAGE(ip,MESSAGE) \
      { \
	long long bcg,bmg,bcl,bml; \
	bcl = jb::mem.bytes_local_max(); \
	bml = jb::mem.bytes_local_max_high(); \
     fprintf (stderr,"\n%30s: proc %4d mem-curr    %lld %lld\n", MESSAGE, ip, bcl, bcl - MHDPreviousCurrent);\
     fprintf (stderr,"\n%30s: proc %4d mem-high    %lld %lld\n", MESSAGE, ip, bml, bml - MHDPreviousGlobal); \
           MHDPreviousCurrent=bcl;\
           MHDPreviousGlobal=bml;\
           fflush(stdout); \
      }
#endif

#else //JBMEM_LOCAL


#ifndef USE_JBMEM
#  define JBMEM_MESSAGE(ip,MESSAGE) ;
#else
#  include "jbMem.h"
#  define JBMEM_MESSAGE(ip,MESSAGE) \
      { \
	long long bcg,bmg,bcl,bml; \
	bcg = jb::mem.bytes_global(); \
	bmg = jb::mem.bytes_global_high(); \
	bcl = jb::mem.bytes_local_max(); \
	bml = jb::mem.bytes_local_max_high(); \
        if (ip==0) { \
     printf ("%30s: mem-global-curr    %lld %lld\n", MESSAGE, bcg, bcg - MHDPreviousCurrent);\
     printf ("%30s: mem-global-high    %lld %lld\n", MESSAGE, bmg, bmg - MHDPreviousGlobal); \
     printf ("%30s: mem-local-max-curr %lld\n", MESSAGE, bcl); \
     printf ("%30s: mem-local-max-high %lld\n", MESSAGE, bml); \
           MHDPreviousCurrent=bcg;\
           MHDPreviousGlobal=bmg;\
           fflush(stdout); \
        } \
      }
#endif

#endif //JBMEM_Local
//----------------------------------------------------------------------
// Redefine new[] and delete[] as NEW() and DELETE()
//----------------------------------------------------------------------

#define NEW(VAR,TYPE,SIZE) VAR = new TYPE [SIZE];
#define DELETE(VAR) delete [] VAR;

//----------------------------------------------------------------------
// jbPerf
//----------------------------------------------------------------------

#ifndef USE_JBPERF

#   define JBPERF_INIT ;
#   define JBPERF_START(REGION) ;
#   define JBPERF_STOP(REGION) ;
#   define JBPERF_WRITE ;
#   define JBPERF_ADVANCE ;
#   define JBPERF_NEXT(REGION) ;

# else

#   include "jbPerf.h"
//----------------------------------------------------------------------
#   define JBPERF_INIT \
     const int jbLength = 20; \
     char jbRegion[jbLength];
//----------------------------------------------------------------------
#   define JBPERF_START(REGION) \
     jb::perf.start(REGION);
//----------------------------------------------------------------------
#   define JBPERF_STOP(REGION) \
     jb::perf.stop(REGION);
//----------------------------------------------------------------------
#   define JBPERF_WRITE \
       jb::perf.write ();
//----------------------------------------------------------------------
#   define JBPERF_ADVANCE \
       jb::perf.advance ();
//----------------------------------------------------------------------
#   define JBPERF_NEXT(REGION) \
     jb::perf.next(REGION);
//----------------------------------------------------------------------
#endif

//----------------------------------------------------------------------
// jbPerf with levels
//----------------------------------------------------------------------

#if !defined(USE_JBPERF)
#   define JBPERF_START_LEVEL(REGION) ;
#   define JBPERF_STOP_LEVEL(REGION) ;
#   define JBPERF_NEXT_LEVEL(REGION) ;
#else
#   if !defined (USE_JBPERF_LEVEL)
#      define JBPERF_START_LEVEL(REGION) JBPERF_START(REGION) ;
#      define JBPERF_STOP_LEVEL(REGION) JBPERF_STOP(REGION) ;
#      define JBPERF_NEXT_LEVEL(REGION) JBPERF_NEXT(REGION) ;
#   else
#      define JBPERF_START_LEVEL(REGION) \
        if (snprintf (jbRegion,jbLength,REGION".%d",level) >= jbLength) { \
          WARNING_MESSAGE; \
        } \
        jb::perf.start(jbRegion);
//----------------------------------------------------------------------
#      define JBPERF_STOP_LEVEL(REGION) \
        if (snprintf (jbRegion,jbLength,REGION".%d",level) >= jbLength) { \
          WARNING_MESSAGE; \
        } \
        jb::perf.stop(jbRegion);
//----------------------------------------------------------------------
#      define JBPERF_NEXT_LEVEL(REGION) \
        if (snprintf (jbRegion,jbLength,REGION".%d",level) >= jbLength) { \
          WARNING_MESSAGE; \
        } \
        jb::perf.next(jbRegion);
#   endif
#endif

//----------------------------------------------------------------------
// jbPerf && MPI
//----------------------------------------------------------------------

#if !defined(USE_JBPERF) || !defined(USE_MPI)
#   define JBPERF_START_MPI_SEND(REGION,COUNT,TYPE) ;
#   define JBPERF_STOP_MPI_SEND(REGION,COUNT,TYPE) ;
#   define JBPERF_START_MPI_RECV(REGION,COUNT,TYPE) ;
#   define JBPERF_STOP_MPI_RECV(REGION,COUNT,TYPE) ;
#   define JBPERF_START_MPI_SENDRECV(REGION,SCOUNT,STYPE,RCOUNT,RTYPE) ;
#   define JBPERF_STOP_MPI_SENDRECV(REGION,SCOUNT,STYPE,RCOUNT,RTYPE) ;
#   define JBPERF_START_MPI_GATHER(REGION,COUNT,TYPE) ;
#   define JBPERF_STOP_MPI_GATHER(REGION,COUNT,TYPE) ;
#   define JBPERF_START_MPI_REDUCE(REGION,COUNT,TYPE) ;
#   define JBPERF_STOP_MPI_REDUCE(REGION,COUNT,TYPE) ;
#   define JBPERF_START_MPI_BARRIER(REGION) ;
#   define JBPERF_STOP_MPI_BARRIER(REGION) ;
#   define JBPERF_COUNT_MPI(REGION,COUNT,TYPE) ;
#else
   namespace jb {
     extern int mpi_send_start,mpi_send_stop;
     extern int mpi_recv_start,mpi_recv_stop;
     extern int mpi_sendrecv_start,mpi_sendrecv_stop;
     extern int mpi_barrier_start,mpi_barrier_stop;
     extern int mpi_gather_start,mpi_gather_stop;
     extern int mpi_reduce_start,mpi_reduce_stop;
   };
//----------------------------------------------------------------------
#   define JBPERF_START_MPI_SEND(REGION,COUNT,TYPE) \
     JBPERF_MPE(jb::mpi_send_start,0,"send-start"); \
     JBPERF_MPI_START(REGION);
//----------------------------------------------------------------------
#   define JBPERF_STOP_MPI_SEND(REGION,COUNT,TYPE) \
     JBPERF_COUNT_MPI("mpi-send",COUNT,TYPE) \
     JBPERF_MPE(jb::mpi_send_stop,0,"send-stop"); \
     JBPERF_MPI_STOP(REGION);
//----------------------------------------------------------------------
#   define JBPERF_START_MPI_RECV(REGION,COUNT,TYPE) \
     JBPERF_MPE(jb::mpi_recv_start,0,"recv-start"); \
     JBPERF_MPI_START(REGION);
//----------------------------------------------------------------------
#   define JBPERF_STOP_MPI_RECV(REGION,COUNT,TYPE) \
     JBPERF_COUNT_MPI("mpi-recv",COUNT,TYPE) \
     JBPERF_MPE(jb::mpi_recv_stop,0,"recv-stop"); \
     JBPERF_MPI_STOP(REGION);
//----------------------------------------------------------------------
#   define JBPERF_START_MPI_SENDRECV(REGION,SCOUNT,STYPE,RCOUNT,RTYPE) \
     JBPERF_MPE(jb::mpi_sendrecv_start,0,"sendrecv-start"); \
     JBPERF_MPI_START(REGION);
//----------------------------------------------------------------------
#   define JBPERF_STOP_MPI_SENDRECV(REGION,SCOUNT,STYPE,RCOUNT,RTYPE) \
     JBPERF_COUNT_MPI("mpi-send",SCOUNT,STYPE) \
     JBPERF_COUNT_MPI("mpi-recv",RCOUNT,RTYPE) \
     JBPERF_MPE(jb::mpi_sendrecv_stop,0,"sendrecv-stop"); \
     JBPERF_MPI_STOP(REGION);
//----------------------------------------------------------------------
#   define JBPERF_START_MPI_GATHER(REGION,COUNT,TYPE) \
     JBPERF_MPE(jb::mpi_gather_start,0,"gather-start"); \
     JBPERF_MPI_START(REGION);
//----------------------------------------------------------------------
#   define JBPERF_STOP_MPI_GATHER(REGION,COUNT,TYPE) \
     JBPERF_COUNT_MPI("mpi-gather",COUNT,TYPE) \
     JBPERF_MPE(jb::mpi_gather_stop,0,"gather-stop"); \
     JBPERF_MPI_STOP(REGION);
//----------------------------------------------------------------------
#   define JBPERF_START_MPI_REDUCE(REGION,COUNT,TYPE) \
     JBPERF_MPE(jb::mpi_reduce_start,0,"reduce-start"); \
     JBPERF_MPI_START(REGION);
//----------------------------------------------------------------------
#   define JBPERF_STOP_MPI_REDUCE(REGION,COUNT,TYPE) \
     JBPERF_COUNT_MPI("mpi-reduce",COUNT,TYPE) \
     JBPERF_MPE(jb::mpi_reduce_stop,0,"reduce-stop"); \
     JBPERF_MPI_STOP(REGION);
//----------------------------------------------------------------------
#   define JBPERF_START_MPI_BARRIER(REGION) \
     JBPERF_MPE(jb::mpi_barrier_start,0,"barrier-start"); \
     JBPERF_MPI_START(REGION);
//----------------------------------------------------------------------
#   define JBPERF_STOP_MPI_BARRIER(REGION) \
     JBPERF_MPE(jb::mpi_barrier_stop,0,"barrier-stop"); \
     JBPERF_MPI_STOP(REGION);
//----------------------------------------------------------------------
#   define JBPERF_COUNT_MPI(REGION,COUNT,TYPE) \
     { \
       int bytes; \
       MPI_Type_size(TYPE,&bytes); \
       jb::perf.increment(REGION,COUNT*bytes); \
     }
//----------------------------------------------------------------------
#endif

//----------------------------------------------------------------------
// USE_JBPERF_MPI
//----------------------------------------------------------------------

#ifdef USE_JBPERF_MPI
#  define JBPERF_MPI_START(REGION) jb::perf.start(REGION);
#  define JBPERF_MPI_STOP(REGION)  jb::perf.stop(REGION);
#else
#  define JBPERF_MPI_START(REGION) ;
#  define JBPERF_MPI_STOP(REGION)  ;
#endif

//----------------------------------------------------------------------
// jbPerf && MPI && MPE
//----------------------------------------------------------------------

#if defined(USE_JBPERF) && defined(USE_MPI) && defined(USE_MPE)
#   define JBPERF_MPE(EVENT,NUM,NAME) MPE_Log_event(EVENT,NUM,NAME);
#else
#   define JBPERF_MPE(EVENT,NUM,NAME) ;
#endif

//----------------------------------------------------------------------
// jbPerf && HDF5
//----------------------------------------------------------------------

// NOT IMPLEMENTED YET!

#ifdef USE_HDF5
#   define JBPERF_COUNT_READ(DATASET,TYPE,SPACE) ;
#   define JBPERF_COUNT_WRITE(DATASET,TYPE,SPACE) ;
#else
#   define JBPERF_COUNT_READ(DATASET,TYPE,SPACE) ;
#   define JBPERF_COUNT_WRITE(DATASET,TYPE,SPACE) ;
#endif

#endif
