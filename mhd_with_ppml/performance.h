/*****************************************************************************
 *                                                                           *
 * Copyright 2006 James Bordner
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Board of Trustees of the University of Illinois            *
 * Copyright 2006 Regents of the University of California                    *
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
// To use jbPerf, compile with -DUSE_JBPERF and link with -ljbperf
// To use PAPI with jbPerf, compile also with -DUSE_PAPI and link with -lpapi
// To instrument HDF5 read/writes with jbPerf, compile with -DUSE_JBPERF_HDF5
// To instrument MPI send/recvs, link with -lpmpi
//
#ifdef PPML
// To use jbMem: -DUSE_JBMEM.  Prints memory usage periodically.
#endif //PPML
//----------------------------------------------------------------------
//
// James Bordner (jbordner@cosmos.ucsd.edu)
// 2003-06-20
//
//======================================================================

#ifdef PPML
//----------------------------------------------------------------------
// jbMem
//----------------------------------------------------------------------
#ifndef JBMEM_LOCAL
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
     printf ("%30s: mem-global-curr    %lld %lld\n", MESSAGE, bcg, bcg - MEMPreviousCurrent);\
     printf ("%30s: mem-global-high    %lld %lld\n", MESSAGE, bmg, bmg - MEMPreviousGlobal); \
     printf ("%30s: mem-local-max-curr %lld\n", MESSAGE, bcl); \
     printf ("%30s: mem-local-max-high %lld\n", MESSAGE, bml); \
           MEMPreviousCurrent=bcg;\
           MEMPreviousGlobal=bmg;\
           fflush(stdout); \
        } \
      }
#endif

#else //JBMEM_LOCAL
//<dbg dcc>
//The Only Local version, for detailed debugging of parallel jobs.
//	bcl = jb::mem.bytes_local_max(); \
//	bml = jb::mem.bytes_local_max_high(); \

#ifndef USE_JBMEM
#  define JBMEM_MESSAGE(ip,MESSAGE) ;
#else

#  include "jbMem.h"
#  define JBMEM_MESSAGE(ip,MESSAGE) \
      { \
	long long bcl,bml; \
	bcl = jb::mem.bytes_local();\
	bml = jb::mem.bytes_local_high();\
     fprintf(stderr,"%30s xp%d: mem-curr %lld %lld\n", MESSAGE, ip, bcl, bcl - MEMPreviousCurrent); \
     fprintf(stderr,"%30s xp%d: mem-high %lld %lld\n", MESSAGE, ip,  bml, bml - MEMPreviousGlobal); \
           MEMPreviousCurrent=bcl;\
           MEMPreviousGlobal=bml;\
           fflush(stderr); \
      }
#endif

#endif //JBMEM_LOCAL
//</dbg dcc>

#endif //PPML

// jbPerf
//----------------------------------------------------------------------

#ifdef USE_JBPERF

#   include "jbPerf.h"

void jbPerfInitialize (int max_level);

#define JB_ITER_PER_SEGMENT 1  /* How frequently to dump data to files. */
#define JB_PERF_EL             /* Whether to include EL?? regions-- */
                               /*   Undefine if jbPerf overhead is too much. */
#define JB_PERF_LEVELS          /* Whether to accumulate jbPerf data for */
                               /*   individual levels--undefine if jbPerf */
                               /*   overhead or memory usage is too much. */

#endif /* USE_JBPERF */

//----------------------------------------------------------------------
// jbPerf && HDF5
//----------------------------------------------------------------------

// Whether to instrument HDF5 read/writes

// NOT IMPLEMENTED YET

#if defined (USE_JBPERF) && defined(USE_JBPERF_HDF5)
#   define JBPERF_HDF5_READ(DATASET) \
    jbPerf.increment ("hdf5-read-calls",1); \
    jbPerf.increment ("hdf5-read-bytes",H5Dget_storage_size (DATASET));
#   define JBPERF_HDF5_WRITE(DATASET) \
    jbPerf.increment ("hdf5-write-calls",1); \
    jbPerf.increment ("hdf5-write-bytes",H5Dget_storage_size (DATASET));
#else
#   define JBPERF_HDF5_READ(DATASET) ;
#   define JBPERF_HDF5_WRITE(DATASET) ;
#endif

//----------------------------------------------------------------------
// jbPerf levels: use lower levels if overhead is too big
//----------------------------------------------------------------------

#ifdef USE_JBPERF

#   if defined (JBPERF_LEVEL_1)   /* Only instrument EvolveLevel as a whole */
#      define JBPERF_START(region)     ;
#      define JBPERF_STOP(region)      ;
#      define JBPERF_START_LOW(region) ;
#      define JBPERF_STOP_LOW(region)  ;
#   endif

#   if defined (JBPERF_LEVEL_2)   /* Instrument major EvolveLevel regions */
#      define JBPERF_START(region)     jbPerf.start (region);
#      define JBPERF_STOP(region)      jbPerf.stop (region);
#      define JBPERF_START_LOW(region) 
#      define JBPERF_STOP_LOW(region)  
#   endif

#   if defined (JBPERF_LEVEL_3)   /* Instrument all EvolveLevel regions */
#      define JBPERF_START(region)     jbPerf.start (region);
#      define JBPERF_STOP(region)      jbPerf.stop (region);
#      define JBPERF_START_LOW(region) jbPerf.start (region);
#      define JBPERF_STOP_LOW(region)  jbPerf.stop (region);
#   endif

#else

#   define JBPERF_START(region)     ;
#   define JBPERF_STOP(region)      ;
#   define JBPERF_START_LOW(region) ;
#   define JBPERF_STOP_LOW(region)  ;

#endif /* USE_JBPERF */

#endif /* PERFORMANCE_H */

