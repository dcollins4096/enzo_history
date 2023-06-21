/*****************************************************************************
 *                                                                           *
 * Copyright 2004 James Bordner                                              *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
#ifndef JBMEM_H
#define JBMEM_H

//======================================================================
//
// File:        jbMem.h
//
// Description: Simple memory statistics for new / delete
// Description: ( NOTE: overrides system new / delete functions )
//
//----------------------------------------------------------------------
//
// Classes:     jbMem
//
//----------------------------------------------------------------------
//
// Functions: + jbMem ()
// Functions: + ~jbMem ()
//
// Functions: + bytes_global ()         Current bytes on all cpus
// Functions: + bytes_global_high ()    High-water bytes on all cpus
// Functions: + bytes_local ()          Current bytes on local cpu
// Functions: + bytes_local_high ()     High-water bytes on individual cpu
// Functions: + bytes_local_max ()      Maximum bytes_local() of all cpus
// Functions: + bytes_local_max_high () Maximum bytes_local_high() of all cpus
// 
// Functions: + new_count ()
// Functions: + delete_count ()
//
// Functions: + new() (overrides system)
// Functions: + delete() (overrides system)
//
//----------------------------------------------------------------------
//
// Defines: USE_MPI
//
//----------------------------------------------------------------------
//
// James Bordner
// UCSD
//
//======================================================================

#ifdef USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

namespace jb {

  class jbMem {

  public:
    jbMem() : bytes_(0), bytesHigh_(0), newCalls_(0), deleteCalls_(0) {};
    ~jbMem()  {};

    //----------------------------------------------------
    long long bytes_global ()
      {
	long long bytesSum;
#ifdef USE_MPI	
	MPI_Allreduce ((void *)(&bytes_),
		       (void *)(&bytesSum),
		       1,MPI_LONG_LONG, MPI_SUM,MPI_COMM_WORLD);
#else
	bytesSum = bytes_;  
#endif
	return bytesSum;
      };
    //----------------------------------------------------
    long long bytes_global_high ()
      {
	long long bytesHighSum;
#ifdef USE_MPI	
	MPI_Allreduce ((void *)(&bytesHigh_),
		       (void *)(&bytesHighSum),
		       1,MPI_LONG_LONG, MPI_SUM,MPI_COMM_WORLD);
#else
	bytesHighSum = bytesHigh_;  
#endif
	return bytesHighSum;
      };
    long long bytes_local () { return bytes_; };
    //----------------------------------------------------
    long long bytes_local_high () { return bytesHigh_;  };
    //----------------------------------------------------
    long long bytes_local_max () 
      { 
	long long bytesMax;
#ifdef USE_MPI	
	MPI_Allreduce ((void *)(&bytes_),
		       (void *)(&bytesMax),
		       1,MPI_LONG_LONG, MPI_MAX,MPI_COMM_WORLD);
#else
	bytesMax = bytes_;  
#endif
	return bytesMax;
      };
    //----------------------------------------------------
    long long bytes_local_max_high () 
      { 
	long long bytesMax;
#ifdef USE_MPI	
	MPI_Allreduce ((void *)(&bytesHigh_),
		       (void *)(&bytesMax),
		       1,MPI_LONG_LONG, MPI_MAX,MPI_COMM_WORLD);
#else
	bytesMax = bytesHigh_;  
#endif
	return bytesMax;
      };
    //----------------------------------------------------
    long long new_count () { return newCalls_;  };
    long long delete_count () { return deleteCalls_; };

  public:

    long long bytes_;
    long long bytesHigh_;
    long long newCalls_;
    long long deleteCalls_;

  };

  extern jbMem mem;

  // *********************************************************************
}; // namespace jb
// *********************************************************************

#endif
