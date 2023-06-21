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
// Functions: + bytes ()
// Functions: + bytes_max ()
// Functions: + bytes_local ()
// Functions: + bytes_max_local ()
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
    jbMem() : bytes_(0), bytesMax_(0), newCalls_(0), deleteCalls_(0) {};
    ~jbMem()  {};

    long long bytes_local () { return bytes_; };
    long long bytes_max_local () { return bytesMax_;  };
    long long bytes ()
      {
	long long bytesGlobal;
#ifdef USE_MPI	
	MPI_Allreduce ((void *)(&bytes_),
		       (void *)(&bytesGlobal),
		       1,MPI_LONG_LONG, MPI_SUM,MPI_COMM_WORLD);
#else
	bytesGlobal = bytes_;  
#endif
	return bytesGlobal;
      };
    long long bytes_max ()
      {
	long long bytesMaxGlobal;
#ifdef USE_MPI	
	MPI_Allreduce ((void *)(&bytesMax_),
		       (void *)(&bytesMaxGlobal),
		       1,MPI_LONG_LONG, MPI_SUM,MPI_COMM_WORLD);
#else
	bytesMaxGlobal = bytesMax_;  
#endif
	return bytesMaxGlobal;
      };

    long long new_count () { return newCalls_;  };
    long long delete_count () { return deleteCalls_; };

  public:

    long long bytes_;
    long long bytesMax_;
    long long newCalls_;
    long long deleteCalls_;

  };

  extern jbMem mem;

  // *********************************************************************
}; // namespace jb
// *********************************************************************

#endif
