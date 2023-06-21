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
#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "jbMem.h"
#include "message.h"

// #define TRACE

//----------------------------------------------------------------------

void *operator new (size_t bytes) throw (std::bad_alloc)
  // WARNING!!! OVERRIDES new () in <new>
{
#ifdef TRACE
  printf ("jbMem new ()\n");
#endif
  
  if (bytes==0) return NULL;
  jb::mem.bytes_ += bytes;
  jb::mem.bytesHigh_ = (jb::mem.bytesHigh_ > jb::mem.bytes_)?jb::mem.bytesHigh_ : jb::mem.bytes_;
  ++jb::mem.newCalls_;
  long long *pi = (long long *)(malloc(bytes+sizeof(long long)));
  if (pi==0) {
    fprintf (stderr,"jbMem error: out of memory!\n");
    ERROR_MESSAGE;
  }
  pi[0] = bytes;
  return (void *)(pi+1);
}

//----------------------------------------------------------------------

void operator delete (void *p) throw ()
  // WARNING!!! OVERRIDES delete () in <new>
{
#ifdef TRACE
  printf ("jbMem delete ()\n");
#endif
  if (p==0) return;
  long long *pi = (long long *)(p)-1;
  jb::mem.bytes_ -= pi[0];
  ++jb::mem.deleteCalls_;
  free(pi);
}

//----------------------------------------------------------------------

void *operator new [] (size_t bytes) throw (std::bad_alloc)
  // WARNING!!! OVERRIDES new () in <new>
{
#ifdef TRACE
  printf ("jbMem new []()\n");
#endif
  if (bytes==0) return NULL;
  jb::mem.bytes_ += bytes;
  jb::mem.bytesHigh_ = (jb::mem.bytesHigh_ > jb::mem.bytes_)?jb::mem.bytesHigh_ : jb::mem.bytes_;
  ++jb::mem.newCalls_;
  long long *pi = (long long *)(malloc(bytes+sizeof(long long)));
  if (pi==0) {
    fprintf (stderr,"jbMem error: out of memory!\n");
    ERROR_MESSAGE;
  }
  pi[0] = bytes;
  return (void *)(pi+1);
}

//----------------------------------------------------------------------

void operator delete [] (void *p) throw ()
  // WARNING!!! OVERRIDES delete () in <new>
{
#ifdef TRACE
  printf ("jbMem delete []()\n");
#endif
  if (p==0) return;
  long long *pi = (long long *)(p)-1;
  jb::mem.bytes_ -= pi[0];
  ++jb::mem.deleteCalls_;
  free(pi);
}

namespace jb {

  jbMem mem;

};
