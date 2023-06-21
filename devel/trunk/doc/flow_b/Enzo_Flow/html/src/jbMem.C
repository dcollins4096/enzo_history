#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "jbMem.h"

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
  jb::mem.bytesMax_ = (jb::mem.bytesMax_ > jb::mem.bytes_)?jb::mem.bytesMax_ : jb::mem.bytes_;
  ++jb::mem.newCalls_;
  long long *pi = (long long *)(malloc(bytes+sizeof(long long)));
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
  jb::mem.bytesMax_ = (jb::mem.bytesMax_ > jb::mem.bytes_)?jb::mem.bytesMax_ : jb::mem.bytes_;
  ++jb::mem.newCalls_;
  long long *pi = (long long *)(malloc(bytes+sizeof(long long)));
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
