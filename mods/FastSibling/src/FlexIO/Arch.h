#ifndef _MACHDEFS_H_
#define _MACHDEFS_H_

#ifdef ANSI
#define PROTO(x) x
#else
#define PROTO(x) ()
#endif

#if defined(SGI) || defined(CM5) || defined (DEC) 
#define F77NAME(a,b,c) a
#else
#ifdef HP
#define  F77NAME(a,b,c) b
#else
#ifdef CRAY
#define F77NAME(a,b,c) c
#endif
#endif
#endif

#ifndef F77NAME
#define F77NAME(a,b,c) a
#endif


#ifndef SGI
typedef long Long8;
#else /* default SGI behavior.  long longs=8 bytes long=4 bytes*/
typedef long long Long8;
#endif

#ifdef T3E
typedef short Int; /* T3E uses 8-byte integers as default */
typedef short Long; /* this is wierd, but its a T3E... */
#else
typedef int Int;  /* every other sane design uses 4-byte ints */
typedef int Long; /* for now Long *MUST* be 32bit integer for PC compatability */
/* if we run into problems later we'll just need have a separate rev of the file format for PC's */
#endif

#ifdef WIN32 /* this aint a happenin thang yet... */
/* #include <sys/types.h> ... bahh!!! Win32 doesn't have this! */
typedef unsigned int IOFile; /* cast to integer for pointers :( */
union Integer8 {
  long l;
  int i[2];
  char c[8];
}; /* what can be said about the byte order though? */
#else /* its a Mac or Unix box probably */
#include <sys/types.h> 
typedef caddr_t IOFile; /* use caddr_t to store pointers */
#endif

#ifdef HP
#define NO_RECURSIVE_INLINE
#endif

#ifndef CONST_BROKEN
#define CONST const
#else
#ifdef CONST
#undef CONST
#endif
#endif

#endif
