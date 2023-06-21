/***********************************************************************
/
/  REPORT MEMORY USAGE (MALLOC, ETC)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/* Defines */

#define NO_OVERLOAD_NEW
#define NO_MALLOC_IRIS4
#define NO_LOCAL_MALLOC

#if defined(MALLOC_IRIS4)
#include <sys/types.h>
#include <malloc.h>
#endif /* MALLOC_IRIS4 */

#ifdef LOCAL_MALLOC
#include "malloc.h"
#endif /* LOCAL_MALLOC */

#ifdef OVERLOAD_NEW
extern int CurrentMemoryUsage;   // in words
extern int MaximumMemoryUsage;
extern int NumberOfCalls     ;
extern void *FirstAddress;
extern void *LargestAddress;
#endif /* OVERLOAD_NEW */

/* Start routine. */

int ReportMemoryUsage(char *header = NULL)
{

#if (defined(MALLOC_IRIS4) || defined(LOCAL_MALLOC))

  static long int MaximumMemoryUsed = 0;

#ifdef MALLOC_IRIS4
  if (MaximumMemoryUsed == 0)
    mallopt(M_MXCHK, 10000);
#endif /* MALLOC_IRIS4 */

  if (header != NULL)
    printf("== %s == ", header);

  /* Use SGI's routines. */
  /* Call mallinfo. */

  struct mallinfo proc;
  proc = mallinfo();
  MaximumMemoryUsed = max(MaximumMemoryUsed, proc.uordblks);
  //  printf("Total space in arena               = %ld\n", proc.arena);
//  printf("Number of ordinary blocks          = %d\n", proc.ordblks);
//  printf("Number of small blocks             = %d\n", proc.smblks);
//  printf("Number of holding blocks           = %d\n", proc.hblkhd);
//  printf("Space in holding block headers     = %d\n", proc.hblkhd);
//  printf("Space in free small blocks         = %d\n", proc.fsmblks);
//  printf("Space in free ordinary blocks      = %d\n", proc.fordblks);
//  printf("Space in used small blocks         = %d\n", proc.usmblks);
  //  printf("Space in used ordinary blocks      = %ld (%f, max = %f)\n", 
  //	 proc.uordblks, float(proc.uordblks)/float(proc.arena),
  //	 float(MaximumMemoryUsed)/float(proc.arena));
  printf("P(%d): Space used, total alloc'ed  = %d %d (%f, max = %f)\n", 
	 MyProcessorNumber, proc.uordblks, proc.arena,
	 float(proc.uordblks)/float(proc.arena),
	 float(MaximumMemoryUsed)/float(proc.arena));
//  printf("Space penalty if keep option used  = %d\n", proc.keepcost);

#endif /* MALLOC_IRIS4 */

  /* ----------------------------------------------------- */

#ifdef OVERLOAD_NEW
  
  //  if (header != NULL)
  //    printf("====== %s ======\n", header);

  float Current = float(CurrentMemoryUsage),
        Maximum = float(MaximumMemoryUsage),
        Arena   = float(((int *) LargestAddress - 
			 (int *) FirstAddress))*sizeof(int);
  printf("%s: P(%d):CurrentMemoryUsage = %g  Max = %g  Arena = %g (%g%%)\n",
	 header, MyProcessorNumber, Current, Maximum, Arena, 
	 Maximum/Arena*100.0);

#endif /* OVERLOAD_NEW */

  return SUCCESS;
}
