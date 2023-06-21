/***********************************************************************
/
/  RETURN THE CURRENT ELAPSED CPU TIME
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>

/* Defines */

#if defined(IRIS4)
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <limits.h>
#endif /* IRIS4 */
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "macros_and_parameters.h"

/* Start routine. */

float ReturnCPUTime()
{
  float result = 0;

#ifdef USE_MPI

  result = MPI_Wtime();

#else /* USE_MPI */

#if defined(IRIS4)

  struct tms buffer;
  times(&buffer);
  result = float(buffer.tms_utime)/float(CLK_TCK);

#endif /* IRIS4 */

#endif /* USE_MPI */

  return result;
}

double ReturnWallTime()
{
#ifdef USE_MPI

  return MPI_Wtime();

#else /* USE_MPI */

  return 0.0;

#endif /* USE_MPI */
}
