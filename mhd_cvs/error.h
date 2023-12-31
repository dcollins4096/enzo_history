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
#ifndef ERROR_H
#define ERROR_H

#ifdef USE_MPI
#   include <mpi.h>
#   include "message.h"
    static char mpi_error_string[MPI_MAX_ERROR_STRING];
    static int mpi_resultlen;
#   define CHECK_MPI_ERROR(IERR) \
    { \
      int ierr=IERR; \
      if (ierr != MPI_SUCCESS) { \
	MPI_Error_string (ierr,mpi_error_string,&mpi_resultlen); \
	WARNING_MESSAGE; \
	printf ("MPI Error: %s\n",mpi_error_string); \
	fflush(stdout); \
      } \
    }
#else /* USE_MPI */
#   define CHECK_MPI_ERROR(IERR) IERR
#endif

       
#endif

