//======================================================================
//
// File:        c_message.C
//
// Description: Display warning and error messages
//
//----------------------------------------------------------------------
//
// James Bordner
// UCSD
//
//======================================================================
 
#include<stdio.h>
#include<stdlib.h>
 
#ifdef USE_MPI
#include <mpi.h>
#ifdef USE_MPE
#include <mpe.h>
#endif /* USE_MPE */
#endif /* USE_MPI */
 
//----------------------------------------------------------------------
 
int static warning_count = 0;
 
//----------------------------------------------------------------------
 
void c_error (char *sourcefile, int linenumber)
{
#ifdef USE_MPI
  int  ierr;
#endif
  int error_code;
  int id;
 
#ifdef USE_MPI
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &id);
#else
  id = 0;
#endif
 
  printf ("==================\n");
  printf ("=== ENZO ERROR ===   %s: %d   node %d\n",
	  sourcefile,linenumber,id);
  printf ("==================\n");
  fflush(stdout);
 
  error_code = -1;
#ifdef USE_MPI
  ierr = MPI_Abort( MPI_COMM_WORLD, error_code);
#else
  exit(error_code);
#endif
}
 
//----------------------------------------------------------------------
 
void c_warning (char *sourcefile, int linenumber)
{
#ifdef USE_MPI
  int  ierr;
#endif
  int id;
 
  ++ warning_count;
 
#ifdef USE_MPI
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &id);
#else
  id = 0;
#endif
 
  printf ("--- ENZO WARNING #%d ---   %s: %d   node %d\n",
	  warning_count,sourcefile,linenumber,id);
  fflush(stdout);
 
}
