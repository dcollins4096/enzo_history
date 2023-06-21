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
#include <stdlib.h>
#include "macros_and_parameters.h"
/***********************************************************************
/
/  f77rand
/
/  Portable random number generator for use by F77 programs
/
/  Calls rand() and srand() in ANSI Standard C library 
/
/  void f77srand (int *): seed the random number generator
/  double f77rand ():     return random number 0 <= x < 1
/
/  James Bordner
/
/  2003-05-21  jb  Created
/
***********************************************************************/

extern "C" double FORTRAN_NAME(f77rand) ();
extern "C" void FORTRAN_NAME(f77srand) ();

double FORTRAN_NAME(f77rand) ()
{
  return (double) (rand()) / (RAND_MAX+1);
}

void FORTRAN_NAME(f77srand) (int *seed)
{
  srand(*seed);
}
