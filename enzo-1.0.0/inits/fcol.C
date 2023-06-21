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
#include <stdio.h>

#include "macros_and_parameters.h"

void fcol(float *x, int n, int m, FILE *log_fptr)
{

  int nrow,mrow;
  int i,j;

  nrow = n/m;
  mrow = n - nrow * m;

  if( mrow > 0 )
  {
    nrow = nrow+1;
  }

  fprintf(log_fptr, "\n");

  for(j=0;j<n;j=j+m)
  {
    for(i=j;i<min(j+m,n);i++)
    {
      fprintf(log_fptr, "%12.4e", x[i]);
    }
    fprintf(log_fptr, "\n");
  }

}
