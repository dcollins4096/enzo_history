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

/* trace flag and trace_level */

//extern int trace;
//extern int trace_level;
//extern FILE *trace_fptr;

#include "flowdefs.h"


void flow1( const char *name )
{
  void print_trace( const char *io, const char *name );

  trace_level = trace_level + 1;
  print_trace("> ",name);
}

void flow2( const char *name )
{
  void print_trace( const char *io, const char *name );

  print_trace("< ",name);
  trace_level = trace_level - 1;
}


void print_trace( const char *io, const char *name )
{

/*
char line[132];
const char pad[2]=".";
int ll;

  if (strlen(name)+strlen(io)+trace_level < 132)
  {

    line[0]=NULL;

    for (ll=0; ll<trace_level; ll++)
    {
      strcat(line,pad);
    }

    strcat(line,io);
    strcat(line,name);

    fprintf(trace_fptr,"%8d : %s\n", trace_level, line);
  }
*/

fprintf(trace_fptr,"%8d %s\n", trace_level, name);

}
