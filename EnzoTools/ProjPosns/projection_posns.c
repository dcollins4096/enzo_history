#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int main(){

  double x,y,z,xfix,yfix,zfix;

  double halfdelta, Nproj, Nrootgrid;

  char *dirstring = ".";
  char *namestring = "../DataDump0033";
  char *fileextension = "project";
  char *datadump = "DataDump0033";

  int level,xaxis,yaxis,zaxis;

  x = 0.494212743390;
  y = 0.497036805755;
  z = 0.494722934236;

  Nproj = 512.0;
  Nrootgrid = 128.0;

  xaxis = 0;
  yaxis = 0;
  zaxis = 1;

  printf("#!/bin/csh\n\n");

  xfix = yfix = zfix = 0.500;

  /* make 512^2 images over a box centered at various levels! */

#define DO_THIS
#ifdef DO_THIS  
  for(level = 2; level <= 6; level += 1)
    {

      halfdelta = Nproj / Nrootgrid * pow( 2.0, -1.0 * ((double) level) - 1.0 );
      
      if(xaxis)
	printf("enzoproj -p x -l %i -b %.20lf %.20lf %.20lf -f %.20lf %.20lf %.20lf -o %s.x.l%i.%s %s\n\n",
	       level,xfix-halfdelta,yfix-halfdelta,zfix-halfdelta,xfix+halfdelta,yfix+halfdelta,zfix+halfdelta,
	       namestring,level,fileextension,datadump);

      if(yaxis)
	printf("enzoproj -p y -l %i -b %.20lf %.20lf %.20lf -f %.20lf %.20lf %.20lf -o %s.y.l%i.%s %s\n\n",
	       level,xfix-halfdelta,yfix-halfdelta,zfix-halfdelta,xfix+halfdelta,yfix+halfdelta,zfix+halfdelta,
	       namestring,level,fileextension,datadump);

      if(zaxis)
	printf("enzoproj -p z -l %i -b %.20lf %.20lf %.20lf -f %.20lf %.20lf %.20lf -o %s.z.l%i.%s %s\n\n",
	       level,xfix-halfdelta,yfix-halfdelta,zfix-halfdelta,xfix+halfdelta,yfix+halfdelta,zfix+halfdelta,
	       namestring,level,fileextension,datadump);

      printf("\n\n");

    }
#endif

  /* make 512^2 images over a box centered at various levels! */
  
  for(level = 6; level <= 28; level += 1)
    {

      halfdelta = Nproj / Nrootgrid * pow( 2.0, -1.0 * ((double) level) - 1.0 );
      
      if(xaxis)
	printf("enzoproj -p x -l %i -b %.20lf %.20lf %.20lf -f %.20lf %.20lf %.20lf -o %s.x.l%i.%s %s\n\n",
	       level,x-halfdelta,y-halfdelta,z-halfdelta,x+halfdelta,y+halfdelta,z+halfdelta,
	       namestring,level,fileextension,datadump);
      if(yaxis)
	printf("enzoproj -p y -l %i -b %.20lf %.20lf %.20lf -f %.20lf %.20lf %.20lf -o %s.y.l%i.%s %s\n\n",
	       level,x-halfdelta,y-halfdelta,z-halfdelta,x+halfdelta,y+halfdelta,z+halfdelta,
	       namestring,level,fileextension,datadump);

      if(zaxis)
	printf("enzoproj -p z -l %i -b %.20lf %.20lf %.20lf -f %.20lf %.20lf %.20lf -o %s.z.l%i.%s %s\n\n",
	       level,x-halfdelta,y-halfdelta,z-halfdelta,x+halfdelta,y+halfdelta,z+halfdelta,
	       namestring,level,fileextension,datadump);

      printf("\n\n");

    }

  return 0;
}
