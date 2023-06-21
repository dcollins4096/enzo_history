#include <stdio.h>
#include <math.h>

int main(){

  double x,y,z;

  x=0.534603;
  y=0.511068;
  z=0.490500;

  double halfdelta;

  int level;

  printf("#!/bin/csh\n\n");

  /* make 512^2 images over a box centered at various levels! */
  
  for(level = 6; level <= 10; level+=2)
    {
      halfdelta = 2.0 * pow( 2.0, -1.0 * ((double) level) );
      
      printf("~/EnzoProj/enzoproj -p x -l %i -b %.16lf %.16lf %.16lf -b %.16lf %.16lf %.16lf star_0022\n",
	     level,x-halfdelta/4.0,y-halfdelta,z-halfdelta,x+halfdelta/4.0,y+halfdelta,z+halfdelta);
      printf("mv enzo.project enzo.x.l%i.project\n\n",level);

      /*
      printf("~/EnzoProj/enzoproj -p y -l %i -b %.16lf %.16lf %.16lf -b %.16lf %.16lf %.16lf star_0022\n",
	     level,x-halfdelta,y-halfdelta,z-halfdelta,x+halfdelta,y+halfdelta,z+halfdelta);
      printf("mv enzo.project enzo.y.l%i.project\n\n",level);

      printf("~/EnzoProj/enzoproj -p z -l %i -b %.16lf %.16lf %.16lf -b %.16lf %.16lf %.16lf star_0022\n",
	     level,x-halfdelta,y-halfdelta,z-halfdelta,x+halfdelta,y+halfdelta,z+halfdelta);
      printf("mv enzo.project enzo.z.l%i.project\n\n\n",level);
      */



    }

  return 0;
}
