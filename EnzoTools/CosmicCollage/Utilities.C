#include "standard_includes.h"
#include "global_variables.h"
#include "subroutines.h"

/* -------------------------------------------------------------------- *
   This is just a wrapper for the integer random number generator used to
   set the axis which we project down for a given projection.
 * -------------------------------------------------------------------- */
int chooseaxis(int proj_number){
  if(debug){ printf("chooseaxis: entering\n"); fflush(stdout); }

  int axis;
  
  if(read_in_random_shifts == 1){

    axis = random_axes[proj_number];

  } else {

    axis = rand()%3;

  }

  if(debug){ printf("chooseaxis: the axis is:  %d\n",axis);  fflush(stdout); }

  if(debug){ printf("chooseaxis: exiting\n"); fflush(stdout); }
  return axis;
}


/* -------------------------------------------------------------------- *
   This is the wrapper for the double-precision random number generator
   we use to calculate the offset for the read-in projection.
 * -------------------------------------------------------------------- */
void chooseoffset(double *xoffset, double *yoffset, int proj_number){
  if(debug){ printf("chooseoffset: entering\n"); fflush(stdout); }


  if(read_in_random_shifts == 1){

    *xoffset = xshifts[proj_number];
    *yoffset = yshifts[proj_number];
    
  } else {

    *xoffset = drand48();
    *yoffset = drand48();

  }

  if(debug){ printf("chooseoffset: exiting\n"); fflush(stdout); }
  return;
}
