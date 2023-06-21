#ifndef __AMR_GRID_H_
#define __AMR_GRID_H_
/*
  Definition of the AMR grid structure for C
  Defined so that the structure will be flat.
*/
typedef struct AmrGrid {
  int level;
  int maxlevel;
  int maxtime,timestep,persistence; /* at finest-level timestep */
  int rank,dims[3]; /* hardcode to 3D maximum dims */
  double delta[3],origin[3]; /* stored in file as delta and origin */
  //
  double minext[3],maxext[3];
  int timerefinement;
  int nbytes;
  int dataveclen;
  int datatype;
  void *data;
} AmrGrid;

#endif
