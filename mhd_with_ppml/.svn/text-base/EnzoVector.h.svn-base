/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Enzo Vector Class, used within linear and nonlinear solvers
/
/  written by: Daniel Reynolds
/  date:       May, 2006
/  modified1:  
/
/  PURPOSE: This class defines a general Enzo vector (one level only).  
/  In addition to containing general 3D local domain data, it contains 
/  a data pointer array.  For each species included in the vector, 
/  this data pointer array contains a pointer to a 3D dataset, 
/  containing the actual vector entries.  Each 3D data set is defined 
/  as generally as possible, such that in each dimension the local 
/  array length is set to be, e.g. Nx0 + Ng0l + Ng0r, where Nx0 
/  represents the number of active local entries, Nx0l represents the 
/  number of ghost cells at the left x0 face, and Nx0r represents the 
/  number of ghost cells at the right x0 face.
/
************************************************************************/

#ifndef ENZO_VECTOR_DEFINED__
#define ENZO_VECTOR_DEFINED__

/* #include <stdio.h> */
/* #include <stdlib.h> */
/* #include <string.h> */
/* #include <math.h> */
/* #include <mpi.h> */

/* #include "macros_and_parameters.h" */
/* #include "typedefs.h" */
/* #include "global_data.h" */

class EnzoVector {
  
 private:

  // local lengths
  int Nx0;        // x0 local active length
  int Nx1;        // x1 local active length
  int Nx2;        // x2 local active length
  int Ng0l;       // x0 left face ghost cells
  int Ng0r;       // x0 right face ghost cells
  int Ng1l;       // x1 left face ghost cells
  int Ng1r;       // x1 right face ghost cells
  int Ng2l;       // x2 left face ghost cells
  int Ng2r;       // x2 right face ghost cells
  int Nspecies;   // number of species in this vector
  int Nglobal;    // total active cells (all procs)
  float **data;   // array of data arrays
  bool owndata;   // flag denoting whether vector or calling routine owns data
  
  // parallelism information
  int Nbors[3][2];  // neighbor procs in each direction 
                    // (set to MPI_PROC_NULL for none)

 public:

  // Constructor
  EnzoVector(int Nx0, int Nx1, int Nx2, 
	     int Ng0l, int Ng0r, int Ng1l, 
	     int Ng1r, int Ng2l, int Ng2r, int Ns,
	     int NBx0L, int NBx0R, int NBx1L, 
	     int NBx1R, int NBx2L, int NBx2R);
  
  // Constructor (overloaded -- sets data to pre-allocated space)
  EnzoVector(int Nx0, int Nx1, int Nx2, 
	     int Ng0l, int Ng0r, int Ng1l, 
	     int Ng1r, int Ng2l, int Ng2r, int Ns,
	     int NBx0L, int NBx0R, int NBx1L, 
	     int NBx1R, int NBx2L, int NBx2R,
	     float **userdata);
  
  // Constructor (overloaded -- sets data arrays to NULL)
  EnzoVector(int Nx0, int Nx1, int Nx2, 
	     int Ng0l, int Ng0r, int Ng1l, 
	     int Ng1r, int Ng2l, int Ng2r, int Ns,
	     int NBx0L, int NBx0R, int NBx1L, 
	     int NBx1R, int NBx2L, int NBx2R, int Empty);
  
  // Destructor
  ~EnzoVector();
  

  // Vector operations

  //   Clone a vector, creating new data arrays
  EnzoVector* clone() const;

  //   Clone a vector, using provided data arrays
  EnzoVector* clone(float **userdata) const;

  //   Writes a given species to file (no ghosts)
  int write(char *outfile, int species) const;

  //   Writes a given species to file (with ghosts)
  int writeall(char *outfile, int species) const;

  //   Communicates ghost cells with neighbors
  int exchange();

  //   Set/Get data array for a given species
  float* GetData(int species);
  int SetData(int species, float *NewArray);

  //   Returns data values by location in 3-space and by species
  float operator()(int i, int j, int k, int s) {
    return data[s][(k*(Nx1+Ng1l+Ng1r)+j)*(Nx0+Ng0l+Ng0r)+i];
  };

  //   Returns dimensional size
  int size(int *n0, int *n1, int *n2, int *ns, int *g0l, 
	   int *g0r, int *g1l, int *g1r, int *g2l, int *g2r);

  //   Copies the values from a vector x (including ghost cells)
  int copy(EnzoVector *x);

  //   Vector linear sum operation, this = a*x + b*y
  int linearsum(float a, EnzoVector *x, float b, EnzoVector *y);

  //   Vector axpy operation, this += a*x
  int axpy(float a, EnzoVector *x);

  //   Vector axpy operation (single component), this += a*x
  int axpy_component(float a, EnzoVector *x, int c);

  //   Vector scale operation, this *= a
  int scale(float a);

  //   Vector constant operation, this(i) = a
  int constant(float a);

  //   Vector add const operation, this(i) += a
  int addconst(float a);

  //   Vector absolute value,  this(i) = |x(i)|
  int abs(EnzoVector *x);

  //   Vector product operation, this(i) = x(i)*y(i)
  int product(EnzoVector *x, EnzoVector *y);

  //   Vector quotient operation, this(i) = x(i)/y(i)
  //   [assumes y(i) != 0]
  int quotient(EnzoVector *x, EnzoVector *y);

  //   Vector minquotient operation, returns
  //   min(this(i)/y(i)) over all y(i)!=0
  float minquotient(EnzoVector *y);

  //   Vector inverse operation, this(i) = 1.0/x(i)
  //   [assumes x(i) != 0]
  int inverse(EnzoVector *x);

  //   Vector constraint checking operation, 
  //   this(i) = 0.0 where constraints true, 1.0 where false
  //   Test:  if c[i] =  2.0, then x[i] must be >  0.0
  //          if c[i] =  1.0, then x[i] must be >= 0.0
  //          if c[i] = -1.0, then x[i] must be <= 0.0
  //          if c[i] = -2.0, then x[i] must be <  0.0
  //   if all constraints satisfied, returns true
  bool constraintcheck(EnzoVector *c, EnzoVector *x);

  //   Vector dot-product,  dot(this,x)
  float dot(EnzoVector *x) const;

  //   Vector RMS norm,  sqrt(dot(this,this)/Nglobal)
  float rmsnorm() const;

  //   Vector weighted RMS norm,  sqrt(dot(this*w,this*w)/Nglobal)
  float wrmsnorm(EnzoVector *w) const;

  //   Vector weighted L-2 norm, sqrt(dot(this*w,this*w))
  float wl2norm(EnzoVector *w) const;

  //   Vector L-1 norm,  sum(abs(this))
  float l1norm() const;

  //   Vector infinity (max) norm,  max(abs(this))
  float infnorm() const;
  
  //   Vector minimum value,  min(this)
  float minval() const;
  
  //   Vector test routine
  int test();

};
  
#endif
