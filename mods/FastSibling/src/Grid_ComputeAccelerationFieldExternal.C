/***********************************************************************
/
/  GRID CLASS (ADD FIXED ACCELERATION TO ACCEL FIELDS)
/
/  written by: Greg Bryan
/  date:       June, 1996
/  modified1:
/
/  PURPOSE: Certain problems required external acceleration fields.
/    This routine adds them to the existing self-gravitating fields, or
/     creates new fields if necessary.
/    There is currently support for:
/      1) A uniform gravitational field in one of the orthogonal directions
/         for grid cells only
/      2) A 3D point source field for grid cells only
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

#define Mpc (3.0856e24)         //Mpc [cm] 
#define SolarMass (1.989e33)    //Solar Mass [g]
#define GravConst (6.67e-8)     //Gravitational Constant [cm3g-1s-2]
#define pi (3.14159)
#define mh (1.67e-24)           //Mass of Hydrogen [g]
#define kboltz (1.381e-16)      //Boltzmann's Constant [ergK-1]

/* function prototypes */

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int grid::ComputeAccelerationFieldExternal()
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int dim, i, j, k, size = 1;

  /* Compute field size (in floats). */

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Check if acceleration field exists.  If not create it and zero it. */

  if (AccelerationField[0] == NULL)
    for (dim = 0; dim < GridRank; dim++) {
      AccelerationField[dim] = new float[size];
      for (i = 0; i < size; i++)
	AccelerationField[dim][i] = 0;
    }

  /* -----------------------------------------------------------------
     Point Source gravity
     ----------------------------------------------------------------- */

  if (PointSourceGravity) {

    FLOAT a = 1, accel, dadt, radius, rcubed, xpos, ypos = 0, zpos = 0, rcore;

    /* Compute adot/a at time = t+1/2dt (time-centered). */

    float DensityUnits = 1 , LengthUnits = 1, TemperatureUnits = 1, 
          TimeUnits = 1, VelocityUnits = 1, AccelUnits = 1;
    if (ComovingCoordinates) {

      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	  == FAIL) {
	fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
	return FAIL;
      }

      CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			&TimeUnits, &VelocityUnits, Time);
      AccelUnits = LengthUnits/TimeUnits/TimeUnits;
    } // end: if comoving coordinates

    /* Loop over grid, adding acceleration to field. */

    for (dim = 0; dim < GridRank; dim++) {
      int n = 0;

      for (k = 0; k < GridDimension[2]; k++) {
	if (GridRank > 2) 
	  zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - 
	    PointSourceGravityPosition[2];
	if (dim == 2 && HydroMethod == Zeus_Hydro) 
	  zpos -= 0.5*CellWidth[2][k];

	for (j = 0; j < GridDimension[1]; j++) {
	  if (GridRank > 1) 
	    ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] -
	      PointSourceGravityPosition[1];
	  if (dim == 1 && HydroMethod == Zeus_Hydro) 
	    ypos -= 0.5*CellWidth[1][j];

	  for (i = 0; i < GridDimension[0]; i++, n++) {
	    xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] -
	      PointSourceGravityPosition[0];
	    if (dim == 0 && HydroMethod == Zeus_Hydro) 
	      xpos -= 0.5*CellWidth[0][i];
	    
	    /* Compute distance from center. */

	    radius = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
	    rcubed = pow(radius, 3);

	    /* (1) is a real (softened) point-source, (2) is NFW profile */

	    if (PointSourceGravity == 1) {

	      /* (1) Point Source:
		 Multiply by a(t) to offset the 1/a(t) in ComovingAccelTerm(). 
		 (i.e. 1/a^2 * a = 1/a). */

	      rcore = max(0.1*CellWidth[0][0], PointSourceGravityCoreRadius);
	      rcubed += rcore*rcore*rcore;
	      accel = PointSourceGravityConstant/(rcubed*a);
	    }
	    else {

	      /* (2) NFW Profile: assume CoreRadius is rs in cm and Constant
		 is mass within rs in g. */

	      rcore = PointSourceGravityCoreRadius/LengthUnits;
	      FLOAT x = radius/rcore;
	      accel = GravConst*PointSourceGravityConstant*
		      ((log(1.0+x  )-x  /(1.0+x  )) /
		       (log(1.0+1.0)-1.0/(1.0+1.0))) / 
   		      pow(radius*LengthUnits, 2) / AccelUnits;
	      accel = accel/radius;  // this radius normalizes the multiplication by xpos,ypos,zpos done below
	      //	      printf("%g %g %g a=%g units=%g r=%g\n", xpos, ypos, zpos, 
	      //		     accel, AccelUnits, radius);
	    }
	    
	    /* Apply force. */
	    
	    if (dim == 0)
	      AccelerationField[0][n] -= accel*xpos;
	    if (dim == 1)
	      AccelerationField[1][n] -= accel*ypos;
	    if (dim == 2)
	      AccelerationField[2][n] -= accel*zpos;

	  }
	}
      } // end: loop over grid
    } // end: loop over dims

    /* DO PARTICLES HERE! */

    if (NumberOfParticles > 0 && GridRank != 3) {
      fprintf(stderr, "PointSourceGravity assumes 3D\n");
      return FAIL;
    }
    
    for (i = 0; i < NumberOfParticles; i++) {

      /* Compute vector between particle (advanced by 1/2 step) and
	 gravity center. */
      
      xpos = ParticlePosition[0][i] + 0.5*dtFixed/a - 
	PointSourceGravityPosition[0];
      ypos = ParticlePosition[1][i] + 0.5*dtFixed/a - 
	PointSourceGravityPosition[1];
      zpos = ParticlePosition[2][i] + 0.5*dtFixed/a - 
	PointSourceGravityPosition[2];
	    
      /* Compute distance from center. */

      radius = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
      rcubed = pow(radius, 3);

      /* (1) is a real (softened) point-source, (2) is NFW profile */

      if (PointSourceGravity == 1) {

	/* (1) Point Source:
	   Multiply by a(t) to offset the 1/a(t) in ComovingAccelTerm(). 
	   (i.e. 1/a^2 * a = 1/a). */

	rcore = max(0.1*CellWidth[0][0], PointSourceGravityCoreRadius);
	rcubed += rcore*rcore*rcore;
	accel = PointSourceGravityConstant/(rcubed*a);
      }
      else {

	/* (2) NFW Profile: assume CoreRadius is rs in cm and Constant
	   is mass within rs in g. */

	rcore = PointSourceGravityCoreRadius/LengthUnits;
	FLOAT x = radius/rcore;
	accel = GravConst*PointSourceGravityConstant*
		      ((log(1.0+x  )-x  /(1.0+x  )) /
		       (log(1.0+1.0)-1.0/(1.0+1.0))) / 
   		      pow(radius*LengthUnits, 2) / AccelUnits;
	accel = accel/radius;  // this radius normalizes the multiplication by xpos,ypos,zpos done below
	      //	      printf("%g %g %g a=%g units=%g r=%g\n", xpos, ypos, zpos, 
	      //		     accel, AccelUnits, radius);
	    }
	    
      /* Apply force. */
	    
      ParticleAcceleration[0][i] -= accel*xpos;
      ParticleAcceleration[1][i] -= accel*ypos;
      ParticleAcceleration[2][i] -= accel*zpos;

    } // end: loop over number of particles

  } // end: if (PointSourceGravity)

  /* -----------------------------------------------------------------
     Uniform gravity field
     ----------------------------------------------------------------- */

  if (UniformGravity) {

//    if (debug) 
//      printf("grid->ComputeAccelerationFieldExternal: dir, g = %d, %g\n",
//	     UniformGravityDirection, UniformGravityConstant);
		      
    for (dim = 0; dim < GridRank; dim++) {

      /* Set constant for this dimension. */

      float Constant = 0.0;
      if (dim == UniformGravityDirection)
	Constant = UniformGravityConstant;
	
      /* Set field. */

      for (i = 0; i < size; i++)
	AccelerationField[dim][i] = Constant;

    } // loop over dims

  } // end: if (UniformGravity)

  return SUCCESS;
}

