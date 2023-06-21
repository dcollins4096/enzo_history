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
/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID TO A UNIFORM POOL OF GAS)
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
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

extern "C" void FORTRAN_NAME(curl_of_e)(float *bx, float *by, float *bz,
                                        float *ex, float *ey, float *ez,
                                        float *dx, float *dy, float *dz,
                                        int *idim, int *jdim, int *kdim,
                                        int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
                                        float *dt, int *method);


#ifdef HAOXU
int grid::InitializeGridwithSphere( int NumberofSpheres, int SphereType[MAX_SPHERES],
                                float SpherePosition[MAX_SPHERES][MAX_DIMENSION],
                                float CoreDensity[MAX_SPHERES], 
                                float SphereRadius[MAX_SPHERES],
                                float SphereCoreRadius[MAX_SPHERES],
                                float Sphere_n[MAX_SPHERES],
                                float UniformDensity,
				float UniformPressure,
				float UniformVelocity[MAX_DIMENSION],
				float *UniformMagneticField) // default = NULL
{
  /* declarations */

  int dim, i,j,k, size, field, sphere, index;
  float x1;
 
  float One=1;
  int  Two=2;
  float Pi = 3.141592654;

  float Scale[3];

  /* create fields */
  
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* compute size of fields */

  size = 1;
  for (dim = 0; dim < GridRank; dim++){
    size *= GridDimension[dim];
    Scale[dim]=(GridRightEdge[dim]-GridLeftEdge[dim])/(GridDimension[dim]-2*DEFAULT_GHOST_ZONES);
    }

  /* allocate fields */

  this->AllocateGrids();

   /* Loop over the mesh. */

  float density, dens1, Velocity[MAX_DIMENSION];
  FLOAT r, x=0, y = 0, z = 0;
  int n = 0;

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {

        /* Compute position */

        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        if (GridRank > 1)
          y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        if (GridRank > 2)
          z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

        /* Loop over spheres. */

        density = UniformDensity;
        for (dim = 0; dim < MAX_DIMENSION; dim++)
          Velocity[dim] = UniformVelocity[dim];
        for (sphere = 0; sphere < NumberofSpheres; sphere++) {

          /* Find distance from center. */

          r = sqrt(POW(fabs(x-SpherePosition[sphere][0]), 2) +
                   POW(fabs(y-SpherePosition[sphere][1]), 2)
    +                  POW(fabs(z-SpherePosition[sphere][2]), 2) );
//          r = max(r, 0.3*CellWidth[0][0]);
        if(SphereType[sphere]==8) r=sqrt(POW(fabs(x-SpherePosition[sphere][0]), 2) +
                   POW(fabs(y-SpherePosition[sphere][1]), 2) );

          if (r < SphereRadius[sphere]||1) {

            /* 1) Uniform */

            if (SphereType[sphere] == 1)
              dens1 = CoreDensity[sphere];

            /* 2) r^-2 power law */

            if (SphereType[sphere] == 2)
              dens1 = CoreDensity[sphere]*POW(r/SphereRadius[sphere], -2);

            /* 3) Kornreich & Scalo 2000's setting*/

            if (SphereType[sphere] == 3) {
              x1 = r/SphereRadius[sphere];
              dens1 = UniformDensity+(CoreDensity[sphere]-UniformDensity)/(1+POW(x1,Sphere_n[sphere]));
            }

           /* 4) Gaussian */

            if (SphereType[sphere] == 4) {
              dens1 = CoreDensity[sphere]*
                      exp(-0.5*POW(r/SphereRadius[sphere], 2));
            }

            /* 5) r^-2 power law with core radius */

/*            if (SphereType[sphere] == 5) {
              if (r < SphereCoreRadius[sphere])
                dens1 = CoreDensity[sphere]*POW(SphereCoreRadius[sphere]/
                                                  CoreRadius[sphere], -2);
              else
                dens1 = CoreDensity[sphere]*POW(r/CoreRadius[sphere], -2);
            }
*/
             /* 2D setting*/
              if(SphereType[sphere] == 8){
              x1 = r/SphereRadius[sphere];
              dens1 = UniformDensity+(CoreDensity[sphere]-UniformDensity)/(1+POW(x1,Sphere_n[sphere]));
            }

             /* 10) disk (ok, it's not a sphere, so shoot me) */

            if (SphereType[sphere] == 10) {
              float   SphereVelocity[1][1];
              
              FLOAT xpos, ypos, zpos, xpos1, ypos1, zpos1, zheight, drad;
              FLOAT ScaleHeightz = SphereCoreRadius[sphere]/6.0,
                    ScaleHeightR = SphereCoreRadius[sphere];

              /* Loop over dims if using Zeus (since vel's face-centered). */

              for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
                   dim++) {

                /* Compute position. */

                xpos = x-SpherePosition[sphere][0] -
                  (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
                ypos = y-SpherePosition[sphere][1] -
                  (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
                zpos = z-SpherePosition[sphere][2] -
                  (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);

                  /* Compute z and r_perp (SphereVelocity is angular momentum
                   and must have unit length). */

                zheight = SphereVelocity[sphere][0]*xpos +
                          SphereVelocity[sphere][1]*ypos +
                          SphereVelocity[sphere][2]*zpos;
                xpos1 = xpos - zheight*SphereVelocity[sphere][0];
                ypos1 = ypos - zheight*SphereVelocity[sphere][1];
                zpos1 = zpos - zheight*SphereVelocity[sphere][2];
                drad = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1);

                /* If we're above the disk, then exit. */

//              if (zheight > max(5.0*ScaleHeightz, 2.0*CellWidth[0][0]))
//                continue;

                /* Compute density (Kruit & Searle 1982). */

                if (dim == 0)
                  dens1 = CoreDensity[sphere]*exp(-drad/ScaleHeightR)/
                    POW(cosh(zheight/max(ScaleHeightz, CellWidth[0][0])), 2);

                if (dens1 < density)
                  break;

                 /* Compute velocity magnitude (divided by drad).
                   This assumes PointSourceGravityPosition and Sphere center
                   are the same.  This should be fixed to use the disk mass
                   as well, but that's a bit tricky. */

//              float vel = sqrt(PointSourceGravityConstant/drad)/drad;

                float accel = 0;
                if (PointSourceGravity == 1)
                  accel = PointSourceGravityConstant/
                    (POW(drad,3) + POW(PointSourceGravityCoreRadius, 3));
                if (PointSourceGravity == 2) {
                  FLOAT x = drad/PointSourceGravityCoreRadius;
                  accel = PointSourceGravityConstant*(log(1+x)-x/(1+x))/
                           POW(drad, 3);
                }

                float vel = sqrt(accel);

                /* Compute velocty: L x r_perp. */

                if (dim == 0 || dim == 1)
                  Velocity[0] = vel*(SphereVelocity[sphere][1]*zpos1 -
                                     SphereVelocity[sphere][2]*ypos1);
                if (dim == 0 || dim == 2)
                  Velocity[1] = vel*(SphereVelocity[sphere][2]*xpos1 +
                                     SphereVelocity[sphere][0]*zpos1);
                if (dim == 0 || dim == 3)
                  Velocity[2] = vel*(SphereVelocity[sphere][0]*ypos1 -
                                     SphereVelocity[sphere][1]*xpos1);

              } // end: loop over dims

            } // end: disk

 
            /* If the density is larger than the background (or the previous
               sphere), then set the velocity. */

            if (dens1 > density) 
              density = dens1;

          } // end: if (r < SphereRadius)
        } // end: loop over spheres


  /* set density */

    BaryonField[0][n] = density;

  /* set velocities */

  for (dim = 0; dim < GridRank; dim++)
      BaryonField[vel+dim][n] = UniformVelocity[dim];


  } //Loop over meshs

//Initialize magnetic fields as curl of vector potnetial
   sphere = 0;
  if(MHD_Used) {
    float r_temp;
    float B0=1.0e-10; //0.001;
    for(field=0;field<3;field++)
     for( k=0;k<ElectricDims[field][2];k++)
    for( j=0;j<ElectricDims[field][1];j++)
      for( i=0;i<ElectricDims[field][0];i++){
        index=i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k);

         /* Compute position */

        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        if (GridRank > 1)
          y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        if (GridRank > 2)
          z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
       
    r_temp=sqrt(pow(x-SpherePosition[sphere][0],2)+pow(y-SpherePosition[sphere][1],2)+pow(z-SpherePosition[sphere][2],2));
    if(field==0) ElectricField[field][index]=-B0*(x-SpherePosition[sphere][0])*exp(-r_temp*r_temp);
    if(field==1) ElectricField[field][index]=B0*(y-SpherePosition[sphere][1])*exp(-r_temp*r_temp);
    if(field==2) ElectricField[field][index]=0.0;  
  }


  //
  //Curl is B=Curl(A)
  //


  FORTRAN_NAME(curl_of_e)(MagneticField[0], MagneticField[1], MagneticField[2],
                          ElectricField[0], ElectricField[1], ElectricField[2],
                          CellWidth[0], CellWidth[1], CellWidth[2],
                          GridDimension, GridDimension +1, GridDimension +2,
                          GridStartIndex, GridEndIndex,
                          GridStartIndex+1, GridEndIndex+1,
                          GridStartIndex+2, GridEndIndex+2,
                          &One, &Two);


  if( this->CenterMagneticField() == FAIL ) {
    fprintf(stderr," error with CenterMagneticField\n");
    return FAIL;
  }

  }//MHD_Used

/*
  if( MHD_Used ){
    for(field=0;field<3;field++){
      for(i=0;i<MagneticSize[field];i++)
	MagneticField[field][i] = 
	  ((UniformMagneticField != NULL)?UniformMagneticField[field]:1.0e-10);
      for(i=0;i<size;i++)
	CenteredB[field][i]= 
	  ((UniformMagneticField != NULL)?UniformMagneticField[field]:1.0e-10);
    }      
  }
*/
  /* Set Total Energy */

 if(MHD_Used){
   for (i=0; i < size; i++){
     BaryonField[1][i] =  (UniformPressure- 0.5*(CenteredB[0][i]*CenteredB[0][i]+CenteredB[1][i]*CenteredB[1][i]+CenteredB[2][i]*CenteredB[2][i]))/(Gamma-1.0) + 0.5*BaryonField[0][i]*
        (BaryonField[vel+0][i]*BaryonField[vel+0][i]+BaryonField[vel+1][i]*BaryonField[vel+1][i]+BaryonField[vel+2][i]*BaryonField[vel+2][i])
       + 0.5*(CenteredB[0][i]*CenteredB[0][i]+CenteredB[1][i]*CenteredB[1][i]+CenteredB[2][i]*CenteredB[2][i]);

   /* Set internal energy if necessary. */

  if (DualEnergyFormalism)
      BaryonField[2][i] = (UniformPressure- 0.5*(CenteredB[0][i]*CenteredB[0][i]+CenteredB[1][i]*CenteredB[1][i]+CenteredB[2][i]*CenteredB[2][i]))/(BaryonField[0][i]*(Gamma-1.0));

}//i
}else{
  for(i=0;i<size;i++){
  BaryonField[1][i] = UniformPressure/(BaryonField[0][i]*(Gamma-1.0)) + 0.5*
        (BaryonField[vel+0][i]*BaryonField[vel+0][i]+BaryonField[vel+1][i]*BaryonField[vel+1][i]+BaryonField[vel+2][i]*BaryonField[vel+2][i]);

   /* Set internal energy if necessary. */

  if (DualEnergyFormalism)
      BaryonField[2][i] = UniformPressure/(BaryonField[0][i]*(Gamma-1.0));
  }
 }


sphere=0;
n=0;
 if(PointSourceGravity == 1){
for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++,n++) {
          x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        if (GridRank > 1)
          y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        if (GridRank > 2)
          z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];   

          /* Find distance from center. */

          r = sqrt(POW(fabs(x-SpherePosition[sphere][0]), 2) +
                   POW(fabs(y-SpherePosition[sphere][1]), 2));
//   +                 POW(fabs(z-SpherePosition[sphere][2]), 2) );
          r = max(r, 1.0*CellWidth[0][0]);

              if(Sphere_n[sphere]==2.0)
              BaryonField[1][n] =PointSourceGravityConstant *
                  (CoreDensity[sphere]/r  + (UniformDensity-CoreDensity[sphere])*atan(SphereRadius[sphere]/r)/SphereRadius[sphere])/((Gamma-1.0)*BaryonField[0][n]);
    if(Sphere_n[sphere]==8.0)
    BaryonField[1][n] = PointSourceGravityConstant *
   (1.0/8)*(8*CoreDensity[sphere]/r -
((UniformDensity-CoreDensity[sphere])*cos(Pi/8)*log(1 + pow(SphereRadius[sphere]/r,2) - 2*SphereRadius[sphere]/r*cos(Pi/8)))/SphereRadius[sphere] -
((UniformDensity-CoreDensity[sphere])*cos((3*Pi)/8)*log(1 + pow(SphereRadius[sphere]/r,2) - 2*SphereRadius[sphere]/r*cos((3*Pi)/8)))/SphereRadius[sphere] - 
   ((UniformDensity-CoreDensity[sphere])*cos((5*Pi)/8)*log(1 + pow(SphereRadius[sphere]/r,2) - 2*SphereRadius[sphere]/r*cos((5*Pi)/8)))/SphereRadius[sphere] -
((UniformDensity-CoreDensity[sphere])*cos((7*Pi)/8)*log(1 + pow(SphereRadius[sphere]/r,2) - 2*SphereRadius[sphere]/r*cos((7*Pi)/8)))/SphereRadius[sphere] - 
   (2*(UniformDensity-CoreDensity[sphere])*atan(1/tan(Pi/8) - SphereRadius[sphere]/r/sin(Pi/8))*sin(Pi/8))/SphereRadius[sphere] +
(2*(UniformDensity-CoreDensity[sphere])*atan((SphereRadius[sphere]/r - cos((3*Pi)/8))/sin((3*Pi)/8))*sin((3*Pi)/8))/SphereRadius[sphere] + 
   (2*(UniformDensity-CoreDensity[sphere])*atan((SphereRadius[sphere]/r - cos((5*Pi)/8))/sin((5*Pi)/8))*sin((5*Pi)/8))/SphereRadius[sphere] + 
   (2*(UniformDensity-CoreDensity[sphere])*atan((SphereRadius[sphere]/r - cos((7*Pi)/8))/sin((7*Pi)/8))*sin((7*Pi)/8))/SphereRadius[sphere]);

    if (DualEnergyFormalism)
   BaryonField[2][n] = BaryonField[1][n]/((Gamma-1.0)*BaryonField[0][n]);
 
       BaryonField[1][n] = BaryonField[1][n]/(Gamma-1.0)+0.5*BaryonField[0][n]*
        (BaryonField[vel+0][n]*BaryonField[vel+0][n]+BaryonField[vel+1][n]*BaryonField[vel+1][n]+BaryonField[vel+2][n]*BaryonField[vel+2][n])
       + 0.5*(CenteredB[0][n]*CenteredB[0][n]+CenteredB[1][n]*CenteredB[1][n]+CenteredB[2][n]*CenteredB[2][n]);

}//i
}//PointSourceGravity

  

  return SUCCESS;
}
#endif //HAOXU
