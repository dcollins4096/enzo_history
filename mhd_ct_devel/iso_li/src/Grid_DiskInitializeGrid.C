#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "pout.h"


int grid::DiskInitializeGrid(float DiskDensity,
			     float DiskPressure,
			     float DiskTheta,
			     float DiskWaveNumber,
			     float DiskRotation,
			     float DiskCompression,
			     float DiskTranslation[],
			     float DiskMagneticField[],
			     float DiskCenter[],
			     float DiskRotationExp,
			     float DiskCompressionExp){
  
#ifndef ATHENA   
   int  EquationOfState = 0 ;
   float IsothermalSoundSpeed = 0;
#endif
  
  
  if( EquationOfState == 0 ){
    NumberOfBaryonFields = 5;
    if(DualEnergyFormalism) NumberOfBaryonFields = 6;
    fprintf(stderr,"GridDim %d %d %d\n",GridDimension[0],GridDimension[1],GridDimension[2]);
    FieldType[0] = Density;
    FieldType[1] = TotalEnergy;
    FieldType[2] = Velocity1;
    FieldType[3] = Velocity2;
    FieldType[4] = Velocity3;
    if(DualEnergyFormalism) FieldType[5] = InternalEnergy;
  }else if( EquationOfState == 1){
    NumberOfBaryonFields = 4;
    if(DualEnergyFormalism){
      NumberOfBaryonFields = 6;
      fprintf(stderr," WARNING!!!! Dual Energy Formalism and Isothermal EOS not checked.\n");
    }
    FieldType[0] = Density;
    FieldType[1] = Velocity1;
    FieldType[2] = Velocity2; 
    FieldType[3] = Velocity3;
    if(DualEnergyFormalism) FieldType[4] = InternalEnergy;
  }
  
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i,j,k, index;
  float xt,yt,zt; //temp coordinates
  float xp,yp,zp; //Rotated coordinates, relative to disk center.
  float vx,vy,vz;
  float ExpZ, ExpP;
  float Pi = 3.14159265;

  int Eeng, Eden, Ev[3], Egas; 
  if (this->IdentifyPhysicalQuantities(Eden, Egas, Ev[0], Ev[1], 
				       Ev[2], Eeng) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  this->AllocateGrids();

  DiskWaveNumber *= 2 * Pi;


  for(k=0;k<GridDimension[2]-1; k++)
  for(j=0;j<GridDimension[1]-1; j++)
    for(i=0;i<GridDimension[0]-1; i++){
      index = i + GridDimension[0]*(j + GridDimension[1]*k);

      BaryonField[ Eden ][index] = DiskDensity;

      
      if( MHD_Used == TRUE ){
	CenteredB[0][index] = DiskMagneticField[0];
	CenteredB[1][index] = DiskMagneticField[1];
	CenteredB[2][index] = DiskMagneticField[2];
      }

      //And now, the velocity.
      xt = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - DiskCenter[0];
      yt = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - DiskCenter[1];
      zt = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - DiskCenter[2];
      xp = xt * cos( DiskTheta ) + zt* sin( DiskTheta );
      yp = yt;
      zp = xt * sin( DiskTheta ) + zt* cos( DiskTheta );

      ExpZ = exp( -1*DiskRotationExp*zp*zp );
      ExpP = exp( -1*DiskCompressionExp*(xp*xp + yp*yp + zp*zp ) );

      //Note the 2pi above on DiskWaveNumber
      vx= 
	(DiskRotation * cos(DiskTheta) * sin( DiskWaveNumber * yp )* ExpZ
	   -DiskCompression * sin( DiskWaveNumber * xp ) * ExpP
	   +DiskTranslation[0]);
      
      vy= (DiskRotation * sin( DiskWaveNumber * xp )* ExpZ
	   -DiskCompression * sin( DiskWaveNumber * yp ) * ExpP
	   +DiskTranslation[1]);
      
      vz= (DiskRotation * sin(DiskTheta) * sin( DiskWaveNumber * yp ) *ExpZ
	   -DiskCompression * sin( DiskWaveNumber * zp ) * ExpP
	   +DiskTranslation[2]);
      
      BaryonField[ Ev[0] ][index] = vx;
      BaryonField[ Ev[1] ][index] = vy;
      BaryonField[ Ev[2] ][index] = vz;
      if( EquationOfState == 0 ) {
	if( MHD_Used == TRUE ){
	  BaryonField[ Eeng ][ index ] = DiskPressure/(Gamma - 1.0) 
	    +0.5 * BaryonField[ Eden ][index] * (vx*vx + vy*vy + vz*vz)
	    +0.5 * (CenteredB[0][index]*CenteredB[0][index] + 
		    CenteredB[1][index]*CenteredB[1][index] +
		    CenteredB[2][index]*CenteredB[2][index]);
	}else{
	  BaryonField[ Eeng ][ index ] = DiskPressure/( (Gamma - 1.0) * Density )+
	    0.5 * (vx*vx + vy*vy + vz*vz);
	}
      }//eos


    }
  return SUCCESS;
  
}
