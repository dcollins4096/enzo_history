/***********************************************************************
/
/  GRID CLASS (COMPUTE MHD SOURCE TERMS)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

#define USE
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int field, int farray[], int numfields);


int grid::MHDSourceTerms(float **dU)
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
    B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
				   TENum, B1Num, B2Num, B3Num, PhiNum);




#ifdef USE
  /* Dedner MHD formulation source terms */

  FLOAT dtdx = dtFixed/CellWidth[0][0],
    dtdy = (GridRank > 1) ? dtFixed/CellWidth[1][0] : 0.0,
    dtdz = (GridRank > 2) ? dtFixed/CellWidth[2][0] : 0.0;  
  float Bx, By, Bz;
  float coeff = 1.;

  //  if (EOSType == 3)  coeff = 0.; // turn of adding dissipated B-field to Etot if isothermal (worth a try ...)

  int n = 0, igrid, igridyp1, igridym1, igridzp1, igridzm1;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	igrid = i+(j+k*GridDimension[1])*GridDimension[0];

	Bx = BaryonField[B1Num][igrid];
	By = BaryonField[B2Num][igrid];
	Bz = BaryonField[B3Num][igrid];
	/*
	igridyp1 = i+(j+1+k*GridDimension[1])*GridDimension[0];
	igridym1 = i+(j-1+k*GridDimension[1])*GridDimension[0];
	igridzp1 = i+(j+(k+1)*GridDimension[1])*GridDimension[0];
	igridzm1 = i+(j+(k-1)*GridDimension[1])*GridDimension[0];
	divB[n] = 0.5*(BaryonField[B1Num][igrid+1]-BaryonField[B1Num][igrid-1])*dtdx +
	  0.5*(BaryonField[B2Num][igridyp1]-BaryonField[B2Num][igridym1])*dtdy +
	  0.5*(BaryonField[B3Num][igridzp1]-BaryonField[B3Num][igridzm1])*dtdz; 
	gradPhi[0][n] = 0.5*(BaryonField[PhiNum][igrid+1]-BaryonField[PhiNum][igrid-1])*dtdx;
	gradPhi[1][n] = 0.5*(BaryonField[PhiNum][igridyp1]-BaryonField[PhiNum][igridym1])*dtdy;
	gradPhi[2][n] = 0.5*(BaryonField[PhiNum][igridzp1]-BaryonField[PhiNum][igridzm1])*dtdz;
	*/ 
	dU[iS1  ][n] -= divB[n]*Bx * coeff;
	dU[iS2  ][n] -= divB[n]*By * coeff;
	dU[iS3  ][n] -= divB[n]*Bz * coeff;
	dU[iEtot][n] -= coeff * (Bx*gradPhi[0][n] + By*gradPhi[1][n] + Bz*gradPhi[2][n]);


      }
    }
  }
#endif

  if (DualEnergyFormalism) {
    int igrid, ip1, im1, jp1, jm1, kp1, km1;
    FLOAT dtdx = 0.5*dtFixed/CellWidth[0][0],
      dtdy = (GridRank > 1) ? 0.5*dtFixed/CellWidth[1][0] : 0.0,
      dtdz = (GridRank > 2) ? 0.5*dtFixed/CellWidth[2][0] : 0.0;
    float min_coeff = 0.0;
    if (UseMinimumPressureSupport) {
      min_coeff = MinimumPressureSupportParameter*
	0.32*pow(CellWidth[0][0],2)/(Gamma*(Gamma-1.0));
    }
    float rho, eint, p, divVdt, h, cs, dpdrho, dpde;
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  ip1 = igrid + 1;
	  im1 = igrid - 1;
	  jp1 = (GridRank > 1) ? i + GridDimension[0]*(j + 1 + k*GridDimension[1]) : 0;
	  jm1 = (GridRank > 1) ? i + GridDimension[0]*(j - 1 + k*GridDimension[1]) : 0;
	  kp1 = (GridRank > 2) ? i + GridDimension[0]*(j + (k+1)*GridDimension[1]) : 0;
	  km1 = (GridRank > 2) ? i + GridDimension[0]*(j + (k-1)*GridDimension[1]) : 0;
	  divVdt = dtdx*(BaryonField[Vel1Num][ip1] - BaryonField[Vel1Num][im1]) +
	    dtdy*(BaryonField[Vel2Num][jp1] - BaryonField[Vel2Num][jm1]) +
	    dtdz*(BaryonField[Vel3Num][kp1] - BaryonField[Vel3Num][km1]);
	  rho = BaryonField[DensNum][igrid];
	  eint = BaryonField[GENum][igrid];
	  eint = max(eint, min_coeff*rho);
	  EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
	  dU[iEint][n] -= p*divVdt;

	}
      }
    }
  }

  if (Coordinate == Cylindrical) {
    float rho, etot, eint, vx, vy, vz, v2, e, h, cs, p, 
      dpdrho, dpde, coty, Bx, By, Bz, B2;
    FLOAT x, dtxinv;
    int n = 0, igrid;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
          igrid = i+(j+k*GridDimension[1])*GridDimension[0];
          x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];

          rho = BaryonField[DensNum][igrid];
          vx  = BaryonField[Vel1Num][igrid];
          vy  = BaryonField[Vel2Num][igrid];
          vz  = BaryonField[Vel3Num][igrid];
	  Bx  = BaryonField[B1Num  ][igrid];
	  By  = BaryonField[B2Num  ][igrid];
	  Bz  = BaryonField[B3Num  ][igrid];
	  if (DualEnergyFormalism) {
	    eint = BaryonField[ieint][igrid];
	  }
	  else {
	    etot = BaryonField[TENum][igrid];
	    v2 = vx*vx + vy*vy + vz*vz;
	    B2 = Bx*Bx + By*By + Bz*Bz;
	    eint = etot - 0.5*v2 - 0.5*B2/rho;
	  }
                  
          EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
         
          dtxinv = dtFixed/x;
          dU[iS1][n]  += dtxinv*(p + rho*vz*vz);
          dU[iS3][n]  += -dtxinv*rho*vx*vz;

	
        }
      }
    }
  }




  if (UseConstantAcceleration) {
    int igrid;
    float rho, gx, gy, gz;
    float vx, vy, vz, vx_old, vy_old, vz_old;
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  rho = BaryonField[DensNum][igrid];
	  
	  gx = ConstantAcceleration[0];
	  gy = ConstantAcceleration[1];
	  gz = ConstantAcceleration[2];
	  vx = BaryonField[Vel1Num][igrid];
	  vy = BaryonField[Vel2Num][igrid];
	  vz = BaryonField[Vel3Num][igrid];
	  
	  dU[iS1][n] += dtFixed*gx*rho;
	  dU[iS2][n] += dtFixed*gy*rho;
	  dU[iS3][n] += dtFixed*gz*rho;
	  dU[iEtot][n] += dtFixed*rho*(gx*vx + gy*vy + gz*vz);

	if (i==3 && j==3 && k==4 && GridLeftEdge[0]==0.0 && GridLeftEdge[1]==1.0)
	  printf("StermStart4 old %"GSYM" \n", dU[iS2][n])  ;
	}
      }
    }
  }

   
  if ((SelfGravity && GridRank == 3) || ExternalGravity) {
    int igrid;
    float rho, gx, gy, gz;
    float vx, vy, vz, vx_old, vy_old, vz_old;
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  rho = BaryonField[DensNum][igrid];
	  
	  gx = AccelerationField[0][igrid];
	  gy = AccelerationField[1][igrid];
	  gz = AccelerationField[2][igrid];
	  vx = BaryonField[Vel1Num][igrid];
	  vy = BaryonField[Vel2Num][igrid];
	  vz = BaryonField[Vel3Num][igrid];
	  
	  dU[iS1  ][n] += dtFixed*gx*rho;
	  dU[iS2  ][n] += dtFixed*gy*rho;
	  dU[iS3  ][n] += dtFixed*gz*rho;
	  dU[iEtot][n] += dtFixed*rho*(gx*vx + gy*vy + gz*vz);

	
	}
      }
    }
  }

  if ((ComovingCoordinates == 1)) { // add cosmological expansion terms here

    int igrid;
    float rho, coef=0.;
    FLOAT a, dadt;
    int n = 0;
    CosmologyComputeExpansionFactor(0.5*(Time+OldTime), &a, &dadt);
    coef = -0.5*dadt/a;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  rho = BaryonField[DensNum][igrid];
	  
	  
	  dU[iBx  ][n] += dtFixed*coef*BaryonField[B1Num][n];
	  dU[iBy  ][n] += dtFixed*coef*BaryonField[B2Num][n];
	  dU[iBz  ][n] += dtFixed*coef*BaryonField[B3Num][n];
	  dU[iEtot][n] -= dtFixed*coef*(BaryonField[B1Num][n]*BaryonField[B1Num][n]+
					BaryonField[B2Num][n]*BaryonField[B2Num][n]+
					BaryonField[B3Num][n]*BaryonField[B3Num][n]);
	  dU[iPhi][n] += 0.0; // Add correct Phi term here .....


	}
      }
    }
  }

  /* Apply external driving force */

  if (UseDrivingField) {

    float rhou, lenu, tempu, tu, velu;
    GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);

    int Drive1Num, Drive2Num, Drive3Num;
    if (IdentifyDrivingFields(Drive1Num, Drive2Num, Drive3Num) == FAIL) {
      printf("grid::SourceTerms: canot identify driving fields.\n");
      return FAIL;
    }
    int igrid;
    float drivex, drivey, drivez, vx, vy, vz, rho;
    int n = 0;
    FLOAT x, y, z, r;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];

	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	  r = sqrt(pow(x-0.5,2) + pow(y-0.5,2) + pow(z-0.5,2));
	

	  rho    = BaryonField[DensNum     ][igrid];
	  drivex = BaryonField[Drive1Num][igrid];
	  drivey = BaryonField[Drive2Num][igrid];
	  drivez = BaryonField[Drive3Num][igrid];
	  vx     = BaryonField[Vel1Num      ][igrid];
	  vy     = BaryonField[Vel2Num      ][igrid];
	  vz     = BaryonField[Vel3Num      ][igrid];

	  float eint = BaryonField[TENum][igrid] - 0.5*sqrt(vx*vx + vy*vy + vz*vz)
	    -0.5*sqrt(pow(BaryonField[B1Num][igrid],2)+pow(BaryonField[B2Num][igrid],2)
		      +pow(BaryonField[B3Num][igrid],2))/rho;

	  float T = eint*Mu*(Gamma-1.0)*tempu;

	  FLOAT R = sqrt(pow(x-0.5,2) + pow(y-0.5,2));

	  //if (T > 90.0) {
	  //if (z < 0.1 || z > 0.9) {
	  //if (r > 0.3) {
	  //if (R > 0.45) {
	    dU[iS1  ][n] += dtFixed*rho*drivex*DrivingEfficiency;
	    dU[iS2  ][n] += dtFixed*rho*drivey*DrivingEfficiency;
	    dU[iS3  ][n] += dtFixed*rho*drivez*DrivingEfficiency;
	    dU[iEtot][n] += dtFixed*rho*(drivex*vx + drivey*vy + drivez*vz)*DrivingEfficiency;


	    //}
	}
      }
    }
  }


  
  /* Add centrifugal force for the shearing box */

  if ((ProblemType == 35 || ProblemType == 36 ||ProblemType == 37) && ShearingBoxProblemType !=0) {


 int igrid;
    float rho, gx, gy, gz;
    FLOAT xPos[3];
    float vels[3]; 
    int n = 0;

    int iden=FindField(Density, FieldType, NumberOfBaryonFields);
    int ivx=FindField(Velocity1, FieldType, NumberOfBaryonFields);
    int ivy=FindField(Velocity2, FieldType, NumberOfBaryonFields);
    int ivz;
    if (GridRank==3)  ivz=FindField(Velocity3, FieldType, NumberOfBaryonFields);
 
    int indexNumbers[3]={iS1,iS2,iS3};

    float A[3]={0,0,0};//Omega
    A[ShearingOtherDirection]=AngularVelocity;
    
    float lengthx=DomainRightEdge[0]-DomainLeftEdge[0]; 
    float lengthy=DomainRightEdge[1]-DomainLeftEdge[1];
    float lengthz;
    if (GridRank==3) lengthz=DomainRightEdge[2]-DomainLeftEdge[2];
    else lengthz-0.0;
    

    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {

	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  rho = BaryonField[iden][igrid];
	  xPos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i]-lengthx/2.0;
	  xPos[1] = CellLeftEdge[1][i] + 0.5*CellWidth[1][i]-lengthy/2.0;
	  if (GridRank==3) xPos[2] = CellLeftEdge[2][i] + 0.5*CellWidth[2][i]-lengthz/2.0;
	  else xPos[2]=0;
	  
	  vels[0] = BaryonField[ivx][igrid];
	  vels[1] = BaryonField[ivy][igrid];
	  if (GridRank==3) vels[2] = BaryonField[ivz][igrid];
	  else vels[2]=0;

	  //Omega cross V

	  dU[indexNumbers[0]][n] -= dtFixed*2.0*rho*(A[1]*vels[2]-A[2]*vels[1]);
	  dU[indexNumbers[1]][n] -= dtFixed*2.0*rho*(A[2]*vels[0]-A[0]*vels[2]);
	  if (GridRank==3) dU[indexNumbers[2]][n] -= dtFixed*2.0*rho*(A[0]*vels[1]-A[1]*vels[0]);
	

	  dU[indexNumbers[ShearingBoundaryDirection]][n] += dtFixed*2.0*rho*VelocityGradient*AngularVelocity*AngularVelocity*xPos[ShearingBoundaryDirection];
	  
	  
	  dU[iEtot][n] +=  dtFixed*2.0*rho*VelocityGradient*AngularVelocity*AngularVelocity*xPos[ShearingBoundaryDirection]*vels[ShearingBoundaryDirection];
	



	  
 	}
      }
    }
  }


  return SUCCESS;
}
