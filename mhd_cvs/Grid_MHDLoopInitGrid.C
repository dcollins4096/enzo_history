
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
//start


int grid::MHDLoopInitGrid(float LoopDensity,float Pressure, float Vx, float Vy, float Vz, float B0, float R0, 
			  float Center[], int CurrentAxis){ 

  //Every processor needs to know this for every grid,
  //WHETHER OR NOT IT HAS THE DATA.

  NumberOfBaryonFields = 5;
#ifdef HAOXU
        if (DualEnergyFormalism){
        NumberOfBaryonFields = 6;
        FieldType[5] = InternalEnergy;
        }
#endif

  fprintf(stderr,"GridDim %d %d %d\n",GridDimension[0],GridDimension[1],GridDimension[2]);
  int field = 0;

#ifdef ATHENA
    NumberOfBaryonFields = 5;
    if(DualEnergyFormalism) NumberOfBaryonFields++;
    if(EquationOfState == 1 ) NumberOfBaryonFields--;

    FieldType[field++] = Density;
    if( EquationOfState == 0 ) FieldType[field++] = TotalEnergy;
    FieldType[field++] = Velocity1;
    FieldType[field++] = Velocity2;
    FieldType[field++] = Velocity3;
    if(DualEnergyFormalism) FieldType[field++] = InternalEnergy;

#else //Athena

  FieldType[0] = Density;
  FieldType[1] = TotalEnergy;
  FieldType[2] = Velocity1;
  FieldType[3] = Velocity2;
  FieldType[4] = Velocity3;
#endif //ATHENA


  //Variable names.
  //I know, I know, I just set Field Type.  This is for generality and consistency with
  //other 
  int Eeng, Eden, Ev[3], Egas; 
  if (this->IdentifyPhysicalQuantities(Eden, Egas, Ev[0], Ev[1], 
				       Ev[2], Eeng) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }


  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"======== Loop ===================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");



  float Pi = 3.14159265, One=1.0;
  float R,X, Y, GasEnergy=Pressure/(Gamma-1), LoopTotalEnergy=0; 
  int index, size=1, i,j,k, Three=3,TENum=1; 
  float Scale[3];

  this->AllocateGrids();  

  for(i=0;i<GridRank;i++){
    size*=GridDimension[i];
    Scale[i]=(GridRightEdge[i]-GridLeftEdge[i])/(GridDimension[i]-2*DEFAULT_GHOST_ZONES);

  }
  for( i=GridRank; i<3; i++){
    Scale[i] = 0;
    Center[i] = 0;
  }

  fprintf(stderr,"Density %f Pressure %f (vx,vy) (%f,%f) B0 %f \n", LoopDensity,Pressure,Vx,Vy,B0);
  fprintf(stderr,"Scale: %f %f %f\n", Scale[0],Scale[1],Scale[2]);
  fprintf(stderr,"Center; %f %f %f\n", Center[0], Center[1], Center[2]);

  //Vector Potential. 
  //Due to the similarity in centering, and lack of foresight in naming,
  //I'm using the Electric Field as a Vector Potential to initialize the Magnetic Field.  

  for(field=0;field<3;field++)
    for( k=0;k<ElectricDims[field][2];k++)
      for( j=0;j<ElectricDims[field][1];j++)
	for( i=0;i<ElectricDims[field][0];i++){

	  index=i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k);

	  if( field== CurrentAxis ){

	    switch( CurrentAxis ){
	    case 0:
	      X=(j-GridStartIndex[1])*Scale[1]-(Center[1]-DomainLeftEdge[1]);
	      Y=(k-GridStartIndex[2])*Scale[2]-(Center[2]-DomainLeftEdge[2]);
	      break;
	    case 1:
	      X=(k-GridStartIndex[2])*Scale[2]-(Center[2]-DomainLeftEdge[2]);
	      Y=(i-GridStartIndex[0])*Scale[0]-(Center[0]-DomainLeftEdge[0]);
	      break;
	    case 2:
	      X=(i-GridStartIndex[0])*Scale[0]-(Center[0]-DomainLeftEdge[0]);
	      Y=(j-GridStartIndex[1])*Scale[1]-(Center[1]-DomainLeftEdge[1]);
	      break;
	    default:
	      fprintf(stderr," Hey, Jerk:  MHDLoopCurrentAxis = %d not defined. Try again.\n",
		      CurrentAxis);
	      break;
	    }//Current Axis switch.
	    R=sqrt(X*X+Y*Y);
	    ElectricField[field][index]=(R<R0)? B0*(R0-R):0.0;

	  }else{
	    ElectricField[field][index]=0.0;

	  }

	}
  

  //
  //Curl is B=Curl(A)
  //


  //the last argument indicates that this isn't a time update.  See the source for details.
  if( this->MHD_Curl(GridStartIndex, GridEndIndex, 0) == FAIL )
    {fprintf(stderr," error occored in MHD_Curl\n"); return FAIL;}


  if( this->CenterMagneticField() == FAIL ) 
    {fprintf(stderr," error with CenterMagneticField , second call\n");return FAIL;}

  //Set the rest of the fields.
  for(k=0;k<GridDimension[2];k++)
    for(j=0;j<GridDimension[1];j++)
      for(i=0;i<GridDimension[0];i++){

	index=index0(i,j,k);
	X=(i-GridStartIndex[0])*Scale[0];
	Y=(j-GridStartIndex[1])*Scale[1];
	
	LoopTotalEnergy=GasEnergy + 0.5*LoopDensity*(Vx*Vx+Vy*Vy + Vz*Vz)
	  +0.5*(CenteredB[0][index]*CenteredB[0][index]+
		CenteredB[1][index]*CenteredB[1][index]+
		CenteredB[2][index]*CenteredB[2][index]);

	BaryonField[0][index]=Density;
	BaryonField[1][index]=TotalEnergy;
#ifdef HAOXU
        if (DualEnergyFormalism)
        BaryonField[5][index]=GasEnergy/Density;
#endif
	BaryonField[2][index]=Vx;
	BaryonField[3][index]=Vy;
	BaryonField[4][index]=0.0;

#ifdef ATHENA
	BaryonField[Eden][index]=LoopDensity;
	if( EquationOfState == 0 ) BaryonField[Eeng][index]=LoopTotalEnergy;
	BaryonField[ Ev[0] ][index]=Vx;
	BaryonField[ Ev[1] ][index]=Vy;
	BaryonField[ Ev[2] ][index]=Vz;
#endif 
      }
  

  return SUCCESS;

}
