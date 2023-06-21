#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

float nrg(float Pressure, float Gamma, float density, float v1, float v2, float v3, float b1, float b2, float b3)
{
  
  return (Pressure)/((Gamma-1))
    + 0.5*(v1*v1+v2*v2+v3*v3)*density
    + 0.5*(b1*b1+b2*b2+b3*b3);
  
}


int grid::MHDShockInitializeGrid(int ShockTubeDirection)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;


  
  //int nxt = GridDimension[0], nyt = GridDimension[1], nzt = GridDimension[2];
  //t nb = DEFAULT_GHOST_ZONES;
  // int nx = nxt-2*nb,  ny = nyt-2*nb, nz = nzt-2*nb;
  int i,j,k, field, index, size = 1, dim;
  float pi = 3.14159265;

  //The spatial extents of the initialization.
  int InitStartIndex[3]={0,0,0};
  int InitEndIndex[3] = {GridDimension[0]-1,GridDimension[1]-1,GridDimension[2]-1};
  
  //Initial conditions taken from Ryu et all 1998
  float Density0 = 1.0, 
    Density1 = 1.0,
    Energy0 = 20.0,  
    Energy1 = 1.0, 
    Velocity10 = 10.0,
    Velocity11 = -10.0,
    Velocity20 = 0.0,
    Velocity21 = 0.0,
    Velocity30 = 0.0,
    Velocity31 = 0.0,
    B0[3] ={0.000001,0.0,1.41},
    B1[3] ={0.000001,0.0,1.41},
    Pressure0 = 460.894,
    Pressure1 = 46.0950;

  if(ShockTubeDirection == 1){
  //Numbers used for Kim et al Isothermal shock
    Density0 = 5.99924; 
    Density1 = 5.99242;
    
    Velocity10 = 19.5975;//1.2; 
    Velocity11 = -6.19633;
    Velocity20 = 0.0;//originally 0.01
    Velocity21 = 0.0;
    Velocity30 = 0.0; //originally 0.5
    Velocity31 = 0.0;
    B0[0] = 0.0;0.846296855141 ;
    B1[0] = 0.0;0.846296855141 ;
    B0[1] = 0.0;1.41049475857;
    B1[1] = 0.0;0.564197903428 ;
    B0[2] = 0.0;
    B1[2] = 0.0;

  }
  
  if(ShockTubeDirection == 2){
    Density0 = 1.08; 
    Density1 = 1.0;
    Energy0 = 0.95;  
    Energy1 = 1.0;
#ifdef HAOXU
    Pressure0 = 0.95;
    Pressure1 =1.0;
#endif 
    Velocity10 = 0.5;
    Velocity11 = 0.0;
    Velocity20 = 1.2;
    Velocity21 = 0.0;
    Velocity30 = 0.01;
    Velocity31 = 0.0;
    B0[0] = 2.0/sqrt(4.0*pi);
    B1[0] = 2.0/sqrt(4.0*pi);
    B0[1] = 2.0/sqrt(4.0*pi);
    B1[1] = 2.0/sqrt(4.0*pi);
    B0[2] = 3.6/sqrt(4.0*pi);
    B1[2] = 4.0/sqrt(4.0*pi);
    
  }
  if(ShockTubeDirection == 3){
    Density0 = 1.08; 
    Density1 = 1.0;
    Energy0 = 0.95;  
    Energy1 = 1.0; 
    Velocity10 = 0.01;
    Velocity11 = 0.0;
    Velocity20 = 0.5;
    Velocity21 = 0.0;
    Velocity30 = 1.2;
    Velocity31 = 0.0;
    B0[0] = 3.6/sqrt(4.0*pi);
    B1[0] = 4.0/sqrt(4.0*pi);
    B0[1] = 2.0/sqrt(4.0*pi);
    B1[1] = 2.0/sqrt(4.0*pi);
    B0[2] = 2.0/sqrt(4.0*pi);
    B1[2] = 2.0/sqrt(4.0*pi);
    
  }
  //4 is a diagonal shock tube along the 1,1,1 line.
  if(ShockTubeDirection == 4){
    
    Density0 = 1.08; 
    Density1 = 1.0;
    Energy0 = 0.95;  
    Energy1 = 1.0; 
    Velocity10 = 1.2; 
    Velocity11 = 0.0;
    Velocity20 = 0.01;
    Velocity21 = 0.0;
    Velocity30 = 0.5;
    Velocity31 = 0.0;
    B0[0] = 2.0/sqrt(4.0*pi);
    B1[0] = 2.0/sqrt(4.0*pi);
    B0[1] = 3.6/sqrt(4.0*pi);
    B1[1] = 4.0/sqrt(4.0*pi);
    B0[2] = 2.0/sqrt(4.0*pi);
    B1[2] = 2.0/sqrt(4.0*pi);
  }
  
  
  //Initialize Magnetic and Electric Fields
  
  for(field = 0;field<3;field++){
    MagneticField[field] = new float[MagneticSize[field]];
    ElectricField[field] = new float[ElectricSize[field]];
  }
  
  //Initialize Barryon Fields

  NumberOfBaryonFields = 5;
#ifdef HAOXU
if(DualEnergyFormalism == 1){
  NumberOfBaryonFields = 6;
  FieldType[5] = InternalEnergy;
}
#endif //HAOXU
  FieldType[0] = Density;
  FieldType[1] = TotalEnergy;
  FieldType[2] = Velocity1;
  FieldType[3] = Velocity2;
  FieldType[4] = Velocity3;

  size = 1;
  for(i = 0; i < GridRank; i++) size *= GridDimension[i];

  for(i = 0; i < NumberOfBaryonFields; i++) BaryonField[i] = new float[size];

  for(field=0;field< NumberOfBaryonFields; field++)
    for(i=0;i<size;i++) 
      BaryonField[field][i] = 0.0;

  for(i=0;i<3;i++){
    CenteredB[i] = new float[size];
    Current[i] = new float[size];
  }
  DivB = new float[size];
  

  
  for(i=0;i<size;i++){
    CenteredB[0][i] = 0.0;
    CenteredB[1][i] = 0.0;
    CenteredB[2][i] = 0.0;
    DivB[i] = 0.0;
    Current[1][i] = 0.0;
    Current[2][i] = 0.0;
    Current[0][i] = 0.0;
  }

  for(dim = 0; dim<3; dim++) {
    for(i=0; i<ElectricSize[dim]; i++) ElectricField[dim][i] = 0.0;
    for(i=0; i<MagneticSize[dim]; i++) MagneticField[dim][i] = 0.0;
  }



  //
  // Set values for shock 
  //
  fprintf(stderr,"(GridDimension[0]/2 - 1) or not %d %d\n",
	  (GridDimension[0]/2 - 1),(GridDimension[0]/2));
  if( ShockTubeDirection == 1){
    
    
    fprintf(stderr, "Energy Kludge!!!! \n");
    for(k=InitStartIndex[2];k<=InitEndIndex[2]; k++)
      for(j=InitStartIndex[1];j<=InitEndIndex[1];j++){
	
	for(i=InitStartIndex[0];i<(GridDimension[0]/2+1);i++){
	  index = i+GridDimension[0]*(j+GridDimension[1]*k);
	  //oi(i,j,k);
	  //cout << index << " ";
	  BaryonField[0][index] = Density0; 
	  BaryonField[2][index] = Velocity10;
	  BaryonField[3][index] = Velocity20;
	  BaryonField[4][index] = Velocity30;
	  CenteredB[0][index] = B0[0];
	  CenteredB[1][index] = B0[1];
	  CenteredB[2][index] = B0[2];
	  //	  BaryonField[1][index] = (Pressure0)/((Gamma-1)*Density0)
	  BaryonField[1][index] = nrg(Pressure0, Gamma, Density0, 
				      Velocity10, Velocity20, Velocity30, 
				      B0[0],B0[1], B0[2]);
#ifdef HAOXU
      if(DualEnergyFormalism == 1){
         BaryonField[5][index] = (Pressure0)/((Gamma-1)*Density0);
       }
#endif //HAOXU

	}
	//cout << "\n--------------\n";
	for(i=GridDimension[0]/2+1;i<=InitEndIndex[0];i++){
	  index = i+GridDimension[0]*(j+GridDimension[1]*k);
	  //oi(i,j,k);
	  //cout << index << " ";
	  BaryonField[0][index] = Density1; 
	  BaryonField[1][index] = Energy1;
	  BaryonField[2][index] = Velocity11;
	  BaryonField[3][index] = Velocity21;
	  BaryonField[4][index] = Velocity31;
	  CenteredB[0][index] = B1[0];
	  CenteredB[1][index] = B1[1];
	  CenteredB[2][index] = B1[2];
	  BaryonField[1][index] = nrg(Pressure1, Gamma, Density1, 
				      Velocity11, Velocity21, Velocity31, 
				      B1[0],B1[1], B1[2]);
#ifdef HAOXU
      if(DualEnergyFormalism == 1){
         BaryonField[5][index] = (Pressure1)/((Gamma-1)*Density1);
       }
#endif //HAOXU

	}
	//cout << "\n++++++++++++++\n";
      }
    
    
    
    for(field=0;field<3;field++)

      for(k=0; k<MagneticDims[field][2]; k++)
	for(j=0; j<MagneticDims[field][1]; j++){
	  
	  for(i=0; i<MagneticDims[field][0]/2+1;i++){
	    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
	    MagneticField[field][index] = B0[field];
	  }
	  
	  for(i=MagneticDims[field][0]/2+1; i<MagneticDims[field][0];i++){
	    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
	    MagneticField[field][index] = B1[field];
	  }
	  
	}
    
    
  }     
  else if( ShockTubeDirection == 2 ){
    
    for(k=InitStartIndex[2];k<=InitEndIndex[2]; k++)
      for(i=InitStartIndex[0];i<=InitEndIndex[0];i++){
	
	for(j=InitStartIndex[1];j<GridDimension[1]/2;j++){
	  
	  index = i+GridDimension[0]*(j+GridDimension[1]*k);
	  //oi(i,j,k);
	  //cout << index << " ";
	  BaryonField[0][index] = Density0; 
	  BaryonField[2][index] = Velocity10;
	  BaryonField[3][index] = Velocity20;
	  BaryonField[4][index] = Velocity30;
	  CenteredB[0][index] = B0[0];
	  CenteredB[1][index] = B0[1];
	  CenteredB[2][index] = B0[2];
	  BaryonField[1][index] = nrg(Pressure0, Gamma, Density0, 
				      Velocity10, Velocity20, Velocity30, 
				      B0[0],B0[1], B0[2]);
#ifdef HAOXU
      if(DualEnergyFormalism == 1){
         BaryonField[5][index] = (Pressure0)/((Gamma-1)*Density0);
       }
#endif //HAOXU
	  
	}
	//cout << "\n--------------\n";
	for(j=GridDimension[1]/2;j<=InitEndIndex[1];j++){
	  
	  index = i+GridDimension[0]*(j+GridDimension[1]*k);
	  //oi(i,j,k);
	  //cout << index << " ";
	  BaryonField[0][index] = Density1; 
	  BaryonField[1][index] = Energy1;
	  BaryonField[2][index] = Velocity11;
	  BaryonField[3][index] = Velocity21;
	  BaryonField[4][index] = Velocity31;
	  CenteredB[0][index] = B1[0];
	  CenteredB[1][index] = B1[1];
	  CenteredB[2][index] = B1[2];
	  BaryonField[1][index] = nrg(Pressure1, Gamma, Density1, 
				      Velocity11, Velocity21, Velocity31, 
				      B1[0],B1[1], B1[2]);
#ifdef HAOXU
      if(DualEnergyFormalism == 1){
         BaryonField[5][index] = (Pressure1)/((Gamma-1)*Density1);
       }
#endif //HAOXU

	}
	//cout << "\n++++++++++++++\n";
	
      }
    //cout << "\n++++++++++++++\n";
    
    
    for(field=0;field<3;field++)
      for(k=MHDStartIndex[field][2]; k<=MHDEndIndex[field][2]; k++)
	for(i=MHDStartIndex[field][0];i<=MHDEndIndex[field][0];i++){
	  
	  for(j=MHDStartIndex[field][1];j<MagneticDims[field][1]/2;j++){
	    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
	    MagneticField[field][index] = B0[field];
	  }
	  
	  for(j=MagneticDims[field][1]/2;j<=MHDEndIndex[field][1];j++){
	    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
	    MagneticField[field][index] = B1[field];
	  }
	  
	}
    
  } //end Shock Tube == 2
  else if(ShockTubeDirection == 3 ) {
    for(k=InitStartIndex[2];k<GridDimension[2]/2; k++)
      for(j=InitStartIndex[1];j<=InitEndIndex[1];j++){
	for(i=InitStartIndex[0];i<=InitEndIndex[0];i++){
	  
	  index = i+GridDimension[0]*(j+GridDimension[1]*k);
	  //oi(i,j,k);
	  //cout << index << " ";
	  BaryonField[0][index] = Density0; 
	  BaryonField[2][index] = Velocity10;
	  BaryonField[3][index] = Velocity20;
	  BaryonField[4][index] = Velocity30;
	  CenteredB[0][index] = B0[0];
	  CenteredB[1][index] = B0[1];
	  CenteredB[2][index] = B0[2];
	  BaryonField[1][index] = nrg(Pressure0, Gamma, Density0, 
				      Velocity10, Velocity20, Velocity30, 
				      B0[0],B0[1], B0[2]);
#ifdef HAOXU
      if(DualEnergyFormalism == 1){
         BaryonField[5][index] = (Pressure0)/((Gamma-1)*Density0);
       }
#endif //HAOXU
	  
	}
      }
    for(k=GridDimension[2]/2;k<=InitEndIndex[2];k++)
      for(j=InitStartIndex[1];j<=InitEndIndex[1];j++)
	for(i=InitStartIndex[0];i<=InitEndIndex[0];i++){
	  
	  index = i+GridDimension[0]*(j+GridDimension[1]*k);
	  //oi(i,j,k);
	  //cout << index << " ";
	  BaryonField[0][index] = Density1; 
	  BaryonField[1][index] = Energy1;
	  BaryonField[2][index] = Velocity11;
	  BaryonField[3][index] = Velocity21;
	  BaryonField[4][index] = Velocity31;
	  CenteredB[0][index] = B1[0];
	  CenteredB[1][index] = B1[1];
	  CenteredB[2][index] = B1[2];
	  BaryonField[1][index] = nrg(Pressure1, Gamma, Density1, 
				      Velocity11, Velocity21, Velocity31, 
				      B1[0], B1[1], B1[2]);
#ifdef HAOXU
      if(DualEnergyFormalism == 1){
         BaryonField[5][index] = (Pressure1)/((Gamma-1)*Density1);
       }
#endif //HAOXU
	}
    
    
    for(field=0;field<3;field++)
      for(k=MHDStartIndex[field][2];k<MagneticDims[field][2]/2;k++)
	for(j=MHDStartIndex[field][1];j<=MHDEndIndex[field][1];j++)
	  for(i=MHDStartIndex[field][0];i<=MHDEndIndex[field][0];i++){
	    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
	    MagneticField[field][index] = B0[field];
	  }
    for(field=0;field<3;field++)
      for(k=MagneticDims[field][2]/2;k<=MHDEndIndex[field][2];k++)
	for(j=MHDStartIndex[field][1];j<=MHDEndIndex[field][1];j++)
	  for(i=MHDStartIndex[field][0];i<=MHDEndIndex[field][0];i++){
	    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
	    MagneticField[field][index] = B1[field];
	  }
    
    
  } //end Shock Tube == 3
  if( 0 == 1 ){
    for(field=0;field<NumberOfBaryonFields;field++)
      for(k=0; k<GridDimension[2]; k++)
	for(j=0;j<GridDimension[1]; j++)
	  for(i=0;i<GridDimension[0]; i++)
	    {
	      index = i+GridDimension[0]*(j+GridDimension[1]*k);
	      BaryonField[field][index] = 1000*(k+1000*(j+1000*i)) + field+1; 
	    }
  }

  return SUCCESS;
}

