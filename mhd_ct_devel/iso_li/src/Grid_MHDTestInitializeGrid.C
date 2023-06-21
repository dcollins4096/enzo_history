#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
//start

int grid::MHDTestInitializeGrid(float Density0, float Density1,
				 float Energy0,  float Energy1,
				 float Velocity10, float Velocity11,
				 float Velocity20, float Velocity21,
				 float Velocity30, float Velocity31,
				 float B10,float B20,float B30, float Radius)

{

 if (ProcessorNumber != MyProcessorNumber)
   return SUCCESS;

  int nxt = GridDimension[0], nyt = GridDimension[1], nzt = GridDimension[2];
  int nb = DEFAULT_GHOST_ZONES;
  int nx = nxt-2*nb,  ny = nyt-2*nb, nz = nzt-2*nb;
  int i,j,k, field, index, size = 1, dim;

  int LongAxis = GridDimension[0];
  LongAxis = LongAxis > GridDimension[1] ? LongAxis : GridDimension[1];
  LongAxis = LongAxis > GridDimension[2] ? LongAxis : GridDimension[2];
  LongAxis = LongAxis == GridDimension[0] ? 0 : LongAxis;
  LongAxis = LongAxis == GridDimension[1] ? 1 : LongAxis;
  LongAxis = LongAxis == GridDimension[2] ? 2 : LongAxis;


  //LongAxis = 0;
  //Initialize Magnetic and Electric Fields
  
  /* later change this to a loop and use Magnetic and Electric Size */
  MagneticField[0] = new float[(nxt+1)*(nyt)*(nzt)];
  MagneticField[1] = new float[(nxt)*(nyt+1)*(nzt)];
  MagneticField[2] = new float[(nxt)*(nyt)*(nzt+1)];
  ElectricField[0] = new float[(nxt)*(nyt+1)*(nzt+1)];
  ElectricField[1] = new float[(nxt+1)*(nyt)*(nzt+1)];
  ElectricField[2] = new float[(nxt+1)*(nyt+1)*(nzt)];

  size = 1;
  for(i=0;i<3;i++)
    size *= GridDimension[i];

  for(i=0;i<3;i++)
    CenteredB[i] = new float[size];
  
  for(i=0;i<size;i++){
    CenteredB[0][i] = 0.0;
    CenteredB[1][i] = 0.0;
    CenteredB[2][i] = 0.0;
  }
      
  //Initialize Barryon Fields

  NumberOfBaryonFields = 5;
  FieldType[0] = Density;
  FieldType[1] = TotalEnergy;
  FieldType[2] = Velocity1;
  FieldType[3] = Velocity2;
  FieldType[4] = Velocity3;

  size = 1;
  for(i = 0; i < GridRank; i++) size *= GridDimension[i];

  for(i = 0; i < NumberOfBaryonFields; i++) BaryonField[i] = new float[size];

  B10 = 0.0;
  B20 = 0.0;
  B30 = 0.0;
  

  for(i = 0; i<size; i++)
    {
      BaryonField[0][i] = 0.0;
      BaryonField[1][i] = 0.0;
      BaryonField[2][i] = 0.0;
      BaryonField[3][i] = 0.0;
      BaryonField[4][i] = 0.0;
    }

  for(dim = 0; dim<3; dim++) 
    for(i=0; i<ElectricSize[dim]; i++) ElectricField[dim][i] = 0.0;

  for(dim=0;dim<3;dim++)
    for(i=0; i<MagneticSize[dim]; i++) MagneticField[dim][i] = 0.0;
  
  //
  // Set Magnetic Field to test values.
  // Note: these are totally unphysical quantities, to be used only for 
  // debugging of data handling routines (IO, communication, etc.)
  //

  if( 1 == 0 ){
    
    for(field=0;field<3;field++)
      for(k=MHDStartIndex[field][2]; k<=MHDEndIndex[field][2]; k++)
	for(j=MHDStartIndex[field][1]; j<=MHDEndIndex[field][1]; j++)
	  for(i=MHDStartIndex[field][0]; i<=MHDEndIndex[field][0];i++){
	    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
	    if(LongAxis == 0)
	      MagneticField[field][index] = 10*(k+10*(j+10*i));//field + 6;
	    else if(LongAxis == 1)
	      MagneticField[field][index] = 10*(i+10*(k+10*j));//field + 6;
	    else if(LongAxis == 2)
	      MagneticField[field][index] = 10*(j+10*(i+10*k));//field + 6;
	  }
    /*
      for(field=0;field<NumberOfBaryonFields;field++)
      for(k=GridStartIndex[2];k<=GridEndIndex[2]; k++)
      for(j=GridStartIndex[1];j<=GridEndIndex[1];j++)
      for(i=GridStartIndex[0]; i<=GridEndIndex[0];i++)
      {
      index = i+GridDimension[0]*(j+GridDimension[1]*k);
      BaryonField[field][index] = 10*(k+10*(j+10*i)) + field+1; 
      }
    */
    for(field=0;field<NumberOfBaryonFields;field++)
      for(k=0; k<GridDimension[2]; k++)
	for(j=0;j<GridDimension[1]; j++)
	  for(i=0;i<GridDimension[0]; i++)
	    {
	      index = i+GridDimension[0]*(j+GridDimension[1]*k);
	      if(LongAxis == 0)
		BaryonField[field][index] = 10*(k+10*(j+10*i)) + field+1; 
	      if(LongAxis == 1)
		BaryonField[field][index] = 10*(i+10*(k+10*j)) + field+1; 
	      if(LongAxis == 2)
		BaryonField[field][index] = 10*(j+10*(i+10*k)) + field+1; 
	    }
    for(field=0;field<3;field++)
      for(k=GridStartIndex[2];k<=GridEndIndex[2]; k++)
	for(j=GridStartIndex[1];j<=GridEndIndex[1];j++)
	  for(i=GridStartIndex[0]; i<=GridEndIndex[0];i++)
	    {
	      index = i+GridDimension[0]*(j+GridDimension[1]*k);
	      CenteredB[field][index] = 10*(k+10*(j+10*i)) + field+7; 
	    }
  }// End Kludgie Skip.

  /*
    
    for(field=0;field<3;field++)
    for(k=0;k<MagneticDims[field][2]; k++)
    for(j=0;j<MagneticDims[field][1];j++)
    for(i=0;i<MagneticDims[field][0]; i++){
    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
    MagneticField[field][index] = 10*(k+10*(j+10*i)) + field;
    //MagneticField[field][index] = 0.0;
    }
  */
  
  //
  // Set Electric Field
  //

  field = 0;
  for(k=MHDeStartIndex[field][2];k<=MHDeEndIndex[field][2];k++)
    for(j=MHDeStartIndex[field][1];j<=MHDeEndIndex[field][1];j++)
      for(i=MHDeStartIndex[field][0];i<=MHDeEndIndex[field][0];i++){
	index = i+ElectricDims[field][0]*(j+ElectricDims[field][1]*(k));
	ElectricField[field][index] = 0.0;
      }

  field = 1;
  for(k=MHDeStartIndex[field][2];k<=MHDeEndIndex[field][2];k++)
    for(j=MHDeStartIndex[field][1];j<=MHDeEndIndex[field][1];j++)
      for(i=MHDeStartIndex[field][0];i<=MHDeEndIndex[field][0];i++){
	index = i+ElectricDims[field][0]*(j+ElectricDims[field][1]*(k));
	ElectricField[field][index] = 0.0;
      }
  
  field = 2;
  for(k=MHDeStartIndex[field][2];k<=MHDeEndIndex[field][2];k++)
    for(j=MHDeStartIndex[field][1];j<=MHDeEndIndex[field][1];j++)
      for(i=MHDeStartIndex[field][0];i<=MHDeEndIndex[field][0];i++){
	index = i+ElectricDims[field][0]*(j+ElectricDims[field][1]*(k));
	ElectricField[field][index] = 0.0;
      }
  
  return SUCCESS;
}

