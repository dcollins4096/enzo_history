
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

#ifdef PPML
#include "PPML.h"
#endif //PPML 

#ifndef MHD
//This is an old version, easier to install than the proper grid method.
extern "C" void FORTRAN_NAME(curl_of_e)(float *bx, float *by, float *bz,
                                        float *ex, float *ey, float *ez,
                                        float *dx, float *dy, float *dz,
                                        int *idim, int *jdim, int *kdim,
                                        int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
                                        float *dt, int *method);
extern "C" void FORTRAN_NAME(center_magnetic_field)(float *bxf, float *byf, float *bzf,
                                    float *bxc, float *byc, float *bzc,
                                    float *energy,
                                    FLOAT *dx, FLOAT *dy, FLOAT *dz,
                                    int *idim, int *jdim, int *kdim,
                                    int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
                                    int *method);

#endif //!MHD


int grid::MHDOrszagTangInitGrid(float DensityIn,float Pressure, float V0, float B0 ){ 

  //Every processor needs to know this for every grid,
  //WHETHER OR NOT IT HAS THE DATA.
  

  NumberOfBaryonFields = 0;
  fprintf(stderr,"GridDim %d %d %d\n",GridDimension[0],GridDimension[1],GridDimension[2]);
  FieldType[NumberOfBaryonFields++] = Density;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  if( EquationOfState == 0 )
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;

  if( MHD_Used == TRUE ){
    FieldType[NumberOfBaryonFields++] = Magnetic1;
    FieldType[NumberOfBaryonFields++] = Magnetic2;
    FieldType[NumberOfBaryonFields++] = Magnetic3;
  }

#ifdef PPML
  PPML_InterfacePointerBundle  *Face;
  if( this->PPML_InitInterfaceTypesAndLabels() == FAIL ){
    fprintf(stderr," Grid_DiscontInit...: Failure in PPML_InitInterface...\n");
    return FAIL;
  }
#endif //PPML

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"======== tang ===================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");

  float Pi = 3.14159265, One=1.0;
  float X, Y, Vx, Vy, GasEnergy=Pressure/(Gamma-1), TotalEnergy=0; 
  int index, size=1, i,j,k, field;
  float Scale[3];
  this->AllocateGrids();  

#ifdef PPML
  if( HydroMethod == PPM_Local )
    Face = new PPML_InterfacePointerBundle( this );
#endif //PPML

  for(i=0;i<GridRank;i++){
    Scale[i]=(GridRightEdge[i]-GridLeftEdge[i])/(GridDimension[i]-2*DEFAULT_GHOST_ZONES);
  }

  IndexPointerMap ind;
  if( this->IdentifyPhysicalQuantities_2(ind) == FAIL ){
    fprintf(stderr," IdentifyPhysicalQuantities_2 failed\n");return FAIL;}

  fprintf(stderr,"monkey d %d  e %d x %d y %d z %d x %d y %d z %d\n", 
	  ind.D, ind.TE, ind.VX, ind.VY, ind.VZ, ind.BX, ind.BY, ind.BZ);

  fprintf(stderr,"Density %f Pressure %f V0 %f B0 %f \n", DensityIn,Pressure,V0,B0);
  fprintf(stderr,"Scale: %f %f %f\n", Scale[0],Scale[1],Scale[2]);

#ifndef MHDF  
  int ElectricDims[3][3], MagneticDims[3][3];
  float *ElectricField[3], *MagneticField[3];
  size = 0;
  //Set up things that will exist in the full MHD implementation.
  for(field=0;field<3;field++){
    for(int dim=0; dim<3; dim++){
      ElectricDims[field][dim] = GridDimension[dim] + ( (field == dim ) ? 0 : 1 );
      MagneticDims[field][dim] = GridDimension[dim] + ( (field == dim ) ? 1 : 0 );
    }
    size = ElectricDims[field][0]* ElectricDims[field][1]* ElectricDims[field][2];
    ElectricField[field] = new float[ size ];
    size = MagneticDims[field][0]* MagneticDims[field][1]* MagneticDims[field][2];
    MagneticField[field] = new float[ size ];
    
  }
#else
  MHD_Allocate();
#endif

  //Vector Potential. 
  //Due to the similarity in centering, and lack of foresight in naming,
  //I'm using the Electric Field as a Vector Potential to initialize the Magnetic Field.  
  
  field=2;
  for( k=0;k<ElectricDims[field][2];k++)
    for( j=0;j<ElectricDims[field][1];j++)
      for( i=0;i<ElectricDims[field][0];i++){
	index=i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k);
	X=(i-GridStartIndex[0])*Scale[0];
	Y=(j-GridStartIndex[1])*Scale[1];
	
	ElectricField[field][index]=B0*( cos(4*Pi*X)/(4*Pi) + cos(2*Pi*Y)/(2*Pi) );
      }
  
  
  //
  //Curl is B=Curl(A)
  //
  
#ifdef MHDF
  //the last argument indicates that this isn't a time update.
  if( this->MHD_Curl(GridStartIndex, GridEndIndex, 0) == FAIL )
    {fprintf(stderr," error occored in MHD_Curl\n"); return FAIL;}
  
  if( this->CenterMagneticField() == FAIL ) {
    fprintf(stderr," error with CenterMagneticField\n");
    return FAIL;
  }
#else
  
  float dtUsed = -1.0;  //-1 because B = B - dt*(curlE)
  int MHD_CenteringMethod = 2;
  
  FORTRAN_NAME(curl_of_e)(MagneticField[0], MagneticField[1], MagneticField[2],
			  ElectricField[0], ElectricField[1], ElectricField[2],
			  CellWidth[0], CellWidth[1], CellWidth[2],
			  GridDimension, GridDimension +1, GridDimension +2,
			  GridStartIndex, GridEndIndex,
			  GridStartIndex+1, GridEndIndex+1,
			  GridStartIndex+2, GridEndIndex+2,
			  &dtUsed, &MHD_CenteringMethod);
  
  FORTRAN_NAME(center_magnetic_field)
    (MagneticField[0], MagneticField[1], MagneticField[2],
     BaryonField[ind.BX],BaryonField[ind.BY],BaryonField[ind.BZ], BaryonField[ind.TE],
     CellWidth[0], CellWidth[1], CellWidth[2],
     GridDimension, GridDimension + 1, GridDimension +2,
     GridStartIndex, GridEndIndex,
     GridStartIndex+1, GridEndIndex+1,
     GridStartIndex+2, GridEndIndex+2,
     &MHD_CenteringMethod);


#endif //!MHD

  for(k=0;k<GridDimension[2];k++)
    for(j=0;j<GridDimension[1];j++)
      for(i=0;i<GridDimension[0];i++){
	index=i + GridDimension[0]*(j+GridDimension[1]*k);
	X=(i-GridStartIndex[0]+0.5)*Scale[0];
	Y=(j-GridStartIndex[1]+0.5)*Scale[1];

	Vx=-V0*sin(2*Pi*Y);
	Vy=V0*sin(2*Pi*X);


	BaryonField[ ind.D ][index]= DensityIn;

	if( EquationOfState == 0 ){
	  TotalEnergy=GasEnergy + 0.5*DensityIn*(Vx*Vx+Vy*Vy)
	    +0.5*(BaryonField[ind.BX][index]*BaryonField[ind.BX][index]+
		  BaryonField[ind.BY][index]*BaryonField[ind.BY][index]+
		  BaryonField[ind.BZ][index]*BaryonField[ind.BZ][index]);
	  
	  BaryonField[ ind.TE ][index]=TotalEnergy;
	}
	if( DualEnergyFormalism )
	  BaryonField[ ind.GE][index]=GasEnergy;
	BaryonField[ ind.VX ][index]=Vx;
	BaryonField[ ind.VY ][index]=Vy;
	BaryonField[ ind.VZ ][index]=0.0;


      }
  
#ifdef MHDF
  if( HydroMethod == PPM_Local){
    PPML_InitInterfaceMethod = 2;
    PPML_InitInterfaceDataGrid();
    PPML_InitInterfaceMethod = 0;
    PPML_MHD_Clone(0);
  }
  MHD_Deallocate();
#endif
  return SUCCESS;
}



