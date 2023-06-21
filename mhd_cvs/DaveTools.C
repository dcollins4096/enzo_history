#include <assert.h>
#include "hdf5.h"
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "extern_hdf5.h"

#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#define DCCFORM "5.2f"
#define DCCFORM2 "10.8f"

#include "extern_hdf5.h"
//hey:  this ACRU thing doesn't work, but it's good to remind me what the acuracy is.

#include "fortran.def"

//returns the wall time.
extern FILE * wall_ptr;
void wall_time_start(){
  char filename[20];
  sprintf(filename,"data_time_%04d",MyProcessorNumber);
  wall_ptr = fopen(filename,"w");
}
void wall_time_flush(){
  fclose(wall_ptr);
  char filename[20];
  sprintf(filename,"data_time_%04d",MyProcessorNumber);
  wall_ptr = fopen(filename,"a");

}
void wall_time_stop(){
  fclose(wall_ptr);
}  
void wall_time (char * string)
{
  //if( MyProcessorNumber != ROOT_PROCESSOR)
  //return;
  // get current time in seconds (specifically the number of 
  // seconds since 00:00:00 UTC January 1 1970)              
  
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv,&tz);
  fprintf(wall_ptr,"%d TIME %s %f \n", MyProcessorNumber, string, tv.tv_sec + 1e-6*tv.tv_usec);
  fflush(wall_ptr);
}

extern "C" void
FORTRAN_NAME(tvtoolf)(int * YesProblem, int * nx, int * ny, int * nz,
		      float *d, float *e, float *vx, float *vy, float *vz,
		      float * bcx,float * bcy,float * bcz,
		      float *bx, float *by, float *bz);


//This routine checks the list WriteInThisA for the given flag value.
int WriteInThisF(int flag){
  
  int WriteInThisLocal = FALSE;
  for( int SnapperJoe = 0; SnapperJoe < N_DbgWrites; SnapperJoe++){
    if( WriteInThisA[SnapperJoe] == flag ) {
      WriteInThisLocal = TRUE;
    }
  }
  
  return WriteInThisLocal;
}

//Writes an HDF5 cube.  Quick and dirty.  
// WriteCube(array, [nx,ny,nz], "ID string", dNum, gNum)
// prints out ID string to std err, along with file name.
// Filename is data111(dNum).grid(gNum)
// datasets are all called "monkey"
// Flow:
// 1.) create filename, print message
// 2.) define size of float 
// 3.) create file
// 3.5) invert array dimensions
// 4.) create dataspace, set
// 5.) close file,space,set.
void WriteCube(float * array, int Dims[], char* string, int dNum, int gNum){
  
  hid_t       file_id, dataset_id, dataspace_id, float_type_id;
  herr_t      status, h5_status, h5_error = -1;
  int FieldRankOut = 3;
  hsize_t     DimsInv[FieldRankOut];
  
  char filename[20];
  
  sprintf(filename, "data111%4.4d.grid%4.4d",dNum,gNum);
  fprintf(stderr,"GPFS WriteCube: %s %s [%d,%d,%d]\n", string, filename, Dims[0],Dims[1],Dims[2]);
  
#define floatdcc double  
  int jj = sizeof(floatdcc);
  switch(jj)
    {
    case 0:
      float_type_id = H5T_NATIVE_INT;
      break;
    case 4:
      float_type_id = HDF5_R4;
      break;
    case 8:
      float_type_id = HDF5_R8;
      break;
    case 16:
      float_type_id = HDF5_R16;
      break;
    default:
      printf("INCORRECT FLOAT DEFINITION\n");
    }
  
  
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  for(int moo=0;moo<FieldRankOut;moo++)
    DimsInv[moo]=Dims[FieldRankOut-moo-1];
  
  
  //Create Dataspace 
  dataspace_id=H5Screate_simple(FieldRankOut, DimsInv, NULL);
  
  //create set
  //                       duh, name,      datatype,  shape of data, Something I dont get
  dataset_id = H5Dcreate(file_id, "monkey", float_type_id, dataspace_id, H5P_DEFAULT);
  
  //Write the Data Set
  //                (set, memory type, mem. space, file space, transfer shit, actual data)
  fprintf(stderr,"Writing set\n");
  status = H5Dwrite(dataset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		    array);
  
  
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  
  
}


//a totalview/ debugging tool.
int grid::TVtool(char * Label){
  
  if( MyProcessorNumber != ProcessorNumber )
    return SUCCESS;
  
  return SUCCESS;
  int YesProblem = 0;
  
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
  
  //grid identification
  if(0==1)
    if( fabs( CellWidth[0][0] - 1.0/16 ) < 1e-7 )
      if( fabs( GridLeftEdge[0] - 0.0 ) < 1e-7 )
	if( fabs( GridLeftEdge[1] - 0.0 ) < 1e-7 )
	  if( fabs( GridLeftEdge[2] - 0.0 ) < 1e-7 ){
	    fprintf(stderr,"this is the grid you want to monitor\n");
	  }      
  
  FORTRAN_NAME(tvtoolf)(&YesProblem, &GridDimension[0],&GridDimension[1],&GridDimension[2],
			BaryonField[DensNum], BaryonField[TENum],
			BaryonField[Vel1Num],  BaryonField[Vel2Num], BaryonField[Vel3Num],
			CenteredB[0],CenteredB[1],CenteredB[2],
			MagneticField[0],MagneticField[1],MagneticField[2]);
  
  if( YesProblem == 1 ){
    fprintf(stderr,"TVtool; Width*16=%f,corner = %f %f %f\n",
	    CellWidth[0][0]*16, GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
    return FAIL;
  }
  
  return SUCCESS;
}

int grid::CreateLevelNumberField(int level){
  fprintf(stderr,"Creating DataLabelField. NumberOfBF %d, level %d\n",
	  NumberOfBaryonFields, level);
  
  int i,j,k, size=1,index;
  NumberOfBaryonFields += 1;
  DataLabel[NumberOfBaryonFields-1] = "LevelNumberField";
  DataUnits[NumberOfBaryonFields-1] = "Dimensionless";
  
  for( i=0;i<GridRank; i++){
    size *= GridDimension[i];
  }
  
  BaryonField[NumberOfBaryonFields-1] = new float[size];
  for(i=0;i<size;i++){
    BaryonField[NumberOfBaryonFields-1][i] = (float) level;
  }
  
  for(k=GridStartIndex[2];k<=GridEndIndex[2];k++)
    for(j=GridStartIndex[1];j<=GridEndIndex[1];j++)
      for(i=GridStartIndex[0];i<=GridEndIndex[0];i++){
	index=index0(i,j,k);
	if(i==GridStartIndex[0] || i==GridEndIndex[0] || 
	   j == GridStartIndex[1] || j==GridEndIndex[1] ||
	   k == GridStartIndex[2] || k== GridEndIndex[2])
	  BaryonField[NumberOfBaryonFields-1][index] += 0.1*level;
      }
  
  return SUCCESS;
}

int grid::CheckForNans(char * Label){
  return SUCCESS;
  if( TVtool(Label) == FAIL ) 
    return FAIL;
  int field, i,j,k, index, FoundOne = FALSE;
  float Why;
  int Zero=0;
  if( MyProcessorNumber != ProcessorNumber )
    return SUCCESS;
  
  Zero++;
  Zero--;
  
  for(field=0;field<3;field++){
    //for(i=0; i<MagneticSize[field];i++){
    for(k=MHDStartIndex[field][2]; k<=MHDEndIndex[field][2]; k++)
      for(j=MHDStartIndex[field][1]; j<=MHDEndIndex[field][1]; j++)
	for(i=MHDStartIndex[field][0]; i<=MHDEndIndex[field][0];i++){
	  index = indexba(i,j,k,field);
	  Why = MagneticField[field][index];
	  //look for nans
	  if( Why != Why ){
	    fprintf(stderr,"   nnann, %d %d", i, field);
	    FoundOne = TRUE;
	  }
	  
	  //This section will search for specific values of the magnetic field, and fail on contact.
	  //Only useful for debugging memory management, grid manipulation.
	  if( 0 == 1 ) {
	    if( Why > 2.0 && 0 == 1 ) {
	      fprintf(stderr, "CFN2 %d %d %f %s %f \n" , i, field, MagneticField[field][index], Label, CellWidth[0][0]*16);
	      return FAIL;
	    }
	    //look for inappropriate projection, fill, alter.  
	    if( fabs(CellWidth[0][0] -1.0/32 ) > 1e-6 ) {
	      if( MagneticField[field][i] < 10.1 && MagneticField[field][i] > 6.9 ) 
		fprintf(stderr, " wtf?\n");
	      
	      if( fabs( MagneticField[field][i] -1.0 ) > 1e-6 ) {
		if( fabs( MagneticField[field][i]  ) > 1e-6 ) {
		  fprintf(stderr, "CFN %d %d %f %s %f \n" , i, field, MagneticField[field][index], Label, CellWidth[0][0]*16);
		  return FAIL;
		}
	      }
	    }
	  }//0==1
	}//i
  }//field
  
  if( FoundOne == TRUE ){
    fprintf(stderr,"CFN y : %s \n", Label);
    return FAIL;
  }
  return SUCCESS;        
}
void grid::TotalMass(char * Label){
  int i,j,k;
  if( 0==1){
    //PreviousMass is relative to the first time this gets called.  (well, as of now.)
    int Grundel;
    float TotalMass=0, Max = -1, Min = 1e20;
    int i,j,k, size=GridDimension[0]*GridDimension[1]*GridDimension[2];
    int index;
    int ncells = 0;
    /*
      for(i=0;i<size;i++)
      TotalMass += BaryonField[0][i];
    */
    for(i=GridStartIndex[0];i<=GridEndIndex[0];i++)
      for(j=GridStartIndex[1];j<=GridEndIndex[1];j++)
	for(k=GridStartIndex[2];k<=GridEndIndex[2];k++){
	  index = index0(i,j,k);
	  TotalMass += BaryonField[0][index];
	  if( BaryonField[0][i] > Max ) Max = BaryonField[0][index];
	  if( BaryonField[0][i] < Min ) Min = BaryonField[0][index];
	  
	  ncells++;
	}
    fprintf(stderr,"TotalMass %f dM/M %e %s\n",
	    TotalMass, (TotalMass - PreviousMass)/TotalMass, Label);
    
    PreviousMass = TotalMass;
  }
}

int grid::TestForAccelerationField(){//int MaybeIWantAnInt, float MaybeIWantAFloat){
  fprintf(stderr,"\n");
  for(int field=0;field<GridRank;field++){
    if( AccelerationField[field] != NULL )
      fprintf(stderr, " AccelerationField[%d] Exists\n", field);
    else
      fprintf(stderr, " AccelerationField[%d] Does Not Exist\n", field);
  }
  return SUCCESS;
}

int grid::IsItShit(char * stringn){    

  int size=0, field=0, index1, index2, i,j,k;
  i=4;j=4;k=4;

  index1=index0(i,j,k);
  index2=index0(i,j-1,k-1);
  if( fabs(BaryonField[1][index1]-BaryonField[1][index2]) > 1e-5 ){
    fprintf(stderr,"Yes, it is shit.  In %s\n", stringn);
    //return FAIL;
  }
  return SUCCESS;
}

void commas(char *str, int numin){

  char temp[100]="";
  int num = numin;
  sprintf(str, "");
  for(int i=10; i>=0; i--){
    int dig = (int) pow((double) 10,(double) i);
    if((int) numin/dig < 1 ) continue;
    sprintf(temp, "%i", num/dig);
    strcat(str, temp);
    if( i % 3 == 0 && i != 0 ) strcat(str, ",");
    num=num-(num/dig)*dig;
  }

}


void PoutF( char * string, float f1 = -12345.6,float f2 = -12345.6,float f3 = -12345.6,
float f4 = -12345.6,float f5 = -12345.6,float f6 = -12345.6 ){
  
  if( MHD_pout != TRUE ) return;

  char filename[50];
  sprintf(filename, "file.%d", MyProcessorNumber);
  
  FILE * fptr = fopen(filename, "a");
  
  fprintf( fptr, "%d ", MyProcessorNumber );
  
  fprintf( fptr, string );

  if( f1 !=  -12345.6 ) fprintf( fptr, " %"DCCFORM2" ", f1);
  if( f2 !=  -12345.6 ) fprintf( fptr, " %"DCCFORM2" ", f2);
  if( f3 !=  -12345.6 ) fprintf( fptr, " %"DCCFORM2" ", f3);
  if( f4 !=  -12345.6 ) fprintf( fptr, " %"DCCFORM2" ", f4);
  if( f5 !=  -12345.6 ) fprintf( fptr, " %"DCCFORM2" ", f5);
  if( f6 !=  -12345.6 ) fprintf( fptr, " %"DCCFORM2" ", f6);

  fprintf( fptr, "\n");

  fclose( fptr );
}

void Pout( char * string, int i1 = -12345, int i2 = -12345,
           int i3 = -12345, int i4 = -12345, int i5 = -12345)
{
  if( MHD_pout != TRUE ) return;

  char filename[50];
  sprintf(filename, "file.%d", MyProcessorNumber);
  
  FILE * fptr = fopen(filename, "a");
  
  fprintf( fptr, "%d ", MyProcessorNumber );
  
  fprintf( fptr, string );

  if( i1 != -12345 ) fprintf( fptr, " %d ", i1 );
  if( i2 != -12345 ) fprintf( fptr, " %d ", i2 );
  if( i3 != -12345 ) fprintf( fptr, " %d ", i3 );
  if( i4 != -12345 ) fprintf( fptr, " %d ", i4 );
  if( i5 != -12345 ) fprintf( fptr, " %d ", i5 );
  
  fprintf( fptr, "\n");

  fclose( fptr );
}

