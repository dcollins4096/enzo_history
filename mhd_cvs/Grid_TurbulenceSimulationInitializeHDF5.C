/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A TURBULENCE SIMULATION)
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

//#define USE_HDF5 done by command line. dcc.

#ifdef USE_HDF5

#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

//#include "performance.h"
#include "hdf4.h"  /* HDF4 int32, float32 and VOIDP definitions */

#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "flowdefs.h"
#include "error.h"

#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif

// HDF5 function prototypes

#include "extern_hdf5.h"

// function prototypes

int ReadFileHDF5(char *name, int Rank, int Dim[], int StartIndex[],
                  int EndIndex[], int BufferOffset[], float *buffer,
		double **tempbuffer, int Part, int Npart);

int grid::TurbulenceSimulationInitializeGridHDF5(
                          float TurbulenceSimulationInitialDensity,
                          float TurbulenceSimulationInitialTemperature,
                          float TurbulenceSimulationInitialMagneticField[],
			  char *TurbulenceSimulationMagneticNames[],
                          char *TurbulenceSimulationDensityName,
                          char *TurbulenceSimulationTotalEnergyName,
                          char *TurbulenceSimulationGasEnergyName,
                          char *TurbulenceSimulationVelocityNames[],
                          char *TurbulenceSimulationRandomForcingNames[],
                          int   TurbulenceSimulationSubgridsAreStatic,
                          int   TotalRefinement)
{


  /* declarations */


  int idim, dim, i, j, vel;
  int DeNum;

  int ExtraField[2];

  inits_type *tempbuffer = NULL;

  FILE *log_fptr;

#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

  /*
  if ( NumberOfProcessors > 64 )
    if (ParallelRootGridIO != TRUE) {
      fprintf(stderr, "ParallelRootGridIO MUST be set for > 64 cpus!\n");
      return FAIL;
    }
  */

  char pid[5];
  sprintf(pid, "%4.4d", MyProcessorNumber);

  char *logname = new char[strlen(pid)+6+1];
  strcpy(logname, "TSlog.");
  strcat(logname,pid);

  if (io_log) {
    log_fptr = fopen(logname, "a");
    fprintf(log_fptr, "\n");
    fprintf(log_fptr, "TSIG ParallelRootGridIO = %d\n", ParallelRootGridIO);
    fprintf(log_fptr, "Processor %d, Target processor %d\n", 
	    MyProcessorNumber, ProcessorNumber);
    fprintf(log_fptr, "TotalRefinement = %d\n", TotalRefinement);
  }

  /* Determine if the data should be loaded in or not. */

  int ReadData = TRUE, Offset[] = {0,0,0};

  if (ParallelRootGridIO == TRUE && TotalRefinement == 1)
    ReadData = FALSE;

  /*
  if (TurbulenceSimulationVelocityNames[0] == NULL) 
    ReadData = FALSE;
  */

  if (io_log) fprintf(log_fptr, "ReadData = %d\n", ReadData);

  /* Calculate buffer Offset (same as Grid unless doing ParallelRootGridIO 
     (TotalRefinement = -1 if used as a signal that we should really load 
     in the data regardless of the value of ParallelRootGridIO). */

  if (ParallelRootGridIO == TRUE && TotalRefinement == -1)
    for (dim = 0; dim < GridRank; dim++)
      Offset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/
			 CellWidth[dim][0]);

  /*----------------------------------------------------*/
  /* Create baryon fields. */

  NumberOfBaryonFields = 0;

  //if (TurbulenceSimulationVelocityNames[0] != NULL) { dcc kludge: removed this check.
    FieldType[NumberOfBaryonFields++] = Density;
#ifdef ATHENA
    if( EquationOfState == 0 )
      FieldType[NumberOfBaryonFields++] = TotalEnergy;
#else //ATHENA
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
#endif //ATHENA
    if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;
    FieldType[NumberOfBaryonFields++] = Velocity1;
    vel = NumberOfBaryonFields - 1;
    if (GridRank > 1)
      FieldType[NumberOfBaryonFields++] = Velocity2;
    if (GridRank > 2)
      FieldType[NumberOfBaryonFields++] = Velocity3;
    //}


  /* Set the subgrid static flag. */

  SubgridsAreStatic = TurbulenceSimulationSubgridsAreStatic;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber == MyProcessorNumber) {

  /* Skip following if NumberOfBaryonFields == 0. */

  if (NumberOfBaryonFields > 0) {

  /* Determine the size of the fields. */

  int size = 1;

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Allocate space for the fields. */

  if (ReadData == TRUE) {
    if (io_log) fprintf(log_fptr, "Allocate %d fields, %d floats per field\n", 
			NumberOfBaryonFields, size);
    for (dim = 0; dim < GridRank; dim++)
      if (io_log) fprintf(log_fptr, "  Field dim %d size %d\n", dim, GridDimension[dim]);
    printf("Allocating %d baryon fields of size %d\n", NumberOfBaryonFields, size);
    for (int field = 0; field < NumberOfBaryonFields; field++)
      BaryonField[field] = new float[size];
    if( MHD_Used == TRUE ){
      for (int field = 0; field < 3; field++ ){
	MagneticField[field] = new float[ MagneticSize[field] ];
	CenteredB[field] = new float[size];
	for( i=0;i<size;i++)
	  CenteredB[field][size] = 0.0;
      }
    }
  }

  /* If set, allocate space for RandomForcing fields. */

  if (RandomForcing == TRUE && ReadData == TRUE)
      for (int dim = 0; dim < GridRank; dim++)
        if (RandomForcingField[dim] == NULL)
          RandomForcingField[dim] = new float[size];

  /* Read the density field. */

  if (TurbulenceSimulationDensityName != NULL && ReadData)
    if (ReadFileHDF5(TurbulenceSimulationDensityName, GridRank, GridDimension,
//  if (ReadFileHDF5(TurbulenceSimulationDensityName, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[0],
              &tempbuffer, 0, 1) == FAIL) {
      fprintf(stderr, "Error reading density field.\n");
      return FAIL;
    }

  /* Read the total energy field. */

  if (TurbulenceSimulationTotalEnergyName != NULL  && ReadData
#ifdef ATHENA
      && EquationOfState == 0
#endif //ATHENA
      )
    if (ReadFileHDF5(TurbulenceSimulationTotalEnergyName, GridRank, 
		    GridDimension, GridStartIndex, GridEndIndex, Offset,
		    BaryonField[1], &tempbuffer, 0, 1) == FAIL) {
      fprintf(stderr, "Error reading total energy field.\n");
      return FAIL;
    }

  /* Read the gas energy field. */

  if (TurbulenceSimulationGasEnergyName != NULL && DualEnergyFormalism && 
      ReadData)
    if (ReadFileHDF5(TurbulenceSimulationGasEnergyName, GridRank, GridDimension,
		     GridStartIndex, GridEndIndex, Offset, BaryonField[2],
		     &tempbuffer, 0, 1) == FAIL) {
      fprintf(stderr, "Error reading gas energy field.\n");
      return FAIL;
    }

  /* Read the velocity fields. */

  //0,1 for 3 different files
  //dim,3 for one 3 component file
  int part=0, npart = 0;
  if (TurbulenceSimulationVelocityNames[0] != NULL && ReadData)
    for (dim = 0; dim < GridRank; dim++){
  if( strcmp( TurbulenceSimulationVelocityNames[0] , TurbulenceSimulationVelocityNames[1] ) == 0 ){
    part = dim;
    npart = 3;
  }else{
    part = 0;
    npart = 1;
  }

      if (ReadFileHDF5(TurbulenceSimulationVelocityNames[dim], GridRank,
		       GridDimension, GridStartIndex, GridEndIndex, Offset,
		       BaryonField[vel+dim], &tempbuffer, part, npart) == FAIL) {
	fprintf(stderr, "Error reading velocity field %d.\n", dim);
	return FAIL;
      }
    }
  /* Get RandomForcing data */

  if (RandomForcing == TRUE && TurbulenceSimulationRandomForcingNames[0] != NULL) {
    //if (RandomForcing == TRUE && TurbulenceSimulationVelocityNames[0] != NULL) {

  /* Read the random forcing fields. */

    if (ReadData)
      for (dim = 0; dim < GridRank; dim++)
	if (ReadFileHDF5(TurbulenceSimulationRandomForcingNames[dim], GridRank,
			  GridDimension, GridStartIndex, GridEndIndex, Offset,
			 RandomForcingField[dim], &tempbuffer, 0, 1) == FAIL) {
			 //RandomForcingField[dim], &tempbuffer, dim, 3) == FAIL) {
	  fprintf(stderr, "Error reading RandomForcing field %d.\n", dim);
	  return FAIL;
      }
  }

  else {

  /* OR: copy random forcing fields from initial velocities. */

    if( ReadData ){ //&& RandomForcing == TRUE){
	if (TurbulenceSimulationVelocityNames[0] == NULL) {
	  //OR compute a brand new one!
	  if( this->ComputeRandomForcingFields(0) == FAIL ){
	    fprintf(stderr," Error in ComputeRandomForcingFields\n"); return FAIL;}
	  
	}
      }//read
      
      if (ReadData == TRUE && RandomForcing == TRUE)
	for (dim = 0; dim < GridRank; dim++)
	  for (i = 0; i < size; i++)
	    RandomForcingField[dim][i] = BaryonField[vel+dim][i];
      

  }
  
  if (TurbulenceSimulationMagneticNames[0] != NULL && ReadData && MHD_Used){
    for (int field = 0; field < GridRank; field++){
      //If all 3 names match, then all 3 velocity components are in one file.
      //Otherwise, not.
      if( strcmp( TurbulenceSimulationMagneticNames[0] , TurbulenceSimulationMagneticNames[1] ) == 0 ){
	part = field;
	npart = 3; 
      }else{
	part = 0;
	npart = 1;
	
	//Interface method = 0 means you'll bring the goods, so they'd better be in the file.
      }
      if (ReadFileHDF5(TurbulenceSimulationMagneticNames[field], GridRank,
			MagneticDims[field], MHDStartIndex[field], MHDEndIndex[field], Offset,
			MagneticField[field], &tempbuffer, part, npart) == FAIL) {
	fprintf(stderr, "Error reading magnetic field %d.\n", field);
	return FAIL;
      }
      //fprintf(stderr,"MagNames! %d %s\n", field, TurbulenceSimulationMagneticNames[field]);
      //fprintf(stderr,"MagDims[%d] = %d %d %d\n", field
      //, MagneticDims[field][0], MagneticDims[field][1], MagneticDims[field][2]);
    }//dim
    this->CenterMagneticField();
  }//magnetic names, mhd_used

   /* If they were not read in above, set the total & gas energy fields now. */
  
  //if (TurbulenceSimulationVelocityNames[0] != NULL && ReadData) {  //dcc: removed for in house generator.
  if ( ReadData) {
    if (TurbulenceSimulationDensityName == NULL)
      for (i = 0; i < size; i++)
        BaryonField[0][i] = TurbulenceSimulationInitialDensity;

    
    if (TurbulenceSimulationTotalEnergyName == NULL
#ifdef ATHENA
	&& EquationOfState == 0
#endif //ATHENA
	)
      if( MHD_Used ){
	for (i = 0; i < size; i++)
	  BaryonField[1][i] = 
	    TurbulenceSimulationInitialTemperature*TurbulenceSimulationInitialTemperature*
	    BaryonField[0][i]/(Gamma*(Gamma-1.));
      }else{
	for (i = 0; i < size; i++)
	  BaryonField[1][i] = TurbulenceSimulationInitialTemperature/(Gamma-1.);
      }
      
    if (TurbulenceSimulationGasEnergyName == NULL && DualEnergyFormalism)
      for (i = 0; i < size; i++)
        BaryonField[2][i] = BaryonField[1][i];
    
    if( MHD_Used == TRUE )
      if( TurbulenceSimulationMagneticNames[0] == NULL ){
	
	for(int field=0;field<3;field++){
	  
	  if( MagneticSize[field] <0 ){
	    fprintf(stderr,"Shit!  Improperly defined magnetic size.\n");
	    return FAIL;}
	  CenteredB[field] = new float[size];
	  for(i=0;i<size;i++)
	    CenteredB[field][i]=TurbulenceSimulationInitialMagneticField[field];
	  MagneticField[field] = new float[MagneticSize[field]];
	  for(i=0;i<MagneticSize[field];i++)
	    MagneticField[field][i]=TurbulenceSimulationInitialMagneticField[field];
	  
	}//field
      }
    
    if (TurbulenceSimulationTotalEnergyName == NULL &&
        HydroMethod != Zeus_Hydro 
#ifdef ATHENA
	&& EquationOfState == 0
#endif //ATHENA
	)

      if( MHD_Used == TRUE ){
	for( dim = 0; dim < GridRank; dim++){
	  for(i=0;i<size;i++){
	    BaryonField[1][i] +=
	      0.5*BaryonField[0][i]*(BaryonField[vel+dim][i] * BaryonField[vel+dim][i])
	      +0.5*CenteredB[dim][i]*CenteredB[dim][i];

	  }
	}
      }else{
	
	for (dim = 0; dim < GridRank; dim++)
	  for (i = 0; i < size; i++){
	    BaryonField[1][i] +=
	      0.5 * BaryonField[vel+dim][i] * BaryonField[vel+dim][i];
	  }
	}

  }

  //dcc kludge: initialize data fields by hand, here.
  
  if( TurbulenceSimulationVelocityNames[0] == NULL && 0 == 1){

    int size = 1, i,j,k, index;
    float X[3], Scale[3], Center[3];
    fprintf(stderr,"InitializeByHand\n");
    
    
    this->AllocateGrids();
    
    for (dim = 0; dim < GridRank; dim++){
    size *= GridDimension[dim];
    Scale[dim]=(GridRightEdge[dim]-GridLeftEdge[dim])/
    (GridDimension[dim]-2*DEFAULT_GHOST_ZONES);
    Center[dim]=(GridRightEdge[dim]-GridLeftEdge[dim])/2.0;
    
    }
    for(dim=0;dim<GridRank;dim++)
    RandomForcingField[dim]= new float[size];
    
    float Vrms    = RandomForcingMachNumber/
    sqrt(TurbulenceSimulationInitialTemperature);
    
    for(k=0;k<GridDimension[2];k++)
    for(j=0;j<GridDimension[1];j++)
    for(i=0;i<GridDimension[0];i++){
    //Might use these later...
    X[0]=(i-GridStartIndex[0])*Scale[0]-Center[0];
    X[1]=(j-GridStartIndex[1])*Scale[1]-Center[1];
    X[2]=(k-GridStartIndex[2])*Scale[2]-Center[2];
    index = i+GridDimension[0]*(j+GridDimension[1]*k);
    
    for (dim = 0; dim < GridRank; dim++) {
    
    RandomForcingField[dim][index] = Vrms* ((float)rand()/(float)(RAND_MAX)- .5);
    
    }
    BaryonField[0][index]=TurbulenceSimulationInitialDensity;
    
    #ifdef ATHENA
    if( EquationOfState == 0 )
    #endif //ATHENA
    BaryonField[1][index]=
    TurbulenceSimulationInitialTemperature/((Gamma-1));
    
    BaryonField[vel + 0][index]=0.0;
    BaryonField[vel + 1][index]=0.0;
    BaryonField[vel + 2][index]=0;
    if( MHD_Used == TRUE ){
    CenteredB[0][index]=TurbulenceSimulationInitialMagneticField[0];
    CenteredB[1][index]=TurbulenceSimulationInitialMagneticField[1];
    CenteredB[2][index]=TurbulenceSimulationInitialMagneticField[2];
    #ifdef ATHENA
    if( EquationOfState == 0 )
    #endif //ATHENA
    BaryonField[1][index]+= 0.5*(CenteredB[0][index]*CenteredB[0][index]+
    CenteredB[1][index]*CenteredB[1][index]+
    CenteredB[2][index]*CenteredB[2][index]);
    
    
    
    }
    
    }
    if( MHD_Used == TRUE ){
    for(int field=0;field<3;field++)
    for(k=0; k<MagneticDims[field][2]; k++)
    for(j=0; j<MagneticDims[field][1]; j++)
    for(i=0; i<MagneticDims[field][0];i++){
    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
    
    
    
    MagneticField[field][index]=TurbulenceSimulationInitialMagneticField[field];
    
    }
    
    }//MHD_Used
  } // Initialize By Hand.

  } // end: if (NumberOfBaryonFields > 0)

  } // end: if (ProcessorNumber == MyProcessorNumber)

  OldTime = Time;

  if (io_log) fclose(log_fptr);

  return SUCCESS;
}

#endif
