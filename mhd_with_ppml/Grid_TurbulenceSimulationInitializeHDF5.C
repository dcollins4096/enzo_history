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
 
#define USE_HDF5
 
#ifdef USE_HDF5
 
#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
 
#include "hdf4.h"
 
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

#ifdef PPML
#include "PPML.h"
#endif //PPML 
// HDF5 function prototypes
 
#include "extern_hdf5.h"
 
// function prototypes
 
int HDF5_ReadFile(char *name, int Rank, int Dim[], int StartIndex[],
                  int EndIndex[], int BufferOffset[], float *buffer,
                  inits_type **tempbuffer, int Part, int Npart);
 
int grid::TurbulenceSimulationInitializeGridHDF5(
                          float TurbulenceSimulationInitialDensity,
                          float TurbulenceSimulationInitialTemperature,
#ifdef PPML
			  float TurbulenceSimulationInitialMagneticField[],
                          char *TurbulenceSimulationMagneticNames[],			  
#endif //PPML
                          char *TurbulenceSimulationDensityName,
                          char *TurbulenceSimulationTotalEnergyName,
                          char *TurbulenceSimulationGasEnergyName,
                          char *TurbulenceSimulationVelocityNames[],
                          char *TurbulenceSimulationRandomForcingNames[],
                          int   TurbulenceSimulationSubgridsAreStatic,
                          int   TotalRefinement)
{
 
  /* declarations */
 
#ifdef PPML
  //dcc: future note, pass this in as an arguement so it can be read

  PPML_InterfacePointerBundle  *Face =NULL; //filled after InitInterfaceTypesAndLabels
  int part=0, npart = 0;  //for reading data.

#endif //PPML  
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
  sprintf(pid, "%4.4"ISYM, MyProcessorNumber);
 
  char *logname = new char[strlen(pid)+6+1];
  strcpy(logname, "TSlog.");
  strcat(logname,pid);
 
  if (io_log) {
    log_fptr = fopen(logname, "a");
    fprintf(log_fptr, "\n");
    fprintf(log_fptr, "TSIG ParallelRootGridIO = %"ISYM"\n", ParallelRootGridIO);
    fprintf(log_fptr, "Processor %"ISYM", Target processor %"ISYM"\n",
	    MyProcessorNumber, ProcessorNumber);
    fprintf(log_fptr, "TotalRefinement = %"ISYM"\n", TotalRefinement);
  }
#ifdef SHOULD_BE
  delete logname;
#endif //SHOULD_BE
  /* Determine if the data should be loaded in or not. */

  int ReadData = TRUE, Offset[] = {0,0,0};
 
  if (ParallelRootGridIO == TRUE && TotalRefinement == 1)
    ReadData = FALSE;
 
  if (io_log) fprintf(log_fptr, "ReadData = %"ISYM"\n", ReadData);
 
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
#ifdef PPML
  int BNum[3];
#endif //PPML

#ifndef PPML
  if (TurbulenceSimulationVelocityNames[0] != NULL) {
#endif //PPML.  We want to have the option of auto-generated fields.
    FieldType[NumberOfBaryonFields++] = Density;
#ifdef PPML
    if( EquationOfState == 0 )
#endif //PPML
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;
    FieldType[NumberOfBaryonFields++] = Velocity1;
    vel = NumberOfBaryonFields - 1;
    if (GridRank > 1 OR_MHD)
      FieldType[NumberOfBaryonFields++] = Velocity2;
    if (GridRank > 2 OR_MHD)
      FieldType[NumberOfBaryonFields++] = Velocity3;
#ifdef PPML
    if( MHD_Used == TRUE ){
      BNum[0] = NumberOfBaryonFields;
      FieldType[NumberOfBaryonFields++] = Magnetic1;
      BNum[1] = NumberOfBaryonFields;
      FieldType[NumberOfBaryonFields++] = Magnetic2;
      BNum[2] = NumberOfBaryonFields;
      FieldType[NumberOfBaryonFields++] = Magnetic3;
    }

    if( HydroMethod == PPM_Local ){
      if( this->PPML_InitInterfaceTypesAndLabels() == FAIL ){
	fprintf(stderr," Grid_DiscontInit...: Failure in PPML_InitInterface...\n");
	return FAIL;
      }

    }
  
    IndexPointerMap ind;
    if( this->IdentifyPhysicalQuantities_2(ind) == FAIL ){
      fprintf(stderr, "TurbSimInitGrid: error in identify physical quantities.\n");
    }
#endif //PPML
#ifndef PPML
  }
#endif //PPML.  This closes the VelocityField name requirement.
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
    if (io_log) fprintf(log_fptr, "Allocate %"ISYM" fields, %"ISYM" floats per field\n",
			NumberOfBaryonFields, size);
    for (dim = 0; dim < GridRank; dim++)
      if (io_log) fprintf(log_fptr, "  Field dim %"ISYM" size %"ISYM"\n", dim, GridDimension[dim]);
    printf("Allocating %"ISYM" baryon fields of size %"ISYM"\n", NumberOfBaryonFields, size);
    for (int field = 0; field < NumberOfBaryonFields; field++)
      BaryonField[field] = new float[size];

#ifdef PPML
    if( HydroMethod == PPM_Local )
      Face = new PPML_InterfacePointerBundle( this );

    
    //error check that there are enough possible baryon fields
    if( NumberOfBaryonFields + RandomForcingNumberOfFields > MAX_NUMBER_OF_BARYON_FIELDS ){
      fprintf(stderr, "Fatal error: AppendForcingToBaryonFields: MAX_NUMBER_OF_BARYON_FIELDS not large enough.\n");
      fprintf(stderr, "    Must be at least as large as NumberOfBaryonFields (%"ISYM") + RandomForcingNumberOfFields (%d)\n",
	      NumberOfBaryonFields, RandomForcingNumberOfFields);
      return FAIL;
    }
#endif //PPML
  }

  /* If set, allocate space for RandomForcing fields. */

#ifdef PPML
  if (RandomForcing == TRUE && ReadData == TRUE)
      for (int dim = 0; dim < RandomForcingNumberOfFields; dim++)
        if (RandomForcingField[dim] == NULL)
          RandomForcingField[dim] = new float[size];
#else //PPML
  if (RandomForcing == TRUE && ReadData == TRUE)
      for (int dim = 0; dim < GridRank; dim++)
        if (RandomForcingField[dim] == NULL)
          RandomForcingField[dim] = new float[size];
#endif //PPML 

  /* Read the density field. */
 
#ifdef PPML
  if( HydroMethod == PPM_Local && PPML_InitInterfaceMethod == 0 )
    npart = PPML_NFaces +1; //(+1 for the centered field)
  part = 0;
  if (TurbulenceSimulationDensityName != NULL && ReadData){

    if (HDF5_ReadFile(TurbulenceSimulationDensityName, GridRank, GridDimension,
//  if (ReadFileHDF5(TurbulenceSimulationDensityName, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[0],
              &tempbuffer, part, npart) == FAIL) {
      fprintf(stderr, "Error reading density field.\n");
      return FAIL;
    }
    
    if( HydroMethod == PPM_Local && PPML_InitInterfaceMethod == 0 ){
      for( int face_part = 0; face_part < PPML_NFaces; face_part++ ){
	part++;
	fprintf(stderr,"Read that bastard\n");
	if (HDF5_ReadFile(TurbulenceSimulationDensityName, GridRank,
			  GridDimension, GridStartIndex, GridEndIndex, Offset,
			  Face->All[face_part][ind.D], &tempbuffer, part, npart) == FAIL) {
	  fprintf(stderr, "Error reading velocity field %d, face %d %d.\n", dim, face_part); 
	  return FAIL; 
	}
      }//loop
    }// Read faces.
  }

#else //PPML
  if (TurbulenceSimulationDensityName != NULL && ReadData)
    if (HDF5_ReadFile(TurbulenceSimulationDensityName, GridRank, GridDimension,
//  if (ReadFileHDF5(TurbulenceSimulationDensityName, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[0],
              &tempbuffer, 0, 1) == FAIL) {
      fprintf(stderr, "Error reading density field.\n");
      return FAIL;
    }

#endif //PPML
 
  /* Read the total energy field. */
 
  if (TurbulenceSimulationTotalEnergyName != NULL && ReadData )
#ifdef PPML
  if( EquationOfState == 0 )
#endif //PPML
    if (HDF5_ReadFile(TurbulenceSimulationTotalEnergyName, GridRank,
		    GridDimension, GridStartIndex, GridEndIndex, Offset,
		    BaryonField[1], &tempbuffer, 0, 1) == FAIL) {
      fprintf(stderr, "Error reading total energy field.\n");
      return FAIL;
    }
 
  /* Read the gas energy field. */
 
  if (TurbulenceSimulationGasEnergyName != NULL && DualEnergyFormalism &&
      ReadData)
    if (HDF5_ReadFile(TurbulenceSimulationGasEnergyName, GridRank, GridDimension,
		     GridStartIndex, GridEndIndex, Offset, BaryonField[2],
		     &tempbuffer, 0, 1) == FAIL) {
      fprintf(stderr, "Error reading gas energy field.\n");
      return FAIL;
    }
 
  /* Read the velocity fields. */
#ifdef SHOULD_BE
  //0,1 for 3 different files
  //dim,3 for one 3 component file

  if (TurbulenceSimulationVelocityNames[0] != NULL && ReadData)
    for (dim = 0; dim < GridRank; dim++){

      //If all 3 names match, then all 3 velocity components are in one file.
      //Otherwise, not.
  if( strcmp( TurbulenceSimulationVelocityNames[0] , TurbulenceSimulationVelocityNames[1] ) == 0 ){
    part = dim;
    npart = 3; 
  }else{
    part = 0;
    npart = 1;

#ifdef PPML
    //Interface method = 0 means you'll bring the goods, so they'd better be in the file.
    if( HydroMethod == PPM_Local && PPML_InitInterfaceMethod == 0 )
      npart = PPML_NFaces +1; //(+1 for the centered field)
#endif //PPML
  }
  if (HDF5_ReadFile(TurbulenceSimulationVelocityNames[dim], GridRank,
		   GridDimension, GridStartIndex, GridEndIndex, Offset,
		   BaryonField[vel+dim], &tempbuffer, part, npart) == FAIL) {
    fprintf(stderr, "Error reading velocity field %d.\n", dim);
    return FAIL;
  }
#ifdef PPML
  if( HydroMethod == PPM_Local && PPML_InitInterfaceMethod == 0 ){
    for( int face_part = 0; face_part < PPML_NFaces; face_part++ ){
      part++;
      fprintf(stderr,"stoogie: %"ISYM" %"ISYM" vel %d dim %d\n", face_part, part, vel, dim);
      if (HDF5_ReadFile(TurbulenceSimulationVelocityNames[dim], GridRank,
			GridDimension, GridStartIndex, GridEndIndex, Offset,
			Face->All[face_part][vel+dim], &tempbuffer, part, npart) == FAIL) {
	fprintf(stderr, "Error reading velocity field %d, face %d %d.\n", dim, face_part); 
	return FAIL; 
      }
    }//loop
  }// Read faces.
#endif //PPML
    }//dim
#else //SHOULD_BE
  if (TurbulenceSimulationVelocityNames[0] != NULL && ReadData)
    for (dim = 0; dim < GridRank; dim++)
      if (HDF5_ReadFile(TurbulenceSimulationVelocityNames[dim], GridRank,
		   GridDimension, GridStartIndex, GridEndIndex, Offset,
                  BaryonField[vel+dim], &tempbuffer, 0, 1) == FAIL) {
       //         BaryonField[vel+dim], &tempbuffer, dim, 3) == FAIL) {
	fprintf(stderr, "Error reading velocity field %"ISYM".\n", dim);
	return FAIL;
      }
#endif //SHOULD_BE
  /* Get RandomForcing data */
  
#ifdef SHOULD_BE
  if (RandomForcing == TRUE ){
    if( TurbulenceSimulationRandomForcingNames[0] != NULL ){
#else //SHOULD_BE
      //if (RandomForcing == TRUE && TurbulenceSimulationVelocityNames[0] != NULL ){
#endif //SHOULD_BE
      
      
      /* Read the random forcing fields. */
      
      if (ReadData)
	for (dim = 0; dim < GridRank; dim++)
	  if (HDF5_ReadFile(TurbulenceSimulationRandomForcingNames[dim], GridRank,
			    GridDimension, GridStartIndex, GridEndIndex, Offset,
			    RandomForcingField[dim], &tempbuffer, 0, 1) == FAIL) {
	    //		  RandomForcingField[dim], &tempbuffer, dim, 3) == FAIL) {
	    fprintf(stderr, "Error reading RandomForcing field %"ISYM".\n", dim);
	    return FAIL;
	  }
    } else {
      /* OR: copy random forcing fields from initial velocities. */
      if( ReadData ){
	if (TurbulenceSimulationVelocityNames[0] == NULL) {
	  //OR compute a brand new one!
	  if( this->ComputeRandomForcingFields(0) == FAIL ){
	    fprintf(stderr," Error in ComputeRandomForcingFields\n"); return FAIL;}
	  
	}
      }//read
      
#ifdef PPML
      //This will get re-coppied in PPML_InitInterfaceGrid.
    
      if (ReadData == TRUE){
	int counter = 0;
	for(dim=0; dim < GridRank; dim++){
	  for (i = 0; i < size; i++)
	    RandomForcingField[counter][i] = BaryonField[vel+dim][i];
	  counter++;
	}
	if( HydroMethod == PPM_Local) {
	  for( int face_part = 0; face_part < PPML_NFaces; face_part++ )
	    for(dim=0; dim < GridRank; dim++){
	      for (i = 0; i < size; i++)
		RandomForcingField[counter][i] = Face->All[face_part][vel+dim][i];
	      counter++;
	    }
	}//ppm_local
	
      }


#else //PPML
  if (ReadData == TRUE)
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < size; i++)
        RandomForcingField[dim][i] = BaryonField[vel+dim][i];
#endif  //PPML
    }
  }//Random Forcing 

#ifdef PPML

  if (TurbulenceSimulationMagneticNames[0] != NULL && ReadData){
    for (dim = 0; dim < GridRank; dim++){
      fprintf(stderr,"Steak for breakfast\n");
      //If all 3 names match, then all 3 velocity components are in one file.
      //Otherwise, not.
      if( strcmp( TurbulenceSimulationMagneticNames[0] , TurbulenceSimulationMagneticNames[1] ) == 0 ){
	part = dim;
	npart = 3; 
      }else{
	part = 0;
	npart = 1;
	
	//Interface method = 0 means you'll bring the goods, so they'd better be in the file.
      }
      if( HydroMethod == PPM_Local && PPML_InitInterfaceMethod == 0 )
	npart = PPML_NFaces +1; //(+1 for the centered field)

      if (HDF5_ReadFile(TurbulenceSimulationMagneticNames[dim], GridRank,
			GridDimension, GridStartIndex, GridEndIndex, Offset,
			BaryonField[ind.B[dim]], &tempbuffer, part, npart) == FAIL) {
	fprintf(stderr, "Error reading magnetic field %d.\n", dim);
	return FAIL;
      }
      fprintf(stderr,"Magnetic Field Names %d %s\n", dim, TurbulenceSimulationMagneticNames[dim]);
      if( HydroMethod == PPM_Local && PPML_InitInterfaceMethod == 0 ){
	for( int face_part = 0; face_part < PPML_NFaces; face_part++ ){
	  part++;
	  fprintf(stderr,"stoogie: %"ISYM" %"ISYM" magnetic %d dim %d\n", face_part, part, vel, dim);
	  if (HDF5_ReadFile(TurbulenceSimulationMagneticNames[dim], GridRank,
			    GridDimension, GridStartIndex, GridEndIndex, Offset,
			    Face->All[face_part][ind.B[dim]], &tempbuffer, part, npart) == FAIL) {
	    fprintf(stderr, "Error reading magneticfield %d, face %d %d.\n", dim, face_part); 
	    return FAIL; 
	  }
	}//loop
      }// Read faces.
    }//dim
  }//magnetic names
#endif //PPML ok

  /* If they were not read in above, set the total & gas energy fields now. */
 
  if (/*TurbulenceSimulationVelocityNames[0] != NULL &&*/ ReadData) {
#ifndef PPML
    if (TurbulenceSimulationDensityName == NULL)
      for (i = 0; i < size; i++)
        BaryonField[0][i] = TurbulenceSimulationInitialDensity;
#else 
    if (TurbulenceSimulationDensityName == NULL){
      for (i = 0; i < size; i++)
        BaryonField[0][i] = TurbulenceSimulationInitialDensity;
      if( HydroMethod == PPM_Local ){
	for( int face_part = 0; face_part < PPML_NFaces; face_part++)
	for (i = 0; i < size; i++) Face->All[face_part][0][i] = TurbulenceSimulationInitialDensity;
      }//ppml
    }
#endif //PPML
    if (TurbulenceSimulationTotalEnergyName == NULL)
#ifdef PPML
    if( EquationOfState == 0 )
#endif //PPML
	
      for (i = 0; i < size; i++)
        BaryonField[1][i] = TurbulenceSimulationInitialTemperature/(Gamma-1.);
 
    if (TurbulenceSimulationGasEnergyName == NULL && DualEnergyFormalism)
      for (i = 0; i < size; i++)
        BaryonField[2][i] = BaryonField[1][i];
 
    if (TurbulenceSimulationTotalEnergyName == NULL &&
        HydroMethod != Zeus_Hydro )
#ifdef PPML
    if( EquationOfState == 0 )
#endif //PPML

      for (dim = 0; dim < GridRank; dim++)
        for (i = 0; i < size; i++)
          BaryonField[1][i] +=
            0.5 * BaryonField[vel+dim][i] * BaryonField[vel+dim][i];
  }
#ifdef PPML
    if( TurbulenceSimulationMagneticNames[0] == NULL && MHD_Used == TRUE && ReadData == TRUE){
      for(int field=0;field<3;field++){
	if( BaryonField[ BNum[field] ] == NULL )
	  BaryonField[ BNum[field] ] = new float[size];
	for(i=0;i<size;i++)
	  BaryonField[ BNum[field] ][i]=TurbulenceSimulationInitialMagneticField[field];
	if( HydroMethod == PPM_Local ){
	  for( int face_part = 0; face_part < PPML_NFaces; face_part++){
	    for (i = 0; i < size; i++) Face->All[ face_part ][ BNum[field] ][i] = 
					 TurbulenceSimulationInitialMagneticField[field];
	  }
	}//ppm
#ifdef MHD

	MagneticField[field] = new float[MagneticSize[field]];
	for(i=0;i<MagneticSize[field];i++)
	  MagneticField[field][i]=TurbulenceSimulationInitialMagneticField[field];
#endif //MHD
      }//field
    }//MHD_Used
#endif //PPML
 
  } // end: if (NumberOfBaryonFields > 0)
 
  } // end: if (ProcessorNumber == MyProcessorNumber)
 
  OldTime = Time;
 
  if (io_log) fclose(log_fptr);

#ifdef PPML
  if( Face != NULL )
    delete Face;
#endif //PPML  

  return SUCCESS;
}
 
#endif
