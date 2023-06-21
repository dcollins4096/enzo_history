/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COSMOLOGY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  James Bordner, June 2003   added USE_HDF4
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#ifdef USE_HDF4

#include <stdio.h>
#include <math.h>
#include <df.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "fortran.def"
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif
#ifdef USE_FLEXIO
#include "IO.hh"
#include "IEEEIO.hh"
#endif

/* Specify the I/O format for reading (HDF unless it's FlexIO) */

#ifdef WRITEGRID_FLEXIO
#define READFILE ReadFlexIOFile
@@@ FLEXIO WONT WORK: PARAMETER LISTS FOR ReadFlexIOFile and ReadFileHDF
@@@ HAVE DIVERGED!
#else
#define READFILE ReadHDF4File
#endif

/* function prototypes */

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
int ReadHDF4File(char *name, int Rank, int Dims[], int StartIndex[],
		int EndIndex[], int BufferOffset[], float *buffer,
		float32 **tempbuffer);
int ReadFlexIOFile(char *name, int Rank, int Dims[], int StartIndex[],
		   int EndIndex[], int BufferOffset[], float *buffer);
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

void pcol32(float32 *x, int n, int m, FILE *log_fptr);
void fcol(float *x, int n, int m, FILE *log_fptr);
void pcol(FLOAT *x, int n, int m, FILE *log_fptr);
void icol(int *x, int n, int m, FILE *log_fptr);


int grid::CosmologySimulationInitializeGridHDF4(
			  float CosmologySimulationOmegaBaryonNow,
			  float CosmologySimulationOmegaCDMNow,
			  float CosmologySimulationInitialTemperature,
			  char *CosmologySimulationDensityName,
			  char *CosmologySimulationTotalEnergyName,
			  char *CosmologySimulationGasEnergyName,
			  char *CosmologySimulationVelocityNames[],
			  char *CosmologySimulationParticlePositionName,
			  char *CosmologySimulationParticleVelocityName,
			  char *CosmologySimulationParticleMassName,
			  int   CosmologySimulationSubgridsAreStatic,
			  int   TotalRefinement,
			  float CosmologySimulationInitialFractionHII,
			  float CosmologySimulationInitialFractionHeII,
			  float CosmologySimulationInitialFractionHeIII,
			  float CosmologySimulationInitialFractionHM,
			  float CosmologySimulationInitialFractionH2I,
			  float CosmologySimulationInitialFractionH2II,
			  int   UseMetallicityField,
			  int  &CurrentParticleNumber)
{
  /* declarations */

  int  idim, dim, i, j, vel;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, MetalNum;

  int ExtraField[2];

  float32 *tempbuffer = NULL;

  FILE *log_fptr;

  char pid[5];
  sprintf(pid, "%4.4d", MyProcessorNumber);

  char *logname = new char[strlen(pid)+6+1];
  strcpy(logname, "CSlog.");
  strcat(logname,pid);
  log_fptr = fopen(logname, "a");

  fprintf(log_fptr,"\n");
  fprintf(log_fptr,"CSIG ParallelRootGridIO = %d\n",ParallelRootGridIO);
  fprintf(log_fptr,"Processor %d, Target processor %d\n",MyProcessorNumber,ProcessorNumber);
  fprintf(log_fptr,"TotalRefinement = %d\n",TotalRefinement);

  /* Determine if the data should be loaded in or not. */

  int ReadData = TRUE, Offset[] = {0,0,0};

  if (ParallelRootGridIO == TRUE && TotalRefinement == 1)
    ReadData = FALSE;

  fprintf(log_fptr,"ReadData = %d\n",ReadData);

  /* Calculat buffer Offset (same as Grid unless doing ParallelRootGridIO 
     (TotalRefinement = -1 if used as a signal that we should really load 
     in the data regardless of the value of ParallelRootGridIO). */

  if (ParallelRootGridIO == TRUE && TotalRefinement == -1)
    for (dim = 0; dim < GridRank; dim++)
      Offset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/
			 CellWidth[dim][0]);

  /*----------------------------------------------------*/
  /* Create baryon fields (unless they are not needed). */

  NumberOfBaryonFields = 0;
  if (CosmologySimulationDensityName != NULL) {
    FieldType[NumberOfBaryonFields++] = Density;
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;
    FieldType[NumberOfBaryonFields++] = Velocity1;
    vel = NumberOfBaryonFields - 1;
    if (GridRank > 1)
      FieldType[NumberOfBaryonFields++] = Velocity2;
    if (GridRank > 2)
      FieldType[NumberOfBaryonFields++] = Velocity3;
    if (MultiSpecies) {
      FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
      FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
      FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
      FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
      FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
      FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
      if (MultiSpecies > 1) {
	FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
	FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
	FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
      }
      if (MultiSpecies > 2) {
	FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
	FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
	FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
      }
    }
    if (UseMetallicityField) {
      FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity;
      FieldType[ExtraField[0] = NumberOfBaryonFields++] = ExtraType0;
      FieldType[ExtraField[1] = NumberOfBaryonFields++] = ExtraType1;
    }
  }

  /* Set the subgrid static flag. */

  SubgridsAreStatic = CosmologySimulationSubgridsAreStatic;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber == MyProcessorNumber) {

  /* Skip following if NumberOfBaryonFields == 0. */

  if (NumberOfBaryonFields > 0) {

  /* Get the cosmology units so we can convert temperature later. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
        VelocityUnits;
  if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		&TimeUnits, &VelocityUnits, InitialTimeInCodeUnits) == FAIL) {
    fprintf(stderr, "Error in CosmologyGetUnits.\n");
    return FAIL;
  }

  /* Determine the size of the fields. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Allocate space for the fields. */

  if (ReadData == TRUE)
  {
    fprintf(log_fptr,"Allocate %d fields, %d floats per field\n",NumberOfBaryonFields,size);
    for (dim=0; dim < GridRank; dim++)
    {
      fprintf(log_fptr,"  Field dim %d size %d\n",dim,GridDimension[dim]);
    }
  }

  if (ReadData == TRUE)
    for (int field = 0; field < NumberOfBaryonFields; field++)
      BaryonField[field] = new float[size];

  /* Read the density field. */

  if (CosmologySimulationDensityName != NULL && ReadData)
  {
    if (READFILE(CosmologySimulationDensityName, GridRank, GridDimension,
                 GridStartIndex, GridEndIndex, Offset, BaryonField[0], 
		 &tempbuffer) == FAIL) {
      fprintf(stderr, "Error reading density field.\n");
      return FAIL;
    }
//  fcol(BaryonField[0],size,10,log_fptr);
  }

  /* Read the total energy field. */

  if (CosmologySimulationTotalEnergyName != NULL && ReadData)
    if (READFILE(CosmologySimulationTotalEnergyName, GridRank, 
		    GridDimension, GridStartIndex, GridEndIndex, Offset,
		    BaryonField[1], &tempbuffer) == FAIL) {
      fprintf(stderr, "Error reading total energy field.\n");
      return FAIL;
    }

  /* Read the gas energy field. */

  if (CosmologySimulationGasEnergyName != NULL && DualEnergyFormalism && 
      ReadData)
    if (READFILE(CosmologySimulationGasEnergyName, GridRank, GridDimension,
		 GridStartIndex, GridEndIndex, Offset, BaryonField[2],
		 &tempbuffer) == FAIL) {
      fprintf(stderr, "Error reading gas energy field.\n");
      return FAIL;
    }

  /* Read the velocity fields. */

  if (CosmologySimulationVelocityNames[0] != NULL && ReadData)
  {
    for (dim = 0; dim < GridRank; dim++)
      if (READFILE(CosmologySimulationVelocityNames[dim], GridRank,
		   GridDimension, GridStartIndex, GridEndIndex, Offset,
		   BaryonField[vel+dim], &tempbuffer) == FAIL) {
	fprintf(stderr, "Error reading velocity field %d.\n", dim);
	return FAIL;
      }
//  for (dim = 0; dim < GridRank; dim++)
//    fcol(BaryonField[vel+dim],size,10,log_fptr);
  }

  /* If using multi-species, set the fields. */

  if (MultiSpecies && ReadData)
    for (i = 0; i < size; i++) {

      BaryonField[HIINum][i] = CosmologySimulationInitialFractionHII *
	CoolData.HydrogenFractionByMass * BaryonField[0][i] *
	sqrt(OmegaMatterNow)/
	(CosmologySimulationOmegaBaryonNow*HubbleConstantNow);
      
      BaryonField[HeIINum][i] = CosmologySimulationInitialFractionHeII*
	BaryonField[0][i] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
      BaryonField[HeIIINum][i] = CosmologySimulationInitialFractionHeIII*
	BaryonField[0][i] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
      BaryonField[HeINum][i] = 
	(1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][i] -
	BaryonField[HeIINum][i] - BaryonField[HeIIINum][i];

      if (MultiSpecies > 1) {
	BaryonField[HMNum][i] = CosmologySimulationInitialFractionHM*
	  BaryonField[HIINum][i]*
	  POW(CosmologySimulationInitialTemperature,float(0.88));
	BaryonField[H2IINum][i] = CosmologySimulationInitialFractionH2II*2.0*
	  BaryonField[HIINum][i]*
	  POW(CosmologySimulationInitialTemperature,float(1.8));
	BaryonField[H2INum][i] = CosmologySimulationInitialFractionH2I*
	  BaryonField[0][i]*CoolData.HydrogenFractionByMass*POW(301.0,5.1)*
	  POW(OmegaMatterNow, float(1.5))/
	  CosmologySimulationOmegaBaryonNow/HubbleConstantNow*2.0;
      }

      BaryonField[HINum][i] = CoolData.HydrogenFractionByMass*BaryonField[0][i]
	- BaryonField[HIINum][i];
      if (MultiSpecies > 1)
	BaryonField[HINum][i] -= BaryonField[HMNum][i]
	  + BaryonField[H2IINum][i]
	  + BaryonField[H2INum][i];

      BaryonField[DeNum][i] = BaryonField[HIINum][i] + 
	0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];
      if (MultiSpecies > 1)
	BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] - 
	  BaryonField[HMNum][i];

      /* Set Deuterium species (assumed to be negligible). */

      if (MultiSpecies > 2) {
	BaryonField[DINum][i] = CoolData.DeuteriumToHydrogenRatio*
	                             BaryonField[HINum][i];
	BaryonField[DIINum][i] = CoolData.DeuteriumToHydrogenRatio*
	                             BaryonField[HIINum][i];
	BaryonField[HDINum][i] = CoolData.DeuteriumToHydrogenRatio*
	                             BaryonField[H2INum][i];
      }

    }
 
  /* If using metallicity, set the field. */
 
  if (UseMetallicityField && ReadData)
    for (i = 0; i < size; i++) {
      BaryonField[MetalNum][i] = 1.0e-10 * BaryonField[0][i];
      BaryonField[ExtraField[0]][i] = 1.0e-10 * BaryonField[0][i];
      BaryonField[ExtraField[1]][i] = 1.0e-10 * BaryonField[0][i];
    }

  /* If they were not read in above, set the total & gas energy fields now. */

  if (CosmologySimulationDensityName != NULL && ReadData) {
    if (CosmologySimulationTotalEnergyName == NULL)
      for (i = 0; i < size; i++)
	BaryonField[1][i] = CosmologySimulationInitialTemperature/
	                      TemperatureUnits;
/*          * POW(BaryonField[0][i]/CosmologySimulationOmegaBaryonNow,Gamma-1) 
	                    / (Gamma-1); */

    if (CosmologySimulationGasEnergyName == NULL && DualEnergyFormalism)
      for (i = 0; i < size; i++)
	BaryonField[2][i] = BaryonField[1][i];

    if (CosmologySimulationTotalEnergyName == NULL &&
	HydroMethod != Zeus_Hydro)
      for (dim = 0; dim < GridRank; dim++)
	for (i = 0; i < size; i++)
	  BaryonField[1][i] += 
	    0.5 * BaryonField[vel+dim][i] * BaryonField[vel+dim][i];
  }

  } // end: if (NumberOfBaryonFields > 0)

  /*----------------------------------------------------*/
  /* Read particle fields if the was a name specifieid. */

  if (CosmologySimulationParticlePositionName != NULL && ReadData) {

    /* Get the number of particles by quizing the file. */

    int TempInt, Dim[1], Start[1] = {0}, End[1], Zero[1] = {0};
#ifdef WRITEGRID_FLEXIO
    int TempIntArray[MAX_DIMENSION];
    IObase::DataType InputDataType;
    IObase *FlexIOreader = new IEEEIO(CosmologySimulationParticlePositionName,
				      IObase::Read);
    if (!FlexIOreader->isValid()) {
      fprintf(stderr, "Error opening FlexIO file %s.\n",
	      CosmologySimulationParticlePositionName);
      return FAIL;
    }
    FlexIOreader->readInfo(InputDataType, TempInt, TempIntArray);
    delete FlexIOreader;
#else
    int32 TempIntArray[MAX_DIMENSION];
    if (DFSDgetdims(CosmologySimulationParticlePositionName, 
		    &TempInt, TempIntArray, MAX_DIMENSION) == HDF_FAIL) {
      fprintf(stderr, "Error reading dims from %s.\n", 
	      CosmologySimulationParticlePositionName);
      return FAIL;
    }
    DFSDrestart();
#endif

    /* Error check */

    if (TempInt != 1) {
      fprintf(stderr, "Rank (%d) is not one in file %s.\n", TempInt,
	      CosmologySimulationParticlePositionName);
      return FAIL;
    }

    /* If doing parallel root grid IO then read in the full list of particle
       positions and construct a mask of those within this grid. */

    int *Mask = NULL, TotalParticleCount = TempIntArray[0];
    if (ParallelRootGridIO == TRUE && TotalRefinement == -1) {

      /* Generate a mask indicating if the particle is inside or outside. */

      Mask = new int[TotalParticleCount];
      for (i = 0; i < TotalParticleCount; i++)
	Mask[i] = TRUE;

      float32 *TempField = new float32[TotalParticleCount];
      End[0] = TotalParticleCount - 1;
      Dim[0] = TotalParticleCount;
      for (dim = 0; dim < GridRank; dim++) {
	//	if (READFILE(CosmologySimulationParticlePositionName, 1, Dim,
	//		     Start, End, Zero, TempField) == FAIL) {
	DFSDgetdims(CosmologySimulationParticlePositionName, &TempInt,
		    TempIntArray, MAX_DIMENSION);
	if (DFSDgetdata(CosmologySimulationParticlePositionName, TempInt,
			TempIntArray, (VOIDP) TempField) == HDF_FAIL) {
	  fprintf(stderr, "Error reading particle position %d.\n", dim);
	  return FAIL;
	}

	for (i = 0; i < TotalParticleCount; i++)
	  if (TempField[i] < GridLeftEdge[dim] ||
	      TempField[i] > GridRightEdge[dim] )
	    Mask[i] = FALSE;
      }

      /* Compute the number of particles left. */

      NumberOfParticles = 0;
      for (i = 0; i < TotalParticleCount; i++)
	NumberOfParticles += Mask[i];

      fprintf(log_fptr,"NumberOfParticles %d\n",NumberOfParticles);

      /* Allocate space for the particles. */

      this->AllocateNewParticles(NumberOfParticles);

      /* Read particle positions using mask and delete temp positions. */

      int n;
#ifdef WRITEGRID_FLEXIO
      return FAIL;
#else
      DFSDrestart();
#endif
      for (dim = 0; dim < GridRank; dim++) {
	  //	if (READFILE(CosmologySimulationParticlePositionName, 1, Dim,
	  //		     Start, End, Zero, TempField) == FAIL) {
	DFSDgetdims(CosmologySimulationParticlePositionName, &TempInt,
		    TempIntArray, MAX_DIMENSION);
	if (DFSDgetdata(CosmologySimulationParticlePositionName, TempInt,
			TempIntArray, (VOIDP) TempField) == HDF_FAIL) {
	  fprintf(stderr, "Error reading particle position %d.\n", dim);
	  return FAIL;
	}
	for (n = 0, i = 0; i < TotalParticleCount; i++)
	  if (Mask[i])
	    ParticlePosition[dim][n++] = TempField[i];
      }

      /* Read particle velocities using mask. */

      if (CosmologySimulationParticleVelocityName != NULL)
	for (dim = 0; dim < GridRank; dim++) {
	  //	  if (READFILE(CosmologySimulationParticleVelocityName, 1, Dim,
	  //		       Start, End, Zero, TempField) == FAIL) {
	DFSDgetdims(CosmologySimulationParticleVelocityName, &TempInt,
		    TempIntArray, MAX_DIMENSION);
	if (DFSDgetdata(CosmologySimulationParticleVelocityName, TempInt,
			TempIntArray, (VOIDP) TempField) == HDF_FAIL) {
	    fprintf(stderr, "Error reading particle velocity %d.\n", dim);
	    return FAIL;
	  }
	  for (n = 0, i = 0; i < TotalParticleCount; i++)
	    if (Mask[i])
	      ParticleVelocity[dim][n++] = TempField[i];
	}
      
      /* Read particle mass. */

      if (CosmologySimulationParticleMassName != NULL) {
	//if (READFILE(CosmologySimulationParticleMassName, 1, Dim, Start, End,
	//		     Zero, TempField) == FAIL) {
	DFSDgetdims(CosmologySimulationParticleMassName, &TempInt,
		    TempIntArray, MAX_DIMENSION);
	if (DFSDgetdata(CosmologySimulationParticleMassName, TempInt,
			TempIntArray, (VOIDP) TempField) == HDF_FAIL) {
	  fprintf(stderr, "Error reading particle masses.\n");
	  return FAIL;
	}
	for (n = 0, i = 0; i < TotalParticleCount; i++)
	  if (Mask[i])
	    ParticleMass[n++] = TempField[i];

      }  /* End of loop over particle mass dimension */

/*
      for (idim = 0; idim < GridRank; idim++)
      {
        pcol(ParticlePosition[idim],NumberOfParticles,10,log_fptr);
      }
      for (idim = 0; idim < GridRank; idim++)
      {
        fcol(ParticleVelocity[idim],NumberOfParticles,10,log_fptr);
      }
      fcol(ParticleMass,NumberOfParticles,10,log_fptr);
*/

      /* delete mask and temp field. */

      fprintf(log_fptr,"De-Allocate Mask and TempField\n");

      delete [] Mask;
      delete [] TempField;
  
    } else {

      /* For regular IO (non-parallel) */

      fprintf(log_fptr,"CSIG Serial I/O\n");

      /* Set Number Of particles. */

      NumberOfParticles = TempIntArray[0];

      End[0] = NumberOfParticles - 1;
      Dim[0] = NumberOfParticles;

      /* Allocate space for the particles. */

      this->AllocateNewParticles(NumberOfParticles);
  
      /* Read particle positions. */

      for (dim = 0; dim < GridRank; dim++) {
	if (READFILE(CosmologySimulationParticlePositionName, 1, Dim,
		     Start, End, Zero, NULL, &tempbuffer) == FAIL) {
	  fprintf(stderr, "Error reading particle position %d.\n", dim);
	  return FAIL;
	}
	for (i = Start[0]; i <= End[0]; i++)
	  ParticlePosition[dim][i] = FLOAT(tempbuffer[i]);

//      pcol32(tempbuffer,NumberOfParticles,10,log_fptr);

        fprintf(log_fptr,"De-Allocate [] tempbuffer for dim = %d\n",dim);

	delete [] tempbuffer;
      }

      /* Read particle velocities. */

      if (CosmologySimulationParticleVelocityName != NULL)
	for (dim = 0; dim < GridRank; dim++) {
	  if (READFILE(CosmologySimulationParticleVelocityName, 1, Dim,
	     Start, End, Zero, ParticleVelocity[dim], &tempbuffer) == FAIL) {
	    fprintf(stderr, "Error reading particle velocity %d.\n", dim);
	    return FAIL;
	  }
//        fcol(ParticleVelocity[dim],NumberOfParticles,10,log_fptr);
	}
  
      /* Read particle mass. */

      if (CosmologySimulationParticleMassName != NULL) {
	if (READFILE(CosmologySimulationParticleMassName, 1, Dim, Start, End,
		     Zero, ParticleMass, &tempbuffer) == FAIL) {
	  fprintf(stderr, "Error reading particle masses.\n");
	  return FAIL;
	}
//      fcol(ParticleMass,NumberOfParticles,10,log_fptr);
      }

    } // end: if ParallelRootGridIO == TRUE

    /* If there are particles, but no velocity file, then set the velocities
       to zero. */
      
    if (NumberOfParticles > 0 && 
	CosmologySimulationParticleVelocityName == NULL) {
      fprintf(stderr, "Error -- no velocity field specified.\n");
      return FAIL;
      //  printf("CosmologySimulation warning: setting velocities to zero.\n");
      //      for (dim = 0; dim < GridRank; dim++)
      //	for (i = 0; i < NumberOfParticles; i++)
      //	  ParticleVelocity[dim][i] = 0.0;
    }

    /* If there are particles, but a mass file name wasn't specified, make
       up the particles mass. */

    if (NumberOfParticles > 0 && CosmologySimulationParticleMassName == NULL) {

      /* Compute mass of each particle. */

      float UniformParticleMass = CosmologySimulationOmegaCDMNow/
	                          OmegaMatterNow;

      /* If there are exactly 1/8 as many particles as cells, then set the
	 particle mass to 8 times the usual.  */

      int NumberOfActiveCells = (GridEndIndex[0]-GridStartIndex[0]+1)*
	                        (GridEndIndex[1]-GridStartIndex[1]+1)*
	                        (GridEndIndex[2]-GridStartIndex[2]+1);
      if (NumberOfParticles*8 == NumberOfActiveCells)
	UniformParticleMass *= 8;
      if (NumberOfParticles == NumberOfActiveCells*8)
        UniformParticleMass /= 8;

      //      UniformParticleMass *= float(POW(TotalRefinement, GridRank));
	/*      for (dim = 0; dim < GridRank; dim++)
	UniformParticleMass *= float(GridEndIndex[dim]-GridStartIndex[dim]+1);
      UniformParticleMass /= float(NumberOfParticles); */

      /* Set all particle masses to this value. */

      if (debug) printf("CosmologySimulationInitializeGrid: particle mass = %g\n",
			UniformParticleMass);
      for (i = 0; i < NumberOfParticles; i++)
	ParticleMass[i] = UniformParticleMass;
    
    }

    /* Set Particle attributes to FLOAT_UNDEFINED. */

    for (j = 0; j < NumberOfParticleAttributes; j++)
      for (i = 0; i < NumberOfParticles; i++)
	ParticleAttribute[j][i] = FLOAT_UNDEFINED;

    /* Assign particles unique identifiers. */

    if (ParallelRootGridIO == TRUE)
      CurrentParticleNumber = 0;  /* unique ID's will be added later */
    for (i = 0; i < NumberOfParticles; i++)
      ParticleNumber[i] = CurrentParticleNumber + i;
    CurrentParticleNumber += NumberOfParticles;

  } // end: if (CosmologySimulationParticleName != NULL)

  } // end: if (ProcessorNumber == MyProcessorNumber)

  /* Share number of particles amoung processors. */

  if (ReadData == TRUE)
    CommunicationBroadcastValue(&NumberOfParticles, ProcessorNumber);

  OldTime = Time;

  fclose(log_fptr);
  
  return SUCCESS;
}

#else 

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "message.h"

// HDF4 is not used, so CosmologySimulationInitializeGridHDF4 should
//   not be called

int grid::CosmologySimulationInitializeGridHDF4(
			  float CosmologySimulationOmegaBaryonNow,
			  float CosmologySimulationOmegaCDMNow,
			  float CosmologySimulationInitialTemperature,
			  char *CosmologySimulationDensityName,
			  char *CosmologySimulationTotalEnergyName,
			  char *CosmologySimulationGasEnergyName,
			  char *CosmologySimulationVelocityNames[],
			  char *CosmologySimulationParticlePositionName,
			  char *CosmologySimulationParticleVelocityName,
			  char *CosmologySimulationParticleMassName,
			  int   CosmologySimulationSubgridsAreStatic,
			  int   TotalRefinement,
			  float CosmologySimulationInitialFractionHII,
			  float CosmologySimulationInitialFractionHeII,
			  float CosmologySimulationInitialFractionHeIII,
			  float CosmologySimulationInitialFractionHM,
			  float CosmologySimulationInitialFractionH2I,
			  float CosmologySimulationInitialFractionH2II,
			  int   UseMetallicityField,
			  int  &CurrentParticleNumber)
{ 
  ERROR_MESSAGE;
  return FAIL;
}
#endif
