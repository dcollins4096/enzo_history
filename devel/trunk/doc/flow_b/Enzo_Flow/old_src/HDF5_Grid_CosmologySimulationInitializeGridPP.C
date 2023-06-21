/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COSMOLOGY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness, July 2002
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "hdf4.h"  /* HDF4 int32, float32 and VOIDP definitions */

#include "mpi.h"

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

#define READFILE HDF5_ReadFile

// HDF5 function prototypes

#include "extern_hdf5.h"

// function prototypes

int HDF5_ReadFile(char *name, int Rank, int Dims[], int StartIndex[],
                  int EndIndex[], int BufferOffset[], float *buffer,
                  float32 **tempbuffer, int Part, int Npart);

int HDF5_ReadAttr(char *Fname, int *Rank, int Dims[], int *NSeg, int *LSeg, FILE *log_fptr);

void pcol32(float32 *x, int n, int m, FILE *log_fptr);
void fcol(float *x, int n, int m, FILE *log_fptr);
void pcol(FLOAT *x, int n, int m, FILE *log_fptr);
void icol(int *x, int n, int m, FILE *log_fptr);

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);




int grid::CosmologySimulationInitializeGrid(
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

  char *subroutine="CosmologySimulationInitializeGrid";

  int idim, dim, i, j, vel;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, MetalNum;

  float32 *tempbuffer = NULL;

  hid_t       file_id, dset_id, attr_id, type_id;
  hid_t       mem_dsp_id, file_dsp_id, attr_dsp_id;

  herr_t      h5_status;
  herr_t      h5_error = -1;

  hsize_t     xfer_size;
  hsize_t     Slab_Dims[4];
  int         Slab_Rank;

  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     slab_stride[4], slab_count[4], slab_block[4];
  hsize_t     attr_count;

  hssize_t    mem_offset;
  hssize_t    slab_offset[4];
  
  int NSeg, LSeg, Part;

  int component_rank_attr;
  int component_size_attr;
  int field_rank_attr;
  int field_dims_attr[3];

  FILE *log_fptr;

#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

//  printf("Enter %s on processor %d\n", subroutine, MyProcessorNumber);

  char pid[5];
  sprintf(pid, "%4.4d", MyProcessorNumber);

  char *logname = new char[strlen(pid)+6+1];
  strcpy(logname, "CSlog.");
  strcat(logname,pid);

  if (io_log) log_fptr = fopen(logname, "a");

  if (io_log) fprintf(log_fptr, "\n");
  if (io_log) fprintf(log_fptr, "CSIG ParallelRootGridIO = %d\n", ParallelRootGridIO);
  if (io_log) fprintf(log_fptr, "Processor %d, Target processor %d\n", MyProcessorNumber, ProcessorNumber);
  if (io_log) fprintf(log_fptr, "TotalRefinement = %d\n", TotalRefinement);

  /* Determine if the data should be loaded in or not. */

  int ReadData = TRUE, Offset[] = {0,0,0};

  if (ParallelRootGridIO == TRUE && TotalRefinement == 1)
    ReadData = FALSE;

  if (io_log) fprintf(log_fptr, "ReadData = %d\n", ReadData);

  /* Calculate buffer Offset (same as Grid unless doing ParallelRootGridIO 
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
    if (UseMetallicityField)
      FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity;
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
    if (io_log) fprintf(log_fptr, "Allocate %d fields, %d floats per field\n", NumberOfBaryonFields, size);
    for (dim=0; dim < GridRank; dim++)
    {
      if (io_log) fprintf(log_fptr, "  Field dim %d size %d\n", dim, GridDimension[dim]);
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
              &tempbuffer, 0, 1) == FAIL) {
      fprintf(stderr, "Error reading density field.\n");
      return FAIL;
    }
//  fcol(BaryonField[0], size, 10, log_fptr);
  }

  /* Read the total energy field. */

  if (CosmologySimulationTotalEnergyName != NULL && ReadData)
    if (READFILE(CosmologySimulationTotalEnergyName, GridRank, 
		    GridDimension, GridStartIndex, GridEndIndex, Offset,
		    BaryonField[1], &tempbuffer, 0, 1) == FAIL) {
      fprintf(stderr, "Error reading total energy field.\n");
      return FAIL;
    }

  /* Read the gas energy field. */

  if (CosmologySimulationGasEnergyName != NULL && DualEnergyFormalism && 
      ReadData)
    if (READFILE(CosmologySimulationGasEnergyName, GridRank, GridDimension,
		 GridStartIndex, GridEndIndex, Offset, BaryonField[2],
		 &tempbuffer, 0, 1) == FAIL) {
      fprintf(stderr, "Error reading gas energy field.\n");
      return FAIL;
    }

  /* Read the velocity fields. */

  if (CosmologySimulationVelocityNames[0] != NULL && ReadData)
  {
    for (dim = 0; dim < GridRank; dim++)
      if (READFILE(CosmologySimulationVelocityNames[dim], GridRank,
		   GridDimension, GridStartIndex, GridEndIndex, Offset,
		   BaryonField[vel+dim], &tempbuffer, dim, 3) == FAIL) {
	fprintf(stderr, "Error reading velocity field %d.\n", dim);
	return FAIL;
      }
//  for (dim = 0; dim < GridRank; dim++)
//    fcol(BaryonField[vel+dim], size, 10, log_fptr);
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
    for (i = 0; i < size; i++)
      BaryonField[MetalNum][i] = 1.0e-10;

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

    int TempIntArray[MAX_DIMENSION];

    HDF5_ReadAttr(CosmologySimulationParticlePositionName,
                  &TempInt, TempIntArray, &NSeg, &LSeg, log_fptr);

    /* Error check */

    if (TempInt != 1) {
      fprintf(stderr, "Rank (%d) is not one in file %s.\n", TempInt,
	      CosmologySimulationParticlePositionName);
      return FAIL;
    }

    /* If doing parallel root grid IO then read in the full list of particle
       positions and construct a mask of those within this grid. */

//    int *Mask = NULL;
//    int *XMask = NULL;

    int TotalParticleCount = TempIntArray[0];

    unsigned int *BitMask = NULL;
    unsigned int BitsPerInt;
    unsigned int BitMaskSize;
    unsigned int BitMaskTrue;
    unsigned int *BitArray = NULL;
    unsigned int TestBit;
    int XMask;

    BitsPerInt = 8 * sizeof(BitsPerInt);
    BitMaskSize = (TotalParticleCount/BitsPerInt)+1;

    if (ParallelRootGridIO == TRUE && TotalRefinement == -1) {

      for(i=0; i<GridRank;i++)
      {
         fprintf(log_fptr,"Proc %d Left %10.4f Right %10.4f\n",MyProcessorNumber,GridLeftEdge[i],GridRightEdge[i]);
      }

      if (io_log) fprintf(log_fptr, "CSIG ParallelRootGridIO\n");
      if (io_log) fprintf(log_fptr, "TotalParticleCount %d\n", TotalParticleCount);

      /* Generate a mask indicating if the particle is inside or outside. */

      if (io_log) fprintf(log_fptr, "Allocate %d ints for Mask\n", TotalParticleCount);

//      Mask = new int[TotalParticleCount];
//      XMask = new int[TotalParticleCount];

      BitMask = new unsigned int[BitMaskSize];

      BitArray = new unsigned int[BitsPerInt];

//      for (i = 0; i < TotalParticleCount; i++)
//	Mask[i] = TRUE;
//      for (i = 0; i < TotalParticleCount; i++)
//        XMask[i] = TRUE;

      for (i = 0; i < BitMaskSize; i++)
        BitMask[i] = 0;

      BitArray[BitsPerInt-1] = 1;
      BitMaskTrue = 1;
      for (i=BitsPerInt-1; i>0; i--)
      {
        BitArray[i-1] = 2 * BitArray[i];
        BitMaskTrue = (BitMaskTrue | BitArray[i-1]);
      }

//      fprintf(log_fptr, "All true %16o\n", BitMaskTrue);

      for (i = 0; i < BitMaskSize; i++)
        BitMask[i] = BitMaskTrue;

      unsigned int MaskAddr;
      unsigned int WordAddr;
      unsigned int BitAddr;

      if (io_log) fprintf(log_fptr, "Allocate %d float32s for TempField\n", TotalParticleCount);

      float32 *TempField = new float32[TotalParticleCount];

      End[0] = TotalParticleCount - 1;
      Dim[0] = TotalParticleCount;

      int ii = sizeof(float32);

      switch(ii)
      {

        case 4:
          type_id = HDF5_R4;
          break;

        case 8:
          type_id = HDF5_R8;
          break;

        default:
          type_id = HDF5_R4;

      }

      int ppbuffsize, pp1, pp2, ppoffset, ppcount;

      ppbuffsize = 512*512;
//    printf("ppbuffsize = %d\n", ppbuffsize);

//    if (io_log) fprintf(log_fptr, "Allocate %d float32s for TempField\n", ppbuffsize);
//    float32 *TempField = new float32[ppbuffsize];

      for (pp1 = 0; pp1 < TotalParticleCount; pp1 = pp1 + ppbuffsize)
      {
         pp2 = min(pp1+ppbuffsize-1,TotalParticleCount-1);
         ppoffset = pp1;
         ppcount = pp2-pp1+1;
//       printf("PP1 = %d, PP2 = %d, PPoffset = %d, PPCount = %d\n", pp1, pp2, ppoffset, ppcount);
      }

      

      for (dim = 0; dim < GridRank; dim++) {

        Part = dim;

        if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", CosmologySimulationParticlePositionName);

        file_id = H5Fopen(CosmologySimulationParticlePositionName, H5F_ACC_RDONLY, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Fopen id: %d\n", file_id);
          assert( file_id != h5_error );

        if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", CosmologySimulationParticlePositionName);

        dset_id = H5Dopen(file_id, CosmologySimulationParticlePositionName);
          if (io_log) fprintf(log_fptr, "H5Dopen id: %d\n", dset_id);
          assert( dset_id != h5_error );


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Rank\n");

        attr_id = H5Aopen_name(dset_id, "Component_Rank");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, &component_rank_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        if (io_log) fprintf(log_fptr, "COMPONENT_RANK %d\n", component_rank_attr);


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Size\n");

        attr_id = H5Aopen_name(dset_id, "Component_Size");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, &component_size_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        if (io_log) fprintf(log_fptr, "COMPONENT_SIZE %d\n", component_size_attr);


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Rank\n");

        attr_id = H5Aopen_name(dset_id, "Rank");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, &field_rank_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        if (io_log) fprintf(log_fptr, "RANK %d\n", field_rank_attr);


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Dimensions\n");

        attr_id = H5Aopen_name(dset_id, "Dimensions");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        attr_count = field_rank_attr;

        attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
          if (io_log) fprintf(log_fptr, "Attr_dsp_id: %d\n", attr_dsp_id);
          assert( attr_dsp_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, field_dims_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(attr_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        for (idim = 0; idim < field_rank_attr; idim++)
        {
          if (io_log) fprintf(log_fptr, "DIMS %d:  %d\n", idim, (int) field_dims_attr[idim]);
        }

        xfer_size = 1;

        for ( idim = 0; idim < field_rank_attr; idim++ )
        {
          xfer_size = xfer_size * field_dims_attr[idim];
        }

        if (io_log) fprintf(log_fptr, "  Grid Elements %d\n", (int) xfer_size);

        Slab_Rank = field_rank_attr+1;

        Slab_Dims[0] = component_rank_attr;

        for ( idim = 1; idim < Slab_Rank; idim++ )
        {
          Slab_Dims[idim] = field_dims_attr[idim-1];
        }

        if (io_log) fprintf(log_fptr, "  Extended Rank %d\n", (int) Slab_Rank);
        for ( idim = 0; idim < Slab_Rank; idim++ )
        {
          if (io_log) fprintf(log_fptr, "    %d:  %d\n", idim, (int) Slab_Dims[idim]);
        }

        slab_offset[0] = Part;
        slab_stride[0] = 1;
        slab_count[0] = 1;
        slab_block[0] = 1;

        xfer_size = 1;

        for ( idim=1; idim < Slab_Rank; idim++)
        {
          slab_offset[idim] = 0;
          slab_stride[idim] = 1;
          slab_count[idim] = field_dims_attr[idim-1];
          slab_block[idim] = 1;
          xfer_size = xfer_size * slab_count[idim];
        }

        for ( idim=0; idim < Slab_Rank; idim++)
        {
          if (io_log) fprintf(log_fptr, "Slab [%d] Offset %d, Count %d, Stride %d, Block %d\n",
                  idim, (int) slab_offset[idim], (int) slab_count[idim],
                       (int) slab_stride[idim], (int) slab_block[idim]);
        }
        if (io_log) fprintf(log_fptr, "Xfer_size %d\n", (int) xfer_size);
        if (io_log) fprintf(log_fptr, "Comp_size %d\n", (int) component_size_attr);

//        xfer_size = component_size_attr;

// Data in memory is considered 1D, stride 1, with zero offset

        mem_stride = 1;
        mem_count = xfer_size;
        mem_offset = 0;
        mem_block = 1;

// 1D memory model

        mem_dsp_id = H5Screate_simple(1, &xfer_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %d\n", mem_dsp_id);
          assert( mem_dsp_id != h5_error );

        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, &mem_block);
          if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %d\n", h5_status);
          assert( h5_status != h5_error );

        file_dsp_id = H5Screate_simple(Slab_Rank, Slab_Dims, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %d\n", file_dsp_id);
          assert( file_dsp_id != h5_error );

        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
          if (io_log) fprintf(log_fptr, "H5Sselect file slab: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) TempField);
          if (io_log) fprintf(log_fptr, "H5Dread: %d\n", (int) h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(mem_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(file_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Dclose(dset_id);
          if (io_log) fprintf(log_fptr, "H5Dclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Fclose(file_id);
          if (io_log) fprintf(log_fptr, "H5Fclose: %d\n", h5_status);
          assert( h5_status != h5_error );

//      j = 0;
//      for (i = pp1; i < pp1 + ppcount; i++)
//      j = j+1;
//      TempField[j]

	for (i = 0; i < TotalParticleCount; i++)
        {

//          if (TempField[i] < GridLeftEdge[dim] || TempField[i] > GridRightEdge[dim] )
//          {
//	    Mask[i] = FALSE;
//          }

          if (TempField[i] < GridLeftEdge[dim] || TempField[i] > GridRightEdge[dim] )
          { 
            MaskAddr = i;
            WordAddr = MaskAddr/BitsPerInt;
            BitAddr  = MaskAddr%BitsPerInt;
            TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
            if ( TestBit != 0 )
              BitMask[WordAddr] = (BitMask[WordAddr] ^ BitArray[BitAddr]);
//          fprintf(log_fptr, "%4d  %4u  %4u  %16o %16o\n", i, WordAddr, BitAddr, BitMask[WordAddr], BitArray[BitAddr]);
          }
        }

      }  /* End of loop over particle position dimensions */

      for (i=0; i<BitMaskSize; i++)
        fprintf(log_fptr, "%4d  %16o\n", i, BitMask[i]);

//    icol(Mask, TotalParticleCount, 128, log_fptr);

      NumberOfParticles = 0;

      for (i = 0; i < TotalParticleCount; i++)
      {
        XMask = 1;
        MaskAddr = i;
        WordAddr = MaskAddr/BitsPerInt;
        BitAddr  = MaskAddr%BitsPerInt;
        TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
        if ( TestBit == 0 )
          XMask = 0;
        NumberOfParticles += XMask;
//      fprintf(log_fptr, "%4d    %2d %2d    %16o  %16o  %16o\n", i, Mask[i], XMask[i], BitMask[WordAddr], BitArray[BitAddr], TestBit); 
      }

//    icol(XMask, TotalParticleCount, 128, log_fptr);

      /* Compute the number of particles left. */

//      NumberOfParticles = 0;
//      for (i = 0; i < TotalParticleCount; i++)
//	NumberOfParticles += Mask[i];

      if (io_log) fprintf(log_fptr, "NumberOfParticles %d\n", NumberOfParticles);

      /* Allocate space for the particles. */

      this->AllocateNewParticles(NumberOfParticles);

      /* Read particle positions using mask and delete temp positions. */

      int n;

      for (dim = 0; dim < GridRank; dim++) {

        Part = dim;

        if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", CosmologySimulationParticlePositionName);

        file_id = H5Fopen(CosmologySimulationParticlePositionName, H5F_ACC_RDONLY, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Fopen id: %d\n", file_id);
          assert( file_id != h5_error );

        if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", CosmologySimulationParticlePositionName);

        dset_id = H5Dopen(file_id, CosmologySimulationParticlePositionName);
          if (io_log) fprintf(log_fptr, "H5Dopen id: %d\n", dset_id);
          assert( dset_id != h5_error );

        slab_offset[0] = Part;
        slab_stride[0] = 1;
        slab_count[0] = 1;
        slab_block[0] = 1;

        xfer_size = 1;

        for ( idim=1; idim < Slab_Rank; idim++)
        {
          slab_offset[idim] = 0;
          slab_stride[idim] = 1;
          slab_count[idim] = field_dims_attr[idim-1];
          slab_block[idim] = 1;
          xfer_size = xfer_size * slab_count[idim];
        }

        for ( idim=0; idim < Slab_Rank; idim++)
        {
          if (io_log) fprintf(log_fptr, "Slab [%d] Offset %d, Count %d, Stride %d, Block %d\n",
                  idim, (int) slab_offset[idim], (int) slab_count[idim],
                       (int) slab_stride[idim], (int) slab_block[idim]);
        }
        if (io_log) fprintf(log_fptr, "Xfer_size %d\n", (int) xfer_size);
        if (io_log) fprintf(log_fptr, "Comp_size %d\n", (int) component_size_attr);

// Data in memory is considered 1D, stride 1, with zero offset

        mem_stride = 1;
        mem_count = xfer_size;
        mem_offset = 0;
        mem_block = 1;

// 1D memory model

        mem_dsp_id = H5Screate_simple(1, &xfer_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %d\n", mem_dsp_id);
          assert( mem_dsp_id != h5_error );

        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, &mem_block);
          if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %d\n", h5_status);
          assert( h5_status != h5_error );

        file_dsp_id = H5Screate_simple(Slab_Rank, Slab_Dims, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %d\n", file_dsp_id);
          assert( file_dsp_id != h5_error );

        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
          if (io_log) fprintf(log_fptr, "H5Sselect file slab: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) TempField);
          if (io_log) fprintf(log_fptr, "H5Dread: %d\n", (int) h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(mem_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(file_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Dclose(dset_id);
          if (io_log) fprintf(log_fptr, "H5Dclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Fclose(file_id);
          if (io_log) fprintf(log_fptr, "H5Fclose: %d\n", h5_status);
          assert( h5_status != h5_error );

	for (n = 0, i = 0; i < TotalParticleCount; i++)
        {
          XMask = TRUE;
          MaskAddr = i;
          WordAddr = MaskAddr/BitsPerInt;
          BitAddr  = MaskAddr%BitsPerInt;
          TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);

          if ( TestBit == 0 )
            XMask = FALSE;

	  if (XMask)
	    ParticlePosition[dim][n++] = TempField[i];
        }

      }  /* End of loop over particle position dimensions */

      /* Read particle velocities using mask. */

      if (CosmologySimulationParticleVelocityName != NULL)
	for (dim = 0; dim < GridRank; dim++) {

        Part = dim;

        if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", CosmologySimulationParticleVelocityName);

        file_id = H5Fopen(CosmologySimulationParticleVelocityName, H5F_ACC_RDONLY, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Fopen id: %d\n", file_id);
          assert( file_id != h5_error );

        if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", CosmologySimulationParticleVelocityName);

        dset_id = H5Dopen(file_id, CosmologySimulationParticleVelocityName);
          if (io_log) fprintf(log_fptr, "H5Dopen id: %d\n", dset_id);
          assert( dset_id != h5_error );


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Rank\n");

        attr_id = H5Aopen_name(dset_id, "Component_Rank");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, &component_rank_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr," H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        if (io_log) fprintf(log_fptr, "COMPONENT_RANK %d\n", component_rank_attr);


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Size\n");

        attr_id = H5Aopen_name(dset_id, "Component_Size");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, &component_size_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        if (io_log) fprintf(log_fptr, "COMPONENT_SIZE %d\n", component_size_attr);


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Rank\n");

        attr_id = H5Aopen_name(dset_id, "Rank");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, &field_rank_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        if (io_log) fprintf(log_fptr, "RANK %d\n", field_rank_attr);


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Dimensions\n");

        attr_id = H5Aopen_name(dset_id, "Dimensions");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        attr_count = field_rank_attr;

        attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
          if (io_log) fprintf(log_fptr, "Attr_dsp_id: %d\n", attr_dsp_id);
          assert( attr_dsp_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, field_dims_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(attr_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose: %d\n", h5_status);
          assert( h5_status != h5_error );


        for (idim = 0; idim < field_rank_attr; idim++)
        {
          if (io_log) fprintf(log_fptr, "DIMS %d:  %d\n", idim, (int) field_dims_attr[idim]);
        }

        xfer_size = 1;

        for ( idim = 0; idim < field_rank_attr; idim++ )
        {
          xfer_size = xfer_size * field_dims_attr[idim];
        }

        if (io_log) fprintf(log_fptr, "  Grid Elements %d\n", (int) xfer_size);

        Slab_Rank = field_rank_attr+1;

        Slab_Dims[0] = component_rank_attr;

        for ( idim = 1; idim < Slab_Rank; idim++ )
        {
          Slab_Dims[idim] = field_dims_attr[idim-1];
        }

        if (io_log) fprintf(log_fptr, "  Extended Rank %d\n", (int) Slab_Rank);
        for ( idim = 0; idim < Slab_Rank; idim++ )
        {
          if (io_log) fprintf(log_fptr, "    %d:  %d\n", idim, (int) Slab_Dims[idim]);
        }

        slab_offset[0] = Part;
        slab_stride[0] = 1;
        slab_count[0] = 1;
        slab_block[0] = 1;

        xfer_size = 1;

        for ( idim=1; idim < Slab_Rank; idim++)
        {
          slab_offset[idim] = 0;
          slab_stride[idim] = 1;
          slab_count[idim] = field_dims_attr[idim-1];
          slab_block[idim] = 1;
          xfer_size = xfer_size * slab_count[idim];
        }

        for ( idim=0; idim < Slab_Rank; idim++)
        {
          if (io_log) fprintf(log_fptr, "Slab [%d] Offset %d, Count %d, Stride %d, Block %d\n",
                  idim, (int) slab_offset[idim], (int) slab_count[idim],
                       (int) slab_stride[idim], (int) slab_block[idim]);
        }
        if (io_log) fprintf(log_fptr, "Xfer_size %d\n", (int) xfer_size);
        if (io_log) fprintf(log_fptr, "Comp_size %d\n", (int) component_size_attr);

// Data in memory is considered 1D, stride 1, with zero offset

        mem_stride = 1;
        mem_count = xfer_size;
        mem_offset = 0;
        mem_block = 1;

// 1D memory model

        mem_dsp_id = H5Screate_simple(1, &xfer_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %d\n", mem_dsp_id);
          assert( mem_dsp_id != h5_error );

        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, &mem_block);
          if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %d\n", h5_status);
          assert( h5_status != h5_error );

        file_dsp_id = H5Screate_simple(Slab_Rank, Slab_Dims, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %d\n", file_dsp_id);
          assert( file_dsp_id != h5_error );

        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
          if (io_log) fprintf(log_fptr, "H5Sselect file slab: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) TempField);
          if (io_log) fprintf(log_fptr, "H5Dread: %d\n", (int) h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(mem_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(file_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Dclose(dset_id);
          if (io_log) fprintf(log_fptr, "H5Dclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Fclose(file_id);
          if (io_log) fprintf(log_fptr, "H5Fclose: %d\n", h5_status);
          assert( h5_status != h5_error );

	  for (n = 0, i = 0; i < TotalParticleCount; i++)
          {
            XMask = TRUE;
            MaskAddr = i;
            WordAddr = MaskAddr/BitsPerInt;
            BitAddr  = MaskAddr%BitsPerInt;
            TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);

            if ( TestBit == 0 )
              XMask = FALSE;

	    if (XMask)
	      ParticleVelocity[dim][n++] = TempField[i];
          }

        }  /* End of loop over particle velocity dimensions */
      
      /* Read particle mass. */

      if (CosmologySimulationParticleMassName != NULL) {

        Part = 0;

        if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", CosmologySimulationParticleMassName);

        file_id = H5Fopen(CosmologySimulationParticleMassName, H5F_ACC_RDONLY, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Fopen id: %d\n", file_id);
          assert( file_id != h5_error );

        if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", CosmologySimulationParticleMassName);

        dset_id = H5Dopen(file_id, CosmologySimulationParticleMassName);
          if (io_log) fprintf(log_fptr, "H5Dopen id: %d\n", dset_id);
          assert( dset_id != h5_error );


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Rank\n");

        attr_id = H5Aopen_name(dset_id, "Component_Rank");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, &component_rank_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        if (io_log) fprintf(log_fptr, "COMPONENT_RANK %d\n", component_rank_attr);


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Size\n");

        attr_id = H5Aopen_name(dset_id, "Component_Size");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, &component_size_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        if (io_log) fprintf(log_fptr, "COMPONENT_SIZE %d\n", component_size_attr);


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Rank\n");

        attr_id = H5Aopen_name(dset_id, "Rank");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, &field_rank_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        if (io_log) fprintf(log_fptr, "RANK %d\n", field_rank_attr);


        if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Dimensions\n");

        attr_id = H5Aopen_name(dset_id, "Dimensions");
          if (io_log) fprintf(log_fptr, "H5Aopen_name id: %d\n", attr_id);
          assert( attr_id != h5_error );

        attr_count = field_rank_attr;

        attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
          if (io_log) fprintf(log_fptr, "Attr_dsp_id: %d\n", attr_dsp_id);
          assert( attr_dsp_id != h5_error );

        h5_status = H5Aread(attr_id, HDF5_I4, field_dims_attr);
          if (io_log) fprintf(log_fptr, "H5Aread: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        for (idim = 0; idim < field_rank_attr; idim++)
        {
          if (io_log) fprintf(log_fptr, "DIMS %d:  %d\n", idim, (int) field_dims_attr[idim]);
        }

        xfer_size = 1;

        for ( idim = 0; idim < field_rank_attr; idim++ )
        {
          xfer_size = xfer_size * field_dims_attr[idim];
        }

        if (io_log) fprintf(log_fptr, "  Grid Elements %d\n", (int) xfer_size);

        Slab_Rank = field_rank_attr+1;

        Slab_Dims[0] = component_rank_attr;

        for ( idim = 1; idim < Slab_Rank; idim++ )
        {
          Slab_Dims[idim] = field_dims_attr[idim-1];
        }

        if (io_log) fprintf(log_fptr, "  Extended Rank %d\n", (int) Slab_Rank);
        for ( idim = 0; idim < Slab_Rank; idim++ )
        {
          if (io_log) fprintf(log_fptr, "    %d:  %d\n", idim, (int) Slab_Dims[idim]);
        }

        slab_offset[0] = Part;
        slab_stride[0] = 1;
        slab_count[0] = 1;
        slab_block[0] = 1;

        xfer_size = 1;

        for ( idim=1; idim < Slab_Rank; idim++)
        {
          slab_offset[idim] = 0;
          slab_stride[idim] = 1;
          slab_count[idim] = field_dims_attr[idim-1];
          slab_block[idim] = 1;
          xfer_size = xfer_size * slab_count[idim];
        }

        for ( idim=0; idim < Slab_Rank; idim++)
        {
          if (io_log) fprintf(log_fptr, "Slab [%d] Offset %d, Count %d, Stride %d, Block %d\n",
                  idim, (int) slab_offset[idim], (int) slab_count[idim],
                       (int) slab_stride[idim], (int) slab_block[idim]);
        }
        if (io_log) fprintf(log_fptr, "Xfer_size %d\n", (int) xfer_size);
        if (io_log) fprintf(log_fptr, "Comp_size %d\n", (int) component_size_attr);

// Data in memory is considered 1D, stride 1, with zero offset

        mem_stride = 1;
        mem_count = xfer_size;
        mem_offset = 0;
        mem_block = 1;

// 1D memory model

        mem_dsp_id = H5Screate_simple(1, &xfer_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %d\n", mem_dsp_id);
          assert( mem_dsp_id != h5_error );

        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, &mem_block);
          if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %d\n", h5_status);
          assert( h5_status != h5_error );

        file_dsp_id = H5Screate_simple(Slab_Rank, Slab_Dims, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %d\n", file_dsp_id);
          assert( file_dsp_id != h5_error );

        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
          if (io_log) fprintf(log_fptr, "H5Sselect file slab: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) TempField);
          if (io_log) fprintf(log_fptr, "H5Dread: %d\n", (int) h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(mem_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(file_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Dclose(dset_id);
          if (io_log) fprintf(log_fptr, "H5Dclose: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Fclose(file_id);
          if (io_log) fprintf(log_fptr, "H5Fclose: %d\n", h5_status);
          assert( h5_status != h5_error );

	for (n = 0, i = 0; i < TotalParticleCount; i++)
        {
          XMask = TRUE;
          MaskAddr = i;
          WordAddr = MaskAddr/BitsPerInt;
          BitAddr  = MaskAddr%BitsPerInt;
          TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);

          if ( TestBit == 0 )
            XMask = FALSE;

	  if (XMask)
	    ParticleMass[n++] = TempField[i];
        }

      }  /* End of loop over particle mass dimension */

/*
      for (idim = 0; idim < GridRank; idim++)
      {
        pcol(ParticlePosition[idim], NumberOfParticles, 10, log_fptr);
      }
      for (idim = 0; idim < GridRank; idim++)
      {
        fcol(ParticleVelocity[idim], NumberOfParticles, 10, log_fptr);
      }
      fcol(ParticleMass, NumberOfParticles, 10, log_fptr);
*/

      /* delete mask and temp field. */

      if (io_log) fprintf(log_fptr, "De-Allocate Mask and TempField\n");

//      delete [] Mask;
      delete [] TempField;
      delete [] BitMask;
      delete [] BitArray;


    } else {


      /* For regular IO (non-parallel) */

      if (io_log) fprintf(log_fptr, "CSIG Serial I/O\n");

      /* Set Number Of particles. */

      NumberOfParticles = TempIntArray[0];

      End[0] = NumberOfParticles - 1;
      Dim[0] = NumberOfParticles;

      /* Allocate space for the particles. */

      this->AllocateNewParticles(NumberOfParticles);
  
      /* Read particle positions. */

      for (dim = 0; dim < GridRank; dim++) {
	if (READFILE(CosmologySimulationParticlePositionName, 1, Dim,
		     Start, End, Zero, NULL, &tempbuffer, dim, 3) == FAIL) {
	  fprintf(stderr, "Error reading particle position %d.\n", dim);
	  return FAIL;
	}
	for (i = Start[0]; i <= End[0]; i++)
	  ParticlePosition[dim][i] = FLOAT(tempbuffer[i]);

//      pcol32(tempbuffer, NumberOfParticles, 10, log_fptr);

        if (io_log) fprintf(log_fptr, "De-Allocate [] tempbuffer for dim = %d\n", dim);

	delete [] tempbuffer;
      }

      /* Read particle velocities. */

      if (CosmologySimulationParticleVelocityName != NULL)
	for (dim = 0; dim < GridRank; dim++) {
	  if (READFILE(CosmologySimulationParticleVelocityName, 1, Dim,
	     Start, End, Zero, ParticleVelocity[dim], &tempbuffer, dim, 3) == FAIL) {
	    fprintf(stderr, "Error reading particle velocity %d.\n", dim);
	    return FAIL;
	  }
//        fcol(ParticleVelocity[dim], NumberOfParticles, 10, log_fptr);
	}
  
      /* Read particle mass. */

      if (CosmologySimulationParticleMassName != NULL) {
	if (READFILE(CosmologySimulationParticleMassName, 1, Dim, Start, End,
		     Zero, ParticleMass, &tempbuffer, 0, 1) == FAIL) {
	  fprintf(stderr, "Error reading particle masses.\n");
	  return FAIL;
	}
//      fcol(ParticleMass, NumberOfParticles, 10, log_fptr);
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

  if (io_log) fclose(log_fptr);
 
  return SUCCESS;
}
