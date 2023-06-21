/***********************************************************************
/
/  GRID CLASS (DEPOSIT PARTICLE POSITIONS ONTO THE SPECIFIED FIELD)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/     This routine deposits the particle living in this grid into either
/       the GravitatingMassField or the GravitatingMassFieldParticles of
/       the TargetGrid, depending on the value of DepositField.
/     It also moves the particles in this grid so they are at the correct
/       time and adjusts their 'mass' to be appropriate for TargetGrid
/       (since particle 'mass' changed from level to level).
/
/  NOTE: 
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "communication.h"

/* This controls the maximum particle mass which will be deposited in
   the MASS_FLAGGING_FIELD.  Only set in Grid_SetFlaggingField. */

float DepositParticleMaximumParticleMass = 0;

/* This variable is only set in ReadParameterFile. */

float DepositPositionsParticleSmoothRadius = 0;

/* function prototypes */

extern "C" void FORTRAN_NAME(cic_deposit)(FLOAT *posx, FLOAT *posy, 
			FLOAT *posz, int *ndim, int *npositions, 
                        float *densfield, float *field, FLOAT *leftedge, 
                        int *dim1, int *dim2, int *dim3, float *cellsize);

extern "C" void FORTRAN_NAME(smooth_deposit)(FLOAT *posx, FLOAT *posy, 
			FLOAT *posz, int *ndim, int *npositions, 
                        float *densfield, float *field, FLOAT *leftedge, 
                        int *dim1, int *dim2, int *dim3, float *cellsize,
			       float *rsmooth);
float ReturnCPUTime();
double ReturnWallTime();
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
                              int Target, int Tag, MPI_Comm CommWorld, 
			      int BufferSize);
#endif /* USE_MPI */

int grid::DepositParticlePositions(grid *TargetGrid, FLOAT DepositTime,
				   int DepositField)
{

  /* Return if this doesn't concern us or if no particles. */

  if (TargetGrid->CommunicationMethodShouldExit(this) ||
      NumberOfParticles == 0)
    return SUCCESS;

  /* Declarations. */

  int dim, i, j, k, size, index1, index2;
  int Dimension[] = {1,1,1}, OriginalDimension[] = {1,1,1},
    Offset[] = {0,0,0};
  float MassFactor = 1.0, *ParticleMassTemp, *ParticleMassPointer,
    TimeDifference = 0;
  float *DepositFieldPointer, *OriginalDepositFieldPointer, CellSize;
  FLOAT LeftEdge[MAX_DIMENSION], OriginalLeftEdge[MAX_DIMENSION];

  /* DepositField specifies where the particles should go.  Set LeftEdge, 
     Dimension, CellSize, DepositFieldPointer according to it. */

  /* 1) GravitatingMassField. */

  if (DepositField == GRAVITATING_MASS_FIELD) {
    if (TargetGrid->GravitatingMassFieldCellSize <= 0)
      TargetGrid->InitializeGravitatingMassField(RefineBy);
    DepositFieldPointer = TargetGrid->GravitatingMassField;
    CellSize            = float(TargetGrid->GravitatingMassFieldCellSize);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = TargetGrid->GravitatingMassFieldLeftEdge[dim];
      Dimension[dim] = TargetGrid->GravitatingMassFieldDimension[dim];
    }
  }

  /* 2) GravitatingMassFieldParticles. */

  else if (DepositField == GRAVITATING_MASS_FIELD_PARTICLES) {
    if (TargetGrid->GravitatingMassFieldParticlesCellSize <= 0)
      TargetGrid->InitializeGravitatingMassFieldParticles(RefineBy);
    DepositFieldPointer = TargetGrid->GravitatingMassFieldParticles;
    CellSize       = float(TargetGrid->GravitatingMassFieldParticlesCellSize);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = TargetGrid->GravitatingMassFieldParticlesLeftEdge[dim];
      Dimension[dim] = TargetGrid->GravitatingMassFieldParticlesDimension[dim];
    }
  }

  /* 3) MassFlaggingField */

  else if (DepositField == MASS_FLAGGING_FIELD) {
    DepositFieldPointer = TargetGrid->MassFlaggingField;
    CellSize            = float(TargetGrid->CellWidth[0][0]);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = TargetGrid->CellLeftEdge[dim][0];
      Dimension[dim] = TargetGrid->GridDimension[dim];
    }
  }

  /* 4) error */

  else {
    fprintf(stderr, "DepositField = %d not recognized.\n", DepositField);
    return FAIL;
  }

  /* If on different processors, generate a temporary field to hold the
     density. */

  if (ProcessorNumber != TargetGrid->ProcessorNumber) {

    /* If this is the target grid processor, then record the orginal
       field characteristics so we can add it in when the data arrives. */

    for (dim = 0; dim < GridRank; dim++) {
      OriginalLeftEdge[dim] = LeftEdge[dim];
      OriginalDimension[dim] = Dimension[dim];
    }
    OriginalDepositFieldPointer = DepositFieldPointer;

    /* Resize the deposit region so it is just big enough to contain the
       grid where the particles reside. */

    size = 1;
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim] = (int(GridLeftEdge[dim]/CellSize)-2)*CellSize;
      Offset[dim] = nint((LeftEdge[dim] - OriginalLeftEdge[dim])/CellSize);
      if (Offset[dim] < 0) {
    fprintf(stderr, "P(%d): dx=%g/%g  %g %g %g   %g %g %g - %g %g %g  %g %g %g   %d %d %d  %d %d %d\n", 
	    MyProcessorNumber, CellSize, CellWidth[0][0],
	    OriginalLeftEdge[0], OriginalLeftEdge[1], OriginalLeftEdge[2],
	    GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2],
	    GridRightEdge[0], GridRightEdge[1], GridRightEdge[2],
	    LeftEdge[0], LeftEdge[1], LeftEdge[2],
	    Offset[0], Offset[1], Offset[2],
	    Dimension[0], Dimension[1], Dimension[2]);
	fprintf(stderr, "Offset[%d] = %d < 0\n", dim, Offset[dim]);
	return FAIL;
      }
      Dimension[dim] = int((GridRightEdge[dim] - LeftEdge[dim])/CellSize) + 3;
      size *= Dimension[dim];
    }

    /* Allocate buffer to communicate deposit region (unless in receive-mode
       in which case the buffer was already allocated in post-receive mode). */

#ifdef USE_MPI
    if (CommunicationDirection == COMMUNICATION_RECEIVE)
      DepositFieldPointer = 
	CommunicationReceiveBuffer[CommunicationReceiveIndex];
    else {
      DepositFieldPointer = new float[size];
      if (MyProcessorNumber == ProcessorNumber)
	for (i = 0; i < size; i++)
	  DepositFieldPointer[i] = 0;
    }
#endif /* USE_MPI */

  }

  if (MyProcessorNumber == ProcessorNumber) {

    /* If the Target is this grid and the DepositField is MassFlaggingField,
       then multiply the Particle density by the volume to get the mass. */

    if (this == TargetGrid && DepositField == MASS_FLAGGING_FIELD)
      for (dim = 0; dim < GridRank; dim++)
	MassFactor *= CellWidth[dim][0];

    /* If the DepositGrid and this grid are not the same, we must adjust the
       particle 'mass'. */

    if (this != TargetGrid) {

      /* Find the difference in resolution between this grid and TargetGrid. */

      float RefinementFactors[MAX_DIMENSION];
      this->ComputeRefinementFactorsFloat(TargetGrid, RefinementFactors);

      /* Compute the implied difference in 'mass' between particles in this
	 grid and those in TargetGrid. */

      for (dim = 0; dim < GridRank; dim++)
	MassFactor *= RefinementFactors[dim];

    }

    /* If required, Change the mass of particles in this grid. */

    if (MassFactor != 1.0) {
      ParticleMassTemp = new float[NumberOfParticles];
      for (i = 0; i < NumberOfParticles; i++)
	ParticleMassTemp[i] = ParticleMass[i]*MassFactor;
      ParticleMassPointer = ParticleMassTemp;
    } else
      ParticleMassPointer = ParticleMass;

    /* If target field is MASS_FLAGGING_FIELD, then set masses of particles
       which are too large to zero (to prevent run-away refinement). */

    if (DepositField == MASS_FLAGGING_FIELD && 
	DepositParticleMaximumParticleMass > 0 && MassFactor != 1.0)
      for (i = 0; i < NumberOfParticles; i++)
	ParticleMassPointer[i] = min(DepositParticleMaximumParticleMass,
				   ParticleMassPointer[i]);

    /* Compute difference between current time and DepositTime. */

    TimeDifference = DepositTime - Time;

    /* Move particles to positions at Time + TimeDifference. */

    this->UpdateParticlePosition(TimeDifference);

    /* Deposit particles. */

    if (DepositPositionsParticleSmoothRadius < CellSize)

      /* Deposit to field using CIC. */

      FORTRAN_NAME(cic_deposit)(
           ParticlePosition[0], ParticlePosition[1], ParticlePosition[2], 
	   &GridRank, &NumberOfParticles, ParticleMassPointer, DepositFieldPointer, 
	   LeftEdge, Dimension, Dimension+1, Dimension+2, &CellSize);

    else

      /* Deposit to field using large-spherical CIC, with radius of
	 DepositPositionsPartaicleSmoothRadius . */

      FORTRAN_NAME(smooth_deposit)(
           ParticlePosition[0], ParticlePosition[1], ParticlePosition[2], 
	   &GridRank, &NumberOfParticles, ParticleMassPointer, DepositFieldPointer, 
	   LeftEdge, Dimension, Dimension+1, Dimension+2,
	   &CellSize, &DepositPositionsParticleSmoothRadius);

  } // end: if (MyProcessorNumber == ProcessorNumber)

  /* If on different processors, copy deposited field back to the target
     grid and add to the correct field. */

  if (ProcessorNumber != TargetGrid->ProcessorNumber) {

#ifdef USE_MPI

    /* If posting a receive, then record details of call. */

    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
      CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = TargetGrid;
      CommunicationReceiveCallType[CommunicationReceiveIndex] = 3;
      CommunicationReceiveArgument[0][CommunicationReceiveIndex] = DepositTime;
      CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] =
	                                                         DepositField;
    }

    MPI_Status status;
    float time1 = ReturnCPUTime();

    if (MyProcessorNumber == ProcessorNumber)
      CommunicationBufferedSend(DepositFieldPointer, size, MPI_FLOAT, 
			     TargetGrid->ProcessorNumber, MPI_SENDREGION_TAG, 
			     MPI_COMM_WORLD, BUFFER_IN_PLACE);

    if (MyProcessorNumber == TargetGrid->ProcessorNumber &&
	CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
      MPI_Recv(DepositFieldPointer, size, MPI_FLOAT, ProcessorNumber, 
	       MPI_SENDREGION_TAG, MPI_COMM_WORLD, &status);

    if (MyProcessorNumber == TargetGrid->ProcessorNumber &&
	CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      MPI_Irecv(DepositFieldPointer, size, MPI_FLOAT, ProcessorNumber, 
	        MPI_SENDREGION_TAG, MPI_COMM_WORLD, 
	        CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
      CommunicationReceiveBuffer[CommunicationReceiveIndex] = 
	                                                  DepositFieldPointer;
      CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
	  CommunicationReceiveCurrentDependsOn;
      CommunicationReceiveIndex++;
    }      

    CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

    if (MyProcessorNumber == TargetGrid->ProcessorNumber &&
	CommunicationDirection != COMMUNICATION_POST_RECEIVE) {

      index1 = 0;
      for (k = 0; k < Dimension[2]; k++)
	for (j = 0; j < Dimension[1]; j++) {
	  index2 = ((k+Offset[2])*OriginalDimension[1] + j + Offset[1])*
	           OriginalDimension[0] + Offset[0];
	  for (i = 0; i < Dimension[0]; i++)
	    OriginalDepositFieldPointer[index2++] += 
	                                      DepositFieldPointer[index1++];
	}

      delete [] DepositFieldPointer;

    } // end: if (MyProcessorNumber == TargetGrid->ProcessorNumber)

  } // end: If (ProcessorNumber != TargetGrid->ProcessorNumber)

  if (MyProcessorNumber == ProcessorNumber) {

    /* If necessary, delete the particle mass temporary. */

    if (MassFactor != 1.0)
      delete ParticleMassTemp;

    /* Return particles to positions at Time. */

    this->UpdateParticlePosition(-TimeDifference);

  }

  return SUCCESS;
}
