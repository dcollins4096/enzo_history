/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1: James Bordner 2003-06  Added USE_HDF4 and USE_HDF5
/
/  PURPOSE:
/
************************************************************************/

class ExternalBoundary
{
 private:
  int  BoundaryRank;                      // This is the rank and dimension 
  int  BoundaryDimension[MAX_DIMENSION];  //  of the grid to which the boundary
                                          //  values apply
  int  NumberOfBaryonFields;              // Number of boundary fields
  int  BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];           
                                          // Field types

  boundary_type *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
                                          /* Type of boundary value:
                                              1 - reflecting
					      2 - outflow
					      3 - inflow
					      4 - periodic   */

  boundary_type ParticleBoundaryType;

  float *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];  
					  // boundary values for inflow (3)

  friend class grid;						    

 public:
//
// Constructor (set rank and dims from TopGrid)
//
  ExternalBoundary();
//
// Prepare External Boundary (set dims, etc.) based on TopGrid in argument.
//
  int Prepare(class grid *TopGrid);
//
// Checks if the External Boundary has been prepared
//
  int AmIPrepared() {return (BoundaryRank > 0) ? TRUE : FALSE;};
//
// Set one face of external boundaries to a constant value 
//  (Note: this is not suitable for setting inflow conditions as
//         BoundaryValue is not set).
//  Returns SUCCESS or FAIL.
//
  int InitializeExternalBoundaryFace(int Dimension, 
			      boundary_type LeftBoundaryType,
			      boundary_type RightBoundaryType, 
			      float LeftBoundaryValue[],
			      float RightBoundaryValue[]);
//
//  Initialize particle boundary conditions
//
  void InitializeExternalBoundaryParticles(boundary_type ParticleBoundary)
    {ParticleBoundaryType = ParticleBoundary;};
//
// Read an external boundary
//

  // See HDF_functions.C for the implementation: calls HDF4* or HDF5*

  inline int ReadExternalBoundary (FILE *fptr);  

  // The following are private since they should only be called by 
  // ReadExternalBoundary()

 private:
  int ReadExternalBoundaryHDF4(FILE *fptr);
  int ReadExternalBoundaryHDF5(FILE *fptr);
 public:


//
// Write an external boundary
//

  // See HDF_functions.C for the implementation: calls HDF4* or HDF5*

  inline int WriteExternalBoundary (FILE *fptr, char *hdfname);

  // The following are private since they should only be called by
  // WriteExternalBoundary()

 private:
  int WriteExternalBoundaryHDF4(FILE *fptr, char *hdfname);
  int WriteExternalBoundaryHDF5(FILE *fptr, char *hdfname);
 public:

//
// Given a pointer to a field and its field type, find the equivalent
//   field type in the list of boundary's and apply that boundary value/type.
//   Returns: 0 on failure
//

  int SetExternalBoundary(int FieldRank, int GridDims[], int GridOffset[],
                          int StartIndex[], int EndIndex[],
                          float *Field, int FieldType);
//
// This routine handle the boundary conditions for particles.  The conditions
//   are assumed to be the same as the mass field.
//
  int SetExternalBoundaryParticles(int FieldRank, int NumberOfParticles,
                                   FLOAT *Position[], float *Velocity[]);
//
// Finds and returns the indexes to commonly used physical quantities.
//
  int IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, 
                                 int &Vel2Num, int &Vel3Num, int &TENum);
//
/************************************************************************/
//
// WavePool test problem:
//  This routine sets up the inflow boundary conditions to model an inflowing
//   linear wave (from the left boundary).  See also WavePoolGlobalData.h.
//
  int SetWavePoolBoundary(FLOAT time);
//
// ShockPool test problem:
//  This routine sets up the inflow boundary conditions to model an inflowing
//   shock wave (from the left boundary).  See also ShockPoolGlobalData.h.
//
  int SetShockPoolBoundary(FLOAT time);
//
// DoubleMach problem:
//  This routine sets up the necessary inflow boundary conditions.
//
  int SetDoubleMachBoundary(FLOAT time, FLOAT CellLeftEdge[], 
                            FLOAT CellWidth[]);

};

#include "message.h"

/***********************************************************************/

inline int ExternalBoundary::ReadExternalBoundary (FILE *fptr)
{
#if defined (USE_HDF4)
  return ReadExternalBoundaryHDF4 (fptr);
#elif defined (USE_HDF5)
  return ReadExternalBoundaryHDF5 (fptr);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif
}

/***********************************************************************/

inline int ExternalBoundary::WriteExternalBoundary (FILE *fptr, char *hdfname)
{
#if defined (USE_HDF4)
  return WriteExternalBoundaryHDF4 (fptr, hdfname);
#elif defined (USE_HDF5)
  return WriteExternalBoundaryHDF5 (fptr, hdfname);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif
}
