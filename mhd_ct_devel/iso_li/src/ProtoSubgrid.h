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
/  PROTO SUBGRID CLASS
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/


class ProtoSubgrid
{
 private:

  int GridRank;
  int GridDimension[MAX_DIMENSION];

  FLOAT GridLeftEdge[MAX_DIMENSION];
  FLOAT GridRightEdge[MAX_DIMENSION];

  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];

  int NumberFlagged;

  int *GridFlaggingField;
  int *Signature[MAX_DIMENSION];

 public:

  ProtoSubgrid();
  ~ProtoSubgrid();

  int AcceptableSubgrid();
  int ReturnNthLongestDimension(int n);
  int ComputeSignature(int dim);
  int FindGridsByZeroSignature(int dim, int &NumberOfNewGrids, 
			       int GridEnds[MAX_NUMBER_OF_SUBGRIDS][2]);
  int CopyToNewSubgrid(int dim, int GridStart, int GridEnd, 
		       ProtoSubgrid *NewGrid);
  int ComputeSecondDerivative(int dim, int &ZeroCrossStrength, 
			      int GridEnds[2][2]);
  int CopyFlaggedZonesFromGrid(grid *Grid);
  int ShrinkToMinimumSize();
  int CleanUp();

  int ReturnGridRank() {return GridRank;};
  int *ReturnGridDimension() {return GridDimension;};
  FLOAT *ReturnGridLeftEdge() {return GridLeftEdge;};
  FLOAT *ReturnGridRightEdge() {return GridRightEdge;};

};
