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
/  GRID CLASS (INTERPOLATE FIELD VALUES FROM PARENT GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/    This function interpolates boundary values from the parent grid
/    (specified in the argument) to the current grid.  The interpolation
/    used should be monotonic and conservative (and preferably third-
/    order accurate).  The values are also interpolated linearly in time,
/    using the OldBaryonField and BaryonField of the parent.  Make sure
/    that both of these fields have intact boundary values.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
/* function prototypes */

int FindField(int f, int farray[], int n);
extern "C" void FORTRAN_NAME(interpolate)
                             (int *rank, float *pfield, int pdim[],
			      int pis[], int pie[], int r[], 
			      float *field, int dim[], int is[], float *work,
			      interpolation_type *imethod, int *posflag);

extern "C" void FORTRAN_NAME(mhd_interpolate)
                               (float *parentx, float *parenty, float *parentz,
				int parentdim[], int refine[],
				float * childx, float * childy, float *childz,
				int childdim[], int childstart[], int refinedim[],
				float * otherbx, float * otherby, float * otherbz,
				int otherdim[], int otherstart[],
				float* DyBx, float* DzBx, float* DyzBx,
				float* DxBy, float* DzBy, float* DxzBy,
				float* DxBz, float* DyBz, float* DxyBz,
				int* DBxFlag,int* DByFlag,int* DBzFlag,
				FLOAT *cellsizex, FLOAT *cellsizey, FLOAT *cellsizez,
                                MHD_interpolation_type *method, int * step, int * counter);



extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest, 
				     int *sdim1, int *sdim2, int *sdim3, 
				     int *ddim1, int *ddim2, int *ddim3,
				     int *sstart1, int *sstart2, int *sstart3, 
				     int *dstart1, int *dstart2, int *dstart3);
extern "C" void FORTRAN_NAME(mult3d)(float *source, float *dest, 
				     int *sdim1, int *sdim2, int *sdim3, 
				     int *ddim1, int *ddim2, int *ddim3,
				     int *sstart1, int *sstart2, int *sstart3, 
				     int *dstart1, int *dstart2, int *dstart3);
extern "C" void FORTRAN_NAME(div3d)(float *source, float *dest, 
				    int *sdim1, int *sdim2, int *sdim3, 
				    int *ddim1, int *ddim2, int *ddim3,
				    int *sstart1, int *sstart2, int *sstart3, 
				    int *dstart1, int *dstart2, int *dstart3,
                                    int *rstart1, int *rstart2, int *rstart3,
                                    int *rend1, int *rend2, int *rend3);
extern "C" void FORTRAN_NAME(combine3d)(
               float *source1, float *weight1, float *source2, float *weight2,
	       float *dest, int *sdim1, int *sdim2, int *sdim3, 
	       int *ddim1, int *ddim2, int *ddim3,
	       int *sstart1, int *sstart2, int *sstart3, 
	       int *dstart1, int *dstart2, int *dstart3,
	       int *ivel_flag, int *irefine);

/* Interpolate Boundary From Parent function */

int grid::InterpolateFieldValues(grid *ParentGrid, LevelHierarchyEntry * OldFineLevel)
{

  /* set grid time to the parent grid time */

  Time = ParentGrid->Time;

  /* Return if this doesn't involve us. */

  if (ProcessorNumber != MyProcessorNumber &&
      ParentGrid->ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int StartIndex[MAX_DIMENSION], EndIndex[MAX_DIMENSION];
  int ParentStartIndex[MAX_DIMENSION], ParentTempEndIndex[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION], Offset[MAX_DIMENSION];
  int ZeroVector[MAX_DIMENSION], ParentTempStartIndex[MAX_DIMENSION];
  int ParentTempDim[MAX_DIMENSION], TempDim[MAX_DIMENSION];
  int ParentDim[MAX_DIMENSION];
  int ParentTempSize, WorkSize, TempSize, GridSize, One = 1, Zero = 0;
  int dim, field;
  float *TemporaryField, *TemporaryDensityField, *Work, 
        *ParentTemp[MAX_NUMBER_OF_BARYON_FIELDS], *FieldPointer;
  
  if (NumberOfBaryonFields > 0) {

    /* Compute refinement factors and set zero. */

    ParentGrid->ComputeRefinementFactors(this, Refinement);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ZeroVector[dim] = 0;

    /* Find density field */

    int densfield;
    if ((densfield=FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "No density field!\n");
      return FAIL;
    }

    /* Set up array of flags if we are using SecondOrderB interpolation
       method.  These flags indicate if the quantity must always by > 0. 
       For Zeus interpolation, they are overloaded as cell vs. face flags. */

    int SecondOrderBFlag[MAX_NUMBER_OF_BARYON_FIELDS];
    if (InterpolationMethod == SecondOrderB)
      for (field = 0; field < NumberOfBaryonFields; field++) {
	if (FieldType[field] == TotalEnergy || FieldType[field] == Pressure ||
	    FieldType[field] == InternalEnergy)
	  SecondOrderBFlag[field] = 2;   // enforce monotonicity
	else
	  SecondOrderBFlag[field] = 2;   // enforce only positivity
	if (FieldType[field] >= Velocity1 && FieldType[field] <= Velocity3)
	  SecondOrderBFlag[field] = 2;   //  no positivity for velocity
      }
    if (HydroMethod == Zeus_Hydro)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (FieldType[field] >= Velocity1 && FieldType[field] <= Velocity3)
	  SecondOrderBFlag[field] = FieldType[field] - Velocity1 + 1;
	else
	  SecondOrderBFlag[field] = 0;

    /* Compute the start and end indicies (in this grid units) of this grid 
       within it's parent */
    /* StartIndex = cells from left edge of parent active region
	                      to left edge of grid total region
		  + boundary cells of parent grid (current grid units). */
    /*   EndIndex = cells from left  edge of parent active region
	                      to right edge of grid total region
		  + boundary cells of parent grid (current grid units). */
    /* Ee adjust StartIndex and EndIndex to start and end at a parental
       cell edge - this means interpolating a larger region than is
       necessary. */

    ParentTempSize = TempSize = WorkSize = GridSize = 1;
    for (dim = 0; dim < GridRank; dim++) {

      StartIndex[dim] = nint((CellLeftEdge[dim][0] - 
			      ParentGrid->CellLeftEdge[dim][0])/
			     CellWidth[dim][0]);
      EndIndex[dim] = nint((CellLeftEdge[dim][GridDimension[dim]-1] +
			    CellWidth[dim][GridDimension[dim]-1]
			    - ParentGrid->CellLeftEdge[dim][0])/
			   CellWidth[dim][0])
		      - 1;

      /* Record the offset between the real StartIndex and the one adjusted
	 to conform with the parent. */

      Offset[dim] = StartIndex[dim] -
	            int(StartIndex[dim]/Refinement[dim])*Refinement[dim];

      /* Adjust Start/end index to conform with parental cell edges. */

      StartIndex[dim] = int(StartIndex[dim]/Refinement[dim])*Refinement[dim];
      EndIndex[dim] = (int(EndIndex[dim]/Refinement[dim])+1)*Refinement[dim]-1;

      /* Compute values for the ParentTemp fields. */

      ParentStartIndex[dim]     = StartIndex[dim] / Refinement[dim];
      ParentTempDim[dim]        = (EndIndex[dim] - StartIndex[dim] + 1) / 
	                          Refinement[dim];
      ParentDim[dim]            = ParentGrid->GridDimension[dim];

      /* Set the start and and end indicies to pass to interpolate. */

      ParentTempStartIndex[dim] = Refinement[dim];
      ParentTempEndIndex[dim]   = Refinement[dim]*(ParentTempDim[dim]+1) - 1;

      /* Add an extra cell to each side of the ParentTemp field since the
         interpolator will require it (and error check). */

      ParentStartIndex[dim]     -= 1;
      ParentTempDim[dim]        += 2;
      if (ParentStartIndex[dim] < 0 || 
          ParentStartIndex[dim]+ParentTempDim[dim] > 
          ParentGrid->GridDimension[dim]) {
        fprintf(stderr, "Parent grid not big enough for interpolation.\n");
        fprintf(stderr, " ParentStartIndex[%d] = %d  ParentTempDim = %d\n",
                dim, ParentStartIndex[dim], ParentTempDim[dim]);
        return FAIL;
      }

      /* Compute the dimensions of the current grid temporary field. */

      TempDim[dim]            = EndIndex[dim] - StartIndex[dim] + 1;

      /* Compute size of current grid (in floats) for the temp fields. */

      ParentTempSize *= ParentTempDim[dim];
      TempSize       *= TempDim[dim];
      WorkSize       *= (TempDim[dim]/Refinement[dim] + 1);
      GridSize       *= GridDimension[dim];
    }
    //fprintf(stderr,"0: EndIndex %d %d %d \n", EndIndex[0], EndIndex[1], EndIndex[2]);
    //fprintf(stderr,"0: StartIndex %d %d %d \n", StartIndex[0], StartIndex[1], StartIndex[2]);

    /* Fill out the rest if dim < MAX_DIMENSION. */

    for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
      ParentTempDim[dim]      = 1;
      ParentStartIndex[dim]   = 0;
      ParentTempEndIndex[dim] = 0;
      ParentDim[dim]          = 1;
      TempDim[dim]            = 1;
      Offset[dim]             = 0;
    }
    
    /* Copy data from other processor if needed (modify ParentDim and
       ParentStartIndex to reflect the fact that we are only coping part of
       the grid. */

    if (ProcessorNumber != ParentGrid->ProcessorNumber) {
      ParentGrid->CommunicationSendRegion(ParentGrid, ProcessorNumber, 
			ALL_FIELDS, NEW_ONLY, ParentStartIndex, ParentTempDim);
      for (dim = 0; dim < GridRank; dim++) {
	ParentDim[dim] = ParentTempDim[dim];
	ParentStartIndex[dim] = 0;
      }
    }

    /* Return if this is not our concern. */

    if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;
    
    /* Allocate temporary space. */

    TemporaryField        = new float[TempSize];
    TemporaryDensityField = new float[TempSize];
    Work                  = new float[WorkSize];
    for (field = 0; field < NumberOfBaryonFields; field++)
      ParentTemp[field]     = new float[ParentTempSize];

    /* Copy just the required section from the parent fields to the temp
       space. */

    //    if (HydroMethod != Zeus_Hydro)
      for (field = 0; field < NumberOfBaryonFields; field++)
	FORTRAN_NAME(copy3d)(ParentGrid->BaryonField[field], ParentTemp[field],
			     ParentDim, ParentDim+1, ParentDim+2,
			     ParentTempDim, ParentTempDim+1, ParentTempDim+2,
			     &Zero, &Zero, &Zero,
			     ParentStartIndex, ParentStartIndex+1, 
			     ParentStartIndex+2);
/*
    if (HydroMethod == Zeus_Hydro)
      for (field = 0; field < NumberOfBaryonFields; field++) {
	float FloatOne = 1.0, FloatZero = 0.0;
	int VelocityShiftFlag = 0;
	if (FieldType[field] >= Velocity1 && FieldType[field] <= Velocity3)
	  VelocityShiftFlag = FieldType[field] - Velocity1 + 1;
	FORTRAN_NAME(combine3d)(ParentGrid->BaryonField[field], &FloatOne,
			       ParentGrid->BaryonField[field], &FloatZero,
			       ParentTemp[field],
			       ParentDim, ParentDim+1, ParentDim+2,
			       ParentTempDim, ParentTempDim+1, ParentTempDim+2,
			       &Zero, &Zero, &Zero,
			       ParentStartIndex, ParentStartIndex+1, 
			       ParentStartIndex+2, 
			       &VelocityShiftFlag, Refinement);
      }
*/
    /* Multiply ParentTemp fields by their own density to get conserved
       quantities. */

    if (ConservativeInterpolation)
      for (field = 0; field < NumberOfBaryonFields; field++){
	if ( FieldType[field] == TotalEnergy &&
	     ( HydroMethod == MHD_Harten ||
	       HydroMethod == Athena ||
	       HydroMethod == MHD_None ||
	       HydroMethod == MHD_Li ||
	       HydroMethod == PPM_Local ) )
	  continue;
	
	if (FieldTypeIsDensity(FieldType[field]) == FALSE)
	  FORTRAN_NAME(mult3d)(ParentTemp[densfield], ParentTemp[field],
                               &ParentTempSize, &One, &One, 
			       &ParentTempSize, &One, &One,
                               &Zero, &Zero, &Zero, &Zero, &Zero, &Zero);
      }
    /* Do the interpolation for the density field. */

    if (HydroMethod == Zeus_Hydro)
      InterpolationMethod = (SecondOrderBFlag[densfield] == 0) ? 
	SecondOrderA : SecondOrderC;

    FORTRAN_NAME(interpolate)(&GridRank, 
			      ParentTemp[densfield], ParentTempDim,
			      ParentTempStartIndex, ParentTempEndIndex, 
                                 Refinement,
			      TemporaryDensityField, TempDim, ZeroVector, Work,
			      &InterpolationMethod, 
			      &SecondOrderBFlag[densfield]);

    /* Loop over all the fields. */
    
    for (field = 0; field < NumberOfBaryonFields; field++) {

    /* Interpolating from the ParentTemp field to a Temporary field.  This
       is done for the entire current grid, not just it's boundaries.
       (skip density since we did it already) */

      if (HydroMethod == Zeus_Hydro)
	InterpolationMethod = (SecondOrderBFlag[field] == 0) ? 
	  SecondOrderA : SecondOrderC;

      if (FieldType[field] != Density)
	FORTRAN_NAME(interpolate)(&GridRank, 
				  ParentTemp[field], ParentTempDim,
				  ParentTempStartIndex, ParentTempEndIndex, 
                                     Refinement,
				  TemporaryField, TempDim, ZeroVector, Work,
				  &InterpolationMethod, 
				  &SecondOrderBFlag[field]);

      /* Divide by density field to convert from conserved to physical 
         variables (skipping density). */

      if (ConservativeInterpolation)
	if (FieldTypeIsDensity(FieldType[field]) == FALSE)
	  if ( ! ( FieldType[field] == TotalEnergy &&
		   ( HydroMethod == MHD_Harten ||
		     HydroMethod == Athena ||
		     HydroMethod == MHD_None ||
		     HydroMethod == MHD_Li ||
		     HydroMethod == PPM_Local ) ) )
	    FORTRAN_NAME(div3d)(TemporaryDensityField, TemporaryField,
				&TempSize, &One, &One, 
				&TempSize, &One, &One,
				&Zero, &Zero, &Zero, &Zero, &Zero, &Zero,
				&Zero, &Zero, &Zero, &TempSize, &Zero, &Zero);
      
      /* Set FieldPointer to either the correct field (density or the one we
	 just interpolated to). */

      if (FieldType[field] == Density)
	FieldPointer = TemporaryDensityField;
      else
	FieldPointer = TemporaryField;

      /* Copy needed portion of temp field to current grid. */

      if (BaryonField[field] == NULL)
	BaryonField[field] = new float[GridSize];
      if (BaryonField[field] == NULL) {
	fprintf(stderr, "malloc error (out of memory?)\n");
	return FAIL;
      }
    

      FORTRAN_NAME(copy3d)(FieldPointer, BaryonField[field],
			   TempDim, TempDim+1, TempDim+2,
			   GridDimension, GridDimension+1, GridDimension+2,
			   &Zero, &Zero, &Zero,
			   Offset, Offset+1, Offset+2);

    } // end loop over fields
    
    //      fprintf(stderr,"[[[[[[[[[[[[[[[[ yo slippery ]]]]]]]]]]]]]]]]\n");

    //      fprintf(stderr,"EndIndex %d %d %d \n", EndIndex[0], EndIndex[1], EndIndex[2]);
    //fprintf(stderr,"StartIndex %d %d %d \n", StartIndex[0], StartIndex[1], StartIndex[2]);

    if( MHD_Used == 1 ){
      //This is used for the averaging of the Electric field.

      //It's probably not necessary to make the TempSize an array, but fuckit.
      //float * MHDParentTemp[3], *MHDChildTemp[3];
      int MHDParentTempDims[3][3], MHDChildTempDims[3][3];

      float *MHDChildTemp[3], *dummy = new float;

      
      int MHDParentTempSize[3]={1,1,1}, MHDChildTempSize[3]={1,1,1}, One[3] = {1,1,1};
      int i;
      int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;      
      
      *dummy = 1;
      for( field=0;field<3;field++ ){

	//save these values for use in the prolongation,
	//called from Grid_CopyZonesFromGrid
	MHDParentTempPermanent[field] = ParentTempDim[field];
	MHDRefinementFactors[field] = Refinement[field];
        ParentDx =ParentGrid->CellWidth[0][0];
	ParentDy =ParentGrid->CellWidth[1][0];
	ParentDz =ParentGrid->CellWidth[2][0];


	for(dim=0;dim<3;dim++){
	  MHDParentTempDims[field][dim] = ParentTempDim[dim]+MHDAdd[field][dim];
	  MHDChildTempDims[field][dim] = TempDim[dim] + MHDAdd[field][dim];
	  MHDParentTempSize[field] *= MHDParentTempDims[field][dim];
	  MHDChildTempSize[field] *= MHDChildTempDims[field][dim];
	  
	}
	//	fprintf(stderr, "Parent Temp Size %d \n ", MHDParentTempSize[field] );
	
	if(MHDParentTemp[field] != NULL ) 
	  delete MHDParentTemp[field];
	
	MHDParentTemp[field] = new float[ MHDParentTempSize[field] ];
	MHDChildTemp[field] = new float[ MHDChildTempSize[field] ];

	if( CenteredB[field] == NULL ) 
	  CenteredB[field] = new float[GridSize];

	if( MagneticField[field] == NULL )
	  MagneticField[field] = new float[ MagneticSize[field] ];
	

	for(i=0;i<GridSize;i++) CenteredB[field][i] = 1.0;
	for(i=0;i<MagneticSize[field]; i++) MagneticField[field][i] = 1.0;
	for(i=0;i<MHDParentTempSize[field];i++) MHDParentTemp[field][i] = 1.0;
	for(i=0;i<MHDChildTempSize[field];i++) MHDChildTemp[field][i] = 1.0;
	/*
	fprintf(stderr,"MHDParentTempDim (field %d) %d %d %d\n"
		, field, MHDParentTempDims[field][0],
		MHDParentTempDims[field][1],MHDParentTempDims[field][2]);
	fprintf(stderr,"MHDChildTempDim (field %d) %d %d %d\n"
		, field, MHDChildTempDims[field][0],
		MHDChildTempDims[field][1],MHDChildTempDims[field][2]);

	fprintf(stderr, "Cell Width %"PSYM" %"PSYM" %"PSYM, ParentGrid->CellWidth[0][0], 
				    ParentGrid->CellWidth[1][0], 
		ParentGrid->CellWidth[2][0]);
	*/
	FORTRAN_NAME(copy3d)(ParentGrid->MagneticField[field], MHDParentTemp[field],
			     ParentGrid->MagneticDims[field], 
			     ParentGrid->MagneticDims[field]+1, 
			     ParentGrid->MagneticDims[field]+2,
			     MHDParentTempDims[field], 
			     MHDParentTempDims[field]+1, 
			     MHDParentTempDims[field]+2,
			     &Zero,&Zero,&Zero,
			     ParentStartIndex, ParentStartIndex+1,
			     ParentStartIndex+2);
      }//field loop
    
      //Dimensions passed in refer to cell centered quantities.  This is done to 
      //minimize variables called. 
      /*
      fprintf(stderr,"ParentDim %d %d %d\n", ParentDim[0],ParentDim[1],ParentDim[2]);
      fprintf(stderr,"ParentStartIndex %d %d %d\n"
	      , ParentStartIndex[0],ParentStartIndex[1],ParentStartIndex[2]);
      fprintf(stderr,"tempdim %d %d %d\n", TempDim[0], TempDim[1], TempDim[2]);

      fprintf(stderr,"[[[[[[[[[[[[[[[[ yo slippery ]]]]]]]]]]]]]]]]\n");
      */

      //Allocates space for DxyBz, etc, derivatives.  The derivatives are currently not
      //used in this stage, but space should be allocated anyway.  They will be used in 
      //the prolongation crap later.

      //fprintf(stderr,"Heycf1p\n");
      MHD_ProlongAllocate(TempDim);
      MHD_DCheck(TempDim, "one");

      //fprintf(stderr,"Heycf1 here\n");


      //mhd_interpolate used to take derivatives and do the rest of the interpolation 
      //at the same time, right here.  I then realized that I needed all derivatives from all
      //old fine grids before the new fine grid is actually filled.  
      //So now the code is called twice, with the flag "step" to split this one routine into two.
      //I will re-arrange things in the future, but for now, watch that flag.
      //Please forgive the incomprehensibility of this code.  All should be fixed shortly,
      //but if you're reading this comment, obviously it hasn't.

      //Allocate derivatives on the parent grid.
      int Step = 1;
      //fprintf(stderr,"IVF\n");
      FORTRAN_NAME(mhd_interpolate)(MHDParentTemp[0], MHDParentTemp[1], MHDParentTemp[2],
				    ParentTempDim,    Refinement,
				    MHDChildTemp[0], MHDChildTemp[1], MHDChildTemp[2],
				    TempDim, ZeroVector, TempDim, 
				    dummy, dummy, dummy, One, One,
				      DyBx, DzBx, DyzBx,
				      DxBy, DzBy, DxzBy,
				      DxBz, DyBz, DxyBz,
				    DBxFlag,DByFlag,DBzFlag,
				    ParentGrid->CellWidth[0], 
				    ParentGrid->CellWidth[1], 
				    ParentGrid->CellWidth[2],
				    &MHD_InterpolationMethod, &Step, &Step);

      
      if( MHD_DCheck(TempDim, "1") == FAIL ){
	//fprintf(stderr, "dcheck1\n");
	return FAIL;
      }



      //calculate derivatives for the other, old fine grids.  
      //This is equivalent to the Prolongation
      //described in Balsara, but here we simply calculate all 
      //fine grid derivatives and do a single
      //interpolation.  This avoids any by-hand mucking about with grid 
      //ownership of face centered
      //variables, one of the primary concerns of the MHD 
      //implimentor. (he now realizes as he's finishing.)
      //int grid::MHD_CID(LevelHierarchyEntry * OldFineLevel, int Offset[], int TempDim[], 
      //int Refinement[], float * MHDParentTemp[], int ParentTempDim[]){

      //MHDparentTemp and ParentTempDim aren't actually needed by this routine, 
      //but they're needed
      //for the call to mhd_interpolate.  They aren't actually needed byt that routine, either,
      //but I'm making sure all this shit works before going about the process of seperating the
      //mhd_interpolate into two routines.
      CommunicationDirection = COMMUNICATION_RECEIVE;

      if( MHD_CID(OldFineLevel, Offset, TempDim, Refinement) == FAIL ){
	fprintf(stderr, "Shit!  MHD_CID failed.\n");
	return FAIL;
      }

      if( MHD_DCheck(TempDim, "2") == FAIL ){
	//fprintf(stderr, "dcheck2\n");
	return FAIL;
      }


      CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
      //Now do the actual interpolation with all derivatives found from all old fine grids and 
	//the parent grid.

      Step = 2;
      for(field=0;field<3;field++)
	for(i=0;i<MHDParentTempSize[field];i++){
	  if( MHDParentTemp[field][i] != MHDParentTemp[field][i]  ){
	    fprintf(stderr,"MHDParentTemp nan \n");
	  }
	}


      for(field=0;field<3;field++)
	for(i=0;i<MHDChildTempSize[field];i++)
	  if( MHDChildTemp[field][i] != MHDChildTemp[field][i] )
	    fprintf(stderr,"shit.  nan.\n");

      FORTRAN_NAME(mhd_interpolate)(MHDParentTemp[0], MHDParentTemp[1], MHDParentTemp[2],
				    ParentTempDim,    Refinement,
				    MHDChildTemp[0], MHDChildTemp[1], MHDChildTemp[2],
				    TempDim, ZeroVector, TempDim, 
				    dummy, dummy, dummy, One, One,
				      DyBx, DzBx, DyzBx,
				      DxBy, DzBy, DxzBy,
				      DxBz, DyBz, DxyBz,
				    DBxFlag,DByFlag,DBzFlag,
				    ParentGrid->CellWidth[0], 
				    ParentGrid->CellWidth[1], 
				    ParentGrid->CellWidth[2],
				    &MHD_InterpolationMethod, &Step, &Step);

      
      for(field=0;field<3;field++){
	
	for(i=0;i<MHDChildTempSize[field];i++)
	  if( MHDChildTemp[field][i] != MHDChildTemp[field][i] )
	    fprintf(stderr,"shit.  nan after.\n");
      }//field
      if( MHD_DCheck(TempDim, "I can't believe this thing doesnt work") == FAIL ){
	//fprintf(stderr, "dcheck3\n");
	return FAIL;
      }
      //free memory allocated for Interpolation Derivatives.
      MHD_ProlongFree();

      for(field=0;field<3;field++){
	
	FORTRAN_NAME(copy3d)(MHDChildTemp[field], this->MagneticField[field],
		  MHDChildTempDims[field],MHDChildTempDims[field]+1,MHDChildTempDims[field]+2,
		  MagneticDims[field],MagneticDims[field]+1,MagneticDims[field]+2,
		  &Zero, &Zero, &Zero, 
		  Offset, Offset+1, Offset+2);

	delete [] MHDParentTemp[field];
	MHDParentTemp[field] = NULL;
	delete [] MHDChildTemp[field];



      }//field
      if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					   Vel3Num, TENum) == FAIL) {
	fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
	return FAIL;
      }


      if( this->CenterMagneticField() == FAIL ) {
	fprintf(stderr," error with CenterMagneticField\n");
	return FAIL;
      }
      
    }// MHD_Used

    
    
    
    delete Work;
    delete TemporaryField;
    delete TemporaryDensityField;
    for (field = 0; field < NumberOfBaryonFields; field++)
      delete ParentTemp[field];
    
    /* If using the dual energy formalism, then modify the total energy field
       to maintain consistency between the total and internal energy fields.
       This is necessary because the interpolation introduces small
       descrepancies between the two fields which are normally kept in sync. */
    
    if (DualEnergyFormalism)
      if (this->RestoreEnergyConsistency(ENTIRE_REGION) == FAIL) {
	fprintf(stderr, "Error in grid->RestoreEnergyConsisitency.\n");
	return FAIL;
      }
    //      if (this->RestoreEnergyConsistency(ONLY_BOUNDARY) == FAIL) {
    
  } // end: if (NumberOfBaryonFields > 0)
  
  this->DebugCheck("InterpolateFieldValues (after)");

  return SUCCESS;
  
}
