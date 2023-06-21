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
/  GRID CLASS (PROJECT SOLUTION IN CURRENT GRID TO PARENT GRID
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  NOTE: This routine assumes that the parent and current grids have the
/        same baryon fields.
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "pout.h"
/* function prototypes */

int FindField(int f, int farray[], int n);
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


int grid::ProjectSolutionToParentGrid(grid &ParentGrid, int MHD_Extraction)
{
  /* Return if this doesn't involve us. */
  //ParentGrid.MHD_GROWTH("PSTPG start Parent");
  //this->MHD_GROWTH("PSTPG start this");
  //fprintf(stderr, "kludge: no project solution to parent\n");
  //return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber && 
      MyProcessorNumber != ParentGrid.ProcessorNumber)
    return SUCCESS;
 
  // What this actually says: Unwrapping all this Negative Logic
  // If SEND: Only run if a.) the proceesors are different and b.) I have the Subgrid.
  // if Receive: run if the processors are the same OR I have the Parent Grid.

  if (CommunicationDirection == COMMUNICATION_SEND &&
      (MyProcessorNumber == ParentGrid.ProcessorNumber || 
       ProcessorNumber == ParentGrid.ProcessorNumber))
    return SUCCESS;

  if (CommunicationDirection == COMMUNICATION_RECEIVE &&
      MyProcessorNumber != ParentGrid.ProcessorNumber &&
      ProcessorNumber != ParentGrid.ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  char basename[30];

  if( WriteInThisF(40) == TRUE) {
    //                                      cycle          subgrid            grid Incremented in ELRO
    sprintf(basename, "data40%d%d%d.grid",dccCounter8, dccCounter6, dccCounter5);

    FILE *dummy = fopen(basename, "a");
    if( ParentGrid.WriteGrid(dummy, basename, MyProcessorNumber) == FAIL ){
      fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
      return FAIL;
    }
    fclose(dummy);

  }

  

  this->DebugCheck("ProjectSolutionToParentGrid (before)");

  /* declarations */
    
  int i, j, k, dim, field, One = 1, Zero = 0, skipi, skipj, skipk;
  int ParentStartIndex[MAX_DIMENSION], ParentDim[MAX_DIMENSION],
      ParentEndIndex[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION], Dim[MAX_DIMENSION];

  /* compute size of current grid fields */
    
  int ParentSize = 1, Size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    Size *= GridDimension[dim];
    ParentSize *= ParentGrid.GridDimension[dim];
    ParentDim[dim] = ParentGrid.GridDimension[dim];
    Dim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
  }

  /* set values that are needed for triply-nested loops. */

  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    ParentDim[dim]        = 1;
    ParentStartIndex[dim] = 0;
    ParentEndIndex[dim]   = 0;
    Dim[dim]              = 1;
  }
    
  /* compute refinement factor */
   
  ParentGrid.ComputeRefinementFactors(this, Refinement);
  
  /* compute the offset (in parent grid index units) from the edge of the
     parent grid to the beginning of the active region in this grid. */

  for (dim = 0; dim < GridRank; dim++) {
    if (GridLeftEdge[dim] >= ParentGrid.GridRightEdge[dim] ||
	GridRightEdge[dim] <= ParentGrid.GridLeftEdge[dim])
      return SUCCESS;
    ParentStartIndex[dim] = nint((GridLeftEdge[dim] - 
				  ParentGrid.GridLeftEdge[dim])/
				 (ParentGrid.CellWidth[dim][0]))
                          + ParentGrid.GridStartIndex[dim];
    ParentEndIndex[dim] = ParentStartIndex[dim] + Dim[dim]/Refinement[dim] - 1;
  }
    
  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  /* Compute the ratio Volume[ThisGridCell]/Volume[ParentCell]. */

  float RelativeVolume = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    RelativeVolume /= float(Refinement[dim]);

  /* Multiply all fields by the density to get conserved quantities. */

  if (ProcessorNumber == MyProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++){
      if ( FieldType[field] == TotalEnergy &&
	   ( HydroMethod == MHD_Harten ||
	     HydroMethod == Athena ||
	     HydroMethod == MHD_None ||
	     HydroMethod == MHD_Li ||
	     HydroMethod == PPM_Local ) )
	continue;
	   
      if (FieldTypeIsDensity(FieldType[field]) == FALSE && 
	  ((FieldType[field] < Velocity1 || FieldType[field] > Velocity3)
	   || HydroMethod != Zeus_Hydro      )       ){
	FORTRAN_NAME(mult3d)(BaryonField[DensNum], BaryonField[field],
			     &Size, &One, &One, &Size, &One, &One,
			     &Zero, &Zero, &Zero, &Zero, &Zero, &Zero);

      }
    }
  /* Allocate Parent temps if it is in another processor. */


  if (ParentGrid.ProcessorNumber != MyProcessorNumber) {
    int ParentSize = 1;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      ParentStartIndex[dim] = 0;
      ParentDim[dim]        = Dim[dim]/Refinement[dim];
      ParentEndIndex[dim]   = ParentDim[dim] - 1;
      ParentSize *= ParentDim[dim];
    }
    for (field = 0; field < NumberOfBaryonFields; field++) {
      if(ParentGrid.BaryonField[field] != NULL )
	delete ParentGrid.BaryonField[field];
      ParentGrid.BaryonField[field] = new float[ParentSize];
    }
  }

  /* For each field, zero the appropriate parental zones. */

  int i1, j1, k1, pindex, gindex;
  if (ProcessorNumber == MyProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++)
      for (k = ParentStartIndex[2]; k <= ParentEndIndex[2]; k++)
	for (j = ParentStartIndex[1]; j <= ParentEndIndex[1]; j++) {
	  pindex = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, pindex++)
	    ParentGrid.BaryonField[field][pindex] = 0.0;
	}


  /* For each field, accumulate it's conserved quantities in the parent 
     grid. */

  if (ProcessorNumber == MyProcessorNumber){

    for (field = 0; field < NumberOfBaryonFields; field++) {
      skipi = skipj = skipk = 1;
      float weight = RelativeVolume;
      if (HydroMethod == Zeus_Hydro) {
	if (FieldType[field] == Velocity1) skipi = Refinement[0];
	if (FieldType[field] == Velocity2) skipj = Refinement[1];
	if (FieldType[field] == Velocity3) skipk = Refinement[2];
      }
      if (skipi*skipj*skipk != 1)
	weight *= float(skipi*skipj*skipk);
      for (k = 0; k < Dim[2]; k += skipk) {
	k1 = k/Refinement[2];
	for (j = 0; j < Dim[1]; j += skipj) {
	  j1 = j/Refinement[1];
	  
	  pindex = (0  + ParentStartIndex[0])                            + 
	    (j1 + ParentStartIndex[1])*ParentDim[0]               +
	    (k1 + ParentStartIndex[2])*ParentDim[0]*ParentDim[1];
	  
	  gindex = 0+GridStartIndex[0]                                      + 
	    (j+GridStartIndex[1])*GridDimension[0]                    +
	    (k+GridStartIndex[2])*GridDimension[0]*GridDimension[1];
	  
	  for (i = 0; i < Dim[0]; i += skipi) { 
	    i1 = i/Refinement[0];
	    ParentGrid.BaryonField[field][pindex+i1] += 
	      BaryonField[field][gindex+i]*weight;
	  }
	}
      }
    }//field
 
    //if(MHD_Used == TRUE && 0==1){ // && MHD_Extraction == TRUE){ dcc 10/14/05 removed 
    if(MHD_Used == TRUE ){

      if(MHD_ProjectE == TRUE ){
       
	int MHDeDim[3][3], MHDeParentDim[3][3], MHDeParentSize[3]={1,1,1};
	float RefineInv;
	for(field=0;field<3;field++){
	  for(dim=0;dim<3;dim++){
	    MHDeDim[field][dim]=Dim[dim]+((field==dim)?0:1);
	    MHDeParentDim[field][dim]=ParentDim[dim]+((field==dim)?0:1);
	    MHDeParentSize[field]*=MHDeParentDim[field][dim];
	  }
	  
	  if(ParentGrid.ProcessorNumber != MyProcessorNumber ){
	    if(ParentGrid.ElectricField[field] != NULL ){
	      fprintf(stderr,"ProjectSolution: ElectricField not null where it should be.\n");
	      fprintf(stderr,"Find out why, and where.\n");
	    }
	      ParentGrid.ElectricField[field]=new float[ MHDeParentSize[field] ];

	  }
	  
	  for(k=ParentStartIndex[2];k<=ParentEndIndex[2]+((field==2)?0:1); k++)
	    for(j=ParentStartIndex[1];j<=ParentEndIndex[1]+((field==1)?0:1);j++)
	      for(i=ParentStartIndex[0];i<=ParentEndIndex[0]+((field==0)?0:1);i++){
		pindex=i+MHDeParentDim[field][0]*(j+MHDeParentDim[field][1]*k);
		ParentGrid.ElectricField[field][pindex]=0.0;
	      }
	}//field
	
	  //Now do the actual projection
	  //Since the Parent and Subgrid electric fields are co-located along one axis,
	  //we skip the interspacing points when doing the projection.
      
      
	  for(field=0;field<3;field++){
	    RefineInv=1.0/Refinement[field];
	    for(k=0;k<MHDeDim[field][2];k+=((field==2)?1:Refinement[2]) ){
	      k1=k/Refinement[2];
	      for(j=0;j<MHDeDim[field][1];j+=((field==1)?1:Refinement[1])){
		j1=j/Refinement[1];
	
		pindex= 0+ParentStartIndex[0]
		  +(j1+ParentStartIndex[1])*MHDeParentDim[field][0]
		  +(k1+ParentStartIndex[2])*MHDeParentDim[field][1]*MHDeParentDim[field][0];
	
		gindex = 0 + GridStartIndex[0]
		  +(j+GridStartIndex[1])*ElectricDims[field][0]
		  +(k+GridStartIndex[2])*ElectricDims[field][1]*ElectricDims[field][0];
                

		//Note that we use AvgElectricField on the subgrid, but ElectricField on the 
		//Parent.  This is because Parent.ElectricField must reflect the time structure
		//of the subgrid advance.
               
		for(i=0;i<MHDeDim[field][0];i+=((field==0)?1:Refinement[0])){
                  i1=i/Refinement[0];
		  ParentGrid.ElectricField[field][pindex+i1] += 
		    AvgElectricField[field][gindex+i]*RefineInv;
		}//i
              	
	      }//j
	    }//k  
	  }//field
        
      }//MHD_ProjectE

      if(MHD_ProjectB == TRUE){


	  fprintf(stderr, "PROJB my proc %d parent proc %d\n", MyProcessorNumber, ParentGrid.ReturnProcessorNumber());
	/*
	  fprintf(stderr,"warning: projb should not be used, for simulations, only extractions.\n");
	  fprintf(stderr,"         Flux correction not acurate enough.\n");
	  fprintf(stderr,"         Also:  It should not be used in parallel, since the proper communication isn't in\n");
	*/
	int MHDDim[3][3], MHDParentDim[3][3], MHDParentSize[3]={1,1,1};
	
	for(field=0;field<3;field++){
	  for(dim=0;dim<3;dim++){
	    MHDDim[field][dim] = Dim[dim]+MHDAdd[field][dim];
	    MHDParentDim[field][dim] = ParentDim[dim]+MHDAdd[field][dim];
	    MHDParentSize[field] *= MHDParentDim[field][dim];
	  }
	  if( ParentGrid.ProcessorNumber != MyProcessorNumber) {
	    fprintf(stderr,"00000000000000000 allocating magnetic field\n");
	    delete ParentGrid.MagneticField[field];
	    ParentGrid.MagneticField[field] = new float[MHDParentSize[field]];
	    
	    
	  }
	  
	  
	}//field
	
	for (field = 0; field < 3; field++)
	  for (k = ParentStartIndex[2]; k <= ParentEndIndex[2]+MHDAdd[field][2]; k++)
	    for (j = ParentStartIndex[1]; j <= ParentEndIndex[1]+MHDAdd[field][1]; j++) 
	      for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]+MHDAdd[field][0]; i++){
		pindex = i+(k*MHDParentDim[field][1] + j)*MHDParentDim[field][0];
		ParentGrid.MagneticField[field][pindex] = 0.0;
		
	      }
	
	for (field = 0; field < 3; field++) {
	  skipi = skipj = skipk = 1;
	  float weight = RelativeVolume;
	  
	  if (field == 0) skipi = Refinement[0];
	  if (field == 1) skipj = Refinement[1];
	  if (field == 2) skipk = Refinement[2];
	  
	  weight *= float(skipi*skipj*skipk);
	  
	  for (k = 0; k < MHDDim[field][2]; k += skipk) {
	    k1 = k/Refinement[2];
	    for (j = 0; j < MHDDim[field][1]; j += skipj) {
	      j1 = j/Refinement[1];
	      
	      pindex = (0  + ParentStartIndex[0])                            + 
		(j1 + ParentStartIndex[1])*MHDParentDim[field][0]     +
		(k1 + ParentStartIndex[2])*MHDParentDim[field][0]*MHDParentDim[field][1];
	      
	      gindex = 0+GridStartIndex[0]                                      + 
		(j+GridStartIndex[1])*MagneticDims[field][0]              +
		(k+GridStartIndex[2])*MagneticDims[field][0]*MagneticDims[field][1];
	      
	      for (i = 0; i < MHDDim[field][0]; i += skipi) { 
		i1 = i/Refinement[0];
		ParentGrid.MagneticField[field][pindex+i1]
		  +=MagneticField[field][gindex+i]*weight;
		
		
	      }
	    }
	  }
	  
	}//field
      }//Proj B
      
      
  }//mhd_used
   
 
}//If MyProcessor==ThisProcessor

  /* If necessary, copy the projected field from the 'fake' ParentGrid to
     the real one. */

  //(MHD_ProjectB
  int FieldToSend = JUST_BARYONS; // When MHD_ProjectFace was doing all the projection,this was the prefered flag.

  if( MHD_Used == TRUE ){

    if( MHD_ProjectB == TRUE ){
      FieldToSend = BARYONS_MAGNETIC; 
    }
    if( MHD_ProjectE == TRUE ){
      FieldToSend =  BARYONS_ELECTRIC;
    }
    
  }

  //int FieldToSend = ALL_FIELDS;     // With the fast SUBling search, this is the prefered method.
  
  if( this->CheckForNans("PSTP: Before Comm SG") == FAIL )
    return FAIL;
  if( ParentGrid.CheckForNans("PSTP: Before Comm PG") == FAIL )
    return FAIL;

  if (ProcessorNumber != ParentGrid.ProcessorNumber) {
    int ParentRegionDim[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ParentRegionDim[dim] = ParentEndIndex[dim] - ParentStartIndex[dim] + 1;
    ParentGrid.CommunicationReceiveRegion(&ParentGrid, ProcessorNumber,
					  FieldToSend, NEW_ONLY, ParentStartIndex, ParentRegionDim, TRUE);
  }

  if( this->CheckForNans("PSTP: after Comm SG") == FAIL )
    return FAIL;
  if( ParentGrid.CheckForNans("PSTP: after Comm PG") == FAIL )
    return FAIL;

  /* Divide all fields by mass to return to original quantity. */
    
  for (field = 0; field < NumberOfBaryonFields; field++){
    if ( FieldType[field] == TotalEnergy &&
	 ( HydroMethod == MHD_Harten ||
	   HydroMethod == Athena ||
	   HydroMethod == MHD_None ||
	   HydroMethod == MHD_Li ||
	   HydroMethod == PPM_Local ) )
      continue;
    
    if (FieldTypeIsDensity(FieldType[field]) == FALSE && 
	((FieldType[field] < Velocity1 || FieldType[field] > Velocity3)	 
          || HydroMethod != Zeus_Hydro      )       ) {
      if (ProcessorNumber == MyProcessorNumber)
	FORTRAN_NAME(div3d)(BaryonField[DensNum], BaryonField[field],
			    &Size, &One, &One, &Size, &One, &One,
			    &Zero, &Zero, &Zero, &Zero, &Zero, &Zero,
			    &Zero, &Zero, &Zero, &Size, &Zero, &Zero);
      if (ParentGrid.ProcessorNumber == MyProcessorNumber)
	FORTRAN_NAME(div3d)(ParentGrid.BaryonField[DensNum],
			    ParentGrid.BaryonField[field],
			    ParentGrid.GridDimension, 
			      ParentGrid.GridDimension+1,
			      ParentGrid.GridDimension+2,
			    ParentGrid.GridDimension, 
			      ParentGrid.GridDimension+1,
			      ParentGrid.GridDimension+2,
			    &Zero, &Zero, &Zero, &Zero, &Zero, &Zero,
			    ParentStartIndex, ParentStartIndex+1, 
                             ParentStartIndex+2,
			    ParentEndIndex, ParentEndIndex+1, ParentEndIndex+2);
			  
    }
  }
    
  /* If appropriate, restore consistency between total and internal
     energy in projected regions. */

  if( this->CheckForNans("PSTP: after a/rho SG") == FAIL )
    return FAIL;
  if( ParentGrid.CheckForNans("PSTP: after a/rho PG") == FAIL )
    return FAIL;

  if (ParentGrid.ProcessorNumber == MyProcessorNumber)
   if (DualEnergyFormalism)
    for (k = ParentStartIndex[2]; k <= ParentEndIndex[2]; k++)
      for (j = ParentStartIndex[1]; j <= ParentEndIndex[1]; j++) {

	i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];

	for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++)
	  ParentGrid.BaryonField[TENum][i1] = 
	    ParentGrid.BaryonField[GENum][i1] + 0.5*
	    ParentGrid.BaryonField[Vel1Num][i1] * 
	    ParentGrid.BaryonField[Vel1Num][i1];

	i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];

	if (GridRank > 1)
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++)
	    ParentGrid.BaryonField[TENum][i1] += 0.5*
	      ParentGrid.BaryonField[Vel2Num][i1] *
	      ParentGrid.BaryonField[Vel2Num][i1];

	i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];

	if (GridRank > 2)
	  for (i = ParentStartIndex[0]; i <= ParentEndIndex[0]; i++, i1++)
	    ParentGrid.BaryonField[TENum][i1] += 0.5*
	      ParentGrid.BaryonField[Vel3Num][i1] *
	      ParentGrid.BaryonField[Vel3Num][i1];
        
#ifdef HAOXU
        if(MHD_Used==TRUE){
           i1 = (k*ParentDim[1] + j)*ParentDim[0] + ParentStartIndex[0];
          for(i = ParentStartIndex[0];i <= ParentEndIndex[0];i++,i1++){
            ParentGrid.BaryonField[TENum][i1]*=ParentGrid.BaryonField[DensNum][i1];
            ParentGrid.BaryonField[TENum][i1]+=0.5*(pow(ParentGrid.CenteredB[0][i1],2)+pow(ParentGrid.CenteredB[1][i1],2)
               +pow(ParentGrid.CenteredB[2][i1],2));
          }
         }
#endif
		  
      } // end: loop over faces

  /* Clean up the fake ParentGrid. */

  if (ParentGrid.ProcessorNumber != MyProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++) {
      delete ParentGrid.BaryonField[field];
      ParentGrid.BaryonField[field] = NULL;
    }


  if( WriteInThisF(41) == TRUE) {
    fprintf(stderr,"moo data41 %d %d %d .grid",dccCounter8, dccCounter6, dccCounter5);
    //                                      cycle          subgrid            grid Incremented in ELRO
    sprintf(basename, "data41%d%d%d.grid",dccCounter8, dccCounter6, dccCounter5);

    FILE *dummy = fopen(basename, "a");
    if( ParentGrid.WriteGrid(dummy, basename, MyProcessorNumber) == FAIL ){
      fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
      return FAIL;
    }
    fclose(dummy);

  }

  if( this->CheckForNans("PSTP: after Comm this") == FAIL )
    return FAIL;
  if( ParentGrid.CheckForNans("PSTP: after Comm pg") == FAIL )
    return FAIL;

  ParentGrid.DebugCheck("ProjectSolutionToParentGrid (Parent, after)");
  //if( MHDAnis(" PSTPG: End" ) == FAIL ) return FAIL;
  //if( ParentGrid.MHDAnis(" PSTPG: Parent End" ) == FAIL ) return FAIL;
 
  //ParentGrid.MHD_GROWTH("PSTPG end Parent");
  //this->MHD_GROWTH("PSTPG end this");

  return SUCCESS;
}
