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
/  GRID CLASS (SET EXTERNAL BOUNDARY VALUES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

// Copy the current baryon fields to the old baryon fields/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
void dump(float * A,int nx,int ny,int nz,int LongAxis, char * filename);


int grid::SetExternalBoundaryValues(ExternalBoundary *Exterior)
{
  int dim, field;
  //dcc Shit.
  int LongAxis = GridDimension[0];
  LongAxis = LongAxis > GridDimension[1] ? LongAxis : GridDimension[1];
  LongAxis = LongAxis > GridDimension[2] ? LongAxis : GridDimension[2];
  LongAxis = LongAxis == GridDimension[0] ? 0 : LongAxis;
  LongAxis = LongAxis == GridDimension[1] ? 1 : LongAxis;
  LongAxis = LongAxis == GridDimension[2] ? 2 : LongAxis;


  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  //if( this->MHDAnis(" +++ beginning of SetExtBdryVal ") == FAIL ) return FAIL;

  /* For Wave Pool problem, compute the new inflow boundary conditions. */

  if (ProblemType == 2)
    if (Exterior->SetWavePoolBoundary(Time) == FAIL) {
      fprintf(stderr, "Error in exterior->SetWavePoolBoundary.\n");
      return FAIL;
    }

  /* For Shock Pool problem, compute the new inflow boundary conditions. */

  if (ProblemType == 3)
    if (Exterior->SetShockPoolBoundary(Time) == FAIL) {
      fprintf(stderr, "Error in exterior->SetShockPoolBoundary.\n");
      return FAIL;
    }
   
  /* For the DoubleMach problem, set the bew inflow boundary conditions. */

  if (ProblemType == 4)
    if (Exterior->SetDoubleMachBoundary(Time, CellLeftEdge[0], CellWidth[0]) 
	== FAIL) {
      fprintf(stderr, "Error in exterior->SetDoubleMachBoundary.\n");
      return FAIL;
    }

  /* Compute offset from corner of domain. */

  int GridOffset[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    if (dim < GridRank)
      GridOffset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/
			     CellWidth[dim][0]);
    else
      GridOffset[dim] = 0;

  /* loop through fields, setting each */


  for (field = 0; field < NumberOfBaryonFields; field++) {
    
    if (BaryonField[field] == NULL) {
      fprintf(stderr, "Baryon field missing.\n");
      return FAIL;
    }

#ifdef OOC_BOUNDARY
    ExternalBoundaryField = field;
#endif
    if (Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
				      GridStartIndex, GridEndIndex,
				      BaryonField[field], FieldType[field])
	== FAIL) {
      fprintf(stderr, "Error in Exterior->SetExternalBoundary.\n");
      return FAIL;
    }


  }

#ifdef EMF_BOUNDARY
     if(EMF_Boundary==1){
     if(MHD_Used==1)
       for(field=0;field<3;field++)
       if (Exterior->SetExternalBoundary_EMF(GridRank, GridDimension, GridOffset,
                                      GridStartIndex, GridEndIndex,
                                      ElectricField[field], field)
        == FAIL) {
      fprintf(stderr, "Error in Exterior->SetExternalBoundary_EMF.\n");
      return FAIL;
    }
      }else{
        /*If there's a magnetic field, set it as well.
    It is still unclear if this is a valid way to do things:
    The Pseudo-Vecor nature of B poses a problem that I haven't sorted out all the way.
    Currently, it's being set as if it were a Plain vector.*/
  //if( this->MHDAnis(" +++ SEBV: After BarryonFields ") == FAIL ) return FAIL;

  if(MHD_Used && NumberOfBaryonFields > 0)
    {

      //A switch over the centering method.  If face centerd variables (i.e. balsara method) the SetMagneticboundary should be used.
      // If cell centered, Exterior->SetExternalBoundary should be used, with magnetic field treated as the corresponding Velocity.
      // Again, this isn't entirely correct.  If you see this comment, contact David Collins.

      switch( MHD_DivB){

      case MHD_DivB_none:

      case MHD_DivB_Balsara:



        if(Exterior->SetMagneticBoundary(GridRank, GridDimension, GridOffset,
                                         GridStartIndex, GridEndIndex,
                                         MagneticField[0],MagneticField1type) == FAIL){
          fprintf(stderr, "Error in Exterior->SetMagneticBoundary, B1\n");
            return FAIL;
        }
        //if(MHD_Diagnose(" Set External: After Bx")==FAIL) return FAIL;
        if(Exterior->SetMagneticBoundary(GridRank, GridDimension, GridOffset,
                                         GridStartIndex, GridEndIndex,
                                         MagneticField[1],MagneticField2type) == FAIL){
          fprintf(stderr, "Error in Exterior->SetMagneticBoundary, B2\n");
            return FAIL;
        }
        //if(MHD_Diagnose(" Set External: After By")==FAIL) return FAIL;
        if(Exterior->SetMagneticBoundary(GridRank, GridDimension, GridOffset,
                                         GridStartIndex, GridEndIndex,
                                         MagneticField[2],MagneticField3type) == FAIL){
          fprintf(stderr, "Error in Exterior->SetMagneticBoundary, B3\n");
          return FAIL;
           }
        //if(MHD_Diagnose(" Set External: After Bz")==FAIL) return FAIL;
        //if( this->MHDAnis(" +++ SEBV: After Bf ") == FAIL ) return FAIL;
        break;


      case MHD_DivB_Poisson:
        fprintf(stderr, "You need to install the poisson cleaner for MHD. I haven't done it yet.\n");
        return FAIL;
        break;

        //case MHD_DivB_none:
        //this has been moved so that the boundary set is called for debugging.
        break;
      } // divb switch

      if( Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
                                        GridStartIndex, GridEndIndex,
                                        CenteredB[0], Velocity1) == FAIL){
        fprintf( stderr, "Shit!  Something's wrong with the CenteredB[0] boundary\n");
        return FAIL;
      }

      if( Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
                                        GridStartIndex, GridEndIndex,
                                        CenteredB[1], Velocity2) == FAIL){
        fprintf( stderr, "Shit!  Something's wrong with the CenteredB[1] boundary\n");
        return FAIL;
      }

      if( Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
                                        GridStartIndex, GridEndIndex,
                                        CenteredB[2], Velocity3) == FAIL){
        fprintf( stderr, "Shit!  Something's wrong with the CenteredB[2] boundary\n");
        return FAIL;
      }
      //if( this->MHDAnis(" +++ SEBV: After Bc ") == FAIL ) return FAIL;

       }// if(MHDUsed)
    
      }// EMF_Boundary

#else // EMF_BOUNDARY
  /*If there's a magnetic field, set it as well.
    It is still unclear if this is a valid way to do things: 
    The Pseudo-Vecor nature of B poses a problem that I haven't sorted out all the way.
    Currently, it's being set as if it were a Plain vector.*/
  //if( this->MHDAnis(" +++ SEBV: After BarryonFields ") == FAIL ) return FAIL;

  if(MHD_Used && NumberOfBaryonFields > 0)
    {

      //A switch over the centering method.  If face centerd variables (i.e. balsara method) the SetMagneticboundary should be used.
      // If cell centered, Exterior->SetExternalBoundary should be used, with magnetic field treated as the corresponding Velocity.
      // Again, this isn't entirely correct.  If you see this comment, contact David Collins.

      switch( MHD_DivB){

      case MHD_DivB_none:

      case MHD_DivB_Balsara:


      
	if(Exterior->SetMagneticBoundary(GridRank, GridDimension, GridOffset,
					 GridStartIndex, GridEndIndex,
					 MagneticField[0],MagneticField1type) == FAIL){
	  fprintf(stderr, "Error in Exterior->SetMagneticBoundary, B1\n");
	    return FAIL;
	}
	//if(MHD_Diagnose(" Set External: After Bx")==FAIL) return FAIL;
	if(Exterior->SetMagneticBoundary(GridRank, GridDimension, GridOffset,
					 GridStartIndex, GridEndIndex,
					 MagneticField[1],MagneticField2type) == FAIL){
	  fprintf(stderr, "Error in Exterior->SetMagneticBoundary, B2\n");
	    return FAIL;
	}
	//if(MHD_Diagnose(" Set External: After By")==FAIL) return FAIL;
	if(Exterior->SetMagneticBoundary(GridRank, GridDimension, GridOffset,
					 GridStartIndex, GridEndIndex,
					 MagneticField[2],MagneticField3type) == FAIL){
	  fprintf(stderr, "Error in Exterior->SetMagneticBoundary, B3\n");
	  return FAIL;
	}
	//if(MHD_Diagnose(" Set External: After Bz")==FAIL) return FAIL; 
	//if( this->MHDAnis(" +++ SEBV: After Bf ") == FAIL ) return FAIL;
	break;

	  
      case MHD_DivB_Poisson:
	fprintf(stderr, "You need to install the poisson cleaner for MHD. I haven't done it yet.\n");
	return FAIL;
	break;
	
	//case MHD_DivB_none:
	//this has been moved so that the boundary set is called for debugging.
	break;
      } // divb switch

      if( Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
					GridStartIndex, GridEndIndex,
					CenteredB[0], Velocity1) == FAIL){
	fprintf( stderr, "Shit!  Something's wrong with the CenteredB[0] boundary\n");
	return FAIL;
      }
      
      if( Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
					GridStartIndex, GridEndIndex,
					CenteredB[1], Velocity2) == FAIL){
	fprintf( stderr, "Shit!  Something's wrong with the CenteredB[1] boundary\n");
	return FAIL;
      }

      if( Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
					GridStartIndex, GridEndIndex,
					CenteredB[2], Velocity3) == FAIL){
	fprintf( stderr, "Shit!  Something's wrong with the CenteredB[2] boundary\n");
	return FAIL;
      }
      //if( this->MHDAnis(" +++ SEBV: After Bc ") == FAIL ) return FAIL;
          
    }// if(MHDUsed)
#endif //EMF_BOUNDARY
  /* Now we handle the particles (if any). */



  if (NumberOfParticles > 0)

    if (Exterior->SetExternalBoundaryParticles(GridRank, NumberOfParticles,
					       ParticlePosition,
					       ParticleVelocity) == FAIL) {
      fprintf(stderr, "Error in Exterior->SetExternalBoundaryParticles.\n");
      return FAIL;
    }



  //if( this->MHDAnis(" +++ SEBV: end ") == FAIL ) return FAIL;
  
  return SUCCESS;

}
