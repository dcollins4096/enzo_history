#ifdef MHDF
//
// C wrapper to Magnetic Field centering routines.
//
// Uses MagneticField[3] to calculate CenteredB[3].
// Switches between one of 4 methods.  Only two of them are actually useful.
// Method determined by the Enzo global variable MHD_CenteringMethod.

//
// Options for MHD_CenteringMethod
// 1: Directionally split averaging. Bxc(i,j,k) = 0.5(Bxf(i+1,j,k) + Bxf(i,j,k)) etc.
//    This method also applies a correction to the Energy, where the Magnetic poriton of                      
//    the total energy is replaced by the newly centered total energy.  This option causes total energy
//    conservation to be violated.  Its usefullness hasn't been fully examined as of this writing.
// 2: Same as above without Energy correction.  Conserves energy.   (this is what you should use.)
// 3: Volumetric average.  The volume integral of a 3d quadratic reconstruction is used.
//    See that source code for more details.  Alternatively, see Balsara (2001) (refernce below.)
// 4: stores only dB/dt in the centered field.  A debugging option only.  If this hasn't been removed
//    by the time the MHD version of the code has been released, notify David Collins.
//
// 5: None.  
//
// 6: Non3d. This option is for rank != 3 simulations, whereing Bcy ==
//    Bfy (for 1d,) Bcz == Bfz (for 2d and 1d.)  This is done directly
//    in this routine (while the others are all done in FOrtran calls).
//    In this routine, rank!= 3 changes are taken care of by shifting
//    the offset along the 'homogenous' axis to zero, so those fields
//    are averaged with themselves.
//    
// Note that if MHD_DivB != MHD_DivB_Balsara, MHD_CenteringMethod is forced to 5 = MHD_none.
// This is due to the fact that the cell centered field is the only field used in this case, so
// it's evolution needs to NOT be overwritten by the face centered field

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "PPML.h"

extern "C" void 
FORTRAN_NAME(center_magnetic_field)(float *bxf, float *byf, float *bzf, 
				    float *bxc, float *byc, float *bzc, 
				    float *energy,
				    FLOAT *dx, FLOAT *dy, FLOAT *dz,
				    int *idim, int *jdim, int *kdim,
				    int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
				    int *method);


extern "C" void 
FORTRAN_NAME(mhd_center_volumetric)(float *bxface, float *byface, float *bzface, 
				    float *bxc, float *byc, float *bzc, float *energy,
				    FLOAT *dx, FLOAT *dy, FLOAT *dz,
				    int *idim, int *jdim, int *kdim,
				    int *method);



int grid::CenterMagneticField(int * Start, int * End){

  //Setting default start & end indicies.  If there's a slicker way to do this, I'm all ears, but I think
  // default setting in C++ must be a static variable.

  if( Start == NULL ) Start = this->GridStartIndex;
  if( End   == NULL )   End = this->GridEndIndex;

  //If this processor (MyProcessorNumber) doesn't have the data belonging to this grid (ProcessorNumber) 
  //leave quietly.
  if( ProcessorNumber != MyProcessorNumber )
    return SUCCESS;

  //CenterMagneticField is called from Set Boundary Conditions.  Set Boundary Conditions is sometimes called from
  //SetAccelerationBoundary, which does some pointer juggling in order to apply the proper boundary conditions 
  //to the acceleration field.

  if( AccelerationHack == TRUE )
    return SUCCESS;

  IndexPointerMap ind;
  IdentifyPhysicalQuantities_2( ind );

  //
  // Ensure that the centering method is appropriate:
  //  For non-CT runs, it needs to be OFF (because you evolve CenteredB
  //  For GridRank < 3, use directy copy on the flat dimensions, simple averaging on non-flat.

  if( MHD_DivB == NoDivB ){
    if( MHD_CenteringMethod != MHD_none ){
      fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      fprintf(stderr," !!! Warning: Incompatable DivB (%d) and Centering (%d) method.  Will cause errors. !!!\n",
	      (int) MHD_DivB, (int) MHD_CenteringMethod);
      fprintf(stderr," !!! Chaning the parameter MHD_CenteringMethod to MHD_none (5)\n");
      fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");  
      fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");      
      fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      MHD_CenteringMethod = MHD_none;
    }
  }//divb
  if(  GridRank < 3 && MHD_CenteringMethod != MHD_Non3d && MHD_CenteringMethod != MHD_none ){
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," Warning: GridRank < 3.  Using Non3d Centering method. (6)(Simple average on non-flat axis)\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    MHD_CenteringMethod = MHD_Non3d;
  }

  int method;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;

  //Variables for the Non-3d centering.
  int i,j,k,dim, indexC, indexB1, indexB2;
  int Offset[3]={1,MagneticDims[1][0], MagneticDims[2][1]*MagneticDims[2][0]};

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
  
  //This switch came about because the volumetric average method came much later than the others,
  //and is significantly more complicated.
  method = 3;
  switch(MHD_CenteringMethod){

  case MHD_SplitWithEng:
    method = 1;
    FORTRAN_NAME(center_magnetic_field)
      (MagneticField[0], MagneticField[1], MagneticField[2],
       ind.CenteredB[0], ind.CenteredB[1], ind.CenteredB[2], BaryonField[TENum],
       CellWidth[0], CellWidth[1], CellWidth[2],
       GridDimension, GridDimension + 1, GridDimension +2,
       Start, End,
       Start+1, End+1,
       Start+2, End+2,
       &method);
    break;

  case MHD_dtOnly:
    //This method is for development purposes only.
    method = 3;
    FORTRAN_NAME(center_magnetic_field)
      (MagneticField[0], MagneticField[1], MagneticField[2],
       ind.CenteredB[0], ind.CenteredB[1], ind.CenteredB[2], BaryonField[TENum],
       CellWidth[0], CellWidth[1], CellWidth[2],
       GridDimension, GridDimension + 1, GridDimension +2,
       Start, End,
       Start+1, End+1,
       Start+2, End+2,
       &method );
    break;

  case MHD_Split:
    //fprintf(stderr," = Old Mehtod ==========================================================================\n");
    method = 2;
    FORTRAN_NAME(center_magnetic_field)
      (MagneticField[0], MagneticField[1], MagneticField[2],
       ind.CenteredB[0], ind.CenteredB[1], ind.CenteredB[2], BaryonField[TENum],
       CellWidth[0], CellWidth[1], CellWidth[2],
       GridDimension, GridDimension + 1, GridDimension +2,
       Start, End,
       Start+1, End+1,
       Start+2, End+2,
       &method);
    break;

  case MHD_Volumetric:
    method = 666;
    //fprintf(stderr,"= the clap ==========================================================================\n");
    FORTRAN_NAME(mhd_center_volumetric)
      (MagneticField[0], MagneticField[1], MagneticField[2],
       ind.CenteredB[0], ind.CenteredB[1], ind.CenteredB[2], BaryonField[TENum],
       CellWidth[0], CellWidth[1], CellWidth[2],
       GridDimension, GridDimension + 1, GridDimension +2,
       &method);

      
    break;

    //You don't *have* to use a centering method, especially if you're not doing CT.
  case MHD_none:
    fprintf(stderr,"No centering method\n");
    break;

  case MHD_Non3d:
    
    //Offset is the distance in memory between i&i+1, j&j+1, or k&k+1, depending on the magnetic field in question
    //For 1d or 2d runs, j+1 or k+1 may not be updated for the magnetic field, so a direct copy is in order.
    //For ease of coding, this is done by averaging Bx[i,j,k] with itself.

    if( GridRank < 3 ) Offset[2] = 0;
    if( GridRank < 2 ) Offset[1] = 0;

    for( k=Start[2]; k<=End[2]; k++)
      for( j=Start[1]; j<=End[1]; j++)
	for( i=Start[0]; i<= End[0]; i++)
	  for(dim=0;dim<3;dim++){
	    indexC = i + GridDimension[0]*(j + GridDimension[1]*k);
	    indexB1= i + MagneticDims[dim][0]*(j+MagneticDims[dim][1]*k);
	    indexB2= i + MagneticDims[dim][0]*(j+MagneticDims[dim][1]*k) + Offset[dim];
	    ind.CenteredB[dim][indexC] =  0.5 *(MagneticField[dim][indexB1] + MagneticField[dim][indexB2]);

	  }

    break;
    
  default:
    fprintf(stderr, "Please select MHD Centering Method: \n");
    fprintf(stderr, "      MHD_Centering = 1: Directionally split, with Energy Correciton (probably bad)\n");
    fprintf(stderr, "                    = 2: Directionally split, no energy correction (better)\n");
    fprintf(stderr, "                    = 3: Volumetric average (best)\n");
    fprintf(stderr, "      For more info on each method, See Balsara, J. Comp. Phys, 174, 614-648 (2001) [method 2&3] \n");
    fprintf(stderr, "      and Balsara & Spicer, J. Comp Phys, 149, 270-292 (1999) \n");
    return FAIL;
    break;

  }//switch

  return SUCCESS;
}
#endif //! MHDF
