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
#include "CosmologyParameters.h"

#ifdef HAOXU
 int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
#endif

extern "C" void FORTRAN_NAME(curl_of_e)(float *bx, float *by, float *bz,
                                        float *ex, float *ey, float *ez,
                                        float *dx, float *dy, float *dz,
                                        int *idim, int *jdim, int *kdim,
                                        int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
                                        float *dt, MHD_Centering *method);


int grid::MHD_UpdateMagneticField(int level, LevelHierarchyEntry * NextLevel){
  //fprintf(stderr,"kludge: probably need to fix the ars for MHD_UpdateMagneticField\n");

   if(MyProcessorNumber != ProcessorNumber || MHD_Used != TRUE)
    return SUCCESS;

  //if we're projecting the magnetic field, this routine isn't necessary.
  if( MHD_ProjectE != TRUE )
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, field;
  float dtUsed;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
  
  //dtUsed is a hack-- in the electric field is NOT E, but E/dt.  This is to ensure
  // that the field is properly averaged over subgrid timesteps of potentially uneven lengths.
  // I left in the option of dt != 1 in order to maintain the versitility of the code.
  // And, because I was lazy.

  dtUsed = 1.0;

  //Magnetic Field needs to be re-updated with the New And Improved (projected from subgrid)
  //Electric field.  THis step is equivalent to Flux Correction, but Enzo doesn't deal with
  //face centered fields well enough to work properly.
  //My choices were to re-write a big chunk of code, or to hack.  I chose to hack.

#ifdef HAOXU
  float * PlaceHolder;
  for(field=0;field<3; field++){
  if(OldMagneticField[field]!=NULL){
    PlaceHolder=MagneticField[field];
    MagneticField[field]=OldMagneticField[field];
    OldMagneticField[field]=PlaceHolder;
  }else{
    PlaceHolder=MagneticField[field];
    MagneticField[field]= new float[MagneticSize[field]];
    OldMagneticField[field]=PlaceHolder;
  }
  }
#else

#ifndef ATHENA  
  float * PlaceHolder;
  for(field=0;field<3; field++){
    PlaceHolder=MagneticField[field];
    MagneticField[field]=OldMagneticField[field];
    OldMagneticField[field]=PlaceHolder;
  }
#endif //ATHENA

#endif //HAOXU

  /* dcc 9/1/06: This has all be updated by the new MHD_Curl routine, which actuall does the proper thing.
     As of 9/1/06, it was un-tested, so watch it.
     You'll know if it worked if you can do a several level run without divergence.
     The only real open question is the fate of OldMagneticField from this point.  I believe that it isn't used
     and that's why I got away with this juggle.  
     //Magnetic Field needs to be re-updated with the New And Improved (projected from subgrid)
     //Electric field.  THis step is equivalent to Flux Correction, but Enzo doesn't deal with
     //face centered fields well enough to work properly.
     //My choices were to re-write a big chunk of code, or to hack.  I chose to hack.
     
     float * PlaceHolder;
     for(field=0;field<3; field++){
     PlaceHolder=MagneticField[field];
     MagneticField[field]=OldMagneticField[field];
     OldMagneticField[field]=PlaceHolder;
     }
  */

  
  MHD_divbmethod AcceptableMethods[6] = 
    {MHD_DivB_Balsara, MHD_DivB_RJ, MHD_DivB_Athena,
     MHD_DivB_Athena_Switch, MHD_DivB_Athena_LF, MHD_DivB_Athena_Balsara};

  int ReUpdate = FALSE;
  if( HydroMethod != MHD_None )
    for(int meth=0;meth<6;meth++)
      if( MHD_DivB == AcceptableMethods[meth] ) ReUpdate = TRUE;

  if( ReUpdate == TRUE ){

    Pout("   UMF Enter UpdateMagneticField");

#ifdef HAOXU
   /* If using comoving coordinates, multiply dx by a(n+1/2). */
   // no longer used, multiply dx by a(n+1/2). inside solver 
    FLOAT a=1.0,dadt;
 
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time-0.5*dtFixed, &a, &dadt)
          == FAIL) {
        fprintf(stderr, "Error in CsomologyComputeExpansionFactors.\n");
        return FAIL;
      }

    /* Create a cell width array to pass (and convert to absolute coords). */
    //this is no longer using, now the flux is divided by a inside solve    
    float *CellWidthTemp[MAX_DIMENSION];
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      CellWidthTemp[dim] = new float[GridDimension[dim]];
      for (int i = 0; i < GridDimension[dim]; i++)
        if (dim < GridRank)
          CellWidthTemp[dim][i] = float(a*CellWidth[dim][i]);
        else
          CellWidthTemp[dim][i] = 1.0;
    }

     FORTRAN_NAME(curl_of_e)(MagneticField[0], MagneticField[1], MagneticField[2],
                          ElectricField[0], ElectricField[1], ElectricField[2],
                          CellWidth[0], CellWidth[1], CellWidth[2],
                          GridDimension, GridDimension +1, GridDimension +2,
                          GridStartIndex, GridEndIndex,
                          GridStartIndex+1, GridEndIndex+1,
                          GridStartIndex+2, GridEndIndex+2,
                          &dtUsed, &MHD_CenteringMethod);


#else

#ifdef ATHENA
     //<dbg> test, db21
     //MHD_Curl( GridStartIndex, GridEndIndex, 1);
     int CurlStart[3] = {0,0,0}, CurlEnd[3] = {GridDimension[0]-1,GridDimension[1]-1,GridDimension[2]-1};
     MHD_Curl( CurlStart,CurlEnd, 1);
     //</dbg>

      //if( this->MHD_Curl(GridStartIndex, GridEndIndex,2) == FAIL )
      //{fprintf(stderr," error occored in MHD_Curl\n"); return FAIL;}
#else

    FORTRAN_NAME(curl_of_e)(MagneticField[0], MagneticField[1], MagneticField[2],
			  ElectricField[0], ElectricField[1], ElectricField[2],
			  CellWidth[0], CellWidth[1], CellWidth[2],
			  GridDimension, GridDimension +1, GridDimension +2,
			  GridStartIndex, GridEndIndex,
			  GridStartIndex+1, GridEndIndex+1,
			  GridStartIndex+2, GridEndIndex+2,
			  &dtUsed, &MHD_CenteringMethod);
#endif //ATHENA
#endif

#ifdef HAOXU
    if (ComovingCoordinates && MHD_Equation == 2){
     int i;
     for(field=0;field<3;field++) 
     for(i=0;i<MagneticSize[field];i++)
     MagneticField[field][i]
      *=(1.0-0.25*dtFixed*dadt/a)/(1.0+0.25*dtFixed*dadt/a);
  }
#endif //HAOXU 
  
    if( this->CenterMagneticField() == FAIL ) {
      fprintf(stderr," error with UpdateMagneticField , second call\n");
      return FAIL;
    }

#ifdef HAOXU
    int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  if(MHD_Used == 1 && DualEnergyFormalism ==1){

     /* Compute the field size. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];


    for(int i=0;i<size;i++)
       BaryonField[TENum][i] = BaryonField[GENum][i]*BaryonField[DensNum][i]
 + 0.5*BaryonField[DensNum][i]*
             (BaryonField[Vel1Num][i]*BaryonField[Vel1Num][i]
              +BaryonField[Vel2Num][i]*BaryonField[Vel2Num][i]
             +BaryonField[Vel3Num][i]*BaryonField[Vel3Num][i])
             + 0.5*(CenteredB[0][i]*CenteredB[0][i]+CenteredB[1][i]
             *CenteredB[1][i]+CenteredB[2][i]*CenteredB[2][i]);
   }
#endif
  
  
  //Also update AvgElectricField, if this is a subgrid.
  if( level > 0 ){
    
    
    for(int field=0;field<3;field++){
      
      if(AvgElectricField[field] == NULL){
	fprintf(stderr,"AvgField null\n");
	return FAIL;
      }
      for(int i=0;i<ElectricSize[field];i++){
	AvgElectricField[field][i] += ElectricField[field][i];
      }
    }//field
  }//level>0
  }//DivB == Balsara

  Pout("   End of UpdateMagneticField");
  return SUCCESS;
}

