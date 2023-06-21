#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "pout.h"

#ifdef ATHENA

#ifdef HAOXU
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
#endif

#define DEFINE_STORAGE  //used to switch between declaration and storage vs. external linking.
#include "MHD_Athena.h"


int grid::MHD_Athena(int CycleNumber, int level, int grid,
		     int NumberOfSubgrids, fluxes * SubgridFluxes[],
		     float * RandomForcingNormalization, float TopGridTimeStep){

  // <dbg>
  int dbgOut = 1;

  //<dbg> I'm not sure the best way to install the CellWidth necessary for the 
  //      the dissipation in the electric field routine.  I've
  //      algebra-ed it out for now, but each derivative should get
  //      its own CellWidth.
  if( MHD_ElectricRecon == 1 )
    if( ( fabs( (CellWidth[0][0] -CellWidth[1][0])/CellWidth[0][0] ) >0.01 ) ||
	( fabs( (CellWidth[0][0] -CellWidth[2][0])/CellWidth[2][0] ) >0.01 && GridRank == 3) ){
      fprintf(stderr,"You can't currently run an MHD Enzo Athena simulation with differeing cell widths\n");
      fprintf(stderr,"It comes from the implimentation of the FastestWaveSpeed, necessary for the\n");
      fprintf(stderr,"lax-Fredrichs formalism in the electric field computation.\n");
      fprintf(stderr,"fix it or consult with mr collins.\n");
      fprintf(stderr,"Widths: %f %f %f\n", CellWidth[0][0], CellWidth[1][0], CellWidth[2][0]);
      return FAIL;
    }
  
  //JBMEM_MESSAGE(MyProcessorNumber,"jb: Enter Ath");
  // Loop & other ancillary Variables
  int size = 1, FluxSize,field,dim,cell, index,index2, i,j,k, SubCycle;
  int result;


#ifdef HAOXU
  FLOAT a[4]={1.0,0.0,1.0,0.0};
  if(ComovingCoordinates==1){
    if (CosmologyComputeExpansionFactor(Time, &a[0],&a[1])
        == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
    if (CosmologyComputeExpansionFactor(Time+(FLOAT)0.5*dtFixed, &a[2],&a[3])
        == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
}
 
 float temp_sqrt_a=sqrt(a[0]);
#endif


  float * Fluxes[3];
  float dT = dtFixed; 
  float SubRandomForcingNormalization;  //norm for the subcycle.

  //<dbg> variables for debugging.
  float t1, t2;
  char basename[20];
  char FlatteningName[3][20];
  for( dim=0; dim<3; dim++){
    sprintf(FlatteningName[dim], "Flattening%d", dim);
  }

  //<dbg>
  char SolverProblemName[3][20];
  for( dim=0;dim<3;dim++)
    sprintf(SolverProblemName[dim], "SolverProblem%d",dim);
  //</dbg>

  int WriteInThisNumber, MyGridNumber = 1;


#ifdef HAOXU  // B=B/sqrt(a)
if(ComovingCoordinates==1){
for( i=0; i<size; i++){
BaryonField[Eeng][i] -= 0.5*(CenteredB[0][i]*CenteredB[0][i]+
                                    CenteredB[1][i]*CenteredB[1][i]+
                                    CenteredB[2][i]*CenteredB[2][i]);
CenteredB[0][i] /= temp_sqrt_a;
CenteredB[1][i] /= temp_sqrt_a;
CenteredB[2][i] /= temp_sqrt_a;
BaryonField[Eeng][i] += 0.5*(CenteredB[0][i]*CenteredB[0][i]+
                                    CenteredB[1][i]*CenteredB[1][i]+
                                    CenteredB[2][i]*CenteredB[2][i]);
                                                                                                                                                             
OldBaryonField[Eeng][i] -= 0.5*(OldCenteredB[0][i]*OldCenteredB[0][i]+
                                    OldCenteredB[1][i]*OldCenteredB[1][i]+
                                    OldCenteredB[2][i]*OldCenteredB[2][i]);
OldCenteredB[0][i] /= temp_sqrt_a;
OldCenteredB[1][i] /= temp_sqrt_a;
OldCenteredB[2][i] /= temp_sqrt_a;
OldBaryonField[Eeng][i] += 0.5*(OldCenteredB[0][i]*OldCenteredB[0][i]+
                                    OldCenteredB[1][i]*OldCenteredB[1][i]+
                                    OldCenteredB[2][i]*OldCenteredB[2][i]);
}
for(field=0;field<3;field++)
for( i=0; i<MagneticSize[field]; i++){
MagneticField[field][i] /= temp_sqrt_a;
OldMagneticField[field][i] /= temp_sqrt_a;
}

                                                                                                                                                             
}
                                                                                                                                                             
#endif


  //</dbg>
  //
  // Allocate Solver Fluxes.
  //

  //  To Do:
  //  1.) Set up parameters.
  //  2.) Allocate Fluxes
  //  3.) Allocate Electric Field
  //  4.) Allocate Flattening Field

  // 1.)

  if( this->MHD_AthenaSetup() == FAIL )
    {fprintf(stderr,"MHD_Athena: Athena Setup failed.\n"); return FAIL;}

  // 2.) 
  size = 1;
  for( dim=0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  FluxSize = 1;
  for( dim=0;dim<GridRank;dim++)
    FluxSize *= GridDimension[dim];

  FluxSize = size * NumberOfMHDFluxes;

  //note that this is intentionally 3, not grid rank--this is MHD.
  for( dim=0; dim < 3; dim++){
    Fluxes[dim] = new float[FluxSize];
    for( cell=0; cell< FluxSize; cell++)
      Fluxes[dim][cell] = 0.0;
  }

  // 3.)
  for( field=0;field<3;field++){
    if( ElectricField[field] == NULL )
      ElectricField[field] = new float[ ElectricSize[field] ];
    //IBM doesn't play nice with memory problems, so we should double check
    if( ElectricField[field] == NULL ){
      fprintf(stderr,"Grid_MHD_Athena.C: Memory error.\n");
      return FAIL;
    }
  }//field


  // 4.) 
  switch( MHD_Flattening ){

  case 1:
  case 3:
    FlatteningField = new float[size*GridRank];
    ShockDirection = new int[size*GridRank];
    for( i=0;i<size*GridRank;i++)
      FlatteningField[i] = 30000 ;//1.2345; //<dbg> should be 0
    break;
  case 0:
  default:
    //If off, don't allocate.
    break;
  }

  //<dbg>
  //keeping track of how badly the roe solver sucks.
  CheckForSolverProblems = FALSE;
  if( CheckForSolverProblems == TRUE ){
    for( dim=0; dim<3; dim++){
      if( SolverProblems[dim] != NULL ){
	delete[] SolverProblems[dim];
	fprintf(stderr,"Grid_MHD_Athena: SolverProblems != NULL.  Examine!\n");
      }
      SolverProblems[dim] = new float[size];
      for(i=0;i<size;i++)
	SolverProblems[dim][i] = 3;
    }
  }//
  

  // <dbg>
  // Something for debugging
  //

  for( dim=0; dim<GridRank; dim++){

  for(k=GridStartIndex[2]; k<=GridEndIndex[2]; k++)
    for(j=GridStartIndex[1]; j<=GridEndIndex[1]; j++)
      for(i=GridStartIndex[0]; i<=GridEndIndex[0]; i++){
	index = i + GridDimension[0]*( j + GridDimension[1]*k);

	for(field = 0; field < NumberOfMHDFluxes; field++){
	  index = index0(i,j,k)*NumberOfMHDFluxes + field;
	  Fluxes[dim][ index ] = 12345.666;//10*i + field;
	  if( index >= FluxSize ){ fprintf(stderr,"shit.\n"); return FAIL;}
	}//field
      }//cell
  }//dim

  // </dbg>

  //
  //
  // MagneticField /= sqrt(a)
  // OldMagneticField /= sqrt(a)
  // 

  //
  //If momentum is interpolated (instead of velocity) make the velocity fields into momentum.
  //

  if( MHD_ReconField[0] == 1 ){
    for( field = 0; field<3; field++){
      for( i = 0; i<size; i++){
	BaryonField[ Ev[field] ][i] *= BaryonField[ Eden ][i];
	OldBaryonField[ Ev[field] ][i] *= OldBaryonField[ Eden ][i];
	
      }
    }}//momentum recon

  //
  // Set up the State field.
  //

  State[Sden] = BaryonField[Eden];
  if( EquationOfState == 0)
    State[Seng] = BaryonField[Eeng];

  //
  // INTEGRATE!!!
  //
  
  //Number of subcycles selected
  int nSubCycle =  (MHD_Recon[1] == -1 ) ? ( (MHD_Recon[0] == -1 ) ? 0 : 1 ):2;
  if( nSubCycle == 1 ) fprintf(stderr, "Hey!:  only one sub-cycle.\n");
  if( nSubCycle == 0 ) fprintf(stderr, "KLUDGE!!! NO DAMNED SOLVER!!!\n");
  for(SubCycle = 0; SubCycle < nSubCycle; SubCycle++){


    if( dbgOut > 0 ){ fprintf(stderr,"subcycle %d \n", SubCycle);}

    //Set up cycle dependant parameters.

    dT = ( ( SubCycle == 0 && nSubCycle != 1 ) ) ? dtFixed*0.5: dtFixed; 

    if( TopGridTimeStep > tiny_number )
      SubRandomForcingNormalization = (*RandomForcingNormalization) * dT / TopGridTimeStep;
    else
      SubRandomForcingNormalization = 0.0;
    for(dim=0;dim<GridRank;dim++){
      AthStart[dim] = 0;
      AthEnd[dim] = GridDimension[dim]-1;
      for(int i=0;i<=SubCycle; i++){
	AthStart[dim] += GhostZones[i];
	AthEnd[dim]   -= GhostZones[i];
      }

    }
    if( dbgOut > 2 ){
      fprintf(stderr, "ath start: %d %d %d\n",AthStart[0],AthStart[1],AthStart[2]);
      fprintf(stderr, "ath end  : %d %d %d\n",AthEnd[0],AthEnd[1],AthEnd[2]);
      fprintf(stderr, "ath gdim : %d %d %d\n",
	      GridDimension[0],GridDimension[1],GridDimension[2]);
    }



    //
    // Reconstruction options: Energy (0), Pressure(1), TotalPressure(2), Enthalpy (3).
    // Note-- since pressure isn't evolved, it must be recomputed every step.
    //


    //
    // actually do the solver steps:
    // -1.0) Gravitational Pre-Step
    // 0.) Compute pressure, if needed
    // 1.) Set flattening field
    // 1.5.) Compute total pressure or enthalpy, as needed.
    // 2.) Compute Fluxes
    // 3.) Difference fluxes
    // 4.) Compute Electric Field
    // 5.) dB/dt = Curl( E )
    // 6.) Center magnetic field
    // 7.) Fill SubgridFluxes arrays
    // 8.) IO for debugging: Fluxes
    // 9.) IO for debugging: Fields

    //-1.0.) Note that this is a pre-step only, so that the gravitational acceleration
    //       is incorporated in the computation of the Left and Right states.
    //       Thus, it doesn't need to be added to OldBaryonField, only BaryonField

    if( 0== 0){
      if (SelfGravity || UniformGravity || PointSourceGravity) {

	for( i=0; i<size; i++ ){
	  //<dbg>
	  //for( dim=0; dim<GridRank; dim++)
	  //AccelerationField[dim][size] = 0;
	  //</dbg>
	  if( EquationOfState == 0 )
	    for(dim=0;dim<GridRank;dim++)
	      BaryonField[ Eeng ][i] += 0.5 * dT * 
		AccelerationField[dim][i] * BaryonField[ Ev[dim] ][i]* BaryonField[ Eden ][i];
	  
	  if( MHD_ReconField[0] == 0 ){
	    for( dim=0; dim<GridRank; dim++)
	      BaryonField[ Ev[dim] ][i] += 0.5*dT*AccelerationField[dim][i];
	  }else if( MHD_ReconField[0] == 1 ){
	    for( dim=0; dim<GridRank; dim++)
	      BaryonField[ Ev[dim] ][i] += 0.5*dT*AccelerationField[dim][i]*BaryonField[ Eden ][i];
	  }
	  
	}//i < size
      }//gravity
    }//hack out gravity kludge
    else{fprintf(stderr,"kludge! no gravity in pre-pre step.\n");}

    // 0.)
    if( MHD_ReconField[1] > 0 || MHD_Flattening != 0 ){
      if( MHD_Pressure == NULL ) MHD_Pressure = new float[size];
      
      if( EquationOfState == 0 )
	State[Seng] = MHD_Pressure;




      if (DualEnergyFormalism)
	result = this->ComputePressureDualEnergyFormalism(Time, MHD_Pressure);
      else
	result = this->ComputePressure(Time, MHD_Pressure);
      
      if (result == FAIL) 
	{fprintf(stderr, "Error in grid->ComputePressure, called in MHD_Athena\n");return FAIL;}
    }

    // 1.) 
    // Only done for non-constant reconstruction.



    if( MHD_Recon[SubCycle] != 0 && MHD_Flattening != 0 ){
      if( MHD_SetFlattening() == FAIL )
	{ fprintf(stderr, "ERROR: Set Flattening\n"); return FAIL;}
    }


    // 1.5)
    if( EquationOfState == 0 && MHD_ReconField[1] > 1){
	
      for( i=0; i<size; i++){
	MHD_Pressure[i] += 0.5*(CenteredB[0][i]*CenteredB[0][i]+
				CenteredB[1][i]*CenteredB[1][i]+
				CenteredB[2][i]*CenteredB[2][i]);
	
	//Enthalpy
	if( MHD_ReconField[1] > 2 ){
	  MHD_Pressure[i] += BaryonField[Eeng][i];
	  MHD_Pressure[i] /= BaryonField[Eden][i];
	}
      }//i
      if( MHD_ReconField[1] > 3 ){
	fprintf(stderr," Grid_MHD_Athena.C: MHD_ReconField[1] = %d isn't valid.\n", MHD_ReconField[1]);
	return FAIL;
      }

    }//EOS





    // 2.)
    if( dbgOut > 1 ) fprintf(stderr," Athena Fluxes: %d\n", MHD_Recon[SubCycle]);
    if( MHD_Fluxes( Fluxes, MHD_Recon[SubCycle], MHD_Riemann[SubCycle], MHD_DiffusionMethod[SubCycle], dT) == FAIL ) 
      { fprintf(stderr,"ERROR: FLuxesLinear\n"); return FAIL;}

    // 3.)
    if( dbgOut > 1 ) fprintf(stderr," Flux Difference\n");
    if( MHD_FluxDifference( dT, Fluxes, SubRandomForcingNormalization ) == FAIL) 
      {fprintf(stderr,"Error: FluxDifference"); return FAIL;}
    
    switch( MHD_DivB){
    case MHD_DivB_Balsara:
    case MHD_DivB_Athena:
    case MHD_DivB_RJ:
      
      // 4.)
      if( MHD_AthenaElectric( dT, Fluxes ) == FAIL )
	{fprintf(stderr,"Error in Athena Fluxes.\n"); return FAIL;}

#ifdef HAOXU  // ElectricField=ElectricField*/a
if(ComovingCoordinates==1){
   for(field=0;field<3;field++)
     for(i=0;i<ElectricSize[field];i++)
   ElectricField[field][i] /= a[2];  
}
#endif

      // 5.)
      //the last argument indicates MagneticField = OldMagneticField - Curl( ElectricField );
      if( MHD_Curl(AthStart,AthEnd, 2) == FAIL )
	{fprintf(stderr,"Error in Athena Curl.\n"); return FAIL;}

      // 6.)
      if( this->CenterMagneticField(AthStart, AthEnd) == FAIL ) 
	{fprintf(stderr," error with CenterMagneticField, first call \n");return FAIL;}
      break;
    default:
      break;

    }//divb switch

    // 7.) dcc: move this outside loop.
    //On the last step, save the fluxes for AMR
    if( SubCycle == nSubCycle-1 ){
      if( this->FillFluxes(NumberOfSubgrids, SubgridFluxes,Fluxes, dT) == FAIL)
	{fprintf(stderr," error with FillFluxes in MHD_Athena.\n"); return FAIL;}
    }//SubgridFluxes 

    //<dbg>
    MHD_Diagnose("cycle.");
    //</dbg>
    // 8.)

    //Write FLUXes.  
    //Kind of icky pointer juggling.
    // 31 32 33, 35 36 37

    int dbgPPML = 1; //for debugging PPML: shifts all fluxes 'ahead' by 1.
    for( dim=0; dim<3; dim++){
      WriteInThisNumber = 31 + dim + SubCycle*4;
      if( WriteInThisF(WriteInThisNumber) == TRUE) {
	fprintf(stderr, " dim %d SubCycle %d, WriteInThisNumber %d\n", dim, SubCycle, WriteInThisNumber);
	
	for(field=0; field<NumberOfBaryonFields; field++){
	  ActualBaryonField[field] = BaryonField[field];
	  BaryonField[field] = new float[size];
	  if( dbgPPML == 0 ){
	    for( i=0;i<size; i++)
	      BaryonField[field][i] = Fluxes[dim][ i * NumberOfMHDFluxes + MapEtoS[dim][field] ];
	  }else{
	    //For comparison with PPML, shift x fluxes ahead by 1.
	    for( i=0;i<size - NumberOfMHDFluxes; i++)
	      BaryonField[field][i] = Fluxes[dim][ (i+1) * NumberOfMHDFluxes + MapEtoS[dim][field] ];
	  }
	}
	for(field=0;field<3;field++){
	  ActualBaryonField[field + NumberOfBaryonFields] = CenteredB[field];
	  CenteredB[field] = new float[size];
	}
	if( dbgPPML == 0 ){
	  for( i=0;i<size;i++){
	    CenteredB[ BNum[dim][0] ][i]  = Fluxes[dim][ i*NumberOfMHDFluxes +  Sb[0] ];
	    CenteredB[ BNum[dim][1] ][i]  = Fluxes[dim][ i*NumberOfMHDFluxes +  Sb[1] ];
	  }
	}else{
	  for( i=0;i<size - NumberOfMHDFluxes;i++){
	    CenteredB[ BNum[dim][0] ][i]  = Fluxes[dim][ (i+1)*NumberOfMHDFluxes +  Sb[0] ];
	    CenteredB[ BNum[dim][1] ][i]  = Fluxes[dim][ (i+1)*NumberOfMHDFluxes +  Sb[1] ];
	  }
	}//ppml
      
      
	
	sprintf(basename, "data%d%d%d%d.grid",WriteInThisNumber,CycleNumber,level, grid);
		/*WriteInThisNumber, CycleNumber, level, grid);*/
	FILE *dummy = fopen(basename, "a");    
	if( this->WriteGrid(dummy, basename, MyGridNumber) == FAIL ){
	  fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
	  return FAIL;
	}
	fclose(dummy);
	for(field=0;field<NumberOfBaryonFields;field++){
	  delete BaryonField[field];
	  BaryonField[field] = ActualBaryonField[field];
	}
	for(field=0;field<3;field++){
	  delete CenteredB[field];
	  CenteredB[field] = ActualBaryonField[ field + NumberOfBaryonFields ];
	}
      }//if write in this
    }//dim
    
    // 9.)
    WriteInThisNumber = 34 + SubCycle*4;
    if( WriteInThisF(WriteInThisNumber) == TRUE ){
      sprintf(basename, "data%d%d%d%d.grid",WriteInThisNumber,CycleNumber, level, grid);

      FILE *dummy = fopen(basename, "a");    

      //<dbg>
      if( CheckForSolverProblems == TRUE ){
	for( dim =0; dim<3;dim++){
	  BaryonField[NumberOfBaryonFields] = SolverProblems[dim];
	  DataLabel[NumberOfBaryonFields] = SolverProblemName[dim];
	  DataUnits[NumberOfBaryonFields] = "choke on it.";
	  NumberOfBaryonFields += 1;
	}
      }
      //</dbg>


      //attach flattening field:
      switch( MHD_Flattening ){
      case 0: default: break;
      case 1:
      case 3:
	for( dim=0;dim<GridRank; dim++){

	  BaryonField[NumberOfBaryonFields] = new float[ size ];
	  DataLabel[NumberOfBaryonFields] = FlatteningName[dim];
	  DataUnits[NumberOfBaryonFields]  = "none";
	  for(k=0;k<GridDimension[2]; k++)
	    for(j=0;j<GridDimension[1]; j++)
	      for(i=0;i<GridDimension[0]; i++){
		index2 = dim + GridRank*(i + GridDimension[0]*(j + GridDimension[1]*k));
		index = i + GridDimension[0]*(j + GridDimension[1]*k);
		BaryonField[NumberOfBaryonFields][index] = FlatteningField[ index2 ];
	      }//i,j,k
	  NumberOfBaryonFields++;
	}//dim
      }//switch

      if( this->WriteGrid(dummy, basename, MyGridNumber) == FAIL ){
	fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
	return FAIL;
      }

      switch( MHD_Flattening ){
      case 0: default: break;  
      case 1:
      case 3:
	for( dim=0;dim<GridRank; dim++){
	  delete BaryonField[NumberOfBaryonFields--];
	}//dim
      }//switch

      //<dbg>
      //Just detach SolverProblems from BaryonField by decreasing NumberOfBaryonFields.
      //The field will be deleted later.
      if( CheckForSolverProblems == TRUE ){
	for(dim=0;dim<3;dim++){
	  if( SolverProblems[dim] != NULL ) {
	    NumberOfBaryonFields-- ;
	  }
	}
      }
      //</dbg>

      fclose(dummy);
    }  

  }//SubCycle.
  
  //If momentum is interpolated (instead of velocity) make the velocity fields into momentum.


  if( MHD_ReconField[0] == 1 ){
    for( field = 0; field<3; field++){
      for( i = 0; i<size; i++){
	BaryonField[ Ev[field] ][i] /= BaryonField[ Eden ][i];
	OldBaryonField[ Ev[field] ][i] /= BaryonField[ Eden ][i];

      }
    }    
  }//momentum recon


#ifdef HAOXU  // B=B*sqrt(a)
if(ComovingCoordinates==1){
for( i=0; i<size; i++){
BaryonField[Eeng][i] -= 0.5*(CenteredB[0][i]*CenteredB[0][i]+
                                    CenteredB[1][i]*CenteredB[1][i]+
                                    CenteredB[2][i]*CenteredB[2][i]);
CenteredB[0][i] *= temp_sqrt_a;
CenteredB[1][i] *= temp_sqrt_a;
CenteredB[2][i] *= temp_sqrt_a;
BaryonField[Eeng][i] += 0.5*(CenteredB[0][i]*CenteredB[0][i]+
                                    CenteredB[1][i]*CenteredB[1][i]+
                                    CenteredB[2][i]*CenteredB[2][i]);

OldBaryonField[Eeng][i] -= 0.5*(OldCenteredB[0][i]*OldCenteredB[0][i]+
                                    OldCenteredB[1][i]*OldCenteredB[1][i]+
                                    OldCenteredB[2][i]*OldCenteredB[2][i]);
OldCenteredB[0][i] *= temp_sqrt_a;
OldCenteredB[1][i] *= temp_sqrt_a;
OldCenteredB[2][i] *= temp_sqrt_a;                                                                                                                                                             
OldBaryonField[Eeng][i] += 0.5*(OldCenteredB[0][i]*OldCenteredB[0][i]+
                                    OldCenteredB[1][i]*OldCenteredB[1][i]+
                                    OldCenteredB[2][i]*OldCenteredB[2][i]);
}
for(field=0;field<3;field++)
for( i=0; i<MagneticSize[field]; i++){
MagneticField[field][i] *= temp_sqrt_a;
OldMagneticField[field][i] *= temp_sqrt_a;
}

 // ElectricField=ElectricField*sqrt(a)
if(ComovingCoordinates==1){
   for(field=0;field<3;field++)
     for(i=0;i<ElectricSize[field];i++)
   ElectricField[field][i] *= temp_sqrt_a;
}
                                                                                                                                                             
}
                                                                                                                                                             
#endif
                                                                                                                                                            
  //
  //
  // MagneticField *= sqrt(a)  
  // OldMagneticField *= sqrt(a)  
  //

  //<dbg>
  //fprintf(stderr," NumberOfBaryonFields; %d\n",NumberOfBaryonFields);
  //</dbg>

  //
  // Delete solver fluxes, pressure field.
  //

  //<dbg>
  for( dim=0;dim<3; dim++)
    if( SolverProblems != NULL ){
      delete [] SolverProblems[dim];
      SolverProblems[dim] = NULL;
    }
  //</dbg>
  
  if( MHD_Pressure != NULL ){
    delete [] MHD_Pressure;
    MHD_Pressure = NULL;
  }
  for( dim=0; dim<3; dim++)
    delete [] Fluxes[dim];

  if( FlatteningField != NULL )
    delete [] FlatteningField;

  //JBMEM_MESSAGE(MyProcessorNumber,"jb: Leave Ath");
  return SUCCESS;
}


#endif //ATHENA







