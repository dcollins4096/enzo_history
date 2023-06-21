//
// Computes fluxes, takes differences.  PPML.
// Nov. 18 2007: added Athena CT. Tested against DaveThena, identical results.
//
#include <stdio.h>
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "PPML.h"
#include "Athena.h"
#include "DaveTools.h"


extern "C" void FORTRAN_NAME(ppml_flux)( int * nx, float * dx, float * dt,     
					 float * ql, float * qr, float * qm, 
					 float * udy, float * udz, float * f,
					 float * source,
					 float * a1, float * a2, 
					 float * courant,float * omloc, int * kstep, 
					 int* Riemann, int * Reconstruction,
					 int* DirectionalSplitting,
					 int* EigenSystem);
/*
extern "C" void FORTRAN_NAME(potoky)( int * nx, float * dx, float * dt,     
				      float * ql, float * qr, float * qm, 
				      float * udx, float * udz, float * f,
				      float * source,
				      float * a1, float * a2, 
				      float * courant,float * omloc, int * kstep, 
				      int* Riemann, int * Reconstruction,
				      int* DirectionalSplitting);

extern "C" void FORTRAN_NAME(potokz)( int * nx, float * dx, float * dt,     
				      float * ql, float * qr, float * qm, 
				      float * udx, float * udy, float * f,
				      float * source,
				      float * a1, float * a2, 
				      float * courant,float * omloc, int * kstep,
				      int* Riemann, int * Reconstruction,
				      int* DirectionalSplitting);
*/
extern "C" void FORTRAN_NAME(ppml_ct)(int * nx,int * ny,int * nz, int* Ngz,
				      float * dt, float * dx, float * dy, float * dz, 
				      float * F, float * G,float * H, 
				      float * Bx,  float * By, float * Bz, int * rank);


//PPML and Enzo use different orderings for the variables.
//code at the bottom of this routine.
void GeneratePPMLtoEnzoMap(int ** PPMLtoEnzoMap, int NumberOfFluidQuantities, IndexPointerMap ind); 


//Source Term:  Will need a better home in the future.
//Dim is currently unused-- might be used later.
float grid::ReturnSourceTerm(int field, IndexPointerMap ind, int index, int dim){
  float Output = 0;
  if( PointSourceGravity || UniformGravity || SelfGravity ){
    if( field == ind.VX ) Output += AccelerationField[0][ index ]*dtFixed*0.5;
    if( GridRank > 1 && field == ind.VY ) Output += AccelerationField[1][ index ]*dtFixed*0.5;
    if( GridRank > 2 && field == ind.VZ ) Output += AccelerationField[2][ index ]*dtFixed*0.5;
  }
  return Output;
}

int grid::PPML_Update(int CycleNumber, int level, int grid){


  if( ProcessorNumber != MyProcessorNumber )
    return SUCCESS;
#ifdef MHDF
  if( MHD_DivB != BalsaraToth ){
    MHD_Allocate();
    PPML_MHD_Clone(1);
  }
#endif //MHDF

  if( MidWayDumpCheck( 30 ) == TRUE ){
    char basename[30];
    sprintf(basename, "data30%d%d.grid",CycleNumber, level);
    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, grid+1) == FAIL ){
      fprintf(stderr, "Error in Write Grid, PPML_Update.\n");
      return FAIL;
    }
    fclose(dummy);
  }

  //fprintf(stderr," PPM Update\n");

  int i,j,k,field, index, *PPMLtoEnzoMap[3];
  int EnzoIndex, PpmlIndex, FluxIndex;
  int B_ijk, F_ijk, F_im1, F_jm1, F_km1;
  int size = 1;
  int dim;
  float * Qpr, * Qpl,* Qpm;        //Right, Left and Center stripes.
  float * Qdx, *Qdy, *Qdz;         //transverse differences
  float * Fm, *Gm, *Hm;            // Temp 1d fluxes
  float * FluxX, * FluxY, * FluxZ; //Full 3d fluxes
  float * Source;                  //Source terms.  Currently gravity.
  int RiemannLocal;  //for 2.5d and testing.
  int SplittingLoop, DoX, DoY, DoZ; //For directional splitting.  

  float a1 = IsothermalSoundSpeed;
  float a2 = IsothermalSoundSpeed*IsothermalSoundSpeed;

  //Parameters for use in the MUSTA riemann solver.
  //Currently not called, but the hook is there for future experimentation.
  float omloc = 1.0/(1.0 + CourantSafetyNumber );
  int kstep = 2;
  //
  // Layer between Physical quantities and Data structures.
  // Not in the slightest bit pretty.
  //

  PPML_InterfacePointerBundle Face( this );
  if( Face.Error == TRUE )
    ERROR_MESSAGE;
#ifdef trash
  if( this->ReturnInterfacePointers( Face ) == FAIL ){
    fprintf(stderr," ReturnInterfacePointers failed.\n"); return FAIL;}
#endif //trash
  IndexPointerMap ind;
  if( this->IdentifyPhysicalQuantities_2(ind) == FAIL ){
    fprintf(stderr," IdentifyPhysicalQuantities_2 failed\n");return FAIL;}
  
  //
  // Artificial Viscosity from Athena.
  // The Athena method attaches itself to the Grid, and modifies the BaryonField data.
  //
  
  Athena ATH(this);
  if( ATH.Error == TRUE ){
    fprintf(stderr,"Grid_PPML_Update: Error with Athena object setup.\n");
    return FALSE;
  }
  

  //
  // To ease our coding, we also make a map between PPML indexing and Enzo indexing.
  // It takes the PPML index, and returns the Enzo index.

  /*
    fprintf(stderr,"monkey d %d  e %d x %d y %d z %d x %d y %d z %d\n", 
    ind.D, ind.TE, ind.VX, ind.VY, ind.VZ, ind.BX, ind.BY, ind.BZ);
  */

  GeneratePPMLtoEnzoMap( PPMLtoEnzoMap ,NumberOfFluidQuantities, ind);
  for( i=0; i<GridRank; i++)
    size *= GridDimension[i];


  float ** TranX = new float * [NumberOfFluidQuantities];
  float ** TranY = new float * [NumberOfFluidQuantities];
  float ** TranZ = new float * [NumberOfFluidQuantities];
  for( field=0; field< NumberOfFluidQuantities; field++){
    TranX[field] = new float[size];
    TranY[field] = new float[size];
    TranZ[field] = new float[size];
    for( i=0;i<size; i++){
      TranX[field][i] = Face.X_R[ field ][i] - Face.X_L[field][i];
      TranY[field][i] = Face.Y_R[ field ][i] - Face.Y_L[field][i];
      TranZ[field][i] = Face.Z_R[ field ][i] - Face.Z_L[field][i];
    }
  }

  //
  // Directional Sweep loop-- Still somewhat experimental.
  //
  int Ndims = 3; //Always 3.
  int StrangCycleStart = ( (MHD_DirectionalSplitting == 1 ) ? CycleNumber % Ndims : 0 );
  int StrangCycleEnd = ( (MHD_DirectionalSplitting == 1 ) ? (StrangCycleStart+Ndims) : 1 );
  for( int StrangLoop = StrangCycleStart; StrangLoop < StrangCycleEnd; StrangLoop++ ){   
    SplittingLoop = StrangLoop % Ndims;
    if( MHD_DirectionalSplitting == 1 ){
      DoX = ( (SplittingLoop  == 0 ) ? 1 : 0 );
      DoY = ( (SplittingLoop  == 1 ) ? 1 : 0 );
      DoZ = ( (SplittingLoop  == 2 ) ? 1 : 0 );
    }else{
      DoX = DoY = DoZ = 1;
    }

    if( MyProcessorNumber == ROOT_PROCESSOR) 
      fprintf(stderr,"Splitting Loop %d X %d Y %d Z %d\n", SplittingLoop, DoX, DoY, DoZ);

    if( DoX == 1 ){

      // 
      // X flux
      //
      dim = 0;
      RiemannLocal = MHD_Riemann[0];       
      Qpr = new float[ GridDimension[0] * NumberOfFluidQuantities];
      Qpl = new float[ GridDimension[0] * NumberOfFluidQuantities];
      Qpm = new float[ GridDimension[0] * NumberOfFluidQuantities];
      Qdy = new float[ GridDimension[0] * NumberOfFluidQuantities];
      Qdz = new float[ GridDimension[0] * NumberOfFluidQuantities];
      Fm  = new float[ GridDimension[0] * NumberOfFluidQuantities];
      Source = new float[ GridDimension[0] * NumberOfFluidQuantities];
      FluxX = new float[ size * NumberOfFluidQuantities];
      
      if( dccCounter09++ == -1 ){
	a1 = -1;
	RiemannLocal = MHD_Riemann[0];
	//RiemannLocal = -13;// MHD_Riemann[0];
	
	fprintf(stderr,"Riemann Hack, Bitches!\n");
	fprintf(stderr,"Riemann Hack, Bitches!\n");
	fprintf(stderr,"Riemann Hack, Bitches!\n");
      }else{
	RiemannLocal = MHD_Riemann[0];
      }
      for(k=0; k<GridDimension[2]; k++){
	for(j=0; j<GridDimension[1]; j++){
	  
	  for(field=0;field<NumberOfFluidQuantities; field++)
	    for(i=0; i<GridDimension[0]; i++){
	      EnzoIndex = i + GridDimension[0]*(j + GridDimension[1]*k);
	      PpmlIndex = i + GridDimension[0]*field;
	      Qpr[ PpmlIndex ] = Face.X_R[ PPMLtoEnzoMap[dim][ field ] ][ EnzoIndex ];
	      Qpl[ PpmlIndex ] = Face.X_L[ PPMLtoEnzoMap[dim][ field ] ][ EnzoIndex ];
	      Qpm[ PpmlIndex ] = BaryonField[ PPMLtoEnzoMap[dim][ field ] ][ EnzoIndex ];
	      Qdy[ PpmlIndex ] = TranY[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ];
	      Qdz[ PpmlIndex ] = TranZ[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ];
	      Source[ PpmlIndex] = ReturnSourceTerm( PPMLtoEnzoMap[dim][field], ind, EnzoIndex, 0);
	      //<dbg>
	      //Fm[ PpmlIndex ] =  Face.X_R[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ];
	      //</dbg>
	      
	    }//field, i
	  
	  FORTRAN_NAME(ppml_flux)( GridDimension , CellWidth[0] , &dtFixed, 
				Qpl, Qpr, Qpm, Qdy, Qdz, Fm, Source, &a1, &a2,
				&CourantSafetyNumber,&omloc, &kstep, &RiemannLocal, &(MHD_Recon[0]),
				&MHD_DirectionalSplitting, &PPML_EigenSystem);
	  for(field=0;field<NumberOfFluidQuantities; field++)
	    for(i=0; i<GridDimension[0]; i++){
	      FluxIndex = i + GridDimension[0]*(j + GridDimension[1]*
						(k + GridDimension[2]*PPMLtoEnzoMap[dim][field]));
	      EnzoIndex = i + GridDimension[0]*(j + GridDimension[1]*k);
	      PpmlIndex = i + GridDimension[0]*field;
	      Face.X_L[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ] = Qpl[ PpmlIndex ];
	      Face.X_R[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ] = Qpr[ PpmlIndex ];
	      FluxX[ FluxIndex ] =  Fm[ PpmlIndex ];
	      //FluxX[ FluxIndex ] = CourantSafetyNumber;
	    }//field,i
	}//j
      }//k
      
      delete Qpr ;
      delete Qpl ;
      delete Qpm; 
      delete Qdy ;
      delete Qdz ;
      delete Fm  ;
      delete Source;
    }//doX

    if( DoY == 1 ){

      // 
      // Y flux
      //
      RiemannLocal = ( (GridRank > 1 ) ? MHD_Riemann[0] : 3 );       
      dim = 1;
      Qpr = new float[ GridDimension[1] * NumberOfFluidQuantities];
      Qpl = new float[ GridDimension[1] * NumberOfFluidQuantities];
      Qpm = new float[ GridDimension[1] * NumberOfFluidQuantities];
      Qdx = new float[ GridDimension[1] * NumberOfFluidQuantities];
      Qdz = new float[ GridDimension[1] * NumberOfFluidQuantities];
      Gm  = new float[ GridDimension[1] * NumberOfFluidQuantities];
      FluxY = new float[ size * NumberOfFluidQuantities];
      Source = new float[ GridDimension[1] * NumberOfFluidQuantities];
      
      for(k=0; k<GridDimension[2]; k++){
	for(i=0; i<GridDimension[0]; i++){
	  
	  for(field=0;field<NumberOfFluidQuantities; field++)
	    for(j=0; j<GridDimension[1]; j++){
	      
	      EnzoIndex = i + GridDimension[0]*(j + GridDimension[1]*k);
	      PpmlIndex = j + GridDimension[1]*field;
	      Qpr[ PpmlIndex ] = Face.Y_R[ PPMLtoEnzoMap[dim][ field ] ][ EnzoIndex ];
	      Qpl[ PpmlIndex ] = Face.Y_L[ PPMLtoEnzoMap[dim][ field ] ][ EnzoIndex ];
	      Qpm[ PpmlIndex ] = BaryonField[ PPMLtoEnzoMap[dim][ field ] ][ EnzoIndex ];
	      Qdx[ PpmlIndex ] = TranX[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ];
	      Qdz[ PpmlIndex ] = TranZ[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ];
	      Source[ PpmlIndex] = ReturnSourceTerm( PPMLtoEnzoMap[dim][field], ind, EnzoIndex, 1);
	      
	    }//field, i
	  
	  
	  FORTRAN_NAME(ppml_flux)( GridDimension+1 , CellWidth[1] , &dtFixed, 
				   Qpl, Qpr, Qpm, Qdx, Qdz, Gm, Source, &a1, &a2,
				   &CourantSafetyNumber,&omloc, &kstep, &RiemannLocal, &(MHD_Recon[0]),
				   &MHD_DirectionalSplitting, &PPML_EigenSystem);

	  for(field=0;field<NumberOfFluidQuantities; field++)
	    for(j=0; j<GridDimension[1]; j++){
	      FluxIndex = i + GridDimension[0]*(j + GridDimension[1]*
						(k + GridDimension[2]*PPMLtoEnzoMap[dim][field]));
	      EnzoIndex = i + GridDimension[0]*(j + GridDimension[1]*k);
	      PpmlIndex = j + GridDimension[1]*field;
	      Face.Y_L[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ] = Qpl[ PpmlIndex ];
	      Face.Y_R[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ] = Qpr[ PpmlIndex ];
	      FluxY[ FluxIndex ] =  Gm[ PpmlIndex ];
	    }//field,i
	}//j
      }//k
      
      delete Qpr ;
      delete Qpl ;
      delete Qpm; 
      delete Qdx ;
      delete Qdz ;
      delete Gm  ;
      delete Source;

    }//DoY
    if( DoZ == 1 ){
      // 
      // Z flux
      //
      

      //locally sensitive switch for riemann solver: for 2.5 d problems, switch to the first order flux.
      //Numbering consistent with the DaveThena numbering scheme.
      dim = 2;
      RiemannLocal = ( (GridRank > 2 ) ? MHD_Riemann[0] : 3 ); 
      
      Qpr = new float[ GridDimension[2] * NumberOfFluidQuantities];
      Qpl = new float[ GridDimension[2] * NumberOfFluidQuantities];
      Qpm = new float[ GridDimension[2] * NumberOfFluidQuantities];
      Qdx = new float[ GridDimension[2] * NumberOfFluidQuantities];
      Qdy = new float[ GridDimension[2] * NumberOfFluidQuantities];
      Hm  = new float[ GridDimension[2] * NumberOfFluidQuantities];
      FluxZ = new float[ size * NumberOfFluidQuantities];
      Source = new float[ GridDimension[2] * NumberOfFluidQuantities];
      for(j=0; j<GridDimension[1]; j++){
	for(i=0; i<GridDimension[0]; i++){
	  
	  for(field=0;field<NumberOfFluidQuantities; field++)
	    for(k=0; k<GridDimension[2]; k++){
	      EnzoIndex = i + GridDimension[0]*(j + GridDimension[1]*k);
	      PpmlIndex = k + GridDimension[2]*field;
	      Qpr[ PpmlIndex ] = Face.Z_R[ PPMLtoEnzoMap[dim][ field ] ][ EnzoIndex ];
	      Qpl[ PpmlIndex ] = Face.Z_L[ PPMLtoEnzoMap[dim][ field ] ][ EnzoIndex ];
	      Qpm[ PpmlIndex ] = BaryonField[ PPMLtoEnzoMap[dim][ field ] ][ EnzoIndex ];
	      Qdx[ PpmlIndex ] = TranX[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ];
	      Qdy[ PpmlIndex ] = TranY[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ];
	      Source[ PpmlIndex] = ReturnSourceTerm( PPMLtoEnzoMap[dim][field], ind, EnzoIndex, 0);
	      
	    }//field, i
	  
	  
	  FORTRAN_NAME(ppml_flux)( GridDimension+2 , CellWidth[2] , &dtFixed, 
				   Qpl, Qpr, Qpm, Qdx, Qdy, Hm, Source,&a1, &a2,
				   &CourantSafetyNumber,&omloc, &kstep, &RiemannLocal, &(MHD_Recon[0]),
				   &MHD_DirectionalSplitting, &PPML_EigenSystem);
	  
	  for(field=0;field<NumberOfFluidQuantities; field++)
	    for(k=0; k<GridDimension[2]; k++){
	      FluxIndex = i + GridDimension[0]*(j + GridDimension[1]*
						(k + GridDimension[2]*PPMLtoEnzoMap[dim][field]));
	      EnzoIndex = i + GridDimension[0]*(j + GridDimension[1]*k);
	      PpmlIndex = k + GridDimension[2]*field;
	      Face.Z_L[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ] = Qpl[ PpmlIndex ];
	      Face.Z_R[ PPMLtoEnzoMap[dim][field] ][ EnzoIndex ] = Qpr[ PpmlIndex ];
	      FluxZ[ FluxIndex ] =  Hm[ PpmlIndex ];
	    }//field,i
	}//j
      }//k
      
      delete Qpr ;
      delete Qpl ;
      delete Qpm; 
      delete Qdx ;
      delete Qdy ;
      delete Hm  ;
      delete Source;
    }//doZ
    
    if( MHD_DiffusionMethod[0] > 0 ){
      if( ATH.PPMLViscosity( FluxX, FluxY, FluxZ, SplittingLoop ) == FAIL )
	{fprintf(stderr," PPML Update: PPMLViscosity failure\n"); return FAIL;}
    }
    
    //convert velocity to momentum.
    //Cumbersome.  Should be temporary.

    for( i=0; i<size; i++ ){
      BaryonField[ ind.VX ][ i ] *= BaryonField[ ind.D ][ i ];
      BaryonField[ ind.VY ][ i ] *= BaryonField[ ind.D ][ i ];
      BaryonField[ ind.VZ ][ i ] *= BaryonField[ ind.D ][ i ];
    }
    
    float dtdx = dtFixed/CellWidth[0][0];
    float dtdy = dtFixed/CellWidth[0][1];
    float dtdz = dtFixed/CellWidth[0][2];
    
    //Density and velocities at the present time step. 
    //Used for gravity step.
    float d_n, v_n[3];
    
    //Booleans to determine which differences to take.
    int diffX, diffY, diffZ;
    //Extents of the flux difference loop.
    //   For non BalsaraToth (MHD_DivB == 4) CT methods, an extra layer of differences
    //   is needed to update the data used for CT.  
    //   BalsaraSpicer might not need this extra layer-- this hasn't been examined/tested.
    //
    int ExtraX, ExtraY, ExtraZ;
    ExtraX = ( ( MHD_DivB == BalsaraToth ) ? 0 : 1 );
    ExtraY = ( ( MHD_DivB == BalsaraToth || GridRank <= 1 ) ? 0 : 1 );
    ExtraZ = ( ( MHD_DivB == BalsaraToth || GridRank <= 2 ) ? 0 : 1 );
    for(k=GridStartIndex[2] -ExtraZ ; k<=GridEndIndex[2] + ExtraZ; k++)
      for(j=GridStartIndex[1] -ExtraY ; j<=GridEndIndex[1] +ExtraY; j++)
	for(i=GridStartIndex[0] -ExtraX; i<=GridEndIndex[0] +ExtraX; i++){
	  
	  
	  //BaryonField index.
	  B_ijk = (i  ) + GridDimension[0]*((j  ) + GridDimension[1]*(k  ) );
	  
	  //Save ThisStep for gravity solve.  See Colella & Woodward JCP 1984 eqn (3.8)
	  if( SelfGravity || UniformGravity || PointSourceGravity ){
	    d_n = BaryonField[ ind.D  ][ B_ijk ];
	    v_n[0]= BaryonField[ ind.VX ][ B_ijk ];
	    v_n[1]= BaryonField[ ind.VY ][ B_ijk ];
	    v_n[2]= BaryonField[ ind.VZ ][ B_ijk ];
	  }
	  for( field=0;field< NumberOfFluidQuantities; field++){
	    //Magnetic field updated w/ ct, not flux differences:
#ifndef MHDF	     
	    if( MHD_Used == TRUE ){ 
	      if( field == ind.BX || field == ind.BY || field == ind.BZ ){
		continue;}}
#else
	    if( MHD_Used == TRUE || MHD_DivB == BalsaraToth){ 
	      if( field == ind.BX || field == ind.BY || field == ind.BZ ){
		continue;}}
	    
#endif //! MHDF
	    //Flux arrays are indexed such that
	    //Fx[i,j,k] maps to {i+1/2,j,k}
	    //Fy[i,j,k] maps to {i,j+1/2,k}
	    //Fz[i,j,k] maps to {i,j,k+1/2}
	    F_ijk =  i    + GridDimension[0]*( j    + GridDimension[1]*(k +   GridDimension[2]*field));
	    F_im1 = (i-1) + GridDimension[0]*((j  ) + GridDimension[1]*(k   + GridDimension[2]*field) );
	    F_jm1 = (i  ) + GridDimension[0]*((j-1) + GridDimension[1]*(k   + GridDimension[2]*field) );
	    F_km1 = (i  ) + GridDimension[0]*((j  ) + GridDimension[1]*(k-1 + GridDimension[2]*field) );
	    
#ifdef MHDF
	    diffX = 1;
	    diffY = ( (GridRank > 1 ) ? 1 : 0 );
	    diffZ = ( (GridRank > 2 ) ? 1 : 0 );
	    //Magnetic quantities only want the cross differences.  
	    //There might be junk in the longitudinal fluxes that we want to ignore.
	    diffX *= ( (field == ind.BX ) ? 0 : 1 );
	    diffY *= ( (field == ind.BY ) ? 0 : 1 );
	    diffZ *= ( (field == ind.BZ ) ? 0 : 1 );
	    //For dimensional splitting
	    diffX *= DoX;
	    diffY *= DoY;
	    diffZ *= DoZ;

	    BaryonField[ field ][B_ijk] =  (BaryonField[ field][B_ijk] 
					    - (( 1 == diffX ) ? dtdx*( FluxX[F_ijk] - FluxX[F_im1]) : 0 )
					    - (( 1 == diffY ) ? dtdy*( FluxY[F_ijk] - FluxY[F_jm1]) : 0 )
					    - (( 1 == diffZ ) ? dtdz*( FluxZ[F_ijk] - FluxZ[F_km1]) : 0 ) );
	    
#else
	    BaryonField[ field ][B_ijk] =  (BaryonField[ field][B_ijk] 
					    -dtdx*( FluxX[F_ijk] - FluxX[F_im1])
					    - (( GridRank > 1)? dtdy*( FluxY[F_ijk] - FluxY[F_jm1]) : 0 )
					    - (( GridRank > 2)? dtdz*( FluxZ[F_ijk] - FluxZ[F_km1]) : 0 ) );
	    
	    
#endif //MHDF	  
	  }//field
	  
	  //Gravitational difference.
	  if( SelfGravity || UniformGravity || PointSourceGravity ){
	    for( dim=0;dim<GridRank;dim++){
	      if( dim== 0 && DoX == 0 ) continue;
	      if( dim== 1 && DoY == 0 ) continue;
	      if( dim== 2 && DoZ == 0 ) continue;

	      BaryonField[ ind.V[dim] ][B_ijk] += dtFixed * AccelerationField[ dim ][ B_ijk ]*
		0.5*(BaryonField[ ind.D ][ B_ijk ] + d_n );
	      
	      //Currently, only Isothermal EOS works with PPML; this is here for ease in the future.
	      if( EquationOfState == 0 ){
		BaryonField[ ind.TE ][ B_ijk ] +=dtFixed * AccelerationField[ dim ][ B_ijk ]*
		  0.5*( BaryonField[ ind.D ][ B_ijk ]*BaryonField[ ind.V[dim] ][ B_ijk ] +
			d_n*v_n[dim] );
	      }//EOS
	    }//dim loop
	  }//gravity
	}//i,j,k
    
    for( i=0; i<size; i++ ){
      BaryonField[ ind.VX ][ i ] /= BaryonField[ ind.D ][ i ];
      BaryonField[ ind.VY ][ i ] /= BaryonField[ ind.D ][ i ];
      BaryonField[ ind.VZ ][ i ] /= BaryonField[ ind.D ][ i ];
    }
  }//DirectionalSplitting loop
  // Other cleanup
  for( dim=0; dim<3; dim++)
    delete PPMLtoEnzoMap[dim];

  int ngz = DEFAULT_GHOST_ZONES;
  
  //JBMEM_MESSAGE(MyProcessorNumber, "MHD_Allocate");  
  //JBMEM_MESSAGE(MyProcessorNumber, "MHD_Middle");  
  
#ifdef MHDF
  //Some might say it was a mistake to not have designed it this way the first time.
  float* Fluxes[3] = {FluxX, FluxY, FluxZ};
  switch( MHD_DivB){
  case BalsaraSpicer: //1
  case Athena_LF: //2
  case Athena_Switch: //3
    
    ATH.ComputeElectricField(dtFixed, Fluxes);
    MHD_Curl( GridStartIndex, GridEndIndex, 1);
    CenterMagneticField();
    //Copy the newly created Face Centered Field to the Pseudo Face Centered
    //field stored for PPML.
    PPML_MHD_Clone(0);
    break;
  case BalsaraToth: //4
    FORTRAN_NAME(ppml_ct)(GridDimension,GridDimension+1,GridDimension+2, &ngz,
			  &dtFixed, CellWidth[0], CellWidth[1], CellWidth[2],
			  FluxX, FluxY, FluxZ, 
			  BaryonField[ ind.BX ], BaryonField[ ind.BY ], BaryonField[ ind.BZ ], 
			  &GridRank);
    break;
  default:
    fprintf(stderr,"Unrecognized DivB selected (%"ISYM"). I hope you know what your doing.\n",
	    MHD_DivB);
    
  }//divB switch

#else
  
  FORTRAN_NAME(ppml_ct)(GridDimension,GridDimension+1,GridDimension+2, &ngz,
			&dtFixed, CellWidth[0], CellWidth[1], CellWidth[2],
			FluxX, FluxY, FluxZ, 
			BaryonField[ ind.BX ], BaryonField[ ind.BY ], BaryonField[ ind.BZ ], 
			&GridRank);
#endif //MHDF
  //JBMEM_MESSAGE(MyProcessorNumber, "MHD_Deallocate");
  //
  // Press Box.  Transmit the Data to the people.
  // Some bad pointer juggling: save Fluxes in Centered BaryonFields.
  //

  //X flux.
  if( MidWayDumpCheck( 31 ) == TRUE ){
    char basename[30];
    float **ActualBaryonField = new float *[NumberOfFluidQuantities];
    for( field=0;field<NumberOfFluidQuantities; field++){
      ActualBaryonField[ field ] = BaryonField[ field ];
      BaryonField[ field ] = FluxX + size*field;
    }
    sprintf(basename, "data31%d%d.grid",CycleNumber, level);
    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, grid+1) == FAIL ){
      fprintf(stderr, "Error in Write Grid, PPML_Update.\n");
      return FAIL;
    }
    for( field=0;field<NumberOfFluidQuantities; field++){
      BaryonField[ field ]=ActualBaryonField[ field ];
    }

    fclose(dummy);
  }


  //y flux.
  if( MidWayDumpCheck( 32 ) == TRUE ){
    char basename[30];
    float **ActualBaryonField = new float *[NumberOfFluidQuantities];
    for( field=0;field<NumberOfFluidQuantities; field++){
      ActualBaryonField[ field ] = BaryonField[ field ];
      BaryonField[ field ] = FluxY + size*field;
    }
    sprintf(basename, "data32%d%d.grid",CycleNumber, level);
    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, grid+1) == FAIL ){
      fprintf(stderr, "Error in Write Grid, PPML_Update.\n");
      return FAIL;
    }
    for( field=0;field<NumberOfFluidQuantities; field++){
      BaryonField[ field ]=ActualBaryonField[ field ];
    }

    fclose(dummy);
  }


  //Z flux.
  if( MidWayDumpCheck( 33 ) == TRUE ){
    char basename[30];
    float **ActualBaryonField = new float *[NumberOfFluidQuantities];
    for( field=0;field<NumberOfFluidQuantities; field++){
      ActualBaryonField[ field ] = BaryonField[ field ];
      BaryonField[ field ] = FluxZ + size*field;
    }
    sprintf(basename, "data33%d%d.grid",CycleNumber, level);
    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, grid+1) == FAIL ){
      fprintf(stderr, "Error in Write Grid, PPML_Update.\n");
      return FAIL;
    }
    for( field=0;field<NumberOfFluidQuantities; field++){
      BaryonField[ field ]=ActualBaryonField[ field ];
    }

    fclose(dummy);
  }


  if( MidWayDumpCheck( 39 ) == TRUE ){
    char basename[30];
    sprintf(basename, "data39%d%d.grid",CycleNumber, level);
    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, grid+1) == FAIL ){
      fprintf(stderr, "Error in Write Grid, PPML_Update.\n");
      return FAIL;
    }
    fclose(dummy);
  }
  delete FluxX;
  delete FluxY;
  delete FluxZ;

  for( field=0; field< NumberOfFluidQuantities; field++){
    delete [] TranX[field];
    delete [] TranY[field];
    delete [] TranZ[field];
  }
  delete TranX;
  delete TranY;
  delete TranZ;

#ifdef trash
  if( this->DeleteInterfacePointers( Face ) == FAIL ){
    fprintf(stderr,"Grid_PPML_Mono1: error with Delete InterfacePointers");
    return FAIL;
  }
#endif //trash
#ifdef MHDF
  if( MHD_DivB != BalsaraToth ){
    MHD_Deallocate();
  }
#endif //MHDF
  return SUCCESS;
}

void GeneratePPMLtoEnzoMap(int ** PPMLtoEnzoMap, int NumberOfFluidQuantities, IndexPointerMap ind )
{
  int i = 0, dim=0;  
  for( dim=0; dim<3; dim++)
    PPMLtoEnzoMap[dim] = new int[ NumberOfFluidQuantities] ;
  
  // X
  dim = 0;
  i = 0;
  //The treatment of the total energy variable is, for now, a place holder.
  PPMLtoEnzoMap[dim]  [ i++ ] = ind.D;
  if( EquationOfState == 0 )
    {PPMLtoEnzoMap[dim]  [ i++ ] = ind.TE;
    fprintf(stderr," Mono1: PPML doesn't yet support TotalEnergy, so you'd better not have it in the map.\n");}
  PPMLtoEnzoMap[dim][ i++ ] = ind.VX;
  PPMLtoEnzoMap[dim][ i++ ] = ind.VY;
  PPMLtoEnzoMap[dim][ i++ ] = ind.VZ;
  if( MHD_Used == TRUE ){
    PPMLtoEnzoMap[dim][ i++ ] = ind.BX;
    PPMLtoEnzoMap[dim][ i++ ] = ind.BY;
    PPMLtoEnzoMap[dim][ i++ ] = ind.BZ;
  }

  // Y
  dim = 1;
  i = 0;
  //The treatment of the total energy variable is, for now, a place holder.
  PPMLtoEnzoMap[dim]  [ i++ ] = ind.D;
  if( EquationOfState == 0 )
    {PPMLtoEnzoMap[dim]  [ i++ ] = ind.TE;
    fprintf(stderr," Mono1: PPML doesn't yet support TotalEnergy, so you'd better not have it in the map.\n");}
  PPMLtoEnzoMap[dim][ i++ ] = ind.VY;
  PPMLtoEnzoMap[dim][ i++ ] = ind.VZ;
  PPMLtoEnzoMap[dim][ i++ ] = ind.VX;
  if( MHD_Used == TRUE ){
    PPMLtoEnzoMap[dim][ i++ ] = ind.BY;
    PPMLtoEnzoMap[dim][ i++ ] = ind.BZ;
    PPMLtoEnzoMap[dim][ i++ ] = ind.BX;
  }

  // Z
  dim = 2;
  i = 0;
  //The treatment of the total energy variable is, for now, a place holder.
  PPMLtoEnzoMap[dim]  [ i++ ] = ind.D;
  if( EquationOfState == 0 )
    {PPMLtoEnzoMap[dim]  [ i++ ] = ind.TE;
    fprintf(stderr," Mono1: PPML doesn't yet support TotalEnergy, so you'd better not have it in the map.\n");}
  PPMLtoEnzoMap[dim][ i++ ] = ind.VZ;
  PPMLtoEnzoMap[dim][ i++ ] = ind.VX;
  PPMLtoEnzoMap[dim][ i++ ] = ind.VY;
  if( MHD_Used == TRUE ){
    PPMLtoEnzoMap[dim][ i++ ] = ind.BZ;
    PPMLtoEnzoMap[dim][ i++ ] = ind.BX;
    PPMLtoEnzoMap[dim][ i++ ] = ind.BY;
  }

}
