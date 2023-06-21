
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Athena.h"



#ifdef HAOXU
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
#endif 

//
// Nov. 7 2007: Objectified the entire thing.
// Only built up as far as I needed to get the Viscosity running.
// Will need to be refactored when the rest of the Athena solver is rolled in.
// Some things to consider upon said refactor:
// 1.) Needs a big 'ol sed to replase E.D, S.D w/ E.D, S.D (etc.)
// 2.) Probably should have a pint of whiskey ready.

//
// Sets up parameters for athena:
// Offset[3]                : distance in memory between neighboring points in the {x,y,z} direction.
// dxI[3]                   : 1/d{x,y,z}
// GhostZones[2]            : # ghost zones for the first and second pass.
// AthStart[3], AthEnd[3]   : start and end indicies.  Changed for 1st and second pass.
// NumberOfMHDFluxes        : 7 for adiabatic, 6 for isothermal.
// E.D, Egas, E.V[3], E.TE  : indicies for BaryonField.
// S.D, S.V[3], S.TE        : indicies for Solver.
// MHD_PLM_SlopeLocal       : maps the input parameter to the useful parameter.
// MapEtoS                  : BaryonField[ field ] = State[  MapEtoS[field] ]
// BNum[dim][transvers0,1]  : transvers B component indecies for each face.
// AllCellCentered          : flag to denote that only the cell centered field is used.
// NumberReconstructed      : Number of fluxes, +1 if  B_{face} isn't used.
//

//Dumb Constructor
Athena::Athena(){
}
Athena::~Athena(){
  fprintf(stderr,"Delete Athena.\n");
}
Athena::Athena(grid * Grid){

  int dim, field;
  Error = TRUE; //Assume the worst.

  //
  // Clone parameters from Grid.
  //

  GridRank = Grid->GridRank;
  for( dim=0; dim< MAX_DIMENSION; dim++){
    GridDimension[dim] = Grid->GridDimension[dim];
    GridStartIndex[dim] = Grid->GridStartIndex[dim];
    GridEndIndex[dim] = Grid->GridEndIndex[dim];
    CellWidth[dim] = Grid->CellWidth[dim];
#ifdef MHDF
  NFaces = Grid->NFaces;
  NEdges = Grid->NEdges;

    for(field=0; field<NFaces; field++){
      MagneticDims[field][dim] = Grid->MagneticDims[field][dim];
      ElectricDims[field][dim] = Grid->ElectricDims[field][dim];
    }
#endif //MHDF

    AthStart[dim] = Grid->GridStartIndex[dim];
    AthEnd[dim] = Grid->GridEndIndex[dim];
  }

  //Enzo index map E

  if( Grid->IdentifyPhysicalQuantities_2(E) == FAIL ){
    fprintf(stderr," IdentifyPhysicalQuantities_2 failed\n");Error = TRUE;}

  dtFixed = Grid->dtFixed;
  NumberOfBaryonFields = Grid->NumberOfBaryonFields;
  NumberOfFluidQuantities = Grid->NumberOfFluidQuantities;

  // Coppies for maximum confusion.  
  // Unclear if this is the right way to go.  Might not be.  
  for( field = 0; field< NumberOfBaryonFields; field++){
    BaryonField[field] = Grid->BaryonField[field];
    OldBaryonField[field] = Grid->OldBaryonField[field];
  }

#ifdef MHDF
  for( field=0; field<NFaces; field++){
    CenteredB[field] = Grid->BaryonField[ E.B[field] ];
    MagneticField[field] = Grid->MagneticField[field];
    ElectricField[field] = Grid->ElectricField[field];
  }
#endif //MHDF

  for( field = 0; field<RandomForcingNumberOfFields; field++){
    RandomForcingField[field] = RandomForcingField[field];
  }

  for( field = 0; field<MAX_NUMBER_OF_BARYON_FIELDS + 3; field++)
    State[field] = NULL;

  //
  // Parameters derived from cloned grid members
  //


  for(dim=0;dim<GridRank; dim++){
    dxI[dim] = 1.0/CellWidth[dim][0]; 
  }

#ifdef HAOXU
if(ComovingCoordinates==1){
 FLOAT a, dadt;
    if (CosmologyComputeExpansionFactor(Time+(FLOAT)0.5*dtFixed, &a, &dadt)
        == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      Error = TRUE;
    }
 for(dim=0;dim<GridRank; dim++)
    dxI[dim] /= a;
}
#endif 

  Offset[0] = 1;
  Offset[1] = GridDimension[0];
  Offset[2] = GridDimension[0]*GridDimension[1];
  
  //Set up ghost zone array.  (each step has a different number of ghost zones.)
  //Note that this is relative to the zone that is fully updated!!!
  //<dbg>
  GhostZones[0] = ( (MHD_Recon[0] == 0 ) ? 1 : 3 );
  GhostZones[1] = 5;

  /*  should This be somewhere else?

  //Check that the number of boundary zones is consistent with 
  //what's needed for oblique runs.
  //Note that the automatic ghost zone checker isn't really done,
  //so this assumes 7 ghost zones.
  int GhostNeeded = 7;
  int ShiftGhostZones;
  
  if( MHDBlastNormal[0] != 0 ){
  for( dim=1;dim<GridRank;dim++){
  ShiftGhostZones = GhostNeeded + (DomainRightEdge[dim] - DomainLeftEdge[dim]) 
  *(MHDBlastNormal[dim]/(MHDBlastNormal[0]*CellWidth[dim][0]) );
  
  if( DEFAULT_GHOST_ZONES < ShiftGhostZones ){
  fprintf(stderr,"Athena Setup ERROR: Not enough ghost zones.  To accomidate ShiftedPeriodic, and to avoid\n");
  fprintf(stderr,"       any more MPI routines (costly to code) there must be enough ghost zones\n");
  fprintf(stderr,"       to accomidate both the solver and the shift. (also note, if not using\n");
  fprintf(stderr,"       the athena solver, bad things might happen anyways.)\n");
  fprintf(stderr,"       In your case, you need %d, but have %d\n", ShiftGhostZones,
  DEFAULT_GHOST_ZONES);
  Error = TRUE;
  }
  }//dim
  }//blast normal ghost zone check.
  if( DEFAULT_GHOST_ZONES < GhostNeeded ){
  fprintf(stderr,"Athena Setup Error: not enough ghost zones.  Please recompile.\n");
  Error = TRUE;
  
  }
  */


  NumberOfMHDFluxes = ( EquationOfState == 0 ) ? 7 : 6;


  //<dbg>
  ath_counter = 0;
  //</dbg>

  //
  //Set up the map between physics and code arrays.
  //

  //State 
  //Don't change these, things will crash.  They have a necessary relation
  //that relates to the Athena Eigensystem code.

  //int D,TE, VX, VY, VZ, BX,BY,BZ, GE;
  //int V[3], B[3];
  
  if( EquationOfState == 0 ){
    S.TE = 4;
    S.D = 0;
    S.V[0]  = 1;
    S.V[1]  = 2;
    S.V[2]  = 3;
    S.B[0]  = 5;
    S.B[1]  = 6;
    S.B[2]  = 7;
  }else if( EquationOfState == 1 ){
    //S.TE = 4; no energy variable needed.
    S.TE = -1;
    S.D = 0;
    S.V[0]  = 1;
    S.V[1]  = 2;
    S.V[2]  = 3;
    S.B[0]  = 4;
    S.B[1]  = 5;
    S.B[2]  = 6;
  }
  /* Won't need this for now.
  MHD_PLM_SlopeLocal[ S.D ]  = MHD_PLM_Slope[ 0 ];
  MHD_PLM_SlopeLocal[ S.V[0] ] = MHD_PLM_Slope[1];
  MHD_PLM_SlopeLocal[ S.V[1] ] = MHD_PLM_Slope[1];
  MHD_PLM_SlopeLocal[ S.V[2] ] = MHD_PLM_Slope[1];
  MHD_PLM_SlopeLocal[ S.B[0] ] = MHD_PLM_Slope[2];
  MHD_PLM_SlopeLocal[ S.B[1] ] = MHD_PLM_Slope[2];
  MHD_PLM_SlopeLocal[ S.B[2] ] = MHD_PLM_Slope[2];
  if( EquationOfState == 0 ) MHD_PLM_SlopeLocal[ S.TE ] = MHD_PLM_Slope[3];
  
  if( MHD_PLM_Slope[4] != -1 ){
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," Too many arguments given to MHD_PLM_Slope in parameter file.\n");
    fprintf(stderr," Need to update that parameter with the new version.  {d,v,b,e}, all v's and b's treated the same.\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!\n");
      Error = TRUE;

  }
  */
  // takes an enzo index, returns a solver index.


  MapEtoS[0][E.D] = S.D;
  MapEtoS[0][E.V[0]] =  S.V[0];
  MapEtoS[0][E.V[1]] =  S.V[1];
  MapEtoS[0][E.V[2]] =  S.V[2];
  if( EquationOfState == 0 )MapEtoS[0][E.TE] = S.TE;
  //This is wrong.  At
  MapEtoS[0][E.B[0]] = S.B[0];
  MapEtoS[0][E.B[1]] = S.B[1];
  MapEtoS[0][E.B[2]] = S.B[2];


  MapEtoS[1][E.D] = S.D;
  MapEtoS[1][E.V[0]] =  S.V[2];
  MapEtoS[1][E.V[1]] =  S.V[0];
  MapEtoS[1][E.V[2]] =  S.V[1];
  if( EquationOfState == 0 )MapEtoS[1][E.TE  ] = S.TE;
  MapEtoS[1][E.B[0]] = S.B[2];
  MapEtoS[1][E.B[1]] = S.B[0];
  MapEtoS[1][E.B[2]] = S.B[1];

  MapEtoS[2][E.D] = S.D;
  MapEtoS[2][E.V[0]] =  S.V[1];
  MapEtoS[2][E.V[1]] =  S.V[2];
  MapEtoS[2][E.V[2]] =  S.V[0];
  if( EquationOfState == 0 )MapEtoS[2][E.TE  ] = S.TE;
  MapEtoS[2][E.B[0]] = S.B[1];
  MapEtoS[2][E.B[1]] = S.B[2];
  MapEtoS[2][E.B[2]] = S.B[0];
  

  //Map physical vector fields to the one the solver sees.  
  //Cyclic permutaion.  I try to use cyclic permutation whenever I
  //map from three to one  (denoted B2 and B3, repsectively.)
  // BNum[dim][MagneticField]
  
  BNum[0][0] = 1;
  BNum[0][1] = 2;
  BNum[1][0] = 2;
  BNum[1][1] = 0;
  BNum[2][0] = 0;
  BNum[2][1] = 1;


  /* I should also clean up the DivB model.
    switch( MHD_DivB){
    
    default:
    case MHD_DivB_Poisson:    
    case MHD_DivB_none:
    //If the selected method isn't a CT method, it must be a cell
    //centered method, in which case the FaceCentered field is
    //ignored, and the CellCentered Longitudinal field is reconstructed (so NumberOfFluxes increases)
    
    AllCellCentered = 1;
    NumberReconstructed = NumberOfMHDFluxes + 1;
    break;
    case MHD_DivB_Balsara:
    case MHD_DivB_Athena:
    case MHD_DivB_RJ:
    AllCellCentered = 0;
    NumberReconstructed = NumberOfMHDFluxes;
    
    break;
    
    }
  */
  FlatteningField = NULL;


  Error = FALSE;
}

