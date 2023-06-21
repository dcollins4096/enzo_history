
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "MHD_Athena.h"

#ifdef ATHENA

#ifdef HAOXU
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
#endif 

//
// Sets up parameters for athena:
// Offset[3]                : distance in memory between neighboring points in the {x,y,z} direction.
// dxI[3]                   : 1/d{x,y,z}
// GhostZones[2]            : # ghost zones for the first and second pass.
// AthStart[3], AthEnd[3]   : start and end indicies.  Changed for 1st and second pass.
// NumberOfMHDFluxes        : 7 for adiabatic, 6 for isothermal.
// Eden, Egas, Ev[3], Eeng  : indicies for BaryonField.
// Sden, Sv[3], Seng        : indicies for Solver.
// MHD_PLM_SlopeLocal       : maps the input parameter to the useful parameter.
// MapEtoS                  : BaryonField[ field ] = State[  MapEtoS[field] ]
// BNum[dim][transvers0,1]  : transvers B component indecies for each face.
// AllCellCentered          : flag to denote that only the cell centered field is used.
// NumberReconstructed      : Number of fluxes, +1 if  B_{face} isn't used.
//

int grid::MHD_AthenaSetup(){
  
  //
  // Set up parameters for Athena.
  //

  //
  // For isothermal EOS, ReconField[1] isn't necessary, since there's no internal
  // energy field to reconstruct on.  So as a safety check, set it to -1.
  // I need a better way, in the long run, to deal with this.  
  //
  if( EquationOfState == 1 )
    MHD_ReconField[1] = -1;

  int dim;

  for(dim=0;dim<GridRank; dim++)

    dxI[dim] = 1.0/CellWidth[dim][0]; // *1/a

#ifdef HAOXU
if(ComovingCoordinates==1){
 FLOAT a, dadt;
    if (CosmologyComputeExpansionFactor(Time+(FLOAT)0.5*dtFixed, &a, &dadt)
        == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
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
	return FAIL;
      }
    }//dim
  }//blast normal ghost zone check.
  if( DEFAULT_GHOST_ZONES < GhostNeeded ){
    fprintf(stderr,"Athena Setup Error: not enough ghost zones.  Please recompile.\n");
    return FAIL;
  }
  //
  for(dim=0;dim<3;dim++){
    AthStart[dim] = GridStartIndex[dim];
    AthEnd[dim]   = GridEndIndex[dim];
  }

  NumberOfMHDFluxes = ( EquationOfState == 0 ) ? 7 : 6;

  //<dbg>
  ath_counter = 0;
  //</dbg>

  //
  //Set up the map between physics and code arrays.
  //

  //Enzo raw numbers
  if (this->IdentifyPhysicalQuantities(Eden, Egas, Ev[0], Ev[1], 
				       Ev[2], Eeng) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }


  //Don't change these, things will crash.

  if( EquationOfState == 0 ){
    Seng = 4;
    Sden = 0;
    Sv[0]  = 1;
    Sv[1]  = 2;
    Sv[2]  = 3;
    Sb[0]  = 5;
    Sb[1]  = 6;
    Sb[2]  = 7;
  }else if( EquationOfState == 1 ){
    //Seng = 4; no energy variable needed.
    Sden = 0;
    Sv[0]  = 1;
    Sv[1]  = 2;
    Sv[2]  = 3;
    Sb[0]  = 4;
    Sb[1]  = 5;
    Sb[2]  = 6;
  }
  MHD_PLM_SlopeLocal[ Sden ]  = MHD_PLM_Slope[ 0 ];
  MHD_PLM_SlopeLocal[ Sv[0] ] = MHD_PLM_Slope[1];
  MHD_PLM_SlopeLocal[ Sv[1] ] = MHD_PLM_Slope[1];
  MHD_PLM_SlopeLocal[ Sv[2] ] = MHD_PLM_Slope[1];
  MHD_PLM_SlopeLocal[ Sb[0] ] = MHD_PLM_Slope[2];
  MHD_PLM_SlopeLocal[ Sb[1] ] = MHD_PLM_Slope[2];
  MHD_PLM_SlopeLocal[ Sb[2] ] = MHD_PLM_Slope[2];
  if( EquationOfState == 0 ) MHD_PLM_SlopeLocal[ Seng ] = MHD_PLM_Slope[3];
  
  if( MHD_PLM_Slope[4] != -1 ){
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," Too many arguments given to MHD_PLM_Slope in parameter file.\n");
    fprintf(stderr," Need to update that parameter with the new version.  {d,v,b,e}, all v's and b's treated the same.\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr," !!!!!!!!!!!!!!!!!!!!!!!\n");
    return FAIL;
  }
  /* <dbg>
      fprintf(stderr,"stew ");
      for( dim=0; dim<4; dim++)
      fprintf(stderr," %d ", MHD_PLM_Slope[dim]);
      fprintf(stderr,"\n");
      fprintf(stderr,"stew ");
      for( dim=0; dim<8; dim++)
      fprintf(stderr," (%d,%d) ", MHD_PLM_SlopeLocal[dim], MHD_PLM_Slope!![dim]);
      fprintf(stderr,"\n");

      return FAIL;

      stew  1  3  1  1 
      WriteAllData: writing file data0001.
      stew  (1,1)  (3,1)  (3,3)  (3,3)  (1,3)  (1,1)  (1,1)  (1,1) 
      MHD_Athena: Athena Setup failed.
      </dbg>
  */

  // takes an enzo index, returns a solver index.


  MapEtoS[0][Eden] = Sden;
  MapEtoS[0][Ev[0]] =  Sv[0];
  MapEtoS[0][Ev[1]] =  Sv[1];
  MapEtoS[0][Ev[2]] =  Sv[2];
  if( EquationOfState == 0 )MapEtoS[0][Eeng] = Seng;

  MapEtoS[1][Eden] = Sden;
  MapEtoS[1][Ev[0]] =  Sv[2];
  MapEtoS[1][Ev[1]] =  Sv[0];
  MapEtoS[1][Ev[2]] =  Sv[1];
  if( EquationOfState == 0 )MapEtoS[1][Eeng  ] = Seng;

  MapEtoS[2][Eden] = Sden;
  MapEtoS[2][Ev[0]] =  Sv[1];
  MapEtoS[2][Ev[1]] =  Sv[2];
  MapEtoS[2][Ev[2]] =  Sv[0];
  if( EquationOfState == 0 )MapEtoS[2][Eeng  ] = Seng;
  
  

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

  //<dbg> Useful for debugging PPML
  //AllCellCentered = 1;
  //NumberReconstructed = NumberOfMHDFluxes + 1;
  //fprintf(stderr,"Reconstruction Changes for PPML debugging.\n");
  //</dbg>

  //Check that Hancock1 uses conserved reconstruction.
  if( MHD_Recon[0] == 2 || MHD_Recon[1] == 2 ){
    if( MHD_ReconField[0] != 1 ){
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING MHD_ReconField[0] must be 1 for MHD_Recon = 2 (Hancock1)\n");
      fprintf(stderr,"WARNING Forcing.\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      MHD_ReconField[0] = 1;
    }
    if(EquationOfState == 0 &&  MHD_ReconField[1] != 0 ){
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING MHD_ReconField[1] must be 0 for MHD_Recon = 2 (Hancock1)\n");
      fprintf(stderr,"WARNING Forcing this to be true.  It was %d\n",MHD_ReconField[1]);
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      MHD_ReconField[1] = 0;
    }
  }

  //Check that Hancock2 uses primitive reconstruction.
  if( MHD_Recon[0] == 3 || MHD_Recon[1] == 3 ){
    if( MHD_ReconField[0] != 0 ){
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING MHD_ReconField[0] must be 0 for MHD_Recon = 3 (Hancock2)\n");
      fprintf(stderr,"WARNING Forcing.\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      MHD_ReconField[0] = 0;
    }
    if(EquationOfState == 0 &&  MHD_ReconField[1] != 1 ){
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING MHD_ReconField[1] must be 1 for MHD_Recon = 3 (Hancock2)\n");
      fprintf(stderr,"WARNING Forcing this to be true.  It was %d\n",MHD_ReconField[1]);
      fprintf(stderr,"WARNING\n");
      fprintf(stderr,"WARNING\n");
      MHD_ReconField[1] = 1;
    }
  }

  FlatteningField = NULL;

  return SUCCESS;
}
#endif //ATHENA
