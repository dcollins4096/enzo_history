
//
// Compute the Electric Field using the algorithm presented by GS04, with some modification.
// In brief, this routine computes the Electric Field from fluxes returned from the solver.
// Additionally it employs a higher order, possibly upwinded (depending on parameter choices) corrections.

// Note that in Enzo, the field called ElectricField is actually
// dT*E.  Putting the timestep here, instead of in the curl, allows
// for our MHD analogue of flux correction.  

// I exploit the cyclic permutation symmetry of the MHD equations pretty heavily in this code:  only one 
// electric field is computed, but written in a generic enough way as to be looped over.  



// Variables:
// -------------------
//
// Cyclic variables:
//    dimX, dimY, dimZ
//   These variables are cyclicly permuted.  dimX is looped over and
//   indexed the electric field.  
//

// Field index variables:  These variables relate the position in the
//   algorithm to position in the 1d array.  
//   There are 3 different array structures used:  Cell Centered, Face
//   centered, and a second type of face centered.
 
//    C1,2,3,4:  cell centered variables:  velocity and Centered
//               Magnetic Field
//    F1,2,3,4:  Fluxes: face centered.
//    B1,2,3,4:  Face centere Magnetic Field.  Not that the fluxes and
//               the face centered magnetic field use different data
//               structure shapes, so require different indicies.
//    The diagram below shows the spatial location of each index,
//    relative to the electric field being calculated (at the center)
//  
//  
//    ----------------- 
//    |       |       | 
//    |       |       | 
//    |  C2   F2  C1  | 
//    |       |       | 
//    |       |       | 
//    ---F3---E---F1---   
//    |       |       | 
//    |       |       | 
//    |  C4   F4  C3  | 
//    |       |       | 
//    |       |       | 
//    ----------------- 

// Offset Variables:
//   Db[dim][direction], Dc[dim][direction]
//    direction are the other two variables, in cyclic permutation.
//    The one dimensional nature of the actual data structures used in
//    Enzo break the cyclic symitry of the MHD equations, so Dc and Db
//    restore that symmetry, so I can loop over the field. Dc is for
//    the centered varaibles, Db is for the Face Centered Magnetic Field.
//    Note that (just to make things difficult), the Flux is stored in an
//    array the same shape as the Centered Variables, so offsets are
//    computed like the Cell Centered Field, not the face centered field.
//    To compute these, index maps for (i,j,k) and (i-1,j-1,k) were
//    written down, and the difference was taken.  For instance, the
//    difference in memory between (i,j,k) and (i,j,k-1) is
//    i+nx(j+ny k) - i+nx(j+ny(k-1) = nx*ny
//      Note that the face centered magetic offset variables only index the
//    Magnetic Field compontents that are in a given plane.

// Macros and Pointers:
//    Enzo uses irritatingly long names for things.  These pointers
//    make my life as a programmer easier.
//

// Weight parameters:
//  As,Bs,Cs:
//   Needs some kind of upwinding.  I'm thinking wind direction, but not
//   clear what's best yet.

// EndX, EndY, EndZ
//   For Rank<3, the size of the cube looped over needs to be modified.

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
#include "MHD_Athena.h"


#ifdef ATHENA


inline void SSSe(float * BsL, float * BsR, float WindDirection){       
  
  float Tolerance = 1e-10;

  if( WindDirection > Tolerance ){
    *BsR = 0; *BsL=2;
  }else if( WindDirection < -1*Tolerance){
    *BsR = 2; *BsL=0;
  }else{
    *BsR = 1; *BsL=1;
  }

}

int grid::MHD_AthenaElectric(float dT, float * Fluxes[]){


  //loop variables.
  int i,j,k,dimX, dimY, dimZ, flag;
  int EndX=1, EndY=1,EndZ=1;

  //Permutation variables.  
  int C1, C2, C3, C4;                   // Cell centered fields.
  int B1, B2, B3, B4;                   // Face centered Mag. Field.
  int F1, F2, F3, F4;                   // Face centered flux.
  int Edex;                             // Corner centered field being updated.

  //inverst timestep.
  float dTi = 1.0/dT;
   // Cell Centered Offsets.  Used for cell centered quantities as well as fluxes.
   int Dc[3][2] =       
     {{GridDimension[0],GridDimension[0]*GridDimension[1]}, 
      {GridDimension[0]*GridDimension[1],1},
      {1,GridDimension[0]}};    

   // Face centered offsets.  Used for the face centered magnetic field
   int Db[3][2] =        
     {{GridDimension[0],GridDimension[0]*(GridDimension[1] + 1)}, 
      {(GridDimension[0]+1)*GridDimension[1],1},
      {1,GridDimension[0]+1}};    


   //Pointers (and one macro) for brevity.
   //note that Ev[] is defined in Grid_MHD_Athena, and is equal to the
   //velocity indicies.
   //Also note that dimY and dimZ (used in Ec) are changed during the
   //relevant loop.
   float * Bc[3] = {CenteredB[0], CenteredB[1],CenteredB[2]};
   float * Bf[3] = {MagneticField[0], MagneticField[1], MagneticField[2]};
   float * Vel[3] = {BaryonField[ Ev[0] ], BaryonField[ Ev[1] ], BaryonField[ Ev[2] ] };
 #define Ec(index) (Bc[ dimY ][index]*Vel[ dimZ ][index] - Bc[ dimZ ][index]*Vel[ dimY ][index] )

   // Weights. Still working on the best way to go about this.  I'm
   // thinking some kind of convex combination of the 4 choices based
   // on wind speed and direction.  

   // 7/20/6: all these need to be positive for convex combination.
   // So something like term = (a[i] > 0 ) ? a[i]*stuff : 0;, or would that be slow as hell?
   float As[4] = {1.0,1.0,1.0,1.0};  
   float Bs[8] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
   float Cs[4] = {1.0,1.0,1.0,1.0};

   //
   // Alterations for 2.5, 1.5d
   // These are really just fail-safes, any fields that actually acces
   // differences along the 'homogenous' axii are removed
   //

   if( GridRank < 3 ){
     Dc[0][1] = 0;
     Dc[1][0] = 0;
     Db[0][1] = 0;
     Db[1][0] = 0;
     EndZ = 0;
   }
   if( GridRank < 2){
     Dc[0][0] = 0;
     Dc[2][1] = 0;
     Db[0][0] = 0;
     Db[2][1] = 0;

     EndY = 0;
   }


   //loop; space, dim.

   for(k=AthStart[2]; k<=AthEnd[2]+EndZ; k++)
   for(j=AthStart[1]; j<=AthEnd[1]+EndY; j++)
   for(i=AthStart[0]; i<=AthEnd[0]+EndX; i++){

     for( dimX=0; dimX< 3; dimX++){

       //Double check that only the used variables get looped over.

       if( (dimX == 0 && i == AthEnd[0] +1 )||
	   (dimX == 1 && j == AthEnd[1] +1 )||
	   (dimX == 2 && k == AthEnd[2] +1 ) ) 
	 continue;

       //the cyclic variables.
       dimY = (dimX == 0 ) ? 1 : (dimX == 1 ) ? 2: 0;
       dimZ = (dimX == 0 ) ? 2 : (dimX == 1 ) ? 0: 1;

       //The index variables.
       Edex = i + ElectricDims[dimX][0] * (j + ElectricDims[dimX][1]*k);
       C1 = i + GridDimension[0] * (j + GridDimension[1]*k);
       C2 = C1 - Dc[dimX][0];
       C3 = C1 - Dc[dimX][1];
       C4 = C1 - (Dc[dimX][0] + Dc[dimX][1]);

       F1 = C1 * NumberOfMHDFluxes + Sb[1];
       F2 = C1 * NumberOfMHDFluxes + Sb[0];
       F3 = (C1 - Dc[dimX][0])*NumberOfMHDFluxes + Sb[1];
       F4 = (C1 - Dc[dimX][1])*NumberOfMHDFluxes + Sb[0];

       B1 = i + MagneticDims[ dimZ ][0]*(j + MagneticDims[ dimZ ][1]*k);
       B2 = i + MagneticDims[ dimY ][0]*(j + MagneticDims[ dimY ][1]*k);
       B3 = B1 - Db[dimX][0];
       B4 = B2 - Db[dimX][1];

       //
       // Solution/Parameter Sensitive Switches:
       // The flags are for, in order,  -d/dy, +d/dy, -d/dz, +d/dz
       //     terms in the reconstruction.

       for(flag=0; flag<4;flag++){
	 As[flag]=1;
	 switch(MHD_ElectricRecon){
	 case 0:
	   Cs[flag]= 0;
	   Bs[2*flag]  = 0;
	   Bs[2*flag+1]= 0;
	   break;
	 case 1:
	   Cs[flag]=1;
	   Bs[2*flag]  = 1;
	   Bs[2*flag+1]= 1;
	   break;
	 case 2:
	   Cs[flag]= 0;

	   //The Bs for MHD_ElectricRecon = 2 is below the loop-- it just 
	   //wasn't smooth to have it in this loop.
	   break;
	 }
       }//flag loop

       // Set the electric component to use based on wind direction.

       if( MHD_ElectricRecon == 2 ){
	 // -d/dy
	 SSSe(Bs+1,Bs+0,Vel[dimZ][C1] + Vel[dimZ][C3]);
	 // +d/dy
	 SSSe(Bs+3,Bs+2,Vel[dimZ][C2] + Vel[dimZ][C4]);
	 // -d/dz
	 SSSe(Bs+5,Bs+4,Vel[dimY][C1] + Vel[dimY][C2]);
	 // +d/dz
	 SSSe(Bs+7,Bs+6,Vel[dimY][C3] + Vel[dimY][C4]);

       }
       
       //
       //The bulk.
       //
       
       ElectricField[dimX][Edex] = 0.25*( As[0]*Fluxes[ dimZ ][F1]
					  +As[1]*Fluxes[ dimZ ][F3]
					  -As[2]*Fluxes[ dimY ][F2]
					  -As[3]*Fluxes[ dimY ][F4] );
       
       //The electric linear correction
       ElectricField[dimX][Edex] +=0.125*(
					  
					  // dimY directed derivative
					  //Note the (-) on the fluxes[dimY]. See documentation for the reason.
					  
					  -Bs[0]*(Ec(C1) + Fluxes[dimY][F2])  
					  -Bs[1]*(Ec(C3) + Fluxes[dimY][F4]) 
					  +Bs[2]*(-Fluxes[dimY][F2] - Ec(C2) )
					  +Bs[3]*(-Fluxes[dimY][F4] - Ec(C4) )

					  // dimZ derivative.
					  -Bs[4]*(Ec(C1) - Fluxes[dimZ][F1] )
					  -Bs[5]*(Ec(C2) - Fluxes[dimZ][F3])
					  +Bs[6]*(Fluxes[dimZ][F1]- Ec(C3) )
					  +Bs[7]*(Fluxes[dimZ][F3] - Ec(C4) )
					  
					  );
       
       //The dissipative derivative correction.
       ElectricField[dimX][Edex] += 
	 0.125*MHD_DivBparam*CellWidth[0][0]*dTi*(
						 //Dissipation term for -dimY derivative
						 +Cs[0]*(Bc[dimY][C1] - Bf[dimY][B2]
							 -(Bc[dimY][C3] - Bf[dimY][B4]))
						 //for +dimY
						 -Cs[1]*(Bf[dimY][B2] - Bc[dimY][C2]
							 -(Bf[dimY][B4] - Bc[dimY][C4]))
						 //for -dimZ (the sign change comes out of the method.)
						 -Cs[2]*(Bc[dimZ][C1] - Bf[dimZ][B1] 
							 -(Bc[dimZ][C2]-Bf[dimZ][B3]) )
						 //for +dimZ
						 +Cs[3]*(Bf[dimZ][B1] - Bc[dimZ][C3]
							 -(Bf[dimZ][B3] -Bc[dimZ][C4]))
						 );
       
       // Finally, add the dT.  While conceptually it doesn't belong
       // attached to the electric field, its  here for AMR concerns.
       // (the flux correction)
       
       //if( dccCounter9++ < 5000 ) fprintf(stderr,"kludge: no dt\n");
       ElectricField[dimX][Edex] *= dT; // * 1/a
       
     }//dim
   }//i,j,k
   
   return SUCCESS;
}

#endif //ATHENA
