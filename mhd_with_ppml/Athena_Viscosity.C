#ifdef PPML
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


#ifndef  PPML
#include "pout.h"
#endif //PPML

#include "Athena.h"

// 
// Diffusion mechanism.
// Used initially for Athena-- Stolen for PPM-L because I'm lazy.
//
// F => F - \nu * (U_{i+1} - U{i})
// \nu = -min(\Div V, 0)
// Returns Success,
//      or -1, (if failure)
//      or \nu (hense the -1: often \nu = 0)
// Alters the Flux array.  If Flux passed in is null,
// then \nu is returned, instead.
//
// Note that this is a multi use code: check all instances before modifying.
// In particular, in MHD_dtVisc, the internal energy variable isn't set, so
// dont' go using that here without fixing it.

// map of indicies:   
//     space:         | i-2 | i-1 | i | i+1 | i+2 | 
//     index[];       |  3  |  1  | 0 |  2  |  4  |
//                    Swap i with j or k as needed.
//                    This swap is controlled by the Offset variable.
//                    cell(i  ,j+1,k  ) = i + nx*(j+1) + nx*ny(k)  = i + nx*j + ny*k + nx
//                                      = index[0] + nx
//                    cell(i  ,j  ,k+2) = i + nx*j + nx*ny*(k+2) 
//                                      = index[0] + 2*nx*ny 
//                    etc.
//    This map was chosen for extensibility:  if a new solver with a more complicated stencil is 
//    written, index[] can be added to without altering the structure of existing code.

//    ----------------- 
//    |       |       | 
//    |       |       | 
//    |  Y3   |   Y1  |  j+1
//    |       |       | 
//    |       |       | 
//    -----------------   
//    |       |       | 
//    |       |       | 
//    |    Flux (i,j) |  j
//    |       |       | 
//    |       |       | 
//    -----------------   
//    |       |       | 
//    |       |       | 
//    |  Y4   |   Y2  |  j-1
//    |       |       | 
//    |       |       | 
//    ----------------- 
//       i-1      i
//
// Offset Variables:
//   Dc[dim][direction]:  For the [dim] axis, gives the memory offset of one
//                        index in the next two cyclic directions.
//                        So for the x face, Dc[0][0] is the memory offset of y+1 (nx)
//                           for the z face, Dc[2][0] is the memory offset of x+1 (1)


//
// NOTICE!!!!!!!
// NOTICE!!!!!!!
// NOTICE!!!!!!!
// 
// This has been made multi purpose code.  If you're goin going to change anything, 
// it would be wise to also ensure that its consistent with the call in Grid_MHD_dtVisc.C

//Slope limiters: see MHD_PLM for details.
inline float MinMod(float dL,float dR){ 
  return ( dL*dR < 0 ) ? 0  : ( ( fabs( dL ) < fabs( dR ) ) ? dL : dR ) ;
}


float Athena::MHD_Viscosity(float * Fluxes, float * Lhs, float * Rhs, int * index, 
			  int dim, int MHD_DiffusionMethodInput, float dT){

  // Index variables.  Offsets, indicies for cross-velocities.
  int Dc[3][2] =       
    {{GridDimension[0],GridDimension[0]*GridDimension[1]}, 
     {GridDimension[0]*GridDimension[1],1},
     {1,GridDimension[0]}};    
  
  //This is some counter intuitive index naming.
  //It's kind of upside down: increasing numbers indicates decreasing indicies.  Hmmm...

  int 
    Y1 = index[0] + Dc[dim][0], // (i  ,j+1,k  )
    Y2 = index[0] - Dc[dim][0], // (i  ,j-1,k  )
    Y3 = index[1] + Dc[dim][0], // (i-1,j+1,k  )
    Y4 = index[1] - Dc[dim][0], // (i-1,j-1,k  )
    Z1 = index[0] + Dc[dim][1], // (i  ,j  ,k+1)
    Z2 = index[0] - Dc[dim][1], // (i  ,j  ,k-1)
    Z3 = index[1] + Dc[dim][1], // (i-1,j  ,k+1)
    Z4 = index[1] - Dc[dim][1], // (i-1,j  ,k-1) (and cyclic permutations therin.)
    //And the diagonal terms
    D11= index[0] + Dc[dim][0] + Dc[dim][1], // (i  ,j+1,k+1)
    D12= index[0] - Dc[dim][0] + Dc[dim][1], // (i  ,j-1,k+1)n
    D13= index[0] + Dc[dim][0] - Dc[dim][1], // (i  ,j+1,k-1)n
    D14= index[0] - Dc[dim][0] - Dc[dim][1], // (i  ,j-1,k-1)n
    D21= index[1] + Dc[dim][0] + Dc[dim][1], // (i-1,j+1,k+1)n
    D22= index[1] - Dc[dim][0] + Dc[dim][1], // (i-1,j-1,k+1)n
    D23= index[1] + Dc[dim][0] - Dc[dim][1], // (i-1,j+1,k-1)n
    D24= index[1] - Dc[dim][0] - Dc[dim][1]; // (i-1,j-1,k-1)n
  
  //Which transerse coordinate to choose.
  //For 1d, only x gradients exist.
  //For 2d, 'transverse' gradents are in the y direction for F_x, in the X direciton for F_y.
  //For 3d, all 3 fluxes have 2 transverse directions.
  //Sticking with the cyclic permutation used through the New MHD code, transverse 1 is first cyclic, 2 is the second.

  int Transverse1 = 1;
  int Transverse2 = 1;
  if( GridRank == 1 ){
    Transverse1 = 0;
    Transverse2 = 0;
  }
  if( GridRank == 2 ){
    Transverse1 = ( (dim == 0 ) ? 1 : 0 );
    Transverse2 = ( (dim == 0 ) ? 0 : 1 );
  }
  //select transverse axii
  int CyclicAxis[3][3]= { {0,1,2},{1,2,0}, {2,0,1}};
  int Taxis1 = CyclicAxis[dim][1];
  int Taxis2 = CyclicAxis[dim][2];
  int field;
  //The actual coefficient of diffusion:
  float nu = 0.0;

  // If this has been called from MHD_dtVisc, it doens't have LHS or RHS computed.

  // So check that only method 3 is being passed in (for now.)
  if( Lhs == NULL || Rhs == NULL )
    if(  MHD_DiffusionMethodInput != 0 && 
	 MHD_DiffusionMethodInput != 3 &&
	 MHD_DiffusionMethodInput != 7){
      fprintf(stderr, " Error: Call to MHD_Diffusion with Null Left and Right states, ");
      fprintf(stderr," and improper diffusion method (%"ISYM").\n",MHD_DiffusionMethodInput);
      return -1;
    }
  

  switch( MHD_DiffusionMethodInput ){

  case 0:
    return SUCCESS;
    break;

  case 1:
    nu = Lhs[S.V[0]] - Rhs[S.V[0] ]
      + ( (Transverse1 == 1 ) ? 
	  (dxI[Taxis1] * CellWidth[dim][0] * 0.25 * ( ( State[S.V[1]][ Y2 ] - State[ S.V[1] ][ Y1 ] )
						      +(State[S.V[1]][ Y4 ] - State[ S.V[1] ][ Y3 ] ) ) )
	  : 0 )
      + ( (Transverse2 == 1 ) ? 
	  (dxI[Taxis2] * CellWidth[dim][0] * 0.25 * ( ( State[S.V[2]][ Z2 ] - State[ S.V[2] ][ Z1 ] )
						      +(State[S.V[2]][ Z4 ] - State[ S.V[2] ][ Z3 ] ) ) )
	  : 0 );


    nu = max(nu,0);
    break;
  case 2:
    nu = 1.0;
    break;

  case 3:
    nu = State[ S.V[0] ][ index[1] ] - State[ S.V[0] ][ index[0] ]
      + ( (Transverse1 == 1 ) ? 
	  (dxI[Taxis1] * CellWidth[dim][0] * 0.25 * ( ( State[S.V[1]][ Y2 ] - State[ S.V[1] ][ Y1 ] )
						      +(State[S.V[1]][ Y4 ] - State[ S.V[1] ][ Y3 ] ) ) )
	  : 0 )
      + ( (Transverse2 == 1 ) ? 
	  (dxI[Taxis2] * CellWidth[dim][0] * 0.25 * ( ( State[S.V[2]][ Z2 ] - State[ S.V[2] ][ Z1 ] )
						      +(State[S.V[2]][ Z4 ] - State[ S.V[2] ][ Z3 ] ) ) )
	  : 0 );


    nu = max(nu,0);
    break;

  case 4:
    nu = Lhs[S.V[0]] - Rhs[S.V[0] ]
      + ( (Transverse1 == 1 ) ? 
	  (dxI[Taxis1] * CellWidth[dim][0] * 0.25 * ( ( State[S.V[1]][ Y2 ] - State[ S.V[1] ][ Y1 ] )
						      +(State[S.V[1]][ Y4 ] - State[ S.V[1] ][ Y3 ] ) ) )
	  : 0 )
      + ( (Transverse2 == 1 ) ? 
	  (dxI[Taxis2] * CellWidth[dim][0] * 0.25 * ( ( State[S.V[2]][ Z2 ] - State[ S.V[2] ][ Z1 ] )
						      +(State[S.V[2]][ Z4 ] - State[ S.V[2] ][ Z3 ] ) ) )
	  : 0 );


    nu = fabs(nu);
    break;


  case 5:
    nu = Lhs[S.V[0]] - Rhs[S.V[0] ]
      + ( (Transverse1 == 1 ) ? 
	  (dxI[Taxis1] * CellWidth[dim][0] * 0.5 * MinMod( ( State[S.V[1]][ Y2 ] - State[ S.V[1] ][ Y1 ] ),
							   (State[S.V[1]][ Y4 ] - State[ S.V[1] ][ Y3 ] ) ) )
	  : 0 )
      + ( (Transverse2 == 1 ) ? 
	  (dxI[Taxis2] * CellWidth[dim][0] * 0.5 * MinMod( ( State[S.V[2]][ Z2 ] - State[ S.V[2] ][ Z1 ] ),
							   (State[S.V[2]][ Z4 ] - State[ S.V[2] ][ Z3 ] ) ) )
	  : 0 );


    nu = fabs(nu);
    break;

  case 6:
    nu = Lhs[S.V[0]] - Rhs[S.V[0] ];    
    nu = max(nu,0);
    break;

  case 7:
    //   From J. Saltzman, J. Comp Phys, 115, 153-168 (1994)
    // p. 159.  Same basic formulation of the viscosity as in Colella Woodward 1984 (the ppm paper)
    // F_{i+1/2} -> F_{i+1/2} - nu_{i+1/2}( U_{i} - U_{i-1} )
    //   Here the viscosity coefficient at the face, nu_{i+1/2}, (called C_{i+1/2} in the paper)
    // is the average of the divercence computed at the corners of the face (only the negative terms)
    // and the corner centered divergence is computed by the differences abutting the corner 
    //   I believe there to be a typo in the paper, as the way it's written there gives
    // C_{i+1/2, j+1/2, k+1/2} involving Vx_{i} - Vx{i-1}, which isn't centered at i+1/2. 
    // I believe the easiest way around this descrepancy is reworking the finite difference
    // \Delta_i to be a forward difference instead of backwards.  I have done that here.
    // I have also algebraically reduced the differences.

    nu =
      // d(Vx)/dx 
      dxI[dim]*
      (4*(State[ S.V[0] ][ index[0] ]- State[ S.V[0] ][ index[1] ])
       +2*( ( State[ S.V[0] ][ Y1 ] + State[ S.V[0] ][ Y2 ] + State[ S.V[0] ][ Z1 ] + State[ S.V[0] ][ Z2 ])
	    -(State[ S.V[0] ][ Y3 ] + State[ S.V[0] ][ Y4 ] + State[ S.V[0] ][ Z3 ] + State[ S.V[0] ][ Z4 ]))
       +( ( State[ S.V[0] ][ D11 ] + State[ S.V[0] ][ D12 ] + State[ S.V[0] ][ D13 ] + State[ S.V[0] ][ D14 ])
	  -(State[ S.V[0] ][ D21 ] + State[ S.V[0] ][ D22 ] + State[ S.V[0] ][ D23 ] + State[ S.V[0] ][ D24 ]) )
       )

    // d(Vy)/dy

      +dxI[Taxis1]*
      (2*( State[ S.V[1] ][ Y1 ] - State[ S.V[1] ][ Y2 ] )
       +( State[ S.V[1] ][ D11 ] - State[ S.V[1] ][ D12 ] )
       +( State[ S.V[1] ][ D13 ] - State[ S.V[1] ][ D14 ] )
       +2*( State[ S.V[1] ][ Y3 ] - State[ S.V[1] ][ Y4 ] )
       +( State[ S.V[1] ][ D21 ] - State[ S.V[1] ][ D22 ] )
       +( State[ S.V[1] ][ D23 ] - State[ S.V[1] ][ D24 ] )
       )

      // d(Vz)/dz
      +dxI[Taxis2]*
      (2*( State[ S.V[2] ][ Z1 ] - State[ S.V[2] ][ Z2 ] )
       +( State[ S.V[2] ][ D11 ] - State[ S.V[2] ][ D13 ] )
       +( State[ S.V[2] ][ D12 ] - State[ S.V[2] ][ D14 ] )
       +2*( State[ S.V[2] ][ Z3 ] - State[ S.V[2] ][ Z4 ] )
       +( State[ S.V[2] ][ D21 ] - State[ S.V[2] ][ D23 ] )
       +( State[ S.V[2] ][ D22 ] - State[ S.V[2] ][ D24 ] )
       )
      ;

    
    // Renormalize:
    // The CellWidth will be devided out again during the Flux Difference.
    // The 16 is from 4 corners, and 4 differences averaged for each corner.
    
    nu *= CellWidth[dim][0]/16.0;  
    nu = -1*min( nu, 0 );  
    
  }//switch

  //<new>
  
  
  if( nu > 0 ) 
    if( dT > CellWidth[dim][0]/(State[ S.D ][ index[0] ] * MHD_DiffusionParameter * 2 * nu ) )
      fprintf(stderr," visc: dt/stuff= %e > 1\n",
	      dT/(CellWidth[dim][0]/(State[ S.D ][ index[0] ] * MHD_DiffusionParameter * 2 * nu ) ) );
  
  //</new>
  int FluxIndex;
  if( Fluxes != NULL ){
    if( HydroMethod != PPM_Local ){
      for(field = 0; field< NumberOfBaryonFields; field++)
	Fluxes[ MapEtoS[dim][field] ] -= MHD_DiffusionParameter * nu *
	  (State[ MapEtoS[dim][field] ][ index[0] ] - State[ MapEtoS[dim][field] ][ index[1] ] );
      
    }else{ //HydroMethod == PPM_Local
      //Cache Viscous.  
      for( field = 0; field < NumberOfFluidQuantities - ((MHD_Used) ? 3 : 0); field++){
	FluxIndex = index[1] + GridDimension[0]*GridDimension[1]*GridDimension[2] * field;
	Fluxes[ FluxIndex ]  -= MHD_DiffusionParameter * nu *
	  (State[ MapEtoS[dim][field] ][ index[0] ] - State[ MapEtoS[dim][field] ][ index[1] ] );
      }
    }

  }else{ //Fluxes == NULL
    return nu * MHD_DiffusionParameter;
  }

  return SUCCESS;
}

#endif //PPML
