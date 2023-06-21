#include <math.h>
#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

//
//  Routines to facilitate the Oblique Shock generation.
//
//
int SetupNormal(float Normal[], float MHDBlastCenter[3], TopGridData & MetaData);
void RotateY(float * Vector, double Angle);
void RotateZ(float * Vector, double Angle);
float AreaWrapper(int i,int j,int k,int dim, int rank, float dx ,float dy, float dz,float Normal[]);
float VolumeWrapper(int i,int j,int k, int rank, float dx ,float dy, float dz,float Normal[]);
float AreaBelow(float pA[], float pB[], float pC[], float pD[], float Line[]);
float Triangle(float Base, float Height);
float Quadralateral(float pA[], float pB[], float pC[], float pD[]);


// Rotates a vector so the X component is along the normal.  
// Done with a rotation first (clockwise) about Z, then (counter clockwise) about the new Y axis.  
// This is because the X unit vector is, conceptually, normal to the plane
// Why? Easier to draw.
// Also note: as of 12/14/06, only single angle rotations are allowed.
// The general formulation is implimented for ease of updating later.  

void RotateVector( float* Vector, float* Normal){

  double SinTheta = Normal[2]/sqrt( Normal[0]*Normal[0]+Normal[1]*Normal[1] + Normal[2]*Normal[2]);
  double Theta = asin( SinTheta );
  double SinPhi   = Normal[1]/sqrt( Normal[0]*Normal[0]+Normal[1]*Normal[1]);
  double Phi   = asin( SinPhi );

  RotateZ( Vector, Phi);
  RotateY( Vector, Theta);
}
void RotateY( float * Vector, double Angle){

  double ct = cos( Angle );
  double st = sin( Angle );
  double matrix[3][3] = { { ct, 0,-st },
			  {  0, 1,  0 },
			  { st, 0, ct } };
  double temp[3] = {0,0,0};
  for( int i=0; i<3; i++)
    for( int j=0; j<3; j++)
      temp[i] += Vector[j] * matrix[i][j];
  for(int i = 0;i<3;i++)
    Vector[i] = temp[i];
}
void RotateZ( float * Vector, double Angle){

  double ct = cos( Angle );
  double st = sin( Angle );
  double matrix[3][3] = { { ct,-st, 0 },
			  { st, ct, 0 },
			  {  0, 0,  1 } };

  double temp[3] = {0,0,0};
  for( int i=0; i<3; i++)
    for( int j=0; j<3; j++)
      temp[i] += Vector[j] * matrix[i][j];
  for(int i = 0;i<3;i++)
    Vector[i] = temp[i];
}

//Area of triangle and quadralateral.
float Triangle( float Base, float Height ){
  return 0.5*Base*Height;
}
float Quadralateral( float pA[], float pB[], float pC[], float pD[]){
  return 0.5 * fabs( (pD[0] - pB[0])*(pC[1] - pA[1]) - (pD[1] - pB[1])*(pC[0] - pA[0]) );
}

// Given a line (described by Line[]) and a rectangle (defined by pA, pB, pC, pD, clockwise from lower left)
// return the fraction of the rectangle that lies below the line.
// The line is defined by Line[0]*X + Line[1]*Y + Line[2] = 0
// There are six cases to check: this was written with the assumption that the slope of the line 
// is negative.

float AreaBelow(float pA[], float pB[], float pC[], float pD[], float Line[]){

  float Area;
  float TotalArea = (pB[0]-pA[0])*(pC[1]-pA[1]);

  //Intersections of the Line and the line going through A[0], B[0], A[1], C[1]
  //For LAX,LBX the value is the Y coordinate of the intersection
  //LAY, LCY is the X coordinate of the intersection of the Y cell edge and the line.

  //For flat lines, watch out:

  float LAX[2], LBX[2], LAY[2], LCY[2];
  LAX[0] =  pA[0];
  LAX[1] = ( (Line[1] != 0 ) ? -(Line[0]*pA[0] + Line[2])/Line[1]: 0 );
  LBX[0] = pB[0];
  LBX[1] = ( (Line[1] != 0 ) ? -(Line[0]*pB[0] + Line[2])/Line[1]: 0 );
  LAY[0] = ( (Line[0] != 0 ) ? -(Line[1]*pA[1] + Line[2])/Line[0]: 0 );
  LAY[1] = pA[1];
  LCY[0] = ( (Line[0] != 0 ) ? -(Line[1]*pC[1] + Line[2])/Line[0]: 0 );
  LCY[1] = pC[1];

  
  if( Line[0] == 0 ){
    //
    // Horizontal Lines
    // 
    if( LAX[1] <= pA[1] ){
      Area = 0;
    }else{
      if( LAX[1] >= pC[1] ){
	Area = TotalArea;
      }else{
	Area = Quadralateral( pA, LAX, LBX, pB);
      }
    }
  }else if( Line[1] == 0 ){
    //
    // Vertical Lines
    //

    if( LAY[0] <= pA[0] ){
      Area = 0;
    }else{
      if( LAY[0] >= pB[0] ){
	Area = TotalArea;
      }else{
	Area = Quadralateral( pA, pC, LCY, LAY );
      }
    }
  }else if( Line[0] != 0 && Line[1] != 0 ){
    //
    // Slant Lines
    //

    if( LAX[1]  <= pA[1] ){ 
      //case 1
      Area = 0; //line starts below cell
    }else if ( LAX[1] > pA[1] ) {

      if( LAX[1] <= pC[1] ){
	if( LBX[1] <= pB[1] ){
	  //case 2
	  Area = Triangle( LAY[0] - pA[0], LAX[1] - pA[1]);
	}else if( LBX[1] > pB[1] ){
	  //case 3
	  Area = Quadralateral( pA, LAX, LBX, pB);
	} //lbx < by
	
      }else if( LAX[1] > pC[1] ) { 
	
	if( LBX[1] < pB[1] ){
	  //case 4
	  Area = Quadralateral( pA, pC, LCY, LAY);
	}else if( LBX[1] >= pB[1] ){
	  if( LBX[1] < pC[1] ){
	    //case 5
	    Area = TotalArea - Triangle( pD[0] - LCY[0], pD[1] - LBX[1] );
	  }else if ( LBX[1] >= pC[1] ){
	    //case 6
	    Area = TotalArea;
	  }
	}
      }//lax > cy
    }// lax<pa[1]
  }// lines[0] !=0, Line[1] != 0

  //<dbg>
  //fprintf(stderr,"area: %f total %f area/total %f\n", Area, TotalArea, Area/TotalArea);
  //</dbg>
  return Area / TotalArea;
}
float AreaWrapper(int i,int j,int k,int dim, int rank, float dx ,float dy, float dz,float Normal[]){

  // Passes the cell information to the bit that actuall computes the area.
  // Also converts i,j,k from "Cells From Data Wall" to "Distance From Domain Wall,"
  // since that's how the BlastCenter is defined.

  float pA[2], pB[2], pC[2], pD[2], Line[3];

  // convert to Cells From Domain Wall
  i = i - DEFAULT_GHOST_ZONES;
  j = j - DEFAULT_GHOST_ZONES;
  k = k - DEFAULT_GHOST_ZONES;
  
  //Just to make sure we don't get any extraneous zeros
  if( rank == 2 )
    dz = 1;
  switch( dim ){
  case 0:
    pA[0] = (j     ) * dy;
    pA[1] = (k     ) * dz;
    pB[0] = (j +1.0) * dy;
    pB[1] = (k     ) * dz;
    pC[0] = (j     ) * dy;
    pC[1] = (k +1.0) * dz;
    pD[0] = (j +1.0) * dy;
    pD[1] = (k +1.0) * dz;
    Line[0] = Normal[1];
    Line[1] = Normal[2];
    Line[2] = Normal[0] * (i      ) * dx + Normal[3];
    break;
  case 1:
    pA[0] = (i     ) * dx;
    pA[1] = (k     ) * dz;
    pB[0] = (i +1.0) * dx;
    pB[1] = (k     ) * dz;
    pC[0] = (i     ) * dx;
    pC[1] = (k +1.0) * dz;
    pD[0] = (i +1.0) * dx;
    pD[1] = (k +1.0) * dz;
    Line[0] = Normal[0];
    Line[1] = Normal[2];
    Line[2] = Normal[1] * (j      ) * dy + Normal[3];
    break;
  case 2:
    pA[0] = (i      ) * dx;
    pA[1] = (j      ) * dy;
    pB[0] = (i + 1.0) * dx;
    pB[1] = (j      ) * dy;
    pC[0] = (i      ) * dx;
    pC[1] = (j + 1.0) * dy;
    pD[0] = (i + 1.0) * dx;
    pD[1] = (j + 1.0) * dy;
    Line[0] = Normal[0];
    Line[1] = Normal[1];
    Line[2] = Normal[2] * k * dz + Normal[3];

    break;
  }

  return AreaBelow(pA,pB,pC,pD,Line);

}
float VolumeWrapper(int i,int j,int k, int rank, float dx ,float dy, float dz,float Normal[]){


  // Currently we restrict our attention to rotations about a single axis. This is due to a lack of 
  // general formulation of the volume of a cube below a plane. Currently, the zoology of shapes
  // leads to some non-trivial volumes to integrate, and I just haven't done it yet.
  // SO one of the normal components has to be zero.
  // This leaves the volume fraction below the plane the same as the area fraction on the 
  // interesting side.

  
  int LocalDim = -1;
  float SideFraction = 0;

  //Select side through which plane will cross.
  if( Normal[0] == 0 ){ LocalDim = 0;}
  else if( Normal[1] == 0 ){ LocalDim = 1;}
  else if( Normal[2] == 0 ){ LocalDim = 2;}
  else{
    fprintf(stderr," Fatal Error: Currently, one component of the normal vector must be zero\n");
    return -1;
  }
  
  //Modify for 2d
  if( rank == 2 ){
    LocalDim = 2;
  }
  if(rank == 1 ){
    fprintf(stderr,"Grid_MHDBlastInitializeGrid: Don't use InitStyle 10 with a 1d problem.\n");
    return -1;
  }
  SideFraction = AreaWrapper(i,j,k,LocalDim, rank, dx,dy,dz,Normal);

  return SideFraction;
}
    
int SetupNormal(float Normal[], float MHDBlastCenter[3], TopGridData & MetaData){

  int dim;

  //Without loss of generality;
  //1.) Normal components must all be positive.
  //2.) Normal[0] must be the longest 
  //4.) Plane must not cut X faces (well, a little generality, but YOU figure out those boundaries.)
  //5.) The plane agree with the shifted periodic paradigm, i.e. N_zones mod N_shift = 0
  //6.) One of the normal components must be zero.
  //    The only reason (that I can think of right now, I'm sure there's more bugs)
  //    is that the volume average for the setup isn't consistent.  
  //7.) Number of ghost zones MUST accomidated shift.  Here we assume 3, if you do something different,
  //    check.

  //1.)
  for( dim=0;dim<3;dim++)
    if( Normal[dim] < 0 ){
      fprintf(stderr, "ERROR: Normal[%d] = %f < 0.  All must be positive. (rotate your problem.)\n"
	      , dim, Normal[dim]);
      return FAIL;
    }

  //2.)
  if( Normal[0] < Normal[1] || Normal[0] < Normal[2] ){
    fprintf(stderr,"ERROR: Normal[0] must be the longest. (just rotate your problem.)\n");
    return FAIL;
  }


  //0.) Compute the offset.
  Normal[3] = -(Normal[0]*(MHDBlastCenter[0]-DomainLeftEdge[0]) + 
			Normal[1]*(MHDBlastCenter[1]-DomainLeftEdge[1]) + 
			Normal[2]*(MHDBlastCenter[2]-DomainLeftEdge[2])); 


  //4.) plane shouldn't cross x face (makes the boundary conditions less general, less sucky.)
  if( ( DomainLeftEdge[0] * Normal[0] + 
      DomainRightEdge[1]*Normal[1] + 
	DomainRightEdge[1]*Normal[2] + Normal[3] >0 ) ||
      ( DomainRightEdge[0] * Normal[0] + 
	DomainLeftEdge[1]*Normal[1] + 
	DomainLeftEdge[1]*Normal[2] + Normal[3] < 0 ) ){
    fprintf(stderr,"ERROR: Bisection plane cuts through X face.  Either increase the number of zones");
    fprintf(stderr," in that direction, or change the angle.\n");
    return FAIL;
  }

  //5.) Right number of zones.
  //<dbg>
  fprintf(stderr,"normal 2: %f %f %f\n", Normal[0], Normal[1], Normal[2]);
  //</dbg>

  for(dim=1;dim<MetaData.TopGridRank;dim++){

    //\Delta X *nx + \Delta Y * ny == 0
    //\Delta X needs to be an integer number of zones, when \Delta Y is the
    //width of the domain.

    float dX = (DomainRightEdge[0] - DomainLeftEdge[0])/MetaData.TopGridDims[0];
    float Xshift = (DomainRightEdge[dim]-DomainLeftEdge[dim])*Normal[dim]/(dX*Normal[0]);
    if( fabs( nint(Xshift) - Xshift) > 1e-10 ){
      fprintf(stderr,"Parameter Error: Line is not periodic.\n");
      fprintf(stderr,"   DomainSize[%d] * Normal[%d]/(CellWidth[0]*Normal[0]) = %f. Must be an integer.\n",
	      dim,dim, Xshift);
	      
      return FAIL;
    }

    
  }//
  
  //6.) All positive.
  if( Normal[0] * Normal[1] * Normal[2] != 0 ) {
    fprintf(stderr,"Error: As of now, one of the components of Normal must be zero.\n");
    fprintf(stderr," (formulation of the totally general area fractions hasn't been finished.)\n");
    return FAIL;
  }
  

  //7) Ghost Zones.
  float CellWidth[3];
  for(dim=0;dim<MetaData.TopGridRank;dim++)
    CellWidth[dim] = (DomainRightEdge[dim] - DomainLeftEdge[dim])/(MetaData.TopGridDims[dim]);

  for( dim=1;dim<MetaData.TopGridRank;dim++){

    if( DEFAULT_GHOST_ZONES <
	3 + (DomainRightEdge[dim] - DomainLeftEdge[dim])*(Normal[dim]/Normal[0] )  ){
      fprintf(stderr,"NormalSetup: Not enough ghost zones.  To accomidate ShiftedPeriodic, and to avoid\n");
      fprintf(stderr,"       any more MPI routines (costly to code) there must be enough ghost zones\n");
      fprintf(stderr,"       to accomidate both the solver and the shift. (also note, if not using\n");
      fprintf(stderr,"       the athena solver, bad things might happen anyways.)\n");
      fprintf(stderr,"       In your case, you need %f\n",
	      3 + (DomainRightEdge[dim] - DomainLeftEdge[dim]) *(Normal[dim]/(Normal[0]*CellWidth[dim]) ));


      return FAIL;
    }
  }

  return SUCCESS;
}
