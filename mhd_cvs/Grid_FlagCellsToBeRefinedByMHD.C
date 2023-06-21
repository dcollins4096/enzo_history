#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

//dcc 12/30/05 Updated this code to be slightly more usable.  Requires #define DC_NEW_FLAG_8
//The new version now takes floating, global coordinates.

//THe old version can be found after the new code.
int grid::FlagCellsToBeRefinedByMHD(){

  //
#define DC_NEW_FLAG_8
#ifdef DC_NEW_FLAG_8

  //An error check.  You never know when these things will be usefull
  //fprintf(stderr,"FlagByMHD, just checking that we're here\n");
  //return -1;

  
  //A more practical error check
  if( FlaggingField == NULL ) {
    fprintf(stderr, "++++++++++ FlaggCellsMHD: No FlaggingField ++++++++++\n");
    return -1;
  }

  //Initialize  some variables.
  //This is the maximum number of allowed subgrids in this routine.
  //Lefts0, Lefts1, and Lefts2 are the left extents of the new refinement region. (same for rights.)
  //Not that this works on strict containments-- it's the edges of the cell that get compared to the
  //edges of the box.
  //nSg is the number of USED subgrids.
  //nStart is position in Lefts(etc) to start (should be less that nSg)  This is used for dynamic
  //       subgrid changes.
  //CW is the cell width, its just short hand.
  //EPS is a small ammount.  Direct comparison of floats is trouble.

#define MaxSubgridsHere 5
  int NumberOfFlaggedCells = 0;
  int size=1,i,j,k,index;
  int nSg, n, nStart = 0; 
  FLOAT CW[3] = {CellWidth[0][0], CellWidth[1][0],CellWidth[2][0]};
  FLOAT EPS[3] = {CW[0]/10.0,CW[1]/10.0,CW[2]/10.0};

  int nb = DEFAULT_GHOST_ZONES;
  //Really, man, this should be a friggin grid variable.
  for(int dim=0;dim<3;dim++){
    size *= GridDimension[dim];
  }

  //For MHD FC fixing, I've been using [2:12,2:12,2:7] or [2:12,2:12,8:12] (cells) on 16^3
  //For 256^3 driving failure and SibSUB, [64:192]^3
  //for FC fast, 8->2, for Sedov, 25
  //for FCfast, 12->7, for sedov 38


  //Note-- when manipulating these values, make sure you know what 'NumberOfBufferZones' is.
  //It will change things.
  FLOAT aaa =  0.601562, bbb = 0.699219;
  nSg = 1;
  nStart = 0;
  //If you wan a specific cell NUMBER, don't forget to add 1 to the Right Edge.

  FLOAT Lefts0[MaxSubgridsHere] ={aaa , (32)*CW[0],0.949219, 0.0, 0.0};
  FLOAT Lefts1[MaxSubgridsHere] ={aaa , (32)*CW[1],0.949219, 0.0, 0.0};
  FLOAT Lefts2[MaxSubgridsHere] ={aaa , (32)*CW[2],0.949219, 0.0, 0.0};
  FLOAT Rights0[MaxSubgridsHere]={bbb , (32)*CW[0], 1.0, 0.0, 0.0};
  FLOAT Rights1[MaxSubgridsHere]={bbb , (32)*CW[1], 1.0, 0.0, 0.0};
  FLOAT Rights2[MaxSubgridsHere]={bbb , (32)*CW[2], 1.0, 0.0, 0.0};
  
  //For time dependant by-hand gridding, move to the second column (or more)
  //in the Left/Right arrays.

  if(0==1)
    if(dccCounter11++>=2 ){
      nStart++;
      //nSg++; not necessary anymore.
    }

  for(i=nStart;i<nStart+nSg; i++){                              
    fprintf( stderr, "FMHD FLAGGING [%f:%f, %f:%f, %f:%f], \nFMHD     This [%f:%f, %f:%f, %f:%f]\n",
	     Lefts0[i],Rights0[i],Lefts1[i],Rights1[i],Lefts2[i],Rights2[i],
	     GridLeftEdge[0],GridRightEdge[0],GridLeftEdge[1],GridRightEdge[1],
	     GridLeftEdge[2],GridRightEdge[2]);
  }

  for(k=GridStartIndex[2];k<=GridEndIndex[2];k++)
    for(j=GridStartIndex[1];j<=GridEndIndex[1];j++)
      for(i=GridStartIndex[0];i<=GridEndIndex[0];i++)
	{
	  index=index0(i,j,k);
	  
	  for(int n=nStart;n<nStart+nSg;n++){
	    
	    if((CellLeftEdge[0][i] >= Lefts0[n]-EPS[0]) && (CellLeftEdge[0][i]+CW[0] <= Rights0[n]+EPS[0]) &&
	       (CellLeftEdge[1][j] >= Lefts1[n]-EPS[1]) && (CellLeftEdge[1][j]+CW[1] <= Rights1[n]+EPS[1]) &&
	       (CellLeftEdge[2][k] >= Lefts2[n]-EPS[2]) && (CellLeftEdge[2][k]+CW[2] <= Rights2[n]+EPS[2]) ){
	      FlaggingField[index0(i,j,k)]++;
	      NumberOfFlaggedCells++;
	      
	    }
	    
	  }//nSg	    
	  
	  if( nSg == -1){
	    return -1;
	    if( BaryonField[0][index] > 0.011 ){
	      
	      FlaggingField[index0(i,j,k)]++;
	      NumberOfFlaggedCells++;
	    }
	  }//nSg == -1
	}//i,j,k

  return NumberOfFlaggedCells;
  
#endif /*DC_NEW_FLAG_8*/

  //
  //The Old Stuff
  //

#ifndef DC_NEW_FLAG_8

  //fprintf(stderr,"shit-- why no FlagByMHD??\n");
  //return -1;
  int NumberOfFlaggedCells = 0;
  int size=1,i,j,k;
  int index;
  if( FlaggingField == NULL ) {
    fprintf(stderr, "++++++++++ FlaggCellsMHD: No FlaggingField ++++++++++\n");
    return -1;
  }
  size = 1;
  
  for(int dim=0;dim<3;dim++){
    size *= GridDimension[dim];
  }

#define MaxSubgridsHere 20

  int nSg = -2, n, nStart = 0; //nSg will be set later.

  int min[MaxSubgridsHere][3];
  int max[MaxSubgridsHere][3];

  int gs[3];//Global start index.  hard wired for 16^3 problems.

  //gs is for "global start."  it's the start index of this grid 
  //with respect to the total active volume.
  //It is used to re-scale the i,j,k indicies so this 
  //by-hand refinement is processor count independant.
  //The -GridStartIndex[i] is for agreement between TotalActive 
  //coordinates and IncludingGhost coordinates.

  for(i=0;i<3;i++){
    gs[i] = GridLeftEdge[i]*16-GridStartIndex[i];
  }
  fprintf(stderr,"doodie: %d %d %d\n",gs[0],gs[1],gs[2]);


  //Initialize
  for(n=0;n<MaxSubgridsHere;n++){
    for(i=0;i<3;i++){
      min[n][i]=GridDimension[i];
      max[n][i]=0;
    }}

  //fake: 0= 3d block expansion for basic prolong test
  //      1= hard corner problem, 2 step
  //      2= hard corner problem, 1 step
  //      3 = just a square
  //      4 = just a different square.   
  //      5 = two squares in the enter
  //Specify all indicies wrt active volume.
  
  nSg = 1;

//number of possible slots in these arrays.  DOES NOT reflect the actual number of SG
#define Npos 5 
  //int aaa=13, bbb=15;//test case 6 (see notes)
  //aaa=5;bbb=7;//test case 4
  //aaa=8;bbb=10;//test case 3
  //aaa=0; bbb=3; //test 5

  //For MHD FC fixing, I've been using [2:12,2:12,2:7] or [2:12,2:12,8:12]
  int Lefts0[Npos] ={85,0,3,2,12};
  int Lefts1[Npos] ={85,4,3,2,2};
  int Lefts2[Npos] ={85,2,0,2,2};//for FC fast, 8->2, for Sedov, 25
  int Rights0[Npos]={170,GridDimension[0],13,11,15};
  int Rights1[Npos]={170,23,13,15,15};
  int Rights2[Npos]={170,23,1,15,15};//for FCfast, 12->7, for sedov 38

  //For time dependant by-hand gridding, move to the second column (or more)
  //in the Left/Right arrays.

  if(0==1)
    if(dccCounter11++>=2 ){
      nStart++;
      nSg++;
    }

  for(n=nStart;n<nSg;n++){
    min[n][0]=Lefts0[n];
    min[n][1]=Lefts1[n];
    min[n][2]=Lefts2[n];
    max[n][0]=Rights0[n];
    max[n][1]=Rights1[n];
    max[n][2]=Rights2[n];
  }
  
  //declare what's going on, for all the world to see.

  for(i=nStart;i<nSg; i++){
    fprintf( stderr, " FLAGGING [%d:%d, %d:%d, %d:%d], ThisGrid [%d:%d, %d:%d, %d:%d]\n",
	     min[i][0], max[i][0], min[i][1], max[i][1], min[i][2], max[i][2], 
	     GridStartIndex[0],GridEndIndex[0],GridStartIndex[1],
	     GridEndIndex[1],GridStartIndex[2],GridEndIndex[2]);
  }

  for(k=GridStartIndex[2];k<=GridEndIndex[2];k++)
    for(j=GridStartIndex[1];j<=GridEndIndex[1];j++)
      for(i=GridStartIndex[0];i<=GridEndIndex[0];i++)
	{
	  index=index0(i,j,k);
	  for(int n=nStart;n<nSg;n++){
	    /*
	    fprintf(stderr,"moo, i %d j %d k %d gs(%d %d %d) n %d (%d %d)\n",
		    i,j,k,gs[0],gs[1],gs[2], n, min[n][0], max[n][0]);
	    */
	    if( i+gs[0] >= min[n][0]
		&& i+gs[0] <= max[n][0]
		)if( j+gs[1] >= min[n][1]
		     && j+gs[1] <= max[n][1]
		     )if( k+gs[2] >= min[n][2]
			  && k+gs[2] <= max[n][2]
			  ){
		       FlaggingField[index]++;
		       NumberOfFlaggedCells++;

		     }//if in extents
	  }//nSg	    
	  
	  if( nSg == -1){

	    if( BaryonField[0][index] > 0.011 ){

		FlaggingField[index0(i,j,k)]++;
		NumberOfFlaggedCells++;
	    }
	  }//nSg == -1
	}//i,j,k
  
  return NumberOfFlaggedCells;
#endif /*DC_NEW_FLAG_8*/
}
