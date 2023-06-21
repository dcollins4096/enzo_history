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

/* This initializes genearlized discontinuities.  Good for blast waves, Kelving Helmholtz, shock tubes.  Probably others.

Code flow:
1.) Declare/Define Defaults for  parameters.
2.) Read parameters from file.
3.) Calculate TotalEnergy from quantity given (options are Total Energy, Gas Energy, Pressure.)
4.) Set up data labels, units.  
5.) Declare Hierarchy Object.
6.) Define linked list
7.) Call initializer on each level.
8.) Project to parent.

*/

//in MHD_ObliqueRoutines.C
void RotateVector( float * Vector, float * Normal);
int SetupNormal(float Normal[], float MHDDiscontCenter[3], TopGridData & MetaData);

int DiscontInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData, ExternalBoundary &Exterior)
{

  fprintf(stderr, "====================================== \n");
  fprintf(stderr, "====================================== \n");
  fprintf(stderr, "========= DiscontInitialize ========= \n");  
  fprintf(stderr, "====================================== \n");
  fprintf(stderr, "======== MyProcessorNumber %"ISYM"  ========= \n",
	  MyProcessorNumber);
  fprintf(stderr, "====================================== \n");
  fprintf(stderr, "====================================== \n");


  //
  //
  // Labels and Units.  (For IO.)
  //  
  int counter = 0;
  char *DensName = "Density";
  char *TEName = "Total_Energy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BXName = "MagneticField_C_1";
  char *BYName = "MagneticField_C_2";
  char *BZName = "MagneticField_C_3";
#ifdef HAOXU
  if(DualEnergyFormalism ){
    char *GEName = "GasEnergy";
    DataLabel[5] = GEName;
    DataUnits[5] = NULL;   
  }
#endif
  

  DataLabel[counter++] = DensName;
  if( EquationOfState == 0 )
    DataLabel[counter++] = TEName;
  DataLabel[counter++] = Vel1Name;
  DataLabel[counter++] = Vel2Name;
  DataLabel[counter++] = Vel3Name;
  
  if( MHD_Used == TRUE ){
    DataLabel[counter++] = BXName;
    DataLabel[counter++] = BYName;
    DataLabel[counter++] = BZName;
  }
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  
#ifdef MHD
  /* Obsolete.
     MHDcLabel[0] = "MagneticField_C_1";
     MHDcLabel[1] = "MagneticField_C_2";
     MHDcLabel[2] = "MagneticField_C_3";
  */
  
  MHDLabel[0] = "MagneticField_F_1";
  MHDLabel[1] = "MagneticField_F_2";
  MHDLabel[2] = "MagneticField_F_3";
  
  MHDeLabel[0] = "ElectricField_1";
  MHDeLabel[1] = "ElectricField_2";
  MHDeLabel[2] = "ElectricField_3";
  
  MHDcUnits[0] = "FourPiGauss";
  MHDcUnits[1] = "FourPiGauss";
  MHDcUnits[2] = "FourPiGauss";
  
  MHDUnits[0] = "FourPiGauss";
  MHDUnits[1] = "FourPiGauss";
  MHDUnits[2] = "FourPiGauss";
  
  MHDeUnits[0] = "FourPiGauss";
  MHDeUnits[1] = "FourPiGauss";
  MHDeUnits[2] = "FourPiGauss";
  
  CurrentLabel[0] = "Current_1";
  CurrentLabel[1] = "Current_2";
  CurrentLabel[2] = "Current_3";
#endif //MHD

#ifdef PPML
  //
  // Units and Labels for PPML are now set in the Grid Initialization routine PPML_InitializeTypesAndLabels.
  //
#endif //PPML
  //
  // General controll variable
  //

  int dim;

  //
  // Parameters and their defaults.
  // 


  char line[MAX_LINE_LENGTH];
  int ret = 0, GasFlag = 0, Pflag=0, TotalFlag=0;
  int ObsFlag = 0;
  int RefineOnStartup;
  // Or 1.0
  float fpi = 1.0;//1/(sqrt(4.0*3.14159265)*4.0*3.14159265);
  float DensityA = 1.0666,
    DensityB = 1.0,
    GasEnergyA = 3.666,
    GasEnergyB = 1000.666,
    TotalEnergyA = 1.0,
    TotalEnergyB = 1.0;


  float PressureA, PressureB;
  float VelocityA[3] = {0.666, 0.666, 0.666};
  float VelocityB[3] = {0.3666, 0.3666, 0.3666};
  float BA[3]  = {0.5, 0.0, 0.0};
  float BB[3]  = {0.5, 0.0, 0.0};
  float EnergyA, EnergyB;
  float Radius = 4.0;
  int InitStyle = 0, PerturbMethod = -1;
  float DiscontCenter[3] = {0.5,0.5,0.5};
  float PerturbAmplitude = 0.0;
  float PerturbWavelength[3] = {0.0, 0.0, 0.0};
  //now this is a global value, because it's needed by the periodic boundary conditions.
  //float DiscontNormal[4] = {0.0, 0.0, 0.0, 0.0}; //3 cmpts of normal for rotation + offset.
  FLOAT DiscontSubgridLeft[3]  = {DomainLeftEdge[0] ,DomainLeftEdge[1] ,DomainLeftEdge[2]};
  FLOAT DiscontSubgridRight[3] = {DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2]};
  
  //Obsolete variable names.
  float Pressure0, Pressure1;
  float B0[3],B1[3],Energy0, Energy1;
  float Density0,Density1, GasEnergy0, GasEnergy1, TotalEnergy0,TotalEnergy1;
    
  //
  // Set up dummy RandomForcing info.
  //

  if( RandomForcing == TRUE ){
    if( RandomForcingEdot < 0 ) {
      fprintf(stderr, "DiscontInitialize: RandomForcing must be explicitly defined.\n");
      return FAIL;
    }
    RandomForcingNumberOfFields = ( ( HydroMethod != PPM_Local ) ? MetaData.TopGridRank : MetaData.TopGridRank* 7);
  }

  //
  // Read Parameter File.
  //

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    //I changed some nominclature, to make things easier on myself
    //This checks for old nominclature.
    ObsFlag = 0;
    
    ret += sscanf(line, "DiscontDA = %"PSYM, &DensityA);
    ret += sscanf(line, "DiscontDB = %"PSYM, &DensityB);

    ret += sscanf(line, "DiscontBA = %"PSYM" %"PSYM" %"PSYM, BA, BA+1, BA+2);
    ret += sscanf(line, "DiscontBB = %"PSYM" %"PSYM" %"PSYM, BB, BB+1, BB+2);

    ret += sscanf(line, "DiscontVelocityA = %"PSYM" %"PSYM" %"PSYM,
		  VelocityA, VelocityA+1,VelocityA+2);
    ret += sscanf(line, "DiscontVelocityB = %"PSYM" %"PSYM" %"PSYM,
		  VelocityB, VelocityB+1,VelocityB+2);
    
    Pflag += sscanf(line, "DiscontPA = %"PSYM, &Pressure0);
    Pflag += sscanf(line, "DiscontPB = %"PSYM, &Pressure1);

    GasFlag += sscanf(line, "DiscontGasEnergyA = %"PSYM, &GasEnergyA);
    GasFlag += sscanf(line, "DiscontGasEnergyB = %"PSYM, &GasEnergyB);

    TotalFlag += sscanf(line, "DiscontTotalEnergyA = %"PSYM, &TotalEnergyA);
    TotalFlag += sscanf(line, "DiscontTotalEnergyB = %"PSYM, &TotalEnergyB);
    
    if( sscanf(line, "PerturbAmplitude      = %"PSYM, &PerturbAmplitude) != 0 ||
	sscanf(line, "PerturbMethod     = %"ISYM, &PerturbMethod) != 0           ||
	sscanf(line, "PerturbWavelength = %"PSYM,PerturbWavelength)      != 0 ){
      fprintf(stderr," Parameter renamed: PerturbAmplitude(,Method,Wavelength) -> DiscontPerturbAmplitude(etc).  Fix it.\n");
      return FAIL;
    }

    ////

    ret += sscanf(line, "DiscontRadius = %"PSYM, &Radius);

    ret += sscanf(line, "DiscontInitStyle = %"ISYM, &InitStyle);

    ret += sscanf(line, "DiscontCenter = %"PSYM" %"PSYM" %"PSYM,
		  DiscontCenter, DiscontCenter+1,DiscontCenter+2);

    ret += sscanf(line, "DiscontSubgridLeft  = %"PSYM" %"PSYM" %"PSYM,
		  DiscontSubgridLeft, DiscontSubgridLeft +1 , DiscontSubgridLeft +2);
    ret += sscanf(line, "DiscontSubgridRight = %"PSYM" %"PSYM" %"PSYM,
		  DiscontSubgridRight, DiscontSubgridRight +1 , DiscontSubgridRight +2);

    ret += sscanf(line, "DiscontPerturbAmplitude      = %"PSYM, &PerturbAmplitude);
    ret += sscanf(line, "DiscontPerturbMethod     = %"ISYM, &PerturbMethod);
    ret += sscanf(line, "DiscontPerturbWavelength = %"PSYM" %"PSYM" %"PSYM,
                  PerturbWavelength,PerturbWavelength+1,PerturbWavelength+2);

    ret += sscanf(line, "DiscontRefineOnStartup  = %"ISYM, &RefineOnStartup);

    ret += sscanf(line, "DiscontNormal = %"PSYM" %"PSYM" %"PSYM,
		  DiscontNormal, DiscontNormal+1,DiscontNormal+2);

  }//line loop

  if( RandomForcing == TRUE ) 
    RandomForcingNumberOfFields = ( ( HydroMethod != PPM_Local ) ? MetaData.TopGridRank : MetaData.TopGridRank* 7);

  //
  //Re scale the subgrid edges to line up with the parent grid.
  // nCellsL and nCellsR are the number of cells from the domain left edge.
  //

  int nCellsL[3],nCellsR[3];
  int nCells[3] = {0,0,0};  
  for( dim = 0; dim < 3; dim++){
    nCellsL[dim]= nint(( DiscontSubgridLeft[dim] - DomainLeftEdge[dim] )/
		     (DomainRightEdge[dim]-DomainLeftEdge[dim])*MetaData.TopGridDims[dim]);


    DiscontSubgridLeft[dim]=max( nCellsL[dim]*(DomainRightEdge[dim]-DomainLeftEdge[dim])/MetaData.TopGridDims[dim],
				  DomainLeftEdge[dim]);
    
    nCellsR[dim] = nint(( DiscontSubgridRight[dim] - DomainLeftEdge[dim] )/
			(DomainRightEdge[dim]-DomainLeftEdge[dim])*MetaData.TopGridDims[dim]);
    
    DiscontSubgridRight[dim] = min( nCellsR[dim]*(DomainRightEdge[dim]-DomainLeftEdge[dim])/MetaData.TopGridDims[dim],
				     DomainRightEdge[dim]);
    nCells[dim] =  nint( (DiscontSubgridRight[dim]-DiscontSubgridLeft[dim])/
      (DomainRightEdge[dim]-DomainLeftEdge[dim])*MetaData.TopGridDims[dim] );
    
    
  }

  if( RefineOnStartup == 1 ){
    fprintf(stderr,"Subgrid Left %f %f %f\n", DiscontSubgridLeft[0], DiscontSubgridLeft[1], DiscontSubgridLeft[2]);
    fprintf(stderr,"Subgrid Right %f %f %f\n",DiscontSubgridRight[0],DiscontSubgridRight[1],DiscontSubgridRight[2]);
    fprintf(stderr,"nCells %"ISYM" %"ISYM" %"ISYM"\n", nCells[0], nCells[1], nCells[2]);
  }
  // Long Dimension is used to conver the radius from Physical units to Grid Units;
  // We want the axis, though, so figure out which is the longest edge (in Grid Units) 
  // then figure out which one it is.  A more elegant solution would be welcomed.

  int LongDimension = 0;
  LongDimension = (nCells[0] > nCells[1] ) ? nCells[0] : nCells[1];
  LongDimension = (LongDimension > nCells[2] ) ? LongDimension : nCells[2];
  for( dim=0; dim<3; dim++)
    if( LongDimension == nCells[dim] ){
      LongDimension = dim;
      break;
    }


  //
  // Calculate Total Energy.
  //


  if( Pflag > 0 ) {
    GasEnergyA = Pressure0/(Gamma-1);
    GasEnergyB = Pressure1/(Gamma-1);

  }
  /*
  for(int i=0;i<3;i++){
    BA[i] *= fpi;
    BB[i] *= fpi;
  }
  */
  //The variable stored is Gas+Kinetic+Magnetic Energy.
  if( GasFlag > 0 || Pflag > 0){
  Energy0 = GasEnergyA + 
    0.5*DensityA*(VelocityA[0]*VelocityA[0] + VelocityA[1]*VelocityA[1] + VelocityA[2]*VelocityA[2])
    +0.5*(BA[0]*BA[0]+BA[1]*BA[1]+BA[2]*BA[2]);

  Energy1 = GasEnergyB + 
    0.5*DensityB*(VelocityB[0]*VelocityB[0] + VelocityB[1]*VelocityB[1] + VelocityB[2]*VelocityB[2])
    +0.5*(BB[0]*BB[0]+BB[1]*BB[1]+BB[2]*BB[2]);
  }

  if( TotalFlag > 0){
    Energy0=TotalEnergyA;
    Energy1=TotalEnergyB;
    
  }


  if( InitStyle == 10 ){
    //
    // Set up normal vector.  
    // Normal[0] * x + Normal[1]*y + Normal[2]*z + Normal[3] = 0
    // x,y,z in Code Units from Domain Wall.
    // The plane is defined to go throught the DiscontCenter, so

    //<dbg>
    fprintf(stderr,"normal before %f %f %f\n", DiscontNormal[0], DiscontNormal[1], DiscontNormal[2]);
    //</dbg>


    //in MHD_ObliqueRoutines
    if( SetupNormal(DiscontNormal, DiscontCenter, MetaData) == FAIL ){
      fprintf(stderr,"Failure in DiscontNormalSetup.\n");
      return FAIL;
    }

    //<dbg>
    fprintf(stderr,"normal after %f %f %f %f\n", 
	    DiscontNormal[0], DiscontNormal[1], DiscontNormal[2], DiscontNormal[3]);
    //</dbg>

    //The initial conditions get rotated, assuming they were initially alligned 
    //with the X axis, in the cyclic order (so VelocityA = Vx, Vy, Vz)
    //The Vx axis is now aligned with the normal component.
    //No rotation about the normal is performed.

    RotateVector( VelocityA, DiscontNormal);
    RotateVector( VelocityB, DiscontNormal);
    RotateVector( BA, DiscontNormal);
    RotateVector( BB, DiscontNormal);
    
  }//rotated vector.

  //
  // Initialize the top grid.  Cant' decide if I want a uniform grid here or DiscontInitialize.
  //
  if( TopGrid.GridData->DiscontInitializeGrid(DensityA, DensityB,
					       Energy0,  Energy1,
					       VelocityA, VelocityB,
					       BA, BB, 
					       Radius, DiscontCenter, LongDimension,
					       PerturbAmplitude, PerturbMethod,PerturbWavelength,
					       InitStyle, DiscontNormal) == FAIL )
    {
      fprintf(stderr, "DiscontInitialize:  Error in DiscontInitializeGrid.\n");
      return FAIL;
    }


  //
  // Generate Hierarchy.
  //
  if( RefineOnStartup == 1 ){
    //Create as many subgrids as there are refinement levels 
    //needed to resolve the initial explosion region upon the start-up. 
    
    HierarchyEntry ** Subgrid;
    if (MaximumRefinementLevel > 0) 
      Subgrid   = new HierarchyEntry*[MaximumRefinementLevel];
    
    //
    //Create new HierarchyEntries.  Note that 'lev' loops are only for the SUBGRIDS.
    //
    
    int lev;
    int NumberOfSubgridZones[3], SubgridDims[3];
    
    for (lev = 0; lev < MaximumRefinementLevel; lev++) 
      Subgrid[lev] = new HierarchyEntry;
    
    for (lev = 0; lev < MaximumRefinementLevel; lev++) {
      
      //Calculate number of cells on this level.
      
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	NumberOfSubgridZones[dim] = nCells[dim]*POW(RefineBy, lev + 1);
      
      fprintf(stderr,"uncle Discont:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1, 
	      NumberOfSubgridZones[0]);
      
      if (NumberOfSubgridZones[0] > 0) {
	
	// fill them out 
	
	if (lev == 0)
	  TopGrid.NextGridNextLevel  = Subgrid[0];
	Subgrid[lev]->NextGridThisLevel = NULL;
	if (lev == MaximumRefinementLevel-1)
	  Subgrid[lev]->NextGridNextLevel = NULL;
	else
	  Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];
	if (lev == 0)
	  Subgrid[lev]->ParentGrid        = &TopGrid;
	else
	  Subgrid[lev]->ParentGrid        = Subgrid[lev-1];
	
	//  compute the dimensions and left/right edges for the subgrid 
	
	for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	  SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*DEFAULT_GHOST_ZONES;
	}
	
	// create a new subgrid and initialize it 
	
	Subgrid[lev]->GridData = new grid;
	Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
	Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
					    DiscontSubgridLeft,DiscontSubgridRight, 0);
	
	
	if( Subgrid[lev]->GridData->DiscontInitializeGrid(DensityA, DensityB,
							   Energy0,  Energy1,
							   VelocityA, VelocityB,
							   BA, BB, 
							   Radius, DiscontCenter, LongDimension,
							   PerturbAmplitude, PerturbMethod,PerturbWavelength,
							   InitStyle, DiscontNormal) == FAIL )
	  {
	    fprintf(stderr, "DiscontInitialize: Error in DiscontInitializeGrid.\n");
	    return FAIL;
	  }
	
	
	
      }//NumberOfSubgridZones > 0
      else{
	printf("SedovDiscont: single grid start-up.\n");
      }
      
    }//level
    
    // Make sure each grid has the best data with respect to the finer grids.
    // This projection juggle is to ensure that, regardless of how the hierarchy is evolved, the field gets projected
    // properly here.
    
#ifdef MHD
    int MHD_ProjectEtmp = MHD_ProjectE;
    int MHD_ProjectBtmp = MHD_ProjectB;
    MHD_ProjectE=FALSE;
    MHD_ProjectB=TRUE;
#endif //MHD
    for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
      if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
							      *(Subgrid[lev-1]->GridData))
	  == FAIL) {
	fprintf(stderr, "Error in ProjectSolutionToParentGrid.\n");
	return FAIL;
      }
    
    // set up the root grid 
    
    if (MaximumRefinementLevel > 0) {
      if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	  == FAIL) {
	fprintf(stderr, "Error in ProjectSolutionToParentGrid.\n");
	return FAIL;
      }
    }
    
    
    /* This, I believe, doens't need to be here.
       else
       if( Subgrid[lev]->GridData->DiscontInitializeGrid(DensityA, DensityB,
       Energy0,  Energy1,
       VelocityA, VelocityB,
       BA, B1, 
       Radius, DiscontCenter, LongDimension
       PerturbAmplitude, PerturbMethod,PerturbWavelength,
       InitStyle) == FAIL )
       {
       fprintf(stderr, "Error in DiscontInitializeGrid.\n");
       return FAIL;
       }
    */
    
#ifdef MHD
    // Put the projection options back to the inital.
    MHD_ProjectE = MHD_ProjectEtmp;
    MHD_ProjectB = MHD_ProjectBtmp;
#endif //MHD
    
    /*
      if( TopGrid.GridData->DiscontInitializeGrid(DensityA, DensityB,
      Energy0,  Energy1,
      VelocityA, VelocityB,
      BA, B1, 
      Radius, DiscontCenter, LongDimension
      PerturbAmplitude, PerturbMethod,PerturbWavelength,
      InitStyle) == FAIL )
      {
      fprintf(stderr, "Error in DiscontInitializeGrid.\n");
      return FAIL;
      }
      
    */
  }//RefineOnStartup

  return SUCCESS;
}








