/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRID CLASS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifndef GRID_DEFINED__
#define GRID_DEFINED__


#include <assert.h>
#include "ProtoSubgrid.h"
#include "ListOfParticles.h"
#include "region.h"

#ifdef ANALYSIS_TOOLS
#include "AnalyzeClusters.h"
#endif

//<dcc> try always including this.
//#ifdef JB_OPT_FLUXES_FIX
#include "TopGridData.h"
//#endif

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
#  include "FastSiblingLocator.h"
#endif

struct LevelHierarchyEntry;
struct HierarchyEntry;

//dcc things I like
struct IndexPointerMap;
class GlobalStats;
#include "EnzoArray.h"
class grid
{
 private:
//
//  General grid class data
//
  int GridRank;                        // number of dimensions
  int GridDimension[MAX_DIMENSION];    // total dimensions of all grids
  int GridStartIndex[MAX_DIMENSION];   // starting index of the active region
                                       //   (zero based)
  int GridEndIndex[MAX_DIMENSION];     // stoping index of the active region
public:                                       //   (zero based)
  FLOAT GridLeftEdge[MAX_DIMENSION];   // starting pos (active problem space)
  FLOAT GridRightEdge[MAX_DIMENSION];  // ending pos (active problem space)
private:
  float dtFixed;                       // current (fixed) timestep
  FLOAT Time;                          // current problem time
  FLOAT OldTime;                       // time corresponding to OldBaryonField
  int   SubgridsAreStatic;             // 
//
//  Baryon grid data
//
  int    NumberOfBaryonFields;                        // active baryon fields
  float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    // pointers to arrays
  float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
  float *RandomForcingField[MAX_DIMENSION];           // pointers to arrays //AK
  int    FieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  FLOAT *CellLeftEdge[MAX_DIMENSION];  //Of the TOTAL VOLUME. Set in PrepareGridDerivedQuan.
  FLOAT *CellWidth[MAX_DIMENSION];
  fluxes BoundaryFluxes;

  float  CourantSafetyNumber;                       // Hydro parameter
  int    PPMFlatteningParameter;                    // PPM parameter
  int    PPMDiffusionParameter;                     // PPM parameter
  int    PPMSteepeningParameter;                    // PPM parameter

//
//  Particle data
//
  int    NumberOfParticles;
  FLOAT *ParticlePosition[MAX_DIMENSION];  // pointers to position arrays
  float *ParticleVelocity[MAX_DIMENSION];  // pointers to velocity arrays
  float *ParticleAcceleration[MAX_DIMENSION+1];  // 
  float *ParticleMass;                     // pointer to mass array
  int   *ParticleNumber;                   // unique identifier
  float *ParticleAttribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
//
//  Gravity data
// 
  float *PotentialField;
  float *AccelerationField[MAX_DIMENSION]; // cell cntr acceleration at n+1/2
  float *GravitatingMassField;
  FLOAT  GravitatingMassFieldLeftEdge[MAX_DIMENSION];
  int    GravitatingMassFieldDimension[MAX_DIMENSION];
  FLOAT  GravitatingMassFieldCellSize;     // all dimensions must be the same
  float *GravitatingMassFieldParticles;     // for particles only
  FLOAT  GravitatingMassFieldParticlesLeftEdge[MAX_DIMENSION];
  FLOAT  GravitatingMassFieldParticlesCellSize;
  int    GravitatingMassFieldParticlesDimension[MAX_DIMENSION];
  gravity_boundary_type GravityBoundaryType;
  float  PotentialSum;
//
//  Rebuild Hierarchy Temporaries
//
  int *FlaggingField;              // Boolean flagging field (for refinement)
  float *MassFlaggingField;        // Used by mass flagging criterion
//
//  Parallel Information
//
  int ProcessorNumber;
//
// Friends
//
  friend int ExternalBoundary::Prepare(grid *TopGrid);
  friend int ProtoSubgrid::CopyFlaggedZonesFromGrid(grid *Grid);

 public:

// -------------------------------------------------------------------------
//  Main hydro/AMR functions
//




/* Grid constructor (Set all data to null/default state). */

   grid();

/* Grid deconstructor (free up memory usage) */

   ~grid();

/* Read grid data from a file (returns: success/failure) */

   inline int ReadGrid(FILE *main_file_pointer);
 private:
  // The following are private since they should only be called by ReadGrid()
   int ReadGridHDF4(FILE *main_file_pointer);
   int ReadGridHDF5(FILE *main_file_pointer);
 public:

/* Get field or particle data based on name or integer 
   defined in typedefs.h. Details are in Grid_CreateFieldArray.C. */

  EnzoArrayInt *CreateFieldArrayInt(field_type_int field);
  EnzoArrayInt *CreateFieldArrayInt(char *field_name);
  
  EnzoArrayFloat *CreateFieldArrayFloat(field_type_int field);
  EnzoArrayFloat *CreateFieldArrayFloat(char *field_name);
  
  EnzoArrayFLOAT *CreateFieldArrayFLOAT(field_type_int field);
  EnzoArrayFLOAT *CreateFieldArrayFLOAT(char *field_name);


/* Write grid data to a file (returns: success/failure) */

   inline int WriteGrid(FILE *main_file_pointer, char *base_name, int grid_id);
  // The following are private since they should only be called by ReadGrid()
   //  (see macros_and_parameters.h for descriptions)
 private:
   int WriteGridFortran(FILE *main_file_pointer, char *base_name, int grid_id);
   int WriteGridRaw(FILE *main_file_pointer, char *base_name, int grid_id);
   int WriteGridHDFDFSD(FILE *main_file_pointer, char *base_name, int grid_id);
   int WriteGridHDFSD(FILE *main_file_pointer, char *base_name, int grid_id);
   int WriteGridSD(FILE *main_file_pointer, char *base_name, int grid_id);
   int WriteGridHDF5(FILE *main_file_pointer, char *base_name, int grid_id);
 public:

   
/* Write grid data to separate files (returns: success/failure) */

   int WriteGridX(FILE *main_file_pointer, char *base_name, int grid_id);

/* Interpolate to specified time and write grid data to a file
   (returns: success/failure). */

   int WriteGridInterpolate(FLOAT WriteTime, FILE *main_file_pointer, 
			    char *base_name, int grid_id);

/* Compute the timestep constraint for this grid
    (for steps #3 and #4) */

//level is here only as a kludge.
   float ComputeTimeStep(int level);

/* Set the timestep in this grid to the timestep in the argument
    (for step #3) */

   void SetTimeStep(float dt) {dtFixed = dt;};

/* Check timestep (dtFixed) against argument (return fail if dtFixed > dt).
    (for step #4) */

   int CheckTimeStep(float dt) {return ((dtFixed > dt) ? FAIL : SUCCESS);};

/* Return time, timestep */

   FLOAT ReturnTime() {return Time;};
   float ReturnTimeStep() {return dtFixed;};

/* Baryons: Interpolate (parental) grid in argument to current grid.
            (returns success or fail).
    (for step #16) */

   int InterpolateBoundaryFromParent(grid *ParentGrid, LevelHierarchyEntry * Level);

/* Baryons: Copy current solution to Old solution (returns success/fail)
    (for step #16) */

   int CopyBaryonFieldToOldBaryonField();

/* Baryons: Update boundary according to the external boundary values
    (for step #16) */

   int SetExternalBoundaryValues(ExternalBoundary *Exterior);

/* Baryons: solve hydro equations in this grid (returns: the fluxes of the
           subgrids in the argument).  Returns SUCCESS or FAIL.
    (for step #16) */

   int SolveHydroEquations(int CycleNumber, int NumberOfSubgrids, 
			   fluxes *SubgridFluxes[], int level, int grid);

/* Baryons: return pointer to the BoundaryFluxes of this grid */

   int ReturnFluxDims(fluxes &f, int RefinementFactors[]);

/* Baryons: clear the accumulated boundary fluxes for this grid.
    (for step #16) */

   void ClearBoundaryFluxes();

/* Baryons: projected solution in current grid to the grid in the 
           argument which must have a lower resolution (i.e. downsample 
           the current grid to the appropriate level).
    (for step #18) */

   //The extraction flag causes the magnetic and electric fields to be projected by this routine
   //ONLY for grid extractions.  See the documentation for why this is done.
   int ProjectSolutionToParentGrid(grid &ParentGrid, int MHD_Extraction=FALSE);

/* Baryons: return boundary fluxes from this grid.  Downsample them to
           the refinement factors specified in the argument.
	   Returns FAIL or SUCCESS.
    (for step #19) */

   int GetProjectedBoundaryFluxes(grid *ParentGrid, fluxes &ProjectedFluxes);

/* Return the refinement factors as compared to the grid in the argument
   (integer version) (for step #19) */

   void ComputeRefinementFactors(grid *SubGrid, int RefinementFactors[]) {
     int dim;
     for (dim = 0; dim < GridRank; dim++) RefinementFactors[dim] = 
	 int( CellWidth[dim][0] / SubGrid->CellWidth[dim][0] + 0.5);
     for (dim = GridRank; dim < MAX_DIMENSION; dim++)
       RefinementFactors[dim] = 1;
   };

/* Return the refinement factors as compared to the grid in the argument
   (float version) (for step #19) */

   void ComputeRefinementFactorsFloat(grid *SubGrid, float Factors[]) {
     int dim;
     for (dim = 0; dim < GridRank; dim++) Factors[dim] = 
       (*CellWidth[dim]) / (*(SubGrid->CellWidth[dim]));;
     for (dim = GridRank; dim < MAX_DIMENSION; dim++)
       Factors[dim] = 1.0;
   };

/* Baryons: Search for redundant overlap between two sets of fluxes (other
            and refined).  If found, set the refined fluxes equal to the
	    initial fluxes so there will be no double corrections.(step #19) */

   void CorrectRedundantFluxes(fluxes *OtherFluxes, fluxes *InitialFluxes, 
                               fluxes *RefinedFluxes);

/* Baryons: correct for better flux estimates produced by subgrids
           (i.e given the initial flux estimates and the subgrid flux 
	   estimates, correct the grid to account for the subgrid 
	   flux estimates).  Returns SUCCESS or FAIL.
    (for step #19) */

   int CorrectForRefinedFluxes(fluxes *InitialFluxes, fluxes *RefinedFluxes,
			       fluxes *BoundaryFluxesThisTimeStep
#ifdef JB_OPT_FLUXES_FIX
			       , int SUBlingGrid,
			       TopGridData *MetaData
#endif
);

/* Baryons: add the fluxes pointed to by the argument to the boundary fluxes
            of this grid (sort of for step #16).  Note that the two fluxes
	    must have the same size. */

   int AddToBoundaryFluxes(fluxes *BoundaryFluxesToBeAdded);

/* set new time (time += dt)
    (step #21) */

   void SetTimeNextTimestep() {Time += dtFixed;};

/* set time of this grid (used in setup) */

   void SetTime(FLOAT NewTime) {Time = NewTime;};

/* set hydro parameters (used in setup) */

   void SetHydroParameters(float co, int p1, int p2, int p3) 
     {
       CourantSafetyNumber    = co;
       PPMFlatteningParameter = p1;
       PPMDiffusionParameter  = p2;
       PPMSteepeningParameter = p3;

     }

/* Baryons: compute the pressure at the requested time. */

   int ComputePressure(FLOAT time, float *pressure);

/* Baryons: compute the pressure at the requested time using the dual energy
            formalism. */

   int ComputePressureDualEnergyFormalism(FLOAT time, float *pressure);

/* Baryons: compute the temperature. */

   int ComputeTemperatureField(float *temperature);

/* Baryons: compute X-ray emissivity in specified band. */

   int ComputeXrayEmissivity(float *temperature,
			     float *xray_emissivity, float keV1, float keV2,
			     char *XrayFileName);

/* Baryons: compute number density of ionized elements (just O7 and O8). */

   int ComputeElementalDensity(float *temperature, float *elemental_density,
			       int Type);

/* Baryons: compute the ratio of specific heats. */

   int ComputeGammaField(float *GammaField);

/* Baryons: compute the cooling time. */

   int ComputeCoolingTime(float *cooling_time);

/* Baryons & DualEnergyFormalism: Restore consistency between total and
                                  internal energy fields. */

   int RestoreEnergyConsistency(int Region);

/* Returns some grid info. */

   int ReturnGridInfo(int *Rank, int Dims[], FLOAT Left[], FLOAT Right[]);

   void ReturnFieldType(int OutType[]){
     for(int i=0;i<MAX_NUMBER_OF_BARYON_FIELDS; i++) OutType[i] = FieldType[i];}
/* Subtracts kinetic component from total energy. */

   int ConvertTotalEnergyToGasEnergy();

/* Sets the energy to provide Jean's level support (Zeus: returns coeff). */
   
   int SetMinimumSupport(float &MinimumSupportEnergyCoefficient);

/* Debugging support. */

   int DebugCheck(char *message = "Debug");

// -------------------------------------------------------------------------
// Functions used for analysis
//

/* Calculate the angular momentum of a grid (given center). */

   int CalculateAngularMomentum(FLOAT Center[], float AngularMomentum[],
				float MeanVelocity[], float DMVelocity[],
				FLOAT CenterOfMass[], FLOAT DMCofM[]);

/* Find and track density peaks. */

   int AnalyzeTrackPeaks(int level, int ReportLevel);

/* Project some of the fields to a plane. */

   int ProjectToPlane(FLOAT ProjectedFieldLeftEdge[], 
		      FLOAT ProjectedFieldRightEdge[],
		      int ProjectedFieldDims[], float *ProjectedField[], 
		      int ProjectionDimension, int ProjectionSmooth,
                      int NumberOfProjectedFields, int level,
		      int XrayUseLookupTable, float XrayLowerCutoffkeV,
		      float XrayUpperCutoffkeV, char *XrayFileName);

/* Set the fields to zero under the active region of the specified subgrid. */

   int ZeroSolutionUnderSubgrid(grid *Subgrid, int FieldsToZero, 
                                float Value = 1.0);

/* Convert the grid data to particle data for output. */

   int OutputAsParticleData(FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[],
                           ListOfParticles *ParticleList[2], float BaseRadius);

/* Output star particles to a binary file */

   int OutputStarParticleInformation(FILE *StarFile);

/* Return some information about the grid. */

   int CollectGridInformation(int &GridMemory, float &GridVolume, 
                              int &CellsTotal, float &AxialRatio,
			      int &CellsActive);

/* Output grid information (for movie generation). */

   int OutputGridMovieData(FILE *Gridfptr, FILE *DMfptr, FILE *Starfptr,
			   FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[],
			   FLOAT WriteOutTime, int NumberOfPoints[3],
			   int NumberOfValuesPerPoint[3],
			   char *PointValueNames[3][20], float BaseRadius);

// -------------------------------------------------------------------------
// Functions for radiative cooling and multi-species rate equations
//

/* Solve the radiative cooling/heating equations  */

   int SolveRadiativeCooling();

/* Solve the rate equations. */

   int SolveRateEquations();

/* Compute densities of various species for RadiationFieldUpdate. */

   int RadiationComputeDensities(int level);

// -------------------------------------------------------------------------
// Functions for grid (re)generation.
//

/* Remove un-needed arrays before rebuilding. */

   void CleanUp();

/* Delete all the fields, but leave other grid data. */

   void DeleteAllFields();

/* Clear mass flagging field (gg #1) */

   void ClearMassFlaggingField();

/* Clear boolean flagging field (gg #0) */

   void ClearFlaggingField();

/* Set boolean flagging field */

   int SetFlaggingField(int &NumberOfFlaggedCells, int level);

/* Set flagging field from static regions */

   int SetFlaggingFieldStaticRegions(int level, int &NumberOfFlaggedCells);

/* Delete flagging field */

   void DeleteFlaggingField();

/* Particles: deposit particles living in this grid into the Mass Flagging
             field (gg #2) */

   void DepositParticlesToMassFlaggingField() {};

/* baryons: add baryon density to mass flaggin field (so the mass flagging
            field contains the mass in the cell (not the density) 
            (gg #3) */

   int AddFieldMassToMassFlaggingField();

/* Flag all points that require refining  (and delete Mass Flagging Field).
     Returns the number of flagged cells.  Returns the number of flagged cells
     (gg #4) */

   int FlagCellsToBeRefinedByMass(int level, int method);

/* Flag all points that require refining by their slope.
     Returns the number of flagged cells.  Returns the number of flagged cells
     (gg #4) */

   int FlagCellsToBeRefinedBySlope();

/* Flag all points that require refinging by the presence of shocks.
     Returns the number of flagged cells.  Returns the number of flagged cells
     (gg #4) */

   int FlagCellsToBeRefinedByShocks();

/* Flag all points that require refining by the Jean's length criterion. */

   int FlagCellsToBeRefinedByJeansLength();

/* Flag all points that require refining by Shear. */

   int FlagCellsToBeRefinedByShear();

/* Flag all cells for which tcool < dx/sound_speed. */

   int FlagCellsToBeRefinedByCoolingTime();

/* Flagging all cell adjacent to a previous flagged cell.  Also, remove all
   Flagged cells in the boundary zones and within one zone of the boundary. */

   int FlagBufferZones();

/* Identify new subgrids for this grid (and prove Fermat's last theorem too)
   (gg #5) */

   void IdentifyNewSubgrids(GridList &list);

/* Identify new subgrids for this grid (1x1x1 subgrids).
   (gg #5) */

   void IdentifyNewSubgridsSmall(GridList &list);

/* Coalesce neighbouring subgrids */

   void CoalesceSubgrids(GridList &list);

/* Inherit properties (rank, baryon field types, etc.) from ParentGrid
   (gg # 5,6) */

   void InheritProperties(grid *ParentGrid);

/* set the grid dimensions, left, right edges and cell quantities based
   on arguments (gg #5,6) */

   void PrepareGrid(int Rank, int Dimensions[], 
		    FLOAT LeftEdge[], FLOAT RightEdge[], int NumParticles);

/* Allocates space for grids (dims and NumberOfBaryonFields must be set). */

   void AllocateGrids();

/* set the grid derived quantites (CellLeftEdge, CellWidth & BoundaryFluxes) */

   void PrepareGridDerivedQuantities();

/* baryons: interpolate field values from the Parent Grid (gg #6).
            Returns SUCCESS or FAIL. */


   //dcc 01/04/04 Changed this to accept LevelHierarchyEntry
   int InterpolateFieldValues(grid *ParentGrid, LevelHierarchyEntry * Level = NULL);

/* baryons: check for coincident zones between grids & copy if found.
            (correctly includes periodic boundary conditions). */

   int CheckForOverlap(grid *OtherGrid,
		       boundary_type LeftFaceBoundaryCondition[],
		       boundary_type RightFaceBoundaryCondition[],
		       int (grid::*CopyFunction)(grid *OtherGrid,
						 FLOAT EdgeOffset[]));

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
   int CheckForPossibleOverlap(grid *OtherGrid,
			       boundary_type LeftFaceBoundaryCondition[],
			       boundary_type RightFaceBoundaryCondition[]);
   int CheckForPossibleOverlapHelper(grid *OtherGrid,
				     FLOAT EdgeOffset[MAX_DIMENSION]);

   /* These two routines add grids to the chaining mesh used in the
   FastSiblingLocator method and use the chaining mesh to find
   possible siblings. */

   int FastSiblingLocatorAddGrid(ChainingMeshStructure *mesh );
				 
   int FastSiblingLocatorFindSiblings(ChainingMeshStructure *mesh,
				      SiblingGridList *list,
				      boundary_type LeftBoundaryCondition[],
				      boundary_type RightBoundaryCondition[]);

#endif

#ifdef JB_OPT_FLUXES_FIX

   int CheckForSharedFace(grid *OtherGrid,
			       boundary_type LeftFaceBoundaryCondition[],
			       boundary_type RightFaceBoundaryCondition[]);
   int CheckForSharedFaceHelper(grid *OtherGrid,
				     FLOAT EdgeOffset[MAX_DIMENSION]);
#endif
/* baryons: copy coincident zone from the (old) grid in the argument
            (gg #7).  Return SUCCESS or FAIL. */

   int CopyZonesFromGrid(grid *GridOnSameLevel, 
			 FLOAT EdgeOffset[MAX_DIMENSION]);

/* gravity: copy coincident potential field zones from grid in the argument
            (gg #7).  Return SUCCESS or FAIL. */

   int CopyPotentialField(grid *GridOnSameLevel, 
			  FLOAT EdgeOffset[MAX_DIMENSION]);

/* baryons: check for coincident zone from the (old) grid in the argument
            (gg #7).  Return SUCCESS or FAIL. */

   int CopyZonesFromGridCountOnly(grid *GridOnSameLevel, int &Overlap);

/* Returns whether or not the subgrids of this grid are static. */

   int AreSubgridsStatic() {return SubgridsAreStatic;};

/* Check the energy conservation. */

   int ComputeEnergy(float EnergySum[]);

   /* hack: add density squared field to grid (used in ExtractSection). */

   void CreateDensitySquaredField() {
     int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
     BaryonField[NumberOfBaryonFields] = new float[size];
     for (int i = 0; i < size; i++)
       BaryonField[NumberOfBaryonFields][i] = 
	 BaryonField[0][i]*BaryonField[0][i];
     FieldType[NumberOfBaryonFields++] = Density;
   };

#ifdef HAOXU
  float SetESSpeed(); 
#endif

// -------------------------------------------------------------------------
// Functions for use with gravity.
//

/* Set the gravity boundary type of a grid. */

   void SetGravityParameters(gravity_boundary_type Boundary) {
     GravityBoundaryType = Boundary;};
   gravity_boundary_type ReturnGravityBoundaryType() 
     {return GravityBoundaryType;};

/* Gravity: Initialize, the gravitating Mass Field
    (for steps #5, #6). */

   int InitializeGravitatingMassField(int RefinementFactor);

/* Gravity: Initialize, the particle component of the mass field. */

   int InitializeGravitatingMassFieldParticles(int RefinementFactor);

/* Gravity: allocate & clear the GravitatingMassField. */

   int ClearGravitatingMassField();

/* Gravity & baryons: Copy the parent density field to the extra boundary
      region of GravitatingMassField (if any). */

   int CopyParentToGravitatingFieldBoundary(grid *ParentGrid);

/* Gravity & Particles: allocate & clear the GravitatingMassFieldParticles. */

   int ClearGravitatingMassFieldParticles();

/* Baryons: add the baryon mass to the GravitatingMassField. */

   int AddBaryonsToGravitatingMassField();

/* Generic deposit particles/grids to grid (either GravitatingMassField or
   GravitatingMassFieldParticles depending on the value of DepositField). */

   int DepositPositions(FLOAT *Positions[], float *Mass, int Number, 
			int DepositField);

/* deposit particles/grids to grid (if they are on the grid). */

/* int DepositPositionsEdgeOff(float *Positions[], float *Mass, int Number);*/

/* Gravity: Difference potential to get acceleration field. */

   int ComputeAccelerationField(int DifferenceType, int level);

/* Gravity: Interpolate accelerations from other grid. */

   int InterpolateAccelerations(grid *FromGrid);

/* Gravity: Compute particle and grid accelerations. */

   int ComputeAccelerations(int level);

/* Particles: add overlapping ParticleMassField to Target's 
   GravitatingMassField. */

   int CopyOverlappingMassField(grid *TargetGrid, 
				FLOAT EdgeOffset[MAX_DIMENSION]);

/* Gravity: Allocate and make initial guess for PotentialField. */

   int PreparePotentialField(grid *ParentGrid);

/* Gravity: Allocate and make initial guess for PotentialField. */

   //dcc HEY! Here's an example of Default Value setting.
   int SolveForPotential(int &Done, int level, FLOAT PotentialTime = -1);

/* Gravity: Prepare the Greens Function. */

   int PrepareGreensFunction();
   int PreparePeriodicGreensFunction(region *GreensRegion);

/* Gravity: Copy potential/density into/out of FFT regions. */

   int PrepareFFT(region *InitialRegion, int Field, int DomainDim[]);
   int FinishFFT(region *InitialRegion, int Field, int DomainDim[]);

/* Gravity: Set the external acceleration fields. */

   int ComputeAccelerationFieldExternal();

/* Particles + Gravity: Clear ParticleAccleration. */

   int ClearParticleAccelerations();

/* Baryons + Gravity: Interpolate the AccelerationField in FromGrid to
             AccelerationFieldForCells at the GridPositions in this grid. */

   int InterpolateGridPositions(grid *FromGrid);

/* Particles + Gravity: Interpolate the AccelerationField in FromGrid to
             ParticleAcceleration at the ParticlePositions in this grid. */

   int InterpolateParticlePositions(grid *FromGrid, int DifferenceType);

/* Generic routine for interpolating particles/grid. */

   int InterpolatePositions(FLOAT *Positions[], int dim, float *Field, 
			    int Number);

/* Gravity: Delete GravitatingMassField. */

   void DeleteGravitatingMassField() {
     delete GravitatingMassField; 
     GravitatingMassField = NULL;
   };

/* Gravity: Delete AccelerationField. */

   void DeleteAccelerationField() {
     for (int dim = 0; dim < GridRank; dim++) {
       delete AccelerationField[dim];
       AccelerationField[dim] = NULL;
     }
   };

/* Gravity: Add fixed, external acceleration to baryons & particles. */

   int AddExternalAcceleration();

/* Gravity: deposit baryons into target GravitatingMassField. */

   int DepositBaryons(grid *TargetGrid, FLOAT DepositTime);




// -------------------------------------------------------------------------
// Functions for accessing various grid-based information
//
   int GetGridRank() {return GridRank;}
   int GetGridDimension(int Dimension) {return GridDimension[Dimension];}
   int GetGridStartIndex(int Dimension) {return GridStartIndex[Dimension];}
   int GetGridEndIndex(int Dimension) {return GridEndIndex[Dimension];}
   FLOAT GetGridLeftEdge(int Dimension) {return GridLeftEdge[Dimension];}
   FLOAT GetGridRightEdge(int Dimension) {return GridRightEdge[Dimension];}
   int GetGravitatingMassFieldDimension(int Dimension) {
     return GravitatingMassFieldDimension[Dimension];}
   FLOAT GetGravitatingMassFieldLeftEdge(int Dimension) {
     return GravitatingMassFieldLeftEdge[Dimension];}
   FLOAT GetGravitatingMassFieldCellSize() {
     return GravitatingMassFieldCellSize;}
   //int GetGravityGhostZones() {return GravityGhostZones;}

// -------------------------------------------------------------------------
// Functions for use with particles.
//

/* Particles: Deposit particles in the specified field (DepositField) of the
              TargetGrid at the given time. */

   int DepositParticlePositions(grid *TargetGrid, FLOAT DepositTime, 
				int DepositField);

/* Particles: add overlapping ParticleMassField to Target's 
   GravitatingMassField. */

   int AddOverlappingParticleMassField(grid *TargetGrid, 
				       FLOAT EdgeOffset[MAX_DIMENSION]);

/* Particles: Apply particle acceleration to velocity for particles in this 
              grid
    (for step #9) */

   int UpdateParticleVelocity(float TimeStep);

/* Particles: Update particle positions (push)
    (for step #13) */

   int UpdateParticlePosition(float TimeStep);

/* Particles: Move particles from TargetGrid to this grid. */

   int MoveAllParticles(int NumberOfGrids, grid* TargetGrids[]);

/* Particles: Move particles that lie within this grid from the TargetGrid
              to this grid. */

   int MoveSubgridParticles(grid *TargetGrid,
                            int *Counter,
                            int *Number,
                            float *Mass,
                            FLOAT *Position[],
                            float *Velocity[],
                            float *Attribute[]);

/* Particles: same as above, but a version that is much more efficient. */

   int MoveSubgridParticlesFast(int NumberOfSubgrids, grid *ToGrids[],
				int AllLocal);

/* Particles: Clean up moved particles (from MoveSubgridParticles). */

   int CleanUpMovedParticles();

/* Particles: delete accleration fields. */

   void DeleteParticleAcceleration() {
     for (int dim = 0; dim < GridRank+ComputePotential; dim++) {
       delete ParticleAcceleration[dim];
       ParticleAcceleration[dim] = NULL;
     }
   };

/* Particles & Gravity: Delete GravitatingMassField. */

   void DeleteGravitatingMassFieldParticles() {
     delete GravitatingMassFieldParticles; 
     GravitatingMassFieldParticles = NULL;
     GravitatingMassFieldParticlesCellSize = FLOAT_UNDEFINED;
   };

/* Particles: return number of particles. */

   int ReturnNumberOfParticles() {return NumberOfParticles;};

/* Particles: set number of particles. */

   void SetNumberOfParticles(int num) {NumberOfParticles = num;};

/* Particles: delete particle fields and set null. */

   void DeleteParticles() {
     delete ParticleMass;
     delete ParticleNumber;
     ParticleMass = NULL;
     ParticleNumber = NULL;
     for (int dim = 0; dim < GridRank; dim++) {
       delete ParticlePosition[dim];
       delete ParticleVelocity[dim];
       ParticlePosition[dim] = NULL;
       ParticleVelocity[dim] = NULL;
     }
     for (int i = 0; i < NumberOfParticleAttributes; i++) {
       delete ParticleAttribute[i];
       ParticleAttribute[i] = NULL;
     }   
   };

/* Particles: allocate new particle fields. */

   void AllocateNewParticles(int NumberOfNewParticles) {
     ParticleMass = new float[NumberOfNewParticles];
     ParticleNumber = new int[NumberOfNewParticles];
     for (int dim = 0; dim < GridRank; dim++) {
       ParticlePosition[dim] = new FLOAT[NumberOfNewParticles];
       ParticleVelocity[dim] = new float[NumberOfNewParticles];
     }
     for (int i = 0; i < NumberOfParticleAttributes; i++)
       ParticleAttribute[i] = new float[NumberOfNewParticles];
   };

/* Particles: Copy pointers passed into into grid. */

   void SetParticlePointers(float *Mass, int *Number, FLOAT *Position[], 
			    float *Velocity[], float *Attribute[]) {
     assert (ParticleMass==NULL);
    ParticleMass   = Mass;
     assert (ParticleNumber==NULL);
    ParticleNumber = Number;
    for (int dim = 0; dim < GridRank; dim++) {
      assert (ParticlePosition[dim]==NULL);
      ParticlePosition[dim] = Position[dim];
     assert (ParticleVelocity[dim]==NULL);
      ParticleVelocity[dim] = Velocity[dim];
    }
    for (int i = 0; i < NumberOfParticleAttributes; i++) {
     assert (ParticleAttribute[i]==NULL);
      ParticleAttribute[i] = Attribute[i];
    }
   };

/* Particles: Set new star particle index. */

   void SetNewParticleIndex(int &NumberCount, int BaseNumber) {
    for (int n = 0; n < NumberOfParticles; n++)
      if (ParticleNumber[n] == INT_UNDEFINED)
	ParticleNumber[n] = BaseNumber + NumberCount++;
   };

/* Particles: Add given number to particle index. */

   void AddToParticleNumber(int *Count) {
     if (MyProcessorNumber == ProcessorNumber)
       for (int n = 0; n < NumberOfParticles; n++)
	 ParticleNumber[n] += *Count;
     *Count += NumberOfParticles;
   }

/* Particles: sort particle data in ascending order by number (id). */

void SortParticlesByNumber();

// -------------------------------------------------------------------------
// Communication functions
//

/* Set grid's processor number. */

  void SetProcessorNumber(int Proc) {
    ProcessorNumber = Proc;
  };

/* Return grid's processor number. */

  int ReturnProcessorNumber() {
    return ProcessorNumber;
  }

/* Send a region from a real grid to a 'fake' grid on another processor. */

  int CommunicationSendRegion(grid *ToGrid, int ToProcessor, int SendField, 
			     int NewOrOld, int RegionStart[], int RegionDim[]);

/* Send a region from a 'fake' grid to a real grid on another processor. */

  int CommunicationReceiveRegion(grid *ToGrid, int ToProcessor, 
				 int SendField, int NewOrOld, 
				 int RegionStart[], int RegionDim[],
				 int IncludeBoundary);

/* Move a grid from one processor to another. */

  int CommunicationMoveGrid(int ToProcessor);

/* Send particles from one grid to another. */

  int CommunicationSendParticles(grid *ToGrid, int ToProcessor, 
				int FromStart, int FromNumber, int ToStart);

/* Transfer particle amount level 0 grids. */

  int CommunicationTransferParticles(grid* Grids[], int NumberOfSubgrids, 
		 int ToGrid[6], int NumberToMove[6], 
		 float_int *ParticleData[6], int CopyDirection);

// -------------------------------------------------------------------------
// Helper functions (should be made private)
//

/* Baryons: find certain commonly used variables from the list of fields. */

  int IdentifyPhysicalQuantities(int &DensNum, int &GENum,   int &Vel1Num, 
				 int &Vel2Num, int &Vel3Num, int &TENum);

/* Identify Multi-speciesi fields. */

  int IdentifySpeciesFields(int &DeNum, int &HINum, int &HIINum, 
			    int &HeINum, int &HeIINum, int &HeIIINum,
			    int &HMNum, int &H2INum, int &H2IINum,
                            int &DINum, int &DIINum, int &HDINum);

  //--------------------------------------------------------------------
  // MHD
  //

 public:




  //dcc This debugging tool simply passes the variables into a fortran array
  //This is useful for debugging using Totalview.  Setting watchpoints or 
  //examining the fluid variables from the C end is a pain.
  int TVtool(char * Label );

  //dcc Need a grid that contains level information, for extracton purposes.
  int CreateLevelNumberField(int level);

  int CheckForNans(char * label);
  int ComputeGlobalStatsGrid(GlobalStats *stat);
  int IdentifyPhysicalQuantities_2( IndexPointerMap &index);

  //Mass Conservation check.  For debuggin purposes-- not parallel, doesn't work in amr.

  float PreviousMass;

  void TotalMass(char * label);

  //Used for boundary condition set of AccelerationField.

  //keeps track of the acceleration boundary trick.
  int AccelerationHack;

  int ActualNumberOfBaryonFields;
  int AttachAcceleration();
  int DetachAcceleration();

  float *ActualBaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    
  float *ActualOldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    

  int    ActualFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  
  //CenteredB is used in the Riemann solver (SolveMHDequations) and the timestep (dtMagnetic)
  //MagneticField is the face centered magnetic field, and is the quantity ultimately updated by the 
  //CT style algorithm.

  float *MagneticField[3]; 
  float *CenteredB[3]; 
  float *ElectricField[3];
  float *AvgElectricField[3];
  float *OldMagneticField[3];
  float *OldElectricField[3];
  float *OldCenteredB[3];
  float *OldAccelerationField[3];


  float dtParent; //used for the Electric Field update.

  //float ShitDummyForMemoryProblems[30][30];

  //MHDParentTemp is used for the interpolation.
  //It's a member of the grid class because it's also used in CopyZonesFromGrid,
  //for the prolongation (see Balsara's AMR MHD paper), which has no knowledge of the parent.


  float *MHDParentTemp[3];
  int MHDParentTempPermanent[3];
  int    MHDRefinementFactors[3];
  FLOAT ParentDx, ParentDy, ParentDz;


  float *DyBx, *DzBx, *DyzBx;
  float *DxBy, *DzBy, *DxzBy;
  float *DxBz, *DyBz, *DxyBz;
  int * DBxFlag, *DByFlag, *DBzFlag;
  int MHD_ProlongAllocate(int * ChildDim);
  int MHD_DCheck(int * ChildDim, char * mess);
  int MHD_ProlongFree();


  //these two routines need to go in favor of the 3rd.
  int MHD_ProlongWrapper(LevelHierarchyEntry * Level);
  int MHD_ProlongFineGrid(grid * OtherGrid, int Step, char * Label);

  int MHD_SendOldFineGrids(LevelHierarchyEntry * OldFineLevel, grid *ParentGrid);

  int MHD_CID(LevelHierarchyEntry * OldFineLevel, int Offset[], int TempDim[], int Refinement[]);
  int MHD_CIDWorker(grid* OtherGrid, FLOAT EdgeOffset[MAX_DIMENSION]);

  int MHD_ProjectFace(grid &ParentGrid,
		  boundary_type LeftFaceBoundaryCondition[],
		      boundary_type RightFaceBoundaryCondition[]);

  int WriteFlaggingField(int cycle, int level, int gridnumber);
  
  int ClearAvgElectricField();
  
  void MHDCleanUpTemp(){
    
    MHD_ProlongFree();

    for(int field=0; field<3; field++) {
      if(MHDParentTemp[field] != NULL ){
	delete [] MHDParentTemp[field];
	MHDParentTemp[field] = NULL;
      }

    }};

  int IsItShit(char * stringn);
					      
  
  float *DivB;
  float *Current[3];


  int MagneticSize[3];
  int ElectricSize[3];

  friend class AthenaObj;


#ifdef BIERMANN
 float *BiermannTerm[3];
 int ComputeBiermannTerms();
 int DeleteBiermannTerms();  
#endif

  //Magnetic dimensions: MagneticDims[field][axis]
  int MagneticDims[3][3], ElectricDims[3][3];
  int MHDStartIndex[3][3], MHDEndIndex[3][3];//For the MagneticField
  int MHDeStartIndex[3][3], MHDeEndIndex[3][3];//Electric Field
  int MHDAdd[3][3]; //How much to add to Barryon Dimensions.

  //ElectricFlag is a field that prevents double counting in the Flux Correction step.
  //The check needs to happen for each electric-magnetic field pair.
  bool * ElectricFlag[3][3];

  inline void NewElectricFlag(){ 
    if( MyProcessorNumber != ProcessorNumber || MHD_Used != TRUE)
      return;
    for(int efield=0; efield<3; efield++)
      for(int bfield=0;bfield<3;bfield++){
	if(efield == bfield) continue; 
	ElectricFlag[efield][bfield] = new bool[ ElectricSize[efield] ];
	
	for( int j=0;j<ElectricSize[efield]; j++)
	  ElectricFlag[efield][bfield][j] = 0;
      }
  }
  inline void DeleteElectricFlag(){
    if( MyProcessorNumber != ProcessorNumber || MHD_Used != TRUE )
      return;
    for(int efield=0; efield<3; efield++)
      for(int bfield=0;bfield<3;bfield++){
	if(efield==bfield) continue;
	delete ElectricFlag[efield][bfield];
      }
  }


  //Calulate CenterdB from MagneticField.
  int CenterMagneticField(int * Start = NULL, int * End = NULL);
  int SolveMHDEquations(int CycleNumber, int NumberOfSubgrids, fluxes *SubgridFluxes[],
			ExternalBoundary *Exterior, int level, int grid);


  int MHD_UpdateMagneticField(int level, LevelHierarchyEntry * Level);

  //
  // New Solver Routines
  //
  int MHD_Curl( int * Start, int * End, int Method);
#ifdef ATHENA
  int ZeroAcceleration(){
    float size = 1;
    int dim, i;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
    for (dim = 0; dim < GridRank; dim++) {
      if( AccelerationField[dim] != NULL ) 
	for (i = 0; i < size; i++)
	  AccelerationField[dim][i] = 0;
    }
    return SUCCESS;
  }//zero

  int MHD_dtVisc(float & dTvisc );
  float MHD_Diffusion
    (float * Fluxes, float * Lhs, float * Rhs, int * index, int dim, int MHD_DiffusionMethodInput, float dT);
  int MHD_SetFlattening();
  int MHD_Flatten(float * Lhs, float * Rhs, int * index, int dim);
  int NewSMHD(int CycleNumber, int NumberOfSubgrids, fluxes *SubgridFluxes[],
	      ExternalBoundary *Exterior, int level, int grid,
	      float *RandomForcingNormalization, float TopGridTimeStep);

  int MHD_Athena(int CycleNumber, int level, int grid,
		 int NumberOfSubgrids, fluxes *SubgridFluxes[],
		 float * RandomForcingNormalization, float TopGridTimeStep);

  int MHD_AthenaSetup();
  int MHD_Fluxes(float * Fluxes[], int MHD_ReconInput, int MHD_RiemannInput, int MHD_DiffusionMethodInput, float dT );
  int MHD_FluxDifference(float dT, float * Fluxes[], float RandomForcingNormalization);
  int AllocateFluxes( int NumberOfSubgrids, fluxes *SubgridFluxes[] );
  int MHD_PLM(float * Lhs, float * Rhs, int * index, int dim);
  int MHD_AthenaElectric( float dT, float * Fluxes[]);

  //For reconstruction methods.  Might actually store Total Pressure, might be enthalpy.
  float * MHD_Pressure;

  // the solution sensitive switch for reconstruction.  Define in Grid_MHD_Fluxes.C
  inline void MHD_SSSr(int dim, int* index, int MHD_ReconInput, int MHD_RiemannInput);

  //For filling the subgrid fluxes
  int FillFluxes( int NumberOfSubgrids, fluxes *SubgridFluxes[], float * SolverFluxes[], float dT );

#endif //ATHENA

  //
  //
  // Initialization Routines.
  //
  //


  int MHDBlastInitializeGrid(float Density0, float Density1,
			     float Energy0,  float Energy1,
			     float Velocity0[], float Velocity1[],
			     float B0[], float B1[], 
			     float Radius, float MHDBlastCenter[], int LongDimension,
			     float PerturbAmplitude, int PerturbMethod, float PerturbWavelength[],
			     int InitStyle, float Normal[]);

  int MHDTestInitializeGrid(float Density0, float Density1,
                                 float Energy0,  float Energy1,
                                 float Velocity10, float Velocity11,
                                 float Velocity20, float Velocity21,
                                 float Velocity30, float Velocity31,
                                 float B10,float B20,float B30, float Radius);

  int DiskInitializeGrid(float DiskDensity,
			 float DiskPressure,
			 float DiskTheta,
			 float DiskWaveNumber,
			 float DiskRotation,
			 float DiskCompression,
			 float DiskTranslation[],
			 float DiskMagneticField[],
			 float DiskCenter[],
			 float DiskRotationExp,
			 float DiskCompressionExp);

  int MHDShockInitializeGrid(int ShockTubeDirection);
  int MHDOrszagTangInitGrid(float Density,float Pressure, float V0, float B0 );
  int MHDLoopInitGrid(float Density,float Pressure, float Vx, float Vy, float Vz, float B0, float R0,
		      float Center[], int CurrentAxis); 
  int MHD_Diagnose(char * label);
  int ComputeRandomForcingFields(int mode);
//#ifdef HAOXU
  // hx
  int MHDZeldovichPancakeInitializeGrid(int ZeldovichPancakeDirection,
				      float ZeldovichPancakeCentralOffset,
				      float ZeldovichPancakeOmegaBaryonNow,
				      float ZeldovichPancakeOmegaCDMNow,
				      float ZeldovichPancakeCollapseRedshift,
				      float ZeldovichPancakeInitialTemperature,
			              float ZeldovichPancakeMagneticfieldx,
				      float ZeldovichPancakeMagneticfieldy,
				      float ZeldovichPancakeMagneticfieldz);

  int MHDSphericalInfallInitializeGrid(float InitialPerturbation, 
					int UseBaryons, 
					float SphericalInfallOmegaBaryonNow,
					float SphericalInfallOmegaCDMNow,
					int SubgridIsStatic,
					float SphericalInfallMagneticfieldx,
					float SphericalInfallMagneticfieldy,
					float SphericalInfallMagneticfieldz
					  );
 
  int MHDSphericalInfallGetProfile(int level, int ReportLevel);

  int MHDCausticInitializeGrid(float CausticVelocityAmplitude,
			       float CausticInitialPressure,
			       float CausticInitialDensity,
			       float CausticMagneticfieldx,
			       float CausticMagneticfieldy,
			       float CausticMagneticfieldz);

                               


//#endif /* HAOXU */

   #define MAX_SPHERES 10
                                        
  int MHDCollapseTestInitializeGrid(int NumberOfSpheres,
                             FLOAT SphereRadius[MAX_SPHERES],
                             FLOAT SphereCoreRadius[MAX_SPHERES],
			     float CollapseTestBackgroundDensity,
                             float SphereDensity[MAX_SPHERES],
                             float SphereTemperature[MAX_SPHERES],
                             FLOAT SpherePosition[MAX_SPHERES][MAX_DIMENSION],
                             float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
                             int   SphereType[MAX_SPHERES],
                             int   SphereUseParticles,
                             float UniformVelocity[MAX_DIMENSION],
                             int   SphereUseColour,
                             float InitialTemperature, int level,
                             float MHDCollapseTestMagneticFieldx,
                             float MHDCollapseTestMagneticFieldy,    
                             float MHDCollapseTestMagneticFieldz);


  inline int index0(int i, int j, int k)
    {return i+GridDimension[0]*(j+GridDimension[1]*k);}

  inline int indexb1(int i, int j, int k)
    {return i+MagneticDims[0][0]*(j+MagneticDims[0][1]*k);}

  inline int indexb2(int i, int j, int k)
    {return i+MagneticDims[1][0]*(j+MagneticDims[1][1]*k);}

  inline int indexb3(int i, int j, int k)
    {return i+MagneticDims[2][0]*(j+MagneticDims[2][1]*k);}

  inline int indexba(int i, int j, int k, int field)
    {
      if(field == 0)
	return i+MagneticDims[0][0]*(j+MagneticDims[0][1]*k);

      if(field == 1)
	return i+MagneticDims[1][0]*(j+MagneticDims[1][1]*k);
      if(field == 2)
	  return i+MagneticDims[2][0]*(j+MagneticDims[2][1]*k);
      
      return -666;
    }
  
  
  int FlagCellsToBeRefinedByMHD();
  int FlagCellsToBeRefinedByHand(int level);
  int MHDAnis(char * string);
  
  
  
  void MHDMax(){
      float big = -3000, little = 3000;
      
      for(int field=0;field<3;field++){
	  for(int i=MHDStartIndex[field][0];i<=MHDEndIndex[field][0];i++)
	      for(int j=MHDStartIndex[field][1];j<=MHDEndIndex[field][1];j++)
		  for(int k=MHDStartIndex[field][2];k<=MHDEndIndex[field][2];k++){
		      int index=indexba(i,j,k,field);
		      
		      big = MagneticField[field][index] >big? MagneticField[field][index]:big;
		      little = MagneticField[field][index] <little? MagneticField[field][index]:little;
		      
		  }
	  fprintf(stderr, "MAX AND SHIT: %f %f\n",big, little);
      }
  }
  
  
  void ShowOff(char * basename){
      
      char filename[50];
      sprintf(filename,"%s.meta",basename);
      
      FILE * fptr = fopen(filename, "a");
      fprintf(fptr, ",,,,,,,,,,,,,,,,,,,\n");
      fprintf(fptr, "Dims  %d %d %d\n",GridDimension[0],GridDimension[1],GridDimension[2]);
      fprintf(fptr, "LeftEdge  %"PSYM" %"PSYM" %"PSYM"\n",
	      GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2]);
      fprintf(fptr, "RightEdge  %"PSYM" %"PSYM" %"PSYM"\n",
	      GridRightEdge[0],GridRightEdge[1],GridRightEdge[2]);
      fclose(fptr);
      
  }
  
  inline int CheckOutMHDf(){ 
      if( ProcessorNumber == MyProcessorNumber )
	  for( int field=0;field<3;field++) 
	      if( MagneticField[field] == NULL ) 
		  return field; 
	      else return -1;
      else return -2;
      return 0;
  }
  inline int CheckOutMHDc(){ 
      if( ProcessorNumber == MyProcessorNumber )
	  for( int field=0;field<3;field++) 
	      if( CenteredB[field] == NULL ) 
		  return field; 
	else return -1;
    else return -2;
    return 0;
  }

  inline int CheckOutMHDb(){ 
    if( ProcessorNumber == MyProcessorNumber )
      for( int field=0;field<NumberOfBaryonFields;field++) 
	if( BaryonField[field] == NULL ) 
	  return field; 
	else return -1;
    else return -2;
    return 0;
  }
  int TestForAccelerationField();

  friend int MHDExit(grid* TopGrid, TopGridData * MetaData);



// -------------------------------------------------------------------------
// Functions for Specific problems (usually problem generator functions).
//

/* ShockTube problem: initialize grid (returns SUCCESS or FAIL) */

  int ShockTubeInitializeGrid(int Direction, float Boundary, float Density[],
			      float Pressure[], float Velocity[]);

/* Implosion problem (Kritsuk): initialize grid (returns SUCCESS or FAIL) */

  int ImplosionInitializeGrid(float densityInner, float totalEnergyInner,
			      float densityOuter, float totalEnergyOuter);

  /* Initialize a grid for Sedov Explosion */

  int SedovBlastInitializeGrid(float dr,
			       float SedovBlastInnerTotalEnergy,
			       grid *TopGrid);

  int ProtostellarCollapseInitializeGrid(
				     float ProtostellarCollapseCoreDensity,
				     float ProtostellarCollapseCoreEnergy,
				     float ProtostellarCollapseCoreRadius,
				     float ProtostellarCollapseAnguarVelocity);

/* Initialize for a uniform grid (returns SUCCESS or FAIL) */

  int InitializeUniformGrid(float UniformDensity, float UniformTotalEnergy,
			    float UniformGasEnergy, float UniformVelocity[],
			    float *UniformMagneticField = NULL);

#ifdef HAOXU
/* initialize for a grid with spheres wanted (return SUCCESS of FAIL) */ 

#define MAX_SPHERES 10

  int InitializeGridwithSphere( int NumberofSpheres,
                                int SphereType[MAX_SPHERES],
                                float SpherePosition[MAX_SPHERES][MAX_DIMENSION],
                                float CoreDensity[MAX_SPHERES],
                                float SphereRadius[MAX_SPHERES],
                                float SphereCoreRadius[MAX_SPHERES],
                                float Sphere_n[MAX_SPHERES],
                            float UniformDensity,
                            float UniformPressure, float UniformVelocity[MAX_DIMENSION],
                            float *UniformMagneticField = NULL);
#endif

/* Initialize a grid for the Double Mach reflection problem. */

  int DoubleMachInitializeGrid(float d0, float e0, float u0,float v0,float w0);

/* Zeldovich Pancake: initial grid (returns SUCCESS or FAIL). */

  int ZeldovichPancakeInitializeGrid(int   ZeldovichPancakeDirection,
				     float ZeldovichPancakeCentralOffset,
				     float ZeldovichPancakeOmegaBaryonNow,
				     float ZeldovichPancakeOmegaCDMNow,
				     float ZeldovichPancakeCollapseRedshift,
				     float ZeldovichPancakeInitialTemperature);

/* 1D Pressureless Collapse: initialize grid. */

  int PressurelessCollapseInitializeGrid(int PressurelessCollapseDirection,
				    float PressurelessCollapseInitialDensity,
				      int PressurelessCollapseNumberOfCells);

/* Gravity Test: initialize grid. */

  int TestGravityInitializeGrid(float CentralDensity, 
				int NumberOfNewParticles, int UseBaryons);
/* Gravity Test: check results. */

  int TestGravityCheckResults(FILE *fptr, grid *TopGrid);

/* Gravity Test Motion: initialize grid. */

  int TestGravityMotionInitializeGrid(float InitialVelocity);

/* Gravity Test (Sphere): initialize grid. */

  int TestGravitySphereInitializeGrid(float InteriorDensity, 
				      float ExteriorDensity,
				      float SphereRadius, 
				      int SphereType, int UseBaryons,
				      FLOAT SphereCenter[]);

/* Gravity Test (Sphere): check results. */

  int TestGravitySphereCheckResults(FILE *fptr);

/* Spherical Infall Test: initialize grid. */

  int SphericalInfallInitializeGrid(float InitialPerturbation, int UseBaryons,
				    float SphericalInfallOmegaBaryonNow,
				    float SphericalInfallOmegaCDMNow,
				    int SubgridIsStatic);


/* Spherical Infall Test: get the profile from the center. */

  int SphericalInfallGetProfile(int level, int ReportLevel);

/* GravityEquilibriumTest: initialize grid. */

  int GravityEquilibriumTestInitializeGrid(float ScaleHeight);

/* CollapseTest: initialize grid. */

#define MAX_SPHERES 10

  int CollapseTestInitializeGrid(int NumberOfSpheres,
			     FLOAT SphereRadius[MAX_SPHERES],
			     FLOAT SphereCoreRadius[MAX_SPHERES],
			     float CollapseTestBackgroundDensity,
			     float SphereDensity[MAX_SPHERES],
			     float SphereTemperature[MAX_SPHERES],
			     FLOAT SpherePosition[MAX_SPHERES][MAX_DIMENSION],
			     float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
			     int   SphereType[MAX_SPHERES],
			     int   SphereUseParticles, 
			     float UniformVelocity[MAX_DIMENSION],
			     int   SphereUseColour,
			     float InitialTemperature, int level);

/* CosmologySimulation: initialize grid. */

#define COSMOLOGY_INIT_PARAMETERS_DECL \
    float CosmologySimulationOmegaBaryonNow, \
     float CosmologySimulationOmegaCDMNow, \
     float CosmologySimulationInitialTemperature, \
     char *CosmologySimulationDensityName, \
     char *CosmologySimulationTotalEnergyName, \
     char *CosmologySimulationGasEnergyName, \
     char *CosmologySimulationVelocityNames[], \
     char *CosmologySimulationParticlePositionName, \
     char *CosmologySimulationParticleVelocityName, \
     char *CosmologySimulationParticleMassName, \
     int   CosmologySimulationSubgridsAreStatic, \
     int   TotalRefinement, \
     float CosmologySimulationInitialFrctionHII, \
     float CosmologySimulationInitialFrctionHeII, \
     float CosmologySimulationInitialFrctionHeIII, \
     float CosmologySimulationInitialFrctionHM, \
     float CosmologySimulationInitialFrctionH2I, \
     float CosmologySimulationInitialFrctionH2II, \
     int   CosmologySimulationUseMetallicityField, \
     int  &CurrentNumberOfParticles 

#define COSMOLOGY_INIT_PARAMETERS \
     CosmologySimulationOmegaBaryonNow, \
     CosmologySimulationOmegaCDMNow, \
     CosmologySimulationInitialTemperature, \
     CosmologySimulationDensityName, \
     CosmologySimulationTotalEnergyName, \
     CosmologySimulationGasEnergyName, \
     CosmologySimulationVelocityNames, \
     CosmologySimulationParticlePositionName, \
     CosmologySimulationParticleVelocityName, \
     CosmologySimulationParticleMassName, \
     CosmologySimulationSubgridsAreStatic, \
     TotalRefinement, \
     CosmologySimulationInitialFrctionHII, \
     CosmologySimulationInitialFrctionHeII, \
     CosmologySimulationInitialFrctionHeIII, \
     CosmologySimulationInitialFrctionHM, \
     CosmologySimulationInitialFrctionH2I, \
     CosmologySimulationInitialFrctionH2II, \
     CosmologySimulationUseMetallicityField, \
     CurrentNumberOfParticles 

  inline int CosmologySimulationInitializeGrid (COSMOLOGY_INIT_PARAMETERS_DECL);

  // The following are private since they should only be called by
  // CosmologySimulationInitializeGrid()

 private:
  int CosmologySimulationInitializeGridHDF4(COSMOLOGY_INIT_PARAMETERS_DECL);
  int CosmologySimulationInitializeGridHDF5(COSMOLOGY_INIT_PARAMETERS_DECL);
 public:


#ifdef HAOXU  


#define MHDCOSMOLOGY_INIT_PARAMETERS_DECL \
    float CosmologySimulationOmegaBaryonNow, \
     float CosmologySimulationOmegaCDMNow, \
     float CosmologySimulationInitialTemperature, \
     char *CosmologySimulationDensityName, \
     char *CosmologySimulationTotalEnergyName, \
     char *CosmologySimulationGasEnergyName, \
     char *CosmologySimulationVelocityNames[], \
     char *CosmologySimulationParticlePositionName, \
     char *CosmologySimulationParticleVelocityName, \
     char *CosmologySimulationParticleMassName, \
     int   CosmologySimulationSubgridsAreStatic, \
     int   TotalRefinement, \
     float CosmologySimulationInitialFrctionHII, \
     float CosmologySimulationInitialFrctionHeII, \
     float CosmologySimulationInitialFrctionHeIII, \
     float CosmologySimulationInitialFrctionHM, \
     float CosmologySimulationInitialFrctionH2I, \
     float CosmologySimulationInitialFrctionH2II, \
     int   CosmologySimulationUseMetallicityField, \
     int  &CurrentNumberOfParticles, \
     float CosmologySimulationInitialMagneticField[], \
     char *CosmologySimulationMagneticFieldNames[]

  
#define MHDCOSMOLOGY_INIT_PARAMETERS \
     CosmologySimulationOmegaBaryonNow, \
     CosmologySimulationOmegaCDMNow, \
     CosmologySimulationInitialTemperature, \
     CosmologySimulationDensityName, \
     CosmologySimulationTotalEnergyName, \
     CosmologySimulationGasEnergyName, \
     CosmologySimulationVelocityNames, \
     CosmologySimulationParticlePositionName, \
     CosmologySimulationParticleVelocityName, \
     CosmologySimulationParticleMassName, \
     CosmologySimulationSubgridsAreStatic, \
     TotalRefinement, \
     CosmologySimulationInitialFrctionHII, \
     CosmologySimulationInitialFrctionHeII, \
     CosmologySimulationInitialFrctionHeIII, \
     CosmologySimulationInitialFrctionHM, \
     CosmologySimulationInitialFrctionH2I, \
     CosmologySimulationInitialFrctionH2II, \
     CosmologySimulationUseMetallicityField, \
     CurrentNumberOfParticles, \
     CosmologySimulationInitialMagneticField, \
     CosmologySimulationMagneticFieldNames

  inline int MHDCosmologySimulationInitializeGrid(MHDCOSMOLOGY_INIT_PARAMETERS_DECL);

  inline int BiermannBatteryInitializeGrid(MHDCOSMOLOGY_INIT_PARAMETERS_DECL);

  // The following are private since they should only be called by
  // MHDCosmologySimulationInitializeGrid()
                               
 private:
  int MHDCosmologySimulationInitializeGridHDF5(MHDCOSMOLOGY_INIT_PARAMETERS_DECL);

  int BiermannBatteryInitializeGridHDF5(MHDCOSMOLOGY_INIT_PARAMETERS_DECL);
 public:
   

#endif /*HAOXU*/


/* Supernova restart initialize grid. */

  int SupernovaRestartInitialize(float EjectaDensity,
				 float EjectaMetalDensity, float EjectaRadius,
				 float EjectaThermalEnergy, 
				 FLOAT EjectaCenter[3], int ColourField,
				 int *NumberOfCellsSet);

/* Tricks for Random Forcing. */

  int ReturnNumberOfBaryonFields(){return NumberOfBaryonFields;};
  int SetNumberOfBaryonFields(int &number)
    {NumberOfBaryonFields = number; return 0;};
  int AppendForcingToBaryonFields();
  int DetachForcingFromBaryonFields();
  int RemoveForcingFromBaryonFields();
  int AddRandomForcing(float * norm, float * bulkMomentum, float dtTopGrid);
  int PrepareRandomForcingNormalization(float * GlobVal, int GlobNum);
  inline int ReadRandomForcingFields(FILE *main_file_pointer);
 private:
  int ReadRandomForcingFieldsHDF4(FILE *main_file_pointer);
  int ReadRandomForcingFieldsHDF5(FILE *main_file_pointer);
 public:

/* TurbulenceSimulation: initialize grid. */

#define TURBULENCE_INIT_PARAMETERS_DECL \
     float TurbulenceSimulationInitialDensity, \
     float TurbulenceSimulationInitialTemperature, \
     float TurbulenceSimulationInitialMagneticField[], \
     char *TurbulenceSimulationMagneticNames[], \
     char *TurbulenceSimulationDensityName, \
     char *TurbulenceSimulationTotalEnergyName, \
     char *TurbulenceSimulationGasEnergyName, \
     char *TurbulenceSimulationVelocityNames[], \
     char *TurbulenceSimulationRandomForcingNames[], \
     int   TurbulenceSimulationSubgridsAreStatic, \
     int   TotalRefinement


#define TURBULENCE_INIT_PARAMETERS \
     TurbulenceSimulationInitialDensity, \
     TurbulenceSimulationInitialTemperature, \
     TurbulenceSimulationInitialMagneticField, \
     TurbulenceSimulationMagneticNames, \
     TurbulenceSimulationDensityName, \
     TurbulenceSimulationTotalEnergyName, \
     TurbulenceSimulationGasEnergyName, \
     TurbulenceSimulationVelocityNames, \
     TurbulenceSimulationRandomForcingNames, \
     TurbulenceSimulationSubgridsAreStatic, \
     TotalRefinement


inline int TurbulenceSimulationInitializeGrid(TURBULENCE_INIT_PARAMETERS_DECL);

  // The following are private since they should only be called by
  // TurbulenceSimulationInitializeGrid()

 private:
  int TurbulenceSimulationInitializeGridHDF4(TURBULENCE_INIT_PARAMETERS_DECL);
  int TurbulenceSimulationInitializeGridHDF5(TURBULENCE_INIT_PARAMETERS_DECL);
 public:

/* Comoving coordinate expansion terms. */

  int ComovingExpansionTerms();

/* Adjust the gravity source terms for comoving coordinates. */

  int ComovingGravitySourceTerm();

/* Star Particle handler routine. */

  int StarParticleHandler(int level);

/* Apply a time-action to a grid. */

  int ApplyTimeAction(int Type, float Parameter);

/* Include functions for analysis tools */

#ifdef ANALYSIS_TOOLS
#include "Grid_AnalyzeClusters.h"
#endif


// -------------------------------------------------------------------------
// Analysis functions for AnalysisBaseClass and it's derivatives.
//
 public:
/* Flag cells that overlap a subgrid (used for analysis). */
int FlagRefinedCells(grid *Subgrid);

/* Check to see if the grid overlaps a volume. */
/* Doesn't care about periodicity. */
   inline int IsInVolume( FLOAT LeftEdge[], FLOAT RightEdge[] ){
     for( int i = 0; i < GridRank; i++ ){
       if( (GridLeftEdge[i] >= RightEdge[i]) ||
	   (GridRightEdge[i] <= LeftEdge[i]) ){
	 return FALSE;
       }
     }
     return TRUE;
   }

   /* Check to see if a point is in the grid. */
   inline int PointInGrid( float *point ){
     for( int i = 0; i < GridRank; i++ ){
       if( in_range(point[i], GridLeftEdge[i], GridRightEdge[i]) == FALSE )
	 return FALSE;
     }
     return TRUE;
   }

   /* Flags a 3D array where the grid overlaps. */
   void FlagGridArray( HierarchyEntry ***GridArray, int *dx,
		       float *cell_width, HierarchyEntry *my_HE );


   void BinPDF( void *dens_func,void *en_func );
		

  /* Includes for analysis tools */


};

/***********************************************************************/

#include "message.h"

inline int grid::ReadGrid (FILE *main_file_pointer)
{
#if defined (USE_HDF4)
  return ReadGridHDF4 (main_file_pointer);
#elif defined (USE_HDF5)
  return ReadGridHDF5 (main_file_pointer);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif
}

/***********************************************************************/

inline int 
grid::CosmologySimulationInitializeGrid (COSMOLOGY_INIT_PARAMETERS_DECL)
{
#if defined (USE_HDF4)
  return CosmologySimulationInitializeGridHDF4 (COSMOLOGY_INIT_PARAMETERS);
#elif defined (USE_HDF5)
  return CosmologySimulationInitializeGridHDF5 (COSMOLOGY_INIT_PARAMETERS);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif
}

#ifdef HAOXU
 /***********************************************************************/
    
inline int
grid::MHDCosmologySimulationInitializeGrid (MHDCOSMOLOGY_INIT_PARAMETERS_DECL)
{
#if defined (USE_HDF5)
  return MHDCosmologySimulationInitializeGridHDF5 (MHDCOSMOLOGY_INIT_PARAMETERS);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif
}

inline int
grid::BiermannBatteryInitializeGrid (MHDCOSMOLOGY_INIT_PARAMETERS_DECL)
{
#if defined (USE_HDF5)
  return BiermannBatteryInitializeGridHDF5 (MHDCOSMOLOGY_INIT_PARAMETERS);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif
}

#endif

/***********************************************************************/

inline int 
grid::WriteGrid(FILE *main_file_pointer, char *base_name, int grid_id)
{
#if defined (WRITEGRID_HDF_DFSD)
  return WriteGridHDFDFSD(main_file_pointer, base_name, grid_id);
#elif defined (WRITEGRID_HDF_SD)
#  if defined (USE_HDF4)
     return WriteGridHDFSD(main_file_pointer, base_name, grid_id);
#  elif defined (USE_HDF5)
     return WriteGridHDF5(main_file_pointer, base_name, grid_id);
#  else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#  endif
#elif defined (WRITEGRID_RAW)
  return WriteGridRaw(main_file_pointer, base_name, grid_id);
#elif defined (WRITEGRID_FLEXIO)
  WARNING_MESSAGE;
#elif defined (WRITEGRID_FORTRAN)
  return WriteGridFortran(main_file_pointer, base_name, grid_id);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if WRITEGRID_* not set 
#endif
}


/***********************************************************************/
inline int grid::ReadRandomForcingFields (FILE *main_file_pointer)
{
#if defined (USE_HDF4)
  return ReadRandomForcingFieldsHDF4 (main_file_pointer);
#elif defined (USE_HDF5)
  return ReadRandomForcingFieldsHDF5 (main_file_pointer);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif


}



inline int 
grid::TurbulenceSimulationInitializeGrid (TURBULENCE_INIT_PARAMETERS_DECL)
{
#if defined (USE_HDF4)
  return TurbulenceSimulationInitializeGridHDF4 (TURBULENCE_INIT_PARAMETERS);
#elif defined (USE_HDF5)
  return TurbulenceSimulationInitializeGridHDF5 (TURBULENCE_INIT_PARAMETERS);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif
}

#endif // GRID_DEFINED__
