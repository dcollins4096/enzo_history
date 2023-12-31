/***********************************************************************
/
/  GRID CLASS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Many times by AK, DC, RH, JB, DR...
/
/  PURPOSE:
/
************************************************************************/

#ifndef GRID_DEFINED__
#define GRID_DEFINED__

#include "ProtoSubgrid.h"
#include "ListOfParticles.h"
#include "region.h"
#include "FastSiblingLocator.h"
#include "StarParticleData.h"
#include "AMRH5writer.h"
#include "Star.h"
#include "FOF_allvars.h"
#include "MemoryPool.h"

#ifdef FLUX_FIX
#include "TopGridData.h"
#endif

struct HierarchyEntry;

#include "EnzoArray.h"

//#ifdef ANALYSIS_TOOLS
#include "../anyl/AnalyzeClusters.h"
//#endif

#ifdef TRANSFER
#include "PhotonPackage.h"
#include "ListOfPhotonsToMove.h"
#endif /* TRANSFER */

//extern int CommunicationDirection;

//struct ParticleEntry {
//  FLOAT Position[3];
//  float Mass;
//  float Velocity[3];
//  float Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
//  int Number;
// int Type;
//};

extern int CommunicationDirection;
int FindField(int f, int farray[], int n);
struct LevelHierarchyEntry;

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
                                       //   (zero based)
  FLOAT GridLeftEdge[MAX_DIMENSION];   // starting pos (active problem space)
  FLOAT GridRightEdge[MAX_DIMENSION];  // ending pos (active problem space)
  float dtFixed;                       // current (fixed) timestep
  FLOAT Time;                          // current problem time
  FLOAT OldTime;                       // time corresponding to OldBaryonField
  int   SubgridsAreStatic;             // 
  int   ID;                            // Grid ID Number
//
//  Baryon grid data
//
  int    NumberOfBaryonFields;                        // active baryon fields
  float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    // pointers to arrays
  float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
  float *InterpolatedField[MAX_NUMBER_OF_BARYON_FIELDS]; // For RT and movies
  float *RandomForcingField[MAX_DIMENSION];           // pointers to arrays //AK
  int    FieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  FLOAT *CellLeftEdge[MAX_DIMENSION];
  FLOAT *CellWidth[MAX_DIMENSION];
  fluxes *BoundaryFluxes;

  // For restart dumps

  int NumberOfSubgrids;
  fluxes **SubgridFluxStorage;

  // MHD data
  float *divB;
  float *gradPhi[MAX_DIMENSION];

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
  PINT  *ParticleNumber;                   // unique identifier
  int   *ParticleType;                     // type of particle
  float *ParticleAttribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
//
//  Star particle data
//
  int NumberOfStars;
  Star *Stars;
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
//  Top grid parallelism (for implicit solvers)
//
  int ProcLayout[MAX_DIMENSION];
  int ProcLocation[MAX_DIMENSION];
  int ProcNeighbors[MAX_DIMENSION][2];
//
//  Radiation data
//
#ifdef TRANSFER
  float MaxRadiationDt;
#endif
//
//  Rebuild Hierarchy Temporaries
//
  int   *FlaggingField;             // Boolean flagging field (for refinement)
  float *MassFlaggingField;         // Used by mass flagging criteria
  float *ParticleMassFlaggingField; // Used by particle mass flagging criteria
//
//  Parallel Information
//
  int ProcessorNumber;
//
// Movie Data Format
//
  int TimestepsSinceCreation; 	// Not really since creation anymore... 
  				// resets everytime the grid outputs


//
// Friends
//
  friend int ExternalBoundary::Prepare(grid *TopGrid);
  friend int ProtoSubgrid::CopyFlaggedZonesFromGrid(grid *Grid);
  friend class Star;

#ifdef TRANSFER
#include "PhotonGrid_Variables.h"
#endif

 public:

// -------------------------------------------------------------------------
//  Main hydro/AMR functions
//

/* Grid constructor (Set all data to null/default state). */

   grid();

/* Grid deconstructor (free up memory usage) */

   ~grid();

/* Read grid data from a file (returns: success/failure) */

   int ReadGrid(FILE *main_file_pointer, int GridID, 
		int ReadText = TRUE, int ReadData = TRUE);

/* Read grid data from a group file (returns: success/failure) */

#ifndef NEW_GRID_IO
  int Group_ReadGrid(FILE *fptr, int GridID, HDF5_hid_t file_id, 
			 int ReadText, int ReadData);
#else
   int Group_ReadGrid(FILE *main_file_pointer, int GridID, HDF5_hid_t file_id, 
		      int ReadText = TRUE, int ReadData = TRUE, int ReadEverything = FALSE);
#endif

/* Get field or particle data based on name or integer 
   defined in typedefs.h. Details are in Grid_CreateFieldArray.C. */

   EnzoArrayBool *CreateFieldArrayBool(field_type field);
   EnzoArrayBool *CreateFieldArrayBool(char *field_name);

   EnzoArrayInt *CreateFieldArrayInt(field_type field);
   EnzoArrayInt *CreateFieldArrayInt(char *field_name);
  
   EnzoArrayFloat *CreateFieldArrayFloat(field_type field);
   EnzoArrayFloat *CreateFieldArrayFloat(char *field_name);
  
   EnzoArrayFLOAT *CreateFieldArrayFLOAT(field_type field);
   EnzoArrayFLOAT *CreateFieldArrayFLOAT(char *field_name);

   EnzoArrayPINT *CreateFieldArrayPINT(field_type field);
   EnzoArrayPINT *CreateFieldArrayPINT(char *field_name);

/* Write unigrid cubes to a file (returns: success/failure) */

   int WriteCube(char *base_name, int grid_id, int TGdims[]);

/* Write grid data to a file (returns: success/failure) */

   int WriteGrid(FILE *main_file_pointer, char *base_name, int grid_id);

/* Write grid data to a group file (returns: success/failure) */

#ifndef NEW_GRID_IO
   int Group_WriteGrid(FILE *fptr, char *base_name, int grid_id, HDF5_hid_t file_id);
#else
   int Group_WriteGrid(FILE *main_file_pointer, char *base_name, int grid_id, HDF5_hid_t file_id, int WriteEverything = FALSE);
#endif

   int WriteAllFluxes(hid_t grid_node);
   int WriteFluxGroup(hid_t top_group, fluxes *fluxgroup);

   int ReadAllFluxes(hid_t grid_node);
   int ReadFluxGroup(hid_t top_group, fluxes *fluxgroup);

   int FillFluxesFromStorage(int *ThisNumberOfSubgrids,
        fluxes ***fluxgroup) {
        *ThisNumberOfSubgrids = this->NumberOfSubgrids;
        *fluxgroup = this->SubgridFluxStorage;
        this->SubgridFluxStorage = NULL;
        if(ProcessorNumber != MyProcessorNumber) return -1;
        return 0;
    }

/* Write grid data to separate files (returns: success/failure) */

   int WriteGridX(FILE *main_file_pointer, char *base_name, int grid_id);

/* Write task memory map */

   int WriteMemoryMap(FILE *file_pointer, char *base_name, int grid_id);

/* Write grid-to-task map */

   int WriteTaskMap(FILE *file_pointer, char *base_name, int grid_id);

/* Write grid hierarchy only */

   int WriteStuff(FILE *main_file_pointer, char *base_name, int grid_id);

/* Interpolate to specified time and write unigrid cube data to a file
   (returns: success/failure). */

   int WriteCubeInterpolate(FLOAT WriteTime, char *base_name, int grid_id, int TGdims[]);

/* Interpolate to specified time and write grid data to a file
   (returns: success/failure). */

   int WriteGridInterpolate(FLOAT WriteTime, FILE *main_file_pointer, 
			    char *base_name, int grid_id);

/* Interpolate to specified time and write grid data to a group file
   (returns: success/failure). */

   int Group_WriteGridInterpolate(FLOAT WriteTime, FILE *main_file_pointer,
                            char *base_name, int grid_id, HDF5_hid_t file_id);

   int ComputeVectorAnalysisFields(field_type fx, field_type fy, field_type fz,
                                   float* &curl_x, float* &curl_y, float* &curl_z,
                                   float* &div);

private:
   int write_dataset(int ndims, hsize_t *dims, char *name, hid_t group, 
       hid_t data_type, void *data, int active_only = TRUE,
       float *temp=NULL);
   int read_dataset(int ndims, hsize_t *dims, char *name, hid_t group,
       hid_t data_type, void *read_to, int copy_back_active=FALSE,
       float *copy_to=NULL, int *active_dims=NULL);
   int ReadExtraFields(hid_t group_id);
public:

/* Compute the timestep constraint for this grid
    (for steps #3 and #4) */

   float ComputeTimeStep();

/* Set the timestep in this grid to the timestep in the argument
    (for step #3) */

   void SetTimeStep(float dt) {dtFixed = dt;};

/* Check timestep (dtFixed) against argument (return fail if dtFixed > dt).
    (for step #4) */

   int CheckTimeStep(float dt) {return ((dtFixed > dt) ? FAIL : SUCCESS);};

/* Return time, timestep */

   FLOAT ReturnTime() {return Time;};
   float ReturnTimeStep() {return dtFixed;};

  /* Return, set grid ID */

  void SetGridID(int id) { ID = id; };
  int GetGridID(void) { return ID; };
   
  /* Baryons: return field types. */

  int ReturnFieldType(int type[]) 
  {
    for (int i = 0; i < NumberOfBaryonFields; i++) type[i] = FieldType[i];
    return SUCCESS;
  };

/* Baryons: Interpolate (parental) grid in argument to current grid.
            (returns success or fail).
    (for step #16) */

   int InterpolateBoundaryFromParent(grid *ParentGrid);

/* Baryons: Copy current solution to Old solution (returns success/fail)
    (for step #16) */

   int CopyBaryonFieldToOldBaryonField();
   int CopyOldBaryonFieldToBaryonField();


/* Copy potential field to baryon potential for output purposes. */

   int CopyPotentialToBaryonField();

/* Baryons: Update boundary according to the external boundary values
    (for step #16) */

   int SetExternalBoundaryValues(ExternalBoundary *Exterior);

/* Baryons: solve hydro equations in this grid (returns: the fluxes of the
           subgrids in the argument).  Returns SUCCESS or FAIL.
    (for step #16) */

   int SolveHydroEquations(int CycleNumber, int NumberOfSubgrids, 
			   fluxes *SubgridFluxes[], int level);

/* Baryons: return pointer to the BoundaryFluxes of this grid */

   int ReturnFluxDims(fluxes &f, int RefinementFactors[]);

/* Baryons: prepare and clear the accumulated boundary fluxes for this grid.
    (for step #16) */

   void PrepareBoundaryFluxes();
   void ClearBoundaryFluxes();

/* Baryons: projected solution in current grid to the grid in the 
           argument which must have a lower resolution (i.e. downsample 
           the current grid to the appropriate level).
    (for step #18) */

   int ProjectSolutionToParentGrid(grid &ParentGrid);

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

#ifdef FLUX_FIX
   int CorrectForRefinedFluxes(fluxes *InitialFluxes, fluxes *RefinedFluxes,
			       fluxes *BoundaryFluxesThisTimeStep,
			       int SUBlingGrid,
			       TopGridData *MetaData);
#else
   int CorrectForRefinedFluxes(fluxes *InitialFluxes, fluxes *RefinedFluxes,
			       fluxes *BoundaryFluxesThisTimeStep);
#endif

/* Baryons: add the fluxes pointed to by the argument to the boundary fluxes
            of this grid (sort of for step #16).  Note that the two fluxes
	    must have the same size. */

   int AddToBoundaryFluxes(fluxes *BoundaryFluxesToBeAdded);

/* set new time (time += dt)
    (step #21) */

   void SetTimeNextTimestep() {Time += dtFixed;};
   void SetTimePreviousTimestep() {Time -= dtFixed;};

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

/* Subtracts kinetic component from total energy. */

   int ConvertTotalEnergyToGasEnergy();

/* Sets the energy to provide Jean's level support (Zeus: returns coeff). */
   
   int SetMinimumSupport(float &MinimumSupportEnergyCoefficient);

/* Debugging support. */

   int DebugCheck(char *message = "Debug");

#ifdef EMISSIVITY
   /* define function prototype as a grid member function */
   int ClearEmissivity();
   int CheckEmissivity();
#endif

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

   int ProjectToPlane2(FLOAT ProjectedFieldLeftEdge[], 
		       FLOAT ProjectedFieldRightEdge[],
		       int ProjectedFieldDims[], float *ProjectedField[], 
		       int ProjectionDimension, int ProjectionSmooth,
		       int NumberOfProjectedFields, int level,
		       int MetalLinesUseLookupTable, char *MetalLinesFilename);

/* Set the fields to zero under the active region of the specified subgrid. */

   int ZeroSolutionUnderSubgrid(grid *Subgrid, int FieldsToZero, 
                                float Value = 1.0, int AllProcessors = FALSE,
				int IncludeGhostZones = FALSE);

/* Convert the grid data to particle data for output. */

   int OutputAsParticleData(FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[],
                           ListOfParticles *ParticleList[NUM_PARTICLE_TYPES],
                           float BaseRadius);

/* Output star particles to a binary file */

   int OutputStarParticleInformation(FILE *StarFile);

/* Return some information about the grid. */

   int CollectGridInformation(int &GridMemory, float &GridVolume, 
                              int &NumberOfCells, float &AxialRatio,
                              int &CellsTotal, int &Particles);

/* Output grid information (for movie generation). */

   int OutputGridMovieData(FILE *Gridfptr, FILE *DMfptr, FILE *Starfptr,
			   FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[],
			   FLOAT WriteOutTime, int NumberOfPoints[3],
			   int NumberOfValuesPerPoint[3],
			   char *PointValueNames[3][20], float BaseRadius);

/* Output movie data (sequential format) */

   int WriteNewMovieData(FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[], 
			 int RootResolution, FLOAT StopTime, 
			 AMRHDF5Writer &AmiraGrid,
			 int lastMovieStep, int TopGridCycle, 
			 int WriteMe, int MovieTimestepCounter, int open, 
			 FLOAT WriteTime);

   int ReturnMovieTimestep() { return TimestepsSinceCreation; };

/* Output tracer particle information (interpolated from baryon grid). */

   int TracerParticleOutputData(FILE *ptr, FLOAT WriteOutTime);

// -------------------------------------------------------------------------
// Functions for radiative cooling and multi-species rate equations
//

/* Handle the selection of cooling and chemistry modules */

   int MultiSpeciesHandler();

/* Handle the selection of shock finding algorithm */

   int ShocksHandler();

/* Solve the radiative cooling/heating equations  */

   int SolveRadiativeCooling();

/* Solve the rate equations. */

   int SolveRateEquations();

/* Solve the joint rate and radiative cooling/heating equations  */

   int SolveRateAndCoolEquations();

/* Solve the joint rate and radiative cooling/heating equations using MTurk's Solver */

   int SolveHighDensityPrimordialChemistry();

/* Compute densities of various species for RadiationFieldUpdate. */

   int RadiationComputeDensities(int level);

// -------------------------------------------------------------------------
// Functions for grid (re)generation.
//

/* Remove un-needed arrays before rebuilding. */

   void CleanUp();

/* Delete all the fields, but leave other grid data. */

   void DeleteAllFields();

/* Delete all the fields except for the particle data */

   void DeleteAllButParticles();

/* Sum particle mass flagging fields into ProcessorNumber if particles
   aren't local. */

   int SetParticleMassFlaggingField(int StartProc=0, int EndProc=0, int level=-1, 
				    int ParticleMassMethod=-1, int *SendProcs=NULL, 
				    int NumberOfSends=0);
   int CollectParticleMassFlaggingField(void);
   void ClearParticleMassFlaggingField(void);

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

/* Particles: deposit particles to particle mass flagging field. */

   int DepositMustRefineParticles(int pmethod, int level);

/* baryons: add baryon density to mass flaggin field (so the mass flagging
            field contains the mass in the cell (not the density) 
            (gg #3) */

   int AddFieldMassToMassFlaggingField();

/* Flag all points where we are forbidding refinement from a color field */

   int FlagCellsToAvoidRefinement();

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

/* Flag all points that require refining by the Resistive Scale length criterion. 
   abs(B)/abs(curl(B)) should be larger than cell size*/

   int FlagCellsToBeRefinedByResistiveLength();

/* Flag all points that require refining by Shear. */

   int FlagCellsToBeRefinedByShear();

/* Flag all cells for which tcool < dx/sound_speed. */

   int FlagCellsToBeRefinedByCoolingTime();

/* Flag all cells which are near a must-refine particle. */

   int FlagCellsToBeRefinedByMustRefineParticles();

/* Flag all cells which are within a user-specified refinement region. */

   int FlagCellsToBeRefinedByMustRefineRegion(int level);

/* Flag all cells which are above a user-specified metallicity. */

   int FlagCellsToBeRefinedByMetallicity(int level);


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

   // void CoalesceSubgrids(GridList &list);

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

   int InterpolateFieldValues(grid *ParentGrid);

/* Interpolate one radiation field.  Based on InterpolateFieldValues
   but removed all of the conservative stuff. */   

   int InterpolateRadiationFromParent(grid *ParentGrid, int Field);

/* baryons: check for coincident zones between grids & copy if found.
            (correctly includes periodic boundary conditions). */

   int CheckForOverlap(grid *OtherGrid,
		       boundary_type LeftFaceBoundaryCondition[],
		       boundary_type RightFaceBoundaryCondition[],
		       int (grid::*CopyFunction)(grid *OtherGrid,
						 FLOAT EdgeOffset[]));

 



/* baryons: check for subgrids adjacent to external boundary with reflecting BCs. */

   int CheckForExternalReflections(
				   boundary_type LeftFaceBoundaryCondition[],
				   boundary_type RightFaceBoundaryCondition[]);

/* David Collins flux correction - July 2005 */
#ifdef FLUX_FIX
   int CheckForSharedFace(grid *OtherGrid,
			       boundary_type LeftFaceBoundaryCondition[],
			       boundary_type RightFaceBoundaryCondition[]);

   int CheckForSharedFaceHelper(grid *OtherGrid,
				     FLOAT EdgeOffset[MAX_DIMENSION]);
#endif

/* baryons: check for overlap between grids & return TRUE if it exists
            (correctly includes periodic boundary conditions). */

   int CheckForPossibleOverlap(grid *OtherGrid,
                       boundary_type LeftFaceBoundaryCondition[],
                       boundary_type RightFaceBoundaryCondition[]);
   int CheckForPossibleOverlapHelper(grid *OtherGrid,
                                        FLOAT EdgeOffset[MAX_DIMENSION]);

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

/* These two routines add grids to the chaining mesh used in the
   FastSiblingLocator method and use the chaining mesh to find
   possible siblings. */

   int FastSiblingLocatorAddGrid(ChainingMeshStructure *mesh);

   int FastSiblingLocatorFindSiblings(ChainingMeshStructure *mesh,
                          SiblingGridList *list,
                          boundary_type LeftBoundaryCondition[],
                          boundary_type RightBoundaryCondition[]);

   /* hack: add density squared field to grid (used in ExtractSection). */

   void CreateDensitySquaredField() {
     int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
     BaryonField[NumberOfBaryonFields] = new float[size];
     for (int i = 0; i < size; i++)
       BaryonField[NumberOfBaryonFields][i] = 
	 BaryonField[0][i]*BaryonField[0][i];
     FieldType[NumberOfBaryonFields++] = Density;
   };

   void PrintBaryonFieldValues(int field, int index)
     {fprintf(stdout, "Baryonfield[field = %"ISYM"][index = %"ISYM"] = %g\n", 
	      field, index, BaryonField[field][index]);};

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

   int SolveForPotential(int level, FLOAT PotentialTime = -1);

/* Gravity: Prepare the Greens Function. */

   int PrepareGreensFunction();
   int PreparePeriodicGreensFunction(region *GreensRegion);

/* Gravity: Copy potential/density into/out of FFT regions. */

   int PrepareFFT(region *InitialRegion, int Field, int DomainDim[]);
   int FinishFFT(region *InitialRegion, int Field, int DomainDim[]);

/* Gravity: set the potential boundary for isolated BC's */

   int SetIsolatedPotentialBoundary();

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
     delete [] GravitatingMassField; 
     GravitatingMassField = NULL;
   };

/* Gravity: Delete AccelerationField. */

   void DeleteAccelerationField() {
     if (!((SelfGravity || UniformGravity || PointSourceGravity))) return;
     for (int dim = 0; dim < GridRank; dim++) {
       delete [] AccelerationField[dim];
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


#ifdef TRANSFER
// -------------------------------------------------------------------------
// Functions for use with coupled radiation-hydrodynamics solver.
//
   void SetMaxRadiationDt(float MaxRadDt) {MaxRadiationDt = MaxRadDt;}
#endif


// -------------------------------------------------------------------------
// Functions for accessing specific baryon fields 
// (all sources combined in the file Grid_AccessBaryonFields.C)
//
   float* AccessDensity();
   float* AccessTotalEnergy();
   float* AccessGasEnergy();
   float* AccessVelocity1();
   float* AccessVelocity2();
   float* AccessVelocity3();
   float* AccessElectronDensity();
   float* AccessHIDensity();
   float* AccessHIIDensity();
   float* AccessHeIDensity();
   float* AccessHeIIDensity();
   float* AccessHeIIIDensity();
   float* AccessHMDensity();
   float* AccessH2IDensity();
   float* AccessH2IIDensity();
   float* AccessDIDensity();
   float* AccessDIIDensity();
   float* AccessHDIDensity();
   float* AccessSNColour();
   float* AccessMetallicity();
   float* AccessExtraType0();
   float* AccessExtraType1();
   float* AccessKPhHI();
   float* AccessPhotoGamma();
   float* AccessKPhHeI();
   float* AccessGammaHeI();
   float* AccessKPhHeII();
   float* AccessGammaHeII();
   float* AccessKDissH2I();
   float* AccessGravPotential();
   float* AccessAcceleration0();
   float* AccessAcceleration1();
   float* AccessAcceleration2();
   float* AccessRadPressure0();
   float* AccessRadPressure1();
   float* AccessRadPressure2();
   float* AccessEmissivity0();
   float* AccessRadiationFrequency0();
   float* AccessRadiationFrequency1();
   float* AccessRadiationFrequency2();
   float* AccessRadiationFrequency3();
   float* AccessRadiationFrequency4();
   float* AccessRadiationFrequency5();
   float* AccessRadiationFrequency6();
   float* AccessRadiationFrequency7();
   float* AccessRadiationFrequency8();
   float* AccessRadiationFrequency9();


// -------------------------------------------------------------------------
// Functions for accessing top-grid parallelism information
// (note: information only available/valid for this level)
//

/* Processor layout: get and set the number of procs in each 
   dim within the cartesian processor grid
   (1-based, i.e. {1 1 1} defines a single-processor layout) */ 
   int GetProcessorLayout(int Dimension) {return ProcLayout[Dimension];}
   void SetProcessorLayout(int Dimension, int procs) {
     if (Dimension < 0 || Dimension > MAX_DIMENSION)
       fprintf(stderr,"SetProcessorLayout: invalid Dimension.\n");
     else
       if (procs > 0)  ProcLayout[Dimension] = procs; 
       else fprintf(stderr,"SetProcessorLayout: invalid procs value.\n");
   }

/* Processor location: get and set the location of this grid's proc
   within the cartesian processor grid defined in ProcLayout
   (0-based, i.e. {0 0 0} defines the 1st proc in each dimension) */
   int GetProcessorLocation(int Dimension) {return ProcLocation[Dimension];}
   void SetProcessorLocation(int Dimension, int location) {
     if (Dimension < 0 || Dimension > MAX_DIMENSION)
       fprintf(stderr,"SetProcessorLocation: invalid Dimension.\n");
     else
       if (location >= 0)  ProcLocation[Dimension] = location; 
       else fprintf(stderr,"SetProcessorLocation: invalid location.\n");     
   }

/* Processor neighbors: get and set the grid IDs (not MPI process IDs) of this
   grid's neighbors within the cartesian processor grid defined in ProcLayout. 
     Get... returns the {left=0,right=1} neighbor grid ID in a given dim
     Set... provides access to set neighbor information into the grid */
   int GetProcessorNeighbors(int Dimension, int LR) {
     return ProcNeighbors[Dimension][LR];}
   void SetProcessorNeighbors(int Dimension, int LR, int NBid) { 
     if (Dimension < 0 || Dimension > MAX_DIMENSION)
       fprintf(stderr,"SetProcessorNeighbors: invalid Dimension.\n");
     else
       if (LR < 0 || LR > 1) 
	 fprintf(stderr,"SetProcessorNeighbors: invalid neighbor.\n");    
       else
	 if (NBid >= 0)  ProcNeighbors[Dimension][LR] = NBid; 
	 else fprintf(stderr,"SetProcessorNeighbors: invalid grid ID.\n");    
   }


// -------------------------------------------------------------------------
// Functions for use with particles.
//

/* Particles: Deposit particles in the specified field (DepositField) of the
              TargetGrid at the given time. */

   int DepositParticlePositions(grid *TargetGrid, FLOAT DepositTime, 
				int DepositField);

   int DepositParticlePositionsLocal(FLOAT DepositTime, int DepositField);

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

   int UpdateParticlePosition(float TimeStep, int OffProcessorUpdate = FALSE);

/* Particles: Move particles from TargetGrid to this grid. */

   int MoveAllParticles(int NumberOfGrids, grid* TargetGrids[]);

   int MoveAllParticlesOld(int NumberOfGrids, grid* TargetGrids[]);

/* Particles: Move particles that lie within this grid from the TargetGrid
              to this grid. */

//   int MoveSubgridParticles(grid *TargetGrid);


   int MoveSubgridParticles(grid *TargetGrid,
                            int *Counter,
                            PINT *Number,
                            int *Type,
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
     if (!((SelfGravity || UniformGravity || PointSourceGravity))) return;
     for (int dim = 0; dim < GridRank+ComputePotential; dim++) {
       delete [] ParticleAcceleration[dim];
       ParticleAcceleration[dim] = NULL;
     }
   };

/* Particles & Gravity: Delete GravitatingMassField. */

   void DeleteGravitatingMassFieldParticles() {
     delete [] GravitatingMassFieldParticles; 
     GravitatingMassFieldParticles = NULL;
     GravitatingMassFieldParticlesCellSize = FLOAT_UNDEFINED;
   };

/* Particles: return number of particles. */

   int ReturnNumberOfParticles() {return NumberOfParticles;};

   int ReturnNumberOfStarParticles(void);

/* Particles: set number of particles. */

   void SetNumberOfParticles(int num) {NumberOfParticles = num;};

/* Particles: delete particle fields and set null. */

   void DeleteParticles() {
     delete [] ParticleMass;
     delete [] ParticleNumber;
     delete [] ParticleType;
     ParticleMass = NULL;
     ParticleNumber = NULL;
     ParticleType = NULL;
     for (int dim = 0; dim < GridRank; dim++) {
       delete [] ParticlePosition[dim];
       delete [] ParticleVelocity[dim];
       ParticlePosition[dim] = NULL;
       ParticleVelocity[dim] = NULL;
     }
     for (int i = 0; i < NumberOfParticleAttributes; i++) {
       delete [] ParticleAttribute[i];
       ParticleAttribute[i] = NULL;
     }   
   };

/* Particles: allocate new particle fields. */

   void AllocateNewParticles(int NumberOfNewParticles) {
     ParticleMass = new float[NumberOfNewParticles];
     ParticleNumber = new PINT[NumberOfNewParticles];
     ParticleType = new int[NumberOfNewParticles];
     for (int dim = 0; dim < GridRank; dim++) {
       ParticlePosition[dim] = new FLOAT[NumberOfNewParticles];
       ParticleVelocity[dim] = new float[NumberOfNewParticles];
     }
     for (int i = 0; i < NumberOfParticleAttributes; i++)
       ParticleAttribute[i] = new float[NumberOfNewParticles];
   };

/* Particles: Copy pointers passed into into grid. */

   void SetParticlePointers(float *Mass, PINT *Number, int *Type,
                            FLOAT *Position[], 
			    float *Velocity[], float *Attribute[]) {
    ParticleMass   = Mass;
    ParticleNumber = Number;
    ParticleType   = Type;
    for (int dim = 0; dim < GridRank; dim++) {
      ParticlePosition[dim] = Position[dim];
      ParticleVelocity[dim] = Velocity[dim];
    }
    for (int i = 0; i < NumberOfParticleAttributes; i++)
      ParticleAttribute[i] = Attribute[i];
   };

/* Particles: Set new star particle index. */

   void SetNewParticleIndex(int &NumberCount1, PINT &NumberCount2);

/* Particles: Set new star particle index. - Old version */

   void SetNewParticleIndexOld(int &NumberCount, int BaseNumber) {
     for (int n = 0; n < NumberOfParticles; n++) 
      if (ParticleNumber[n] == INT_UNDEFINED)
	ParticleNumber[n] = BaseNumber + NumberCount++;
   };

/* Particles: Add given number to particle index. */

   void AddToParticleNumber(PINT *Count) {
     if (MyProcessorNumber == ProcessorNumber)
       for (int n = 0; n < NumberOfParticles; n++)
	 ParticleNumber[n] += *Count;
     *Count += NumberOfParticles;
   }

  float ReturnTotalSinkMass() {
    float total = 0;
    double dx3 = CellWidth[0][0] * CellWidth[0][0] * CellWidth[0][0];
    if (MyProcessorNumber == ProcessorNumber)
      for (int n = 0; n < NumberOfParticles; n++)
	if (ParticleType[n] == PARTICLE_TYPE_MUST_REFINE)
	  total += ParticleMass[n] * dx3;
    return total;
  };

/* Particles: return particle type (1 - dm, 2 - star).   Note that this
   has now been superceded by a real particle type field. */

   int ReturnParticleType(int index) {
     if (NumberOfParticleAttributes > 0 && StarParticleCreation > 0)
       if (ParticleAttribute[0][index] > 0)
         return PARTICLE_TYPE_STAR;
     return PARTICLE_TYPE_DARK_MATTER;
   }

/* Particles: return particle information in structure array */

   int ReturnParticleEntry(ParticleEntry *ParticleList);

/* Particles: set mass of merged particles to be -1 */

   void RemoveMergedParticles(ParticleEntry *List, const int &Size, int *Flag);

/* Particles: append particles belonging to this grid from a list */

   int AddParticlesFromList(ParticleEntry *List, const int &Size, int *AddedNewParticleNumber);
   int CheckGridBoundaries(FLOAT *Position);

/* Particles: sort particle data in ascending order by number (id) or type. */

void SortParticlesByNumber();
void SortParticlesByType();

int CreateParticleTypeGrouping(hid_t ptype_dset,
                               hid_t ptype_dspace,
                               hid_t parent_group,
                               hid_t file_id);

 int ChangeParticleTypeBeforeSN(int _type, int level, 
				int *ParticleBufferSize=NULL);

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

  int CommunicationMoveGrid(int ToProcessor, int MoveParticles = TRUE);

/* Send particles from one grid to another. */

  int CommunicationSendParticles(grid *ToGrid, int ToProcessor, 
				int FromStart, int FromNumber, int ToStart);

/* Transfer particle amount level 0 grids. */

#ifdef OPTIMIZED_CTP
  int CommunicationTransferParticles(grid* Grids[], int NumberOfGrids,
				     int ThisGridNum, int *&NumberToMove, 
				     int StartIndex, int EndIndex, 
				     particle_data *&List,
				     int *Layout, int *GridMap, 
				     int CopyDirection);
  int CommunicationTransferStars(grid* Grids[], int NumberOfGrids,
				 int ThisGridNum, int *&NumberToMove, 
				 int StartIndex, int EndIndex, 
				 star_data *&List,
				 int *Layout, int *GridMap, 
				 int CopyDirection);
#else
  int CommunicationTransferParticles(grid* Grids[], int NumberOfGrids,
			     int ToGrid[6], int NumberToMove[6],
			     float_int *ParticleData[6], int CopyDirection);
  int CommunicationTransferStars(grid* Grids[], int NumberOfGrids,
			 int ToGrid[6], int NumberToMove[6],
			 StarBuffer *StarData[6], int CopyDirection);
#endif

  int CollectParticles(int GridNum, int* &NumberToMove, 
		       int &StartIndex, int &EndIndex, 
		       particle_data* &List, int CopyDirection);

  int CollectStars(int GridNum, int* &NumberToMove, 
		   int &StartIndex, int &EndIndex, 
		   star_data* &List, int CopyDirection);

  // Only used for static hierarchies
  int MoveSubgridStars(int NumberOfSubgrids, grid* ToGrids[],
		       int AllLocal);

  int TransferSubgridParticles(grid* Subgrids[], int NumberOfSubgrids, 
			       int* &NumberToMove, int StartIndex, 
			       int EndIndex, particle_data* &List, 
			       bool KeepLocal, bool ParticlesAreLocal,
			       int CopyDirection,
			       int IncludeGhostZones = FALSE);

  int TransferSubgridStars(grid* Subgrids[], int NumberOfSubgrids, 
			   int* &NumberToMove, int StartIndex, 
			   int EndIndex, star_data* &List, 
			   bool KeepLocal, bool ParticlesAreLocal,
			   int CopyDirection,
			   int IncludeGhostZones = FALSE);

// -------------------------------------------------------------------------
// Helper functions (should be made private)
//

  /* This is a simple helper function which determines if the method
     should return immediately because of communication-mode reasons.
     Assumes that info is being sent from the "other-grid" processor to
     the "this-grid" processor (i.e. the one that holds this object). */

  int CommunicationMethodShouldExit(grid *OtherGrid) {

    /* Return if neither grid lives on this processor. */
    //    if (NumberOfProcessors == 1) return SUCCESS;

    if (MyProcessorNumber != ProcessorNumber && 
        MyProcessorNumber != OtherGrid->ProcessorNumber)
      return SUCCESS;

    /* If the two grids are on the same processor then return if
       either in post-receive or receive modes to avoid duplicating method
       (i.e. action is only carried out if in send mode (or send-receive)). */

    if (ProcessorNumber == OtherGrid->ProcessorNumber)
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	  CommunicationDirection == COMMUNICATION_RECEIVE)
        return SUCCESS;

    /* If in send mode then exit if the send grid is not on this processor. */

    if (CommunicationDirection == COMMUNICATION_SEND &&
        MyProcessorNumber != OtherGrid->ProcessorNumber)
      return SUCCESS;

    /* If in either receive phase then exit if receive grid is not on this 
       processor. */

    if ((CommunicationDirection == COMMUNICATION_RECEIVE ||
         CommunicationDirection == COMMUNICATION_POST_RECEIVE) &&
        MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    return FAIL; /* i.e. method should not exit immediately. */
  }

/* Baryons: find certain commonly used variables from the list of fields. */

  int IdentifyPhysicalQuantities(int &DensNum, int &GENum,   int &Vel1Num, 
				 int &Vel2Num, int &Vel3Num, int &TENum);

  int IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, 
				 int &Vel2Num, int &Vel3Num, int &TENum,
				 int &B1Num, int &B2Num, int &B3Num);

  int IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, 
				 int &Vel2Num, int &Vel3Num, int &TENum,
				 int &B1Num, int &B2Num, int &B3Num, int &PhiNum);

  /* Identify driving fields */

  int IdentifyDrivingFields(int &Drive1Num, int &Drive2Num, int &Drive3Num);

  /* Identify potential field */

  int IdentifyPotentialField(int &PotenNum, int &Acce1Num, int &Acce2Num, int &Acce3Num);

  /* Identify colour field */

  int IdentifyColourFields(int &SNColourNum, int &MetalNum, int &MBHColourNum,
			   int &Galaxy1ColourNum, int &Galaxy2ColourNum);

  /* Identify Multi-species fields. */

  int IdentifySpeciesFields(int &DeNum, int &HINum, int &HIINum, 
			    int &HeINum, int &HeIINum, int &HeIIINum,
			    int &HMNum, int &H2INum, int &H2IINum,
                            int &DINum, int &DIINum, int &HDINum);

  // Identify Simon Glover Species Fields
  int IdentifyGloverSpeciesFields(int &HIINum,int &HINum,int &H2INum,
				  int &DINum,int &DIINum,int &HDINum,
				  int &HeINum,int &HeIINum,int &HeIIINum,
				  int &CINum,int &CIINum,int &OINum,
				  int &OIINum,int &SiINum,int &SiIINum,
				  int &SiIIINum,int &CHINum,int &CH2INum,
				  int &CH3IINum,int &C2INum,int &COINum,
				  int &HCOIINum,int &OHINum,int &H2OINum,
				  int &O2INum);

/* Zeus Solver. */

  int ZeusSolver(float *gamma, int igamfield, int nhy, 
		 float dx[], float dy[], float dz[], 
		 int gravity, int NumberOfSubgrids, long_int GridGlobalStart[],
		 fluxes *SubgridFluxes[],
		 int NumberOfColours, int colnum[], int bottom,
		 float minsupecoef);

/* PPM Direct Euler Solver. */

int SolvePPM_DE(int CycleNumber, int NumberOfSubgrids, 
		fluxes *SubgridFluxes[], float *CellWidthTemp[], 
		Elong_int GridGlobalStart[], int GravityOn, 
		int NumberOfColours, int colnum[]);

int xEulerSweep(int k, int NumberOfSubgrids, fluxes *SubgridFluxes[], 
		Elong_int GridGlobalStart[], float *CellWidthTemp[], 
		int GravityOn, int NumberOfColours, int colnum[]);

int yEulerSweep(int i, int NumberOfSubgrids, fluxes *SubgridFluxes[], 
		Elong_int GridGlobalStart[], float *CellWidthTemp[], 
		int GravityOn, int NumberOfColours, int colnum[]);

int zEulerSweep(int j, int NumberOfSubgrids, fluxes *SubgridFluxes[], 
		Elong_int GridGlobalStart[], float *CellWidthTemp[], 
		int GravityOn, int NumberOfColours, int colnum[]);

// AccelerationHack

  int AccelerationHack;

#ifdef SAB
  //These should be moved later.
  //Used for boundary condition set of AccelerationField.
  int ActualNumberOfBaryonFields;
  int AttachAcceleration();
  int DetachAcceleration();

  int    ActualFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  float *ActualBaryonField[MAX_NUMBER_OF_BARYON_FIELDS];
  float *ActualOldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS];
  float *OldAccelerationField[3];
#endif

// -------------------------------------------------------------------------
// Functions for Specific problems (usually problem generator functions).
//

/* Generalized Extra Field Grid Initializer (returns NumberOfBaryonFields) */
  int InitializeTestProblemGrid(int field_counter);

/* Protostellar Collapse problem: initialize grid (returns SUCCESS or FAIL) */
  int ProtostellarCollapseInitializeGrid(float CoreDensity,
					 float CoreEnergy,
					 float CoreRadius,
					 float AngularVelocity);

/* HydroShockTubes problems: initialize grid (returns SUCCESS or FAIL) */

  int HydroShockTubesInitializeGrid(float InitialDiscontinuity,
				    float LeftDensity, float RightDensity, 
				    float LeftVelocityX, float RightVelocityX,
				    float LeftVelocityY, float RightVelocityY,
				    float LeftVelocityZ, float RightVelocityZ,
				    float LeftPressure, float RightPressure);

/* Initialize for a uniform grid (returns SUCCESS or FAIL) */

  int InitializeUniformGrid(float UniformDensity, float UniformTotalEnergy,
			    float UniformGasEnergy, float UniformVelocity[], 
			    float UniformBField[]);


/* Initialize a grid for the Double Mach reflection problem. */

  int DoubleMachInitializeGrid(float d0, float e0, float u0,float v0,float w0);


/* Initialize a grid for Implosion test problem */

  int ImplosionInitializeGrid(float ImplosionDiamondDensity,
			      float ImplosionDiamondTotalEnergy);


/* Initialize a grid for Sedov Explosion */

  int SedovBlastInitializeGrid(float SedovBlastInitialRadius,
                               float SedovBlastInnerTotalEnergy);

  int SedovBlastInitializeGrid3D(char * fname);

/* Initialize a grid for RadiatingShock (Sedov+Cooling) Explosion */

  int RadiatingShockInitializeGrid(FLOAT RadiatingShockInitialRadius,
				   float RadiatingShockInnerDensity,
				   float RadiatingShockInnerTotalEnergy,
				   int RadiatingShockUseDensityFluctuations,
				   int RadiatingShockRandomSeed,
				   float RadiatingShockDensityFluctuationLevel,
				   int RadiatingShockInitializeWithKE,
				   int RadiatingShockUseSedovProfile,
				   FLOAT RadiatingShockSedovBlastRadius,
				   float RadiatingShockEnergy,
				   float RadiatingShockPressure,
				   float RadiatingShockKineticEnergyFraction,
				   float RadiatingShockRhoZero,
				   float RadiatingShockVelocityZero,
				   int RadiatingShockRandomSeedInitialize,
				   FLOAT RadiatingShockCenterPosition[MAX_DIMENSION]);

  /* Initialize a grid for a rotating cylinder collapse */
  int RotatingCylinderInitializeGrid(FLOAT RotatingCylinderRadius,
				     FLOAT RotatingCylinderCenterPosition[MAX_DIMENSION],
				     float RotatingCylinderLambda,
				     float RotatingCylinderOverdensity);

  /* Initialize a grid for the KH instability problem. */

  int KHInitializeGrid(float KHInnerDensity,
                       float KHInnerInternalEnergy,
                       float KHOuterInternalEnergy,
                       float KHPerturbationAmplitude,
                       float KHInnerVx, float KHOuterVx);

  /* Initialize a grid and set boundary for the 2D/3D Noh problem. */

  int NohInitializeGrid(float d0, float p0, float u0);
  int ComputeExternalNohBoundary();

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

/* Gravity Test: particles in isolated boundaries */
  int TestOrbitInitializeGrid(int NumberOfTestParticles,
			      FLOAT TestRadius,
			      float CentralMass,
			      float TestMass,
			      int UseBaryons);

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
			     float SphereDensity[MAX_SPHERES],
			     float SphereTemperature[MAX_SPHERES],
			     float SphereMetallicity[MAX_SPHERES],
			     FLOAT SpherePosition[MAX_SPHERES][MAX_DIMENSION],
			     float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
			     float SphereFracKeplarianRot[MAX_SPHERES],
			     float SphereTurbulence[MAX_SPHERES],
			     float SphereDispersion[MAX_SPHERES],
			     float SphereCutOff[MAX_SPHERES],
			     float SphereAng1[MAX_SPHERES],
			     float SphereAng2[MAX_SPHERES],
			     int   SphereNumShells[MAX_SPHERES],
			     int   SphereType[MAX_SPHERES],
			     int   SphereUseParticles,
			     float ParticleMeanDensity,
			     float UniformVelocity[MAX_DIMENSION],
			     int   SphereUseColour,
			     int   SphereUseMetals,
			     float InitialTemperature, 
			     float InitialDensity, int level);

  /* CosmologySimulation: initialize grid. */
  int CosmologySimulationInitializeGrid(
			  int   InitialGridNumber,
			  float CosmologySimulationOmegaBaryonNow,
			  float CosmologySimulationOmegaCDMNow,
			  float CosmologySimulationInitialTemperature,
			  char *CosmologySimulationDensityName,
			  char *CosmologySimulationTotalEnergyName,
			  char *CosmologySimulationGasEnergyName,
			  char *CosmologySimulationVelocityNames[],
			  char *CosmologySimulationParticlePositionName,
			  char *CosmologySimulationParticleVelocityName,
			  char *CosmologySimulationParticleMassName,
			  char *CosmologySimulationParticleTypeName,
			  char *CosmologySimulationParticleVelocityNames[],
			  int   CosmologySimulationSubgridsAreStatic,
			  int   TotalRefinement,
			  float CosmologySimulationInitialFrctionHII,
			  float CosmologySimulationInitialFrctionHeII,
			  float CosmologySimulationInitialFrctionHeIII,
			  float CosmologySimulationInitialFrctionHM,
			  float CosmologySimulationInitialFrctionH2I,
			  float CosmologySimulationInitialFrctionH2II,
#ifdef TRANSFER
			  float RadHydroInitialRadiationEnergy,
#endif
			  int   CosmologySimulationUseMetallicityField,
			  PINT &CurrentNumberOfParticles,
			  int CosmologySimulationManuallySetParticleMassRatio,
			  float CosmologySimulationManualParticleMassRatio,
			  int CosmologySimulationCalculatePositions);

  int CosmologyInitializeParticles(
		   char *CosmologySimulationParticleVelocityName,
		   char *CosmologySimulationParticleMassName,
		   char *CosmologySimulationParticleTypeName,
		   char *CosmologySimulationParticleVelocityNames[],
		   float CosmologySimulationOmegaBaryonNow,
		   int *Offset, int level);

  /* CosmologySimulation: initialize partitioned nested grids. */
  int NestedCosmologySimulationInitializeGrid(
			  int   InitialGridNumber,
			  float CosmologySimulationOmegaBaryonNow,
			  float CosmologySimulationOmegaCDMNow,
			  float CosmologySimulationInitialTemperature,
			  char *CosmologySimulationDensityName,
			  char *CosmologySimulationTotalEnergyName,
			  char *CosmologySimulationGasEnergyName,
			  char *CosmologySimulationVelocityNames[],
			  char *CosmologySimulationParticlePositionName,
			  char *CosmologySimulationParticleVelocityName,
			  char *CosmologySimulationParticleMassName,
			  char *CosmologySimulationParticleTypeName,
			  char *CosmologySimulationParticleVelocityNames[],
			  int   CosmologySimulationSubgridsAreStatic,
			  int   TotalRefinement,
			  float CosmologySimulationInitialFrctionHII,
			  float CosmologySimulationInitialFrctionHeII,
			  float CosmologySimulationInitialFrctionHeIII,
			  float CosmologySimulationInitialFrctionHM,
			  float CosmologySimulationInitialFrctionH2I,
			  float CosmologySimulationInitialFrctionH2II,
			  int   CosmologySimulationUseMetallicityField,
			  PINT &CurrentNumberOfParticles,
			  int CosmologySimulationManuallySetParticleMassRatio,
			  float CosmologySimulationManualParticleMassRatio,
			  int CosmologySimulationCalculatePositions);


  /* Initialization for isolated galaxy sims */
  int GalaxySimulationInitializeGrid(
				     FLOAT DiskRadius,
				     float GalaxyMass,
				     float GasMass,
				     FLOAT DiskPosition[MAX_DIMENSION], 
				     FLOAT ScaleHeightz,
				     FLOAT ScaleHeightR, 
				     float DMConcentration,
				     float DiskTemperature,
				     float InitialTemperature,
				     float AngularMomentum[MAX_DIMENSION],
				     float UniformVelocity[MAX_DIMENSION], 
				     int UseMetallicityField, 
				     float GalaxySimulationInflowTime,
				     float GalaxySimulationInflowDensity,
				     int level);

  /* Free expansion test */
  int FreeExpansionInitializeGrid(int FreeExpansionFullBox,
				  float FreeExpansionDensity,
				  double FreeExpansionEnergy,
				  float FreeExpansionMaxVelocity,
				  float FreeExpansionMass,
				  float FreeExpansionRadius,
				  float DensityUnits, float VelocityUnits,
				  float LengthUnits, float TimeUnits);

  /* Supernova restart initialize grid. */
  int SupernovaRestartInitialize(float EjectaDensity, float EjectaRadius,
				 float EjectaThermalEnergy, 
				 FLOAT EjectaCenter[3], int ColourField,
				 int *NumberOfCellsSet);

  /* Put Sink restart initialize grid. */
  int PutSinkRestartInitialize(int level ,int *NumberOfCellsSet);

  /* Free-streaming radiation test problem: initialize grid (SUCCESS or FAIL) */
  int FSMultiSourceInitializeGrid(float DensityConst, float V0Const, 
				  float V1Const, float V2Const, float TEConst, 
				  float RadConst, int local);

  /* FLD Radiation test problem: initialize grid (SUCCESS or FAIL) */
  int RadHydroConstTestInitializeGrid(int NumChem, float DensityConst, 
				      float V0Const, float V1Const, 
				      float V2Const, float IEConst, 
				      float EgConst, float HMassFrac, 
				      float InitFracHII, float InitFracHeII, 
				      float InitFracHeIII, int local);

  /* FLD Radiation ionization test problem: initialize grid (SUCCESS or FAIL) */
  int RHIonizationTestInitializeGrid(int NumChem, float DensityConst, 
				     float V0Const, float V1Const, 
				     float V2Const, float IEConst, 
				     float EgConst, float HMassFrac, 
				     float InitFracHII, float InitFracHeII, 
				     float InitFracHeIII, int local);

  /* FLD Radiation clump ionization problem: initialize grid (SUCCESS or FAIL) */
  int RHIonizationClumpInitializeGrid(int NumChem, float NumDensityIn, 
				      float NumDensityOut, float V0Const,
				      float V1Const, float V2Const,
				      float IEConstIn, float IEConstOut, 
				      float EgConst, float HMassFrac, 
				      float InitFracHII, float ClumpCenterX0, 
				      float ClumpCenterX1, float ClumpCenterX2, 
				      float ClumpRadius, int local);

  /* FLD Rad r^{-2} density ionization problem: initialize grid (SUCCESS or FAIL) */
  int RHIonizationSteepInitializeGrid(int NumChem, float NumDensity, 
				      float DensityRadius, float DensityCenter0, 
				      float DensityCenter1, float DensityCenter2, 
				      float V0Const, float V1Const, 
				      float V2Const, float IEConst, 
				      float EgConst, float InitFracHII, int local);

  /* FLD Radiation test problem: cosmological HII ioniztion (SUCCESS or FAIL) */
  int CosmoIonizationInitializeGrid(int NumChem, float VxConst, float VyConst, 
				    float VzConst, float IEConst, 
				    float EgConst, float InitFracHII, 
				    float OmegaBaryonNow, int local);

  /* FLD Radiation test problem: stream test (SUCCESS or FAIL) */
  int RadHydroStreamTestInitializeGrid(float DensityConst, float EgConst,
				       int RadStreamDim, int RadStreamDir,
				       int local);

  /* FLD Radiation test problem: pulse test (SUCCESS or FAIL) */
  int RadHydroPulseTestInitializeGrid(float DensityConst, float EgConst,
				      int RadPulseDim, int local);

  /* FLD Radiation test problem: grey Marshak wave test (SUCCESS or FAIL) */
  int RadHydroGreyMarshakWaveInitializeGrid(float DensityConst, float IEConst, 
			      	            float EgConst, int GreyMarshDir,
					    int local);

  /* FLD Radiation test problem: radiating shock test (SUCCESS or FAIL) */
  int RadHydroRadShockInitializeGrid(float DensityConst, float TEConst, 
			      	     float REConst, float VelocityConst,
                                     int ShockDir, int local);

  /* Cooling test initialization */
  int CoolingTestInitializeGrid();

  /* Reset internal energy to initial values for cooling test. */
  int CoolingTestResetEnergies();

/* Tricks for Random Forcing. */

  int ReturnNumberOfBaryonFields(){return NumberOfBaryonFields;};
  int SetNumberOfBaryonFields(int &number)
    {NumberOfBaryonFields = number; return 0;};
  int AppendForcingToBaryonFields();
  int DetachForcingFromBaryonFields();
  int RemoveForcingFromBaryonFields();
  int AddRandomForcing(float * norm, float dtTopGrid);
  int PrepareRandomForcingNormalization(float * GlobVal, int GlobNum);
  int ReadRandomForcingFields(FILE *main_file_pointer);

  int AddFields(int TypesToAdd[], int NumberOfFields);
  int DeleteObsoleteFields(int *ObsoleteFields, 
			   int NumberOfObsoleteFields);
 
  inline bool isLocal () {return MyProcessorNumber == ProcessorNumber; };

 private:
//  int ReadRandomForcingFields(FILE *main_file_pointer);
 public:

/* TurbulenceSimulation: initialize grid. */

#define TURBULENCE_INIT_PARAMETERS_DECL \
     float TurbulenceSimulationInitialDensity, \
     float TurbulenceSimulationInitialTemperature, \
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
     TurbulenceSimulationDensityName, \
     TurbulenceSimulationTotalEnergyName, \
     TurbulenceSimulationGasEnergyName, \
     TurbulenceSimulationVelocityNames, \
     TurbulenceSimulationRandomForcingNames, \
     TurbulenceSimulationSubgridsAreStatic, \
     TotalRefinement


 int TurbulenceSimulationInitializeGrid(TURBULENCE_INIT_PARAMETERS_DECL);

  // The following are private since they should only be called by
  // TurbulenceSimulationInitializeGrid()

 private:
//  int TurbulenceSimulationInitializeGrid(TURBULENCE_INIT_PARAMETERS_DECL);
 public:


/* Comoving coordinate expansion terms. */

  int ComovingExpansionTerms();

/* Adjust the gravity source terms for comoving coordinates. */

  int ComovingGravitySourceTerm();

/* Star Particle handler routine. */

  int StarParticleHandler(HierarchyEntry* SubgridPointer, int level
#ifdef EMISSIVITY
			  // pass in dtLevelAbove for calculation of Geoffrey's Emissivity0 baryon field 
			  , float dtLevelAbove
#endif
                         );

/* Particle splitter routine. */

  int ParticleSplitter(int level);

/* Magnetic field resetting routine. */

  int MagneticFieldResetter(int level);

/* Apply a time-action to a grid. */

  int ApplyTimeAction(int Type, float Parameter);

/* Routine to set the tracer particle velocities from the grid velocity. */

  int TracerParticleSetVelocity();

/* Creates tracer particles in this grid. */

  int TracerParticleCreateParticles(FLOAT LeftEdge[], FLOAT RightEdge[],
                                    FLOAT Spacing, PINT &TotalParticleCount);


/* ShearingBox: initialize grid. */

  int ShearingBoxInitializeGrid(float ThermalMagneticRatio, float fraction, 
				float ShearingGeometry, 
				int InitialMagneticFieldConfiguration);

  int ShearingBox2DInitializeGrid(float ThermalMagneticRatio, float fraction, 
				float ShearingGeometry, 
			       int InitialMagneticFieldConfiguration);

  int ShearingBoxStratifiedInitializeGrid(float ThermalMagneticRatio, float fraction, 
				float ShearingGeometry, 
				int InitialMagneticFieldConfiguration);
// -------------------------------------------------------------------------
// Analysis functions for AnalysisBaseClass and it's derivatives.
//

  // Flag cells that overlap a subgrid (used for analysis).
  int FlagRefinedCells(grid *Subgrid);
  
  inline int IsInVolume( Eflt32 *LeftEdge, Eflt32 *RightEdge ){
    for( int i = 0; i < GridRank; i++ ){
      if( (GridLeftEdge[i] >= RightEdge[i]) ||
	  (GridRightEdge[i] <= LeftEdge[i]) ){
	return FALSE;
      }
    }
     return TRUE;
  }
  inline int IsInVolume( Eflt64 *LeftEdge, Eflt64 *RightEdge ){
    for( int i = 0; i < GridRank; i++ ){
      if( (GridLeftEdge[i] >= RightEdge[i]) ||
	  (GridRightEdge[i] <= LeftEdge[i]) ){
	return FALSE;
      }
    }
     return TRUE;
  }
  inline int IsInVolume( Eflt128 *LeftEdge, Eflt128 *RightEdge ){
    for( int i = 0; i < GridRank; i++ ){
      if( (GridLeftEdge[i] >= RightEdge[i]) ||
	  (GridRightEdge[i] <= LeftEdge[i]) ){
	return FALSE;
      }
    }
     return TRUE;
  }



  // Check to see if a FLOAT point is in the grid.
  inline int PointInGrid( Eflt32 *point ){
    for( int i = 0; i < GridRank; i++ ){
      if( ((point[i] >= GridLeftEdge[i]) &&
	  (point[i] < GridRightEdge[i])) == FALSE )
	return FALSE;
    }
    return TRUE;
  }
  inline int PointInGrid( Eflt64 *point ){
    for( int i = 0; i < GridRank; i++ ){
      if( ((point[i] >= GridLeftEdge[i]) &&
	  (point[i] < GridRightEdge[i])) == FALSE )
	return FALSE;
    }
    return TRUE;
  }
  inline int PointInGrid( Eflt128 *point ){
    for( int i = 0; i < GridRank; i++ ){
      if( ((point[i] >= GridLeftEdge[i]) &&
	  (point[i] < GridRightEdge[i])) == FALSE )
	return FALSE;
    }
    return TRUE;
  }

  // Check to see if a FLOAT point is in the grid (excluding boundaries)
  inline int PointInGridNB( Eflt32 *point ){
    for( int i = 0; i < GridRank; i++ ){
      if( ((point[i] > GridLeftEdge[i]) &&
	  (point[i] < GridRightEdge[i])) == FALSE )
	return FALSE;
    }
    return TRUE;
  }
  inline int PointInGridNB( Eflt64 *point ){
    for( int i = 0; i < GridRank; i++ ){
      if( ((point[i] > GridLeftEdge[i]) &&
	  (point[i] < GridRightEdge[i])) == FALSE )
	return FALSE;
    }
    return TRUE;
  }
  inline int PointInGridNB( Eflt128 *point ){
    for( int i = 0; i < GridRank; i++ ){
      if( ((point[i] > GridLeftEdge[i]) &&
	  (point[i] < GridRightEdge[i])) == FALSE )
	return FALSE;
    }
    return TRUE;
  }

  // Flags a 3D array where the grid overlaps.
  // Very similar to the FastSib stuff. (I think.)
  void FlagGridArray( HierarchyEntry ***GridArray, int *dx,
		      float *cell_width, HierarchyEntry *my_HE );

  /* Includes for analysis tools */

//#ifdef ANALYSIS_TOOLS
#include "../anyl/Grid_AnalyzeClusters.h"
//#endif

#ifdef USE_PYTHON
    void ConvertToNumpy(int GridID, PyArrayObject *container[],
                        int ParentID, int level, FLOAT WriteTime);
#endif
//------------------------------------------------------------------------
// Methods for star formation
//------------------------------------------------------------------------

  Star *ReturnStarPointer(void) { return Stars; };
  int ReturnNumberOfStars(void) { return NumberOfStars; };
  void SetNumberOfStars(int value) { NumberOfStars = value; };

  /* Calculate enclosed mass within a radius */

  int GetEnclosedMass(Star *star, float radius, float &mass,
		      float &metallicity, float &coldgas_mass, 
		      float AvgVelocity[]);
  int GetEnclosedMass(FLOAT star_pos[], float radius, float &mass,
		      float &metallicity, float &coldgas_mass, 
		      float AvgVelocity[]);

  int RemoveParticle(int ID, bool disable=false);

  int AddFeedbackSphere(Star *cstar, int level, float radius, float DensityUnits,
			float LengthUnits, float VelocityUnits, 
			float TemperatureUnits, float TimeUnits, double EjectaDensity, 
			double EjectaMetalDensity, double EjectaThermalEnergy,
			int &CellsModified);

  int MoveAllStars(int NumberOfGrids, grid* FromGrid[], int TopGridDimension);

  int MoveAllStarsOld(int NumberOfGrids, grid* FromGrid[], int TopGridDimension);

  int CommunicationSendStars(grid *ToGrid, int ToProcessor);

  int TransferSubgridStars(int NumberOfSubgrids, grid* ToGrids[], int AllLocal);
  
  int FindNewStarParticles(int level);

  int FindAllStarParticles(int level);

  int MirrorStarParticles(void);

  int UpdateStarParticles(int level);

  int AddH2Dissociation(Star *AllStars);

  int ReturnStarStatistics(int &Number, float &minLife);

//------------------------------------------------------------------------
// Radiative transfer methods that don't fit in the TRANSFER define
//------------------------------------------------------------------------

  int IdentifyRadiativeTransferFields(int &kphHINum, int &gammaNum,
				      int &kphHeINum, int &kphHeIINum, 
				      int &kdissH2INum);

#ifdef TRANSFER
#include "PhotonGrid_Methods.h"
#endif /* TRANSFER */

//------------------------------------------------------------------------
// Vertex interpolation (for radiative transfer and streaming data)
//------------------------------------------------------------------------

  int ComputeVertexCenteredField(int Num);
  int ComputeCellCenteredField(int Num);
  
  float ComputeInterpolatedValue(int Num, int vci, int vcj, int vck, 
				 float mx, float my, float mz);
  
  int NeighborIndices(int index, int vi[]);
  int NeighborVertexIndices(int index, int vi[]);

  void DeleteInterpolatedFields() {
    int i;
    for (i = 0; i < NumberOfBaryonFields; i++)
      if (InterpolatedField[i] != NULL) {
	delete [] InterpolatedField[i];
	InterpolatedField[i] = NULL;
      }
  }

//-----------------------------------------------------------------------
//  Returns radiative cooling by component
//-----------------------------------------------------------------------

  int ComputeLuminosity(float *luminosity, int NumberOfLuminosityFields);
  int ComputeMetalLineLuminosity(float *total_luminosity, float *all_emis, 
				 float *temperature);


//------------------------------------------------------------------------
//  Shearing Boundary Conditions
//------------------------------------------------------------------------

  bool isTopGrid(){
    for(int i=0; i<GridRank; i++){
      if (CellWidth[i][0]<TopGridDx[i]*0.95) return FALSE;
      if (CellWidth[i][0]>TopGridDx[i]*1.05) return FALSE;
    }
    return TRUE;};
  
//------------------------------------------------------------------------
//  Misc.
//------------------------------------------------------------------------

  int FindMinimumParticleMass(float &min_mass, int level);

  int FindMassiveParticles(float min_mass, int level, FLOAT *pos[], int &npart,
			   int CountOnly);

//------------------------------------------------------------------------
//  Inline FOF halo finder and particle interpolation using a tree
//------------------------------------------------------------------------

  int MoveParticlesFOF(int level, FOF_particle_data* &P, 
		       int &Index, FOFData AllVars, float VelocityUnits, 
		       double MassUnits, int CopyDirection);

  int InterpolateParticlesToGrid(FOFData *D);

//------------------------------------------------------------------------
//  Grid star particles onto the AMR mesh
//------------------------------------------------------------------------

  int InterpolateStarParticlesToGrid(int NumberOfSPFields);  

//------------------------------------------------------------------------
// new hydro & MHD routines
//------------------------------------------------------------------------

  int SetNumberOfColours(void);
  int SaveSubgridFluxes(fluxes *SubgridFluxes[], int NumberOfSubgrids,
                        float *Flux3D[], int flux, float fluxcoef, float dt);
  void ZeroFluxes(fluxes *SubgridFluxes[], int NumberOfSubgrids);
  int RungeKutta2_1stStep(fluxes *SubgridFluxes[],
                          int NumberOfSubgrids, int level,
                          ExternalBoundary *Exterior);
  int RungeKutta2_2ndStep(fluxes *SubgridFluxes[],
                          int NumberOfSubgrids, int level,
                          ExternalBoundary *Exterior);
  int ReturnHydroRKPointers(float **Prim, bool ReturnMassFractions = true);
  int ReturnOldHydroRKPointers(float **Prim, bool ReturnMassFractions = true);
  int UpdateElectronDensity(void);
  int UpdatePrim(float **dU, float c1, float c2);
  int Hydro3D(float **Prim, float **dU, float dt,
	      fluxes *SubgridFluxes[], int NumberOfSubgrids,
	      float fluxcoef, int fallback);
  int TurbulenceInitializeGrid(float CloudDensity, float CloudSoundSpeed, FLOAT CloudRadius, 
			       float CloudMachNumber, float CloudAngularVelocity, float InitialBField,
			       int SetTurbulence, int CloudType, int TurbulenceSeed, int PutSink, 
			       int level, int SetBaryonFields);
  int Collapse3DInitializeGrid(int n_sphere,
			       FLOAT r_sphere[MAX_SPHERES],
			       FLOAT rc_sphere[MAX_SPHERES],
			       float rho_sphere[MAX_SPHERES],
			       float p_sphere[MAX_SPHERES],
			       float cs_sphere[MAX_SPHERES],
			       FLOAT sphere_position[MAX_SPHERES][MAX_DIMENSION],
			       float omega_sphere[MAX_SPHERES],
			       int   sphere_type[MAX_SPHERES],
			       float rho_medium, float p_medium, int level);
  int Collapse1DInitializeGrid(FLOAT r_sphere,
			       FLOAT rc_sphere,
			       float rho_sphere,
			       float p_sphere,
			       float cs_sphere,
			       float omega_sphere,
			       int   sphere_type,
			       float rho_medium, float p_medium);
  int AddSelfGravity(float coef);
  int SourceTerms(float **dU);
  int MHD1DTestInitializeGrid(float rhol, float rhor,
			      float vxl,  float vxr,
			      float vyl,  float vyr,
			      float vzl,  float vzr,
			      float pl,   float pr,
			      float Bxl,  float Bxr,
			      float Byl,  float Byr,
			      float Bzl,  float Bzr);
  int MHD2DTestInitializeGrid(int MHD2DProblemType, 
			      float RampWidth,
			      float rhol, float rhou,
			      float vxl,  float vxu,
			      float vyl,  float vyu,
			      float pl,   float pu,
			      float Bxl,  float Bxu,
			      float Byl,  float Byu);
  int MHD3DTestInitializeGrid(int MHD3DProblemType,
			      float rhol, float rhou,
			      float vxl,  float vxu,
			      float vyl,  float vyu,
			      float pl,   float pu,
			      float Bxl,  float Bxu,
			      float Byl,  float Byu);
  int CollapseMHD3DInitializeGrid(int n_sphere,
				  FLOAT r_sphere[MAX_SPHERES],
				  FLOAT rc_sphere[MAX_SPHERES],
				  float rho_sphere[MAX_SPHERES],
				  float p_sphere[MAX_SPHERES],
				  float cs_sphere[MAX_SPHERES],
				  FLOAT sphere_position[MAX_SPHERES][MAX_DIMENSION],
				  float omega_sphere[MAX_SPHERES], float Bnaught, float theta_B,
				  int   sphere_type[MAX_SPHERES],
				  float rho_medium, float p_medium, int level);
  int MHDTurbulenceInitializeGrid(float rho_medium, float cs_medium, float mach, 
				  float Bnaught, int seed, int level, int SetBaryonFields);

  int PrepareVelocityNormalization(double *v_rms, double *Volume);
  int NormalizeVelocities(Eflt factor);

  int GalaxyDiskInitializeGrid(int NumberOfHalos,
			       FLOAT HaloRadius[MAX_SPHERES],
			       FLOAT HaloCoreRadius[MAX_SPHERES],
			       float HaloDensity[MAX_SPHERES],
			       float HaloTemperature[MAX_SPHERES],
			       FLOAT HaloPosition[MAX_SPHERES][MAX_DIMENSION],
			       float HaloSpin[MAX_SPHERES],
			       float HaloVelocity[MAX_SPHERES][MAX_DIMENSION],
			       float HaloAngVel[MAX_SPHERES],
			       float HaloMagneticField,
			       FLOAT DiskRadius[MAX_SPHERES],
			       FLOAT DiskHeight[MAX_SPHERES],
			       float DiskDensity[MAX_SPHERES],
			       float DiskTemperature[MAX_SPHERES],
			       int   GalaxyType[MAX_SPHERES],
			       int   UseParticles, int UseGas,
			       float UniformVelocity[MAX_DIMENSION],
			       float MediumTemperature, float MediumDensity, int level);
  int AGNDiskInitializeGrid(float BlackHoleMass,
			    int BlackHoleType,
			    int DiskType,
			    float DiskDensity,
			    float DiskTemperature,
			    FLOAT DiskRadius,
			    FLOAT DiskHeight, 
			    int UseGas, int level);
  int MHDRK2_1stStep(fluxes *SubgridFluxes[], 
		     int NumberOfSubgrids, int level,
		     ExternalBoundary *Exterior);
  int MHDRK2_2ndStep(fluxes *SubgridFluxes[], 
		     int NumberOfSubgrids, int level,
		     ExternalBoundary *Exterior);
  int MHD3D(float **Prim, float **dU, float dt,
	    fluxes *SubgridFluxes[], int NumberOfSubgrids,
	    float fluxcoef, int fallback);
  int MHDSourceTerms(float **dU);
  int UpdateMHDPrim(float **dU, float c1, float c2);
  int SaveMHDSubgridFluxes(fluxes *SubgridFluxes[], int NumberOfSubgrids,
			   float *Flux3D[], int flux, float fluxcoef, float dt);
  int SetFloor();


  /* Poisson clean routines */

  int PoissonSolver(int level);

  int PoissonSolverSOR();
  int PoissonSolverSOR2();
  int PoissonSolverFFT();
  int PoissonSolverMultigrid();

  int PoissonSolverCGA(int difftype, double *divB_p);
  template <typename T> int multA(T* input, T* output,  int *MatrixStartIndex, int *MatrixEndIndex);
  template <typename T> int multA2(T* input, T* output,  int *MatrixStartIndex, int *MatrixEndIndex);
  template <typename T> T dot(T *a, T *b,  int size);
  template <typename T> int setNeumannBC(T* x, int *MatrixStartIndex, int *MatrixEndIndex,int type);
  template <typename T> int setDirichletBC(T* x, int *MatrixStartIndex, int *MatrixEndIndex);

  int PoissonCleanStep(int level);

  int GetIndex(int i, int j, int k) { 
    return i + j*(GridDimension[0])
      + k*(GridDimension[0])*(GridDimension[1]);
  }
  
  int PoissonSolverTestInitializeGrid(int TestType, float GeometryControl);

  

  int PrintToScreenBoundaries(float *field, char *display, int direction, int slice,
			      int check, float diffvalue);  
  int PrintToScreenBoundaries(float *field, char *display, int direction, int slice);
  int PrintToScreenBoundaries(float *field, char *display);
  int PrintToScreenBoundaries();
  int PrintToScreenBoundaries(int field);

  int getField(int i){return FieldType[i];};
  
  int ReduceWindBoundary();

  /* New particle routines */
  int CheckParticle();

  int ReturnMaximumParticleNumber();

  int ReturnNumberOfNewParticles() {
    int np = 0;
    for (int n = 0; n < NumberOfParticles; n++)
      if (ParticleNumber[n] < 0) np++;
    return np;
  };
  
  /* Non-ideal effects */

  int AddViscosity();
  int ComputeViscosity(float *viscosity, int DensNum);

  int AddAmbipolarDiffusion();

  int AddResistivity();
  int ComputeResistivity(float *resistivity, int DensNum);
  /* END OF NEW STANFORD HYDRO/MHD ROUTINES */

};

// inline int grid::ReadRandomForcingFields (FILE *main_file_pointer);

// inline int grid::TurbulenceSimulationInitializeGrid (TURBULENCE_INIT_PARAMETERS_DECL);


#endif
