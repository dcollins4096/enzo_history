/***********************************************************************
/
/  GRID CLASS (SOLVE THE HYDRO EQUATIONS, SAVING SUBGRID FLUXES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

// Solve the hydro equations with the solver, saving the subgrid fluxes
//

#include <stdio.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int f, int farray[], int n);

#ifdef PPM_LR
extern "C" void FORTRAN_NAME(ppm_lr)(
			  float *d, float *E, float *u, float *v, float *w,
			    float *ge,
                          int *grav, float *gr_ax, float *gr_ay, float *gr_az,
			  float *gamma, float *dt, int *cycle_number,
                            float dx[], float dy[], float dz[],
			  int *rank, int *in, int *jn, int *kn,
                            int is[], int ie[],
			  float gridvel[], int *flatten, int *ipresfree,
			  int *diff, int *steepen, int *idual, 
                            float *eta1, float *eta2,
			  int *num_subgrids, int leftface[], int rightface[],
			  int istart[], int iend[], int jstart[], int jend[],
			  float *standard, int dindex[], int Eindex[],
			  int uindex[], int vindex[], int windex[],
			    int geindex[], float *temp,
                          int *ncolour, float *colourpt, int *coloff,
                            int colindex[], int *ifloat_size);
#endif /* PPM_LR */

int grid::SolveHydroEquations(int CycleNumber, int NumberOfSubgrids, 
			      fluxes *SubgridFluxes[], int level)
{

  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || !UseHydro)
    return SUCCESS;

  JBPERF_START("grid_SolveHydroEquations");

  this->DebugCheck("SolveHydroEquations");

  if (NumberOfBaryonFields > 0) {

    /* initialize */

    int dim, i, j, field, size, subgrid, n, colnum[MAX_COLOR];   // MAX_COLOR is defined in fortran.def
    Elong_int GridGlobalStart[MAX_DIMENSION];
    FLOAT a = 1, dadt;

    /* Compute size (in floats) of the current grid. */

    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

    /* If multi-species being used, then treat them as colour variables
       (note: the solver has been modified to treat these as density vars). */

    int NumberOfColours = 0, ColourNum;
    if (MultiSpecies > 0) {
      NumberOfColours = 6 + 3*(MultiSpecies-1);

      if ((ColourNum =
           FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Could not find ElectronDensity.");
      }

      /* Generate an array of field numbers corresponding to the colour fields
	 (here assumed to start with ElectronDensity and continue in order). */

      for (i = 0; i < NumberOfColours; i++)
        colnum[i] = ColourNum+i;

    }

    /* Add metallicity as a colour variable. */

    int MetalNum;

    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1) {
      colnum[NumberOfColours++] = MetalNum;
      if (MultiMetals || TestProblemData.MultiMetals) {
	colnum[NumberOfColours++] = MetalNum+1;
	colnum[NumberOfColours++] = MetalNum+2;
      }
    }

    /* Add SN colour */

    int SNColourNum;
    if ((SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields))
	!= -1) {
      colnum[NumberOfColours++] = SNColourNum;
    }

    /* Add Simon Glover's chemistry species as color fields */

    if(TestProblemData.GloverChemistryModel){

      // Declarations for Simon Glover's cooling.
      int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
	DINum, DIINum, HDINum;

      int CINum,CIINum,OINum,OIINum,SiINum,SiIINum,SiIIINum,CHINum,CH2INum,
	CH3IINum,C2INum,COINum,HCOIINum,OHINum,H2OINum,O2INum;

      int GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

      if (IdentifyGloverSpeciesFields(HIINum,HINum,H2INum,DINum,DIINum,HDINum,
				      HeINum,HeIINum,HeIIINum,CINum,CIINum,OINum,
				      OIINum,SiINum,SiIINum,SiIIINum,CHINum,CH2INum,
				      CH3IINum,C2INum,COINum,HCOIINum,OHINum,H2OINum,
				      O2INum) == FAIL) {
	ENZO_FAIL("Error in IdentifyGloverSpeciesFields.");
      }

      colnum[NumberOfColours++] = HIINum;
      colnum[NumberOfColours++] = HINum;
      colnum[NumberOfColours++] = H2INum;

      if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
	colnum[NumberOfColours++] = DINum;
	colnum[NumberOfColours++] = DIINum;
	colnum[NumberOfColours++] = HDINum;
	colnum[NumberOfColours++] = HeINum;
	colnum[NumberOfColours++] = HeIINum;
	colnum[NumberOfColours++] = HeIIINum;
      }

      if( (GCM==3) || (GCM==5) || (GCM==7) ){
	colnum[NumberOfColours++] = COINum;
      }

      if( (GCM==2) || (GCM==3) || (GCM==7) ){
	colnum[NumberOfColours++] = CINum;
	colnum[NumberOfColours++] = CIINum;
	colnum[NumberOfColours++] = OINum;
	colnum[NumberOfColours++] = OIINum;
      }

      if( (GCM==2) || (GCM==3) ){
	colnum[NumberOfColours++] = SiINum;
	colnum[NumberOfColours++] = SiIINum;
	colnum[NumberOfColours++] = SiIIINum;
      }

      if( (GCM==3) || (GCM==7) ){
	colnum[NumberOfColours++] = CHINum;
	colnum[NumberOfColours++] = CH2INum;
	colnum[NumberOfColours++] = CH3IINum;
	colnum[NumberOfColours++] = C2INum;
	colnum[NumberOfColours++] = HCOIINum;
	colnum[NumberOfColours++] = OHINum;
	colnum[NumberOfColours++] = H2OINum;
	colnum[NumberOfColours++] = O2INum;
      }
      
    } // if(TestProblemData.GloverChemistryModel)


    /* Determine if Gamma should be a scalar or a field. */

    int UseGammaField = FALSE;
    float *GammaField;
    if (HydroMethod == Zeus_Hydro && MultiSpecies > 1) {
      UseGammaField = TRUE;
      GammaField = new float[size];
      if (this->ComputeGammaField(GammaField) == FAIL) {
	ENZO_FAIL("Error in grid->ComputeGammaField.");
      }
    } else {
      GammaField = new float[1];
      GammaField[0] = Gamma;
    }

    /* Set lowest level flag (used on Zeus hydro). */

    int LowestLevel = (level > MaximumRefinementLevel-1) ? TRUE : FALSE;

    /* Set minimum support (used natively in zeus hydro). */

    float MinimumSupportEnergyCoefficient = 0;
    if (UseMinimumPressureSupport == TRUE && level > MaximumRefinementLevel-1)
      if (this->SetMinimumSupport(MinimumSupportEnergyCoefficient) == FAIL) {
	ENZO_FAIL("Error in grid->SetMinimumSupport,");
      }

    /* allocate space for fluxes */

    /* Set up our restart dump fluxes container */
    this->SubgridFluxStorage = SubgridFluxes;
    this->NumberOfSubgrids = NumberOfSubgrids;

    for (i = 0; i < NumberOfSubgrids; i++) {
      for (dim = 0; dim < GridRank; dim++)  {

	/* compute size (in floats) of flux storage */

        size = 1;
        for (j = 0; j < GridRank; j++)
          size *= SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] -
                  SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] + 1;

	/* set unused dims (for the solver, which is hardwired for 3d). */

        for (j = GridRank; j < 3; j++) {
          SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] = 0;
          SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] = 0;
          SubgridFluxes[i]->RightFluxStartGlobalIndex[dim][j] = 0;
          SubgridFluxes[i]->RightFluxEndGlobalIndex[dim][j] = 0;
        }

	/* Allocate space (if necessary). */

        for (field = 0; field < NumberOfBaryonFields; field++) {
	  //
	  if (SubgridFluxes[i]->LeftFluxes[field][dim] == NULL)
	    SubgridFluxes[i]->LeftFluxes[field][dim]  = new float[size];
	  if (SubgridFluxes[i]->RightFluxes[field][dim] == NULL)
	    SubgridFluxes[i]->RightFluxes[field][dim] = new float[size];
	  for (n = 0; n < size; n++) {
	    SubgridFluxes[i]->LeftFluxes[field][dim][n] = 0;
	    SubgridFluxes[i]->RightFluxes[field][dim][n] = 0;
	  }
        }

	for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	     field++) {
          SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
          SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
	}

      }  // next dimension

      /* make things pretty */

      for (dim = GridRank; dim < 3; dim++)
        for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
          SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
          SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
	}

    } // end of loop over subgrids

    /* compute global start index for left edge of entire grid 
       (including boundary zones) */

    for (dim = 0; dim < GridRank; dim++)
      GridGlobalStart[dim] = nlongint((GridLeftEdge[dim]-DomainLeftEdge[dim])/(*(CellWidth[dim]))) -
	GridStartIndex[dim];

    /* fix grid quantities so they are defined to at least 3 dims */

    for (i = GridRank; i < 3; i++) {
      GridDimension[i]   = 1;
      GridStartIndex[i]  = 0;
      GridEndIndex[i]    = 0;
      GridVelocity[i]    = 0.0;
      GridGlobalStart[i] = 0;
    }

    /* If using comoving coordinates, multiply dx by a(n+1/2).
       In one fell swoop, this recasts the equations solved by solver
       in comoving form (except for the expansion terms which are taken
       care of elsewhere). */

    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	  == FAIL) {
	ENZO_FAIL("Error in CsomologyComputeExpansionFactors.");
      }

    /* Create a cell width array to pass (and convert to absolute coords). */

    float *CellWidthTemp[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CellWidthTemp[dim] = new float[GridDimension[dim]];
      for (i = 0; i < GridDimension[dim]; i++)
	if (dim < GridRank)
	  CellWidthTemp[dim][i] = float(a*CellWidth[dim][i]);
	else
	  CellWidthTemp[dim][i] = 1.0;
    }

    /* Prepare Gravity. */

    int GravityOn = 0, FloatSize = sizeof(float);
    if (SelfGravity || UniformGravity || PointSourceGravity)
      GravityOn = 1;

    /* Call Solver on this grid.
       Notice that it is hard-wired for three dimensions, but it does
       the right thing for < 3 dimensions. */
    /* note: Start/EndIndex are zero based */
        
    if (HydroMethod == PPM_DirectEuler)
      this->PPMDirectEuler(CycleNumber, NumberOfSubgrids, 
			   SubgridFluxes, CellWidthTemp, 
			   GridGlobalStart, GravityOn,
			   NumberOfColours, colnum);

    /* PPM LR has been withdrawn. */

    if (HydroMethod == PPM_LagrangeRemap) {
#ifdef PPM_LR
      FORTRAN_NAME(ppm_lr)(
			density, totalenergy, velocity1, velocity2, velocity3,
                          gasenergy,
			&GravityOn, AccelerationField[0],
                           AccelerationField[1],
                           AccelerationField[2],
			&Gamma, &dtFixed, &CycleNumber,
                          CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
			&GridRank, &GridDimension[0], &GridDimension[1],
                           &GridDimension[2], GridStartIndex, GridEndIndex,
			GridVelocity, &PPMFlatteningParameter,
                           &PressureFree,
			&PPMDiffusionParameter, &PPMSteepeningParameter,
                           &DualEnergyFormalism, &DualEnergyFormalismEta1,
                           &DualEnergyFormalismEta2,
			&NumberOfSubgrids, leftface, rightface,
			istart, iend, jstart, jend,
			standard, dindex, Eindex, uindex, vindex, windex,
			  geindex, temp,
                        &NumberOfColours, colourpt, coloff, colindex);
#else /* PPM_LR */
      ENZO_FAIL("PPM LR is not supported.");
#endif /* PPM_LR */
    }

    if (HydroMethod == Zeus_Hydro)
      if (this->ZeusSolver(GammaField, UseGammaField, CycleNumber, 
			   CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
			   GravityOn, NumberOfSubgrids, GridGlobalStart,
			   SubgridFluxes,
			   NumberOfColours, colnum, LowestLevel,
			   MinimumSupportEnergyCoefficient) == FAIL)
	ENZO_FAIL("");
	

    /* Clean up allocated fields. */

    delete [] GammaField;

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      delete [] CellWidthTemp[dim];

  /* If we're supposed to be outputting on Density, we need to update
     the current maximum value of that Density. */

    if(OutputOnDensity == 1){
      int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);
      for(i = 0; i < size; i++)
        max(BaryonField[DensNum][size], CurrentMaximumDensity);
    }

  }  // end: if (NumberOfBaryonFields > 0)

  this->DebugCheck("SolveHydroEquations (after)");

  JBPERF_STOP("grid_SolveHydroEquations");
  return SUCCESS;

}
