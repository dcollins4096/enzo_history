/***********************************************************************
/
/  GRID CLASS (SOLVE THE COOLING/HEATING AND RATE EQUATIONS COUPLED 
               WITH RADIATIVE TRANSFER)
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:  July, 2005 to solve cool and rate equations simultaneously
/  modified2:  January, 2006 to incorporate radiative transfer
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"

/* This parameter controls whether the cooling function recomputes
   the metal cooling rates.  It is reset by RadiationFieldUpdate. */

extern int RadiationFieldRecomputeMetalRates;

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);
double ReturnWallTime();
extern "C" void FORTRAN_NAME(solve_rate_cool)(
	float *d, float *e, float *ge, float *u, float *v, float *w, float *de,
	float *HI, float *HII, float *HeI, float *HeII, float *HeIII,
	int *in, int *jn, int *kn, int *nratec, int *iexpand, 
           hydro_method *imethod,
        int *idual, int *ispecies, int *imetal, int *imcool, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, int *ih2co, 
	   int *ipiht,
	float *dt, float *aye, float *temstart, float *temend,
	float *utem, float *uxyz, float *uaye, float *urho, float *utim,
	float *eta1, float *eta2, float *gamma, float *fh, float *dtoh,
	float *k1a, float *k2a, float *k3a, float *k4a, float *k5a, 
	   float *k6a, float *k7a, float *k8a, float *k9a, float *k10a,
	float *k11a, float *k12a, float *k13a, float *k13dda, float *k14a, 
           float *k15a,
        float *k16a, float *k17a, float *k18a, float *k19a, float *k22a,
	float *k24, float *k25, float *k26, float *k27, float *k28, float *k29,
	   float *k30, float *k31,
	float *k50a, float *k51a, float *k52a, float *k53a, float *k54a,
	   float *k55a, float *k56a,
	float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, 
	   float *ciHeIa, 
	float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a, 
	float *reHeII2a, float *reHeIIIa, float *brema, float *compa,
	float *comp_xraya, float *comp_temp, 
           float *piHI, float *piHeI, float *piHeII,
	float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI,
           float *metal,
	float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
	float *gpldl, float *gphdl, float *HDltea, float *HDlowa,
	float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela,
	float *metala, int *n_xe, float *xe_start, float *xe_end,
	float *inutot, int *iradtype, int *nfreq, int *imetalregen,
	int *iradshield, float *avgsighp, float *avgsighep, float *avgsighe2p,
	int *iradtrans, int *iradcoupled, int *iradstep, int *ierr,
	float *kphHI, float *kphHeI, float *kphHeII, 
	float *kdissH2I, float *gammaHI, float *gammaHeI, float *gammaHeII);

int grid::SolveCoupledRateEquations()
{

  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  this->DebugCheck("SolveRadiativeCooling");

  /* Declarations */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  FLOAT a = 1.0, dadt;
    
  /* Find fields: density, total energy, velocity1-3. */

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }

  /* Find Multi-species fields. */

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      ENZO_FAIL("");
    }

  /* Find photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaHINum, gammaHeINum, gammaHeIINum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				      gammaHeINum, kphHeIINum, gammaHeIINum, 
				      kdissH2INum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyRadiativeTransferFields.\n");
    ENZO_FAIL("");
  }

  /* Get easy to handle pointers for each variable. */

  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];

  /* If using cosmology, compute the expansion factor and get units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  if (ComovingCoordinates) {

    if (CosmologyComputeExpansionFactor(PhotonTime+0.5*dtFixed, &a, &dadt) 
	== FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      ENZO_FAIL("");
    }

    aUnits = 1.0/(1.0 + InitialRedshift);

  }

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, PhotonTime) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    ENZO_FAIL("");
  }

  float afloat = float(a);

  /* Metal cooling codes. */

  int MetalCoolingType = FALSE, MetalNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) == -1)
    MetalNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1);

  // Double check if there's a metal field when we have metal cooling
  if (MetalCooling == JHW_METAL_COOLING) {
    if (MetalNum != -1)
      MetalCoolingType = JHW_METAL_COOLING;
    else {
      fprintf(stderr, 
	      "Warning: No metal field found.  Turning OFF MetalCooling.\n");
      MetalCooling = FALSE;
      MetalNum = 0;
    }
  }
  if (MetalCooling == CEN_METAL_COOLING)
    if (MetalNum != 1)
      MetalCoolingType = CEN_METAL_COOLING;
    else {
      fprintf(stderr, 
	      "Warning: No metal field found.  Turning OFF MetalCooling.\n");
      MetalCooling = FALSE;
      MetalNum = 0;
    }

  /* Calculate the rates due to the radiation field. */

  if (RadiationFieldCalculateRates(PhotonTime+0.5*dtFixed) == FAIL) {
    fprintf(stderr, "Error in RadiationFieldCalculateRates.\n");
    ENZO_FAIL("");
  }

  /* Set up information for rates which depend on the radiation field. */

  int RadiationShield = (RadiationFieldType == 11) ? TRUE : FALSE;

  /* Precompute factors for self shielding (this is the cross section * dx). */

  float HIShieldFactor = RadiationData.HIAveragePhotoHeatingCrossSection * 
                         double(LengthUnits) * CellWidth[0][0];
  float HeIShieldFactor = RadiationData.HeIAveragePhotoHeatingCrossSection * 
                          double(LengthUnits) * CellWidth[0][0];
  float HeIIShieldFactor = RadiationData.HeIIAveragePhotoHeatingCrossSection * 
                           double(LengthUnits) * CellWidth[0][0];

  /* Call the fortran routine to solve cooling equations. */

  int ierr = 0;
  int RTCoupledSolverIntermediateStep = TRUE;

  FORTRAN_NAME(solve_rate_cool)(
    density, totalenergy, gasenergy, velocity1, velocity2, velocity3,
    BaryonField[DeNum], BaryonField[HINum], BaryonField[HIINum], 
       BaryonField[HeINum], BaryonField[HeIINum], BaryonField[HeIIINum], 
    GridDimension, GridDimension+1, GridDimension+2, 
       &CoolData.NumberOfTemperatureBins, &ComovingCoordinates, &HydroMethod, 
    &DualEnergyFormalism, &MultiSpecies, &MetalFieldPresent, &MetalCoolingType, 
       &GridRank, GridStartIndex, GridStartIndex+1, GridStartIndex+2, 
       GridEndIndex, GridEndIndex+1, GridEndIndex+2,
       &CoolData.ih2co, &CoolData.ipiht,
    &dtPhoton, &afloat, &CoolData.TemperatureStart, &CoolData.TemperatureEnd,
    &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
    &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
       &CoolData.HydrogenFractionByMass, &CoolData.DeuteriumToHydrogenRatio,
    RateData.k1, RateData.k2, RateData.k3, RateData.k4, RateData.k5, 
       RateData.k6, RateData.k7, RateData.k8, RateData.k9, RateData.k10,
    RateData.k11, RateData.k12, RateData.k13, RateData.k13dd, RateData.k14, 
       RateData.k15, RateData.k16,
    RateData.k17, RateData.k18, RateData.k19, RateData.k22,
    &RateData.k24, &RateData.k25, &RateData.k26, &RateData.k27,
       &RateData.k28, &RateData.k29, &RateData.k30, &RateData.k31,
    RateData.k50, RateData.k51, RateData.k52, RateData.k53,
       RateData.k54, RateData.k55, RateData.k56,
    CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI,
       CoolData.ciHeI, 
    CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, 
    CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, &CoolData.comp,
    &CoolData.comp_xray, &CoolData.temp_xray,
       &CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
    BaryonField[HMNum], BaryonField[H2INum], BaryonField[H2IINum],
       BaryonField[DINum], BaryonField[DIINum], BaryonField[HDINum],
       BaryonField[MetalNum],
    CoolData.hyd01k, CoolData.h2k01, CoolData.vibh, CoolData.roth,CoolData.rotl,
    CoolData.GP99LowDensityLimit, CoolData.GP99HighDensityLimit, 
       CoolData.HDlte, CoolData.HDlow,
    CoolData.GAHI, CoolData.GAH2, CoolData.GAHe, CoolData.GAHp,
    CoolData.GAel,
       CoolData.metals, &CoolData.NumberOfElectronFracBins, 
       &CoolData.ElectronFracStart, &CoolData.ElectronFracEnd,
    RadiationData.Spectrum[0], &RadiationFieldType, 
          &RadiationData.NumberOfFrequencyBins, 
          &RadiationFieldRecomputeMetalRates,
    &RadiationShield, &HIShieldFactor, &HeIShieldFactor, &HeIIShieldFactor,
    &RadiativeTransfer, &RadiativeTransferCoupledRateSolver,
    &RTCoupledSolverIntermediateStep, &ierr,
    BaryonField[kphHINum], BaryonField[kphHeINum], BaryonField[kphHeIINum], 
    BaryonField[kdissH2INum], BaryonField[gammaHINum], BaryonField[gammaHeINum], 
    BaryonField[gammaHeIINum]);

  if (ierr) {
      fprintf(stdout, "Error in FORTRAN rate/cool solver\n");
      fprintf(stdout, "GridLeftEdge = %"FSYM" %"FSYM" %"FSYM"\n",
	      GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
      fprintf(stdout, "GridRightEdge = %"FSYM" %"FSYM" %"FSYM"\n",
	      GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
      fprintf(stdout, "GridDimension = %"ISYM" %"ISYM" %"ISYM"\n",
	      GridDimension[0], GridDimension[1], GridDimension[2]);
      ENZO_FAIL("");
  }

  return SUCCESS;

}
