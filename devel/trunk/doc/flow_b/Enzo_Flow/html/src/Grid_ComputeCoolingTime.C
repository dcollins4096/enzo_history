/***********************************************************************
/
/  GRID CLASS (COMPUTE THE COOLING TIME FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/

// Compute the cooling time

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);
extern "C" void FORTRAN_NAME(cool_multi_time)(
	float *d, float *e, float *ge, float *u, float *v, float *w, float *de,
	   float *HI, float *HII, float *HeI, float *HeII, float *HeIII,
           float *cooltime,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
           hydro_method *imethod,
        int *idual, int *ispecies, int *imetal, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, int *ih2co, 
	   int *ipiht,
	float *dt, float *aye, float *temstart, float *temend,
	float *utem, float *uxyz, float *uaye, float *urho, float *utim,
	float *eta1, float *eta2, float *gamma,
	float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, 
	   float *ciHeIa, 
	float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a, 
	float *reHeII2a, float *reHeIIIa, float *brema, float *compa,
	float *comp_xraya, float *comp_temp,
           float *piHI, float *piHeI, float *piHeII,
	float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI,
           float *metal,
	float *hyd01ka, float *h2k01a, float *vibha, float *rotha, 
	   float *rotla,
	float *gpldl, float *gphdl, float *HDltea, float *HDlowa,
	float *inutot, int *iradfield, int *nfreq, int *imetalregen,
	int *iradshield, float *avgsighp, float *avgsighep, float *avgsighe2p);

extern "C" void FORTRAN_NAME(cool_time)(
	float *d, float *e, float *ge, float *u, float *v, float *w,
           float *cooltime,
	int *in, int *jn, int *kn, int *nratec, int *iexpand, 
           hydro_method *imethod, int *idual, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
	float *dt, float *aye, float *temstart, float *temend,
	   float *fh, float *utem,
	float *eta1, float *eta2, float *gamma, float *coola);


int grid::ComputeCoolingTime(float *cooling_time)
{

  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  /* Compute the size of the fields. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find fields: density, total energy, velocity1-3. */

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  /* Find Multi-species fields. */

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      return FAIL;
    }

  /* Get easy to handle pointers for each variable. */

  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];

  /* Compute the cooling time. */

  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;    
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
    CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		      &TimeUnits, &VelocityUnits, Time);

    aUnits = 1.0/(1.0 + InitialRedshift);
  }
  float afloat = float(a);

  /* Renyue's metal cooling code. */

  int MetallicityField = FALSE, MetalNum = 0;
#ifdef CEN_METALS
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) 
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;
#endif /* CEN_METALS */

  /* Calculate the rates due to the radiation field. */

  if (RadiationFieldCalculateRates(Time+0.5*dtFixed) == FAIL) {
    fprintf(stderr, "Error in RadiationFieldCalculateRates.\n");
    return FAIL;
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

  /* Call the appropriate FORTRAN routine to do the work. */

  if (MultiSpecies)
    FORTRAN_NAME(cool_multi_time)(
       density, totalenergy, gasenergy, velocity1, velocity2, velocity3,
       BaryonField[DeNum], BaryonField[HINum], BaryonField[HIINum], 
       BaryonField[HeINum], BaryonField[HeIINum], BaryonField[HeIIINum], 
       cooling_time,
       GridDimension, GridDimension+1, GridDimension+2, 
          &CoolData.NumberOfTemperatureBins, &ComovingCoordinates,
          &HydroMethod, 
       &DualEnergyFormalism, &MultiSpecies, &MetallicityField, &GridRank,
       GridStartIndex, GridStartIndex+1, GridStartIndex+2,
          GridEndIndex, GridEndIndex+1, GridEndIndex+2,
          &CoolData.ih2co, &CoolData.ipiht,
       &dtFixed, &afloat, &CoolData.TemperatureStart,
          &CoolData.TemperatureEnd,
       &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
       &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
       CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI,
          CoolData.ciHeI, 
       CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII, 
          CoolData.reHeII1, 
       CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, &CoolData.comp,
       &CoolData.comp_xray, &CoolData.temp_xray,
          &CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
       BaryonField[HMNum], BaryonField[H2INum], BaryonField[H2IINum],
       BaryonField[DINum], BaryonField[DIINum], BaryonField[HDINum],
          BaryonField[MetalNum],
       CoolData.hyd01k, CoolData.h2k01, CoolData.vibh, 
          CoolData.roth, CoolData.rotl,
       CoolData.GP99LowDensityLimit, CoolData.GP99HighDensityLimit, 
          CoolData.HDlte, CoolData.HDlow,
       RadiationData.Spectrum[0], &RadiationFieldType, 
          &RadiationData.NumberOfFrequencyBins, 
          &RadiationFieldRecomputeMetalRates,
       &RadiationShield, &HIShieldFactor, &HeIShieldFactor, &HeIIShieldFactor);
  else
    FORTRAN_NAME(cool_time)(
       BaryonField[DensNum], BaryonField[TENum], BaryonField[GENum], 
          BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], 
          cooling_time,
       GridDimension, GridDimension+1, GridDimension+2, 
          &CoolData.NumberOfTemperatureBins, &ComovingCoordinates,
          &HydroMethod,
       &DualEnergyFormalism, &GridRank,
       GridStartIndex, GridStartIndex+1, GridStartIndex+2,
          GridEndIndex, GridEndIndex+1, GridEndIndex+2,
       &dtFixed, &afloat, &CoolData.TemperatureStart,
          &CoolData.TemperatureEnd, &CoolData.HydrogenFractionByMass,
          &TemperatureUnits,
       &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
          CoolData.EquilibriumRate);

  return SUCCESS;
}
