/***********************************************************************
/
/  GRID CLASS (COMPUTE THE PRESSURE FIELD AT THE GIVEN TIME)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the pressure at the requested time.  The pressure here is
//   just the ideal-gas equation-of-state.
 
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
 
/* function prototypes */
 
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
 
int grid::ComputePressure(FLOAT time, float *pressure)
{
 
  /* declarations */
 
  float density, gas_energy, total_energy;
  float velocity1, velocity2 = 0, velocity3 = 0;
  int i, size = 1;
 
  /* Error Check */
 
  if (time < OldTime || time > Time) {
    fprintf(stderr, "requested time is outside available range.\n");
    return FAIL;
  }
 
  /* Compute interpolation coefficients. */
 
  float coef, coefold;
  if (Time - OldTime > 0)
    coef    = (time - OldTime)/(Time - OldTime);
  else
    coef    = 1;
 
  coefold = 1 - coef;
 
  /* Compute the size of the grid. */
 
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
#ifdef PPML
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, BXNum, BYNum, BZNum;
  if (this->IdentifyMHDQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum,BXNum,BYNum, BZNum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
#else //PPML
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
#endif //PPML 
  /* If using Zeus_Hydro, then TotalEnergy is really GasEnergy so don't
     subtract the kinetic energy term. */
 
  float OneHalf = 0.5;
  if (HydroMethod == Zeus_Hydro)
    OneHalf = 0.0;
 
  /* Loop over the grid, compute the thermal energy, then the pressure,
     the timestep and finally the implied timestep. */
 
  /* special loop for no interpolate. */

  if (time == Time)
    
#ifdef PPML
  switch( EquationOfState ){
      case 0: //Adiabatic
#endif //PPML
    for (i = 0; i < size; i++) {
      
      total_energy  = BaryonField[TENum][i];
      density       = BaryonField[DensNum][i];
      velocity1     = BaryonField[Vel1Num][i];
      if (GridRank > 1 OR_MHD ) // see macros_and_parameters for OR_MHD
	velocity2   = BaryonField[Vel2Num][i];
      if (GridRank > 2 OR_MHD )
	velocity3   = BaryonField[Vel3Num][i];
      
#ifdef PPML
      
      if( MHD_Used == 1 ){
	fprintf(stderr,"warning: pressure computation in PPML MHD not inspected.\n");
	gas_energy = total_energy - (OneHalf*density*(velocity1*velocity1 +
						   velocity2*velocity2 +
						   velocity3*velocity3)+
				     BaryonField[BXNum][i]*BaryonField[BXNum][i] +
				     BaryonField[BYNum][i]*BaryonField[BYNum][i] +
				     BaryonField[BZNum][i]*BaryonField[BZNum][i]);
	pressure[i] = (Gamma - 1.0 )*gas_energy;
	
      }else{
#endif //PPML
	
      /* gas energy = E - 1/2 v^2. */
	
      gas_energy    = total_energy - OneHalf*(velocity1*velocity1 +
					      velocity2*velocity2 +
					      velocity3*velocity3);
      
      pressure[i] = (Gamma - 1.0)*density*gas_energy;
      
#ifdef PPML
      }//mhd_used
#endif //PPML
      if (pressure[i] < tiny_number)
	pressure[i] = tiny_number;
      
    } // end of loop
#ifdef PPML
    break;
  case 1: //EquationOfState == 1, isothermal
    for(i=0;i<size;i++)
      pressure[i] =  IsothermalSoundSpeed *IsothermalSoundSpeed * BaryonField[DensNum][i];
    break;
  default:
    fprintf(stderr,"ComputePressure: %d is an invalid EquationOfState choice.\n", EquationOfState);
    return FAIL;
  }//EOS switch
#endif //PPML
  else
  
    /* general case: */
#ifdef PPML
    switch( EquationOfState ){
    case 0: //adiabatic
#endif //PPML
      
    for (i = 0; i < size; i++) {
      
      total_energy  = coef   *   BaryonField[TENum][i] +
  	              coefold*OldBaryonField[TENum][i];
      density       = coef   *   BaryonField[DensNum][i] +
  	              coefold*OldBaryonField[DensNum][i];
      velocity1     = coef   *   BaryonField[Vel1Num][i] +
	              coefold*OldBaryonField[Vel1Num][i];
      
      if (GridRank > 1 OR_MHD)
	velocity2   = coef   *   BaryonField[Vel2Num][i] +
 	              coefold*OldBaryonField[Vel2Num][i];
      if (GridRank > 2 OR_MHD)
	velocity3   = coef   *   BaryonField[Vel3Num][i] +
	              coefold*OldBaryonField[Vel3Num][i];
      
#ifdef PPML
      if( HydroMethod == PPM_Local ){
	fprintf(stderr,"Compute Pressure: Adiabatic PPML not yet installed\n");
	return FAIL;
      }
#endif //PPML
      /* gas energy = E - 1/2 v^2. */
      
      gas_energy    = total_energy - OneHalf*(velocity1*velocity1 +
					  velocity2*velocity2 +
					  velocity3*velocity3);
      
      pressure[i] = (Gamma - 1.0)*density*gas_energy;
      
      if (pressure[i] < tiny_number)
	pressure[i] = tiny_number;
      
    }
#ifdef PPML
      break;
    case 1: //EquationOfState == Isothermal
      for( i=0; i<size; i++)
	pressure[i] = IsothermalSoundSpeed *IsothermalSoundSpeed * (coef*BaryonField[DensNum][i]+
								    coefold*OldBaryonField[DensNum][i]);
      break;
    default:
      fprintf(stderr, "ComputePressure: %d not a defined equation of state\n", EquationOfState);
      return FAIL;
    }//eos switch
#endif //PPML

  /* Correct for Gamma from H2. */
  
  if (MultiSpecies > 1) {
 
    float TemperatureUnits = 1, number_density, nH2, GammaH2Inverse,
      GammaInverse = 1.0/(Gamma-1.0), x, Gamma1, temp;
    float DensityUnits, LengthUnits, VelocityUnits, TimeUnits;
 
    /* Find Multi-species fields. */
 
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum,
        H2IINum, DINum, DIINum, HDINum;
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      return FAIL;
    }
 
    /* Find the temperature units if we are using comoving coordinates. */
 
    if (ComovingCoordinates)
      if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			    &TimeUnits, &VelocityUnits, Time) == FAIL) {
	fprintf(stderr, "Error in CosmologyGetUnits.\n");
	return FAIL;
      }
 
    for (i = 0; i < size; i++) {
 
      number_density =
	  0.25*(BaryonField[HeINum][i]  + BaryonField[HeIINum][i] +
		BaryonField[HeIIINum][i]                        ) +
	        BaryonField[HINum][i]   + BaryonField[HIINum][i]  +
                BaryonField[DeNum][i];
 
      nH2 = 0.5*(BaryonField[H2INum][i]  + BaryonField[H2IINum][i]);
 
      /* First, approximate temperature. */
 
      if (number_density == 0)
	number_density = tiny_number;
      temp = max(TemperatureUnits*pressure[i]/(number_density + nH2), 1);
 
      /* Only do full computation if there is a reasonable amount of H2.
	 The second term in GammaH2Inverse accounts for the vibrational
	 degrees of freedom. */
 
      GammaH2Inverse = 0.5*5.0;
      if (nH2/number_density > 1e-3) {
	x = temp/6100.0;
	if (x < 10.0)
	  GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/POW(exp(x)-1.0,2));
      }
 
      Gamma1 = 1.0 + (nH2 + number_density) /
	             (nH2*GammaH2Inverse + number_density * GammaInverse);
	
      /* Correct pressure with improved Gamma. */
 
      pressure[i] *= (Gamma1 - 1.0)/(Gamma - 1.0);
 
    } // end: loop over i
 
  } // end: if (MultiSpecies > 1)
 
  /* To emulate the opacity limit in turbulent star formation 
     simulations */
  
  float Gamma1 = Gamma;
  if ((ProblemType == 60 || ProblemType == 61) && SelfGravity == 1)
    for (i=0; i<size; i++) {
      Gamma1 = min(Gamma + (log10(BaryonField[DensNum][i])-8.0)*0.3999/2.5, 1.4);
      pressure[i] *= (Gamma1 - 1.0)/(Gamma - 1.0);
    }


  return SUCCESS;
}
