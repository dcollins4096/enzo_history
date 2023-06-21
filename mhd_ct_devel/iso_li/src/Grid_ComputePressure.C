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
  
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
  
  /* If using Zeus_Hydro, then TotalEnergy is really GasEnergy so don't
     subtract the kinetic energy term. */
  
  float OneHalf = 0.5;
  if (HydroMethod == Zeus_Hydro)
    OneHalf = 0.0;
  
  /* Loop over the grid, compute the thermal energy, then the pressure,
     the timestep and finally the implied timestep. */
  
  /* special loop for no interpolate. */
  
  /* Note: In enzo, energy is SPECIFIC energy; Total Energy Per unit Mass. dcc */

  
  if (time == Time){
    
#ifdef ATHENA
    if( EquationOfState == 0 ){
#endif //ATHENA

      for (i = 0; i < size; i++) {
	
	total_energy  = BaryonField[TENum][i];
	density       = BaryonField[DensNum][i];
	velocity1     = BaryonField[Vel1Num][i];
	if (GridRank > 1 || MHD_Used == TRUE)  //MHD needs 2.5d.
	  velocity2   = BaryonField[Vel2Num][i];
	if (GridRank > 2 || MHD_Used == TRUE)
	  velocity3   = BaryonField[Vel3Num][i];
	
	/* gas energy = E - 1/2 v^2. */
	
	
	if(MHD_Used){
	  
	  //hey, man-- if you change the Energy Variable, make sure
	  //you change the harten flag in CorrectForRefinedFluxes.
	  gas_energy = total_energy -density*OneHalf*(velocity1*velocity1 +
						      velocity2*velocity2 +
						      velocity3*velocity3) -
	    OneHalf*(CenteredB[0][i]*CenteredB[0][i]+
		     CenteredB[1][i]*CenteredB[1][i]+
		     CenteredB[2][i]*CenteredB[2][i]);

	  pressure[i] = (Gamma - 1.0)*gas_energy;



	  if( pressure[i] != pressure[i] ){
	    TVtool("bad pressure?");
	    fprintf(stderr,"Shit!  Bad Pressure!\n");
	    return FAIL;
	  }
	  
	}else{
	  gas_energy    = total_energy - OneHalf*(velocity1*velocity1 +
						  velocity2*velocity2 +
						  velocity3*velocity3);
	  
	  //pressure[i] = ( (Gamma - 1.0 ) < 0.01 )? density :(Gamma - 1.0)*gas_energy;
	  pressure[i] = (Gamma - 1.0)*density*gas_energy;
	  
	}				
	
#ifdef HAOXU
	if (pressure[i] < tiny_pressure)
	  pressure[i] = tiny_pressure;
#else      
	if (pressure[i] < tiny_number)
	  pressure[i] = tiny_number;
#endif /* HAOXU */
	
      } // end of loop
      
#ifdef ATHENA
    }else if( EquationOfState == 1 ){ 
      for( i=0;i<size;i++){
	pressure[i] = IsothermalSoundSpeed *IsothermalSoundSpeed * BaryonField[DensNum][i];
      }
    }else{ //EOS
      fprintf(stderr,"ComputePressure: %d is an invalid EquationOfState choice.\n", EquationOfState);
      return FAIL;
    }// end EOS
#endif //ATHENA
  }else{ //time==Time
    
      /* general case: */
      /* Note: In enzo, energy is SPECIFIC energy; Total Energy Per unit Mass */
      
#ifdef ATHENA 
      if( EquationOfState == 0 ){
#endif //ATHENA
      for (i = 0; i < size; i++) {
	total_energy  = coef   *   BaryonField[TENum][i] + 
	  coefold*OldBaryonField[TENum][i];
	density       = coef   *   BaryonField[DensNum][i] + 
	  coefold*OldBaryonField[DensNum][i];
	velocity1     = coef   *   BaryonField[Vel1Num][i] + 
	  coefold*OldBaryonField[Vel1Num][i];
	
	if (GridRank > 1)
	  velocity2   = coef   *   BaryonField[Vel2Num][i] + 
	    coefold*OldBaryonField[Vel2Num][i];
	if (GridRank > 2)
	  velocity3   = coef   *   BaryonField[Vel3Num][i] + 
	    coefold*OldBaryonField[Vel3Num][i];
	
	/* gas energy = E - 1/2 v^2. */
	
	if(MHD_Used){
	  gas_energy    = total_energy - density * OneHalf*(velocity1*velocity1 +
							    velocity2*velocity2 +
							    velocity3*velocity3);
	  
	  gas_energy  -= OneHalf*(coef*(CenteredB[0][i]*CenteredB[0][i]+
					CenteredB[1][i]*CenteredB[1][i]+
					CenteredB[2][i]*CenteredB[2][i]) +
				  coefold*(OldCenteredB[0][i]*OldCenteredB[0][i]+
					   OldCenteredB[1][i]*OldCenteredB[1][i]+
					   OldCenteredB[2][i]*OldCenteredB[2][i]));
	  
	  pressure[i] = density*gas_energy;//(Gamma - 1.0)*gas_energy;
	}else{ //MHD
	  gas_energy    = total_energy - OneHalf*(velocity1*velocity1 +
						  velocity2*velocity2 +
						  velocity3*velocity3);
	  
	  
	  pressure[i] = (Gamma - 1.0)*density*gas_energy;
	}
	/*
	  if(MHD_Used){
	  pressure[i]  -= OneHalf*(coef*(CenteredB[0][i]*CenteredB[0][i]+
	  CenteredB[1][i]*CenteredB[1][i]+
	  CenteredB[2][i]*CenteredB[2][i]) +
	  coefold*(OldCenteredB[0][i]*OldCenteredB[0][i]+
	  OldCenteredB[1][i]*OldCenteredB[1][i]+
	  OldCenteredB[2][i]*OldCenteredB[2][i]));
	  }
	*/
      	
#ifdef HAOXU
	if (pressure[i] < tiny_pressure)
	  pressure[i] = tiny_pressure;
#else 
	if (pressure[i] < tiny_number)
	  pressure[i] = tiny_number;
#endif /* HAOXU */     
      }//i
#ifdef ATHENA
    }else if( EquationOfState == 1 ){
      for( i=0;i<size;i++){
	pressure[i] = IsothermalSoundSpeed *IsothermalSoundSpeed * BaryonField[DensNum][i];
      }
    }else{  //EOS
      fprintf(stderr,"ComputePressure: %d is an invalid EquationOfState choice.\n", EquationOfState);
      return FAIL;
    }// end EOS
#endif //ATHENA
  }// time
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
  
  return SUCCESS;
}


