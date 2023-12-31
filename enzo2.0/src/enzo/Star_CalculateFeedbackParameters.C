/***********************************************************************
/
/  CALCULATE FEEDBACK SPHERE PARAMETERS
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: Ji-hoon Kim
/             July, 2009
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

void Star::CalculateFeedbackParameters(float &Radius, 
				       float RootCellWidth,
				       float SNe_dt, double &EjectaDensity,
				       double &EjectaThermalEnergy,
				       double &EjectaMetalDensity,
				       float DensityUnits, float LengthUnits, 
				       float TemperatureUnits, float TimeUnits,
				       float VelocityUnits, float dtForThisStar)
{

  // Parameters for the Stroemgen sphere in Whalen et al. (2004)
  const float	BirthRadius	  = 50;		// pc
  const float	WhalenTemperature = 20000;	// K
  const float	WhalenDensity	  = 1;	        // cm^-3
  const float	WhalenMaxVelocity = 35;		// km/s

  const double pc = 3.086e18, Msun = 1.989e33, Grav = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13, 
    k_b = 1.38e-16, m_h = 1.673e-24, c = 3.0e10, sigma_T = 6.65e-25, h=0.70;

  float StarLevelCellWidth;
  double EjectaVolume, SNEnergy, HeliumCoreMass, Delta_SF;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum;

  int igrid[MAX_DIMENSION], dim, index;
  int size=1;
  float mdot;

  Radius = 0.0;
  EjectaDensity = 0.0;
  EjectaThermalEnergy = 0.0;
  EjectaMetalDensity = 0.0;
  StarLevelCellWidth = RootCellWidth / powf(float(RefineBy), float(this->level));

  switch (this->FeedbackFlag) {
  case SUPERNOVA:  // pair-instability SNe
    Radius = PopIIISupernovaRadius * pc / LengthUnits;
    Radius = max(Radius, 3.5*StarLevelCellWidth);
    EjectaVolume = 4.0/3.0 * 3.14159 * pow(PopIIISupernovaRadius*pc, 3);
    EjectaDensity = Mass * Msun / EjectaVolume / DensityUnits;
    HeliumCoreMass = (13./24.) * (Mass - 20);
    SNEnergy = (5.0 + 1.304 * (HeliumCoreMass - 64)) * 1e51;
    EjectaMetalDensity = HeliumCoreMass * Msun / EjectaVolume / 
      DensityUnits;
    EjectaThermalEnergy = SNEnergy / (Mass * Msun) / VelocityUnits /
      VelocityUnits;

    // Exaggerate influence radius because the blastwave will enter
    // into some of the surrounding parent grids within the next
    // timestep if we inject the energy into a small radius.
    Radius *= 1.0;
    break;

  case STROEMGREN:
    Radius = BirthRadius * pc / LengthUnits;
    Radius = max(Radius, 1.5*StarLevelCellWidth);
    EjectaDensity = WhalenDensity * m_h / DensityUnits;
    EjectaThermalEnergy =
      WhalenTemperature / (TemperatureUnits * (Gamma-1.0) * 0.6);
    break;

  case FORMATION:
    Radius = 0;
    break;

  case CONT_SUPERNOVA:
    // Inject energy into a sphere
    Radius = StarClusterSNRadius * pc / LengthUnits;
    Radius = max(Radius, 2*StarLevelCellWidth);

    // Release SNe energy constantly over 16 Myr (t = 4-20 Myr), which is defined in Star_SetFeedbackFlag.C.
    //Delta_SF = StarMassEjectionFraction * Mass * SNe_dt * TimeUnits / (16.0*Myr);
    Delta_SF = StarMassEjectionFraction * Mass * dtForThisStar * 
      TimeUnits / (16.0*Myr);
    EjectaVolume = 4.0/3.0 * 3.14159 * pow(Radius*LengthUnits, 3);   
    EjectaDensity = Delta_SF * Msun / EjectaVolume / DensityUnits;   
    EjectaMetalDensity = EjectaDensity * StarMetalYield;
    EjectaThermalEnergy = StarClusterSNEnergy / Msun /   
      (VelocityUnits * VelocityUnits);
    break;

  } // ENDSWITCH FeedbackFlag
  
//    fprintf(stdout, "star::CFP:  EjectaThermalEnergy = %g, EjectaDensity = %g, 
//                Radius = %g, mdot = %g, dtForThisStar = %g\n", 
//    	    EjectaThermalEnergy, EjectaDensity, Radius, mdot, dtForThisStar);  

  return;
}




