/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE STAR PARTICLES
/
/  written by: Greg Bryan
/  date:       February, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#ifndef STAR_PARTICLE_DATA_DEFINED__
#define STAR_PARTICLE_DATA_DEFINED__

#ifdef DEFINE_STORAGE
# define SPEXTERN
#else /* DEFINE_STORAGE */
# define SPEXTERN extern
#endif /* DEFINE_STORAGE */

#define STAR_PARTICLE_NUMBER_START 1000000000

/* #include "macros_and_parameters.h" */

struct ParticleEntry {
  FLOAT Position[3];
  float Mass;
  float Velocity[3];
  float Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  PINT Number;
  int Type;
};


/* Number of Star particles. */

SPEXTERN int NumberOfStarParticles;
SPEXTERN int NumberOfDeletedParticles;
SPEXTERN PINT NumberOfOtherParticles; //all the particles other than type=2
SPEXTERN int G_TotalNumberOfStars;

/* Star particle parameters. */

SPEXTERN int StarFeedbackType;
SPEXTERN float StarMakerOverDensityThreshold;
SPEXTERN float StarMakerSHDensityThreshold;
SPEXTERN float StarMakerMassEfficiency;
SPEXTERN float StarMakerMinimumMass;
SPEXTERN float StarMakerMinimumDynamicalTime;
SPEXTERN float StarMassEjectionFraction;
SPEXTERN float StarMetalYield;
SPEXTERN float StarEnergyToThermalFeedback;
SPEXTERN float StarEnergyFeedbackRate;
SPEXTERN float StarEnergyToStellarUV;
SPEXTERN float StarEnergyToQuasarUV;

SPEXTERN float PopIIIStarMass;
SPEXTERN int   PopIIIBlackHoles;
SPEXTERN float PopIIIBHLuminosityEfficiency;
SPEXTERN float PopIIIOverDensityThreshold;
SPEXTERN float PopIIIH2CriticalFraction;
SPEXTERN float PopIIIMetalCriticalFraction;
SPEXTERN float PopIIISupernovaRadius;
SPEXTERN int   PopIIISupernovaUseColour;
SPEXTERN int   PopIIISupernovaMustRefine;
SPEXTERN int   PopIIISupernovaMustRefineResolution;
SPEXTERN float PopIIIColorDensityThreshold;
SPEXTERN float PopIIIColorMass;

SPEXTERN int    StarClusterUseMetalField;
SPEXTERN int    StarClusterHeliumIonization;
SPEXTERN float  StarClusterMinDynamicalTime;
SPEXTERN double StarClusterIonizingLuminosity;
SPEXTERN double StarClusterSNEnergy;
SPEXTERN float  StarClusterSNRadius;
SPEXTERN float  StarClusterFormEfficiency;
SPEXTERN float  StarClusterMinimumMass;
SPEXTERN float  StarClusterCombineRadius;
SPEXTERN float  StarClusterRegionLeftEdge[3];
SPEXTERN float  StarClusterRegionRightEdge[3];

SPEXTERN float minStarLifetime;
SPEXTERN FLOAT LastSupernovaTime;

#endif
