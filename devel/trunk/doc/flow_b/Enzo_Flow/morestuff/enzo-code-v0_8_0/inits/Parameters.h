/***********************************************************************
/
/  STRUCTURE FOR PARAMETERS
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

struct parmstruct {

  /* Size parameters */

  int Rank;
  int GridDims[3];
  int ParticleDims[3];
  int MaxDims[3];
  int NewCenter[3];
  int StartIndex[3];
  int GridRefinement;
  int ParticleRefinement;

  /* Temporary parameters (used to set other parameters, in
     ReadParameterFile and then not used after). */

  float NewCenterFloat[3];
  int StartIndexInNewCenterTopGridSystem[3];
  int EndIndexInNewCenterTopGridSystem[3];
  int RootGridDims[3];

  /* Boolean flags. */

  int InitializeParticles;
  int InitializeGrids;

  /* Names. */

  char *ParticlePositionName;
  char *ParticleVelocityName;
  char *ParticleMassName;
  char *GridDensityName;
  char *GridVelocityName;

  /* Power spectrum. */

  int WaveNumberCutoff;

};
