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
