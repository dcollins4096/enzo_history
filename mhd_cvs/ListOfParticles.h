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
/* ListOfParticles declarations */

struct ListOfParticles {
  int NumberOfParticles;
  int NumberOfValues;
  float *ParticlePosition[MAX_DIMENSION];
  float *ParticleVelocity[MAX_DIMENSION];
  int   *ParticleIndex;
  float *ParticleRadius;
  float *ParticleValue[MAX_NUMBER_OF_BARYON_FIELDS];
  ListOfParticles *NextList;
};
