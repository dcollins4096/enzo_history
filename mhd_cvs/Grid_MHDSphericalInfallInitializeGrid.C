//  GRID CLASS (INITIALIZE THE GRID FOR A MHD SPHERICAL INFALL TEST)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "CosmologyParameters.h"
#include "Grid.h"
#include "SphericalInfall.h"


int grid::MHDSphericalInfallInitializeGrid(float InitialPerturbation, 
					int UseBaryons, 
					float SphericalInfallOmegaBaryonNow,
					float SphericalInfallOmegaCDMNow,
					int SubgridIsStatic,
					float SphericalInfallMagneticfieldx,
					float SphericalInfallMagneticfieldy,
					float SphericalInfallMagneticfieldz
					  )
{
  /* declarations */

  int dim, i, j, k, n, size, field, vel;
  int ParticleDimension[MAX_DIMENSION], ParticleCenter[MAX_DIMENSION];;
  float DelParticle[MAX_DIMENSION], DelCenter[MAX_DIMENSION];;

  /* Error check. */

  if (fabs(SphericalInfallOmegaBaryonNow + SphericalInfallOmegaCDMNow
	   - 1.0) > 1.0e-4 ||
      fabs(OmegaMatterNow - 1.0) > 1.0e-4) {
    fprintf(stderr, "SphericalInfall only works for Omega = 1");
    return FAIL;
  }

  if (UseBaryons) {

    /* create fields */

    NumberOfBaryonFields = 0;
    FieldType[NumberOfBaryonFields++] = Density;
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;
    vel = NumberOfBaryonFields;
    FieldType[NumberOfBaryonFields++] = Velocity1;
    if (GridRank > 1) 
      FieldType[NumberOfBaryonFields++] = Velocity2;
    if (GridRank > 2)
      FieldType[NumberOfBaryonFields++] = Velocity3;
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Set particles. */

  if (NumberOfParticles > 0) {

    /* Allocate space. */

    ParticleMass = new float[NumberOfParticles];
    ParticleNumber = new int[NumberOfParticles];
    for (dim = 0; dim < GridRank; dim++) {
      ParticlePosition[dim] = new FLOAT[NumberOfParticles];
      ParticleVelocity[dim] = new float[NumberOfParticles];
    }

    /* Error check number of particles. */

    if (POW(nint(POW(NumberOfParticles, 1.0/float(GridRank))), 
	    float(GridRank)) != NumberOfParticles) {
      fprintf(stderr, "NumberOfParticles must be N^%d.\n", GridRank);
      return FAIL;
    }

    /* Set ParticleDimension to be the number of particle per dim. */

    for (dim = 0; dim < GridRank; dim++)
      ParticleDimension[dim] = 
	nint(POW(NumberOfParticles, 1.0/float(GridRank)));
    for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
      ParticleDimension[dim] = 1;
      ParticleCenter[dim] = 0;
    }

#ifdef UNUSED
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ParticleDimension[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;

    if (ParticleDimension[0]*ParticleDimension[1]*ParticleDimension[2] != 
	NumberOfParticles) {
      fprintf(stderr, "NumberOfParticles must equal active cells.\n");
      return FAIL;
    }
#endif /* UNUSED */

    /* Compute: DelParticle - the interval between particles.
                ParticleCenter - index of the particle at the center of pert.
		DelCenter   - the amount to shift all particles to insure
		              the center particle is exactly at center. */

    for (dim = 0; dim < GridRank; dim++) {
      DelParticle[dim]    = (GridRightEdge[dim]-GridLeftEdge[dim])/
	float(ParticleDimension[dim]);
      ParticleCenter[dim] = int((SphericalInfallCenter[dim]/
				 DelParticle[dim] - 0.5     ) );
      DelCenter[dim]  = ((SphericalInfallCenter[dim]/DelParticle[dim] - 0.5) -
			  float(ParticleCenter[dim]))*DelParticle[dim];
    }
    if (debug) {
      printf("SphericalInfallInitialize: DelCenter = %f %f %f\n", 
	     DelCenter[0], DelCenter[1], DelCenter[2]);
      printf("SphericalInfallInitialize: ParticleCenter = %d %d %d\n",
	     ParticleCenter[0], ParticleCenter[1], ParticleCenter[2]);
    }

    /* Compute the number of active cells. */

    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridEndIndex[dim] - GridStartIndex[dim] + 1;

    /* Set particle positions, velocities and masses. */

    n = 0;
    for (k = 0; k < ParticleDimension[2]; k++)
      for (j = 0; j < ParticleDimension[1]; j++)
	for (i = 0; i < ParticleDimension[0]; i++, n++) {

	  ParticlePosition[0][n] = (float(i) + 0.5)*DelParticle[0] +
	    GridLeftEdge[0] + DelCenter[0];
	  if (GridRank > 1)
	    ParticlePosition[1][n] = (float(j) + 0.5)*DelParticle[1] +
	      GridLeftEdge[1] + DelCenter[1];
	  if (GridRank > 2)
	    ParticlePosition[2][n] = (float(k) + 0.5)*DelParticle[2] +
	      GridLeftEdge[2] + DelCenter[2];

	  for (dim = 0; dim < GridRank; dim++)
	    ParticleVelocity[dim][n] = 0.0;

	  ParticleMass[n] = SphericalInfallOmegaCDMNow*float(size)/
	    float(NumberOfParticles);

	  ParticleNumber[n] = n;
	}

    /* Give the central particle density it's perturbation. */


    ParticleMass[  (ParticleCenter[2])*ParticleDimension[0] 
                                      *ParticleDimension[1]
                 + (ParticleCenter[1])*ParticleDimension[0]
                 + (ParticleCenter[0])] += 
		       InitialPerturbation*SphericalInfallOmegaCDMNow;

  }

  /* If requested, set up the baryon field. */

  if (UseBaryons) {

    /* compute size of fields */

    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

    /* allocate fields */

    for (field = 0; field < NumberOfBaryonFields; field++)
      if (BaryonField[field] == NULL)
	BaryonField[field] = new float[size];
    
    /* set density to cosmic mean and total energy to near zero. */

    for (i = 0; i < size; i++) {
      BaryonField[0][i] = SphericalInfallOmegaBaryonNow;
      BaryonField[1][i] = 1.0e-10;
    }

    /* Set central density. */

    BaryonField[0][  (GridDimension[2]/2)*GridDimension[0]*GridDimension[1]
                   + (GridDimension[1]/2)*GridDimension[0]
                   + (GridDimension[0]/2)] += 
		       InitialPerturbation*SphericalInfallOmegaBaryonNow;

    /* set velocities */
    
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < size; i++)
	BaryonField[vel+dim][i] = 0.0;

    /* If open boundary conditions, make the density into a spherical
       top hat. */
    
    float Center = 0.5*(DomainRightEdge[0] - DomainLeftEdge[0]) + 
      DomainLeftEdge[0] + 0.5*CellWidth[0][0];
    printf("SphericalInfall: Center = %g\n", Center);
    n = 0;
#if 0
    float xpos, ypos, zpos, radius;
    if (GravityBoundaryType == TopGridIsolated)
      for (k = 0; k < GridDimension[2]; k++) {
	zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	for (j = 0; j < GridDimension[1]; j++) {
	  ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  for (i = 0; i < GridDimension[0]; i++, n++) {
	    xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	    radius = sqrt(POW((zpos - Center), 2) + POW((ypos - Center),2) +
			  POW((xpos - Center), 2));
	    if (radius > Center - DomainLeftEdge[0] - 4.1*CellWidth[0][0]) {
	      BaryonField[0][n] *= 1.0e-5;
	      BaryonField[1][n] *= 1.0e+5; // to keep pressure balance
	    }
	  }
	}
      }
#endif


	   /* Determine the size of the fields. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  /* Allocate space for the Magnetic fields */
    CenteredB[0] = new float[size];
    Current[0] = new float[size];
    MagneticField[0] = new float[MagneticSize[0]];
    ElectricField[0] = new float[ElectricSize[0]];
    CenteredB[1] = new float[size]; 
	Current[1] = new float[size];
	MagneticField[1] = new float[MagneticSize[1]];
    ElectricField[1] = new float[ElectricSize[1]]; 
	CenteredB[2] = new float[size]; 
	Current[2] = new float[size]; 
	MagneticField[2] = new float[MagneticSize[2]];
    ElectricField[2] = new float[ElectricSize[2]];
    DivB = new float[size];

     //Initialize Magnetic Fields

  for(i=0;i < size; i++) {
	  CenteredB[0][i] = SphericalInfallMagneticfieldx;
	  CenteredB[1][i] = SphericalInfallMagneticfieldy;
	  CenteredB[2][i] = SphericalInfallMagneticfieldz;
	  
  }
  for(i=0;i<MagneticSize[0];i++) MagneticField[0][i] = SphericalInfallMagneticfieldx;
  for(i=0;i<MagneticSize[1];i++) MagneticField[1][i] = SphericalInfallMagneticfieldy;
  for(i=0;i<MagneticSize[2];i++) MagneticField[2][i] = SphericalInfallMagneticfieldz;
  

    
    /* Set internal energy if necessary. */

    if (DualEnergyFormalism)
		for (i = 0; i < size; i++) {
	BaryonField[2][i] = BaryonField[1][i];

    /* Add thermal part to kinetic energy. */

	BaryonField[1][i] = BaryonField[2][i];

    for (dim = 0; dim < GridRank; dim++)
	BaryonField[1][i] +=  0.5*BaryonField[vel+dim][i]*BaryonField[vel+dim][i]*BaryonField[0][i]
	        + 0.5*CenteredB[dim][i]*CenteredB[dim][i];
		}
  }

  /* Set static subgrid flag. */

  SubgridsAreStatic = SubgridIsStatic;

  return SUCCESS;
}
