/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GALAXY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, Feb, 2004
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

#define Mpc (3.0856e24)         //Mpc [cm] 
#define SolarMass (1.989e33)    //Solar Mass [g]
#define GravConst (6.67e-8)     //Gravitational Constant [cm3g-1s-2]
#define pi (3.14159)
#define mh (1.67e-24)           //Mass of Hydrogen [g]
#define kboltz (1.381e-16)      //Boltzmann's Constant [ergK-1]


int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
float gasdev();
float gasvel(float radius, float DiskDensity, float ExpansionFactor, float GalaxyMass, float ScaleHeightR, float ScaleHeightz, float DMConcentration);
float av_den(float r, float DiskDensity, float ScaleHeightR, float ScaleHeightz, float cellwidth, float z,float xpos, float ypos, float zpos);


static float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits;


static int GalaxySimulationParticleCount = 0;

static float CosmologySimulationInitialFractionHII   = 1.2e-5;
static float CosmologySimulationInitialFractionHeII  = 1.0e-14;
static float CosmologySimulationInitialFractionHeIII = 1.0e-17;
static float CosmologySimulationInitialFractionHM    = 2.0e-9;
static float CosmologySimulationInitialFractionH2I   = 2.0e-20;
static float CosmologySimulationInitialFractionH2II  = 3.0e-14;

int grid::GalaxySimulationInitializeGrid(int NumberOfDisks, 
					 float DiskRadius[MAX_SPHERES], 
					 float GalaxyMass[MAX_SPHERES],
					 float GasMass[MAX_SPHERES], 
					 float DiskPosition[MAX_SPHERES][MAX_DIMENSION], 
					 float ScaleHeightz[MAX_SPHERES],
					 float ScaleHeightR[MAX_SPHERES], 
					 float DMConcentration[MAX_SPHERES],
					 float DiskTemperature[MAX_SPHERES],
					 float InitialTemperature,
					 int DiskUseParticles,
					 int DiskUseColour, 
					 float AngularMomentum[MAX_DIMENSION][MAX_DIMENSION],
					 float UniformVelocity[MAX_DIMENSION], int level)

{
  /* declarations */

  int dim, i, j, k, m, field, disk, size;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  float DiskDensity, DiskVelocityMag,dens2;


  
  //  FILE *fout = fopen("densities.dat", "w");
  
  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }
  int ColourNum = NumberOfBaryonFields;
  if (DiskUseColour)
    FieldType[NumberOfBaryonFields++] = Metallicity; /* fake it with metals */

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    NumberOfParticles = (DiskUseParticles > 0) ? 1 : 0;
    for (dim = 0; dim < GridRank; dim++)
      NumberOfParticles *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    return SUCCESS;
  }

  /* Set various units. */

  // float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, 
  // VelocityUnits, CriticalDensity = 1, BoxLength = 1, mu = 0.6;
  float CriticalDensity = 1, BoxLength = 1, mu = 0.6;
  float a, dadt, ExpansionFactor = 1;
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    ExpansionFactor = a/(1.0+InitialRedshift);
    CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		      &TimeUnits, &VelocityUnits, Time);
    CriticalDensity = 2.78e11*pow(HubbleConstantNow, 2); // in Msolar/Mpc^3
    BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
  }


  /* Loop over the set-up twice, once to count the particles, the second
     time to initialize them. */

  int SetupLoopCount, npart = 0;
  for (SetupLoopCount = 0; SetupLoopCount < 1+min(DiskUseParticles, 1);
       SetupLoopCount++) {

  /* Set densities */

  float BaryonMeanDensity = DiskUseParticles ? 0.1 : 1.0;
  if (DiskUseParticles == 2) BaryonMeanDensity = 0.9;
  float ParticleMeanDensity = 1.0 - BaryonMeanDensity, ParticleCount = 0;

  /* Set particles. */

  if (DiskUseParticles > 0 && SetupLoopCount > 0) {

    /* If particles already exist (coarse particles), then delete. */

    if (NumberOfParticles > 0) {
      delete ParticleMass;
      delete ParticleNumber;
      for (dim = 0; dim < GridRank; dim++) {
	delete ParticlePosition[dim];
	delete ParticleVelocity[dim];
      }
    }

    /* Use count from previous loop to set particle number. */

    NumberOfParticles = npart;
    npart = 0;

    /* Allocate space. */

    ParticleMass = new float[NumberOfParticles];
    ParticleNumber = new int[NumberOfParticles];
    for (dim = 0; dim < GridRank; dim++) {
      ParticlePosition[dim] = new FLOAT[NumberOfParticles];
      ParticleVelocity[dim] = new float[NumberOfParticles];
    }

    /* Particle values will be set below. */

  } // end: particle initialization

  /* Set up the baryon field. */

  /* compute size of fields */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* allocate fields */

  if (SetupLoopCount == 0)
    for (field = 0; field < NumberOfBaryonFields; field++)
      if (BaryonField[field] == NULL)
	BaryonField[field] = new float[size];

  /* Loop over the mesh. */

  float density, dens1, Velocity[MAX_DIMENSION],
    temperature, temp1, sigma, sigma1, colour;
  FLOAT r, x, y = 0, z = 0;
  int n = 0;

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {

	/* Compute position */

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	/* Loop over disks. */

	density = 1.0;
	temperature = temp1 = InitialTemperature;
	sigma = sigma1 = 0;
	colour = 1.0e-10;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Velocity[dim] = 0;
	for (disk = 0; disk < NumberOfDisks; disk++) {

	  /* Find distance from center. */

	  r = sqrt(pow(fabs(x-DiskPosition[disk][0]), 2) +
		   pow(fabs(y-DiskPosition[disk][1]), 2) +
		   pow(fabs(z-DiskPosition[disk][2]), 2) );
	  r = max(r, 0.1*CellWidth[0][0]);

	  if (r < DiskRadius[disk]) {

	      FLOAT xpos, ypos, zpos, xpos1, ypos1, zpos1, zheight, drad;

	      /* Loop over dims if using Zeus (since vel's face-centered). */

	      for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
		   dim++) {

		/* Compute position. */

		xpos = x-DiskPosition[disk][0] - 
		  (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
		ypos = y-DiskPosition[disk][1] -
		  (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
		zpos = z-DiskPosition[disk][2] -
		  (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);

		/* Compute z and r_perp (AngularMomentum is angular momentum 
		   and must have unit length). */    

		/* magnitude of z = r.L in L direction */

		zheight = AngularMomentum[disk][0]*xpos + 
		          AngularMomentum[disk][1]*ypos +
		          AngularMomentum[disk][2]*zpos;

		/* position in plane of disk */

		xpos1 = xpos - zheight*AngularMomentum[disk][0];
		ypos1 = ypos - zheight*AngularMomentum[disk][1];
		zpos1 = zpos - zheight*AngularMomentum[disk][2];
		drad = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1);


		/* Normalize the vector r_perp = unit vector pointing along plane of disk */

		xpos1 = xpos1/drad;
		ypos1 = ypos1/drad;
		zpos1 = zpos1/drad;



		/* If we're above the disk, then exit. */

		//		if (fabs(zheight) > max(5.0*ScaleHeightz[disk]*Mpc/LengthUnits, 
		//					2.0*CellWidth[0][0]))
		//		  break;


		DiskDensity = (GasMass[disk]*SolarMass/(8.0*pi*ScaleHeightz[disk]*Mpc*pow(ScaleHeightR[disk]*Mpc,2.0)))/DensityUnits;   //Code units (rho_0)

		DiskVelocityMag = gasvel(drad, DiskDensity, ExpansionFactor, GalaxyMass[disk], ScaleHeightR[disk], ScaleHeightz[disk], DMConcentration[disk]);

		if (dim == 0)

		  if (ScaleHeightz[disk]*Mpc/LengthUnits > CellWidth[0][0])
		    {

		      dens1=av_den(drad*LengthUnits, DiskDensity*DensityUnits, ScaleHeightR[disk]*Mpc, ScaleHeightz[disk]*Mpc, CellWidth[0][0]*LengthUnits, zheight*LengthUnits,xpos1*drad*LengthUnits,ypos1*drad*LengthUnits,zpos1*drad*LengthUnits);

		    }
		  else
		    {

		      dens1 = DiskDensity*exp(-drad/(ScaleHeightR[disk]*Mpc/LengthUnits))/pow(cosh(zheight/CellWidth[0][0]), 2);

		    }	  

		 dens2 = DiskDensity*exp(-drad/(ScaleHeightR[disk]*Mpc/LengthUnits))/pow(cosh(zheight/max((ScaleHeightz[disk]*Mpc/LengthUnits), CellWidth[0][0])), 2);

    //		fprintf(fout, "den1 %g\t den2  %g \t cellwdth %g \t drad %g\n", dens1, dens2, CellWidth[0][0],drad);
	    
		if (dens2 < density)
		  break;

		/* Compute velocity magnitude (divided by drad). 
		   This assumes PointSourceGravityPosition and Disk center 
		   are the same. */

		/* Compute velocty: L x r_perp. */

		if (dim == 0 || dim == 1)
		  Velocity[0] = DiskVelocityMag*(AngularMomentum[disk][1]*zpos1 -
				     AngularMomentum[disk][2]*ypos1);
		if (dim == 0 || dim == 2)
		  Velocity[1] = DiskVelocityMag*(AngularMomentum[disk][2]*xpos1 -
				     AngularMomentum[disk][0]*zpos1);
		if (dim == 0 || dim == 3)
		  Velocity[2] = DiskVelocityMag*(AngularMomentum[disk][0]*ypos1 -
				     AngularMomentum[disk][1]*xpos1);
		//				printf("pos = %g %g %g   velocity = %g %g %g  DiskVelocityMag = %g (%g cm/s)\n",
		//				       xpos, ypos, zpos, Velocity[0], Velocity[1], Velocity[2],
		//				       DiskVelocityMag, DiskVelocityMag*VelocityUnits);

	      } // end: loop over dims

	      //	    } // end: disk
	    
	    /* If the density is larger than the background (or the previous
	       disk), then set the velocity. */

	    if (dens2 > density) {
	      density = dens1;
	      if (temp1 == InitialTemperature)
		temp1 = DiskTemperature[disk];
	      temperature = temp1;
	      sigma = sigma1;
	      if (disk == 0)
		colour = dens1; /* only mark first disk */
	    }

	  } // end: if (r < DiskRadius)
	} // end: loop over disks

	/* Set density. */

	BaryonField[0][n] = density*BaryonMeanDensity;
	//	printf("pos = %g %g %g    dens = %g   temp = %g\n",
	//	       x, y, z, density, temperature);

	/* If doing multi-species (HI, etc.), set these. */

	if (MultiSpecies > 0) {
	  
	  BaryonField[HIINum][n] = CosmologySimulationInitialFractionHII *
	    CoolData.HydrogenFractionByMass * BaryonField[0][n] *
	    sqrt(OmegaMatterNow)/
	    (OmegaMatterNow*BaryonMeanDensity*HubbleConstantNow);
	    //	    (CosmologySimulationOmegaBaryonNow*HubbleConstantNow);
      
	  BaryonField[HeIINum][n] = CosmologySimulationInitialFractionHeII*
	    BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	  BaryonField[HeIIINum][n] = CosmologySimulationInitialFractionHeIII*
	    BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	  BaryonField[HeINum][n] = 
	    (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][n] -
	    BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];

	  if (MultiSpecies > 1) {
	    BaryonField[HMNum][n] = CosmologySimulationInitialFractionHM*
	      BaryonField[HIINum][n]* pow(temperature,float(0.88));
	    BaryonField[H2IINum][n] = CosmologySimulationInitialFractionH2II*
	      2.0*BaryonField[HIINum][n]* pow(temperature,float(1.8));
	    BaryonField[H2INum][n] = CosmologySimulationInitialFractionH2I*
	      BaryonField[0][n]*CoolData.HydrogenFractionByMass*pow(301.0,5.1)*
	      pow(OmegaMatterNow, float(1.5))/
	      (OmegaMatterNow*BaryonMeanDensity)/
	    //	      CosmologySimulationOmegaBaryonNow/
	      HubbleConstantNow*2.0;
	  }

	  BaryonField[HINum][n] = 
	    CoolData.HydrogenFractionByMass*BaryonField[0][n]
	    - BaryonField[HIINum][n];
	  if (MultiSpecies > 1)
	    BaryonField[HINum][n] -= BaryonField[HMNum][n]
	      + BaryonField[H2IINum][n]
	      + BaryonField[H2INum][n];

	  BaryonField[DeNum][n] = BaryonField[HIINum][n] + 
	    0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
	  if (MultiSpecies > 1)
	    BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] - 
	      BaryonField[HMNum][n];

	  /* Set Deuterium species (assumed to be negligible). */

	  if (MultiSpecies > 2) {
	    BaryonField[DINum][n] = CoolData.DeuteriumToHydrogenRatio*
	                              BaryonField[HINum][n];
	    BaryonField[DIINum][n] = CoolData.DeuteriumToHydrogenRatio*
	                             BaryonField[HIINum][n];
	    BaryonField[HDINum][n] = CoolData.DeuteriumToHydrogenRatio*
	                             BaryonField[H2INum][n];
	  }
	}

	/* If there is a colour field, set it. */

	if (DiskUseColour)
	  BaryonField[ColourNum][n] = colour;

	/* Set Velocities. */

	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[vel+dim][n] = Velocity[dim] + UniformVelocity[dim];

	/* Set energy (thermal and then total if necessary). */

	BaryonField[1][n] = temperature/TemperatureUnits/
                            ((Gamma-1.0)*mu);

	if (DualEnergyFormalism)
	  BaryonField[2][n] = BaryonField[1][n];

	if (HydroMethod != Zeus_Hydro)
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[1][n] += 0.5*pow(BaryonField[vel+dim][n], 2);

	if (BaryonField[1][n] <= 0)
	  printf("n = %d  temp = %g   e = %g\n", 0, temperature, 
	       BaryonField[1][0]);

	/* Set particles if being used (generate a number of particle
	   proportional to density). */

	if (DiskUseParticles) {
	  if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
	      j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
	      k >= GridStartIndex[2] && k <= GridEndIndex[2]  ) {
	    ParticleCount += density/pow(float(RefineBy), GridRank*level);
	    while (ParticleCount > 1) {
	      if (SetupLoopCount > 0) {
		ParticleMass[npart] = ParticleMeanDensity*
		                      pow(float(RefineBy), GridRank*level);
		ParticleNumber[npart] = GalaxySimulationParticleCount++;

		/* Set random position within cell. */

		ParticlePosition[0][npart] = x + 
		        CellWidth[0][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		ParticlePosition[1][npart] = y +
		        CellWidth[1][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		ParticlePosition[2][npart] = z +
		        CellWidth[2][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);

		/* Set bulk velocity. */

		for (dim = 0; dim < GridRank; dim++)
		  ParticleVelocity[dim][npart] = 
		    Velocity[dim]+UniformVelocity[dim];

		/* Add random velocity; */

		if (sigma != 0)
		  for (dim = 0; dim < GridRank; dim++)
		    ParticleVelocity[dim][npart] += 
		      gasdev()*sigma/VelocityUnits;
		
	      }
	      npart++;
	      ParticleCount -= 1.0;
	    }
	  } // end: if statement
	} // end: if DiskUseParticles
  
      } // end loop over grid

  } // end loop SetupLoopCount

  if (DiskUseParticles && debug)
    printf("GalaxySimulationInitialize: NumberOfParticles = %d\n", 
	   NumberOfParticles);

  //     fclose(fout);

  return SUCCESS;
}


float gasvel(float radius, float DiskDensity, float ExpansionFactor, float GalaxyMass, float ScaleHeightR, float ScaleHeightz, float DMConcentration)
{
  //Elizabeth Tasker 20/03/04

  double OMEGA=OmegaLambdaNow+OmegaMatterNow;                 //Flat Universe

  double r = radius*LengthUnits/100;    // Radius [m]

  double M_200 = GalaxyMass*SolarMass/1000.0;      // Virial Mass [kg]

  double H = sqrt(HubbleConstantNow*100*HubbleConstantNow*100*(OmegaLambdaNow+OmegaMatterNow*pow(ExpansionFactor,-3)-(OMEGA-1.)*pow(ExpansionFactor,-2)));                                

  double r_200 = (1.63e-2*pow(GalaxyMass,1.0/3.0)*pow((OmegaLambdaNow+OmegaMatterNow*pow(ExpansionFactor, -3)-(OMEGA-1.0)*pow(ExpansionFactor,-2)),-1.0/3.0)*ExpansionFactor*pow(H,-2.0/3.0)*pow(100,2.0/3.0))*Mpc/1.0e5;
  //virial radius [m]: M_200/M_Solar = GalaxyMass

  double M_gas, M_DM, M_Tot, Acc, V_Circ;
  double f_C = log(1.0+DMConcentration)-DMConcentration/(1.0+DMConcentration);
  double r_s = r_200/DMConcentration;  //[m]
    
  // Mass of gas disk and DM at given radius

      M_gas=8.0*M_PI*ScaleHeightz*Mpc/100*ScaleHeightR*Mpc/100*ScaleHeightR*Mpc/100*DiskDensity*DensityUnits*1000*exp(-r/(ScaleHeightR*Mpc/100))*(exp(r/(ScaleHeightR*Mpc/100))-r/(ScaleHeightR*Mpc/100)-1.0);
     
      M_DM=(M_200/f_C)*(log(1.0+r/r_s)-(r/r_s)/(1.0+r/r_s));

      if (SelfGravity==1){
	M_Tot=M_DM+M_gas;
      }
      else{
	M_Tot=M_DM;
      }
  // Set the point source gravity parameters.  This is the DM mass (in g)
  //   within rs.  Also set the core radius to rs in cm.

      PointSourceGravityConstant = (M_200/f_C)*(log(1.0+1.0)-1.0/(1.0+1.0))*1000.0;
      PointSourceGravityCoreRadius = r_s*100.0;
    
  // Force per unit mass on disk (i.e. acceleration) [ms-2]

      Acc=((GravConst/1000.0)*M_Tot)/(r*r);

  // Magnitude of Circular Velocity of disk 

      V_Circ = sqrt(r*Acc)*100;       //cms-1

      /*      printf("r = %g  M_Tot = %g  Acc = %g  M_DM = %g  M_gas = %g  f_C = %g\n",
	     r, M_Tot, Acc, M_DM, M_gas, f_C);
      printf("r_s = %g  DMConcentration = %g  r_200 = %g  r/r_s = %g\n",
	     r_s, DMConcentration, r_200, r/r_s);
      printf("EF = %g  H = %g  OMEGA = %g\n", ExpansionFactor, H, OMEGA);
      printf("radius = %g  v_circ = %g\n", radius, V_Circ);  */
      
      return (V_Circ/VelocityUnits);  //code units
}

float av_den(float r, float DiskDensity, float ScaleHeightR, float ScaleHeightz, float cellwidth, float z, float xpos, float ypos, float zpos)
{
  // routine to return the average gas density in a grid cell
  // Routine samples density in r plane of grid and averages
  // Assumes all input units are CGS!!

  int i,points;
  double den,r1,nx,ny,nz;

  points = 100;
  den = DiskDensity*exp(-r/ScaleHeightR)/(pow(cosh(z/(2.0*ScaleHeightz)),2));

  for (i=0;i<points;i++)
    {
      nx = drand48()*cellwidth-cellwidth/2.0;
      ny = drand48()*cellwidth-cellwidth/2.0;
      nz = drand48()*cellwidth-cellwidth/2.0;
      r1 = sqrt(pow((xpos+nx),2)+pow((ypos+ny),2)+pow((zpos+nz),2)); 
      den = den+DiskDensity*exp(-r1/ScaleHeightR)/(pow(cosh(z/(2.0*ScaleHeightz)),2));
    }

  double av_den = den/points;

  return (av_den/DensityUnits); //code unites

}
