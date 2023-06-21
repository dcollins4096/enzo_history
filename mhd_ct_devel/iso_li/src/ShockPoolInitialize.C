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
/  INITIALIZE A SHOCK POOL SIMULATION
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:
/
/  PURPOSE:
/    The shock pool sets up a system which introduces a shock from the
/    left boundary.  The initial active region is completely uniform, 
/    and wave enters via inflow boundary conditions.   
/
/    In the frame in which the shock is stationary, the definitions are:
/                     |
/       Velocity2     |   Velocity1
/       Density2      |   Density1
/       Pressure2     |   Pressure1
/                     |
/
/    The laboratory frame (in which LabVelocity1 = 0) is related to this
/      frame through:
/                       Velocity1 = LabVelocity1 - ShockVelocity
/                       Velocity2 = LabVelocity2 - ShockVelocity
/
/    The MachNumber = |     Velocity1 / SoundSpeed1 |
/                   = | ShockVelocity / SoundSpeed1 |
/                     (if LabVelocity1 = 0)
/
/    See also Mihalas & Mihalas (Foundations of Radiation Hydrodynamics,
/        p. 236) eq. 56-40
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#define DEFINE_STORAGE
#include "ShockPoolGlobalData.h"
#undef DEFINE_STORAGE

#ifdef HAOXU
  void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
#endif

int ShockPoolInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
#ifdef HAOXU
   char *GEName   = "GasEnergy";
#endif
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";

  /* parameter declarations */
#ifndef HAOXU
  float ShockPoolMachNumber;
#endif
  FLOAT ShockPoolSubgridLeft, ShockPoolSubgridRight;

  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION], 
       SubgridDims[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  float MachSquared, SoundSpeed1;
#ifdef HAOXU
  float ShockPoolShockVel;
#else
  float ShockPoolShockVel, ShockPoolPressure, ShockPoolShockPressure;
#endif
  const float TwoPi = 6.283185;

  /* set default parameters */

  ShockPoolAngle         = 0.0;    // x direction
  ShockPoolMachNumber    = 10.0;    // Velocity1 / SoundSpeed1

  ShockPoolDensity       = 1.0;    // Density in region 1 (preshock)
  ShockPoolPressure      = 1.0;    // Pressure in region 1
  ShockPoolVelocity[0]   = 0.0;    // LabVelocity in region 1
  ShockPoolVelocity[1]   = 0.0;    //  Note: these should all be zero for
  ShockPoolVelocity[2]   = 0.0;    //        the MachNumber to be correct

  ShockPoolSubgridLeft   = 0.0;    // start of subgrid
  ShockPoolSubgridRight  = 0.0;    // end of subgrid

#ifdef HAOXU
  int sphere;
  int NumberofSphere = 1;
  int SphereType[MAX_SPHERES];
  float SpherePosition[MAX_SPHERES][MAX_DIMENSION];
  float CoreDensity[MAX_SPHERES];
  float SphereRadius[MAX_SPHERES];
  float SphereCoreRadius[MAX_SPHERES];
  float Sphere_n[MAX_SPHERES];
   
    for (sphere = 0; sphere < MAX_SPHERES; sphere++) {
    SphereRadius[sphere]     = 1.0;
    SphereCoreRadius[sphere] = 0.1;
    CoreDensity[sphere]    = 1.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      SpherePosition[sphere][dim] = 0.5;
    }
      SphereType[sphere]       = 0;
      Sphere_n[sphere] = 0.0;
  }

   	
#endif //HAOXU

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "ShockPoolAngle = %"FSYM, &ShockPoolAngle);
    ret += sscanf(line, "ShockPoolMachNumber = %"FSYM, &ShockPoolMachNumber);

    ret += sscanf(line, "ShockPoolDensity = %"FSYM, &ShockPoolDensity);
    ret += sscanf(line, "ShockPoolPressure = %"FSYM, &ShockPoolPressure);
    ret += sscanf(line, "ShockPoolVelocity1 = %"FSYM, &ShockPoolVelocity[0]);
    ret += sscanf(line, "ShockPoolVelocity2 = %"FSYM, &ShockPoolVelocity[1]);
    ret += sscanf(line, "ShockPoolVelocity3 = %"FSYM, &ShockPoolVelocity[2]);

    ret += sscanf(line, "ShockPoolSubgridLeft = %"PSYM, &ShockPoolSubgridLeft);
    ret += sscanf(line, "ShockPoolSubgridRight = %"PSYM, &ShockPoolSubgridRight);
   
#ifdef HAOXU
    ret +=sscanf(line, "NumberofSphere = %d", &NumberofSphere);
    if (sscanf(line, "SphereType[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereType[%d] = %d", &sphere,
                    &SphereType[sphere]);
    if (sscanf(line, "CoreDensity[%d]", &sphere) > 0)
      ret += sscanf(line, "CoreDensity[%d] = %"FSYM, &sphere,
                    &CoreDensity[sphere]);
    if (sscanf(line, "SphereRadius[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereRadius[%d] = %"FSYM, &sphere,
                   &SphereRadius[sphere]);
    if (sscanf(line, "SphereCoreRadius[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereCoreRadius[%d] = %"FSYM, &sphere,
                   &SphereCoreRadius[sphere]);
    if (sscanf(line, "SpherePosition[%d]", &sphere) > 0)
      ret += sscanf(line, "SpherePosition[%d] = %"FSYM" %"FSYM" %"FSYM,
                    &sphere, &SpherePosition[sphere][0],
                    &SpherePosition[sphere][1],
                    &SpherePosition[sphere][2]);
      if (sscanf(line, "Sphere_n[%d]", &sphere) > 0)
      ret += sscanf(line, "Sphere_n[%d] = %"FSYM, &sphere,
                   &Sphere_n[sphere]);

   
#endif //HAOXU
    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "ShockPool") &&
        line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* Compute the physical variables in the postshock region */

  MachSquared            = ShockPoolMachNumber * ShockPoolMachNumber;
  ShockPoolShockDensity  = ShockPoolDensity *
                           ((Gamma + 1.0) * MachSquared      ) /
                           ((Gamma - 1.0) * MachSquared + 2.0);
  ShockPoolShockPressure = ShockPoolPressure *
                           (2.0 * Gamma * MachSquared - (Gamma - 1.0)) /
                           (Gamma + 1.0);
  SoundSpeed1 = sqrt(Gamma * ShockPoolPressure / ShockPoolDensity);
  ShockPoolShockVel = SoundSpeed1 * ShockPoolMachNumber * 
                      (1.0 - ShockPoolDensity / ShockPoolShockDensity);
  ShockPoolShockVelocity[0] = cos(ShockPoolAngle*TwoPi/360.)*ShockPoolShockVel;
  ShockPoolShockVelocity[1] = sin(ShockPoolAngle*TwoPi/360.)*ShockPoolShockVel;
  ShockPoolShockVelocity[2] = 0.0;

  /* Compute total energies */

  ShockPoolTotalEnergy = ShockPoolPressure/((Gamma - 1.0)*ShockPoolDensity);
  ShockPoolShockTotalEnergy = ShockPoolShockPressure/((Gamma - 1.0)*
						      ShockPoolShockDensity);

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    ShockPoolTotalEnergy      += 0.5*POW(ShockPoolVelocity[dim], 2);
    ShockPoolShockTotalEnergy += 0.5*POW(ShockPoolShockVelocity[dim], 2);
  }

#ifdef HAOXU
// in MHD case , the definition of Energy are different
  if(MHD_Used){
   ShockPoolTotalEnergy = ShockPoolPressure/(Gamma - 1.0);
  ShockPoolShockTotalEnergy = ShockPoolShockPressure/(Gamma - 1.0);

   for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    ShockPoolTotalEnergy      += 0.5*POW(ShockPoolVelocity[dim], 2)*ShockPoolDensity;
    ShockPoolShockTotalEnergy += 0.5*POW(ShockPoolShockVelocity[dim], 2)*ShockPoolShockDensity;
  }
  }

#endif


  /* Compute the speed of the shock itself. */

  ShockPoolShockSpeed = SoundSpeed1 * ShockPoolMachNumber;

  /* set the inflow boundary on the left, otherwise leave things alone. */

  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    MetaData.LeftFaceBoundaryCondition[dim] = inflow;

#ifdef HAOXU
    // shock is only from x direction
    MetaData.LeftFaceBoundaryCondition[1] = outflow;
    MetaData.LeftFaceBoundaryCondition[2] = outflow;
#endif

#ifdef HAOXU

    if (TopGrid.GridData->InitializeGridwithSphere(NumberofSphere,
                                                 SphereType,
                                                 SpherePosition,
                                                 CoreDensity, 
                                                 SphereRadius,
                                                 SphereCoreRadius,
                                                 Sphere_n,
                                                 ShockPoolDensity,
                                                 ShockPoolPressure,
                                                 ShockPoolVelocity) == FAIL) {
      fprintf(stderr, "Error in InitializeGridwithSphere (subgrid).\n");
      return FAIL;
    }
  


#else


  /* set up grid */

  if (TopGrid.GridData->InitializeUniformGrid(ShockPoolDensity, 
					      ShockPoolTotalEnergy,
					      ShockPoolTotalEnergy,
					      ShockPoolVelocity) == FAIL) {
    fprintf(stderr, "Error in InitializeUniformGrid.\n");
    return FAIL;
  }
 
#endif //HAOXU

  /* If requested, create a subgrid */

  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] = 
      nint((ShockPoolSubgridRight - ShockPoolSubgridLeft)/
	   ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	    float(MetaData.TopGridDims[dim])))
	*RefineBy;

  if (NumberOfSubgridZones[0] > 0) {

    /* create a new HierarchyEntry, attach to the top grid and fill it out */

    HierarchyEntry *Subgrid    = new HierarchyEntry;
    TopGrid.NextGridNextLevel  = Subgrid;
    Subgrid->NextGridNextLevel = NULL;
    Subgrid->NextGridThisLevel = NULL;
    Subgrid->ParentGrid        = &TopGrid;

    /* compute the dimensions and left/right edges for the subgrid */

    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*DEFAULT_GHOST_ZONES;
      LeftEdge[dim]    = ShockPoolSubgridLeft;
      RightEdge[dim]   = ShockPoolSubgridRight;
    }

#ifdef HAOXU
   
     /* create a new subgrid and initialize it */

    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
                                   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->InitializeGridwithSphere(NumberofSphere, SphereType,
                                                 SpherePosition,
                                                 CoreDensity, 
                                                 SphereRadius,
                                                 SphereCoreRadius,
                                                 Sphere_n,
                                                 ShockPoolDensity,
                                                 ShockPoolPressure,
                                                 ShockPoolVelocity) == FAIL) {
      fprintf(stderr, "Error in InitializeGridwithSphere (subgrid).\n");
      return FAIL;
    }
  }


#else 

    /* create a new subgrid and initialize it */

    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->InitializeUniformGrid(ShockPoolDensity, 
						 ShockPoolTotalEnergy,
						 ShockPoolTotalEnergy,
						 ShockPoolVelocity) == FAIL) {
      fprintf(stderr, "Error in InitializeUniformGrid (subgrid).\n");
      return FAIL;
    }			   
  }

#endif //HAOXU


#ifdef HAOXU
 /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;

  for (int i = 0; i < count; i++)
    DataUnits[i] = NULL;

  if(MHD_Used){
  MHDcLabel[0] = "MagneticField_C_1";
  MHDcLabel[1] = "MagneticField_C_2";
  MHDcLabel[2] = "MagneticField_C_3";

  MHDLabel[0] = "MagneticField_F_1";
  MHDLabel[1] = "MagneticField_F_2";
  MHDLabel[2] = "MagneticField_F_3";

  MHDeLabel[0] = "ElectricField_1";
  MHDeLabel[1] = "ElectricField_2";
  MHDeLabel[2] = "ElectricField_3";

  CurrentLabel[0] = "Current_1";
  CurrentLabel[1] = "Current_2";
  CurrentLabel[2] = "Current_3";
 }

#else
  /* set up field names and units */

  DataLabel[0] = DensName;
  DataLabel[1] = TEName;
  DataLabel[2] = Vel1Name;
  DataLabel[3] = Vel2Name;
  DataLabel[4] = Vel3Name;

  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
#endif
  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "ShockPoolAngle        = %f\n"  , ShockPoolAngle);
    fprintf(Outfptr, "ShockPoolMachNumber   = %f\n\n", ShockPoolMachNumber);

    fprintf(Outfptr, "ShockPoolDensity      = %f\n"  , ShockPoolDensity);
    fprintf(Outfptr, "ShockPoolPressure     = %f\n"  , ShockPoolPressure);
    fprintf(Outfptr, "ShockPoolVelocity1    = %f\n"  , ShockPoolVelocity[0]);
    fprintf(Outfptr, "ShockPoolVelocity2    = %f\n"  , ShockPoolVelocity[1]);
    fprintf(Outfptr, "ShockPoolVelocity3    = %f\n\n", ShockPoolVelocity[2]);

    fprintf(Outfptr, "ShockPoolSubgridLeft  = %"GOUTSYM"\n"  , ShockPoolSubgridLeft);
    fprintf(Outfptr, "ShockPoolSubgridRight = %"GOUTSYM"\n\n", ShockPoolSubgridRight);

#ifdef HAOXU
   for (sphere = 0; sphere < NumberofSphere; sphere++) {
      fprintf(Outfptr, "SphereType[%d] = %d\n", sphere,
              SphereType[sphere]);
      fprintf(Outfptr, "SphereRadius[%d] = %"GOUTSYM"\n", sphere,
              SphereRadius[sphere]);
      fprintf(Outfptr, "SphereCoreRadius[%d] = %"GOUTSYM"\n", sphere,
              SphereCoreRadius[sphere]);
      fprintf(Outfptr, "CoreDensity[%d] = %f\n", sphere,
              CoreDensity[sphere]);
      fprintf(Outfptr, "Sphere_n[%d] = %"GOUTSYM"\n", sphere,
              Sphere_n[sphere]);
      fprintf(Outfptr, "SpherePosition[%d] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
                        SpherePosition[sphere]);
//      fprintf(Outfptr, "SphereVelocity[%d] = ", sphere);
//      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
//                        SphereVelocity[sphere]);
    }
#endif //HAOXU

  }

  /* For Zeus solver, subtract kinetic component from TotalEnergy. */

  if (HydroMethod == Zeus_Hydro)
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      ShockPoolTotalEnergy      -= 0.5*POW(ShockPoolVelocity[dim], 2);
      ShockPoolShockTotalEnergy -= 0.5*POW(ShockPoolShockVelocity[dim], 2);
    }

  return SUCCESS;

}
