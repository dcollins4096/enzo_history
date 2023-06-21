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
/  Restart A SHOCK POOL SIMULATION
/
/  written by: Hao Xu
/  date:       July, 2006
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
#include "ShockPoolGlobalData.h"

int ShockPoolRestart(FILE *fptr,  HierarchyEntry *TopGrid,
		       TopGridData &MetaData)
{
  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  float MachSquared, SoundSpeed1;
  float ShockPoolShockVel;
  const float TwoPi = 6.283185;

  /* set default parameters */

  ShockPoolAngle         = 0.0;    // x direction
  ShockPoolMachNumber    = 10.0;    // Velocity1 / SoundSpeed1

  ShockPoolDensity       = 1.0;    // Density in region 1 (preshock)
  ShockPoolPressure      = 1.0;    // Pressure in region 1
  ShockPoolVelocity[0]   = 0.0;    // LabVelocity in region 1
  ShockPoolVelocity[1]   = 0.0;    //  Note: these should all be zero for
  ShockPoolVelocity[2]   = 0.0;    //        the MachNumber to be correct


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

  /* Compute the speed of the shock itself. */

  ShockPoolShockSpeed = SoundSpeed1 * ShockPoolMachNumber;

  /* set the inflow boundary on the left, otherwise leave things alone. */

  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    MetaData.LeftFaceBoundaryCondition[dim] = inflow;


  /* For Zeus solver, subtract kinetic component from TotalEnergy. */

  if (HydroMethod == Zeus_Hydro)
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      ShockPoolTotalEnergy      -= 0.5*POW(ShockPoolVelocity[dim], 2);
      ShockPoolShockTotalEnergy -= 0.5*POW(ShockPoolShockVelocity[dim], 2);
    }

  return SUCCESS;

}
