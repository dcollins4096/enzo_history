/***********************************************************************
/
/  ENZO Version and Options in Effect
/
/  written by: Robert Harkness
/  date:       January, 2006
/
************************************************************************/
 
#include <stdio.h>
 
#ifdef USE_MPI
#include <mpi.h>
#endif
 
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
#include "CosmologyParameters.h"
#include "svn_version.def"
 
#ifdef MEM_TRACE
Eint64 mused(void);
#endif
 
 
int ENZO_OptionsinEffect(void) 
{

  FILE *opf;

  if (MyProcessorNumber == 0) {

    opf = fopen("Enzo_Options", "w");

    fprintf(opf,"=========================\n");
    fprintf(opf,"Enzo SVN Branch   %s\n",ENZO_SVN_BRANCH);
    fprintf(opf,"Enzo SVN Revision %d\n",ENZO_SVN_REVISION);
    fprintf(opf,"=========================\n");

#ifdef SMALL_INTS
    fprintf(opf, " 32 bit Integer version\n");
#endif

#ifdef LARGE_INTS
    fprintf(opf, " 64 bit Integer version\n");
#endif

#ifdef INITS32
    fprintf(opf, " 32 bit Integer initial conditions\n");
#endif

#ifdef INITS64
    fprintf(opf, " 64 bit Integer initial conditions\n");
#endif

#ifdef r4
    fprintf(opf, " Float precision is 32 bits\n");
#endif

#ifdef r8
    fprintf(opf, " Float precision is 64 bits\n");
#endif

#ifdef p4
    fprintf(opf, " Position and time precision is 32 bits - NOT SUPPORTED!\n");
#endif

#ifdef p8
    fprintf(opf, " Position and time precision is 64 bits\n");
#endif

#ifdef p16
    fprintf(opf, " Position and time precision is 128 bits\n");
#endif

#ifdef HDF5_BE
    fprintf(opf, " HDF5 native i/o type is big-endian\n");
#endif

#ifdef HDF5_LE
    fprintf(opf, " HDF5 native i/o type is little-endian\n");
#endif


    fprintf(opf, "\n");
    fprintf(opf, " Optimizations in Effect\n");

#ifdef OOC_BOUNDARY
    fprintf(opf, " Out-of-core Top Grid boundary conditions\n");
#endif

#ifdef SIB1
    fprintf(opf, " Fast Sibling Locator 1\n");
#endif

#ifdef SIB2
    fprintf(opf, " Fast Sibling Locator 2\n");
#endif

#ifdef SIB3
    fprintf(opf, " Fast Sibling Locator 3\n");
#endif

#ifdef SIB4
    fprintf(opf, " Fast Sibling Locator 4\n");
#endif

#ifdef FLUX_FIX
    fprintf(opf, " New Flux Correction scheme by Collins & Wagner\n");
#endif

#ifdef SAB
    fprintf(opf, " AccelerationHack by Collins\n");
#endif

#ifdef UNIGRID
    fprintf(opf, " Minimum memory start-up => non-apative mesh only\n");
#endif

#ifdef UNIGRID_TRANSPOSE
    fprintf(opf, " Fast book-keeping for FFT TopGrid neighbours\n");
#endif

    fprintf(opf, "\n");
    fprintf(opf, "Macro and Parameter Definitions\n");

    fprintf(opf, "  MAX_NUMBER_OF_TASKS                 %8d\n", MAX_NUMBER_OF_TASKS);
    fprintf(opf, "  MAX_NUMBER_OF_BARYON_FIELDS (>=6)   %8d\n", MAX_NUMBER_OF_BARYON_FIELDS);
    fprintf(opf, "  MAX_NUMBER_OF_SUBGRIDS              %8d\n", MAX_NUMBER_OF_SUBGRIDS);
    fprintf(opf, "  MAX_DEPTH_OF_HIERARCHY              %8d\n", MAX_DEPTH_OF_HIERARCHY);
    fprintf(opf, "  MAX_LINE_LENGTH                     %8d\n", MAX_LINE_LENGTH);
    fprintf(opf, "  DEFAULT_GHOST_ZONES (>=3)           %8d\n", DEFAULT_GHOST_ZONES);
    fprintf(opf, "  MAX_NUMBER_OF_OUTPUT_REDSHIFTS      %8d\n", MAX_NUMBER_OF_OUTPUT_REDSHIFTS);
    fprintf(opf, "  GRAVITY_BUFFER_SIZE                 %8d\n", GRAVITY_BUFFER_SIZE);
    fprintf(opf, "  MAX_FLAGGING_METHODS                %8d\n", MAX_FLAGGING_METHODS);
    fprintf(opf, "  MAX_STATIC_REGIONS                  %8d\n", MAX_STATIC_REGIONS);
    fprintf(opf, "  MAX_NUMBER_OF_PARTICLE_ATTRIBUTES   %8d\n", MAX_NUMBER_OF_PARTICLE_ATTRIBUTES);
    fprintf(opf, "  MAX_TIME_ACTIONS                    %8d\n", MAX_TIME_ACTIONS);
    fprintf(opf, "  MAX_CUBE_DUMPS                      %8d\n", MAX_CUBE_DUMPS);
    fprintf(opf, "  MAX_POTENTIAL_ITERATIONS            %8d\n", MAX_POTENTIAL_ITERATIONS);

    fprintf(opf, "\n");

#ifdef SUN
    fprintf(opf, "Processor type is SUN\n");
#endif

#ifdef SUN_OLD
    fprintf(opf, "Processor type is SUN_OLD\n");
#endif

#ifdef SPP
    fprintf(opf, "Processor type is SPP\n");
#endif

#ifdef SP2
    fprintf(opf, "Processor type is SP2\n");
#endif

#ifdef BGL
    fprintf(opf, "Processor type is BGL\n");
#endif

#ifdef IRIS4
    fprintf(opf, "Processor type is IRIS4\n");
#endif

#ifdef COMPAQ
    fprintf(opf, "Processor type is COMPAQ\n");
#endif

#ifdef CONVEX
    fprintf(opf, "Processor type is CONVEX\n");
#endif

#ifdef LINUX
    fprintf(opf, "Processor type is LINUX\n");
#endif

#ifdef IA64
    fprintf(opf, "Processor type is IA64\n");
#endif

#ifdef CRAYX1
    fprintf(opf, "Processor type is CRAYX1\n");
#endif

#ifdef MEM_TRACE
    fprintf(opf, "Memory tracing enabled\n");
#endif

#ifdef TP_VELOCITY
    fprintf(opf, "Tracer particle velocity output is enabled\n");
#else
    fprintf(opf, "Tracer particle velocity output is disabled\n");
#endif

#ifdef ISO_GRAV
    fprintf(opf, "Gravity isolating B.C.s are enabled\n");
#else
    fprintf(opf, "Gravity isolating B.C.s are disabled\n");
#endif

#ifdef RAD_HYDRO
    fprintf(opf, "RHD is enabled\n");
#else
    fprintf(opf, "RHD is disabled\n");
#endif

    fclose(opf);

  } // processor zero only

  return SUCCESS;

}
