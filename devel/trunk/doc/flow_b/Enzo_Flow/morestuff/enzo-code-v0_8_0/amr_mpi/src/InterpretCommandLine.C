/***********************************************************************
/
/  INTERPRET COMMAND-LINE OPTIONS
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "performance.h"
#include <ctype.h>
#include "macros_and_parameters.h"

/* function prototypes */

void PrintUsage(char *myname);
void my_exit(int status);

int InterpretCommandLine(int argc, char *argv[], char *myname,
			 int &restart, int &debug, int &extract,
			 int &InformationOutput,
			 int &OutputAsParticleData,
			 int &project, int &ProjectionDimension, 
			 int &ProjectionSmooth,
			 char *ParameterFile[],
			 int RegionStart[], int RegionEnd[], 
			 FLOAT RegionStartCoordinate[], 
			 FLOAT RegionEndCoordinate[], 
			 int &RegionLevel, int MyProcessorNumber)
{

  int dim;

  /* Initialize Region start/end. */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    RegionStart[dim] = INT_UNDEFINED;
    RegionEnd[dim]   = INT_UNDEFINED;
    RegionStartCoordinate[dim] = FLOAT_UNDEFINED;
    RegionEndCoordinate[dim]   = FLOAT_UNDEFINED;
  }
  RegionLevel = INT_UNDEFINED;

  /* Interpret command-line arguments */

  char c;
  bool event_initialized = false; // Whether arguments include -P event <event>

  while (--argc > 0 && (*++argv)[0] == '-')
    while ((c = *++argv[0]))
      switch (c) {

	/* get beginning of region selection (float coordinate). */

      case 'b':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"PSYM, &RegionStartCoordinate[dim++]) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	      fprintf(stderr, "%s: error reading Begin coordinates\n", myname);
	    return FAIL;
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* debug */

      case 'd':
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  debug = TRUE;
	break;

	/* get end of region selection (integer index). */

      case 'e':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%d", &RegionEnd[dim++]) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	      fprintf(stderr, "%s: error reading End indexes.\n", myname);
	    return FAIL;
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* get finish of region selection (float coordinate). */

      case 'f':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"PSYM, &RegionEndCoordinate[dim++]) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	      fprintf(stderr, "%s: error reading Finish coordinates\n",myname);
	    return FAIL;
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* help */

      case 'h':
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  PrintUsage(myname);
	my_exit(EXIT_SUCCESS);
	break;

	/* Information output */

      case 'i':
	InformationOutput = TRUE;
	break;

	/* level of region selection. */

      case 'l':
	if (--argc > 0) {
	  if (sscanf((*++argv), "%d", &RegionLevel) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	      fprintf(stderr, "%s: error reading level.\n", myname);
	    return FAIL;
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	else {
	  if (MyProcessorNumber == ROOT_PROCESSOR)
	    fprintf(stderr, "%s: Need to specify level.\n", myname);
	  return FAIL;
	}
	break;

	/* Smooth projection */

      case 'm':
	ProjectionSmooth = TRUE;
	break;

	/* Output as particle data */

      case 'o':
	OutputAsParticleData = TRUE;
	break;

	/* Project to a plane. */

      case 'p':
	if (--argc > 0) {
	  if (sscanf((*++argv), "%d", &ProjectionDimension) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	      fprintf(stderr, "%s: error reading ProjectionDimension.\n",
		      myname);
	    return FAIL;
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	else {
	  if (MyProcessorNumber == ROOT_PROCESSOR)
	    fprintf(stderr, "%s: Need to specify level.\n", myname);
	  return FAIL;
	}
	project = TRUE;
	break;

#ifdef USE_JBPERF
	//
	// -P mode <m>    m=1 next(), 2 advance(), 4 category(), 8 trace()
	// -P event <e>   e = xx-yyy for PAPI_XX_YYY, eg fp-ins for PAPI_FP_INS
	// -P dir <d>     where d is directory to dump to.  default JBPERF
	//
      case 'P':
	if (--argc >= 2) { // -P commands have two arguments not one
	  --argc;
	  ++argv;
	  if (strcmp(*argv,"mode")==0) {
	    jb::perf.mode(atoi(*++argv));
	  } else if (strcmp(*argv,"event")==0) {
	    jb::perf.event(*++argv);
	    event_initialized = true;
	  } else if (strcmp(*argv,"dir")==0) {
	    printf ("WARNING: -P dir not implemented yet--ignoring\n");
	    *++argv;
	  } else {
	    fprintf(stderr, "%s: Unrecognize argument \"%s\" for option -P.\n",
		    myname,*argv);
	    return FAIL;
	  }
	} else {
	  if (MyProcessorNumber == ROOT_PROCESSOR)
	    fprintf(stderr, "%s: Need to specify performance argument.\n", 
		    myname);
	  return FAIL;
	}
	while (*(argv[0]+1))
	  ++argv[0];
	break;
#endif
	/* restart file. */

      case 'r':
	restart = TRUE;
	break;

	/* get start of region selection. */

      case 's':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%d", &RegionStart[dim++]) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	      fprintf(stderr, "%s: error reading Start indexes.\n", myname);
	    return FAIL;
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* Extract section. */

      case 'x':
	extract = TRUE;
	break;

	/* Unknown */

      default:
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  fprintf(stderr, "%s: unknown command-line option: -%c.\n",myname,c);
	return FAIL;
	
      } // end of switch(c)

  /* Use "fp-ins" as default for "-P event" switch if no event specified */

#ifdef USE_JBPERF

  if (! event_initialized) {
    jb::perf.event("fp-ins");
  }

#endif
  
  /* Error check for number of parameters, and set parameter file. */
  
  if (argc != 1) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      PrintUsage(myname);
    return FAIL;
  }
  *ParameterFile = argv[0];

  return SUCCESS;
}


/* Explain how to run the program. */

void PrintUsage(char *myname)
{
  fprintf(stderr, 
	  "usage: %s [options] <param_file>\n"
	  "\n"
	  "   general options:\n"
	  "      -d                            display debug information\n"
	  "      -r                            restart\n"
	  "      -x                            extract\n"
	  "      -l <level>                    level of extract\n"
	  "      -p <dimension>                project to plane\n"
	  "      -m                            smooth projection\n"
	  "      -o                            output as particle data\n"
	  "      -h                            help\n"
	  "      -i                            information output\n"
	  "      -s <dim0> [<dim1> [<dim2>]]   start index region\n"
	  "      -e <dim0> [<dim1> [<dim2>]]   end index region\n"
	  "      -b <dim0> [<dim1> [<dim2>]]   begin coordinates\n"
	  "      -f <dim0> [<dim1> [<dim2>]]   finish coordinate region\n"
	  "\n"
	  "   performance options:\n"
	  "      -P mode <modeval>             set jbPerf mode\n"
	  "      -P event <eventname>          set jbPerf event\n"
	  "      -P dir <directory>            set jbPerf directory\n"
          ,myname);
}
