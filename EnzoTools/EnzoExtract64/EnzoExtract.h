#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#include "extern_hdf5.h"  // hdf 5 prototypes


// ---------------- THE USER CAN MODIFY THESE!  ----------------------
//
#define FIELD_VALUES_DOUBLE  /*  This controls the precision of the input
				 and output datasets (which can be 32 or
				 64 bit).  Acceptable values are:
 
				 FIELD_VALUES_FLOAT      (32 bit)
				 FIELD_VALUES_DOUBLE     (64 bit)
			     */


#define PACK_AMR  /*  This controls whether the IO is packed (Harkness
		      IO) vs. non-packed (original IO).  Acceptable
		      values are:

		      PACK_AMR          (packed IO on)
		      NO_PACK_AMR       (packed IO off)
		  */

#define DEBUG         1   // debug flag:  1 on, 0 off
#define VERBOSEDEBUG  1   // verbose debug flat:  1 on, 0 off -- very noisy!

#define HIER_RICK  /* This controls whether the hierarchy file type is 
		      the one needed by Rick Wagner's packed IO modes
		      or the "standard" one.  Acceptable values are:

		      HIER_RICK         (rick's version)
		      HIER_NORMAL       (normal version)
		   */

// 
// ------------------------------------------------------------------

#define MAX_LINE_LENGTH 512

#define HDF5_FILE_I4 H5T_IEEE_F32BE  
#define HDF5_FILE_I8 H5T_IEEE_F64BE  

#define MAX_EXTRACTION_FIELDS 20
#define MAX_NAME_LENGTH 256
#define HUGENUMBER 1.0e+20

// return calls
#define SUCCESS 1
#define FAILURE 0

#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))

#ifdef ORIGINAL_DEFINE
# define EXTERN
#else /* ORIGINAL_DEFINE */
# define EXTERN extern
#endif /* ORIGINAL_DEFINE */


// misc. basic info about the simulation
EXTERN int numberofgrids, rootgridsize, maxlevel, refineby;

// pointers to arrays of info about the individual grids
EXTERN int *gridlevel,*griddx,*griddy,*griddz;
EXTERN double *gridlex,*gridley,*gridlez,
  *gridrex,*gridrey,*gridrez;

// values taken in from the command line
EXTERN int outputrowmajor, outputmaxlevel;
EXTERN double xstart, ystart, zstart, xend, yend, zend;

// file names (hierarchy files, parameter files)
EXTERN char **gridfilenames,*parameterfile;

// information about the size of the extraction volume
EXTERN int xextractnumcells,yextractnumcells,zextractnumcells,extractionarraysize;


// pointers to extraction field name arrays, field value arrays
EXTERN char **ext_fieldnames;

#ifdef FIELD_VALUES_DOUBLE
EXTERN double **ext_fieldvalues;
#else
EXTERN float **ext_fieldvalues;
#endif

// values corresponding to the various extraction fields
// and is used for ext_fieldnames and ext_fieldvalues
EXTERN int densitynum,metalnum,electronnum,temperaturenum;

// this is the number of extraction fields and is used for
// various operations relating to ext_fieldnames and
// ext_fieldvalues
EXTERN int numberofextractionfields;

// grid flagging buffer
EXTERN int *flagbuffer,numberofgridcells;

// buffers full of grid datasets
#ifdef FIELD_VALUES_DOUBLE
EXTERN double **gridbuffers;
#else
EXTERN float **gridbuffers;
#endif

// in file EnzoExtract_GridStuff.C
void GetGridInfo(char *hierfilename);

// in file EnzoExtract_ParseInputs.C
void ParseInputs(int numarg,char *arguments[]);
void ReadParameterFile(char *filename);
void CheckBoundaryValues(void);
void HelpMe(void);

// in file EnzoExtract_Extract.C
int Extract(int gridnum);

// in file EnzoExtract_Misc.C
void CreateExtractionArrays(void);
void DeleteExtractionArrays(void);
void WriteExtractionArrays(void);
void SwitchArraysToRowOrder(void);
