#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <hdf5.h>

#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"


ProjectionParameters *DefaultProjectionParameters(){

  ProjectionParameters *new_parameters
    = new ProjectionParameters;

  strcpy(new_parameters->OutputName, "AnalysisOutput");
  
  for(int p=0;p<MAX_PROJECTIONS; p++){
    new_parameters->ProjectionAxis[p] = INT_UNDEFINED;
    //<dbg>
    //fprintf(stderr,"b Projection[%d] Axis = %d (iu = %d)\n"
    //,p, new_parameters->ProjectionAxis[p],
    //INT_UNDEFINED);
    //</dbg>
    for(int i=0;i<MAX_DIMENSION;i++){
      new_parameters->ProjectionLeftEdge[p][i] = 0;
      new_parameters->ProjectionRightEdge[p][i] = 0;
    }
    new_parameters->ProjectionTitle[p] = NULL;
    
    for(int i=0;i<MAX_PROJ_FIELDS;i++)
      new_parameters->fields[p][i] = FieldUndefined;
    
    new_parameters->ProjectionLevel[p] = -1;
    new_parameters->number_of_fields[p] = 0;
    new_parameters->output_number = 0;
  }
  return new_parameters;
  
}

PDFGeneratorParameters *DefaultPDFGeneratorParameters(){
  PDFGeneratorParameters *new_parameters
    = new PDFGeneratorParameters;
  strcpy(new_parameters->OutputName, "AnalysisOutput");
  new_parameters->max_level = 0;
  new_parameters->num_pdf_bins = 500;
  new_parameters->output_number = 0;
  new_parameters->gplot = TRUE;

  return new_parameters;
}
