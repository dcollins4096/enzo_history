#ifndef _ANALYSIS_PARAMETERS_H_
#define _ANALYSIS_PARAMETERS_H_

struct ProjectionParameters{
  char OutputName[MAX_LINE_LENGTH];
  // WARN WARN WARN
  // Hard code value, made up for no good reason.

  int ProjectionAxis[MAX_PROJECTIONS];
  float ProjectionLeftEdge[MAX_PROJECTIONS][MAX_DIMENSION];
  float ProjectionRightEdge[MAX_PROJECTIONS][MAX_DIMENSION];
  //float ProjectionLevel[MAX_PROJECTIONS];
  int ProjectionLevel[MAX_PROJECTIONS];
  char * ProjectionTitle[MAX_PROJECTIONS];
  int fields[MAX_PROJECTIONS][MAX_PROJ_FIELDS];
  int number_of_fields[MAX_PROJECTIONS];
  int output_number; // number appended to end of projection name


};

struct PDFGeneratorParameters{

  char OutputName[MAX_LINE_LENGTH];
  int max_level;
  int gplot;
  int num_pdf_bins;
  int output_number; // number appended to end of structure functions name
};


// Function prototypes

ProjectionParameters *DefaultProjectionParameters();
PDFGeneratorParameters *DefaultPDFGeneratorParameters();
#endif //#ifndef _ANALYSIS_PARAMTERS_H_
