

//
// void ReadAnalysisParameterFile(char * Filename, FILE * fptr)
// 
// Reads parameters for the Inline Analysis routines.
// If Filename isn't found or not given (NULL) then the (assumed to be open) fptr is used.
//


#include <stdio.h>
#include <string.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"

void ReadIntsFromString(char * string, int * array, int size_of_array);
void ReadFloatFromString(char * string, float * array, int size_of_array);
void ReadTokensFromString(char * string, char ** array, char * delim,int size_of_array);
field_type_int get_field_id(char *field_name);

void ReadAnalysisParameterFile(char * Filename, FILE * fptr,TopGridData &MetaData){

  FILE * fptrB;
  if( Filename != NULL ){
    fptrB = fopen(Filename,"r");
    fprintf(stderr, "Opening AnalysisParameter %s\n",Filename);
    if( fptrB == NULL ){
      fprintf(stderr,
	      "Error reading AnalysisParameterFile '%s'.  Reading AnalysisParameters from primary file.\n",
	      Filename);
      fptrB = fptr;
    }
  }else{
    fptrB = fptr;
  }

  int ret = 0;
  int proj, this_id;
  char line[MAX_LINE_LENGTH];
  char *dummy = new char[MAX_LINE_LENGTH];
  int int_dummy;

  char * ProjectionFieldList_Char[MAX_PROJ_FIELDS];
  for( proj=0;proj<MAX_PROJ_FIELDS; proj++)
    ProjectionFieldList_Char[proj] = NULL;

  ProjectionParameters *proj_params = DefaultProjectionParameters();
  int proj_analysis_on = FALSE;

  PDFGeneratorParameters *pdf_params = DefaultPDFGeneratorParameters();
  int pdf_analysis_on = FALSE;


  while (fgets(line, MAX_LINE_LENGTH, fptrB) != NULL) {
    //fprintf(stderr,"L:%s!!!",line);
    ret += sscanf(line, "TimeLastAnalysis = %"PSYM,
                  &MetaData.TimeLastAnalysis);

    ret += sscanf(line, "dtAnalysis       = %"PSYM,
		  &MetaData.dtAnalysis);


    if (sscanf(line,  "AnalysisBaseName               = %s", dummy) == 1)
      MetaData.AnalysisName = dummy;


    //
    // Projection Parameters
    //

    if(sscanf(line, "ProjectionAnalysis   = %"ISYM,
	      &int_dummy) == 1){
      ret++;
      if(int_dummy)
	proj_analysis_on = TRUE;
    };

    if( sscanf(line,"ProjectionAxis[%"ISYM"]", &proj) != 0 ){
      ReadIntsFromString( strstr(line,"=") +1 , &(proj_params->ProjectionAxis[proj]), 1);
    }
    if( sscanf(line,"ProjectionLevel[%"ISYM"]", &proj) != 0 ){
      ReadIntsFromString( strstr(line,"=") +1 , &(proj_params->ProjectionLevel[proj]), 1);
    }


    if( sscanf(line,"ProjectionRightEdge[%"ISYM"]", &proj) != 0 ){
      ReadFloatFromString( strstr(line,"=") +1 , proj_params->ProjectionRightEdge[proj], 
			    MetaData.TopGridRank);
    }
    if( sscanf(line,"ProjectionLeftEdge[%"ISYM"]", &proj) != 0 ){
      ReadFloatFromString( strstr(line,"=") +1 , proj_params->ProjectionLeftEdge[proj], 
			    MetaData.TopGridRank);
    }

    ret += sscanf(line, "ProjectionOutputNumber = %"ISYM, &proj_params->output_number);

    if( sscanf(line,"ProjectionFieldList[%"ISYM"]", &proj) != 0 ){
      ReadTokensFromString( strstr(line,"=")+1, ProjectionFieldList_Char," ",MAX_PROJ_FIELDS);
      for( int q=0; q<MAX_PROJ_FIELDS; q++){
	if( ProjectionFieldList_Char[q] != NULL ){

	  this_id =  get_field_id(ProjectionFieldList_Char[q]);
	  if( this_id >= 0 ){
	    proj_params->fields[proj][q] = this_id;
	    proj_params->number_of_fields[proj]++;
	  }

	}
      }
    }

    if( sscanf(line,"ProjectionTitle[%"ISYM"]", &proj) != 0 ){
      ReadTokensFromString( strstr(line,"=")+1, &(proj_params->ProjectionTitle[proj])," ",1);
    }


    //
    // PDF Parameters
    //



    if(sscanf(line, "PDFAnalysis                      = %"ISYM,
	      &int_dummy) == 1){
      ret++;
      if(int_dummy)
	pdf_analysis_on = TRUE;
    };

    ret += sscanf(line,"PDFGeneratorMaximumLevel       = %"ISYM"\n", 
		  &(pdf_params->max_level) );
    ret += sscanf(line,"PDFGeneratorNumberOfPDFBins    = %"ISYM"\n", 
		  &(pdf_params->num_pdf_bins) );
    ret += sscanf(line,"PDFGeneratorCurrentOutputNumber = %"ISYM"\n", 
		  &(pdf_params->output_number) );    

    /* If the dummy char space was used, then make another. */
    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }


    
  }

  delete dummy;
  rewind(fptrB);

  if ( proj_analysis_on ){
    MetaData.ProjectionParams = proj_params;

  }else{
    delete proj_params;
  }

  if ( pdf_analysis_on ){
    MetaData.PDFGeneratorParams = pdf_params;
  }else{
    delete pdf_params;
  }



  if( Filename != NULL ){
    fclose(fptrB);
  }


}
