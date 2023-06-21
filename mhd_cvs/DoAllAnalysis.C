#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <hdf5.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#ifdef ISO_GRAV
#include "GravityPotentialBoundary.h"
#endif
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
 
// HDF5 function prototypes
 
#include "extern_hdf5.h"
#include "ealInt.h"
#include "ealFloat.h"

#include "AnalysisBaseClass.h"
#include "Projection.h"
#include "ProbabilityDistributionFunction.h"
#include "PDFGenerator.h"

void DoAllAnalysis(HierarchyEntry *TopGrid, TopGridData &MetaData) {

  char buffer[10];
  char outputname[MAX_LINE_LENGTH];
  int i;

  // =====================================================================
  // Projection
  // =====================================================================
  
  if(MetaData.ProjectionParams){
    for( int p = 0; p<MAX_PROJECTIONS;p++){
      if( MetaData.ProjectionParams->ProjectionAxis[p] != INT_UNDEFINED ){
	Projection *proj;
	sprintf(outputname, "%s%4.4"ISYM, MetaData.AnalysisName, MetaData.ProjectionParams->output_number);

	if( debug && 0 == 1){
	  fprintf(stderr,"Projection[%d] Axis = %d\n",p, MetaData.ProjectionParams->ProjectionAxis[p]);
	  fprintf(stderr,"Projection[%d] Right = %f %f %f\n",
		  p, MetaData.ProjectionParams->ProjectionRightEdge[p][0],
		  MetaData.ProjectionParams->ProjectionRightEdge[p][1],
		  MetaData.ProjectionParams->ProjectionRightEdge[p][2]);
	  fprintf(stderr,"Projection[%d] Left = %f %f %f\n",
		  p, MetaData.ProjectionParams->ProjectionLeftEdge[p][0],
		  MetaData.ProjectionParams->ProjectionLeftEdge[p][1],
		  MetaData.ProjectionParams->ProjectionLeftEdge[p][2]);
	  
	  fprintf(stderr,"ProjectionFieldList[%d] = ");
	  for( int q=0;q<MAX_PROJ_FIELDS;q++)
	    if( MetaData.ProjectionParams->fields[p][q] != FieldUndefined )
	      fprintf(stderr," %d", MetaData.ProjectionParams->fields[p][q]);
	  fprintf(stderr,"\n");
	  
	  if( MetaData.ProjectionParams->ProjectionTitle[p] != NULL )
	    fprintf(stderr, "Projection[%d] Title = %s\n",p,
		    MetaData.ProjectionParams->ProjectionTitle[p] );
	  
	}


	proj = new Projection(&MetaData,TopGrid,
			      MetaData.ProjectionParams->number_of_fields[p],
			      MetaData.ProjectionParams->fields[p],
			      MetaData.ProjectionParams->ProjectionAxis[p],
			      MetaData.ProjectionParams->ProjectionLevel[p],
			      MetaData.ProjectionParams->ProjectionTitle[p],
			      MetaData.ProjectionParams->ProjectionLeftEdge[p],
			      MetaData.ProjectionParams->ProjectionRightEdge[p]);

	proj->Project();
	if(MyProcessorNumber == ROOT_PROCESSOR)
	  proj->Write(outputname);

	delete proj;

      }

	

    }
    /*
      Projection *proj;
      // full box projections
      for(i = 0; i < 3; i++){
      sprintf(outputname, "%s%4.4"ISYM".%"ISYM".proj", MetaData.AnalysisName,
      MetaData.ProjectionParams->output_number, i);
      proj = new Projection(&MetaData, TopGrid, 
      MetaData.ProjectionParams->number_of_fields,
      MetaData.ProjectionParams->fields, i, -1);
      proj->Project();
      if(MyProcessorNumber == ROOT_PROCESSOR)
      proj->Write(outputname);
      
      delete proj;
      }
    */
    /* Removed to loop the above business.
    // slices
    for(i = 0; i < 3; i++){
      FLOAT left[] = {0., 0., 0.};
      FLOAT right[] = {1., 1., 1.};
      right[i] = 1./(MetaData.TopGridDims[i]);

      sprintf(outputname, "%s%4.4"ISYM".%"ISYM".slice", MetaData.AnalysisName,
	      MetaData.ProjectionParams->output_number, i);
      
      proj = new Projection(&MetaData, TopGrid, 
			    MetaData.ProjectionParams->number_of_fields,
			    MetaData.ProjectionParams->fields, i, -1,
			    left, right);
      proj->Project();
      if(MyProcessorNumber == ROOT_PROCESSOR)
	proj->Write(outputname);
      
      delete proj;
    }
    */
    MetaData.ProjectionParams->output_number++;
  }


  // =====================================================================
  // Probability Distribution Functions
  // =====================================================================
  if(MetaData.PDFGeneratorParams){
    strcpy(MetaData.PDFGeneratorParams->OutputName, MetaData.AnalysisName);
    sprintf(buffer, "%4.4"ISYM, MetaData.PDFGeneratorParams->output_number);
    strcat(MetaData.PDFGeneratorParams->OutputName, buffer);
  
    MetaData.PDFGeneratorParams->output_number++;

    PDFGenerator pdfgen(&MetaData, TopGrid);
    pdfgen.GeneratePDF();
  }

}
