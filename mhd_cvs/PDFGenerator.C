#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <hdf5.h>
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
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "ealInt.h"
#include "ealFloat.h"
#include "ProbabilityDistributionFunction.h"
#include "PDFGenerator.h"

void my_exit(int status);

PDFGenerator::PDFGenerator( TopGridData *mdata, 
			  HierarchyEntry *topgrid): 
  AnalysisBaseClass( mdata, topgrid, 
		     mdata->PDFGeneratorParams->max_level,
		     NULL, NULL){

  this->MyParameters = mdata->PDFGeneratorParams;

  this->dens_func = new ProbabilityDistributionFunction( this->MyParameters, 
							 PDF_DENSITY_LOG, "log_density", "Log Density");
  
  this->en_func = new ProbabilityDistributionFunction( this->MyParameters, 
						       PDF_TOTAL_ENERGY, "energy", "TotalEnergy");
}

PDFGenerator::~PDFGenerator(void){

  /***** Delete PDFs *******/   
  delete this->dens_func;
  delete this->en_func;
}

void PDFGenerator::GeneratePDF(){
  LoopGeneratePDF(TopGrid, 0);

  this->dens_func->Reduce();
  this->en_func->Reduce();

  this->dens_func->Write(MetaData->Time);
  this->en_func->Write(MetaData->Time);
}

void PDFGenerator::LoopGeneratePDF( HierarchyEntry *Grid,
				    int current_level){

  while( Grid ){
    if(Grid->NextGridNextLevel)
      LoopGeneratePDF(Grid->NextGridNextLevel, current_level + 1);
    if(Grid->GridData){
      if ( Grid->GridData->IsInVolume( AnalysisRegionLeftEdge,
				       AnalysisRegionRightEdge)){
	
	if(Grid->GridData->ReturnProcessorNumber() == MyProcessorNumber){
	  if(current_level < MaximumAnalysisLevel){
	    FlagGridCells( Grid );
	  }else{
	    Grid->GridData->ClearFlaggingField();
	  }
#ifdef StandAloneTools
	  if(LoadGridDataAtStart == FALSE){
#ifndef UNIT_TESTING
	    if(Grid->GridData->ReadGrid(NULL, -1, FALSE, TRUE) == FAIL){
	      fprintf(stderr, "Failed to read in grid data, exiting. Grid information follows.\n");
	      my_exit(EXIT_FAILURE);
	    }
#else
            float dens = 1.25, en = 1049.0, inten = 3.25, vel[3] = {10., 0.0, 0.0};
	    Grid->GridData->InitializeUniformGrid(dens, en, inten, vel);
#endif
	  }
#endif // StandAloneTools

	  if(current_level <= MaximumAnalysisLevel){
	    Grid->GridData->BinPDF((void *)(this->dens_func), (void *)(this->en_func));  
	  }

	  Grid->GridData->DeleteFlaggingField();
#ifdef StandAloneTools
	  if(LoadGridDataAtStart == FALSE){
	    Grid->GridData->DeleteAllFields();
	  }
#endif // StandAloneTools
	}
      }
    }
    Grid = Grid->NextGridThisLevel;
  }
}

