#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <hdf5.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_GNUPLOT
#include "gnuplot_i.h"
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

ProbabilityDistributionFunction::ProbabilityDistributionFunction( PDFGeneratorParameters *params,
								  int type, char *fname, 
								  char *descript ){
  
  this->MyParameters = params;
  this->num_pdf_bins = MyParameters->num_pdf_bins;
  this->func_type = type;

  NextFunction = NULL;

  filename = new char[MAX_LINE_LENGTH];
  strcpy(this->filename, fname);
  description = new char[MAX_LINE_LENGTH];
  strcpy(this->description, descript);

  /* WARN - Hard coded pdf_min and pdf_max */

  pdf_min = 1E-4;

  switch(func_type){
  case PDF_DENSITY_LOG: // log10 density
    pdf_max = 1.0E4;
    break;
  case PDF_TOTAL_ENERGY: // log10 total energy
    pdf_min = 1.0E2;
    pdf_max = 1.0E4;
    break;
  }

  this->over = new ealFloat(1);
  this->under = new ealFloat(1);
  this->pdf = new ealFloat(this->num_pdf_bins);

  log_bin_width = log10((pdf_max)/pdf_min)/((float)num_pdf_bins);
}

ProbabilityDistributionFunction::~ProbabilityDistributionFunction(){
  delete[] this->filename;
  delete[] this->description;

  delete over;
  delete under;
  delete pdf;
}

void ProbabilityDistributionFunction::Reduce( ){

  float v_actual = 1.0;
  int i;


  // WARN WARN WARN - currently uses domain vs analysis region
  for( i = 0; i < MAX_DIMENSION; i++ )
    v_actual *= (DomainRightEdge[i] - DomainLeftEdge[i]);

  under->ReduceSum();
  over->ReduceSum();
  pdf->ReduceSum();

  v_actual -= (under->Array[0] + over->Array[0]);

  // normalize by actual volume measured
  if( MyProcessorNumber == ROOT_PROCESSOR ){
    if(v_actual > 0){
      *(pdf)/= v_actual;
    }
  }
}

void ProbabilityDistributionFunction::Write(FLOAT Time){

  int i, j;

  float *centers = log_bin_centers( );

  char *pdf_out_name =  new char[MAX_LINE_LENGTH];

  strcpy(pdf_out_name, MyParameters->OutputName);  
  strcat(pdf_out_name, ".PDF." );
  strcat(pdf_out_name, filename);

  FILE *pdf_out_file = fopen( pdf_out_name, "w" );

  fprintf( pdf_out_file,     "# Observable          = %s\n",
	   description );

  WriteHeader( pdf_out_file, Time );

  fprintf( pdf_out_file, "# PDFMinimum           = %le\n", pdf_min);  
  fprintf( pdf_out_file, "# PDFMaximum           = %le\n", pdf_max);  
  fprintf( pdf_out_file, "# PDFLogBinWidth       = %le\n", log_bin_width);

  fprintf( pdf_out_file, "\n# Outliers (Volume Weighted)\n");  
  fprintf( pdf_out_file, "# %12.8e < PDFMinimum\n", under->Array[0]);
  fprintf( pdf_out_file, "# %12.8e > PDFMaximum\n", over->Array[0]);

  fprintf( pdf_out_file, "\n" );
  fprintf( pdf_out_file, "# bin, center, pdf\n" );
   
  for( i = 0; i < num_pdf_bins; i++ ){
    fprintf( pdf_out_file, "%"ISYM", %12.8e, %12.8e\n", i+1, centers[i], 
	       pdf->Array[i] );    
  }
   
  fclose( pdf_out_file );


#ifdef HAVE_GNUPLOT
  if(MyParameters->gplot){
   
    char eps_out_name[MAX_LINE_LENGTH];
    char *plot_cmd = new char[MAX_LINE_LENGTH*10];

    sprintf(plot_cmd, "plot '%s' using 2:3 with lines ", pdf_out_name);

    strcpy(eps_out_name, pdf_out_name);
    strcat(eps_out_name, ".eps" );

    gnuplot_ctrl * g = gnuplot_init();
    if(g){
      gnuplot_cmd(g, "set terminal postscript color");
      
      gnuplot_cmd(g, "set logscale x");
      gnuplot_cmd(g, "set logscale y");
      gnuplot_cmd(g, "set output '%s'", eps_out_name);
      gnuplot_cmd(g, "set title '%s'", description);
      gnuplot_cmd(g, plot_cmd);  
      
      gnuplot_close(g);
    }
    delete [] plot_cmd;
  }
#endif // HAVE_GNUPLOT


  delete[] centers;
  //delete[] out_name;
  delete[] pdf_out_name;
   
}

void ProbabilityDistributionFunction::WriteHeader( FILE *OutFile, FLOAT Time ){

  fprintf( OutFile, "# NumberOfProcessors  = %"ISYM"\n", NumberOfProcessors);
  fprintf( OutFile, "# NumberOfPDFBins     = %"ISYM"\n", MyParameters->num_pdf_bins);
  fprintf( OutFile, "# Time                = %"GOUTSYM"\n", Time);
  fprintf( OutFile, "\n");

  return;
}

float *ProbabilityDistributionFunction::log_bin_centers(  ){

  int i;
  float *centers = new float[num_pdf_bins];
  float log_bin_width = log10((pdf_max)/(pdf_min))/((float)num_pdf_bins);
  float scale = pow(10.0, log_bin_width);

  centers[0] = (pdf_min)*pow(10.0, log_bin_width/2.0);

  for( i=1; i<num_pdf_bins; i++ )
    centers[i] = centers[i-1]*scale;

  return centers;
}
