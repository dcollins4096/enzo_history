#ifndef PROB_DIST_FUNC_H
#define PROB_DIST_FUNC_H

#define PDF_DENSITY_LOG 0
#define PDF_TOTAL_ENERGY 1

class ProbabilityDistributionFunction {

 public:
  ProbabilityDistributionFunction( PDFGeneratorParameters *params,
		     int type, char *filename, 
		     char *description );

  ~ProbabilityDistributionFunction();
  ProbabilityDistributionFunction* NextFunction;

  //   void BinPDF(float density, 
  // 	      float total_energy, 
  // 	      float vx,
  // 	      float vy,
  // 	      float vz);

  void Reduce();
  void Write(FLOAT Time);

  // Making everything public so the grids can
  // access them.

  PDFGeneratorParameters *MyParameters;
  int func_type;

  char *description;
  char *filename;

  float pdf_min;
  float pdf_max;

  int num_pdf_bins;

  float log_bin_width;

  ealFloat *over;
  ealFloat *under;
  ealFloat *pdf;

  float *log_bin_centers();
  void WriteHeader(FILE *OutFile, FLOAT Time);
};

#endif /* PROB_DIST_FUNC_H */
