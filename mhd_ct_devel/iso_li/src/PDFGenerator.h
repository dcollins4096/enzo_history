#ifndef PDFGENERATOR_H
#define PDFGENERATOR_H
/***********************************************************************
/  
/ PROBABILITY DISTRIBUTION FUNCTION GENERATOR CLASS, 
/ SUBCLASS OF ENZO ANALYSIS BASE CLASS
/
************************************************************************/
#include "AnalysisBaseClass.h"

class PDFGenerator : public AnalysisBaseClass {
  
 public:
  PDFGenerator( TopGridData *mdata, 
		HierarchyEntry *topgrid );

  ~PDFGenerator(void);
  void GeneratePDF();
  void PrintParameters();

#ifdef UNIT_TESTING
 public:
#else
 private:
#endif

  void LoopGeneratePDF(HierarchyEntry *Grid,
		       int current_level);
  
  PDFGeneratorParameters *MyParameters;
  ProbabilityDistributionFunction *dens_func;
  ProbabilityDistributionFunction *en_func;
};

#endif
