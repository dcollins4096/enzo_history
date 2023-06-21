#ifndef __WriteEvolution_h_
#define __WriteEvolution_h_
#ifdef FAIL /* fix inconsistency between HDF version of FAIL vs. ENZO def */
#undef FAIL
#endif
#include <FlexArrayTmpl.H>
#include <IEEEIO.hh>
#include <HDFIO.hh>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "TopGridData.h"

class WriteEvolution {
  TopGridData &metadata;
  FlexArray<IObase*> jadfiles;
  IObase *starfile,*dmfile;
public:
  WriteEvolution(TopGridData &mdata);
  ~WriteEvolution();
};

#ifdef FAIL /* restore FAIL definition to what enzo expects. */
#undef FAIL
#endif
#define FAIL 0

#endif
