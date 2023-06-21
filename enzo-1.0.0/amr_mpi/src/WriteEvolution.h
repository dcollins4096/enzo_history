/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
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
