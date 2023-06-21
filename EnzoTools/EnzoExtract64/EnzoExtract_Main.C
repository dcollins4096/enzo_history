/*--------------------------------------------------------------------
 *
 *                               EnzoExtract 
 *
 *
 *  EnzoExtract is a tool that extracts rectangular solid sections of
 *  Enzo hdf 5 datasets and outputs them as single monolithic files.
 *  This is useful for vizualization and some sorts of data analysis
 *  (spectrum generation, etc.).  It is intended to fully replace the
 *  extraction feature built into the code, ie, anything that is done
 *  with the enzo -x command.
 *  
 *
 *  Compilation Instructions:
 *    compile with gnu make (gmake on some systems):  
 *    just type "make".
 *    you may have to tweak the makefile in order to find the hdf
 *    libraries.  Also, you have to change the top line of the 
 *    makefile to reflect the site you're at - this will hopefully
 *    change soon.
 * 
 *  Usage Instructions:
 *    After compilation type 'enzoextract -h' or 'enzoextract -help'
 *    (without the quotes), or see the enzo user's manual or cookbook.
 *    You must have the entire enzo hierarchy (all grid files, the
 *    hierarchy file and the parameter file) in the directory that
 *    you run enzoextract in (though of course the executable doesn't
 *    have to be in that directory).  The output files also appear in
 *    that directory.
 *
 *  Written By:  
 *      Brian O'Shea, bwoshea@cosmos.ucsd.edu
 *      July 30, 2003
 *
 *  Revision history:
 *
 *      June 2004:  updated EnzoExtract so that all grid and cell locations
 *        are calculated using 64 bits of precision, in order to work well
 *        for high-resolution simulations.  (BWO) 
 *
 *      June 12 2006:   added Alexei Kritsuk's modifications so that we can
 *        use the Harkness packed IO format and also have arbitrary
 *        RefineBy ratios.  Note that the packed IO stuff is controlled
 *        in EnzoExtract.h using the PACK_AMR define.  I also cleaned up
 *        some of the misc. printf statements that make this code so darned
 *        noisy.  (BWO)
 *
 *      June 14, 2006:  made a few modifications - now the code can handle
 *        fully 64-bit float IO (controlled by a define statement in the
 *        EnzoExtract.h header file), and a few remaining incidences of "float"
 *        in position information, etc. were removed.  Some mods to the 
 *        HelpMe routine were also added.  A define statement now controls
 *        the hierarchy file format read in - is it Rick's modified hierarchy
 *        file, or the standard file?  (BWO)
 *
 *  Notes for modifications:
 *    The routine CreateExtractionArrays() in EnzoExtract_Misc.C is the
 *    file that must be modified to add more extraction arrays.
 *
 *--------------------------------------------------------------------*/
#define ORIGINAL_DEFINE
#include "EnzoExtract.h"
#undef ORIGINAL_DEFINE

int main(int argc, char *argv[]){
  if(DEBUG) fprintf(stderr,"main:  entering\n");
  int returnflag,i;

  if(argc < 2){
    HelpMe();
    exit(SUCCESS);
  }

  // read in inputs
  ParseInputs(argc,argv);

  // read in the parameter file for stuff we need 
  ReadParameterFile(parameterfile);

  // check boundary values for sanity - report errors
  // this needs stuff from ReadParameterFile, so don't move it!
  CheckBoundaryValues();

  // get hierarchy file name
  int inlen = (int) strlen(parameterfile);
  char *hierfile = new char[inlen+6];
  strcpy(hierfile,parameterfile);
  hierfile=strcat(hierfile,".hierarchy");

  // read hierarchy file and get grid information
  GetGridInfo(hierfile);

  // generate extraction arrays
  CreateExtractionArrays();

  // extract desired quantities to grid  -- arrays in column-order
  for(i=0; i<numberofgrids; i++) Extract(i);

  // if the input flag is switched on, change extraction arrays to row-order
  if(outputrowmajor) SwitchArraysToRowOrder();

  // write extraction arrays to file
  WriteExtractionArrays();

  // erase extraction arrays
  DeleteExtractionArrays();
 
  if(DEBUG) fprintf(stderr,"main:  exiting\n");
  exit(SUCCESS);
}
