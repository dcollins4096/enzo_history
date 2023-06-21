#ifdef USE_FLEXIO
#include "WriteEvolution.h"

WriteEvolution::WriteEvolution(TopGridData &mdata):metadata(mdata){
  // Must have metadata intialized
  // before this class is instantiated
  if(metadata.commandline.jad.enabled){
    // setup the JAD stuff
    jadfiles.setSize(1);
    if(metadata.commandline.jad.filetype==Commandline::JadOpt::IEEE)
      jadfiles[0] = new IEEEIO(metadata.commandline.stars.filename,IObase::Write);
    else if(metadata.commandline.jad.filetype==Commandline::JadOpt::HDF4)
      jadfiles[0] = new HDFIO(metadata.commandline.stars.filename,IObase::Write);
    else{
      fprintf(stderr,"WriteEvolution: Unsupported filetype for JAD output\n");
      metadata.commandline.jad.enabled = 0; // disable Jad output
    }
  }
  if(metadata.commandline.stars.enabled){
    // setup the Star output stuff
     if(metadata.commandline.stars.filetype==Commandline::StarOpt::IEEE)
      starfile = new IEEEIO(metadata.commandline.stars.filename,IObase::Write);
    else if(metadata.commandline.stars.filetype==Commandline::StarOpt::HDF4)
      starfile = new HDFIO(metadata.commandline.stars.filename,IObase::Write);
    else{
      fprintf(stderr,"WriteEvolution: Unsupported filetype for STARS output\n");
      metadata.commandline.stars.enabled = 0; // disable Jad output
    }
  }
  if(metadata.commandline.dm.enabled){
    // setup the DM output stuff
    if(metadata.commandline.dm.filetype==Commandline::DMOpt::IEEE)
      dmfile = new IEEEIO(metadata.commandline.dm.filename,IObase::Write);
    else if(metadata.commandline.dm.filetype==Commandline::DMOpt::HDF4)
      dmfile = new HDFIO(metadata.commandline.dm.filename,IObase::Write);
    else{
      fprintf(stderr,"WriteEvolution: Unsupported filetype for DM output\n");
      metadata.commandline.dm.enabled = 0; // disable Jad output
    }
  }
}

WriteEvolution::~WriteEvolution(){
	for(int i=0;i<jadfiles.getSize();i++)
		if(jadfiles[i])
			delete jadfiles[i];
	if(starfile) delete starfile;
	if(dmfile)   delete dmfile;
}

#endif /* USE_FLEXIO */
