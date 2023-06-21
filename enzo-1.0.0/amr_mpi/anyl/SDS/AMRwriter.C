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
// AMRwriter.cc
#include <stdio.h>
#include <stdlib.h>
#include "AMRwriter.hh"

void AMRwriter::writeAMRattributes(){
  double origin[5],delta[5],ext[5],currenttime;
  // OK, now we get to do the attribs thang...
  // first compute the coorect origin and delta for this level.
  //int i;
  int level_timestep,persistence;
  
  /*****
	Must compute the floating point time from the integer timestep.
	The timestep is based on the finest computational level.
	Should compute from level_timestep for less numerical error
  *****/
  currenttime = (double)currentstep*basetimestep/(double)(levels[levels.getSize()-1].trefine);
  // compute persistence
  
  persistence = levels[levels.getSize()-currentlevel-1].trefine;
  level_timestep = currentstep/persistence;
  /******
	 Generate delta, max_ext, and origin from the integer origin and 
	 spatial refinement.
  ******/
  for(int i=0;i<drank;i++){
    delta[i]=ddelta[i]/(double)levels[currentlevel].srefine[i];
    // is the origin relative to the parent or is it relative to
    // the overall grid structure?  I'm choosing relative to the 
    // global origin.
    origin[i] = dorigin[i] + (double)(iorigin[i])/(double)levels[currentlevel].prefine[i];
    // compute from min-ext...
    // y-know... there is more point density near 0 in the float
    // representation... so having sizes instead of extents would
    // be smarter.... But viz systems generally use the extent... :(
		
    // this is VERY inaccurate... Have to think about this more...
    // Its OK for small subgrids, but inaccuracy is unbearable past 
    // 128 elements.  Need to compute relative to basegrid extents!, 
    // but that isn't possible without dims for the basegrid!
    ext[i] = origin[i] + (double)(ddims[i]) * delta[i];
  }
  // IObase::setOrigin(origin);// this will set the IObase:: basegrid origin :(
  // IObase::setDelta(delta);
  // Storage is in order of liklihood of that param being accessed.
  // Less-frequently accessed params go to the end.
  file.writeAttribute("level",IObase::Int32,1,&currentlevel);
  file.writeAttribute("origin",IObase::Float64,drank,origin);
  file.writeAttribute("delta",IObase::Float64,drank,delta);
  file.writeAttribute("min_ext",IObase::Float64,drank,origin);
  file.writeAttribute("max_ext",IObase::Float64,drank,ext);
  file.writeAttribute("time",IObase::Float64,1,&currenttime);

  // writeBounds();
  file.writeAttribute("timestep",IObase::Int32,1,&currentstep);
  file.writeAttribute("level_timestep",IObase::Int32,1,&level_timestep);
  file.writeAttribute("persistence",IObase::Int32,1,&persistence);
  file.writeAttribute("time_refinement",IObase::Int32,1,&(levels[currentlevel].trefine));

  // need timestep vs. leveltimestep (must define time-is-local vs. time-is-global)
  file.writeAttribute("spatial_refinement",IObase::Int32,drank,levels[currentlevel].srefine);
  file.writeAttribute("grid_placement_refinement",IObase::Int32,1,&(levels[currentlevel].prefine));
  file.writeAttribute("iorigin",IObase::Int32,drank,iorigin);
  // I guess thats it for integer parameters.
  // These params are only for informational purposes
  // Most viz systems will probably find the double precision 
  // parameters most suitable
}
	
AMRwriter::AMRwriter(IObase &descriptor):
  Writer(descriptor),basetimestep(1.0),levels(1),
  currentlevel(0),currentstep(0){
    for(int i=0;i<5;i++){
      levels[0].srefine[i]=levels[0].prefine[i]=1;
      iorigin[i]=0;
    }
    levels[0].trefine=0;
}

AMRwriter::~AMRwriter(){
  // destroy whatever....
  // ... nothing to destroy now (went all static)
}

void AMRwriter::setToplevelParameters(double *origin,
				      double *delta, double timestep){
  Writer::setOrigin(origin);
  Writer::setDelta(delta);
  /*
    for(int i=0;i<rank;i++){
    gorigin[i]=origin[i];
    gdelta[i]=delta[i];
    }*/
  basetimestep=timestep;
}

// Could also make placement refinement FINESTLEVEL=-1000 and then you could
// say things like gridplacement = FINESTLEVEL -1 or FINESTLEVEL +1
// or maybe even FINESTLEVELREFINEMENT*2 or FINESTLEVELREFINEMENT/2
// to get the correct value?  Would require float support?
// This would generate a lot of IF statements to monitor state though...
void AMRwriter::setLevelParameters(int level,int timerefinement,
				   int *spatialrefinement,
				   int *gridplacementrefinement){
  if(level>=levels.getSize()) levels.setSize(level+1);
  levels[level].trefine=timerefinement;
  for(int i=0;i<drank;i++) {
    (levels[level]).srefine[i]=spatialrefinement[i];
    if(gridplacementrefinement)
      (levels[level]).prefine[i]=gridplacementrefinement[i];
    else
      (levels[level]).prefine[i]=spatialrefinement[i];
  }
}

void AMRwriter::setLevelParameters(int level,int timerefinement,
				   int spatialrefinement,int gridplacementrefinement){
  if(level>=levels.getSize()) levels.setSize(level+1);
  levels[level].trefine=timerefinement;
  printf("\tsetLevelParameters trefine=%u\n",levels[level].trefine);
  for(int i=0;i<drank;i++) {
    (levels[level]).srefine[i]=spatialrefinement;
    if(gridplacementrefinement)
      (levels[level]).prefine[i]=gridplacementrefinement;
    else
      (levels[level]).prefine[i]=spatialrefinement;
  }
}// using scalar values)



void AMRwriter::setOrigin(int *origin){
  // now this is the integer origin
  for(int i=0;i<drank;i++) iorigin[i]=origin[i];
}

void AMRwriter::setGrid(int level,int timestep,int *origin,int *dims){
  setLevel(level);
  setDims(dims);
  setTime(timestep);
  setOrigin(origin);
}

void AMRwriter::write(void *data){
  // write data
  // then write AMR attributes (which will include
  // the regular complement of attribs)
  // IObase::write(data); // should call virtual writebounds?
  file.write(dtypeID,drank,ddims,data);
  writeAMRattributes();
}

void FrameworkAMRwriter::setParameters(int rank,
				    double *origin,
				    double *delta, 
				    double timestep,
				    int interlevelRefinementFactor,
				    int numLevels){
  register int i,ref;
  drank = rank;
  nlevels=numLevels;
  refinement = interlevelRefinementFactor; // class member for FrameworkAMRwriter
  for(i=0,maxrefinement=1;i<nlevels-1;i++) maxrefinement*=interlevelRefinementFactor;
  AMRwriter::setToplevelParameters(origin,delta,timestep);
  AMRwriter::setLevelParameters(0,1,1,maxrefinement);
  for(i=1,ref=interlevelRefinementFactor;i<nlevels;i++,ref*=refinement)
     AMRwriter::setLevelParameters(i,ref,ref,maxrefinement);
}
void FrameworkAMRwriter::setParameters(IObase::DataType dt,
				       int rank,
				       double *origin,
				       double *delta, 
				       double timestep,
				       int interlevelRefinementFactor,
				       int numLevels){
  register int i,ref;
  drank = rank;
  nlevels=numLevels;
  Writer::setType(dt);
  refinement = interlevelRefinementFactor; // class member for FrameworkAMRwriter 
  for(i=1,maxrefinement=1;i<nlevels;i++) maxrefinement*=interlevelRefinementFactor;
  printf("FrameworkAMRwriter maxrefinement=%u\n",maxrefinement);
  AMRwriter::setToplevelParameters(origin,delta,timestep);
  AMRwriter::setLevelParameters(0,1,1,maxrefinement);
  for(i=1,ref=interlevelRefinementFactor;i<nlevels;i++,ref*=interlevelRefinementFactor){
     AMRwriter::setLevelParameters(i,ref,ref,maxrefinement);
     printf("setLevelParameters[%u](level=%u,timeref=%u,Spaceref=%u,gridref=%u)\n",
  	    i,i,ref,ref,maxrefinement);
  }
}

//===========C Interface======================================
// should have an RTTI interface and inherit everything from
// IOobject which contains the RTTI isOfType() information.
// isOfType() should propagate recursively to determin type info
// match.  must grab typeID from floating point pool.  And then
// we need a static initializer for everything.

// How does performer/inventor do this?


AMRFile AMRbeginFile(IOFile *descriptor){
  IObase *io = (IObase*)descriptor;
  return (AMRFile)(new AMRwriter(*io));
}

void AMRendFile(AMRFile afile){
  AMRwriter *w = (AMRwriter*)afile;
  delete w;
}

void AMRsetRank(AMRFile afile,int rank){
  AMRwriter *w = (AMRwriter*)afile;
  w->setRank(rank);
}

void AMRsetType(AMRFile afile,int numbertype){
  AMRwriter *w = (AMRwriter*)afile;
  w->setType(IObase::Int2DataType(numbertype));
}

void AMRsetToplevelParameters(AMRFile afile,double *origin,
			      double *delta, double timestep){
  AMRwriter *w = (AMRwriter*)afile;
  w->setToplevelParameters(origin,delta,timestep);
}


void AMRsetLevelParameters(AMRFile afile,int level,
			   int timerefinement,
			   int *spatialrefinement,
			   int *gridplacementrefinement){
  AMRwriter *w = (AMRwriter*)afile;
  w->setLevelParameters(level,timerefinement,spatialrefinement,gridplacementrefinement);
}


/* using scalar values) */
void AMRsetScalarLevelParameters(AMRFile afile,int level,
			   int timerefinement,
			   int spatialrefinement,
			   int gridplacementrefinement){
  AMRwriter *w = (AMRwriter*)afile;
  w->setLevelParameters(level,timerefinement,spatialrefinement,gridplacementrefinement);
}


void AMRsetDims(AMRFile afile,int *dims){
  AMRwriter *w = (AMRwriter*)afile;
  w->setDims(dims);
}


void AMRsetRankDims(AMRFile afile,int rank, int *dims){
  AMRwriter *w = (AMRwriter*)afile;
  w->setDims(rank,dims);
}


void AMRsetOrigin(AMRFile afile,int *origin){
  AMRwriter *w = (AMRwriter*)afile;
  w->setOrigin(origin);
}


void AMRwrite(AMRFile afile,void *data){
  AMRwriter *w = (AMRwriter*)afile;
  w->write(data);
}

fAMRFile fAMRbeginFile(IOFile *descriptor){
  IObase *io = (IObase*)descriptor;
  return (fAMRFile)(new FrameworkAMRwriter(*io));
}

void fAMRendFile(fAMRFile afile){
  FrameworkAMRwriter *w = (FrameworkAMRwriter*)afile;
  delete w;
}
void fAMRsetLevelParameters(fAMRFile afile,
				   int datatype,
				   int rank,
				   double *origin,
				   double *delta,
				   double timestep,
				   int interlevelRefinementRatio,
			    int nlevels){
  FrameworkAMRwriter *w = (FrameworkAMRwriter*)afile;
  IObase::DataType dt = IObase::Int2DataType(datatype);
  w->setParameters(dt,rank,origin,delta,timestep,interlevelRefinementRatio,nlevels);
}
void fAMRwrite(fAMRFile afile,
	       int level,
	       int globaltimestep,
	       int *origin,
	       int *dims,
	       void *data){
  FrameworkAMRwriter *w = (FrameworkAMRwriter*)afile;
  w->write(level,globaltimestep,origin,dims,data);
}
