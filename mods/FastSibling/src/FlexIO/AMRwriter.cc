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
  // compute from min-ext...
  // y-know... there is more point density near 0 in the float
  // representation... so having sizes instead of extents would
  // be smarter.... But viz systems generally use the extent... :(		
  // this is VERY inaccurate... Have to think about this more...
  // Its OK for small subgrids, but inaccuracy is unbearable past 
  // 128 elements.  Need to compute relative to basegrid extents!, 
  // but that isn't possible without dims for the basegrid!
  // IObase::setOrigin(origin);// this will set the IObase:: basegrid origin :(
  // IObase::setDelta(delta);
  // Storage is in order of liklihood of that param being accessed.
  // Less-frequently accessed params go to the end.
  for(int i=0;i<drank;i++) {
    currentdelta[i] = ddelta[i]/(levels[currentlevel]).srefine[i];
    ext[i] = currentorigin[i] + 
      (double)((ddims[i]+1) * currentdelta[i]);
  }
  
  file.writeAttribute("level",IObase::Int32,1,&currentlevel);
  file.writeAttribute("origin",IObase::Float64,drank,currentorigin);
  file.writeAttribute("delta",IObase::Float64,drank,currentdelta);
  file.writeAttribute("min_ext",IObase::Float64,drank,currentorigin);
  file.writeAttribute("max_ext",IObase::Float64,drank,ext);
  file.writeAttribute("time",IObase::Float64,1,&currenttime);

  // writeBounds();
  file.writeAttribute("timestep",IObase::Int32,1,&currentstep);
  file.writeAttribute("level_timestep",IObase::Int32,1,&level_timestep);
  file.writeAttribute("persistence",IObase::Int32,1,&persistence);
  file.writeAttribute("time_refinement",IObase::Int32,1,&(levels[currentlevel].trefine));

  // need timestep vs. leveltimestep (must define time-is-local vs. time-is-global)
  file.writeAttribute("spatial_refinement",IObase::Int32,drank,levels[currentlevel].srefine);
  file.writeAttribute("grid_placement_refinement",IObase::Int32,drank,levels[currentlevel].prefine);
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

void AMRwriter::setTopLevelParameters(int rank,double *origin,
				      double *delta, double timestep,int maxdepth){
  setRank(rank);
  Writer::setOrigin(origin);
  Writer::setDelta(delta);
  levels.setSize(maxdepth);
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
void AMRwriter::setRefinement(int timerefinement,
			      int spatialrefinement,
			      int gridplacementrefinement){
  int tref=1,sref=1,gref=gridplacementrefinement;
  int maxdepth=levels.getSize()+1;
  // if gref==-1, then we need to multiply by maxdepth
  for(int level=0;level<maxdepth;level++){
    setLevelRefinement(level,tref,sref,gref);
    tref*=timerefinement;
    sref*=spatialrefinement;
    gref*=spatialrefinement;
  }
  if(gridplacementrefinement<0){
    gref = -(levels[maxdepth-1]).prefine[0];
    for(int level=0;level<maxdepth;level++) 
      for(int i=0;i<3;i++) (levels[level]).prefine[i] = gref;
  }
}
void AMRwriter::setRefinement(int timerefinement,
			      int *spatialrefinement,
			      int *gridplacementrefinement){
  int maxdepth = 1+levels.getSize();
  int tref = 1;
  int *sref = new int[drank];
  int *gref = new int[drank];
  for(int i=0;i<drank;i++) 
  {     sref[i]=gref[i]=1;
        if(gridplacementrefinement) gref[i]+=gridplacementrefinement[i];
  }
  for(int level=0;level<maxdepth;level++){
    setLevelRefinement(level,tref,sref,gref);
    for(int i=0;i<drank;i++) {
      sref[i]*=spatialrefinement[i];
      gref[i]*=spatialrefinement[i];
    }
  }
  int isnegative=0,negmask[3]={0,0,0};
  {for(int i=0;i<drank && i<3 ;i++)
    if(gridplacementrefinement[i]<0)
      negmask[i]=isnegative=1;
  }
  if(isnegative){ // set max depth for all levels
    for(int i=0;i<drank;i++) gref[i]=-(levels[maxdepth-1]).prefine[i];

    for(int level=0;level<maxdepth;level++)
    {
      for(int i=0;i<drank;i++)
	if(negmask[i])
	  (levels[level]).prefine[i]=gref[i];
    }
  }
  delete sref;
  delete gref;
}

void AMRwriter::setLevelRefinement(int level,int timerefinement,
				   int *spatialrefinement,
				   int *gridplacementrefinement){
  if(level>=levels.getSize()) levels.setSize(level+1);
  levels[level].trefine=timerefinement;
  for(int i=0;i<drank;i++) {
    (levels[level]).srefine[i]=spatialrefinement[i];
    if(gridplacementrefinement>0)
      (levels[level]).prefine[i]=gridplacementrefinement[i];
    else
      (levels[level]).prefine[i]=spatialrefinement[i];
  }
}

void AMRwriter::setLevelRefinement(int level,int timerefinement,
				   int spatialrefinement,
				   int gridplacementrefinement){
  if(level>=levels.getSize()) levels.setSize(level+1);
  levels[level].trefine=timerefinement;
  for(int i=0;i<drank;i++) {
    (levels[level]).srefine[i]=spatialrefinement;
    if(gridplacementrefinement>0)
      (levels[level]).prefine[i]=gridplacementrefinement;
    else
      (levels[level]).prefine[i]=spatialrefinement;
  }
}

void AMRwriter::setOrigin(int *origin){
  for(int i=0;i<drank;i++) {
    double dx = ddelta[i]/(levels[currentlevel]).prefine[i];
    iorigin[i]=origin[i]; // integer origin (logical origin)
    currentorigin[i]=(origin[i]*dx+dorigin[i]); // float origin (real origin)
  }
}

void AMRwriter::setOrigin(float *origin){
  for(int i=0;i<drank;i++) {
    double dx = ddelta[i]/(levels[currentlevel]).prefine[i];
    iorigin[i]=(origin[i]-dorigin[i])/dx;
    // set floating point origin
    currentorigin[i]=origin[i];
  }
}

void AMRwriter::setOrigin(double *origin){
  for(int i=0;i<drank;i++) {
    double dx = ddelta[i]/(levels[currentlevel]).prefine[i];
    iorigin[i]=(origin[i]-dorigin[i])/dx;
    // set floating point origin
    currentorigin[i]=origin[i];
  }
}

void AMRwriter::write(void *data){
  // write data
  // then write AMR attributes (which will include
  // the regular complement of attribs)
  // IObase::write(data); // should call virtual writebounds?
  file.write(dtypeID,drank,ddims,data);
  writeAMRattributes();
}

void FrameworkAMRwriter::init(int rank,
			      double *origin,
			      double *delta, 
			      double timestep,
			      int interlevelRefinementFactor,
			      int numLevels){
  register int i,ref;
  nlevels=numLevels;
  refinement = interlevelRefinementFactor; // class member for FrameworkAMRwriter
  for(i=0,maxrefinement=1;i<nlevels-1;i++) maxrefinement*=interlevelRefinementFactor;
  AMRwriter::setTopLevelParameters(rank,origin,delta,timestep,numLevels);
  AMRwriter::setRefinement(numLevels,interlevelRefinementFactor,interlevelRefinementFactor);
}
void FrameworkAMRwriter::init(IObase::DataType dt,
			      int rank,
			      double *origin,
			      double *delta, 
			      double timestep,
			      int interlevelRefinementFactor,
			      int numLevels){
  Writer::setType(dt);
  init(rank,origin,delta,timestep,interlevelRefinementFactor,numLevels);
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

void AMRsetType(AMRFile afile,int numbertype){
  AMRwriter *w = (AMRwriter*)afile;
  w->setType(IObase::Int2DataType(numbertype));
}

void AMRsetToplevelParameters(AMRFile afile,int rank, double *origin,
			      double *delta, double timestep,int maxdepth){
  AMRwriter *w = (AMRwriter*)afile;
  w->setTopLevelParameters(rank,origin,delta,timestep,maxdepth);
}

void AMRsetRefinement(AMRFile afile,
		      int timerefinement,
		      int *spatialrefinement,
		      int *gridplacementrefinement){
  AMRwriter *w = (AMRwriter*)afile;
  w->setRefinement(timerefinement,spatialrefinement,gridplacementrefinement);
}
void AMRsetScalarRefinement(AMRFile afile,
			    int timerefinement,
			    int spatialrefinement,
			    int gridplacementrefinement){
  AMRwriter *w = (AMRwriter*)afile;
  w->setRefinement(timerefinement,spatialrefinement,gridplacementrefinement);
}
void AMRsetLevelRefinement(AMRFile afile,int level,
			   int timerefinement,
			   int *spatialrefinement,
			   int *gridplacementrefinement){
  AMRwriter *w = (AMRwriter*)afile;
  w->setLevelRefinement(level,timerefinement,spatialrefinement,gridplacementrefinement);
}


/* using scalar values) */
void AMRsetScalarLevelRefinement(AMRFile afile,int level,
				 int timerefinement,
				 int spatialrefinement,
				 int gridplacementrefinement){
  AMRwriter *w = (AMRwriter*)afile;
  w->setLevelRefinement(level,timerefinement,spatialrefinement,gridplacementrefinement);
}

void AMRlevel(AMRFile afile,int level){
  AMRwriter *w = (AMRwriter*)afile;
  w->setLevel(level);
}

void AMRtime(AMRFile afile,int time){
  AMRwriter *w = (AMRwriter*)afile;
  w->setTime(time);
}

void AMRincrementTime(AMRFile afile){
  AMRwriter *w = (AMRwriter*)afile;
  w->incrementTime();
}

void AMRwrite(AMRFile afile,int *origin,int *dims,void *data){
  AMRwriter *w = (AMRwriter*)afile;
  w->write(origin,dims,data);
}

void AMRwriteFloat(AMRFile afile,float *origin,int *dims,void *data){
  AMRwriter *w = (AMRwriter*)afile;
  w->write(origin,dims,data);
}

void AMRwriteDouble(AMRFile afile,double *origin,int *dims,void *data){
  AMRwriter *w = (AMRwriter*)afile;
  w->write(origin,dims,data);
}

fAMRFile fAMRbeginFile(IOFile *descriptor){
  IObase *io = (IObase*)descriptor;
  return (fAMRFile)(new FrameworkAMRwriter(*io));
}

void fAMRendFile(fAMRFile afile){
  FrameworkAMRwriter *w = (FrameworkAMRwriter*)afile;
  delete w;
}

void fAMRinit(fAMRFile afile,
	      int datatype,
	      int rank,
	      double *origin,
	      double *delta,
	      double timestep,
	      int interlevelRefinementRatio,
	      int nlevels){
  FrameworkAMRwriter *w = (FrameworkAMRwriter*)afile;
  IObase::DataType dt = IObase::Int2DataType(datatype);
  w->init(dt,rank,origin,delta,timestep,interlevelRefinementRatio,nlevels);
}
void fAMRsetParameters(fAMRFile afile,
		       int datatype,
		       int rank,
		       double *origin,
		       double *delta,
		       double timestep,
		       int interlevelRefinementRatio,
		       int nlevels){
  fAMRinit(afile,datatype,rank,origin,delta,timestep,interlevelRefinementRatio,nlevels);
}
inline void fAMRsetLevelParameters(fAMRFile afile,
				   int datatype,
				   int rank,
				   double *origin,
				   double *delta,
				   double timestep,
				   int interlevelRefinementRatio,
				   int nlevels){
  fAMRinit(afile,datatype,rank,origin,delta,timestep,interlevelRefinementRatio,nlevels);
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

//===========F77/F90 Interface==================================

Long8 f_amr_begin(Long8 *descriptor){
  IObase *io = (IObase*)descriptor;
  return (Long8)(new AMRwriter(*io));
}

int f_amr_end(Long8 *afile){
  AMRwriter *w = (AMRwriter*)afile;
  delete w;
  return 1;
}

int f_amr_settype(Long8 *afile,int *numbertype){
  AMRwriter *w = (AMRwriter*)afile;
  w->setType(IObase::Int2DataType(*numbertype));
  return 1;
}

int f_amr_setparams(Long8 *afile,int *rank,double *origin,
		    double *delta, double *timestep,int *maxdepth){
  AMRwriter *w = (AMRwriter*)afile;
  w->setTopLevelParameters(*rank,origin,delta,*timestep,*maxdepth);
  return 1;
}

int f_amr_setref(Long8 *afile,int *veclen,
		 int *timerefinement,
		 int *spatialrefinement,
		 int *gridplacementrefinement){
  AMRwriter *w = (AMRwriter*)afile;
  if(*veclen==1)
    w->setRefinement(*timerefinement,
		     *spatialrefinement,*gridplacementrefinement);
  else
    w->setRefinement(*timerefinement,
		     spatialrefinement,gridplacementrefinement);
  return 1;
}

int f_amr_setlref(Long8 *afile,int *level,int *veclen,
		  int *timerefinement,
		  int *spatialrefinement,
		  int *gridplacementrefinement){
  AMRwriter *w = (AMRwriter*)afile;
  if(*veclen==1)
    w->setLevelRefinement(*level,*timerefinement,
			  *spatialrefinement,*gridplacementrefinement);
  else
    w->setLevelRefinement(*level,*timerefinement,
			  spatialrefinement,gridplacementrefinement);
  return 1;
}

int f_amr_setlevel (Long8 *afile,int *level){
  AMRwriter *w = (AMRwriter*)afile;
  w->setLevel(*level);
  return 1;
}

int f_amr_settime (Long8 *afile,int *timestep){
  AMRwriter *w = (AMRwriter*)afile;
  w->setTime(*timestep);
  return 1;
}

int f_amr_incrtime(Long8 *afile){
  AMRwriter *w = (AMRwriter*)afile;
  w->incrementTime();
  return 1;
}

int f_amr_write(Long8 *afile,int *origin, int *dims, void *data){
  AMRwriter *w = (AMRwriter*)afile;
  w->write(origin,dims,data);
  return 1;
}
int f_amr_writef(Long8 *afile,float *origin, int *dims, void *data){
  AMRwriter *w = (AMRwriter*)afile;
  w->write(origin,dims,data);
  return 1;
}
int f_amr_writed(Long8 *afile,double *origin, int *dims, void *data){
  AMRwriter *w = (AMRwriter*)afile;
  w->write(origin,dims,data);
  return 1;
}
