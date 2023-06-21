// AMRwriter
#ifndef __AMRWRITER_HH_
#define __AMRWRITER_HH_

#include "Arch.h"
#include "IO.hh"
#include "Writer.hh"
#include "FlexArrayTmpl.H"

class AMRwriter : public Writer{
protected:
  struct LevelParams {
    int srefine[3]; // grid refinement
    int prefine[3]; // placement refinement
    int trefine; // time refinement
  };
  double basetimestep,currentorigin[3],currentdelta[3];
  FlexArray<LevelParams> levels;
  int iorigin[3];
  int currentlevel,currentstep;
	
  void writeAMRattributes();
private:
  virtual void setOrigin(int *origin); // logical origin
  virtual void setOrigin(float *origin); // real origin
  virtual void setOrigin(double *origin); // real origin
  // virtual void setDelta(double *delta) { Writer::setOrigin(delta); }
  virtual void setDims(int *dims) { Writer::setDims(dims);}
  virtual void setDims(int rank, int *dims) { Writer::setDims(rank,dims); }
  virtual void write(void *data);
  virtual void setRank(int rank) {Writer::setRank(rank); } 
public: //====================================================
  enum Flags {MaxDepth=-1};
  AMRwriter(IObase &descriptor);
  virtual ~AMRwriter();
  //------------Initialization Methods------------------
  virtual void setType(IObase::DataType numbertype) { Writer::setType(numbertype); }
  virtual void setTopLevelParameters(int rank,double *origin,
				     double *delta,double timestep,int maxdepth);
  virtual void setRefinement(int timerefinement,
			     int *spatialrefinement,
			     int *gridplacementrefinement=0);  
  virtual void setRefinement(int timerefinement,
			     int spatialrefinement,
			     int gridplacementrefinement=1);
  virtual void setLevelRefinement(int level,
				  int timerefinement,
				  int *spatialrefinement,
				  int *gridplacementrefinement=0);
  virtual void setLevelRefinement(int level,
				  int timerefinement,
				  int spatialrefinement,
				  int gridplacementrefinement=1);
  //-----------Stepping Parameters------------------
  virtual void setLevel(int level) { currentlevel=level; }
  virtual void setTime(int timestep) { currentstep=timestep;}
  virtual void incrementTime() {currentstep++;}
  // virtual void setDeltaTime(double dt) {deltatime=dt;}
  virtual void write(int *origin,int *dims,void *data){
    setOrigin(origin);
    setDims(dims);
    write(data);
  } 
  virtual void write(float *origin,int *dims,void *data){
    setOrigin(origin);
    setDims(dims);
    write(data);
  }
  virtual void write(double *origin,int *dims,void *data){
    // create iorigin from origin...
    setOrigin(origin);
    setDims(dims);
    write(data);
  }
};


class FrameworkAMRwriter : protected AMRwriter {
  int nlevels;
  int maxrefinement;
  int refinement;
public:
   /*@@
   @routine    FrameworkAMRwriter::FrameworkAMRwriter
   @date       Tue Apr 15 14:23:32 1997
   @author     John Shalf
   @desc 
   Constructor for the FrameworkAMRwriter.  It makes many assumptions about the
   way that users of the framework would like to store their data.  This eliminates
   or hides most of the calls that AMRwriter:: would need to describe the AMR grid
   heirarchy.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

  FrameworkAMRwriter(IObase &descriptor):AMRwriter(descriptor),
    nlevels(1),maxrefinement(1),refinement(2){}
  virtual ~FrameworkAMRwriter() {}
  /*
    This is set ONCE right after you open your AMR file.  It is not
    built into the constructor because you may not have the information
    necessary to do this at the time of construction...
  */
   /*@@
   @routine    FrameworkAMRwriter::setParameters
   @date       Tue Apr 15 14:00:40 1997
   @author     John Shalf
   @desc   This is set ONCE right after you open your AMR file.  It is not
    built into the constructor because you may not have the information
    necessary to do this at the time of construction...
    This sets parameters necessary do define the AMR heirarchy, namely
    the refinement between levels and the floating-point parameters that
    form the basis for locating the grids in 3-space in a visualization system.
   @enddesc 
      @par     rank
   @pdesc   The number of dimensions for the dataset
   @ptype   int
   @pvalues 1-5
   @pcomment 
   This calls setRank() in the base Writer() class
   @endpar 

      @par     origin
   @pdesc   The origin of the toplevel (coarsest) grid in REAL coordinates
   @ptype   double*
   @endpar 
   
      @par     delta
   @pdesc   The inter-grid spacing at the toplevel (coarsest level) in REAL coordinates.
   @ptype   double
   @pcomment 
   The inter-grid spacing is divided by the refinement factor for a particular level
   to determine the real grid spacing for that nested grid.
   This assumes you have a uniform grid delta (same delta in all spatial directions).
   @endpar 

      @par     timestep
   @pdesc   The REAL size of a toplevel (coarsest) timestep
   @ptype   double
   @pvalues 
   @endpar 

      @par     interlevelRefinementRatio
   @pdesc   The integer ratio of time refinement and spatial refinement between levels.
   @ptype   int
   @pvalues Any positive integer > 0
   @pcomment 
   This ratio is used to determine the spatial refinement for a level.
   So spatialrefinement_on_this_level = toplevel_delta / interlevelRefinementRatio^level
   Where the levels are numbered from 0:nlevels-1.
   This assumes that you have the same refinement factor in all spatial directions as
   well as in time.
   @endpar 

      @par     numLevels
   @pdesc   The maximum number of levels in the AMR heirarchy.
   @ptype   int
   @pvalues Any positive integer > 0
   @pcomment 
   This is necessary to find out the resolution of the finest-grid for the 
   purposes of grid-placement in the AMR file.  Otherwise we wouldn't be able
   to interpret the above parameters.
   @endpar 

   @calls     AMRwriter::setTopLevelParameters
   @calledby   
   @history 
   @endhistory 

@@*/

  virtual void init(int rank, // number of dimensions in the dataset
		    double *origin, // REAL origin or coarsest level grid
		    double *delta, // float grid spacing at coarsest level
		    double timestep, // float timestep at coarsest level
		    int interlevelRefinementRatio, //refinement ratio with
		    // between levels.  This covers refinement of time stepping
		    // as well
		    int numLevels); // the maximum depth of the AMR grid heirarchy
  virtual void init(  IObase::DataType dt, // type for all data
		      int rank, // number of dimensions in the dataset
		      double *origin, // REAL origin or coarsest level grid
		      double *delta, // float grid spacing at coarsest level
		      double timestep, // float timestep at coarsest level
		      int interlevelRefinementRatio, //refinement ratio with
		      // between levels.  This covers refinement of time stepping
		      // as well
		      int numLevels);
  // This is called on every write
  /*@@
    @routine    FrameworkAMRwriter::write
    @date       Tue Apr 15 14:14:48 1997
   @author     John Shalf
   @desc 
   This is called on every write to set the parameters for a particular grid.
   @enddesc 
   @calls  AMRwriter::setGrid AMRwriter::write 
   @calledby
      @par     level
   @pdesc   The current AMR level for this grid
   @ptype   int
   @pvalues 0 to nlevels-1
   @endpar 

      @par     globaltimestep
   @pdesc   The current global timestep (stepping at the finest resolution level)
   @ptype   int
   @pvalues Any positive integer
   @pcomment 
   So we step at finest time resolution instead of resolution relative to this level.
   @endpar 

      @par     origin
   @pdesc   The integer origin of the grid using coordinates relative to the finest resolution grid in the heirarchy.
   @ptype   int*
   @pvalues Any positive integer > 0
   @pcomment 
   So grid placement is with respect to the finest level integer coordinates.
   @endpar 

      @par     dims
   @pdesc   The dimensions of the array (dims of the grid).
   @ptype   int*
   @endpar 

      @par     data
   @pdesc   Pointer to the data array for the grid.
   @ptype   void*
   @endpar 


   @history 
 
   @endhistory 

@@*/

  virtual void write(int level,
	     int globaltimestep,
	     int *origin,
	     int *dims,
	     void *data) {
    AMRwriter::setLevel(level);
    AMRwriter::setTime(globaltimestep);
    AMRwriter::write(origin,dims,data);
  }
  
};

#define f_amr_begin F77NAME(amr_begin_,amr_begin,AMR_BEGIN)
#define f_amr_end F77NAME(amr_end_,amr_end,AMR_END)
#define f_amr_settype F77NAME(amr_settype_,amr_settype,AMR_SETTYPE)
#define f_amr_setparams F77NAME(amr_setparams_,amr_setparams,AMR_SETPARAMS)
#define f_amr_setref F77NAME(amr_setref_,amr_setref,AMR_SETREF)
#define f_amr_setlref F77NAME(amr_setlref_,amr_setlref,AMR_SETLREF)
#define f_amr_setdims F77NAME(amr_setdims_,amr_setdims,AMR_SETDIMS)
#define f_amr_setlevel F77NAME(amr_setlevel_,amr_setlevel,AMR_SETLEVEL)
#define f_amr_settime F77NAME(amr_settime_,amr_settime,AMR_SETTIME)
#define f_amr_incrtime F77NAME(amr_incrtime_,amr_incrtime,AMR_INCRTIME)
#define f_amr_write F77NAME(amr_write_,amr_write,AMR_WRITE)

extern "C" {
#include "AMRwriter.h"
Long8 f_amr_begin (Long8 *descriptor);
int f_amr_end (Long8 *afile);
int f_amr_settype (Long8 *afile,int *numbertype);
int f_amr_setparams (Long8 *afile,int *rank,double *origin,
				   double *delta, double *timestep,int maxdepth);
int f_amr_setref (Long8 *afile,int veclen,
		      int *timerefinement,
		      int *spatialrefinement,int *gridplacementrefinement);
int f_amr_setlref (Long8 *afile,int *level,int *veclen,
			 int *timerefinement,
			 int *spatialrefinement,int *gridplacementrefinement);
int f_amr_setlevel (Long8 *afile,int *level);
int f_amr_settime (Long8 *afile,int *timestep);
int f_amr_incrtime (Long8 *afile);
int f_amr_write (Long8 *afile,int *origin,int *dims,void *data);
  int f_amr_writef (Long8 *afile,float *origin,int *dims,void *data);
  int f_amr_writed (Long8 *afile,double *origin,int *dims,void *data);
}

#endif
