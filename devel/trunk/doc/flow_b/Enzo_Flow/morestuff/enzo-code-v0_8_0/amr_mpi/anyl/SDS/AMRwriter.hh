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
    int srefine[5]; // grid refinement
    int prefine[5]; // placement refinement
    int trefine; // time refinement
  };
  double basetimestep;
  FlexArray<LevelParams> levels;
  int iorigin[5];
  int currentlevel,currentstep;
	
  void writeAMRattributes();
  virtual void setOrigin(double *origin) { Writer::setOrigin(origin); }
  // virtual void setDelta(double *delta) { Writer::setOrigin(delta); }
public:
  AMRwriter(IObase &descriptor);
  virtual ~AMRwriter();	
  virtual void setRank(int rank) {Writer::setRank(rank); } 
  virtual void setType(IObase::DataType numbertype) { Writer::setType(numbertype); }
  virtual void setToplevelParameters(double *origin,
				     double *delta, 
				     double timestep);
  virtual void setLevelParameters(int level,
				  int timerefinement,
				  int *spatialrefinement,
				  int *gridplacementrefinement=0);
  // using scalar values)
  virtual void setLevelParameters(int level,
				  int timerefinement,
				  int spatialrefinement,
				  int gridplacementrefinement=0);
  virtual void setDims(int *dims) { Writer::setDims(dims);}
  virtual void setDims(int rank, int *dims) { Writer::setDims(rank,dims); }
  virtual void setOrigin(int *origin);
  virtual void setLevel(int level) { currentlevel=level; }
  virtual void setTime(int timestep) { currentstep=timestep;}
  virtual void setGrid(int level,
		       int timestep,
		       int *origin,
		       int *dims);
  virtual void write(void *data);
};


class FrameworkAMRwriter : protected AMRwriter {
  int nlevels;
  int maxrefinement;
  int refinement;
  virtual void write(void *data) {AMRwriter::write(data);}
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

  virtual void setParameters(int rank, // number of dimensions in the dataset
			     double *origin, // REAL origin or coarsest level grid
			     double *delta, // float grid spacing at coarsest level
			     double timestep, // float timestep at coarsest level
			     int interlevelRefinementRatio, //refinement ratio with
			     // between levels.  This covers refinement of time stepping
			     // as well
			     int numLevels); // the maximum depth of the AMR grid heirarchy
  virtual void setParameters(  IObase::DataType dt, // type for all data
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
    AMRwriter::setGrid(level,globaltimestep,origin,dims);
    AMRwriter::write(data);
  }
  // virtual void write(void *data) { AMRwriter::write(data); }
};

extern "C" {
#include "AMRwriter.h"
}

#endif
