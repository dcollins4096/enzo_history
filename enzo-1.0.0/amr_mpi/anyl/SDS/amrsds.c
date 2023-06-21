/*****************************************************************************
 *                                                                           *
 * Copyright 1993-2004, Laboratory for Computational Astrophysics            *
 * Copyright 1993-2000, Board of Trustees of the University of Illinois      *
 * Copyright 2000-2004, Regents of the University of California              *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf.h>

#include "amrsds.h"
#include "sds.h"

/*@@
   @file      amrsds.c
   @date      Wed Feb 26 20:05:00 1996
   @author   John Shalf 
   @desc 
   	The AMR routines are used to store Heirarchial Adaptive Mesh
   	Refinement data in an HDF-compatible format.  It is implemented
   	using the <a href="http://bach.ncsa.uiuc.edu/IO/SDSlibrary.html">
   	SDS wrapper routines</a>. (seealso DAGH)<p>

	Heirarchial Adaptive Mesh Refinement can be very difficult to 
	represent in a flattened format.  There are many different ways
	to represent the grid information.  The choices made here are
	intended to be applicable to AMR systems which use rectilinear
	grids with an unbounded number levels that can vary during the 
	computation.  It also allows for the resolution of the grids
	to be independent in each dimension (including time) and
	to change upon each regridding.  The grids are uniquely 
	identified by their level, gridID, and timestep.  However, the
	processor ID is not stored on multiprocessor systems. Processor
	information can be stored in a separate NetCDF attribute using
	SDSwriteAttrib().
	<hr>
	
<h4>Taxonomy</h4>
The AMR heirarchy resembles a tree.  The root of the tree is level 0
and the lower levels are numbered consecutively.  The nodes of this
tree are "grids".  A "grid" is a rectilinear array of data at a
particular AMR level.  At a particular level there can be multiple
"steps" of the level as it evolves.  The grid can be of a different
size and location on every timestep.<p>

<h4>Explanation of AMR-SDS Parameters</h4>
<i>(the same names are used for all function calls)</i>
<dl>
<dt><b>level:</b>  <dd>The level in the amr heirarchy, starting from 0.
<dt><b>gridID:</b> <dd>At a particular level, there can be multiple grids.
        This is an integer (starting from 0) that differentiates
	grids at the same level.  There is no other implied meaning
	for this number.
<dt><b>step:</b>   <dd>This is the current integer step for a particular 
        level.  It need not be associated with time so that the
	AMR file can store more than one component as successive
	steps in the same graph if they share the same tree structure.

<dt><b>rank:</b>   <dd>As with regular HDF, this is the number of dimensions of 
        the data

<dt><b>dims:</b>   <dd>An n-element array of the dimensions of the data where
        n=the rank of the data.

<dt><b>origin:</b> <dd>An n-element array locating the origin of the current
        grid with respect to the origin of the toplevel grid
	(n = rank of the data at top level).  The resolution
	of the coordinates is with respect to the current level
        is with respect to the toplevel.  ie. A 1D grid that is
	offset by 2 top-level grid points with respect to
	the toplevel origin, and has a refinement factor of 5
	has its origin[0]=10.

<dt><b>resolution:</b> <dd>This is the integer refinement factor by which the 
       current level grid has been refined with respect to the 
       toplevel grid (ie. a refinement factor of 3 means  
       3 grid points at this level fit in the space of one grid point 
       at the toplevel.  It is an (n+1) element array where n=rank 
       of the toplevel grid.  So each dimension can be of a different
       refinement factor.  The last element is the refinement
       factor for time with respect to the toplevel.

<dt><b>realTime:</b> <dd>Since this may not be derivable from any combination of
        the time resolution and level/gridID/step, then we'd best
	store this explicitly.  Since it is a floating point number,
	it can't be relied upon for searching or seeking through
	the heirarchy, but it does provide useful information
	once the grid is found.
</dl>	
   @enddesc 
   @includes amrsds.h sds.h
   @seeheader amrsds.h sds.h
   @version 0.9a
 @@*/

/*------------------------------------------------------*/
#define DEFAULT_ANNOTATION_SZ 512
int AmrDebug=0; /* set this flag to TRUE to turn on verbose error messages */
static AMRinterpolator AMRdefaultInterp=PatchReplace;
/* parses an $amrsds annotation 
   The annotation has the form
   $amrsds=level:gridID:step:rank:origin[0],origin[1],..origin[N]:
           resolution[0],resultion[1],..resolution[N],resolution[N+1]:
	   realtime
   resolution[N+1] is the time-resolution for the grid.
   all items are integers except for "realtime" which is a float.
*/
/*@@
  @routine AMRgetLocation
  @date      Wed Feb 26 20:05:00 1996
  @author   John Shalf
  @desc Reads the minimum subset of the AMR attributes necessary to 
  uniquely identify the current grid; the level, gridID, and timestep.
  @enddesc
  @calls SDSwriteAttrib
  @calledby AMRwriteData
  @seeroutine AMRreadData AMRreadAttribs 

  @par     filename
  @pdesc   Name of the HDF file to open or create
  @ptype   char*
  @endpar 

  @par	   level
  @pdesc   The level in the amr heirarchy, starting from 0.
  @ptype   int32
  @endpar
  
  @par	   gridID
  @pdesc   At a particular level, there can be multiple grids.
        This is an integer (starting from 0) that differentiates
	grids at the same level at a particular timestep.  There is no 
	other implied meaning for this number (ie. it has no meaning
	between timesteps).
  @ptype   int32
  @endpar
  
  @par	   step
  @pdesc   This is the current integer step for a particular 
        level.  It need not be associated with time so that the
	AMR file can store more than one component as successive
	steps in the same graph if they share the same tree structure.
  @ptype   int32
  @endpar

  @par 
  @retval
  @ptype int
  @pvalues 0 on success, -1 on failure
  @endpar

   @history 
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf 
   @hdesc Initial Version   
   @endhistory 

@@*/
int AMRgetLocation(filename,level,gridID,step)
char *filename;
int32 *level;
int32 *gridID;
int32 *step;
{
  SDSreadAttrib(filename,"AMRlevel",1,level);
  SDSreadAttrib(filename,"AMRgridID",1,gridID);
  SDSreadAttrib(filename,"AMRtimestep",1,step);
  return 0;
}

/*@@
  @routine AMRwriteAttribs
  @date      Wed Feb 26 20:05:00 1996
  @author   John Shalf
  @desc Used internally by AMRwriteData to write the attibutes
  	of the AMR grid into the HDF file.  This uses SDSwriteAttrib
  	to attach NetCDF attributes to the HDF dataset.
  @enddesc
  @calls SDSwriteAttrib
  @calledby AMRwriteData
  @seeroutine AMRreadData AMRreadAttribs 

  @par     filename
  @pdesc   Name of the HDF file to open or create
  @ptype   char*
  @endpar 

  @par	   level
  @pdesc   The level in the amr heirarchy, starting from 0.
  @ptype   int32
  @endpar
  
  @par	   gridID
  @pdesc   At a particular level, there can be multiple grids.
        This is an integer (starting from 0) that differentiates
	grids at the same level at a particular timestep.  There is no 
	other implied meaning for this number (ie. it has no meaning
	between timesteps).
  @ptype   int32
  @endpar
  
  @par	   step
  @pdesc   This is the current integer step for a particular 
        level.  It need not be associated with time so that the
	AMR file can store more than one component as successive
	steps in the same graph if they share the same tree structure.
  @ptype   int32
  @endpar
   
  @par	   realtime
  @pdesc   Since this may not be derivable from any combination of
        the time resolution and level/gridID/step, then we'd best
	store this explicitly.  Since it is a floating point number,
	it can't be relied upon for searching or seeking through
	the heirarchy, but it does provide useful information
	once the grid is found.  You can use AMRseekTime to select
	grids based on a particular floating point time value.
  @ptype   float64
  @endpar
  
  @par	   rank
  @pdesc   The number of dimensions of the data (this particular grid).
  @ptype   int32
  @endpar
  
  @par	   dims
  @pdesc   The dimensions of the data (for this particular grid).
  @ptype   int32*
  @endpar
  
  @par	   origin
  @pdesc   An array locating the origin of the current
        grid with respect to the origin of the toplevel grid
	(where the number of elements in the array = rank of the data at top level).  
	The resolution of the coordinates is with respect to the current level
    is with respect to the toplevel.  ie. A 1D grid that is
	offset by 2 top-level grid points with respect to
	the toplevel origin, and has a refinement factor of 5
	has its origin[0]=10.
  @ptype   int32*
  @endpar

  @par	   resolution
  @pdesc   This is the integer refinement factor by which the 
       current level grid has been refined with respect to the 
       toplevel grid (ie. a refinement factor of 3 means  
       3 grid points at this level fit in the space of one grid point 
       at the toplevel.  It is an (rank+1) element array.
      The first 'rank' elements are the refinement factors for each
      of the spatial dimensions.  The last element is the refinement
       factor for time with respect to the toplevel.
  @ptype   int32*
  @endpar

  @par 
  @retval
  @ptype int
  @pvalues 0 on success, -1 on failure
  @endpar
  
  @history 
  @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf 
  @hdesc Initial Version   
  @endhistory 
   
   
@@*/
/* Internal Only */
int AMRwriteAttribs(filename,
		     level,gridID,step,realTime,
		     rank,origin,resolution)
char *filename;
int32 level;
int32 gridID;
int32 step;
float64 realTime;
int32 rank;
int32 *origin;
int32 *resolution;
{
  SDSwriteAttrib(filename,"AMRlevel",DFNT_INT32,1,&level);
  SDSwriteAttrib(filename,"AMRgridID",DFNT_INT32,1,&gridID);
  SDSwriteAttrib(filename,"AMRtimestep",DFNT_INT32,1,&step);
  SDSwriteAttrib(filename,"AMRrank",DFNT_INT32,1,&rank);
  SDSwriteAttrib(filename,"AMRorigin",DFNT_INT32,rank,origin);
  SDSwriteAttrib(filename,"AMRresolution",DFNT_INT32,rank+1,resolution);
  SDSwriteAttrib(filename,"AMRrealtime",DFNT_FLOAT64,1,&realTime);
}

/*@@ 
  @routine AMRwriteData
  @date      Wed Feb 26 20:05:00 1996
  @author   John Shalf
  @desc Works just like DFSDadddata(), except it includes
         all of the AMR attributes as well.
  @enddesc
  @calls AMRwriteAttribs SDSwriteData
  @seeroutine AMRreadData AMRreadAttribs

  @par     filename
  @pdesc   Name of the HDF file to open or create
  @ptype   char*
  @endpar 
  
  @par	   dataname
  @pdesc   Name of the dataset.  Can use this to identify the
  	name of the component this grid refers to.
  @ptype   char*
  @endpar
  
  @par	   level
  @pdesc   The level in the amr heirarchy, starting from 0.
  @ptype   int32
  @endpar
  
  @par	   gridID
  @pdesc   At a particular level, there can be multiple grids.
        This is an integer (starting from 0) that differentiates
	grids at the same level at a particular timestep.  There is no 
	other implied meaning for this number (ie. it has no meaning
	between timesteps).
  @ptype   int32
  @endpar
  
  @par	   step
  @pdesc   This is the current integer step for a particular 
        level.  It need not be associated with time so that the
	AMR file can store more than one component as successive
	steps in the same graph if they share the same tree structure.
  @ptype   int32
  @endpar
   
  @par	   realtime
  @pdesc   Since this may not be derivable from any combination of
        the time resolution and level/gridID/step, then we'd best
	store this explicitly.  Since it is a floating point number,
	it can't be relied upon for searching or seeking through
	the heirarchy, but it does provide useful information
	once the grid is found.  You can use AMRseekTime to select
	grids based on a particular floating point time value.
  @ptype   float64
  @endpar
  
  @par	   rank
  @pdesc   The number of dimensions of the data (this particular grid).
  @ptype   int32
  @endpar
  
  @par	   dims
  @pdesc   The dimensions of the data (for this particular grid).
  @ptype   int32*
  @endpar
  
  @par	   origin
  @pdesc   An array locating the origin of the current
        grid with respect to the origin of the toplevel grid
	(where the number of elements in the array = rank of the data at top level).  
	The resolution of the coordinates is with respect to the current level
    is with respect to the toplevel.  ie. A 1D grid that is
	offset by 2 top-level grid points with respect to
	the toplevel origin, and has a refinement factor of 5
	has its origin[0]=10.
  @ptype   int32*
  @endpar

  @par	   resolution
  @pdesc   This is the integer refinement factor by which the 
       current level grid has been refined with respect to the 
       toplevel grid (ie. a refinement factor of 3 means  
       3 grid points at this level fit in the space of one grid point 
       at the toplevel.  It is an (rank+1) element array.
      The first 'rank' elements are the refinement factors for each
      of the spatial dimensions.  The last element is the refinement
       factor for time with respect to the toplevel.
  @ptype   int32*
  @endpar

  @par data
  @pdesc	The data for the grid to be stored in the HDF file.
  @ptype VOIDP
  @endpar

  @par 
  @retval
  @ptype int
  @pvalues 0 on success, -1 on failure
  @endpar
  
   @history 
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf 
   @hdesc Initial Version   
   @endhistory   
 @@*/
int AMRwriteData (filename,dataname,level,gridID,step,realTime,
		rank,dims,origin,resolution,data)
char *filename;
char *dataname;
int32 level;
int32 gridID;
int32 step;
float64 realTime;
int32 rank;
int32 *dims;
int32 *origin;
int32 *resolution;
void *data;
{
  if(SDSwriteData(filename,dataname,rank,dims,data)<0)
  {
    if(AmrDebug) fprintf(stderr,"AMRaddData: DFSDaddData() failed for %s\n",filename);
    return -1;
  }
  AMRwriteAttribs(filename,level,gridID,step,realTime,rank,origin,resolution);
  return 0;
}

/*@@
  @routine AMRreadAttribs
  @date      Wed Feb 26 20:05:00 1996
  @author   John Shalf
  @desc Works just like DFSDadddata(), except it includes
         all of the AMR attributes as well.
  @enddesc
  @calls SDSgetDims SDSreadAttrib
  @seeroutine AMRwriteData

  @par     filename
  @pdesc   Name of the HDF file to open or create
  @ptype   char*
  @endpar 

  @par     dataname
  @pdesc   Array in which to store the name of the dataset
  @ptype   char*
  @pvalues Preallocated character array or NULL
  @endpar  

  @par	   level
  @pdesc   The level in the amr heirarchy, starting from 0.
  @ptype   int32*
  @endpar
  
  @par	   gridID
  @pdesc   At a particular level, there can be multiple grids.
        This is an integer (starting from 0) that differentiates
	grids at the same level at a particular timestep.  There is no 
	other implied meaning for this number (ie. it has no meaning
	between timesteps).
  @ptype   int32*
  @endpar
  
  @par	   step
  @pdesc   This is the current integer step for a particular 
        level.  It need not be associated with time so that the
	AMR file can store more than one component as successive
	steps in the same graph if they share the same tree structure.
  @ptype   int32*
  @endpar
   
  @par	   realtime
  @pdesc   Since this may not be derivable from any combination of
        the time resolution and level/gridID/step, then we'd best
	store this explicitly.  Since it is a floating point number,
	it can't be relied upon for searching or seeking through
	the heirarchy, but it does provide useful information
	once the grid is found.  You can use AMRseekTime to select
	grids based on a particular floating point time value.
  @ptype   float64*
  @endpar
  
  @par	   rank
  @pdesc   The number of dimensions of the data (this particular grid).
  @ptype   int32*
  @endpar
  
  @par	   dims
  @pdesc   The dimensions of the data (for this particular grid).
  @ptype   int32*
  @endpar
  
  @par	   origin
  @pdesc   An array locating the origin of the current
        grid with respect to the origin of the toplevel grid
	(where the number of elements in the array = rank of the data at top level).  
	The resolution of the coordinates is with respect to the current level
    is with respect to the toplevel.  ie. A 1D grid that is
	offset by 2 top-level grid points with respect to
	the toplevel origin, and has a refinement factor of 5
	has its origin[0]=10.
  @ptype   int32*
  @endpar

  @par	   resolution
  @pdesc   This is the integer refinement factor by which the 
       current level grid has been refined with respect to the 
       toplevel grid (ie. a refinement factor of 3 means  
       3 grid points at this level fit in the space of one grid point 
       at the toplevel.  It is an (rank+1) element array.
      The first 'rank' elements are the refinement factors for each
      of the spatial dimensions.  The last element is the refinement
       factor for time with respect to the toplevel.
  @ptype   int32*
  @endpar

  @par 
  @retval
  @ptype int
  @pvalues 0 on success, -1 on failure
  @endpar

   @history 
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf 
   @hdesc Initial Version   
   @endhistory   

 @@*/
int AMRreadAttribs(filename,dataname,
		     level,gridID,step,realTime,
		     rank,dims,origin,resolution)
char *filename;
char *dataname;
int32 *level;
int32 *gridID;
int32 *step;
float64 *realTime;
int32 *dims;
int32 *rank;
int32 *origin;
int32 *resolution;
{
  SDSgetDims(filename,dataname,rank,dims);
  SDSreadAttrib(filename,"AMRlevel",1,level);
  SDSreadAttrib(filename,"AMRgridID",1,gridID);
  SDSreadAttrib(filename,"AMRtimestep",1,step);
  SDSreadAttrib(filename,"AMRrank",1,rank);
  SDSreadAttrib(filename,"AMRorigin",*rank,origin);
  SDSreadAttrib(filename,"AMRresolution",*rank+1,resolution);
  SDSreadAttrib(filename,"AMRrealtime",1,realTime);
}

/*@@
   @routine    AMRreadData
   @date       Wed Feb 28 15:37:43 1996
   @desc
          Reads data from an HDF file.
	  This must be preceded by an AMRreadAttribs() which
	  returns the dimensions of the data so that the user
	  can allocate the space to read in the data as well as
	  all other relevant AMR parameters.
   @enddesc
   @seeroutine SDSgetDims AMRreadAttribs
   @calls      SDSreadData

   @par     filename
   @pdesc   Name of the HDF file
   @ptype   char*
   @endpar 

   @par     rank
   @pdesc   Number of dimension in the dataset
   @ptype   int32
   @endpar 

   @par     dims
   @pdesc   Dimensions of the dataset
   @ptype   int32*
   @endpar 

   @par     data
   @pdesc   Pointer to the actual dataset (preallocated by the user)
   @ptype   VOIDP
   @endpar 

  @par 
  @retval
  @ptype int
  @pvalues 0 on success, -1 on failure
  @endpar

   @history
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf
   @hdesc Initial Version
   @endhistory
 @@*/
int AMRreadData (filename,rank,dims,data)
char *filename;
int32 rank;
int32 *dims;
void *data;
{
  return SDSreadData(filename,rank,dims,data);
}



/*@@ ------------------------------------------------------
  @routine AMRseekTime
  @date      Wed Feb 26 20:05:00 1996
  @author   John Shalf
  @desc 	Seeks to the nearest grid with matching
  	level, gridID, and timestep.  The value
  	AMR_ANY can be used as a wildcard that matches
  	any value.
  @enddesc
  @calls SDSreadAttrib SDSseek
  @seeroutine AMRseek

  @par     filename
  @pdesc   Name of the HDF file to open or create
  @ptype   char*
  @endpar 

  @par	   realtime
  @pdesc   The floating point global "real time" value for the
  	AMR evolution. 
  @ptype   float64
  @endpar
  
  @par	   matchtype
  @pdesc  Selects the criterion used for matching a particular
  	timevalue.  For example, <b>AMRseekTime("file.hdf",6.125,"&lt=")</b>
  	matches the next grid (stored linearly in the file) that has
  	a "realtime" attribute which is less than or equal to 6.125.
  	Likewise, <b>AMRseekTime("file.hdf",6.125,"&gt=")</b>
  	matches the next grid greater than or equal to 6.125.
  	The search will wrap around to the beginning
  	of the HDF file and will continue as far as the grid that
  	the search started on.  However, it will not select the
  	current grid as a successful match under any circumstance.
  @ptype   int32
  @pvalues integer or AMR_ANY for wildcard
  @endpar

  @par 
  @retval
  @ptype int
  @pvalues 0 on success, -1 on failure
  @endpar

   @history 
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf 
   @hdesc Initial Version   
   @endhistory 
 @@*/
int AMRseekTime(char *filename,float64 realtime,char *matchtype)
{
  /* seek based on the RealTime.
   Return -1 if you end up seeking back to the start (even on success) */  
  float64 rt;
  int32 lref,i,rank;
  int32 dims[8]; /* actually superflouous */
  int32 tlevel,tgridID,tstep;
  char annotation[DEFAULT_ANNOTATION_SZ];
  int32 nsds=SDSgetNumDatasets(filename)+1; /* limit the search to one cycle */
  int32 start=SDSgetIndex(filename);

  for(i=start;i<=nsds;i++){
    SDSseek(filename,i);
    if(SDSisCoord(filename)) 
      continue;
    if(SDSreadAttrib(filename,"AMRrealtime",1,&rt)<0)
    {
      if(AmrDebug) fprintf(stderr,
	"AMRseek: () AMRparseLocation failed for annotation %s\n",annotation);
      return -1;
    }
    if(!matchtype){ /* use default */
      if(realtime>=rt) /* default matchtype */
	return 0; /* a match (sort-of) */
    }
    if(matchtype[0]=='='){ /* exact match */
      if(rt==realtime)
	return 0;
    }
    else if(matchtype[0]=='<'){
      if(matchtype[1]=='='){
	if(rt<=realtime)
	  return 0;
      }
      else { /* just < */
	if(rt<realtime)
	  return 0;
      }
    }
    else if(matchtype[0]=='>'){
      if(matchtype[1]=='='){
	if(rt>=realtime)
	  return 0;
      }
      else { /* just > */
	if(rt>realtime)
	  return 0;
      }
    }
    else {
      fprintf(stderr,"AMRseekRealTime: Don't understand matchtype of %s\n",
	      matchtype);
      break;
    }
  }
  for(i=0;i<=(start-1);i++) /* if start is only place with this time, then fail */
  {
    SDSseek(filename,i);
    if(SDSisCoord(filename)) 
      continue;
    if(SDSreadAttrib(filename,"AMRrealtime",1,&rt)<0)
    {
      if(AmrDebug) fprintf(stderr,
	"AMRseek: () AMRparseLocation failed for annotation %s\n",annotation);
      return -1;
    }if(SDSreadAttrib(filename,"AMRrealtime",1,&rt)<0)
    {
      if(AmrDebug) fprintf(stderr,
	"AMRseek: () AMRparseLocation failed for annotation %s\n",annotation);
      return -1;
    }
    if(!matchtype){ /* use default */
      if(realtime>=rt) /* default matchtype */
	return 0; /* a match (sort-of) */
    }
    if(matchtype[0]=='='){ /* exact match */
      if(rt==realtime)
	return 0;
    }
    else if(matchtype[0]=='<'){
      if(matchtype[1]=='='){
	if(rt<=realtime)
	  return 0;
      }
      else { /* just < */
	if(rt<realtime)
	  return 0;
      }
    }
    else if(matchtype[0]=='>'){
      if(matchtype[1]=='='){
	if(rt>=realtime)
	  return 0;
      }
      else { /* just > */
	if(rt>realtime)
	  return 0;
      }
    }
    else {
      fprintf(stderr,"AMRseekRealTime: Don't understand matchtype of %s\n",
	      matchtype);
      break;
    }
  }
  SDSseek(filename,start); /* seek back to the original */
  return -1; /* no match found */
}


/*@@ 
  @routine AMRseek
  @date      Wed Feb 26 20:05:00 1996
  @author   John Shalf
  @desc 	Seeks to the nearest grid with matching
  	level, gridID, and timestep.  The value
  	AMR_ANY can be used as a wildcard that matches
  	any value.  The search will wrap around to the beginning
  	of the HDF file and will continue as far as the grid that
  	the search started on.  However, it will not select the
  	current grid as a successful match under any circumstance.
  @enddesc
  @calls AMRgetLocation SDSseek
  @seeroutine AMRseekTime

  @par     filename
  @pdesc   Name of the HDF file to open or create
  @ptype   char*
  @endpar 

  @par	   level
  @pdesc   The level in the amr heirarchy, starting from 0.
  @ptype   int32
  @pvalues	integer or AMR_ANY for wildcard
  @endpar
  
  @par	   gridID
  @pdesc   At a particular level, there can be multiple grids.
        This is an integer (starting from 0) that differentiates
	grids at the same level at a particular timestep.  There is no 
	other implied meaning for this number (ie. it has no meaning
	between timesteps).
  @ptype   int32
  @pvalues integer or AMR_ANY for wildcard
  @endpar
  
  @par	   step
  @pdesc   This is the current integer step for a particular 
        level.  It need not be associated with time so that the
	AMR file can store more than one component as successive
	steps in the same graph if they share the same tree structure.
  @ptype   int32
  @pvalues integer or AMR_ANY for wildcard
  @endpar

  @par 
  @retval
  @ptype int
  @pvalues 0 on success, -1 on failure
  @endpar

   @history 
   @hdate Wed Feb 28 20:38:57 1996 @hauthor John Shalf 
   @hdesc Initial Version   
   @endhistory 
 @@*/
int AMRseek(filename,level,gridID,step)
char *filename;
int32 level;
int32 gridID;
int32 step;
{
  
  int lref,i,rank;
  int32 dims[8]; /* actually superflouous */
  int tlevel,tgridID,tstep;
  char annotation[DEFAULT_ANNOTATION_SZ];
  int nsds=SDSgetNumDatasets(filename)+1; /* limit the search to one cycle */
  int start=SDSgetIndex(filename);
  
  for(i=start;i<=nsds;i++){
    SDSseek(filename,i);
    if(SDSisCoord(filename)) 
      continue;
    SDSgetDataName(filename,annotation);
    if(AMRgetLocation(annotation,&tlevel,&tgridID,&tstep))
    {
      if(AmrDebug) fprintf(stderr,
	"AMRseek: () AMRparseLocation failed for annotation %s\n",annotation);
      return -1;
    }
    if((level==AMR_ANYVAL || tlevel==level) && 
       (gridID==AMR_ANYVAL || tgridID==gridID) && 
       (step==AMR_ANYVAL || tstep==step))
      return 0; /* a match */
  }
  for(i=0;i<=(start-1);i++) /* if start is only place with this time, then fail */
  {
    SDSseek(filename,i);
    if(SDSisCoord(filename)) 
      continue;
    SDSgetDataName(filename,annotation);
    if(AMRgetLocation(annotation,&tlevel,&tgridID,&tstep))
    {
      if(AmrDebug) fprintf(stderr,
	"AMRseek: () AMRparseLocation failed for annotation %s\n",annotation);
      return -1;
    }
    if((level==AMR_ANYVAL || tlevel==level) && 
       (gridID==AMR_ANYVAL || tgridID==gridID) && 
       (step==AMR_ANYVAL || tstep==step))
      return 0; /* a match */
  }
  SDSseek(filename,start); /* seek back to the original */
  return -1; /* no match found */
}

int AMRgetRegion(filename,level,step,rank,lb,ub,dims,data)
char *filename;
int32 level; /* get at level of detail of this level */
int32 step;  /* get at this step at specified level */
int32 rank;
int32 *dims,*lb,*ub;
void *data;
{
  return 0;
}

int RegionOverlap(rank,
		  sourceDims,sourceOrigin,sourceResolution,
		  destDims,destOrigin,destResolution)
int *rank;
int32 *sourceDims;
int32 *sourceOrigin;
int32 *sourceResolution;
int32 *destDims;
int32 *destOrigin;
int32 *destResolution;
{
  /* must compute normalized dimensions (normalize to a single resolution)
     And then check for overlap on ALL dimensions */
  int i;
  for(i=0;i<*rank;i++)
  {
    float normSrc,normDst,nStSrc,nStDst,nEndSrc,nEndDst;
    float maxres=(float)((sourceResolution[i]>destResolution[i])?
      (sourceResolution[i]):(destResolution[i]));
    normSrc = (float)(sourceResolution[i]) / maxres;
    normDst = (float)(destResolution[i])   / maxres;
    nStSrc =   normSrc * (float)(sourceOrigin[i]);
    nStDst =   normDst * (float)(destOrigin[i]);
    nEndSrc  = nStSrc  + (float)(sourceDims[i]) * normSrc;
    nEndDst  = nStDst  + (float)(destDims[i])   * normDst;
    if(nStSrc>nEndDst || nStDst>nEndSrc) /* if any dim doesn't overlap */
      return 0;  /* then nothing overlaps */
  }
  return 1;
}

void AMRsetInterpolator (interp)
AMRinterpolator interp;
{
  AMRdefaultInterp=interp;
}

void PatchReplace(rank,
		  sourceDims,sourceOrigin,sourceResolution,sourcedata,
		  destdims,destOrigin,destResolution,destdata)
int32 *rank; /*pointer for f77 compat */
int32 *sourceDims;
int32 *sourceOrigin;
int32 *sourceResolution;
void *sourcedata;
int32 *destdims;
int32 *destOrigin;
int32 *destResolution;
void *destdata;
{
  float *sNorm=(float *)malloc(sizeof(float)**rank);
  int i;
  if(!sNorm) {perror("PatchReplace: malloc failed"); return;}
  /* compute normalization factor for given resolution */
  for(i=0;i<*rank;i++)
    sNorm[i]=(float)(destResolution[i])/(float)(sourceResolution[i]);
  for(i=0;i<*rank;i++)
  {}
  free(sNorm);
}
