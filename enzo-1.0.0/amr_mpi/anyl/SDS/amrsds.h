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
#ifndef __AMRSDS_H_
#define __AMRSDS_H_

/*-----------------------------
CPROTO's for people stuck with K&R compilers (ie. Sun's)
------------------------------*/
#ifndef ANSI
#if defined(__STDC__) ||defined(__cplusplus)||defined(_LANGUAGE_C_PLUS_PLUS)
#define ANSI
#endif
#endif

#ifdef ANSI
#define CPROTO(x) x
#else
#define CPROTO(x) ()
#endif

#include <hdf.h>
extern int AmrDebug; /* set this flag to 1 to get verbose error messages */

/*@@
  @header amrsds.h
  @date Wed Feb 28 20:30:57 1996 @author John Shalf
  @version 0.9a
@@*/


/*========================================================================
Explanation of Parameters

Taxonomy
The AMR heirarchy resembles a tree.  The root of the tree is level 0
and the lower levels are numbered consecutively.  The nodes of this
tree are "grids".  A "grid" is a rectilinear array of data at a
particular AMR level.  At a particular level there can be multiple
"steps" of the level as it evolves.  The grid can be of a different
size and location on every timestep.

AMR Parameters  (the same names are used for all function calls)
level:  The level in the amr heirarchy, starting from 0.
gridID: At a particular level, there can be multiple grids.
        This is an integer (starting from 0) that differentiates
	grids at the same level.  There is no other implied meaning
	for this number.
step:   This is the current integer step for a particular 
        level.  It need not be associated with time so that the
	AMR file can store more than one component as successive
	steps in the same graph if they share the same tree structure.

rank:   As with regular HDF, this is the number of dimensions of 
        the data

dims:   An n-element array of the dimensions of the data where
        n=the rank of the data.

origin: An n-element array locating the origin of the current
        grid with respect to the origin of the toplevel grid
	(n = rank of the data at top level).  The resolution
	of the coordinates is with respect to the current level
        is with respect to the toplevel.  ie. A 1D grid that is
	offset by 2 top-level grid points with respect to
	the toplevel origin, and has a refinement factor of 5
	has its origin[0]=10.

resolution: This is the integer refinement factor by which the 
       current level grid has been refined with respect to the 
       toplevel grid (ie. a refinement factor of 3 means  
       3 grid points at this level fit in the space of one grid point 
       at the toplevel.  It is an (n+1) element array where n=rank 
       of the toplevel grid.  So each dimension can be of a different
       refinement factor.  The last element is the refinement
       factor for time with respect to the toplevel.

realTime: Since this may not be derivable from any combination of
        the time resolution and level/gridID/step, then we'd best
	store this explicitly.  Since it is a floating point number,
	it can't be relied upon for searching or seeking through
	the heirarchy, but it does provide useful information
	once the grid is found.
====================================================================*/


/*------------------------------------------------------
Name: AMRwriteData
Purpose: Works just like DFSDadddata(), except it includes
         all of the AMR attributes as well.
returns: 0 on success
         -1 on failure
-------------------------------------------------------*/
int AMRwriteData CPROTO((char *filename,
			 char *dataname,
			 int32 level,
			 int32 gridID,
			 int32 step,
			 float64 realTime,
			 int32 rank,
			 int32 *dims,
			 int32 *origin,
			 int32 *resolution,
			 void *data));

/*------------------------------------------------------
Name: AMRreadAttribs
Purpose: Analogous to DFSDgetdims() but gets amr attribs
         as well.  Can be used to step through an AMR/HDF file
returns: 0 on success
         -1 on failure
-------------------------------------------------------*/
int AMRreadAttribs CPROTO((char *filename,
			   char *dataname,
			   int32 *level,
			   int32 *gridID,
			   int32 *step,
			   float64 *realtime,
			   int32 *rank,
			   int32 *dims,
			   int32 *origin,
			   int32 *resolution));

/*------------------------------------------------------
Name: AMRreadData
Purpose: Identical in function to DFSDgetdata() or
         SDread() or SDSreadData()
returns: 0 on success
         -1 on failure
-------------------------------------------------------*/
int AMRreadData CPROTO((char *filename,
		      int32 rank,
		      int32 *dims,
		      void *data));

/*------------------------------------------------------
Name: AMRseek
Purpose: A generally inefficient seek function for
         AMR-SDS files.  It searches linearly through
	 the file (hence the inefficiency) to find
	 an SDS of the specified level,gridID, and
	 step.  You can use AMR_ANYVAL as a wildcard
	 for any of the level, grid & step.
returns: 0 on success
         -1 on failure
Note: This can be made more efficient in the future
      by building a table of all grids & their associated
      attributes and SDS ref-numbers in memory.  Then 
      the SDSseek() can go directly to the proper SDS instead
      of searching linearly through the file.
-------------------------------------------------------*/
int AMRseek CPROTO((char *filename,
		    int32 level,
		    int32 gridID,
		    int32 step));

int AMRseekTime CPROTO((char *filename,
			float64 realtime,
			char *matchtype));

/* Interpolates a rectilinear dataset of a specified resolution
   out of the heirarchial data */
int AMRgetRegion CPROTO((char *filename,
			 int32 level,/* get at level of detail of this level */
			 int32 step, /* get at this step at specified level */
			 int32 rank,
			 int32 *dims,
			 int32 *lb,
			 int32 *ub,
			 void *data));

/* returns TRUE if the source & destination regions overlap */
int RegionOverlap CPROTO((int *rank,
			 int32 *sourceDims,
			 int32 *sourceOrigin,
			 int32 *sourceResolution,
			 int32 *destDims, 
			 int32 *destOrigin,
			 int32 *destResolution));

/* 
   type definition for 
   function interpolates source 
   grid into destination grid 
*/
typedef void (*AMRinterpolator) CPROTO((int32 *rank,/*pointer for f77 compat */
					int32 *sourceDims,
					int32 *sourceOrigin,
					int32 *sourceResolution,
					void *sourcedata,
					int32 *destdims, 
					int32 *destOrigin,
					int32 *destResolution,
					void *destdata));

/* Select an interpolation function.  Default=PatchReplace */
void AMRsetInterpolator CPROTO((AMRinterpolator interp));

/* 
   A Simple Patch-Replacement heirarchial interpolator
   It linearly interpolates data from the source Grid into
   the resolution of the destination grid and overwrites
   patches in the destination grid that it replaces.  A
   more sophisticated interpolator would make some attempt
   at flux normalization (ie. interpolate the source data
   into the destination grid so as to minimize mismatches)
*/
void PatchReplace CPROTO((int32 *rank,/*pointer for f77 compat */
			  int32 *sourceDims,
			  int32 *sourceOrigin,
			  int32 *sourceResolution,
			  void *sourcedata,
			  int32 *destdims, 
			  int32 *destOrigin,
			  int32 *destResolution,
			  void *destdata));

#define AMR_ANYVAL -1
#endif /* __AMRSDS_H_ */
