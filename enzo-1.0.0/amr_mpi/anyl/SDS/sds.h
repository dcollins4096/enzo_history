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
#ifndef __SDS_H_
#define __SDS_H_

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
#include <netcdf.h>

/*@@
  @header sds.h
  @date Wed Feb 28 20:30:57 1996 @author John Shalf
  @version 0.9a
@@*/

/**************** FileID juggling Operations ********************/
/* make file handle inactive (use if shortage of descriptors) */
void SDSdeactivate CPROTO((char *filename));
void SDSreactivate CPROTO((char *filename));


/******************** Explicit file open & close ops ******************/
/* use to specify file opening and closing, but not 
   really needed since all routines will do this automatically */
void SDSopen CPROTO((char *filename,char *access));
void SDSclose CPROTO((char *filename));
void SDScloseAll CPROTO((void));
void SDSflush CPROTO((char *filename));
void SDSpurge CPROTO((char *filename));


/********************  Random Access Operations ************************/
/* need to make sure that seek leaves everything
   in the proper state so that a second getDims
   after the seek will not advance one more sds */
void SDSseek CPROTO((char *filename,int setnum));
int32 SDSseekName CPROTO((char *filename,char *name));
#define SDSgetLocation(x) SDSgetIndex(x)
int32 SDSgetIndex CPROTO((char *filename));
int32 SDSisCoord CPROTO((char *filename));


/********************** File and Dataset Info **********************/
/* analogue to DFSDgetDims() */
/* should have set coordinates routine which can also handle coords
   that are n-dimensional instead of just edge coords.  So it can
   handle Uniform/Rectilinear/Irregular dataset types */
void SDSgetFileInfo CPROTO((char *filename,int32 *nsds, int32 *nattrib));
int32 SDSgetDataName CPROTO((char *filename,char *name));
int32 SDSgetNumDatasets CPROTO((char *filename));
long SDSgetDims CPROTO((char *filename,char *dataname,
		       int32 *rank,int32 *dims));
int32 SDSgetNT CPROTO((char *filename));
void SDSsetNT CPROTO((char *filename,int32 numbertype));
long SDScomputeSize CPROTO((int32 rank,int32 *dims,int32 numbertype));

/* seems rank should be a pointer */
/* SD and DFSD are rather different here.  In SD, you are selecting
   how much of each dimension to read from the file.  In DFSD, you
   are requesting that the dims of the dataset you read be returned
   to you. */
int SDSreadData CPROTO((char *filename,
		     int32 rank,int32 *dims,VOIDP data));
int SDSreadChunk CPROTO((char *filename,int32 rank,int32 *dims,
		 int32 *origin,int32 *stride,VOIDP data));
/* Can simply write data or write the name with the data.
   The SD interface kind-of requires that the name for the dataset
   be written here, so it would have been inconsistent with the DFSD
   interface that this library is attempting to replace. Hence the #defines */
int32 SDSallocateDataset CPROTO((char *filename,char *dataname,
			  int32 rank,int32 *dims));
int SDSwriteData CPROTO((char *filename,char *dataname,
		      int32 rank,int32 *dims,VOIDP data));
int SDSwriteChunk CPROTO((char *filename,
		       int32 rank,int32 *dims,
		       int32 *origin,int32 *stride,VOIDP data));

/************ Managing Annotations *******************/
int32 SDSaddAnnotation CPROTO((char *filename,char *annotation));
int32 SDSgetAnnotationSize CPROTO((char *filename));
int32 SDSgetAnnotation CPROTO((char *filename,char *annotation,int32 maxlen));


/************ Managing Attributes ******************/
int32 SDSgetNumAttribs CPROTO((char *filename));
int32 SDSwriteAttrib CPROTO((char *filename,char *attribname,
			   int32 numbertype,int32 nelements,VOIDP buffer));
/* 0=Noattrib -1=error +<n>=Number of elements read */
/* nelements= max number of elements the buffer can hold */
int32 SDSreadAttrib CPROTO((char *filename,char *attribname,
			    int32 nelements,VOIDP buffer));
int32 SDSfindAttribInfo CPROTO((char *filename,char *attribname,
			     int32 *numbertype, int32 *nelements));
int32 SDSgetAttribInfo CPROTO((char *filename,int32 index,char *attribname,
			     int32 *numbertype, int32 *nelements));
/* ** Standard String Attributes ** */
int32 SDSgetDataStrs CPROTO((char *filename,char *label,char *units,
			    char *format, char *coordsys, int len));
#endif /* __SDS_H_ */
