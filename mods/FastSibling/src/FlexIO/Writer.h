#ifndef __WRITER_H_
#define __WRITER_H_

#include "Arch.h"

typedef IOFile WRFile;/* its the same, but it is a different object underneath */
WRFile WRbeginFile PROTO((IOFile *descriptor));
void WRendFile PROTO((WRFile afile));

void WRsetRank PROTO((WRFile afile,int rank));
void WRsetType PROTO((WRFile afile,int numbertype));
void WRsetParams PROTO((WRFile afile,
			int rank,int *dims,int type,
			double *origin,double *delta));
void WRsetDims PROTO((WRFile afile,int *dims));
void WRsetRankDims PROTO((WRFile afile,int rank, int *dims));
void WRsetOrigin PROTO((WRFile afile,double *origin));
void WRsetDelta PROTO((WRFile afile,double *delta));
void WRwrite PROTO((WRFile afile,void *data));
void WRreserveChunk PROTO((WRFile afile));
void WRwriteChunk PROTO((WRFile afile,int *dims,int *origin,void *data));
#endif

