#ifndef __AMRWRITER_H_
#define __AMRWRITER_H_

#include "Arch.h"

typedef IOFile AMRFile; /* its the same, but it is a different object underneath */
AMRFile AMRbeginFile PROTO((IOFile descriptor));
void AMRendFile PROTO((AMRFile afile));
void AMRsetType PROTO((AMRFile afile,int numbertype));
void AMRsetTopLevelParameters PROTO((AMRFile afile,int rank,double *origin,
				   double *delta, double timestep,int maxdepth));
void AMRsetRefinement PROTO((AMRFile afile,int timerefinement,
				int *spatialrefinement,int *gridplacementrefinement));
void AMRsetScalarRefinement PROTO((AMRFile afile,int timerefinement,
				int spatialrefinement,int gridplacementrefinement));
void AMRsetLevelRefinement PROTO((AMRFile afile,int level,int timerefinement,
				int *spatialrefinement,int *gridplacementrefinement));
void AMRsetScalarLevelRefinement PROTO((AMRFile afile,int level,int timerefinement,
				int spatialrefinement,int gridplacementrefinement));
/* Stepping Methods */
void AMRsetLevel PROTO((AMRFile afile,int level));
void AMRsetTime PROTO((AMRFile afile,int timestep));
void AMRincrementTime PROTO((AMRFile afile));
void AMRwrite PROTO((AMRFile afile,int *origin, int *dims, void *data));
void AMRwriteFloat PROTO((AMRFile afile,float *origin, int *dims, void *data));
void AMRwriteDouble PROTO((AMRFile afile,double *origin, int *dims, void *data));

/*----------For the Framework AMR---------------*/
typedef IOFile fAMRFile;
fAMRFile fAMRbeginFile PROTO((IOFile descriptor));
void fAMRendFile PROTO((fAMRFile afile));
void fAMRsetParameters PROTO((fAMRFile afile,
				   int datatype,
				   int rank,
				   double *origin,
				   double *delta,
				   double timestep,
				   int interlevelRefinementRatio,
				   int nlevels));
void fAMRwrite PROTO((fAMRFile afile,
		      int level,
		      int globaltimestep,
		      int *origin,
		      int *dims,
		      void *data));
#endif

