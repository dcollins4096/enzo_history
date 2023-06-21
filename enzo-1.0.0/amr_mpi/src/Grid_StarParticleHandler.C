/***********************************************************************
/
/  GRID CLASS (HANDLE THE CREATION AND FEEDBACK OF STAR PARTICLES)
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "StarParticleData.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

#define NO_STAR1

#ifdef STAR1
extern "C" void FORTRAN_NAME(star_maker1)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *dx, FLOAT *t, float *z, int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     		 int *ibuff, hydro_method *imethod,
             float *odthresh, float *massff, float *smthrest, int *level,
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp);
#endif /* STAR1 */

extern "C" void FORTRAN_NAME(star_maker2)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *dx, FLOAT *t, float *z, 
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     		 int *ibuff, 
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf);

extern "C" void FORTRAN_NAME(star_maker3)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *zfield1, float *zfield2,
             float *dx, FLOAT *t, float *z, 
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     		 int *ibuff, 
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf);

extern "C" void FORTRAN_NAME(star_maker4)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *zfield1, float *zfield2,
             float *dx, FLOAT *t, float *z, 
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     		 int *ibuff, 
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf);


#ifdef STAR1
extern "C" void FORTRAN_NAME(star_feedback1)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, 
		       float *w, float *dt, float *r, float *dx, 
                       FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     				 int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *te, float *ge, 
		       int *idual);
#endif /* STAR1 */

extern "C" void FORTRAN_NAME(star_feedback2)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v, 
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *r, float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, 
			float *justburn);

extern "C" void FORTRAN_NAME(star_feedback3)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v, 
		       float *w, float *metal, float *zfield1, float *zfield2,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *r, float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, 
			float *justburn);

extern "C" void FORTRAN_NAME(star_feedback4)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v, 
		       float *w, float *metal, float *zfield1, float *zfield2,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *r, float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, 
			float *justburn);

extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest, 
                                   int *sdim1, int *sdim2, int *sdim3, 
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3, 
                                   int *dstart1, int *dstart2, int *dststart3);



int grid::StarParticleHandler(int level)
{
  return FAIL;
}
