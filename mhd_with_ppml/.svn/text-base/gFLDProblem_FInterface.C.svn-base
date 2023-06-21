/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Implicit Problem Class Fortran 
/  interfaces.
/
/  written by: Daniel Reynolds
/  date:       November, 2006
/  modified1:  
/
/  PURPOSE: provides C++ interfaces to the relevant Fortran 
/           computational kernels
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"

/* Fortran function prototypes */
extern "C" void FORTRAN_NAME(gfldproblem_matrixentries)(
   double *matentries, float *Eg, float *Eg0, float *Temp, 
   float *sigmaA, float *sigmaS, float *adjvec, int *LimImp, 
   float *dt, FLOAT *a, float *theta, float *dx, float *dy, 
   float *dz, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, 
   int *BCxL, int *BCxR, int *BCyL, int *BCyR, int *BCzL, 
   int *BCzR, int *xlface, int *xrface, int *ylface, 
   int *yrface, int *zlface, int *zrface, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_diffrhs)(
   float *drhs, float *Eg, float *Eg0, float *Temp, float *sigmaA, 
   float *sigmaS, int *LimImp, FLOAT *a, float *dx, float *dy, 
   float *dz, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);
  
extern "C" void FORTRAN_NAME(gfldproblem_locecrhs)(
   float *ecrhs, float *vx, float *vy, float *vz, float *rho, 
   float *ec, float *Eg, float *nHI, float *nHeI, float *nHeII, 
   float *ne, float *Temp, FLOAT *a, FLOAT *adot, float *gamma, 
   int *model, int *NTbins, float *Tlo, float *Thi, float *ceHI, 
   float *ceHeI, float *ceHeII, float *ciHI, float *ciHeI, 
   float *ciHeIS, float *ciHeII, float *reHII, float *reHeII1, 
   float *reHeII2, float *reHeIII, float *brem, float *CompA, 
   float *Comp_xray, float *Comp_temp, float *IsE, float *IsEsHI,
   float *IsEsHInu, float *IsEsHeI, float *IsEsHeInu, 
   float *IsEsHeII, float *IsEsHeIInu, float *aunits, 
   float *rhounits, float *timeunits, float *lenunits,  
   int *Nchem, float *dx, float *dy, float *dz, int *Nx, 
   int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, 
   int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_locegrhs)(
   float *Egrhs, float *Eg, float *Temperature, float *n_HI, 
   float *rho, float *ne, float *Kappa, FLOAT *a, FLOAT *adot, 
   int *Nchem, int *Model, float *IsE, int *Nx, int *Ny, 
   int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, 
   int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_locnirhs)(
   float *rhs_HI, float *rhs_HeI, float *rhs_HeII, int *Nchem, 
   float *n_HI, float *n_HeI, float *n_HeII, float *Eg, 
   float *Temperature, float *rho, float *ne, int *Model, 
   FLOAT *a, FLOAT *adot, float *IsE, float *IsEsHI, float *IsEsHInu, 
   float *IsEsHeI, float *IsEsHeInu, float *IsEsHeII, 
   float *IsEsHeIInu, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_blocksolve)(
   float *Amat, float *xvec, float *bvec, int *N, int *M, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_opacity)(
   float *Kappa, float *n_HI, float *n_HeI, float *n_HeII, FLOAT *a, 
   int *Model, float *IsE, float *IsEsHI, float *IsEsHInu, 
   float *IsEsHeI, float *IsEsHeInu, float *IsEsHeII, 
   float *IsEsHeIInu, int *Nchem, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);



/* C++ interface wrappers */
/* These are each designed to extract as much relevant information 
   as possible from the gFLDProblem class, so that the argument 
   lists are simplified to only the relevant input/output arguments. */

/********/
int gFLDProblem::MatrixEntries(double *matentries, float *Eg, float *Eg0, 
			       float *Temperature, float *sigA, 
			       float *sigS, float *adjvec) 
{
  int xlface = (OnBdry[0][0]) ? 1 : 0;
  int xrface = (OnBdry[0][1]) ? 1 : 0;
  int ylface = (OnBdry[1][0]) ? 1 : 0;
  int yrface = (OnBdry[1][1]) ? 1 : 0;
  int zlface = (OnBdry[2][0]) ? 1 : 0;
  int zrface = (OnBdry[2][1]) ? 1 : 0;
  int ier;
  FORTRAN_NAME(gfldproblem_matrixentries)
    (matentries, Eg, Eg0, Temperature, sigA, sigS, adjvec, 
     &LimImp, &dt, &anew, &theta, &dx[0], &dx[1], &dx[2], &LocDims[0], 
     &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], 
     &BdryType[0][0], &BdryType[0][1], &BdryType[1][0], 
     &BdryType[1][1], &BdryType[2][0], &BdryType[2][1], 
     &xlface, &xrface, &ylface, &yrface, &zlface, &zrface, &ier);
  return(ier);
}

/********/
int gFLDProblem::DiffRHS(float *drhs, float *Eg, float *Eg0, 
			 float *Temperature, float *sigA, 
			 float *sigS, FLOAT *a)
{
  int ier;
  FORTRAN_NAME(gfldproblem_diffrhs)
    (drhs, Eg, Eg0, Temperature, sigA, sigS, &LimImp, a, 
     &dx[0], &dx[1], &dx[2], &LocDims[0], &LocDims[1], &LocDims[2], 
     &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
     &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}

  
/********/
int gFLDProblem::LocEcRHS(float *ecrhs, float *ec, float *Eg, 
			  float *Temperature, float *nHI, float *nHeI, 
			  float *nHeII, FLOAT *a, FLOAT *adot, 
			  float *aUnits, float *DensityUnits,
			  float *TimeUnits, float *LenUnits)
{
  int ier;
  FORTRAN_NAME(gfldproblem_locecrhs)
    (ecrhs, vx, vy, vz, rho, ec, Eg, nHI, nHeI, nHeII, ne, 
     Temperature, a, adot, &Gamma, &Model, 
     &CoolData.NumberOfTemperatureBins, &CoolData.TemperatureStart, 
     &CoolData.TemperatureEnd, CoolData.ceHI, CoolData.ceHeI, 
     CoolData.ceHeII, CoolData.ciHI, CoolData.ciHeI, CoolData.ciHeIS, 
     CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, 
     CoolData.reHeII2, CoolData.reHeIII, CoolData.brem,
     &CoolData.comp, &CoolData.comp_xray, &CoolData.temp_xray, 
     &intSigE, &intSigESigHI, &intSigESigHInu, &intSigESigHeI, 
     &intSigESigHeInu, &intSigESigHeII, &intSigESigHeIInu, 
     aUnits, DensityUnits, TimeUnits, LenUnits, &Nchem, &dx[0], &dx[1], 
     &dx[2], &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], 
     &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], 
     &GhDims[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::LocEgRHS(float *Egrhs, float *Eg, float *Temperature, 
			  float *nHI, float *Kappa, FLOAT *a, FLOAT *adot)
{
  int ier;
  FORTRAN_NAME(gfldproblem_locegrhs)
    (Egrhs, Eg, Temperature, nHI, rho, ne, Kappa, a, 
     adot, &Nchem, &Model, &intSigE, &LocDims[0], &LocDims[1], 
     &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::LocNiRHS(float *rhs_HI, float *rhs_HeI, float *rhs_HeII, 
			  float *n_HI, float *n_HeI, float *n_HeII, 
			  float *Eg, float *Temperature, FLOAT *a, FLOAT *adot)
{
  int ier;
  FORTRAN_NAME(gfldproblem_locnirhs)
    (rhs_HI, rhs_HeI, rhs_HeII, &Nchem, n_HI, n_HeI, n_HeII, Eg, 
     Temperature, rho, ne, &Model, a, adot, &intSigE, &intSigESigHI, 
     &intSigESigHInu, &intSigESigHeI, &intSigESigHeInu, 
     &intSigESigHeII, &intSigESigHeIInu, &LocDims[0], &LocDims[1], 
     &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::BlockSolve(float *Amat, float *xvec, float *bvec, 
			    int *N, int *M)
{
  int ier;
  FORTRAN_NAME(gfldproblem_blocksolve)(Amat, xvec, bvec, N, M, &ier);
  return(ier);
}


/********/
int gFLDProblem::Opacity(float *Kappa, float *n_HI, float *n_HeI, 
			 float *n_HeII, FLOAT *a)
{
  int ier;
  FORTRAN_NAME(gfldproblem_opacity)
    (Kappa, n_HI, n_HeI, n_HeII, a, &Model, &intSigE, &intSigESigHI, 
     &intSigESigHInu, &intSigESigHeI, &intSigESigHeInu, 
     &intSigESigHeII, &intSigESigHeIInu, &Nchem, &LocDims[0], 
     &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}

#endif
