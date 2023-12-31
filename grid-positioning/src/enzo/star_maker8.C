/***********************************************************************
  CREATES AND ADDS TO SINK PARTICLES

  written by: Greg Bryan
  date:       November, 2002
  modified1: Peng Wang, May, 2008
             Added Bondi-Hoyle accretion and stellar wind feedback. 
             Modified particle merge.
  modified2: Elizabeth Harper-Clark, June 2009
             significantly changed isotropic wind as written by Peng to make more isotropic!
             SW = 1 - magnetic field protostellar jets
             SW = 2 - random direction protostellar jets
             SW = 3 - isotropic main sequence stellar wind
  INPUTS:
    d     - density field
    u,v,w - velocity fields
    r     - refinement field (non-zero if zone is further refined)
    dt    - current timestep
    dx    - zone size (code units)
    t     - current time
    z     - current redshift
    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
    nx,ny,nz - dimensions of field arrays
    ibuff    - number of buffer zones at each end of grid
    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
    massthresh - maximum allowed mass in a cell
    level - current level of refinement
    procnum - processor number (for output)
    npold - number of sink particles already on grid
    xp-zpold - current sink particle positions
    up-wpold - current sink particle velocities
    mpold - curernt sink particle masses
    tcp   - star particle creation time (0 for dm particles)
    tdp   - star particle dynamic time (not used here)
    nmax     - particle array size specified by calling routine
    x/y/z start - starting position of grid origin
    jlrefine - Jeans length refinement (<0 for none)
    temp - temperature

  OUTPUTS:
    np   - number of particles created
    xp,yp,zp - positions of created particles
    up,vp,wp - velocities of created particles
    mp       - mass of new particles
***********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#define USE

/* function prototypes */
int  GetUnits(float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MAssUnits, FLOAT Time);

int star_maker8(int *nx, int *ny, int *nz, int *size, float *d, float *te, float *ge, 
		float *u, float *v, float *w, float *bx, float *by, float *bz,
		float *dt, float *r, float *dx, FLOAT *t, 
		float *z, int *procnum, float *d1, float *x1, float *v1, 
		float *t1, int *nmax, FLOAT *xstart, FLOAT *ystart, 
		FLOAT *zstart, int *ibuff, int *imethod, int *idual,
		float *massthresh, int *level, int *np, FLOAT *xp, FLOAT *yp, 
		FLOAT *zp, float *up, float *vp, float *wp, float *mp, 
		float *tcp, float *tdp, float *dm, 
		int *type, int *npold, FLOAT *xpold, 
		FLOAT *ypold, FLOAT *zpold, float *upold, float *vpold, 
		float *wpold, float *mpold, float *tcpold, float *tdpold, float *dmold,
		float *nx_jet, float *ny_jet, float *nz_jet,
		int *typeold, int *idold, int *ctype, float *jlrefine, float *temp,
		float *gamma, float *mu, int *nproc, int *nstar)
{

  int		i, j, k, index, ii, inew, n, bb, cc, nsinks, closest;
  int           xo, yo, zo;
  float		densthresh, maxdens, adddens, ugrid, vgrid, wgrid;
  double	jeansthresh, jlsquared, dx2, dist2, total_density, nearestdx2;
  FLOAT		xpos, ypos, zpos, delx, dely, delz;
  double        Pi = 3.1415926;

  /* Compute Units. */
 
  /* float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
    VelocityUnits = 1;
  double MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits,  TimeUnits) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
    }*/

  /* Convert mass threshold to density */

  densthresh = *massthresh / pow(*dx,3);
  //densthresh = 1e-12/(*d1);
  dx2 = (*dx)*(*dx);

  if (*jlrefine > 0) {
    jlsquared = ((double)((*gamma) * 3.14159 * 1.38e-16 / 6.673e-08) / 
		 ((double)(*d1) * 1.673e-24)) / pow(*x1,2) / (*mu) / pow((*jlrefine),2);
  }

  /* Set new particle index to number of created star particles */

  ii = *np;

  /* Look for existing sink particles */

  nsinks = 0;
  int *sink_index = new int[50000];
  for (n = 0; n < *npold; n++) {
    if (typeold[n] == *ctype) {
      sink_index[nsinks++] = n;
    }
  }
  printf("star_maker8: nsinks = %"ISYM"\n", nsinks);

  /* Merge any sink particles that are close enough to each other */

  double mfrac_b, mfrac_c, total_mass, mi, mj;
  double msun = 1.989e33;
  double umass = (*d1)*pow(*x1,3)/msun;
  float SinkMergeMass = 0.001/umass;
  //printf("star_maker8: SinkMergeDistance = %"FSYM"\n",SinkMergeDistance );
  if (*level == MaximumRefinementLevel && SinkMergeDistance > 0.0) {
    for (i = 0; i < nsinks-1; i++) {
      printf("star_maker8: Merging alogrithm called\n");
      bb = sink_index[i];
      mi = mpold[bb]*pow(*dx,3);
      
      if (mi <= 0.0) continue;
      
      for (j = i+1; j < nsinks; j++) {
	
	cc = sink_index[j];
	mj = mpold[cc]*pow(*dx,3);
	
	if (mj <= 0.0 || cc == bb) continue;
	if (mi > SinkMergeMass && mj > SinkMergeMass) continue;
		
	delx = xpold[bb] - xpold[cc];
	dely = ypold[bb] - ypold[cc];
	delz = zpold[bb] - zpold[cc];
	dist2 = delx*delx + dely*dely + delz*delz;
	
	if (dist2 > pow(SinkMergeDistance,2)) continue;
	
	/* Do the merging */

	if (mj < SinkMergeMass) {

	  /* Merge j to i */

	  total_mass = mpold[bb] + mpold[cc];
	  mfrac_b = mpold[bb] / total_mass;
	  mfrac_c = mpold[cc] / total_mass;
	  xpold[bb] = xpold[bb]*mfrac_b + xpold[cc]*mfrac_c;
	  ypold[bb] = ypold[bb]*mfrac_b + ypold[cc]*mfrac_c;
	  zpold[bb] = zpold[bb]*mfrac_b + zpold[cc]*mfrac_c;
      
	  upold[bb] = upold[bb]*mfrac_b + upold[cc]*mfrac_c;
	  vpold[bb] = vpold[bb]*mfrac_b + vpold[cc]*mfrac_c;
	  wpold[bb] = wpold[bb]*mfrac_b + wpold[cc]*mfrac_c;
	  mpold[bb] = total_mass;
	  dmold[bb] += dmold[cc];
    
	  // Set second particle to be ignored (no mass)
	  tcpold[cc] = 0.0;
	  dmold[cc] = 0.0;
	  upold[cc] = vpold[cc] = wpold[cc] = 0.0;
	  mpold[cc] = FLOAT_UNDEFINED;
	  
	} else {

	  /* Merge i to j */

	  total_mass = mpold[bb] + mpold[cc];
	  mfrac_b = mpold[bb] / total_mass;
	  mfrac_c = mpold[cc] / total_mass;
	  xpold[cc] = xpold[bb]*mfrac_b + xpold[cc]*mfrac_c;
	  ypold[cc] = ypold[bb]*mfrac_b + ypold[cc]*mfrac_c;
	  zpold[cc] = zpold[bb]*mfrac_b + zpold[cc]*mfrac_c;
      
	  upold[cc] = upold[bb]*mfrac_b + upold[cc]*mfrac_c;
	  vpold[cc] = vpold[bb]*mfrac_b + vpold[cc]*mfrac_c;
	  wpold[cc] = wpold[bb]*mfrac_b + wpold[cc]*mfrac_c;
	  mpold[cc] = total_mass;
	  dmold[cc] += dmold[bb];
    
	  // Set second particle to be ignored (no mass)
	  tcpold[bb] = 0.0;
	  dmold[bb] = 0.0;
	  upold[bb] = vpold[bb] = wpold[bb] = 0.0;
	  mpold[bb] = FLOAT_UNDEFINED;

	  /* Now we are done with the ith sink */
	  break;

	}	  
	  
      }  // ENDIF merge particle 

    } // ENDFOR first old particle

    /* Remove deleted particle from sink particle list */

    int nRemoved = 0;
    for (n = 0; n < nsinks; n++) {
      if (mpold[sink_index[n]] < 0.0) {
	for (bb = n+1; bb < nsinks; bb++) {
	  sink_index[bb-1] = sink_index[bb];
	}
	nRemoved++;
      }
    }
    nsinks -= nRemoved;

    if (nRemoved > 0) 
      printf("star_maker8[remove]: Ignoring %"ISYM" sink particles.\n", nRemoved);

  } // if (level == maxlevel)	  





  /* sink particle accretes gas from parent cell according to modified Bondi-Hoyle formula. 
   Reference: M. Ruffert, ApJ (1994) 427 342 */

  double csgrid2;
  float densgrid, tempgrid, msink, drho,
    usink, vsink, wsink, vrel2, mdot, e, de;
  float pi = 4.0*atan(1.0);
  FLOAT r_bh;
  for (n = 0; n < nsinks; n++) {

    bb = sink_index[n];

    if (mpold[bb] < 0.0) continue;

    i = (xpold[bb] - *xstart) / (*dx);
    j = (ypold[bb] - *ystart) / (*dx);
    k = (zpold[bb] - *zstart) / (*dx);
    index = i + (j + k * (*ny)) * (*nx);
    
    densgrid  = d[index];
    ugrid     = u[index];
    vgrid     = v[index];
    wgrid     = w[index];
    tempgrid  = temp[index];    
    msink     = mpold[bb] * pow(*dx,3);
    usink     = upold[bb];
    vsink     = vpold[bb];
    wsink     = wpold[bb];

    csgrid2 = 1.38e-16 * tempgrid / 1.673e-24 / pow(*v1,2);
    vrel2 = pow(ugrid-usink,2) + pow(vgrid-vsink,2) + pow(wgrid-wsink,2);
    r_bh = msink / (csgrid2 + vrel2);
    densgrid *= min(pow((*dx)/r_bh, 1.5), 1.0);
    mdot = 4.0 * pi * densgrid * pow(r_bh, 2) * sqrt(1.2544*csgrid2 + vrel2);
    drho = min(mdot * (*dt) / pow(*dx,3), 0.25 * d[index]);
    
    /*maxdens = jlsquared * temp[index] / dx2;
      drho = max(0.0, d[index] - maxdens);*/        
    //printf("star_maker8: Accretion routine, mass added = %"FSYM"\n",drho*pow(*dx,3)*umass);

    upold[bb] = (mpold[bb]*usink + drho*ugrid) / (mpold[bb] + drho);
    vpold[bb] = (mpold[bb]*vsink + drho*vgrid) / (mpold[bb] + drho);
    wpold[bb] = (mpold[bb]*wsink + drho*wgrid) / (mpold[bb] + drho);
    mpold[bb] += drho;
    dmold[bb] += drho*pow(*dx,3);
    d[index]  -= drho;

  }

  /* For a 3D->1D index, put x+1, y+1, z+1 indices into nice variables */

  xo = 1;
  yo = *nx;
  zo = (*nx) * (*ny);

  /* Add stellar wind feedback */

  //  double msun = 1.989e33;
  //  double umass = (*d1)*pow(*x1,3)/msun;
  float StellarWindThresholdMass = 0.1;
  float StellarWindMomentumPerStellarMass = 5e6;
  float StellarWindEjectionFraction = 0.2;
  float fe = StellarWindEjectionFraction;
  float p_star = StellarWindMomentumPerStellarMass/(*v1);
  float m_wind = StellarWindThresholdMass/umass;
  float p_wind, m_cell;
#define MAX_SUPERCELL_NUMBER 1000
  int n_cell, ind_cell[MAX_SUPERCELL_NUMBER];
  float nx_cell[MAX_SUPERCELL_NUMBER], 
    ny_cell[MAX_SUPERCELL_NUMBER], nz_cell[MAX_SUPERCELL_NUMBER];
  float b_c, bx_c, by_c, bz_c, costheta = cos(3.1415926/3.9), v_wind, rho_wind,
    nx_b, ny_b, nz_b;
  FLOAT r_cell, x_cell, y_cell, z_cell;

  if (StellarWindFeedback == 1 && bx != NULL) { /*protostellar jets by B direction*/
    for (n = 0; n < nsinks; n++) {

      bb = sink_index[n];

      if (mpold[bb] < 0.0) continue;
      
      i = (xpold[bb] - *xstart)/(*dx);
      j = (ypold[bb] - *ystart)/(*dx);
      k = (zpold[bb] - *zstart)/(*dx);
      
      index = i + (j + k*(*ny))*(*nx);

      /* Decide whether mass increase reached the ejection threshold */

      if (dmold[bb] < m_wind) continue;
      if (mpold[bb]*pow(*dx,3)*umass < StellarWindTurnOnMass && (*t - tcpold[bb])*(*t1) < 1e5*3.1557e7) continue;

      int first = 0;
      if (dmold[bb] > 0.99*mpold[bb]*pow(*dx,3)) first = 1;
      

      /* Decide whether the current grid contains the whole supercell */

      if (i < *ibuff+2 || i > *nx-*ibuff-3 || 
	  j < *ibuff+2 || j > *ny-*ibuff-3 ||
	  k < *ibuff+2 || k > *nz-*ibuff-3) continue;


      /* Cacluate the total wind momentum */

      p_wind = p_star*sqrt(mpold[bb]*pow(*dx,3)*umass/10.0)*(1-fe)*dmold[bb];

      n_cell = 0;      
      m_cell = 0;
      
      /* Find the local B field direction */

      if (first == 1) {
	bx_c = bx[index], by_c = by[index], bz_c = bz[index];
	b_c = sqrt(pow(bx_c,2) + pow(by_c,2) + pow(bz_c,2));
	nx_b = bx_c/b_c, ny_b = by_c/b_c, nz_b = bz_c/b_c;
	printf("StellarWindFeedback = 1 l219, nx_b = %f, ny_b = %f, nz_b = %f\n",nx_b,ny_b,nz_b);	
	nx_jet[bb] = nx_b;
	ny_jet[bb] = ny_b;
	nz_jet[bb] = nz_b;
      } else {
	nx_b = nx_jet[bb];
	ny_b = ny_jet[bb];
	nz_b = nz_jet[bb];
      }

      /* Find the supercell and caclualte its total mass */

      for (int kk = -2; kk <= 2; kk++) {
	for (int jj = -2; jj <= 2; jj++) {
	  for (int ii = -2; ii <= 2; ii++) {
	    if (fabs(ii) != 2 && fabs(jj) != 2 && fabs(kk) != 2) continue;

	    x_cell = ii*(*dx), y_cell = jj*(*dx), z_cell = kk*(*dx);
	    r_cell = sqrt(pow(x_cell,2) + pow(y_cell,2) + pow(z_cell,2));	    

	    if (fabs((x_cell*nx_b + y_cell*ny_b + z_cell*nz_b)/r_cell) > costheta) { 
	      ind_cell[n_cell] = i+ii+(j+jj+(k+kk)*(*ny))*(*nx);
	      nx_cell[n_cell]  = x_cell / r_cell;
	      ny_cell[n_cell]  = y_cell / r_cell;
	      nz_cell[n_cell]  = z_cell / r_cell;
	      m_cell += d[ind_cell[n_cell]] * pow(*dx,3);
	      n_cell++;
	    }
	  }
	}
      }

      /* Calculate the feedback velocity */
      v_wind = p_wind / (m_cell + fe * dmold[bb]);

      /* Calculate the jet density */
      rho_wind = (m_cell + fe * dmold[bb]) / (n_cell * pow(*dx,3));

      /*printf("Wind injected: id=%"ISYM", vwind=%g, n_cell=%"ISYM", x=(%g, %g, %g), n=(%g,%g,%g,), ",
	     idold[bb], v_wind*(*v1), n_cell, xpold[bb], ypold[bb], zpold[bb], 
             nx_b, ny_b, nz_b);
      printf(" m_cell=%g, dm=%g, rho_wind=%g, p_wind=%g\n",
      m_cell*umass, dmold[bb]*umass, rho_wind*(*d1), p_wind);*/

      if (v_wind*(*v1) > 1e9) return FAIL;

      /* Do the feedback */
      for (int ic = 0; ic < n_cell; ic++) {

	d[ind_cell[ic]] = rho_wind;

	/* Fan-like jet */

	/*u[ind_cell[ic]] += nx_cell[ic]*v_wind;
	v[ind_cell[ic]] += ny_cell[ic]*v_wind;
	w[ind_cell[ic]] += nz_cell[ic]*v_wind;*/

	/* Parallel jet */
	
	int temp = sign(nx_cell[ic]*nx_b + ny_cell[ic]*ny_b + nz_cell[ic]*nz_b);
	te[ind_cell[ic]] -= 0.5*(pow(u[ind_cell[ic]],2) + pow(v[ind_cell[ic]],2) + pow(w[ind_cell[ic]],2));
	u[ind_cell[ic]] = temp*nx_b*v_wind;
	v[ind_cell[ic]] = temp*ny_b*v_wind;
	w[ind_cell[ic]] = temp*nz_b*v_wind;
	te[ind_cell[ic]] += 0.5*(pow(u[ind_cell[ic]],2) + pow(v[ind_cell[ic]],2) + pow(w[ind_cell[ic]],2));

      }

      /* Substract the ejected mass and set dm to be zero */

      dmold[bb] = 0.0;
      mpold[bb] -= fe*dmold[bb]/pow(*dx,3);
    }
  }

  if (StellarWindFeedback == 2 && bx == NULL) { /*protostellar jets by random direction*/
    for (n = 0; n < nsinks; n++) {
      //printf("StellarWindFeedback = 2 running\n");
      bb = sink_index[n];

      if (mpold[bb] < 0.0) continue;
      
      i = (xpold[bb] - *xstart)/(*dx);
      j = (ypold[bb] - *ystart)/(*dx);
      k = (zpold[bb] - *zstart)/(*dx);
      
      index = i + (j + k*(*ny))*(*nx);

      /* Decide whether mass increase reached the ejection threshold */

      if (dmold[bb] < m_wind)	continue;
      if (mpold[bb]*pow(*dx,3)*umass < StellarWindTurnOnMass && (*t - tcpold[bb])*(*t1) < 1e5*3.1557e7)	continue;
  
      
      int first = 0;
      if (dmold[bb] > 0.99*mpold[bb]*pow(*dx,3)) first = 1;
     
      
      /* Decide whether the current grid contains the whole supercell */

      if (i < *ibuff+2 || i > *nx-*ibuff-3 || 
	  j < *ibuff+2 || j > *ny-*ibuff-3 ||
	  k < *ibuff+2 || k > *nz-*ibuff-3) continue;


      /* Cacluate the total wind momentum */

      p_wind = p_star*sqrt(mpold[bb]*pow(*dx,3)*umass/10.0)*(1-fe)*dmold[bb];

      n_cell = 0;      
      m_cell = 0;
      
      /* Find the jet direction */

      if (first == 1) {
	//srand((unsigned)time(NULL));
	
	bx_c = rand()/332767.0, by_c = rand()/332767.0 , bz_c = rand()/332767.0 ;
	b_c = sqrt(pow(bx_c,2) + pow(by_c,2) + pow(bz_c,2));
	nx_b = bx_c/b_c, ny_b = by_c/b_c, nz_b = bz_c/b_c;
	nx_jet[bb] = nx_b;
	ny_jet[bb] = ny_b;
	nz_jet[bb] = nz_b;
      } else {
	nx_b = nx_jet[bb];
	ny_b = ny_jet[bb];
	nz_b = nz_jet[bb];
      }
      //printf("StellarWindFeedback = 2 l351\n");
      /* Find the supercell and caclualte its total mass */
  printf("%f\n",*nx_jet);
      for (int kk = -2; kk <= 2; kk++) {
	for (int jj = -2; jj <= 2; jj++) {
	  for (int ii = -2; ii <= 2; ii++) {
	    if (fabs(ii) != 2 && fabs(jj) != 2 && fabs(kk) != 2) continue;

	    x_cell = ii*(*dx), y_cell = jj*(*dx), z_cell = kk*(*dx);
	    r_cell = sqrt(pow(x_cell,2) + pow(y_cell,2) + pow(z_cell,2));	    

	    if (fabs((x_cell*nx_b + y_cell*ny_b + z_cell*nz_b)/r_cell) > costheta) { 
	      ind_cell[n_cell] = i+ii+(j+jj+(k+kk)*(*ny))*(*nx);
	      nx_cell[n_cell]  = x_cell / r_cell;
	      ny_cell[n_cell]  = y_cell / r_cell;
	      nz_cell[n_cell]  = z_cell / r_cell;
	      m_cell += d[ind_cell[n_cell]] * pow(*dx,3);
	      n_cell++;
	    }
	  }
	}
      }

      /* Calculate the feedback velocity */
      v_wind = p_wind / (m_cell + fe * dmold[bb]);

      /* Calculate the jet density */
      rho_wind = (m_cell + fe * dmold[bb]) / (n_cell * pow(*dx,3));

      /*printf("Wind injected: id=%"ISYM", vwind=%g, n_cell=%"ISYM", x=(%g, %g, %g), n=(%g,%g,%g,), ",
	     idold[bb], v_wind*(*v1), n_cell, xpold[bb], ypold[bb], zpold[bb], 
             nx_b, ny_b, nz_b);
      printf(" m_cell=%g, dm=%g, rho_wind=%g, p_wind=%g\n",
      m_cell*umass, dmold[bb]*umass, rho_wind*(*d1), p_wind);*/

      if (v_wind*(*v1) > 1e9) return FAIL;

      /* Do the feedback */
      for (int ic = 0; ic < n_cell; ic++) {

	d[ind_cell[ic]] = rho_wind;

	/* Fan-like jet */

	/*u[ind_cell[ic]] += nx_cell[ic]*v_wind;
	v[ind_cell[ic]] += ny_cell[ic]*v_wind;
	w[ind_cell[ic]] += nz_cell[ic]*v_wind;*/

	/* Parallel jet */
	
	int temp = sign(nx_cell[ic]*nx_b + ny_cell[ic]*ny_b + nz_cell[ic]*nz_b);
	te[ind_cell[ic]] -= 0.5*(pow(u[ind_cell[ic]],2) + pow(v[ind_cell[ic]],2) + pow(w[ind_cell[ic]],2));
	u[ind_cell[ic]] = temp*nx_b*v_wind;
	v[ind_cell[ic]] = temp*ny_b*v_wind;
	w[ind_cell[ic]] = temp*nz_b*v_wind;
	te[ind_cell[ic]] += 0.5*(pow(u[ind_cell[ic]],2) + pow(v[ind_cell[ic]],2) + pow(w[ind_cell[ic]],2));

      }

      /* Substract the ejected mass and set dm to be zero */

      dmold[bb] = 0.0;
      mpold[bb] -= fe*dmold[bb]/pow(*dx,3);
    }
  }


  /* StellarWindFeedback 3: Isotropic wind */

  float mdot_wind = 1e-5*(*dt)*(*t1)/(3.1557e7*umass);  /* 10^-5 solar mases per year - this is in code units: density x length^3*/
  v_wind = 2.0e6/(*v1);
  mdot_wind = mdot_wind/(4.0*Pi); /* mass Per solid angle */
  //printf("Adding Stellar wind 3: dt =%e, mdot =%e, Vwind =%e, rho_wind =%e \n",dt,mdot_wind*umass/(*t1),v_wind*(*v1),rho_wind*(*d1));
  FLOAT radius_cell[MAX_SUPERCELL_NUMBER]; 
  FLOAT radius2_cell[MAX_SUPERCELL_NUMBER];
  float SolidAngle;
  if (StellarWindFeedback == 3) {
  printf("mdotwind = %e\n",mdot_wind);
    // printf("STELLAR WIND FEEDBACK = 3\n");
    for (n = 0; n < nsinks; n++) {
      
      bb = sink_index[n];
      
      i = (xpold[bb] - *xstart)/(*dx);
      j = (ypold[bb] - *ystart)/(*dx);
      k = (zpold[bb] - *zstart)/(*dx);
      
      index = i + (j + k*(*ny))*(*nx);

      /* Decide whether the current grid contains the whole supercell */

      /*if (i < *ibuff+2 || i > *nx-*ibuff-3 || 
	  j < *ibuff+2 || j > *ny-*ibuff-3 ||
	  k < *ibuff+2 || k > *nz-*ibuff-3) continue;*/

      n_cell = 0;      
      m_cell = 0;
      
      /* Caclualte the total mass of the supercell */

      for (int kk = -3; kk <= 3; kk++) {
	for (int jj = -3; jj <= 3; jj++) {
	  for (int ii = -3; ii <= 3; ii++) {
	    if (ii == 0 && jj == 0 && kk == 0) continue;
	    /*if ((fabs(kk) == 2 && fabs(jj) == 2) ||
		(fabs(kk) == 2 && fabs(ii) == 2) ||
		(fabs(jj) == 2 && fabs(ii) == 2))
		continue;*/
	    if (fabs(ii) != 3 && fabs(jj) != 3 && fabs(kk) != 3) continue;
	    
	    x_cell = ii*(*dx), y_cell = jj*(*dx), z_cell = kk*(*dx);
	    double radius2 = pow(x_cell,2) + pow(y_cell,2) + pow(z_cell,2);
	    double radius = sqrt(radius2);
	    radius_cell[n_cell] = radius;
	    radius2_cell[n_cell] = pow(ii,2)+pow(jj,2)+pow(kk,2);
	    // printf("Radius Squared = %f\n",radius2_cell[n_cell]);
	    ind_cell[n_cell] = i+ii+(j+jj+(k+kk)*(*ny))*(*nx);
	    nx_cell[n_cell]  = x_cell/radius_cell[n_cell];
	    ny_cell[n_cell]  = y_cell/radius_cell[n_cell];
	    nz_cell[n_cell]  = z_cell/radius_cell[n_cell];
	    m_cell += d[ind_cell[n_cell]]*pow(*dx,3);
	    n_cell++;
	  
	  }
	}
      }

      /* Calculate the wind mass */

      // rho_wind = rho_wind/(n_cell*pow((*dx),3)) /* Density per solid angle in units of density code units */


      /* Do the feedback */
      float m_wind = 0.0;
      float cells_volume = 0.0;
      for (int ic = 0; ic < n_cell; ic++) {
	//v_wind = mdot_wind/(4.0*Pi*pow(radius_cell[ic],2)*rho_wind);
	//rho_wind = mdot_wind/(4.0*Pi*pow(radius_cell[ic],2)*v_wind);
	// m_wind += rho_wind*pow(*dx,3);;
	float u1 = u[ind_cell[ic]], 
	  v1 = v[ind_cell[ic]],
	  w1 = v[ind_cell[ic]];
	if (radius2_cell[ic] == 9.0) SolidAngle = 0.101796321; /*requires -3 -> 3 for suppercell */
	else if (radius2_cell[ic] == 10.0) SolidAngle = 0.0886801635;
	else if (radius2_cell[ic] == 11.0) SolidAngle = 0.0783382419;
	else if (radius2_cell[ic] == 13.0) SolidAngle = 0.0622518433;
	else if (radius2_cell[ic] == 14.0) SolidAngle = 0.0567666590;
	else if (radius2_cell[ic] == 17.0) SolidAngle = 0.0442920207;
	else if (radius2_cell[ic] == 18.0) SolidAngle = 0.0454334545;
	else if (radius2_cell[ic] == 19.0) SolidAngle = 0.0426303472;
	else if (radius2_cell[ic] == 22.0) SolidAngle = 0.0356113280;
	else if (radius2_cell[ic] == 27.0) SolidAngle = 0.0302870901;
	 else { 
	   SolidAngle = 4.*3.1415926/n_cell; 
	   //printf("star_maker8.C line 373: Radius squared is wrong?!? radius =%f, n_cell = %i\n",radius2_cell[ic],n_cell); 
	 }
	rho_wind = mdot_wind*SolidAngle/(pow((*dx),3));
	cells_volume += pow((*dx),3);
	//printf("rho_wind = %e cgs\n",rho_wind*(*d1));
	m_wind += rho_wind*pow(*dx,3);
	d[ind_cell[ic]] = rho_wind;
	u[ind_cell[ic]] = nx_cell[ic]*v_wind;
	v[ind_cell[ic]] = ny_cell[ic]*v_wind;
	w[ind_cell[ic]] = nz_cell[ic]*v_wind;

	te[ind_cell[ic]] -= 0.5*(u1*u1+v1*v1+w1*w1);
	te[ind_cell[ic]] += 0.5*v_wind*v_wind;
      }

      /* Substract the ejected mass and set dm to be zero */
      /*printf("Iso-wind injected: dt = %e s, vwind=%g, n_cell=%"ISYM", xp=(%g, %g, %g),m_star = %e, m_cell=%e, m_wind=%e, rho=%e, umass = %e\n", 
	     (*dt)*(*t1),v_wind*(*v1), n_cell, xpold[bb], ypold[bb], zpold[bb],mpold[bb]*umass , m_cell*umass,m_wind*umass, rho_wind*(*d1),umass);
      printf("Iso-wind injected: volume = %e code, %e cgs\n",cells_volume, cells_volume*pow((*x1),3));*/


      mpold[bb] += (m_cell - m_wind)/pow(*dx,3);

      /*dmold[bb] = 0.0;
	mpold[bb] -= fe*dmold[bb]/pow(*dx,3);*/

    }
  }


  /* Loop over grid looking for a cell with mass larger than massthres */


  if (*level == MaximumRefinementLevel) {
    float oldrho;
    float SinkCollapseDistance = SinkMergeDistance;
    for (k = *ibuff; k < *nz-*ibuff; k++) {
      for (j = *ibuff; j < *ny-*ibuff; j++) {
	index = (k * (*ny) + j) * (*nx) + (*ibuff);
	for (i = *ibuff; i < *nx-*ibuff; i++, index++) {

	  /* Finest level of refinement and density greater than threshold? */
	  
	  if (*jlrefine > 0)
	    jeansthresh = jlsquared * temp[index] / d[index];


	  if (r[index] == 0 && (d[index] > densthresh ||
				(*jlrefine > 0 && dx2 > jeansthresh))) {
	  printf("star_maker8: density above thresh-hold - will now make a new star?!\n");
	    
	    xpos = *xstart + ((float) i - 0.5)*(*dx);
	    ypos = *ystart + ((float) j - 0.5)*(*dx);
	    zpos = *zstart + ((float) k - 0.5)*(*dx);
	  
	    /* Calculate change in density */

	    if (*jlrefine > 0)
	      maxdens = min(jlsquared * temp[index] / dx2, densthresh);
	    else
	      maxdens = densthresh;
	    oldrho = d[index];
	    adddens = d[index] - maxdens;
	    
	    /* Remove mass from grid */
	    
	    d[index] = maxdens;

	    /* Get velocity at grid center (averaged over cells ... ONLY
	       for PPM) */

	    if (*imethod == 2) {
	      ugrid = 0.5*(u[index] + u[index+xo]);
	      vgrid = 0.5*(v[index] + v[index+yo]);
	      wgrid = 0.5*(w[index] + w[index+zo]);
	    } else {
	      total_density = d[index] + d[index-xo] + d[index+xo] + d[index-yo] +
		d[index+yo] + d[index-zo] + d[index+zo];

	      ugrid = (u[index]*d[index] + 
		       u[index-xo]*d[index-xo] + u[index+xo]*d[index+xo] +
		       u[index-yo]*d[index-yo] + u[index+yo]*d[index+yo] +
		       u[index-zo]*d[index-zo] + u[index+zo]*d[index+zo]) /
		total_density;

	      vgrid = (v[index]*d[index] + 
		       v[index-xo]*d[index-xo] + v[index+xo]*d[index+xo] +
		       v[index-yo]*d[index-yo] + v[index+yo]*d[index+yo] +
		       v[index-zo]*d[index-zo] + v[index+zo]*d[index+zo]) /
		total_density;

	      wgrid = (w[index]*d[index] + 
		       w[index-xo]*d[index-xo] + w[index+xo]*d[index+xo] +
		       w[index-yo]*d[index-yo] + w[index+yo]*d[index+yo] +
		       w[index-zo]*d[index-zo] + w[index+zo]*d[index+zo]) /
		total_density;
	    }

	    /* Look for a nearby OLD sink particle to add the mass to */
	    
	    inew = 1;
	    nearestdx2 = 1e20;
	    for (cc = 0; cc < nsinks; cc++) {
	      
	      n = sink_index[cc];
	      
	      delx = xpos - xpold[n];
	      dely = ypos - ypold[n];
	      delz = zpos - zpold[n];
	      dist2 = delx*delx + dely*dely + delz*delz;

	      /* If sink is within 5 cells and closest one, then add to it */
	      //printf("star_maker8:  distance: dist=%"FSYM", SCD =%"FSYM" \n",pow(dist2,0.5),SinkCollapseDistance );	      
	      if (dist2 < pow(SinkCollapseDistance,2) && dist2 < nearestdx2) {
		nearestdx2 = dist2;
		closest = n;
	      }

	    } // ENDFOR old particles
	    printf("star_maker8: nearest old star = %"FSYM"\n",pow(nearestdx2,0.5) );

	    /* Add momentum and mass to nearest OLD sink */

	    if (nearestdx2 < 1) {
	    
	      upold[closest] = (upold[closest] * mpold[closest] + ugrid*adddens) /
		(mpold[closest] + adddens);
	      vpold[closest] = (vpold[closest] * mpold[closest] + vgrid*adddens) /
		(mpold[closest] + adddens);
	      wpold[closest] = (wpold[closest] * mpold[closest] + wgrid*adddens) /
		(mpold[closest] + adddens);
	      mpold[closest] = mpold[closest] + adddens;
	      dmold[closest] += adddens*pow(*dx,3);
	    
	      /* Record that a new particle is not needed */
	      
	      inew = 0;
	      printf("star_maker8:  new star not needed \n" );	      
	    }  // ENDIF add to particle

	    /* Now look for nearby NEW sinks */

	    nearestdx2 = 1e20;
	    for (n = 0; n < ii; n++) {
	      
	      delx = xpos - xp[n];
	      dely = ypos - yp[n];
	      delz = zpos - zp[n];
	      dist2 = delx*delx + dely*dely + delz*delz;

	      /* If sink is within SinkCollapseDistance, then add to it */
	      
	      if (dist2 < pow(SinkCollapseDistance,2) && dist2 < nearestdx2) {
		nearestdx2 = dist2;
		closest = n;
	      }
	      
	    } // ENDFOR new particles
	  
	    /* Add momentum and then mass to NEW sink*/

	    printf("star_maker8: nearest new star = %"FSYM"\n",pow(nearestdx2,0.5) );
	    if (nearestdx2 < 1) {

	      up[closest] = (up[closest] * mp[closest] + ugrid*adddens) /
		(mp[closest] + adddens);
	      vp[closest] = (vp[closest] * mp[closest] + vgrid*adddens) /
		(mp[closest] + adddens);
	      wp[closest] = (wp[closest] * mp[closest] + wgrid*adddens) /
		(mp[closest] + adddens);
	      mp[closest] = mp[closest] + adddens;
	      dm[closest] += adddens*pow(*dx,3);
	    
	      /* Record that a new particle is not needed */
	      
	      inew = 0;
	      printf("star_maker8:  new star not needed \n" );
	    } // ENDIF add to new particle

	    /* Create a new sink particle if necessary and if there's room */
	    
	    if (inew == 1 && ii < *nmax) {
	      
	    printf("star_maker8: making new star\n" );
	      mp[ii] = adddens;
	      type[ii] = *ctype;
	      
	      /* Set positions and velocities */
	    
	      xp[ii] = xpos;
	      yp[ii] = ypos;
	      zp[ii] = zpos;
	      up[ii] = ugrid;
	      vp[ii] = vgrid;
	      wp[ii] = wgrid;

	      /* Set creation time */
	      
	      tcp[ii] = (float) *t;
	      tdp[ii] = 0.0;
	      dm[ii]  = adddens*pow(*dx,3);

	      ii++;

	    } // ENDIF create a new sink
	    
	  } // ENDIF make sink particle

	} // ENDFOR i
      } // ENDFOR j
    } // ENDFOR k

  } // if (level == maxlevel)

  if (ii > 0)
    printf("P(%"ISYM"): star_maker8[add]: %"ISYM" new sink particles\n", *nproc, ii);

  if (ii >= *nmax) {
    fprintf(stdout, "star_maker8: reached max new particle count");
    return FAIL;
  }

  delete sink_index;

  *np = ii;
  return SUCCESS;

}
