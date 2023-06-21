
#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h" 
void esys_roe_adb_mhd(const float d, const float v1, const float v2, const float v3,
		      const float h, const float b1, const float b2, const float b3, 
		      const float x, const float y,
		      float eigenvalues[],
		      float right_eigenmatrix[][7], float left_eigenmatrix[][7])
{
  
  float Gamma_1 = Gamma - 1, Gamma_2 = Gamma - 2, ONE_OVER_SQRT2 = 1.0/sqrt(2.0);
  
  float vsq,btsq,bt_starsq,vaxsq,hp,twid_asq,q_starsq,cfsq,cf,cssq,cs;
  float bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,vbet,alpha_f,alpha_s;
  float sqrtd,s,twid_a,qf,qs,af_prime,as_prime,afpbb,aspbb,vax;
  float norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
  vsq = v1*v1 + v2*v2 + v3*v3;
  btsq = b2*b2 + b3*b3;
  bt_starsq = (Gamma_1 - Gamma_2*y)*btsq;

/* Compute fast- and slow-magnetosonic speeds (eqs. ) */

  vaxsq = b1*b1/d;
  hp = h - (vaxsq + btsq/d);
  twid_asq = max((Gamma_1*(hp-0.5*vsq)-Gamma_2*x), tiny_number); 
  q_starsq  = twid_asq + (vaxsq + bt_starsq/d);
  cfsq = 0.5*(q_starsq + sqrt(q_starsq*q_starsq - 4.0*twid_asq*vaxsq));
  cf  = sqrt(cfsq);
  cssq = 0.5*(q_starsq - sqrt(q_starsq*q_starsq - 4.0*twid_asq*vaxsq));
  if (cssq <= 0.0) {
    cssq = 0.0;
  }
  cs  = sqrt(cssq);

/* Compute beta(s) (eqs. ) */

  bt = sqrt(btsq);
  bt_star = sqrt(bt_starsq);
  if (bt == 0.0) {
    bet2 = ONE_OVER_SQRT2 ;
    bet3 = ONE_OVER_SQRT2 ;
    bet2_star = ONE_OVER_SQRT2 ;
    bet3_star = ONE_OVER_SQRT2 ;
  } else {
    bet2 = b2/bt;
    bet3 = b3/bt;
    bet2_star = b2/bt_star;
    bet3_star = b3/bt_star;
  }
  bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;
  vbet = v2*bet2_star + v3*bet3_star;

/* Compute alpha(s) (eq. ) */

  if ((cfsq-cssq) == 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else if ( (twid_asq - cssq) <= 0.0) {
    alpha_f = 0.0;
    alpha_s = 1.0;
  } else if ( (cfsq - twid_asq) <= 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else {
    alpha_f = sqrt((twid_asq - cssq)/(cfsq - cssq));
    alpha_s = sqrt((cfsq - twid_asq)/(cfsq - cssq));
  }

/* Compute Q(s) and A(s) (eq. ), etc. */

  sqrtd = sqrt(d);
  s = sign(b1);
  twid_a = sqrt(twid_asq);
  qf = cf*alpha_f*s;
  qs = cs*alpha_s*s;
  af_prime = twid_a*alpha_f/sqrtd;
  as_prime = twid_a*alpha_s/sqrtd;
  afpbb = af_prime*bt_star*bet_starsq;
  aspbb = as_prime*bt_star*bet_starsq;

/* Compute eigenvalues (eq. ) */

  vax = sqrt(vaxsq);
  eigenvalues[0] = v1 - cf;
  eigenvalues[1] = v1 - vax;
  eigenvalues[2] = v1 - cs;
  eigenvalues[3] = v1;
  eigenvalues[4] = v1 + cs;
  eigenvalues[5] = v1 + vax;
  eigenvalues[6] = v1 + cf;
  if (right_eigenmatrix == NULL || left_eigenmatrix == NULL) return;

/* Right-eigenvectors, stored as COLUMNS (eq. ) */

  right_eigenmatrix[0][0] = alpha_f;
  right_eigenmatrix[1][0] = alpha_f*(v1 - cf);
  right_eigenmatrix[2][0] = alpha_f*v2 + qs*bet2_star;
  right_eigenmatrix[3][0] = alpha_f*v3 + qs*bet3_star;
  right_eigenmatrix[4][0] = alpha_f*(hp - v1*cf) + qs*vbet + aspbb;
  right_eigenmatrix[5][0] = as_prime*bet2_star;
  right_eigenmatrix[6][0] = as_prime*bet3_star;

/*right_eigenmatrix[0][1] = 0.0; */
/*right_eigenmatrix[1][1] = 0.0; */
  right_eigenmatrix[2][1] = -bet3;
  right_eigenmatrix[3][1] = bet2;
  right_eigenmatrix[4][1] = -(v2*bet3 - v3*bet2);
  right_eigenmatrix[5][1] = -bet3*s/sqrtd;
  right_eigenmatrix[6][1] = bet2*s/sqrtd;

  right_eigenmatrix[0][2] = alpha_s;
  right_eigenmatrix[1][2] = alpha_s*(v1 - cs);
  right_eigenmatrix[2][2] = alpha_s*v2 - qf*bet2_star;
  right_eigenmatrix[3][2] = alpha_s*v3 - qf*bet3_star;
  right_eigenmatrix[4][2] = alpha_s*(hp - v1*cs) - qf*vbet - afpbb;
  right_eigenmatrix[5][2] = -af_prime*bet2_star;
  right_eigenmatrix[6][2] = -af_prime*bet3_star;

  right_eigenmatrix[0][3] = 1.0;
  right_eigenmatrix[1][3] = v1;
  right_eigenmatrix[2][3] = v2;
  right_eigenmatrix[3][3] = v3;
  right_eigenmatrix[4][3] = 0.5*vsq + Gamma_2*x/Gamma_1;
/*right_eigenmatrix[5][3] = 0.0; */
/*right_eigenmatrix[6][3] = 0.0; */

  right_eigenmatrix[0][4] = alpha_s;
  right_eigenmatrix[1][4] = alpha_s*(v1 + cs);
  right_eigenmatrix[2][4] = alpha_s*v2 + qf*bet2_star;
  right_eigenmatrix[3][4] = alpha_s*v3 + qf*bet3_star;
  right_eigenmatrix[4][4] = alpha_s*(hp + v1*cs) + qf*vbet - afpbb;
  right_eigenmatrix[5][4] = right_eigenmatrix[5][2];
  right_eigenmatrix[6][4] = right_eigenmatrix[6][2];

/*right_eigenmatrix[0][5] = 0.0; */
/*right_eigenmatrix[1][5] = 0.0; */
  right_eigenmatrix[2][5] = bet3;
  right_eigenmatrix[3][5] = -bet2;
  right_eigenmatrix[4][5] = -right_eigenmatrix[4][1];
  right_eigenmatrix[5][5] = right_eigenmatrix[5][1];
  right_eigenmatrix[6][5] = right_eigenmatrix[6][1];

  right_eigenmatrix[0][6] = alpha_f;
  right_eigenmatrix[1][6] = alpha_f*(v1 + cf);
  right_eigenmatrix[2][6] = alpha_f*v2 - qs*bet2_star;
  right_eigenmatrix[3][6] = alpha_f*v3 - qs*bet3_star;
  right_eigenmatrix[4][6] = alpha_f*(hp + v1*cf) - qs*vbet + aspbb;
  right_eigenmatrix[5][6] = right_eigenmatrix[5][0];
  right_eigenmatrix[6][6] = right_eigenmatrix[6][0];

/* Left-eigenvectors, stored as ROWS (eq. ) */

/* Normalize by 1/2a^{2}: quantities denoted by \hat{f} */
  norm = 0.5/twid_asq;
  cff = norm*alpha_f*cf;
  css = norm*alpha_s*cs;
  qf *= norm;
  qs *= norm;
  af = norm*af_prime*d;
  as = norm*as_prime*d;
  afpb = norm*af_prime*bt_star;
  aspb = norm*as_prime*bt_star;

/* Normalize by (gamma-1)/2a^{2}: quantities denoted by \bar{f} */
  norm *= Gamma_1;
  alpha_f *= norm;
  alpha_s *= norm;
  q2_star = bet2_star/bet_starsq;
  q3_star = bet3_star/bet_starsq;
  vqstr = (v2*q2_star + v3*q3_star);
  norm *= 2.0;
  //kludge: fucking with these vectors.
  left_eigenmatrix[0][0] = alpha_f*(vsq-hp) + cff*(cf+v1) - qs*vqstr - aspb;
  left_eigenmatrix[0][1] = -alpha_f*v1 - cff;
  left_eigenmatrix[0][2] = -alpha_f*v2 + qs*q2_star;
  left_eigenmatrix[0][3] = -alpha_f*v3 + qs*q3_star;
  left_eigenmatrix[0][4] = alpha_f;
  left_eigenmatrix[0][5] = as*q2_star - alpha_f*b2;
  left_eigenmatrix[0][6] = as*q3_star - alpha_f*b3;

  left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
/*left_eigenmatrix[1][1] = 0.0; */
  left_eigenmatrix[1][2] = -0.5*bet3;
  left_eigenmatrix[1][3] = 0.5*bet2;
/*left_eigenmatrix[1][4] = 0.0; */
  left_eigenmatrix[1][5] = -0.5*sqrtd*bet3*s;
  left_eigenmatrix[1][6] = 0.5*sqrtd*bet2*s;

  left_eigenmatrix[2][0] = alpha_s*(vsq-hp) + css*(cs+v1) + qf*vqstr + afpb;
  left_eigenmatrix[2][1] = -alpha_s*v1 - css;
  left_eigenmatrix[2][2] = -alpha_s*v2 - qf*q2_star;
  left_eigenmatrix[2][3] = -alpha_s*v3 - qf*q3_star;
  left_eigenmatrix[2][4] = alpha_s;
  left_eigenmatrix[2][5] = -af*q2_star - alpha_s*b2;
  left_eigenmatrix[2][6] = -af*q3_star - alpha_s*b3;

  left_eigenmatrix[3][0] = 1.0 - norm*(0.5*vsq - Gamma_2*x/Gamma_1); 
  left_eigenmatrix[3][1] = norm*v1;
  left_eigenmatrix[3][2] = norm*v2;
  left_eigenmatrix[3][3] = norm*v3;
  left_eigenmatrix[3][4] = -norm;
  left_eigenmatrix[3][5] = norm*b2;
  left_eigenmatrix[3][6] = norm*b3;

  left_eigenmatrix[4][0] = alpha_s*(vsq-hp) + css*(cs-v1) - qf*vqstr + afpb;
  left_eigenmatrix[4][1] = -alpha_s*v1 + css;
  left_eigenmatrix[4][2] = -alpha_s*v2 + qf*q2_star;
  left_eigenmatrix[4][3] = -alpha_s*v3 + qf*q3_star;
  left_eigenmatrix[4][4] = alpha_s;
  left_eigenmatrix[4][5] = left_eigenmatrix[2][5];
  left_eigenmatrix[4][6] = left_eigenmatrix[2][6];

  left_eigenmatrix[5][0] = -left_eigenmatrix[1][0];
/*left_eigenmatrix[5][1] = 0.0; */
  left_eigenmatrix[5][2] = -left_eigenmatrix[1][2];
  left_eigenmatrix[5][3] = -left_eigenmatrix[1][3];
/*left_eigenmatrix[5][4] = 0.0; */
  left_eigenmatrix[5][5] = left_eigenmatrix[1][5];
  left_eigenmatrix[5][6] = left_eigenmatrix[1][6];

  left_eigenmatrix[6][0] = alpha_f*(vsq-hp) + cff*(cf-v1) + qs*vqstr - aspb;
  left_eigenmatrix[6][1] = -alpha_f*v1 + cff;
  left_eigenmatrix[6][2] = -alpha_f*v2 - qs*q2_star;
  left_eigenmatrix[6][3] = -alpha_f*v3 - qs*q3_star;
  left_eigenmatrix[6][4] = alpha_f;
  left_eigenmatrix[6][5] = left_eigenmatrix[0][5];
  left_eigenmatrix[6][6] = left_eigenmatrix[0][6];
}
 
 
 
/*------------------------------ ISOTHERMAL MHD --------------------------------
 *
 * Input: d,v1,v2,v3,b1,b2,b3 = Roe averaged density, velocities, and B field
 *        x,y = numerical factors (eqs. )
 * Output: eigenvalues[6], right_eigenmatrix[6,6], left_eigenmatrix[6,6];
 */
 
void esys_roe_iso_mhd(const float d, const float v1, const float v2, const float v3,
		      const float Iso_csound2, const float b1, const float b2, const float b3, 
		      const float x, const float y, 
		      float eigenvalues[],
		      float right_eigenmatrix[][7], float left_eigenmatrix[][7])
   {
  float ONE_OVER_SQRT2 = 1.0/sqrt(2.0);
  float btsq,bt_starsq,vaxsq,twid_csq,q_starsq,cfsq,cf,cssq,cs;
  float bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,alpha_f,alpha_s;
  float sqrtd,s,twid_c,qf,qs,af_prime,as_prime,vax;
  float norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
  btsq = b2*b2 + b3*b3;
  bt_starsq = btsq*y;

/* Compute fast- and slow-magnetosonic speeds (eq. ) */

  vaxsq = b1*b1/d;
  twid_csq = Iso_csound2 + x;
  q_starsq  = twid_csq + (vaxsq + bt_starsq/d);
  cfsq = 0.5*(q_starsq + sqrt(q_starsq*q_starsq - 4.0*twid_csq*vaxsq));
  cf  = sqrt(cfsq);
  cssq = 0.5*(q_starsq - sqrt(q_starsq*q_starsq - 4.0*twid_csq*vaxsq));
  if (cssq <= 0.0) {
    cssq = 0.0;
  }
  cs  = sqrt(cssq);

/* Compute beta's (eqs. ) */

  bt = sqrt(btsq);
  bt_star = sqrt(bt_starsq);
  if (bt == 0.0) {
    bet2 = ONE_OVER_SQRT2 ;
    bet3 = ONE_OVER_SQRT2 ;
    bet2_star = ONE_OVER_SQRT2 ;
    bet3_star = ONE_OVER_SQRT2 ;
  } 
  else {
    bet2 = b2/bt;
    bet3 = b3/bt;
    bet2_star = b2/bt_star;
    bet3_star = b3/bt_star;
  }
  bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;

/* Compute alpha's (eq. ) */

  if ((cfsq-cssq) == 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else if ((twid_csq - cssq) <= 0.0) {
    alpha_f = 0.0;
    alpha_s = 1.0;
  } else if ((cfsq - twid_csq) <= 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else {
    alpha_f = sqrt((twid_csq - cssq)/(cfsq - cssq));
    alpha_s = sqrt((cfsq - twid_csq)/(cfsq - cssq));
  }

/* Compute Q's (eq. ), etc. */

  sqrtd = sqrt(d);
  s = sign(b1);
  twid_c = sqrt(twid_csq);
  qf = cf*alpha_f*s;
  qs = cs*alpha_s*s;
  af_prime = twid_c*alpha_f/sqrtd;
  as_prime = twid_c*alpha_s/sqrtd;

/* Compute eigenvalues (eq. ) */

  vax  = sqrt(vaxsq);
  eigenvalues[0] = v1 - cf;
  eigenvalues[1] = v1 - vax;
  eigenvalues[2] = v1 - cs;
  eigenvalues[3] = v1 + cs;
  eigenvalues[4] = v1 + vax;
  eigenvalues[5] = v1 + cf;
  if (right_eigenmatrix == NULL || left_eigenmatrix == NULL) return;

/* Right-eigenvectors, stored as COLUMNS (eq. ) */

  right_eigenmatrix[0][0] = alpha_f;
  right_eigenmatrix[1][0] = alpha_f*(v1 - cf);
  right_eigenmatrix[2][0] = alpha_f*v2 + qs*bet2_star;
  right_eigenmatrix[3][0] = alpha_f*v3 + qs*bet3_star;
  right_eigenmatrix[4][0] = as_prime*bet2_star;
  right_eigenmatrix[5][0] = as_prime*bet3_star;

/*right_eigenmatrix[0][1] = 0.0; */
/*right_eigenmatrix[1][1] = 0.0; */
  right_eigenmatrix[2][1] = -bet3;
  right_eigenmatrix[3][1] = bet2;
  right_eigenmatrix[4][1] = -bet3*s/sqrtd;
  right_eigenmatrix[5][1] = bet2*s/sqrtd;

  right_eigenmatrix[0][2] = alpha_s;
  right_eigenmatrix[1][2] = alpha_s*(v1 - cs);
  right_eigenmatrix[2][2] = alpha_s*v2 - qf*bet2_star;
  right_eigenmatrix[3][2] = alpha_s*v3 - qf*bet3_star;
  right_eigenmatrix[4][2] = -af_prime*bet2_star;
  right_eigenmatrix[5][2] = -af_prime*bet3_star;

  right_eigenmatrix[0][3] = alpha_s;
  right_eigenmatrix[1][3] = alpha_s*(v1 + cs);
  right_eigenmatrix[2][3] = alpha_s*v2 + qf*bet2_star;
  right_eigenmatrix[3][3] = alpha_s*v3 + qf*bet3_star;
  right_eigenmatrix[4][3] = right_eigenmatrix[4][2];
  right_eigenmatrix[5][3] = right_eigenmatrix[5][2];

/*right_eigenmatrix[0][4] = 0.0; */
/*right_eigenmatrix[1][4] = 0.0; */
  right_eigenmatrix[2][4] = bet3;
  right_eigenmatrix[3][4] = -bet2;
  right_eigenmatrix[4][4] = right_eigenmatrix[4][1];
  right_eigenmatrix[5][4] = right_eigenmatrix[5][1];

  right_eigenmatrix[0][5] = alpha_f;
  right_eigenmatrix[1][5] = alpha_f*(v1 + cf);
  right_eigenmatrix[2][5] = alpha_f*v2 - qs*bet2_star;
  right_eigenmatrix[3][5] = alpha_f*v3 - qs*bet3_star;
  right_eigenmatrix[4][5] = right_eigenmatrix[4][0];
  right_eigenmatrix[5][5] = right_eigenmatrix[5][0];

/* Left-eigenvectors, stored as ROWS (eq. ) */
/* Normalize by 1/2a^{2}: quantities denoted by \hat{f} */

  norm = 0.5/twid_csq;
  cff = norm*alpha_f*cf;
  css = norm*alpha_s*cs;
  qf *= norm;
  qs *= norm;
  af = norm*af_prime*d;
  as = norm*as_prime*d;
  afpb = norm*af_prime*bt_star;
  aspb = norm*as_prime*bt_star;

  q2_star = bet2_star/bet_starsq;
  q3_star = bet3_star/bet_starsq;
  vqstr = (v2*q2_star + v3*q3_star);

  left_eigenmatrix[0][0] = cff*(cf+v1) - qs*vqstr - aspb;
  left_eigenmatrix[0][1] = -cff;
  left_eigenmatrix[0][2] = qs*q2_star;
  left_eigenmatrix[0][3] = qs*q3_star;
  left_eigenmatrix[0][4] = as*q2_star;
  left_eigenmatrix[0][5] = as*q3_star;

  left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
/*left_eigenmatrix[1][1] = 0.0; */
  left_eigenmatrix[1][2] = -0.5*bet3;
  left_eigenmatrix[1][3] = 0.5*bet2;
  left_eigenmatrix[1][4] = -0.5*sqrtd*bet3*s;
  left_eigenmatrix[1][5] = 0.5*sqrtd*bet2*s;

  left_eigenmatrix[2][0] = css*(cs+v1) + qf*vqstr + afpb;
  left_eigenmatrix[2][1] = -css;
  left_eigenmatrix[2][2] = -qf*q2_star;
  left_eigenmatrix[2][3] = -qf*q3_star;
  left_eigenmatrix[2][4] = -af*q2_star;
  left_eigenmatrix[2][5] = -af*q3_star;

  left_eigenmatrix[3][0] = css*(cs-v1) - qf*vqstr + afpb;
  left_eigenmatrix[3][1] = css;
  left_eigenmatrix[3][2] = -left_eigenmatrix[2][2];
  left_eigenmatrix[3][3] = -left_eigenmatrix[2][3];
  left_eigenmatrix[3][4] = left_eigenmatrix[2][4];
  left_eigenmatrix[3][5] = left_eigenmatrix[2][5];

  left_eigenmatrix[4][0] = -left_eigenmatrix[1][0];
/*left_eigenmatrix[4][1] = 0.0; */
  left_eigenmatrix[4][2] = -left_eigenmatrix[1][2];
  left_eigenmatrix[4][3] = -left_eigenmatrix[1][3];
  left_eigenmatrix[4][4] = left_eigenmatrix[1][4];
  left_eigenmatrix[4][5] = left_eigenmatrix[1][5];

  left_eigenmatrix[5][0] = cff*(cf-v1) + qs*vqstr - aspb;
  left_eigenmatrix[5][1] = cff;
  left_eigenmatrix[5][2] = -left_eigenmatrix[0][2];
  left_eigenmatrix[5][3] = -left_eigenmatrix[0][3];
  left_eigenmatrix[5][4] = left_eigenmatrix[0][4];
  left_eigenmatrix[5][5] = left_eigenmatrix[0][5];
}
