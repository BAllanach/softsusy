/* Functions needed for 1- and 2-loop scalar self energies */

/* Conventions: 
   If function does *not* require TSIL evaluation, then it takes
   x,y,z,... args, plus s and qq, and returns the function value.  If
   TSIL is required, it takes one or more TSIL_RESULTs as args and
   passes a collection of results back via pointers.  In such cases
   the first TSIL_RESUL in the arg list contains the function
   parameters in the correct order. What is appropriate for the others
   depends on circumstances.
*/

#include "supermodel.h"
#include "self_scalar.h"

/* ---------------------------------------------------------------- */
/* Some functions and wrappers (notations as in hep-ph/0312092)     */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq. 3.2 */

SUMO_COMPLEX A_S (SUMO_REAL x, TSIL_REAL qq)
{
  return TSIL_A (x, qq);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq. 3.3 */

SUMO_COMPLEX B_SS (SUMO_REAL x, SUMO_REAL y, SUMO_REAL s, SUMO_REAL qq)
{
  return -TSIL_B (x, y, s, qq);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Derivative of B_SS with respect to first mass^2 argument x.  */

SUMO_COMPLEX B_SpS (SUMO_REAL x, SUMO_REAL y, SUMO_REAL s, SUMO_REAL qq)
{
  return -SUMO_Bp (x, y, s, qq, interp);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs. 3.5 and 3.5 */

int bFF (TSIL_REAL x, TSIL_REAL y, TSIL_REAL s, TSIL_REAL qq, 
	 TSIL_COMPLEX *BFF, TSIL_COMPLEX *Bff) 
{
  TSIL_COMPLEX Bxys; 

  Bxys = TSIL_B(x,y,s,qq);

  *BFF = (x + y - s)*Bxys - TSIL_A(x,qq) - TSIL_A(y,qq);
  *Bff = 2.0L * Bxys;

  return 0;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Derivative of BFF and Bff with respect to first mass^2 argument
   x: */

int bFpF (TSIL_REAL x, TSIL_REAL y, TSIL_REAL s, TSIL_REAL qq, 
	 TSIL_COMPLEX *BFpF, TSIL_COMPLEX *Bfpf) 
{
  TSIL_COMPLEX Bxys, Bxpys; 

  Bxys = TSIL_B(x,y,s,qq);
  Bxpys = SUMO_Bp(x,y,s,qq,interp);

  *BFpF = (x + y - s)*Bxpys + Bxys - TSIL_Ap(x,qq);
  *Bfpf = 2.0L * Bxpys;

  return 0;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq. 3.6 */

SUMO_COMPLEX A_V (SUMO_REAL x, TSIL_REAL qq)
{
  return 3.0L * TSIL_A (x, qq);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Derivative of A_V with respect to mass^2 argument.               */

SUMO_COMPLEX A_Vp (SUMO_REAL x, TSIL_REAL qq)
{
  return 3.0L * TSIL_Ap (x, qq);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq. 3.7 
   NOTE: the imaginary part of B(x,0,s,qq) is removed when y>0, 
   because it corresponds to an absorptive part that should be 
   cancelled by the corresponding Goldstone contributions in B_SS
   and B_VV. The real part of those integrals is also taken. That 
   way the unphysical absorptive parts are dropped in all parts of 
   the calculation. Otherwise, imperfect cancellation would ensue.  
*/

SUMO_COMPLEX B_SV (SUMO_REAL x,SUMO_REAL y,TSIL_REAL s,TSIL_REAL qq)
{
  SUMO_COMPLEX result;

  if (TSIL_FABS(s/(x+y+s+qq)) < TSIL_TOL) return 0.0L;

  if (TSIL_FABS(y/(x+y+s)) < TSIL_TOL) {
    result = 3.0L * ((x+s) * TSIL_B(x,0,s,qq) + TSIL_A(x,qq)) - 2.0L * s; 
  } else {
    result = (2.0L*x - y + 2.0L*s)*TSIL_B(x,y,s,qq) +
      (x - s)*(x - s)*(TSIL_CREAL(TSIL_B(x,0,s,qq)) - TSIL_B(x,y,s,qq))/y
      + TSIL_A(x,qq) + (x-y-s)*TSIL_A(y,qq)/y;
  }

  return result;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Derivative of B_SV with respect to first argument x.             */

SUMO_COMPLEX B_SpV (SUMO_REAL x,SUMO_REAL y,TSIL_REAL s,TSIL_REAL qq)
{
  if (TSIL_FABS(s/(x+y+s+qq)) < 100.0*TSIL_TOL) return(0.0L);
  
  if (TSIL_FABS(y/(x+y+s)) < 100.0*TSIL_TOL) 
    return (3.0L * ((x+s) * SUMO_Bp(x,0,s,qq,interp) + TSIL_B(x,0,s,qq) + 
             TSIL_Ap(x,qq))); 
  
  if (TSIL_FABS((x-s)/(x+y+s)) < 100.0*TSIL_TOL)
    return ( (4.0L*x - y)*SUMO_Bp(x,y,s,qq,interp) + 2.0L*TSIL_B(x,y,s,qq)
             + TSIL_Ap(x,qq) + TSIL_A(y,qq)/y);

  return ((2.0L*x - y + 2.0L*s)*SUMO_Bp(x,y,s,qq,interp) + 2.0L*TSIL_B(x,y,s,qq) 
	  + (x - s)*(x - s)*(TSIL_CREAL(SUMO_Bp(x,0,s,qq,interp)) - SUMO_Bp(x,y,s,qq,interp))/y
             + 2.0L*(x - s)*(TSIL_CREAL(TSIL_B(x,0,s,qq)) - TSIL_B(x,y,s,qq))/y
             + TSIL_Ap(x,qq) + TSIL_A(y,qq)/y);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Derivative of B_SV with respect to second argument.              */

SUMO_COMPLEX B_SVp (SUMO_REAL x,SUMO_REAL y,TSIL_REAL s,TSIL_REAL qq)
{
  SUMO_COMPLEX result;

  if (TSIL_FABS(s/(x+y+s+qq)) < TSIL_TOL) result = 0.0L;
  
  result = (2.0L*x - y + 2.0L*s - (x - s)*(x - s)/y)*SUMO_Bp(y,x,s,qq,interp) 
           + ((x - s)*(x - s)/(y*y) - 1.0) * TSIL_B(x,y,s,qq) 
           - ((x - s)*(x - s)/(y*y)) * TSIL_CREAL(TSIL_B(x,0,s,qq))
           + (x-s)/y - TSIL_Ap (y,qq);

  return result;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq. 3.8 (Landau gauge) 
   NOTE: the imaginary parts of B(0,0,s,qq), B(x,0,s,qq), and
   B(0,y,s,qq) are removed, because they correspond to absorptive 
   parts that should be cancelled by the corresponding Goldstone 
   contributions in B_SS and B_SV. The real part of those integrals 
   is also taken. That way the unphysical absorptive parts are 
   dropped in all parts of the calculation. Otherwise, imperfect 
   cancellation of the unphysical absorptive parts would ensue.  
*/

SUMO_COMPLEX B_VV (SUMO_REAL x, SUMO_REAL y, TSIL_REAL s, TSIL_REAL qq)
{
  SUMO_COMPLEX result;

  if (TSIL_FABS(s/(x+y+s+qq)) < TSIL_TOL) {
    if (TSIL_FABS((x-y)/(x+y+qq)) < TSIL_TOL) {
      return (3.0L * (TSIL_A(x,qq)/x + 1.0L));
    } else
    return (3.0L * (TSIL_A(x,qq) - TSIL_A(y,qq))/(x-y));
  }

  if ((TSIL_FABS(x/(x+y+s)) > 100.0 * TSIL_TOL) && 
      (TSIL_FABS(y/(x+y+s)) > 100.0 * TSIL_TOL)) 
  {
    result = (y * TSIL_A(x,qq) + x * TSIL_A(y,qq) 
             - s * s * TSIL_CREAL(TSIL_B(0,0,s,qq))
             + (s-y) * (s-y) * TSIL_CREAL(TSIL_B(0,y,s,qq))
             + (s-x) * (s-x) * TSIL_CREAL(TSIL_B(0,x,s,qq))
             + (-s*s + 2*s*x - x*x + 2*s*y - 10*x*y - y*y) * TSIL_B(x,y,s,qq)
             )/(4.0L * x * y);
  } 
  else if (TSIL_FABS(x/(x+y+s)) > 100.0 * TSIL_TOL) {
    result = 0.75L * (TSIL_A(x,qq)- s * TSIL_CREAL(TSIL_B(0,0,s,qq)) +
            (s - 3.0L * x) * TSIL_B(0,x,s,qq))/x;
  } 
  else if (TSIL_FABS(y/(x+y+s)) > 100.0 * TSIL_TOL) {
    result = 0.75L * (TSIL_A(y,qq)- s * TSIL_CREAL(TSIL_B(0,0,s,qq)) +
            (s - 3.0L * y) * TSIL_B(0,y,s,qq))/y;
  } 
  else result = 1.5L - 3.0L * TSIL_B(0,0,s,qq);

  return (result);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Derivative of previous function with respect to first mass^2
   argument, x. */

SUMO_COMPLEX B_VpV (SUMO_REAL x, SUMO_REAL y, TSIL_REAL s, TSIL_REAL qq)
{
  SUMO_COMPLEX result;

  if (TSIL_FABS(s/(x+y+s+qq)) < TSIL_TOL) {
    if (TSIL_FABS((x-y)/(x+y+qq)) < TSIL_TOL) {
      return (1.5L/x);
    } else
    return (3.0L * (y*(TSIL_Ap(y,qq) - TSIL_Ap(x,qq))/((x-y)*(x-y)) + 
                     1.0L/(x-y)) );
  }

  if ((TSIL_FABS(x/(x+y+s)) > 100.0 * TSIL_TOL) && 
      (TSIL_FABS(y/(x+y+s)) > 100.0 * TSIL_TOL)) 
  {
    result = (y * TSIL_Ap (x,qq) + TSIL_A(y,qq)
              -2.0 * (s-x) * TSIL_CREAL(TSIL_B(0,x,s,qq))
              + (s-x) * (s-x) * TSIL_CREAL(SUMO_Bp(x,0,s,qq,interp))
              + (2.0 * s - 2.0 * x - 10.0 * y) * TSIL_B(x,y,s,qq)
              + (-s*s + 2*s*x - x*x + 2*s*y - 10*x*y - y*y) 
	      * SUMO_Bp(x,y,s,qq,interp))/(4.0L * x * y)
              - (y * TSIL_A(x,qq) + x * TSIL_A(y,qq) 
             - s * s * TSIL_CREAL(TSIL_B(0,0,s,qq))
             + (s-y) * (s-y) * TSIL_CREAL(TSIL_B(0,y,s,qq))
             + (s-x) * (s-x) * TSIL_CREAL(TSIL_B(0,x,s,qq))
             + (-s*s + 2*s*x - x*x + 2*s*y - 10*x*y - y*y) * TSIL_B(x,y,s,qq)
             )/(4.0L * x * x * y);
  } 
  else if (TSIL_FABS(x/(x+y+s)) > 100.0 * TSIL_TOL) {
    result = 0.75L * (x + s * (TSIL_CREAL(TSIL_B(0,0,s,qq)) 
			       - TSIL_B(0,x,s,qq)) +
		      x * (s - 3.0 * x) * SUMO_Bp(x,0,s,qq,interp))/(x*x);
  } 
  else result = 1.0L/0.0L; /* Don't have x=0 case, probably don't need. */

  return result;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Functions needed for the two-loop self energy (notations as in
   hep-ph/0405022). */
/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/*  hep-ph/0312092 eqs. (5.32) and (5.33).
    For G_FF (x,y) the TSIL_RESULT supplied should have args
    (x,x,y,y,0). */

void G_FF (TSIL_RESULT *r, TSIL_COMPLEX *gFF, TSIL_COMPLEX *gff)
{
  SUMO_COMPLEX Mxxyy0;
  SUMO_COMPLEX Uxyy0, Uyxx0;
  SUMO_COMPLEX Vxyy0, Vyxx0;
  SUMO_COMPLEX S0xy, Tx0y, Ty0x, Tbarx0y;
  SUMO_COMPLEX Bxy, Bxpy, Bxyp, Ax, Ay;
  SUMO_COMPLEX Mplus, Mminus;
  TSIL_REAL x, y, s, qq;
  TSIL_REAL delta;
  char funcname[] = "G_FF";
  int check;

  check = SUMO_FPCompare (r->x, r->y) 
    * SUMO_FPCompare (r->z, r->u) 
    * SUMO_FPCompare (r->v, 0.0L);

  if (check == 0)
    TSIL_Error ("G_FF","TSIL data has wrong arguments.",23);

  x = r->x;
  y = r->z;
  s = r->s;
  qq = r->qq;

  /* Test for threshold case and if yes, do minimal interpolation: */
  delta = s/TSIL_POW(TSIL_SQRT(x)+TSIL_SQRT(y),2) - 1.0L;

  if (TSIL_FABS(delta) < THRESH_TOL) {
    TSIL_Manalytic (x,x,y,y,0,s*(1.0L - THRESH_TOL), &Mminus);
    TSIL_Manalytic (x,x,y,y,0,s*(1.0L + THRESH_TOL), &Mplus);
    Mxxyy0 = 0.5L*(1.0L + delta/THRESH_TOL)*Mplus + 
             0.5L*(1.0L - delta/THRESH_TOL)*Mminus;
  } else
    Mxxyy0 = r->M;

  Uxyy0 = r->U[xzuv];
  Uyxx0 = r->U[uyxv];
  Vxyy0 = SUMO_GetFunctionR (r, "Vxzuv", interp);
  Vyxx0 = SUMO_GetFunctionR (r, "Vuyxv", interp);

  S0xy = r->S[uxv];
  Tx0y = r->T[xuv];
  Ty0x = r->T[uxv];
  Tbarx0y = r->TBAR[xuv];

  Bxy = r->B[xz];
  Ax = TSIL_A (x, qq);
  Ay = TSIL_A (y, qq);
  Bxpy = SUMO_Bp (x, y, s, qq, interp);
  Bxyp = SUMO_Bp (y, x, s, qq, interp);

  if ((TSIL_FABS(x) > TSIL_TOL) && (TSIL_FABS(y) > TSIL_TOL)) {
    /* x and y are both nonzero: */

    *gFF = -2.0L * TSIL_POW(x + y - s, 2) * Mxxyy0
      + 4.0L * (x + y - s) * (Uxyy0 - y*Vxyy0 + Uyxx0 - x*Vyxx0)
      + 2.0L * S0xy
      + (3.0L*x + y - s) * Tx0y
      + (3.0L*y + x - s) * Ty0x
      + 2.0L * (x + y) * Bxy * Bxy
      + (y - x - s)*Ax*Bxy/x + (x - y - s)*Ay*Bxy/y
      + 2.0L*(s - x - y)*(Ax*Bxpy + Ay*Bxyp)
      - 5.0L*(TSIL_I2(0,x,x,qq) + TSIL_I2(0,y,y,qq))
      + 6.0L*(Ax + Ay) - 4.0L*(x + y) - s/2.0L;

    *gff = 4.0L*(s - x - y)*Mxxyy0 + 8.0L*(Uxyy0 + Uyxx0)
      - 8.0L*(y*Vxyy0 + x*Vyxx0) + 4.0L*Bxy*Bxy - 4.0L*(Ax*Bxpy + Ay*Bxyp);
  }
  else if ((TSIL_FABS(x) < TSIL_TOL) && (TSIL_FABS(y) > TSIL_TOL)) {
    /* x = 0 and y != 0 */

    *gFF = (-29.0L*s + 32.0L*y)/4.0L + 7.0L*Ay - Ay*Ay/y + 4.0L*(-3.0L*s + 4.0L*y)*Bxy + 
      ((s - 9.0L*y)*Ay*Bxy)/y + 2.0L*y*Bxy*Bxy - 
      2.0L*(-s + y)*(-s+y)*Mxxyy0 + (s - y)*Ty0x - 
      5.0L*(-s + y)*Tbarx0y;

    /* This works since the term multiplying gff will contain a zero
       (quark mass). */
    *gff = 0.0L;

    /* If needed someday... */
    /* if (TSIL_FABS (y-s) > 10.0*TSIL_TOL) { */
      /* y != s */

    /*   *gff = (8*(-2*s + y))/(-s + y) - 20*Ay/(s - y) - 12*Ay*Ay/(y*(-s + y)) + */
    /*   	4*(-7*s + 9*y)*Bxy/(-s + y) - 4*(-s + 4*y)*Ay*Bxy/(y*(-s + y)) + */
    /*   	4*Bxy*Bxy - 4*(-s + y)*Mxxyy0 - 4*Ty0x - 12*Tbarx0y; */
    /* } */
    /* else if (TSIL_FABS(y) > TSIL_TOL) { */
      /* y = s but y and s != 0 */

      /* *gff = -40*Ay/y + 8*Ay*Ay/(y*y) - 4*(Ty0x + 3*(-4 + Tbarx0y)); */
    /* } */
    /* else { */
      /* y = s = 0; not yet implemented */
    /*   TSIL_Error (funcname, "This can't happen; case not implemented.", 321); */
    /* } */

  }
  else if ((TSIL_FABS(x) < TSIL_TOL) && (TSIL_FABS(y) < TSIL_TOL)) {
    /* x and y are both 0 */

    *gFF = (-29.0L*s)/4.0L - 13.0L*s*Bxy - 2.0L*s*s*Mxxyy0 + 6.0L*s*Tbarx0y;

    *gff = 0.0L;

    /* If needed someday... */
    /* *gff = 4*(4 + 8*Bxy + Bxy*Bxy + s*Mxxyy0 - 4*Tbarx0y); */
  }
  else {
    /* x !=0 and y = 0; should never happen */
    TSIL_Error (funcname, "Can't happen: x !=0 and y = 0.", 314);
  }

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq. (5.31). 
   The TSIL_RESULT supplied should have args (x,x,y,y,0). */

void G_SS (TSIL_RESULT *r,TSIL_COMPLEX *gSS)
{
  SUMO_COMPLEX Mxxyy0;
  SUMO_COMPLEX Uxyy0, Uyxx0;
  SUMO_COMPLEX Vxyy0, Vyxx0;
  /* SUMO_COMPLEX S0xy, Tx0y, Ty0x; */
  SUMO_COMPLEX Bxy, Bxpy, Bxyp, Ax, Ay;
  TSIL_REAL x, y, s, qq;
  int check;

  check = SUMO_FPCompare (r->x, r->y)
    * SUMO_FPCompare (r->z, r->u) 
    * SUMO_FPCompare (r->v, 0.0L);

  if (check == 0)
    TSIL_Error ("G_SS", "TSIL data has wrong arguments.", 23);

  x = r->x;
  y = r->z;
  s = r->s;
  qq = r->qq;

  Mxxyy0 = r->M;
  Uxyy0 = r->U[xzuv];
  Uyxx0 = r->U[uyxv];
  Vxyy0 = SUMO_GetFunctionR (r, "Vxzuv", interp);
  Vyxx0 = SUMO_GetFunctionR (r, "Vuyxv", interp);
  /* Vxyy0 = r->V[xzuv]; */
  /* Vyxx0 = r->V[uyxv]; */

  Bxy = r->B[xz];
  Ax = TSIL_A (x, qq);
  Ay = TSIL_A (y, qq);
  Bxpy = SUMO_Bp (x, y, s, qq, interp);
  Bxyp = SUMO_Bp (y, x, s, qq, interp);

  *gSS = 2.0L*(x + y - s)*Mxxyy0 + 4.0L*y*Vxyy0 + 4.0L*x*Vyxx0
    - 4.0L*Uxyy0 - 4.0L*Uyxx0 - Ay*Bxyp - Ax*Bxpy + Bxy*Bxy;

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Special form for massless vector bosons (hep-ph/0312092
   eq. (5.1)). Needs to be expanded eventually to cover the general
   case. */

TSIL_COMPLEX W_SSSV (TSIL_REAL x, TSIL_REAL qq)
{
  return 3.0L*TSIL_I2 (0,x,x,qq) - TSIL_A(x,qq) + 2.0L*x;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs. 4.25 to 4.27 
   These are functions of x,y,z,u,v. The TSIL_RESULT supplied should
   have these same arguments. */

void M_SFSFF (TSIL_RESULT *r, 
	      TSIL_COMPLEX *mSfSff, 
	      TSIL_COMPLEX *mSFSfF,
	      TSIL_COMPLEX *mSFSFf)
{
  TSIL_REAL x, y, u, v, s;

  x = r->x;
  y = r->y;
  u = r->u;
  v = r->v;
  s = r->s;

  *mSfSff = 2.0L* (r->M);

  *mSFSfF = (v - x + y)*(r->M) + r->U[yuzv] - r->U[xzuv]
    - (r->B[xz]) * (r->B[yu]);

  *mSFSFf = (y + u - s)*(r->M) - r->U[xzuv] - r->U[zxyv];

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs A.11 through A.13
   This is Vbar (x,0,u,v)
   The TSIL_RESULT provided should have params x,*,0,u,v,s,qq */

TSIL_COMPLEX Vbar (TSIL_RESULT *r)
{
  TSIL_REAL x, u, v, s, qq;
  TSIL_COMPLEX Ux0uv, xTxuv, uTuxv, vTvxu, TBAR00x, Sxuv;
  TSIL_COMPLEX Ax, Au, Av, B0x, I0uu;
  TSIL_COMPLEX res;

  if (r->z > TSIL_TOL)
    TSIL_Error ("Vbar", "Wrong data arguments!", 13);

  x = r->x;
  u = r->u;
  v = r->v;
  s = r->s;
  qq = r->qq;
  
  if (TSIL_FABS(u-v) > TSIL_TOL) {
    /* u != v */

    Ux0uv = r->U[xzuv];
    if (x > TSIL_TOL) xTxuv = x*(r->T[xuv]); else xTxuv = 0.0L;
    if (u > TSIL_TOL) uTuxv = u*(r->T[uxv]); else uTuxv = 0.0L;
    if (v > TSIL_TOL) vTvxu = v*(r->T[vxu]); else vTvxu = 0.0L;
    Sxuv = r->S[uxv];
    Ax = TSIL_A(x,qq);
    Au = TSIL_A(u,qq);
    Av = TSIL_A(v,qq);
    B0x = r->B[xz];

    res = -x*Ux0uv/((s-x)*(s-x)) + ((u+v)/((u-v)*(u-v)*(s-x)) - 
       1.0L/((s-x)*(s-x)))*xTxuv + (v*(u-v-x+s)/((u-v)*(u-v)*(s-x)) - 
				    2.0L*x/((u-v)*(s-x)*(s-x)))*uTuxv
      + (u*(u-v+x-s)/(TSIL_POW(u-v,3)*(s-x)) +
	 2.0L*x/((u-v)*(s-x)*(s-x)))*vTvxu
      + ((u+v)*(Sxuv-Ax/2.0L+x+u+v-5.0L*s/8.0L)+Au*Av)/((u-v)*(u-v)*(s-x))
      + B0x*((2.0L*x+s)*(Av-Au)/((s-x)*(s-x)*(u-v)) +
	     (v*Au - u*Av)/TSIL_POW(u-v,3) -
	     (u+v)/(2.0L*(u-v)*(u-v))) + Ax*(Av-Au)/((u-v)*(s-x)*(s-x)) +
      ((x*u+x*v-2.0L*s*v)*Au + 
       (x*u+x*v-2.0L*s*u)*Av)/((u-v)*(u-v)*(s-x)*(s-x));
  }
  else if (u > TSIL_TOL) {
    /* u = v but both are nonzero */

    if (u > TSIL_TOL) uTuxv = u*(r->T[uxv]); else uTuxv = 0.0L;
    if (x > TSIL_TOL) xTxuv = x*(r->T[xuv]); else xTxuv = 0.0L;
    Sxuv = r->S[uxv];
    Ax = TSIL_A(x,qq);
    Au = TSIL_A(u,qq);
    I0uu = TSIL_I2 (0,u,u,qq);
    B0x = r->B[xz];

    res = 1.0L/(3.0L*u*(s-x)*(s-x))*
      ((4.0L*u-3.0L*x-s)*uTuxv + (s-x)*xTxuv + (4.0L*u-x+s)*Sxuv
       -2.0L*u*I0uu + 2.0L*Au*(x-2.0L*u-Ax) + Ax*(x/2.0L-4.0L*u-s/2.0L)
       +4.0L*u*u - u*x - x*x + 13.0L*s*x/8.0L + s*u/2.0L - 5.0L*s*s/8.0L)
      - B0x*((Au/u + 1.0L)*(s+x)/((s-x)*(s-x)) + 1.0L/(6.0L*u));

  }
  else {
    /* u = v = 0 */

    TBAR00x = r->TBAR[xuv];
    Sxuv = r->S[uxv];
    B0x = r->B[xz];
    Ax = TSIL_A(x,qq);

    res = (-(s+x)*TBAR00x + 2.0L*Sxuv - (s+x)*B0x - 3.0L*Ax + x - s/4.0L)/
      ((s-x)*(s-x));
  }

  return res;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs 4.33 through 4.43 */
/* Note that r1 is the tsil struct with the specified arguments for
   these functions; r2 is a tsil result with the same arguments but
   with y <-> z. */

/* Tested for the y = z != 0 case so far!! */

void V_FFFFS (TSIL_RESULT  *r1, TSIL_RESULT  *r2, 
	      TSIL_COMPLEX *vFFFFS, TSIL_COMPLEX *vFffFS,
	      TSIL_COMPLEX *vfFfFS, TSIL_COMPLEX *vfFFfS,
	      TSIL_COMPLEX *vFFffS, TSIL_COMPLEX *vffffS)
{
  TSIL_REAL x, y, z, u, v, s, qq;
  TSIL_COMPLEX Uxzuv, Uxyuv, Vxyuv;
  TSIL_COMPLEX xTxuv, uTuvx, vTvxu;
  TSIL_COMPLEX Sxuv;
  TSIL_COMPLEX Ax, Ay, Az, Av, Au, Bxy, Bxz;
  TSIL_COMPLEX Iyuv, Izuv;
  TSIL_COMPLEX Ipyuv, Bxyp;
  int check;

  /* Check arguments for sanity: x, u, v, s, qq should all agree. */
  check = SUMO_FPCompare(r1->x, r2->x)
    * SUMO_FPCompare(r1->u, r2->u)
    * SUMO_FPCompare(r1->v, r2->v)
    * SUMO_FPCompare(r1->s, r2->s)
    * SUMO_FPCompare(r1->qq, r2->qq);

  if (check == 0)
    TSIL_Error ("V_FFFFS", "Wrong data arguments", 12);

  /* Following is disabled so that trickeration may be used in eq 4.13
     of Pi2nonQCD. */

  /* Also check that y and z are exchanged between the two
     TSIL_RESULTs (separate from the above to that a separete warning
     message can be issued). */
/*   mustbe0 = TSIL_FABS(r1->z - r2->y) */
/*           + TSIL_FABS(r1->y - r2->z); */

/*   if (mustbe0 > TSIL_TOL) */
/*     TSIL_Error ("V_FFFFS", */
/* 		"Second TSIL_RESULT must be first with y <-> z", 13); */

  /* For convenience: */
  x = r1->x;
  y = r2->z;
  z = r1->z;
  u = r1->u;
  v = r1->v;
  s = r1->s;
  qq = r1->qq;

  Uxyuv = r2->U[xzuv];    
  if (x > TSIL_TOL) xTxuv = x*(r1->T[xuv]); else xTxuv = 0.0L;
  if (u > TSIL_TOL) uTuvx = u*(r1->T[uxv]); else uTuvx = 0.0L;
  if (v > TSIL_TOL) vTvxu = v*(r1->T[vxu]); else vTvxu = 0.0L;
  Sxuv = r1->S[uxv];
  Ax = TSIL_A (x,qq);
  Ay = TSIL_A (y,qq);
  Au = TSIL_A (u,qq);
  Av = TSIL_A (v,qq);
  Iyuv = TSIL_I2 (y,u,v,qq);
  Bxy = r2->B[xz];

  /* DGR ADDED 1/16/15: */
  Bxz = r1->B[xz];

  if (TSIL_FABS((y - z)/(x + y + z + u + v)) > TSIL_TOL) {
    /* Generic case with y != z */

    Uxzuv = r1->U[xzuv];
    Az = TSIL_A (z,qq);
    Izuv = TSIL_I2 (z,u,v,qq);

    *vffffS = -2.0L*(Uxyuv - Uxzuv)/(y - z);

    *vfFfFS = ((v - y - u)*Uxyuv + (Av - Au)*Bxy)/(y - z)
            + ((v - z - u)*Uxzuv + (Av - Au)*Bxz)/(z - y);

    *vfFFfS = 2.0L*(z*Uxzuv - y*Uxyuv)/(y - z);

    *vFffFS = ((y + u - v)*((s - x - y)*Uxyuv + Iyuv) +
	       (Au - Av)*((s - x - y)*Bxy + Ay))/(2.0L*y*(y - z)) +
              ((z + u - v)*((s - x - z)*Uxzuv + Izuv) +
               (Au - Av)*((s - x - z)*Bxz + Az))/(2.0L*z*(z - y)) +
      ((u - v)*(2.0L*xTxuv + 2.0L*Sxuv - Ax - Au - Av + x + u + v - s/4.0L)
       + (x + u - v - s)*uTuvx - (x + v - u - s)*vTvxu)/(2.0L*y*z);

    *vFFffS = ((s - x - y)*Uxyuv + Iyuv)/(y - z)
            + ((s - x - z)*Uxzuv + Izuv)/(z - y);

    *vFFFFS = ((y + u - v)*((s - x - y)*Uxyuv + Iyuv) + y*Sxuv +
	       (Au - Av)*((s - x - y)*Bxy + Ay))/(2.0L*(y - z)) +
              ((z + u - v)*((s - x - z)*Uxzuv + Izuv) + z*Sxuv +
               (Au - Av)*((s - x - z)*Bxz + Az))/(2.0L*(z - y));
  }
  else {
    /* y = z (y = z = 0 has some special features below) */

    Vxyuv = SUMO_GetFunctionR (r2, "Vxzuv", interp);
    /* Vxyuv = r2->V[xzuv]; */
    Ipyuv = TSIL_I2p (y,u,v,qq);
    Bxyp = SUMO_Bp (y,x,s,qq,interp);

    *vffffS = -2.0L*(-Vxyuv);
    *vfFfFS = (y + u - v)*Vxyuv - Uxyuv + (Av - Au)*Bxyp;
    *vfFFfS = -2.0L*Uxyuv;

    if (y > TSIL_TOL) 
      *vfFFfS += 2.0L*y*Vxyuv;

    *vFffFS = (y*(y + u - v)*((x + y - s)*Vxyuv + Ipyuv) +
	       ((u - v)*(x - s) - y*y)*Uxyuv
	       + (u - v)*(2.0L*xTxuv + 2.0L*Sxuv - Iyuv
			  - Ax - Au - Av + x + u + v - s/4.0L)
	       + (x + u - v - s)*uTuvx - (x + v - u - s)*vTvxu
	       + (Au - Av)*((s - x - y)*y*Bxyp + (x - s)*Bxy+y))/(2.0L*y*y);

    *vFFffS = (x + y - s)*Vxyuv - Uxyuv + Ipyuv;

    if (y > TSIL_TOL) {

      *vFFFFS = ((y+u-v)*((x+y-s)*Vxyuv + Ipyuv) +
		 (s-x-2.0L*y-u+v)*Uxyuv + Sxuv + Iyuv +
		 (Au-Av)*((s-x-y)*Bxyp - Bxy + TSIL_CLOG(y/qq)))/2.0L;
    }
    else {
      Bxyp = SUMO_Bp (x,0,s,qq,interp);
      *vFFFFS = ((u-v)*(x-s)*Vbar (r1) + (s-x-u+v)*Uxyuv + Sxuv
		 + Iyuv)/2.0L + (Au-Av)*(s*Bxyp + Ax/x);
    }
  }

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs 4.16 and 4.17 */
/* Requires no TSIL evaluation... */

TSIL_COMPLEX W_SSFF (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z,
		     TSIL_REAL u, TSIL_REAL qq)
{
  if (TSIL_FABS(y - x) < TSIL_TOL)
    return (z+u-x)*TSIL_I2p (x,z,u,qq) - TSIL_I2 (x,z,u,qq)
      - (TSIL_A(u,qq) + TSIL_A(z,qq))*TSIL_LOG(x/qq);
  else
    return ((z+u-x)*TSIL_I2 (x,z,u,qq) - (z+u-y)*TSIL_I2 (y,z,u,qq)
      - (TSIL_A(x,qq) - TSIL_A(y,qq)) * (TSIL_A(z,qq) + TSIL_A(u,qq)))/(x-y);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs 4.2 and 4.10 */
/* Requires no TSIL evaluation... */

TSIL_COMPLEX W_SSSS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, 
		     TSIL_REAL u, TSIL_REAL qq)
{
  if (TSIL_FABS(y - x) < TSIL_TOL)
    return -TSIL_I2p (x,z,u,qq);
  else
    return (TSIL_I2(x,z,u,qq) - TSIL_I2(y,z,u,qq))/(y - x);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq 4.15 */
/* Requires no TSIL evaluation... */

TSIL_COMPLEX W_SSff (TSIL_REAL x, TSIL_REAL y,TSIL_REAL z,
		     TSIL_REAL u, TSIL_REAL qq)
{
  return -2.0L * W_SSSS (x,y,z,u,qq);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs 4.8 and 4.13 */

/* For function args (x,y,z,u,v) r1 should have the same args, r2
   should have y <-> z. */

void V_SSSSS (TSIL_RESULT *r1, 
	      TSIL_RESULT *r2, 
	      TSIL_COMPLEX *vSSSSS)
{
  TSIL_REAL y, z;
  int check;

  y = r2->z;
  z = r1->z;

  /* Check arguments for sanity: x, u, v, s, qq should all agree. */
  check = SUMO_FPCompare(r1->x, r2->x)
    * SUMO_FPCompare(r1->u, r2->u)
    * SUMO_FPCompare(r1->v, r2->v)
    * SUMO_FPCompare(r1->s, r2->s)
    * SUMO_FPCompare(r1->qq, r2->qq);

  if (check == 0)
    TSIL_Error ("V_SSSSS", "Wrong data arguments", 12);

  if (TSIL_FABS(y - z) < TSIL_TOL) 
    *vSSSSS = -SUMO_GetFunctionR (r1, "Vxzuv", interp);
  else
    *vSSSSS = (r2->U[xzuv] - (r1->U[xzuv]))/(y - z);

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs 4.29 through 4.31 */

/* For function args (x,y,z,u,v) r1 should have the same args, r2
   should have y <-> z. */

void V_SSSFF (TSIL_RESULT *r1,TSIL_RESULT *r2,
	      TSIL_COMPLEX *vSSSFF,TSIL_COMPLEX *vSSSff)
{
  TSIL_REAL x,y,z,u,v,qq;
  TSIL_COMPLEX s;
  /* TSIL_COMPLEX vSSSSS; */
  TSIL_COMPLEX Uxyuv, Uxzuv;
  TSIL_COMPLEX Vxyuv;
  TSIL_COMPLEX Bxy, Bxz, Bxyp, Au, Av;
  int check;

  /* Check arguments for sanity: x, u, v, s, qq should all agree. */
  check = SUMO_FPCompare(r1->x, r2->x)
    * SUMO_FPCompare(r1->u, r2->u)
    * SUMO_FPCompare(r1->v, r2->v)
    * SUMO_FPCompare(r1->s, r2->s)
    * SUMO_FPCompare(r1->qq, r2->qq);

  if (check == 0)
    TSIL_Error ("V_SSSSS", "Wrong data arguments", 12);

  x = r1->x;
  y = r2->z;
  z = r1->z;
  u = r1->u;
  v = r1->v;
  s = r1->s;
  qq = r1->qq;

  Au = TSIL_A (u,qq);
  Av = TSIL_A (v,qq);
  Bxy = r2->B[xz];
  Uxyuv = r2->U[xzuv];

  if (TSIL_FABS (y-z) > TSIL_TOL) {
    Uxzuv = r1->U[xzuv];
    Bxz = r1->B[xz];
    *vSSSFF = ((y-u-v)*Uxyuv + (Au+Av)*Bxy
              -(z-u-v)*Uxzuv - (Au+Av)*Bxz)/(y-z);
    *vSSSff = -2.0 * (Uxyuv - Uxzuv)/(y-z);
  }
  else {

    Vxyuv = SUMO_GetFunctionR (r2, "Vxzuv", interp);
    /* Vxyuv = r2->V[xzuv]; */
    Bxyp = SUMO_Bp(z,x,s,qq,interp);

    *vSSSFF = (u+v-z) * Vxyuv + Uxyuv + (Au+Av) * Bxyp;
    *vSSSff = 2.0L * Vxyuv;
  }

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs 4.3 and 4.11 */

TSIL_COMPLEX X_SSS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq)
{
  if (TSIL_FABS((x - y)/(x + y + z)) < TSIL_TOL)
    return TSIL_A(z, qq) * TSIL_LOG(x/qq);
  else 
    return TSIL_A(z, qq) * (TSIL_A(x, qq) - TSIL_A(y, qq))/(x - y);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs 4.4 and 4.12 */

TSIL_COMPLEX Y_SSSS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z,
		     TSIL_REAL u, TSIL_REAL s, TSIL_REAL qq)
{
  if (TSIL_FABS(z - y) < TSIL_TOL)
    return -TSIL_A(u, qq)*SUMO_Bp(y,x,s,qq,interp);
  else
    return TSIL_A(u, qq)*(TSIL_B(x,z,s,qq)-TSIL_B(x,y,s,qq))/(y - z);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq 4.5 */

TSIL_COMPLEX Z_SSSS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z,
		     TSIL_REAL u,TSIL_REAL s, TSIL_REAL qq)
{
  return TSIL_B(x,y,s,qq)*TSIL_B(z,u,s,qq);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq 4.6 */
/* For input data args x,y,z,u,v returns -S(v,y,z). */

void S_SSS (TSIL_RESULT *r, TSIL_COMPLEX *sSSS)
{
  *sSSS = -(r->S[vyz]);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq 4.9 */

void M_SSSSS (TSIL_RESULT *r, TSIL_COMPLEX *mSSSSS)
{
  *mSSSSS = -(r->M);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq 4.7 */
/* For input data args x,y,z,u,v returns U(x,z,u,v) */

void U_SSSS (TSIL_RESULT *r, TSIL_COMPLEX *uSSSS)
{
  *uSSSS = r->U[xzuv];
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs 4.19 through 4.23 */

void M_FFFFS (TSIL_RESULT *r,
	      TSIL_COMPLEX *mFFFFS, TSIL_COMPLEX *mFFffS,
	      TSIL_COMPLEX *mFfFfS, TSIL_COMPLEX *mFffFS,
	      TSIL_COMPLEX *mffffS)
{
  TSIL_REAL x,y,z,u,v,qq;
  TSIL_REAL s;
  TSIL_COMPLEX M, Uzxyv, Uuyxv, Uxzuv, Uyuzv;
  TSIL_COMPLEX BxzByu, Sxuv, Syzv;

  x = r->x;
  y = r->y;
  z = r->z;
  u = r->u;
  v = r->v;
  s = r->s;
  qq = r->qq;

  M      = r->M;
  Uzxyv  = r->U[zxyv];
  Uuyxv  = r->U[uyxv];
  Uxzuv  = r->U[xzuv];
  Uyuzv  = r->U[yuzv];
  BxzByu = (r->B[xz])*(r->B[yu]);
  Sxuv   = r->S[uxv];
  Syzv   = r->S[vyz];

  *mffffS = 2.0L*M;

  *mFffFS = (y+z-v-s)*M - Uxzuv - Uuyxv + BxzByu;

  *mFfFfS = (x+z-s)*M - Uyuzv - Uuyxv;

  *mFFffS = (x+y-v)*M - Uxzuv - Uyuzv + BxzByu;

  *mFFFFS = (x*u + y*z - v*s)*M - x*Uzxyv - z*Uxzuv - u*Uyuzv - y*Uuyxv
    + Sxuv + Syzv + s*BxzByu;

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eqs 4.19 through 4.23 */

/* This is a bit of a hack to accomodate hep-ph/0405022 eq 4.14. It
   computes four functions normally, then does a PermuteTSIL and
   computes two more.  Functions mFffFS and mFfFfSp have the permuted
   args. */

void M_FFFFS2 (TSIL_RESULT *r,
	       TSIL_COMPLEX *mFFFFS, TSIL_COMPLEX *mFFffS,
	       TSIL_COMPLEX *mFfFfS, TSIL_COMPLEX *mFffFS,
	       TSIL_COMPLEX *mffffS, TSIL_COMPLEX *mFfFfSp)
{
  TSIL_REAL x,y,z,u,v,qq;
  TSIL_REAL s;
  TSIL_COMPLEX M, Uzxyv, Uuyxv, Uxzuv, Uyuzv;
  TSIL_COMPLEX BxzByu, Sxuv, Syzv;
  TSIL_RESULT gaak;

  x = r->x;
  y = r->y;
  z = r->z;
  u = r->u;
  v = r->v;
  s = r->s;
  qq = r->qq;

  M = r->M;
  Uzxyv = r->U[zxyv];
  Uuyxv = r->U[uyxv];
  Uxzuv = r->U[xzuv];
  Uyuzv = r->U[yuzv];
  BxzByu = (r->B[xz])*(r->B[yu]);
  Sxuv = r->S[uxv];
  Syzv = r->S[vyz];

  *mFFFFS = (x*u+y*z-v*s)*M - x*Uzxyv - z*Uxzuv - u*Uyuzv - y*Uuyxv
    + Sxuv + Syzv + s*BxzByu;
  *mFFffS = (x+y-v)*M - Uxzuv - Uyuzv + BxzByu;
  *mFfFfS = (x+z-s)*M - Uyuzv - Uuyxv;

  /* This one computed with swapped args here! */
/*   *mFffFS = (y+z-v-s)*M - Uxzuv - Uuyxv + BxzByu; */

  *mffffS = 2.0L*M;

  /* Now swap args and evaluate the remaining two functions */  
  TSIL_PermuteResult (r, XYandZU, &gaak);

  x = gaak.x;
  y = gaak.y;
  z = gaak.z;
  v = gaak.v;

  Uuyxv = gaak.U[uyxv];
  Uxzuv = gaak.U[xzuv];
  Uyuzv = gaak.U[yuzv];
  BxzByu = (gaak.B[xz])*(gaak.B[yu]);

  *mFfFfSp = (x+z-s)*M - Uyuzv - Uuyxv;
  *mFffFS = (y+z-v-s)*M - Uxzuv - Uuyxv + BxzByu;

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq. (5.7). Supplied TSIL_RESULT should have args
   (x,y,y,y,0). */

void V_SSSSV (TSIL_RESULT *foo, TSIL_COMPLEX *vSSSSV)
{
  TSIL_COMPLEX vxyy0,uxyy0,tbar0xy;
  TSIL_REAL x,y,s,qq;

  x = foo->x;
  y = foo->y;
  s = foo->s;
  qq = foo->qq;

  vxyy0 = SUMO_GetFunctionR (foo, "Vxzuv", interp);
  /* vxyy0 = foo->V[xzuv]; */
  uxyy0 = foo->U[xzuv];
  tbar0xy = foo->TBAR[vxu];

  *vSSSSV = 4.0L*y*vxyy0 - 2.0L*uxyy0
    - TSIL_A(y,qq)*SUMO_Bp(y,x,s,qq,interp) + tbar0xy;

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

TSIL_COMPLEX V_VSSSS (TSIL_REAL y, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v,
		      TSIL_REAL s, TSIL_REAL qq)
{
  TSIL_DATA foo;
  TSIL_COMPLEX u0yuv,u0zuv,v0yuv;

  if (TSIL_FABS(z-y) > TSIL_TOL) {
    TSIL_SetParameters (&foo, 0.0L, y, z, u, v, qq);
    TSIL_Evaluate (&foo, s);
    u0zuv = foo.U[xzuv].value;
    TSIL_SetParameters (&foo, 0.0L, z, y, u, v, qq);
    TSIL_Evaluate (&foo, s);
    u0yuv = foo.U[xzuv].value;

    return -3.0L*((y + s)*u0yuv + TSIL_I2(y,u,v,qq))/(y - z)
      - 3.0L*((z + s)*u0zuv + TSIL_I2(z,u,v,qq))/(z - y);
  }
  else {
    TSIL_SetParameters (&foo, 0.0L, y, y, u, v, qq);
    TSIL_Evaluate (&foo, s);
    u0yuv = foo.U[xzuv].value;
    v0yuv = SUMO_GetFunction (&foo, "Vxzuv", interp);
    /* v0yuv = foo.V[xzuv].value; */

    return 3.0L*((y + s)*v0yuv - u0yuv) - 3.0L*TSIL_I2p(y,u,v,qq);
  }
}


/* V_VSSSS (0,y,z,u,v)

   hep-ph/0312092 eq. (5.8). First supplied TSIL_RESULT should have
   args (0,y,z,u,v), second should have (0,z,y,u,v).  Assumed \xi = 1
   (Feynman gauge). */

/* void V_VSSSS (TSIL_RESULT *foo, TSIL_RESULT *bar, */
/* 	      TSIL_COMPLEX *vVSSSS) */
/* { */
/*   TSIL_COMPLEX u0yuv,u0zuv,v0yuv; */
/*   TSIL_REAL y,z,u,v,s,qq; */

/*   y = foo->y; */
/*   z = foo->z; */
/*   u = foo->u; */
/*   v = foo->v; */
/*   s = foo->s; */
/*   qq = foo->qq; */

/*   if (TSIL_FABS(y-z) < TSIL_TOL) { */
/*     u0yuv = bar->U[xzuv]; */
/*     v0yuv = bar->V[xzuv]; */
/*     *vVSSSS = 2.0L*((y + s)*v0yuv - u0yuv) - TSIL_I2p(y,u,v,qq); */
/*   } */
/*   else { */
/*     u0yuv = bar->U[xzuv]; */
/*     u0zuv = foo->U[xzuv]; */
/*     *vVSSSS = (2.0L*(y + s)*u0yuv - TSIL_I2(y,u,v,qq))/(y - z) */
/*       + (2.0L*(z + s)*u0zuv - TSIL_I2(z,u,v,qq))/(z - y) */
/*       ; */
/*   } */
/* } */

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq. (5.12), in Landau gauge (xi = 0). 
   Supplied TSIL_RESULT should have args (0,y,z,u,y). */

void M_VSSSS (TSIL_RESULT *foo, TSIL_COMPLEX *mVSSSS)
{
  TSIL_COMPLEX m0yzuy,uyuyz,uuy0y,u0zyu,uz0yy,ty0u,tu0y,s0yu,byu,b0z;
  TSIL_REAL y,z,u,s,qq;

  y = foo->y;
  z = foo->z;
  u = foo->u;
  s = foo->s;
  qq = foo->qq;

  m0yzuy = foo->M;
  uyuyz = foo->U[yuzv];
  uuy0y = foo->U[uyxv];
  u0zyu = foo->U[xzuv];
  uz0yy = foo->U[zxyv];
  ty0u = foo->T[vxu];
  tu0y = foo->T[uxv];
  s0yu = foo->S[uxv];
  byu = TSIL_B(y,u,s,qq);
  b0z = TSIL_B(0.0L,z,s,qq);

  *mVSSSS = (2.0L*y + z - 2.0L*u + s)*m0yzuy
    + uyuyz - uuy0y - u0zyu + 2.0L*uz0yy - b0z*byu
    + ((z+s)*(u0zyu - b0z*byu)
       - 2.0L*y*ty0u - 2.0L*u*tu0y - 4.0L*s0yu
       + 2.0L*TSIL_I2(y,z,u,qq) + 2.0L*(s - TSIL_A(z,qq))*byu
       + 2.0L*TSIL_A(u,qq) + 2.0L*TSIL_A(y,qq)
       - 2.0L*u - 2.0L*y + 0.5L*s)/(z - s)
    ;

}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0312092 eq. (5.34). Here we depart from the general
   convention -- for this function one supplies arguments
   (x,y,z,u,s,qq) and receives a return value despite that TSIL
   evaluation is required. All TSIL calls are carried out inside the
   routine. */

TSIL_COMPLEX G_SSSS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u,
		     TSIL_REAL s, TSIL_REAL qq)
{
  TSIL_COMPLEX result;
  TSIL_COMPLEX vSSSSV, vVSSSS, mVSSSS1, mVSSSS2;
  TSIL_DATA foo;
  TSIL_RESULT r1;

  /* V_SSSSV term: */
  TSIL_SetParameters (&foo, z, u, u, u, 0.0L, qq);
  TSIL_Evaluate (&foo, s);
  TSIL_CopyResult (&foo, &r1);
  V_SSSSV (&r1, &vSSSSV);

  /* V_VSSSS term: */
  /* TSIL_SetParameters (&foo, 0, x, y, z, u, qq); */
  /* TSIL_Evaluate (&foo, s); */
  /* TSIL_CopyResult (&foo, &r1); */
  /* TSIL_SetParameters (&foo, 0, y, x, z, u, qq); */
  /* TSIL_Evaluate (&foo, s); */
  /* TSIL_CopyResult (&foo, &r2); */
  /* V_VSSSS (&r1, &r2, &vVSSSS); */

  vVSSSS = V_VSSSS (x, y, z, u, s, qq);

  /* M_VSSSS term 1: */
  TSIL_SetParameters (&foo, 0.0L, u, x, z, u, qq);
  TSIL_Evaluate (&foo, s);
  TSIL_CopyResult (&foo, &r1);
  M_VSSSS (&r1, &mVSSSS1);

  /* M_VSSSS term 2: */
  TSIL_SetParameters (&foo, 0.0L, u, y, z, u, qq);
  TSIL_Evaluate (&foo, s);
  TSIL_CopyResult (&foo, &r1);
  M_VSSSS (&r1, &mVSSSS2);

  result = vSSSSV + vVSSSS + mVSSSS1 + mVSSSS2;

  return result;
}

/* ---------------------------------------------------------------- */
/* Following are defined in hep-ph/0502168, needed for squarks.     */
/* ---------------------------------------------------------------- */
/* hep-ph/0502168 eq. (3.40). Note repeated argument is suppressed. */

TSIL_COMPLEX Gtilde_SSS (TSIL_REAL x, TSIL_REAL z, TSIL_REAL qq)
{
  return (4.0L - 3.0L*TSIL_LOG(x/qq))*TSIL_A (z, qq);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/0502168 eq. (3.41). For Gtilde_SSSS (x,x,y,z) the supplied
   TSIL_RESULT object should have args (0,z,x,y,z). */

void Gtilde_SSSS (TSIL_RESULT *foo,
		  TSIL_COMPLEX *result)
{
  TSIL_COMPLEX m0zxyz,ux0zz,uzyzx,tbar0yz,s0yz,ty0z,tz0y;
  TSIL_REAL x,y,z,s,qq,lnbarx,lnbary,lnbarz;
  
  x = foo->z;
  y = foo->u;
  z = foo->y;
  s = foo->s;
  qq = foo->qq;

  m0zxyz  = foo->M;
  ux0zz   = foo->U[zxyv];
  uzyzx   = foo->U[yuzv];
  tbar0yz = foo->TBAR[xuv];
  s0yz    = foo->S[uxv];
  ty0z    = foo->T[uxv];
  tz0y    = foo->T[vxu];

  lnbarx = TSIL_LOG (x/qq);
  lnbary = TSIL_LOG (y/qq);
  lnbarz = TSIL_LOG (z/qq);

  *result = 4.0L*(x - y + z)*m0zxyz
    + 4.0L*ux0zz + 2.0L*uzyzx + 4.0L*tbar0yz - 4.0L*s0yz/x
    - 2.0L*y*ty0z/x + (2.0L - 2.0L*z/x)*tz0y
    + 2.0L*TSIL_I2(x,y,z,qq)/x + 3.0L*TSIL_I2p(x,y,z,qq)
    + z*(3.0L*lnbarz - 7.0L)*SUMO_Bp(z,y,s,qq,interp)
    + 2.0L*(lnbarx + lnbarz - 8.0L)*TSIL_B(y,z,s,qq)
    + (2.0L*y*lnbary + 2.0L*z*lnbarz - 7.L*x/2.L - 4.L*y - 4.L*z)/x
    ;

  /* printf("Ready to return from Gtilde_SSSS: g~ssss = ");TSIL_cprintf(*result);printf("\n"); */
  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/502168 eq. (3.33) (in the DRbar' scheme, of course). */

TSIL_REAL G_S (TSIL_REAL y, TSIL_REAL qq)
{
  TSIL_REAL lnbary = TSIL_LOG(y/qq);

  return y*(-12.0L + 11.0L*lnbary - 3.0L*lnbary*lnbary);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/502168 eqs. (3.42) and (3.43).  
   For Gtilde (x,x,y,z), the TSIL_RESULT supplied should have args
   (0,z,x,y,z) */

void Gtilde_SSFF (TSIL_RESULT *foo, TSIL_COMPLEX *gSSFF, TSIL_COMPLEX *gSSff)
{
  TSIL_REAL x,y,z,s,qq,lnbarx,lnbary,lnbarz;
  TSIL_COMPLEX m0zxyz,ux0zz,uzyzx,tbar0yz,tz0y,ty0z,s0yz,sxzz,gtildeSSSS;
  TSIL_COMPLEX Ax,Ay,Az,b0z,b0y;
  /* char funcname[] = "Gtilde_SSFF"; */

  x = foo->z;
  y = foo->u;
  z = foo->y;
  s = foo->s;
  qq = foo->qq;

  m0zxyz  = foo->M;
  ux0zz   = foo->U[zxyv];
  uzyzx   = foo->U[yuzv];
  tbar0yz = foo->TBAR[xuv];
  tz0y    = foo->T[vxu];
  ty0z    = foo->T[uxv];
  s0yz    = foo->S[uxv];
  sxzz    = foo->S[vyz];

  if (z < TSIL_TOL) {

    Ax = TSIL_A (x,qq);
    Ay = TSIL_A (y,qq);
    b0y = foo->B[yu];

    *gSSFF =
      (150*x + 8*PI2*x - 63*y - 16*PI2*y)/12. - (2*(5*x - 4*y)*Ax)/x +
      ((x - 2*y)*Ax*Ax)/(x*x) - ((3*x + y)*Ay)/x + (6*Ax*Ay)/x -
      (9*x*x - 14*x*y + y*y)*b0y/x - 2*(x + y)*Ax*b0y/x +
      (2*(2*x - y)*TSIL_I2(0,x,y,qq))/x + 4*(x - y)*(x - y)*m0zxyz -
      (x - 2*y)*s0yz/x + 5*(x - y)*tbar0yz - 4*y*uzyzx
      ;

  }
  else if (y > TSIL_TOL) {

    lnbarx = TSIL_LOG (x/qq);
    lnbary = TSIL_LOG (y/qq);
    lnbarz = TSIL_LOG (z/qq);

    /* The TSIL_RESULT foo has the correct args already for this call */
    Gtilde_SSSS (foo, &gtildeSSSS);

    *gSSff = -2.0L*gtildeSSSS + 6.0L*z*SUMO_Bp(z,y,s,qq,interp)*(1.0L - lnbarz);

    *gSSFF = 4.0L*((x - y)*(x - y) - z*z)*m0zxyz
      + 4.0L*(x - y - z)*ux0zz 
      - 4.0L*y*uzyzx
      + 4.0L*(x - y - z)*tbar0yz 
      + (x - 2.0L*z)*(x - y - z)*tz0y/x
      + 2.0L*y*(y + z)*ty0z/x 
      + (4.0L*(y + z)/x - 1.0L)*s0yz
      + 2.0L*sxzz + 3.0L*(x - y - z)*TSIL_I2p(x,y,z,qq)
      + (1.0L - 2.0L*(y + z)/x)*TSIL_I2(x,y,z,qq) 
      + 2.0L*(x - y - z)*z*(3.0L*lnbarz - 5.0L)*SUMO_Bp(z,y,s,qq,interp)
      + (15.0L*y + 17.0L*z - 7.0L*x - 2.0L*(x + y + z)*lnbarx 
	 + (x - y - 3.0L*z)*lnbarz)*TSIL_B(y,z,s,qq)
      + 3.0L*lnbarx*(y*lnbary - y + z*lnbarz - z)
      - y*(3.0L + 2.0L*(y + z)/x)*lnbary
      + z*(5.0L*lnbarz - 15.0L - 2.0L*(y + z)/x)*lnbarz
      + 4.0L*(y + z)*(y + z)/x - 15.0L*x/4.0L + 11.0L*y/2.0L + 33.0L*z/2.0L
      ;
  }
  else {     /* y = 0 and z != 0 */

    Ax = TSIL_A (x,qq);
    Az = TSIL_A (z,qq);
    b0z = foo->B[yu];

    /* This works because terms multiplying this will include a zero
       (the quark mass). */
    *gSSff = 0.0L;

    *gSSFF = (5*x*x - 29*x*z + 78*z*z)/(8.*z) - 3*Ax -
      (-x*x - 7*x*z + 2*z*z)*Az/(2.*x*z) + 6*Ax*Az/x -
      Az*Az/z + (x - 2*z)*(x*x - 18*x*z + z*z)*b0z/(2.*x*z) -
      2*(x + z)*Ax*b0z/x + (x - 9*z)*Az*b0z/z + 2*(2*x - z)*TSIL_I2(0,x,z,qq)/x +
      4*(x - z)*(x + z)*m0zxyz -
      (x*x - 2*x*z - 2*z*z)*s0yz/(x*z) + 2*sxzz +
      4*(x - z)*tbar0yz + 4*(x - z)*ux0zz
      ;
  }
  
  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/502168 eq. (3.20). Only valid when s = x! */

TSIL_COMPLEX F1tilde (TSIL_REAL x, TSIL_REAL qq)
{
  TSIL_REAL lnbarx = TSIL_LOG(x/qq);

  return x*(12.0L*PI2 - 16.0*PI2*ln2 - 11.0L/8.0L + 24.0L*Zeta3
	    - 39.0L*lnbarx/2.0L + 15.0L*lnbarx*lnbarx/2.0L);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/502168 eq. (3.21). Only valid when s = x! */

TSIL_COMPLEX F2tilde (TSIL_REAL x, TSIL_REAL qq)
{
  TSIL_REAL lnbarx = TSIL_LOG(x/qq);

  return x*(1147.0L/16.0L - 10.0L*PI2/3.0L + 8.0*PI2*ln2 - 12.0L*Zeta3
	    - 409.0L*lnbarx/12.0L + 19.0L*lnbarx*lnbarx/4.0L);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/502168 eq. (3.24), for convenience */

TSIL_COMPLEX SPMf (TSIL_REAL z)
{
  return z*(TSIL_Dilog((1.0L-z)/(1.0L+z)) - TSIL_Dilog((z-1.0L)/(1.0L+z))
	    + PI2/4.0L);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/502168 eqs.(3.22) and (3.25) and (3.27) and (3.29).
   Only valid when s = x! */

TSIL_COMPLEX F3tilde (TSIL_REAL x, TSIL_REAL y, TSIL_REAL qq)
{
  TSIL_REAL lnbarx = TSIL_LOG(x/qq);
  TSIL_REAL lnbary = TSIL_LOG(y/qq);

  if (y < TSIL_TOL && x > TSIL_TOL)
    return x*(19.0L*lnbarx/3.0L - lnbarx*lnbarx - 49.0L/4.0L - 2.0L*PI2/3.0L);
  else if (x < TSIL_TOL && y > TSIL_TOL)
    return y*(4.0L - 12.0L*lnbary);
  else if (TSIL_FABS(x-y) < TSIL_TOL)
    return x*(23.0L/4.0L - 17.0L*lnbarx/3.0L - lnbarx*lnbarx);
  else
    return 2.0L*(x - y*y/x)*TSIL_Dilog(1.0L - x/y) 
      + 8.0L*(x - y)*SPMf(TSIL_SQRT(y/x)) 
      + y*(18.0L + PI2*y/(3.0L*x) - 6.0L*(lnbarx + lnbary))
      + x*(-49.0L/4.0L - PI2/3.0L + 19.0L*lnbarx/3.0L
	   - 2.0L*lnbarx*lnbary + lnbary*lnbary);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* hep-ph/502168 eqs.(3.23) and (3.26) and (3.28) and (3.30).
   Only valid when s = x! */

TSIL_COMPLEX F4tilde (TSIL_REAL x, TSIL_REAL y, TSIL_REAL qq)
{
  TSIL_REAL lnbarx = TSIL_LOG(x/qq);
  TSIL_REAL lnbary = TSIL_LOG(y/qq);

  if (y < TSIL_TOL && x > TSIL_TOL)
    return x*(25.0L*lnbarx/6.0L - lnbarx*lnbarx/2.0L - 75.0L/8.0L - PI2/3.0L);
  else if (x < TSIL_TOL && y > TSIL_TOL)
    return y*(11.0L + 3.0L*lnbary*lnbary);
  else if (TSIL_FABS(x-y) < TSIL_TOL)
    return x*(8.0L*PI2/3.0L - 107.0L/8.0L + 25.0L*lnbarx/6.0L 
	      + 5.0L*lnbarx*lnbarx/2.0L);
  else
    return (x + 6.0L*y + y*y/x)*TSIL_Dilog(1.0L - x/y) 
      + 8.0L*(x + y)*SPMf(TSIL_SQRT(y/x)) 
      + y*(-4.0L - PI2*(1.0L + y/(6.0L*x)) 
	   + 7.0L*TSIL_LOG(x/y) + 3.0L*lnbary*lnbary)
      + x*(-75.0L/8.0L - PI2/6.0L + 25.0L*lnbarx/6.0L
	   - lnbarx*lnbary + lnbary*lnbary/2.0L);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

TSIL_COMPLEX BpFF (TSIL_REAL x, TSIL_REAL y, TSIL_REAL s, TSIL_REAL qq)
{
  return (x + y - s)*SUMO_dBds(x,y,s,qq,interp) - TSIL_B(x,y,s,qq);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

TSIL_COMPLEX Bpff (TSIL_REAL x, TSIL_REAL y, TSIL_REAL s, TSIL_REAL qq)
{
  return 2.0L*SUMO_dBds(x,y,s,qq,interp);
}
