/*
  Implements one-loop and two-loop vacuum energy functions given in
  hep-ph/0111209.
*/

#include "supermodel.h"

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (3.13) */

TSIL_COMPLEX SUMO_OneLoopVacuum (TSIL_REAL x, TSIL_REAL QQ)
{
  return (0.25L * x * (TSIL_A (x,QQ) - 0.5L * x));
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (4.12) */

TSIL_COMPLEX SUMO_Fvac_SSS (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z, 
			    TSIL_REAL QQ)
{
  return -TSIL_I2 (x, y, z, QQ);
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (4.13) */

TSIL_COMPLEX SUMO_Fvac_SS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL QQ)
{
  return (TSIL_A (x, QQ) * TSIL_A (y, QQ));
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (4.14) */

TSIL_COMPLEX SUMO_Fvac_FFS (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z, 
			    TSIL_REAL QQ)
{
  TSIL_COMPLEX Ax, Ay, Az;

  Ax = TSIL_A (x, QQ);
  Ay = TSIL_A (y, QQ);
  Az = TSIL_A (z, QQ);

  return (Ax * (Ay - Az) - Ay * Az + (x + y - z) * TSIL_I2 (x, y, z, QQ));
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (4.15) NOTE: MASS INSERTIONS ARE INCLUDED!!!! */

TSIL_COMPLEX SUMO_Fvac_ffS (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z, 
			    TSIL_REAL QQ)
{
  return (2.0L * TSIL_CSQRT(x * y) * TSIL_I2 (x, y, z, QQ));
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (4.16), (4.24) */

TSIL_COMPLEX SUMO_Fvac_SSV (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z, 
			    TSIL_REAL QQ)
{
  TSIL_COMPLEX Ax, Ay, Az;

  Ax = TSIL_A (x, QQ);
  Ay = TSIL_A (y, QQ);

  if (TSIL_FABS (z) < TSIL_TOL) 
  {
    return (3.L * (Ax * Ay + (x + y) * TSIL_I2 (x, y, 0.L, QQ))
      - 2.L * (x * Ax + y * Ay) + (x + y)*(x + y));
  }

  Az = TSIL_A (z, QQ);

  return (( 
    (2.L * (x*y + x*z + y*z) - x*x - y*y - z*z) * TSIL_I2 (x, y, z, QQ)
    + (x-y)*(x-y) * TSIL_I2 (x, y, 0.L, QQ)
    + (y - x - z) * Ax * Az + (x - y - z) * Ay * Az 
    + z * Ax * Ay + 2.L * z * (x + y - (z/3.L)) * Az
    )/z);
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (6.21) */

TSIL_COMPLEX SUMO_Fvac_VS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL QQ)
{
  return (3.L * TSIL_A (x, QQ) * TSIL_A (y, QQ));
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (6.22), (6.26) and (6.27) */

TSIL_COMPLEX SUMO_Fvac_VVS (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z, 
			    TSIL_REAL QQ)
{
  TSIL_REAL tmp;
  TSIL_COMPLEX Ax, Ay, Az;

  if (TSIL_FABS (x) < TSIL_FABS (y)) {tmp = y; y = x; x = tmp;}

  Az = TSIL_A (z, QQ);

  if (TSIL_FABS (x) < TSIL_TOL) 
    {
      return 1.5L * Az - 0.25L * z - 3.L * TSIL_I2 (0.L, 0.L, z, QQ);
    }
  
  Ax = TSIL_A (x, QQ);

  if (TSIL_FABS (y) < TSIL_TOL) 
    {
      return ((x + 2.L * z + 3.L * Ax * Az/x
	       + (3.L * z/x - 9.L) * TSIL_I2 (0.L, x, z, QQ)
	       - 3.L * (z/x) * TSIL_I2 (0.L, 0.L, z, QQ))/4.L);
    }
  
  Ay = TSIL_A (y, QQ);
  
  return (((-(x + y - z) * (x + y - z) -8.L * x * y) * TSIL_I2 (x, y, z, QQ)
	   + (y-z) * (y-z) * TSIL_I2 (0.L, y, z, QQ)
	   + (x-z) * (x-z) * TSIL_I2 (0.L, x, z, QQ)
	   - z * z * TSIL_I2 (0.L, 0.L, z, QQ)
	   + (z - x - y) * Ax * Ay + (y * Ax  + x * Ay) * Az
	   + 2.L * x * y * (Ax + Ay)
	   )/(4.L * x * y));
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (6.23), (6.28) */

TSIL_COMPLEX SUMO_Fvac_FFV (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z, 
			    TSIL_REAL QQ)
{
  TSIL_COMPLEX Ax, Ay, Az;

  Ax = TSIL_A (x, QQ);
  Ay = TSIL_A (y, QQ);

  if (TSIL_FABS (z) < TSIL_TOL) 
    return (2.L * (x * Ax + y * Ay) - (x + y)*(x + y));
  
  Az = TSIL_A (z, QQ);

  return ((((x-y)*(x-y) - 2.L*z*z + (x+y)*z) * TSIL_I2 (x, y, z, QQ)
	   -(x-y)*(x-y) * TSIL_I2 (x, y, 0.L, QQ)
	   + (x - y - 2.L * z) * Ax * Az + (y - x - 2.L * z) * Ay * Az 
	   + 2.L * z * Ax * Ay + 2.L * z * ((z/3.L) - x - y) * Az
	   )/z);
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (6.24) NOTE: MASS INSERTIONS ARE INCLUDED!!!! */

TSIL_COMPLEX SUMO_Fvac_ffV (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z, 
			    TSIL_REAL QQ)
{
  return (6.L * TSIL_CSQRT(x * y) * TSIL_I2 (x, y, z, QQ));
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* hep-ph/0111209 eq. (6.25), (6.29) and (6.30) */

TSIL_COMPLEX SUMO_Fvac_gauge (TSIL_REAL x,
			      TSIL_REAL y,
			      TSIL_REAL z, 
			      TSIL_REAL QQ)
{
  TSIL_REAL tmp;
  TSIL_COMPLEX Ax, Ay, Az;

  if (TSIL_FABS(x) < TSIL_FABS(y)) {tmp = y; y = x; x = tmp;}
  if (TSIL_FABS(y) < TSIL_FABS(z)) {tmp = z; z = y; y = tmp;}
  if (TSIL_FABS(x) < TSIL_FABS(y)) {tmp = y; y = x; x = tmp;}

  Ax = TSIL_A (x, QQ);

  if (TSIL_FABS(y) < TSIL_TOL) 
    {
      return (x * (13.L * TSIL_I2 (0.L, 0.L, x, QQ) - 71.L * Ax/6.L + 4.75L * x));
    }

  Ay = TSIL_A (y, QQ);

  if (TSIL_FABS(z) < TSIL_TOL) 
    {
      return (3.L*x*x + 5.5L*x*y + 3.L*y*y + 5.L*(y*Ax + x*Ay) 
	      - 25.L*(x*Ax + y*Ay)/3.L  
	      + (x*x*(7.L*x + 2.L*y)*TSIL_I2 (0.L, 0.L, x, QQ) 
		 +y*y*(2.L*x + 7.L*y)*TSIL_I2 (0.L, 0.L, y, QQ) 
		 - (7.L*x*x*x - 43.L*x*x*y - 43.L*x*y*y + 7.L*y*y*y)*
		 TSIL_I2 (0.L, x, y, QQ)
		 - (7.L*x*x - 34.L*x*y + 7.L*y*y)*Ax*Ay)/(4.L*x*y));
    }

  Az = TSIL_A (z, QQ);

  return ((-TSIL_I2 (x, y, z, QQ) * (x*x + y*y + z*z - 2.L*(x*y + x*z + y*z))*
	   (x*x + y*y + z*z + 10.L*(x*y + x*z + y*z))
	   + TSIL_I2 (0.L, y, z, QQ) * (y - z) * (y - z) *
	   (y*y + z*z + 10.L*y*z)
	   + TSIL_I2 (0.L, x, z, QQ) * (x - z) * (x - z) *
	   (x*x + z*z + 10.L*x*z)
	   + TSIL_I2 (0.L, x, y, QQ) * (x - y) * (x - y) *
	   (x*x + y*y + 10.L*x*y)
	   + TSIL_I2 (0.L, 0.L, x, QQ) * (x*x*(2.L*y*z - x*x))
	   + TSIL_I2 (0.L, 0.L, y, QQ) * (y*y*(2.L*x*z - y*y))
	   + TSIL_I2 (0.L, 0.L, z, QQ) * (z*z*(2.L*x*y - z*z))
	   + x*(x*x + 14.L * y*z + 9.L * (x*(y + z) - y*y - z*z))*Ay*Az
	   + y*(y*y + 14.L * x*z + 9.L * (y*(x + z) - x*x - z*z))*Ax*Az
	   + z*(z*z + 14.L * x*y + 9.L * (z*(x + y) - x*x - y*y))*Ax*Ay    
	   )/(4.L * x * y * z) 
	  + Ax * (5.5L * (y+z) - 10.L * x/3)
	  + Ay * (5.5L * (x+z) - 10.L * y/3)
	  + Az * (5.5L * (x+y) - 10.L * z/3));
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
