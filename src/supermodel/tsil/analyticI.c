/* Contains the two-loop vacuum integral function I(x,y,z) (eq. 2.20
   in hep-ph/0307101), and derivatives of it. For negative arguments,
   x,y,z have a small negative imaginary part. 
*/

#include "internal.h"

/* ***************************************************************** */
/* Fixed in v1.32 to work properly with negative squared mass args.  */

TSIL_COMPLEX TSIL_I2 (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL Z, TSIL_REAL QQ)
{
  TSIL_REAL tmp, absX, absY, absZ;
  TSIL_COMPLEX DeltaXYZ, sqDeltaXYZ, rp, rm, xiXYZ;
  TSIL_COMPLEX lnbarX, lnbarY, lnbarZ;
  TSIL_COMPLEX lnrp, lnrm, dilogrp, dilogrm;
  TSIL_COMPLEX x, y, z, result;

  if (X < Y) {tmp = Y; Y = X; X = tmp;}
  if (Y < Z) {tmp = Z; Z = Y; Y = tmp;}
  if (X < Y) {tmp = Y; Y = X; X = tmp;}
 
  absX = TSIL_FABS(X);
  if (absX < TSIL_TOL) return TSIL_I20xy (Y, Z, QQ);

  absY = TSIL_FABS(Y);
  if (absY < TSIL_TOL) return TSIL_I20xy (X, Z, QQ);

  absZ = TSIL_FABS(Z);
  if (absZ < TSIL_TOL) return TSIL_I20xy (X, Y, QQ);

  lnbarX = TSIL_LOG(absX/QQ);
  lnbarY = TSIL_LOG(absY/QQ);
  lnbarZ = TSIL_LOG(absZ/QQ);

  x = X; y = Y; z = Z;

  if (X<0) {lnbarX -= I*PI; x -= I*TSIL_EPSILON;}
  if (Y<0) {lnbarY -= I*PI; y -= I*TSIL_EPSILON;}
  if (Z<0) {lnbarZ -= I*PI; z -= I*TSIL_EPSILON;}

  DeltaXYZ = x*x + y*y + z*z - 2.0L*(x*y + x*z + y*z);
  sqDeltaXYZ = TSIL_CSQRT(DeltaXYZ);
  rp = (x + z - y - sqDeltaXYZ)/(2.0L * x);
  rm = (x + y - z - sqDeltaXYZ)/(2.0L * x);

  lnrp = TSIL_CLOG(rp);
  lnrm = TSIL_CLOG(rm);
  dilogrp = TSIL_Dilog(rp);
  dilogrm = TSIL_Dilog(rm);

  xiXYZ = sqDeltaXYZ*(2.0L * (lnrp * lnrm - dilogrp - dilogrm + Zeta2)
                - (lnbarY - lnbarX) * (lnbarZ - lnbarX));
  
  result = 0.5L*((X - Y - Z) * lnbarY * lnbarZ +
                 (Y - X - Z) * lnbarX * lnbarZ +
                 (Z - X - Y) * lnbarX * lnbarY -xiXYZ)
          + 2.0L * (X * lnbarX + Y * lnbarY + Z * lnbarZ)
          - 2.5L * (X + Y + Z);

  if (Z>0) result = TSIL_CREAL(result);

  return(result);
}

/* ***************************************************************** */
/* Fixed in v1.32 to work properly with negative squared mass args.  */

TSIL_COMPLEX TSIL_I20xy (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL QQ)
{
  TSIL_COMPLEX lnbarX, lnbarY, result;
  TSIL_REAL tmp, absX, absY;
  
  if (X < Y) {tmp = Y; Y = X; X = tmp;}

  absX = TSIL_FABS(X);
  if (absX < TSIL_TOL) return TSIL_I200x (Y,QQ);

  absY = TSIL_FABS(Y);
  if (absY < TSIL_TOL) return TSIL_I200x (X,QQ);

  lnbarX = TSIL_LOG(absX/QQ);
  if (X<0) lnbarX -= I*PI;

  if (TSIL_FABS(X-Y) < TSIL_TOL)
    return (X * (-lnbarX*lnbarX + 4.0L * lnbarX - 5.0L));
 
  lnbarY = TSIL_LOG(absY/QQ);
  if (Y<0) lnbarY -= I*PI;

  if (X*Y > 0) {
    result = (Y-X) * (TSIL_Dilog(1 - Y/X) + 0.5L * lnbarX * lnbarX) 
      - Y * lnbarX * lnbarY + 2.0L * (X * lnbarX + Y * lnbarY) 
      - 2.5L * (X + Y);   
  } else {
    result = ((X - Y) * (TSIL_Dilog(Y/X) 
      + (lnbarY - lnbarX) * TSIL_CLOG((X - Y)/QQ)
      + 0.5L * lnbarX * lnbarX - Zeta2) - 2.5L * (X+Y)
      + 2.0L * (X * lnbarX +  Y * lnbarY) - X * lnbarX * lnbarY);
  }

  return (result);
}

/* ***************************************************************** */
/* Fixed in v1.32 to work properly with negative squared mass args.  */

TSIL_COMPLEX TSIL_I200x (TSIL_REAL X,  TSIL_REAL QQ)
{
  TSIL_COMPLEX lnbarX;
  TSIL_REAL absX;
  
  absX = TSIL_FABS(X);
  if (absX < TSIL_TOL) return 0.0L + I*0.0L;
  lnbarX = TSIL_LOG(absX/QQ);
  if (X < 0) lnbarX -= I*PI;

  return (X * (-0.5L*lnbarX*lnbarX + 2.0L * lnbarX - 2.5L - Zeta2));
}

/* ************************************************************** */
/* I2p is the derivative of I2 with respect to its first argument:
   I(x',y,z) */

TSIL_COMPLEX TSIL_I2p (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL QQ)
{
  TSIL_REAL Deltaxyz, tmp, onemyox, onemyox2, onemyox3, onemyox4, onemyox5;
  TSIL_COMPLEX Ax, Ay, Az, X, Y, sqrtx, sqrty, sqrtz;

  if (y < z) {tmp = z; z = y; y = tmp;}

  if (TSIL_FABS(x) < TSIL_TOL) {
    TSIL_Warn("I2p", "I(x',y,z) is undefined for x = 0.");
    return TSIL_Infinity;
  }

  Ax = TSIL_A(x,QQ);
  Ay = TSIL_A(y,QQ);

  if (TSIL_FABS(z) < TSIL_TOL * (TSIL_FABS(x) + TSIL_FABS(y))) {
    onemyox = 1.0L - y/x;
    if (TSIL_FABS(onemyox) > 0.01) 
      return ((TSIL_I20xy(x,y,QQ) + x + y - Ax - Ay + Ax*Ay/x)/(x - y));

    onemyox2 = onemyox*onemyox;
    onemyox3 = onemyox*onemyox2;
    onemyox4 = onemyox2*onemyox2;
    onemyox5 = onemyox2*onemyox3;
   
    return (-0.5L*Ax*Ax/(x*x) - onemyox - onemyox2/4.0L - onemyox3/9.0L 
	    - onemyox4/16.0L - onemyox5/25.0L - onemyox5*onemyox/36.0L 
	    - onemyox5*onemyox2/49.0L - onemyox5*onemyox3/64.0L 
	    - onemyox5*onemyox4/81 - onemyox5*onemyox5/100.0L);
  }

  Az = TSIL_A(z,QQ);

  Deltaxyz = TSIL_Delta(x,y,z);

  if (TSIL_FABS(Deltaxyz/(x*x + y*y)) < TSIL_TOL) {
    X = x; Y = y;
    if (x < 0) X -= I * TSIL_EPSILON;
    if (y < 0) Y -= I * TSIL_EPSILON;
    sqrtx = TSIL_CSQRT(X);
    sqrty = TSIL_CSQRT(Y);
    sqrtz = sqrtx - sqrty;

    /* By design, sqrtz is not necessarily positive. */

    return (Ay/(sqrtx*sqrty) + Az/(sqrtx*sqrtz) - Ax/x 
	    +0.5L*(x*Ay*Az - y*Ax*Az - z*Ax*Ay)/(x*y*z));
  }
  else
    return (((x-y-z)*TSIL_I2(x,y,z,QQ) + (x-y+z)*Ax*Ay/x 
	     + (x+y-z)*Ax*Az/x -2.0L*Ay*Az 
	     + (y+z-x)*(Ax+Ay+Az) + x*x - (y+z)*(y+z))/Deltaxyz);
}

/* ****************************************************************** */
/* I2p2 gives the second derivative of I2 with respect to its first   */
/* argument: I(x'',y,z).                                              */

TSIL_COMPLEX TSIL_I2p2 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL QQ)
{
  TSIL_REAL tmp;

  if (TSIL_FABS(x) < TSIL_TOL) {
    TSIL_Warn("I2p2", "I2(x'',y,z) is undefined for x=0.");
    return TSIL_Infinity;
  }
  
  if (y < z) {tmp = y; y = z; z = tmp;}
  if (y > x)
    return TSIL_xI2p2(y,x,z,QQ)/x;
  else
    return TSIL_xI2p2(x,y,z,QQ)/x;
}


/* ****************************************************************** */
/* The following function is useful because it is symmetric under     */
/* interchange of x,y,z.                                              */

TSIL_COMPLEX TSIL_xI2p2 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL QQ)
{
  TSIL_REAL tmp, Deltaxyz, onemyox;
  TSIL_COMPLEX Ax, Ay, Az, X, Y, Z, sqrtx, sqrty, sqrtz;
  TSIL_REAL onemyox2, onemyox3, onemyox4, onemyox5;
  TSIL_REAL x2, x3, y2, y3, z2, z3;

  if (x < y) {tmp = x; x = y; y = tmp;}
  if (x < z) {tmp = x; x = z; z = tmp;}
  if (y < z) {tmp = y; y = z; z = tmp;}

  if (TSIL_FABS(x) < TSIL_TOL) {
    TSIL_Warn("TSIL_xI2p2", "x I(x'',y,z) is undefined for x = y = z = 0.");
    return TSIL_Infinity;
  }

  if (TSIL_FABS(z/x) < TSIL_TOL) {
    onemyox = 1.0L - y/x;
    if (TSIL_FABS(onemyox) > 0.01)
      return (TSIL_A(x,QQ)-TSIL_A(y,QQ))/(y-x);

    onemyox2 = onemyox*onemyox;
    onemyox3 = onemyox*onemyox2;
    onemyox4 = onemyox2*onemyox2;
    onemyox5 = onemyox2*onemyox3;

    return (-TSIL_Ap(x,QQ) + onemyox/2.0L + onemyox2/6.0L + onemyox3/12.0L 
	    + onemyox4/20.0L + onemyox5/30.0L + onemyox5*onemyox/42.0L 
	    + onemyox5*onemyox2/56.0L + onemyox5*onemyox3/72.0L 
	    + onemyox5*onemyox4/90.0L + onemyox5*onemyox5/110.0L); 
  }

  Ax = TSIL_A(x,QQ);
  Ay = TSIL_A(y,QQ);
  Az = TSIL_A(z,QQ);
  Deltaxyz = TSIL_Delta(x,y,z);

  if (TSIL_FABS(Deltaxyz/(x*x)) < TSIL_TOL) {
    X = x; Y = y; Z = z;
    if (x < 0) X -= I * TSIL_EPSILON;
    if (y < 0) Y -= I * TSIL_EPSILON;
    if (z < 0) Z -= I * TSIL_EPSILON;
    sqrtx = TSIL_CSQRT(X);
    sqrty = TSIL_CSQRT(Y);
    sqrtz = TSIL_CSQRT(Z);
    return ((Az/(sqrtx*sqrty) + Ay/(sqrtx*sqrtz) - Ax/(sqrty*sqrtz) 
     - 1.0L)/3.0L);
  }

  x2 = x*x;
  x3 = x*x2;
  y2 = y*y;
  y3 = y*y2;
  z2 = z*z;
  z3 = z*z2;

  return ((-4.L*x2*y*z - 4.L*x*y2*z + Ay*Az*(2.L*x2 - 2.L*x*y - 2.L*x*z) +
	   Ax*Az*(-2.L*x*y + 2.L*y2 - 2.L*y*z) - 4.L*x*y*z2 +
	   Ax*Ay*(-2.L*x*z - 2.L*y*z + 2.L*z2) +
	   Az*(x3 - x2*y - x*y2 + y3 - 3.L*x2*z + 2.L*x*y*z 
	       - 3.L*y2*z + 3.L*x*z2 + 3.L*y*z2 - z3) + 
	   Ay*(x3 - 3.L*x2*y + 3.L*x*y2 - y3 - x2*z + 2.L*x*y*z 
	       + 3.L*y2*z - x*z2 - 3.L*y*z2 + z3) +
	   Ax*(-x3 + 3.L*x2*y - 3.L*x*y2 + y3 + 3.L*x2*z + 2.L*x*y*z 
	       - y2*z - 3.L*x*z2 - y*z2 + z3) - 4.L*x*y*z*TSIL_I2(x, y, z, QQ))/
          (Deltaxyz*Deltaxyz));
}

/* ****************************************************************** */
/* I2pp is the derivative of I2p with respect to its second argument: */
/* I(x',y',z)                                                         */

TSIL_COMPLEX TSIL_I2pp (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL QQ)
{
  TSIL_COMPLEX Ax, Ay, Az, X, Y, Z, sqrtx, sqrty, sqrtz;
  TSIL_REAL tmp, Deltaxyz, onemyox;
  TSIL_REAL onemyox2, onemyox3, onemyox4, onemyox5;
  TSIL_REAL x2, x3, y2, y3, z2, z3, xpymz;

  if (x < y) {tmp = x; x = y; y = tmp;}

  if (TSIL_FABS(y) < TSIL_TOL) {
    TSIL_Warn("I2pp", "I(x',y',z) is undefined for x = 0.");
    return TSIL_Infinity;
  }
  
  if (TSIL_FABS(z) < TSIL_TOL) {
    onemyox = 1.0L - y/x;
    if (TSIL_FABS(onemyox) > 0.01) 
      return (TSIL_Ap(x,QQ) - TSIL_Ap(y,QQ))/(x-y);

    onemyox2 = onemyox*onemyox;
    onemyox3 = onemyox*onemyox2;
    onemyox4 = onemyox2*onemyox2;
    onemyox5 = onemyox2*onemyox3;

    return (1.0L + onemyox/2.0L + onemyox2/3.0L + onemyox3/4.0L 
	    + onemyox4/5.0L + onemyox5/6.0L + onemyox5*onemyox/7.0L 
	    + onemyox5*onemyox2/8.0L + onemyox5*onemyox3/9.0L 
	    + onemyox5*onemyox4/10.0L + onemyox5*onemyox5/11.0L)/x;  
  }

  Ax = TSIL_A(x,QQ);
  Ay = TSIL_A(y,QQ);
  Az = TSIL_A(z,QQ);
  Deltaxyz = TSIL_Delta(x,y,z);
  
  if (TSIL_FABS(Deltaxyz/(x*x + z*z)) < TSIL_TOL) {
    X = x; Y = y; Z = z;
    if (x < 0) X -= I * TSIL_EPSILON;
    if (y < 0) Y -= I * TSIL_EPSILON;
    if (z < 0) Z -= I * TSIL_EPSILON;
    sqrtx = TSIL_CSQRT(X);
    sqrty = TSIL_CSQRT(Y);
    sqrtz = TSIL_CSQRT(Z);
    return ((1.0L + 0.5L*(Az-Ax-3.0L*Ay)/y + (Ax-Ay)/(sqrty*sqrtz) 
	     +(sqrtz/sqrty))/(3.0L*x));
  }
   
  x2 = x*x;   
  x3 = x*x2;
  y2 = y*y; 
  y3 = y*y2;
  z2 = z*z;    
  z3 = z*z2;
  xpymz = x + y - z;

  return ((Ay*Az*x*xpymz*(y + z - x) + Ax*Az*xpymz*y*(x - y + z) 
	   + Ax*Ay*xpymz*xpymz*z +  2.0L*x*xpymz*y*z*(x + y + z) 
	   - 2.0L*Az*x*y*(x2 - 2*x*y + y2 - x*z - y*z) 
	   - Ay*x*(x3 - 3.0L*x2*y + 3.0L*x*y2 - y3 - 3.0L*x2*z + 4.0L*x*y*z 
		   + 3.0L*y2*z + 3.0L*x*z2 - y*z2 - z3) 
	   + Ax*y*(x3 - 3.0L*x2*y + 3.0L*x*y2 - y3 - 3.0L*x2*z - 4.0L*x*y*z 
		   + 3.0L*y2*z + x*z2 - 3.0L*y*z2 + z3) 
	   + 2.0L*x*xpymz*y*z*TSIL_I2(x, y, z, QQ))/(x*y*Deltaxyz*Deltaxyz));
}


/* ****************************************************************** */
/* Third derivative of I2 with respect to its first argument:         */
/*   I(x''',y,z)                                                      */

TSIL_COMPLEX TSIL_I2p3 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL QQ)
{
  TSIL_COMPLEX Ax, Ay, Az, X, Y, Z, sqrtx, sqrty, sqrtz;
  TSIL_REAL tmp, Deltaxyz;
  TSIL_REAL onemyox, onemyox2, onemyox3, onemyox4, onemyox5;
  TSIL_REAL xm1, xm2, x2, x3, x4;
  TSIL_REAL y2, y3, y4, y5, z2, z3, z4, z5;

  if (TSIL_FABS(x) < TSIL_TOL) {
    TSIL_Warn("I2p3", "I2(x''',y,z) is undefined for x=0.");
    return TSIL_Infinity;
  }

  if (y < z) {tmp = y; y = z; z = tmp;}

  Ax = TSIL_A(x,QQ);
  Ay = TSIL_A(y,QQ);

  if (TSIL_FABS(z/x) < TSIL_TOL) {
    onemyox = 1.0L - y/x;
    if (TSIL_FABS(onemyox) > 0.01) 
      return ((x*(y - x) + x*Ax + (y - 2.L*x)*Ay)/(x*x*(x-y)*(x-y)));

    onemyox2 = onemyox*onemyox;
    onemyox3 = onemyox*onemyox2;
    onemyox4 = onemyox2*onemyox2;
    onemyox5 = onemyox2*onemyox3;

    return ((TSIL_Ap(x,QQ) -0.5L -2.L*onemyox/3.L - onemyox2/4.L 
	     - 2.L*onemyox3/15.L - onemyox4/12.L -  2.L*onemyox5/35.L 
	     - onemyox5*onemyox/24.L - 2.L*onemyox5*onemyox2/63.L 
	     - onemyox5*onemyox3/40.L - 2.L*onemyox5*onemyox4/99.L 
	     - onemyox5*onemyox5/60.L)/(x*x));
  }

  Az = TSIL_A(z,QQ);
  Deltaxyz = TSIL_Delta(x,y,z);

  if (TSIL_FABS(Deltaxyz/(x*x + y*y)) < TSIL_TOL) {
    X = x; Y = y; 
    if (x < 0) X -= I * TSIL_EPSILON;
    if (y < 0) Y -= I * TSIL_EPSILON;
    sqrtx = TSIL_CSQRT(X);
    sqrty = TSIL_CSQRT(Y);
    sqrtz = sqrtx - sqrty;

    /* By design, sqrtz is not necessarily positive. */          

    return ( (Ax/(y*z) - 2.L*(sqrty*sqrtz + y + z)/(sqrty*sqrtz*x) 
	      - Az*(5.L*sqrty + sqrtz)/(sqrtx*x*y) 
	      - Ay*(sqrty + 5.L*sqrtz)/(sqrtx*x*z) )/(10.L*x) );
  }

  xm1 = 1/x;
  xm2 = 1/(x*x);
  x2 = x*x;
  x3 = x*x2;
  x4 = x2*x2;
  y2 = y*y;
  y3 = y*y2;
  y4 = y2*y2;
  y5 = y2*y3;
  z2 = z*z;
  z3 = z*z2;
  z4 = z2*z2;
  z5 = z2*z3;

  return ( (-x4 - 10.L*x2*y2 + 10.L*x*y3 - 5.L*y4 + xm1*y5 
	    + 6.L*x*y2*z - 8.L*y3*z - 3.L*xm1*y4*z + 5.L*x3*(y + z) 
	    - 10.L*x2*z2 + 6.L*x*y*z2 - 22.L*y2*z2 + 2.L*xm1*y3*z2 
	    - 6.L*Ay*Az*(x2 + y2 + 2.L*y*z - 2.L*x*(y + z) + z2) 
	    + 12.L*TSIL_I2(x,y,z,QQ)*(x*y*z - y2*z - y*z2) 
	    + 6.L*Ax*Az*(x*y - 2.L*y2 + xm1*(y3 - y*z2)) 
	    + 10.L*x*z3 - 8.L*y*z3 + 2.L*xm1*y2*z3 
	    + 6.L*Ax*Ay*(x*z - 2.L*z2 + xm1*(-y2*z + z3)) 
	    - 5.L*z4 - 3.L*xm1*y*z4 + Ax*(x3 + 6.L*x*y2 - 4.L*y3 
            + xm1*y4 + 4.L*x*y*z - 8.L*y2*z + 8.L*xm1*y3*z 
            - 4.L*x2*(y + z) + 6.L*x*z2 - 8.L*y*z2 - 18.L*xm1*y2*z2 
            - 4.L*z3 + 8.L*xm1*y*z3 + xm1*z4) + xm1*z5 
	    + Ay*(-2.L*x3 - 16.L*x*y2 + 14.L*y3 - 6.L*xm1*y4 + xm2*y5 
            - 14.L*x*y*z + 8.L*y2*z + 10.L*xm1*y3*z - 5.L*xm2*y4*z 
            + x2*(9.L*y + z) + 10.L*x*z2 + 18.L*y*z2 + 6.L*xm1*y2*z2 
            + 10.L*xm2*y3*z2 - 16.L*z3 - 18.L*xm1*y*z3 - 10.L*xm2*y2*z3 
            + 8.L*xm1*z4 + 5.L*xm2*y*z4 - xm2*z5) +  Az*(-2.L*x3 
            + 10.L*x*y2 - 16.L*y3 + 8.L*xm1*y4 - xm2*y5 - 14.L*x*y*z 
            + 18.L*y2*z - 18.L*xm1*y3*z + 5.L*xm2*y4*z + x2*(y + 9.L*z) 
            - 16.L*x*z2 + 8.L*y*z2 + 6.L*xm1*y2*z2 - 10.L*xm2*y3*z2 
            + 14.L*z3 + 10.L*xm1*y*z3 + 10.L*xm2*y2*z3 - 6.L*xm1*z4 
            - 5.L*xm2*y*z4 + xm2*z5))/(Deltaxyz*Deltaxyz*Deltaxyz));
}
