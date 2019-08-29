
/** \file numerics.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "numerics.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_LOOPTOOLS
#include "clooptools.h"
#endif

int cubicRoots(double alpha, double beta, double gamma, double delta,
	       DoubleVector & ans) {
  if (fpclassify(alpha) != FP_ZERO) {
    int a =
      cubicRootsInside(beta / alpha, gamma / alpha, delta / alpha, ans);
    if (a > 2) {
      /// order the roots in increasing absolute order
      if (fabs(ans(1)) > fabs(ans(2))) ans.swap(1, 2);
      if (fabs(ans(2)) > fabs(ans(3))) ans.swap(2, 3);
      if (fabs(ans(1)) > fabs(ans(2))) ans.swap(1, 2);      
    }
    return a;
  }
  if (fpclassify(beta) == FP_ZERO) {
    /// It's a linear equation
    if (fpclassify(gamma) != FP_ZERO) {
	ans.setEnd(1);
	ans(1) = -delta / gamma;
	return 1;
      }
  }
  double sgn_delta = 1.;
  if (delta < 0.) sgn_delta = -1.;
  /// it's a quadratic really
  double q = -0.5 *
    (delta + sgn_delta * sqrt(sqr(delta) - 4.0 * beta * gamma));
  ans.setEnd(2);
  ans(1) = q / beta;
  ans(2) = gamma / q;
  
  return ans.displayEnd();
}

int cubicRootsInside(double a, double b, double c, 
	       DoubleVector & ans) {
  ans.setEnd(3);
  double Q = (sqr(a) - 3.0 * b) / 9.0;
  double R = (2.0 * sqr(a) * a - 9.0 * a * b + 27.0 * c) / 54.;
  double Qcub = Q * Q * Q;
  double Rsq  = R * R;
  double rtQ  = sqrt(Q);

  if (Rsq <= Qcub) {
    /// Three real roots
    double theta = acos(R / sqrt(Qcub));
    ans(1) = -2.0 * rtQ * cos(theta / 3.0) - a / 3.0;
    ans(2) = -2.0 * rtQ * cos((theta + 2.0 * PI) / 3.0)
      - a / 3.0;
    ans(3) = -2.0 * rtQ * cos((theta - 2.0 * PI) / 3.0)
      - a / 3.0;
  } else {
    /// only one real root
    double sgn_R = 1;
    if (R < 0.) sgn_R = -1.;
    double A = -sgn_R * cbrt(abs(R) + sqrt(Rsq - Qcub));
    double B = 1.;
    if (fpclassify(A) != FP_ZERO) B = Q / A;
    ans.setEnd(1);
    ans(1) = A + B - a / 3.0;
  }
  return ans.displayEnd();
}

double dgauss(double (*f)(double x), double a, double b, double eps) {
  static DoubleVector w(12), x(12);
  /// set the initial values if they are uninitialised
  if (w(1) < 1.0e-5) {
    w(1)  =    0.1012285362903762591525313543;
    w(2)  =    0.2223810344533744705443559944;
    w(3)  =    0.3137066458778872873379622020;
    w(4)  =    0.3626837833783619829651504493;
    w(5)  =    0.0271524594117540948517805725;
    w(6)  =    0.0622535239386478928628438370;
    w(7)  =    0.0951585116824927848099251076;
    w(8)  =    0.1246289712555338720524762822;
    w(9)  =    0.1495959888165767320815017305;
    w(10) =    0.1691565193950025381893120790;
    w(11) =    0.1826034150449235888667636680;
    w(12) =    0.1894506104550684962853967232;
    x(1)  =    0.9602898564975362316835608686;
    x(2)  =    0.7966664774136267395915539365;
    x(3)  =    0.5255324099163289858177390492;
    x(4)  =    0.1834346424956498049394761424;
    x(5)  =    0.9894009349916499325961541735;
    x(6)  =    0.9445750230732325760779884155;
    x(7)  =    0.8656312023878317438804678977;
    x(8)  =    0.7554044083550030338951011948;
    x(9)  =    0.6178762444026437484466717640; 
    x(10)  =   0.4580167776572273863424194430;
    x(11)  =   0.2816035507792589132304605015;
    x(12)  =   0.0950125098376374401853193354;
  }
  const double constant = 1.0e-25;
  double delta = constant * fabs(a - b);
  double dgauss = 0.;
  double aa = a;
 lab5: double y = b - aa;
  if (fabs(y) <= delta) return dgauss;
 lab2: double bb = aa + y;
  double c1 = 0.5 * (aa + bb);
  double c2 = c1 - aa;
  double s8 = 0.;
  double s16 = 0.;
  for (int i=1; i<=4; i++) {
    double u = x(i) * c2;
    s8 = s8 + w(i) * (f(c1 + u) + f(c1 - u));
  }
  for (int i=5; i<=12; i++) {
    double u = x(i) * c2;
    s16 = s16 + w(i) * (f(c1 + u) + f(c1 - u));
  }
  s8 = s8 * c2;
  s16 = s16 * c2;
  if (fabs(s16 - s8) > eps * (1. + fabs(s16))) goto lab4;
  dgauss = dgauss + s16;
  aa = bb;
  goto lab5;
 lab4: y = 0.5 * y;
  if (fabs(y) > delta) goto lab2;
  cout << "Too high accuracy required in numerics.cpp:dgauss\n";
  dgauss = 0.;
  throw("Too high accuracy required in numerics.cpp:dgauss\n");
  return dgauss;
}

double dgauss(double (*f)(double x, const DoubleVector & v),
	      const DoubleVector & v, double a, double b, double eps) {
  static DoubleVector w(12), x(12);
  /// set the initial values if they are uninitialised
  if (w(1) < 1.0e-5) {
    w(1)  =    0.1012285362903762591525313543;
    w(2)  =    0.2223810344533744705443559944;
    w(3)  =    0.3137066458778872873379622020;
    w(4)  =    0.3626837833783619829651504493;
    w(5)  =    0.0271524594117540948517805725;
    w(6)  =    0.0622535239386478928628438370;
    w(7)  =    0.0951585116824927848099251076;
    w(8)  =    0.1246289712555338720524762822;
    w(9)  =    0.1495959888165767320815017305;
    w(10) =    0.1691565193950025381893120790;
    w(11) =    0.1826034150449235888667636680;
    w(12) =    0.1894506104550684962853967232;
    x(1)  =    0.9602898564975362316835608686;
    x(2)  =    0.7966664774136267395915539365;
    x(3)  =    0.5255324099163289858177390492;
    x(4)  =    0.1834346424956498049394761424;
    x(5)  =    0.9894009349916499325961541735;
    x(6)  =    0.9445750230732325760779884155;
    x(7)  =    0.8656312023878317438804678977;
    x(8)  =    0.7554044083550030338951011948;
    x(9)  =    0.6178762444026437484466717640; 
    x(10)  =   0.4580167776572273863424194430;
    x(11)  =   0.2816035507792589132304605015;
    x(12)  =   0.0950125098376374401853193354;
  }
  const double constant = 1.0e-25;
  double delta = constant * fabs(a - b);
  double dgauss = 0.;
  double aa = a;
 lab5: double y = b - aa;
  if (fabs(y) <= delta) return dgauss;
 lab2: double bb = aa + y;
  double c1 = 0.5 * (aa + bb);
  double c2 = c1 - aa;
  double s8 = 0.;
  double s16 = 0.;
  for (int i=1; i<=4; i++) {
    double u = x(i) * c2;
    s8 = s8 + w(i) * (f(c1 + u, v) + f(c1 - u, v));
  }
  for (int i=5; i<=12; i++) {
    double u = x(i) * c2;
    s16 = s16 + w(i) * (f(c1 + u, v) + f(c1 - u, v));
  }
  s8 = s8 * c2;
  s16 = s16 * c2;
  if (fabs(s16 - s8) > eps * (1. + fabs(s16))) goto lab4;
  dgauss = dgauss + s16;
  aa = bb;
  goto lab5;
 lab4: y = 0.5 * y;
  if (fabs(y) > delta) goto lab2;
  cout << "Too high accuracy required in numerics.cpp:dgauss\n";
  dgauss = 0.;
  throw("Too high accuracy required in numerics.cpp:dgauss\n");
  return dgauss;
}


double accurateSqrt1Plusx(double x) {
  if (x > 1.) return sqrt(1.0 + x);
  
  /// calculate approximate number of terms we're going to need to get 16
  /// digit precision
  int num = -16.0 / log10(x); 

  cout << "Number of terms: " << num << endl;

  double factor = 0.5 * x;
  double ans = 1. + factor;
  double pow = -0.5;

  for (int i=2; i<=num; i++) {
    factor *= x * pow / double(i);
    pow -= 1.;
    ans += factor;
  }
  return ans;
}

double log1minusx(double x) {
  if (x > 1.) return 0.; 
  else if (close(1., x, EPSTOL)) return numberOfTheBeast;
  else if (fabs(x) > 0.125) return log(1. - x);
  else if (fpclassify(x) == FP_ZERO) return 0.;
  double test = -x; int i = 1; double l1mx = -x;
  /// Find largest power that we need from the expansion
  do { 
    i++;
    test = l1mx; 
    l1mx = l1mx - pow(x, i) / double(i);
    if (close(l1mx, test, 1.0e-15)) break;
  } while (!close(l1mx, test, 1.0e-15));
  l1mx = 0.;
  int nMax = i;
  for (i=nMax; i>=1; i--) {
    l1mx -= pow(x, i) / double(i);
  }
  return l1mx;
}

// returns >0 if there's a problem:
int integrateOdes(DoubleVector & ystart, double from, double to, double eps,
	      double h1, double hmin, 
	      DoubleVector (*derivs)(double, const DoubleVector &),
	      int (*rkqs)(DoubleVector & y, const DoubleVector & dydx, double
			   *x, double htry, double eps, DoubleVector & yscal,
			   double *hdid, double *hnext, 
			   DoubleVector (*derivs)(double, const DoubleVector
						  &)) ) {  
  int nvar =  ystart.displayEnd();
  int nstp, i;
  double x, hnext, hdid, h;
  DoubleVector yscal(nvar), y(ystart), dydx(nvar);
  
  x = from;
  h = sign(h1, to - from);
  
  const int MAXSTP = 400;
  const double TINY = 1.0e-16;

  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    dydx = (*derivs)(x, y);
    for (i = 1; i <= nvar; i++)
      yscal(i) = fabs(y(i)) + fabs(dydx(i) * h) + TINY;
    if ((x + h - to) * (x + h - from) > 0.0) h = to - x;
    int smallStep = (*rkqs)(y, dydx, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if (smallStep) return 1;

    if ((x - to) * (to - from) >= 0.0) {
      for (i = 1; i<= nvar; i++) ystart(i) = y(i);
      return 0;
    }
      
    if (fabs(hnext) <= hmin) {
      nstp = MAXSTP; // bail out
      if (PRINTOUT > 1) {
	cout << "Step size too small in numerics.cpp:integrateOdes\n";
	cout << "**********x = " << x << "*********\n";
	for (i = 1;i<= nvar;i++) 
	  cout << "y(" << i << ") = " << y(i) << " dydx(" << i <<
	    ") = " << dydx(i) << endl;
	cout.flush();
      }
    }
    
    h = hnext;
  }
  
  if (PRINTOUT > 1) {
    cout << "Bailed out of numerics.cpp:too many steps in integrateOdes\n";
    cout << "**********x = " << x << "*********\n";
    for (i = 1;i<= nvar;i++) 
      cout << "y(" << i << ") = " << y(i) << " dydx(" << i <<
	") = " << dydx(i) << endl;
    cout.flush();
  }
  
  return 1;
}

double lnLPoisson(unsigned k, double lambda) {
  if (lambda < 6.66e-66) throw("In lnPoisson: lambda not positive\n");
  double lnL = -lambda + double(k) * log(lambda);
  for (unsigned i=2; i<=k; i++) lnL -= log(double(i));

  return lnL;
}

double LPoisson(unsigned k, double lambda) {
  if (lambda < 6.66e-66) throw("In lnPoisson: lambda not positive\n");
  double L = exp(-lambda) * pow(lambda,int(k));
  for (unsigned i=2; i<=k; i++) L /= double(i);

  return L;
}

int odeStepper(DoubleVector & y, const DoubleVector & dydx, double *x, double
		htry, double eps, DoubleVector & yscal, double *hdid, 
		double *hnext,		
		DoubleVector (*derivs)(double, const DoubleVector &))
{
  const double SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;

  int i, n = y.displayEnd();
  double errmax, h, htemp, xnew;
  
  DoubleVector yerr(n), ytemp(n);
  h = htry;
  for (;;) {
    rungeKuttaStep(y, dydx, *x, h, ytemp, yerr, derivs);
    errmax = 0.0;

    for (i = 1; i<= n;i++) errmax = maximum(errmax, fabs(yerr(i) / yscal(i)));
    errmax  /= eps;

    if (errmax <= 1.0) break;
    htemp = SAFETY * h * pow(errmax, PSHRNK);
    h = (h >= 0.0 ? maximum(htemp ,0.1 * h) : minimum(htemp, 0.1 * h));
    xnew = (*x) + h;
    if (xnew == *x) {
	if (PRINTOUT) {
	cout << "At x = " << *x;
	cout << ",stepsize underflow in odeStepper" << flush << endl;
	}
	return 1;
      }
  }
  if (errmax > ERRCON) *hnext = SAFETY * h * pow(errmax,PGROW);
  else *hnext = 5.0 * h;
  *x += (*hdid = h);
  y = ytemp;
  return 0;
}

void rungeKuttaStep(const DoubleVector & y, const DoubleVector & dydx, 
	     double x, double h, DoubleVector & yout, DoubleVector & yerr, 
	     DoubleVector (*derivs)(double, const DoubleVector &)) {
  int i;
  const static double a2 = 0.2,a3 = 0.3,a4 = 0.6,a5 = 1.0,a6 = 0.875,b21 =
    0.2,b31 = 3.0 / 40.0,b32 = 9.0 / 40.0,b41 = 0.3,b42 = -0.9,b43 = 1.2,
    b51 = -11.0 / 54.0, b52 = 2.5,b53 = -70.0 / 27.0,b54 = 35.0 / 27.0,
    b61 = 1631.0 / 55296.0,b62 = 175.0 / 512.0,b63 = 575.0 / 13824.0,
    b64 = 44275.0 / 110592.0,b65 = 253.0 / 4096.0,c1 = 37.0 / 378.0,
    c3 = 250.0 / 621.0,c4 = 125.0 / 594.0,c6 = 512.0 / 1771.0,
    dc5 = -277.00 / 14336.0;
  const double dc1 = c1-2825.0 / 27648.0,dc3 = c3-18575.0 / 48384.0,
    dc4 = c4-13525.0 / 55296.0,dc6 = c6-0.25;
  
  int n = y.displayEnd();

  DoubleVector ytemp(b21 * h * dydx.display() + y.display());
  DoubleVector ak2((*derivs)(x + a2 * h, ytemp));

  // Allowing piece-wise calculating of ytemp for speed reasons
  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b31 * dydx.display(i) + b32 * ak2(i));
  DoubleVector ak3((*derivs)(x + a3 * h, ytemp));

  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b41 * dydx.display(i) + b42 * ak2(i) + b43
				   * ak3(i));
  DoubleVector ak4((*derivs)(x+a4*h,ytemp));

  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b51 * dydx.display(i) + b52 * ak2(i) + b53
				   * ak3(i) + b54 * ak4(i));
  DoubleVector ak5((*derivs)(x + a5 * h, ytemp));

  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b61 * dydx.display(i) + b62 * ak2(i) + b63
				   * ak3(i) + b64 * ak4(i) + b65 * ak5(i));
  DoubleVector ak6((*derivs)(x + a6 * h, ytemp));

  for (i = 1; i<= n; i++)
    yout(i) = y.display(i) + h * (c1 * dydx.display(i) + c3 * ak3(i) + c4 *
				  ak4(i) + c6 * ak6(i));
  for (i = 1; i<= n; i++) 
    yerr(i) = h * (dc1 * dydx.display(i) + dc3 * ak3(i) + 
		   dc4 * ak4(i) + dc5 * ak5(i) + dc6 * ak6(i));
}

double calcDerivative(double (*func)(double), double x, double h, double
		      *err){
  const double CON = 1.4, CON2 = CON * CON, BIG = 1.0e30, 
    SAFE = 2.0; 
  const int NTAB = 10;
  
  int i, j;
  double errt, fac, hh, ans = 0.0;
  
  if (h == 0.0) throw "h must be nonzero in numerics.cpp:calcDerivative";


  DoubleMatrix a(NTAB, NTAB);
  hh = h;
  a(1, 1) = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
  *err = BIG;
  for (i=2; i<=NTAB; i++) {
    hh /= CON;
    a(1, i) = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
    fac = CON2;
    for (j=2; j<=i; j++) {
      a(j, i) = (a(j-1, i) * fac - a(j-1, i-1)) / (fac - 1.0);
      fac = CON2 * fac;
      errt = maximum(fabs(a(j, i) - a(j-1, i)), fabs(a(j, i) - a(j-1, i-1)));
      if (errt <= *err) {
	*err = errt;
	ans = a(j, i);
      }
    }
    if (fabs(a(i, i) - a(i-1, i-1)) >= SAFE * (*err)) break;
  }

  return ans;
}

double calcDerivative(double (*func)(double, void*), double x, double h, double
		      *err, void* params){
  const double CON = 1.4, CON2 = CON * CON, BIG = 1.0e30, 
    SAFE = 2.0; 
  const int NTAB = 10;
  
  int i, j;
  double errt, fac, hh, ans = 0.0;
  
  if (h == 0.0) throw "h must be nonzero in numerics.cpp:calcDerivative";


  DoubleMatrix a(NTAB, NTAB);
  hh = h;
  a(1, 1) = ((*func)(x + hh, params) - (*func)(x - hh, params)) / (2.0 * hh);
  *err = BIG;
  for (i=2; i<=NTAB; i++) {
    hh /= CON;
    a(1, i) = ((*func)(x + hh, params) - (*func)(x - hh, params)) / (2.0 * hh);
    fac = CON2;
    for (j=2; j<=i; j++) {
      a(j, i) = (a(j-1, i) * fac - a(j-1, i-1)) / (fac - 1.0);
      fac = CON2 * fac;
      errt = maximum(fabs(a(j, i) - a(j-1, i)), fabs(a(j, i) - a(j-1, i-1)));
      if (errt <= *err) {
	*err = errt;
	ans = a(j, i);
      }
    }
    if (fabs(a(i, i) - a(i-1, i-1)) >= SAFE * (*err)) break;
  }

  return ans;
}

inline void shft2(double & a, double & b, double c) { a = b; b = c; }

inline void shft3(double & a, double & b, double & c, double d) { 
  a = b; b = c; c = d;
}

double findMinimum(double ax, double bx, double cx, double (*f)(double),
		   double tol, double &xmin)
{
  const double R = 0.61803399, C = 1.0 - R;
  double f1, f2, x0, x1, x2, x3;
  
  x0 = ax; 
  x3 = cx; 
  if (fabs(cx - bx) > fabs(bx - ax)) {
    x1 = bx; 
    x2 = bx + C * (cx - bx); 
  } else {
    x2 = bx; 
    x1 = bx - C * (bx - ax); 
  }
  f1 = (*f)(x1); 
  f2 = (*f)(x2); 
  while (fabs(x3 - x0) > tol * (fabs(x1) + fabs(x2))) {
    if (f2 < f1) {
      shft3(x0, x1, x2, R * x1 + C * x3);
      shft2(f1, f2, (*f)(x2));
    } else {
      shft3(x3, x2, x1, R * x2 + C * x0);
      shft2(f2, f1, (*f)(x1));
	}
  }
  if (f1 < f2) {
    xmin = x1; 
    return f1; 
  } else {
    xmin = x2; 
    return f2; 
  }
}



DoubleVector dd(double x, const DoubleVector & /* y */) {
  DoubleVector dydx(1);
  dydx(1) = -integrandThreshbnr(x);
  return dydx;
}

double integrandThreshbnr(double x) {
  return fnfn(x).real();
}

// Integration routine needs these variables
static double mtInt, pInt, m1Int, m2Int;
static int nInt;

Complex fnfn(double x) {
  const static Complex iEpsilon(0.0, TOLERANCE * 1.0e-20);
  
  double xn = 1.0;
  int i; for (i=1; i<=nInt; i++) xn *= x;
  return xn * 
    log( ((1 - x) * sqr(m1Int) + x * sqr(m2Int) - x * (1 - x) *
	  sqr(pInt) - iEpsilon)
	 / sqr(mtInt));
}

double bIntegral(int n1, double p, double m1, double m2, double mt) {
  // Set global variables so that integration function can access them
  nInt = n1; pInt = p; m1Int = m1; m2Int = m2; mtInt = mt;
  double from = 0.0, to = 1.0, hmin = TOLERANCE * 1.0e-5;

  double ans = dgauss(integrandThreshbnr, from, to, hmin);
  return -ans;
}

double fB(const Complex & a) {
  /// First, special cases at problematic points
  double x = a.real();
  if (fabs(x) < EPSTOL) {
    double ans = -1. - x + sqr(x) * 0.5;
    return ans;
  }
  if (close(x, 1., EPSTOL)) return -1.;

  Complex ans = log(1. - a) - 1. - a * log(1.0 - 1.0 / a);

  return ans.real();
}

Complex fBc(const Complex & a) {
  /// First, special cases at problematic points
  double x = a.real();
  if (fabs(x) < EPSTOL) {
    double ans = -1. - x + sqr(x) * 0.5;
    return ans;
  }
  if (close(x, 1., EPSTOL)) return -1.;

  Complex ans = log(1. - a) - 1. - a * log(1.0 - 1.0 / a);

  return ans;
}
  
/*
  Analytic expressions follow for above integrals: sometimes useful!
  From hep-ph/9606211
  Note it returns the REAL PART ONLY. 
*/
double b0(double p, double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  if (fabs(m1) < 1.0e-30) m1 = 1.0e-30;
  if (fabs(m2) < 1.0e-30) m2 = 1.0e-30;  
  double b0l = B0(p*p, m1*m1, m2*m2).real();
  return B0(p*p, m1*m1, m2*m2).real();
#endif

  /// Avoids IR infinities
  if (close(p, 0.0, EPSTOL) && close(m1, 0.0, EPSTOL)
      && close(m2, 0.0, EPSTOL))
    return 0.0;

  double ans  = 0.;
  double mMin = minimum(fabs(m1), fabs(m2));
  double mMax = maximum(fabs(m1), fabs(m2));

  double pSq = sqr(p), mMinSq = sqr(mMin), mMaxSq = sqr(mMax);
  double s = 0.;
  /// Try to increase the accuracy of s
  double dmSq = mMaxSq - mMinSq;
  s = pSq + dmSq;

  double pTest = sqr(p) / sqr(mMax);
  /// Decides level at which one switches to p=0 limit of calculations
  const double pTolerance = 1.0e-6; 

  /// p is not 0  
  if (pTest > pTolerance) {  
    Complex iEpsilon(0.0, EPSTOL * sqr(mMax));
    
    Complex xPlus, xMinus;

    /// alternative form: should be more accurate
    Complex oneiEpsilon(1.0, EPSTOL);
    Complex qq =
      -0.5 * (-s + sqrt(sqr(s) - 4.0 * pSq * sqr(mMax) * oneiEpsilon));
    xMinus  = qq / pSq;
    xPlus = sqr(mMax) * oneiEpsilon / qq;
    
    /*    xPlus = (s + sqrt(sqr(s) - 4. * sqr(p) * (sqr(mMax) - iEpsilon))) /
      (2. * sqr(p));
    xMinus = 2. * (sqr(mMax) - iEpsilon) / 
    (s + sqrt(sqr(s) - 4. * sqr(p) * (sqr(mMax) - iEpsilon)));*/
    //    cout << "DEBUG2: " << xPlus << " " << xMinus << endl;
    ans = -2.0 * log(p / q) - fB(xPlus) - fB(xMinus);
  } else {
    if (close(m1, m2, EPSTOL)) {
      //      	cout << "DEBUG5\n";
      ans = - log(sqr(m1 / q));
    } else {
      double Mmax2 = sqr(mMax), Mmin2 = sqr(mMin); 
      if (Mmin2 < 1.e-30) {
	//	cout << "DEBUG3\n";
	ans = 1.0 - log(Mmax2 / sqr(q));
      } else {
	//	cout << "DEBUG4\n";
	ans = 1.0 - log(Mmax2 / sqr(q)) + Mmin2 * log(Mmax2 / Mmin2) 
	  / (Mmin2 - Mmax2);
      }
    }
  }   

#ifdef USE_LOOPTOOLS
  if (!close(b0l, ans, 1.0e-3)) {
    cout << "DEBUG Err: DB0(" << p << ", " << m1 << ", " << m2 
	 << ", "  << q << ")=" << 1.-b0l/ans << endl;
    cout << "SOFTSUSY  B0=" << ans << endl;
    cout << "LOOPTOOLS B0=" << b0l << endl;
  }
#endif

  return ans;
}

/*
  Analytic expressions follow for above integrals: sometimes useful!
  From hep-ph/9606211
  Note it returns the REAL PART ONLY. 
*/
Complex b0c(double p, double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  if (fabs(m1) < 1.0e-30) m1 = 1.0e-30;
  if (fabs(m2) < 1.0e-30) m2 = 1.0e-30;  
  Complex b0l = B0(p*p, m1*m1, m2*m2);
  return B0(p*p, m1*m1, m2*m2);
#endif

  /// Avoids IR infinities
  if (close(p, 0.0, EPSTOL) && close(m1, 0.0, EPSTOL)
      && close(m2, 0.0, EPSTOL))
    return 0.0;

  Complex ans  = 0.;
  double mMin = minimum(fabs(m1), fabs(m2));
  double mMax = maximum(fabs(m1), fabs(m2));

  double pSq = sqr(p), mMinSq = sqr(mMin), mMaxSq = sqr(mMax);
  double s = 0.;
  /// Try to increase the accuracy of s
  double dmSq = mMaxSq - mMinSq;
  s = pSq + dmSq;

  double pTest = sqr(p) / sqr(mMax);
  /// Decides level at which one switches to p=0 limit of calculations
  const double pTolerance = 1.0e-6; 

  /// p is not 0  
  if (pTest > pTolerance) {  
    //    Complex iEpsilon(0.0, EPSTOL * sqr(mMax));
    Complex xPlus, xMinus;
    
    /// alternative form: should be more accurate
    Complex oneiEpsilon(1.0, EPSTOL);
    Complex qq =
      -0.5 * (-s + sqrt(sqr(s) - 4.0 * pSq * sqr(mMax) * oneiEpsilon));
    xMinus  = qq / pSq;
    xPlus   = sqr(mMax) * oneiEpsilon / qq;
     
    /*    xPlus = (s + sqrt(sqr(s) - 4. * sqr(p) * (sqr(mMax) - iEpsilon))) /
      (2. * sqr(p));
    xMinus = 2. * (sqr(mMax) - iEpsilon) / 
    (s + sqrt(sqr(s) - 4. * sqr(p) * (sqr(mMax) - iEpsilon)));*/

    /*cout << "DEBUG: " << pSq << " " << mMax << " " << s << endl;
    cout << "Comparison: " << xPlus << "=" << xPlus1 << " " << xMinus << "=" << xMinus1 << endl;
    exit(0);*/
    ans = -2.0 * log(p / q) - fBc(xPlus) - fBc(xMinus);
  } else {
    if (close(m1, m2, EPSTOL)) {
      ans = - log(sqr(m1 / q));
    } else {
      double Mmax2 = sqr(mMax), Mmin2 = sqr(mMin); 
      if (Mmin2 < 1.e-30) {
	ans = 1.0 - log(Mmax2 / sqr(q));
      } else {
	ans = 1.0 - log(Mmax2 / sqr(q)) + Mmin2 * log(Mmax2 / Mmin2) 
	  / (Mmin2 - Mmax2);
      }
    }
  }   

#ifdef USE_LOOPTOOLS
  if (!close(b0l.real(), ans.real(), 1.0e-3) || 
      !close(b0l.imag(), ans.imag(), 1.0e-3)) {
    if (b0l.imag() < EPSTOL && ans.imag() < 5.0e-8) return ans;
    else {
	cout << "DEBUG Err: DB0(" << p << ", " << m1 << ", " << m2 
	     << ", "  << q << ")=" << 1.-b0l/ans << endl;
	cout << "SOFTSUSY  B0=" << ans << endl;
	cout << "LOOPTOOLS B0=" << b0l << endl;
      }
  }
#endif

  return ans;
}


/// Note that b1 is NOT symmetric in m1 <-> m2!!!
double b1(double p, double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  if (fabs(m1) < 1.0e-30) m1 = 1.0e-30;
  if (fabs(m2) < 1.0e-30) m2 = 1.0e-30;  
  double b1l = -B1(p*p, m1*m1, m2*m2).real();
  return b1l;
#endif

  double ans = 0.;
  double pTest = sqr(p) / maximum(sqr(m1), sqr(m2));

  /// Decides level at which one switches to p=0 limit of calculations
  const double pTolerance = 1.0e-4; 

  if (pTest > pTolerance) {
    ans = (a0(m2, q) - a0(m1, q) + (sqr(p) + sqr(m1) - sqr(m2)) 
	   * b0(p, m1, m2, q)) / (2.0 * sqr(p));
  } else if (fabs(m1) > 1.0e-15 && fabs(m2) > 1.0e-15) { ///< checked
    const double m12 = sqr(m1), m22 = sqr(m2);
    const double m14 = sqr(m12), m24 = sqr(m22);
    const double m16 = m12*m14 , m26 = m22*m24;
    const double m18 = sqr(m14), m28 = sqr(m24);
    const double p2 = sqr(p), p4 = sqr(p2), q2 = sqr(q);
    if (fabs(m12 - m22) < pTolerance * maximum(m12, m22)) {
      ans = 0.08333333333333333*p2/m22
	+ 0.008333333333333333*p4/m24
	+ sqr(m12 - m22)*(0.041666666666666664/m24 +
			  0.016666666666666666*p2/m26 +
			  0.005357142857142856*p4/m28)
	+ (m12 - m22)*(-0.16666666666666666/m22 -
		       0.03333333333333333*p2/m24 -
		       0.007142857142857142*p4/m26)
	- 0.5*log(m22/q2);
    } else {
      ans = (3*m14 - 4*m12*m22 + m24 - 2*m14*log(m12/m22))/(4.*sqr(m12 - m22))
	+ (p2*(4*pow(m12 - m22,3)*
	       (2*m14 + 5*m12*m22 - m24) +
	       (3*m18 + 44*m16*m22 - 36*m14*m24 - 12*m12*m26 + m28)*p2
	       - 12*m14*m22*(2*sqr(m12 - m22) + (2*m12 + 3*m22)*p2)*
	       log(m12/m22)))/
	(24.*pow(m12 - m22,6)) - 0.5*log(m22/q2);
    }
  } else {
    ans = bIntegral(1, p, m1, m2, q); 
  }

#ifdef USE_LOOPTOOLS
  if (!close(b1l, ans, 1.0e-3)) {
    cout << " Test=" << pTest << " ";
    cout << "DEBUG Err: Db1(" << p << ", " << m1 << ", " << m2
    	 << ", "  << q << ")=" << 1.-b1l/ans << endl;
    cout << "SOFTSUSY  B1=" << ans << endl;
    cout << "LOOPTOOLS B1=" << b1l << endl;
  }
#endif
  
  return ans;
}

/// Note that b1 is NOT symmetric in m1 <-> m2!!!
Complex b1c(double p, double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  if (fabs(m1) < 1.0e-30) m1 = 1.0e-30;
  if (fabs(m2) < 1.0e-30) m2 = 1.0e-30;    
  Complex b1l = -B1(p*p, m1*m1, m2*m2);
  return b1l;
#endif

  Complex ans = 0.;
  double pTest = sqr(p) / maximum(sqr(m1), sqr(m2));

  /// Decides level at which one switches to p=0 limit of calculations
  const double pTolerance = 1.0e-4; 

  if (pTest > pTolerance) {
    ans = (a0(m2, q) - a0(m1, q) + (sqr(p) + sqr(m1) - sqr(m2)) 
	   * b0c(p, m1, m2, q)) / (2.0 * sqr(p)); 
  } else if (fabs(m1) > 1.0e-15 && !close(m1, m2, EPSTOL) 
	     && fabs(m2) > 1.0e-15) { ///< checked
    double Mmax2 = maximum(sqr(m1) , sqr(m2)), x = sqr(m2 / m1);
    ans = 0.5 * (-log(Mmax2 / sqr(q)) + 0.5 + 1.0 / (1.0 - x) + log(x) /
		 sqr(1.0 - x) - theta(1.0 - x) * log(x)); ///< checked
    ans = 0.5 * (1. + log(sqr(q) / sqr(m2)) + 
		 sqr(sqr(m1) / (sqr(m1) - sqr(m2))) * log(sqr(m2) / sqr(m1)) +
		 0.5 * (sqr(m1) + sqr(m2)) / (sqr(m1) - sqr(m2))
		 );
  } else {
    ans = bIntegral(1, p, m1, m2, q); 
  }
#
#ifdef USE_LOOPTOOLS
  if (!close(b1l, ans, 1.0e-3) && 
      (b1l.imag() > 1.0e-6 && ans.imag() > 1.0e-6)) {

    cout << "DEBUG Err: Db1(" << p << ", " << m1 << ", " << m2 
	 << ", "  << q << ")=" << (ans - b1l) / ans.mod() << endl;
    cout << "SOFTSUSY  B1=" << ans << endl;
    cout << "LOOPTOOLS B1=" << b1l << endl;
  }
#endif

return ans;
}

double b22(double p,  double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  if (fabs(m1) < 1.0e-30) m1 = 1.0e-30;
  if (fabs(m2) < 1.0e-30) m2 = 1.0e-30;  
  double b22l = B00(p*p, m1*m1, m2*m2).real();
  return b22l;
#endif

  double answer = 0.;
  
  /// Decides level at which one switches to p=0 limit of calculations
  const double pTolerance = 1.0e-6; 

  if (sqr(p) < pTolerance * maximum(sqr(m1), sqr(m2)) ) {
    // m1 == m2 with good accuracy
    if (close(m1, m2, EPSTOL)) {

      answer = -sqr(m1) * log(sqr(m1 / q)) * 0.5 + sqr(m1) * 0.5;
    }
    else
      /// This zero p limit is good
      if (fabs(m1) > EPSTOL && fabs(m2) > EPSTOL) {
	answer = 0.375 * (sqr(m1) + sqr(m2)) - 0.25 * 
	  (sqr(sqr(m2)) * log(sqr(m2 / q)) - sqr(sqr(m1)) * 
	   log(sqr(m1 / q))) / (sqr(m2) - sqr(m1)); 
      }
      else
	if (fabs(m1) < EPSTOL) {
	  answer = 0.375 * sqr(m2) - 0.25 * sqr(m2) * log(sqr(m2 / q));
	}
	else {
	  answer = 0.375 * sqr(m1) - 0.25 * sqr(m1) * log(sqr(m1 / q));
	}
  }
  else {// checked
    double b0Save = b0(p, m1, m2, q);

    answer = 1.0 / 6.0 * 
      (0.5 * (a0(m1, q) + a0(m2, q)) + (sqr(m1) + sqr(m2) - 0.5 * sqr(p))
       * b0Save + (sqr(m2) - sqr(m1)) / (2.0 * sqr(p)) *
       (a0(m2, q) - a0(m1, q) - (sqr(m2) - sqr(m1)) * b0Save) +
       sqr(m1) + sqr(m2) - sqr(p) / 3.0);
  }

#ifdef USE_LOOPTOOLS
  if (!close(b22l, answer, 1.0e-3)) {
    cout << " DEBUG Err: Db22(" << p << ", " << m1 << ", " << m2 
	 << ", "  << q << ")=" << 1.-b22l/answer << endl;
    cout << "SOFTSUSY  B22=" << answer << " B0=" << b0(p, m1, m2, q) << endl;
    cout << "LOOPTOOLS B22=" << b22l << " B0=" << B0(p*p, m1*m1, m2*m2).real() 
	 << endl;
  }
#endif

  return answer;
}

Complex b22c(double p,  double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  if (fabs(m1) < 1.0e-30) m1 = 1.0e-30;
  if (fabs(m2) < 1.0e-30) m2 = 1.0e-30;  
  Complex b22l = B00(p*p, m1*m1, m2*m2);
  return b22l;
#endif

  Complex answer = 0.;
  
  /// Decides level at which one switches to p=0 limit of calculations
  const double pTolerance = 1.0e-6; 

  if (sqr(p) < pTolerance * maximum(sqr(m1), sqr(m2)) ) {
    // m1 == m2 with good accuracy
    if (close(m1, m2, EPSTOL)) {

      answer = -sqr(m1) * log(sqr(m1 / q)) * 0.5 + sqr(m1) * 0.5;
    }
    else
      /// This zero p limit is good
      if (fabs(m1) > EPSTOL && fabs(m2) > EPSTOL) {
	answer = 0.375 * (sqr(m1) + sqr(m2)) - 0.25 * 
	  (sqr(sqr(m2)) * log(sqr(m2 / q)) - sqr(sqr(m1)) * 
	   log(sqr(m1 / q))) / (sqr(m2) - sqr(m1)); 
      }
      else
	if (fabs(m1) < EPSTOL) {
	  answer = 0.375 * sqr(m2) - 0.25 * sqr(m2) * log(sqr(m2 / q));
	}
	else {
	  answer = 0.375 * sqr(m1) - 0.25 * sqr(m1) * log(sqr(m1 / q));
	}
  }
  else {// checked
    Complex b0Save = b0c(p, m1, m2, q);

    answer = 1.0 / 6.0 * 
      (0.5 * (a0(m1, q) + a0(m2, q)) + (sqr(m1) + sqr(m2) - 0.5 * sqr(p))
       * b0Save + (sqr(m2) - sqr(m1)) / (2.0 * sqr(p)) *
       (a0(m2, q) - a0(m1, q) - (sqr(m2) - sqr(m1)) * b0Save) +
       sqr(m1) + sqr(m2) - sqr(p) / 3.0);
  }

#ifdef USE_LOOPTOOLS
  if (!close(b22l, answer, 1.0e-3)) {
    /// mitigating against when the imaginary part is supposed to be zero
    if (!(b22l.imag() / b22l.real() < 1.0e-9 && 
	  answer.imag() / answer.mod() < 1.0e-6)) {
      cout << " DEBUG Err: Db22(" << p << ", " << m1 << ", " << m2 
	   << ", "  << q << ")=" << 1.-b22l/answer << endl;
      cout << "SOFTSUSY  B22=" << answer << " B0=" << b0(p, m1, m2, q) << endl;
      cout << "LOOPTOOLS B22=" << b22l << " B0=" 
	   << B0(p*p, m1*m1, m2*m2).real() 
	   << endl;
    }
  }
#endif

  return answer;
}

// debugged 23.01.07 - thanks to Shindou Tetsuo 
double d0(double m1, double m2, double m3, double m4) {
  if (close(m1, m2, EPSTOL)) {
    double m2sq = sqr(m2), m3sq = sqr(m3), m4sq = sqr(m4);

    if (close(m2,0.,EPSTOL)) {
       // d0 is undefined for m1 == m2 == 0
       return 0.;
    } else if (close(m3,0.,EPSTOL)) {
       return (-sqr(m2) + sqr(m4) - sqr(m2) * log(sqr(m4/m2)))/
          sqr(m2 * sqr(m2) - m2 * sqr(m4));
    } else if (close(m4,0.,EPSTOL)) {
       return (-sqr(m2) + sqr(m3) - sqr(m2) * log(sqr(m3/m2)))/
          sqr(m2 * sqr(m2) - m2 * sqr(m3));
    } else if (close(m2, m3, EPSTOL) && close(m2, m4, EPSTOL)) {
      return 1.0 / (6.0 * sqr(m2sq));
    } else if (close(m2, m3, EPSTOL)) {
      return (sqr(m2sq) - sqr(m4sq) + 2.0 * m4sq * m2sq * log(m4sq / m2sq)) / 
	(2.0 * m2sq * sqr(m2sq - m4sq) * (m2sq - m4sq));
    } else if (close(m2, m4, EPSTOL)) {
      return (sqr(m2sq) - sqr(m3sq) + 2.0 * m3sq * m2sq * log(m3sq / m2sq)) / 
	(2.0 * m2sq * sqr(m2sq - m3sq) * (m2sq - m3sq));
    } else if (close(m3, m4, EPSTOL)) {
      return -1.0 / sqr(m2sq - m3sq) * 
	((m2sq + m3sq) / (m2sq - m3sq) * log(m3sq / m2sq) + 2.0);
    }

    return 
      (m4sq / sqr(m2sq - m4sq) * log(m4sq / m2sq) + 
       m4sq / (m2sq * (m2sq - m4sq)) -
       m3sq / sqr(m2sq - m3sq) * log(m3sq / m2sq) -
       m3sq / (m2sq * (m2sq - m3sq))) / (m3sq - m4sq);
  }
  return (c0(m1, m3, m4) - c0(m2, m3, m4)) / (sqr(m1) - sqr(m2));
}

double d27(double m1, double m2, double m3, double m4) {// checked

  if (close(m1, m2, EPSTOL)) {
    const double m1n = m1 + TOLERANCE * 0.01;
    return (sqr(m1n) * c0(m1n, m3, m4) - sqr(m2) * c0(m2, m3, m4)) 
      / (4.0 * (sqr(m1n) - sqr(m2)));
  }
  return (sqr(m1) * c0(m1, m3, m4) - sqr(m2) * c0(m2, m3, m4)) 
    / (4.0 * (sqr(m1) - sqr(m2)));
}

// Bug-fixed 14.10.02 by T. Watari and collaborators - many thanks!
double c0(double m1, double m2, double m3) {

#ifdef USE_LOOPTOOLS
  double q = 100.;
  setmudim(q*q); 
  double psq = 0.;
  double c0l = C0(psq, psq, psq, m1*m1, m2*m2, m3*m3).real();
#endif

  double ans = 0.;

  if (close(m1,0.,EPSTOL) && close(m2,0.,EPSTOL) && close(m3,0.,EPSTOL)) {
     // c0 is undefined for m1 == m2 == m3 == 0
     ans = 0.;
  } else if (close(m2,0.,EPSTOL) && close(m3,0.,EPSTOL)) {
     // c0 is undefined for m2 == m3 == 0
     ans = 0.;
  } else if (close(m1,0.,EPSTOL) && close(m3,0.,EPSTOL)) {
     // c0 is undefined for m1 == m3 == 0
     ans = 0.;
  } else if (close(m1,0.,EPSTOL) && close(m2,0.,EPSTOL)) {
     // c0 is undefined for m1 == m2 == 0
     ans= 0.;
  } else if (close(m1,0.,EPSTOL)) {
     if (close(m2,m3,EPSTOL)) {
        ans = -1./sqr(m2);
     } else {
        ans = (-log(sqr(m2)) + log(sqr(m3)))/(sqr(m2) - sqr(m3));
     }
  } else if (close(m2,0.,EPSTOL)) {
     if (close(m1,m3,EPSTOL)) {
        ans = -1./sqr(m1);
     } else {
        ans = log(sqr(m3/m1))/(sqr(m1) - sqr(m3));
     }
  } else if (close(m3,0.,EPSTOL)) {
     if (close(m1,m2,EPSTOL)) {
        ans = -1./sqr(m1);
     } else {
        ans = log(sqr(m2/m1))/(sqr(m1) - sqr(m2));
     }
  } else if (close(m2, m3, EPSTOL)) {
    if (close(m1, m2, EPSTOL)) {
      ans = ( - 0.5 / sqr(m2) ); // checked 14.10.02
    } else {
      ans = ( sqr(m1) / sqr(sqr(m1)-sqr(m2) ) * log(sqr(m2)/sqr(m1))
               + 1.0 / (sqr(m1) - sqr(m2)) ) ; // checked 14.10.02
    }
  } else if (close(m1, m2, EPSTOL)) {
     ans = ( - ( 1.0 + sqr(m3) / (sqr(m2)-sqr(m3)) * log(sqr(m3)/sqr(m2)) )
             / (sqr(m2)-sqr(m3)) ) ; // checked 14.10.02
  } else if (close(m1, m3, EPSTOL)) {
     ans = ( - (1.0 + sqr(m2) / (sqr(m3)-sqr(m2)) * log(sqr(m2)/sqr(m3)))
             / (sqr(m3)-sqr(m2)) ); // checked 14.10.02
  } else {
     ans = (1.0 / (sqr(m2) - sqr(m3)) *
            (sqr(m2) / (sqr(m1) - sqr(m2)) *
             log(sqr(m2) / sqr(m1)) -
             sqr(m3) / (sqr(m1) - sqr(m3)) *
             log(sqr(m3) / sqr(m1))) );
  }

#ifdef USE_LOOPTOOLS
  if (!close(c0l, ans, 1.0e-3)) {
    cout << " DEBUG Err: C0" << m1 << ", " << m2 
	 << ", "  << m3 << ")=" << 1.-c0l/ans << endl;
    cout << "SOFTSUSY  C0=" << ans << endl;
    cout << "LOOPTOOLS C0=" << c0l << endl;
  }
#endif

  return ans;
}

double truncGaussWidthHalf(long & idum) {
  for (;;) {
    double a = ran1(idum);
    if (a < 0.5) return a * 2.0;
	else return fabs(gasdev(idum)) + 1.;
  }
}

double gasdev(long & idum) {
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

double ran1(long & idum) {
  const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, NTAB = 32, 
    NDIV = 1+(IM-1)/NTAB;
  const double AM = 1.0 / double(IM), EPS = 1.2e-15, RNMX = 1.0 - EPS;
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (idum <= 0 || !iy) {
    if (-(idum) < 1) idum=1;
    else idum = -(idum);
    for (j = NTAB+7; j >= 0; j--) {
      k = (idum) / IQ;
      idum = IA * (idum - k * IQ) - IR * k;
      if (idum < 0) idum += IM;
      if (j < NTAB) iv[j] = idum;
    }
    iy=iv[0];
  }
  k = (idum) / IQ;
  idum = IA * (idum - k * IQ) - IR * k;
  if (idum < 0) idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = idum;
  if ((temp = AM * iy) > RNMX) return RNMX;
  else return temp;
}

double cauchyRan(long & idum) {
  double x = ran1(idum) - 0.5;
  double unNormalised = tan(x * PI);
  return unNormalised;
}

int bin(double data, double start, double end, int numBins) {
  double range = end - start;
  double binSize = range / double(numBins);
  return int((data - start) / binSize + 1.);
}

double logOfSum(double a, double b) {
  double max = maximum(a, b);
  double min = minimum(a, b);

  if (max + min < 0. || max == 0.) return asin(1.0);
  double ans = log(max);
  ans = ans + log (1.0 + min / max);

  return ans;
}

double sumOfExp(double a, double b) {
  double max = maximum(a, b);
  double min = minimum(a, b);
  
  /// this underflow determines whether it's actually worth adding the two
  /// figures... any more than 15 digits and it's just going to be the maximum
  /// one. 
  const double underflowTol = 15.0 * log(10.);
  if ((max - min) > underflowTol) return exp(max);
  else return exp(max) * (1.0 + exp(min - max));
}

double dilog(double x) {
   // The DiLogarithm function
   // Code translated by R.Brun from CERNLIB DILOG function C332

   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {0.42996693560813697, 0.40975987533077105,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};
   
   double T,H,Y,S,A,ALFA,B1,B2,B0;
   
   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i=19;i>=0;i--){
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

Complex dilog(const Complex & x) {
  double a = x.real(), b = x.imag();
  double ansreal = 0., ansimag = 0.;

  dilogwrap_(&a, &b, &ansreal, &ansimag);
  return Complex(ansreal, ansimag);
}

namespace {
  long double sqr(long double a) { return a * a; }
  std::complex<long double> sqr(const std::complex<long double>& a) { return a * a; }
} // anonymous namespace

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
std::complex<long double> dilogcl(const std::complex<long double>& z) {
   using std::log;
   const long double PI = 3.1415926535897932384626433832795l;
   std::complex<long double> cy, cz;
   int jsgn, ipi12;
   static const int N = 20;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 19}]
   const long double bf[N] = {
      - 1.l/4.l,
      + 1.l/36.l,
      - 1.l/36.e2l,
      + 1.l/21168.e1l,
      - 1.l/108864.e2l,
      + 1.l/52690176.e1l,
      - 4.0647616451442255268059093862919666745470571274397078e-11l,
      + 8.9216910204564525552179873167527488515142836130490451e-13l,
      - 1.9939295860721075687236443477937897056306947496538801e-14l,
      + 4.5189800296199181916504765528555932283968190144666184e-16l,
      - 1.0356517612181247014483411542218656665960912381686505e-17l,
      + 2.3952186210261867457402837430009803816789490019429743e-19l,
      - 5.5817858743250093362830745056254199055670546676443981e-21l,
      + 1.3091507554183212858123073991865923017498498387833038e-22l,
      - 3.0874198024267402932422797648664624315955652561327457e-24l,
      + 7.315975652702203420357905609252148591033401063690875e-26l,
      - 1.7408456572340007409890551477597025453408414217542713e-27l,
      + 4.1576356446138997196178996207752266734882541595115639e-29l,
      - 9.9621484882846221031940067024558388498548600173944888e-31l,
      + 2.3940344248961653005211679878937495629342791569329158e-32l,
   };

   const long double rz = std::real(z);
   const long double iz = std::imag(z);
   const long double az = std::sqrt(sqr(rz) + sqr(iz));

   // special cases
   if (iz == 0.l) {
      if (rz <= 1.l)
         return std::complex<long double>(dilog(rz), 0.l);
      if (rz > 1.l)
         return std::complex<long double>(dilog(rz), -PI*std::log(rz));
   } else if (az < std::numeric_limits<long double>::epsilon()) {
      return z;
   }

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5l) {
      if (az > 1.l) {
         cy = -0.5l * sqr(log(-z));
         cz = -log(1.l - 1.l / z);
         jsgn = -1;
         ipi12 = -2;
      } else { // (az <= 1.)
         cy = 0;
         cz = -log(1.l - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5
      if (az <= std::sqrt(2*rz)) {
         cz = -log(z);
         cy = cz * log(1.l - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // (az > sqrt(2*rz))
         cy = -0.5l * sqr(log(-z));
         cz = -log(1.l - 1.l / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   // the dilogarithm
   const std::complex<long double> cz2(sqr(cz));
   std::complex<long double> sumC;

   for (int i1 = 2; i1 < N; i1++)
      sumC = cz2 * (sumC + bf[N + 1 - i1]);

   // lowest order terms w/ different powers
   sumC = cz + cz2 * (bf[0] + cz * (bf[1] + sumC));

   const std::complex<long double> result
      = static_cast<long double>(jsgn) * sumC + cy + ipi12 * PI * PI / 12.l;

   return result;
}

DoubleVector getRandomDirection(int n, int & numChanged,
				long & idum) {
  DoubleVector del(n);

  int i; 

  double length = 0.;
  do {
    numChanged = 0;
    for (i=1; i<=n; i++) {
	del(i) = gasdev(idum);
	numChanged++;
    }
    length = sqrt(del.dot(del));
  }
  while (length < EPSTOL);

  del = 1.0 / sqrt(del.dot(del)) * del;
  return del;
}


double calcCL(double cl, const DoubleVector & l) {
  DoubleVector s(l.sort());

  double fr = 0., likelihood = 0., tot = 0.;
  int i; 
  for (i=s.displayStart(); i<=s.displayEnd(); i++) tot = tot + s.display(i);
  for (i=s.displayStart(); i<=s.displayEnd(); i++) {
    fr = fr + s.display(i) / tot;
    if (fr > (1.0 - cl)) {
      likelihood = s.display(i);
      // bail out of loop
      i = s.displayEnd() + 1;
    }
  }
  return likelihood;
}

int *ivector(long nl, long nh) {
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
	if (!v) throw("allocation failure in ivector()\n");
	return v-nl+1;
}

/* free an int vector allocated with ivector() */
void free_ivector(int *v, long nl, long /* nh */) {
	free((char *) (v+nl-1));
}

double calc1dFraction(double y, const DoubleVector & l) {
  double n = 0., lMax = l.max();
  double tot = 0.;
  int i; for (i=1; i<=l.displayEnd(); i++) {
    if (l.display(i) <= y * lMax) n = n + l.display(i);
    tot = tot + l.display(i);
  }

  return n / tot;
}

double fps(double z) {
  if (z < 0.25) {
    double y = sqrt(1.0 - 4.0 * z);
    return 2.0 * z / y * 
      (dilog(1.0 - (1.0 - y) / (2.0 * z)) - 
       dilog(1.0 - (1.0 + y) / (2.0 * z)));
  }
  Complex zz(z);
  Complex y = sqrt(1.0 - 4.0 * zz);

  Complex ans = 2.0 * zz / y * 
    (dilog(1.0 - (1.0 - y) / (2.0 * zz)) - 
     dilog(1.0 - (1.0 + y) / (2.0 * zz)));

  /// answer should always be real
  if (ans.imag() > EPSTOL) throw("Error in fps\n");
  return ans.real();
}

double fs(double z) {
  return (2.0 * z - 1) * fps(z) - 2.0 * z * (2.0 + log(z));
}

double ffbar(double z) {
  return z * 0.5 * (2.0 + log(z) - fps(z));
}

#define JMAX 20

#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double), double a, double b, int n, double /* EPS */) {
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

double qtrap(double (*func)(double), double a, double b, double EPS) {
  //	double midpnt(double (*func)(double), double a, double b, int n);
  //	void nrerror(char error_text[]);

	int j;
	double s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=midpnt(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		if (s == 0.0 && olds == 0.0 && j > 6) return s;
		olds=s;
	}
	throw("Too many steps in routine qtrap\n");
	return 0.0;
}
#undef JMAX

double midpnt(double (*func)(double), double a, double b, int n) {
	double x,tnm,sum,del,ddel;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}

double edgefn(double x, double y, double z) {
  return (x - y) * (y - z) / y;
}
double mllMax(double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double expression = edgefn(xi, l, chi);
  if (mChi2 < mSlep) expression = sqr(mChi2-mChi1);
  if (expression < 0.) return -numberOfTheBeast; 
  return sqrt(expression);
}

double mllq(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double q     = sqr(mSq);
  if (sqr(l) < q * chi && q * chi < sqr(xi) && sqr(xi) * chi < q * sqr(l))
    return (mSq - mChi1);
  double exp1  = edgefn(q, xi, chi);
  if (exp1 < 0.) exp1 = -numberOfTheBeast; 
  double exp2  = edgefn(q, l, chi);
  if (exp2 < 0.) exp2 = -numberOfTheBeast; 
  double exp3  = (q * l - xi * chi) * (xi - l) / (xi * l);
  if (exp3 < 0.) exp3 = -numberOfTheBeast; 
  double expression = maximum(maximum(exp1, exp2), exp3);
  if (mChi2 < mSlep  && xi > sqr(q) * sqr(chi) ) return (mSq - mChi1);
  if (mChi2 < mSlep  && xi < sqr(q) * sqr(chi) )  expression = edgefn(q,xi,chi);  
  return sqrt(expression);
}

double llqThresh(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double q     = sqr(mSq);
  double ans   = 2.0 * l * (q - xi) * (xi - chi) + 
    (q + xi) * (xi - l) * (l - chi) - (q - xi) * 
    sqrt(sqr(xi + l) * sqr(l + chi) - 16. * xi * sqr(l) * chi); 
  if (ans < 0.) return -numberOfTheBeast;
  return sqrt(ans / (4.0 * l * xi));
}

double lqnear(double mSq, double /* mChi1 */, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double q     = sqr(mSq);
  double expression = (q - xi) * (xi - l) / xi;
  if (expression < 0.) return -numberOfTheBeast;
  return sqrt(expression);
}

double lqfar(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double q     = sqr(mSq);
  double expression = (q - xi) * (l - chi) / l;
  if (expression < 0.) return -numberOfTheBeast;
  return sqrt(expression);
}

double lqhigh(double mSq, double mChi1, double mSlep, double mChi2) {
  double mlqnear = lqnear(mSq, mChi1, mSlep, mChi2);
  double mlqfar  = lqfar(mSq, mChi1, mSlep, mChi2);
if (mChi2 < mSlep) {return sqrt(edgefn(sqr(mSq),sqr(mChi2),sqr(mChi1)));}
  return maximum(mlqnear, mlqfar);
}

double lqlow(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double q     = sqr(mSq);
  double mlqnear = lqnear(mSq, mChi1, mSlep, mChi2);
  double mlqother  = (q - xi) * (l - chi) / (2.0 * l - chi);
  if (mlqother < 0.) mlqother = -numberOfTheBeast;
  else mlqother = sqrt(mlqother);

  return minimum(mlqnear, mlqother);
}



#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

double qromb(double (*func)(double), double a, double b, double EPS) {
  double ss = 0.,dss = 0.;
  DoubleVector s(JMAXP+1),h(JMAXP+1);
  int j;
  
  h(1) = 1.0;
  for (j=1; j<=JMAX; j++) {
    s(j)=trapzd(func,a,b,j);
    if (j >= K) {
      DoubleVector h2(K), s2(K);
      for (unsigned int ii=1; ii<=K; ii++) {
	h2(ii) = h(j-K+ii); 
	s2(ii) = s(j-K+ii);
      }	
      polint(h2,s2,0.0, ss, dss);
      if (fabs(dss) <= EPS*fabs(ss)) return ss;
    }
    s(j+1)=s(j);
    h(j+1)=0.25 * h(j);
  }
  throw("Too many steps in routine qromb\n");
  return 0.0;
}

double qromb2(double (*func)(double), double a, double b, double EPS) {
  double ss = 0.,dss = 0.;
  DoubleVector s(JMAXP+1),h(JMAXP+1);
  int j;
  
  h(1) = 1.0;
  for (j=1; j<=JMAX; j++) {
    s(j)=trapzd(func,a,b,j);
    if (j >= K) {
      DoubleVector h2(K), s2(K);
      for (unsigned int ii=1; ii<=K; ii++) {
	h2(ii) = h(j-K+ii); 
	s2(ii) = s(j-K+ii);
      }	
      polint(h2,s2,0.0, ss, dss);
      if (fabs(dss) <= EPS*fabs(ss)) return ss;
    }
    s(j+1)=s(j);
    h(j+1)=0.25 * h(j);
  }
  throw("Too many steps in routine qromb\n");
  return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K

void polint(const DoubleVector&  xa, const DoubleVector & ya, double x, 
	    double & y, double & dy) {
  if (xa.displayEnd() != ya.displayEnd()) throw("Incompatible vector lengths in numerics.cpp:polint\n");
  
  int n = xa.displayEnd();
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  DoubleVector c(n), d(n);

  dif=fabs(x-xa.display(1));
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa.display(i))) < dif) {
      ns=i;
      dif=dift;
    }
    c(i)=ya.display(i);
    d(i)=ya.display(i);
  }
  y=ya.display(ns--);
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa.display(i)-x;
      hp=xa.display(i+m)-x;
      w=c(i+1)-d(i);
      if ( (den=ho-hp) == 0.0) throw("Error in routine polint\n");
      den=w/den;
      d(i)=hp*den;
      c(i)=ho*den;
    }
    y += (dy=(2*ns < (n-m) ? c(ns+1) : d(ns--)));
  }
}

DoubleMatrix display3x3RealMixing(double theta12, double theta13, 
				  double theta23, double d) {
  Complex eID(cos(d), sin(d));
  double s12 = sin(theta12);
  double s13 = sin(theta13);
  double s23 = sin(theta23);
  
  double c12 = cos(theta12);
  double c13 = cos(theta13);
  double c23 = cos(theta23);

  DoubleMatrix ckmMatrix(3, 3);
  ckmMatrix(1, 1) = c12 * c13;      
  ckmMatrix(1, 2) = s12 * c13; 

  /// phase factor e^i delta: we'll set it to + or - 1 depending on the sign
  /// of s13
  int pf = 1;
  if (s13 < 0.) pf = -1;
  
  ckmMatrix(1, 3) = pf * s13;
  ckmMatrix(2, 1) = (-s12 * c23 - pf * c12 * s23 * s13);
  ckmMatrix(2, 2) = (c12 * c23 - pf * s12 * s23 * s13);
  ckmMatrix(2, 3) = s23 * c13; 
  ckmMatrix(3, 1) = (s12 * s23 - pf * c12 * c23 * s13); 
  ckmMatrix(3, 2) = (-c12 * s23 - pf * s12 * c23 * s13); 
  ckmMatrix(3, 3) = c23 * c13;

  return ckmMatrix;
  
}

void getAngles(const DoubleMatrix & v, double & t12, double & t13, 
	       double & t23, double & d) {
  int pf = 1;
  t23 = atan(v.display(2, 3) / v.display(3, 3));
  t13 = asin(v.display(1, 3));
  if (t13 < 0.) {
    pf = -1; 
    t13 = -t13;
  }
  d = acos(double(pf));
  t12 = asin(v.display(1, 2) / cos(t13));
}

bool midPtStep(DoubleVector & xi, 
	      DoubleVector (*derivs)(double t, const DoubleVector & v), 
	      double tInitial, double tStep) {
  /// initial guess
  DoubleVector xiPlus1Old(xi + tStep * derivs(tInitial + tStep, xi));
  DoubleVector xiPlus1New(xi.displayEnd());

  double diff = numberOfTheBeast, delta = 0.;
  int count = 0;
  while (diff > 1.0e-15 && count < 100) {
    count++;
    xiPlus1New = xi + tStep * derivs(tInitial + tStep, 0.5*(xi + xiPlus1Old));
    double max = -1.; /// difference between iterations
      int i; for (i=1; i<=xi.displayEnd(); i++) {
	if (fabs(xiPlus1New(i)) < 1.0e-12) 
	delta = fabs(xiPlus1New(i) - xiPlus1Old(i));
      else delta = fabs(1.- xiPlus1Old(i) / xiPlus1New(i));
      if (delta > max) max = delta;
      }

    diff = max;

    xiPlus1Old = xiPlus1New;
  }
  xi = xiPlus1New;
  if (count < 100) return false;
  return true;
}

/// Do a fixed number of steps of approximately reversible integration: hoping
/// it will help convergence in difficult cases. It is really *slow* though.
bool integrateReversibly(DoubleVector & xi, 
			 DoubleVector (*derivs)(double t, 
						const DoubleVector & v), 
			 double tInitial, double tFinal, int numSteps) {
  bool err = false;
  double tStep = (tFinal - tInitial) / double(numSteps);

  double t = tInitial;
  int i; for (i=0; i<=numSteps; i++) {

    if (midPtStep(xi, derivs, t, tStep)) err = true;
    t = tStep * double(i) + tInitial;
  }
  return err;
}

double den(double a, int b) {
  double aa = a;
  int i; for (i=1; i<b; i++) aa = aa * a;
  return 1. / aa;
}

double fin(double mm1, double mm2) {
  if (mm1>mm2)
    return (-3.5 - (7.*mm2)/(2.*mm1) +
	    dilog(mm2/mm1) - (mm2*dilog(mm2/mm1))/mm1 -
	    (3.*mm2*log(mm1))/mm1 - log(mm1)*log(mm1 - mm2) +
	    (mm2*log(mm1)*log(mm1 - mm2))/mm1 + (3.*mm2*log(mm2))/mm1 -
	    log(mm1)*log(mm2) + (2.*mm2*log(mm1)*log(mm2))/mm1 +
	    log(mm1 - mm2)*log(mm2) - (mm2*log(mm1 - mm2)*log(mm2))/mm1 - 
	    sqr(PI) * 0.25 +
	    (mm2*sqr(PI))/(12.*mm1) + sqr(log(mm1)) - 
	    (3.*mm2*sqr(log(mm1)))/(2.*mm1) -
	    (mm2*sqr(log(mm2)))/(2.*mm1));
  else if (mm1<mm2)
    return (-3.5 - (7.*mm2)/(2.*mm1) -
	    dilog(mm1/mm2) + (mm2*dilog(mm1/mm2))/mm1 -
	    (3*mm2*log(mm1))/mm1 + (3*mm2*log(mm2))/mm1 +
	    (mm2*log(mm1)*log(mm2))/mm1 - log(mm1)*
	    log(-mm1 + mm2)+ (mm2*log(mm1)*log(-mm1 + mm2))/mm1 
	    + log(mm2)*log(-mm1 + mm2) -
	    (mm2*log(mm2)*log(-mm1 + mm2))/mm1 + sqr(PI) / 12. -
	    (mm2*sqr(PI))/(4.*mm1) + sqr(log(mm1)) * 0.5 - 
	    (mm2*sqr(log(mm1)))/mm1 -
	    sqr(log(mm2))/2);
  else return 7.-sqr(PI)/6.;  
}

double zriddr(double (*func)(double), double x1, double x2, double xacc) {
  const int MAXIT = 60;
  const double UNUSED = -1.11e30;
  int j;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
  
  fl=(*func)(x1);
  fh=(*func)(x2);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl=x1;
    xh=x2;
    ans=UNUSED;
    for (j=1;j<=MAXIT;j++) {
      xm=0.5*(xl+xh);
      fm=(*func)(xm);
      s=sqrt(fm*fm-fl*fh);
      if (s == 0.0) return ans;
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xnew-ans) <= xacc) return ans;
      ans=xnew;
      fnew=(*func)(ans);
      if (fnew == 0.0) return ans;
      if (SIGN(fm,fnew) != fm) {
	xl=xm;
	fl=fm;
	xh=ans;
	fh=fnew;
      } else if (SIGN(fl,fnew) != fl) {
	xh=ans;
	fh=fnew;
      } else if (SIGN(fh,fnew) != fh) {
	xl=ans;
	fl=fnew;
      } else throw("In zriddr: never get here.");
      if (fabs(xh-xl) <= xacc) return ans;
    }
    throw("zriddr exceeded maximum iterations\n");
  }
  else {
    if (close(fl, 0.0, EPSTOL)) return x1;
    if (close(fh, 0.0, EPSTOL)) return x2;
    ostringstream s;
    s << "Root must be bracketed in zriddr. f(" << x1 << ")=" << fl 
      << ", fh(" << x2 << ")=" << fh << endl;
    throw(s.str());
  }
  return 0.0;
}

/// You will need to clear this lot up....
DoubleMatrix fdjac(int n, DoubleVector x, const DoubleVector & fvec,
           int (*vecfunc)(const DoubleVector &, void*, DoubleVector &),
           void* params) {
  double EPS = maximum(TOLERANCE, 1.0e-4);
  int i,j;
  double h,temp;
  
  DoubleVector f(1, n);
  DoubleMatrix df(n, n);
  for (j=1; j<=n; j++) {
    temp = x(j);
    h = EPS * fabs(temp);
    if (h == 0.0) h = EPS;
    x(j) = temp + h;
    h = x(j) - temp;
    (*vecfunc)(x, params, f);
    x(j) = temp;
    for (i=1; i<=n; i++) df(i, j) = (f(i) - fvec.display(i)) / h;
  }
  return df;
}

bool lnsrch(const DoubleVector & xold, double fold, const DoubleVector & g, 
	    DoubleVector & p, 
	    DoubleVector & x, double & f, double stpmax, 
	    int (*vecfunc)(const DoubleVector &, void*, DoubleVector &),
	    DoubleVector & fvec, void* params) {
  double ALF = TOLERANCE;
  double TOLX = TOLERANCE * 1.0e-3;
  
  int i;
  double a,alam,alam2 = 0.,alamin,b,disc,f2 = 0.,fold2 = 0.,rhs1,rhs2,slope,sum,temp,
    test,tmplam;
  
  bool err = false;
  sum = sqrt(p.dot(p));
  if (sum > stpmax) p = (stpmax / sum) * p;

  for (slope=0.0, i=1; i<=xold.displayEnd(); i++)
    slope += g(i) * p(i);
  test = 0.0;
  for (i=1; i<=xold.displayEnd(); i++) {
    temp=fabs(p(i)) / maximum(fabs(xold(i)), 1.0);
    if (temp > test) test = temp;
  }
  alamin = TOLX / test;
  alam = 1.0;
  for (;;) {
    x = xold + alam * p;
    vecfunc(x, params, fvec);
    f = fvec.dot(fvec);
    if (alam < alamin) {
      x = xold;
      err = true;
      return err;
    } else if (f <= fold + ALF * alam * slope) return err;
    else {
      if (alam == 1.0)
	tmplam = -slope / (2.0 * (f - fold - slope));
      else {
	rhs1 = f - fold - alam * slope;
	rhs2 = f2 - fold2 - alam2 * slope;
	a=(rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
	b=(-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / 
	  (alam - alam2);
	if (a == 0.0) tmplam = -slope / (2.0 * b);
	else {
	  disc = b * b - 3.0 * a * slope;
	  if (disc < 0.0) throw("Roundoff problem in lnsrch.\n");
	  else tmplam = (-b + sqrt(disc)) / (3.0 * a);
	}
	if (tmplam > 0.5 * alam)
	  tmplam = 0.5 * alam;
      }
    }
    alam2 = alam;
    f2 = f;
    fold2 = fold;
    alam = maximum(tmplam, 0.1 * alam);
  }
  return err;
}

/// Get rid of int n
void lubksb(const DoubleMatrix & a, int n, int *indx, DoubleVector & b) {
  int i,ii=0,ip,j;
  float sum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b(ip);
    b(ip)=b(i);
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a(i, j) * b(j);
    else if (sum) ii=i;
    b(i)=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b(i);
    for (j=i+1;j<=n;j++) sum -= a(i, j)*b(j);
    b(i)=sum/a(i, i);
  }
}


/// Get rid of int n
void ludcmp(DoubleMatrix & a, int n, int *indx, double & d) {
  const double TINY = 1.0e-20;
  int i, imax = 0, j, k;
  float big,dum,sum,temp;
  DoubleVector vv(n);

  d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a(i, j))) > big) big=temp;
    if (big == 0.0) {
      ostringstream ii; ii.setf(ios::scientific, ios::floatfield);
      ii.precision(6);
      ii << "Singular matrix in routine ludcmp: ";
      ii << a;
      throw (ii.str());
    }
    vv(i)=1.0 / big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a(i, j);
      for (k=1;k<i;k++) sum -= a(i, k)*a(k, j);
      a(i, j)=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a(i, j);
      for (k=1;k<j;k++)
	sum -= a(i, k)*a(k, j);
      a(i, j)=sum;
      if ( (dum=vv(i)*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a(imax, k);
	a(imax, k)=a(j, k);
	a(j, k)=dum;
      }
      d = -(d);
      vv(imax)=vv(j);
    }
    indx[j]=imax;
    if (a(j, j) == 0.0) a(j, j) = TINY;
    if (j != n) {
      dum=1.0/(a(j, j));
      for (i=j+1;i<=n;i++) a(i, j) *= dum;
    }
  }
}


/// More work can be done on this: get rid of int n and in subfunctions too
bool newt(DoubleVector & x, 
	  int (*vecfunc)(const DoubleVector &, void*, DoubleVector &),
          void* params) {
  bool err = false; 
  const int MAXITS = 200;    ///< max iterations
  double TOLF   = TOLERANCE; ///< convergence on function values
  double TOLMIN = TOLF * 1.e-2; ///< spurious convergence to min of fmin
  double TOLX   = TOLF * 1.e-3; ///< maximum dx convergence criterion
  const double STPMX  = 100.0; ///< maximum step length allowed in line searches
   
  int n = x.displayEnd();
  int i,its,j,*indx;
  double d,den,f,fold,stpmax,sum = 0., temp, test;
  
  indx = ivector(1, n);
  DoubleMatrix fjac(n, n);
  DoubleVector g(n), p(n), xold(n);
  DoubleVector fvec(n);

  vecfunc(x, params, fvec);
  f = 0.5 * fvec.dot(fvec);
  test = 0.0;
  test = fvec.apply(fabs).max();
  if (test < 0.01 * TOLF) {
    err = false;
    free_ivector(indx, 1, n); return err;
  }
  sum += x.dot(x);
  stpmax = STPMX * maximum(sqrt(sum), (double) n);
  for (its=1; its<=MAXITS; its++) {
    //    cout << its << endl; ///< DEBUG
    fjac = fdjac(n, x, fvec, vecfunc, params);
    for (i=1;i<=n;i++) {
      for (sum=0.0, j=1; j<=n; j++) sum += fjac(j, i) * fvec(j);
      g(i) = sum;
    }
    for (i=1; i<=n; i++) xold(i) = x(i);
    fold = f;
    for (i=1; i<=n; i++) p(i) = -fvec(i);
    ludcmp(fjac, n, indx, d);
    lubksb(fjac, n, indx, p);
    err = lnsrch(xold, fold, g, p, x, f, stpmax, vecfunc, fvec, params);
    test = 0.0;
    for (i=1; i<=n; i++)
      if (fabs(fvec(i)) > test) test = fabs(fvec(i));
    if (test < TOLF) {
      err = false;
      free_ivector(indx, 1, n);
      return err;
    }
    if (err) {
      test = 0.0;
      den = maximum(f, 0.5 * n);
      for (i=1; i<=n; i++) {
	temp = fabs(g(i)) * maximum(fabs(x(i)), 1.0) / den;
	if (temp > test) test = temp;
      }
      err = (test < TOLMIN ? true : false);
      free_ivector(indx, 1, n);
      return err;
    }
    /// Check points aren't getting too close together
    test = 0.0;
    for (i=1; i<=n; i++) {
      temp = (fabs(x(i) - xold(i))) / maximum(fabs(x(i)), 1.0);
      if (temp > test) test = temp;
    }
    if (test < TOLX) { 
      free_ivector(indx, 1, n); 
      return err; 
    }
  }

  return true;
}

DoubleVector testDerivs(double /* x */, const DoubleVector & y) {
  DoubleVector dydx(3);
  dydx(1) = y(1) * y(2) * y(2);
  dydx(2) = y(2) * y(1) * y(3);
  dydx(3) = y(1);
  return dydx;
}

/*void shoot(const DoubleVector & v, DoubleVector & f) {
  double h1, hmin = 0.0;

  const double EPS = 1.0e-6;
  
  DoubleVector y(3);
  
  double x1 = 1., x2 = 2.;
  /// Initial stepsize guess for integration
  h1 = (x2 - x1) * 0.01;

  /// set initial BCs: y1(1)=1
  y(1) = 1.;
  y(2) = v.display(1);
  y(3) = v.display(2);

  /// integrate up from x1 to x2
  int err = integrateOdes(y, x1, x2, EPS, h1, hmin, testDerivs, odeStepper);

  if (err) 
  
  /// now, determine a vector showing how far (WITH SIGN) the solution is from
  /// the second boundary condition: y2(2)=1.
  f(1) = y(2) - 1.;
  f(2) = y(3) - y(1);

  return;
  }*/
  
void broydn(DoubleVector x, int & check, 
	    int (*vecfunc)(const DoubleVector &, void*, DoubleVector &),
            void* params) {
  const int MAXITS = 200;
  double TOLF =  TOLERANCE;
  const double EPS = TOLF * 1.0e-3;
  const double TOLX = EPS;
  const double STPMX = 100.0;
  double TOLMIN = TOLERANCE * 1.0e-2;

  int n = x.displayEnd();
 
  /*	void fdjac(int n, float x[], float fvec[], float **df,
		void (*vecfunc)(int, float [], float []));
	float fmin(float x[]);
	void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
		 float *f, float stpmax, int *check, float (*func)(float []));
	void qrdcmp(float **a, int n, float *c, float *d, int *sing);
	void qrupdt(float **r, float **qt, int n, float u[], float v[]);
	void rsolv(float **a, int n, float d[], float b[]);
	int i,its,j,k,restrt,sing,skip;
	float den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
	float *g,*p,**qt,**r,*s,*t,*w,*xold; */

  DoubleVector c(n), d(n), fvcold(n), g(n), p(n), s(n), t(n), w(n), xold(n), 
    fvec(n);
  DoubleMatrix qt(n, n), r(n, n);

  vecfunc(x, params, fvec);
  double f = 0.5 * fvec.dot(fvec);
  
  double test = fvec.apply(fabs).max();
  if (test < 0.01 * TOLF) {
    check = 0;
    return;
  }
  double sum = x.dot(x);
  double stpmax = STPMX * maximum(sqrt(sum), (double)n);
  int restrt = 1, sing;
  for (int its=1; its<=MAXITS; its++) {
    if (restrt) {
       r = fdjac(n, x, fvec, vecfunc, params);
      qrdcmp(r, n, c, d, sing);
      if (sing) throw("singular Jacobian in broydn\n");
      for (int i=1; i<=n; i++) {
	for (int j=1; j<=n; j++) qt(i, j) = 0.0;
	qt(i, i) = 1.0;
      }
      for (int k=1; k<n; k++) {
	if (c(k)) {
	  for (int j=1; j<=n; j++) {
	    sum = 0.0;
	    for (int i=k; i<=n; i++)
	      sum += r(i, k) * qt(i, j);
	    sum /= c(k);
	    for (int i=k; i<=n; i++)
	      qt(i, j) -= sum*r(i, k);
	  }
	}
      }
      for (int i=1; i<=n; i++) {
	r(i, i) = d(i);
	for (int j=1; j<i; j++) r(i, j) = 0.0;
      }
    } else {
      s = x - xold;
      int j;
      for (int i=1; i<=n; i++) {
	for (sum=0.0, j=i; j<=n; j++) sum += r(i, j) * s(j);
	t(i) = sum;
      }
      int skip = 1;
      for (int i=1; i<=n; i++) {
	for (sum=0.0, j=1; j<=n; j++) sum += qt(j, i) * t(j);
	w(i) = fvec(i) - fvcold(i) - sum;
	if (fabs(w(i)) >= EPS * (fabs(fvec(i))+fabs(fvcold(i)))) skip=0;
	else w(i)=0.0;
      }
      if (!skip) {
	for (int i=1; i<=n; i++) {
	  for (sum=0.0, j=1; j<=n; j++) sum += qt(i, j) * w(j);
	  t(i) = sum;
	}
	double den = 0.; den += s.dot(s);
	s = s * (1. / den);
	qrupdt(r,qt,n,t,s);
	for (int i=1;i<=n;i++) {
	  if (r(i, i) == 0.0) throw("r singular in broydn\n");
	  d(i) = r(i, i);
	}
      }
    }
    int j;
    for (int i=1; i<=n; i++) {
      for (sum=0., j=1 ; j<=n; j++) sum += qt(i, j) * fvec(j);
      g(i) = sum;
    }
    for (int i=n; i>=1; i--) {
      for (sum=0.0, j=1; j<=i; j++) sum += r(j, i) * g(j);
      g(i) = sum;
    }
    xold = x; 
    fvcold = fvec;
    double fold = f;
    for (int i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += qt(i, j) * fvec(j);
      p(i) = -sum;
    }
    rsolv(r, n, d, p); 
    DoubleVector fvec(n);
    check = lnsrch(xold, fold, g, p, x, f, stpmax, vecfunc, fvec, params);
    test = 0.0;
    for (int i=1; i<=n; i++)
      if (fabs(fvec(i)) > test) test = fabs(fvec(i));
    if (test < TOLF) {
      check=0;
      return;
    }
    if (check) {
      if (restrt) return;
      else {
	test = 0.0;
	double den = maximum(f, 0.5 * n);
	for (int i=1;i<=n;i++) {
	  double temp = fabs(g(i)) * maximum(fabs(x(i)), 1.0) / den;
	  if (temp > test) test = temp;
	}
	if (test < TOLMIN) return;
	else restrt = 1;
      }
    } else {
      restrt = 0;
      test = 0.0;
      for (int i=1; i<=n; i++) {
	double temp = (fabs(x(i) - xold(i))) / maximum(fabs(x(i)), 1.0);
	if (temp > test) test = temp;
      }
      if (test < TOLX) return;
    }
  }
  throw("MAXITS exceeded in broydn\n");
  return;
}

void qrdcmp(DoubleMatrix & a, int n, DoubleVector & c, DoubleVector & d, 
	    int & sing) {
  int i,j,k;
  double scale,sigma,sum,tau;
  
  sing = 0;
  for (k=1; k<n; k++) {
    scale = 0.0;
    for (i=k; i<=n; i++) scale = maximum(scale, fabs(a(i, k)));
    if (scale == 0.0) {
      sing = 1;
      c(k) = d(k) = 0.0;
    } else {
      for (i=k;i<=n;i++) a(i, k) /= scale;
      for (sum=0.0,i=k;i<=n;i++) sum += sqr(a(i, k));
      sigma = sign(sqrt(sum), a(k, k));
      a(k, k) += sigma;
      c(k) = sigma * a(k, k);
      d(k) = -scale*sigma;
      for (j=k+1; j<=n; j++) {
	for (sum=0.0,i=k;i<=n;i++) sum += a(i, k) * a(i, j);
	tau = sum / c(k);
	for (i=k; i<=n; i++) a(i, j) -= tau * a(i, k);
      }
    }
  }
  d(n) = a(n, n);
  if (d(n) == 0.0) sing = 1;
}


void qrupdt(DoubleMatrix & r, DoubleMatrix & qt, int n, 
	    DoubleVector & u, DoubleVector & v) {
  //	void rotate(float **r, float **qt, int n, int i, float a, float b);
	int i,j,k;

	for (k=n; k>=1; k--) {
	  if (u(k)) break;
	}
	if (k < 1) k = 1;
	for (i=k-1; i>=1; i--) {
	  rotate(r, qt, n, i, u(i), -u(i+1));
	  if (u(i) == 0.0) u(i) = fabs(u(i+1));
		else if (fabs(u(i)) > fabs(u(i+1)))
			u(i) = fabs(u(i)) * sqrt(1.0 + sqr(u(i+1) / u(i)));
		else u(i) = fabs(u(i+1)) * sqrt(1.0 + sqr(u(i) / u(i+1)));
	}
	for (j=1; j<=n; j++) r(1, j) += u(1) * v(j);
	for (i=1; i<k; i++)
	  rotate(r, qt, n, i, r(i, i), -r(i+1, i));
}

void rotate(DoubleMatrix & r, DoubleMatrix & qt, int n, int i, float a, 
	    float b) {
  int j;
  double c,fact,s,w,y;
  
  if (a == 0.0) {
    c = 0.0;
    s = (b >= 0.0 ? 1.0 : -1.0);
  } else if (fabs(a) > fabs(b)) {
    fact = b / a;
    c = sign(1.0/sqrt(1.0+(fact*fact)),a);
    s = fact*c;
  } else {
    fact = a / b;
    s = sign(1.0 / sqrt(1.0 + (fact * fact)), b);
    c = fact * s;
  }
  for (j=i; j<=n; j++) {
    y = r(i, j);
    w = r(i+1, j);
    r(i, j) = c * y - s * w;
    r(i+1, j) = s * y + c * w;
  }
  for (j=1; j<=n; j++) {
    y = qt(i, j);
    w = qt(i+1, j);
    qt(i, j) = c * y - s * w;
    qt(i+1, j) = s * y + c * w;
  }
}

void rsolv(const DoubleMatrix & a, int n, const DoubleVector & d, 
	   DoubleVector & b) {
	int i,j;
	double sum;

	b(n) /= d(n);
	for (i=n-1; i>=1; i--) {
	  for (sum=0.0, j=i+1; j<=n; j++) sum += a(i, j) * b(j);
	  b(i)=(b(i) - sum) / d(i);
	}
}

double ccbSqrt(double f) { return sqrt(fabs(f)); }

double signedSqrt(double f){ return f<0 ? -ccbSqrt(double(f)) : ccbSqrt(double(f)); }

double signedSqr(double f){ if (f > 0.) return sqr(f); 
  else return -sqr(f); }

  /// returns the square root of the absolute value of the argument
  // returns sqrt(f) for f>0 
  double ccbSqrt(double f);
  /// returns the square root of the absolute value of the argument
  // returns sqrt(f) for f>0 or -sqrt(|f|) for f<0
  double signedSqrt(double f);
  /// returns f * f * sign(f)
  double signedSqr(double f);

double kinFn(double m1, double m2, double m3) {
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    cout << "m1 < m2 + m3 DEBUG\n";
    return 0.;
  }
  double squareplus  = pow(m1, 2) - pow((m2 + m3), 2);
  double squareminus = pow(m1, 2) - pow((m2 - m3), 2);

  return 0.5 * sqrt(squareplus * squareminus) / m1;
}

double lambda(double a, double b, double c) {
  if (fabs(a) < fabs(b) + fabs(c)) {
    throw("a < b + c in numerics.cpp:lambda\n");
  }
  return sqr(a + b - c) - 4.0 * a * b;
}
