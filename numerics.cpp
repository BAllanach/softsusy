
/** \file numerics.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "numerics.h"

double log1minusx(double x) {
  if (x > 1.) return 0.; 
  else if (close(1., x, EPSTOL)) return 6.66e66;
  else if (fabs(x) > 0.125) return log(1. - x);
  else if (x < 1.e-200) return 0.;
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
    if (xnew == *x) 
      {
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

inline void shft2(double & a, double & b, double c) { a = b; b = c; }

inline void shft3(double & a, double & b, double & c, double d) { 
  a = b; b = c; c = d;
}

double findMinimum(double ax, double bx, double cx, double (*f)(double),
		   double tol, double *xmin)
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
    *xmin = x1; 
    return f1; 
  } else {
    *xmin = x2; 
    return f2; 
  }
}



DoubleVector dd(double x, const DoubleVector & y) {
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
  int i; for (i=1; i<=nInt; i++) xn = xn * x;
  return xn * 
    log( ((1 - x) * sqr(m1Int) + x * sqr(m2Int) - x * (1 - x) *
	  sqr(pInt) - iEpsilon)
	 / sqr(mtInt));
}

DoubleVector dilogarg(double t, const DoubleVector & y) {

  const double eps = TOLERANCE * 1.0e-20;

  DoubleVector dydx(1);
  dydx(1) = -log(fabs(1 - t + eps)) / (t + eps);

  return dydx;
}

/*
double dilog(double x) {
  // Set global variables so that integration function can access them
  double from = 0.0, to = x, guess = 0.1, hmin = TOLERANCE * 1.0e-5;

  DoubleVector v(1); 
  double eps = TOLERANCE * 1.0e-5;
  v(1) = 1.0; 

  // Runge-Kutta, f(b) = int^b0 I(x) dx, I is integrand => d f / db = I(b)
  // odeint has a problem at f(0): therefore, define f'(b)=f(b)+1
  integrateOdes(v, from, to, eps, guess, hmin, dilogarg, odeStepper); 
  
  return v(1) - 1.0;
}
*/

// Returns real part of integral
double bIntegral(int n1, double p, double m1, double m2, double mt) {
  // Set global variables so that integration function can access them
  nInt = n1; pInt = p; m1Int = m1; m2Int = m2; mtInt = mt;
  double from = 0.0, to = 1.0, guess = 0.1, hmin = TOLERANCE * 1.0e-5;
  
  DoubleVector v(1); double eps = TOLERANCE * 1.0e-3;
  v(1) = 1.0; 

  // Runge-Kutta, f(b) = int^b0 I(x) dx, I is integrand => d f / db = I(b)
  // odeint has a problem at f(0): therefore, define f'(b)=f(b)+1
  integrateOdes(v, from, to, eps, guess, hmin, dd, odeStepper); 
  
  return v(1) - 1.0;
}

/// Decides level at which one switches to p=0 limit of calculations
const double pTolerance = 1.0e-6; 

double fB(const Complex & a) {
  /// First, special cases at problematic points
  double x = a.real(), y = a.imag();
  if (fabs(x) < pTolerance) {
    double ans = -1. - x + sqr(x) * 0.5;
    return ans;
  }
  if (close(x, 1., pTolerance)) { 
    if (close(x, 1., EPSTOL)) return -1.;
    double eps = x - 1.;
    double ans = -1. + eps + sqr(eps) * 0.5;
    return ans;
  }
  
  Complex ans = log(1. - a) -1. - a * log(1.0 - 1.0 / a);
  //double ans = log1minusx(a.real()) - 1. - a.real() * log1minusx(1.0 / a.real());
  return ans.real();
}
  
/*
  Analytic expressions follow for above integrals: sometimes useful!
  From hep-ph/9606211
  Note it returns the REAL PART ONLY. 
*/
double b0(double p, double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  double b0l = B0(p*p, m1*m1, m2*m2).real();
  //  return B0(p*p, m1*m1, m2*m2).real();
#endif

  double ans  = 0.;
  double mMin = minimum(fabs(m1), fabs(m2));
  double mMax = maximum(fabs(m1), fabs(m2));

  double pSq = sqr(p), mMinSq = sqr(mMin), mMaxSq = sqr(mMax);
  double s = 0.;
  /// Try to increase the accuracy of s
  double dm = mMaxSq - mMinSq;
  s = pSq + dm;

  char * methodId = (char *) "";

  double pTest = sqr(p) / sqr(mMax);
  /// p is not 0  
  if (pTest > pTolerance) {  
    methodId = (char *) "B0A";
    
    Complex iEpsilon(0.0, EPSTOL * sqr(mMax));
    
    Complex xPlus, xMinus;
    
    xPlus = (s + sqrt(sqr(s) - 4. * sqr(p) * (sqr(mMax) - iEpsilon))) /
      (2. * sqr(p));
    xMinus = 2. * (sqr(mMax) - iEpsilon) / 
      (s + sqrt(sqr(s) - 4. * sqr(p) * (sqr(mMax) - iEpsilon)));
    
    ans = -2.0 * log(p / q) - fB(xPlus) - fB(xMinus);
  } else {
    if (close(m1, m2, EPSTOL)) {
      methodId = (char *) "B0B";
      
      ans = - log(sqr(m1 / q));
    } else {
      double Mmax2 = sqr(mMax), Mmin2 = sqr(mMin); 
      if (Mmin2 < sqr(TOLERANCE)) {
	methodId = (char *) "B0C";

	ans = 1.0 - log(Mmax2 / sqr(q));
      } else {
	methodId = (char *) "B0D";

	ans = 1.0 - log(Mmax2 / sqr(q)) + Mmin2 * log(Mmax2 / Mmin2) 
	  / (Mmin2 - Mmax2);
      }
    }
  }   
  
  if (!close(b0l, ans, 1.0e-3)) {
    cout << methodId << " ";
    cout << "DEBUG Err: DB0(" << p << ", " << m1 << ", " << m2 
	 << ", "  << q << ")=" << 1.-b0l/ans << endl;
    cout << "SOFTSUSY  B0=" << ans << endl;
    cout << "LOOPTOOLS B0=" << b0l << endl;
  }
  return ans;
}

/// Note that b1 is NOT symmetric in m1 <-> m2!!!
double b1(double p, double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  double b1l = -B1(p*p, m1*m1, m2*m2).real();
  //    return b1l;
#endif

  double ans = 0.;
  double pTest = sqr(p) / maximum(sqr(m1), sqr(m2));

  //  cout << "pT=" << pTest << endl; ///< DEBUG

  char * methodId = (char *) "";

  if (pTest > pTolerance * 1.0e3) {
    methodId = (char *) "B1A";

    ans = (a0(m2, q) - a0(m1, q) + (sqr(p) + sqr(m1) - sqr(m2)) 
	   * b0(p, m1, m2, q)) / (2.0 * sqr(p)); 
  } else if (fabs(m1) > EPSTOL && !close(m1, m2, EPSTOL) 
	     && fabs(m2) > EPSTOL) { ///< checked
    methodId = (char *) "B1B";

    double Mmax2 = maximum(sqr(m1) , sqr(m2)), x = sqr(m2 / m1);
    ans = 0.5 * (-log(Mmax2 / sqr(q)) + 0.5 + 1.0 / (1.0 - x) + log(x) /
		 sqr(1.0 - x) - theta(1.0 - x) * log(x)); ///< checked
  } else {
    methodId = (char *) "B1C";

    ans = bIntegral(1, p, m1, m2, q); 
  }

  if (!close(b1l, ans, 5.0e-2)) {
    cout << methodId << " Test=" << pTest << " ";
    cout << "DEBUG Err: Db1(" << p << ", " << m1 << ", " << m2 
	 << ", "  << q << ")=" << 1.-b1l/ans << endl;
    cout << "SOFTSUSY  B1=" << ans << " B0=" << b0(p, m1, m2, q) << endl;
    cout << "LOOPTOOLS B1=" << b1l << " B0=" << B0(p*p, m1*m1, m2*m2).real() 
	 << endl;
  }
  return ans;
}

double b22(double p,  double m1, double m2, double q) {

#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  //  return B00(p*p, m1*m1, m2*m2).real();
#endif

  double answer;
  
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
    double m2p = m2, m1p = m1, pp = p, qp = q;    

    double ans = 1.0 / 6.0 * 
      (0.5 * (a0(m1p, qp) + a0(m2p, qp)) + (sqr(m1p) + sqr(m2p) - 0.5 * sqr(pp))
       * b0Save + (sqr(m2p) - sqr(m1p)) / (2.0 * sqr(pp)) *
       (a0(m2p, qp) - a0(m1p, qp) - (sqr(m2p) - sqr(m1p)) * b0Save) +
       sqr(m1p) + sqr(m2p) - sqr(pp) / 3.0);

    answer = 1.0 / 6.0 * 
      (0.5 * (a0(m1, q) + a0(m2, q)) + (sqr(m1) + sqr(m2) - 0.5 * sqr(p))
       * b0Save + (sqr(m2) - sqr(m1)) / (2.0 * sqr(p)) *
       (a0(m2, q) - a0(m1, q) - (sqr(m2) - sqr(m1)) * b0Save) +
       sqr(m1) + sqr(m2) - sqr(p) / 3.0);

    answer = ans;
  }

  return answer;
}

// debugged 23.01.07 - thanks to Shindou Tetsuo 
double d0(double m1, double m2, double m3, double m4) {
  if (close(m1, m2, EPSTOL)) {
    double m2sq = sqr(m2), m3sq = sqr(m3), m4sq = sqr(m4);

    if (close(m2, m3, EPSTOL) && close(m2, m4, EPSTOL)) 
      return 1.0 / (6.0 * sqr(m2sq));
    else if (close(m2, m3, EPSTOL)) {
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
    double m1n = m1 + TOLERANCE * 0.01;
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
  //  return C0(psq, psq, psq, m1*m1, m2*m2, m3*m3).real();
#endif

  if (close(m2, m3, EPSTOL)) {
    if (close(m1, m2, EPSTOL)) {
      return ( - 0.5 / sqr(m2) ); // checked 14.10.02
    }
    else {
      return ( sqr(m1) / sqr(sqr(m1)-sqr(m2) ) * log(sqr(m2)/sqr(m1))
               + 1.0 / (sqr(m1) - sqr(m2)) ) ; // checked 14.10.02
    }
  }
  else
    if (close(m1, m2, EPSTOL)) {
      return ( - ( 1.0 + sqr(m3) / (sqr(m2)-sqr(m3)) * log(sqr(m3)/sqr(m2)) )
               / (sqr(m2)-sqr(m3)) ) ; // checked 14.10.02
    }
    else
      if (close(m1, m3, EPSTOL)) {
        return ( - (1.0 + sqr(m2) / (sqr(m3)-sqr(m2)) * log(sqr(m2)/sqr(m3))) 
                 / (sqr(m3)-sqr(m2)) ); // checked 14.10.02
      }
      else return (1.0 / (sqr(m2) - sqr(m3)) * 
		   (sqr(m2) / (sqr(m1) - sqr(m2)) *
		    log(sqr(m2) / sqr(m1)) -
		    sqr(m3) / (sqr(m1) - sqr(m3)) *
		    log(sqr(m3) / sqr(m1))) );
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
    for (j=NTAB+7;j>=0;j--) {
      k=(idum)/IQ;
      idum=IA*(idum-k*IQ)-IR*k;
      if (idum < 0) idum += IM;
      if (j < NTAB) iv[j] = idum;
    }
    iy=iv[0];
  }
  k=(idum)/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
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

double trapzd(double (*func)(double), double a, double b, int n, double EPS) {
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
  if (expression < 0.) return -6.66e66; 
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
  if (exp1 < 0.) exp1 = -6.66e66; 
  double exp2  = edgefn(q, l, chi);
  if (exp2 < 0.) exp2 = -6.66e66; 
  double exp3  = (q * l - xi * chi) * (xi - l) / (xi * l);
  if (exp3 < 0.) exp3 = -6.66e66; 
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
  if (ans < 0.) return -6.66e66;
  return sqrt(ans / (4.0 * l * xi));
}

double lqnear(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double q     = sqr(mSq);
  double expression = (q - xi) * (xi - l) / xi;
  if (expression < 0.) return -6.66e66;
  return sqrt(expression);
}

double lqfar(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double q     = sqr(mSq);
  double expression = (q - xi) * (l - chi) / l;
  if (expression < 0.) return -6.66e66;
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
  if (mlqother < 0.) mlqother = -6.66e66;
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

  double diff = 6.66e66, delta = 0.;
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
