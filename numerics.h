
/** \file numerics.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief numerical routines - differential equation solver, differentiator
   and function minimiser for instance
*/

#ifndef NUMERICS_H
#define NUMERICS_H

#include "utils.h"
#include "mycomplex.h"
#include <iostream>
#include "def.h"
#include "linalg.h"
#include "twoloophiggs.h"
//#include "clooptools.h"
using namespace softsusy;

/// Comment if you want default softsusy behaviour
//#define USE_LOOPTOOLS

/// A single step of Runge Kutta (5th order), input: 
/// y and dydx (derivative of y), x is independent variable. yout is value
/// after step. derivs is a user-supplied function
void rungeKuttaStep(const DoubleVector & y, const DoubleVector & dydx, 
	     double x, double h, DoubleVector & yout, DoubleVector & yerr, 
	     DoubleVector (*derivs)(double, const DoubleVector &));

/// organises the variable step-size for Runge-Kutta evolution
int odeStepper(DoubleVector & y, const DoubleVector & dydx, double *x, double
		htry, double eps, DoubleVector & yscal, double *hdid, 
		double *hnext,		
		DoubleVector (*derivs)(double, const DoubleVector &));

/// Calculates log likelihood of a Poisson with k observed events, expecting
/// lambda>0. 
double lnLPoisson(unsigned k, double lambda);

/// Calculates likelihood of a Poisson with k observed events, expecting
/// lambda>0. 
double LPoisson(unsigned k, double lambda);

/// Organises integration of 1st order system of ODEs
int integrateOdes(DoubleVector & ystart, double x1, double x2, double eps,
		  double h1, double hmin, 
		  DoubleVector (*derivs)(double, const DoubleVector &),
		  int (*rkqs)
		  (DoubleVector & y, const DoubleVector & dydx, double *x,
		   double htry, double eps, DoubleVector & yscal, double
		   *hdid, double *hnext, 
		   DoubleVector (*derivs)(double, const DoubleVector &)));

/// func is user-supplied, h is an estimate of what step-size to start with
/// and err returns error flags
double calcDerivative(double (*func)(double), 
		     double x, double h, double *err);

/// f is user-defined function, minimum value returned in xmin. Based on a
/// golden section search
double findMinimum(double ax, double bx, double cx, double (*f)(double),
		   double tol, double *xmin);

void shft2(double & a, double & b, double & c); ///< a=b and b=c
/// a=b, b=c and c=d
void shft3(double & a, double & b, double & c, double & d); 

/// For calculation of PV functions
double integrandThreshbn(double x);
/// Returns real part of b function, less accurate than analytic expressions
double bIntegral(int n, double p, double m1, double m2, double mt);
DoubleVector dd(double x, const DoubleVector & y);

/// Passarino-Veltman function definition
double b0(double p, double m1, double m2, double q);
/// Passarino-Veltman function definition
double b1(double p, double m1, double m2, double q);
/// Passarino-Veltman function definition
double b22(double p,  double m1, double m2, double q);
/// Passarino-Veltman function definition
double c0(double m1, double m2, double m3);
/// Passarino-Veltman function definition
double d27(double m1, double m2, double m3, double m4);
/// Passarino-Veltman function definition
double d0(double m1, double m2, double m3, double m4);

// inlined PV functions
inline double a0(double m, double q) {
  if (m == 0.0) return 0.0;
  return sqr(m) * (1.0 - log(sqr(m / q)));
}

inline double ffn(double p, double m1, double m2, double q) {
  return a0(m1, q) - 2.0 * a0(m2, q) - 
    (2.0 * sqr(p) + 2.0 * sqr(m1) - sqr(m2)) * 
    b0(p, m1, m2, q);
}

inline double gfn(double p, double m1, double m2, double q) {
  return (sqr(p) - sqr(m1) - sqr(m2)) * b0(p, m1, m2, q) - a0(m1, q) 
    - a0(m2, q); 
}

inline double hfn(double p, double m1, double m2, double q) {
  return 4.0 * b22(p, m1, m2, q) + gfn(p, m1, m2, q);
}

inline double b22bar(double p, double m1, double m2, double q) {
  return b22(p, m1, m2, q) - 0.25 * a0(m1, q) - 0.25 * a0(m2, q);
}

inline double fB(const Complex & x) {
  /// special case at 0
  //  if (x.real() < EPSTOL * 10. && x.imag() < EPSTOL * 10.) return -1.;
  return (log(Complex(1.0 - x)) - x * log(Complex(1.0 - 1.0 / x)) 
	  - Complex(1.0)).real();
  } 

double dilogarg(double t);
double dilog(double x);

double integrandThreshbnr(double x);
Complex fnfn(double x);

/// Gaussian deviated random number, mean 0 variance 1. Don't re-set idum once
/// you've initially set it. Initialise with a NEGATIVE integer
double gasdev(long & idum);
/// Normally distributed random number between 0 and 1. Don't re-set idum once
/// you've initially set it. Initialise with a NEGATIVE integer
double ran1(long & idum);
/// Cauchy distribution ie 1 / [ pi gamma (1 + x^2/gamma^2) ]. For a width,
/// you must multiply the x coming out by the width. 
double cauchyRan(long & idum);

/// Returns the number of a bin that the data is in: from 1 to numBins in the
/// range (bins other than this range are also possible - you must deal with
/// them outside the function...)
int bin(double data, double start, double end, int numBins);

/// Adds logs of two numbers in a more careful way that avoids underflow
double logOfSum(double a, double b);

/// returns a random direction: total number of dimensions=n
DoubleVector getRandomDirection(int n, int & numChanged, long & idum); 

/// Calculates the vertical level of confidence level required to contain 
/// a fraction cl of area of the histogrammed binned likelihood l.
/// If err is true, a satisfactory answer could not be found for some reason.
double calcCL(double cl, const DoubleVector & l);

/// given a normalised binned likelihood vector l, calculates the fraction of
/// bins with likelihood less than or equal to y as a % of maximum:
/// approximates area OUTSIDE confidence level y
double calc1dFraction(double y, const DoubleVector & l);

/// These three functions are for the calculation of 2-loop log pieces of g-2
/// of the muon
double fps(double z);
double fs(double z);
double ffbar(double z);

Complex dilog(const Complex & x);

/// Trapezoidal function integration to accuracy tol
double trapzd(double (*func)(double), double a, double b, int n, 
	      double tol = 1.0e-3);
/// Driver for integration
double qtrap(double (*func)(double), double a, double b, double tol);
/// Mid point function integration
double midpnt(double (*func)(double), double a, double b, int n);
/// Integrate a function *func which is analytic to precision EPS between a
/// and b 
double qromb(double (*func)(double), double a, double b, double EPS);
/// Identical copy to facilitate double integral
double qromb2(double (*func)(double), double a, double b, double EPS);

/// This sums two exponentials nof the arguments in a way that tries to avoid
/// underflows. 
double sumOfExp(double a, double b);

/// Returns a value of the polynomial (y) and an error estimate (dy), given 
/// a bunch of points (xa) and their y values (ya).
void polint(const DoubleVector&  xa, const DoubleVector & ya, double x, 
	    double & y, double & dy);

/// Returns a 3 by 3 real mixing matrix. Input angles are standard CKM
/// parameterisation. If the phase d is not zero, the result is only an
/// approximation to the full complex matrix: see SOFTSUSY manual for
/// details. 
DoubleMatrix display3x3RealMixing(double theta12, double theta13, 
				  double theta23, double d);

/// Given a matrix v, determines which angles gave it
void getAngles(const DoubleMatrix & v, double & t12, double & t13, 
	       double & t23, double & d);
/// Evolves the dependent variables xi by one reversible mid-point step of
/// length tStep. Returns true if there is an error.
bool midPtStep(DoubleVector & xi, 
	      DoubleVector (*derivs)(double t, const DoubleVector & v), 
	      double tInitial, double tStep);
/// Reversible integrate with numSteps steps between tInitial and tFinal. 
/// Returns true if there's an error
bool integrateReversibly(DoubleVector & xi, 
			 DoubleVector (*derivs)(double t, 
						const DoubleVector & v), 
			 double tInitial, double tFinal, int numSteps);

/// useful for 2-loop mb/mt corrections
double fin(double mm1, double mm2);
double den(double a, int b); /// 1/a^b
#endif

