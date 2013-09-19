
/** \file utils.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

 */

#include "utils.h"
#include "physpars.h"

double frexp(const Complex & c, int * i) {
  int a, b;
  frexp(c.real(), &a);
  frexp(c.imag(), &b);
  if (fabs(a) > fabs(b)) return *i = a; 
  *i = b;
  return 0.;
}

void FPE_ExceptionHandler(int nErrType) {
  throw "SIGFPE"; ///< This reverts back to softsusy code!
}

int theta(double a) {
  int temp = 0;
  if (a > 0.0) temp = 1;
  return temp;
}

// Just sets precision and format of outputs
void outputCharacteristics(int n) {
  cin.setf(ios::scientific, ios::floatfield);
  cin.precision(n);
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(n);
  cerr.setf(ios::scientific, ios::floatfield);
  cerr.precision(n);
}

// Finds fractional difference between |a| and |b|
double toleranceCheck(double a, double b) {
  double sTin = fabs(a), sTout = fabs(b);
  double maxx = maximum(sTin, sTout);

  const double underflow = 1.0e-20;

  if (maxx < underflow) return 0.0;
  return fabs(1.0 - minimum(sTin, sTout) / maxx);
}

// Outputs a space if greater than zero, a minus otherwise.
// Useful for outputting negative numbers in rows
void printRow(ostream & out, double x) {

  // make it return a character when you've worked out the equivalent of printf

  double underflow = 1.0e-120;
  if (fabs(x) < underflow) x = 0.0; // Traps -0.0
  if (x >= 0.0) out << " " << x;
  else out << x;
}

bool testNan(double f) {
  return (f != f);
}

bool close(double m1, double m2, double tol) {
  double mmax = fabs(maximum(fabs(m1), fabs(m2)));
  double mmin = fabs(minimum(fabs(m1), fabs(m2)));
  double max_tol = tol * mmax;
  if (max_tol == 0.0 && mmax != 0.0 && tol != 0.0)
    return (mmax - mmin <= tol);

  return (mmax - mmin <= tol * mmax);
}

double sTfn(double sTins, double sTouts) {
  double sTin  = fabs(sTins);
  double sTout = fabs(sTouts);
  if (sTin < 1. && sTout < 1.) return fabs(sTin - sTout);
  else return fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
}

/// LCT: Difference between two drBarPars objects
void sumTol(const drBarPars & a, const drBarPars & b, DoubleVector & sT) {
  int k = 1;
  
  sT(k) = sTfn(a.mGluino, b.mGluino); k++;
  int i; for (i=1; i<=a.mh0.displayEnd(); i++) {
    sT(k) = sTfn(a.mh0(i), b.mh0(i)); k++;
  }
  for (i=1; i<=a.mA0.displayEnd(); i++) {
    sT(k) = sTfn(a.mA0(i), b.mA0(i)); k++;
  }
  sT(k) = sTfn(a.mHpm, b.mHpm); k++;
  for (i=1; i<=3; i++) {
    sT(k) = sTfn(a.msnu(i), b.msnu(i)); k++;
  }
  for (i=1; i<=2; i++) {
    sT(k) = sTfn(a.mch(i), b.mch(i)); k++;
  }
  for (i=1; i<=a.mneut.displayEnd(); i++) {
    sT(k) = sTfn(a.mneut(i), b.mneut(i)); k++;
  }
  int j; for (j=1; j<=3; j++)
    for(i=1; i<=2; i++) {
      sT(k) = sTfn(a.mu(i, j), b.mu(i, j)); k++;
      sT(k) = sTfn(a.md(i, j), b.md(i, j)); k++;
      sT(k) = sTfn(a.me(i, j), b.me(i, j)); k++;
    }
}



