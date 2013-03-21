
#ifndef MATHEMATICA_WRAPPERS
#define MATHEMATICA_WRAPPERS

#include <cmath>
#include "linalg.h"
#include "numerics.h"

DoubleMatrix Adj(const DoubleMatrix&);
ComplexMatrix Adj(const ComplexMatrix&);

double Conj(double);
Complex Conj(const Complex&);
DoubleMatrix Conj(const DoubleMatrix&);
ComplexMatrix Conj(const ComplexMatrix&);

double Cos(double);
double Sin(double);

unsigned Delta(unsigned, unsigned);

double Mass2(double);

template <typename Base, typename Exponent>
double Power(Base base, Exponent exp)
{
   return std::pow(base, exp);
}

double Re(double);
double Re(const Complex&);

double Sqrt(double);

int ThetaStep(int, int);

DoubleMatrix Tp(const DoubleMatrix&);
ComplexMatrix Tp(const ComplexMatrix&);

double trace(const DoubleMatrix&);
Complex trace(const ComplexMatrix&);

#define UNITMATRIX(rows) \
   unitMatrix<rows>()

template <int rows>
DoubleMatrix unitMatrix()
{
   DoubleMatrix u(rows,rows);
   for (int i = 1; i <= rows; ++i)
      u(i,i) = 1.0;
   return u;
}

#endif
