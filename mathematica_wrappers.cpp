
#include "mathematica_wrappers.hpp"

DoubleMatrix Adj(const DoubleMatrix& m)
{
   return m.transpose();
}

ComplexMatrix Adj(const ComplexMatrix& m)
{
   return m.hermitianConjugate();
}

double Conj(double a)
{
   return a;
}

Complex Conj(const Complex& a)
{
   return a.conj();
}

DoubleMatrix Conj(const DoubleMatrix& m)
{
   return m;
}

ComplexMatrix Conj(const ComplexMatrix& m)
{
   return m.complexConjugate();
}

double Cos(double x)
{
   return cos(x);
}

double Sin(double x)
{
   return sin(x);
}

unsigned Delta(unsigned i, unsigned j)
{
   return i == j ? 1 : 0;
}

double Mass2(double m)
{
   return m * m;
}

double Re(double x)
{
   return x;
}

double Re(const Complex& x)
{
   return real(x);
}

int ThetaStep(int a, int b)
{
   return a <= b ? 1 : 0;
}

DoubleMatrix Tp(const DoubleMatrix& m)
{
   return m.transpose();
}

ComplexMatrix Tp(const ComplexMatrix& m)
{
   return m.transpose();
}

double trace(const DoubleMatrix& m)
{
   return m.trace();
}

Complex trace(const ComplexMatrix& m)
{
   return m.trace();
}

double Sqrt(double a)
{
   return std::sqrt(a);
}
