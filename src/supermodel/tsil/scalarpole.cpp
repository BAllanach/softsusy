// === scalarpole.cpp ===

// Sample code for calculation of scalar field pole squared mass.

// To compile: 

// g++ -o scalarpole scalarpole.c -L<dir> -ltsil -DTSIL_SIZE_<size>

// where <dir> is the directory containing libtsil.a, and <size> is
// either LONG or DOUBLE, which should agree with the choice made
// when compiling libtsil.a.
 
// Command-line arguments are: 
//        scalar mass squared = x, 
//        cubic coupling = g, 
//        quartic coupling = lambda, 
//        renormalization scale squared = qq. 

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <ctime>

using namespace std;

#include "tsil_cpp.h"  // Required TSIL header file

#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

int main (int argc, char *argv[])
{
  TSIL_DATA       result; // Top-level TSIL data object
  TSIL_REAL       qq;
  TSIL_REAL       x, g, lambda;
  TSIL_COMPLEXCPP pi1, pi1prime, pi2, s1, s2;
  TSIL_REAL       factor = 1.0L/(16.0L*PI*PI);

  // If incorrect number of args, print message on stderr and exit:
  if (argc != 5)
    TSIL_Error("main", "Expected 4 arguments: m^2, g, lambda, and Q^2", 1);

  // Note cast to appropriate floating-point type for safety
  x      = (TSIL_REAL) strtold(argv[1], (char **) NULL); 
  g      = (TSIL_REAL) strtold(argv[2], (char **) NULL);
  lambda = (TSIL_REAL) strtold(argv[3], (char **) NULL); 
  qq     = (TSIL_REAL) strtold(argv[4], (char **) NULL); 

  // All loop integrals have a common squared-mass argument x:
  TSIL_SetParameters (&result, x, x, x, x, x, qq);

  // For the pole mass calculation, evaluate two-loop integrals at s = x:
  TSIL_Evaluate (&result, x);

  // Assemble one- and two-loop mass squared results:
  // Note new wrappers for use in C++ (defined in tsil_cpp.h).
  pi1 = 0.5L*lambda*TSIL_A_ (x,qq) - 0.5L*g*g*TSIL_B_ (x,x,x,qq);

  pi1prime = -0.5L*g*g*TSIL_dBds_ (x, x, x, qq); 

  pi2 = - 0.5L*g*g*g*g*TSIL_GetFunction_ (&result, "M")
    - 0.5L*g*g*g*g*TSIL_GetFunction_ (&result, "Vzxyv")
    + g*g*g*TSIL_GetFunction_ (&result, "Uzxyv")
    - (1.0L/6.0L)*lambda*lambda*TSIL_GetFunction_ (&result, "Svyz")
    + 0.25L*lambda*g*g*pow(TSIL_GetFunction_ (&result, "Bxz"), 2)
    + 0.25L*lambda*lambda*TSIL_A_ (x,qq)*(TSIL_A_ (x,qq)/x + 1.0L)   
    - 0.5L*lambda*g*g*TSIL_A_ (x,qq)*TSIL_Bp_ (x, x, x, qq)
    - 0.25L*lambda*g*g*TSIL_I2p_ (x,x,x,qq);

  s1 = x + factor*pi1;
  s2 = x + factor*pi1 + factor*factor*(pi2 + pi1*pi1prime);

  cout << "Tree-level squared mass:    " << x  << endl; 
  cout << "One-loop pole squared mass: " << real(s1) << endl;
  cout << "Two-loop pole squared mass: " << real(s2) << endl; 

  return 0;
}
