/* 
   This program implements the most general TSIL calculation: it takes
   input values of x, y, z, u, v, s, qq (as command line arguments, in
   that order) and prints the values of all integral functions,
   including "bold" variants.  General and timing information for the
   computation are also printed.

   Compile as, e.g.

   gcc -DTSIL_SIZE_LONG -o tsil basecalc.c -L. -ltsil -lm

   (assuming libtsil.a is in the current directory, and was made using         
   long double floating point data sizes) and run as

   ./tsil 1 2 3 4 5 10 1

*/

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

// long double strtold (const char *, char **);

int main (int argc, char *argv[])
{
  TSIL_DATA result;
  TSIL_REAL x, y, z, u, v, s, qq;
  clock_t t0, t1;

  if (argc != 8)
    TSIL_Error ("main", "Incorrect number of args; should be 7", 1);

  /* Get input parameters */
  x  = (TSIL_REAL) strtold(argv[1], (char **) NULL);
  y  = (TSIL_REAL) strtold(argv[2], (char **) NULL);
  z  = (TSIL_REAL) strtold(argv[3], (char **) NULL);
  u  = (TSIL_REAL) strtold(argv[4], (char **) NULL);
  v  = (TSIL_REAL) strtold(argv[5], (char **) NULL);
  s  = (TSIL_REAL) strtold(argv[6], (char **) NULL);
  qq = (TSIL_REAL) strtold(argv[7], (char **) NULL);

  /* Set parameters in result */
  TSIL_SetParameters_ (&result, x, y, z, u, v, qq);

  /* Evaluate functions */
  t0 = clock ();
  TSIL_Evaluate_ (&result, s);
  t1 = clock ();

  /* Output general information and results */
  TSIL_PrintInfo_ ();
  TSIL_PrintVersion_ ();
  TSIL_PrintStatus_ (&result);
  TSIL_PrintDataM_ (&result);

  printf("\nTotal calculation time (s): %lf\n",
	 difftime(t1, t0)/CLOCKS_PER_SEC);

  return 0;
}
