#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "tsil.h"

/* 
   This program implements a TSIL calculation of the STU integrals: it
   takes input values of x, z, u, v, s, qq (as command line arguments,
   in that order) and prints the values of the integral functions
   U(x,z,u,v), V(x,z,u,v), S(x,u,v), T(x,u,v), T(u,x,v), T(v,x,u),
   Tbar(x,u,v), Tbar(u,x,v), Tbar(v,x,u), and B(x,z), and their "bold"
   variants.  General and timing information for the computation is
   also given.

   Compile as, e.g.

   gcc -DTSIL_SIZE_LONG -o tsilSTU baseSTU.c -L. -ltsil -lm

   (assuming libtsil.a is in the current directory, and was made using         
   long double floating point data sizes) and run as

   ./tsilSTU 1 3 4 5 10 1

*/

long double strtold (const char *, char **);

int main (int argc, char *argv[])
{
  TSIL_DATA result;
  TSIL_REAL x, z, u, v, s, qq;
  clock_t t0, t1;

  if (argc != 7)
    TSIL_Error("main", "Incorrect number of args; should be 6", 1);

  /* Get input parameters */
  x  = (TSIL_REAL) strtold(argv[1], (char **) NULL);
  z  = (TSIL_REAL) strtold(argv[2], (char **) NULL);
  u  = (TSIL_REAL) strtold(argv[3], (char **) NULL);
  v  = (TSIL_REAL) strtold(argv[4], (char **) NULL);
  s  = (TSIL_REAL) strtold(argv[5], (char **) NULL);
  qq = (TSIL_REAL) strtold(argv[6], (char **) NULL);

  /* Set parameters in result */
  TSIL_SetParametersSTU (&result, x, z, u, v, qq);

  /* Evaluate functions */
  t0 = clock ();
  TSIL_Evaluate (&result, s);
  t1 = clock ();

  /* Output general information and results */
  TSIL_PrintInfo ();
  TSIL_PrintVersion ();
  TSIL_PrintStatus (&result);
  TSIL_PrintDataM (&result);

  printf("\nTotal calculation time (s): %lf\n",
	 difftime(t1, t0)/CLOCKS_PER_SEC);

  return 0;
}
