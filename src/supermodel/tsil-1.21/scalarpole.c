/*   === scalarpole.c ===
 *
 *   Sample code for calculation of scalar field pole squared mass.
 *
 *   To compile: 
 *
 *   cc -o scalarpole scalarpole.c -L<dir> -ltsil -lm -DTSIL_SIZE_<size>
 *
 *   where <dir> is the directory containing libtsil.a, and <size> is
 *   either LONG or DOUBLE, which should agree with the choice made
 *   when compiling libtsil.a.
 * 
 *   Command-line arguments are: 
 *          scalar mass squared = x, 
 *          cubic coupling = g, 
 *          quartic coupling = lambda, 
 *          renormalization scale squared = qq. 
 *
 *    Run as, for example: ./scalarpole 1 2 3 1
 */

#include <stdio.h>
#include "tsil.h"  /* Required TSIL header file */

#define PI 4.0L*TSIL_ATAN(1.0L) /* Uses arctan function defined in tsil.h */

long double strtold(const char *, char **);

int main (int argc, char *argv[])
{
  TSIL_DATA    result; /* Top-level TSIL data object */
  TSIL_REAL    qq;     /* Ensures correct basic type; see also TSIL_COMPLEX */ 
  TSIL_REAL    x, g, lambda;
  TSIL_COMPLEX pi1, pi1prime, pi2, s1, s2;
  TSIL_REAL    factor = 1.0L/(16.0L*PI*PI);

  /* If incorrect number of args, print message on stderr and exit: */
  if (argc != 5)
    TSIL_Error("main", "Expected 4 arguments: m^2, g, lambda, and Q^2", 1);

  /* Note cast to appropriate floating-point type for safety */
  x      = (TSIL_REAL) strtold(argv[1], (char **) NULL); 
  g      = (TSIL_REAL) strtold(argv[2], (char **) NULL);
  lambda = (TSIL_REAL) strtold(argv[3], (char **) NULL); 
  qq     = (TSIL_REAL) strtold(argv[4], (char **) NULL); 

  /* All loop integrals have a common squared-mass argument x: */
  TSIL_SetParameters (&result, x, x, x, x, x, qq); 

  /* For the pole mass calculation, evaluate two-loop integrals at s = x: */
  TSIL_Evaluate (&result, x);

  /* Assemble one- and two-loop mass squared results: */
  pi1 = 0.5L*lambda*TSIL_A(x,qq) - 0.5L*g*g*TSIL_B(x,x,x,qq);

  pi1prime = -0.5L*g*g*TSIL_dBds(x, x, x, qq); 

  pi2 = - 0.5L*g*g*g*g*TSIL_GetFunction(&result, "M")
        - 0.5L*g*g*g*g*TSIL_GetFunction(&result, "Vzxyv")
               + g*g*g*TSIL_GetFunction(&result, "Uzxyv")
    - (1.0L/6.0L)*lambda*lambda*TSIL_GetFunction(&result, "Svyz")
    + 0.25L*lambda*g*g*TSIL_POW(TSIL_GetFunction(&result, "Bxz"), 2)
    + 0.25L*lambda*lambda*TSIL_A(x,qq)*(TSIL_A(x,qq)/x + 1.0L)   
    - 0.5L*lambda*g*g*TSIL_A(x,qq)*TSIL_Bp(x, x, x, qq)      
    - 0.25L*lambda*g*g*TSIL_I2p(x,x,x,qq);               

  s1 = x + factor*pi1;
  s2 = x + factor*pi1 + factor*factor*(pi2 + pi1*pi1prime);

  printf("Tree-level squared mass:    %lf\n", (double) x); 
  printf("One-loop pole squared mass: %lf\n", (double) s1); 
  printf("Two-loop pole squared mass: %lf\n", (double) s2); 

  return 0;
}
