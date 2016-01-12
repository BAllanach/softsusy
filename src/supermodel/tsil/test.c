#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "tsil.h"

#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

int main (int argc, char *argv[])
{
  int i;
  TSIL_COMPLEX res, dat[6], bolddat[6][3];
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
  TSIL_SetParameters (&result, x, y, z, u, v, qq);

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

  res = TSIL_GetFunction (&result, "Txuv");
  printf ("result is "); TSIL_cprintf (res); printf("\n");

  TSIL_GetData (&result, "T", dat);
  for (i=0; i<6; i++) {printf("T[%d] = ", i); TSIL_cprintf (dat[i]); printf ("\n"); }

  TSIL_GetBoldData (&result, "T", bolddat);
  for (i=0; i<6; i++) {printf("T[%d][1] = ", i); TSIL_cprintf (bolddat[i][1]); printf ("\n"); }

  TSIL_Vanalytic (x, y, 0.0L, y, s + 0.1L*I, qq, &res);
  printf ("result is "); TSIL_cprintf (res); printf("\n");

  TSIL_Manalytic (0.0L, x, 0.0L, x, x, s, &res);
  printf ("result is "); TSIL_cprintf (res); printf("\n");

  return 0;
}
