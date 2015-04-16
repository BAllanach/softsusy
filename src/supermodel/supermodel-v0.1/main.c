/* Temporary main program for SOFTSUMO testing. */

/* Calculates 1- and 2-loop gluino pole masses for a range of
   renormalization scales (Q = 500 to 1280 GeV). */

#include "supermodel.h"
#define NSTEPS 50

int main ()
{
  int i,loop = 0;
  TSIL_REAL Qi,Qf,Q0,dQ;

  /* Disable annoying warning messages */
  fclose (stderr);

  /* Set benchmark model: */
  SUMO_SetTestModel ();

  /* Minimize 2-loop Veff, just for fun (vevs are already correct in
     benchmark model): */
  SUMO_Minimize_Veff (2);

  /* Set up for stepping in Q */
  Q0 = Q; Qi = 500.0; Qf = 2.0*Q0;
  dQ = (Qf - Qi)/NSTEPS;

  /* This is the gluino running mass at the reference scale Q0: */
  printf("%Lf\n", SUMO_SQRT(m2_gluino));

  /* Now do the scan... */
  SUMO_RGrun (Qi, 2);

  while (Q < Qf) {
    SUMO_GluinoPole (1);
    printf("%Lf\t%Lf\t", Q, SUMO_SQRT(M2_gluino));
    SUMO_GluinoPole (2);
    printf("%Lf\n", SUMO_SQRT(M2_gluino));
    SUMO_RGrun (Q + dQ, 2);
  }

  return 0;
}
