/* Temporary main program for SOFTSUMO testing. */

#include "supermodel.h"
#include "supermodel_defs.h"

#define NSTEPS 50

int main ()
{
  int i,which,loop = 0;
  TSIL_REAL Qi,Qf,Q0,dQ;

  int oneloop, twoloop;

  /* Disable annoying warning messages */
  fclose (stderr);

  SUMO_Initialize ();

  /* Set benchmark model: */
  SUMO_SetTestModel4 ();

  /* printf("Gluino:\n"); */
  /* SUMO_GluinoPole (1); */
  /* printf("%7.2Lf\t%7.2Lf\t%7.2Lf\t", Q, SUMO_SQRT(m2_gluino), SUMO_SQRT(M2_gluino)); */
  /* SUMO_GluinoPole (2); */
  /* printf("%7.2Lf\n", SUMO_SQRT(M2_gluino)); */


  /* for (i=0; i<4; i++) */
  /*   printf("m_phi0[%d] = %7.2Lf\n", i, SUMO_SGNSQRT(m2_phi0[i])); */

  /* Minimize 2-loop Veff, just for fun (vevs are already correct in
     benchmark model): */
  /* SUMO_Minimize_Veff (2); */

  /* Set up for stepping in Q */
  Q0 = Q; Qi = 580.0; Qf = 2.0*Q0;
  dQ = (Qf - Qi)/NSTEPS;

  /* This is the gluino running mass at the reference scale Q0: */
  /* printf("%7.2Lf\n", SUMO_SQRT(m2_gluino)); */

  /* Now do the scan... */
  SUMO_RGrun (Qi, 2);

  /* while (Q < Qf) { */
  /*   SUMO_GluinoPole (1); */
  /*   printf("%7.2Lf\t%7.2Lf\t%7.2Lf\t", Q, SUMO_SQRT(m2_gluino), SUMO_SQRT(M2_gluino)); */
  /*   SUMO_GluinoPole (2); */
  /*   printf("%7.2Lf\n", SUMO_SQRT(M2_gluino)); */
  /*   SUMO_RGrun (Q + dQ, 2); */
  /* } */

  /* exit(0); */

  /* which = 0; */

  /* while (Q < Qf) { */
  /*   SUMO_StopPole (which, 1); */
  /*   printf("%7.2Lf\t%7.2Lf\t%7.2Lf\t", Q, SUMO_SQRT(m2_stop[which]), SUMO_SQRT(M2_stop[which])); */
  /*   SUMO_StopPole (which, 2); */
  /*   printf("%7.2Lf\n", SUMO_SQRT(M2_stop[which])); */
  /*   SUMO_RGrun (Q + dQ, 2); */
  /* } */

  /* exit (0); */

  /* while (Q < Qf) { */
  /*   SUMO_SbotPole (which, 1); */
  /*   printf("%7.2Lf\t%7.2Lf\t%7.2Lf\t", Q, SUMO_SQRT(m2_sbot[which]), SUMO_SQRT(M2_sbot[which])); */
  /*   SUMO_SbotPole (which, 2); */
  /*   printf("%7.2Lf\n", SUMO_SQRT(M2_sbot[which])); */
  /*   SUMO_RGrun (Q + dQ, 2); */
  /* } */

  /* while (Q < Qf) { */
  /*   SUMO_SuLPole (which, 1); */
  /*   printf("%7.2Lf\t%7.2Lf\t%7.2Lf\t", Q, SUMO_SQRT(m2_suL[which]), SUMO_SQRT(M2_suL[which])); */
  /*   SUMO_SuLPole (which, 2); */
  /*   printf("%7.2Lf\n", SUMO_SQRT(M2_suL[which])); */
  /*   SUMO_RGrun (Q + dQ, 2); */
  /* } */



  SUMO_RGrun (1100.0, 2);
  /* SUMO_Minimize_Veff (2); */
  /* SUMO_Find_b_and_mu (2, 1); */

  /* for (i=0; i<4; i++) */
  /*   printf("m_phi0[%d] = %7.2Lf\n", i, SUMO_SGNSQRT(m2_phi0[i])); */

  printf("\nWe are at Q = %7.2Lf\n", Q);
  printf("\t\ttree\t\t1-loop\t\t2-loop\n");
  SUMO_GluinoPole (1);
  printf("Gluino:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_gluino), SUMO_SQRT(M2_gluino));
  SUMO_GluinoPole (2);
  printf("%7.2Lf\n", SUMO_SQRT(M2_gluino));

  SUMO_StopPole (0, 1);
  printf("Stop[0]:\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_stop[0]), SUMO_SQRT(M2_stop[0]));
  SUMO_StopPole (0, 2);
  printf("%7.2Lf\n", SUMO_SQRT(M2_stop[0]));

  SUMO_StopPole (1, 1);
  printf("Stop[1]:\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_stop[1]), SUMO_SQRT(M2_stop[1]));
  SUMO_StopPole (1, 2);
  printf("%7.2Lf\n", SUMO_SQRT(M2_stop[1]));

  SUMO_SbotPole (0, 1);
  printf("Sbot[0]:\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_sbot[0]), SUMO_SQRT(M2_sbot[0]));
  SUMO_SbotPole (0, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_sbot[0]));

  SUMO_SbotPole (1, 1);
  printf("Sbot[1]:\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_sbot[1]), SUMO_SQRT(M2_sbot[1]));
  SUMO_SbotPole (1, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_sbot[1]));

  SUMO_SuLPole (0, 1);
  printf("SuL[0]:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_suL[0]), SUMO_SQRT(M2_suL[0]));
  SUMO_SuLPole (0, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_suL[0]));

  SUMO_SuLPole (1, 1);
  printf("SuL[1]:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_suL[1]), SUMO_SQRT(M2_suL[1]));
  SUMO_SuLPole (1, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_suL[1]));

  SUMO_SuRPole (0, 1);
  printf("SuR[0]:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_suR[0]), SUMO_SQRT(M2_suR[0]));
  SUMO_SuRPole (0, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_suR[0]));

  SUMO_SuRPole (1, 1);
  printf("SuR[1]:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_suR[1]), SUMO_SQRT(M2_suR[1]));
  SUMO_SuRPole (1, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_suR[1]));

  SUMO_SdLPole (0, 1);
  printf("SdL[0]:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_sdL[0]), SUMO_SQRT(M2_sdL[0]));
  SUMO_SdLPole (0, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_sdL[0]));

  SUMO_SdLPole (1, 1);
  printf("SdL[1]:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_sdL[1]), SUMO_SQRT(M2_sdL[1]));
  SUMO_SdLPole (1, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_sdL[1]));

  SUMO_SdRPole (0, 1);
  printf("SdR[0]:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_sdR[0]), SUMO_SQRT(M2_sdR[0]));
  SUMO_SdRPole (0, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_sdR[0]));

  SUMO_SdRPole (1, 1);
  printf("SdR[1]:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_sdR[1]), SUMO_SQRT(M2_sdR[1]));
  SUMO_SdRPole (1, 2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_sdR[1]));

  oneloop = 2; twoloop = 1;

  SUMO_h0Pole (1);
  printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop,
  	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0]));
  SUMO_h0Pole (1.5);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0]));

  oneloop = 2; twoloop = 7;

  SUMO_h0Pole (1);
  printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop,
  	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0]));
  SUMO_h0Pole (2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0]));

  exit(0);

  /* oneloop = 2; twoloop = 2; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 2; twoloop = 3; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 2; twoloop = 4; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 2; twoloop = 5; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 2; twoloop = 6; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  oneloop = 2; twoloop = 7;

  SUMO_h0Pole (1);
  printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop,
  	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0]));
  SUMO_h0Pole (2);
  printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0]));


  /* oneloop = 2; twoloop = 8; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 2; twoloop = 9; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */
  
  /* oneloop = 3; twoloop = 1; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 3; twoloop = 2; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 3; twoloop = 3; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 3; twoloop = 4; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 3; twoloop = 5; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 3; twoloop = 6; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  /* oneloop = 3; twoloop = 7; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */
  
  /* oneloop = 3; twoloop = 8; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */
  
  /* oneloop = 3; twoloop = 9; */

  /* SUMO_M_h0 (oneloop, 0); */
  /* printf("h0 (%d, %d):\t%7.2Lf\t\t%7.2Lf\t\t", oneloop, twoloop, */
  /* 	 SUMO_SGNSQRT(m2_phi0[0]), SUMO_SGNSQRT(M2_phi0[0])); */
  /* SUMO_M_h0 (oneloop, twoloop); */
  /* printf("%7.2Lf\n", SUMO_SGNSQRT(M2_phi0[0])); */

  SUMO_H0Pole (1);
  printf("H0:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_phi0[1]), SUMO_SQRT(M2_phi0[1]));
  SUMO_H0Pole (2);
  printf("%7.2Lf\n", SUMO_SQRT(M2_phi0[1]));

  /* printf("\n"); */
  /* exit(0); */

  SUMO_G0Pole (1);
  printf("G0:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_phi0[2]), SUMO_SGNSQRT(M2_phi0[2]));
  printf("\n");
  /* SUMO_G0Pole (2); */
  /* printf("%7.2Lf\n", SUMO_SQRT(M2_phi0[2])); */

  SUMO_A0Pole (1);
  printf("A0:\t\t%7.2Lf\t\t%7.2Lf\t\t", SUMO_SQRT(m2_phi0[3]), SUMO_SQRT(M2_phi0[3]));
  SUMO_A0Pole (2);
  printf("%7.2Lf\n", SUMO_SQRT(M2_phi0[3]));

  return 0;
}
