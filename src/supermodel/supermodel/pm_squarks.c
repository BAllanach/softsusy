/* Squark pole masses */

#include "supermodel.h"
#include "sumo_squarks.h"

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Implements expansion around tree masses. */

int SUMO_StopPole (int i, int loop_order)
{
  int j = 0;
  char funcname[] = "SUMO_StopPole";

  if (i==0) j = 1;

  if (loop_order > 2) 
    SUMO_Error(funcname, "Invalid loop order specified.", 3);

  CM2_stop[i] = m2_stop[i];

  if (loop_order > 0)
    CM2_stop[i] += SUMO_oneloopfactor * pi1_stop (i, i, m2_stop[i]);

  if (loop_order > 1) {
    CM2_stop[i] += SUMO_twoloopfactor *
      (  pi20tilde_stop (i, m2_stop[i])
       + pi21tilde_stop (i, m2_stop[i])
       + pi22tilde_stop (i, m2_stop[i])
       + pi1_stop (i, j, m2_stop[i]) * pi1_stop(j, i, m2_stop[i])/(m2_stop[i] - m2_stop[j])
       );
  }
  M2_stop[i] = TSIL_CREAL(CM2_stop[i]);

  return 0;    
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

int SUMO_SbotPole (int i, int loop_order)
{
  int j = 0;
  char funcname[] = "SUMO_SbotPole";

  if (i==0) j = 1;

  if (loop_order > 2) 
    SUMO_Error(funcname, "Invalid loop order specified.", 3);

  CM2_sbot[i] = m2_sbot[i];

  if (loop_order > 0)
    CM2_sbot[i] += SUMO_oneloopfactor * pi1_sbot (i, i, m2_sbot[i]);

  if (loop_order > 1) {
    CM2_sbot[i] += SUMO_twoloopfactor *
      (  pi20tilde_sbot (i, m2_sbot[i])
       + pi21tilde_sbot (i, m2_sbot[i])
       + pi22tilde_sbot (i, m2_sbot[i])
       + pi1_sbot (i, j, m2_sbot[i]) * pi1_sbot(j, i, m2_sbot[i])/(m2_sbot[i] - m2_sbot[j])
       );
  }
  M2_sbot[i] = TSIL_CREAL(CM2_sbot[i]);

  return 0;    
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

int SUMO_SuLPole (int i, int loop_order)
{
  int j = 0;
  char funcname[] = "SUMO_SuLPole";

  if (i==0) j = 1;

  if (loop_order > 2) 
    SUMO_Error(funcname, "Invalid loop order specified.", 3);

  CM2_suL[i] = m2_suL[i];

  if (loop_order > 0)
    CM2_suL[i] += SUMO_oneloopfactor * pi1_suL (i, m2_suL[i]);

  if (loop_order > 1) {
    CM2_suL[i] += SUMO_twoloopfactor *
      (  pi20tilde_suL (i, m2_suL[i])
       + pi21tilde_suL (i, m2_suL[i])
       + pi22tilde_suL (i, m2_suL[i])
       /* + pi1_suL (i, j, m2_suL[i]) * pi1_suL(j, i, m2_suL[i])/(m2_suL[i] - m2_suL[j]) */
	 );
  }
  M2_suL[i] = TSIL_CREAL(CM2_suL[i]);

  return 0;    
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

int SUMO_SdLPole (int i, int loop_order)
{
  int j = 0;
  char funcname[] = "SUMO_SdLPole";

  if (i==0) j = 1;

  if (loop_order > 2) 
    SUMO_Error(funcname, "Invalid loop order specified.", 3);

  CM2_sdL[i] = m2_sdL[i];

  if (loop_order > 0)
    CM2_sdL[i] += SUMO_oneloopfactor * pi1_sdL (i, m2_sdL[i]);

  if (loop_order > 1) {
    CM2_sdL[i] += SUMO_twoloopfactor *
      (  pi20tilde_sdL (i, m2_sdL[i])
       + pi21tilde_sdL (i, m2_sdL[i])
       + pi22tilde_sdL (i, m2_sdL[i])
       /* + pi1_sdL (i, j, m2_sdL[i]) * pi1_sdL(j, i, m2_sdL[i])/(m2_sdL[i] - m2_sdL[j]) */
       );
  }
  M2_sdL[i] = TSIL_CREAL(CM2_sdL[i]);

  return 0;    
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

int SUMO_SuRPole (int i, int loop_order)
{
  char funcname[] = "SuRPole";
  int j = 0;

  if (i==0) j = 1;

  if (loop_order > 2) 
    SUMO_Error(funcname, "Invalid loop order specified.", 3);

  CM2_suR[i] = m2_suR[i];

  if (loop_order > 0)
    CM2_suR[i] += SUMO_oneloopfactor * pi1_suR (i, m2_suR[i]);

  if (loop_order > 1) {
    CM2_suR[i] += SUMO_twoloopfactor *
      (  pi20tilde_suR (i, m2_suR[i])
       + pi21tilde_suR (i, m2_suR[i])
       + pi22tilde_suR (i, m2_suR[i])
       /* + pi1_suR (i, j, m2_suR[i]) * pi1_suR(j, i, m2_suR[i])/(m2_suR[i] - m2_suR[j]) */
       );
  }
  M2_suR[i] = TSIL_CREAL(CM2_suR[i]);

  return 0;    
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

int SUMO_SdRPole (int i, int loop_order)
{
  char funcname[] = "SdRPole";
  int j = 0;

  if (i==0) j = 1;

  if (loop_order > 2) 
    SUMO_Error(funcname, "Invalid loop order specified.", 3);

  CM2_sdR[i] = m2_sdR[i];

  if (loop_order > 0)
    CM2_sdR[i] += SUMO_oneloopfactor * pi1_sdR (i, m2_sdR[i]);

  if (loop_order > 1) {
    CM2_sdR[i] += SUMO_twoloopfactor *
      (  pi20tilde_sdR (i, m2_sdR[i])
       + pi21tilde_sdR (i, m2_sdR[i])
       + pi22tilde_sdR (i, m2_sdR[i])
       /* + pi1_sdR (i, j, m2_sdR[i]) * pi1_sdR(j, i, m2_sdR[i])/(m2_sdR[i] - m2_sdR[j]) */
       );
  }
  M2_sdR[i] = TSIL_CREAL(CM2_sdR[i]);

  return 0;    
}
