/*   
   Computes pole mass of gluino at one or two loops, using
   eqs. (5.1)-(5.3) of hep-ph/0509115 
*/

#include "supermodel.h"
#define MAXITERS 10
#define ITERTOL  1.e-5

/* Self energy functions defined in se_gluino.c: */
SUMO_COMPLEX pi1tilde_gluino (void);
SUMO_COMPLEX dpi1tildedM2_gluino (void);
SUMO_COMPLEX dpi1tildedM2_top (void);
SUMO_COMPLEX dpi1tildedM2_stop (int);
SUMO_COMPLEX dpi1tildedM2_sbot (int);
SUMO_COMPLEX dpi1tildedM2_suL (int);
SUMO_COMPLEX dpi1tildedM2_suR (int);
SUMO_COMPLEX dpi1tildedM2_sdL (int);
SUMO_COMPLEX dpi1tildedM2_sdR (int);
SUMO_COMPLEX pi2tilde_gluino (void);

/* Self energy function defined in se_top.c: */
SUMO_COMPLEX pi1tilde_top (void);

/* Self energy functions defined in se1_squarks.c: */
SUMO_COMPLEX pi1_stop (int, int, TSIL_REAL);
SUMO_COMPLEX pi1_sbot (int, int, TSIL_REAL);
SUMO_COMPLEX pi1_suL (int, TSIL_REAL);
SUMO_COMPLEX pi1_suR (int, TSIL_REAL);
SUMO_COMPLEX pi1_sdL (int, TSIL_REAL);
SUMO_COMPLEX pi1_sdR (int, TSIL_REAL);

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

/* Evaluates gluino pole mass at specified loop order.  

   If expand_around_pole == NONE (or 1) then lagrangian masses are
   used everywhere.

   If expand_around_pole == GLUINO (or 2) then the self energy is
   expanded around gluino pole mass only, and iterated until the
   relative change in cabs(m2) from one iteration to the next is less
   than ITERTOL.

   If expand_around_pole == FULL (or 3) then we expand around gluino,
   squark, and top quark pole masses and iterate until the relative
   change in cabs(m2) from one iteration to the next ie less than
   ITERTOL.  Note that squark and top pole masses must be set first!

   The return value is:
     * 0 for evaluation with lagrangian masses.
     * the number of iterations needed in the case of expansion around
       pole masses.
     * -1 if the iteration did not converge in MAXITERS tries.

*/

int SUMO_GluinoPole (int loop_order, int expand_around_pole)
{
  TSIL_COMPLEX cm2_result, current, pi1, dpi1, pi2;
  TSIL_REAL m2_gluino_saf = m2_gluino; /* Save lagrangian mass^2, restore later */
  TSIL_REAL m2_top_saf, m_top_saf, m2_stop_saf[2], m2_sbot_saf[2];
  TSIL_REAL m2_suL_saf[2], m2_suR_saf[2], m2_sdL_saf[2], m2_sdR_saf[2];
  int i, iter, retval = -1;
  /* char funcname[] = "GluinoPole"; */

  if (0 == are_tree_couplings_updated) SUMO_Tree_Couplings ();

  /* If expanding around pole masses, save lagrangian masses, restore
     later... */
  if (expand_around_pole == FULL) {

    for (i=0; i<2; i++) {
      m2_stop_saf[i] = m2_stop[i];
      m2_sbot_saf[i] = m2_sbot[i];
      m2_suL_saf[i] = m2_suL[i];
      m2_suR_saf[i] = m2_suR[i];
      m2_sdL_saf[i] = m2_sdL[i];
      m2_sdR_saf[i] = m2_sdR[i];
    }
    m2_top_saf = m2_top;
    m_top_saf = m_top;
  }

  /* Tree mass^2 is the starting value in any case: */
  cm2_result = m2_gluino;

  if (expand_around_pole == NONE || loop_order < 2) {

    if (loop_order > 0)
      cm2_result += pi1tilde_gluino ();

    if (loop_order > 1)
      cm2_result += pi2tilde_gluino ();

    retval = 0;
  }
  else if (expand_around_pole == GLUINO) {

    for (iter = 0; iter < MAXITERS; iter++) {

      /* Save current value */
      current = cm2_result;

      /* Then tree value is just original lagrangian mass^2: */
      cm2_result = m2_gluino_saf;

      /* These SE functions get evaluated with s = m2_gluino, which is
	 now the re part of the current pole mass: */
      pi1 = pi1tilde_gluino ();
      dpi1 = dpi1tildedM2_gluino();
      pi2 = pi2tilde_gluino ();

      cm2_result += pi1;
      cm2_result += -TSIL_CREAL(pi1)*dpi1;
      cm2_result += pi2;

      /* For next iteration, set tree mass to real part of current pole mass. */
      m2_gluino = TSIL_CREAL (cm2_result);
      m_gluino = TSIL_SQRT (m2_gluino);

      if (TSIL_CABS(cm2_result - current)/TSIL_CABS(current) < ITERTOL) {
	retval = iter+1;
	break;
      }
    }
  }
  else {
    /* Expand around both gluino and squark pole masses... */

    /* Set squark and top masses to be pole masses: */
    for (i=0; i<2; i++) {
      m2_stop[i] = M2_stop[i];
      m2_sbot[i] = M2_sbot[i];
      m2_suL[i] = M2_suL[i];
      m2_suR[i] = M2_suR[i];
      m2_sdL[i] = M2_sdL[i];
      m2_sdR[i] = M2_sdR[i];
    }
    m_top = M_top;
    m2_top = M2_top;

    for (iter = 0; iter < MAXITERS; iter++) {

      /* Save current value */
      current = cm2_result;

      /* Then tree value is just original lagrangian mass^2: */
      cm2_result = m2_gluino_saf;

      /* These SE functions get evaluated with s = m2_gluino, which is
	 now the re part of the current pole mass: */
      pi1 = pi1tilde_gluino ();
      dpi1 = dpi1tildedM2_gluino();
      pi2 = pi2tilde_gluino ();

      cm2_result += pi1;
      cm2_result += -TSIL_CREAL(pi1)*dpi1;
      cm2_result += pi2;

      /* Note that the squark self energy functions do not include
	 SUMO_oneloopfactor!  Should standardize this eventually... */
      for (i=0; i<2; i++)
        cm2_result += SUMO_oneloopfactor * (
         - TSIL_CREAL(pi1_stop (i, i, m2_stop[i])) * dpi1tildedM2_stop (i)
      	 - TSIL_CREAL(pi1_sbot (i, i, m2_sbot[i])) * dpi1tildedM2_sbot (i)
      	 - TSIL_CREAL(pi1_suL (i, m2_suL[i])) * dpi1tildedM2_suL (i)
      	 - TSIL_CREAL(pi1_suR (i, m2_suR[i])) * dpi1tildedM2_suR (i)
      	 - TSIL_CREAL(pi1_sdL (i, m2_sdL[i])) * dpi1tildedM2_sdL (i)
      	 - TSIL_CREAL(pi1_sdR (i, m2_sdR[i])) * dpi1tildedM2_sdR (i)
      					    );

      cm2_result += -TSIL_CREAL(pi1tilde_top()) * dpi1tildedM2_top ();

      /* For next iteration, set tree mass to real part of current pole mass. */
      m2_gluino = TSIL_CREAL (cm2_result);
      m_gluino = TSIL_SQRT (m2_gluino);

      if (TSIL_CABS(cm2_result - current)/TSIL_CABS(current) < ITERTOL) {
	retval = iter+1;
	break;
      }
    }
  }

  /* Set the final values for the gluino mass: */
  CM2_gluino = cm2_result;
  M2_gluino = TSIL_CREAL(cm2_result);
  M_gluino = TSIL_SQRT(M2_gluino);

  /* Put back correct tree values if we changed them: */
  if (expand_around_pole != NONE) {
    m2_gluino = m2_gluino_saf;
    m_gluino = TSIL_SQRT(m2_gluino);
  }
  if (expand_around_pole == FULL) {
    for (i=0; i<2; i++) {
      m2_stop[i] = m2_stop_saf[i];
      m2_sbot[i] = m2_sbot_saf[i];
      m2_suL[i] = m2_suL_saf[i];
      m2_suR[i] = m2_suR_saf[i];
      m2_sdL[i] = m2_sdL_saf[i];
      m2_sdR[i] = m2_sdR_saf[i];
    }
    m2_top = m2_top_saf;
    m_top = m_top_saf;
  }

  return retval;
}
