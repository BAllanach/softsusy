/*   
   Computes pole mass of gluino at one or two loops, using
   eqs. (5.1)-(5.3) of hep-ph/0509115 
*/

#include "supermodel.h"

/* Self energy functions defined in se_gluino.c: */
SUMO_COMPLEX pi1tilde_gluino (int expand_around_pole);
SUMO_COMPLEX pi2tilde_gluino (int expand_around_pole);

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

int SUMO_GluinoPole (int loop_order)
{
  int expand_around_pole = 0; /* Later will hopefully add option 1. */
  TSIL_COMPLEX cm2_result;

  if (0 == are_tree_couplings_updated) SUMO_Tree_Couplings ();

  cm2_result = m2_gluino;

  if (loop_order > 0)
    cm2_result += pi1tilde_gluino (expand_around_pole);

  if (loop_order > 1)
    cm2_result += pi2tilde_gluino (expand_around_pole);

  CM2_gluino = cm2_result;
  M2_gluino = TSIL_CREAL(cm2_result);
  M_gluino = TSIL_SQRT(M2_gluino);

  return 0;
}
