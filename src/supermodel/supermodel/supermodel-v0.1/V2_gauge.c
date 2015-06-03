#include "supermodel.h"

/* Implements hep-ph/0206136 eq. (3.75) */

SUMO_COMPLEX SUMO_V2_gauge (void)
{ 
  TSIL_REAL qq = Q2;

  return (0.5L * SUMO_twoloopfactor 
    * TSIL_CREAL( e2 * SUMO_Fvac_gauge (m2_W, m2_W, 0.L, qq)  
    + (g2 * g2/g2plusgp2) * SUMO_Fvac_gauge (m2_W, m2_W, m2_Z, qq) ));
}  
