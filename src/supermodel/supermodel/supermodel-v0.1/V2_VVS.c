#include "supermodel.h"

/* hep-ph/0206136 eq. (3.60) */

SUMO_COMPLEX SUMO_V2_VVS (void)
{ 
 int i;
 TSIL_COMPLEX V2_VVS = 0.0L;
 TSIL_REAL qq = Q2;

 for (i=0; i<2; i++){
  V2_VVS += (gWgammap[i]) * (gWgammap[i]) * 
    SUMO_Fvac_VVS (m2_W, 0.L, m2_phip[i], qq);

  V2_VVS += (gWZp[i]) * (gWZp[i]) * 
    SUMO_Fvac_VVS (m2_W, m2_Z, m2_phip[i], qq);

  V2_VVS += 0.25L * (gZZ0[i]) * (gZZ0[i]) * 
    SUMO_Fvac_VVS (m2_Z, m2_Z, m2_phi0[i], qq);

  V2_VVS += 0.5L * (gWW0[i]) * (gWW0[i]) * 
    SUMO_Fvac_VVS (m2_W, m2_W, m2_phi0[i], qq);
 }

 return (SUMO_twoloopfactor * V2_VVS);
}  
