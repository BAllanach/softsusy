#include "supermodel.h"

SUMO_COMPLEX SUMO_V2_VS (void)
{ 
  int i;
  TSIL_COMPLEX V2_Wphi = 0.0L;
  TSIL_COMPLEX V2_Zphi = 0.0L;
  TSIL_COMPLEX V2_Wsq = 0.0L;
  TSIL_COMPLEX V2_Wsl = 0.0L;
  TSIL_COMPLEX V2_ZsuR = 0.0L;
  TSIL_COMPLEX V2_ZsuL = 0.0L;
  TSIL_COMPLEX V2_ZsdR = 0.0L;
  TSIL_COMPLEX V2_ZsdL = 0.0L;
  TSIL_COMPLEX V2_Zsnu = 0.0L;
  TSIL_COMPLEX V2_ZseR = 0.0L;
  TSIL_COMPLEX V2_ZseL = 0.0L;
  TSIL_COMPLEX V2_Zsf = 0.0L;
  TSIL_REAL qq = Q2;
  TSIL_COMPLEX tempW, tempZ;

/* -------------------------------------------------------------------- */
/* In the following, use F_VS(x,y) = 3 A(x) A(y) */

  tempW = 3.0L * TSIL_A (m2_W, qq);
  tempZ = 3.0L * TSIL_A (m2_Z, qq);

/* hep-ph/0206136 eq. (3.56), (3.57) */

  for (i=0; i<4; i++){
    V2_Wphi += 0.5L * gWW00[i] * TSIL_A (m2_phi0[i], qq);
    V2_Zphi += 0.25L * gZZ00[i] * TSIL_A (m2_phi0[i], qq);
  }

  for (i=0; i<2; i++){
    V2_Wphi += gWWpm[i] * TSIL_A (m2_phip[i], qq);
    V2_Zphi += 0.5L * gZZpm[i] * TSIL_A (m2_phip[i], qq);
  }

  V2_Wphi *= tempW;
  V2_Zphi *= tempZ;

/* -------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.58) */

  for (i=0; i<2; i++){
    V2_Wsq += TSIL_A (m2_suL[i], qq);
    V2_Wsq += TSIL_A (m2_sdL[i], qq);
    V2_Wsq += Lstop[i] * Lstopc[i] * TSIL_A (m2_stop[i], qq);
    V2_Wsq += Lsbot[i] * Lsbotc[i] * TSIL_A (m2_sbot[i], qq);
  }

  V2_Wsq *= 1.5L * g2 * tempW;

  for (i=0; i<2; i++){
    V2_Wsl += TSIL_A (m2_snu[i], qq);
    V2_Wsl += TSIL_A (m2_seL[i], qq);
    V2_Wsl += Lstau[i] * Lstauc[i] * TSIL_A (m2_stau[i], qq);
  }

  V2_Wsl += TSIL_A (m2_snu[2], qq);

  V2_Wsl *= 0.5L * g2 * tempW;

/* -------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.59) */

  for (i=0; i<2; i++){
    V2_ZsuR += TSIL_A (m2_suR[i], qq);
    V2_ZsuR += Rstop[i] * Rstopc[i] * TSIL_A (m2_stop[i], qq);

    V2_ZsuL += TSIL_A (m2_suL[i], qq);
    V2_ZsuL += Lstop[i] * Lstopc[i] * TSIL_A (m2_stop[i], qq);

    V2_ZsdR += TSIL_A (m2_sdR[i], qq);
    V2_ZsdR += Rsbot[i] * Rsbotc[i] * TSIL_A (m2_sbot[i], qq);

    V2_ZsdL += TSIL_A (m2_sdL[i], qq);
    V2_ZsdL += Lsbot[i] * Lsbotc[i] * TSIL_A (m2_sbot[i], qq);

    V2_Zsnu += TSIL_A (m2_snu[i], qq);

    V2_ZseR += TSIL_A (m2_seR[i], qq);
    V2_ZseR += Rstau[i] * Rstauc[i] * TSIL_A (m2_stau[i], qq);

    V2_ZseL += TSIL_A (m2_seL[i], qq);
    V2_ZseL += Lstau[i] * Lstauc[i] * TSIL_A (m2_stau[i], qq);
  }

  V2_Zsnu += TSIL_A (m2_snu[2], qq);
 
  V2_Zsf = (gp2 * gp2 * (16.0L * V2_ZsuR + 4.0L * V2_ZsdR + 12.0L * V2_ZseR)  
     + (3.0L * g2 - gp2) * (3.0L * g2 - gp2) * V2_ZsuL
     + (3.0L * g2 + gp2) * (3.0L * g2 + gp2) * V2_ZsdL
     + 3.0L * (g2 - gp2) * (g2 - gp2) * V2_ZseL
     + 3.0L * (g2plusgp2) * (g2plusgp2) * V2_Zsnu)/(12.0L * (g2plusgp2));

  V2_Zsf *= tempZ;

/* -------------------------------------------------------------------- */

  return (SUMO_twoloopfactor * (V2_Wphi + V2_Zphi + V2_Wsq + V2_Wsl + V2_Zsf));
}  
