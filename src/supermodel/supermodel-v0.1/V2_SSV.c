#include "supermodel.h"

SUMO_COMPLEX SUMO_V2_SSV (void)
{ 
  int i, j;
  TSIL_COMPLEX V2_sfsfgamma = 0.0L;
  TSIL_COMPLEX V2_sfsfZ = 0.0L;
  TSIL_COMPLEX V2_sfsfW = 0.0L;
  TSIL_COMPLEX V2_pmgamma = 0.0L;
  TSIL_COMPLEX V2_pmZ = 0.0L;
  TSIL_COMPLEX V2_00Z = 0.0L;
  TSIL_COMPLEX V2_0pW = 0.0L;

/* -------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.49) */

  for (i=0; i<2; i++){  
    V2_sfsfgamma += (4.0L/3.0L) * (
        SUMO_Fvac_SSV (m2_suR[i], m2_suR[i], 0.L, Q2)
      + SUMO_Fvac_SSV (m2_suL[i], m2_suL[i], 0.L, Q2)
      + SUMO_Fvac_SSV (m2_stop[i], m2_stop[i], 0.L, Q2));

    V2_sfsfgamma += (1.0L/3.0L) * (
        SUMO_Fvac_SSV (m2_sdR[i], m2_sdR[i], 0.L, Q2)
      + SUMO_Fvac_SSV (m2_sdL[i], m2_sdL[i], 0.L, Q2)
      + SUMO_Fvac_SSV (m2_sbot[i], m2_sbot[i], 0.L, Q2));

    V2_sfsfgamma += 
        SUMO_Fvac_SSV (m2_seR[i], m2_seR[i], 0.L, Q2)
      + SUMO_Fvac_SSV (m2_seL[i], m2_seL[i], 0.L, Q2)
      + SUMO_Fvac_SSV (m2_stau[i], m2_stau[i], 0.L, Q2);
  }

  V2_sfsfgamma *= 0.5L * e2;

/* -------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.50) */

  for (i=0; i<2; i++){  
    V2_sfsfZ += gZsuLsuLc * gZsuLsuLc * 
                SUMO_Fvac_SSV (m2_suL[i], m2_suL[i], m2_Z, Q2);

    V2_sfsfZ += gZsdLsdLc * gZsdLsdLc *
                SUMO_Fvac_SSV (m2_sdL[i], m2_sdL[i], m2_Z, Q2);

    V2_sfsfZ += gZsuRsuRc * gZsuRsuRc * 
                SUMO_Fvac_SSV (m2_suR[i], m2_suR[i], m2_Z, Q2);

    V2_sfsfZ += gZsdRsdRc * gZsdRsdRc *
                SUMO_Fvac_SSV (m2_sdR[i], m2_sdR[i], m2_Z, Q2);

    for (j=0; j<2; j++){  
       V2_sfsfZ += gZstopstopc[i][j] * gZstopstopc[j][i] * 
                   SUMO_Fvac_SSV (m2_stop[i], m2_stop[j], m2_Z, Q2);

       V2_sfsfZ += gZsbotsbotc[i][j] * gZsbotsbotc[j][i] * 
                   SUMO_Fvac_SSV (m2_sbot[i], m2_sbot[j], m2_Z, Q2);
    }
  }

  V2_sfsfZ *= 3.0L;

  for (i=0; i<2; i++){  
    V2_sfsfZ += gZsnusnuc * gZsnusnuc *
                SUMO_Fvac_SSV (m2_snu[i], m2_snu[i], m2_Z, Q2);

    V2_sfsfZ += gZseLseLc * gZseLseLc *
                SUMO_Fvac_SSV (m2_seL[i], m2_seL[i], m2_Z, Q2);

    V2_sfsfZ += gZseRseRc * gZseRseRc *
                SUMO_Fvac_SSV (m2_seR[i], m2_seR[i], m2_Z, Q2);

    for (j=0; j<2; j++){  
       V2_sfsfZ += gZstaustauc[i][j] * gZstaustauc[j][i] * 
                   SUMO_Fvac_SSV (m2_stau[i], m2_stau[j], m2_Z, Q2);
    }
  }

  V2_sfsfZ += gZsnusnuc * gZsnusnuc *
              SUMO_Fvac_SSV (m2_snu[2], m2_snu[2], m2_Z, Q2);

  V2_sfsfZ *= 0.5L;

/* -------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.51) */

  for (i=0; i<2; i++){  
    V2_sfsfW += 3.0L * SUMO_Fvac_SSV (m2_suL[i], m2_sdL[i], m2_W, Q2);

    V2_sfsfW += SUMO_Fvac_SSV (m2_snu[i], m2_seL[i], m2_W, Q2);

    V2_sfsfW += Lstau[i] * Lstauc[i] * SUMO_Fvac_SSV (
                m2_snu[2], m2_stau[i], m2_W, Q2);

    for (j=0; j<2; j++){  
      V2_sfsfW += 3.0L * Lstop[i] * Lstopc[i] *
                  Lsbot[j] * Lsbotc[j] * SUMO_Fvac_SSV (
                  m2_stop[i], m2_sbot[j], m2_W, Q2);
    }
  }

  V2_sfsfW *= 0.5L * g2;

/* -------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.52), (3.53) */

  for (i=0; i<2; i++){  
    V2_pmgamma +=  
      SUMO_Fvac_SSV (m2_phip[i], m2_phip[i], 0.L, Q2);

    V2_pmZ +=  
      SUMO_Fvac_SSV (m2_phip[i], m2_phip[i], m2_Z, Q2);
  }
  
  V2_pmgamma *= 0.5L * e2;
  V2_pmZ *= 0.125L*(g2 - gp2)*(g2 - gp2)/(g2plusgp2);

/* -------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.54) */

  for (i=0; i<2; i++){  
    for (j=2; j<4; j++){  
      V2_00Z += 0.5L * gZ00[i][j] * gZ00[i][j] *  
                SUMO_Fvac_SSV (m2_phi0[i], m2_phi0[j], m2_Z, Q2);
  }}

/* -------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.55) */

  for (i=0; i<4; i++){  
    for (j=0; j<2; j++){  
      V2_0pW += gW0p[i][j] * SUMO_CONJ(gW0p[i][j]) *  
                SUMO_Fvac_SSV (m2_phi0[i], m2_phip[j], m2_W, Q2);
  }}

/* -------------------------------------------------------------------- */

  return (SUMO_twoloopfactor * (V2_sfsfgamma + V2_sfsfZ 
                + V2_sfsfW + V2_pmgamma + V2_pmZ + V2_00Z + V2_0pW));
}  
