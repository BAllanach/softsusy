#include "supermodel.h"

SUMO_COMPLEX SUMO_V2_FFV (void)
{ 
  int i,j;
  TSIL_REAL qq = Q2;  
  TSIL_COMPLEX V2_ffgamma =0.0L;
  TSIL_COMPLEX V2_ffZ = 0.0L;
  TSIL_COMPLEX V2_ffW = 0.0L;
  TSIL_COMPLEX V2_CCgamma = 0.0L;
  TSIL_COMPLEX V2_CCZ = 0.0L;
  TSIL_COMPLEX V2_NNZ = 0.0L;
  TSIL_COMPLEX V2_NCW = 0.0L;

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.62) */

  V2_ffgamma = (e2) * (
    (4.0L/3.0L) * (SUMO_Fvac_FFV (m2_top, m2_top, 0.0L, qq)
      - SUMO_Fvac_ffV (m2_top, m2_top, 0.0L, qq))
    + (1.0L/3.0L) * (SUMO_Fvac_FFV (m2_bot, m2_bot, 0.0L, qq)
      - SUMO_Fvac_ffV (m2_bot, m2_bot, 0.0L, qq))
    + SUMO_Fvac_FFV (m2_tau, m2_tau, 0.0L, qq)
      - SUMO_Fvac_ffV (m2_tau, m2_tau, 0.0L, qq)
			    );

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.63) */

  V2_ffZ = (1.0L/(24.0L * (g2plusgp2))) * (
      (51.0L * g2 * g2 + 6.0L * g2 * gp2 + 83.0L * gp2 * gp2) *
         SUMO_Fvac_FFV (0.0L, 0.0L, m2_Z, qq)
    + (9.0L * g2 * g2 - 6.0L * g2 * gp2 + 17.0L * gp2 * gp2) *
      SUMO_Fvac_FFV (m2_top, m2_top, m2_Z, qq)
    + gp2 * (24.0L * g2 - 8.0L * gp2) *
      SUMO_Fvac_ffV (m2_top, m2_top, m2_Z, qq)
    + (9.0L * g2 * g2 + 6.0L * g2 * gp2 + 5.0L * gp2 * gp2) *
      SUMO_Fvac_FFV (m2_bot, m2_bot, m2_Z, qq)
    + gp2 * (12.0L * g2 + 4.0L * gp2) *
      SUMO_Fvac_ffV (m2_bot, m2_bot, m2_Z, qq)
    + (3.0L * g2 * g2 - 6.0L * g2 * gp2 + 15.0L * gp2 * gp2) *
      SUMO_Fvac_FFV (m2_tau, m2_tau, m2_Z, qq)
    + 12.0L * gp2 * (g2 - gp2) *
      SUMO_Fvac_ffV (m2_tau, m2_tau, m2_Z, qq)
						);

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.64) */

  V2_ffW = 0.5L * g2 * (
      3.0L * SUMO_Fvac_FFV (m2_top, m2_bot, m2_W, qq)
      + SUMO_Fvac_FFV (m2_tau, 0.0L, m2_W, qq)
      + 8.0L * SUMO_Fvac_FFV (0.0L, 0.0L, m2_W, qq));

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.65) */

  for (i=0; i<2; i++) {
    V2_CCgamma += (
       SUMO_Fvac_FFV (m2_Char[i], m2_Char[i], 0.0L, qq)
       -SUMO_Fvac_ffV (m2_Char[i], m2_Char[i], 0.0L, qq));
  }
  
  V2_CCgamma *= (e2);

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.66) */

  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      V2_NNZ += TSIL_CREAL(
        TSIL_CABS(OppL[i][j] * OppL[i][j]) *
     SUMO_Fvac_FFV (m2_Neut[i], m2_Neut[j], m2_Z, qq)
        + OppL[i][j] * OppL[i][j] *
     SUMO_Fvac_ffV (m2_Neut[i], m2_Neut[j], m2_Z, qq));
  }}
  
  V2_NNZ *= 0.5L * (g2plusgp2);

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.67) */

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      V2_CCZ += (TSIL_CABS(OpL[i][j] * OpL[i][j])
               + TSIL_CABS(OpR[i][j] * OpR[i][j])) *
        SUMO_Fvac_FFV (m2_Char[i], m2_Char[j], m2_Z, qq)
        -2.0L * OpL[i][j] * OpR[j][i] *
        SUMO_Fvac_ffV (m2_Char[i], m2_Char[j], m2_Z, qq);
  }}
  
  V2_CCZ *= 0.5L * (g2plusgp2);

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.68) */

  for (i=0; i<4; i++) {
    for (j=0; j<2; j++) {
      V2_NCW += (TSIL_CABS(OL[i][j] * OL[i][j])
               + TSIL_CABS(OR[i][j] * OR[i][j])) *
       SUMO_Fvac_FFV (m2_Neut[i], m2_Char[j], m2_W, qq)
        -2.0L * TSIL_CREAL(OL[i][j] * ORc[i][j]) *
       SUMO_Fvac_ffV (m2_Neut[i], m2_Char[j], m2_W, qq);
  }}
  
  V2_NCW *= g2;

  return (SUMO_twoloopfactor * 
    (V2_ffgamma + V2_ffZ + V2_ffW + V2_CCgamma + V2_CCZ + V2_NNZ + V2_NCW));
}
