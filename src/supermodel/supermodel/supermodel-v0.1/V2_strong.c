#include "supermodel.h"

SUMO_COMPLEX SUMO_V2_strong (void)
{ 
  int i,j;
  TSIL_REAL g3sq = g3*g3;
  TSIL_REAL qq = Q2;
  TSIL_COMPLEX V2_strong_sqsq = 0.0L;
  TSIL_COMPLEX V2_strong_qgluinosq = 0.0L;
  TSIL_COMPLEX V2_strong_sqsqg = 0.0L;
  TSIL_COMPLEX V2_strong_qqg;
  
  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.14) */

  for (i=0; i<2; i++) {
    V2_strong_sqsq += SUMO_Fvac_SS (m2_suL[i], m2_suL[i], qq);
    V2_strong_sqsq += SUMO_Fvac_SS (m2_sdL[i], m2_sdL[i], qq);
    V2_strong_sqsq += SUMO_Fvac_SS (m2_suR[i], m2_suR[i], qq);
    V2_strong_sqsq += SUMO_Fvac_SS (m2_sdR[i], m2_sdR[i], qq);
  }

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      V2_strong_sqsq +=
	SUMO_Fvac_SS (m2_stop[i], m2_stop[j], qq) *
	TSIL_POW(TSIL_CABS((Lstop[i]) * (Lstopc[j])
			   -(Rstop[i]) * (Rstopc[j])),2);

      V2_strong_sqsq +=
	SUMO_Fvac_SS (m2_sbot[i], m2_sbot[j], qq) *
	TSIL_POW(TSIL_CABS((Lsbot[i]) * (Lsbotc[j])
			   -(Rsbot[i]) * (Rsbotc[j])),2);
    }}

  V2_strong_sqsq *= 2.0L * g3sq;

  /* ------------------------------------------------------------------- */
  /* equation (3.25) */

  for (i=0; i<2; i++) {
    V2_strong_qgluinosq +=
      SUMO_Fvac_FFS (0.L, m2_gluino, m2_suL[i], qq);
    V2_strong_qgluinosq +=
      SUMO_Fvac_FFS (0.L, m2_gluino, m2_sdL[i], qq);
    V2_strong_qgluinosq +=
      SUMO_Fvac_FFS (0.L, m2_gluino, m2_suR[i], qq);
    V2_strong_qgluinosq +=
      SUMO_Fvac_FFS (0.L, m2_gluino, m2_sdR[i], qq);
    V2_strong_qgluinosq +=
      SUMO_Fvac_FFS (m2_bot, m2_gluino, m2_sbot[i], qq);
    V2_strong_qgluinosq +=
      SUMO_Fvac_FFS (m2_top, m2_gluino, m2_stop[i], qq);
    V2_strong_qgluinosq +=
      - 2.0L * TSIL_CREAL((Lsbot[i]) * (Rsbotc[i])) *
      SUMO_Fvac_ffS (m2_bot, m2_gluino, m2_sbot[i], qq);
    V2_strong_qgluinosq +=
      - 2.0L * TSIL_CREAL((Lstop[i]) * (Rstopc[i])) *
      SUMO_Fvac_ffS (m2_top, m2_gluino, m2_stop[i], qq);
  }

  V2_strong_qgluinosq *= 8.0L*g3sq;

  /* ------------------------------------------------------------------- */
  /* equation (3.48) */

  for (i=0; i<2; i++) {
    V2_strong_sqsqg +=
      SUMO_Fvac_SSV (m2_suL[i], m2_suL[i], 0.0L, qq);
    V2_strong_sqsqg +=
      SUMO_Fvac_SSV (m2_sdL[i], m2_sdL[i], 0.0L, qq);
    V2_strong_sqsqg +=
      SUMO_Fvac_SSV (m2_suR[i], m2_suR[i], 0.0L, qq);
    V2_strong_sqsqg +=
      SUMO_Fvac_SSV (m2_sdR[i], m2_sdR[i], 0.0L, qq);
    V2_strong_sqsqg +=
      SUMO_Fvac_SSV (m2_stop[i], m2_stop[i], 0.0L, qq);
    V2_strong_sqsqg +=
      SUMO_Fvac_SSV (m2_sbot[i], m2_sbot[i], 0.0L, qq);
  }

  V2_strong_sqsqg *= 2.0L*g3sq;

  /* ------------------------------------------------------------------- */
  /* equation (3.61) */

  V2_strong_qqg = 4.0L*g3sq*(
      SUMO_Fvac_FFV (m2_top, m2_top, 0.0L, qq)
    + SUMO_Fvac_FFV (m2_bot, m2_bot, 0.0L, qq)
    - SUMO_Fvac_ffV (m2_top, m2_top, 0.0L, qq)
    - SUMO_Fvac_ffV (m2_bot, m2_bot, 0.0L, qq)
      );

  /* ------------------------------------------------------------------- */

  return (SUMO_twoloopfactor * 
    (V2_strong_sqsq + V2_strong_qgluinosq + V2_strong_sqsqg + V2_strong_qqg));
}
