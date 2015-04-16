#include "supermodel.h"

SUMO_COMPLEX SUMO_V2_SSS (void)
{ 
  int i,j,k;
  TSIL_REAL qq = Q2;
  TSIL_COMPLEX V2_000 = 0.0L;
  TSIL_COMPLEX V2_0pm = 0.0L;
  TSIL_COMPLEX V2_0sfsf = 0.0L;
  TSIL_COMPLEX V2_psfsf = 0.0L;

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.3) */

  for (i=0; i<2; i++) {
    V2_000 += (lambda000[i][i][i]) * (lambda000[i][i][i]) *
              SUMO_Fvac_SSS (
              m2_phi0[i], m2_phi0[i], m2_phi0[i], qq)/12.0L;

    V2_000 += (lambda000[i][2][2]) * (lambda000[i][2][2]) *
             SUMO_Fvac_SSS (
             m2_phi0[i], m2_phi0[2], m2_phi0[2], qq)/4.0L;

    V2_000 += (lambda000[i][3][3]) * (lambda000[i][3][3]) *
             SUMO_Fvac_SSS (
             m2_phi0[i], m2_phi0[3], m2_phi0[3], qq)/4.0L;

    V2_000 += (lambda000[i][2][3]) * (lambda000[i][2][3]) *
             SUMO_Fvac_SSS (
             m2_phi0[i], m2_phi0[2], m2_phi0[3], qq)/2.0L;
  }

  V2_000 += (lambda000[0][1][1]) * (lambda000[0][1][1]) *
            SUMO_Fvac_SSS (
            m2_phi0[0], m2_phi0[1], m2_phi0[1], qq)/4.0L;

  V2_000 += (lambda000[0][0][1]) * (lambda000[0][0][1]) *
            SUMO_Fvac_SSS (
            m2_phi0[0], m2_phi0[0], m2_phi0[1], qq)/4.0L;

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.4) */

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      V2_0pm += 0.5L * lambda0pm[i][j][j] * lambda0pm[i][j][j] *
                SUMO_Fvac_SSS (
                m2_phi0[i], m2_phip[j], m2_phip[j], qq);
    }}

  for (i=0; i<4; i++) {
    V2_0pm += lambda0pm[i][0][1] * lambda0pm[i][1][0] *
              SUMO_Fvac_SSS (
              m2_phi0[i], m2_phip[0], m2_phip[1], qq);
  }

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.5) */

  for (i=0; i<4; i++) {

    for (j=0; j<2; j++) {
      V2_0sfsf += 1.5L * lambda0suLsuLc[i] * lambda0suLsuLc[i] *
                 SUMO_Fvac_SSS (
                 m2_phi0[i], m2_suL[j], m2_suL[j], qq);

      V2_0sfsf += 1.5L * lambda0sdLsdLc[i] * lambda0sdLsdLc[i] *
                 SUMO_Fvac_SSS (
                 m2_phi0[i], m2_sdL[j], m2_sdL[j], qq);

      V2_0sfsf += 1.5L * lambda0suRsuRc[i] * lambda0suRsuRc[i] *
                 SUMO_Fvac_SSS (
                 m2_phi0[i], m2_suR[j], m2_suR[j], qq);

      V2_0sfsf += 1.5L * lambda0sdRsdRc[i] * lambda0sdRsdRc[i] *
                 SUMO_Fvac_SSS (
                 m2_phi0[i], m2_sdR[j], m2_sdR[j], qq);

      V2_0sfsf += 0.5L * lambda0seLseLc[i] * lambda0seLseLc[i] *
                 SUMO_Fvac_SSS (
                 m2_phi0[i], m2_seL[j], m2_seL[j], qq);

      V2_0sfsf += 0.5L * lambda0seRseRc[i] * lambda0seRseRc[i] *
                 SUMO_Fvac_SSS (
                 m2_phi0[i], m2_seR[j], m2_seR[j], qq);
    }

    for (j=0; j<3; j++) {
      V2_0sfsf += 0.5L * lambda0snusnuc[i] * lambda0snusnuc[i] *
                 SUMO_Fvac_SSS (
                 m2_phi0[i], m2_snu[j], m2_snu[j], qq);
    }

    V2_0sfsf += 3.0L
      * lambda0stopstopc[i][0][1] * lambda0stopstopc[i][1][0]
      * SUMO_Fvac_SSS (m2_phi0[i], m2_stop[0], m2_stop[1], qq);

    V2_0sfsf += 3.0L
      * lambda0sbotsbotc[i][0][1] * lambda0sbotsbotc[i][1][0]
      * SUMO_Fvac_SSS (m2_phi0[i], m2_sbot[0], m2_sbot[1], qq);

    V2_0sfsf += lambda0staustauc[i][0][1] * lambda0staustauc[i][1][0]
      * SUMO_Fvac_SSS (m2_phi0[i], m2_stau[0], m2_stau[1], qq);

    for (j=0; j<2; j++) {

      V2_0sfsf += 1.5L
	* lambda0stopstopc[i][j][j] * lambda0stopstopc[i][j][j]
	* SUMO_Fvac_SSS (m2_phi0[i], m2_stop[j], m2_stop[j], qq);

      V2_0sfsf += 1.5L
	* lambda0sbotsbotc[i][j][j] * lambda0sbotsbotc[i][j][j]
	* SUMO_Fvac_SSS (m2_phi0[i], m2_sbot[j], m2_sbot[j], qq);

      V2_0sfsf += 0.5L
	* lambda0staustauc[i][j][j] * lambda0staustauc[i][j][j]
	* SUMO_Fvac_SSS (m2_phi0[i], m2_stau[j], m2_stau[j], qq);
    }
  }

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.10) */

  for (i=0; i<2; i++){

    for (j=0; j<2; j++){
      V2_psfsf += 3.0L * lambdapsdLsuLc[i] * lambdapsdLsuLc[i]
	* SUMO_Fvac_SSS (m2_phip[i], m2_sdL[j], m2_suL[j], qq);

      V2_psfsf += lambdapseLsnuc[i] * lambdapseLsnuc[i] *
                  SUMO_Fvac_SSS (
                  m2_phip[i], m2_seL[j], m2_snu[j], qq);

      V2_psfsf += SUMO_CONJ(lambdapstausnutauc[i][j]) *
                  lambdapstausnutauc[i][j] * SUMO_Fvac_SSS (
                  m2_phip[i], m2_stau[j], m2_snu[2], qq);
    }

    for (j=0; j<2; j++) {
      for (k=0; k<2; k++) {
	V2_psfsf += 3.0L
	  * SUMO_CONJ(lambdapsbotstopc[i][j][k])
	            * lambdapsbotstopc[i][j][k]
	  * SUMO_Fvac_SSS (m2_phip[i], m2_sbot[j], m2_stop[k], qq);
      }}
  }

  /* ------------------------------------------------------------------- */

  return (SUMO_twoloopfactor * (V2_000 + V2_0pm + V2_0sfsf + V2_psfsf) );
}
