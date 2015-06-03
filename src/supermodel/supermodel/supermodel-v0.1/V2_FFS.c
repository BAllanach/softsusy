#include "supermodel.h"

SUMO_COMPLEX SUMO_V2_FFS (void)
{ 
  int i,j,k;
  TSIL_COMPLEX V2_qNsq = 0.0L;
  TSIL_COMPLEX V2_lNsl = 0.0L;
  TSIL_COMPLEX V2_qCsq = 0.0L;
  TSIL_COMPLEX V2_lCsl = 0.0L;
  TSIL_COMPLEX V2_ff0 = 0.0L;
  TSIL_COMPLEX V2_ffp = 0.0L;
  TSIL_COMPLEX V2_CC0 = 0.0L;
  TSIL_COMPLEX V2_NN0 = 0.0L;
  TSIL_COMPLEX V2_CNp = 0.0L;


/* Note gluino contribution (eq. 3.25) is in V2_strong.c */

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.26) */

  for (i=0; i<4; i++) {
  for (j=0; j<2; j++) {
    V2_qNsq += (YtcNstop[i][j] * conYtcNstop[i][j]
              + YtNstopc[i][j] * conYtNstopc[i][j])
              * SUMO_Fvac_FFS (m2_top, m2_Neut[i], m2_stop[j], Q2);

    V2_qNsq += TSIL_CREAL (YtcNstop[i][j] * YtNstopc[i][j])
      * 2.0L * SUMO_Fvac_ffS (m2_top, m2_Neut[i], m2_stop[j], Q2);

    V2_qNsq += (YbcNsbot[i][j] * conYbcNsbot[i][j]
              + YbNsbotc[i][j] * conYbNsbotc[i][j])
              * SUMO_Fvac_FFS (m2_bot, m2_Neut[i], m2_sbot[j], Q2);

    V2_qNsq += TSIL_CREAL (YbcNsbot[i][j] * YbNsbotc[i][j])
      * 2.0L * SUMO_Fvac_ffS (m2_bot, m2_Neut[i], m2_sbot[j], Q2);

    V2_qNsq += YuNsuLc[i] * conYuNsuLc[i] *
      SUMO_Fvac_FFS (0.L, m2_Neut[i], m2_suL[j], Q2);

    V2_qNsq += YdNsdLc[i] * conYdNsdLc[i] *
      SUMO_Fvac_FFS (0.L, m2_Neut[i], m2_sdL[j], Q2);

    V2_qNsq += YucNsuR[i] * conYucNsuR[i] *
      SUMO_Fvac_FFS (0.L, m2_Neut[i], m2_suR[j], Q2);

    V2_qNsq += YdcNsdR[i] * conYdcNsdR[i] *
      SUMO_Fvac_FFS (0.L, m2_Neut[i], m2_sdR[j], Q2);
   }}

  V2_qNsq *= 3.0L;

  for (i=0; i<4; i++) {
  for (j=0; j<2; j++) {
    V2_lNsl += (YtaucNstau[i][j] * conYtaucNstau[i][j]
              + YtauNstauc[i][j] * conYtauNstauc[i][j])
              * SUMO_Fvac_FFS (m2_tau, m2_Neut[i], m2_stau[j], Q2);

    V2_lNsl += TSIL_CREAL (YtaucNstau[i][j] * YtauNstauc[i][j])
      * 2.0L * SUMO_Fvac_ffS (m2_tau, m2_Neut[i], m2_stau[j], Q2);

    V2_lNsl += YeNseLc[i] * conYeNseLc[i] *
      SUMO_Fvac_FFS (0.L, m2_Neut[i], m2_seL[j], Q2);

    V2_lNsl += YecNseR[i] * conYecNseR[i] *
      SUMO_Fvac_FFS (0.L, m2_Neut[i], m2_seR[j], Q2);
  }}

  for (i=0; i<4; i++) {
   for (j=0; j<3; j++) {
    V2_lNsl += YnuNsnuc[i] * conYnuNsnuc[i] *
      SUMO_Fvac_FFS (0.L, m2_Neut[i], m2_snu[j], Q2);
  }}

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.34) */

  for (i=0; i<2; i++) {
   for (j=0; j<2; j++) {
    V2_qCsq += (YtcCsbot[i][j] * conYtcCsbot[i][j]
               +YtCsbotc[i][j] * conYtCsbotc[i][j]) *
               SUMO_Fvac_FFS (
               m2_top, m2_Char[i], m2_sbot[j], Q2);

    V2_qCsq += TSIL_CREAL (YtcCsbot[i][j] * YtCsbotc[i][j]) *
               2.0L * SUMO_Fvac_ffS (
               m2_top, m2_Char[i], m2_sbot[j], Q2);

    V2_qCsq += (YbcCstop[i][j] * conYbcCstop[i][j]
               +YbCstopc[i][j] * conYbCstopc[i][j]) *
               SUMO_Fvac_FFS (
               m2_bot, m2_Char[i], m2_stop[j], Q2);

    V2_qCsq += TSIL_CREAL (YbcCstop[i][j] * YbCstopc[i][j]) *
               2.0L * SUMO_Fvac_ffS (
               m2_bot, m2_Char[i], m2_stop[j], Q2);

    /* Following two terms appear to have no counterparts in Pi2, so
       leave them commented out when comparing to 4.7,9,11 */

    V2_qCsq += YdCsuLc[i] * conYdCsuLc[i] *
      SUMO_Fvac_FFS (0.L, m2_Char[i], m2_suL[j], Q2);

    V2_qCsq += YuCsdLc[i] * conYuCsdLc[i] *
      SUMO_Fvac_FFS (0.L, m2_Char[i], m2_sdL[j], Q2);
  }}

  V2_qCsq *= 3.0L;

  for (i=0; i<2; i++) {
   for (j=0; j<2; j++) {
    V2_lCsl += YnutauCstauc[i][j] * conYnutauCstauc[i][j] *
      SUMO_Fvac_FFS (0.L, m2_Char[i], m2_stau[j], Q2);

    /* Following two terms appear to have no counterparts in Pi2, so
       leave them commented out when comparing to 4.7,9,11 */

    V2_lCsl += YeCsnuc[i] * conYeCsnuc[i] *
      SUMO_Fvac_FFS (0.L, m2_Char[i], m2_snu[j], Q2);

    V2_lCsl += YnuCseLc[i] * conYnuCseLc[i] *
      SUMO_Fvac_FFS (0.L, m2_Char[i], m2_seL[j], Q2);
  }}

  for (i=0; i<2; i++) {
    V2_lCsl += (YtaucCsnutau[i] * conYtaucCsnutau[i]
               +YtauCsnutauc[i] * conYtauCsnutauc[i]) *
               SUMO_Fvac_FFS (
               m2_tau, m2_Char[i], m2_snu[2], Q2);

    V2_lCsl += TSIL_CREAL (YtaucCsnutau[i] * YtauCsnutauc[i]) *
               2.0L * SUMO_Fvac_ffS (
               m2_tau, m2_Char[i], m2_snu[2], Q2);
  }

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.39) */

  for (i=0; i<4; i++) {
    V2_ff0 += 1.5L * ytop2 * (
              ku[i] * kuc[i] * SUMO_Fvac_FFS (
                         m2_top, m2_top, m2_phi0[i], Q2)
            + ku[i] * ku[i] * SUMO_Fvac_ffS (
                         m2_top, m2_top, m2_phi0[i], Q2));

    V2_ff0 += 1.5L * ybot2 * (
              kd[i] * kdc[i] * SUMO_Fvac_FFS (
                          m2_bot, m2_bot, m2_phi0[i], Q2)
            + kd[i] * kd[i] * SUMO_Fvac_ffS (
                          m2_bot, m2_bot, m2_phi0[i], Q2));

    V2_ff0 += 0.5L * ytau2 * (
              kd[i] * kdc[i] * SUMO_Fvac_FFS (
                          m2_tau, m2_tau, m2_phi0[i], Q2)
            + kd[i] * kd[i] * SUMO_Fvac_ffS (
                          m2_tau, m2_tau, m2_phi0[i], Q2));
  }

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.40) */
/* Checked and correspond with terms in 4.13 of Pi2nonQCD (for
   g=gp=0) */

  for (i=0; i<2; i++) {
    V2_ffp += 3.0L * (ytop2 * kup[i] * kup[i] + ybot2 * kdp[i] * kdp[i])
      * SUMO_Fvac_FFS (m2_top, m2_bot, m2_phip[i], Q2);

    V2_ffp += 6.0L * ytop * ybot * kup[i] * kdp[i]
      * SUMO_Fvac_ffS (m2_top, m2_bot, m2_phip[i], Q2);

    V2_ffp += ytau2 * kdp[i] * kdp[i]
      * SUMO_Fvac_FFS (0.0L, m2_tau, m2_phip[i], Q2);
  }

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.41) */

  for (k=0; k<4; k++) {
   for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
     V2_CC0 += YCC0[i][j][k] * conYCC0[i][j][k] *
               SUMO_Fvac_FFS (
               m2_Char[i], m2_Char[j], m2_phi0[k], Q2);

     V2_CC0 += TSIL_CREAL(YCC0[i][j][k] * YCC0[j][i][k]) *
       SUMO_Fvac_ffS (
       m2_Char[i], m2_Char[j], m2_phi0[k], Q2);
  }}}

/* ------------------------------------------------------------------- */
/* hep-ph/0206136 eq. (3.42) */

  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) {
     V2_NN0 += YNN0[i][i][k] * conYNN0[i][i][k] *
               SUMO_Fvac_FFS (
               m2_Neut[i], m2_Neut[i], m2_phi0[k], Q2);

     V2_NN0 += TSIL_CREAL(YNN0[i][i][k] * YNN0[i][i][k]) *
       SUMO_Fvac_ffS (
       m2_Neut[i], m2_Neut[i], m2_phi0[k], Q2);
  }}

  V2_NN0 = 0.5L * V2_NN0;

  for (i=1; i<4; i++) {
   for (j=0; j<i; j++) {
    for (k=0; k<4; k++) {
     V2_NN0 += YNN0[i][j][k] * conYNN0[i][j][k] *
               SUMO_Fvac_FFS (
               m2_Neut[i], m2_Neut[j], m2_phi0[k], Q2);

     V2_NN0 += TSIL_CREAL(YNN0[i][j][k] * YNN0[i][j][k]) *
       SUMO_Fvac_ffS (
       m2_Neut[i], m2_Neut[j], m2_phi0[k], Q2);
  }}}


/* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.43) */

  for (i=0; i<2; i++) {
   for (j=0; j<4; j++) {
    for (k=0; k<2; k++) {
     V2_CNp += (YCNp[i][j][k] * conYCNp[i][j][k] +
                YCNm[i][j][k] * conYCNm[i][j][k]) *
                SUMO_Fvac_FFS (
                m2_Char[i], m2_Neut[j], m2_phip[k], Q2);

     V2_CNp += 2.0L * TSIL_CREAL(YCNp[i][j][k] * YCNm[i][j][k]) *
       SUMO_Fvac_ffS (
       m2_Char[i], m2_Neut[j], m2_phip[k], Q2);
  }}}

/* ------------------------------------------------------------------- */

  return (SUMO_twoloopfactor * (
         V2_qNsq
       + V2_lNsl
       + V2_qCsq
       + V2_lCsl
       + V2_ff0
       + V2_ffp
       + V2_CC0
       + V2_NN0
       + V2_CNp
       ));
}
