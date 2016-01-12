/* 1-loop self energy functions for squarks */

#include "supermodel.h"
#include "self_scalar.h"

#define Cq (4.0L/3.0L)

/* Lower case (pi1...) use tree masses.
   Upper case (Pi1...) use pole masses (not yet fully
   implemented!). */

TSIL_COMPLEX pi1_squarkQCD (TSIL_REAL);
TSIL_COMPLEX Pi1_squarkQCD (TSIL_REAL);
TSIL_COMPLEX sumX (void);
TSIL_COMPLEX sumXP (void);
TSIL_COMPLEX SumX (void);
TSIL_COMPLEX SumXP (void);

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* This is used for u,d,s,c squarks... */
/* Lagrangian masses used everywhere. */

TSIL_COMPLEX pi1_squarkQCD (TSIL_REAL m2squark) 
{
  TSIL_COMPLEX result, BFF, Bff;

  result = A_S (m2squark, Q2);
  bFF (0.0L, m2_gluino, m2squark, Q2, &BFF, &Bff);
  result += 2.0L * BFF;
  result += B_SV (m2squark, 0, m2squark, Q2);
  result *= Cq * g3 * g3;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* This is used for u,d,s,c squarks... */

TSIL_COMPLEX Pi1_squarkQCD (TSIL_REAL m2squark) 
{
  TSIL_COMPLEX result, BFF, Bff;

  result = A_S (m2squark, Q2);
  bFF (0.0L, m2_gluino, m2squark, Q2, &BFF, &Bff);
  result += 2.0L * BFF;
  result += B_SV (m2squark, 0, m2squark, Q2);
  result *= Cq * g3 * g3;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Lagrangian masses used everywhere. */

TSIL_COMPLEX pi1_stop (int i, int j, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;
  /* TSIL_COMPLEX sumx3 = 0.0L; */
  /* TSIL_COMPLEX sumxp = 0.0L; */
  TSIL_COMPLEX result4;
  /* TSIL_COMPLEX result4X; */

  TSIL_REAL phi02tmp, phip0tmp;

  /* First do QCD part from hep-ph/0502168 eq. (4.7). */
  for (k=0; k<2; k++) {
    QCDresult += (Lstop[i] * Lstopc[k] - Rstop[i] * Rstopc[k]) *
                 (Lstop[k] * Lstopc[j] - Rstop[k] * Rstopc[j]) *  
                 A_S (m2_stop[k], Q2);
  }

  bFF (m2_top, m2_gluino, s, Q2, &BFF, &Bff);

  if (i == j) QCDresult += 2.0L * BFF; 

  QCDresult += -2.0L * (Lstop[i] * Rstopc[j] + Rstop[i] * Lstopc[j]) * 
               m_top * m_gluino * Bff;

  if (i == j) QCDresult += B_SV (m2_stop[i], 0.0L, s, Q2);

  /* QCDresult *= SUMO_oneloopfactor * Cq * g3 * g3; */
  QCDresult *= Cq * g3 * g3;

  /* Rest is from PBMZ... */
  /* photon loop */
  if (i == j) result += (4.0L/9.0L) * e2 * B_SV (m2_stop[i], 0, s, Q2);

  /* Z loops */
  for (k=0; k<2; k++) {
    result += gZstopstopc[i][k] * gZstopstopc[k][j] * 
              B_SV (m2_stop[k], M2_Z, s, Q2);
  }

  result += 0.5L * gZZstopstopc[i][j] * A_V (M2_Z, Q2);

  /* W loops */
  for (k=0; k<2; k++) {
    result += 0.5L * g2 * Lstop[i] * Lstopc[j] * Lsbot[k] * Lsbotc[k] *  
              B_SV (m2_sbot[k], M2_W, s, Q2);
  }

  result += 0.5L * g2 * Lstop[i] * Lstopc[j] * A_V (M2_W, Q2);
  
  /* top-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (m2_top, m2_Neut[k], s, Q2, &BFF, &Bff);
    
    result += (YtNstopc[k][j] * conYtNstopc[k][i] +
               YtcNstop[k][i] * conYtcNstop[k][j]) * BFF;

    result += (YtNstopc[k][j] * YtcNstop[k][i] +
               conYtNstopc[k][i] * conYtcNstop[k][j]) *
               m_top * m_Neut[k] * Bff;
  }

  /* bottom-chargino loop */
  for (k=0; k<2; k++) {
    bFF (m2_bot, m2_Char[k], s, Q2, &BFF, &Bff);
    
    result += (YbCstopc[k][j] * conYbCstopc[k][i] +
               YbcCstop[k][i] * conYbcCstop[k][j]) * BFF;

    result += (YbCstopc[k][j] * YbcCstop[k][i] +
               conYbCstopc[k][i] * conYbcCstop[k][j]) *
               m_bot * m_Char[k] * Bff;
  }

  /* Save Goldstone value and replace with 0 temporarily */
  phi02tmp = m2_phi0[2];
  m2_phi0[2] = 0.0L;

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
  for (k=0; k<2; k++) {
    result += lambda0stopstopc[n][i][k] * lambda0stopstopc[n][k][j] *
              B_SS (m2_phi0[n], m2_stop[k], s, Q2);
  }}

  for (n=0; n<4; n++) {
    result += 0.5L * lambda00stopstopc[n][n][i][j] * A_S (m2_phi0[n], Q2);
  }

  m2_phi0[2] = phi02tmp;

  /* Save Goldstone value and replace with 0 temporarily */
  phip0tmp = m2_phip[0];
  m2_phip[0] = 0.0L;

  /* charged Higgs loops */
  for (n=0; n<2; n++) {
  for (k=0; k<2; k++) {
      result += SUMO_CONJ(lambdapsbotstopc[n][k][i]) *
	lambdapsbotstopc[n][k][j] *
	B_SS (m2_phip[n], m2_sbot[k], s, Q2);
  }}

  for (n=0; n<2; n++) {
    result += lambdapmstopstopc[n][n][i][j] * A_S (m2_phip[n], Q2);
  }

  m2_phip[0] = phip0tmp;

  result4 = 0.0L;

  for (k=0; k<2; k++) {
    result4 += ytop * ytop * (
               3.0L * (Lstop[i] * Rstopc[j] * Lstopc[k] * Rstop[k] +
                       Rstop[i] * Lstopc[j] * Rstopc[k] * Lstop[k]) +
               Lstop[i] * Lstopc[j] * Rstop[k] * Rstopc[k] +
               Rstop[i] * Rstopc[j] * Lstop[k] * Lstopc[k]) * 
               A_S (m2_stop[k], Q2);

    result4 += (ytop * ytop * Rstop[i] * Rstopc[j] * Lsbot[k] * Lsbotc[k]
              + ybot * ybot * Lstop[i] * Lstopc[j] * Rsbot[k] * Rsbotc[k]) *
               A_S (m2_sbot[k], Q2);
  }

  result4 += g2 * 0.25L * Lstop[i] * Lstopc[j] * sumX ();

  result4 += gp2 * 
             (Lstop[i] * Lstopc[j]/6.0L - 2.0L/3.0L * Rstop[i] * Rstopc[j]) * 
             sumXP ();

  for (k=0; k<2; k++) {
    result4 += (g2/4.0L) * Lstop[i] * Lstopc[j] *
               Lstop[k] * Lstopc[k] * A_S (m2_stop[k], Q2);

    result4 += (g2/2.0L) * Lstop[i] * Lstopc[j] * 
               Lsbot[k] * Lsbotc[k] * A_S (m2_sbot[k], Q2);

    result4 += gp2 * 
             (Lstop[i] * Lstopc[k]/6.0L -2.0L/3.0L * Rstop[i] * Rstopc[k]) *
             (Lstop[k] * Lstopc[j]/6.0L -2.0L/3.0L * Rstop[k] * Rstopc[j]) *
             A_S (m2_stop[k], Q2);
  }

  result += result4;
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Lagrangian masses used everywhere. */

TSIL_COMPLEX pi1_sbot (int i, int j, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;
  TSIL_COMPLEX result4;

  TSIL_REAL phi02tmp, phip0tmp;

  /* First do QCD part from hep-ph/0502168 eq. (4.7). */
  for (k=0; k<2; k++) {
    QCDresult += (Lsbot[i] * Lsbotc[k] - Rsbot[i] * Rsbotc[k]) *
                 (Lsbot[k] * Lsbotc[j] - Rsbot[k] * Rsbotc[j]) *  
                 A_S (m2_sbot[k], Q2);
  }

  bFF (m2_bot, m2_gluino, s, Q2, &BFF, &Bff);

  if (i == j) QCDresult += 2.0L * BFF; 

  QCDresult += -2.0L * (Lsbot[i] * Rsbotc[j] + Rsbot[i] * Lsbotc[j]) * 
                m_bot * m_gluino * Bff;

  if (i == j) QCDresult += B_SV (m2_sbot[i], 0.0L, s, Q2);

  QCDresult *= Cq * g3 * g3;

  /* Rest is from PBMZ... */
  /* photon loop */
  if (i == j) result += (1.0L/9.0L) * e2 * B_SV (m2_sbot[i], 0.0L, s, Q2);

  /* Z loops */
  for (k=0; k<2; k++) {
    result += gZsbotsbotc[i][k] * gZsbotsbotc[k][j] * 
              B_SV (m2_sbot[k], M2_Z, s, Q2);
  }

  result += 0.5L * gZZsbotsbotc[i][j] * A_V (M2_Z, Q2);

  /* W loops */
  for (k=0; k<2; k++) {
    result += 0.5L * g2 * Lsbot[i] * Lsbotc[j] * Lstop[k] * Lstopc[k] *  
              B_SV (m2_stop[k], M2_W, s, Q2);
  }

  result += 0.5L * g2 * Lsbot[i] * Lsbotc[j] * A_V (M2_W, Q2);
  
  /* bot-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (m2_bot, m2_Neut[k], s, Q2, &BFF, &Bff);
    
    result += (YbNsbotc[k][j] * conYbNsbotc[k][i] +
               YbcNsbot[k][i] * conYbcNsbot[k][j]) * BFF;

    result += (YbNsbotc[k][j] * YbcNsbot[k][i] +
               conYbNsbotc[k][i] * conYbcNsbot[k][j]) *
               m_bot * m_Neut[k] * Bff;
  }

  /* top-chargino loop */
  for (k=0; k<2; k++) {
    bFF (m2_top, m2_Char[k], s, Q2, &BFF, &Bff);
    
    result += (YtCsbotc[k][j] * conYtCsbotc[k][i] +
               YtcCsbot[k][i] * conYtcCsbot[k][j]) * BFF;

    result += (YtCsbotc[k][j] * YtcCsbot[k][i] +
               conYtCsbotc[k][i] * conYtcCsbot[k][j]) *
               m_top * m_Char[k] * Bff;
  }

  /* Save Goldstone value and replace with 0 temporarily */
  phi02tmp = m2_phi0[2];
  m2_phi0[2] = 0.0L;

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
  for (k=0; k<2; k++) {
    result += lambda0sbotsbotc[n][i][k] * lambda0sbotsbotc[n][k][j] *
              B_SS (m2_phi0[n], m2_sbot[k], s, Q2);
  }}

  for (n=0; n<4; n++) {
    result += 0.5L * lambda00sbotsbotc[n][n][i][j] * A_S (m2_phi0[n], Q2);
  }

  m2_phi0[2] = phi02tmp;

  /* Save Goldstone value and replace with 0 temporarily */
  phip0tmp = m2_phip[0];
  m2_phip[0] = 0.0L;

  /* charged Higgs loops */
  for (n=0; n<2; n++) {
  for (k=0; k<2; k++) {
    result +=
      lambdapsbotstopc[n][i][k] * SUMO_CONJ(lambdapsbotstopc[n][j][k]) *
      B_SS (m2_phip[n], m2_stop[k], s, Q2);
  }}

  for (n=0; n<2; n++) {
    result += lambdapmsbotsbotc[n][n][i][j] * A_S (m2_phip[n], Q2);
  }

  m2_phip[0] = phip0tmp;

  result4 = 0.0L;

  for (k=0; k<2; k++) {
    result4 += ybot * ybot * (
               3.0L * (Lsbot[i] * Rsbotc[j] * Lsbotc[k] * Rsbot[k] +
                       Rsbot[i] * Lsbotc[j] * Rsbotc[k] * Lsbot[k]) +
               Lsbot[i] * Lsbotc[j] * Rsbot[k] * Rsbotc[k] +
               Rsbot[i] * Rsbotc[j] * Lsbot[k] * Lsbotc[k]) * 
               A_S (m2_sbot[k], Q2);

    /* This is the original: (looks plausible based on couplings in 0405022) */
    result4 += (ybot * ybot * Rsbot[i] * Rsbotc[j] * Lstop[k] * Lstopc[k]
               +ytop * ytop * Lsbot[i] * Lsbotc[j] * Rstop[k] * Rstopc[k]) *
               A_S (m2_stop[k], Q2);

    result4 += ybot * ytau * (Lsbot[i] * Rsbotc[j] * Lstauc[k] * Rstau[k] +
                              Rsbot[i] * Lsbotc[j] * Rstauc[k] * Lstau[k]) *
               A_S (m2_stau[k], Q2);
  }

  result4 += -g2 * 0.25L * Lsbot[i] * Lsbotc[j] * sumX ();

  result4 += gp2 * (Lsbot[i] * Lsbotc[j]/6.0L + Rsbot[i] * Rsbotc[j]/3.0) *
             sumXP ();

  for (k=0; k<2; k++) {
    result4 += (g2/4.0L) * Lsbot[i] * Lsbotc[j] *
               Lsbot[k] * Lsbotc[k] * A_S (m2_sbot[k], Q2);

    result4 += (g2/2.0L) * Lsbot[i] * Lsbotc[j] * 
               Lstop[k] * Lstopc[k] * A_S (m2_stop[k], Q2);

    result4 += gp2 *
             (Lsbot[i] * Lsbotc[k]/6.0L + Rsbot[i] * Rsbotc[k]/3.0L) *
             (Lsbot[k] * Lsbotc[j]/6.0L + Rsbot[k] * Rsbotc[j]/3.0L) *
             A_S (m2_sbot[k], Q2);
  }

  result += result4;
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Lagrangian masses everywhere. */

TSIL_COMPLEX pi1_sbot_orig (int i, int j, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;
  TSIL_COMPLEX result4;

  /* First do QCD part from hep-ph/0502168 eq. (4.7).   */
  for (k=0; k<2; k++) {
    QCDresult += (Lsbot[i] * Lsbotc[k] - Rsbot[i] * Rsbotc[k]) *
                 (Lsbot[k] * Lsbotc[j] - Rsbot[k] * Rsbotc[j]) *  
                 A_S (m2_sbot[k], Q2);
  }

  bFF (m2_bot, m2_gluino, s, Q2, &BFF, &Bff);

  if (i == j) QCDresult += 2.0L * BFF; 

  QCDresult += -2.0L * (Lsbot[i] * Rsbotc[j] + Rsbot[i] * Lsbotc[j]) * 
               m_bot * m_gluino * Bff;

  if (i == j) QCDresult += B_SV (m2_sbot[i], 0.0L, s, Q2);

  QCDresult *= Cq * g3 * g3;

  /* photon loop */
  if (i == j) result += (1.0L/9.0L) * e2 * B_SV (m2_sbot[i], 0.0L, s, Q2);

  /* Z loops */
  for (k=0; k<2; k++) {
    result += gZsbotsbotc[i][k] * gZsbotsbotc[k][j] *
              B_SV (m2_sbot[k], M2_Z, s, Q2);
  }

  result += 0.5L * gZZsbotsbotc[i][j] * A_V (M2_Z, Q2);

  /* W loops */
  for (k=0; k<2; k++) {
    result += 0.5L * g2 * Lsbot[i] * Lsbotc[j] * Lstop[k] * Lstopc[k] *
              B_SV (m2_stop[k], M2_W, s, Q2);
  }

  result += 0.5L * g2 * Lsbot[i] * Lsbotc[j] * A_V (M2_W, Q2);

  /* bottom-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (m2_bot, m2_Neut[k], s, Q2, &BFF, &Bff);

    result += (YbNsbotc[k][j] * conYbNsbotc[k][i] +
               YbcNsbot[k][i] * conYbcNsbot[k][j]) * BFF;

    result += (YbNsbotc[k][j] * YbcNsbot[k][i] +
               conYbNsbotc[k][i] * conYbcNsbot[k][j]) *
               m_bot * m_Neut[k] * Bff;
  }

  /* top-chargino loop */
  for (k=0; k<2; k++) {
    bFF (m2_top, m2_Char[k], s, Q2, &BFF, &Bff);

    result += (YtCsbotc[k][j] * conYtCsbotc[k][i] +
               YtcCsbot[k][i] * conYtcCsbot[k][j]) * BFF;

    result += (YtCsbotc[k][j] * YtcCsbot[k][i] +
               conYtCsbotc[k][i] * conYtcCsbot[k][j]) *
               m_top * m_Char[k] * Bff;
  }

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
  for (k=0; k<2; k++) {
    result += lambda0sbotsbotc[n][i][k] * lambda0sbotsbotc[n][k][j] *
              B_SS (m2_phi0[n], m2_sbot[k], s, Q2);
  }}

  for (n=0; n<4; n++) {
    result += 0.5L * lambda00sbotsbotc[n][n][i][j] * A_S (m2_phi0[n], Q2);
  }

  /* charged Higgs loops */
  for (n=0; n<2; n++) {
  for (k=0; k<2; k++) {
    result +=
      lambdapsbotstopc[n][i][k] * SUMO_CONJ(lambdapsbotstopc[n][j][k]) *
      B_SS (m2_phip[n], m2_stop[k], s, Q2)
      ;
  }}

  for (n=0; n<2; n++) {
    result += lambdapmsbotsbotc[n][n][i][j] * A_S (m2_phip[n], Q2);
  }
  
  result4 = 0.0L;

  for (k=0; k<2; k++) {
    result4 += ybot * ybot * (
               3.0L * (Lsbot[i] * Rsbotc[j] * Lsbotc[k] * Rsbot[k] +
                       Rsbot[i] * Lsbotc[j] * Rsbotc[k] * Lsbot[k]) +
               Lsbot[i] * Lsbotc[j] * Rsbot[k] * Rsbotc[k] +
               Rsbot[i] * Rsbotc[j] * Lsbot[k] * Lsbotc[k]) *
               A_S (m2_sbot[k], Q2);

    result4 += (ybot * ybot * Rsbot[i] * Rsbotc[j] * Lstop[k] * Lstopc[k]
               +ytop * ytop * Lsbot[i] * Lsbotc[j] * Rstop[k] * Rstopc[k]) *
               A_S (m2_stop[k], Q2);

    result4 += ybot * ytau * (Lsbot[i] * Rsbotc[j] * Lstauc[k] * Rstau[k] +
                              Rsbot[i] * Lsbotc[j] * Rstauc[k] * Lstau[k]) *
               A_S (m2_stau[k], Q2);
  }     

  result4 += -g2 * 0.25L * Lsbot[i] * Lsbotc[j] * sumX ();

  result4 += gp2 * (Lsbot[i] * Lsbotc[j]/6.0L + Rsbot[i] * Rsbotc[j]/3.0) *
             sumXP ();

  for (k=0; k<2; k++) {
    result4 += (g2/4.0L) * Lsbot[i] * Lsbotc[j] * Lsbot[k] * Lsbotc[k] * A_S (m2_sbot[k], Q2);

    result4 += (g2/2.0L) * Lsbot[i] * Lsbotc[j] * Lstop[k] * Lstopc[k] * A_S (m2_stop[k], Q2);

    result4 += gp2 *
             (Lsbot[i] * Lsbotc[k]/6.0L + Rsbot[i] * Rsbotc[k]/3.0L) *
             (Lsbot[k] * Lsbotc[j]/6.0L + Rsbot[k] * Rsbotc[j]/3.0L) *
             A_S (m2_sbot[k], Q2);
  }

  result += result4;
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX pi1_sdL (int i, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;

  QCDresult = Pi1_squarkQCD (m2_sdL[i]);

  /* photon loop */
  result += (1.0L/9.0L) * e2 * B_SV (m2_sdL[i], 0, s, Q2);

  /* Z loops */
  result += gZsdLsdLc * gZsdLsdLc * B_SV (m2_sdL[i], M2_Z, s, Q2);
  result += 0.5L * gZZsdLsdLc * A_V (M2_Z, Q2);

  /* W loops */
  result += 0.5L * g2 * B_SV (m2_suL[i], M2_W, s, Q2);
  result += 0.5L * g2 * A_V (M2_W, Q2);

  /* down-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (0.0L, m2_Neut[k], s, Q2, &BFF, &Bff);
    result += YdNsdLc[k] * conYdNsdLc[k] * BFF;
  }

  /* up-chargino loop */
  for (k=0; k<2; k++) {
    bFF (0.0L, m2_Char[k], s, Q2, &BFF, &Bff);
    result += YuCsdLc[k] * conYuCsdLc[k] * BFF;
  }

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
    result += lambda0sdLsdLc[n] * lambda0sdLsdLc[n] *
              B_SS (m2_phi0[n], m2_sdL[i], s, Q2);
    result += 0.5L * lambda00sdLsdLc[n][n] * A_S (m2_phi0[n], Q2);
  }

  /* charged Higgs loops */
  for (n=0; n<2; n++) {
    result += lambdapsdLsuLc[n] * SUMO_CONJ(lambdapsdLsuLc[n]) *
              B_SS (m2_phip[n], m2_suL[i], s, Q2);
    result += lambdapmsdLsdLc[n][n] * A_S (m2_phip[n], Q2);
  }
  
  /* 4-sfermion interactions */
  result += -g2 * 0.25L * sumX ();
  result += gp2 * (1.0L/6.0L) * sumXP ();
  result += (g2/4.0L) * A_S (m2_sdL[i], Q2);
  result += (g2/2.0L) * A_S (m2_suL[i], Q2);
  result += gp2 * (1.0L/6.0L) * (1.0L/6.0L) * A_S (m2_sdL[i], Q2);

  /* result *= SUMO_oneloopfactor; */
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX Pi1_sdL (int i, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;

  QCDresult = Pi1_squarkQCD (M2_sdL[i]);

  /* photon loop */
  result += (1.0L/9.0L) * e2 * B_SV (M2_sdL[i], 0, s, Q2);

  /* Z loops */
  result += gZsdLsdLc * gZsdLsdLc * B_SV (M2_sdL[i], M2_Z, s, Q2);
  result += 0.5L * gZZsdLsdLc * A_V (M2_Z, Q2);

  /* W loops */
  result += 0.5L * g2 * B_SV (M2_suL[i], M2_W, s, Q2);
  result += 0.5L * g2 * A_V (M2_W, Q2);

  /* down-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (0.0L, M2_Neut[k], s, Q2, &BFF, &Bff);
    result += YdNsdLc[k] * conYdNsdLc[k] * BFF;
  }

  /* up-chargino loop */
  for (k=0; k<2; k++) {
    bFF (0.0L, M2_Char[k], s, Q2, &BFF, &Bff);
    result += YuCsdLc[k] * conYuCsdLc[k] * BFF;
  }

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
    result += lambda0sdLsdLc[n] * lambda0sdLsdLc[n] *
              B_SS (M2_phi0[n], M2_sdL[i], s, Q2);
    result += 0.5L * lambda00sdLsdLc[n][n] * A_S (M2_phi0[n], Q2);
  }

  /* charged Higgs loops */
  for (n=0; n<2; n++) {
    result += lambdapsdLsuLc[n] * SUMO_CONJ(lambdapsdLsuLc[n]) *
              B_SS (M2_phip[n], M2_suL[i], s, Q2);
    result += lambdapmsdLsdLc[n][n] * A_S (M2_phip[n], Q2);
  }
  
  /* 4-sfermion interactions */
  result += -g2 * 0.25L * SumX ();
  result += gp2 * (1.0L/6.0L) * SumXP ();
  result += (g2/4.0L) * A_S (M2_sdL[i], Q2);
  result += (g2/2.0L) * A_S (M2_suL[i], Q2);
  result += gp2 * (1.0L/6.0L) * (1.0L/6.0L) * A_S (M2_sdL[i], Q2);

  /* result *= SUMO_oneloopfactor; */
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX pi1_sdR (int i, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;

  QCDresult = Pi1_squarkQCD (m2_sdR[i]);

  /* photon loop */
  result += (1.0L/9.0L) * e2 * B_SV (m2_sdR[i], 0, s, Q2);

  /* Z loops */
  result += gZsdRsdRc * gZsdRsdRc * B_SV (m2_sdR[i], M2_Z, s, Q2);
  result += 0.5L * gZZsdRsdRc * A_V (M2_Z, Q2);

  /* down-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (0.0L, m2_Neut[k], s, Q2, &BFF, &Bff);
    result += YdcNsdR[k] * conYdcNsdR[k] * BFF;
  }

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
    result += lambda0sdRsdRc[n] * lambda0sdRsdRc[n] *
              B_SS (m2_phi0[n], m2_sdR[i], s, Q2);
    result += 0.5L * lambda00sdRsdRc[n][n] * A_S (m2_phi0[n], Q2);
  }

  /* charged Higgs loop */
  for (n=0; n<2; n++) {
    result += lambdapmsdRsdRc[n][n] * A_S (m2_phip[n], Q2);
  }
  
  /* 4-sfermion interactions */
  result += gp2 * (1.0L/3.0L) * sumXP ();
  result += gp2 * (1.0L/3.0L) * (1.0L/3.0L) * A_S (m2_sdR[i], Q2);

  /* result *= SUMO_oneloopfactor; */
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX Pi1_sdR (int i, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;

  QCDresult = Pi1_squarkQCD (M2_sdR[i]);

  /* photon loop */
  result += (1.0L/9.0L) * e2 * B_SV (M2_sdR[i], 0, s, Q2);

  /* Z loops */
  result += gZsdRsdRc * gZsdRsdRc * B_SV (M2_sdR[i], M2_Z, s, Q2);
  result += 0.5L * gZZsdRsdRc * A_V (M2_Z, Q2);

  /* down-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (0.0L, M2_Neut[k], s, Q2, &BFF, &Bff);
    result += YdcNsdR[k] * conYdcNsdR[k] * BFF;
  }

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
    result += lambda0sdRsdRc[n] * lambda0sdRsdRc[n] *
              B_SS (M2_phi0[n], M2_sdR[i], s, Q2);
    result += 0.5L * lambda00sdRsdRc[n][n] * A_S (M2_phi0[n], Q2);
  }

  /* charged Higgs loop */
  for (n=0; n<2; n++) {
    result += lambdapmsdRsdRc[n][n] * A_S (M2_phip[n], Q2);
  }
  
  /* 4-sfermion interactions */
  result += gp2 * (1.0L/3.0L) * SumXP ();
  result += gp2 * (1.0L/3.0L) * (1.0L/3.0L) * A_S (M2_sdR[i], Q2);

  /* result *= SUMO_oneloopfactor; */
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX pi1_suL (int i, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;

  QCDresult = Pi1_squarkQCD (m2_suL[i]);

  /* photon loop */
  result += (4.0L/9.0L) * e2 * B_SV (m2_suL[i], 0, s, Q2);

  /* Z loops */
  result += gZsuLsuLc * gZsuLsuLc * B_SV (m2_suL[i], M2_Z, s, Q2);
  result += 0.5L * gZZsuLsuLc * A_V (M2_Z, Q2);

  /* W loops */
  result += 0.5L * g2 * B_SV (m2_sdL[i], M2_W, s, Q2);
  result += 0.5L * g2 * A_V (M2_W, Q2);
  
  /* up-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (0.0L, m2_Neut[k], s, Q2, &BFF, &Bff);
    result += YuNsuLc[k] * conYuNsuLc[k] * BFF;
  }

  /* down-chargino loop */
  for (k=0; k<2; k++) {
    bFF (0.0L, m2_Char[k], s, Q2, &BFF, &Bff);
    result += YdCsuLc[k] * conYdCsuLc[k] * BFF;
  }

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
    result += lambda0suLsuLc[n] * lambda0suLsuLc[n] *
              B_SS (m2_phi0[n], m2_suL[i], s, Q2);
    result += 0.5L * lambda00suLsuLc[n][n] * A_S (m2_phi0[n], Q2);
  }

  /* charged Higgs loops */
  for (n=0; n<2; n++) {
    result += SUMO_CONJ(lambdapsdLsuLc[n]) * lambdapsdLsuLc[n] * 
              B_SS (m2_phip[n], m2_sdL[i], s, Q2);
    result += lambdapmsuLsuLc[n][n] * A_S (m2_phip[n], Q2);
  }

  /* 4-sfermion interactions */
  result += g2 * 0.25L * sumX ();
  result += gp2 * (1.0/6.0L) * sumXP ();
  result += (g2/4.0L) * A_S (m2_suL[i], Q2);
  result += (g2/2.0L) * A_S (m2_sdL[i], Q2);
  result += gp2 * (1.0L/6.0L) * (1.0L/6.0L) * A_S (m2_suL[i], Q2);

  /* result *= SUMO_oneloopfactor; */
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX Pi1_suL (int i, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;

  QCDresult = Pi1_squarkQCD (M2_suL[i]);

  /* photon loop */
  result += (4.0L/9.0L) * e2 * B_SV (M2_suL[i], 0, s, Q2);

  /* Z loops */
  result += gZsuLsuLc * gZsuLsuLc * B_SV (M2_suL[i], M2_Z, s, Q2);
  result += 0.5L * gZZsuLsuLc * A_V (M2_Z, Q2);

  /* W loops */
  result += 0.5L * g2 * B_SV (M2_sdL[i], M2_W, s, Q2);
  result += 0.5L * g2 * A_V (M2_W, Q2);
  
  /* up-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (0.0L, M2_Neut[k], s, Q2, &BFF, &Bff);
    result += YuNsuLc[k] * conYuNsuLc[k] * BFF;
  }

  /* down-chargino loop */
  for (k=0; k<2; k++) {
    bFF (0.0L, M2_Char[k], s, Q2, &BFF, &Bff);
    result += YdCsuLc[k] * conYdCsuLc[k] * BFF;
  }

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
    result += lambda0suLsuLc[n] * lambda0suLsuLc[n] *
              B_SS (M2_phi0[n], M2_suL[i], s, Q2);
    result += 0.5L * lambda00suLsuLc[n][n] * A_S (M2_phi0[n], Q2);
  }

  /* charged Higgs loops */
  for (n=0; n<2; n++) {
    result += SUMO_CONJ(lambdapsdLsuLc[n]) * lambdapsdLsuLc[n] * 
              B_SS (M2_phip[n], M2_sdL[i], s, Q2);
    result += lambdapmsuLsuLc[n][n] * A_S (M2_phip[n], Q2);
  }

  /* 4-sfermion interactions */
  result += g2 * 0.25L * SumX ();
  result += gp2 * (1.0/6.0L) * SumXP ();
  result += (g2/4.0L) * A_S (M2_suL[i], Q2);
  result += (g2/2.0L) * A_S (M2_sdL[i], Q2);
  result += gp2 * (1.0L/6.0L) * (1.0L/6.0L) * A_S (M2_suL[i], Q2);

  /* result *= SUMO_oneloopfactor; */
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX pi1_suR (int i, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;

  QCDresult = Pi1_squarkQCD (m2_suR[i]);

  /* photon loop */
  result += (4.0L/9.0L) * e2 * B_SV (m2_suR[i], 0, s, Q2);

  /* Z loops */
  result += gZsuRsuRc * gZsuRsuRc * B_SV (m2_suR[i], M2_Z, s, Q2);
  result += 0.5L * gZZsuRsuRc * A_V (M2_Z, Q2);

  /* up-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (M2_top, m2_Neut[k], s, Q2, &BFF, &Bff);    
    result += YucNsuR[k] * conYucNsuR[k] * BFF;
  }

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
    result += lambda0sdRsdRc[n] * lambda0sdRsdRc[n] *
              B_SS (m2_phi0[n], m2_sdR[i], s, Q2);
    result += 0.5L * lambda00sdRsdRc[n][n] * A_S (m2_phi0[n], Q2);
  }

  /* charged Higgs loops */
  for (n=0; n<2; n++) {
    result += lambdapmsuRsuRc[n][n] * A_S (m2_phip[n], Q2);
  }

  /* 4-sfermion interactions */
  result += gp2 * (-2.0L/3.0L) * sumXP ();
  result += gp2 * (-2.0L/3.0L) * (-2.0L/3.0L) * A_S (m2_sdR[i], Q2);

  /* result *= SUMO_oneloopfactor; */
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX Pi1_suR (int i, TSIL_REAL s) 
{
  int k, n;
  TSIL_COMPLEX BFF, Bff;
  TSIL_COMPLEX QCDresult = 0.0L;
  TSIL_COMPLEX result = 0.0L;

  QCDresult = Pi1_squarkQCD (M2_suR[i]);

  /* photon loop */
  result += (4.0L/9.0L) * e2 * B_SV (M2_suR[i], 0, s, Q2);

  /* Z loops */
  result += gZsuRsuRc * gZsuRsuRc * B_SV (M2_suR[i], M2_Z, s, Q2);
  result += 0.5L * gZZsuRsuRc * A_V (M2_Z, Q2);

  /* up-neutralino loop */
  for (k=0; k<4; k++) {
    bFF (M2_top, M2_Neut[k], s, Q2, &BFF, &Bff);    
    result += YucNsuR[k] * conYucNsuR[k] * BFF;
  }

  /* neutral Higgs loops */
  for (n=0; n<4; n++) {
    result += lambda0sdRsdRc[n] * lambda0sdRsdRc[n] *
              B_SS (M2_phi0[n], M2_sdR[i], s, Q2);
    result += 0.5L * lambda00sdRsdRc[n][n] * A_S (M2_phi0[n], Q2);
  }

  /* charged Higgs loops */
  for (n=0; n<2; n++) {
    result += lambdapmsuRsuRc[n][n] * A_S (M2_phip[n], Q2);
  }

  /* 4-sfermion interactions */
  result += gp2 * (-2.0L/3.0L) * sumXP ();
  result += gp2 * (-2.0L/3.0L) * (-2.0L/3.0L) * A_S (M2_sdR[i], Q2);

  /* result *= SUMO_oneloopfactor; */
  result += QCDresult;

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX sumX (void)
{
  int k;
  TSIL_COMPLEX result = 0.0L;

  for (k=0; k<2; k++) {
    result +=  3.0L * Lstop[k] * Lstopc[k] * A_S (m2_stop[k], Q2); 
    result += -3.0L * Lsbot[k] * Lsbotc[k] * A_S (m2_sbot[k], Q2); 
    result += -Lstau[k] * Lstauc[k] * A_S (m2_stau[k], Q2); 
    result +=  3.0L * A_S (m2_suL[k], Q2); 
    result += -3.0L * A_S (m2_sdL[k], Q2); 
    result +=  A_S (m2_snu[k], Q2); 
    result += -A_S (m2_seL[k], Q2); 
  }

  result += A_S (m2_snu[2], Q2);

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX SumX (void)
{
  int k;
  TSIL_COMPLEX result = 0.0L;

  for (k=0; k<2; k++) {
    result +=  3.0L * Lstop[k] * Lstopc[k] * A_S (M2_stop[k], Q2); 
    result += -3.0L * Lsbot[k] * Lsbotc[k] * A_S (M2_sbot[k], Q2); 
    result += -Lstau[k] * Lstauc[k] * A_S (M2_stau[k], Q2); 
    result +=  3.0L * A_S (M2_suL[k], Q2); 
    result += -3.0L * A_S (M2_sdL[k], Q2); 
    result +=  A_S (M2_snu[k], Q2); 
    result += -A_S (M2_seL[k], Q2); 
  }

  result += A_S (M2_snu[2], Q2);

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX sumXP (void)
{
  int k;
  TSIL_COMPLEX result = 0.0L;

  result = 0.0L;

  for (k=0; k<2; k++) {
    result +=  3.0L * (Lstop[k] * Lstopc[k]/6.0L - 
               2.0L/3.0L * Rstop[k] * Rstopc[k]) * A_S (m2_stop[k], Q2); 
    result +=  3.0L * (Lsbot[k] * Lsbotc[k]/6.0L + 
               1.0L/3.0L * Rsbot[k] * Rsbotc[k]) * A_S (m2_sbot[k], Q2); 
    result +=  (-Lstau[k] * Lstauc[k]/2.0L + Rstau[k] * Rstauc[k]) * 
               A_S (m2_stau[k], Q2); 

    result +=  3.0L * (1.0L/6.0L) * A_S (m2_suL[k], Q2); 
    result +=  3.0L * (1.0L/6.0L) * A_S (m2_sdL[k], Q2); 
    result +=  3.0L * (-2.0L/3.0L) * A_S (m2_suR[k], Q2); 
    result +=  3.0L * (1.0L/3.0L) * A_S (m2_sdR[k], Q2); 
    result +=  (-1.0L/2.0L) * A_S (m2_snu[k], Q2); 
    result +=  (-1.0L/2.0L) * A_S (m2_seL[k], Q2); 
    result +=  A_S (m2_seR[k], Q2); 
  }

  result +=  -(1.0L/2.0L) * A_S (m2_snu[2], Q2); 

  return (result);
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX SumXP (void)
{
  int k;
  TSIL_COMPLEX result = 0.0L;

  result = 0.0L;

  for (k=0; k<2; k++) {
    result +=  3.0L * (Lstop[k] * Lstopc[k]/6.0L - 
               2.0L/3.0L * Rstop[k] * Rstopc[k]) * A_S (M2_stop[k], Q2); 
    result +=  3.0L * (Lsbot[k] * Lsbotc[k]/6.0L + 
               1.0L/3.0L * Rsbot[k] * Rsbotc[k]) * A_S (M2_sbot[k], Q2); 
    result +=  (-Lstau[k] * Lstauc[k]/2.0L + Rstau[k] * Rstauc[k]) * 
               A_S (M2_stau[k], Q2); 

    result +=  3.0L * (1.0L/6.0L) * A_S (M2_suL[k], Q2); 
    result +=  3.0L * (1.0L/6.0L) * A_S (M2_sdL[k], Q2); 
    result +=  3.0L * (-2.0L/3.0L) * A_S (M2_suR[k], Q2); 
    result +=  3.0L * (1.0L/3.0L) * A_S (M2_sdR[k], Q2); 
    result +=  (-1.0L/2.0L) * A_S (M2_snu[k], Q2); 
    result +=  (-1.0L/2.0L) * A_S (M2_seL[k], Q2); 
    result +=  A_S (M2_seR[k], Q2); 
  }

  result += -(1.0L/2.0L) * A_S (M2_snu[2], Q2); 

  return (result);
}
