/* Neutral Higgs self-energy functions. */

#include "supermodel.h"
#include "self_scalar.h"

#define VERBOSE 0

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

void SUMO_Backup (SUMO_MODEL *foo)
{
  int i;

  foo->gp = gp;
  foo->g = g;
  foo->g3 = g3;
  foo->ytop = ytop;
  foo->ybot = ybot;
  foo->ytau = ytau;
  foo->vu = vu;
  foo->vd = vd;
  foo->Q = Q;
  foo->m2Hu = m2Hu;
  foo->m2Hd = m2Hd;
  foo->mu = mu;
  foo->M1 = M1;
  foo->M2 = M2;
  foo->M3 = M3;
  foo->atop = atop;
  foo->abot = abot;
  foo->atau = atau;
  foo->b = b;
  for (i=0; i<3; i++) {
    foo->m2Q[i] = m2Q[i];
    foo->m2u[i] = m2u[i];
    foo->m2d[i] = m2d[i];
    foo->m2L[i] = m2L[i];
    foo->m2e[i] = m2e[i];
  }

 return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

void SUMO_Restore (SUMO_MODEL *foo)
{
  int i;
  
  gp = foo->gp; 
  g =  foo->g; 
  g3 = foo->g3; 
  ytop = foo->ytop; 
  ybot = foo->ybot; 
  ytau = foo->ytau; 
  vu = foo->vu; 
  vd = foo->vd; 
  Q =  foo->Q;
  m2Hu = foo->m2Hu; 
  m2Hd = foo->m2Hd; 
  mu = foo->mu; 
  M1 = foo->M1; 
  M2 = foo->M2;  
  M3 = foo->M3;  
  atop = foo->atop; 
  abot = foo->abot; 
  atau = foo->atau;
  b = foo->b;
  for (i=0; i<3; i++) {
    m2Q[i] = foo->m2Q[i];
    m2u[i] = foo->m2u[i];
    m2d[i] = foo->m2d[i];
    m2L[i] = foo->m2L[i];
    m2e[i] = foo->m2e[i];
  }

  /* Make sure everything is up to date... */
  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  return;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* One-loop self energy function for neutral Higgs particles
   (hep-ph/0405022, eq. 3.1), with tree masses as arguments of
   kinematic functions.
*/

SUMO_COMPLEX pi1_phi0 (int i, int j, SUMO_REAL s)
{
  int k, n;
  SUMO_COMPLEX BFF, Bff, temp, temp2;
  SUMO_COMPLEX result = 0.0L;
  char funcname[] = "pi1_phi0";

  for (k=0; k<4; k++)
    result += 0.5L * lambda0000[i][j][k][k] * A_S (m2_phi0[k], Q2);

  for (k=0; k<2; k++)
    result += lambda00pm[i][j][k][k] * A_S (m2_phip[k], Q2);

  for (k=0; k<2; k++) {
    result += 3.0L * lambda00stopstopc[i][j][k][k] * A_S (m2_stop[k], Q2);
    result += 3.0L * lambda00sbotsbotc[i][j][k][k] * A_S (m2_sbot[k], Q2);
    result += lambda00staustauc[i][j][k][k] * A_S (m2_stau[k], Q2);
    result += 3.0L * lambda00suLsuLc[i][j] * A_S (m2_suL[k], Q2);
    result += 3.0L * lambda00suRsuRc[i][j] * A_S (m2_suR[k], Q2);
    result += 3.0L * lambda00sdLsdLc[i][j] * A_S (m2_sdL[k], Q2);
    result += 3.0L * lambda00sdRsdRc[i][j] * A_S (m2_sdR[k], Q2);
    result += lambda00seLseLc[i][j] * A_S (m2_seL[k], Q2);
    result += lambda00seRseRc[i][j] * A_S (m2_seR[k], Q2);
    result += lambda00snusnuc[i][j] * A_S (m2_snu[k], Q2);
  }
  result += lambda00snusnuc[i][j] * A_S (m2_snu[2], Q2);

  for (k=0; k<2; k++) {
  for (n=0; n<2; n++) {
    temp = lambda0pm[i][k][n] * lambda0pm[j][n][k];
    if (TSIL_FABS(temp) > TSIL_TOL) {
      if ((m2_phip[k] < 0) && (0 == k)) m2_phip[k] = 1.0;
      if ((m2_phip[n] < 0) && (0 == n)) m2_phip[n] = 1.0;
      temp2 = B_SS (m2_phip[k], m2_phip[n], s, Q2);
      /* Because imaginary part coming from G+ is spurious: */
      if ((0 == k) || (0 == n)) temp2 = SUMO_CREAL(temp2);
      result += temp * temp2;
    }}}

  for (k=0; k<4; k++) {
  for (n=0; n<4; n++) {
    temp = lambda000[i][k][n] * lambda000[j][k][n];
    if (TSIL_FABS(temp) > TSIL_TOL) {
      if ((m2_phi0[k] < 0) && (2 == k)) m2_phi0[k] = 1.0;
      if ((m2_phi0[k] < 0) && (0 == k)) m2_phi0[k] = 10000.0;
      if ((m2_phi0[n] < 0) && (2 == n)) m2_phi0[n] = 1.0;
      if ((m2_phi0[n] < 0) && (0 == n)) m2_phi0[n] = 10000.0;
      temp2 = B_SS (m2_phi0[k], m2_phi0[n], s, Q2);
      /* Because imaginary part coming from G0 is spurious: */
      if((2 == k) || (2 == n)) temp2 = SUMO_CREAL(temp2);
      result += 0.5L * temp * temp2;
    }}}

  for (k=0; k<2; k++) {
  for (n=0; n<2; n++) {
    result += 3.0L * lambda0stopstopc[i][k][n] * lambda0stopstopc[j][n][k] *
              B_SS (m2_stop[k], m2_stop[n], s, Q2);
    result += 3.0L * lambda0sbotsbotc[i][k][n] * lambda0sbotsbotc[j][n][k] *
              B_SS (m2_sbot[k], m2_sbot[n], s, Q2);
    result += lambda0staustauc[i][k][n] * lambda0staustauc[j][n][k] *
              B_SS (m2_stau[k], m2_stau[n], s, Q2);
  }}

  for (k=0; k<2; k++) {
    result += 3.0L * lambda0suLsuLc[i] * lambda0suLsuLc[j] *
              B_SS (m2_suL[k], m2_suL[k], s, Q2);
    result += 3.0L * lambda0suRsuRc[i] * lambda0suRsuRc[j] *
              B_SS (m2_suR[k], m2_suR[k], s, Q2);
    result += 3.0L * lambda0sdLsdLc[i] * lambda0sdLsdLc[j] *
              B_SS (m2_sdL[k], m2_sdL[k], s, Q2);
    result += 3.0L * lambda0sdRsdRc[i] * lambda0sdRsdRc[j] *
              B_SS (m2_sdR[k], m2_sdR[k], s, Q2);
    result += lambda0seLseLc[i] * lambda0seLseLc[j] *
              B_SS (m2_seL[k], m2_seL[k], s, Q2);
    result += lambda0seRseRc[i] * lambda0seRseRc[j] *
              B_SS (m2_seR[k], m2_seR[k], s, Q2);
    result += lambda0snusnuc[i] * lambda0snusnuc[j] *
              B_SS (m2_snu[k], m2_snu[k], s, Q2);
  }

  result += lambda0snusnuc[i] * lambda0snusnuc[j] *
            B_SS (m2_snu[2], m2_snu[2], s, Q2);

  for (k=0; k<2; k++) {
  for (n=0; n<2; n++) {
    bFF (m2_Char[k], m2_Char[n], s, Q2, &BFF, &Bff);
    result += 2.0L * TSIL_CREAL(YCC0[k][n][i] * conYCC0[k][n][j]) * BFF;
    result += 2.0L * TSIL_CREAL(YCC0[k][n][i] * YCC0[n][k][j]) *
              m_Char[k] * m_Char[n] * Bff;
  }}

  for (k=0; k<4; k++) {
  for (n=0; n<4; n++) {
    bFF (m2_Neut[k], m2_Neut[n], s, Q2, &BFF, &Bff);
    result += TSIL_CREAL(YNN0[k][n][i] * conYNN0[k][n][j]) * BFF;
    result += TSIL_CREAL(YNN0[k][n][i] * YNN0[k][n][j]) *
              m_Neut[k] * m_Neut[n] * Bff;
  }}

  bFF (m2_top, m2_top, s, Q2, &BFF, &Bff);
  result += 6.0L * TSIL_CREAL (Yttcphi[i] * conYttcphi[j]) * BFF;
  result += 6.0L * TSIL_CREAL (Yttcphi[i] * Yttcphi[j]) * m2_top * Bff;

  bFF (m2_bot, m2_bot, s, Q2, &BFF, &Bff);
  result += 6.0L * TSIL_CREAL (Ybbcphi[i] * conYbbcphi[j]) * BFF;
  result += 6.0L * TSIL_CREAL (Ybbcphi[i] * Ybbcphi[j]) * m2_bot * Bff;

  bFF (m2_tau, m2_tau, s, Q2, &BFF, &Bff);
  result += 2.0L * TSIL_CREAL (Ytautaucphi[i] * conYtautaucphi[j]) * TSIL_CREAL(BFF);
  result += 2.0L * TSIL_CREAL (Ytautaucphi[i] * Ytautaucphi[j]) * m2_tau * TSIL_CREAL(Bff);

/* NOTE: TO REPRODUCE THE RESULTS OF 0405022 FIGURE 4, COMMENT OUT NEXT THREE LINES!!! */ 
  bFF (m2_tau, m2_tau, s, Q2, &BFF, &Bff);
  result += 2.0L * TSIL_CREAL (Ytautaucphi[i] * conYtautaucphi[j]) * BFF;
  result += 2.0L * TSIL_CREAL (Ytautaucphi[i] * Ytautaucphi[j]) * m2_tau * Bff;
/* NOTE: TO REPRODUCE THE RESULTS OF 0405022 FIGURE 4, COMMENT OUT PREVIOUS THREE LINES!!! */

  for (k=0; k<2; k++) {
    temp = B_SV (m2_phip[k], m2_W, s, Q2);
         /* Because imaginary part coming from Gp is spurious: */
    if (0 == k) temp = SUMO_CREAL(temp);
    result += 2.0L * TSIL_CREAL(gW0p[i][k] * TSIL_CONJ(gW0p[j][k])) * temp; 
  }

  for (k=0; k<4; k++) {
    temp = B_SV (m2_phi0[k], m2_Z, s, Q2);
    /* Because imaginary part coming from G0 is spurious: */
    if (2 == k) temp = SUMO_CREAL(temp);
    result += gZ00[i][k] * gZ00[j][k] * temp;
  }

  if (i==j) {
    result += 0.5L * gZZ00[i] * A_V (m2_Z, Q2);
    result += gWW00[i] * A_V (m2_W, Q2);
  }

  result += 0.5L * gZZ0[i] * gZZ0[j] * B_VV (m2_Z, m2_Z, s, Q2);

  result += gWW0[i] * gWW0[j] * B_VV (m2_W, m2_W, s, Q2);

  return result;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* One-loop self energy function for neutral Higgs particles
   (hep-ph/0405022, eq. 3.1), with POLE MASSES as arguments of
   kinematic functions. Use this if you are doing the Higgs pole mass
   calculation at one-loop order only.
*/

SUMO_COMPLEX Pi1_phi0 (int i, int j, SUMO_REAL s)
{
  SUMO_MODEL temp_model;
  SUMO_COMPLEX result;

  SUMO_Backup (&temp_model);
  SUMO_Tree_masses_from_pole ();
  result = pi1_phi0 (i, j, s);
  SUMO_Restore (&temp_model);

  return result;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* One-loop self energy function for neutral Higgs particles
   (hep-ph/0405022, eq. 3.1), with tree-level masses expanded around
   the corresponding pole masses, EXCEPT for the Goldstones.  Use this
   if you are doing the Higgs pole mass calculation at two-loop order.
*/

SUMO_COMPLEX pi1_phi0_X (int i, int j, SUMO_REAL s)
{
  int k, n;
  SUMO_COMPLEX BFF, Bff, BFpF, Bfpf, temp, temp2;
  char funcname[] = "pi1_phi0_X";
  SUMO_MODEL temp_model;
  SUMO_COMPLEX result = 0.0L;

  SUMO_Backup (&temp_model);
  M2_phi0[0] = m2_phi0[0];
  M2_phi0[2] = m2_phi0[2];
  M2_phip[0] = m2_phip[0];

  for (k=0; k<4; k++) {
    result += 0.5L * lambda0000[i][j][k][k] * (A_S (M2_phi0[k], Q2) +
              (m2_phi0[k] - M2_phi0[k]) * SUMO_Ap (M2_phi0[k], Q2));
  }

  for (k=0; k<2; k++)
    result += lambda00pm[i][j][k][k] * (A_S (M2_phip[k], Q2) +
              (m2_phip[k] - M2_phip[k]) * SUMO_Ap (M2_phip[k], Q2));

  for (k=0; k<2; k++) {
    result += 3.0L * lambda00stopstopc[i][j][k][k] * (A_S (M2_stop[k], Q2) +
              (m2_stop[k] - M2_stop[k]) * SUMO_Ap (M2_stop[k], Q2));
    result += 3.0L * lambda00sbotsbotc[i][j][k][k] * (A_S (M2_sbot[k], Q2) +
              (m2_sbot[k] - M2_sbot[k]) * SUMO_Ap (M2_sbot[k], Q2));
    result += lambda00staustauc[i][j][k][k] * (A_S (M2_stau[k], Q2) +
              (m2_stau[k] - M2_stau[k]) * SUMO_Ap (M2_stau[k], Q2));
    result += 3.0L * lambda00suLsuLc[i][j] * (A_S (M2_suL[k], Q2) +
              (m2_suL[k] - M2_suL[k]) * SUMO_Ap (M2_suL[k], Q2));
    result += 3.0L * lambda00suRsuRc[i][j] * (A_S (M2_suR[k], Q2) +
              (m2_suR[k] - M2_suR[k]) * SUMO_Ap (M2_suR[k], Q2));
    result += 3.0L * lambda00sdLsdLc[i][j] * (A_S (M2_sdL[k], Q2) +
              (m2_sdL[k] - M2_sdL[k]) * SUMO_Ap (M2_sdL[k], Q2));
    result += 3.0L * lambda00sdRsdRc[i][j] * (A_S (M2_sdR[k], Q2) +
              (m2_sdR[k] - M2_sdR[k]) * SUMO_Ap (M2_sdR[k], Q2));
    result += lambda00seLseLc[i][j] * (A_S (M2_seL[k], Q2) +
              (m2_seL[k] - M2_seL[k]) * SUMO_Ap (M2_seL[k], Q2));
    result += lambda00seRseRc[i][j] * (A_S (M2_seR[k], Q2) +
              (m2_seR[k] - M2_seR[k]) * SUMO_Ap (M2_seR[k], Q2));
    result += lambda00snusnuc[i][j] * (A_S (M2_snu[k], Q2) +
              (m2_snu[k] - M2_snu[k]) * SUMO_Ap (M2_snu[k], Q2));
  }

  result += lambda00snusnuc[i][j] * (A_S (M2_snu[2], Q2) +
            (m2_snu[2] - M2_snu[2]) * SUMO_Ap (M2_snu[2], Q2));

  for (k=0; k<2; k++) {
  for (n=0; n<2; n++) {
    temp = lambda0pm[i][k][n] * lambda0pm[j][n][k];
    if (TSIL_FABS(temp) > TSIL_TOL) {
      if ((M2_phip[k] < 0) && (0 == k)) M2_phip[k] = 1.0;
      if ((M2_phip[n] < 0) && (0 == n)) M2_phip[n] = 1.0;
      temp2 = B_SS (M2_phip[k], M2_phip[n], s, Q2) +
           (m2_phip[k] - M2_phip[k]) * B_SpS(M2_phip[k], M2_phip[n], s, Q2) +
           (m2_phip[n] - M2_phip[n]) * B_SpS(M2_phip[n], M2_phip[k], s, Q2);
      /* Because imaginary part coming from G+ is spurious: */
      if ((0 == k) || (0 == n)) temp2 = SUMO_CREAL(temp2);
      result += temp * temp2;
    }
  }}

  for (k=0; k<4; k++) {
  for (n=0; n<4; n++) {
    temp = lambda000[i][k][n] * lambda000[j][k][n];
    if (TSIL_FABS(temp) > TSIL_TOL) {
      if ((M2_phi0[k] < 0) && (2 == k)) M2_phi0[k] = 1.0;
      if ((M2_phi0[k] < 0) && (0 == k)) M2_phi0[k] = 10000.0;
      if ((M2_phi0[n] < 0) && (2 == n)) M2_phi0[n] = 1.0;
      if ((M2_phi0[n] < 0) && (0 == n)) M2_phi0[n] = 10000.0;
      temp2 = B_SS (M2_phi0[k], M2_phi0[n], s, Q2) +
        (m2_phi0[k] - M2_phi0[k]) * B_SpS(M2_phi0[k], M2_phi0[n], s, Q2) +
        (m2_phi0[n] - M2_phi0[n]) * B_SpS(M2_phi0[n], M2_phi0[k], s, Q2);
      /* Because imaginary part coming from G0 is spurious: */
      if ((2 == k) || (2 == n)) temp2 = SUMO_CREAL(temp2);
      result += 0.5L * temp * temp2;
    }
  }}

  for (k=0; k<2; k++) {
  for (n=0; n<2; n++) {
    result += 3.0L * lambda0stopstopc[i][k][n] * lambda0stopstopc[j][n][k] *
              (B_SS (M2_stop[k], M2_stop[n], s, Q2) +
               (m2_stop[k] - M2_stop[k]) * (B_SpS(M2_stop[k], M2_stop[n], s, Q2)) +
               (m2_stop[n] - M2_stop[n]) * (B_SpS(M2_stop[n], M2_stop[k], s, Q2)) );
    result += 3.0L * lambda0sbotsbotc[i][k][n] * lambda0sbotsbotc[j][n][k] *
              (B_SS (M2_sbot[k], M2_sbot[n], s, Q2) +
               (m2_sbot[k] - M2_sbot[k]) * (B_SpS(M2_sbot[k], M2_sbot[n], s, Q2)) +
               (m2_sbot[n] - M2_sbot[n]) * (B_SpS(M2_sbot[n], M2_sbot[k], s, Q2)) );
    result += lambda0staustauc[i][k][n] * lambda0staustauc[j][n][k] *
              (B_SS (M2_stau[k], M2_stau[n], s, Q2) +
               (m2_stau[k] - M2_stau[k]) * (B_SpS(M2_stau[k], M2_stau[n], s, Q2)) +
               (m2_stau[n] - M2_stau[n]) * (B_SpS(M2_stau[n], M2_stau[k], s, Q2)) );
  }}

  for (k=0; k<2; k++) {
    result += 3.0L * lambda0suLsuLc[i] * lambda0suLsuLc[j] *
              (B_SS (M2_suL[k], M2_suL[k], s, Q2) +
              2.0 * (m2_suL[k] - M2_suL[k]) * (B_SpS(M2_suL[k], M2_suL[k], s, Q2)) );
    result += 3.0L * lambda0suRsuRc[i] * lambda0suRsuRc[j] *
              (B_SS (M2_suR[k], M2_suR[k], s, Q2) +
              2.0 * (m2_suR[k] - M2_suR[k]) * (B_SpS(M2_suR[k], M2_suR[k], s, Q2)) );
    result += 3.0L * lambda0sdLsdLc[i] * lambda0sdLsdLc[j] *
              (B_SS (M2_sdL[k], M2_sdL[k], s, Q2) +
              2.0 * (m2_sdL[k] - M2_sdL[k]) * (B_SpS(M2_sdL[k], M2_sdL[k], s, Q2)) );
    result += 3.0L * lambda0sdRsdRc[i] * lambda0sdRsdRc[j] *
              (B_SS (M2_sdR[k], M2_sdR[k], s, Q2) +
              2.0 * (m2_sdR[k] - M2_sdR[k]) * (B_SpS(M2_sdR[k], M2_sdR[k], s, Q2)) );
    result += lambda0seLseLc[i] * lambda0seLseLc[j] *
              (B_SS (M2_seL[k], M2_seL[k], s, Q2) +
              2.0 * (m2_seL[k] - M2_seL[k]) * (B_SpS(M2_seL[k], M2_seL[k], s, Q2)) );
    result += lambda0seRseRc[i] * lambda0seRseRc[j] *
              (B_SS (M2_seR[k], M2_seR[k], s, Q2) +
              2.0 * (m2_seR[k] - M2_seR[k]) * (B_SpS(M2_seR[k], M2_seR[k], s, Q2)) );
    result += lambda0snusnuc[i] * lambda0snusnuc[j] *
              (B_SS (M2_snu[k], M2_snu[k], s, Q2) +
              2.0 * (m2_snu[k] - M2_snu[k]) * (B_SpS(M2_snu[k], M2_snu[k], s, Q2)) );
  }

  result += lambda0snusnuc[i] * lambda0snusnuc[j] *
            (B_SS (M2_snu[2], M2_snu[2], s, Q2) +
            2.0 * (m2_snu[2] - M2_snu[2]) * (B_SpS(M2_snu[2], M2_snu[2], s, Q2)) );

  for (k=0; k<2; k++) {
  for (n=0; n<2; n++) {
    bFF (M2_Char[k], M2_Char[n], s, Q2, &BFF, &Bff);
    bFpF (M2_Char[k], M2_Char[n], s, Q2, &BFpF, &Bfpf);
    result += 2.0L * TSIL_CREAL(YCC0[k][n][i] * conYCC0[k][n][j]) * (BFF +
              (m2_Char[k] - M2_Char[k]) * BFpF);
    result += 2.0L * TSIL_CREAL(YCC0[n][k][i] * conYCC0[n][k][j]) *
              (m2_Char[k] - M2_Char[k]) * BFpF;
    result += 2.0L * TSIL_CREAL(YCC0[k][n][i] * YCC0[n][k][j]) *
              M_Char[k] * M_Char[n] * (Bff +
              (m2_Char[k] - M2_Char[k]) * (Bfpf + Bff/(2.0 * M2_Char[k])));
    result += 2.0L * TSIL_CREAL(YCC0[n][k][i] * YCC0[k][n][j]) *
              M_Char[k] * M_Char[n] *
              (m2_Char[k] - M2_Char[k]) * (Bfpf + Bff/(2.0 * M2_Char[k]));
  }}

  for (k=0; k<4; k++) {
  for (n=0; n<4; n++) {
    bFF (M2_Neut[k], M2_Neut[n], s, Q2, &BFF, &Bff);
    bFpF (M2_Neut[k], M2_Neut[n], s, Q2, &BFpF, &Bfpf);
    result += TSIL_CREAL(YNN0[k][n][i] * conYNN0[k][n][j]) * (BFF +
              (m2_Neut[k] - M2_Neut[k]) * BFpF);
    result += TSIL_CREAL(YNN0[n][k][i] * conYNN0[n][k][j]) *
              (m2_Neut[k] - M2_Neut[k]) * BFpF;
    result += TSIL_CREAL(YNN0[k][n][i] * YNN0[k][n][j]) *
              M_Neut[k] * M_Neut[n] * (Bff +
              (m2_Neut[k] - M2_Neut[k]) * (Bfpf + Bff/(2.0 * M2_Neut[k])));
    result += TSIL_CREAL(YNN0[n][k][i] * YNN0[n][k][j]) *
              M_Neut[k] * M_Neut[n] *
              (m2_Neut[k] - M2_Neut[k]) * (Bfpf + Bff/(2.0 * M2_Neut[k]));
  }}

  bFF (M2_top_phys, M2_top_phys, s, Q2, &BFF, &Bff);
  bFpF (M2_top_phys, M2_top_phys, s, Q2, &BFpF, &Bfpf);
  result += 6.0L * TSIL_CREAL (Yttcphi[i] * conYttcphi[j]) *
            (BFF + 2.0 * (m2_top - M2_top_phys) * BFpF);
  result += 6.0L * TSIL_CREAL (Yttcphi[i] * Yttcphi[j]) * (
            M2_top_phys * Bff + (m2_top - M2_top_phys) * (2.0 * M2_top_phys * Bfpf + Bff));

  bFF (M2_bot_phys, M2_bot_phys, s, Q2, &BFF, &Bff);
  bFpF (M2_bot_phys, M2_bot_phys, s, Q2, &BFpF, &Bfpf);
  result += 6.0L * TSIL_CREAL (Ybbcphi[i] * conYbbcphi[j]) *
            (BFF + 2.0 * (m2_bot - M2_bot_phys) * BFpF);
  result += 6.0L * TSIL_CREAL (Ybbcphi[i] * Ybbcphi[j]) * (
            M2_bot_phys * Bff + (m2_bot - M2_bot_phys) * (2.0 * M2_bot_phys * Bfpf + Bff));

  bFF (M2_tau_phys, M2_tau_phys, s, Q2, &BFF, &Bff);
  bFpF (M2_tau_phys, M2_tau_phys, s, Q2, &BFpF, &Bfpf);
  result += 2.0L * TSIL_CREAL (Ytautaucphi[i] * conYtautaucphi[j]) *
            (BFF + 2.0 * (m2_tau - M2_tau_phys) * BFpF);
  result += 2.0L * TSIL_CREAL (Ytautaucphi[i] * Ytautaucphi[j]) * (
            M2_tau_phys * Bff + (m2_tau - M2_tau_phys) * (2.0 * M2_tau_phys * Bfpf + Bff));

  for (k=0; k<2; k++) {
    temp = B_SV(M2_phip[k], M2_W_phys, s, Q2) +
              (m2_phip[k] - M2_phip[k]) * B_SpV(M2_phip[k], M2_W_phys, s, Q2) +
              (m2_W - M2_W_phys) * B_SVp(M2_phip[k], M2_W_phys, s, Q2);
    /* Because imaginary part coming from Gp is spurious: */
    if (0 == k) temp = SUMO_CREAL(temp);
    result += 2.0L * TSIL_CREAL(gW0p[i][k] * TSIL_CONJ(gW0p[j][k])) * temp;
  }

  for (k=0; k<4; k++) {
    temp = B_SV(M2_phi0[k], M2_Z_phys, s, Q2) +
              (m2_phi0[k] - M2_phi0[k]) * B_SpV(M2_phi0[k], M2_Z_phys, s, Q2) +
              (m2_Z - M2_Z_phys) * B_SVp(M2_phi0[k], M2_Z_phys, s, Q2);
    /* Because imaginary part coming from G0 is spurious: */
    if (2 == k) temp = SUMO_CREAL(temp);
    result += gZ00[i][k] * gZ00[j][k] * temp;
  }

  if (i==j) {
    result += 0.5L * gZZ00[i] * (A_V (M2_Z_phys, Q2) +
              (m2_Z - M2_Z_phys) * A_Vp(M2_Z_phys,Q2));
    result += gWW00[i] * (A_V (M2_W_phys, Q2) +
              (m2_W - M2_W_phys) * A_Vp(M2_W_phys,Q2));
  }

  result += 0.5L * gZZ0[i] * gZZ0[j] * (B_VV(M2_Z_phys, M2_Z_phys, s, Q2) +
            2.0L * (m2_Z - M2_Z_phys) * B_VpV(M2_Z_phys, M2_Z_phys, s, Q2) );

  result += gWW0[i] * gWW0[j] * (B_VV (M2_W_phys, M2_W_phys, s, Q2) +
            2.0L * (m2_W - M2_W_phys) * B_VpV(M2_W_phys, M2_W_phys, s, Q2) );

  SUMO_Restore (&temp_model);

  return result;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Two-loop self energy function for neutral Higgs bosons evaluated at
   the specified squared momentum; QCD contribution (hep-ph/0405022,
   eqs. 4.1-4.5).
*/

SUMO_COMPLEX pi2QCD_phi0 (int i, int j, SUMO_REAL s)
{
  char funcname[] = "Pi2QCD_phi0";
  TSIL_DATA bar;
  TSIL_RESULT gaak;
  TSIL_RESULT stopData[2][2], sbotData[2][2];
  TSIL_COMPLEX vFFFFS, vFffFS, vfFfFS, vfFFfS, vFFffS, vffffS;
  TSIL_COMPLEX gSS, gFF, gff;
  TSIL_COMPLEX mSfSff, mSFSfF, mSFSFf;
  TSIL_COMPLEX vSSSFF, vSSSff;
  int k, l, m, n, p;
  SUMO_COMPLEX topgluonresult = 0; 
  SUMO_COMPLEX botgluonresult = 0; 
  SUMO_COMPLEX stopgluonresult = 0; 
  SUMO_COMPLEX sbotgluonresult = 0; 
  SUMO_COMPLEX squarkgluonresult = 0; 
  SUMO_COMPLEX topstopgluinoresult = 0; 
  SUMO_COMPLEX botsbotgluinoresult = 0; 
  SUMO_COMPLEX quarksquarkgluinoresult = 0; 
  SUMO_COMPLEX stopresult = 0;
  SUMO_COMPLEX sbotresult = 0;
  SUMO_COMPLEX squarkresult = 0;
  SUMO_COMPLEX totalresult;
  SUMO_COMPLEX prefactor = 4.0L * g3 * g3;

#ifdef SUMO_HIGGS_TIMING
  clock_t t0, t1;
  t0 = clock ();
#endif

  /* ===== hep-ph/0405022 EQ. 4.1 ===== */

  TSIL_SetParameters (&bar, m2_top, m2_top, m2_top, m2_top, 0.0L, Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  G_FF (&gaak, &gFF, &gff);

  topgluonresult = prefactor * 2.0L * (SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) * gFF
                   + SUMO_CREAL(Yttcphi[i] * Yttcphi[j]) * m2_top * gff);

  TSIL_SetParameters (&bar, m2_bot, m2_bot, m2_bot, m2_bot, 0.0L, Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  G_FF (&gaak, &gFF, &gff);

  botgluonresult = prefactor * 2.0L * (SUMO_CREAL(Ybbcphi[i] * conYbbcphi[j]) * gFF
                   + SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j]) * m2_bot * gff);

  for (k=0; k<2; k++) {
  for (l=k; l<2; l++) {
    TSIL_SetParameters (&bar,
        m2_stop[k], m2_stop[k], m2_stop[l], m2_stop[l], 0.0L, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &stopData[k][l]);
    if (l != k)
      TSIL_PermuteResult (&stopData[k][l], XZandYU, &stopData[l][k]);
  }}

  for (k=0; k<2; k++) {
  for (l=0; l<2; l++) {
    G_SS (&stopData[k][l], &gSS);
    stopgluonresult += prefactor * lambda0stopstopc[i][k][l] * lambda0stopstopc[j][l][k] * gSS;
  }}

  for (k=0; k<2; k++) {
  for (l=k; l<2; l++) {
    TSIL_SetParameters (&bar,
        m2_sbot[k], m2_sbot[k], m2_sbot[l], m2_sbot[l], 0.0L, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &sbotData[k][l]);
    if (l != k)
      TSIL_PermuteResult (&sbotData[k][l], XZandYU, &sbotData[l][k]);
  }}

  for (k=0; k<2; k++) {
  for (l=0; l<2; l++) {
    G_SS (&sbotData[k][l], &gSS);
    sbotgluonresult += prefactor * lambda0sbotsbotc[i][k][l] * lambda0sbotsbotc[j][l][k] * gSS;
  }}

  for (k=0; k<2; k++) {
    TSIL_SetParameters (&bar,
      m2_suL[k], m2_suL[k], m2_suL[k], m2_suL[k], 0.0L, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);
    G_SS (&gaak, &gSS);
    squarkgluonresult += prefactor * lambda0suLsuLc[i] * lambda0suLsuLc[j] * gSS;

    TSIL_SetParameters (&bar,
      m2_suR[k], m2_suR[k], m2_suR[k], m2_suR[k], 0.0L, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);
    G_SS (&gaak, &gSS);
    squarkgluonresult += prefactor * lambda0suRsuRc[i] * lambda0suRsuRc[j] * gSS;

    TSIL_SetParameters (&bar,
      m2_sdL[k], m2_sdL[k], m2_sdL[k], m2_sdL[k], 0.0L, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);
    G_SS (&gaak, &gSS);
    squarkgluonresult += prefactor * lambda0sdLsdLc[i] * lambda0sdLsdLc[j] * gSS;

    TSIL_SetParameters (&bar,
      m2_sdR[k], m2_sdR[k], m2_sdR[k], m2_sdR[k], 0.0L, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);
    G_SS (&gaak, &gSS);
    squarkgluonresult += prefactor * lambda0sdRsdRc[i] * lambda0sdRsdRc[j] * gSS;
  }

  for (k=0; k<2; k++) {
    stopgluonresult += prefactor * lambda00stopstopc[i][j][k][k] * W_SSSV (m2_stop[k], Q2);
    sbotgluonresult += prefactor * lambda00sbotsbotc[i][j][k][k] * W_SSSV (m2_sbot[k], Q2);
  }

  for (k=0; k<2; k++) {
    squarkgluonresult += prefactor * (lambda00suLsuLc[i][j] * W_SSSV (m2_suL[k], Q2)
                                    + lambda00suRsuRc[i][j] * W_SSSV (m2_suR[k], Q2)
                                    + lambda00sdLsdLc[i][j] * W_SSSV (m2_sdL[k], Q2)
                                    + lambda00sdRsdRc[i][j] * W_SSSV (m2_sdR[k], Q2));
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%4.3lf sec\n", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
#endif

  /* ===== EQ 4.2 ===== */

  /* Perform needed TSIL evaluations: */
  for (k=0; k<2; k++) {
  for (l=k; l<2; l++) {
      TSIL_SetParameters (&bar,
        m2_stop[k], m2_top, m2_stop[l], m2_top, m2_gluino, Q2);
      TSIL_Evaluate (&bar, s);
      TSIL_CopyResult (&bar, &stopData[k][l]);
      if (l != k)
	TSIL_PermuteResult (&stopData[k][l], XZandYU, &stopData[l][k]);
  }}

  for (k=0; k<2; k++) {
  for (l=0; l<2; l++) {

    M_SFSFF (&stopData[k][l], &mSfSff, &mSFSfF, &mSFSFf);

    topstopgluinoresult += prefactor * 4.0 * (SUMO_CREAL (
               (Yttcphi[i] * Lstop[k] * Lstopc[l] + conYttcphi[i] * Rstop[k] * Rstopc[l]) * 
	       lambda0stopstopc[j][l][k]
             + (Yttcphi[j] * Lstop[k] * Lstopc[l] + conYttcphi[j] * Rstop[k] * Rstopc[l]) * 
	       lambda0stopstopc[i][l][k]) * m_top * mSFSfF
          - SUMO_CREAL (Yttcphi[i] * Lstop[l] * Rstopc[k] * lambda0stopstopc[j][k][l] +
               Yttcphi[j] * Lstop[l] * Rstopc[k] * lambda0stopstopc[i][k][l]) * m_gluino * mSFSFf
          - SUMO_CREAL (conYttcphi[i] * Lstop[l] * Rstopc[k] * lambda0stopstopc[j][k][l]
             + conYttcphi[j] * Lstop[l] * Rstopc[k] * lambda0stopstopc[i][k][l]) * 
					      m2_top * m_gluino * mSfSff);
  }}

  for (k=0; k<2; k++) {
  for (l=k; l<2; l++) {
      TSIL_SetParameters (&bar,
        m2_sbot[k], m2_bot, m2_sbot[l], m2_bot, m2_gluino, Q2);
      TSIL_Evaluate (&bar, s);
      TSIL_CopyResult (&bar, &sbotData[k][l]);
      if (l != k)
	TSIL_PermuteResult (&sbotData[k][l], XZandYU, &sbotData[l][k]);
  }}

  for (k=0; k<2; k++) {
  for (l=0; l<2; l++) {

    M_SFSFF (&sbotData[k][l], &mSfSff, &mSFSfF, &mSFSFf);

    botsbotgluinoresult += prefactor * 4.0 * (SUMO_CREAL (
               (Ybbcphi[i] * Lsbot[k] * Lsbotc[l] + conYbbcphi[i] * Rsbot[k] * Rsbotc[l]) *
	       lambda0sbotsbotc[j][l][k]
             + (Ybbcphi[j] * Lsbot[k] * Lsbotc[l] + conYbbcphi[j] * Rsbot[k] * Rsbotc[l]) *
	       lambda0sbotsbotc[i][l][k]) * m_bot * mSFSfF
	  - SUMO_CREAL (Ybbcphi[i] * Lsbot[l] * Rsbotc[k] * lambda0sbotsbotc[j][k][l] +
               Ybbcphi[j]*Lsbot[l]*Rsbotc[k] * lambda0sbotsbotc[i][k][l]) * m_gluino * mSFSFf
	  - SUMO_CREAL (conYbbcphi[i] * Lsbot[l] * Rsbotc[k] * lambda0sbotsbotc[j][k][l]
             + conYbbcphi[j] * Lsbot[l] * Rsbotc[k] * lambda0sbotsbotc[i][k][l])
					      * m2_bot * m_gluino * mSfSff);
  }}

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%4.3lf sec\n", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
#endif

  /* ===== EQ. 4.3 ===== */

  /* (s)top: */
  for (k=0; k<2; k++) {

    /* Reorganize the earlier TSIL results from 4.2 for use in
       V_FFFFS: */
    gaak.x = stopData[k][k].y;
    gaak.y = stopData[k][k].u;
    gaak.z = stopData[k][k].u;
    gaak.u = stopData[k][k].v;
    gaak.v = stopData[k][k].z;
    gaak.s = stopData[k][k].s;
    gaak.qq = stopData[k][k].qq;
    gaak.B[xz] = stopData[k][k].B[yu];
    gaak.U[xzuv] = stopData[k][k].U[yuzv];
    gaak.V[xzuv] = stopData[k][k].V[yuzv];
    gaak.T[xuv] = stopData[k][k].T[uxv];
    gaak.T[uxv] = stopData[k][k].T[vxu];
    gaak.T[vxu] = stopData[k][k].T[xuv];
    gaak.S[uxv] = stopData[k][k].S[uxv];

    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    topstopgluinoresult += prefactor * 4.0 * (SUMO_CREAL (Yttcphi[i]*conYttcphi[j])
            * (vFFFFS + m2_top * vFffFS)
      + 2.0L*SUMO_CREAL (Yttcphi[i]*Yttcphi[j]) * m2_top * vfFfFS
      - 4.0L*SUMO_CREAL (Yttcphi[i]*conYttcphi[j])
            *SUMO_CREAL (Lstop[k]*Rstopc[k])
            * m_top * m_gluino * vFFffS
      - 2.0L*SUMO_CREAL (Yttcphi[i]*Yttcphi[j]*Lstop[k]*Rstopc[k])
            * m_top * m_gluino * vfFFfS
      - 2.0L*SUMO_CREAL (Yttcphi[i]*Yttcphi[j]*Lstopc[k]*Rstop[k])
            * m_top * m2_top * m_gluino * vffffS );
  }

  /* (s)bottom: */
  for (k=0; k<2; k++) {

    /* Reorganize the TSIL results for use in V_FFFFS: */
    gaak.x = sbotData[k][k].y;
    gaak.y = sbotData[k][k].u;
    gaak.z = sbotData[k][k].u;
    gaak.u = sbotData[k][k].v;
    gaak.v = sbotData[k][k].z;
    gaak.s = sbotData[k][k].s;
    gaak.qq = sbotData[k][k].qq;
    gaak.B[xz] = sbotData[k][k].B[yu];
    gaak.U[xzuv] = sbotData[k][k].U[yuzv];
    gaak.V[xzuv] = sbotData[k][k].V[yuzv];
    gaak.T[xuv] = sbotData[k][k].T[uxv];
    gaak.T[uxv] = sbotData[k][k].T[vxu];
    gaak.T[vxu] = sbotData[k][k].T[xuv];
    gaak.S[uxv] = sbotData[k][k].S[uxv];

    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    botsbotgluinoresult += prefactor * 4.0 * (SUMO_CREAL (Ybbcphi[i]*conYbbcphi[j])
            *(vFFFFS + m2_bot * vFffFS)
      + 2.0L*SUMO_CREAL (Ybbcphi[i]*Ybbcphi[j]) * m2_bot * vfFfFS
      - 4.0L*SUMO_CREAL (Ybbcphi[i]*conYbbcphi[j])
            *SUMO_CREAL (Lsbot[k]*Rsbotc[k])
            * m_bot * m_gluino * vFFffS
      - 2.0L*SUMO_CREAL (Ybbcphi[i]*Ybbcphi[j]*Lsbot[k]*Rsbotc[k])
            * m_bot * m_gluino * vfFFfS
      - 2.0L*SUMO_CREAL (Ybbcphi[i]*Ybbcphi[j]*Lsbotc[k]*Rsbot[k])
            * m_bot * m2_bot * m_gluino * vffffS );
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%4.3lf sec\n", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
#endif

  /* ===== EQ. 4.4 ===== */

  for (k=0; k<2; k++) {
    topstopgluinoresult += prefactor * 2.0 * lambda00stopstopc[i][j][k][k] *
                           W_SSFF (m2_stop[k], m2_stop[k], m2_top, m2_gluino, Q2);
  }

  for (k=0; k<2; k++) {
  for (m=0; m<2; m++) {
    topstopgluinoresult += prefactor * (-4.0L) * SUMO_CREAL(lambda00stopstopc[i][j][k][m]*
                           Lstopc[k]*Rstop[m]) * m_top * m_gluino *
                           W_SSff (m2_stop[k], m2_stop[m], m2_top, m2_gluino, Q2);
  }}

  for (k=0; k<2; k++) {
  for (m=0; m<2; m++) {
  for (n=0; n<2; n++) {
    V_SSSFF(&stopData[k][m], &stopData[k][n], &vSSSFF, &vSSSff);

    if (n==m) topstopgluinoresult += prefactor * 4.0L * SUMO_CREAL(lambda0stopstopc[i][k][m] *
                                     lambda0stopstopc[j][m][k]) * vSSSFF;

    topstopgluinoresult += prefactor * (-4.0L) * SUMO_CREAL(
                           lambda0stopstopc[i][k][m] * lambda0stopstopc[j][n][k] *
                           (Lstop[m] * Rstopc[n] + Rstop[m] * Lstopc[n])) *
                           m_top * m_gluino * vSSSff;
  }}}

  for (k=0; k<2; k++) {
    botsbotgluinoresult += prefactor * 2.0L * lambda00sbotsbotc[i][j][k][k] *
                           W_SSFF (m2_sbot[k], m2_sbot[k], m2_bot, m2_gluino, Q2);
  }

  for (k=0; k<2; k++) {
  for (m=0; m<2; m++) {
    botsbotgluinoresult += prefactor * (-4.0L) * SUMO_CREAL(lambda00sbotsbotc[i][j][k][m]*
                           Lsbotc[k]*Rsbot[m]) * m_bot * m_gluino *
	                   W_SSff (m2_sbot[k], m2_sbot[m], m2_bot, m2_gluino, Q2);
  }}


  for (k=0; k<2; k++) {
  for (m=0; m<2; m++) {
  for (n=0; n<2; n++) {
    V_SSSFF(&sbotData[k][m], &sbotData[k][n], &vSSSFF, &vSSSff);

    if (n==m) botsbotgluinoresult += prefactor * 4.0L * SUMO_CREAL(lambda0sbotsbotc[i][k][m] *
                                     lambda0sbotsbotc[j][m][k]) * vSSSFF;

    botsbotgluinoresult += prefactor * (-4.0L) * SUMO_CREAL(
                           lambda0sbotsbotc[i][k][m] * lambda0sbotsbotc[j][n][k] *
                           (Lsbot[m] * Rsbotc[n] + Rsbot[m] * Lsbotc[n])) *
                           m_bot * m_gluino * vSSSff;
  }}}

  for (k=0; k<2; k++) {

    TSIL_SetParametersSTU (&bar, m2_suL[k], m2_suL[k], 0.0L, m2_gluino, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);
    V_SSSFF (&gaak, &gaak, &vSSSFF, &vSSSff);

    quarksquarkgluinoresult += prefactor * 2.0L * 
      (lambda00suLsuLc[i][j] * W_SSFF (m2_suL[k], m2_suL[k], 0.0L, m2_gluino, Q2)
       + 2.0L*SUMO_CREAL(lambda0suLsuLc[i]*lambda0suLsuLc[j])*vSSSFF);

    TSIL_SetParametersSTU (&bar, m2_suR[k], m2_suR[k], 0.0L, m2_gluino, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);

    V_SSSFF (&gaak, &gaak, &vSSSFF, &vSSSff);

    quarksquarkgluinoresult += prefactor * 2.0L *
      (lambda00suRsuRc[i][j] * W_SSFF (m2_suR[k], m2_suR[k], 0.0L, m2_gluino, Q2)
       + 2.0L*SUMO_CREAL(lambda0suRsuRc[i]*lambda0suRsuRc[j])*vSSSFF);

    TSIL_SetParametersSTU (&bar, m2_sdL[k], m2_sdL[k], 0.0L, m2_gluino, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);

    V_SSSFF (&gaak, &gaak, &vSSSFF, &vSSSff);

    quarksquarkgluinoresult += prefactor * 2.0L *
      (lambda00sdLsdLc[i][j] * W_SSFF (m2_sdL[k], m2_sdL[k], 0.0L, m2_gluino, Q2)
       + 2.0L*SUMO_CREAL(lambda0sdLsdLc[i]*lambda0sdLsdLc[j])*vSSSFF);

    TSIL_SetParametersSTU (&bar, m2_sdR[k], m2_sdR[k], 0.0L, m2_gluino, Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);

    V_SSSFF (&gaak, &gaak, &vSSSFF, &vSSSff);

    quarksquarkgluinoresult +=  prefactor * 2.0L *
      (lambda00sdRsdRc[i][j] * W_SSFF (m2_sdR[k], m2_sdR[k], 0.0L, m2_gluino, Q2)
       + 2.0L*SUMO_CREAL(lambda0sdRsdRc[i]*lambda0sdRsdRc[j])*vSSSFF);
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%4.3lf sec\n", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
#endif

  /* ===== EQ. 4.5 =====*/

  /* stops */
  for (k=0; k<2; k++) {
  for (m=0; m<2; m++) {
  for (n=0; n<2; n++) {

     stopresult += prefactor * lambda00stopstopc[i][j][k][m] *
            (Lstopc[k] * Lstop[n] - Rstopc[k] * Rstop[n]) *
            (Lstop[m] * Lstopc[n] - Rstop[m] * Rstopc[n]) *
            X_SSS (m2_stop[k], m2_stop[m], m2_stop[n], Q2);

     for (p=0; p<2; p++) {
       stopresult += prefactor * (2.0L *
          SUMO_CREAL(lambda0stopstopc[i][k][m]*lambda0stopstopc[j][n][k]*
          (Lstop[m]*Lstopc[p] - Rstop[m]*Rstopc[p])*
          (Lstop[p]*Lstopc[n] - Rstop[p]*Rstopc[n]))*
          Y_SSSS(m2_stop[k], m2_stop[m], m2_stop[n], m2_stop[p], s, Q2)
	+ lambda0stopstopc[i][k][m]*lambda0stopstopc[j][n][p]*
	  (Lstop[m]*Lstopc[n] - Rstop[m]*Rstopc[n])*
	  (Lstop[p]*Lstopc[k] - Rstop[p]*Rstopc[k])*
          Z_SSSS(m2_stop[k], m2_stop[m], m2_stop[n], m2_stop[p], s, Q2));
      }
  }}}

  /* sbots */
  for (k=0; k<2; k++) {
  for (m=0; m<2; m++) {
  for (n=0; n<2; n++) {

     sbotresult += prefactor * lambda00sbotsbotc[i][j][k][m]*
          (Lsbotc[k]*Lsbot[n] - Rsbotc[k]*Rsbot[n])*
	  (Lsbot[m]*Lsbotc[n] - Rsbot[m]*Rsbotc[n])*
	  X_SSS (m2_sbot[k], m2_sbot[m], m2_sbot[n], Q2);
     for (p=0; p<2; p++) { 
       sbotresult += prefactor * (2.0L *
	  SUMO_CREAL(lambda0sbotsbotc[i][k][m]*lambda0sbotsbotc[j][n][k]*
          (Lsbot[m]*Lsbotc[p] - Rsbot[m]*Rsbotc[p])*
          (Lsbot[p]*Lsbotc[n] - Rsbot[p]*Rsbotc[n]))*
	  Y_SSSS(m2_sbot[k], m2_sbot[m], m2_sbot[n], m2_sbot[p], s, Q2)
	  + lambda0sbotsbotc[i][k][m]*lambda0sbotsbotc[j][n][p]*
	  (Lsbot[m]*Lsbotc[n] - Rsbot[m]*Rsbotc[n])*
	  (Lsbot[p]*Lsbotc[k] - Rsbot[p]*Rsbotc[k])*
	  Z_SSSS(m2_sbot[k], m2_sbot[m], m2_sbot[n], m2_sbot[p], s, Q2));
     }
  }}}

  /* Last two terms (sums over 1st/2nd gen sqarks) */
  for (k=0; k<2; k++) {
    squarkresult += prefactor * (lambda00suLsuLc[i][j]
            * X_SSS (m2_suL[k], m2_suL[k], m2_suL[k], Q2)
          + lambda0suLsuLc[i]*lambda0suLsuLc[j]*
	 (2.0L*Y_SSSS (m2_suL[k], m2_suL[k], m2_suL[k], m2_suL[k], s, Q2)
	  + Z_SSSS (m2_suL[k], m2_suL[k], m2_suL[k], m2_suL[k], s, Q2)));

    squarkresult += prefactor * (lambda00suRsuRc[i][j]
            * X_SSS (m2_suR[k], m2_suR[k], m2_suR[k], Q2)
          + lambda0suRsuRc[i]*lambda0suRsuRc[j]*
	 (2.0L*Y_SSSS (m2_suR[k], m2_suR[k], m2_suR[k], m2_suR[k], s, Q2)
	  + Z_SSSS (m2_suR[k], m2_suR[k], m2_suR[k], m2_suR[k], s, Q2)));

    squarkresult += prefactor * (lambda00sdLsdLc[i][j]
            * X_SSS (m2_sdL[k], m2_sdL[k], m2_sdL[k], Q2)
          + lambda0sdLsdLc[i]*lambda0sdLsdLc[j]*
	 (2.0L*Y_SSSS (m2_sdL[k], m2_sdL[k], m2_sdL[k], m2_sdL[k], s, Q2)
	  + Z_SSSS (m2_sdL[k], m2_sdL[k], m2_sdL[k], m2_sdL[k], s, Q2)));

    squarkresult += prefactor * (lambda00sdRsdRc[i][j]
            * X_SSS (m2_sdR[k], m2_sdR[k], m2_sdR[k], Q2)
          + lambda0sdRsdRc[i]*lambda0sdRsdRc[j]*
	 (2.0L*Y_SSSS (m2_sdR[k], m2_sdR[k], m2_sdR[k], m2_sdR[k], s, Q2)
	  + Z_SSSS (m2_sdR[k], m2_sdR[k], m2_sdR[k], m2_sdR[k], s, Q2)));
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%4.3lf sec\n", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
#endif

  totalresult = topgluonresult + botgluonresult + 
                stopgluonresult + sbotgluonresult + squarkgluonresult + 
                topstopgluinoresult + botsbotgluinoresult + quarksquarkgluinoresult + 
                stopresult + sbotresult + squarkresult;

  if (1 == VERBOSE) {
  printf("\nTotal QCD part of Higgs-self energy = \n%.10f + %.10f I\n", (double) SUMO_CREAL(totalresult*SUMO_twoloopfactor),(double) SUMO_CIMAG(totalresult*SUMO_twoloopfactor) );
  printf("  Breakdown:\n");
  printf("    top/gluon =            %.10f + %.10fI;\n", (double) SUMO_CREAL(topgluonresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(topgluonresult*SUMO_twoloopfactor));
  printf("    bot/gluon =            %.10f + %.10fI;\n", (double) SUMO_CREAL(botgluonresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(botgluonresult*SUMO_twoloopfactor));
  printf("    stop/gluon =           %.10f + %.10fI;\n", (double) SUMO_CREAL(stopgluonresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(stopgluonresult*SUMO_twoloopfactor));
  printf("    sbot/gluon =           %.10f + %.10fI;\n", (double) SUMO_CREAL(sbotgluonresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(sbotgluonresult*SUMO_twoloopfactor));
  printf("    squarks/gluon =        %.10f + %.10fI;\n", (double) SUMO_CREAL(squarkgluonresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(squarkgluonresult*SUMO_twoloopfactor));
  printf("    top/stop/gluino =      %.10f + %.10fI;\n", (double) SUMO_CREAL(topstopgluinoresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(topstopgluinoresult*SUMO_twoloopfactor));
  printf("    bot/sbot/gluino =      %.10f + %.10fI;\n", (double) SUMO_CREAL(botsbotgluinoresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(botsbotgluinoresult*SUMO_twoloopfactor));
  printf("    quark/squark/gluino =  %.10f + %.10fI;\n", (double) SUMO_CREAL(quarksquarkgluinoresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(quarksquarkgluinoresult*SUMO_twoloopfactor));
  printf("    stops =                %.10f + %.10fI;\n", (double) SUMO_CREAL(stopresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(stopresult*SUMO_twoloopfactor));
  printf("    sbots =                %.10f + %.10fI;\n", (double) SUMO_CREAL(sbotresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(sbotresult*SUMO_twoloopfactor));
  printf("    squarks =              %.10f + %.10fI;\n\n", (double) SUMO_CREAL(squarkresult*SUMO_twoloopfactor), (double) SUMO_CIMAG(squarkresult*SUMO_twoloopfactor));
  }

  return totalresult;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

int printc2 (SUMO_COMPLEX argu) {
  printf("%.10f + %.10f I",
	 (double) SUMO_CREAL(SUMO_twoloopfactor * argu),
	 (double) SUMO_CIMAG(SUMO_twoloopfactor * argu));
  return 0;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Two-loop self energy function for neutral Higgs particles evaluated
   at the specified squared momentum; non-QCD contributions
   (hep-ph/0405022, eqs. 4.6-4.17).
*/

SUMO_COMPLEX pi2nonQCD_phi0 (int i, int j, SUMO_REAL s)
{
  char funcname[] = "pi2nonQCD_phi0";
  TSIL_DATA bar;
  TSIL_RESULT gaak;
  TSIL_RESULT snutauData[2][2], stauData[2][2];
  TSIL_RESULT gaak2;
  TSIL_RESULT topNeutData[2][2][4], botNeutData[2][2][4], tauNeutData[2][2][4];
  TSIL_RESULT topCharData[2][2][4], botCharData[2][2][4], tauCharData[4];
  TSIL_RESULT phipData[2], phipData2[2];
  TSIL_RESULT eq15Data[2][2][2][2], eq16Data[2][2][2][2];
  TSIL_COMPLEX mFFFFS, mFFffS, mFfFfS, mFffFS, mffffS;
  /* TSIL_COMPLEX mFfFfSp; */
  TSIL_COMPLEX fac1, fac2;
  TSIL_COMPLEX sSSS, uSSSS, mSSSSS, vSSSSS;
  TSIL_COMPLEX vFFFFS, vFffFS, vfFfFS, vfFFfS, vFFffS, vffffS;
  TSIL_COMPLEX mSfSff, mSFSfF, mSFSFf;
  TSIL_COMPLEX vSSSFF, vSSSff;
  int k, m, n, p, q;
  TSIL_REAL m2_phi0_temp[4];
  TSIL_REAL m2_phip_temp[2];

  SUMO_COMPLEX topstopNresult = 0.0L;
  SUMO_COMPLEX botsbotNresult = 0.0L;

  SUMO_COMPLEX taustauNresult;
  SUMO_COMPLEX taustauNMresult = 0.0L;
  SUMO_COMPLEX taustauNVFFFFSresult = 0.0L;
  SUMO_COMPLEX taustauNWresult = 0.0L;
  SUMO_COMPLEX taustauNVSSSFFresult = 0.0L;

  SUMO_COMPLEX nutausnutauNresult = 0.0L;
  SUMO_COMPLEX topsbotCresult = 0.0L;  
  SUMO_COMPLEX botstopCresult = 0.0L; 
  SUMO_COMPLEX tausnutauCresult = 0.0L;
  SUMO_COMPLEX nutaustauCresult = 0.0L;
  SUMO_COMPLEX stopphi0result = 0.0L; 
  SUMO_COMPLEX sbotphi0result = 0.0L; 

  SUMO_COMPLEX stauphi0result; 
  SUMO_COMPLEX stauphi0Sresult = 0.0L; 
  SUMO_COMPLEX stauphi0Xresult = 0.0L; 
  SUMO_COMPLEX stauphi0Uresult = 0.0L; 
  SUMO_COMPLEX stauphi0Wresult = 0.0L; 
  SUMO_COMPLEX stauphi0Yresult = 0.0L; 
  SUMO_COMPLEX stauphi0Mresult = 0.0L; 
  SUMO_COMPLEX stauphi0Vresult = 0.0L; 

  SUMO_COMPLEX snutauphi0result = 0.0L;
  SUMO_COMPLEX topphi0result = 0.0L;
  SUMO_COMPLEX botphi0result = 0.0L; 
  SUMO_COMPLEX tauphi0result = 0.0L;

  SUMO_COMPLEX stopsbotphipresult;
  SUMO_COMPLEX sbotstopphipVresult = 0;
  SUMO_COMPLEX stopsbotphipVresult = 0;
  SUMO_COMPLEX stopsbotphipMresult = 0;
  SUMO_COMPLEX stopsbotphipSresult = 0;
  SUMO_COMPLEX stopphipXresult = 0;
  SUMO_COMPLEX sbotphipXresult = 0;
  SUMO_COMPLEX stopsbotphipUresult = 0;
  SUMO_COMPLEX sbotstopphipUresult = 0;
  SUMO_COMPLEX stopsbotphipWresult = 0;
  SUMO_COMPLEX sbotstopphipWresult = 0;
  SUMO_COMPLEX stopphipYresult = 0;
  SUMO_COMPLEX sbotphipYresult = 0;

  SUMO_COMPLEX snutaustauphipresult;
  SUMO_COMPLEX snutaustauphipSresult = 0;
  SUMO_COMPLEX snutauphipXresult = 0;
  SUMO_COMPLEX stauphipXresult = 0;
  SUMO_COMPLEX snutauphipYresult = 0;
  SUMO_COMPLEX stauphipYresult = 0;
  SUMO_COMPLEX snutaustauphipUresult = 0;
  SUMO_COMPLEX stausnutauphipUresult = 0;
  SUMO_COMPLEX snutaustauphipVresult = 0;
  SUMO_COMPLEX stausnutauphipVresult = 0;
  SUMO_COMPLEX snutaustauphipWresult = 0;
  SUMO_COMPLEX stausnutauphipWresult = 0;
  SUMO_COMPLEX snutaustauphipMresult = 0;

  SUMO_COMPLEX topbotphipresult = 0.0L;
  SUMO_COMPLEX nutautauphipresult = 0.0L;

  SUMO_COMPLEX Xresult;
  SUMO_COMPLEX Yresult;
  SUMO_COMPLEX Zresult;
  SUMO_COMPLEX XYZresult;

  SUMO_COMPLEX Xstopresult = 0;
  SUMO_COMPLEX Xstopsbotresult = 0;
  SUMO_COMPLEX Xstopstauresult = 0;
  SUMO_COMPLEX Xstopsnutauresult = 0; 
  SUMO_COMPLEX Xsbotstopresult = 0; 
  SUMO_COMPLEX Xsbotresult = 0; 
  SUMO_COMPLEX Xsbotstauresult = 0; 
  SUMO_COMPLEX Xsbotsnutauresult = 0; 
  SUMO_COMPLEX Xstaustopresult = 0; 
  SUMO_COMPLEX Xstausbotresult = 0; 
  SUMO_COMPLEX Xstauresult  = 0;
  SUMO_COMPLEX Xstausnutauresult = 0; 
  SUMO_COMPLEX Xsnutaustopresult = 0;
  SUMO_COMPLEX Xsnutausbotresult = 0;
  SUMO_COMPLEX Xsnutaustauresult = 0;
  SUMO_COMPLEX Xsnutauresult = 0;

  SUMO_COMPLEX Ystopresult = 0;
  SUMO_COMPLEX Ystopsbotresult = 0;
  SUMO_COMPLEX Ystopstauresult = 0;
  SUMO_COMPLEX Ystopsnutauresult = 0; 
  SUMO_COMPLEX Ysbotstopresult = 0; 
  SUMO_COMPLEX Ysbotresult = 0; 
  SUMO_COMPLEX Ysbotstauresult = 0; 
  SUMO_COMPLEX Ysbotsnutauresult = 0; 
  SUMO_COMPLEX Ystaustopresult = 0; 
  SUMO_COMPLEX Ystausbotresult = 0; 
  SUMO_COMPLEX Ystauresult  = 0;
  SUMO_COMPLEX Ystausnutauresult = 0; 
  SUMO_COMPLEX Ysnutaustopresult = 0;
  SUMO_COMPLEX Ysnutausbotresult = 0;
  SUMO_COMPLEX Ysnutaustauresult = 0;
  SUMO_COMPLEX Ysnutauresult = 0;

  SUMO_COMPLEX Zstopresult = 0;
  SUMO_COMPLEX Zstopsbotresult = 0;
  SUMO_COMPLEX Zstopstauresult = 0;
  SUMO_COMPLEX Zstopsnutauresult = 0; 
  SUMO_COMPLEX Zsbotresult = 0; 
  SUMO_COMPLEX Zsbotstauresult = 0; 
  SUMO_COMPLEX Zsbotsnutauresult = 0; 
  SUMO_COMPLEX Zstauresult  = 0;
  SUMO_COMPLEX Zstausnutauresult = 0; 
  SUMO_COMPLEX Zsnutauresult = 0;

  SUMO_COMPLEX totalresult;

#ifdef SUMO_HIGGS_TIMING
  clock_t t0, t1;
  t0 = clock ();
#endif

  /* Need non-negative versions of the Higgs squared masses. */
  for (k=0; k<4; k++) {
    m2_phi0_temp[k] = m2_phi0[k];
    if (m2_phi0_temp[k] < 0) m2_phi0_temp[k] = 0.0L;
  }

  for (k=0; k<2; k++) {
    m2_phip_temp[k] = m2_phip[k];
    if (m2_phip_temp[k] < 0) m2_phip_temp[k] = 0.0L;
  }

  /* ===== hep-ph/0405022 EQ. 4.6 ===== */

  /* Perform needed TSIL evaluations */
  for (m=0; m<2; m++)
  for (n=m; n<2; n++)
  for (k=0; k<4; k++) {

    TSIL_SetParameters (&bar, 
			m2_stop[m], m2_top, m2_stop[n], m2_top, m2_Neut[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &topNeutData[m][n][k]);

    TSIL_SetParameters (&bar,
			m2_sbot[m], m2_bot, m2_sbot[n], m2_bot, m2_Neut[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &botNeutData[m][n][k]);

    TSIL_SetParameters (&bar,
			m2_stau[m], m2_tau, m2_stau[n], m2_tau, m2_Neut[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &tauNeutData[m][n][k]);

    if (n != m) {
      TSIL_PermuteResult (&topNeutData[m][n][k], XZandYU, &topNeutData[n][m][k]);
      TSIL_PermuteResult (&botNeutData[m][n][k], XZandYU, &botNeutData[n][m][k]);
      TSIL_PermuteResult (&tauNeutData[m][n][k], XZandYU, &tauNeutData[n][m][k]);
    }
  }


  for (m=0; m<2; m++)
  for (n=0; n<2; n++)
  for (k=0; k<4; k++) {

    M_SFSFF (&topNeutData[m][n][k], &mSfSff, &mSFSfF, &mSFSFf);

    topstopNresult += Ncolors * 2.0L * (SUMO_CREAL( (Yttcphi[i] * lambda0stopstopc[j][m][n]
                   + Yttcphi[j] * lambda0stopstopc[i][m][n] )
                 * conYtNstopc[k][n] * conYtcNstop[k][m]) * m_Neut[k] * mSFSFf
      + SUMO_CREAL( (Yttcphi[i] * lambda0stopstopc[j][n][m]
                   + Yttcphi[j] * lambda0stopstopc[i][n][m])
        * (conYtNstopc[k][m] * YtNstopc[k][n] + conYtcNstop[k][n] * YtcNstop[k][m]) ) * m_top * mSFSfF
      + SUMO_CREAL( (Yttcphi[i] * lambda0stopstopc[j][m][n]
                   + Yttcphi[j] * lambda0stopstopc[i][m][n])
                 * YtNstopc[k][m] * YtcNstop[k][n] ) * m2_top * m_Neut[k] * mSfSff);

    /* (s)bottom terms */
    M_SFSFF (&botNeutData[m][n][k], &mSfSff, &mSFSfF, &mSFSFf);

    botsbotNresult += Ncolors * 2.0L * (SUMO_CREAL( (Ybbcphi[i] * lambda0sbotsbotc[j][m][n]
		   + Ybbcphi[j] * lambda0sbotsbotc[i][m][n] )
		 * conYbNsbotc[k][n] * conYbcNsbot[k][m]) * m_Neut[k] * mSFSFf
      + SUMO_CREAL( (Ybbcphi[i] * lambda0sbotsbotc[j][n][m]
		   + Ybbcphi[j] * lambda0sbotsbotc[i][n][m])
	* (conYbNsbotc[k][m] * YbNsbotc[k][n] + conYbcNsbot[k][n] * YbcNsbot[k][m]) ) * m_bot * mSFSfF
      + SUMO_CREAL( (Ybbcphi[i] * lambda0sbotsbotc[j][m][n]
		    + Ybbcphi[j] * lambda0sbotsbotc[i][m][n])
		 * YbNsbotc[k][m] * YbcNsbot[k][n] ) * m2_bot * m_Neut[k] * mSfSff);
  }

  /* (s)tau */
  for (m=0; m<2; m++)
  for (n=0; n<2; n++)
  for (k=0; k<4; k++) {

    M_SFSFF (&tauNeutData[m][n][k], &mSfSff, &mSFSfF, &mSFSFf);

    taustauNMresult += 2.0L * (SUMO_CREAL((Ytautaucphi[i]*lambda0staustauc[j][m][n]
                      + Ytautaucphi[j]*lambda0staustauc[i][m][n])
		       *conYtauNstauc[k][n]*conYtaucNstau[k][m]) *m_Neut[k]*mSFSFf
      + SUMO_CREAL((Ytautaucphi[i]*lambda0staustauc[j][n][m]
		  + Ytautaucphi[j]*lambda0staustauc[i][n][m])
		   *(conYtauNstauc[k][m]*YtauNstauc[k][n]
		     + conYtaucNstau[k][n]*YtaucNstau[k][m])) *m_tau*mSFSfF
      + SUMO_CREAL((Ytautaucphi[i]*lambda0staustauc[j][m][n]
		    + Ytautaucphi[j]*lambda0staustauc[i][m][n])
		   *YtauNstauc[k][m]*YtaucNstau[k][n]) *m2_tau*m_Neut[k]*mSfSff);
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.7 ===== */
  /* Correct; checked 6/7/12 */

  /* TSIL evaluation */
  for (m=0; m<2; m++)
  for (n=m; n<2; n++)
  for (k=0; k<2; k++) {

    TSIL_SetParameters (&bar, 
			m2_sbot[m], m2_top, m2_sbot[n], m2_top, m2_Char[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &topCharData[m][n][k]);

    TSIL_SetParameters (&bar, 
			m2_stop[m], m2_bot, m2_stop[n],	m2_bot, m2_Char[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &botCharData[m][n][k]);

    if (n != m) {
      TSIL_PermuteResult (&topCharData[m][n][k], XZandYU, &topCharData[n][m][k]);
      TSIL_PermuteResult (&botCharData[m][n][k], XZandYU, &botCharData[n][m][k]);
    }}

  for (m=0; m<2; m++)
  for (n=0; n<2; n++)
  for (k=0; k<2; k++) {

    M_SFSFF (&topCharData[m][n][k], &mSfSff, &mSFSfF, &mSFSFf);

    topsbotCresult += Ncolors * 2.0L * (SUMO_CREAL( (Yttcphi[i] * lambda0sbotsbotc[j][m][n]
       	           + Yttcphi[j] * lambda0sbotsbotc[i][m][n])
	     * conYtCsbotc[k][n] * conYtcCsbot[k][m] ) * m_Char[k] * mSFSFf
      + SUMO_CREAL( (Yttcphi[i] * lambda0sbotsbotc[j][n][m]
	           + Yttcphi[j] * lambda0sbotsbotc[i][n][m])
	     * (conYtCsbotc[k][m] * YtCsbotc[k][n] + conYtcCsbot[k][n] * YtcCsbot[k][m]) )
             * m_top * mSFSfF
      + SUMO_CREAL( (Yttcphi[i] * lambda0sbotsbotc[j][m][n]
		   + Yttcphi[j] * lambda0sbotsbotc[i][m][n])
             * YtCsbotc[k][m] * YtcCsbot[k][n] )
             * m2_top * m_Char[k] * mSfSff);

    /* t <-> b terms */
    M_SFSFF (&botCharData[m][n][k], &mSfSff, &mSFSfF, &mSFSFf);

    botstopCresult += Ncolors * 2.0L * (SUMO_CREAL( (Ybbcphi[i] * lambda0stopstopc[j][m][n]
		   + Ybbcphi[j] * lambda0stopstopc[i][m][n])
	     * conYbCstopc[k][n] * conYbcCstop[k][m]) * m_Char[k] * mSFSFf
      + SUMO_CREAL( (Ybbcphi[i]*lambda0stopstopc[j][n][m]
		   + Ybbcphi[j]*lambda0stopstopc[i][n][m])
	      * (conYbCstopc[k][m] * YbCstopc[k][n] + conYbcCstop[k][n] * YbcCstop[k][m]))
              * m_bot * mSFSfF
      + SUMO_CREAL( (Ybbcphi[i]*lambda0stopstopc[j][m][n]
		   + Ybbcphi[j]*lambda0stopstopc[i][m][n])
	      * YbCstopc[k][m] * YbcCstop[k][n]) * m2_bot * m_Char[k] * mSfSff);
  }

  /* Last three lines: */
  for (k=0; k<2; k++) {
    
    TSIL_SetParameters (&bar,
			m2_snu[2], m2_tau, m2_snu[2], m2_tau, m2_Char[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &tauCharData[k]);

    M_SFSFF (&tauCharData[k], &mSfSff, &mSFSfF, &mSFSFf);

    tausnutauCresult += 2.0 * (SUMO_CREAL((Ytautaucphi[i] * lambda0snusnuc[j]
                  + Ytautaucphi[j] * lambda0snusnuc[i])
	      * conYtauCsnutauc[k] * conYtaucCsnutau[k]) * m_Char[k] * mSFSFf
      + SUMO_CREAL((Ytautaucphi[i] * lambda0snusnuc[j]
                  + Ytautaucphi[j] * lambda0snusnuc[i]))
              * (YtauCsnutauc[k] * SUMO_CONJ(YtauCsnutauc[k]) 
	       + YtaucCsnutau[k] * SUMO_CONJ(YtaucCsnutau[k])) * m_tau * mSFSfF
      + SUMO_CREAL((Ytautaucphi[i]*lambda0snusnuc[j]
		  + Ytautaucphi[j]*lambda0snusnuc[i])
	      * YtauCsnutauc[k] * conYtaucCsnutau[k]) * m2_tau * m_Char[k] * mSfSff);
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.8 ===== */
  for (m=0; m<2; m++)
  for (k=0; k<4; k++) {

    /* The same reorganization here as in term 4.3! */
    gaak.z = topNeutData[m][m][k].x;
    gaak.x = topNeutData[m][m][k].y;
    gaak.y = topNeutData[m][m][k].u;
    gaak.u = topNeutData[m][m][k].v;
    gaak.v = topNeutData[m][m][k].z;
    gaak.z = topNeutData[m][m][k].y;
    gaak.s = topNeutData[m][m][k].s;
    gaak.qq = topNeutData[m][m][k].qq;
    gaak.B[xz] = topNeutData[m][m][k].B[yu];
    gaak.U[xzuv] = topNeutData[m][m][k].U[yuzv];
    gaak.V[xzuv] = topNeutData[m][m][k].V[yuzv];
    gaak.T[xuv] = topNeutData[m][m][k].T[uxv];
    gaak.T[uxv] = topNeutData[m][m][k].T[vxu];
    gaak.T[vxu] = topNeutData[m][m][k].T[zyv];
    gaak.S[uxv] = topNeutData[m][m][k].S[uxv];

    /* TRY UNCLEVER TSIL EVAL -- DOES NOT AFFECT RESULTS */
/*     TSIL_SetParameters (&bar, */
/* 			m2_top, m2_top, m2_top, m2_Neut[k], m2_stop[m], Q2); */
/*     TSIL_Evaluate (&bar, s); */
/*     TSIL_CopyResult (&bar, &gaak); */
    /* END OF TEST CODE */

    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    topstopNresult += Ncolors * 2.0L *
      ((YtNstopc[k][m] * conYtNstopc[k][m] + YtcNstop[k][m] * conYtcNstop[k][m])
      * ( SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) * (vFFFFS + m2_top * vFffFS)
      + 2.0L * SUMO_CREAL(Yttcphi[i] * Yttcphi[j]) * m2_top * vfFfFS )
      + 2.0L * m_top * m_Neut[k]
      * ( SUMO_CREAL(conYttcphi[i] * conYttcphi[j] * YtNstopc[k][m] * YtcNstop[k][m]) * vfFFfS
      + 2.0L * SUMO_CREAL(Yttcphi[i] * conYttcphi[j])
	  * SUMO_CREAL(YtNstopc[k][m] * YtcNstop[k][m]) * vFFffS
      + SUMO_CREAL(Yttcphi[i] * Yttcphi[j] * YtNstopc[k][m] * YtcNstop[k][m])
	  * m2_top * vffffS ));

    /* ...and t -> b: */
    gaak.z = botNeutData[m][m][k].x;
    gaak.x = botNeutData[m][m][k].y;
    gaak.y = botNeutData[m][m][k].u;
    gaak.u = botNeutData[m][m][k].v;
    gaak.v = botNeutData[m][m][k].z;
    gaak.z = botNeutData[m][m][k].y;
    gaak.s = botNeutData[m][m][k].s;
    gaak.qq = botNeutData[m][m][k].qq;
    gaak.B[xz] = botNeutData[m][m][k].B[yu];
    gaak.U[xzuv] = botNeutData[m][m][k].U[yuzv];
    gaak.V[xzuv] = botNeutData[m][m][k].V[yuzv];
    gaak.T[xuv] = botNeutData[m][m][k].T[uxv];
    gaak.T[uxv] = botNeutData[m][m][k].T[vxu];
    gaak.T[vxu] = botNeutData[m][m][k].T[zyv];
    gaak.S[uxv] = botNeutData[m][m][k].S[uxv];

    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    botsbotNresult += Ncolors * 2.0L *
      ((YbNsbotc[k][m]*conYbNsbotc[k][m] + YbcNsbot[k][m]*conYbcNsbot[k][m]) *
      (SUMO_CREAL(Ybbcphi[i]*conYbbcphi[j])*(vFFFFS + m2_bot*vFffFS)
       + 2.0L*SUMO_CREAL(Ybbcphi[i]*Ybbcphi[j])*m2_bot*vfFfFS)
      + 2.0L*m_bot*m_Neut[k]*
      (SUMO_CREAL(conYbbcphi[i]*conYbbcphi[j]*YbNsbotc[k][m]*YbcNsbot[k][m])*vfFFfS
       + 2.0L*SUMO_CREAL(Ybbcphi[i]*conYbbcphi[j])*
       SUMO_CREAL(YbNsbotc[k][m]*YbcNsbot[k][m])*vFFffS
       + SUMO_CREAL(Ybbcphi[i]*Ybbcphi[j]*YbNsbotc[k][m]*YbcNsbot[k][m])
       *m2_bot*vffffS));
  }

  /* t -> tau */
  for (m=0; m<2; m++)
  for (k=0; k<4; k++) {

    gaak.z = tauNeutData[m][m][k].x;
    gaak.x = tauNeutData[m][m][k].y;
    gaak.y = tauNeutData[m][m][k].u;
    gaak.u = tauNeutData[m][m][k].v;
    gaak.v = tauNeutData[m][m][k].z;
    gaak.z = tauNeutData[m][m][k].y;
    gaak.s = tauNeutData[m][m][k].s;
    gaak.qq = tauNeutData[m][m][k].qq;
    gaak.B[xz] = tauNeutData[m][m][k].B[yu];
    gaak.U[xzuv] = tauNeutData[m][m][k].U[yuzv];
    gaak.V[xzuv] = tauNeutData[m][m][k].V[yuzv];
    gaak.T[xuv] = tauNeutData[m][m][k].T[uxv];
    gaak.T[uxv] = tauNeutData[m][m][k].T[vxu];
    gaak.T[vxu] = tauNeutData[m][m][k].T[zyv];
    gaak.S[uxv] = tauNeutData[m][m][k].S[uxv];
    
    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    taustauNVFFFFSresult += 2.0L *
      ((YtauNstauc[k][m] * conYtauNstauc[k][m] + YtaucNstau[k][m] * conYtaucNstau[k][m]) * 
      (SUMO_CREAL(Ytautaucphi[i] * conYtautaucphi[j]) * (vFFFFS + m2_tau * vFffFS)
       + 2.0L * SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j]) * m2_tau * vfFfFS)
      + 2.0L * m_tau * m_Neut[k] * (
       SUMO_CREAL(conYtautaucphi[i] * conYtautaucphi[j] * YtauNstauc[k][m] * YtaucNstau[k][m]) * vfFFfS
       + 2.0L * SUMO_CREAL(Ytautaucphi[i] * conYtautaucphi[j]) * 
       SUMO_CREAL(YtauNstauc[k][m] * YtaucNstau[k][m]) * vFFffS
       + SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j] * YtauNstauc[k][m] * YtaucNstau[k][m])
        * m2_tau * vffffS));
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.9 ===== */
  /* Correct; checked 6/7/12 */

  for (m=0; m<2; m++)
  for (k=0; k<2; k++) {

    /* The same reorganization here as in term 4.3! */
    gaak.z = topCharData[m][m][k].x;
    gaak.x = topCharData[m][m][k].y;
    gaak.y = topCharData[m][m][k].u;
    gaak.u = topCharData[m][m][k].v;
    gaak.v = topCharData[m][m][k].z;
    gaak.z = topCharData[m][m][k].y;
    gaak.s = topCharData[m][m][k].s;
    gaak.qq = topCharData[m][m][k].qq;
    gaak.B[xz] = topCharData[m][m][k].B[yu];
    gaak.U[xzuv] = topCharData[m][m][k].U[yuzv];
    gaak.V[xzuv] = topCharData[m][m][k].V[yuzv];
    gaak.T[xuv] = topCharData[m][m][k].T[uxv];
    gaak.T[uxv] = topCharData[m][m][k].T[vxu];
    gaak.T[vxu] = topCharData[m][m][k].T[zyv];
    gaak.S[uxv] = topCharData[m][m][k].S[uxv];
    
    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    topsbotCresult += Ncolors * 2.0L *
      ((YtCsbotc[k][m] * conYtCsbotc[k][m] + YtcCsbot[k][m] * conYtcCsbot[k][m])
      * (SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) * (vFFFFS + m2_top * vFffFS)
      + 2.0L * SUMO_CREAL(Yttcphi[i] * Yttcphi[j]) * m2_top * vfFfFS)
      + 2.0L * m_top * m_Char[k]
      * ( SUMO_CREAL(conYttcphi[i] * conYttcphi[j] * YtCsbotc[k][m] * YtcCsbot[k][m]) * vfFFfS
         + 2.0L * SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) 
	 * SUMO_CREAL(YtCsbotc[k][m] * YtcCsbot[k][m]) * vFFffS
	 + SUMO_CREAL(Yttcphi[i] * Yttcphi[j] * YtCsbotc[k][m] * YtcCsbot[k][m])
	 * m2_top * vffffS ));
    
    /* ...and t -> b: */
    gaak.z = botCharData[m][m][k].x;
    gaak.x = botCharData[m][m][k].y;
    gaak.y = botCharData[m][m][k].u;
    gaak.u = botCharData[m][m][k].v;
    gaak.v = botCharData[m][m][k].z;
    gaak.z = botCharData[m][m][k].y;
    gaak.s = botCharData[m][m][k].s;
    gaak.qq = botCharData[m][m][k].qq;
    gaak.B[xz] = botCharData[m][m][k].B[yu];
    gaak.U[xzuv] = botCharData[m][m][k].U[yuzv];
    gaak.V[xzuv] = botCharData[m][m][k].V[yuzv];
    gaak.T[xuv] = botCharData[m][m][k].T[uxv];
    gaak.T[uxv] = botCharData[m][m][k].T[vxu];
    gaak.T[vxu] = botCharData[m][m][k].T[zyv];
    gaak.S[uxv] = botCharData[m][m][k].S[uxv];
    
    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    botstopCresult += Ncolors * 2.0L *
      ((YbCstopc[k][m] * conYbCstopc[k][m] + YbcCstop[k][m] * conYbcCstop[k][m])
      * (SUMO_CREAL(Ybbcphi[i] * conYbbcphi[j]) * (vFFFFS + m2_bot * vFffFS)
      + 2.0L * SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j]) * m2_bot * vfFfFS)
      + 2.0L * m_bot * m_Char[k]
      * ( SUMO_CREAL(conYbbcphi[i] * conYbbcphi[j] * YbCstopc[k][m] * YbcCstop[k][m])
	  * vfFFfS
	  + 2.0L * SUMO_CREAL(Ybbcphi[i] * conYbbcphi[j])
	  * SUMO_CREAL(YbCstopc[k][m] * YbcCstop[k][m]) * vFFffS
	  + SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j] * YbCstopc[k][m] * YbcCstop[k][m])
	  * m2_bot * vffffS ));
  }
  
  /* t <-> tau */
  for (k=0; k<2; k++) {
    
    gaak.z = tauCharData[k].x;
    gaak.x = tauCharData[k].y;
    gaak.y = tauCharData[k].u;
    gaak.u = tauCharData[k].v;
    gaak.v = tauCharData[k].z;
    gaak.z = tauCharData[k].y;
    gaak.s = tauCharData[k].s;
    gaak.qq = tauCharData[k].qq;
    gaak.B[xz] = tauCharData[k].B[yu];
    gaak.U[xzuv] = tauCharData[k].U[yuzv];
    gaak.V[xzuv] = tauCharData[k].V[yuzv];
    gaak.T[xuv] = tauCharData[k].T[uxv];
    gaak.T[uxv] = tauCharData[k].T[vxu];
    gaak.T[vxu] = tauCharData[k].T[zyv];
    gaak.S[uxv] = tauCharData[k].S[uxv];

    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    tausnutauCresult += 2.0L *
      ((YtauCsnutauc[k] * conYtauCsnutauc[k] + YtaucCsnutau[k] * conYtaucCsnutau[k])
      * (SUMO_CREAL(Ytautaucphi[i] * conYtautaucphi[j]) * (vFFFFS + m2_tau * vFffFS)
      + 2.0L * SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j]) * m2_tau * vfFfFS)
      + 2.0L * m_tau * m_Char[k]
      * ( SUMO_CREAL(conYtautaucphi[i] * conYtautaucphi[j] * YtauCsnutauc[k] * YtaucCsnutau[k])
	  * vfFFfS
	  + 2.0L * SUMO_CREAL(Ytautaucphi[i] * conYtautaucphi[j])
	  * SUMO_CREAL(YtauCsnutauc[k] * YtaucCsnutau[k]) * vFFffS
	  + SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j] * YtauCsnutauc[k] * YtaucCsnutau[k])
	  * m2_tau * vffffS ));
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.10 ===== */
  for (m=0; m<2; m++)
  for (n=0; n<2; n++)
  for (k=0; k<4; k++) {

    topstopNresult += Ncolors * (lambda00stopstopc[i][j][m][n]
      * ( (YtNstopc[k][m]  * conYtNstopc[k][n] + conYtcNstop[k][m] * YtcNstop[k][n])
      * W_SSFF (m2_stop[m], m2_stop[n], m2_top, m2_Neut[k], Q2)
      + (YtNstopc[k][m] * YtcNstop[k][n] + conYtcNstop[k][m] * conYtNstopc[k][n])
      * m_top * m_Neut[k] * W_SSff (m2_stop[m], m2_stop[n], m2_top, m2_Neut[k], Q2) ));

    botsbotNresult += Ncolors * (lambda00sbotsbotc[i][j][m][n] * 
      ((YbNsbotc[k][m] * conYbNsbotc[k][n] + conYbcNsbot[k][m] * YbcNsbot[k][n]) * 
       W_SSFF (m2_sbot[m],m2_sbot[n],m2_bot,m2_Neut[k],Q2)
       + (YbNsbotc[k][m] * YbcNsbot[k][n] + conYbcNsbot[k][m] * conYbNsbotc[k][n]) * 
       m_bot * m_Neut[k] * W_SSff (m2_sbot[m],m2_sbot[n],m2_bot,m2_Neut[k],Q2) ));

    for (p=0; p<2; p++) {

      V_SSSFF(&topNeutData[m][n][k], &topNeutData[m][p][k], &vSSSFF, &vSSSff);

      /* TRY UNCLEVER TSIL EVAL -- MAKES NO DIFFERENCE */
/*       TSIL_SetParameters (&bar, */
/* 			  m2_stop[m], m2_stop[n], m2_stop[p], m2_top, m2_Neut[k], Q2); */
/*       TSIL_Evaluate (&bar, s); */
/*       TSIL_CopyResult (&bar, &gaak); */

/*       TSIL_SetParameters (&bar, */
/* 			  m2_stop[m], m2_stop[p], m2_stop[n], m2_top, m2_Neut[k], Q2); */
/*       TSIL_Evaluate (&bar, s); */
/*       TSIL_CopyResult (&bar, &gaak2); */

/*       V_SSSFF (&gaak, &gaak2, &vSSSFF, &vSSSff); */
      /* END OF TEST CODE */

      topstopNresult += Ncolors *
	(2.0L * SUMO_CREAL( lambda0stopstopc[i][m][n] * lambda0stopstopc[j][p][m]
	* (conYtNstopc[k][n] * YtNstopc[k][p] + YtcNstop[k][n] * conYtcNstop[k][p]) ) * vSSSFF
	+ 2.0L * SUMO_CREAL( lambda0stopstopc[i][m][n] * lambda0stopstopc[j][p][m]
	* (YtcNstop[k][n] * YtNstopc[k][p] + conYtcNstop[k][p] * conYtNstopc[k][n]) )
	* m_top * m_Neut[k] * vSSSff);

      V_SSSFF(&botNeutData[m][n][k], &botNeutData[m][p][k], &vSSSFF, &vSSSff);

      botsbotNresult += Ncolors *
	(2.0L*SUMO_CREAL(lambda0sbotsbotc[i][m][n]*lambda0sbotsbotc[j][p][m]*
			      (conYbNsbotc[k][n]*YbNsbotc[k][p] +
			       YbcNsbot[k][n]*conYbcNsbot[k][p]))*vSSSFF
	+ 2.0L*SUMO_CREAL(lambda0sbotsbotc[i][m][n]*lambda0sbotsbotc[j][p][m]*
			  (YbcNsbot[k][n]*YbNsbotc[k][p] +
			   conYbcNsbot[k][p]*conYbNsbotc[k][n]))* m_bot*m_Neut[k]*vSSSff);
    }
  }

  /* t -> tau terms */
  for (m=0; m<2; m++)
  for (n=0; n<2; n++)
  for (k=0; k<4; k++) {

    taustauNWresult += lambda00staustauc[i][j][m][n] * 
      ((YtauNstauc[k][m] * conYtauNstauc[k][n] + conYtaucNstau[k][m] * YtaucNstau[k][n]) * 
       W_SSFF (m2_stau[m],m2_stau[n],m2_tau,m2_Neut[k],Q2)
       + (YtauNstauc[k][m] * YtaucNstau[k][n] + conYtaucNstau[k][m] * conYtauNstauc[k][n]) * 
       m_tau * m_Neut[k] * W_SSff (m2_stau[m],m2_stau[n],m2_tau,m2_Neut[k],Q2));

    for (p=0; p<2; p++) {

      V_SSSFF(&tauNeutData[m][n][k], &tauNeutData[m][p][k], &vSSSFF, &vSSSff);

      taustauNVSSSFFresult += 2.0L * SUMO_CREAL(lambda0staustauc[i][m][n] * lambda0staustauc[j][p][m] * 
			      (conYtauNstauc[k][n] * YtauNstauc[k][p] +
			       YtaucNstau[k][n] * conYtaucNstau[k][p])) * vSSSFF
	+ 2.0L * SUMO_CREAL(lambda0staustauc[i][m][n] * lambda0staustauc[j][p][m] * 
			  (YtaucNstau[k][n] * YtauNstauc[k][p] +
			   conYtaucNstau[k][p] * conYtauNstauc[k][n])) *  m_tau * m_Neut[k] * vSSSff;
    }
  }

  taustauNresult = taustauNMresult + taustauNVFFFFSresult + taustauNWresult + taustauNVSSSFFresult;

  /* Last line: (NOT INDEPENDENTLY CHECKED -- vanishes for g=gp=0 */
  for (k=0; k<4; k++) {

    TSIL_SetParametersSTU (&bar, m2_snu[2], m2_snu[2], 0.0L, m2_Neut[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);

    V_SSSFF (&gaak, &gaak, &vSSSFF, &vSSSff);

    nutausnutauNresult += YnuNsnuc[k] * conYnuNsnuc[k]
      * (2.0L * lambda0snusnuc[i] * lambda0snusnuc[j] * vSSSFF
	 + lambda00snusnuc[i][j] * W_SSFF (m2_snu[2], m2_snu[2], 0.0L, m2_Neut[k], Q2));

  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.11 ===== */

  for (m=0; m<2; m++)
  for (n=0; n<2; n++)
  for (k=0; k<2; k++) {

    topsbotCresult += Ncolors * (lambda00sbotsbotc[i][j][m][n]
      * ((YtCsbotc[k][m] * conYtCsbotc[k][n] + conYtcCsbot[k][m] * YtcCsbot[k][n])
      * W_SSFF (m2_sbot[m], m2_sbot[n], m2_top, m2_Char[k], Q2)
      + (YtCsbotc[k][m] * YtcCsbot[k][n] + conYtcCsbot[k][m] * conYtCsbotc[k][n])
      * m_top * m_Char[k] * W_SSff (m2_sbot[m], m2_sbot[n], m2_top, m2_Char[k], Q2)));

    botstopCresult += Ncolors * (lambda00stopstopc[i][j][m][n] *
      ((YbCstopc[k][m] * conYbCstopc[k][n] + conYbcCstop[k][m] * YbcCstop[k][n])
      * W_SSFF (m2_stop[m], m2_stop[n], m2_bot, m2_Char[k], Q2)
      + (YbCstopc[k][m] * YbcCstop[k][n] + conYbcCstop[k][m] * conYbCstopc[k][n])
      * m_bot * m_Char[k] * W_SSff (m2_stop[m], m2_stop[n], m2_bot, m2_Char[k], Q2)));

    for (p=0; p<2; p++) {

      V_SSSFF(&topCharData[m][n][k], &topCharData[m][p][k], &vSSSFF, &vSSSff);

      topsbotCresult += Ncolors * (2.0L*SUMO_CREAL(lambda0sbotsbotc[i][m][n] * lambda0sbotsbotc[j][p][m]
        * (conYtCsbotc[k][n] * YtCsbotc[k][p] + YtcCsbot[k][n] * conYtcCsbot[k][p])) * vSSSFF
	+ 2.0L*SUMO_CREAL(lambda0sbotsbotc[i][m][n] * lambda0sbotsbotc[j][p][m]
	* (YtcCsbot[k][n] * YtCsbotc[k][p] + conYtcCsbot[k][p] * conYtCsbotc[k][n]))
	* m_top * m_Char[k] * vSSSff);

      V_SSSFF(&botCharData[m][n][k], &botCharData[m][p][k], &vSSSFF, &vSSSff);
      
      botstopCresult += Ncolors * (2.0L*SUMO_CREAL(lambda0stopstopc[i][m][n] * lambda0stopstopc[j][p][m]
	* (conYbCstopc[k][n] * YbCstopc[k][p] + YbcCstop[k][n] * conYbcCstop[k][p])) * vSSSFF
	+ 2.0L*SUMO_CREAL(lambda0stopstopc[i][m][n] * lambda0stopstopc[j][p][m]
        * (YbcCstop[k][n] * YbCstopc[k][p] + conYbcCstop[k][p] * conYbCstopc[k][n]))
	* m_bot * m_Char[k] * vSSSff);
    }}

  for (k=0; k<2; k++) {

    V_SSSFF(&tauCharData[k], &tauCharData[k], &vSSSFF, &vSSSff);

    tausnutauCresult += 2.0L * lambda0snusnuc[i] * lambda0snusnuc[j]
      * (YtauCsnutauc[k] * conYtauCsnutauc[k] + YtaucCsnutau[k] * conYtaucCsnutau[k]) * vSSSFF
      + 4.0L * lambda0snusnuc[i] * lambda0snusnuc[j]
      * SUMO_CREAL(YtaucCsnutau[k] * YtauCsnutauc[k]) * m_tau * m_Char[k] * vSSSff
      + lambda00snusnuc[i][j]
      * ((YtauCsnutauc[k] * conYtauCsnutauc[k] + YtaucCsnutau[k] * conYtaucCsnutau[k])
      * W_SSFF (m2_snu[2], m2_snu[2], m2_tau, m2_Char[k], Q2)
      + 2.0L * SUMO_CREAL(YtauCsnutauc[k] * YtaucCsnutau[k]) * m_tau * m_Char[k]
      * W_SSff (m2_snu[2], m2_snu[2], m2_tau, m2_Char[k], Q2));

    for (m=0; m<2; m++)
    for (n=0; n<2; n++) {

      /* First do p = n */
      TSIL_SetParametersSTU (&bar, m2_stau[m], m2_stau[n], 0.0L, m2_Char[k], Q2);
      TSIL_Evaluate (&bar, s);
      TSIL_CopyResult (&bar, &gaak);

      p = n;

      V_SSSFF (&gaak, &gaak, &vSSSFF, &vSSSff);

      nutaustauCresult += 2.0L * SUMO_CREAL(lambda0staustauc[i][m][n] * lambda0staustauc[j][p][m]
		    * conYnutauCstauc[k][n] * YnutauCstauc[k][p]) * vSSSFF;

      if (n == 0) p = 1; else p = 0;

      TSIL_SetParametersSTU (&bar, m2_stau[m], m2_stau[p], 0.0L, m2_Char[k], Q2);
      TSIL_Evaluate (&bar, s);
      TSIL_CopyResult (&bar, &gaak2);

      V_SSSFF (&gaak, &gaak2, &vSSSFF, &vSSSff);

      nutaustauCresult += 2.0L * SUMO_CREAL(lambda0staustauc[i][m][n] * lambda0staustauc[j][p][m]
		    * conYnutauCstauc[k][n] * YnutauCstauc[k][p]) * vSSSFF;

      nutaustauCresult += lambda00staustauc[i][j][m][n] * YnutauCstauc[k][m] * conYnutauCstauc[k][n]
	* W_SSFF (m2_stau[m], m2_stau[n], 0.0L, m2_Char[k], Q2);
    }}

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.12 ===== */

  /* These TSIL evaluations are not reused anywhere... */
  for (k=0; k<4; k++) {

    TSIL_SetParameters (&bar, m2_top, m2_top, m2_top, m2_top, m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);

    V_FFFFS (&gaak, &gaak, &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);
    M_FFFFS (&gaak, &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS);

    topphi0result += Ncolors * 4.0L * (
        SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) * Yttcphi[k] * conYttcphi[k] * m2_top * mFffFS
      + 2.0L * SUMO_CREAL(Yttcphi[i] * Yttcphi[j]) * Yttcphi[k] * conYttcphi[k] * m2_top * mFFffS/2.0L
      + SUMO_CREAL(Yttcphi[i] * Yttcphi[j] * conYttcphi[k] * conYttcphi[k]) * mFFFFS/2.0L
      + 2.0L * SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) * SUMO_CREAL(Yttcphi[k] * Yttcphi[k]) *
	m2_top * mFfFfS/2.0L
      + SUMO_CREAL(Yttcphi[i] * Yttcphi[j] * Yttcphi[k] * Yttcphi[k]) * m2_top * m2_top * mffffS/2.0L);

    topphi0result += Ncolors * 4.0L * (
        SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) * Yttcphi[k] * conYttcphi[k] * (vFFFFS + m2_top * vFffFS)
      + 2.0L * SUMO_CREAL(Yttcphi[i] * Yttcphi[j]) * Yttcphi[k] * conYttcphi[k] * m2_top * vfFfFS 
      + SUMO_CREAL(Yttcphi[i] * Yttcphi[j] * conYttcphi[k] * conYttcphi[k]) * m2_top * vfFFfS 
      + 2.0L * SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) * SUMO_CREAL(Yttcphi[k] * Yttcphi[k]) *
	m2_top * vFFffS
      + SUMO_CREAL(Yttcphi[i] * Yttcphi[j] * Yttcphi[k] * Yttcphi[k]) * m2_top * m2_top * vffffS );

    TSIL_SetParameters (&bar, m2_bot, m2_bot, m2_bot, m2_bot, m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);

    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);
    M_FFFFS (&gaak, &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS);

    botphi0result += Ncolors * 4.0L * (
        SUMO_CREAL(Ybbcphi[i] * conYbbcphi[j]) * Ybbcphi[k] * conYbbcphi[k] * m2_bot * mFffFS
      + 2.0L * SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j]) * Ybbcphi[k] * conYbbcphi[k] * m2_bot * mFFffS/2.0L
      + SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j] * conYbbcphi[k] * conYbbcphi[k]) * mFFFFS/2.0L
      + 2.0L * SUMO_CREAL(Ybbcphi[i] * conYbbcphi[j]) * SUMO_CREAL(Ybbcphi[k] * Ybbcphi[k]) *
	m2_bot * mFfFfS/2.0L
      + SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j] * Ybbcphi[k] * Ybbcphi[k]) * m2_bot * m2_bot * mffffS/2.0L);

    botphi0result += Ncolors * 4.0L * (
        SUMO_CREAL(Ybbcphi[i] * conYbbcphi[j]) * Ybbcphi[k] * conYbbcphi[k] * (vFFFFS + m2_bot * vFffFS)
      + 2.0L * SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j]) * Ybbcphi[k] * conYbbcphi[k] * m2_bot * vfFfFS 
      + SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j] * conYbbcphi[k] * conYbbcphi[k]) * m2_bot * vfFFfS
      + 2.0L * SUMO_CREAL(Ybbcphi[i] * conYbbcphi[j]) * SUMO_CREAL(Ybbcphi[k] * Ybbcphi[k]) *
	m2_bot * vFFffS 
      + SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j] * Ybbcphi[k] * Ybbcphi[k]) * m2_bot * m2_bot * vffffS);

  }

  for (k=0; k<4; k++) {

    TSIL_SetParameters (&bar, m2_tau, m2_tau, m2_tau, m2_tau, m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);

    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);
    M_FFFFS (&gaak, &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS);

    tauphi0result += 4.0L * (SUMO_CREAL(Ytautaucphi[i] * conYtautaucphi[j])
			     * Ytautaucphi[k] * conYtautaucphi[k] * m2_tau * mFffFS
      + 2.0L * SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j]) * Ytautaucphi[k] * conYtautaucphi[k] *
			     m2_tau * mFFffS/2.0L
      + SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j] * conYtautaucphi[k] * conYtautaucphi[k]) *
			     mFFFFS/2.0L
      + 2.0L * SUMO_CREAL(Ytautaucphi[i] * conYtautaucphi[j]) *
			     SUMO_CREAL(Ytautaucphi[k] * Ytautaucphi[k]) * m2_tau * mFfFfS/2.0L
      + SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j] * Ytautaucphi[k] * Ytautaucphi[k]) *
			     m2_tau * m2_tau * mffffS/2.0L);

    tauphi0result += 4.0L * (SUMO_CREAL(Ytautaucphi[i] * conYtautaucphi[j]) * Ytautaucphi[k] *
			     conYtautaucphi[k] * (vFFFFS + m2_tau * vFffFS)
      + 2.0L * SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j]) * Ytautaucphi[k] * conYtautaucphi[k] *
			     m2_tau * vfFfFS
      + SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j] * conYtautaucphi[k] * conYtautaucphi[k]) *
			     m2_tau * vfFFfS 
      + 2.0L * SUMO_CREAL(Ytautaucphi[i] * conYtautaucphi[j]) *
			     SUMO_CREAL(Ytautaucphi[k] * Ytautaucphi[k]) * m2_tau * vFFffS 
      + SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j] * Ytautaucphi[k] * Ytautaucphi[k]) *
			     m2_tau * m2_tau * vffffS);
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.14 ===== */

  for (k=0; k<2; k++) {

    TSIL_SetParameters (&bar, m2_top, m2_bot, m2_top, m2_bot, m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &phipData[k]);

    /* These are the tbtb terms... */
    M_FFFFS (&phipData[k],
	     &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS);

    topbotphipresult += Ncolors * 2.0L * (
      (Ytcbphip[k] * Ytcbphip[k] + Ybctphim[k] * Ybctphim[k]) * m_top * m_bot * 
      SUMO_CREAL(Yttcphi[i] * Ybbcphi[j] + Yttcphi[j] * Ybbcphi[i]) * mFFffS
      + Ytcbphip[k] * Ybctphim[k] * (  SUMO_CREAL(Yttcphi[i] * Ybbcphi[j] + Yttcphi[j] * Ybbcphi[i]) *
				       (mFFFFS + m2_top * m2_bot * mffffS)
      + SUMO_CREAL(Yttcphi[i] * conYbbcphi[j] + Yttcphi[j] * conYbbcphi[i]) * m2_bot * mFfFfS));

    /* Now swap to get btbt terms... */
    TSIL_PermuteResult (&phipData[k], XYandZU, &phipData2[k]);

    M_FFFFS (&phipData2[k],
	     &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS);

    topbotphipresult += Ncolors * 2.0L * (
      (Ytcbphip[k] * Ytcbphip[k] + Ybctphim[k] * Ybctphim[k]) * m_top * m_bot *
      SUMO_CREAL(Yttcphi[i] * conYbbcphi[j] + Yttcphi[j] * conYbbcphi[i]) * mFffFS
      + Ytcbphip[k] * Ybctphim[k] *
      SUMO_CREAL(Yttcphi[i] * conYbbcphi[j] + Yttcphi[j] * conYbbcphi[i]) * m2_top * mFfFfS);
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.13 ===== */
  for (k=0; k<2; k++) {

    /* This is for testing: */
    /* TSIL_SetParameters (&bar, m2_top, m2_top, m2_top, m2_bot, m2_phip_temp[k], Q2); */
    /* TSIL_Evaluate (&bar, s); */
    /* TSIL_CopyResult (&bar, &phipData[k]); */

    /* ORIGINAL CALCULATIONS JUST STARTS HERE */
    /* The TSIL results computed in 4.14 can be used as is for the
       terms shown... appears correct (by comparison to calculation
       w/no tricks) */
    V_FFFFS (&phipData[k], &phipData[k],
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    fac1 = Ytcbphip[k] * Ytcbphip[k] + Ybctphim[k] * Ybctphim[k];

    fac2 = 2.0L * Ytcbphip[k] * Ybctphim[k] * m_top * m_bot;

    topbotphipresult += Ncolors * 2.0L *
      (fac1 * (SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) * (vFFFFS + m2_top * vFffFS)
	 + 2.0L * SUMO_CREAL(Yttcphi[i] * Yttcphi[j]) * m2_top * vfFfFS)
      + fac2 * (SUMO_CREAL(Yttcphi[i] * Yttcphi[j]) * (vfFFfS + m2_top * vffffS)
	 + 2.0L * SUMO_CREAL(Yttcphi[i] * conYttcphi[j]) * vFFffS));

    /* THIS SWAP IS PART OF THE ORIGINAL CALC */
    /* For the t <-> b terms we must swap args... */
/*     TSIL_PermuteResult (&phipData[k], XYandZU, &gaak); */

    /* Test of alternate calc */
    TSIL_SetParameters (&bar, m2_bot, m2_bot, m2_bot, m2_top, m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);

    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    topbotphipresult += Ncolors * 2.0L *
      (fac1 * (SUMO_CREAL(Ybbcphi[i] * conYbbcphi[j]) * (vFFFFS + m2_bot * vFffFS)
	 + 2.0L * SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j]) * m2_bot * vfFfFS)
         + fac2 * (SUMO_CREAL(Ybbcphi[i] * Ybbcphi[j]) * (vfFFfS + m2_bot * vffffS)
	 + 2.0L * SUMO_CREAL(Ybbcphi[i] * conYbbcphi[j]) * vFFffS));
  }

  /* Tau terms; STU eval suffices for these. CHECKED: correct 6/4/12 */
  for (k=0; k<2; k++) {

    TSIL_SetParametersSTU (&bar, m2_tau, m2_tau, 0.0L, m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar,s);
    TSIL_CopyResult (&bar, &gaak);

    V_FFFFS (&gaak, &gaak,
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    nutautauphipresult += 2.0L * Ytaucnutauphim[k] * Ytaucnutauphim[k]
      * (SUMO_CREAL(Ytautaucphi[i] * conYtautaucphi[j]) * (vFFFFS + m2_tau * vFffFS)
	 + 2.0L*SUMO_CREAL(Ytautaucphi[i] * Ytautaucphi[j]) * m2_tau * vfFfFS);
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.15 ===== */

  /* stop terms */
  /* Perform a minimal set of TSIL evaluations for each k: */
  for (k=0; k<4; k++) {

    /* (m,n,p,q) = (0,0,0,0) */
    m = n = p = q = 0;
    TSIL_SetParameters (&bar, m2_stop[m], m2_stop[n], m2_stop[p], m2_stop[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);

    /* (0,0,0,1) and three permutations */
    q = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_stop[n], m2_stop[p], m2_stop[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XUandYZ, &eq15Data[q][p][n][m]);

    /* (0,0,1,1) and one permutation */
    p = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_stop[n], m2_stop[p], m2_stop[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);

    /* (0,1,0,1) and one permutation */
    p = 0; n = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_stop[n], m2_stop[p], m2_stop[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);

    /* (1,0,0,1) and one permutation */
    n = 0; m = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_stop[n], m2_stop[p], m2_stop[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);

    /* (1,1,1,0) and three permutations */
    n = p = 1; q = 0;
    TSIL_SetParameters (&bar, m2_stop[m], m2_stop[n], m2_stop[p], m2_stop[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XUandYZ, &eq15Data[q][p][n][m]);

    /* (1,1,1,1) */
    q = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_stop[n], m2_stop[p], m2_stop[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);

    /* TEST -- braindead TSIL evaluation */
/*     for (m=0; m<2; m++) */
/*       for (n=0; n<2; n++) */
/* 	for (p=0; p<2; p++) */
/* 	  for (q=0; q<2; q++) { */

/* 	    TSIL_SetParameters (&bar, m2_stop[m], m2_stop[n], m2_stop[p], m2_stop[q], m2_phi0_temp[k], Q2); */
/* 	    TSIL_Evaluate (&bar, s); */
/* 	    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]); */

/* 	  } */
    /* END OF TEST -- delete to here */

    /* Now compute contributions for this k */
    for (m=0; m<2; m++) {
    for (n=0; n<2; n++) {

      S_SSS (&eq15Data[0][m][n][0], &sSSS);

      stopphi0result += Ncolors * (lambda00stopstopc[i][k][m][n] * lambda00stopstopc[j][k][n][m] * sSSS
        + 0.5L * lambda00stopstopc[i][j][m][n] * lambda00stopstopc[k][k][n][m]
	* X_SSS (m2_stop[m], m2_stop[n], m2_phi0_temp[k], Q2));

      for (p=0; p<2; p++) {

	U_SSSS (&eq15Data[m][0][n][p], &uSSSS);

        stopphi0result += Ncolors *
	  (2.0L * SUMO_CREAL((lambda0stopstopc[i][m][n] * lambda00stopstopc[j][k][p][m]
 			   + lambda0stopstopc[j][m][n] * lambda00stopstopc[i][k][p][m])
			   * lambda0stopstopc[k][n][p]) * uSSSS
	  + lambda00stopstopc[i][j][m][n] * lambda0stopstopc[k][n][p] * lambda0stopstopc[k][p][m]
	  * W_SSSS (m2_stop[m], m2_stop[n], m2_stop[p], m2_phi0_temp[k], Q2)
	  + SUMO_CREAL(lambda0stopstopc[i][m][n] * lambda0stopstopc[j][p][m]
		     * lambda00stopstopc[k][k][n][p]) *
	   Y_SSSS (m2_stop[m], m2_stop[n], m2_stop[p], m2_phi0_temp[k], s, Q2));

	for (q=0; q<2; q++) {

	  M_SSSSS (&eq15Data[m][n][p][q], &mSSSSS);
	  V_SSSSS (&eq15Data[m][n][p][q], &eq15Data[m][p][n][q], &vSSSSS);

          stopphi0result += Ncolors * (lambda0stopstopc[i][m][p] * lambda0stopstopc[j][q][n]
	        * lambda0stopstopc[k][p][q] * lambda0stopstopc[k][n][m] * mSSSSS
             + 2.0L * SUMO_CREAL(lambda0stopstopc[i][m][n] * lambda0stopstopc[j][p][m]
			       * lambda0stopstopc[k][n][q] * lambda0stopstopc[k][q][p]) * vSSSSS);

	}}}}
  }


  /* sbot terms */
  for (k=0; k<4; k++) {

    m2_phi0_temp[k] = m2_phi0[k];
    if (m2_phi0_temp[k] < 0) m2_phi0_temp[k] = 0;

    /* (m,n,p,q) = (0,0,0,0) */
    m = n = p = q = 0;
    TSIL_SetParameters (&bar, m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_sbot[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);

    /* (0,0,0,1) and three permutations */
    q = 1;
    TSIL_SetParameters (&bar, m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_sbot[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XUandYZ, &eq15Data[q][p][n][m]);

    /* (0,0,1,1) and one permutation */
    p = 1;
    TSIL_SetParameters (&bar, m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_sbot[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);

    /* (0,1,0,1) and one permutation */
    p = 0; n = 1;
    TSIL_SetParameters (&bar, m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_sbot[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);

    /* (1,0,0,1) and one permutation */
    n = 0; m = 1;
    TSIL_SetParameters (&bar, m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_sbot[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);

    /* (1,1,1,0) and three permutations */
    n = p = 1; q = 0;
    TSIL_SetParameters (&bar, m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_sbot[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XUandYZ, &eq15Data[q][p][n][m]);

    /* (1,1,1,1) */
    q = 1;
    TSIL_SetParameters (&bar, m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_sbot[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);

    /* Now compute contributions for this k */
    for (m=0; m<2; m++) {
    for (n=0; n<2; n++) {

      S_SSS (&eq15Data[0][m][n][0], &sSSS);

      sbotphi0result += Ncolors * (lambda00sbotsbotc[i][k][m][n]*lambda00sbotsbotc[j][k][n][m]*sSSS
	+ 0.5L * lambda00sbotsbotc[i][j][m][n] * lambda00sbotsbotc[k][k][n][m]
	* X_SSS (m2_sbot[m], m2_sbot[n], m2_phi0_temp[k], Q2));

      for (p=0; p<2; p++) {

	U_SSSS (&eq15Data[m][0][n][p], &uSSSS);

        sbotphi0result += Ncolors *
	  (2.0L * SUMO_CREAL((lambda0sbotsbotc[i][m][n]*lambda00sbotsbotc[j][k][p][m] +
			     lambda0sbotsbotc[j][m][n]*lambda00sbotsbotc[i][k][p][m])
			   * lambda0sbotsbotc[k][n][p]) * uSSSS
	  + lambda00sbotsbotc[i][j][m][n] * lambda0sbotsbotc[k][n][p] * lambda0sbotsbotc[k][p][m]
	  * W_SSSS (m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_phi0_temp[k], Q2)
	  + SUMO_CREAL(lambda0sbotsbotc[i][m][n] * lambda0sbotsbotc[j][p][m]
		     * lambda00sbotsbotc[k][k][n][p])
	  * Y_SSSS (m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_phi0_temp[k], s, Q2));
	
	for (q=0; q<2; q++) {
	    
	  M_SSSSS (&eq15Data[m][n][p][q], &mSSSSS);
	  V_SSSSS (&eq15Data[m][n][p][q], &eq15Data[m][p][n][q], &vSSSSS);
	  
          sbotphi0result += Ncolors * (lambda0sbotsbotc[i][m][p] * lambda0sbotsbotc[j][q][n]
	        * lambda0sbotsbotc[k][p][q] * lambda0sbotsbotc[k][n][m] * mSSSSS
	     + 2.0L * SUMO_CREAL(lambda0sbotsbotc[i][m][n] * lambda0sbotsbotc[j][p][m]
	   	    	       * lambda0sbotsbotc[k][n][q] * lambda0sbotsbotc[k][q][p]) * vSSSSS);
	}}}}
  }

  /* stau terms */
  for (k=0; k<4; k++) {

    /* (m,n,p,q) = (0,0,0,0) */
    m = n = p = q = 0;
    TSIL_SetParameters (&bar, m2_stau[m], m2_stau[n], m2_stau[p], m2_stau[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);

    /* (0,0,0,1) and three permutations */
    q = 1;
    TSIL_SetParameters (&bar, m2_stau[m], m2_stau[n], m2_stau[p], m2_stau[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XUandYZ, &eq15Data[q][p][n][m]);

    /* (0,0,1,1) and one permutation */
    p = 1;
    TSIL_SetParameters (&bar, m2_stau[m], m2_stau[n], m2_stau[p], m2_stau[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);

    /* (0,1,0,1) and one permutation */
    p = 0; n = 1;
    TSIL_SetParameters (&bar, m2_stau[m], m2_stau[n], m2_stau[p], m2_stau[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);

    /* (1,0,0,1) and one permutation */
    n = 0; m = 1;
    TSIL_SetParameters (&bar, m2_stau[m], m2_stau[n], m2_stau[p], m2_stau[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);

    /* (1,1,1,0) and three permutations */
    n = p = 1; q = 0;
    TSIL_SetParameters (&bar, m2_stau[m], m2_stau[n], m2_stau[p], m2_stau[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq15Data[n][m][q][p]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XUandYZ, &eq15Data[q][p][n][m]);

    /* (1,1,1,1) */
    q = 1;
    TSIL_SetParameters (&bar, m2_stau[m], m2_stau[n], m2_stau[p], m2_stau[q], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);

    /* Now compute contributions for this k */
    for (m=0; m<2; m++) {
    for (n=0; n<2; n++) {

      S_SSS (&eq15Data[0][m][n][0], &sSSS);

      stauphi0Sresult += lambda00staustauc[i][k][m][n] * lambda00staustauc[j][k][n][m] * sSSS;

      stauphi0Xresult += 0.5L * lambda00staustauc[i][j][m][n] * lambda00staustauc[k][k][n][m] * 
                         X_SSS (m2_stau[m], m2_stau[n], m2_phi0_temp[k], Q2);

      for (p=0; p<2; p++) {

	U_SSSS (&eq15Data[m][0][n][p], &uSSSS);

        stauphi0Uresult += 2.0L * SUMO_CREAL((lambda0staustauc[i][m][n] * lambda00staustauc[j][k][p][m] 
			   + lambda0staustauc[j][m][n] * lambda00staustauc[i][k][p][m]) *
					     lambda0staustauc[k][n][p]) * uSSSS;

        stauphi0Wresult += lambda00staustauc[i][j][m][n] * lambda0staustauc[k][n][p] *
	  lambda0staustauc[k][p][m] * W_SSSS (m2_stau[m], m2_stau[n], m2_stau[p], m2_phi0_temp[k], Q2);

        stauphi0Yresult += SUMO_CREAL(lambda0staustauc[i][m][n] * lambda0staustauc[j][p][m] *
				      lambda00staustauc[k][k][n][p]) * 
	  Y_SSSS (m2_stau[m], m2_stau[n], m2_stau[p], m2_phi0_temp[k], s, Q2);
	
	for (q=0; q<2; q++) {
	    
	  M_SSSSS (&eq15Data[m][n][p][q], &mSSSSS);
	  V_SSSSS (&eq15Data[m][n][p][q], &eq15Data[m][p][n][q], &vSSSSS);

          stauphi0Mresult += lambda0staustauc[i][m][p] * lambda0staustauc[j][q][n] *
	    lambda0staustauc[k][p][q] * lambda0staustauc[k][n][m] * mSSSSS;

          stauphi0Vresult += 2.0L * SUMO_CREAL(lambda0staustauc[i][m][n] * lambda0staustauc[j][p][m]
			     * lambda0staustauc[k][n][q] * lambda0staustauc[k][q][p]) * vSSSSS;
	}}}}
  }

  stauphi0result = stauphi0Sresult + stauphi0Xresult + stauphi0Uresult + stauphi0Wresult +
                   stauphi0Yresult + stauphi0Mresult + stauphi0Vresult;

  /* tau sneutrino terms -- DEBUGGED BY HAND - DGR */
  for (k=0; k<4; k++) {

    TSIL_SetParameters (&bar, m2_snu[2], m2_snu[2], m2_snu[2], m2_snu[2], m2_phi0_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak);

    S_SSS (&gaak, &sSSS);

    snutauphi0result += lambda00snusnuc[i][k] * lambda00snusnuc[j][k] * sSSS
      + 0.5L * lambda00snusnuc[i][j] * lambda00snusnuc[k][k]
      * X_SSS (m2_snu[2], m2_snu[2], m2_phi0_temp[k],Q2);

    U_SSSS (&gaak, &uSSSS);

    snutauphi0result += 2.0L * SUMO_CREAL((lambda0snusnuc[i] * lambda00snusnuc[j][k] +
			       lambda0snusnuc[j] * lambda00snusnuc[i][k])
			      * lambda0snusnuc[k]) * uSSSS
      + lambda00snusnuc[i][j] * lambda0snusnuc[k] * lambda0snusnuc[k]
      * W_SSSS (m2_snu[2], m2_snu[2], m2_snu[2], m2_phi0_temp[k], Q2)
      + SUMO_CREAL(lambda0snusnuc[i] * lambda0snusnuc[j] * lambda00snusnuc[k][k])
      * Y_SSSS (m2_snu[2], m2_snu[2], m2_snu[2], m2_phi0_temp[k], s, Q2);

    M_SSSSS (&gaak, &mSSSSS);
    V_SSSSS (&gaak, &gaak, &vSSSSS);

    snutauphi0result += lambda0snusnuc[i] * lambda0snusnuc[j] * lambda0snusnuc[k] *
      lambda0snusnuc[k] * mSSSSS
      + 2.0L * SUMO_CREAL(lambda0snusnuc[i] * lambda0snusnuc[j] *
			  lambda0snusnuc[k] * lambda0snusnuc[k]) * vSSSSS;
  }

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.16 ===== */

  /* A minimal set of TSIL evaluations for each k. (Note that in
     generating the initial results only the XZandYU permutation is
     allowed.)  We reuse eq15Data[][][][] to hold the initial results,
     as needed by M_SSSSS, and use eq16Data[][][][] to hold the
     results with an XYandZU permutation (to exchange t <-> b). */

  for (k=0; k<2; k++) {

    /* (m,n,p,q) = (0,0,0,0) */
    m = n = p = q = 0;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);

    /* (0,0,0,1) and one permutation */
    q = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);

    /* (0,0,1,0) and one permutation */
    p = 1; q = 0;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);

    /* (0,0,1,1) and one permutation */
    q = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);

    /* (0,1,0,1) */
    n = 1; p = 0;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);

    /* (1,0,1,0) */
    m = p = 1; n = q = 0;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);

    /* (1,0,0,1) and one permutation */
    m = 1; n = 0; p = 0; q = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);

    /* (1,1,1,0) and one permutation */
    n = p = 1; q = 0;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);

    /* (1,1,0,1) and one permutation */
    p = 0; q = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);
    TSIL_PermuteResult (&eq15Data[m][n][p][q], XZandYU, &eq15Data[p][q][m][n]);

    /* (1,1,1,1) */
    p = 1;
    TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]);


    /* TEST -- braindead TSIL evaluation */
    /* APPEARS TO MAKE NO DIFFERENCE */
/*     for (m=0; m<2; m++) */
/*     for (n=0; n<2; n++) */
/*     for (p=0; p<2; p++) */
/*     for (q=0; q<2; q++) { */

/*       TSIL_SetParameters (&bar, m2_stop[m], m2_sbot[n], m2_stop[p], m2_sbot[q], m2_phi0_temp[k], Q2); */
/*       TSIL_Evaluate (&bar, s); */
/*       TSIL_CopyResult (&bar, &eq15Data[m][n][p][q]); */

/*     } */
    /* END OF TEST -- delete to here */


    /* Generate t <-> b swapped copies of everything: */
    for (m=0; m<2; m++)
    for (n=0; n<2; n++)
    for (p=0; p<2; p++)
    for (q=0; q<2; q++)
      TSIL_PermuteResult (&eq15Data[m][n][p][q], XYandZU, &eq16Data[n][m][q][p]);

    /* Now compute contributions for this k */
    for (m=0; m<2; m++) {
    for (n=0; n<2; n++) {

      S_SSS (&eq15Data[0][m][n][0], &sSSS);

      stopsbotphipSresult += Ncolors * 2.0L *
	SUMO_CREAL(lambda0psbotstopc[i][k][m][n] * SUMO_CONJ(lambda0psbotstopc[j][k][m][n])) * sSSS;
      stopphipXresult += Ncolors * lambda00stopstopc[i][j][m][n] * lambdapmstopstopc[k][k][n][m] *
	X_SSS (m2_stop[m], m2_stop[n], m2_phip_temp[k], Q2);
      sbotphipXresult += Ncolors * lambda00sbotsbotc[i][j][m][n] * lambdapmsbotsbotc[k][k][n][m] *
	X_SSS (m2_sbot[m], m2_sbot[n], m2_phip_temp[k], Q2);

      for (p=0; p<2; p++) {

	U_SSSS (&eq15Data[m][0][n][p], &uSSSS);

        stopsbotphipUresult += Ncolors *
	  (2.0L * SUMO_CREAL((lambda0stopstopc[i][m][n] * lambda0psbotstopc[j][k][p][m] +
			lambda0stopstopc[j][m][n] * lambda0psbotstopc[i][k][p][m]) *
			     SUMO_CONJ(lambdapsbotstopc[k][p][n])) * uSSSS);

	U_SSSS (&eq16Data[m][0][n][p], &uSSSS);

        sbotstopphipUresult += Ncolors *
	  (2.0L * SUMO_CREAL((lambda0sbotsbotc[i][m][n] * SUMO_CONJ(lambda0psbotstopc[j][k][m][p]) +
			lambda0sbotsbotc[j][m][n] * SUMO_CONJ(lambda0psbotstopc[i][k][m][p])) *
			     lambdapsbotstopc[k][n][p]) * uSSSS);

        stopsbotphipWresult += Ncolors *
	  lambda00stopstopc[i][j][m][n] * lambdapsbotstopc[k][p][m] *
	  SUMO_CONJ(lambdapsbotstopc[k][p][n]) * 
	  W_SSSS (m2_stop[m], m2_stop[n], m2_sbot[p], m2_phip_temp[k], Q2);

        sbotstopphipWresult += Ncolors * lambda00sbotsbotc[i][j][m][n] * lambdapsbotstopc[k][n][p] *
	  SUMO_CONJ(lambdapsbotstopc[k][m][p]) * 
	  W_SSSS (m2_sbot[m], m2_sbot[n], m2_stop[p], m2_phip_temp[k], Q2);

        stopphipYresult += Ncolors * 2.0L * SUMO_CREAL(lambda0stopstopc[i][m][n] *
		      lambda0stopstopc[j][p][m] * lambdapmstopstopc[k][k][n][p]) * 
                      Y_SSSS (m2_stop[m], m2_stop[n], m2_stop[p], m2_phip_temp[k], s, Q2);

        sbotphipYresult += Ncolors * 2.0L * SUMO_CREAL(lambda0sbotsbotc[i][m][n] *
		      lambda0sbotsbotc[j][p][m] * lambdapmsbotsbotc[k][k][n][p]) * 
                      Y_SSSS (m2_sbot[m], m2_sbot[n], m2_sbot[p], m2_phip_temp[k], s, Q2);

	for (q=0; q<2; q++) {

	  M_SSSSS (&eq15Data[m][n][p][q], &mSSSSS);
	  V_SSSSS (&eq15Data[m][n][p][q], &eq15Data[m][p][n][q], &vSSSSS);
	  
          stopsbotphipMresult += Ncolors * (lambda0stopstopc[i][m][p] * lambda0sbotsbotc[j][q][n] +
	       lambda0stopstopc[j][m][p] * lambda0sbotsbotc[i][q][n]) *
	       lambdapsbotstopc[k][n][m] * SUMO_CONJ(lambdapsbotstopc[k][q][p]) * mSSSSS;

          stopsbotphipVresult += Ncolors * 2.0L * SUMO_CREAL(lambda0stopstopc[i][m][n] *
	       lambda0stopstopc[j][p][m] * SUMO_CONJ(lambdapsbotstopc[k][q][n]) *
				 lambdapsbotstopc[k][q][p]) * vSSSSS;

	  V_SSSSS (&eq16Data[m][n][p][q], &eq16Data[m][p][n][q], &vSSSSS);
	  
          sbotstopphipVresult += Ncolors * 2.0L *
	    SUMO_CREAL(lambda0sbotsbotc[i][m][n] * lambda0sbotsbotc[j][p][m] *
                      SUMO_CONJ(lambdapsbotstopc[k][p][q]) * lambdapsbotstopc[k][n][q]) * vSSSSS;
	}}}}
  }

  stopsbotphipresult = stopsbotphipSresult + stopphipXresult + sbotphipXresult + 
                       stopphipYresult + sbotphipYresult + 
                       stopsbotphipUresult + sbotstopphipUresult + 
                       stopsbotphipVresult + sbotstopphipVresult + 
                       stopsbotphipWresult +  sbotstopphipWresult + stopsbotphipMresult;

  /* Next are the stop -> snutau and sbot -> stau terms. Checking on these? */
  for (k=0; k<2; k++) {
    
    n = q = 0;
    TSIL_SetParameters (&bar, m2_snu[2], m2_stau[n], m2_snu[2], m2_stau[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &snutauData[n][q]);

    q = 1;
    TSIL_SetParameters (&bar, m2_snu[2], m2_stau[n], m2_snu[2], m2_stau[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &snutauData[n][q]);
    TSIL_PermuteResult (&snutauData[n][q], XZandYU, &snutauData[q][n]);

    n = 1;
    TSIL_SetParameters (&bar, m2_snu[2], m2_stau[n], m2_snu[2], m2_stau[q], m2_phip_temp[k], Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &snutauData[n][q]);

    /* Generate snutau <-> stau swapped copies of everything: */
    for (n=0; n<2; n++)
    for (q=0; q<2; q++)
      TSIL_PermuteResult (&snutauData[n][q], XYandZU, &stauData[n][q]);

    /* Now compute contributions for this k */

      snutauphipXresult += lambda00snusnuc[i][j] * lambdapmsnutausnutauc[k][k] *
	X_SSS (m2_snu[2], m2_snu[2], m2_phip_temp[k], Q2);

      snutauphipYresult += 2.0L *
	SUMO_CREAL(lambda0snusnuc[i] * lambda0snusnuc[j] * lambdapmsnutausnutauc[k][k]) * 
        Y_SSSS (m2_snu[2], m2_snu[2], m2_snu[2], m2_phip_temp[k], s, Q2);

    /* Terms with one more summed index */
      for (p=0; p<2; p++) {

        V_SSSSS (&snutauData[0][p], &snutauData[0][p], &vSSSSS);

        snutaustauphipVresult += 2.0L * SUMO_CREAL(lambda0snusnuc[i] * lambda0snusnuc[j]
                         * SUMO_CONJ(lambdapstausnutauc[k][p]) * lambdapstausnutauc[k][p]) * vSSSSS;

        snutaustauphipWresult += lambda00snusnuc[i][j] * lambdapstausnutauc[k][p] *
	  SUMO_CONJ(lambdapstausnutauc[k][p])
	  * W_SSSS (m2_snu[2], m2_snu[2], m2_stau[p], m2_phip_temp[k], Q2);

        /* Should not depend on 0 vs. 1 in the second index. (Order [p][0] was wrong before?) */
        S_SSS (&stauData[0][p], &sSSS);

        snutaustauphipSresult += 2.0L * 
	  SUMO_CREAL(lambda0pstausnutauc[i][k][p] * SUMO_CONJ(lambda0pstausnutauc[j][k][p])) * sSSS;

        U_SSSS (&snutauData[0][p], &uSSSS);

        snutaustauphipUresult += 2.0L * SUMO_CREAL((lambda0snusnuc[i] * lambda0pstausnutauc[j][k][p] +
				 lambda0snusnuc[j] * lambda0pstausnutauc[i][k][p])
				* SUMO_CONJ(lambdapstausnutauc[k][p])) * uSSSS;

        /* Terms with two summed indices (not including k) */
        for (n=0; n<2; n++) {

	  M_SSSSS (&snutauData[n][p], &mSSSSS);

          snutaustauphipMresult += (lambda0snusnuc[i] * lambda0staustauc[j][p][n] +
		 lambda0snusnuc[j] * lambda0staustauc[i][p][n])
	     * lambdapstausnutauc[k][n] * SUMO_CONJ(lambdapstausnutauc[k][p]) * mSSSSS;

	  U_SSSS (&stauData[p][n], &uSSSS);

          stausnutauphipUresult += 2.0L *
	    SUMO_CREAL((lambda0staustauc[i][p][n] * SUMO_CONJ(lambda0pstausnutauc[j][k][p]) +
			lambda0staustauc[j][p][n] * SUMO_CONJ(lambda0pstausnutauc[i][k][p]))
		       * lambdapstausnutauc[k][n]) * uSSSS;

          stausnutauphipWresult += lambda00staustauc[i][j][p][n] * lambdapstausnutauc[k][n]
	    * SUMO_CONJ(lambdapstausnutauc[k][p])
	    * W_SSSS (m2_stau[p], m2_stau[n], m2_snu[2], m2_phip_temp[k], Q2);

          stauphipXresult += lambda00staustauc[i][j][p][n] * lambdapmstaustauc[k][k][n][p] *
	    X_SSS (m2_stau[p], m2_stau[n], m2_phip_temp[k], Q2);

	  /* Finally, the terms with three summed indices (plus k) */
	  for (m=0; m<2; m++) {
	    V_SSSSS (&stauData[m][n], &stauData[m][p], &vSSSSS);

            stausnutauphipVresult += 2.0L *
	      SUMO_CREAL(lambda0staustauc[i][m][n] * lambda0staustauc[j][p][m]
			 * SUMO_CONJ(lambdapstausnutauc[k][p]) * lambdapstausnutauc[k][n]) * vSSSSS;

            stauphipYresult += 2.0L *
	      SUMO_CREAL(lambda0staustauc[i][m][n] * lambda0staustauc[j][p][m] *
			 lambdapmstaustauc[k][k][n][p])
	      * Y_SSSS (m2_stau[m], m2_stau[n], m2_stau[p], m2_phip_temp[k], s, Q2);

	  }}}}


  snutaustauphipresult = snutaustauphipSresult + snutauphipXresult + stauphipXresult + 
                         snutauphipYresult + stauphipYresult + 
                         snutaustauphipUresult + stausnutauphipUresult + 
                         snutaustauphipVresult + stausnutauphipVresult + 
                         snutaustauphipWresult + stausnutauphipWresult + snutaustauphipMresult;

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  /* ===== EQ. 4.17 ===== */
  /* INCLUDES ONLY 3RD FAMILY SFERMIONS !!! */

  /* === START X TERMS === */
  for (k=0; k<2; k++)
  for (m=0; m<2; m++) {

    Xstopsnutauresult += lambda00stopstopc[i][j][k][m] * 3.0L * lambdastopstopcsnutausnutauc[m][k] *
      X_SSS (m2_stop[k], m2_stop[m], m2_snu[2], Q2);

    for (n=0; n<2; n++) {

      Xstopresult += lambda00stopstopc[i][j][k][m] * 
	(9.0L * lambdastopstopcstopstopc[m][k][n][n] + 3.0L * lambdastopstopcstopstopc[m][n][n][k])
	* X_SSS (m2_stop[k], m2_stop[m], m2_stop[n], Q2);

      Xstopsbotresult += lambda00stopstopc[i][j][k][m] *
	(9.0L * lambdastopstopcsbotsbotc[m][k][n][n] + 3.0L * lambdastopsbotcsbotstopc[m][n][n][k])
	* X_SSS (m2_stop[k], m2_stop[m], m2_sbot[n], Q2);

      Xstopstauresult += lambda00stopstopc[i][j][k][m] * 3.0L * lambdastopstopcstaustauc[m][k][n][n]
	* X_SSS (m2_stop[k], m2_stop[m], m2_stau[n], Q2);
    }

    Xsbotsnutauresult += lambda00sbotsbotc[i][j][k][m] * 3.0L * lambdasbotsbotcsnutausnutauc[m][k] *
      X_SSS (m2_sbot[k], m2_sbot[m], m2_snu[2], Q2);

    for (n=0; n<2; n++) {

      Xsbotstopresult += lambda00sbotsbotc[i][j][k][m] * 
	(9.0L * lambdastopstopcsbotsbotc[n][n][m][k] + 3.0L * lambdasbotstopcstopsbotc[m][n][n][k])
	* X_SSS (m2_sbot[k], m2_sbot[m], m2_stop[n], Q2);

      Xsbotresult += lambda00sbotsbotc[i][j][k][m] *
	(9.0L * lambdasbotsbotcsbotsbotc[m][k][n][n] + 3.0L * lambdasbotsbotcsbotsbotc[m][n][n][k])
	* X_SSS (m2_sbot[k], m2_sbot[m], m2_sbot[n], Q2);

      Xsbotstauresult += lambda00sbotsbotc[i][j][k][m] * 3.0L * lambdasbotsbotcstaustauc[m][k][n][n]
	* X_SSS (m2_sbot[k], m2_sbot[m], m2_stau[n], Q2);
    }

    Xstausnutauresult += lambda00staustauc[i][j][k][m] * (lambdastaustaucsnutausnutauc[m][k] +
							  lambdastausnutaucsnutaustauc[m][k])
      * X_SSS (m2_stau[k], m2_stau[m], m2_snu[2], Q2);
    
    for (n=0; n<2; n++) {

      Xstaustopresult += lambda00staustauc[i][j][k][m] * 3.0L * lambdastopstopcstaustauc[n][n][m][k] * 
	X_SSS (m2_stau[k], m2_stau[m], m2_stop[n], Q2);

      Xstausbotresult += lambda00staustauc[i][j][k][m] * 3.0L * lambdasbotsbotcstaustauc[n][n][m][k] * 
	X_SSS (m2_stau[k], m2_stau[m], m2_sbot[n], Q2);

      Xstauresult += lambda00staustauc[i][j][k][m] * (lambdastaustaucstaustauc[m][k][n][n] +
						      lambdastaustaucstaustauc[m][n][n][k]) * 
	X_SSS (m2_stau[k], m2_stau[m], m2_stau[n], Q2);
    }
  } /* The k,m loop should end here. Before it wrongly included the Xsnutau* terms below. */
  
    /* DGR check this -- is overall factor of two correct (due to identity of the two terms)? */
  Xsnutauresult += lambda00snutausnutauc[i][j] * 2.0L * lambdasnutausnutaucsnutausnutauc *
    X_SSS (m2_snu[2], m2_snu[2], m2_snu[2], Q2);

  for (n=0; n<2; n++) {
    Xsnutaustopresult += lambda00snutausnutauc[i][j] * 3.0L * lambdastopstopcsnutausnutauc[n][n] * 
      X_SSS (m2_snu[2], m2_snu[2], m2_stop[n], Q2);

    Xsnutausbotresult += lambda00snutausnutauc[i][j] * 3.0L * lambdasbotsbotcsnutausnutauc[n][n] * 
      X_SSS (m2_snu[2], m2_snu[2], m2_sbot[n], Q2);

    Xsnutaustauresult += lambda00snutausnutauc[i][j] * (lambdastaustaucsnutausnutauc[n][n] +
							lambdasnutaustaucstausnutauc[n][n]) * 
      X_SSS (m2_snu[2], m2_snu[2], m2_stau[n], Q2);
  }

  Xresult = Xstopresult + Xstopsbotresult + Xstopstauresult + Xstopsnutauresult + 
            Xsbotstopresult + Xsbotresult + Xsbotstauresult + Xsbotsnutauresult + 
            Xstaustopresult + Xstausbotresult + Xstauresult + Xstausnutauresult + 
            Xsnutaustopresult + Xsnutausbotresult + Xsnutaustauresult + Xsnutauresult;

  /* === START Y TERMS === */
  for (k=0; k<2; k++) {
  for (m=0; m<2; m++) {
  for (n=0; n<2; n++) {

    Ystopsnutauresult += 2.0L * SUMO_CREAL(lambda0stopstopc[i][k][m] * lambda0stopstopc[j][n][k] *
					   3.0L * lambdastopstopcsnutausnutauc[m][n]) * 
      Y_SSSS (m2_stop[k], m2_stop[m], m2_stop[n], m2_snu[2], s, Q2);

    Ysbotsnutauresult += 2.0L * SUMO_CREAL(lambda0sbotsbotc[i][k][m] * lambda0sbotsbotc[j][n][k] *
					   3.0L * lambdasbotsbotcsnutausnutauc[m][n]) * 
      Y_SSSS (m2_sbot[k], m2_sbot[m], m2_sbot[n], m2_snu[2], s, Q2);

    Ystausnutauresult += 2.0L * SUMO_CREAL(lambda0staustauc[i][k][m] * lambda0staustauc[j][n][k] * 
                     (lambdastaustaucsnutausnutauc[m][n] + lambdastausnutaucsnutaustauc[m][n])) * 
                     Y_SSSS (m2_stau[k], m2_stau[m], m2_stau[n], m2_snu[2], s, Q2);

    for (p=0; p<2; p++) {

      Ystopresult += 2.0L * SUMO_CREAL(lambda0stopstopc[i][k][m] * lambda0stopstopc[j][n][k] * 
        (9.0L * lambdastopstopcstopstopc[m][n][p][p] + 3.0L * lambdastopstopcstopstopc[m][p][p][n])) * 
	Y_SSSS (m2_stop[k], m2_stop[m], m2_stop[n], m2_stop[p], s, Q2);

      Ystopsbotresult += 2.0L * SUMO_CREAL(lambda0stopstopc[i][k][m]*lambda0stopstopc[j][n][k] * 
        (9.0L * lambdastopstopcsbotsbotc[m][n][p][p] + 3.0L * lambdastopsbotcsbotstopc[m][p][p][n])) * 
	Y_SSSS (m2_stop[k], m2_stop[m], m2_stop[n], m2_sbot[p], s, Q2);

      Ystopstauresult += 2.0L * SUMO_CREAL(lambda0stopstopc[i][k][m] * lambda0stopstopc[j][n][k] *
					   3.0L * lambdastopstopcstaustauc[m][n][p][p]) * 
	Y_SSSS (m2_stop[k], m2_stop[m], m2_stop[n], m2_stau[p], s, Q2);

      Ysbotstopresult += 2.0L * SUMO_CREAL(lambda0sbotsbotc[i][k][m] * lambda0sbotsbotc[j][n][k] * 
        (9.0L * lambdastopstopcsbotsbotc[p][p][m][n] + 3.0L * lambdasbotstopcstopsbotc[m][p][p][n])) * 
	Y_SSSS (m2_sbot[k], m2_sbot[m], m2_sbot[n], m2_stop[p], s, Q2);

      Ysbotresult += 2.0L * SUMO_CREAL(lambda0sbotsbotc[i][k][m] * lambda0sbotsbotc[j][n][k] * 
        (9.0L * lambdasbotsbotcsbotsbotc[m][n][p][p] + 3.0L * lambdasbotsbotcsbotsbotc[m][p][p][n])) * 
	Y_SSSS (m2_sbot[k], m2_sbot[m], m2_sbot[n], m2_sbot[p], s, Q2);

      Ysbotstauresult += 2.0L * SUMO_CREAL(lambda0sbotsbotc[i][k][m] * lambda0sbotsbotc[j][n][k] *
					   3.0L * lambdasbotsbotcstaustauc[m][n][p][p]) * 
	Y_SSSS (m2_sbot[k], m2_sbot[m], m2_sbot[n], m2_stau[p], s, Q2);

      Ystaustopresult += 2.0L * SUMO_CREAL(lambda0staustauc[i][k][m] * lambda0staustauc[j][n][k] *
					   3.0L * lambdastopstopcstaustauc[p][p][m][n]) * 
	Y_SSSS (m2_stau[k], m2_stau[m], m2_stau[n], m2_stop[p], s, Q2);

      Ystausbotresult += 2.0L * SUMO_CREAL(lambda0staustauc[i][k][m] * lambda0staustauc[j][n][k] *
					   3.0L * lambdasbotsbotcstaustauc[p][p][m][n]) * 
	Y_SSSS (m2_stau[k], m2_stau[m], m2_stau[n], m2_sbot[p], s, Q2);

      Ystauresult += 2.0L * SUMO_CREAL(lambda0staustauc[i][k][m] * lambda0staustauc[j][n][k] * 
               (lambdastaustaucstaustauc[m][n][p][p] + lambdastaustaucstaustauc[m][p][p][n])) * 
               Y_SSSS (m2_stau[k], m2_stau[m], m2_stau[n], m2_stau[p], s, Q2);
    }
  }}}

  for (p=0; p<2; p++) {

    Ysnutaustopresult += 2.0L * SUMO_CREAL(lambda0snutausnutauc[i] * lambda0snutausnutauc[j] * 3.0L *
					   lambdastopstopcsnutausnutauc[p][p]) * 
      Y_SSSS (m2_snu[2], m2_snu[2], m2_snu[2], m2_stop[p], s, Q2);

    Ysnutausbotresult += 2.0L * SUMO_CREAL(lambda0snutausnutauc[i] * lambda0snutausnutauc[j] * 3.0L *
					   lambdasbotsbotcsnutausnutauc[p][p]) * 
      Y_SSSS (m2_snu[2], m2_snu[2], m2_snu[2], m2_sbot[p], s, Q2);

    Ysnutaustauresult += 2.0L * SUMO_CREAL(lambda0snutausnutauc[i] * lambda0snutausnutauc[j] * 
                     (lambdastaustaucsnutausnutauc[p][p] + lambdastausnutaucsnutaustauc[p][p])) * 
                     Y_SSSS (m2_snu[2], m2_snu[2], m2_snu[2], m2_stau[p], s, Q2);
  }

  /* DGR - again check factor of 2 for identical terms */
  Ysnutauresult += 2.0L * SUMO_CREAL(lambda0snutausnutauc[i] * lambda0snutausnutauc[j] * 2.0L *
				     lambdasnutausnutaucsnutausnutauc) * 
    Y_SSSS (m2_snu[2], m2_snu[2], m2_snu[2], m2_snu[2], s, Q2);

  Yresult = Ystopresult + Ystopsbotresult + Ystopstauresult + Ystopsnutauresult + 
            Ysbotstopresult + Ysbotresult + Ysbotstauresult + Ysbotsnutauresult + 
            Ystaustopresult + Ystausbotresult + Ystauresult + Ystausnutauresult + 
            Ysnutaustopresult + Ysnutausbotresult + Ysnutaustauresult + Ysnutauresult;

  /* === START Z TERMS === */

  Zsnutauresult += lambda0snutausnutauc[i]*lambda0snutausnutauc[j]  * 2.0L *
    lambdasnutausnutaucsnutausnutauc * Z_SSSS (m2_snu[2], m2_snu[2], m2_snu[2], m2_snu[2], s, Q2);

  for (k=0; k<2; k++) {
  for (m=0; m<2; m++) {
    Zstopsnutauresult += 2.0L * lambda0stopstopc[i][k][m] * lambda0snutausnutauc[j] * 3.0L *
      lambdastopstopcsnutausnutauc[m][k] * 
      Z_SSSS (m2_stop[k], m2_stop[m], m2_snu[2], m2_snu[2], s, Q2);

    Zsbotsnutauresult += 2.0L * lambda0sbotsbotc[i][k][m] * lambda0snutausnutauc[j] * 3.0L *
      lambdasbotsbotcsnutausnutauc[m][k] * 
      Z_SSSS (m2_sbot[k], m2_sbot[m], m2_snu[2], m2_snu[2], s, Q2);

    Zstausnutauresult += 2.0L * lambda0staustauc[i][k][m] * lambda0snutausnutauc[j] *
      (lambdastaustaucsnutausnutauc[m][k] + lambdastausnutaucsnutaustauc[m][k]) * 
      Z_SSSS (m2_stau[k], m2_stau[m], m2_snu[2], m2_snu[2], s, Q2);

    for (n=0; n<2; n++) {
    for (p=0; p<2; p++) {

      Zstopresult += lambda0stopstopc[i][k][m] * lambda0stopstopc[j][n][p] *
	(9.0L * lambdastopstopcstopstopc[m][k][p][n] + 3.0L * lambdastopstopcstopstopc[m][n][p][k]) * 
	Z_SSSS (m2_stop[k], m2_stop[m], m2_stop[n], m2_stop[p], s, Q2); 

      Zstopsbotresult += 2.0L * lambda0stopstopc[i][k][m] * lambda0sbotsbotc[j][n][p] * 
	(9.0L * lambdastopstopcsbotsbotc[m][k][p][n] + 3.0L * lambdastopsbotcsbotstopc[m][n][p][k]) * 
	Z_SSSS (m2_stop[k], m2_stop[m], m2_sbot[n], m2_sbot[p], s, Q2); 

      Zstopstauresult += 2.0L * lambda0stopstopc[i][k][m]*lambda0staustauc[j][n][p] * 3.0L *
	lambdastopstopcstaustauc[m][k][p][n] * 
	Z_SSSS (m2_stop[k], m2_stop[m], m2_stau[n], m2_stau[p], s, Q2);

      Zsbotresult += lambda0sbotsbotc[i][k][m] * lambda0sbotsbotc[j][n][p] *
	(9.0L * lambdasbotsbotcsbotsbotc[m][k][p][n] + 3.0L * lambdasbotsbotcsbotsbotc[m][n][p][k]) * 
	Z_SSSS (m2_sbot[k], m2_sbot[m], m2_sbot[n], m2_sbot[p], s, Q2);

      Zsbotstauresult += 2.0L * lambda0sbotsbotc[i][k][m] * lambda0staustauc[j][n][p] * 3.0L *
	lambdasbotsbotcstaustauc[m][k][p][n] * 
	Z_SSSS (m2_sbot[k], m2_sbot[m], m2_stau[n], m2_stau[p], s, Q2);

      Zstauresult += lambda0staustauc[i][k][m] * lambda0staustauc[j][n][p] *
	(lambdastaustaucstaustauc[m][k][p][n] + lambdastaustaucstaustauc[m][n][p][k]) * 
	Z_SSSS (m2_stau[k], m2_stau[m], m2_stau[n], m2_stau[p], s, Q2);
    }}
  }}

  Zresult = Zstopresult + Zstopsbotresult + Zstopstauresult + Zstopsnutauresult + 
            Zsbotresult + Zsbotstauresult + Zsbotsnutauresult + 
            Zstauresult + Zstausnutauresult + Zsnutauresult;

#ifdef SUMO_HIGGS_TIMING
  t1 = clock ();
  printf("\t%lf secs", difftime(t1,t0)/CLOCKS_PER_SEC);
  t0 = t1;
  printf("\n");
#endif

  XYZresult = Xresult + Yresult + Zresult;

  totalresult = topstopNresult + botsbotNresult + taustauNresult + nutausnutauNresult +
                topsbotCresult + botstopCresult + tausnutauCresult + nutaustauCresult + 
                stopphi0result + sbotphi0result + stauphi0result + snutauphi0result + 
                topphi0result + botphi0result + tauphi0result + 
                stopsbotphipresult + snutaustauphipresult + topbotphipresult + nutautauphipresult +
                XYZresult;
  
  if (1 == VERBOSE) {
  printf("Total nonQCD part of two-loop Higgs-self energy:\n "); printc2(totalresult); printf(";\n");
  printf("  Breakdown:\n");
  printf("    topstopN =      "); printc2(topstopNresult); printf(";\n");
  printf("    botsbotN =      "); printc2(botsbotNresult); printf(";\n");
  printf("    taustauN =      "); printc2(taustauNresult); printf(";\n");
  printf("      Breakdown:\n");
  printf("        taustauN M =      "); printc2(taustauNMresult); printf(";\n");
  printf("        taustauN VFFFFS =      "); printc2(taustauNVFFFFSresult); printf(";\n");
  printf("        taustauN W =      "); printc2(taustauNWresult); printf(";\n");
  printf("        taustauN VSSSFF =      "); printc2(taustauNVSSSFFresult); printf(";\n");
  printf("    nutausnutauN =  "); printc2(nutausnutauNresult); printf(";\n");
 
  printf("    topsbotC =      "); printc2(topsbotCresult); printf(";\n");
  printf("    botstopC =      "); printc2(botstopCresult); printf(";\n");
  printf("    tausnutauC =    "); printc2(tausnutauCresult); printf(";\n");
  printf("    nutaustauC =    "); printc2(nutaustauCresult); printf(";\n");

  printf("    stopphi0 =      "); printc2(stopphi0result); printf(";\n");
  printf("    sbotphi0 =      "); printc2(sbotphi0result); printf(";\n");
  printf("    stauphi0 =      "); printc2(stauphi0result); printf(";\n");
  printf("      Breakdown:\n");
  printf("        stauphi0S =      "); printc2(stauphi0Sresult); printf(";\n");
  printf("        stauphi0X =      "); printc2(stauphi0Xresult); printf(";\n");
  printf("        stauphi0U =      "); printc2(stauphi0Uresult); printf(";\n");
  printf("        stauphi0W =      "); printc2(stauphi0Wresult); printf(";\n");
  printf("        stauphi0Y =      "); printc2(stauphi0Yresult); printf(";\n");
  printf("        stauphi0M =      "); printc2(stauphi0Mresult); printf(";\n");
  printf("        stauphi0V =      "); printc2(stauphi0Vresult); printf(";\n");
  printf("    snutauphi0 =    "); printc2(snutauphi0result); printf(";\n");

  printf("    topphi0 =       "); printc2(topphi0result); printf(";\n");
  printf("    botphi0 =       "); printc2(botphi0result); printf(";\n");
  printf("    tauphi0 =       "); printc2(tauphi0result); printf(";\n");

  printf("    stopsbotphip =  "); printc2(stopsbotphipresult); printf(";\n");
  printf("      Breakdown:\n");
  printf("        stopsbotphipSresult: "); printc2(stopsbotphipSresult); printf("\n");
  printf("        stopphipXresult: "); printc2(stopphipXresult); printf("\n");
  printf("        sbotphipXresult: "); printc2(sbotphipXresult); printf("\n");
  printf("        stopsbotphipUresult: "); printc2(stopsbotphipUresult); printf("\n");
  printf("        sbotstopphipUresult: "); printc2(sbotstopphipUresult); printf("\n");
  printf("        stopsbotphipWresult: "); printc2(stopsbotphipWresult); printf("\n");
  printf("        sbotstopphipWresult: "); printc2(sbotstopphipWresult); printf("\n");
  printf("        stopphipYresult: "); printc2(stopphipYresult); printf("\n");
  printf("        sbotphipYresult: "); printc2(sbotphipYresult); printf("\n");
  printf("        stopsbotphipMresult: "); printc2(stopsbotphipMresult); printf("\n");
  printf("        stopsbotphipVresult: "); printc2(stopsbotphipVresult); printf("\n");
  printf("        sbotstopphipVresult: "); printc2(sbotstopphipVresult); printf("\n");

  printf("    snustauphip =   "); printc2(snutaustauphipresult); printf(";\n");
  printf("      Breakdown:\n");
  printf("        snutaustauphipSresult: "); printc2(snutaustauphipSresult); printf("\n");
  printf("        snutauphipXresult: "); printc2(snutauphipXresult); printf("\n");
  printf("        stauphipXresult: "); printc2(stauphipXresult); printf("\n");
  printf("        snutauphipYresult: "); printc2(snutauphipYresult); printf("\n");
  printf("        stauphipYresult: "); printc2(stauphipYresult); printf("\n");
  printf("        snutaustauphipUresult: "); printc2(snutaustauphipUresult); printf("\n");
  printf("        stausnutauphipUresult: "); printc2(stausnutauphipUresult); printf("\n");
  printf("        snutaustauphipVresult: "); printc2(snutaustauphipVresult); printf("\n");
  printf("        stausnutauphipVresult: "); printc2(stausnutauphipVresult); printf("\n");
  printf("        snutaustauphipWresult: "); printc2(snutaustauphipWresult); printf("\n");
  printf("        stausnutauphipWresult: "); printc2(stausnutauphipWresult); printf("\n");
  printf("        snutaustauphipMresult: "); printc2(snutaustauphipMresult); printf("\n");
  printf("        topbotphip =    "); printc2(topbotphipresult); printf(";\n");
  printf("        nutautauphip =  "); printc2(nutautauphipresult); printf(";\n");

  printf("    sfermion XYZ =  "); printc2(XYZresult); printf(";\n");
  }

  return totalresult;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Two-loop self energy function for neutral Higgs bosons evaluated at
   the specified squared momentum; QCD contribution (hep-ph/0405022,
   eqs. 4.1-4.5). Has pole masses substituted for the tree masses.
*/

SUMO_COMPLEX Pi2QCD_phi0 (int i, int j, SUMO_REAL s)
{
  SUMO_MODEL temp_model;
  SUMO_COMPLEX result;
  
  SUMO_Backup (&temp_model);
  SUMO_Tree_masses_from_pole ();
  result = pi2QCD_phi0 (i, j, s);
  SUMO_Restore (&temp_model);

  return result;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Two-loop self energy function for neutral Higgs particles evaluated
   at the specified squared momentum; non-QCD contributions
   (hep-ph/0405022, eqs. 4.6-4.17). Has pole masses substituted for
   the tree masses.
*/

SUMO_COMPLEX Pi2nonQCD_phi0 (int i, int j, SUMO_REAL s)
{
  SUMO_MODEL temp_model;
  SUMO_COMPLEX result;
  
  SUMO_Backup (&temp_model);
  SUMO_Tree_masses_from_pole ();
  result = pi2nonQCD_phi0 (i, j, s);
  SUMO_Restore (&temp_model);

  return result;
}
