/* 
   Computes tree-level masses and associated mixing matrices for
   particles in the MSSM.
*/

#include "supermodel.h"

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */

int SUMO_Tree_Masses ()
{  
  if (0 == is_updated) SUMO_Update ();

  SUMO_Tree_WZ ();
  SUMO_Tree_Higgs ();
  SUMO_Tree_Neutralinos ();
  SUMO_Tree_Charginos ();
  SUMO_Tree_TopBotTau ();
  SUMO_Tree_Sfermions ();
  SUMO_Tree_Gluino ();
  are_tree_masses_updated = 1;

  return 0;
}

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* hep-ph/0206136 eqs. (2.3), (2.4) */

int SUMO_Tree_WZ ()
{
  m2_W = 0.5L * g2 * v2; 
  m2_Z = 0.5L * g2plusgp2 * v2; 

  return 0;
}

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* hep-ph/0206136 eqs. (2.6) through (2.17) */

int SUMO_Tree_Higgs ()
{
  int i;
  TSIL_REAL m2uu = mu2 + m2Hu;
  TSIL_REAL m2dd = mu2 + m2Hd;

  SUMO_EIGEN_HERM eigenM2GpHp;
  SUMO_EIGEN_HERM eigenM2h0H0;
  SUMO_EIGEN_HERM eigenM2G0A0;

  TSIL_COMPLEX M2h0H0[2][2]; 
  TSIL_COMPLEX M2G0A0[2][2];
  TSIL_COMPLEX M2GpHp[2][2];
 
  M2h0H0[0][0] = m2uu + 0.25L * g2plusgp2 * (3.0L * vu2 - vd2);
  M2h0H0[1][1] = m2dd + 0.25L * g2plusgp2 * (3.0L * vd2 - vu2);
  M2h0H0[0][1] = M2h0H0[1][0] = -b - 0.5L * g2plusgp2 * vuvd;

  M2G0A0[0][0] = m2uu + 0.25L * g2plusgp2 * (vu2 - vd2);
  M2G0A0[1][1] = m2dd + 0.25L * g2plusgp2 * (vd2 - vu2);
  M2G0A0[0][1] = M2G0A0[1][0] = b;

  M2GpHp[0][0] = m2uu + 0.25L * (g2plusgp2 * vu2 + (g2 - gp2) * vd2);
  M2GpHp[1][1] = m2dd + 0.25L * (g2plusgp2 * vd2 + (g2 - gp2) * vu2);
  M2GpHp[0][1] = M2GpHp[1][0] = b + 0.5L * g2 * vuvd;

  SUMO_DiagonalizeHerm (*M2h0H0, 2, &eigenM2h0H0);
  SUMO_DiagonalizeHerm (*M2G0A0, 2, &eigenM2G0A0);
  SUMO_DiagonalizeHerm (*M2GpHp, 2, &eigenM2GpHp);
  
  for (i=0; i<2; i++){
    m2_phi0[i] = eigenM2h0H0.eigenvalues[i];
    ku[i] = eigenM2h0H0.eigenvectors[0][i];
    kd[i] = eigenM2h0H0.eigenvectors[1][i];

    m2_phi0[i+2] = eigenM2G0A0.eigenvalues[i];
    ku[i+2] = I*eigenM2G0A0.eigenvectors[0][i];
    kd[i+2] = I*eigenM2G0A0.eigenvectors[1][i];

    m2_phip[i] = eigenM2GpHp.eigenvalues[i];
    kup[i] = eigenM2GpHp.eigenvectors[0][i];
    kdp[i] = eigenM2GpHp.eigenvectors[1][i];
  }

  calpha = eigenM2h0H0.eigenvectors[0][0];
  salpha = eigenM2h0H0.eigenvectors[0][1];
  sbeta0 = eigenM2G0A0.eigenvectors[0][0];
  cbeta0 = eigenM2G0A0.eigenvectors[0][1];
  sbetapm = kup[0];
  cbetapm = kup[1];

  c2alpha = calpha * calpha - salpha * salpha;
  s2alpha = 2.0L * calpha * salpha;

  c2beta0 = cbeta0 * cbeta0 - sbeta0 * sbeta0;
  s2beta0 = 2.0L * cbeta0 * sbeta0;

  c2betapm = cbetapm * cbetapm - sbetapm * sbetapm;
  s2betapm = 2.0L * cbetapm * sbetapm;

  for (i=0; i<4; i++) {
    kuc[i] = SUMO_CONJ(ku[i]);
    kdc[i] = SUMO_CONJ(kd[i]);
  }
  
  return 0;
}

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* hep-ph/0206136 eqs. (2.19) through (2.23) */

int SUMO_Tree_Neutralinos ()
{
  int i,j,k;
  TSIL_COMPLEX Mneut[4][4]; 
  TSIL_COMPLEX M2neut[4][4]; 
  SUMO_EIGEN_HERM eigenM2neut;
  TSIL_COMPLEX diagcomplex[4];
  TSIL_COMPLEX diagphases[4];
  TSIL_REAL vuosqrt2 = vu/SQRT2;
  TSIL_REAL vdosqrt2 = vd/SQRT2;

  Mneut[0][0] = M1;
  Mneut[1][1] = M2;
  Mneut[2][2] = Mneut[3][3] = Mneut[0][1] = Mneut[1][0] = 0.0L;
  Mneut[2][3] = Mneut[3][2] = -mu;
  Mneut[0][2] = Mneut[2][0] = -gp * vdosqrt2;
  Mneut[0][3] = Mneut[3][0] = gp * vuosqrt2;
  Mneut[1][2] = Mneut[2][1] = g * vdosqrt2;
  Mneut[1][3] = Mneut[3][1] = -g * vuosqrt2;

  for (i=0; i<4; i++){
    for (j=0; j<4; j++){
      M2neut[i][j] = 0.0L;
      for (k=0; k<4; k++){
	M2neut[i][j] += SUMO_CONJ(Mneut[k][i]) * Mneut[k][j];
  }}}

  SUMO_DiagonalizeHerm (*M2neut, 4, &eigenM2neut);

  for (i=0; i<4; i++) {
    m2_Neut[i] = eigenM2neut.eigenvalues[i];
    m_Neut[i] = TSIL_SQRT(m2_Neut[i]);
    diagcomplex[i] = 0.0L;
    for (j=0; j<4; j++) {
      for (k=0; k<4; k++) {
	diagcomplex[i] += Mneut[k][j] * eigenM2neut.eigenvectors[j][i] 
	  * eigenM2neut.eigenvectors[k][i];
      }}

    /* DGR ADDED conditional 7/5/12 to permit setting g=gp=0 */
    if (TSIL_CABS(diagcomplex[i]) > TSIL_TOL)
      diagphases[i] = 
	SUMO_CONJ(TSIL_CSQRT(diagcomplex[i]/TSIL_CABS(diagcomplex[i])));
    else
      diagphases[i] = 1.0L;
  }

  for (i=0; i<4; i++){  
    for (j=0; j<4; j++){
      Nmixc[i][j] = eigenM2neut.eigenvectors[j][i] * diagphases[i];
      Nmix[i][j] = SUMO_CONJ(Nmixc[i][j]);
    }}

  return 0;
}

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* hep-ph/0206136 eqs. (2.24) through (2.26) */

int SUMO_Tree_Charginos ()
{
  int i,j,k;
  TSIL_COMPLEX Mchar[2][2];
  TSIL_COMPLEX M2char[2][2]; 
  SUMO_EIGEN_HERM eigenM2char;

  Mchar[0][0] = M2;
  Mchar[0][1] = g * vu;
  Mchar[1][0] = g * vd;
  Mchar[1][1] = mu;

  for (i=0; i<2; i++){
    for (j=0; j<2; j++){
      M2char[i][j] = 0.0L;
      for (k=0; k<2; k++){
	M2char[i][j] += SUMO_CONJ(Mchar[i][k]) * Mchar[j][k];
      }
    }
  }
  SUMO_DiagonalizeHerm (*M2char, 2, &eigenM2char);
                                     
  for (i=0; i<2; i++){  
    m2_Char[i] = eigenM2char.eigenvalues[i];
    m_Char[i] = TSIL_SQRT(eigenM2char.eigenvalues[i]);
    for (j=0; j<2; j++){
      Umixc[i][j] = eigenM2char.eigenvectors[j][i];
      Umix[i][j] = SUMO_CONJ(Umixc[i][j]);
    }
  }

  for (i=0; i<2; i++){  
    for (j=0; j<2; j++){
      Vmix[i][j] = 0.L;
      for (k=0; k<2; k++){
	Vmix[i][j] += Umixc[i][k] * 
	  Mchar[k][j]/TSIL_CSQRT(m2_Char[i]);
      }
      Vmixc[i][j] = SUMO_CONJ(Vmix[i][j]);
    }
  }

  return 0;
}

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* hep-ph/0206136 eq. (2.28) */
  
int SUMO_Tree_TopBotTau ()
{
  m_top = ytop * vu; 
  m_bot = ybot * vd; 
  m_tau = ytau * vd; 

  m2_top = m_top * m_top; 
  m2_bot = m_bot * m_bot; 
  m2_tau = m_tau * m_tau; 

  return 0;
}

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* hep-ph/0206136 eqs. (2.31) through (2.42) */

int SUMO_Tree_Sfermions ()
{
  int i;
  TSIL_REAL vd2mvu2o2 = 0.5L * (vd2 - vu2);

  TSIL_REAL DsuL = (g2/2.0L - gp2/6.0L) * vd2mvu2o2;
  TSIL_REAL DsdL = (-g2/2.0L - gp2/6.0L) * vd2mvu2o2;
  TSIL_REAL Dsnu = (g2/2.0L + gp2/2.0L) * vd2mvu2o2;
  TSIL_REAL DseL = (-g2/2.0L + gp2/2.0L) * vd2mvu2o2;
  TSIL_REAL DsuR = (2.0L*gp2/3.0L) * vd2mvu2o2;
  TSIL_REAL DsdR = (-gp2/3.0L) * vd2mvu2o2;
  TSIL_REAL DseR = (-gp2) * vd2mvu2o2;

  SUMO_EIGEN_HERM eigenM2stop;
  SUMO_EIGEN_HERM eigenM2sbot;
  SUMO_EIGEN_HERM eigenM2stau;

  TSIL_COMPLEX M2stop[2][2]; 
  TSIL_COMPLEX M2sbot[2][2]; 
  TSIL_COMPLEX M2stau[2][2]; 

  M2stop[0][0] = m2Q[2] + m2_top + DsuL; 
  M2stop[0][1] = vu * atopc - vd * mu * ytop;
  M2stop[1][0] = SUMO_CONJ (M2stop[0][1]);
  M2stop[1][1] = m2u[2] + m2_top + DsuR;

  M2sbot[0][0] = m2Q[2] + m2_bot + DsdL; 
  M2sbot[0][1] = vd * abotc - vu * mu * ybot;
  M2sbot[1][0] = SUMO_CONJ (M2sbot[0][1]);
  M2sbot[1][1] = m2d[2] + m2_bot + DsdR;

  M2stau[0][0] = m2L[2] + m2_tau + DseL; 
  M2stau[0][1] = vd * atauc - vu * mu * ytau;
  M2stau[1][0] = SUMO_CONJ (M2stau[0][1]);
  M2stau[1][1] = m2e[2] + m2_tau + DseR;
  
  SUMO_DiagonalizeHerm (*M2stop, 2, &eigenM2stop);
  SUMO_DiagonalizeHerm (*M2sbot, 2, &eigenM2sbot);
  SUMO_DiagonalizeHerm (*M2stau, 2, &eigenM2stau);

  for (i=0; i<2; i++){
    m2_stop[i] = eigenM2stop.eigenvalues[i];
    Lstop[i] = eigenM2stop.eigenvectors[0][i];
    Rstop[i] = eigenM2stop.eigenvectors[1][i];
    Lstopc[i] = SUMO_CONJ(Lstop[i]);
    Rstopc[i] = SUMO_CONJ(Rstop[i]);

    m2_sbot[i] = eigenM2sbot.eigenvalues[i];
    Lsbot[i] = eigenM2sbot.eigenvectors[0][i];
    Rsbot[i] = eigenM2sbot.eigenvectors[1][i];
    Lsbotc[i] = SUMO_CONJ(Lsbot[i]);
    Rsbotc[i] = SUMO_CONJ(Rsbot[i]);

    m2_stau[i] = eigenM2stau.eigenvalues[i];
    Lstau[i] = eigenM2stau.eigenvectors[0][i];
    Rstau[i] = eigenM2stau.eigenvectors[1][i];
    Lstauc[i] = SUMO_CONJ(Lstau[i]);
    Rstauc[i] = SUMO_CONJ(Rstau[i]);
  }

  m2_snu[2] = m2L[2] + Dsnu;

  for (i=0; i<2; i++){
    m2_suL[i] = m2Q[i] + DsuL;
    m2_sdL[i] = m2Q[i] + DsdL;
    m2_suR[i] = m2u[i] + DsuR;
    m2_sdR[i] = m2d[i] + DsdR;
    m2_snu[i] = m2L[i] + Dsnu;
    m2_seL[i] = m2L[i] + DseL;
    m2_seR[i] = m2e[i] + DseR;
  }

  return 0;
}

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */

int SUMO_Tree_Gluino ()
{
  m2_gluino = TSIL_CABS(M3 * M3); 
  m_gluino = TSIL_CABS(M3); 

  return 0;
}

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
