/* 
   Computes pole masses of top quark at one or two loops, using
   eqs. (5.11)-(5.13) of hep-ph/0509115. (Bottom and tau pole masses
   may be done also, but most likely not.)
*/

#include "supermodel.h"
#include "self_fermion.h"

/* Facts about QCD: */
#define CG 3.L
#define Cq (4.L/3.L)
#define Iq 0.5L

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */

TSIL_COMPLEX pi1tilde_top ()
{
  TSIL_COMPLEX temp = 0.0L;
  TSIL_COMPLEX result = 0.0L;
  TSIL_REAL Top, Bot, Gluino, Stop[2], Sbot[2], Char[2], Phip[2];
  TSIL_REAL Z, W, Neut[4], Phi0[4];
  TSIL_REAL mtop, mgluino; 
  int i,j;
  TSIL_COMPLEX BFS, BfS, BFV, BfV;

  Top = m2_top;
  mtop = m_top;
  Bot = m2_bot;
  Gluino = m2_gluino;
  mgluino = m_gluino;
  for (i=0; i<2; i++) {
    Stop[i] = m2_stop[i];
    Sbot[i] = m2_sbot[i];
    Char[i] = m2_Char[i];
    Phip[i] = m2_phip[i];
/* FOLLOWING IS A CHEESY CHEAT!!!! */
    if (Phip[i] < 0.0L) Phip[i] = 10.0;      
  }
  for (i=0; i<4; i++) {
    Neut[i] = m2_Neut[i];
    Phi0[i] = m2_phi0[i];
/* FOLLOWING IS A CHEESY CHEAT!!!! */
    if (Phi0[i] < 0.0L) Phi0[i] = 10.0;      
  }
  Z = m2_Z;
  W = m2_W;

  temp = Top * (10.L - 6.L * TSIL_LOG(Top/(Q2)));

  for (j=0; j<2; j++) {
    bFS (Gluino, Stop[j], Top, Q2, &BFS, &BfS);

    temp += 2.0L * BFS - 4.0L * TSIL_CREAL(Lstop[j] * Rstopc[j]) * 
            mtop * mgluino * BfS;
  }

  result = g3 * g3 * Cq * temp;
   
  for (j=0; j<4; j++) {
    bFS (Top, Phi0[j], Top, Q2, &BFS, &BfS);

    result += 2.0 * (SUMO_AbsSq(Yttcphi[j]) * BFS +
              TSIL_CREAL (Yttcphi[j] * Yttcphi[j]) *
              mtop * mtop * BfS);
  }

  for (j=0; j<2; j++) {
    bFS (Bot, Phip[j], Top, Q2, &BFS, &BfS);

    result += (SUMO_AbsSq(Ytcbphip[j]) +
               SUMO_AbsSq(Ybctphim[j])) * BFS +
               2.0L * TSIL_CREAL (Ytcbphip[j] * Ybctphim[j]) *
               mtop * m_bot * BfS;
  }

  for (i=0; i<4; i++) {
  for (j=0; j<2; j++) {
    bFS (Neut[i], Stop[j], Top, Q2, &BFS, &BfS);

    result += (SUMO_AbsSq(YtNstopc[i][j]) +
               SUMO_AbsSq(YtcNstop[i][j])) * BFS
        + 2.0L * TSIL_CREAL (YtNstopc[i][j] * YtcNstop[i][j]) * 
          mtop * m_Neut[i] * BfS;
  }}

  for (i=0; i<2; i++) {
  for (j=0; j<2; j++) {
    bFS (Char[i], Sbot[j], Top, Q2, &BFS, &BfS);

    result += (SUMO_AbsSq(YtCsbotc[i][j]) +
               SUMO_AbsSq(YtcCsbot[i][j])) * BFS
        + 2.0L * TSIL_CREAL (YtCsbotc[i][j] * YtcCsbot[i][j]) * 
          mtop * m_Char[i] * BfS;
  }}

  result += (4.0L/9.0L) * e2 * Top * (10.L - 6.L * TSIL_LOG(Top/Q2));

  bFV (Top, Z, Top, Q2, 0.0L, &BFV, &BfV);
 
  result += (g2plusgp2/4.0L - (2.0L/3.0L) * gp2 + 
            (8.0L/9.0L) * gp2 * gp2/g2plusgp2) * BFV;

  result += ((8.0L/9.0L) * gp2 * gp2/g2plusgp2 - (2.0L/3.0L) * gp2) *
            mtop * mtop * BfV; 

  bFV (Bot, W, Top, Q2, 0.0L, &BFV, &BfV);

  result += g2 * BFV/2.0L;

  return (SUMO_oneloopfactor * result);
}
