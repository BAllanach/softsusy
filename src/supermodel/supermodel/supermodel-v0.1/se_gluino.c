/*   
   Gluino self energy functions, given in eqs. (5.1)-(5.3) of
   hep-ph/0509115
*/

#include "supermodel.h"
#include "self_fermion.h"

/* Facts about QCD: */
#define CG 3.L
#define Cq (4.L/3.L)
#define Iq 0.5L

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

SUMO_COMPLEX Pi1tilde_gluino (int expand_around_pole)
{
  int j,n;
  TSIL_COMPLEX pi1Vtilde;
  TSIL_COMPLEX pi1Stilde = 0.L;
  TSIL_COMPLEX BFS, BfS;
  TSIL_REAL Gluino, mgluino, quark[6], mquark[6], squark[6][2];
  TSIL_COMPLEX reLRc[6][2];

  if (0 == expand_around_pole) {
    Gluino = m2_gluino;

    for (n=0; n<4; n++) quark[n] = 0.L;
    quark[4] = m2_bot;
    quark[5] = m2_top;

    mgluino = m_gluino;
    //printf("here gluino %Lg\n", mgluino);
    for (n=0; n<6; n++) {
     mquark[n] = TSIL_SQRT (quark[n]);
     //printf("qua %Lg \n",quark[n]);
    }
    squark[0][0] = m2_sdR[0];
    squark[0][1] = m2_sdL[0];
    squark[1][0] = m2_suR[0];
    squark[1][1] = m2_suL[0];
    squark[2][0] = m2_sdR[1];
    squark[2][1] = m2_sdL[1];
    squark[3][0] = m2_suR[1];
    squark[3][1] = m2_suL[1];
    for (j=0; j<2; j++) {
      squark[4][j] = m2_sbot[j];
      squark[5][j] = m2_stop[j];
    }
     
    
  } else {
    SUMO_Error("Pi1tilde_gluino",
               "Expanding around pole not implemented yet.",3);
  }

  
  for (j=0; j<2; j++) {
    for (n=0; n<4; n++) {
      reLRc[n][j] = 0.0L;
    }
    reLRc[4][j] = TSIL_CREAL(Lsbot[j] * Rsbotc[j]);
    reLRc[5][j] = TSIL_CREAL(Lstop[j] * Rstopc[j]);
  }

  pi1Vtilde = CG * Gluino * (10.L - 6.L * TSIL_LOG(Gluino/Q2));

  for (n=0; n<6; n++){
    for (j=0; j<2; j++){
      bFS (quark[n], squark[n][j], Gluino, Q2, &BFS, &BfS);
      pi1Stilde += BFS - 2.0L * reLRc[n][j] * mgluino * mquark[n] * BfS;
    }}

  return (SUMO_oneloopfactor * g3 * g3 * (pi1Vtilde + 4.L * Iq * pi1Stilde));
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

TSIL_COMPLEX Pi2tilde_gluino (int expand_around_pole)
{
  TSIL_REAL Gluino, gluino, mgluino, quark[6], mquark[6], squark[6][2];
  TSIL_COMPLEX Rsq[6][2], Lsq[6][2], Rsqc[6][2], Lsqc[6][2], reLRc[6][2];
  TSIL_COMPLEX pi1leftover = 0.L;
  TSIL_COMPLEX pi2tildeMSFFSF = 0.L;
  TSIL_COMPLEX pi2tildeVSFFFS = 0.L;
  TSIL_COMPLEX pi2tildeVFSSFF = 0.L;
  TSIL_COMPLEX pi2tildeYFSSS  = 0.L;
  TSIL_COMPLEX pi2tilde1vector = 0.L;
  TSIL_COMPLEX pi2tilde2vector;
  TSIL_COMPLEX BFS[6][2], BfS[6][2], BpFS[6][2], BpfS[6][2];
  TSIL_DATA    int_data[6][2][2];
  TSIL_COMPLEX MSFFSF, MSFFSf, MSFfSF, MSFfSf, MSffSF, MSffSf;
  TSIL_COMPLEX VSFFFS, VSFFfS, VSFfFS, VSFffS, VSffFS, VSfffS;
  TSIL_COMPLEX VFSSFF, VFSSff, VfSSFF, VfSSff;
  TSIL_COMPLEX YFSSS, YfSSS;
  TSIL_COMPLEX Lj, Ljc, Rj, Rjc, Lk, Lkc, Rk, Rkc, Ln, Lnc, Rn, Rnc;
  TSIL_REAL mq, mq2;
  TSIL_COMPLEX f1result, f2result, f3result, f4result, f5result, f6result;
  /* TSIL_COMPLEX BFV, BfV; */
  int i,j,k,n;

  if (0 == expand_around_pole) {
    Gluino = m2_gluino;

    for (n=0; n<4; n++) quark[n] = 0.L;
    quark[4] = m2_bot; 
    quark[5] = m2_top;  

    mgluino = m_gluino;
    gluino = mgluino * mgluino; /* Keep logically distinct from Gluino! */ 
    for (n=0; n<6; n++) mquark[n] = TSIL_SQRT (quark[n]);

    /* d,u,s,c,b,t squark masses and mixings */
    squark[0][0] = m2_sdR[0];
    squark[0][1] = m2_sdL[0];
    squark[1][0] = m2_suR[0];
    squark[1][1] = m2_suL[0];
    squark[2][0] = m2_sdR[1];
    squark[2][1] = m2_sdL[1];
    squark[3][0] = m2_suR[1];
    squark[3][1] = m2_suL[1];
    for (n=0; n<2; n++) {
      squark[4][n] = m2_sbot[n]; 
      squark[5][n] = m2_stop[n]; 
    }
  } else {
    SUMO_Error("Pi2tilde_gluino",
               "Expanding around pole not implemented yet.",3);
  }

  for (n=0; n<4; n++) { 
    Rsq[n][0] = 1.L; 
    Lsq[n][0] = 0.L;
    Rsq[n][1] = 0.L; 
    Lsq[n][1] = 1.L;
  }

  for (n=0; n<2; n++) {
    Rsq[4][n] = Rsbot[n]; 
    Lsq[4][n] = Lsbot[n]; 
    Rsq[5][n] = Rstop[n]; 
    Lsq[5][n] = Lstop[n]; 
  }

  for (n=0; n<6; n++) {
    for (j=0; j<2; j++) {
      Rsqc[n][j] = SUMO_CONJ (Rsq[n][j]);
      Lsqc[n][j] = SUMO_CONJ (Lsq[n][j]);
      reLRc[n][j] = TSIL_CREAL(Lsq[n][j] * Rsqc[n][j]);
    }}

  /* Compute all of the required basis integrals: */
  for (i=0; i<6; i++){
    for (j=0; j<2; j++){
      for (k=j; k<2; k++){
	TSIL_SetParameters (&(int_data[i][j][k]), squark[i][j], quark[i], 
			    quark[i], squark[i][k], Gluino, Q2);
	TSIL_Evaluate (&int_data[i][j][k], Gluino);
      }}}
  
  for (i=0; i<6; i++) {
    mq = mquark[i];
    mq2 = mq * mq;
    for (j=0; j<2; j++) {
      Lj = Lsq[i][j];
      Ljc = Lsqc[i][j];
      Rj = Rsq[i][j];
      Rjc = Rsqc[i][j];
      for (k=0; k<2; k++) {
	Lk = Lsq[i][k];
	Lkc = Lsqc[i][k];
	Rk = Rsq[i][k];
	Rkc = Rsqc[i][k];
	
	if (k < j) {
	  mSFFSF (&int_data[i][k][j], 1, 
		  &MSFFSF, &MSFFSf, &MSFfSF, &MSFfSf, &MSffSF, &MSffSf);
	  vSFFFS (&int_data[i][k][j], &int_data[i][k][j], 1,
		  &VSFFFS, &VSFFfS, &VSFfFS, &VSFffS, &VSffFS, &VSfffS);
	} else {
	  mSFFSF (&int_data[i][j][k], 0, 
		  &MSFFSF, &MSFFSf, &MSFfSF, &MSFfSf, &MSffSF, &MSffSf);
	  vSFFFS (&int_data[i][j][k], &int_data[i][j][k], 0,
		  &VSFFFS, &VSFFfS, &VSFfFS, &VSFffS, &VSffFS, &VSfffS);
	}

	pi2tildeMSFFSF += Lj * Rjc * Lkc * Rk * MSFFSF 
	  - reLRc[i][k] * mgluino * mq * MSFfSF
	  + 0.5L*(Lj * Ljc * Lk * Lkc + Rj * Rjc * Rk * Rkc) * gluino * MSFFSf
	  + Lj * Ljc * Rk * Rkc * mq2 * MSffSF
	  - reLRc[i][j] * mgluino * mq * MSFfSf
	  + TSIL_CREAL(Lj * Rjc * Lk * Rkc) * gluino * mq2 * MSffSf;
	
	pi2tildeVSFFFS += (Lj * Ljc * Lk * Lkc + Rj * Rjc * Rk * Rkc) * VSFFFS
	  - 2.L * reLRc[i][j] * mgluino * mq * VSFfFS
	  + (Lj * Rjc * Lkc * Rk + Ljc * Rj * Lk * Rkc) * gluino * VSFFfS 
	  + (Lj * Ljc * Rk * Rkc + Rj * Rjc * Lk * Lkc) * mq2 * VSffFS
	  - 2.L * reLRc[i][k] * mgluino * mq * VSFffS
	  + (Lj * Rjc * Lk * Rkc + Ljc * Rj * Lkc * Rk) * gluino * mq2 * VSfffS; 
	
	vFSSFF (&int_data[i][j][1], &int_data[i][k][1], 
		&VFSSFF, &VFSSff, &VfSSFF, &VfSSff);

	if (j==k) {
	  pi2tildeVFSSFF += VFSSFF - 2.L * reLRc[i][j] * mgluino * mq * VfSSFF;
	}
	
	pi2tildeVFSSFF += -2.0L * reLRc[i][j] * mgluino * mq * VFSSff 
	  + 2.0L * TSIL_CREAL(Lj * Rkc * (Rjc * Lk + Ljc * Rk)) *
	  gluino * mq2 * VfSSff;
	
	for (n=0; n<2; n++){
	  Ln = Lsq[i][n];
	  Lnc = Lsqc[i][n];
	  Rn = Rsq[i][n];
	  Rnc = Rsqc[i][n];

	  yFSSS(quark[i], squark[i][j], squark[i][k], squark[i][n], 
		Gluino, Q2, &YFSSS, &YfSSS);

	  if (j==k){
	    pi2tildeYFSSS += 0.5L * SUMO_AbsSq (Lj * Lnc - Rj * Rnc) * YFSSS;
	  }
 
	  pi2tildeYFSSS += -TSIL_CREAL(Ljc * Rk * 
        (Lj * Lnc - Rj * Rnc)*(Lkc * Ln - Rkc * Rn)) * mgluino * mq * YfSSS;
	}    

      }}}

  pi2tildeMSFFSF *= 8.L * Iq * (2.L * Cq - CG);
  pi2tildeVSFFFS *= 8.L * Iq * Cq;
  pi2tildeVFSSFF *= 8.L * Iq * Cq;
  pi2tildeYFSSS *= 8.L * Iq * Cq;

  for (n=0; n<6; n++){
    for (j=0; j<2; j++){
      bFS (quark[n], squark[n][j], gluino, Q2, &BFS[n][j], &BfS[n][j]);
      bpFS (quark[n], squark[n][j], gluino, Q2, &BpFS[n][j], &BpfS[n][j]);
    }}

  for (i=0; i<6; i++){
    for (j=0; j<2; j++){
      for (n=0; n<6; n++){
	for (k=0; k<2; k++){	  
	  pi1leftover += 4.0L * (BFS[i][j] - 2.L * reLRc[i][j] * 
				 mgluino * mquark[i] * BfS[i][j]) *
	    (BpFS[n][k] - 2.L * reLRc[n][k] * mgluino * mquark[n] * BpfS[n][k]);

	  pi1leftover += -BFS[i][j] * BFS[n][k]/gluino;

	  pi1leftover += 4.0L * reLRc[i][j] * reLRc[n][k] * 
	    mquark[i] * mquark[n] * BfS[i][j] * BfS[n][k];
	}}}}

  pi1leftover *= 4.0L * Iq * Iq;

  for (i=0; i<6; i++){
    for (j=0; j<2; j++){
      f1f4(Gluino, quark[i], squark[i][j], Q2, &f1result, &f4result);
      pi2tilde1vector += (4.L * Cq - 2.L * CG) * Iq * (f1result 
             - 2.L * reLRc[i][j] * mgluino * mquark[i] * f4result);

      f2f3f5f6(Gluino, quark[i], squark[i][j], Q2, 
	       &f2result, &f3result, &f5result, &f6result);
      pi2tilde1vector += 4.L * CG * Iq * (f3result 
             - 2.L * reLRc[i][j] * mgluino * mquark[i] * f6result);
      /* f2 and f5 are not actually used for anything in this example. */
    }}

  pi2tilde2vector = CG * CG * (F1(Gluino, Q2) + F2(Gluino, Q2) 
			       + F3(Gluino, Gluino, Q2));
  
  for (n=0; n<6; n++) {
    pi2tilde2vector += 2.L * CG * Iq * F3(Gluino, quark[n], Q2);
    for (j=0; j<2; j++) {
      pi2tilde2vector += CG * Iq * F4(Gluino, squark[n][j], Q2);
    }}

  return (SUMO_twoloopfactor*g3*g3*g3*g3* (pi2tildeMSFFSF + 
          pi2tildeVSFFFS + pi2tildeVFSSFF + pi2tildeYFSSS + 
          pi1leftover + pi2tilde1vector + pi2tilde2vector));
}
