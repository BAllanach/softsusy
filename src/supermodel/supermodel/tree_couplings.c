/* 
  Computes (selected) 3-particle and 4-particle couplings in MSSM.
*/

#include "supermodel.h"

#define SQRT8 2.82842712474619009760338L

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

int SUMO_Tree_Couplings (void)
{
  if (NO == are_tree_masses_updated) SUMO_Tree_Masses ();

  SUMO_Tree_SSS_Couplings ();
  SUMO_Tree_SSSS_Couplings ();
  SUMO_Tree_FFS_Couplings ();
  SUMO_Tree_VectorHiggs_Couplings ();
  SUMO_Tree_VectorNC_Couplings ();
  SUMO_Tree_VectorSfermion_Couplings ();
  are_tree_couplings_updated = YES;

  return 0;
}
/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
/* These are the scalar quartic couplings given in hep-ph/0405022
   eq. (2.43) et seq. */

int SUMO_Tree_SSSS_Couplings (void)
{
  int i,j,k,l,n,m;
  TSIL_COMPLEX temp;

  TSIL_COMPLEX xpstopstopc[2][2],xpsbotsbotc[2][2],xpstaustauc[2][2];
  TSIL_COMPLEX x3stopstopc[2][2],x3sbotsbotc[2][2],x3staustauc[2][2];
  TSIL_COMPLEX x3snutausnutauc;
  TSIL_COMPLEX x1stopsbotc[2][2],x1stausnutauc[2];
  TSIL_COMPLEX x1sbotstopc[2][2],x1snutaustauc[2]; /* Complex conjs */
  TSIL_COMPLEX x2stopsbotc[2][2],x2stausnutauc[2];
  TSIL_COMPLEX x2sbotstopc[2][2],x2snutaustauc[2]; /* Complex conjs */


  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.43) */
  for (i=0; i<4; i++){
  for (j=0; j<4; j++){
  for (k=0; k<4; k++){
  for (n=0; n<4; n++){
     lambda0000[i][j][k][n] = g2plusgp2 * 0.25L * (
       TSIL_CREAL(ku[i] * kuc[j] - kd[i] * kdc[j]) *
       TSIL_CREAL(ku[k] * kuc[n] - kd[k] * kdc[n]) +
       TSIL_CREAL(ku[k] * kuc[j] - kd[k] * kdc[j]) *
       TSIL_CREAL(ku[i] * kuc[n] - kd[i] * kdc[n]) +
       TSIL_CREAL(ku[n] * kuc[j] - kd[n] * kdc[j]) *
       TSIL_CREAL(ku[k] * kuc[i] - kd[k] * kdc[i]));
  }}}}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.44) */
  for (i=0; i<4; i++){
  for (j=0; j<4; j++){
  for (k=0; k<2; k++){
  for (n=0; n<2; n++){
    lambda00pm[i][j][k][n] = ku[i] * kd[j] * kdp[k] * kup[n] +
                             kuc[i] * kdc[j] * kup[k] * kdp[n] +
                             ku[j] * kd[i] * kdp[k] * kup[n] +
                             kuc[j] * kdc[i] * kup[k] * kdp[n];
    if ((i==j) && (k==n)) {
      lambda00pm[i][j][k][n] += 1.0L;
    }
    lambda00pm[i][j][k][n] *= 2.0L * g2;
    lambda00pm[i][j][k][n] += gp2 * (
      (ku[i] * kuc[j] - kd[i] * kdc[j] + ku[j] * kuc[i] - kd[j] * kdc[i]) * 
      (kup[k] * kup[n] - kdp[k] * kdp[n]));  
    lambda00pm[i][j][k][n] *= 0.125L;
  }}}}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.45) */
  for (i=0; i<2; i++){
  for (j=0; j<2; j++){
  for (k=0; k<2; k++){
  for (n=0; n<2; n++){
    lambdappmm[i][j][k][n] = g2plusgp2 * 0.25L * (
        2.0L * kup[i] * kup[j] * kup[k] * kup[n] 
      + 2.0L * kdp[i] * kdp[j] * kdp[k] * kdp[n] 
      - kup[i] * kdp[j] * kup[k] * kdp[n]
      - kdp[i] * kup[j] * kdp[k] * kup[n]
      - kdp[i] * kup[j] * kup[k] * kdp[n]
      - kup[i] * kdp[j] * kdp[k] * kup[n]);   
  }}}}  

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022, eq. (2.47) */

  for (i=0; i<4; i++) {
  for (j=0; j<4; j++) {
    temp = SUMO_CREAL(kd[i] * kdc[j] - ku[i] * kuc[j])/2.0L;
    lambda00suLsuLc[i][j] = (g2/2.0L - gp2/6.0L) * temp;
    lambda00sdLsdLc[i][j] = (-g2/2.0L - gp2/6.0L) * temp;
    lambda00snusnuc[i][j] = (g2/2.0L + gp2/2.0L) * temp;
    lambda00seLseLc[i][j] = (-g2/2.0L + gp2/2.0L) * temp;
    lambda00suRsuRc[i][j] = (2.0L * gp2/3.0L) * temp;
    lambda00sdRsdRc[i][j] = (-gp2/3.0L) * temp;
    lambda00seRseRc[i][j] = -gp2 * temp;
    lambda00snutausnutauc[i][j] = lambda00snusnuc[i][j];

/*     lambda00suLsuLc[i][j] = 0.0; */
/*     lambda00sdLsdLc[i][j] = 0.0; */
/*     lambda00snusnuc[i][j] = 0.0; */
/*     lambda00seLseLc[i][j] = 0.0; */
/*     lambda00suRsuRc[i][j] = 0.0; */
/*     lambda00sdRsdRc[i][j] = 0.0; */
/*     lambda00seRseRc[i][j] = 0.0; */
/*     lambda00snutausnutauc[i][j] = 0.0; */
  }}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022, eqs. (2.48), (2.49) and (2.50) */

  for (i=0; i<4; i++) {
  for (j=0; j<4; j++) {

    for (k=0; k<2; k++) {
    for (l=0; l<2; l++) {
      lambda00stopstopc[i][j][k][l] = 
	    lambda00suLsuLc[i][j] * Lstop[k] * Lstopc[l] +
	    lambda00suRsuRc[i][j] * Rstop[k] * Rstopc[l];

      lambda00sbotsbotc[i][j][k][l] =
	    lambda00sdLsdLc[i][j] * Lsbot[k] * Lsbotc[l] +
	    lambda00sdRsdRc[i][j] * Rsbot[k] * Rsbotc[l];

      lambda00staustauc[i][j][k][l] =
	    lambda00seLseLc[i][j] * Lstau[k] * Lstauc[l] +
	    lambda00seRseRc[i][j] * Rstau[k] * Rstauc[l];
    }}
  }}

  for (i=0; i<4; i++) {
  for (j=0; j<4; j++) {
    temp = SUMO_CREAL(ku[i] * kuc[j]) * ytop * ytop;
    for (k=0; k<2; k++)
      lambda00stopstopc[i][j][k][k] += temp;
  }}

  for (i=0; i<4; i++) {
  for (j=0; j<4; j++) {
      temp = SUMO_CREAL(kd[i] * kdc[j])* ybot * ybot;
      for (k=0; k<2; k++)
	lambda00sbotsbotc[i][j][k][k] += temp;
  }}

  for (i=0; i<4; i++) {
  for (j=0; j<4; j++) {
      temp = SUMO_CREAL(kd[i] * kdc[j])* ytau * ytau;
      for (k=0; k<2; k++)
	lambda00staustauc[i][j][k][k] += temp;
  }}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.55). */

  for (i=0; i<2; i++) {
  for (j=0; j<2; j++) {
      temp = (kup[i] * kup[j] - kdp[i] * kdp[j])/2.0L;
      lambdapmsuLsuLc[i][j] = (g2/2.0L + gp2/6.0L) * temp;
      lambdapmsdLsdLc[i][j] = (-g2/2.0L + gp2/6.0L) * temp;
      lambdapmsnusnuc[i][j] = (g2/2.0L - gp2/2.0L) * temp;
      lambdapmseLseLc[i][j] = (-g2/2.0L - gp2/2.0L) * temp;
      lambdapmsuRsuRc[i][j] = (-2.0L * gp2/3.0L) * temp;
      lambdapmsdRsdRc[i][j] = (gp2/3.0L) * temp;
      lambdapmseRseRc[i][j] = gp2 * temp;
  }}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eqs. (2.56), (2.57), (2.58) and (2.59). */

  for (i=0; i<2; i++) {
  for (j=0; j<2; j++) {

      /* Eq. (2.58) */
      lambdapmsnutausnutauc[i][j] = ytau2*kdp[i]*kdp[j] + lambdapmsnusnuc[i][j];

      for (k=0; k<2; k++) {
      for (m=0; m<2; m++) {
	  /* Eq. (2.56) */
	  lambdapmstopstopc[i][j][k][m] =
	    Rstop[k]*Rstopc[m]*(ytop2*kup[i]*kup[j] + lambdapmsuRsuRc[i][j]) +
	    Lstop[k]*Lstopc[m]*(ybot2*kdp[i]*kdp[j] + lambdapmsuLsuLc[i][j]);

	  /* Eq. (2.57) */
	  lambdapmsbotsbotc[i][j][k][m] =
	    Rsbot[k]*Rsbotc[m]*(ybot2*kdp[i]*kdp[j] + lambdapmsdRsdRc[i][j]) +
	    Lsbot[k]*Lsbotc[m]*(ytop2*kup[i]*kup[j] + lambdapmsdLsdLc[i][j]);

	  /* Eq. (2.59) */
	  lambdapmstaustauc[i][j][k][m] =
	    Lstau[k]*Lstauc[m]*lambdapmseLseLc[i][j] +
	    Rstau[k]*Rstauc[m]*(ytau2*kdp[i]*kdp[j] + lambdapmseRseRc[i][j]);
      }}
  }}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.63)-(2.65) */

  for (i=0; i<4; i++) {
  for (j=0; j<2; j++) {
      lambda0psdLsuLc[i][j] = g2*(kuc[i]*kup[j] + kd[i]*kdp[j])/(2.0L*SQRT2);
      lambda0pseLsnuc[i][j] = lambda0psdLsuLc[i][j];

      for (k=0; k<2; k++) {
	lambda0pstausnutauc[i][j][k] = Lstau[k]*(lambda0pseLsnuc[i][j]
						 - ytau2*kd[i]*kdp[j]/SQRT2);
	for (m=0; m<2; m++)
	  lambda0psbotstopc[i][j][k][m] = Lsbot[k]*Lstopc[m]*
	    (lambda0psdLsuLc[i][j] - (ytop2*kuc[i]*kup[j] + ybot2*kd[i]*kdp[j])/SQRT2)
	    - Rsbot[k]*Rstopc[m]*ytop*ybot*(ku[i]*kdp[j] + kdc[i]*kup[j])/SQRT2;
      }
  }}

  /* ------------------------------------------------------------------ */  
  /*  hep-ph/0405022 eqs. (2.66)-(2.89). */

  for (j=0; j<2; j++) {

    x1stausnutauc[j] = Lstau[j]/2.0L;
    x2stausnutauc[j] = -I*Lstau[j]/2.0L;
    /* Complex conjs for simplicity: */
    x1snutaustauc[j] = TSIL_CONJ(x1stausnutauc[j]);
    x2snutaustauc[j] = TSIL_CONJ(x2stausnutauc[j]);

    for (k=0; k<2; k++){
      xpstopstopc[j][k] = Lstop[j]*Lstopc[k]/6.0L - 2.0L*Rstop[j]*Rstopc[k]/3.0L;
      xpsbotsbotc[j][k] = Lsbot[j]*Lsbotc[k]/6.0L + Rsbot[j]*Rsbotc[k]/3.0L;
      xpstaustauc[j][k] = -Lstau[j]*Lstauc[k]/2.0L + Rstau[j]*Rstauc[k];
      x1stopsbotc[j][k] = Lstop[j]*Lsbotc[k]/2.0L;
      x2stopsbotc[j][k] = I*x1stopsbotc[j][k];
      /* Complex conjs for simplicity: */
      x1sbotstopc[k][j] = TSIL_CONJ(x1stopsbotc[j][k]);
      x2sbotstopc[k][j] = TSIL_CONJ(x2stopsbotc[j][k]);

      x3stopstopc[j][k] = Lstop[j]*Lstopc[k]/2.0L;
      x3sbotsbotc[j][k] = -Lsbot[j]*Lsbotc[k]/2.0L;
      x3staustauc[j][k] = -Lstau[j]*Lstauc[k]/2.0L;
      x3snutausnutauc = 0.5L;
    }}

  /* F-term contributions: */
  for (i=0; i<2; i++){
    for (j=0; j<2; j++){
      
      lambdastausnutaucsnutaustauc[i][j] = ytau2*Rstauc[j]*Rstau[i];

      for (k=0; k<2; k++){
	for (m=0; m<2; m++){
	  lambdastopstopcstopstopc[i][j][k][m] =
	    ytop2*(Lstop[i]*Rstopc[j]*Rstop[k]*Lstopc[m] + Rstop[i]*Lstopc[j]*Lstop[k]*Rstopc[m]);

	  lambdasbotsbotcsbotsbotc[i][j][k][m] =
	    ybot2*(Lsbot[i]*Rsbotc[j]*Rsbot[k]*Lsbotc[m] + Rsbot[i]*Lsbotc[j]*Lsbot[k]*Rsbotc[m]);

	  lambdastaustaucstaustauc[i][j][k][m] =
	    ytau2*(Lstau[i]*Rstauc[j]*Rstau[k]*Lstauc[m] + Rstau[i]*Lstauc[j]*Lstau[k]*Rstauc[m]);

	  lambdastopsbotcsbotstopc[i][j][k][m] =
	    ytop2*Rstop[i]*Lsbotc[j]*Lsbot[k]*Rstopc[m] + ybot2*Lstop[i]*Rsbotc[j]*Rsbot[k]*Lstopc[m];

	  lambdasbotsbotcstaustauc[i][j][k][m] =
	    ybot*ytau*(Lsbot[i]*Rsbotc[j]*Rstau[k]*Lstauc[m] + Rsbot[i]*Lsbotc[j]*Lstau[k]*Rstauc[m]);
	}}}}

  /* D-term contributions */
/* COMMENT FOLLOWING TO KILL D-TERMS IN 4-SFERMION COUPLINGS. */

  for (i=0; i<2; i++){
    for (j=0; j<2; j++){
      lambdastopstopcsnutausnutauc[i][j] = g2*x3stopstopc[i][j]*x3snutausnutauc;
      lambdasbotsbotcsnutausnutauc[i][j] = g2*x3sbotsbotc[i][j]*x3snutausnutauc;
      lambdastaustaucsnutausnutauc[i][j] = g2*x3staustauc[i][j]*x3snutausnutauc;

      lambdastausnutaucsnutaustauc[i][j] += 
	g2*(x1stausnutauc[i]*x1snutaustauc[j] + x2stausnutauc[i]*x2snutaustauc[j]);

      for (k=0; k<2; k++){
	for (m=0; m<2; m++){
	  lambdastopstopcstopstopc[i][j][k][m] +=
	    g2*x3stopstopc[i][j]*x3stopstopc[k][m] + gp2*xpstopstopc[i][j]*xpstopstopc[k][m];

	  lambdasbotsbotcsbotsbotc[i][j][k][m] +=
	    g2*x3sbotsbotc[i][j]*x3sbotsbotc[k][m] + gp2*xpsbotsbotc[i][j]*xpsbotsbotc[k][m];

	  lambdastaustaucstaustauc[i][j][k][m] +=
	    g2*x3staustauc[i][j]*x3staustauc[k][m] + gp2*xpstaustauc[i][j]*xpstaustauc[k][m];

	  lambdastopstopcsbotsbotc[i][j][k][m] =
	    g2*x3stopstopc[i][j]*x3sbotsbotc[k][m] + gp2*xpstopstopc[i][j]*xpsbotsbotc[k][m];

	  lambdastopstopcstaustauc[i][j][k][m] =
	    g2*x3stopstopc[i][j]*x3staustauc[k][m] + gp2*xpstopstopc[i][j]*xpstaustauc[k][m];

	  lambdasbotsbotcstaustauc[i][j][k][m] +=
	    g2*x3sbotsbotc[i][j]*x3staustauc[k][m] + gp2*xpsbotsbotc[i][j]*xpstaustauc[k][m];

	  lambdastopsbotcsbotstopc[i][j][k][m] +=
	    g2*(x1stopsbotc[i][j]*x1sbotstopc[k][m] + x2stopsbotc[i][j]*x2sbotstopc[k][m]);
	}}}}

  lambdasnutausnutaucsnutausnutauc = g2*x3snutausnutauc*x3snutausnutauc;
/* COMMENT PREVIOUS TO KILL D-TERMS IN 4-SFERMION COUPLINGS. */


  /* Here are some copies, added for convenience */
  for (i=0; i<2; i++) {
  for (j=0; j<2; j++) {
    lambdasnutaustaucstausnutauc[i][j] = lambdastausnutaucsnutaustauc[j][i];

    for (k=0; k<2; k++) {
    for (m=0; m<2; m++) {
      lambdasbotstopcstopsbotc[i][j][k][m] = lambdastopsbotcsbotstopc[k][m][i][j];

    }}}}

  return 0;
}

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
/* Scalar trilinear couplings found in hep-ph/0405022 */

int SUMO_Tree_SSS_Couplings (void)
{
  int i,j,k;
  TSIL_REAL prelambda000[4][4][4];
  TSIL_COMPLEX temp, tempc;

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.46a) */

  for (i=0; i<4; i++){
    temp = TSIL_CREAL(vu * ku[i] - vd * kd[i]) * g2plusgp2/SQRT8;
    for (j=0; j<4; j++){
      for (k=0; k<4; k++){
	prelambda000[i][j][k] = temp * 
                              TSIL_CREAL(ku[j] * kuc[k] - kd[j] * kdc[k]);
      }}}

  for (i=0; i<4; i++){
    for (j=0; j<4; j++){
      for (k=0; k<4; k++){
	lambda000[i][j][k] = prelambda000[i][j][k] +
	                          prelambda000[j][i][k] +
	                          prelambda000[k][i][j];
      }}}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.46b) */

  for (i=0; i<4; i++){
    temp = g2 * (vd * ku[i] + vu * kd[i])/SQRT8;
    tempc = SUMO_CONJ(temp);
    for (j=0; j<2; j++){
      for (k=0; k<2; k++){
	lambda0pm[i][j][k] = temp * kdp[j] * kup[k]  
                                  + tempc * kup[j] * kdp[k]; 
      }}}

  for (i=0; i<2; i++){
    temp = g2 * (vd * kd[i] + vu * ku[i])/SQRT8;
    for (j=0; j<2; j++){
      lambda0pm[i][j][j] += temp;
    }}

  for (i=0; i<2; i++){
    temp = gp2 * (vd * kd[i] - vu * ku[i])/SQRT8; 
    for (j=0; j<2; j++){
      for (k=0; k<2; k++){
	lambda0pm[i][j][k] += temp * (kdp[j] * kdp[k] - kup[j] * kup[k]);
      }}}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.51). 
     [Note sign differs from hep-ph/0206136 eq. (3.6).] */

  for (i=0; i<4; i++) {
    temp = TSIL_CREAL(vd * kd[i] - vu * ku[i])/SQRT2;
    lambda0suLsuLc[i] = (g2/2.0L - gp2/6.0L) * temp;
    lambda0sdLsdLc[i] = (-g2/2.0L - gp2/6.0L) * temp;
    lambda0suRsuRc[i] = (2.0L * gp2/3.0L) * temp;
    lambda0sdRsdRc[i] = (-gp2/3.0L) * temp;
    lambda0snusnuc[i] = (g2/2.0L + gp2/2.0L) * temp;
    lambda0seLseLc[i] = (-g2/2.0L + gp2/2.0L) * temp;
    lambda0seRseRc[i] = -gp2 * temp;
    lambda0snutausnutauc[i] = lambda0snusnuc[i];
  }

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.52). 
     [Note sign differs from hep-ph/0206136 eq. (3.7).] */

  for (i=0; i<4; i++) {
    temp = (ku[i] * atop - kdc[i] * muc * ytop)/SQRT2;
    tempc = SUMO_CONJ(temp);
    for (j=0; j<2; j++) {
    for (k=0; k<2; k++) {
	lambda0stopstopc[i][j][k] =
	  Lstop[j] * Lstopc[k] * lambda0suLsuLc[i]
        + Rstop[j] * Rstopc[k] * lambda0suRsuRc[i]
        + Lstop[j] * Rstopc[k] * temp + Rstop[j] * Lstopc[k] * tempc;
    }}
  }

  /* i loop goes to i<2 because i=2,3 have imaginary ku; see
     eq. (2.16) of hep-ph/0206136. */
  for (i=0; i<2; i++) {
  for (j=0; j<2; j++) {
    lambda0stopstopc[i][j][j] += SQRT2 * vu * ytop2 * ku[i];
  }} 

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.53). 
     [Note sign differs from hep-ph/0206136 eq. (3.8).] */

  for (i=0; i<4; i++){
    temp = (kd[i] * abot - kuc[i] * muc * ybot)/SQRT2;
    tempc = SUMO_CONJ(temp);
    for (j=0; j<2; j++){
      for (k=0; k<2; k++){
	lambda0sbotsbotc[i][j][k] =
	    Lsbot[j] * Lsbotc[k] * (lambda0sdLsdLc[i]) 
	  + Rsbot[j] * Rsbotc[k] * (lambda0sdRsdRc[i])
	  + Lsbot[j] * Rsbotc[k] * temp + Rsbot[j] * Lsbotc[k] * tempc;
      }}}

  /* i loop goes to i<2 because i=2,3 have imaginary ku;
     see eq. (2.16) of hep-ph/0206136. */
  for (i=0; i<2; i++){
  for (j=0; j<2; j++){
      lambda0sbotsbotc[i][j][j] += SQRT2 * vd * ybot2 * kd[i];
  }} 

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.54). 
     [Note sign differs from hep-ph/0206136 eq. (3.9).] */

  for (i=0; i<4; i++){
    temp = (kd[i] * atau - kuc[i] * muc * ytau)/SQRT2;
    tempc = SUMO_CONJ(temp);
    for (j=0; j<2; j++){
      for (k=0; k<2; k++){
	lambda0staustauc[i][j][k] =
            Lstau[j] * Lstauc[k] * (lambda0seLseLc[i]) 
	  + Rstau[j] * Rstauc[k] * (lambda0seRseRc[i])
	  + Lstau[j] * Rstauc[k] * temp + Rstau[j] * Lstauc[k] * tempc;
      }}}

  /* i loop goes to i<2 because i=2,3 have imaginary ku;
     see eq. (2.16) of hep-ph/0206136. */
  for (i=0; i<2; i++){
  for (j=0; j<2; j++){
      lambda0staustauc[i][j][j] += SQRT2 * vd * ytau2 * kd[i];
  }}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.60). 
     [Note sign differs from hep-ph/0206136 eq. (3.11).] */

  for (i=0; i<2; i++){
    lambdapsdLsuLc[i] = lambdapseLsnuc[i] = 
      0.5L * g2 * (kup[i] * vu + kdp[i] * vd);
  }

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.61). 
     [Note sign differs from hep-ph/0206136 eq. (3.12).] */

  for (i=0; i<2; i++){
  for (j=0; j<2; j++){
  for (k=0; k<2; k++){
	lambdapsbotstopc[i][j][k] =
	  Lsbot[j] * Lstopc[k] *
	  (lambdapsdLsuLc[i] - vu * ytop2 * kup[i] - vd * ybot2 * kdp[i])
	  - Rsbot[j] * Rstopc[k] * 
	  ytop * ybot * (kdp[i] * vu + kup[i] * vd)
	  - Lsbot[j] * Rstopc[k] * 
	  (kup[i] * atop + kdp[i] * muc * ytop)
	  - Rsbot[j] * Lstopc[k] * 
	  (kdp[i] * abotc + kup[i] * mu * ybot);

	/* printf("lambdapsbotstopc[%d][%d][%d] = ", i,j,k); */
	/* TSIL_cprintf(lambdapsbotstopc[i][j][k]); printf("\n"); */
  }}}
  /* printf("\n"); */

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.62). 
     [Note sign differs from hep-ph/0206136 eq. (3.13).] */

  for (i=0; i<2; i++){
  for (j=0; j<2; j++){
      lambdapstausnutauc[i][j] =
	Lstau[j] * (lambdapseLsnuc[i] - vd * ytau2 * kdp[i])
	-Rstau[j] * (kdp[i] * atauc + kup[i] * mu * ytau);
  }}

  return 0;
}
/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

int SUMO_Tree_FFS_Couplings (void)
{
  int i,j,k;
  TSIL_COMPLEX tempi, tempj;
 
  /* ------------------------------------------------------------------ */  
  /* fermion-fermion-Higgs couplings hep-ph/0405022 eqs. (2.1)-(2.6) */

  for (i=0; i<4; i++) {
    Yttcphi[i] = ytop * ku[i] / SQRT2;
    conYttcphi[i] = SUMO_CONJ (Yttcphi[i]);
    Ybbcphi[i] = ybot * kd[i] / SQRT2;
    conYbbcphi[i] = SUMO_CONJ (Ybbcphi[i]);
    Ytautaucphi[i] = ytau * kd[i] / SQRT2;
    conYtautaucphi[i] = SUMO_CONJ (Ytautaucphi[i]);
  }

  for (i=0; i<2; i++) {
    Ytcbphip[i] = -ytop * kup[i];
    conYtcbphip[i] = SUMO_CONJ (Ytcbphip[i]);
    Ybctphim[i] = -ybot * kdp[i];
    conYbctphim[i] = SUMO_CONJ (Ybctphim[i]);
    Ytaucnutauphim[i] = -ytau * kdp[i];
    conYtaucnutauphim[i] = SUMO_CONJ (Ytaucnutauphim[i]);
  }

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.28). 
     [Note sign differs from hep-ph/0206136 eq. (3.44).] */

  for (i=0; i<2; i++) {
  for (j=0; j<2; j++) {
  for (k=0; k<4; k++) {
	conYCC0[i][j][k] = g * (kd[k] * Vmix[i][0] * Umix[j][1]
	                      + ku[k] * Vmix[i][1] * Umix[j][0])/SQRT2;
	
	YCC0[i][j][k] = SUMO_CONJ(conYCC0[i][j][k]);
  }}}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.29). 
     [Note sign differs from hep-ph/0206136 eq. (3.45).] */

  for (i=0; i<4; i++) {
    tempi = 0.5L * (g * Nmix[i][1] - gp * Nmix[i][0]);
    for (j=0; j<4; j++) {
      tempj = 0.5L * (g * Nmix[j][1] - gp * Nmix[j][0]);
      for (k=0; k<4; k++) {
	conYNN0[i][j][k] = tempi *
          (kd[k] * Nmix[j][2] - ku[k] * Nmix[j][3])
	  + tempj * (kd[k] * Nmix[i][2] - ku[k] * Nmix[i][3]);

	YNN0[i][j][k] = SUMO_CONJ(conYNN0[i][j][k]);
      }}}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.30), (2.31). 
     [Note signs same as hep-ph/0206136 eqs. (3.46), (3.47).] */

  for (i=0; i<2; i++) {
  for (j=0; j<4; j++) {
  for (k=0; k<2; k++) {
    conYCNm[i][j][k] = kup[k] * (g * Vmix[i][0] * Nmix[j][3]
             + (g * Nmix[j][1] + gp * Nmix[j][0]) * Vmix[i][1]/SQRT2);
	
    conYCNp[i][j][k] = kdp[k] * (g * Umix[i][0] * Nmix[j][2]
             - (g * Nmix[j][1] + gp * Nmix[j][0]) * Umix[i][1]/SQRT2);

    YCNp[i][j][k] = SUMO_CONJ(conYCNp[i][j][k]);
    YCNm[i][j][k] = SUMO_CONJ(conYCNm[i][j][k]);
  }}}

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.7) through (2.19). 
     [Note signs differ from hep-ph/0206136 eqs. (3.27) through (3.33).] */

  for (i=0; i<4; i++){
    YuNsuLc[i] = (g * Nmixc[i][1] + gp * Nmixc[i][0]/3.0L)/SQRT2;
    YucNsuR[i] = -2.0L * SQRT2 * gp * Nmixc[i][0]/3.0L;
    YdNsdLc[i] = (-g * Nmixc[i][1] + gp * Nmixc[i][0]/3.0L)/SQRT2;
    YdcNsdR[i] = SQRT2 * gp * Nmixc[i][0]/3.0L;
    YeNseLc[i] = (-g * Nmixc[i][1] - gp * Nmixc[i][0])/SQRT2;
    YecNseR[i] = SQRT2 * gp * Nmixc[i][0];
    YnuNsnuc[i] = (g * Nmixc[i][1] - gp * Nmixc[i][0])/SQRT2;

    conYuNsuLc[i] = SUMO_CONJ (YuNsuLc[i]);
    conYucNsuR[i] = SUMO_CONJ (YucNsuR[i]);
    conYdNsdLc[i] = SUMO_CONJ (YdNsdLc[i]);
    conYdcNsdR[i] = SUMO_CONJ (YdcNsdR[i]);
    conYeNseLc[i] = SUMO_CONJ (YeNseLc[i]);
    conYecNseR[i] = SUMO_CONJ (YecNseR[i]);
    conYnuNsnuc[i] = SUMO_CONJ (YnuNsnuc[i]);

    for (j=0; j<2; j++){
      YtNstopc[i][j] = Lstopc[j] * YuNsuLc[i]
	                    + Rstopc[j] * ytop * Nmixc[i][3];

      YtcNstop[i][j] = Rstop[j] * YucNsuR[i]
	                    + Lstop[j] * ytop * Nmixc[i][3];

      YbNsbotc[i][j] = Lsbotc[j] * YdNsdLc[i]
	                    + Rsbotc[j] * ybot * Nmixc[i][2];

      YbcNsbot[i][j] = Rsbot[j] * YdcNsdR[i]
	                    + Lsbot[j] * ybot * Nmixc[i][2];

      YtauNstauc[i][j] = Lstauc[j] * YeNseLc[i]
	                      + Rstauc[j] * ytau * Nmixc[i][2];

      YtaucNstau[i][j] = Rstau[j] * YecNseR[i]
	                      + Lstau[j] * ytau * Nmixc[i][2];

      conYtNstopc[i][j] = SUMO_CONJ (YtNstopc[i][j]);
      conYtcNstop[i][j] = SUMO_CONJ (YtcNstop[i][j]);
      conYbNsbotc[i][j] = SUMO_CONJ (YbNsbotc[i][j]);
      conYbcNsbot[i][j] = SUMO_CONJ (YbcNsbot[i][j]);
      conYtauNstauc[i][j] = SUMO_CONJ (YtauNstauc[i][j]);
      conYtaucNstau[i][j] = SUMO_CONJ (YtaucNstau[i][j]);
    }
  }

  /* ------------------------------------------------------------------ */  
  /* hep-ph/0405022 eq. (2.20) through (2.27). 
     [Note signs differ from hep-ph/0206136 eqs. (3.35) through (3.38).] */

  for (i=0; i<2; i++){
    YdCsuLc[i] = YeCsnuc[i] = YtauCsnutauc[i] 
      = g * Vmixc[i][0];
    YuCsdLc[i] = YnuCseLc[i] = g * Umixc[i][0];

    conYdCsuLc[i] = conYeCsnuc[i] = conYtauCsnutauc[i]
      = g * Vmix[i][0];
    conYuCsdLc[i] = conYnuCseLc[i] = g * Umix[i][0];

    YtaucCsnutau[i] = -Umixc[i][1] * ytau;    
    conYtaucCsnutau[i] = SUMO_CONJ(YtaucCsnutau[i]);

    for (j=0; j<2; j++){
      YbCstopc[i][j] = Lstopc[j] * YdCsuLc[i] 
                            - Rstopc[j] * Vmixc[i][1] * ytop;    
      YbcCstop[i][j] = -Lstop[j] * Umixc[i][1] * ybot;    

      YtCsbotc[i][j] = Lsbotc[j] * YuCsdLc[i] 
                            - Rsbotc[j] * Umixc[i][1] * ybot;    
      YtcCsbot[i][j] = -Lsbot[j] * Vmixc[i][1] * ytop;    

      YnutauCstauc[i][j] = Lstauc[j] * YnuCseLc[i] 
                                - Rstauc[j] * Umixc[i][1] * ytau;    
     
      conYbCstopc[i][j] = SUMO_CONJ(YbCstopc[i][j]);
      conYbcCstop[i][j] = SUMO_CONJ(YbcCstop[i][j]);
      conYtCsbotc[i][j] = SUMO_CONJ(YtCsbotc[i][j]);
      conYtcCsbot[i][j] = SUMO_CONJ(YtcCsbot[i][j]);
      conYnutauCstauc[i][j] = SUMO_CONJ(YnutauCstauc[i][j]);
    }    
  }

  return 0;
}
/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

int SUMO_Tree_VectorHiggs_Couplings (void)
{
  int i,j;

  /* hep-ph/0405022 eqs. (2.33). */ 

  for (i=0; i<4; i++){
    for (j=0; j<4; j++){
      gZ00[i][j] = 0.5L * TSIL_SQRT(g2plusgp2) * 
	TSIL_CIMAG(ku[i] * kuc[j] - kd[i] * kdc[j]);
    }}

  /* hep-ph/0405022 eqs. (2.34). */ 

  for (i=0; i<2; i++)
    gZpm[i] = I * 0.5L * (g2 - gp2)/TSIL_SQRT(g2plusgp2);

  for (i=0; i<4; i++){
    for (j=0; j<2; j++){
      gW0p[i][j] = I * 0.5L * g * (kd[i] * kdp[j] - kuc[i] * kup[j]);
    }}

  /* hep-ph/0405022 eqs. (2.35) through (2.38). */   

  for (i=0; i<4; i++){
    gWW00[i] = 0.5L * g2;
    gZZ00[i] = 0.5L * (g2plusgp2);
  }

  for (i=0; i<2; i++){
    gWWpm[i] = 0.5L * g2;
    gZZpm[i] = 0.5L * (g2 - gp2) * (g2 - gp2)/(g2plusgp2);
  }

  /* hep-ph/0405022 eqs. (2.39), (2.40). */ 

  for (i=0; i<4; i++){
    gWW0[i] = TSIL_CREAL (vu * ku[i] + vd * kd[i])/SQRT2;
    gZZ0[i] = (g2plusgp2) * gWW0[i];
    gWW0[i] *= g2;
  }

  /* hep-ph/0405022 eqs. (2.41), (2.42). */ 

  for (i=0; i<2; i++){
    gWgammap[i] = (g * gp/TSIL_SQRT(g2 + gp2))
      * (vu * kup[i] - vd * kdp[i])/SQRT2;
    gWZp[i] = -gp * (gWgammap[i]);
    gWgammap[i] *= g;
  }

  return 0;
}
/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

int SUMO_Tree_VectorNC_Couplings (void)
{
  int i,j;

  /* hep-ph/0206136 eqs. (3.69) through (3.73) */
  
  for (i=0; i<4; i++){
  for (j=0; j<2; j++){
    OL[i][j] = Nmix[i][1] * Vmixc[j][0] - Nmix[i][3] * Vmixc[j][1]/SQRT2;
    OR[i][j] = Nmixc[i][1] * Umix[j][0] + Nmixc[i][2] * Umix[j][1]/SQRT2;
    OLc[i][j] = SUMO_CONJ(OL[i][j]);
    ORc[i][j] = SUMO_CONJ(OR[i][j]);
  }}

  for (i=0; i<2; i++){
  for (j=0; j<2; j++){
   OpL[i][j] = -Vmix[i][0] * Vmixc[j][0] -0.5L * Vmix[i][1] * Vmixc[j][1];
   OpR[i][j] = -Umixc[i][0] * Umix[j][0] -0.5L * Umixc[i][1] * Umix[j][1];
  }}
  
  for (i=0; i<2; i++){
    OpL[i][i] += gp2/g2plusgp2; 
    OpR[i][i] += gp2/g2plusgp2; 
  }

  for (i=0; i<4; i++){
  for (j=0; j<4; j++){
    OppL[i][j] = 0.5L * (Nmix[i][3] * Nmixc[j][3]
                      - Nmix[i][2] * Nmixc[j][2]);
  }}

  return 0;
}

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
/* No reference, needs checking! */

int SUMO_Tree_VectorSfermion_Couplings (void)
{
  TSIL_REAL denom;
  int i,j;

  denom = 6.0L * TSIL_SQRT(g2plusgp2);
  gZsuLsuLc = (3.0L * g2 - gp2)/denom;
  gZsuRsuRc = -4.0L * gp2/denom;
  gZsdLsdLc = (-3.0L * g2 - gp2)/denom;
  gZsdRsdRc = 2.0L * gp2/denom;
  gZseLseLc = (-3.0L * g2 + 3.0L * gp2)/denom;
  gZseRseRc = 6.0L * gp2/denom;
  gZsnusnuc = 3.0L * g2plusgp2/denom;

  for (i=0; i<2; i++) {
  for (j=0; j<2; j++) {
    gZstopstopc[i][j] = Lstop[i] * Lstopc[j] * gZsuLsuLc +
                        Rstop[i] * Rstopc[j] * gZsuRsuRc;

    gZsbotsbotc[i][j] = Lsbot[i] * Lsbotc[j] * gZsdLsdLc +
                        Rsbot[i] * Rsbotc[j] * gZsdRsdRc;

    gZstaustauc[i][j] = Lstau[i] * Lstauc[j] * gZseLseLc +
                        Rstau[i] * Rstauc[j] * gZseRseRc;
  }}

  gZZsuLsuLc = 2.0L * gZsuLsuLc * gZsuLsuLc;
  gZZsuRsuRc = 2.0L * gZsuRsuRc * gZsuRsuRc;
  gZZsdLsdLc = 2.0L * gZsdLsdLc * gZsdLsdLc;
  gZZsdRsdRc = 2.0L * gZsdRsdRc * gZsdRsdRc;
  gZZseLseLc = 2.0L * gZseLseLc * gZseLseLc;
  gZZseRseRc = 2.0L * gZseRseRc * gZseRseRc;
  gZZsnusnuc = 2.0L * gZsnusnuc * gZsnusnuc;

  for (i=0; i<2; i++) {
  for (j=0; j<2; j++) {
    gZZstopstopc[i][j] = Lstop[i] * Lstopc[j] * gZZsuLsuLc +
                         Rstop[i] * Rstopc[j] * gZZsuRsuRc;
 
    gZZsbotsbotc[i][j] = Lsbot[i] * Lsbotc[j] * gZZsdLsdLc +
                         Rsbot[i] * Rsbotc[j] * gZZsdRsdRc;

    gZZstaustauc[i][j] = Lstau[i] * Lstauc[j] * gZZseLseLc +
                         Rstau[i] * Rstauc[j] * gZZseRseRc;
  }}
  
  return 0;
}
