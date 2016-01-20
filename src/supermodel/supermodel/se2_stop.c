/* Two loop self energy function for stop */

#include "supermodel.h"
#include "self_scalar.h"

/* Facts about QCD: */
#define CG 3.L
#define Cq (4.L/3.L)
#define Iq 0.5L

/* Defined in se_squark.c */
TSIL_COMPLEX pi1_stop (int, int, TSIL_REAL);
TSIL_COMPLEX Pi1_stop (int, int, TSIL_REAL);

/* Used in multiple routines below */
TSIL_COMPLEX Pstop[2][2], Nstop[2][2];
/* int arePandNstopSet = NO; */

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* hep-ph/0502168 eqs. (4.5) and (4.6). */

void SetPandNstop ()
{
  int i,j;
  for (i=0; i<2; i++)
    for (j=0; j<2; j++) {
      Pstop[i][j] = Lstop[i]*Lstopc[j] - Rstop[i]*Rstopc[j];
      Nstop[i][j] = Lstop[i]*Rstopc[j] + Rstop[i]*Lstopc[j];
    }

  /* arePandNstopSet = YES; */
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* 0502168 eq. (4.8) */

TSIL_COMPLEX pi20_stop (int i, int j, TSIL_REAL s) 
{
  int k,m,n,r,t;
  TSIL_COMPLEX result = 0.0L;
  TSIL_COMPLEX term = 0.0L;
  TSIL_DATA bar;
  TSIL_RESULT gaak,gaak2;
  TSIL_COMPLEX vFFFFS, vFffFS, vfFfFS, vfFFfS, vFFffS, vffffS;
  TSIL_COMPLEX mFFFFS, mFFffS, mFfFfS, mFffFS, mffffS;
  TSIL_COMPLEX sSSS;

/* SPM Jan. 14 2016. Apparently this should be NO sometimes, but isn't.
  if (arePandNstopSet == NO) SetPandNstop ();
*/
  SetPandNstop ();

  /* 0502168 eq. (4.8), term by term: */
  for (k=0; k<2; k++)
    for (m=0; m<2; m++)
      for (n=0; n<2; n++) {
  	term += Pstop[i][k] * Pstop[k][m] * Pstop[m][n] * Pstop[n][j] *
  	  X_SSS (m2_stop[k],m2_stop[n],m2_stop[m],Q2);
      }
  result += 0.25L*term;

  term = 0.0L;
  for (k=0; k<2; k++)
    term += Pstop[i][k] * Pstop[k][j] *
      W_SSFF (m2_stop[k],m2_stop[k],m2_top,m2_gluino,Q2);

  for (k=0; k<2; k++)
    for (m=0; m<2; m++) {
      term += -Pstop[i][k] * Nstop[k][m] * Pstop[m][j] * m_top * m_gluino *
  	W_SSff (m2_stop[k],m2_stop[m],m2_top,m2_gluino,Q2);
    }
  result += 0.5*term;

  /* Now term will be the second part, multiplied by (4Cq^2 - 2CG Cq): */
  term = 0.0L;

  for (k=0; k<2; k++) {
    TSIL_SetParameters (&bar,m2_gluino,m2_top,m2_top,m2_gluino,m2_stop[k],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    result += (Lstop[i] * Lstopc[j] * Lstop[k] * Lstopc[k] +
    	       Rstop[i] * Rstopc[j] * Rstop[k] * Rstopc[k]) * vFFFFS
      + (Lstop[i] * Rstopc[j] * Rstop[k] * Lstopc[k] +
    	 Rstop[i] * Lstopc[j] * Lstop[k] * Rstopc[k]) * m2_gluino * vfFFfS
      + (Lstop[i] * Lstopc[j] * Rstop[k] * Rstopc[k] +
    	 Rstop[i] * Rstopc[j] * Lstop[k] * Lstopc[k]) * m2_top * vFffFS
      - Nstop[i][j] * m_top * m_gluino * vfFfFS
      + (Lstop[i] * Rstopc[j] * Lstop[k] * Rstopc[k] +
    	 Rstop[i] * Lstopc[j] * Rstop[k] * Lstopc[k]) * m2_top * m2_gluino * vffffS;

    if (i==j)
      result += -Nstop[k][k] * m_top * m_gluino * vFFffS;

    TSIL_PermuteResult (&gaak, XYandZU, &gaak2);
    M_FFFFS (&gaak2, &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS);

    term += (Lstop[i] * Rstopc[j] * Rstop[k] * Lstopc[k] + 
	     Rstop[i] * Lstopc[j] * Lstop[k] * Rstopc[k]) * mFFFFS
      - Nstop[i][j] * m_top * m_gluino * mFFffS
      + (Lstop[i] * Lstopc[j] * Lstop[k] * Lstopc[k] + 
	 Rstop[i] * Rstopc[j] * Rstop[k] * Rstopc[k]) * m2_gluino * mFffFS
      + (Lstop[i] * Rstopc[j] * Lstop[k] * Rstopc[k] + 
	 Rstop[i] * Lstopc[j] * Rstop[k] * Lstopc[k]) * m2_top * m2_gluino * mffffS;

    if (i==j)
      term += -Nstop[k][k]*m_top*m_gluino*mFfFfS;

    M_FFFFS (&gaak, &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS);
    term += (Lstop[i] * Lstopc[j] * Rstop[k] * Rstopc[k] + 
	     Rstop[i] * Rstopc[j] * Lstop[k] * Lstopc[k]) * m2_top * mFffFS;

    for (m=0; m<2; m++)
      for (n=0; n<2; n++) {
	TSIL_SetParametersST (&bar,m2_stop[k],m2_stop[m],m2_stop[n],Q2);
	TSIL_Evaluate (&bar, s);
	/* For clarity (we just need this one function): */
	sSSS = -bar.S[uxv].value;
	term += 0.25L * Pstop[i][k] * Pstop[k][m] * Pstop[m][n] * Pstop[n][j] * sSSS;
      }
  }
  result *= 4.0L*Cq*Cq;
  result += term*(4.0L*Cq*Cq - 2.0L*CG*Cq);

  /* Next bit is last three lines of eq. (4.8).  These involve summing
     over all 12 quark/squark mass eigenstates. */

  term = 0.0L;

  /* stop */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_top,m2_gluino,m2_gluino,m2_top,m2_stop[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    term += -2.0L * Nstop[i][j] * m_top * m_gluino * vfFfFS
      + Nstop[i][j] * Nstop[r][r] * m_top * m_top * (vfFFfS + m2_gluino*vffffS);

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS - 2.0L*Nstop[r][r]*m_top*m_gluino*vFFffS;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_stop[k],m2_stop[r],m2_stop[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Pstop[i][k] * Pstop[k][j] * Pstop[r][t] * Pstop[t][r] * sSSS;
      }
  }

  /* sbot */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_top,m2_gluino,m2_gluino,0.L,m2_sbot[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Take mbot ~ 0 */
    /* term += -2.0L*Nstop[i][j]*m_top*m_gluino*vfFfFS */
    /*   + Nstop[i][j]*Nstop[r][r]*m_top*m_top*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nstop[r][r]*m_top*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_stop[k],m2_sbot[r],m2_sbot[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Pstop[i][k] * Pstop[k][j] * Pstop[r][t] * Pstop[t][r] * sSSS;
      }
  }

  /* suL */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_top,m2_gluino,m2_gluino,0.L,m2_suL[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*Nstop[i][j]*m_top*m_gluino*vfFfFS */
    /*   + Nstop[i][j]*Nstop[r][r]*m_top*m_top*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nstop[r][r]*m_top*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_stop[k],m2_suL[r],m2_suL[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Pstop[i][k] * Pstop[k][j] * Pstop[r][t] * Pstop[t][r] * sSSS;
      }
  }

  /* suR */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_top,m2_gluino,m2_gluino,0.L,m2_suR[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*Nstop[i][j]*m_top*m_gluino*vfFfFS */
    /*   + Nstop[i][j]*Nstop[r][r]*m_top*m_top*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nstop[r][r]*m_top*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_stop[k],m2_suR[r],m2_suR[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Pstop[i][k] * Pstop[k][j] * Pstop[r][t] * Pstop[t][r] * sSSS;
      }
  }

  /* sdL */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_top,m2_gluino,m2_gluino,0.L,m2_sdL[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*Nstop[i][j]*m_top*m_gluino*vfFfFS */
    /*   + Nstop[i][j]*Nstop[r][r]*m_top*m_top*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nstop[r][r]*m_top*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_stop[k],m2_sdL[r],m2_sdL[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Pstop[i][k] * Pstop[k][j] * Pstop[r][t] * Pstop[t][r] * sSSS;
      }
  }

  /* sdR */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_top,m2_gluino,m2_gluino,0.L,m2_sdR[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*Nstop[i][j]*m_top*m_gluino*vfFfFS */
    /*   + Nstop[i][j]*Nstop[r][r]*m_top*m_top*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nstop[r][r]*m_top*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_stop[k],m2_sdR[r],m2_sdR[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Pstop[i][k] * Pstop[k][j] * Pstop[r][t] * Pstop[t][r] * sSSS;
      }
  }
  result += 4.0L*Cq*Iq*term;

  result *= g3*g3*g3*g3;
  /* End of eq. (4.8) */

  return result;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* 0502168 eq. (4.9) */

TSIL_COMPLEX pi21_stop (int i, int j, TSIL_REAL s) 
{
  int k;
  TSIL_COMPLEX result = 0.0L;
  TSIL_DATA bar;
  TSIL_RESULT gaak;
  TSIL_COMPLEX gFF,gff;
  TSIL_COMPLEX gSSFF,gSSff;

/* SPM Jan. 14 2016. Apparently this should be NO sometimes, but isn't.
  if (arePandNstopSet == NO) SetPandNstop ();
*/
  SetPandNstop ();

  TSIL_SetParameters (&bar,m2_top,m2_top,m2_gluino,m2_gluino,0.0L,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak); 
  G_FF (&gaak, &gFF, &gff);

  TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_stop[i],m2_top,m2_gluino,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  result += Cq*CG*(gFF + gSSFF - Nstop[i][i]*m_top*m_gluino*(gff + gSSff));

  TSIL_SetParameters (&bar,0.0L,m2_top,m2_stop[i],m2_gluino,m2_top,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  result += (2.0L*Cq - CG)*Cq*(gSSFF - Nstop[i][i]*m_top*m_gluino*gSSff);

  for (k=0; k<2; k++)
    result += Cq * Cq * Pstop[i][k] * Pstop[k][i] *
      (G_S(m2_stop[k],Q2) + Gtilde_SSS(m2_stop[i],m2_stop[k],Q2));

  result *= g3*g3*g3*g3;

  return result;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Eq. (4.9), tilde version (i.e., with j=i and tilded versions of
   G_SSS, G_SSFF, and G_SSff). 

   Note s arg technically redundant since this should only ever be
   evaluated for s = m2_stop[i]. */

TSIL_COMPLEX pi21tilde_stop (int i, TSIL_REAL s) 
{
  int k;
  TSIL_COMPLEX result = 0.0L;
  TSIL_DATA bar;
  TSIL_RESULT gaak;
  TSIL_COMPLEX gFF,gff;
  TSIL_COMPLEX gSSFF,gSSff;

/* SPM Jan. 14 2016. Apparently this should be NO sometimes, but isn't.
  if (arePandNstopSet == NO) SetPandNstop ();
*/
  SetPandNstop ();

  TSIL_SetParameters (&bar,m2_top,m2_top,m2_gluino,m2_gluino,0.0L,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak); 
  G_FF (&gaak, &gFF, &gff);

  TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_stop[i],m2_top,m2_gluino,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  result += Cq*CG*(gFF + gSSFF - Nstop[i][i]*m_top*m_gluino*(gff + gSSff));

  TSIL_SetParameters (&bar,0.0L,m2_top,m2_stop[i],m2_gluino,m2_top,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  result += (2.0L*Cq - CG)*Cq*(gSSFF - Nstop[i][i]*m_top*m_gluino*gSSff);

  for (k=0; k<2; k++)
    result += Cq * Cq * Pstop[i][k] * Pstop[k][i] *
      (G_S(m2_stop[k],Q2) + Gtilde_SSS(m2_stop[i],m2_stop[k],Q2));

  result *= g3*g3*g3*g3;

  return result;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Eq. (4.10), tilde version */
/* DGR s arg is useless for same reason as above... */

TSIL_COMPLEX pi22tilde_stop (int i, TSIL_REAL s) 
{
  int r;
  TSIL_COMPLEX result, term;

  result = Cq*F1tilde(m2_stop[i], Q2) 
         + CG*F2tilde(m2_stop[i], Q2)
         + CG*F3tilde(m2_stop[i], m2_gluino, Q2);

  /* Now do final sum over r. First is r = top. Note that the two F3
     terms are the same for r=top, and all remaining F3 terms are
     identical for m_quark = 0. */
  term = 2.0L*F3tilde(m2_stop[i], m2_top, Q2) 
      + 10.0L*F3tilde(m2_stop[i], 0.0L, Q2);

  for (r=0; r<2; r++)
    term += F4tilde(m2_stop[i], m2_stop[r], Q2)
          + F4tilde(m2_stop[i], m2_sbot[r], Q2)
          + F4tilde(m2_stop[i], m2_suL[r], Q2)
          + F4tilde(m2_stop[i], m2_suR[r], Q2)
          + F4tilde(m2_stop[i], m2_sdL[r], Q2)
          + F4tilde(m2_stop[i], m2_sdR[r], Q2);

  /* Final assembly: */
  result += Iq*term;
  result *= g3*g3*g3*g3*Cq;

  return result;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* eq. (4.12); valid for s = m2_squark[i] */

TSIL_COMPLEX pi20tilde_stop (int i, TSIL_REAL s) 
{
/* SPM Jan. 14 2016. Apparently this should be NO sometimes, but isn't.
  if (arePandNstopSet == NO) SetPandNstop ();
*/
  SetPandNstop ();

  return pi20_stop(i, i, s)
    + pi1_stop(i, i, s)*2.0L*g3*g3*Cq*(BpFF(m2_top,m2_gluino,s,Q2)
     - Nstop[i][i]*m_top*m_gluino*Bpff(m2_top,m2_gluino,s,Q2));
}
