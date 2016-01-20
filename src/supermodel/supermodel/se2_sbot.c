/* Two loop self energy function for sbot */

#include "supermodel.h"
#include "self_scalar.h"

/* Facts about QCD: */
#define CG 3.L
#define Cq (4.L/3.L)
#define Iq 0.5L

/* Defined in se_squark.c */
TSIL_COMPLEX pi1_sbot (int, int, TSIL_REAL);
TSIL_COMPLEX Pi1_sbot (int, int, TSIL_REAL);

/* Used in multiple routines below */
TSIL_COMPLEX Psbot[2][2], Nsbot[2][2];
/* int arePandNsbotSet = NO; */

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* hep-ph/0502168 eqs. (4.5) and (4.6). */

void SetPandNsbot () 
{
  int i,j;
  for (i=0; i<2; i++)
    for (j=0; j<2; j++) {
      Psbot[i][j] = Lsbot[i]*Lsbotc[j] - Rsbot[i]*Rsbotc[j];
      Nsbot[i][j] = Lsbot[i]*Rsbotc[j] + Rsbot[i]*Lsbotc[j];
      /* printf("Psbot[%d][%d] = ", i,j);TSIL_cprintf(Psbot[i][j]);printf("\n"); */
      /* printf("Nsbot[%d][%d] = ", i,j);TSIL_cprintf(Nsbot[i][j]);printf("\n"); */
    }
  /* arePandNsbotSet = YES; */
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* 0502168 eq. (4.8) */

TSIL_COMPLEX pi20_sbot (int i, int j, TSIL_REAL s) 
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
  if (arePandNsbotSet == NO) SetPandNsbot ();
*/
  SetPandNsbot ();

  /* 0502168 eq. (4.8), term by term: */
  for (k=0; k<2; k++)
    for (m=0; m<2; m++)
      for (n=0; n<2; n++) {
	term += Psbot[i][k] * Psbot[k][m] * Psbot[m][n] * Psbot[n][j] *
	  X_SSS (m2_sbot[k],m2_sbot[n],m2_sbot[m],Q2);
      }
  result += 0.25L*term;

  term = 0.0L;
  for (k=0; k<2; k++)
    term += Psbot[i][k] * Psbot[k][j]*
      W_SSFF (m2_sbot[k],m2_sbot[k],m2_bot,m2_gluino,Q2);

  /* NOTE m_bot = 0 for bot != top */

  /* for (k=0; k<2; k++) */
  /*   for (m=0; m<2; m++) { */
  /*     term += -Psbot[i][k] * Nsbot[k][m] * Psbot[m][j] * m_bot * m_gluino * */
  /* 	W_SSff (m2_sbot[k],m2_sbot[m],m2_bot,m2_gluino,Q2); */
  /*   } */
  result += 0.5*term;

  /* Now term will be the second part, multiplied by (4Cq^2 - 2CG Cq): */
  term = 0.0L;

  for (k=0; k<2; k++) {
    TSIL_SetParameters (&bar,m2_gluino,m2_bot,m2_bot,m2_gluino,m2_sbot[k],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    result += (Lsbot[i] * Lsbotc[j] * Lsbot[k] * Lsbotc[k] + 
	       Rsbot[i] * Rsbotc[j] * Rsbot[k] * Rsbotc[k]) * vFFFFS
      + (Lsbot[i] * Rsbotc[j] * Rsbot[k] * Lsbotc[k] + 
	 Rsbot[i] * Lsbotc[j] * Lsbot[k] * Rsbotc[k]) * m2_gluino * vfFFfS
      /* + (Lsbot[i] * Lsbotc[j] * Rsbot[k] * Rsbotc[k] +  */
      /* 	 Rsbot[i] * Rsbotc[j] * Lsbot[k] * Lsbotc[k]) * m2_bot * vFffFS */
      /* - Nsbot[i][j] * m_bot * m_gluino * vfFfFS */
      /* + (Lsbot[i] * Rsbotc[j] * Lsbot[k] * Rsbotc[k] +  */
      /* 	 Rsbot[i] * Lsbotc[j] * Rsbot[k] * Lsbotc[k]) * m2_bot * m2_gluino * vffffS */
      ;

    /* if (i==j) */
    /*   result += -Nsbot[k][k] * m_bot * m_gluino * vFFffS; */

    TSIL_PermuteResult (&gaak, XYandZU, &gaak2);
    M_FFFFS (&gaak2, &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS);

    term += (Lsbot[i] * Rsbotc[j] * Rsbot[k] * Lsbotc[k] + 
	     Rsbot[i] * Lsbotc[j] * Lsbot[k] * Rsbotc[k]) * mFFFFS
      /* - Nsbot[i][j] * m_bot * m_gluino * mFFffS */
      + (Lsbot[i] * Lsbotc[j] * Lsbot[k] * Lsbotc[k]
      + Rsbot[i] * Rsbotc[j] * Rsbot[k] * Rsbotc[k]) * m2_gluino * mFffFS
      /* + (Lsbot[i] * Rsbotc[j] * Lsbot[k] * Rsbotc[k] +  */
      /* 	 Rsbot[i] * Lsbotc[j] * Rsbot[k] * Lsbotc[k]) * m2_bot * m2_gluino * mffffS */
      ;

    /* if (i==j) */
    /*   term += -Nsbot[k][k] * m_bot * m_gluino * mFfFfS; */

    /* M_FFFFS (&gaak, &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS); */
    /* term += (Lsbot[i] * Lsbotc[j] * Rsbot[k] * Rsbotc[k] +  */
    /* 	     Rsbot[i] * Rsbotc[j] * Lsbot[k] * Lsbotc[k]) * m2_bot * mFffFS; */

    for (m=0; m<2; m++)
      for (n=0; n<2; n++) {
	TSIL_SetParametersST (&bar,m2_sbot[k],m2_sbot[m],m2_sbot[n],Q2);
	TSIL_Evaluate (&bar, s);
	/* For clarity (we just need this one function): */
	sSSS = -bar.S[uxv].value;
	term += 0.25L * Psbot[i][k] * Psbot[k][m] * Psbot[m][n] * Psbot[n][j] * sSSS;
      }
  }
  result *= 4.0L*Cq*Cq;
  result += term*(4.0L*Cq*Cq - 2.0L*CG*Cq);

  /* Next bit is last three lines of eq. (4.8).  These involve summing
     over all 12 quark/squark mass eigenstates. */

  term = 0.0L;

  /* stop */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_bot,m2_gluino,m2_gluino,m2_top,m2_stop[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    term += -2.0L * Nsbot[i][j] * m_top * m_gluino * vfFfFS
      + Nsbot[i][j] * Nsbot[r][r] * m_bot * m_top * (vfFFfS + m2_gluino*vffffS);

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS - 2.0L*Nsbot[r][r]*m_top*m_gluino*vFFffS;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sbot[k],m2_stop[r],m2_stop[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Psbot[i][k] * Psbot[k][j] * Psbot[r][t] * Psbot[t][r] * sSSS;
      }
  }

  /* sbot */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_bot,m2_gluino,m2_gluino,0.L,m2_sbot[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Take mbot ~ 0 */
    /* term = -2.0L*Nsbot[i][j]*m_bot*m_gluino*vfFfFS */
    /*   + Nsbot[i][j]*Nsbot[r][r]*m_bot*m_bot*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nsbot[r][r]*m_bot*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sbot[k],m2_sbot[r],m2_sbot[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Psbot[i][k] * Psbot[k][j] * Psbot[r][t] * Psbot[t][r] * sSSS;
      }
  }

  /* suL */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_bot,m2_gluino,m2_gluino,0.L,m2_suL[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*Nsbot[i][j]*m_bot*m_gluino*vfFfFS */
    /*   + Nsbot[i][j]*Nsbot[r][r]*m_bot*m_bot*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nsbot[r][r]*m_bot*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sbot[k],m2_suL[r],m2_suL[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Psbot[i][k] * Psbot[k][j] * Psbot[r][t] * Psbot[t][r] * sSSS;
      }
  }

  /* suR */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_bot,m2_gluino,m2_gluino,0.L,m2_suR[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*Nsbot[i][j]*m_bot*m_gluino*vfFfFS */
    /*   + Nsbot[i][j]*Nsbot[r][r]*m_bot*m_bot*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nsbot[r][r]*m_bot*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sbot[k],m2_suR[r],m2_suR[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Psbot[i][k] * Psbot[k][j] * Psbot[r][t] * Psbot[t][r] * sSSS;
      }
  }

  /* sdL */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_bot,m2_gluino,m2_gluino,0.L,m2_sdL[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*Nsbot[i][j]*m_bot*m_gluino*vfFfFS */
    /*   + Nsbot[i][j]*Nsbot[r][r]*m_bot*m_bot*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nsbot[r][r]*m_bot*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sbot[k],m2_sdL[r],m2_sdL[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Psbot[i][k] * Psbot[k][j] * Psbot[r][t] * Psbot[t][r] * sSSS;
      }
  }

  /* sdR */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,m2_bot,m2_gluino,m2_gluino,0.L,m2_sdR[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*Nsbot[i][j]*m_bot*m_gluino*vfFfFS */
    /*   + Nsbot[i][j]*Nsbot[r][r]*m_bot*m_bot*(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*Nsbot[r][r]*m_bot*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sbot[k],m2_sdR[r],m2_sdR[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * Psbot[i][k] * Psbot[k][j] * Psbot[r][t] * Psbot[t][r] * sSSS;
      }
  }
  result += 4.0L*Cq*Iq*term;

  result *= g3*g3*g3*g3;
  /* End of eq. (4.8) */

  return result;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Eq. (4.9), tilde version (i.e., with j=i and tilded versions of
   G_SSS, G_SSFF, and G_SSff). 

   Note s arg technically redundant since this should only ever be
   evaluated for s = m2_sbot[i]. */

TSIL_COMPLEX pi21tilde_sbot (int i, TSIL_REAL s) 
{
  int k;
  TSIL_COMPLEX result = 0.0L;
  TSIL_DATA bar;
  TSIL_RESULT gaak;
  TSIL_COMPLEX gFF,gff;
  TSIL_COMPLEX gSSFF,gSSff;

/* SPM Jan. 14 2016. Apparently this should be NO sometimes, but isn't.
  if (arePandNsbotSet == NO) SetPandNsbot ();
*/
  SetPandNsbot ();

  TSIL_SetParameters (&bar,m2_bot,m2_bot,m2_gluino,m2_gluino,0.0L,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak); 
  G_FF (&gaak, &gFF, &gff);

  TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_sbot[i],m2_bot,m2_gluino,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  /* NOTE m_bot = 0 for bot != top */
  result += Cq*CG*(gFF + gSSFF 
		   - Nsbot[i][i]*m_bot*m_gluino*(gff + gSSff)
		   );

  TSIL_SetParameters (&bar,0.0L,m2_bot,m2_sbot[i],m2_gluino,m2_bot,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  result += (2.0L*Cq - CG)*Cq*(gSSFF
			       - Nsbot[i][i]*m_bot*m_gluino*gSSff
			       );

  for (k=0; k<2; k++)
    result += Cq * Cq * Psbot[i][k] * Psbot[k][i] *
      (G_S(m2_sbot[k],Q2) + Gtilde_SSS(m2_sbot[i],m2_sbot[k],Q2));

  result *= g3*g3*g3*g3;

  return result;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Eq. (4.10), tilde version */
/* DGR s arg is useless for same reason as above... */

TSIL_COMPLEX pi22tilde_sbot (int i, TSIL_REAL s) 
{
  int r;
  TSIL_COMPLEX result, term;

  result = Cq*F1tilde(m2_sbot[i], Q2) 
         + CG*F2tilde(m2_sbot[i], Q2)
         + CG*F3tilde(m2_sbot[i], m2_gluino, Q2);

  /* Now do final sum over r. First is r = top. Note that the two F3
     terms are the same for r=top, and all remaining F3 terms are
     identical for m_quark = 0. */
  term = 2.0L*F3tilde(m2_sbot[i], m2_top, Q2) 
      + 10.0L*F3tilde(m2_sbot[i], 0.0L, Q2);

  for (r=0; r<2; r++)
    term += F4tilde(m2_sbot[i], m2_stop[r], Q2)
          + F4tilde(m2_sbot[i], m2_sbot[r], Q2)
          + F4tilde(m2_sbot[i], m2_suL[r], Q2)
          + F4tilde(m2_sbot[i], m2_suR[r], Q2)
          + F4tilde(m2_sbot[i], m2_sdL[r], Q2)
          + F4tilde(m2_sbot[i], m2_sdR[r], Q2);

  /* Final assembly: */
  result += Iq*term;
  result *= g3*g3*g3*g3*Cq;

  return result;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* eq. (4.12); valid for s = m2_sbot[i] */

TSIL_COMPLEX pi20tilde_sbot (int i, TSIL_REAL s) 
{
/* SPM Jan. 14 2016. Apparently this should be NO sometimes, but isn't.
  if (arePandNsbotSet == NO) SetPandNsbot ();
*/
  SetPandNsbot ();

  return pi20_sbot (i, i, s) 
    + pi1_sbot (i, i, s) * 2.0L * g3 * g3 * Cq * (BpFF(m2_bot,m2_gluino,s,Q2) 
    - Nsbot[i][i] * m_bot * m_gluino * Bpff(m2_bot,m2_gluino,s,Q2)
						  );
}
