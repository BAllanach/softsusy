/* Two loop self energy function for sdL */

#include "supermodel.h"
#include "self_scalar.h"

/* Facts about QCD: */
#define CG 3.L
#define Cq (4.L/3.L)
#define Iq 0.5L

/* Defined in se_squark.c */
TSIL_COMPLEX pi1_sdL (int, int, TSIL_REAL);
TSIL_COMPLEX Pi1_sdL (int, int, TSIL_REAL);

/* Used in multiple routines below */
TSIL_COMPLEX PsdL[2][2], NsdL[2][2];
int arePandNsdLSet = NO;

/* For u,d,c,s, uncomment these lines: */
TSIL_COMPLEX LsdL[2], LsdLc[2], RsdL[2], RsdLc[2];

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* hep-ph/0502168 eqs. (4.5) and (4.6). */

void SetPandNsdL ()
{
  int i,j;

  LsdL[0] = LsdLc[0] = 0.0L;
  LsdL[1] = LsdLc[1] = -1.0L;
  RsdL[0] = RsdLc[0] = 1.0L;
  RsdL[1] = RsdLc[1] = 0.0L;

  for (i=0; i<2; i++)
    for (j=0; j<2; j++) {
      PsdL[i][j] = LsdL[i]*LsdLc[j] - RsdL[i]*RsdLc[j];
      NsdL[i][j] = LsdL[i]*RsdLc[j] + RsdL[i]*LsdLc[j];
    }

  arePandNsdLSet = YES;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* 0502168 eq. (4.8) */

TSIL_COMPLEX pi20_sdL (int i, int j, TSIL_REAL s) 
{
  int k,m,n,r,t;
  TSIL_COMPLEX result = 0.0L;
  TSIL_COMPLEX term = 0.0L;
  TSIL_DATA bar;
  TSIL_RESULT gaak,gaak2;
  TSIL_COMPLEX vFFFFS, vFffFS, vfFfFS, vfFFfS, vFFffS, vffffS;
  TSIL_COMPLEX mFFFFS, mFFffS, mFfFfS, mFffFS, mffffS;
  TSIL_COMPLEX sSSS;

  if (arePandNsdLSet == NO) SetPandNsdL ();

  /* 0502168 eq. (4.8), term by term: */
  for (k=0; k<2; k++)
    for (m=0; m<2; m++)
      for (n=0; n<2; n++)
	term += PsdL[i][k] * PsdL[k][m] * PsdL[m][n] * PsdL[n][j] *
	  X_SSS (m2_sdL[k],m2_sdL[n],m2_sdL[m],Q2);

  result += 0.25L*term;

  term = 0.0L;
  for (k=0; k<2; k++)
    term += PsdL[i][k] * PsdL[k][j] *
      W_SSFF (m2_sdL[k],m2_sdL[k],0.0L,m2_gluino,Q2);

  /* NOTE m_q = 0 for q != top */

  /* for (k=0; k<2; k++) */
  /*   for (m=0; m<2; m++) { */
  /*     term += -PsdL[i][k] * NsdL[k][m] * PsdL[m][j] * 0.0L * m_gluino * */
  /* 	W_SSff (m2_sdL[k],m2_sdL[m],0.0L,m2_gluino,Q2); */
  /*   } */
  result += 0.5*term;

  /* Now term will be the second part, multiplied by (4Cq^2 - 2CG Cq): */
  term = 0.0L;

  for (k=0; k<2; k++) {
    TSIL_SetParameters (&bar,m2_gluino,0.0L,0.0L,m2_gluino,m2_sdL[k],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    result += (LsdL[i] * LsdLc[j] * LsdL[k] * LsdLc[k] + 
	       RsdL[i] * RsdLc[j] * RsdL[k] * RsdLc[k]) * vFFFFS
      + (LsdL[i] * RsdLc[j] * RsdL[k] * LsdLc[k] + 
	 RsdL[i] * LsdLc[j] * LsdL[k] * RsdLc[k]) * m2_gluino * vfFFfS
      /* + (LsdL[i] * LsdLc[j] * RsdL[k] * RsdLc[k] +  */
      /* 	 RsdL[i] * RsdLc[j] * LsdL[k] * LsdLc[k]) * 0.0L * vFffFS */
      /* - NsdL[i][j] * 0.0L * m_gluino * vfFfFS */
      /* + (LsdL[i] * RsdLc[j] * LsdL[k] * RsdLc[k] +  */
      /* 	 RsdL[i] * LsdLc[j] * RsdL[k] * LsdLc[k]) * 0.0L * m2_gluino * vffffS */
      ;

    /* if (i==j) */
    /*   result += -NsdL[k][k] * 0.0L * m_gluino * vFFffS; */

    TSIL_PermuteResult (&gaak, XYandZU, &gaak2);
    M_FFFFS (&gaak2, &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS);

    term += (LsdL[i] * RsdLc[j] * RsdL[k] * LsdLc[k] + 
	     RsdL[i] * LsdLc[j] * LsdL[k] * RsdLc[k]) * mFFFFS
      /* - NsdL[i][j] * 0.0L * m_gluino * mFFffS */
      + (LsdL[i] * LsdLc[j] * LsdL[k] * LsdLc[k]
      + RsdL[i] * RsdLc[j] * RsdL[k] * RsdLc[k]) * m2_gluino * mFffFS
      /* + (LsdL[i] * RsdLc[j] * LsdL[k] * RsdLc[k] +  */
      /* 	 RsdL[i] * LsdLc[j] * RsdL[k] * LsdLc[k]) * 0.0L * m2_gluino * mffffS */
      ;

    /* if (i==j) */
    /*   term += -NsdL[k][k] * 0.0L * m_gluino * mFfFfS; */

    /* M_FFFFS (&gaak, &mFFFFS, &mFFffS, &mFfFfS, &mFffFS, &mffffS); */
    /* term += (LsdL[i] * LsdLc[j] * RsdL[k] * RsdLc[k] +  */
    /* 	     RsdL[i] * RsdLc[j] * LsdL[k] * LsdLc[k]) * 0.0L * mFffFS; */

    for (m=0; m<2; m++)
      for (n=0; n<2; n++) {
	TSIL_SetParametersST (&bar,m2_sdL[k],m2_sdL[m],m2_sdL[n],Q2);
	TSIL_Evaluate (&bar, s);
	/* For clarity (we just need this one function): */
	sSSS = -bar.S[uxv].value;
	term += 0.25L * PsdL[i][k] * PsdL[k][m] * PsdL[m][n] * PsdL[n][j] * sSSS;
      }
  }
  result *= 4.0L*Cq*Cq;
  result += term*(4.0L*Cq*Cq - 2.0L*CG*Cq);

  /* Next bit is last three lines of eq. (4.8).  These involve summing
     over all 12 quark/squark mass eigenstates. */

  term = 0.0L;

  /* stop */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_gluino,m2_top,m2_stop[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    term += -2.0L * NsdL[i][j] * m_top * m_gluino * vfFfFS
      + NsdL[i][j] * NsdL[r][r] * 0.0L * m_top * (vfFFfS + m2_gluino*vffffS);

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS - 2.0L*NsdL[r][r]*m_top*m_gluino*vFFffS;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sdL[k],m2_stop[r],m2_stop[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * PsdL[i][k] * PsdL[k][j] * PsdL[r][t] * PsdL[t][r] * sSSS;
      }
  }

  /* sbot */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_gluino,m2_bot,m2_sbot[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Take mbot ~ 0 */
    /* term += -2.0L*NsdL[i][j]*m_bot*m_gluino*vfFfFS */
    /*   + NsdL[i][j]*NsdL[r][r]*0.0L*m_bot*(vfFFfS + m2_gluino*vffffS) */
    ;

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*NsdL[r][r]*m_bot*m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sdL[k],m2_sbot[r],m2_sbot[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * PsdL[i][k] * PsdL[k][j] * PsdL[r][t] * PsdL[t][r] * sSSS;
      }
  }

  /* sdL */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_gluino,0.0L,m2_sdL[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*NsdL[i][j]* 0.0L *m_gluino*vfFfFS */
    /*   + NsdL[i][j]*NsdL[r][r]*0.0L* 0.0L *(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*NsdL[r][r]* 0.0L *m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sdL[k],m2_sdL[r],m2_sdL[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * PsdL[i][k] * PsdL[k][j] * PsdL[r][t] * PsdL[t][r] * sSSS;
      }
  }

  /* suR */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_gluino,0.0L,m2_suR[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*NsdL[i][j]* 0.0L *m_gluino*vfFfFS */
    /*   + NsdL[i][j]*NsdL[r][r]*0.0L* 0.0L *(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*NsdL[r][r]* 0.0L *m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sdL[k],m2_suR[r],m2_suR[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * PsdL[i][k] * PsdL[k][j] * PsdL[r][t] * PsdL[t][r] * sSSS;
      }
  }

  /* sdL */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_gluino,0.0L,m2_sdL[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*NsdL[i][j]* 0.0L *m_gluino*vfFfFS */
    /*   + NsdL[i][j]*NsdL[r][r]*0.0L* 0.0L *(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*NsdL[r][r]* 0.0L *m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sdL[k],m2_sdL[r],m2_sdL[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * PsdL[i][k] * PsdL[k][j] * PsdL[r][t] * PsdL[t][r] * sSSS;
      }
  }

  /* sdR */
  for (r=0; r<2; r++) {

    TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_gluino,0.0L,m2_sdR[r],Q2);
    TSIL_Evaluate (&bar, s);
    TSIL_CopyResult (&bar, &gaak); 
    V_FFFFS (&gaak, &gaak, 
	     &vFFFFS, &vFffFS, &vfFfFS, &vfFFfS, &vFFffS, &vffffS);

    /* Tese terms are all zero since the fermion mass is zero */
    /* term = -2.0L*NsdL[i][j]* 0.0L *m_gluino*vfFfFS */
    /*   + NsdL[i][j]*NsdL[r][r]*0.0L* 0.0L *(vfFFfS + m2_gluino*vffffS); */

    if (i == j)
      term += vFFFFS + m2_gluino*vFffFS 
	/* - 2.0L*NsdL[r][r]* 0.0L *m_gluino*vFFffS */
	;

    for (k=0; k<2; k++)
      for (t=0; t<2; t++) {
        TSIL_SetParametersST (&bar,m2_sdL[k],m2_sdR[r],m2_sdR[t],Q2);
        TSIL_Evaluate (&bar, s);
        /* For clarity (just need this one function): */
        sSSS = -bar.S[uxv].value;
	term += 0.25L * PsdL[i][k] * PsdL[k][j] * PsdL[r][t] * PsdL[t][r] * sSSS;
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
   evaluated for s = m2_sdL[i]. */

TSIL_COMPLEX pi21tilde_sdL (int i, TSIL_REAL s) 
{
  int k;
  TSIL_COMPLEX result = 0.0L;
  TSIL_DATA bar;
  TSIL_RESULT gaak;
  TSIL_COMPLEX gFF,gff;
  TSIL_COMPLEX gSSFF,gSSff;

  if (arePandNsdLSet == NO) SetPandNsdL ();

  TSIL_SetParameters (&bar,0.0L,0.0L,m2_gluino,m2_gluino,0.0L,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak); 
  G_FF (&gaak, &gFF, &gff);

  TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_sdL[i],0.1L,m2_gluino,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_sdL[i],0.01L,m2_gluino,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  TSIL_SetParameters (&bar,0.0L,m2_gluino,m2_sdL[i],0.0L,m2_gluino,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  /* NOTE m_q = 0 for q != top */
  result += Cq*CG*(gFF + gSSFF 
		   /* - NsdL[i][i]*0.0L*m_gluino*(gff + gSSff) */
		   );

  TSIL_SetParameters (&bar,0.0L,0.0L,m2_sdL[i],m2_gluino,0.0L,Q2);
  TSIL_Evaluate (&bar, s);
  TSIL_CopyResult (&bar, &gaak);
  Gtilde_SSFF (&gaak, &gSSFF, &gSSff);

  result += (2.0L*Cq - CG)*Cq*(gSSFF
			       /* - NsdL[i][i]*0.0L*m_gluino*gSSff */
			       );

  for (k=0; k<2; k++)
    result += Cq * Cq * PsdL[i][k] * PsdL[k][i] *
      (G_S(m2_sdL[k],Q2) + Gtilde_SSS(m2_sdL[i],m2_sdL[k],Q2));

  result *= g3*g3*g3*g3;

  return result;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Eq. (4.10), tilde version */
/* s arg is useless for same reason as above... */

TSIL_COMPLEX pi22tilde_sdL (int i, TSIL_REAL s) 
{
  int r;
  TSIL_COMPLEX result, term;

  result = Cq*F1tilde(m2_sdL[i], Q2) 
         + CG*F2tilde(m2_sdL[i], Q2)
         + CG*F3tilde(m2_sdL[i], m2_gluino, Q2);

  /* Now do final sum over r. First is r = top. Note that the two F3
     terms are the same for r=top, and all remaining F3 terms are
     identical for m_quark = 0. */
  term = 2.0L*F3tilde(m2_sdL[i], m2_top, Q2) 
       + 2.0L*F3tilde(m2_sdL[i], m2_bot, Q2) 
       + 8.0L*F3tilde(m2_sdL[i], 0.0L, Q2);

  for (r=0; r<2; r++)
    term += F4tilde(m2_sdL[i], m2_stop[r], Q2)
          + F4tilde(m2_sdL[i], m2_sbot[r], Q2)
          + F4tilde(m2_sdL[i], m2_sdL[r], Q2)
          + F4tilde(m2_sdL[i], m2_suR[r], Q2)
          + F4tilde(m2_sdL[i], m2_sdL[r], Q2)
          + F4tilde(m2_sdL[i], m2_sdR[r], Q2);

  /* Final assembly: */
  result += Iq*term;
  result *= g3*g3*g3*g3*Cq;

  return result;
}

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* eq. (4.12); valid for s = m2_sdL[i] */

TSIL_COMPLEX pi20tilde_sdL (int i, TSIL_REAL s) 
{
  if (arePandNsdLSet == NO) SetPandNsdL ();

  return pi20_sdL (i, i, s) 
    + pi1_sdL (i, i, s) * 2.0L * g3 * g3 * Cq * (BpFF(0.0L,m2_gluino,s,Q2) 
    /* - NsdL[i][i] * 0.0L * m_gluino * Bpff(0.0L,m2_gluino,s,Q2) */
						  );
}
