/* Supermodel data declarations */

#ifndef _SUPERMODEL_DEFS_H_
#define _SUPERMODEL_DEFS_H_

/* -------------------------------------------------------------------- */
/* Main Supermodel Data */

/* MSBar EM coupling at MZ, 5 flavors */
SUMO_REAL alpha_MZ_MS;

/* Fine structure constant */
SUMO_REAL alpha0;

/* MSBar strong coupling at MZ, 5 flavors */
SUMO_REAL alphaS_MZ_MS;

/* Fermi constant */
SUMO_REAL G_Fermi;

/* MSBar bottom quark mass at Q=mb */
SUMO_REAL mb_mb_MS;

/* Physical masses */
SUMO_REAL M_top_phys, M_bot_phys, M_tau_phys;
SUMO_REAL M_W_phys, M_Z_phys;
SUMO_REAL M2_top_phys, M2_bot_phys, M2_tau_phys;
SUMO_REAL M2_W_phys, M2_Z_phys;

TSIL_REAL Q;          /* Renormalization scale */

/* Lagrangian Parameters: */
TSIL_REAL gp, g, g3;           /* Gauge couplings [U(1)_Y, SU(2),
				  and SU(3), resp.] */
TSIL_REAL ytop, ybot, ytau;    /* 3rd generation Yukawas (we assume
				  others are essentially zero) */
TSIL_COMPLEX mu;               /* mu parameter, including phase */

TSIL_COMPLEX atop, abot, atau; /* 3rd generation trilinear
				  couplings.  We assume the 1st and
				  2nd generation couplings are
				  zero. */
TSIL_COMPLEX b;                /* This can and will be made real by
				  a rotation at any given RG scale,
				  but at intermediate steps it might
				  be complex. */

TSIL_COMPLEX M1;     /* U(1)_Y gaugino mass (Bino) */
TSIL_COMPLEX M2;     /* SU(2)_L gaugino mass (Wino) */
TSIL_COMPLEX M3;     /* SU(3) gaugino mass (gluino) */

/* Note that one might assume 1st/2nd generation universality, in
   which case the following mass parameters are the same for these
   generations, hence 10 independent parameters rather than 15. */

TSIL_REAL m2Q[3]; /* LH squark mass^2 matrix (diagonal els only) */
TSIL_REAL m2L[3]; /* LH slepton mass^2 matrix (diagonal els only) */
TSIL_REAL m2u[3]; /* RH up-type squark mass^2 matrix (diagonal els
		     only) */
TSIL_REAL m2d[3]; /* RH down-type squark mass^2 matrix (diagonal els
		     only) */
TSIL_REAL m2e[3]; /* RH slepton mass^2 matrix (diagonal els only) */

/* Note these can be negative */
TSIL_REAL m2Hu;   /* Up-type Higgs mass squared */
TSIL_REAL m2Hd;   /* Down-type Higgs mass squared */

/* Higgs VEVs */
TSIL_REAL vu, vd;

/* Input value of tanbetaMZ, not to be changed during iterations. */
TSIL_REAL tanbetaMZ;
TSIL_REAL sbetaMZ, cbetaMZ;

/* Convenient combinations of the above.  These are updated by the
   routine SUMO_Update (): */
int          is_updated;
TSIL_REAL    Q2;
TSIL_REAL    g2, gp2, g2plusgp2, e, e2;
TSIL_REAL    vu2, vd2, v2, vuvd;
TSIL_REAL    ytop2, ybot2, ytau2;
TSIL_COMPLEX muc, atopc, abotc, atauc; 
TSIL_REAL    mu2;
TSIL_REAL    tanbeta, sinbeta, cosbeta;

/* Numerical values of beta functions of basic inputs. These are
   computed by SUMO_Betas (): */
TSIL_REAL    beta_gp, beta_g, beta_g3;
TSIL_REAL    beta_ytop, beta_ybot, beta_ytau;
TSIL_COMPLEX beta_mu;
TSIL_COMPLEX beta_M1, beta_M2, beta_M3;
TSIL_COMPLEX beta_atop, beta_abot, beta_atau;
TSIL_COMPLEX beta_b;
TSIL_REAL    beta_vu, beta_vd;
TSIL_REAL    beta_m2Hu, beta_m2Hd;
TSIL_REAL    beta_m2Q12, beta_m2Q3;
TSIL_REAL    beta_m2L12, beta_m2L3;
TSIL_REAL    beta_m2u12, beta_m2u3;
TSIL_REAL    beta_m2d12, beta_m2d3;
TSIL_REAL    beta_m2e12, beta_m2e3;

/* Tree-level squared mass eigenvalues. Computed in
   SUMO_Tree_Masses(): */
int are_tree_masses_updated;
TSIL_REAL m2_W, m2_Z;
TSIL_REAL m2_phi0[4];
TSIL_REAL m2_phip[2];
TSIL_REAL m2_top, m2_bot, m2_tau;
TSIL_REAL m2_gluino;
TSIL_REAL m2_Neut[4], m2_Char[2];

/* Arrays hold the 1,2 mass eigenvalues for 3rd family: */
TSIL_REAL m2_stop[2], m2_sbot[2], m2_stau[2];

/* For these, array elements correspond to families: */
TSIL_REAL m2_suL[2], m2_suR[2];
TSIL_REAL m2_sdL[2], m2_sdR[2];
TSIL_REAL m2_seL[2], m2_seR[2];
TSIL_REAL m2_snu[3];

/* Tree-level fermion mass eigenvalues.  Computed in
   SUMO_Tree_Masses (): */
TSIL_REAL m_top, m_bot, m_tau;
TSIL_REAL m_gluino;
TSIL_REAL m_Neut[4], m_Char[2]; /* Neutralinos, Charginos */

/* (Complex) pole squared mass eigenvalues. */
TSIL_COMPLEX CM2_W, CM2_Z;
TSIL_COMPLEX CM2_phi0[4];
TSIL_COMPLEX CM2_phip[2];
TSIL_COMPLEX CM2_top, CM2_bot, CM2_tau;
TSIL_COMPLEX CM2_gluino;
TSIL_COMPLEX CM2_Neut[4], CM2_Char[2]; 
TSIL_COMPLEX CM2_stop[2], CM2_sbot[2], CM2_stau[2];
TSIL_COMPLEX CM2_suL[2], CM2_suR[2];
TSIL_COMPLEX CM2_sdL[2], CM2_sdR[2];
TSIL_COMPLEX CM2_seL[2], CM2_seR[2];
TSIL_COMPLEX CM2_snu[3];

/* Pole squared masses (real parts). */
/* Entire particle content of the MSSM added, though some are
   useless at the moment; these must be in the same order as the
   list of enums in sumo_pmdata.h. */

TSIL_REAL M2_dn, M2_up, M2_str, M2_cha, M2_bot, M2_top;
TSIL_REAL M2_el, M2_nuel;
TSIL_REAL M2_mu, M2_numu;
TSIL_REAL M2_tau, M2_nutau;
TSIL_REAL M2_gluon, M2_photon;
TSIL_REAL M2_Z, M2_W;
TSIL_REAL M2_phi0[4];
TSIL_REAL M2_phip[2];
TSIL_REAL M2_grav, M2_gravitino;
TSIL_REAL M2_suL[2], M2_sdL[2];
TSIL_REAL M2_suR[2], M2_sdR[2];
TSIL_REAL M2_seL[2], M2_seR[2];
TSIL_REAL M2_snu[3];
TSIL_REAL M2_stop[2], M2_sbot[2], M2_stau[2];
TSIL_REAL M2_gluino;
TSIL_REAL M2_Neut[4], M2_Char[2];

/* Pole masses of fermions (square roots of real parts above). */
TSIL_REAL M_top, M_bot, M_tau;
TSIL_REAL M_gluino;
TSIL_REAL M_Neut[4], M_Char[2]; 

/* Tree-level mixing matrices (notations as in hep-ph/0206136
   section 2). Computed in SUMO_Tree_Couplings(): */
TSIL_COMPLEX Nmix[4][4], Nmixc[4][4];      /* Neutralinos */
TSIL_COMPLEX Umix[2][2], Vmix[2][2];       /* Charginos */
TSIL_COMPLEX Umixc[2][2], Vmixc[2][2]; 
TSIL_COMPLEX ku[4], kd[4], kuc[4], kdc[4]; /* Neutral Higgs
					      scalars */
TSIL_REAL    kup[2], kdp[2];               /* Charged Higgs
					      scalars */

/* In SUMO_Tree_Masses(): */
TSIL_REAL    calpha, salpha, cbeta0, sbeta0, cbetapm, sbetapm;
TSIL_REAL    c2alpha, s2alpha, c2beta0, s2beta0, c2betapm, s2betapm;
TSIL_COMPLEX Rstop[2], Lstop[2], Rstopc[2], Lstopc[2];
TSIL_COMPLEX Rsbot[2], Lsbot[2], Rsbotc[2], Lsbotc[2];
TSIL_COMPLEX Rstau[2], Lstau[2], Rstauc[2], Lstauc[2];

/* Field-independent vacuum energy, only used for checking RG
   invariance. */
TSIL_COMPLEX LambdaVacuum; 

int are_tree_couplings_updated;

/* Higgs^4 couplings (hep-ph/0404055 eqs. (2.43)-(2.45)). Computed in
   SUMO_Tree_Couplings() */
TSIL_REAL    lambda0000[4][4][4][4];
TSIL_COMPLEX lambda00pm[4][4][2][2];
TSIL_COMPLEX lambdappmm[2][2][2][2];

/* Higgs^3 couplings (hep-ph/0404055 eq. (2.46)). Computed in
   SUMO_Tree_Couplings() */
TSIL_REAL    lambda000[4][4][4];
TSIL_COMPLEX lambda0pm[4][2][2];

/* Neutral Higgs-sfermion-sfermion couplings (hep-ph/0405022,
   eqs. (2.51)-(2.54)). */
TSIL_REAL    lambda0suLsuLc[4];
TSIL_REAL    lambda0suRsuRc[4];
TSIL_REAL    lambda0sdLsdLc[4];
TSIL_REAL    lambda0sdRsdRc[4];
TSIL_REAL    lambda0snusnuc[4];
TSIL_REAL    lambda0seLseLc[4];
TSIL_REAL    lambda0seRseRc[4];
TSIL_COMPLEX lambda0stopstopc[4][2][2];
TSIL_COMPLEX lambda0sbotsbotc[4][2][2];
TSIL_COMPLEX lambda0staustauc[4][2][2];
/* DGR this is the same as lambda0snusnuc[] above, but... */
TSIL_COMPLEX lambda0snutausnutauc[4];

/* DGR - added (code in SUMO_Tree_SSSS_Couplings) */
/* Neutral Higgs-neutral Higgs-sfermion-sfermion 4-point couplings
   hep-ph/0405022, eqs. (2.47)-(2.50). */
TSIL_REAL    lambda00suLsuLc[4][4];
TSIL_REAL    lambda00suRsuRc[4][4];
TSIL_REAL    lambda00sdLsdLc[4][4];
TSIL_REAL    lambda00sdRsdRc[4][4];
TSIL_REAL    lambda00snusnuc[4][4];
TSIL_REAL    lambda00seLseLc[4][4];
TSIL_REAL    lambda00seRseRc[4][4];
TSIL_COMPLEX lambda00stopstopc[4][4][2][2];
TSIL_COMPLEX lambda00sbotsbotc[4][4][2][2]; 
TSIL_COMPLEX lambda00staustauc[4][4][2][2];
/* DGR this is the same as lambda00snusnuc[] above, but... */
TSIL_COMPLEX lambda00snutausnutauc[4][4];

/* Charged Higgs-sfermion-sfermion couplings (hep-ph/0405022,
   eqs. (2.60)-2.62)) */
TSIL_REAL    lambdapsdLsuLc[2];
TSIL_REAL    lambdapseLsnuc[2];
TSIL_COMPLEX lambdapsbotstopc[2][2][2];
TSIL_COMPLEX lambdapstausnutauc[2][2];

/* DGR - Added charged Higgs-neutral Higgs-sfermion-antisfermion
     couplings, hep-ph/0405022 eqs. (2.63)-(2.65). */
/* (Code computing these added to SUMO_Tree_SSSS_Couplings.) */
TSIL_COMPLEX lambda0psbotstopc[4][2][2][2];
TSIL_COMPLEX lambda0pstausnutauc[4][2][2];
TSIL_COMPLEX lambda0psdLsuLc[4][2];
TSIL_COMPLEX lambda0pseLsnuc[4][2];
/* The last are equal to the ones above, but perhaps having a
   separate copy will keep things clearer later on... Could perhaps
   also include compex conjugates, as these are occasionally
   needed. */

/* DGR - Added couplings of pairs of charged Higgses to sfermions
   hep-ph/0405022 eqs. (2.55)-(2.59). */
/* (Code computing these added to SUMO_Tree_SSSS_Couplings.) */
TSIL_COMPLEX lambdapmstopstopc[2][2][2][2];
TSIL_COMPLEX lambdapmsbotsbotc[2][2][2][2];
TSIL_COMPLEX lambdapmstaustauc[2][2][2][2];
TSIL_COMPLEX lambdapmsnutausnutauc[2][2];
TSIL_COMPLEX lambdapmsuLsuLc[2][2];
TSIL_COMPLEX lambdapmsuRsuRc[2][2];
TSIL_COMPLEX lambdapmsdLsdLc[2][2];
TSIL_COMPLEX lambdapmsdRsdRc[2][2];
TSIL_COMPLEX lambdapmsnusnuc[2][2];
TSIL_COMPLEX lambdapmseLseLc[2][2];
TSIL_COMPLEX lambdapmseRseRc[2][2];
/* (Again some copies at the end here...) */

/* DGR - Added (code in SUMO_Tree_FFS_Couplings) */
/* Couplings of 3rd-generation fermions to neutral and
   charged Higgs scalars hep-ph/0405022 eqs. (2.1)-(2.6).  */
TSIL_COMPLEX Yttcphi[4], conYttcphi[4];
TSIL_COMPLEX Ybbcphi[4], conYbbcphi[4];
TSIL_COMPLEX Ytautaucphi[4], conYtautaucphi[4];
TSIL_COMPLEX Ytcbphip[2], conYtcbphip[2];
TSIL_COMPLEX Ybctphim[2], conYbctphim[2];
TSIL_COMPLEX Ytaucnutauphim[2], conYtaucnutauphim[2];

/* DGR - a subset only of these couplings added so far. */
/* Sfermion-antisfermion-sfermion-antisfermion couplings,
   hep-ph/0405022 eq. (2.67) et seq. */
TSIL_COMPLEX lambdastopstopcstopstopc[2][2][2][2];
TSIL_COMPLEX lambdasbotsbotcsbotsbotc[2][2][2][2];
TSIL_COMPLEX lambdastaustaucstaustauc[2][2][2][2];
TSIL_COMPLEX lambdasnutausnutaucsnutausnutauc;
TSIL_COMPLEX lambdastopstopcsbotsbotc[2][2][2][2];
TSIL_COMPLEX lambdastopstopcstaustauc[2][2][2][2];
TSIL_COMPLEX lambdastopstopcsnutausnutauc[2][2];
TSIL_COMPLEX lambdasbotsbotcstaustauc[2][2][2][2];
TSIL_COMPLEX lambdasbotsbotcsnutausnutauc[2][2];
TSIL_COMPLEX lambdastaustaucsnutausnutauc[2][2];
TSIL_COMPLEX lambdastopsbotcsbotstopc[2][2][2][2];
TSIL_COMPLEX lambdastausnutaucsnutaustauc[2][2];

/* DGR added 3/12/13 -- last one? */
/*   TSIL_COMPLEX lambdastopsbotcstausnutauc[2][2][2]; */

/* These are related to the above, but perhaps we want separate
   copies? If so there are others we could also include... */

TSIL_COMPLEX lambdasbotstopcstopsbotc[2][2][2][2];
TSIL_COMPLEX lambdasnutaustaucstausnutauc[2][2];
/*   TSIL_COMPLEX staustaucsbtsbotc[2][2][2][2]; */

/* DGR -- end of added items. */

/* Couplings of sfermions to neutralinos and fermions hep-ph/0405022
   eqs. (2.7)-(2.19). */
TSIL_COMPLEX conYuNsuLc[4], YuNsuLc[4];
TSIL_COMPLEX conYucNsuR[4], YucNsuR[4];
TSIL_COMPLEX conYdNsdLc[4], YdNsdLc[4];
TSIL_COMPLEX conYdcNsdR[4], YdcNsdR[4];
TSIL_COMPLEX conYnuNsnuc[4], YnuNsnuc[4];
TSIL_COMPLEX conYeNseLc[4], YeNseLc[4];
TSIL_COMPLEX conYecNseR[4], YecNseR[4];

TSIL_COMPLEX conYtNstopc[4][2], YtNstopc[4][2];
TSIL_COMPLEX conYtcNstop[4][2], YtcNstop[4][2];
TSIL_COMPLEX conYbNsbotc[4][2], YbNsbotc[4][2];
TSIL_COMPLEX conYbcNsbot[4][2], YbcNsbot[4][2];
TSIL_COMPLEX conYtauNstauc[4][2], YtauNstauc[4][2];
TSIL_COMPLEX conYtaucNstau[4][2], YtaucNstau[4][2];

/* Couplings of sfermions to charginos and fermions (hep-ph/0405022
   eqs. (2.20)-(2.27)). */
TSIL_COMPLEX YdCsuLc[2], YeCsnuc[2], YtauCsnutauc[2];
TSIL_COMPLEX YuCsdLc[2], YnuCseLc[2];
TSIL_COMPLEX conYdCsuLc[2], conYeCsnuc[2], conYtauCsnutauc[2];
TSIL_COMPLEX conYuCsdLc[2], conYnuCseLc[2];
TSIL_COMPLEX conYbCstopc[2][2], YbCstopc[2][2];
TSIL_COMPLEX conYbcCstop[2][2], YbcCstop[2][2];
TSIL_COMPLEX conYtCsbotc[2][2], YtCsbotc[2][2];
TSIL_COMPLEX conYtcCsbot[2][2], YtcCsbot[2][2];
TSIL_COMPLEX conYnutauCstauc[2][2], YnutauCstauc[2][2];
TSIL_COMPLEX conYtaucCsnutau[2], YtaucCsnutau[2];

/* Couplings of Higgs scalars to charginos, neutralinos
   (hep-ph/0405022 eqs. (2.28)-(2.31)). */
TSIL_COMPLEX YCC0[2][2][4], conYCC0[2][2][4];
TSIL_COMPLEX YNN0[4][4][4], conYNN0[4][4][4];
TSIL_COMPLEX YCNp[2][4][2], conYCNp[2][4][2];
TSIL_COMPLEX YCNm[2][4][2], conYCNm[2][4][2];

/* Couplings of W,Z to charginos, neutralinos (hep-ph/0206136
   eqs. (3.69)-(3.73); see also Haber-Kane). */
TSIL_COMPLEX OL[4][2], OLc[4][2], OR[4][2], ORc[4][2];
TSIL_COMPLEX OpL[2][2], OpR[2][2];
TSIL_COMPLEX OppL[4][4];

/* Couplings of W,Z to Higgs scalars (hep-ph/0405022
   eqs. (2.32)-(2.42)). */
TSIL_REAL gWW0[4], gZZ0[4];
TSIL_REAL gWZp[2], gWgammap[2];
TSIL_REAL gWW00[4], gZZ00[4], gWWpm[2], gZZpm[2];
TSIL_REAL gZ00[4][4];
TSIL_COMPLEX gZpm[2];
TSIL_COMPLEX gW0p[4][2];

/* Couplings of W,Z to sfermions */
TSIL_COMPLEX gZsuLsuLc, gZsuRsuRc, gZsdLsdLc, gZsdRsdRc;
TSIL_COMPLEX gZseLseLc, gZseRseRc, gZsnusnuc;
TSIL_COMPLEX gZstopstopc[2][2];
TSIL_COMPLEX gZsbotsbotc[2][2];
TSIL_COMPLEX gZstaustauc[2][2];
TSIL_COMPLEX gZZsuLsuLc, gZZsuRsuRc, gZZsdLsdLc, gZZsdRsdRc;
TSIL_COMPLEX gZZseLseLc, gZZseRseRc, gZZsnusnuc;
TSIL_COMPLEX gZZstopstopc[2][2];
TSIL_COMPLEX gZZsbotsbotc[2][2];
TSIL_COMPLEX gZZstaustauc[2][2];

int interp = YES;

#endif /* supermodel_defs.h */
