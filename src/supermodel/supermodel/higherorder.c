/* Interfacing SOFTSUSY with SOFTSUMO program. */

/* Calculates 2-loop gluino pole masses for a range of
   renormalization scales. */

#include "supermodel.h"
#include "supermodel_defs.h"
#include "../../higher_order.h"

void higherorder (supermodel *smodel)
{

  int i;

  /* Disable annoying warning messages */
  TSIL_WarnsOff ();

#ifdef TSIL_SIZE_DOUBLE
  vu = smodel->vu; 
  vd = smodel->vd; 
  gp = smodel->gp; 
  g =  smodel->g; 
  g3 = smodel->g3; 
  ytop = smodel->ytop; 
  ybot = smodel->ybot; 
  ytau = smodel->ytau; 
  M1 = smodel->m1; 
  M2 = smodel->m2;  
  M3 = smodel->m3;  
  atop = smodel->atop; 
  abot = smodel->abot; 
  atau = smodel->atau;
  //
  for(i=0; i<3; i++) {
    m2Q[i] = smodel->m2Q[i];
    m2u[i] = smodel->m2u[i];
    m2d[i] = smodel->m2d[i];
    m2L[i] = smodel->m2L[i];
    m2e[i] = smodel->m2e[i];
  }
  //
  m2Hu = smodel->m2Hu; 
  m2Hd = smodel->m2Hd; 
  mu = smodel->mu; 
  b = smodel->b; 
  //renormalization scale
  Q =  smodel->Q;
  // DGR added:
  M_top = smodel->mtop;
#else
  vu = (long double) smodel->vu; 
  vd = (long double) smodel->vd;
  gp = (long double) smodel->gp; 
  g =  (long double) smodel->g; 
  g3 = (long double) smodel->g3; 
  ytop = (long double) smodel->ytop; 
  ybot = (long double) smodel->ybot; 
  ytau = (long double) smodel->ytau; 
  M1 = (long double) smodel->m1; 
  M2 = (long double) smodel->m2;  
  M3 = (long double) smodel->m3;  
  atop = (long double) smodel->atop; 
  abot = (long double) smodel->abot; 
  atau = (long double) smodel->atau;
  //
  for(i=0; i<3; i++) {
    m2Q[i] = (long double) smodel->m2Q[i];
    m2u[i] = (long double) smodel->m2u[i];
    m2d[i] = (long double) smodel->m2d[i];
    m2L[i] = (long double) smodel->m2L[i];
    m2e[i] = (long double) smodel->m2e[i];
  }
  //
  m2Hu = (long double) smodel->m2Hu; 
  m2Hd = (long double) smodel->m2Hd; 
  mu = (long double) smodel->mu; 
  b = (long double) smodel->b; 
  //renormalization scale  
  Q =  (long double) smodel->Q; 
  // DGR added:
  M_top = (long double) smodel->mtop;
#endif

  tanbeta = vu/vd;

  /// DEBUG: Ben has added for checking
  //  TSIL_REAL m2Ztree, mztree;
  //  TSIL_REAL cos2bet;
  //  TSIL_REAL primer818;
  //  TSIL_REAL primer819;
  // Following uses normalization where VEV is about 175 GeV
  //  vu = vu/sqrt(2.);
  //  vd = vd/sqrt(2.);

  //  printf("Softsusy (vu, vd, tanbeta) = (%Lf, %Lf, %Lf)\n", vu, vd, tanbeta);
  //  SUMO_Update ();
  //  SUMO_Tree_Masses ();
  //  SUMO_Tree_Couplings ();
  //  printf("Supermodel 0-loop (vu, vd) = (%Lf, %Lf)\n", vu, vd);
  //  SUMO_Minimize_Veff (1);
  //  printf("Supermodel 1-loop (vu, vd) = (%Lf, %Lf)\n", vu, vd);
  //  SUMO_Minimize_Veff (2);
  //  printf("Supermodel 2-loop (vu, vd) = (%Lf, %Lf)\n", vu, vd);
  //  m2Ztree = (g*g + gp*gp)*(vu*vu + vd*vd)/2.0;
  //  mztree = sqrt(m2Ztree);
  //  cos2bet = (vd*vd - vu*vu)/(vd*vd + vu*vu);
  //  primer818 = TSIL_CREAL(m2Hu + mu*mu - b*vd/vu - m2Ztree*cos2bet/2.0);
  //  primer819 = TSIL_CREAL(m2Hd + mu*mu - b*vu/vd + m2Ztree*cos2bet/2.0);
  //  printf("Primer eq. (8.1.8) = %Lf\n", primer818);
  //  printf("Primer eq. (8.1.9) = %Lf\n", primer819);
  //  printf("Primer eq. (8.1.8)/mu^2 = %Lf\n", primer818/(mu*mu));
  //  printf("Primer eq. (8.1.9)/mu^2 = %Lf\n", primer819/(mu*mu));
  /// end of DEBUG

  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  /* Minimize 2-loop Veff, to set correct vu,vd,tanbeta: */
  /*  printf("Before minimizing (vu, vd) = (%Lf, %Lf)\n", vu, vd); 
  SUMO_Minimize_Veff (0); 
  printf("After minimizing  (vu, vd) = (%Lf, %Lf)\n", vu, vd); */

  // Top pole mass needs to be set before gluino calculation:
  M2_top = M_top * M_top;

  // DGRDGR Also set lightest Higgs and Goldstones to correct pole
  // masses, in case these have gone negative or are otherwise wonky.
  // Note this has to be done *after* SUMO_Tree_masses () is called!
  // Light Higgs is the pole mass^2 calculated by SOFTSUSY:
  m2_phi0[0] = (smodel->mh0) * (smodel->mh0);

  /* Could also just use this for simplicity: */
  /* m2_phi0[0] = 125.L*125.L; */

  // Exact Goldstone pole masses:
  m2_phi0[2] = m2_phip[0] = 0.0L;
  // DGRDGR end of additions

  int which = 0;  int loops = 2;
  for (which = 0; which < 2; which++) {
    /*    loops = 1;
    SUMO_StopPole (which, loops);
    printf("Loops=%i Which=%i stop mass:     %Lf
    \n",loops,which,SUMO_SQRT(M2_stop[which]));
    loops = 2;
    */
    SUMO_StopPole (which, loops);
    //    printf("Loops=%i Which=%i stop mass:     %Lf \n",loops,which,SUMO_SQRT(M2_stop[which]));

    /*    loops = 1;
    SUMO_SbotPole (which, loops);
    printf("Loops=%i Which=%i sbot mass:     %Lf \n",loops,which,SUMO_SQRT(M2_sbot[which]));
    loops = 2;*/
    SUMO_SbotPole (which, loops);
    //    printf("Loops=%i Which=%i sbot mass:     %Lf \n",loops,which,SUMO_SQRT(M2_sbot[which]));

    /*    loops = 1;
    SUMO_SuLPole (which, loops);
    printf("Loops=%i Which=%i suL mass:     %Lf \n",loops,which,SUMO_SQRT(M2_suL[which]));
    loops = 2;*/
    SUMO_SuLPole (which, loops);
    //    printf("Loops=%i Which=%i suL mass:     %Lf \n",loops,which,SUMO_SQRT(M2_suL[which]));

    /*    loops = 1;
    SUMO_SuRPole (which, loops);
    printf("Loops=%i Which=%i suR mass:     %Lf \n",loops,which,SUMO_SQRT(M2_suR[which]));
    loops = 2; */
    SUMO_SuRPole (which, loops);
    //    printf("Loops=%i Which=%i suR mass:     %Lf \n",loops,which,SUMO_SQRT(M2_suR[which]));

    /*loops = 1;
     SUMO_SdLPole (which, loops);
    printf("Loops=%i Which=%i sdL mass:     %Lf \n",loops,which,SUMO_SQRT(M2_sdL[which]));
    loops = 2; */
    SUMO_SdLPole (which, loops);
    //    printf("Loops=%i Which=%i sdL mass:     %Lf \n",loops,which,SUMO_SQRT(M2_sdL[which]));

    /*    loops = 1;
    SUMO_SdRPole (which, loops);
    printf("Loops=%i Which=%i sdR mass:     %Lf
    \n",loops,which,SUMO_SQRT(M2_sdR[which])); 
    loops = 2;*/
    SUMO_SdRPole (which, loops);
    //    printf("Loops=%i Which=%i sdR mass:     %Lf \n",loops,which,SUMO_SQRT(M2_sdR[which]));
  }

  // DGR Moved this *after* top and squark pole mass calculation
  // (necessary when expanding gluino around pole masses).
  // For debugging purposes, we print
  // SUMO_GluinoPole (0,1);
  // printf("Tree-level gluino mass: %Lf \n",SUMO_SQRT(M2_gluino));
  // SUMO_GluinoPole (1,1);
  // printf("One-loop gluino mass:   %Lf \n",SUMO_SQRT(M2_gluino));
  int eag = smodel->expandAroundGluinoPole;
  SUMO_GluinoPole (2, eag);
  //  printf("Two-loop gluino mass:  %d %Lf \n",eag, SUMO_SQRT(M2_gluino));


  /*  SUMO_h0Pole (loops);
      SUMO_H0Pole (loops);
      SUMO_G0Pole (loops);
      SUMO_A0Pole (loops); */

#ifdef TSIL_SIZE_DOUBLE
  smodel->mgluino = SUMO_SQRT(M2_gluino);
  smodel->mstop1  = SUMO_SQRT(M2_stop[0]);
  smodel->mstop2  = SUMO_SQRT(M2_stop[1]);
  smodel->msbot1  = SUMO_SQRT(M2_sbot[0]);
  smodel->msbot2  = SUMO_SQRT(M2_sbot[1]);
  smodel->muL     = SUMO_SQRT(M2_suL[0]);
  smodel->mcL     = SUMO_SQRT(M2_suL[1]);
  smodel->mdL     = SUMO_SQRT(M2_sdL[0]);
  smodel->msL     = SUMO_SQRT(M2_sdL[1]);
  smodel->muR     = SUMO_SQRT(M2_suR[0]);
  smodel->mcR     = SUMO_SQRT(M2_suR[1]);
  smodel->mdR     = SUMO_SQRT(M2_sdR[0]);
  smodel->msR     = SUMO_SQRT(M2_sdR[1]);
#else
  smodel->mgluino = (double) SUMO_SQRT(M2_gluino);
  smodel->mstop1  = (double) SUMO_SQRT(M2_stop[0]);
  smodel->mstop2  = (double) SUMO_SQRT(M2_stop[1]);
  smodel->msbot1  = (double) SUMO_SQRT(M2_sbot[0]);
  smodel->msbot2  = (double) SUMO_SQRT(M2_sbot[1]);
  smodel->muL     = (double) SUMO_SQRT(M2_suL[0]);
  smodel->mcL     = (double) SUMO_SQRT(M2_suL[1]);
  smodel->mdL     = (double) SUMO_SQRT(M2_sdL[0]);
  smodel->msL     = (double) SUMO_SQRT(M2_sdL[1]);
  smodel->muR     = (double) SUMO_SQRT(M2_suR[0]);
  smodel->mcR     = (double) SUMO_SQRT(M2_suR[1]);
  smodel->mdR     = (double) SUMO_SQRT(M2_sdR[0]);
  smodel->msR     = (double) SUMO_SQRT(M2_sdR[1]);
#endif

}
