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
  fclose (stderr);

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
#endif 


  /* Transfering information */
  /*
  vu = 172.L;
  vd = 17.2L;
  gp = 0.36L;
  g = 0.65L;
  g3 = 1.06L;
  ytop = 0.90L;
  ybot = 0.13L;
  ytau = 0.10L;
  M1 = 150.L;
  M2 = 280.L;
  M3 = 800.L;
  atop = -600.L;
  abot = -150.L;
  atau = -40.L;
  m2Q[1] = m2Q[0] = 780.L*780.L;
  m2u[1] = m2u[0] = 740.L*740.L;
  m2d[1] = m2d[0] = 735.L*735.L;
  m2L[1] = m2L[0] = 280.L*280.L;
  m2e[1] = m2e[0] = 200.L*200.L;
  m2Q[2] = 700.L*700.L;
  m2u[2] = 580.L*580.L;
  m2d[2] = 725.L*725.L;
  m2L[2] = 270.L*270.L;
  m2e[2] = 195.L*195.L;
  m2Hu = -500.L*500.L;
  m2Hd = 270.L*270.L;
  mu = 504.1811202L;
  b = 33937.10367865087L;
  */

  //Check this
  /* DGR can safely be ignored */
  /* LambdaVacuum = 0.L; */

  tanbeta = vu/vd;

  /// DEBUG: Ben has added for checking
  TSIL_REAL m2Ztree, mztree;
  TSIL_REAL cos2bet;
  TSIL_REAL primer818;
  TSIL_REAL primer819;
  // Following uses normalization where VEV is about 175 GeV
  vu = vu/sqrt(2.);
  vd = vd/sqrt(2.);

  m2Ztree = (g*g + gp*gp)*(vu*vu + vd*vd)/2.0;
  mztree = sqrt(m2Ztree);
  cos2bet = (vd*vd - vu*vu)/(vd*vd + vu*vu);
  primer818 = TSIL_CREAL(m2Hu + mu*mu - b*vd/vu - m2Ztree*cos2bet/2.0);
  primer819 = TSIL_CREAL(m2Hd + mu*mu - b*vu/vd + m2Ztree*cos2bet/2.0);
  printf("Primer eq. (8.1.8) = %Lf\n", primer818);
  printf("Primer eq. (8.1.9) = %Lf\n", primer819);
  printf("Primer eq. (8.1.8)/mu^2 = %Lf\n", primer818/(mu*mu));
  printf("Primer eq. (8.1.9)/mu^2 = %Lf\n", primer819/(mu*mu));
  /// end of DEBUG

  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  /* Minimize 2-loop Veff, to set correct vu,vd,tanbeta: */
  /*  printf("Before minimizing (vu, vd) = (%Lf, %Lf)\n", vu, vd); 
  SUMO_Minimize_Veff (0); 
  printf("After minimizing  (vu, vd) = (%Lf, %Lf)\n", vu, vd); */

  // For debugging purposes, we print
  SUMO_GluinoPole (0);
  printf("Tree-level gluino mass: %Lf \n",SUMO_SQRT(M2_gluino));
  SUMO_GluinoPole (1);
  printf("One-loop gluino mass:   %Lf \n",SUMO_SQRT(M2_gluino));
  SUMO_GluinoPole (2);
  printf("Two-loop gluino mass:   %Lf \n", SUMO_SQRT(M2_gluino));

  int which = 1;
  int loops = 2;
  printf("Loops=%i \n",loops);
  SUMO_StopPole (which, loops);
  SUMO_SbotPole (which, loops);
  SUMO_SuLPole (which, loops);
  SUMO_SuRPole (which, loops);
  SUMO_SdLPole (which, loops);
  SUMO_SdRPole (which, loops);
  
  SUMO_h0Pole (loops);
  SUMO_H0Pole (loops);
  SUMO_G0Pole (loops);
  SUMO_A0Pole (loops);

#ifdef TSIL_SIZE_DOUBLE
  smodel->mgluino = SUMO_SQRT(M2_gluino);
#else
  smodel->mgluino = (double) SUMO_SQRT(M2_gluino);
#endif

}
