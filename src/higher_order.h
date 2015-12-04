#ifndef HIGHERORDER_H
#define HIGHERORDER_H

typedef struct{
  //Inputs
  double gp, g, g3, ytop, ybot, ytau;
  double vu, vd, mu, Q;
  double m1, m2, m3, atop, abot, atau, b, m2Hu, m2Hd;
  double m2Q[3], m2u[3], m2d[3], m2L[3], m2e[3]; 
  // DGR added:
  double mtop;
  //Outputs
  double mgluino;
  double mstop1, mstop2, msbot1, msbot2, muL, mdL, mcL, msL, muR, mdR, mcR, msR;
} supermodel;

#ifdef __cplusplus
 extern "C" void higherorder(supermodel &smodel);
#else
 void higherorder(supermodel *smodel); 
#endif

#endif
