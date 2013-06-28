#include <iostream>
#include "nmssmsoftsusy.h"
#include "mycomplex.h"
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "rge.h"
#include "softsusy.h"
#include "nmssmsoftpars.h"
#include "nmssmsusy.h"
#include "softpars.h"
#include "susy.h"
#include "utils.h"
#include "numerics.h"
#include <time.h>

double getRandomNum(double Min, double Max){
   double temp = rand();
   return Min + temp /  RAND_MAX * (Max - Min);
}
int Fillfortest(SoftParsNmssm & r){
   srand(time(NULL));
   double Scale = getRandomNum(0, 300);
   DoubleMatrix yu(3,3), yd(3,3), ye(3,3), mUR(3,3), mDR(3,3), mER(3,3), mQL(3,3), mLL(3,3), ua(3,3), da(3,3), ea(3,3);
   for(int i = 1; i <= 3; i++){
      for(int j = 1; j <= 3; j++){
         yu(i,j) = getRandomNum(0, 1.5);
         yd(i,j) = getRandomNum(0, 1.5);
         ye(i,j) = getRandomNum(0, 1.5);
         mUR(i,j) = getRandomNum(0, 1000 * 1000);
         mDR(i,j) = getRandomNum(0, 1000 * 1000);
         mER(i,j) = getRandomNum(0, 1000 * 1000);
         mQL(i,j) = getRandomNum(0, 1000 * 1000);
         mLL(i,j) = getRandomNum(0, 1000 * 1000);
         mLL(i,j) = getRandomNum(0, 1000 * 1000);
         ua(i,j) = getRandomNum(0, 1000);
         da(i,j) = getRandomNum(0, 1000);
         ea(i,j) = getRandomNum(0, 1000);
      }
   }
   r.setMu(Scale);
   for(int i = 1; i <= 3; i++){
      for(int j = 1; j <= 3; j++){
         r.setYukawaElement(YU, i, j, yu(i,j));
         r.setYukawaElement(YD, i, j, yd(i,j));
         r.setYukawaElement(YE, i, j, ye(i,j));
         r.setSoftMassElement(mUr, i, j, mUR(i,j));
         r.setSoftMassElement(mDr, i, j, mDR(i,j));
         r.setSoftMassElement(mEr, i, j, mER(i,j));
         r.setSoftMassElement(mQl, i, j, mQL(i,j)); 
         r.setSoftMassElement(mLl, i, j, mLL(i,j));
         r.setTrilinearElement(UA, i, j, ua(i,j));
         r.setTrilinearElement(DA, i, j, da(i,j));
         r.setTrilinearElement(EA, i, j, ea(i,j));
      }
   }
   DoubleVector g(3), mG(3);
    for(int j = 1; j <= 3; j++){
         g(j) = getRandomNum(0, 1.5);
         mG(j) = getRandomNum(0, 1000);
    }

    for(int j = 1; j <= 3; j++){
       r.setGaugeCoupling(j, g(j));
       r.setGauginoMass(j, mG(j));
      }
    
   
 double lam = getRandomNum(0, 2.0);
 double kap = getRandomNum(0, 1.9);
 double mupr = getRandomNum(0, 1000);
 double mu = getRandomNum(0, 1000);
 double tb = getRandomNum(1, 15);
 double vev = getRandomNum(10, 1000);
 double svev = getRandomNum(1, 1000);
 double m3sq = getRandomNum(0, 1000 * 1000);
 double mh1sq = getRandomNum(0, 1000 * 1000);
 double mh2sq = getRandomNum(0, 1000 * 1000);
 double mspsq = getRandomNum(0, 1000 * 1000);
 double mssq = getRandomNum(0, 1000 * 1000);
 double al = getRandomNum(0, 1000);
 double ak = getRandomNum(0, 1000);
 double XiS = getRandomNum(0, 1000 * 1000 * 1000);
 double XiF = getRandomNum(0, 1000 * 1000);
 
 cout << "ak = "  << ak << endl;
   r.setLambda(lam);
   //  r.setLambda(0.0);

   r.setKappa(kap);
   r.setMupr(mupr);
   r.setSusyMu(mu);
   r.setTanb(tb);
   r.setHvev(vev);
   r.setSvev(svev);
   r.setXiF(XiF);
   

   r.setM3Squared(m3sq);
   r.setMh1Squared(mh1sq);
   r.setMh2Squared(mh2sq);
   r.setMspSquared(mspsq);
   r.setMsSquared(mssq);
   r.setTrialambda(al);
   r.setTriakappa(ak);
   r.setXiS(XiS);
  
   
   cout << "r.displayTriakappa() = " << r.displayTriakappa() << endl;
  
   
   
  

   return 0;
}

int main() {

  /// Sets format of output: 12 decimal places
  outputCharacteristics(20);
  double eps = 1.0e-4;
  MIXING = 1;  PRINTOUT = 0;  
  bool speak = true;
  
  /// Parameters used: CMSSM parameters
  double m12 = 500., a0 = 0., mGutGuess = 2.0e16, tanb = 10.0, m0 = 125.;
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 10; ///< number of scan points
  //PA: new nmssm parameters to be fed into lowOrg
  double lambda = 0.0001, kappa = 0.00001, s = 1e3, xiF = 100.0, mupr = 10.0;
  xiF = 0; mupr = 0; lambda = -2e-8; kappa =0; s =1;
  
  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters
 
  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMass(mBottom, mbmb);
  oneset.toMz();      ///< Runs SM fermion masses to MZ
  
  NmssmSoftsusy n;
  cout << "n = "  << n << endl;
  Fillfortest(n);
  cout << "n = "  << n << endl;
  n.calcDrBarPars();
  n.physical(0);
  
  DoubleVector pars(3); 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  DoubleVector nmpars(5);
  nmpars(1) = lambda; nmpars(2) = kappa; nmpars(3) = s; 
  nmpars(4) = xiF; nmpars(5) = mupr;
  bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
  n.lowOrg(NmssmMsugraBcs, mGutGuess, pars, nmpars, sgnMu, tanb, oneset, uni);
  
  n.printall();
   
  return 1;
}
