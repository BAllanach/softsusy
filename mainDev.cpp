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
   return Min + temp /  RAND_MAX * Max;
}
int Fillfortest(SoftParsNmssm & r){
   srand(time(NULL));
   cout << "rand() = " << rand() << endl;
   cout << "rand() = " << rand() << endl;
   cout << " RAND_MAX = "  <<  RAND_MAX << endl;
   double Scale = getRandomNum(0, 300);
   cout << "Scale = "  << Scale << endl;
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
         r.setTrilinearElement(DA, i, j, ua(i,j));
         r.setTrilinearElement(EA, i, j, ua(i,j));
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
   double eps = 1.000000e-20;
   MIXING = 1;
   bool speak = true;
 
   NmssmSoftsusy n;
   cout << "n = "  << n << endl;
   Fillfortest(n);
cout << "n = "  << n << endl;
     n.calcDrBarPars();

     n.printall();

     
     
    


    
   
   
   return 1;
}
