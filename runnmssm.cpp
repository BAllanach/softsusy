#include <iostream>
#include "mycomplex.h"
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "rge.h"
#include "nmssmsusy.h"
#include "utils.h"
#include "numerics.h"
#include "nmssmsoftpars.h"
int main() {
 
/// Sets format of output: 6 decimal places
  outputCharacteristics(12);
  cerr << "nmssm SOFTSUSY" 
       << " test program to run  and nmssmsusy and softparsnmssm parameters only, Peter Athron";
  
 
  //Create an empty object
  nMssmSusy n; 
  //show that it's empty
  
  cout << n.display() << endl;
  // cout << "mu_s" << n.displayMu_s() << endl;
  //Fill it with values at the Highscale, Mx

  n.setYukawaElement(YU, 1, 1, 3.2204387e-02);
  // n.setYukawaElement(YU, 1, 2, 5.2098121e-03);
  // n.setYukawaElement(YU, 1, 3, 8.1872365e-03);
  // n.setYukawaElement(YU, 2, 1, 2.4589781e-03); 
  n.setYukawaElement(YU, 2, 2, 1.4639152e-02);
  // n.setYukawaElement(YU, 2, 3, 1.9834821e-01);
  // n.setYukawaElement(YU, 3, 1, 2.3798168e-01);
  // n.setYukawaElement(YU, 3, 2, 9.2832571e-02); 
  n.setYukawaElement(YU, 3, 3, 2.4587251e-01);
  


  n.setYukawaElement(YD, 1, 1, 8.476118e-02);
  // n.setYukawaElement(YD, 1, 2, 2.345788e-05);
  // n.setYukawaElement(YD, 1, 3, 6.321345e-06);
  // n.setYukawaElement(YD, 2, 1, 2.459546e-05); 
  n.setYukawaElement(YD, 2, 2, 3.495862e-03);
  // n.setYukawaElement(YD, 2, 3, 5.302871e-05);
  // n.setYukawaElement(YD, 3, 1, 2.222134e-01);
  // n.setYukawaElement(YD, 3, 2, 5.321345e-05); 
  n.setYukawaElement(YD, 3, 3, 9.342113e-02);

  n.setYukawaElement(YE, 1, 1, 1.222343e-02);
  // n.setYukawaElement(YE, 1, 2, 5.297271e-05);
  // n.setYukawaElement(YE, 1, 3, 8.313451e-01);
  // n.setYukawaElement(YE, 2, 1, 2.111223e-05); 
  n.setYukawaElement(YE, 2, 2, 6.775421e-03);
  // n.setYukawaElement(YE, 2, 3, 3.330733e-05);
  // n.setYukawaElement(YE, 3, 1, 4.687256e-02);
  // n.setYukawaElement(YE, 3, 2, 2.434562e-05); 
  n.setYukawaElement(YE, 3, 3, 1.495761e-01);

  double mx = 1.860292e+16;
  double eps = 1.000000e-05;
    eps = 1.000000e-9;
    MIXING = 0;
   n.setGaugeCoupling(1, 7.134896e-01);
    n.setGaugeCoupling(2, 2.324567e-01);
   n.setGaugeCoupling(3, 1.869472e-01);

    n.setLambda(0.5311);
    //   n.setKappa(0.3486721);
      // n.setKap(0.3486721);
      // n.setMu_s(1);
    // // n.setSusyMuPrime(1);
    // n.setXiF(200);
     n.setZeta(200);
     // n.setSvev(2000);
     //         n.setSusyMu(100);
     n.setTanb(20);
   n.setHvev(249);
  

 double msusy = 8.797252e+02;
 // //cout after filling
 // cout << "n initiallisd to " << n.display() << endl;
 // //run to msusy and output
 // n.run(mx, msusy, eps);
 // cout << "n at msusy = " << n.display() << endl;
 // //run from msusy to MZ and output 
 // n.run(msusy, MZ, eps);
 // cout <<"n at mz = " << n.display() << endl;

 // n.runto(mx, eps);
 // cout << "should be starting values = " << n.display() << endl;


 SoftParsNmssm sn(n);
  cout << "created object from nmssmsusy one should have that data =" << sn.display() << endl;
 // cout << "sn at mx b4 bcs = " << sn.display() << endl;
 //sn.standardsemiSugra(500, 200, 300, 100, 200, 150);
 // cout << "sn at mx after bc's = " << sn.display() << endl;
 // cout << "alambda = "  << sn.displayTrialambda() << endl;

  sn.setM3Squared(1e05);
    sn.setMh1Squared(1e09);
    sn.setMh2Squared(1e08);
  sn.setMspSquared(98);
  sn.setMsSquared(1e10);
     sn.setTrialambda(8.27252e05);
     sn.setTriakappa(2.498275e05);
  // sn.setHLam(8.27252e-02);
  //sn.setHKap(2.498275e03);

     sn.setMu_s(1.487261);
     //sn.setSusyMuPrime(1.487261);
   // sn.setXiS(2.4981);
    sn.setZeta_s(2.4981);
  
 


 sn.setSoftMassElement(mUr, 1, 1, 1.58737e02);
  // sn.setSoftMassElement(mUr, 1, 2, 3.87761e03);
  // sn.setSoftMassElement(mUr, 1, 3, 2.97237e01);
  //sn.setSoftMassElement(mUr, 2, 1, 5.98219e03);
  sn.setSoftMassElement(mUr, 2, 2, 6.09828e01);
  // sn.setSoftMassElement(mUr, 2, 3, 4.98787e02);
  // sn.setSoftMassElement(mUr, 3, 1, 2.24445e04);
  // sn.setSoftMassElement(mUr, 3, 2, 1.11132e01);
  sn.setSoftMassElement(mUr, 3, 3, 5.29287e03);

  sn.setSoftMassElement(mDr, 1, 1, 3.56987e01);
  // sn.setSoftMassElement(mDr, 1, 2, 7.49821e02);
  // sn.setSoftMassElement(mDr, 1, 3, 2.96837e02);
  // sn.setSoftMassElement(mDr, 2, 1, 4.69713e03);
  sn.setSoftMassElement(mDr, 2, 2, 1.28742e04);
  // sn.setSoftMassElement(mDr, 2, 3, 6.97093e01);
  // sn.setSoftMassElement(mDr, 3, 1, 4.87261e02);
  // sn.setSoftMassElement(mDr, 3, 2, 6.98271e02);
  sn.setSoftMassElement(mDr, 3, 3, 2.45987e04);

  sn.setSoftMassElement(mQl, 1, 1, 3.58727e01);
  // sn.setSoftMassElement(mQl, 1, 2, 2.85978e04);
  // sn.setSoftMassElement(mQl, 1, 3, 5.28911e03);
  // sn.setSoftMassElement(mQl, 2, 1, 4.98265e01);
  sn.setSoftMassElement(mQl, 2, 2, 3.79815e01);
  // sn.setSoftMassElement(mQl, 2, 3, 9.76262e04);
  // sn.setSoftMassElement(mQl, 3, 1, 4.29872e02);
  // sn.setSoftMassElement(mQl, 3, 2, 6.19832e03);
  sn.setSoftMassElement(mQl, 3, 3, 3.87763e01);

  sn.setSoftMassElement(mLl, 1, 1, 2.58267e02);
  // sn.setSoftMassElement(mLl, 1, 2, 1.54862e03);
  // sn.setSoftMassElement(mLl, 1, 3, 9.99921e01);
  // sn.setSoftMassElement(mLl, 2, 1, 3.54567e04);
  sn.setSoftMassElement(mLl, 2, 2, 2.48561e01);
  // sn.setSoftMassElement(mLl, 2, 3, 6.22214e02);
  // sn.setSoftMassElement(mLl, 3, 1, 7.21461e02);
  // sn.setSoftMassElement(mLl, 3, 2, 2.89751e01);
  sn.setSoftMassElement(mLl, 3, 3, 5.91919e02);

  sn.setSoftMassElement(mEr, 1, 1, 2.45861e01);
  // sn.setSoftMassElement(mEr, 1, 2, 6.19813e03);
  // sn.setSoftMassElement(mEr, 1, 3, 2.58261e02);
  // sn.setSoftMassElement(mEr, 2, 1, 5.79812e01);
  sn.setSoftMassElement(mEr, 2, 2, 2.98741e03);
  // sn.setSoftMassElement(mEr, 2, 3, 6.29821e04);
  // sn.setSoftMassElement(mEr, 3, 1, 9.86121e01);
  // sn.setSoftMassElement(mEr, 3, 2, 3.58761e01);
  sn.setSoftMassElement(mEr, 3, 3, 1.11234e03);

  sn.setTrilinearElement(UA, 1, 1, 6.28721e01);
  // sn.setTrilinearElement(UA, 1, 2, 2.58731e02);
  // sn.setTrilinearElement(UA, 1, 3, 8.88614e02);
  // sn.setTrilinearElement(UA, 2, 1, 2.47165e03);
  sn.setTrilinearElement(UA, 2, 2, 2.49731e01);
  // sn.setTrilinearElement(UA, 2, 3, 6.29814e01);
  // sn.setTrilinearElement(UA, 3, 1, 7.49822e02);
  // sn.setTrilinearElement(UA, 3, 2, 2.59841e02);
  sn.setTrilinearElement(UA, 3, 3, 5.29713e01);


  sn.setTrilinearElement(DA, 1, 1, 3.68711e01);
  // sn.setTrilinearElement(DA, 1, 2, 2.49875e02);
  // sn.setTrilinearElement(DA, 1, 3, 9.65321e02);
  // sn.setTrilinearElement(DA, 2, 1, 2.38752e03);
   sn.setTrilinearElement(DA, 2, 2, 1.11123e01);
  // sn.setTrilinearElement(DA, 2, 3, 4.29718e01);
  // sn.setTrilinearElement(DA, 3, 1, 5.98913e02);
  // sn.setTrilinearElement(DA, 3, 2, 8.32761e02);
  sn.setTrilinearElement(DA, 3, 3, 2.95872e01);

   sn.setTrilinearElement(EA, 1, 1, 1.7115e01);  
   // sn.setTrilinearElement(EA, 1, 2, 2.5971e01);
   // sn.setTrilinearElement(EA, 1, 3, 5.3919e02);  
   // sn.setTrilinearElement(EA, 2, 1, 3.8711e03);
   sn.setTrilinearElement(EA, 2, 2, 2.2223e02);  
   // sn.setTrilinearElement(EA, 2, 3, 1.3974e01);
   // sn.setTrilinearElement(EA, 3, 1, 5.8713e02);  
   // sn.setTrilinearElement(EA, 3, 2, 6.9821e01);
   sn.setTrilinearElement(EA, 3, 3, 2.5871e02);  
  

  
     sn.setGauginoMass(2, 1e03);
     sn.setGauginoMass(1, 2e03);
     sn.setGauginoMass(3, 3e03);

 // cout << "alambda = "  << sn.displayTrialambda() << endl;
 //  cout << "akappa = "  << sn.displayTriakappa() << endl;
 
     cout << "sn at mx just before run down = " << sn.display() << endl;
   double err =  sn.run(mx, MZ, eps);
     if (err) {
	cout << "itLowsoft gone non-perturbative approaching mgut\n"; 
     }
     cout <<"sn at mz = " << sn.display() << endl;
     
     err = sn.run(MZ, mx, eps);
     if (err) {
        cout << "itLowsoft gone non-perturbative approaching mgut\n"; 
     }
     cout << " sn at MX"  << sn.display() << endl;
     




  // QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  // /// most important Standard Model inputs: you may change these and recompile
  // double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
  // oneset.setAlpha(ALPHAS, alphasMZ);
  // oneset.setPoleMt(mtop);
  // oneset.setMass(mBottom, mbmb);

  // oneset.toMz();      ///< Runs SM fermion masses to MZ
   
 //Create an empty object
//   nMssmSusy n2; 
// cout << n2.display() << endl;
// //Fill it with values at the Highscale, Mx
//  n2.setYukawaElement(YU, 1, 1,  3.220453e-06 ); 
//  n2.setYukawaElement(YU, 2, 2, 1.463900e-03); 
// n2.setYukawaElement(YU, 3, 3,  4.986557e-01); 

// n2.setYukawaElement(YD, 1, 1, 8.476048e-05); 
//  n2.setYukawaElement(YD, 2, 2, 1.855824e-03); 
// n2.setYukawaElement(YD, 3, 3,  8.675091e-02 );

//  n2.setYukawaElement(YE, 1, 1, 3.330709e-05 ); 
//  n2.setYukawaElement(YE, 2, 2,  6.887167e-03  ); 
// n2.setYukawaElement(YE, 3, 3,   1.229572e-01 );

//   mx = 1.860215e+16;
 
// eps = 1.000000e-12;
//  n2.setGaugeCoupling(1,  7.134806e-01 );
// n2.setGaugeCoupling(2, 7.134388e-01 );
// n2.setGaugeCoupling(3, 7.046983e-01);
//   msusy = 8.797263e+02;
//  //cout after filling
//  cout << "n initiallisd to " << n2.display() << endl;
//  //run to msusy and output
//  n2.run(mx, msusy, eps);
//  cout << "n at msusy = " << n2.display() << endl;
//  //run from msusy to MZ and output 
//  n2.run(msusy, MZ, eps);
//  cout <<"n at mz = " << n2.display() << endl;

//  n2.runto(mx, eps);
//  cout << "should be starting values = " << n.display() << endl;
 



  return 0;
}
