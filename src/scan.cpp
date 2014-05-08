
/** 
   Project:     SOFTSUSY 
   File:        main.cpp
   Author:      Ben Allanach 
   Manual:      B.C. Allanach,hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
                B.C. Allanach, M.A. Bernhardt, arXiv:0903.1805, Comp. Phys. 
		Commun. 181 (2010) 232-245
   Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   Description: main calling program example: performs a scan of tan beta 
   (starting at CMSSM10.1.1) and prints out Higgs masses as a result in
   the format
       <tan beta>      <mh0>     <mA0>    <mH0>     <mH+/->
*/

#include <iostream>
#include "mycomplex.h"
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "rge.h"
#include "softsusy.h"
#include "softpars.h"
#include "susy.h"
#include "utils.h"
#include "numerics.h"

void produceSLHAfile(MssmSoftsusy & t, const char * fileName, int sgnMu, 
		    double tanb, const DoubleVector & pars) {
  const char * modelIdent = "sugra"; 
  double qMax = 0.; int numPoints = 1; 
  bool altEwsb = false;
  fstream fout(fileName, ios::out);
  fout.setf(ios::scientific, ios::floatfield);
  fout.precision(10);  
  // new softsusy call
  t.lesHouchesAccordOutput(fout, modelIdent, pars, sgnMu, tanb, qMax, 
			   numPoints, altEwsb);
  fout.close();
}

double doDarkMatter(DoubleVector & pars, double tanb, int sgnMu, 
		    char * fileName) {
  double m0 = pars(1), m12 = pars(2), a0 = pars(3);
  char oFile[500], buff[500];
  sprintf(oFile,"om_%d_%d_%d_%d_%d", int(m0), int(m12), int(a0), 
	  int(tanb), int(sgnMu));
  sprintf(buff,"../../code/micromegas_3.3.13/MSSM/main %s > %s",
	  fileName, oFile);
  int err = system(buff);
  double omega = 0.;
  if (!err) //throw("Problem in micromegas system call: \n");
    { fstream fin2(oFile, ios::in); fin2 >> omega; fin2.close(); }
  
  remove(oFile);
  return omega;
}

/// Returns the object along with omega. Oneset should already be fixed at MZ
void getCmssmAndOmega(MssmSoftsusy & r, DoubleVector & pars, const double tanb, 
			      const int sgnMu, const QedQcd & oneset, 
			      double mGutGuess, bool uni, double & omega, 
			      double & msqAv) {
  double m0=pars.display(1), m12 = pars.display(2), a0 = pars.display(3);
  r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
  r.setData(oneset);

  /// Produces SLHA output file
  char fileName[500]; 
  sprintf(fileName,"lesHout_%d_%d_%d_%d_%d",int(m0),int(m12),int(a0),
	  int(tanb),int(sgnMu));
  produceSLHAfile(r, fileName, sgnMu, tanb, pars);
  if (!r.displayProblem().test()) 
    omega = doDarkMatter(pars, tanb, sgnMu, fileName);

  msqAv = (r.displayPhys().mu(2, 1) + 
		     r.displayPhys().mu(1, 1) +
		     r.displayPhys().md(2, 1) + 
		     r.displayPhys().md(1, 1)) * 
    0.25;
  
  //  remove(fileName);
  return;
}

int main(int argc, char *argv[]) {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  try {
 /// Sets format of output: 6 decimal places
  outputCharacteristics(6);

  /// Parameters used: CMSSM parameters
  double m12 = 500., a0 = 0., mGutGuess = 2.0e16, tanb = 30.0, m0 = 125.;
  int sgnMu = 1;      ///< sign of mu parameter 

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.2, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);
  oneset.toMz();      ///< Runs SM fermion masses to MZ

  if (argc != 6) { exit(-1); }

  m0 = atof(argv[1]);
  m12 = atof(argv[2]);
  a0 = atof(argv[3]);
  tanb = atof(argv[4]);
  sgnMu = atoi(argv[5]);

  /// Preparation for calculation: set up object and input parameters
  DoubleVector pars(3); 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)

  /// Switch off 3-loop RGEs etc
  double omega2=0., msqAv2 = 0.;  
  USE_THREE_LOOP_RGE = false;   USE_TWO_LOOP_THRESHOLD = false;
  MssmSoftsusy twoLoop;
  getCmssmAndOmega(twoLoop, pars, tanb, sgnMu, oneset, mGutGuess, 
		   uni, omega2, msqAv2);

  /// Just 2-loop thresholds for strong coupling constant
  double omegaAs = 0., msqAvAs = 0.; mGutGuess = 2.e16;
  MssmSoftsusy twoLoopAs; 
  twoLoopAs.included_thresholds |= ENABLE_TWO_LOOP_AS_AS_YUK;
  USE_TWO_LOOP_THRESHOLD = true;
  getCmssmAndOmega(twoLoopAs, pars, tanb, sgnMu, oneset, mGutGuess, 
		   uni, omegaAs, msqAvAs); 

  /// Just 2-loop strong thresholds for mt
  USE_TWO_LOOP_THRESHOLD = false;
  double omegaMt = 0., msqAvMt = 0.; mGutGuess = 2.e16;
  MssmSoftsusy twoLoopMt; 
  twoLoopMt.included_thresholds |= ENABLE_TWO_LOOP_MT_AS;
  USE_TWO_LOOP_THRESHOLD = true;
  getCmssmAndOmega(twoLoopMt, pars, tanb, sgnMu, oneset, mGutGuess, 
		   uni, omegaMt, msqAvMt); 

  /// Just 2-loop for mb,mtau
  USE_TWO_LOOP_THRESHOLD = false;
  double omegaMb = 0., msqAvMb = 0.; mGutGuess = 2.e16;
  MssmSoftsusy twoLoopMb; 
  twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MB_AS;
  twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MB_YUK;
  twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MTAU_YUK;
  USE_TWO_LOOP_THRESHOLD = true;
  getCmssmAndOmega(twoLoopMb, pars, tanb, sgnMu, oneset, mGutGuess, 
		   uni, omegaMb, msqAvMb); 

  /// 3-loop etc ON
  double omega3 = 0., msqAv3 = 0.; mGutGuess = 2.0e16;
  USE_THREE_LOOP_RGE = true;
  USE_TWO_LOOP_THRESHOLD = true;
  MssmSoftsusy threeLoop;
  getCmssmAndOmega(threeLoop, pars, tanb, sgnMu, oneset, mGutGuess, 
					  uni, omega3, msqAv3); 

  /// check the point in question is problem free: if so print the output
  //  if (twoLoop.displayProblem().test() ||
  //  threeLoop.displayProblem().test()) cout << "# ";           ///< column
  cout << m0                                                     ///< 1
       << " " << m12                                             ///< 2
       << " " << a0                                              ///< 3
       << " " << tanb                                            ///< 4
       << " " << twoLoop.displayPhys().mh0(1)                    ///< 5
       << " " << twoLoopAs.displayPhys().mh0(1)                  ///< 6
       << " " << twoLoopMt.displayPhys().mh0(1)                  ///< 7
       << " " << twoLoopMb.displayPhys().mh0(1)                  ///< 8
       << " " << threeLoop.displayPhys().mh0(1)                  ///< 9
       << " " << twoLoop.displayPhys().mA0(1)                    ///< 10
       << " " << twoLoopAs.displayPhys().mA0(1)                  ///< 11
       << " " << twoLoopMt.displayPhys().mA0(1)                  ///< 12
       << " " << twoLoopMb.displayPhys().mA0(1)                  ///< 13
       << " " << threeLoop.displayPhys().mA0(1)                  ///< 14
       << " " << twoLoop.displayPhys().mh0(2)                    ///< 15
       << " " << twoLoopAs.displayPhys().mh0(2)                  ///< 16
       << " " << twoLoopMt.displayPhys().mh0(2)                  ///< 17
       << " " << twoLoopMb.displayPhys().mh0(2)                  ///< 18
       << " " << threeLoop.displayPhys().mh0(2)                  ///< 19
       << " " << twoLoop.displayPhys().mHpm                      ///< 20
       << " " << twoLoopAs.displayPhys().mHpm                    ///< 21
       << " " << twoLoopMt.displayPhys().mHpm                    ///< 22
       << " " << twoLoopMb.displayPhys().mHpm                    ///< 23
       << " " << threeLoop.displayPhys().mHpm                    ///< 24
       << " " << twoLoop.displayPhys().mGluino                   ///< 25
       << " " << twoLoopAs.displayPhys().mGluino                 ///< 26
       << " " << twoLoopMt.displayPhys().mGluino                 ///< 27
       << " " << twoLoopMb.displayPhys().mGluino                 ///< 28
       << " " << threeLoop.displayPhys().mGluino                 ///< 29
       << " " << msqAv2                                          ///< 30
       << " " << msqAvAs                                         ///< 31 
       << " " << msqAvMt                                         ///< 32
       << " " << msqAvMb                                         ///< 33
       << " " << msqAv3                                          ///< 34
       << " " << twoLoop.displayPhys().me(1, 1)                  ///< 35 
       << " " << twoLoopAs.displayPhys().me(1, 1)                ///< 36 
       << " " << twoLoopMt.displayPhys().me(1, 1)                ///< 37 
       << " " << twoLoopMb.displayPhys().me(1, 1)                ///< 38 
       << " " << threeLoop.displayPhys().me(1, 1)                ///< 39 
       << " " << twoLoop.displayPhys().me(2, 1)                  ///< 40 
       << " " << twoLoopAs.displayPhys().me(2, 1)                ///< 41 
       << " " << twoLoopMt.displayPhys().me(2, 1)                ///< 42
       << " " << twoLoopMb.displayPhys().me(2, 1)                ///< 43
       << " " << threeLoop.displayPhys().me(2, 1)                ///< 44
       << " " << fabs(twoLoop.displayPhys().mneut(1))            ///< 45 
       << " " << fabs(twoLoopAs.displayPhys().mneut(1))          ///< 46
       << " " << fabs(twoLoopMt.displayPhys().mneut(1))          ///< 47
       << " " << fabs(twoLoopMb.displayPhys().mneut(1))          ///< 48
       << " " << fabs(threeLoop.displayPhys().mneut(1))          ///< 49
       << " " << fabs(twoLoop.displayPhys().mneut(2))            ///< 50
       << " " << fabs(twoLoopAs.displayPhys().mneut(2))          ///< 51
       << " " << fabs(twoLoopMt.displayPhys().mneut(2))          ///< 52
       << " " << fabs(twoLoopMb.displayPhys().mneut(2))          ///< 53
       << " " << fabs(threeLoop.displayPhys().mneut(2))          ///< 54
       << " " << fabs(twoLoop.displayPhys().mneut(3))            ///< 55
       << " " << fabs(twoLoopAs.displayPhys().mneut(3))          ///< 56
       << " " << fabs(twoLoopMt.displayPhys().mneut(3))          ///< 57 
       << " " << fabs(twoLoopMb.displayPhys().mneut(3))          ///< 58
       << " " << fabs(threeLoop.displayPhys().mneut(3))          ///< 59
       << " " << fabs(twoLoop.displayPhys().mneut(4))            ///< 60
       << " " << fabs(twoLoopAs.displayPhys().mneut(4))          ///< 61
       << " " << fabs(twoLoopMt.displayPhys().mneut(4))          ///< 62
       << " " << fabs(twoLoopMb.displayPhys().mneut(4))          ///< 63
       << " " << fabs(threeLoop.displayPhys().mneut(4))          ///< 64
       << " " << twoLoop.displayPhys().mu(1, 3)                  ///< 65
       << " " << twoLoopAs.displayPhys().mu(1, 3)                ///< 66
       << " " << twoLoopMt.displayPhys().mu(1, 3)                ///< 67
       << " " << twoLoopMb.displayPhys().mu(1, 3)                ///< 68
       << " " << threeLoop.displayPhys().mu(1, 3)                ///< 69
       << " " << twoLoop.displayPhys().mu(2, 3)                  ///< 70
       << " " << twoLoopAs.displayPhys().mu(2, 3)                ///< 71
       << " " << twoLoopMt.displayPhys().mu(2, 3)                ///< 72
       << " " << twoLoopMb.displayPhys().mu(2, 3)                ///< 73
       << " " << threeLoop.displayPhys().mu(2, 3)                ///< 74
       << " " << twoLoop.displayPhys().md(1, 3)                  ///< 75
       << " " << twoLoopAs.displayPhys().md(1, 3)                ///< 76
       << " " << twoLoopMt.displayPhys().md(1, 3)                ///< 77
       << " " << twoLoopMb.displayPhys().md(1, 3)                ///< 78
       << " " << threeLoop.displayPhys().md(1, 3)                ///< 79
       << " " << twoLoop.displayPhys().md(2, 3)                  ///< 80
       << " " << twoLoopAs.displayPhys().md(2, 3)                ///< 81
       << " " << twoLoopMt.displayPhys().md(2, 3)                ///< 82
       << " " << twoLoopMb.displayPhys().md(2, 3)                ///< 83
       << " " << threeLoop.displayPhys().md(2, 3)                ///< 84
       << " " << twoLoop.displayPhys().me(1, 3)                  ///< 85
       << " " << twoLoopAs.displayPhys().me(1, 3)                ///< 86
       << " " << twoLoopMt.displayPhys().me(1, 3)                ///< 87
       << " " << twoLoopMb.displayPhys().me(1, 3)                ///< 88
       << " " << threeLoop.displayPhys().me(1, 3)                ///< 89
       << " " << twoLoop.displayPhys().me(2, 3)                  ///< 90
       << " " << twoLoopAs.displayPhys().me(2, 3)                ///< 91
       << " " << twoLoopMt.displayPhys().me(2, 3)                ///< 92
       << " " << twoLoopMb.displayPhys().me(2, 3)                ///< 93
       << " " << threeLoop.displayPhys().me(2, 3)                ///< 94
       << " " << twoLoop.displayYukawaElement(YU, 3, 3)          ///< 95
       << " " << twoLoopAs.displayYukawaElement(YU, 3, 3)        ///< 96
       << " " << twoLoopMt.displayYukawaElement(YU, 3, 3)        ///< 97
       << " " << twoLoopMb.displayYukawaElement(YU, 3, 3)        ///< 98
       << " " << threeLoop.displayYukawaElement(YU, 3, 3)        ///< 99
       << " " << twoLoop.displayYukawaElement(YD, 3, 3)          ///< 100
       << " " << twoLoopAs.displayYukawaElement(YD, 3, 3)        ///< 101
       << " " << twoLoopMt.displayYukawaElement(YD, 3, 3)        ///< 102
       << " " << twoLoopMb.displayYukawaElement(YD, 3, 3)        ///< 103
       << " " << threeLoop.displayYukawaElement(YD, 3, 3)        ///< 104
       << " " << twoLoop.displayYukawaElement(YE, 3, 3)          ///< 105
       << " " << twoLoopAs.displayYukawaElement(YE, 3, 3)        ///< 106
       << " " << twoLoopMt.displayYukawaElement(YE, 3, 3)        ///< 107
       << " " << twoLoopMb.displayYukawaElement(YE, 3, 3)        ///< 108
       << " " << threeLoop.displayYukawaElement(YE, 3, 3);       ///< 109
  int facMusq3= 1., facMusq2 = 1., facMusq2As = 1., facMusq2Mt = 1., 
    facMusq2Mb = 1.;
  if (twoLoop.displayProblem().muSqWrongSign) facMusq2 = -1. ;
  if (twoLoopAs.displayProblem().muSqWrongSign) facMusq2As = -1. ;
  if (twoLoopMt.displayProblem().muSqWrongSign) facMusq2Mt = -1. ;
  if (twoLoopMb.displayProblem().muSqWrongSign) facMusq2Mb = -1. ;
  if (threeLoop.displayProblem().muSqWrongSign) facMusq3 = -1. ;
  cout << " " << twoLoop.displaySusyMu()  * facMusq2             ///< 110
       << " " << twoLoopAs.displaySusyMu()  * facMusq2As         ///< 111
       << " " << twoLoopMt.displaySusyMu()  * facMusq2Mt         ///< 112
       << " " << twoLoopMb.displaySusyMu()  * facMusq2Mb         ///< 113
       << " " << threeLoop.displaySusyMu()  * facMusq3           ///< 114
       << " " << fabs(twoLoop.displayPhys().mch(1))              ///< 115
       << " " << fabs(twoLoopAs.displayPhys().mch(1))            ///< 116
       << " " << fabs(twoLoopMt.displayPhys().mch(1))            ///< 117
       << " " << fabs(twoLoopMb.displayPhys().mch(1))            ///< 118
       << " " << fabs(threeLoop.displayPhys().mch(1))            ///< 119
       << " " << fabs(twoLoop.displayPhys().mch(2))              ///< 120
       << " " << fabs(twoLoopAs.displayPhys().mch(2))            ///< 121
       << " " << fabs(twoLoopMt.displayPhys().mch(2))            ///< 122
       << " " << fabs(twoLoopMb.displayPhys().mch(2))            ///< 123
       << " " << fabs(threeLoop.displayPhys().mch(2))            ///< 124
       << " " << omega2                                          ///< 125
       << " " << omegaAs                                         ///< 126
       << " " << omegaMt                                         ///< 127
       << " " << omegaMb                                         ///< 128
       << " " << omega3                                         ///< 129
       << " " << twoLoop.displayGaugeCoupling(3)                ///< 130
       << " " << twoLoopAs.displayGaugeCoupling(3)              ///< 131
       << " " << twoLoopMt.displayGaugeCoupling(3)              ///< 132
       << " " << twoLoopMb.displayGaugeCoupling(3)              ///< 133
       << " " << threeLoop.displayGaugeCoupling(3);             ///< 134

      if (twoLoop.displayProblem().test()) cout << "# 2-loop problem: " 
					  << twoLoop.displayProblem();
       if (threeLoop.displayProblem().test()) cout << "# 3-loop problem " 
					  << threeLoop.displayProblem(); 
      cout << endl;
    }

  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}
