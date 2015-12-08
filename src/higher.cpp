
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

#ifdef COMPILE_TWO_LOOP_GAUGE_YUKAWA
#ifdef COMPILE_THREE_LOOP_RGE

/// NLLFAST version
void getCrossSection(MssmSoftsusy & r, double m0, double m12, double a0, 
		     double tanb, double & xsGG, double & xsSG, double & xsSS,
		     double & xsSB, double & xsTB) {
  double mg = r.displayPhys().mGluino;
  double mt1 = minimum(r.displayPhys().mu(1, 3), r.displayPhys().mu(2, 3));
  double mq  = (r.displayPhys().mu(1, 1) + r.displayPhys().md(1, 1) +
		r.displayPhys().mu(1, 2) + r.displayPhys().md(1, 2)) * 0.25;

  char buff[500];
  sprintf(buff, "cd /home/bca20/code/nllfast-3.0-13TeV; ./nllfast13 gg mstw %f %f > output 2> err; ./nllfast13 sg mstw %f %f >> output 2>> err; ./nllfast13 ss mstw %f %f >> output 2>> err; ./nllfast13 sb mstw %f %f >> output 2>> err; ./nllfast13 st mstw %f >> output 2>> err",mq,mg,mq,mg,mq,mg,mq,mg,mt1);
  //  cout << buff << endl;
  int err = system(buff);
  xsGG = 0.; xsSG = 0.; xsSS = 0.; xsSB = 0.; xsTB = 0.;

  char c[500];
  char fn[500];
  sprintf(fn, "/home/bca20/code/nllfast-3.0-13TeV/output");
  if (!err) {
    fstream fin(fn, ios::in); 
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> xsGG >> c >> c >> c >> c >> c >> c >> c >> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> xsSG >> c >> c >> c >> c >> c >> c >> c >> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> xsSS >> c >> c >> c >> c >> c >> c >> c >> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> xsSB >> c >> c >> c >> c >> c >> c >> c >> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> xsTB >> c >> c >> c >> c >> c >> c >> c >> c >> c;
    fin.close();
  } else  cout << "CROSS SECTION ERROR\n" << xsGG;
   
  return;
}

/* PROSPINO VERSION
void getCrossSection(MssmSoftsusy & r, double m0, double m12, double a0, 
		     double tanb, double & xsGG, double & xsSG, double & xsSS,
		     double & xsSB, double & xsTB) {
  /// First, make a SLHA file
  DoubleVector pars(3); 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  bool ewsbBCscale = false; int sgnMu = 1;
  const char* modelIdent = "sugra"; 
  double qMax = 0.;
  char fileName[500];
  sprintf(fileName,"/home/bca20/code/prospino2.1/prospino.in.les_houches");
  fstream fout(fileName, ios::out);
  fout.setf(ios::scientific, ios::floatfield);
  r.lesHouchesAccordOutput(fout, modelIdent, pars, sgnMu, tanb, qMax, 
			   0, ewsbBCscale);
  fout.close();

  char buff[500];
  sprintf(buff, "cd /home/bca20/code/prospino2.1; ./prospino_2.run > output 2> err");
  //     cout << buff << endl;
  int err = system(buff);
  xsGG = 0.; xsSG = 0.; xsSS = 0.; xsSB = 0.; xsTB = 0.;

  char c[500];
  char fn[500];
  sprintf(fn, "/home/bca20/code/prospino2.1/prospino.dat");
  if (!err) {
    fstream fin(fn, ios::in); 
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> xsGG >> c >> c >> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> xsSG >> c >> c >> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> xsSS >> c >> c >> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> xsSB >> c >> c >> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> xsTB >> c >> c >> c >> c;
    fin.close();
  } else  cout << "CROSS SECTION ERROR\n" << xsGG;
   
  return;
}
*/
void doScan(double lowRatio, double highRatio, int numPoints) {
    /// Sets format of output: 6 decimal places
    outputCharacteristics(9);

    double m12 = 1000., m0 = 0., m0Overm12 = 0., a0 = 0., tanb = 10.;

    int sgnMu = 1;      ///< sign of mu parameter 
    
    QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

    /// most important Standard Model inputs: you may change these and recompile
    double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
     // double lowRatio = 0.7817552, highRatio = 0.7817581;
    oneset.setAlpha(ALPHAS, alphasMZ);
    oneset.setPoleMt(mtop);
    oneset.setMbMb(mbmb);    
    oneset.toMz();      ///< Runs SM fermion masses to MZ
    

    DoubleVector pars(3); 
    bool uni = true, ewsbBCscale = false; double mGutGuess = 1.e16;
    int i; for (i=0; i<=numPoints; i++) {
      USE_TWO_LOOP_SPARTICLE_MASS = false;
      USE_TWO_LOOP_GAUGE_YUKAWA = false;
      USE_THREE_LOOP_RGE = false;

      MssmSoftsusy r; 
      
      m0Overm12 = lowRatio + 
	(highRatio - lowRatio) / double(numPoints) * double(i);
      m0 = m0Overm12 * m12;
      pars(1) = m0; pars(2) = m12; pars(3) = a0; 
      /// Calculate the 
      r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni, 
	       ewsbBCscale);
      
      MssmSoftsusy ho;
      USE_TWO_LOOP_SPARTICLE_MASS = true;
      ho.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni, 
		ewsbBCscale);
      //      ho.lesHouchesAccordOutput(cout, modelIdent, pars, sgnMu, tanb, qMax, 
      //			       0, ewsbBCscale);

      if (r.displayProblem().test()) cout << "# ";
      cout << m0Overm12 << " "                 // 1
	   << r.displayPhys().mGluino << " "   // 2
	   << ho.displayPhys().mGluino << " "  // 3
	   << r.displayPhys().mu(1, 1) << " "  // 4
	   << ho.displayPhys().mu(1, 1) << " " // 5
	   << r.displayPhys().mu(1, 3) << " "  // 6
	   << ho.displayPhys().mu(1, 3) << " " // 7
	   << r.displayPhys().mu(2, 3) << " "  // 8
	   << ho.displayPhys().mu(2, 3) << " " // 9
	   << r.displayPhys().md(1, 1) << " "  // 10
	   << ho.displayPhys().md(1, 1) << " " // 11
	   << r.displayPhys().md(1, 3) << " "  // 12
	   << ho.displayPhys().md(1, 3) << " " // 13
	   << r.displayPhys().mu(2, 1) << " "  // 14
	   << ho.displayPhys().mu(2, 1) << " " // 15
	   << r.displayPhys().md(2, 1) << " "  // 16
	   << ho.displayPhys().md(2, 1) << " " // 17
	   << r.displayPhys().md(2, 3) << " "  // 18
	   << ho.displayPhys().md(2, 3) << " ";// 19

      r.runto(r.displayMsusy());
      r.calcDrBarPars();
      cout << r.displayDrBarPars().mGluino  << " " // 20
	   << r.displayDrBarPars().mu(1, 1) << " " // 21
	   << r.displayDrBarPars().mu(1, 3) << " " // 22
	   << r.displayDrBarPars().mu(2, 3) << " " // 23
	   << r.displayDrBarPars().md(1, 1) << " " // 24
	   << r.displayDrBarPars().md(1, 3) << " " // 25
	   << r.displayDrBarPars().mu(2, 1) << " " // 26
	   << r.displayDrBarPars().md(2, 1) << " " // 27
	   << r.displayDrBarPars().md(2, 3) << " " // 28
	   << r.displayDrBarPars().mt       << " " // 29
	   << r.displayDrBarPars().mneut(1) << " " // 30
	   << r.displayDrBarPars().mneut(2) << " " // 31
	   << r.displayDrBarPars().mneut(3) << " " // 32
	   << r.displayDrBarPars().mneut(4) << " ";// 33

      /// calculate 13 TeV cross-sections
      double xsGG, xsSG, xsSS, xsSB, xsTB;
      double xsGGho, xsSGho, xsSSho, xsSBho, xsTBho;
      getCrossSection(r, m0, m12, a0, tanb, xsGG, xsSG, xsSS, xsSB, xsTB);
      getCrossSection(ho, m0, m12, a0, tanb, xsGGho, xsSGho, xsSSho, xsSBho, 
		      xsTBho);
      cout 
	<< xsGG << " "   // 34
	<< xsGGho << " " // 35
	<< xsSG << " "   // 36
	<< xsSGho << " " // 37
	<< xsSS << " "   // 38
	<< xsSSho << " " // 39
	<< xsSB << " "   // 40
	<< xsSBho << " " // 41
	<< xsTB << " "   // 42
	<< xsTBho;       // 43

      if (r.displayProblem().test()) cout << " " << r.displayProblem();      
      cout << endl;
    }
}

#endif
#endif


int main(int argc, char *argv[]) {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  TOLERANCE = 1.0e-4;
  try {
#ifdef COMPILE_TWO_LOOP_GAUGE_YUKAWA
#ifdef COMPILE_THREE_LOOP_RGE
    doScan(0.1, 4.5, 20);
    //    doScan(1.96, 1.9786, 10);
    //    doScan(1.9786, 1.9825, 10);
    //    doScan(1.9825, 2.03, 10);
    //    doScan(2.03, 3.2, 10);
#endif 
#endif
  }
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}

