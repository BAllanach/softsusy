/** 
   Project:     SOFTSUSY 
   File:        higher.cpp
   Author:      B.C. Allanach, S. Martin, D. Robertson and R. Ruiz de Austri
   Manual:      B.C. Allanach, S. Martin, D. Robertson and R. Ruiz de Austri, 
                arXiv:1601.06657
   Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   Description: main calling program to be put in src/ in the softsusy
                directory. Running make will produce an executable called
                "higher.x", which takes no arguments. Running "higher.x"
                produces the data given in "twoLoop.dat", which were used to 
                make the figures in the manual, and which illustrate the size 
                of two-loop contributions to gluino and squark masses and the
                resulting corrections to gluino and squark production 
                cross-sections at the 13 TeV LHC. You will need to install 
                NLL-fast (google it) in some directory, and compile it to get 
                an executable which you should call nllfast13. Then change 
                the variable fullPathToNllFast just
                below to point to the directory you've placed nllfast13 in. 
*/

#include "higher.h"

/// Change this variable to the location in your code of nllfast
const char * fullPathToNllFast = "/home/bca20/code/nllfast-3.1-13TeV/";

bool treeLevelGluino = false;

/// NLLFAST version
void getCrossSection(MssmSoftsusy & r, double m0, double m12, double a0, 
		     double tanb, double & xsGG, double & xsSG, double & xsSS,
		     double & xsSB, double & xsTB) {
  double mg = r.displayPhys().mGluino;
  double mt1 = minimum(r.displayPhys().mu(1, 3), r.displayPhys().mu(2, 3));
  double mq  = (r.displayPhys().mu(1, 1) + r.displayPhys().md(1, 1) +
		r.displayPhys().mu(1, 2) + r.displayPhys().md(1, 2)) * 0.25;

  char buff[500];

  sprintf(buff, "cd %s; ./nllfast13 gg mstw %f %f > output 2> err; ./nllfast13 sg mstw %f %f >> output 2>> err; ./nllfast13 ss mstw %f %f >> output 2>> err; ./nllfast13 sb mstw %f %f >> output 2>> err; ./nllfast13 st mstw %f >> output 2>> err",fullPathToNllFast,mq,mg,mq,mg,mq,mg,mq,mg,mt1);

  int err = system(buff);
  xsGG = 0.; xsSG = 0.; xsSS = 0.; xsSB = 0.; xsTB = 0.;

  string c;
  if (!err) {
    char fn[500]; sprintf(fn,"%s/output",fullPathToNllFast);
    fstream fin(fn, ios::in); 
    // Following lines corrected Jan. 15 2016 by SPM, to give NLL+NLO 
    // cross-sections. Previously, xsGG, xsSG, xsSS, and xsSB were 
    // all NLO only, and xsTB was d_mu+.
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c
    	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> xsGG >> c >> c >> c >> c >> c >> c 
	>> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> xsSG >> c >> c >> c >> c >> c >> c 
	>> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> xsSS >> c >> c >> c >> c >> c >> c 
	>> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> xsSB >> c >> c >> c >> c >> c >> c 
	>> c >> c;
    fin >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c >> c 
	>> c >> c >> c >> xsTB >> c >> c >> c >> c >> c >> c >> c >> c
	>> c >> c;
    fin.close();

  } else  cout << "CROSS SECTION ERROR\n" << xsGG;

  /* For debugging purposes; compare with command line output of nllfast. SPM.
  cout << "\n" << "squark, gluino, stop1 masses " << mq << " " << mg << " " << mt1 << "\n";
  cout << "gg " << xsGG << "\n";   
  cout << "sg " << xsSG << "\n";   
  cout << "ss " << xsSS << "\n";   
  cout << "sb " << xsSB << "\n";   
  cout << "tb " << xsTB << "\n\n";   
  exit(0);
  */

  return;
}

void scaleVariation() {
    /// Sets format of output: 6 decimal places
    outputCharacteristics(9);

    double m12 = 1000., m0 = 1000., a0 = 0., tanb = 10.;

    int sgnMu = 1;      ///< sign of mu parameter 
    
    QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

    /// most important Standard Model inputs: you may change these and recompile
    double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
     // double lowRatio = 0.7817552, highRatio = 0.7817581;
    oneset.setAlpha(ALPHAS, alphasMZ);
    oneset.setPoleMt(mtop);
    oneset.setMbMb(mbmb);    
    oneset.toMz();      ///< Runs SM fermion masses to MZ
    
    double qewsbStart = 0.5, qewsbEnd = 2.; int numPoints = 20;

    cout << "# m12 = " << m12 << ", m0 = " << m0 << "\n";
    cout << "# A0 = " << a0 << ", tanBeta = " << tanb << "\n";
    cout << "# alphasMZ = " << alphasMZ << ", Mt = " << mtop << ", mb(mb) = " << mbmb << "\n";

    if (treeLevelGluino) {
      cout << "# Q/MSUSY       Q/GeV           M3(Q)\n";
    } else {
      cout << "# Q/MSUSY       Q/GeV           mgluinopole     msuLpole       	mstop1pole      mstop2pole\n";
    }

    for (int i=0; i<=numPoints; i++) {
      double qewsb = (qewsbEnd - qewsbStart) / double(numPoints) * double(i) + 
	qewsbStart; 
      DoubleVector pars(3); 
      bool uni = true, ewsbBCscale = false; double mGutGuess = 1.e16;
      
      pars(1) = m0; pars(2) = m12; pars(3) = a0; 
      /// Calculate the spectrum
      MssmSoftsusy r; r.setQewsb(qewsb);
      r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni, 
	       ewsbBCscale);
      if (treeLevelGluino) {
	r.runto(r.displayMsusy());
	sPhysical aa = r.displayPhys();
	aa.mGluino = r.displayGaugino(3);
	r.setPhys(aa);
      }
      cout << qewsb << " " << r.displayMsusy() << " "  
	   << r.displayPhys().mGluino << " "   // Gluino pole mass
	   << r.displayPhys().mu(1, 1) << " "  // suL pole mass
	   << r.displayPhys().mu(2, 3) << " "  // stop1 pole mass
	   << r.displayPhys().mu(1, 3)         // stop2 pole mass
           << " # " << r.displayProblem() << endl;
    }    
    cout << endl << endl;
}

void doScan(double lowRatio, double highRatio, int numPoints) {

    /// Sets format of output: 6 decimal places
    outputCharacteristics(9);
    expandAroundGluinoPole = 3;

    double m12 = 1000., m0 = 0., m0Overm12 = 0., a0 = -2000., tanb = 10.;

    int sgnMu = 1;      ///< sign of mu parameter 
    
    QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

    /// most important Standard Model inputs: you may change these and recompile
    double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
    oneset.setAlpha(ALPHAS, alphasMZ);
    oneset.setPoleMt(mtop);
    oneset.setMbMb(mbmb);    
    oneset.toMz();      ///< Runs SM fermion masses to MZ

   cout << "# Figure 3(a) is columns 4, 7, 10, and 13 vs. column 1.\n";
   cout << "# Figure 3(b) is columns 40, 43, 46, 49, and 52 vs. column 1.\n";
   cout << "# m12 = " << m12 << ", A0 = " << a0 << ", tanBeta = " << tanb << "\n";
   cout << "# alphasMZ = " << alphasMZ << ", Mt = " << mtop << ", mb(mb) = " << mbmb << "\n";
   cout << "# m0/M12        m_gl(1-loop)    m_gl(2-loop)    delta(m_gl)     muL(1-loop)     muL(2-loop)     delta(muL)      mt1(1-loop)     mt1(2-loop)     delta(mt1)      mt2(1-loop)     mt2(2-loop)     mdL(1-loop)     mdL(2-loop)     mb1(1-loop)     mb1(2-loop)     muR(1-loop)     muR(2-loop)     muR(1-loop)     muR(2-loop)     mb2(1-loop)     mb2(2-loop)     M3(MSUSY)       muL(MSUSY)      mtL(MSUSY)      mtR(MSUSY)      mdL(MSUSY)      mbL(MSUSY)      muR(MSUSY)      mdR(MSUSY)      mbR(MSUSY)      mt(MSUSY)       mneut1(MSUSY)   mneut2(MSUSY)   mneut3(MSUSY)    mneut4(MSUSY)   gg(1-loop)      gg(2-loop)      delta(gg)        sg(1-loop)      sg(2-loop)      delta(sg)        ss(1-loop)      ss(2-loop)      delta(ss)        sb(1-loop)      sb(2-loop)      delta(sb)        tb(1-loop)      tb(2-loop)      delta(tb)\n";

    DoubleVector pars(3); 
    bool uni = true, ewsbBCscale = false; double mGutGuess = 1.e16;
    int i; for (i=0; i<=numPoints; i++) {
      USE_TWO_LOOP_SPARTICLE_MASS = false;

      MssmSoftsusy r; 
      
      m0Overm12 = lowRatio + 
	(highRatio - lowRatio) / double(numPoints) * double(i);
      m0 = m0Overm12 * m12;
      pars(1) = m0; pars(2) = m12; pars(3) = a0; 
      /// Calculate the spectrum
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
	   << ho.displayPhys().mGluino/r.displayPhys().mGluino -1 << " "  // 4
	   << r.displayPhys().mu(1, 1) << " "  // 5
	   << ho.displayPhys().mu(1, 1) << " " // 6
	   << ho.displayPhys().mu(1, 1)/r.displayPhys().mu(1, 1)-1 << " " // 7
	   << r.displayPhys().mu(2, 3) << " "  // 8
	   << ho.displayPhys().mu(2, 3) << " " // 9
	   << ho.displayPhys().mu(2, 3)/r.displayPhys().mu(2, 3)-1 << " " // 10
	   << r.displayPhys().mu(1, 3) << " "  // 11
	   << ho.displayPhys().mu(1, 3) << " " // 12
	   << ho.displayPhys().mu(1, 3)/r.displayPhys().mu(1, 3)-1 << " " // 13
	   << r.displayPhys().md(1, 1) << " "  // 14
	   << ho.displayPhys().md(1, 1) << " " // 15
	   << r.displayPhys().md(1, 3) << " "  // 16
	   << ho.displayPhys().md(1, 3) << " " // 17
	   << r.displayPhys().mu(2, 1) << " "  // 18
	   << ho.displayPhys().mu(2, 1) << " " // 19
	   << r.displayPhys().md(2, 1) << " "  // 20
	   << ho.displayPhys().md(2, 1) << " " // 21
	   << r.displayPhys().md(2, 3) << " "  // 22
	   << ho.displayPhys().md(2, 3) << " ";// 23

      r.runto(r.displayMsusy());
      r.calcDrBarPars();
      cout << r.displayDrBarPars().mGluino  << " " // 24
	   << r.displayDrBarPars().mu(1, 1) << " " // 25
	   << r.displayDrBarPars().mu(1, 3) << " " // 26
	   << r.displayDrBarPars().mu(2, 3) << " " // 27
	   << r.displayDrBarPars().md(1, 1) << " " // 28
	   << r.displayDrBarPars().md(1, 3) << " " // 29
	   << r.displayDrBarPars().mu(2, 1) << " " // 30
	   << r.displayDrBarPars().md(2, 1) << " " // 31
	   << r.displayDrBarPars().md(2, 3) << " " // 32
	   << r.displayDrBarPars().mt       << " " // 33
	   << r.displayDrBarPars().mneut(1) << " " // 34
	   << r.displayDrBarPars().mneut(2) << " " // 35
	   << r.displayDrBarPars().mneut(3) << " " // 36
	   << r.displayDrBarPars().mneut(4) << " ";// 37

      /// calculate 13 TeV cross-sections
      double xsGG, xsSG, xsSS, xsSB, xsTB;
      double xsGGho, xsSGho, xsSSho, xsSBho, xsTBho;
      getCrossSection(r, m0, m12, a0, tanb, xsGG, xsSG, xsSS, xsSB, xsTB);
      getCrossSection(ho, m0, m12, a0, tanb, xsGGho, xsSGho, xsSSho, xsSBho, 
		      xsTBho);
      cout 
	<< xsGG << " "   // 38
	<< xsGGho << " " // 39
	<< xsGGho/(xsGG + 1.0e-20) - 1.0 << " " // 40
	<< xsSG << " "   // 41
	<< xsSGho << " " // 42
	<< xsSGho/(xsSG + 1.0e-20) - 1.0 << " " // 43
	<< xsSS << " "   // 44
	<< xsSSho << " " // 45
	<< xsSSho/(xsSS + 1.0e-20) - 1.0 << " " // 46
	<< xsSB << " "   // 47
	<< xsSBho << " " // 48
	<< xsSBho/(xsSB + 1.0e-20) - 1.0 << " " // 49
	<< xsTB << " "   // 50
	<< xsTBho << " " // 51
	<< xsTBho/(xsTB + 1.0e-20) - 1.0;       // 52

      if (r.displayProblem().test()) cout << " " << r.displayProblem();      
      cout << endl;
    }
}

int main(int argc, char *argv[]) {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  TOLERANCE = 1.0e-4;
  try {
    cout << "# Data file for plots\n\n";

    doScan(0.1, 4.4, 430); cout << endl << endl;

    treeLevelGluino = false;
    USE_TWO_LOOP_SPARTICLE_MASS = false;
    cout << "\n# Figure 4, 1-loop pole masses\n";
    scaleVariation(); 

    treeLevelGluino = false;
    USE_TWO_LOOP_SPARTICLE_MASS = true;
    expandAroundGluinoPole = 3;
    cout << "\n# Figure 4, 2-loop pole masses, gluino re-expanded about both gluino and squark\n";
    scaleVariation(); 

    treeLevelGluino = false;
    USE_TWO_LOOP_SPARTICLE_MASS = true;
    expandAroundGluinoPole = 2;
    cout << "\n# Figure 4, 2-loop pole masses, gluino re-expanded about gluino only\n";
    scaleVariation(); 

    treeLevelGluino = false;
    USE_TWO_LOOP_SPARTICLE_MASS = true;
    expandAroundGluinoPole = 1;
    cout << "\n# Figure 4, 2-loop pole masses, no re-expansion\n";
    scaleVariation(); 
  }
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}


