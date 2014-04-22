
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

  /// Print out header line
  //  cout << "#       1:m0        2:m12         3:a0       4:tanb         5:mh         6:Dmh         7:mA         8:dmA        9:mH        10:DmH        11:mH+       12:DmH+ ";
   // cout << "      13:Dmg " << "       14:mg " << "     15:Dmsq " 
   //     << "      16:msq " << "      17:meL " << "     18: meL " 
   //     << "     19:DmeR " << "      20:meR " << "  21:Dmneut1 " 
   //     << "   22:mneut1 " << "  23:Dmneut2 " << "  24:mneut2  "
   //     << "  25:Dmneut3 " << "  26: mneut3 " << " 27:Dmneut4  " 
   //     << "   28:mneut4 " << "     29:dmtL " << "      30:mtL " 
   //     << "    31:DmtR  " << "     32: mtR " << "    33:DmbL  " 
   //     << "    34: mbL  " << " 35:DmtauL   " << "   36:mtauL  " 
   //     << "     37:dht  " << "       38:ht " << "   39: dhb   " 
   //     << "      40:hb  " << "   41:dhtau  " << "    42:htau  "
   //     << " 43:Dmu      " << "   44:mu     " 
   //     << " 45:Dmchi+1  " << "   46:mchi+1 " << "  47:Dmchi+2 " 
   //     << "   48:mchi+2 " << endl;

  if (argc != 6) { exit(-1); }

  m0 = atof(argv[1]);
  m12 = atof(argv[2]);
  a0 = atof(argv[3]);
  tanb = atof(argv[4]);
  sgnMu = atoi(argv[5]);

  /// Preparation for calculation: set up object and input parameters
  MssmSoftsusy r; 
  DoubleVector pars(3); 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
  
  /// Switch off 3-loop RGEs etc
#ifdef COMPILE_THREE_LOOP_RGE
  USE_THREE_LOOP_RGE = false;
#endif

#ifdef COMPILE_FULL_SUSY_THRESHOLD
  USE_TWO_LOOP_THRESHOLD = false;
#endif
  /// Calculate the spectrum
  r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
  
  double msq2loop = (r.displayPhys().mu(2, 1) + r.displayPhys().mu(1, 1) +
		     r.displayPhys().md(2, 1) + r.displayPhys().md(1, 1)) * 
    0.25;
  
#ifdef COMPILE_THREE_LOOP_RGE
  USE_THREE_LOOP_RGE = true;
#endif

#ifdef COMPILE_FULL_SUSY_THRESHOLD
  USE_TWO_LOOP_THRESHOLD = true;
  /*
  s.included_thresholds &= ~ENABLE_TWO_LOOP_MB_YUK; 
  s.included_thresholds &= ~ENABLE_TWO_LOOP_MB_AS; 
  s.included_thresholds &= ~ENABLE_TWO_LOOP_MT_AS; 
  s.included_thresholds &= ~ENABLE_TWO_LOOP_MTAU_YUK; 
  s.included_thresholds &= ~ENABLE_TWO_LOOP_AS_AS_YUK; 
  */
#endif

  // AVB: USE_TWO_LOOP_THRESHOLD should be set before declaration of Softsusy object
  
  MssmSoftsusy s; mGutGuess = 2.0e16;
  s.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
  
  double msq3loop = (s.displayPhys().mu(2, 1) + s.displayPhys().mu(1, 1) +
		     s.displayPhys().md(2, 1) + s.displayPhys().md(1, 1)) * 
    0.25;
 /* 
  /// check the point in question is problem free: if so print the output
  if (r.displayProblem().test() || s.displayProblem().test()) cout << "# ";
  cout << m0 << " " << m12 << " " << a0 << " " << tanb 
       << " " << r.displayPhys().mh0(1)
       << " " << (1. - r.displayPhys().mh0(1) / s.displayPhys().mh0(1))
       << " " << r.displayPhys().mA0(1)
       << " " << (1. - r.displayPhys().mA0(1) / s.displayPhys().mA0(1))
       << " " << r.displayPhys().mh0(2)
       << " " << (1. - r.displayPhys().mh0(2) / s.displayPhys().mh0(2))
       << " " << r.displayPhys().mHpm 
       << " " << (1. - r.displayPhys().mHpm / s.displayPhys().mHpm)
       << " " << r.displayPhys().mGluino
       << " " << (1. - r.displayPhys().mGluino / s.displayPhys().mGluino)
       << " " << msq2loop
       << " " << (1. - msq2loop / msq3loop)
       << " " << r.displayPhys().me(1, 1)
	   << " " << (1. - r.displayPhys().me(1, 1) / s.displayPhys().me(1, 1))
	   << " " << r.displayPhys().me(2, 1)
	   << " " << (1. - r.displayPhys().me(2, 1) / s.displayPhys().me(2, 1))
	   << " " << fabs(r.displayPhys().mneut(1))
	   << " " << (1. - fabs(r.displayPhys().mneut(1)) / fabs(s.displayPhys().mneut(1)))
	   << " " << fabs(r.displayPhys().mneut(2))
	   << " " << (1. - fabs(r.displayPhys().mneut(2)) / fabs(s.displayPhys().mneut(2)))
	   << " " << fabs(r.displayPhys().mneut(3))
	   << " " << (1. - fabs(r.displayPhys().mneut(3)) / fabs(s.displayPhys().mneut(3)))
	   << " " << fabs(r.displayPhys().mneut(4))
	   << " " << (1. - fabs(r.displayPhys().mneut(4)) / fabs(s.displayPhys().mneut(4)))
	   << " " << fabs(r.displayPhys().mu(1, 3))
	   << " " << (1. - fabs(r.displayPhys().mu(1, 3)) / fabs(s.displayPhys().mu(1, 3)))
	   << " " << r.displayPhys().mu(2, 3)
	   << " " << (1. - r.displayPhys().mu(2, 3) / s.displayPhys().mu(2, 3))
	   << " " << r.displayPhys().md(1, 3)
	   << " " << (1. - r.displayPhys().md(1, 3) / s.displayPhys().md(1, 3))
	   << " " << r.displayPhys().md(2, 3)
	   << " " << (1. - r.displayPhys().md(2, 3) / s.displayPhys().md(2, 3))
	   << " " << r.displayPhys().me(1, 3)
	   << " " << (1. - r.displayPhys().me(1, 3) / s.displayPhys().me(1, 3))
	   << " " << r.displayPhys().me(2, 3)
	   << " " << (1. - r.displayPhys().me(2, 3) / s.displayPhys().me(2, 3))
	   << " " << r.displayYukawaElement(YU, 3, 3)
	   << " " << (1. - r.displayYukawaElement(YU, 3, 3) / 
		      s.displayYukawaElement(YU, 3, 3))
	   << " " << r.displayYukawaElement(YD, 3, 3)
	   << " " << (1. - r.displayYukawaElement(YD, 3, 3) / 
		      s.displayYukawaElement(YD, 3, 3))
	   << " " << r.displayYukawaElement(YE, 3, 3)
	   << " " << (1. - r.displayYukawaElement(YE, 3, 3) / 
		      s.displayYukawaElement(YE, 3, 3))
	   << " " << r.displaySusyMu() 
	   << " " << (1. - r.displaySusyMu() / s.displaySusyMu())
	   << " " << r.displayPhys().mch(1)
	   << " " << (1. - r.displayPhys().mch(1) / s.displayPhys().mch(1))
	   << " " << r.displayPhys().mch(2)
	   << " " << (1. - r.displayPhys().mch(2) / s.displayPhys().mch(2));
      
      if (r.displayProblem().test()) cout << " 2-loop problem: " 
					  << r.displayProblem();
      if (s.displayProblem().test()) cout << " 3-loop problem " 
					  << s.displayProblem(); 
      cout << endl;
*/
 		ofstream tmpfile2,tmpfile3;                                                                                                                                                            
		std::ostringstream sstream;
		sstream << "m0=" << m0  <<"_m12=" << m12 << "_a0=" << a0 << "_tanb=" << tanb << ".slha";
		std::string filename = sstream.str();

		std::string filename2l = std::string("2loop/") + filename;
		std::string filename3l = std::string("3loop/") + filename;
		tmpfile2.open(filename2l.c_str());
		tmpfile3.open(filename3l.c_str());                                                                                                                                                      
		tmpfile2.setf(ios::scientific,ios::floatfield);
		tmpfile2.precision(8);
		tmpfile3.setf(ios::scientific,ios::floatfield);
		tmpfile3.precision(8);


		r.lesHouchesAccordOutput(tmpfile2,"sugra",pars,sgnMu,tanb,0.,1,mGutGuess);
		s.lesHouchesAccordOutput(tmpfile3,"sugra",pars,sgnMu,tanb,0.,1,mGutGuess);
		tmpfile2.close();
		tmpfile3.close();

    }

  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}
