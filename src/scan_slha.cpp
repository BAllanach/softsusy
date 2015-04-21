
/** \file softpoint.cpp
    - Project:     SOFTSUSY 
   - Authors:     Ben Allanach, Markus Bernhardt 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 

   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: main calling program: command line interface. Reads Les
   - Houches files and command-line inputs and drives the calculation of a point
   - in parameter space.
*/ 

#include "softpoint.h"
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

double getCrossSection(MssmSoftsusy & r, char * fileName, double m0, double m12, double a0, double tanb) {
  double msq = (r.displayPhys().mu(1, 1) + r.displayPhys().mu(2, 1) + 
	 r.displayPhys().md(1, 1) + r.displayPhys().mu(2, 1)) / 4;
  double mg = r.displayPhys().mGluino;
  double mt = r.displayDataSet().displayPoleMt();

  char buff[500];
  char fn[500];
  sprintf(fn, "/home/bca20/code/prospino/output_%d_%d_%d_%d",int(m0),int(m12),int(a0),int(tanb));
  sprintf(buff, "cd /home/bca20/code/prospino; echo \"%10.4f %10.4f %8.2f %d 14000.\" | ./crosssec | grep -v NAN > output_%d_%d_%d_%d",msq,mg,mt,1,int(m0),int(m12),int(a0),int(tanb));
  cout << buff << endl;
  int err = system(buff);
  double xs = 0.;

  if (!err) {
  fstream fin(fn, ios::in); 
  fin >> xs;
  } else  cout << "CROSS SECTION ERROR\n" << xs;
  remove(fn); 
  
  return xs;
}

double doDarkMatter(DoubleVector & pars, double tanb, int sgnMu, 
		    char * fileName) {
  double m0 = pars(1), m12 = pars(2), a0 = pars(3);
  char oFile[500], buff[500];
  sprintf(oFile,"om_%d_%d_%d_%d_%d", int(m0), int(m12), int(a0), 
	  int(tanb), int(sgnMu));
  sprintf(buff,"/home/bca20/code/micromegas_3.6.9.2/MSSM/slha %s > %s",
	  fileName, oFile);
  int err = system(buff);
  double omega = 0.;
  if (!err) //throw("Problem in micromegas system call: \n");
    { fstream fin2(oFile, ios::in); fin2 >> omega; fin2.close(); }
  
  remove(oFile); 

  return omega;
}

void writeTable(MssmSoftsusy & twoLoop, MssmSoftsusy & oddLoop, MssmSoftsusy & twoLoopAs,  
		MssmSoftsusy & twoLoopMt, MssmSoftsusy & twoLoopMb, 
		MssmSoftsusy & threeLoop, double omega2, double omegaOdd, double omega2As, 
		double omega2Mt, double omega2Mb, double omega3, 
		double msqAv2, double msqAvOdd, double msqAv2As, double msqAv2Mt, 
		double msqAv2Mb, double msqAv3, double cs, double csOdd, double csAs, double csMt, double csMb, double cs3) {
  cout << "\\begin{table}\n\\begin{center}\n\\begin{tabular}{|c|ccccccc|}"
       << "\\hline\n  & $m_h$  & $m_{\\tilde g}$ & "
       << "$m_{{\\tilde q}}$ & $m_{\\chi_1^0}$  & $m_{\\chi_2^0}$ & "
       << "$m_{\\chi_3^0}$ & $m_{\\chi_4}^0$ \\\\ \\hline\n$Q$"
       << "               &  ";
  if (omega2 != omega2) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A\\\\\n"); else   printf("%5.1f & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f &%4.0f\\\\\n",
	 twoLoop.displayPhys().mh0(1), 
	 twoLoop.displayPhys().mGluino, 
	 msqAv2, 
	 fabs(twoLoop.displayPhys().mneut(1)), 
	 fabs(twoLoop.displayPhys().mneut(2)), 
	 fabs(twoLoop.displayPhys().mneut(3)), 
	 fabs(twoLoop.displayPhys().mneut(4))
	 );
  cout << "$\\Delta_3$ &  ";
  if (omegaOdd != omegaOdd) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A\\\\\n"); else  printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 oddLoop.displayPhys().mh0(1) - twoLoop.displayPhys().mh0(1), 
	 oddLoop.displayPhys().mGluino - twoLoop.displayPhys().mGluino, 
	 msqAvOdd - msqAv2, 
	 fabs(oddLoop.displayPhys().mneut(1)) - 
	 fabs(twoLoop.displayPhys().mneut(1)), 
	 fabs(oddLoop.displayPhys().mneut(2)) - 
	 fabs(twoLoop.displayPhys().mneut(2)), 
	 fabs(oddLoop.displayPhys().mneut(3)) - 
	 fabs(twoLoop.displayPhys().mneut(3)), 
	 fabs(oddLoop.displayPhys().mneut(4)) - 
	 fabs(twoLoop.displayPhys().mneut(4))
	 ); 
  cout << "$\\Delta \\alpha_s$  &  ";
  if (omega2As!= omega2As) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n"); else   printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 twoLoopAs.displayPhys().mh0(1) - twoLoop.displayPhys().mh0(1), 
	 twoLoopAs.displayPhys().mGluino - twoLoop.displayPhys().mGluino, 
	 msqAv2As - msqAv2, 
	 fabs(twoLoopAs.displayPhys().mneut(1)) - 
	 fabs(twoLoop.displayPhys().mneut(1)), 
	 fabs(twoLoopAs.displayPhys().mneut(2)) - 
	 fabs(twoLoop.displayPhys().mneut(2)), 
	 fabs(twoLoopAs.displayPhys().mneut(3)) - 
	 fabs(twoLoop.displayPhys().mneut(3)), 
	 fabs(twoLoopAs.displayPhys().mneut(4)) - 
	 fabs(twoLoop.displayPhys().mneut(4))
	 );
  cout << "$\\Delta m_t$      & ";
  if (omega2Mt != omega2Mt) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A  \\\\\n"); else   printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 twoLoopMt.displayPhys().mh0(1) - twoLoop.displayPhys().mh0(1), 
	 twoLoopMt.displayPhys().mGluino - twoLoop.displayPhys().mGluino, 
	 msqAv2Mt - msqAv2, 
	 fabs(twoLoopMt.displayPhys().mneut(1)) - 
	 fabs(twoLoop.displayPhys().mneut(1)), 
	 fabs(twoLoopMt.displayPhys().mneut(2)) - 
	 fabs(twoLoop.displayPhys().mneut(2)), 
	 fabs(twoLoopMt.displayPhys().mneut(3)) - 
	 fabs(twoLoop.displayPhys().mneut(3)), 
	 fabs(twoLoopMt.displayPhys().mneut(4)) - 
	 fabs(twoLoop.displayPhys().mneut(4))
	 );
  cout << "$\\Delta m_b, m_\\tau$& ";
  if (omega2Mb != omega2Mb) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A  \\\\\n"); else   printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 twoLoopMb.displayPhys().mh0(1) - twoLoop.displayPhys().mh0(1), 
	 twoLoopMb.displayPhys().mGluino - twoLoop.displayPhys().mGluino, 
	 msqAv2Mb - msqAv2, 
	 fabs(twoLoopMb.displayPhys().mneut(1)) - 
	 fabs(twoLoop.displayPhys().mneut(1)), 
	 fabs(twoLoopMb.displayPhys().mneut(2)) - 
	 fabs(twoLoop.displayPhys().mneut(2)), 
	 fabs(twoLoopMb.displayPhys().mneut(3)) - 
	 fabs(twoLoop.displayPhys().mneut(3)), 
	 fabs(twoLoopMb.displayPhys().mneut(4)) - 
	 fabs(twoLoop.displayPhys().mneut(4))
	 );
  cout << "$\\Delta$ All      & ";
  if (omega3 != omega3) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n"); else   printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 threeLoop.displayPhys().mh0(1) - twoLoop.displayPhys().mh0(1), 
	 threeLoop.displayPhys().mGluino - twoLoop.displayPhys().mGluino, 
	 msqAv3 - msqAv2, 
	 fabs(threeLoop.displayPhys().mneut(1)) - 
	 fabs(twoLoop.displayPhys().mneut(1)), 
	 fabs(threeLoop.displayPhys().mneut(2)) - 
	 fabs(twoLoop.displayPhys().mneut(2)), 
	 fabs(threeLoop.displayPhys().mneut(3)) - 
	 fabs(twoLoop.displayPhys().mneut(3)), 
	 fabs(threeLoop.displayPhys().mneut(4)) - 
	 fabs(twoLoop.displayPhys().mneut(4))
	 );
  cout << "\n%\n\\hline"
       << "& $m_{{\\tilde t}_2}$  & $m_{{\\tilde t}_1}$ &$m_{{\\tilde b}_1}$&"
       << "$m_{{\\tilde b}_2}$&$m_{{\\tilde \\tau}_2}$&$m_{{\\tilde \\tau}_1}$&"
       << "$m_{\\chi_1}^\\pm$ \\\\ \\hline\n"
       << "$Q$             & ";
  if (omega2 != omega2) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n"); else   printf("%4.0f & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f &%4.0f\\\\\n",
	 twoLoop.displayPhys().mu(1, 3), 
	 twoLoop.displayPhys().mu(2, 3), 
	 twoLoop.displayPhys().md(1, 3), 
	 twoLoop.displayPhys().md(2, 3), 
	 twoLoop.displayPhys().me(1, 3), 
	 twoLoop.displayPhys().me(2, 3), 
	 twoLoop.displayPhys().mch(1)
	 );
  cout << "$\\Delta_3$ & ";
  if (omegaOdd != omegaOdd) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n"); else printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 oddLoop.displayPhys().mu(1, 3) - twoLoop.displayPhys().mu(1, 3), 
	 oddLoop.displayPhys().mu(2, 3) - twoLoop.displayPhys().mu(2, 3), 
	 oddLoop.displayPhys().md(1, 3) - twoLoop.displayPhys().md(1, 3), 
	 oddLoop.displayPhys().md(2, 3) - twoLoop.displayPhys().md(2, 3), 
	 oddLoop.displayPhys().me(1, 3) - twoLoop.displayPhys().me(1, 3), 
	 oddLoop.displayPhys().me(2, 3) - twoLoop.displayPhys().me(2, 3), 
	 fabs(oddLoop.displayPhys().mch(1)) - 
	 fabs(twoLoop.displayPhys().mch(1))
	 );

  cout << "$\\Delta \\alpha_s$  & ";
  if (omega2As != omega2As) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n"); else   printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 twoLoopAs.displayPhys().mu(1, 3) - twoLoop.displayPhys().mu(1, 3), 
	 twoLoopAs.displayPhys().mu(2, 3) - twoLoop.displayPhys().mu(2, 3), 
	 twoLoopAs.displayPhys().md(1, 3) - twoLoop.displayPhys().md(1, 3), 
	 twoLoopAs.displayPhys().md(2, 3) - twoLoop.displayPhys().md(2, 3), 
	 twoLoopAs.displayPhys().me(1, 3) - twoLoop.displayPhys().me(1, 3), 
	 twoLoopAs.displayPhys().me(2, 3) - twoLoop.displayPhys().me(2, 3), 
	 fabs(twoLoopAs.displayPhys().mch(1)) - 
	 fabs(twoLoop.displayPhys().mch(1))
	 );
  cout << "$\\Delta m_t$      & ";
  if (omega2Mt != omega2Mt) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n"); else   printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 twoLoopMt.displayPhys().mu(1, 3) - twoLoop.displayPhys().mu(1, 3), 
	 twoLoopMt.displayPhys().mu(2, 3) - twoLoop.displayPhys().mu(2, 3), 
	 twoLoopMt.displayPhys().md(1, 3) - twoLoop.displayPhys().md(1, 3), 
	 twoLoopMt.displayPhys().md(2, 3) - twoLoop.displayPhys().md(2, 3), 
	 twoLoopMt.displayPhys().me(1, 3) - twoLoop.displayPhys().me(1, 3), 
	 twoLoopMt.displayPhys().me(2, 3) - twoLoop.displayPhys().me(2, 3), 
	 fabs(twoLoopMt.displayPhys().mch(1)) - 
	 fabs(twoLoop.displayPhys().mch(1))
	 );
  cout << "$\\Delta m_b, m_\\tau$& ";
  if (omega2Mb != omega2Mb) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n"); else   printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 twoLoopMb.displayPhys().mu(1, 3) - twoLoop.displayPhys().mu(1, 3), 
	 twoLoopMb.displayPhys().mu(2, 3) - twoLoop.displayPhys().mu(2, 3), 
	 twoLoopMb.displayPhys().md(1, 3) - twoLoop.displayPhys().md(1, 3), 
	 twoLoopMb.displayPhys().md(2, 3) - twoLoop.displayPhys().md(2, 3), 
	 twoLoopMb.displayPhys().me(1, 3) - twoLoop.displayPhys().me(1, 3), 
	 twoLoopMb.displayPhys().me(2, 3) - twoLoop.displayPhys().me(2, 3), 
	 fabs(twoLoopMb.displayPhys().mch(1)) - 
	 fabs(twoLoop.displayPhys().mch(1))
	 );
  cout << "$\\Delta$ All      & ";
  if (omega3 != omega3) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n"); else   printf("%+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f & %+5.1f &%+5.1f\\\\\n",
	 threeLoop.displayPhys().mu(1, 3) - twoLoop.displayPhys().mu(1, 3), 
	 threeLoop.displayPhys().mu(2, 3) - twoLoop.displayPhys().mu(2, 3), 
	 threeLoop.displayPhys().md(1, 3) - twoLoop.displayPhys().md(1, 3), 
	 threeLoop.displayPhys().md(2, 3) - twoLoop.displayPhys().md(2, 3), 
	 threeLoop.displayPhys().me(1, 3) - twoLoop.displayPhys().me(1, 3), 
	 threeLoop.displayPhys().me(2, 3) - twoLoop.displayPhys().me(2, 3), 
	 fabs(threeLoop.displayPhys().mch(1)) - 
	 fabs(twoLoop.displayPhys().mch(1))
	 );
  cout << "\n%\n\\hline";
  cout << "      &   $g_3(M_{SUSY})$ & $Y_t(M_{SUSY})$ & "
       << " $Y_b(M_{SUSY})$ & $Y_\\tau(M_{SUSY})$  & $\\mu(M_{SUSY})$"
       << "    & $\\Omega_{CDM} h^2$ & $\\sigma_{SUSY}^{TOT}$\\\\ \\hline\n"
       << " $Q$                   & ";

  if (omega2 != omega2) printf("N/A & N/A & N/A & N/A & N/A & N/A &\\\\\n"); else printf("%5.3f & %5.3f & %5.3f & %5.3f & %4.0f & %5.2f & %5.0f\\\\\n",
	 twoLoop.displayGaugeCoupling(3), 
	 twoLoop.displayYukawaElement(YU, 3, 3), 
	 twoLoop.displayYukawaElement(YD, 3, 3), 
	 twoLoop.displayYukawaElement(YE, 3, 3), 
	 twoLoop.displaySusyMu(),
											 omega2, cs
	 );
  cout << "$\\Delta_3$     &  ";
  if (omegaOdd != omegaOdd) printf("N/A & N/A & N/A & N/A & N/A & N/A &\\\\\n"); else printf("%+5.3f & %+5.3f & %+5.3f & %+5.3f & %+4.0f & %+5.2f & %+5.0f\\\\\n",
	 oddLoop.displayGaugeCoupling(3) - twoLoop.displayGaugeCoupling(3), 
	 oddLoop.displayYukawaElement(YU, 3, 3) - 
	 twoLoop.displayYukawaElement(YU, 3, 3), 
	 oddLoop.displayYukawaElement(YD, 3, 3) -
	 twoLoop.displayYukawaElement(YD, 3, 3), 
	 oddLoop.displayYukawaElement(YE, 3, 3) -
	 twoLoop.displayYukawaElement(YE, 3, 3), 
	 oddLoop.displaySusyMu() - twoLoop.displaySusyMu(),
											     omegaOdd - omega2,
											     csOdd - cs
	 );

  cout << "$\\Delta \\alpha_s$  & ";
  if (omega2As != omega2As) printf("N/A & N/A & N/A & N/A & N/A & N/A &\\\\\n"); else printf("%+5.3f & %+5.3f & %+5.3f & %+5.3f & %+4.0f & %+5.2f & %+5.0f\\\\\n",
	 twoLoopAs.displayGaugeCoupling(3) - twoLoop.displayGaugeCoupling(3), 
	 twoLoopAs.displayYukawaElement(YU, 3, 3) - 
	 twoLoop.displayYukawaElement(YU, 3, 3), 
	 twoLoopAs.displayYukawaElement(YD, 3, 3) -
	 twoLoop.displayYukawaElement(YD, 3, 3), 
	 twoLoopAs.displayYukawaElement(YE, 3, 3) -
	 twoLoop.displayYukawaElement(YE, 3, 3), 
	 twoLoopAs.displaySusyMu() - twoLoop.displaySusyMu(),
											     omega2As - omega2,
											     csAs - cs
	 );
  cout << "$\\Delta m_t$      & ";
    if (omega2Mt != omega2Mt) printf("N/A & N/A & N/A & N/A & N/A & N/A& N/A\\\\\n"); else printf("%+5.3f & %+5.3f & %+5.3f & %+5.3f & %+4.0f & %+5.2f & %+5.0f\\\\\n",
	 twoLoopMt.displayGaugeCoupling(3) - twoLoop.displayGaugeCoupling(3), 
	 twoLoopMt.displayYukawaElement(YU, 3, 3) - 
	 twoLoop.displayYukawaElement(YU, 3, 3), 
	 twoLoopMt.displayYukawaElement(YD, 3, 3) -
	 twoLoop.displayYukawaElement(YD, 3, 3), 
	 twoLoopMt.displayYukawaElement(YE, 3, 3) -
	 twoLoop.displayYukawaElement(YE, 3, 3), 
	 twoLoopMt.displaySusyMu() - twoLoop.displaySusyMu(),
											       omega2Mt - omega2, csMt - cs
	 );
    cout << "$\\Delta m_b, m_\\tau$& ";
    if (omega2Mb != omega2Mb) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A\\\\\n"); else printf("%+5.3f & %+5.3f & %+5.3f & %+5.3f & %+4.0f & %+5.2f & %+5.0f\\\\\n",
	 twoLoopMb.displayGaugeCoupling(3) - twoLoop.displayGaugeCoupling(3), 
	 twoLoopMb.displayYukawaElement(YU, 3, 3) - 
	 twoLoop.displayYukawaElement(YU, 3, 3), 
	 twoLoopMb.displayYukawaElement(YD, 3, 3) -
	 twoLoop.displayYukawaElement(YD, 3, 3), 
	 twoLoopMb.displayYukawaElement(YE, 3, 3) -
	 twoLoop.displayYukawaElement(YE, 3, 3), 
	 twoLoopMb.displaySusyMu() - twoLoop.displaySusyMu(),
												omega2Mb - omega2, csMb - cs
	 );
    cout << "$\\Delta$ All      & ";
    //    if (omega3 != omega3) printf("N/A & N/A & N/A & N/A & N/A & N/A &
    //    \\\\\n"); else 
    printf("%+5.3f & %+5.3f & %+5.3f & %+5.3f & %+4.0f & %+5.2f &%+5.0f\\\\\n",
	 threeLoop.displayGaugeCoupling(3) - twoLoop.displayGaugeCoupling(3), 
	 threeLoop.displayYukawaElement(YU, 3, 3) - 
	 twoLoop.displayYukawaElement(YU, 3, 3), 
	 threeLoop.displayYukawaElement(YD, 3, 3) -
	 twoLoop.displayYukawaElement(YD, 3, 3), 
	 threeLoop.displayYukawaElement(YE, 3, 3) -
	 twoLoop.displayYukawaElement(YE, 3, 3), 
	 threeLoop.displaySusyMu() - twoLoop.displaySusyMu(),
	 omega3 - omega2, 
         cs3-cs
	);
  cout << "\\hline\n\\end{tabular}\n\\end{center}\n";
}

/// Returns the object along with omega. Oneset should already be fixed at MZ
void getMssmAndOmega(MssmSoftsusy & r, DoubleVector & pars, const double tanb, 
		     const int sgnMu, const QedQcd & oneset, 
		     double mGutGuess, bool uni, double & omega, 
		     double & msqAv,
		     void (*boundaryCondition)(MssmSoftsusy &, 
					       const DoubleVector &), 
		     bool ewsbBcScale, double & cs) {
  double m0 = pars(1), m12 = pars(2), a0 = pars(3);
  r.fixedPointIteration(boundaryCondition, mGutGuess, pars, sgnMu, tanb, 
			oneset, uni, ewsbBcScale); 

  r.setData(oneset);

  /// Produces SLHA output file
  char fileName[500]; 
  sprintf(fileName,"lesHout");
  produceSLHAfile(r, fileName, sgnMu, tanb, pars);

  cs = getCrossSection(r, fileName, m0, m12, a0, tanb);

  if (!r.displayProblem().test()) 
    omega = doDarkMatter(pars, tanb, sgnMu, fileName);

  msqAv = (r.displayPhys().mu(2, 1) + 
		     r.displayPhys().mu(1, 1) +
		     r.displayPhys().md(2, 1) + 
		     r.displayPhys().md(1, 1)) * 
    0.25;
  
  remove(fileName); 
  return;
}

// Returns a string with all characters in upper case: very handy
string ToUpper(const string & s) {
        string result;
        unsigned int index;
        for (index = 0; index < s.length(); index++) {
	  char a = s[index];
	  a = toupper(a);
	  result = result + a;
        }
	
	return result;
    }

void errorCall() {
  ostringstream ii;
  ii << "\n\nSOFTSUSY" << SOFTSUSY_VERSION 
     << " called with incorrect arguments. Need to put either:\n";
  ii << "./softpoint.x leshouches < lesHouchesInput\n for SLHA/SLAH2 input, or\n";
  ii << "./softpoint.x sugra [SUGRA parameters] [other options]\n";
  ii << "./softpoint.x amsb [mAMSB parameters] [other options]\n";
  ii << "./softpoint.x gmsb [mGMSB parameters] [other options]\n";
  ii << "./softpoint.x nmssm sugra [NMSSM flags] [NMSSM parameters] [other options]\n\n";
  ii << "[other options]: --mbmb=<value> --mt=<value> --alpha_s=<value> --QEWSB=<value>\n";
  ii << "--alpha_inverse=<value> --tanBeta=<value> --sgnMu=<value>\n";
#ifdef COMPILE_TWO_LOOP_GAUGE_YUKAWA
  if (USE_TWO_LOOP_GAUGE_YUKAWA) ii << "--disable-full_susy_threshold disables the 2-loop SUSY threshold corrections to third generation Yukawa couplings and g3.\n";
#endif //COMPILE_TWO_LOOP_GAUGE_YUKAWA
#ifdef COMPILE_THREE_LOOP_RGE
  if (USE_THREE_LOOP_RGE) ii << "--disable-three_loop disables 3-loop corrections RGEs\n";
#endif //COMPILE_THREE_LOOP_RGE
  ii << "--mgut=unified sets the scale at which SUSY breaking terms are set to the GUT\n";
  ii << "scale where g1=g2. --mgut=<value> sets it to a fixed scale, ";
  ii << "whereas --mgut=msusy\nsets it to MSUSY\n\n";
  ii << "If you want the R-parity violating MSSM calculation, set any of the following:\n";
  ii << "--lambda <i> <j> <k> <coupling>, the word lambda replaceable with lambdaP\nor lambdaPP for LLE, LQD, UDD coupling, respectively.\n\n";
  ii << "[SUGRA parameters]: --m0=<value> --m12=<value> --a0=<value>\n";
  ii << "[mAMSB parameters]: --m0=<value> --m32=<value>\n";
  ii << "[mGMSB parameters]: --n5=<value> --mMess=<value> --LAMBDA=<value> --cgrav=<value>\n\n";
  ii << "Bracketed entries are numerical values, in units of GeV if they are massive.\n";
  ii << "Warning: entries left unspecified will be assumed to be zero for SUSY breaking\nterms, unified (for mgut) or at their default central values for Standard Model parameters\n";
  ii << "\n"
     "[NMSSM flags]:\n"
     "  --lambdaAtMsusy   input lambda at renormalization scale Q = Msusy\n"
     "\n"
     "[NMSSM parameters]:\n"
     "  --m0= , --m12= , --a0= , --tanBeta= , --mHd2= , --mHu2= ,\n"
     "  --mu= , --m3SqrOverCosBetaSinBeta= , --lambda= , --kappa= ,\n"
     "  --Alambda= , --Akappa= , --lambdaS= , --xiF= , --xiS= ,\n"
     "  --muPrime= , --mPrimeS2= , --mS2=\n"
     "\n"
     "  Unset NMSSM parameters are assumed to be zero\n"
     "\n"
     "NMSSM example:\n"
     "  ./softpoint.x nmssm sugra --m0=125 --m12=200 --tanBeta=10 --a0=-300 \\\n"
     "     --lambda=0.1 --lambdaAtMsusy\n";
  throw ii.str();
}

int main(int argc, char *argv[]) {

  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  double lambdaW = 0., aCkm = 0., rhobar = 0., etabar = 0.;
  NMSSM_input nmssm_input; // NMSSM input parameters

  bool flavourViolation = false;

  int numPoints = 1;

  double qMax = 0.;

  // Sets format of output: 4 decimal places
  outputCharacteristics(6);

  void (*boundaryCondition)(MssmSoftsusy &, const DoubleVector &)=sugraBcs;
  void (*nmssmBoundaryCondition)(NmssmSoftsusy&, const DoubleVector&) = NmssmMsugraBcs;

  QedQcd oneset;
  MssmSoftsusy m; FlavourMssmSoftsusy k;
  NmssmSoftsusy nmssm;
  k.setInitialData(oneset);
  MssmSoftsusy * r = &m; 
  RpvNeutrino kw; bool RPVflag = false;
  enum Model_t { MSSM, NMSSM } susy_model = MSSM; // susy model (MODSEL entry 3)
  softsusy::GUTlambda = true;
  softsusy::GUTkappa = true;
  softsusy::GUTmuPrime = true;
  softsusy::GUTxiF = true;
  softsusy::GUTsVev = true;

  try {
  if (argc !=1 && strcmp(argv[1],"leshouches") != 0) {
    cerr << "SOFTSUSY" << SOFTSUSY_VERSION << endl;
    if (!strcmp(argv[1], "-v") || !strcmp(argv[1], "--version")) exit(0);
    cerr << "B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331,";
    cerr << " hep-ph/0104145\n";
    cerr << "For RPV aspects, B.C. Allanach and M.A. Bernhardt, Comput. "
	 << "Phys. Commun. 181 (2010) 232, arXiv:0903.1805.\n\n";
    cerr << "Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
	 << TOLERANCE << endl;
    cerr << "G_F=" << GMU << " GeV^2" << endl;
  }
  
  double mgutGuess = 2.0e16, tanb = 0.;
  int sgnMu = 1;
  bool gaugeUnification = true, ewsbBCscale = false;
  double desiredMh = 0.;

  // If there are no arguments, give error message,
  // or if none of the options are called, then go to error message
  if (argc == 1 || ( (strcmp(argv[1], "sugra") && 
		      strcmp(argv[1], "amsb") &&
		      strcmp(argv[1], "gmsb") && 
		      strcmp(argv[1], "runto") && 
		      strcmp(argv[1], "leshouches") && 
                      strcmp(argv[1], "nmssm") &&
		      strcmp(argv[1], "-v") &&
		      strcmp(argv[1], "--version"))))
    errorCall();
  
  DoubleVector pars(3); 
  
  const char* modelIdent = "";
  
  /// Non model specific options
  if (strcmp(argv[1], "leshouches")) {
      for (int i = 2; i < argc; i++) {
	if (starts_with(argv[i],"--mbmb=")) 
	  oneset.setMass(mBottom, get_value(argv[i], "--mbmb="));
	else if (starts_with(argv[i],"--mt=")) 
	  oneset.setPoleMt(get_value(argv[i], "--mt="));
	else if (starts_with(argv[i],"--alpha_s="))
	  oneset.setAlpha(ALPHAS, get_value(argv[i], "--alpha_s="));      
	else if (starts_with(argv[i],"--alpha_inverse="))
	  oneset.setAlpha(ALPHA, 1.0 / get_value(argv[i],"--alpha_inverse="));
	else if (starts_with(argv[i],"--RPV")) 
	  RPVflag = true;
	else if (starts_with(argv[i], "--tanBeta=")) 
	  tanb = get_value(argv[i], "--tanBeta=");
	else if (starts_with(argv[i], "--sgnMu=")) 
	  sgnMu = get_valuei(argv[i], "--sgnMu=");
	else if (starts_with(argv[i], "--mgut=")) 
	  mgutGuess = mgutCheck(argv[i], gaugeUnification, ewsbBCscale); 
#ifdef COMPILE_TWO_LOOP_GAUGE_YUKAWA
	else if (starts_with(argv[i], "--disable-full_susy_threshold"))
	  USE_TWO_LOOP_GAUGE_YUKAWA = false;
#endif
#ifdef COMPILE_THREE_LOOP_RGE
	else if (starts_with(argv[i], "--disable-three_loop_rge"))
	  USE_THREE_LOOP_RGE = false;
#endif
	else if (starts_with(argv[i], "--QEWSB=")) 
	  QEWSB = get_value(argv[i], "--QEWSB=");
      }
      if (tanb < 1.5 || tanb > 70.) {
	ostringstream ii; 
	ii << "tanBeta=" << tanb 
	   << " in SUGRA input. The point will not yield a sensible answer\n";
	throw ii.str();
      }
      
      /// Pass through to see if there are any RPV options
      if (strcmp(argv[1], "nmssm")) {
        for (int i = 2; i < argc; i++) {
          if (starts_with(argv[i], "--lambda")) {
            if (i + 4 >= argc) {
              throw "ERROR: three indices and one value need to be provided"
                " after --lambda or --lambdaP or --lambdaPP\n";
            }
            RPVflag = true;
            int ii= int(atof(argv[i+1]));
            int j = int(atof(argv[i+2]));
            int k = int(atof(argv[i+3]));
            double d = atof(argv[i+4]);
            if (starts_with(argv[i], "--lambdaPP"))
              kw.setLambda(LU, ii, j, k, d);
            else if (starts_with(argv[i], "--lambdaP"))
              kw.setLambda(LD, k, ii, j, d);
            else if (starts_with(argv[i], "--lambda"))
              kw.setLambda(LE, k, ii, j, d);
          }
          if (starts_with(argv[i], "--kappa")) {
            if (i + 2 >= argc) {
              throw "ERROR: one index and one value need to be provided"
                " after --kappa\n";
            }
            int ii = int(atof(argv[i+1]));
            double d = atof(argv[i+2]);
            kw.setKappa(ii, d);
            RPVflag = true;
          }
        }
      }
      
      /// Model specific options
      if (!strcmp(argv[1], "sugra")) {
	cout << "# SOFTSUSY SUGRA calculation" << endl;
	boundaryCondition = &sugraBcs;
	modelIdent = "sugra";
	double m0 = 0., m12 = 0., a0 = 0.;
	for (int i = 2; i < argc; i++) {
	  if (starts_with(argv[i], "--m0=")) m0 = get_value(argv[i], "--m0=");
	  else if (starts_with(argv[i], "--m12=")) 
	    m12 = get_value(argv[i], "--m12=");
	  else if (starts_with(argv[i], "--a0="))
	    a0  = get_value(argv[i], "--a0=");
	}
	pars(1) = m0; pars(2) = m12; pars(3) = a0;      
	if (m12 < MZ) {
	  ostringstream ii; 
	  ii << "m12=" << m12 
	     << " in SUGRA input. The point will not yield a sensible answer\n";
	  throw ii.str();
	}
	r = &m;
      }
      
      if (!strcmp(argv[1], "amsb")) {
	cout << "# SOFTSUSY mAMSB calculation" << endl;
	boundaryCondition = &amsbBcs;
	modelIdent = "amsb";
	double m0 = 0., m32 = 0.;
	for (int i = 2; i < argc; i++) {
	  if (starts_with(argv[i], "--m0=")) 
	    m0 = get_value(argv[i], "--m0=");
	  else if (starts_with(argv[i], "--m32=")) 
	    m32 = get_value(argv[i], "--m32=");
	}
	pars(1) = m0; pars(2) = m32; 
	if (m32 < 1.0e3) {
	  ostringstream ii; 
	  ii << "m32=" << m32 
	     << " in SUGRA input (too low). The point will not yield a sensible answer\n";
	  throw ii.str();
	}
	r = &m;
      }
      
      if (!strcmp(argv[1], "gmsb")) {
	cout << "# SOFTSUSY mGMSB calculation" << endl;
	boundaryCondition = &gmsbBcs;
	modelIdent = "gmsb";
	double n5 = 0., mMess = 0., LAMBDA = 0., cgrav = 1.;
	for (int i = 2; i < argc; i++) {
	  if (starts_with(argv[i], "--n5=")) n5 = get_value(argv[i], "--n5=");
	  else if (starts_with(argv[i], "--mMess=")) { 
	    gaugeUnification = false; 
	    mMess = get_value(argv[i], "--mMess=");
	    mgutGuess = mMess;
	  }
	  else if (starts_with(argv[i], "--LAMBDA=")) 
	    LAMBDA = get_value(argv[i], "--LAMBDA=");
	  else if (starts_with(argv[i], "--cgrav=")) 
	    cgrav = get_value(argv[i], "--cgrav=");
	}
	
	pars.setEnd(4);
	pars(1) = n5; pars(2) = mMess; pars(3) = LAMBDA; pars(4) = cgrav;
	if (mMess < 1.0e3) {
	  ostringstream ii; 
	  ii << " mMess=" << mMess
	     << " in SUGRA input (too low). The point will not yield a sensible answer\n";
	  throw ii.str();
	}
	
	r = &m;
	if (LAMBDA > mMess) {
	  ostringstream ii;
	  ii << "Input LAMBDA=" << LAMBDA << " should be less than mMess="
	     << mMess << endl;
	  throw ii.str();
	}
	if (cgrav > 1.0) {
	  ostringstream ii;
	  ii << "Input cgrav=" << cgrav << " a real number bigger than or "
	     << " equal to 1 (you can use 1 as a default value).\n";
	  throw ii.str();
	}
      }
    }
    if (!strcmp(argv[1], "nmssm")) {
      susy_model = NMSSM;
      NMSSM_command_line_parser nmssm_parser(&nmssm_input);
      nmssm_parser.parse(argc, argv);
      modelIdent = nmssm_parser.get_modelIdent();
      pars = nmssm_parser.get_pars();
    }
    
  if (!strcmp(argv[1], "leshouches")) {
    cout << "In leshouches portion\n";

	/// SLHA option "leshouches" used.
	outputCharacteristics(8);
	if (argc == 2) {
	  string line, block;
	  int model;
	  
	  while (getline(cin,line)) {
	    //	  mgutGuess = mgutCheck("unified", gaugeUnification);
	    
	    //	cout << line << endl;
	    istringstream input(line); 
	    string word1, word2;
	    input >> word1;
	    
	    if (word1.find("#") == string::npos &&
                !contains_only_whitespace(word1)) {
	      // read in another word if there's no comment
	      input >> word2; 
	      
	      if (ToUpper(word1) == "BLOCK")  { 
		block = ToUpper(word2);
		
	      } else { // ought to be data
		istringstream kk(line);
		if (block == "MODSEL") {
		  int i; kk >> i; 
		  
		  switch(i) {
		  case 1: kk >> model; 
		    switch(model) {
		    case 0: boundaryCondition = &extendedSugraBcs;
		      modelIdent = "nonUniversal"; r=&m;
		      break;
		    case 1: 
		      if (!flavourViolation) {
			pars.setEnd(3); 
			boundaryCondition = &sugraBcs; 
		      }
		      modelIdent = "sugra";
		      break;
		    case 2: 
		      if (!flavourViolation) {
			boundaryCondition = &gmsbBcs; 
			pars.setEnd(4); 
		      } 
		      modelIdent = "gmsb";
		      break;
		    case 3: 		    
		      boundaryCondition = &amsbBcs; 
		      pars.setEnd(2); 
		      modelIdent = "amsb";
		      break;
		    case 4:
		    boundaryCondition = &splitGmsb;
		    pars.setEnd(7); sgnMu = 0; 
		    modelIdent = "splitgmsb";
		    break;
		    default: 
		      ostringstream ii;
		      ii << "SOFTSUSY" << SOFTSUSY_VERSION 
			 << " cannot yet do model " 
			 << model << ": terminal error\n";
		      throw ii.str();
		    }
		    break;
                    // reading entry 3: susy model (MSSM, NMSSM, ...)
                  case 3: { int i; kk >> i;
                      switch(i) {
                      case 0: susy_model = MSSM; // default
                         break;
                      case 1: susy_model = NMSSM;
                         if (flavourViolation) {
                            flavourViolation = false;
                            cout << "# Warning: flavour violation is currtently"
                               " not supported in the NMSSM\n";
                         }
                         break;
                      default:
                         ostringstream ii;
                         ii << "MODSEL 3 choosing silly model switch\n"
                            << "(" << i << ") not a valid switch" << endl;
                         throw ii.str();
                      }
                    }
                    break;
		  case 4: int i; kk >> i;
		    switch(i) {
		    case 0: RPVflag = false;
		      break;
		    case 1: RPVflag = true;
		      break;
		    default:
		      ostringstream ii;
		      ii << "MODSEL 4 choosing silly RPV switch\n"
			 << "(" << i << ") not a valid switch" << endl;
		      throw ii.str();
		    }
		    break;
		  case 6: int j; kk >> j;
                    switch(j) {
                    case 0: flavourViolation = false; break;
                    default:
                       if (susy_model == NMSSM) {
                          flavourViolation = false;
                          cout << "# Warning: flavour violation is currtently"
                             " not supported in the NMSSM\n";
                       } else {
                          r = &k; flavourViolation = true;
                          if (boundaryCondition != & amsbBcs) {
                             pars.setEnd(64); boundaryCondition = &flavourBcs;
                          }
                       }
                    }
                    break;
		  case 11: kk >> numPoints;
		    if (numPoints < 1) {
		      ostringstream ii;
		      ii << "MODSEL 11 selecting silly number of points"
			 << "(" << numPoints << ") to output" << endl;
		      throw ii.str();
		    }
		    break;
		  case 12: double d; kk >> d;
		    if (d < MZ) {
		      ostringstream ii;
		      ii << "MODSEL 12 selecting silly scale Qmax"
			 << "(" << d << ") < MZ to output" << endl;
		      throw ii.str();
		    }
		    qMax = d; break;
		  default:
		    cout << "# WARNING: don't understand first integer " 
			 << word1 << " " << word2 << " in block " << block
			 << ": ignoring it\n";
		    break;
		  }
		}
		else if (block == "MINPAR") {
		  int i; double d; kk >> i >> d; 
		  switch (i) {
		  case 3: tanb = d;
                    nmssm_input.set(NMSSM_input::tanBeta, d);
                    break;
		  case 4: sgnMu = int(d); break;
		  default: 
		    switch(model) {
		    case 0:
		      // SUGRA inputs to fill out the pheno MSSM case
		      switch(i) {
		      case 1: pars(1) = d; break;
		      case 2: pars(2) = d; break;
		      case 5: pars(3) = d; break;
		      default: 
			ostringstream ii;
			ii << "Didn't understand pheno MSSM input " << i << endl;
			break;
		      } break;
		    case 1: // SUGRA inputs
		      switch(i) {
		      case 1: 
			if (flavourViolation) { pars.setEnd(77);
			  double m0 = sqr(d);
			  pars(4) = m0; pars(7) = m0; pars(9) = m0;
			  pars(10) = m0; pars(13) = m0; pars(15) = m0; 
			  pars(16) = m0; pars(19) = m0; pars(21) = m0; 
			  pars(22) = m0; pars(25) = m0; pars(27) = m0; 
			  pars(28) = m0; pars(31) = m0; pars(33) = m0; 
			  pars(63) = m0; pars(64) = m0;
			} else pars(1) = d; 
			break;
		      case 2: 
			if (flavourViolation) {
			  pars(1) = d; pars(2) = d; pars(3) = d;
			} else pars(2) = d; 
			break;
		      case 5: 
			if (flavourViolation) {
			  pars.setEnd(77);
			  pars(62) = d;
			  
			} else pars(3) = d; 
			break;
		      default: 
			ostringstream ii;
			ii << "Didn't understand SUGRA input " << i << endl;
			break;
		      } break;
		    case 2: // GMSB inputs
		      switch(i) {
		      case 1: pars(3) = d; break;
		    case 2: pars(2) = d; mgutGuess = d;
		      gaugeUnification = false; break;
		      case 5: pars(1) = d; break;
		      case 6: pars(4) = d; break;
		      default: 
			ostringstream ii;
			ii << "Didn't understand GMSB input " << i << endl;
			break;
		      } break;
		    case 3: ///< AMSB inputs
		      switch(i) {
		      case 1: pars(2) = d; break;
		      case 2: pars(1) = d; break;
		      default: 
			ostringstream ii;
		      ii << "Didn't understand AMSB input " << i << endl;
		      break;
		      } break;
		    case 4: ///< split GMSB inputs 
		      switch(i) {
		      case 1: pars(2) = d; break;
		      case 2: pars(3) = d; break;
		      case 5: pars(1) = d; break;
		      case 6: pars(7) = d; break;
		      case 7: pars(4) = d; mgutGuess = d; 
			gaugeUnification = false; break;
		      case 8: pars(5) = d; break;
		      case 9: pars(6) = d; m.useAlternativeEwsb(); 
			kw.useAlternativeEwsb();
			break;
		      case 10: desiredMh = d; break;
		      default: 
			ostringstream ii;
			ii << "Didn't understand GMSB input " << i << endl;
			break;
		      } break;
		    default: 
		      ostringstream ii;
		    ii << "Didn't understand model input " << model << endl;
		    break;
		    }
		    break;
		  }
		}
		// Adding non-minimal options. 
		else if (block == "EXTPAR") {
                  int i; double d; kk >> i >> d;

                  // read extra NMSSM input parameters from EXTPAR
                  // (skipping NMSSM parameters if the MSSM was selected)
                  if (susy_model == MSSM) {
                     switch (i) {
                     case 61:
                     case 62:
                     case 63:
                     case 64:
                     case 65:
                     case 66:
                     case 67:
                     case 68:
                     case 69:
                     case 70:
                        cout << "# Warning: NMSSM parameter EXTPAR " << i
                             << " given but MSSM chosen -- ignoring it.\n";
                        continue;
                     }
                  } else if (susy_model == NMSSM) {
                     // read NMSSM susy parameters only and continue
                     switch (i) {
                     case 23: nmssm_input.set(NMSSM_input::mu     , d); 
		       continue;
                     case 61: nmssm_input.set(NMSSM_input::lambda , d); 
		       continue;
                     case 62: nmssm_input.set(NMSSM_input::kappa  , d); 
		       continue;
                     case 65: nmssm_input.set(NMSSM_input::lambdaS, d); 
		       continue;
                     case 66: nmssm_input.set(NMSSM_input::xiF    , d); 
		       continue;
                     case 68: nmssm_input.set(NMSSM_input::muPrime, d);
		       continue;
                     }
                  }
		  /// First, we want to convert our input to EXTPAR if we have
		  /// mSUGRA already
		  if (!strcmp(modelIdent, "sugra")) {
		    modelIdent = "nonUniversal";
		    if (!flavourViolation) {
		      /// We assume mSUGRA BCs with no flavour violation
		      r=&m; 
		      boundaryCondition = &extendedSugraBcs;
		      double m0 = pars(1), m12 = pars(2), a0 = pars(3);
		      pars.setEnd(49);
		      int i; for (i=1; i<=3; i++) pars(i) = m12;
		      for (i=11; i<=13; i++) pars(i) = a0;
		      pars(21) = m0*m0; pars(22) = m0*m0;
		      for (i=31; i<=36; i++) pars(i) = m0;		    
		      for (i=41; i<=49; i++) pars(i) = m0;		    
		      kw.setNumRpcBcs(50); 
                      if (susy_model == NMSSM) {
                         pars.setEnd(56);
                         pars(50) = a0; // Alambda
                         pars(51) = a0; // Akappa
                         pars(52) = 0.; // mS'^2
                         pars(53) = m0*m0; // mS^2
                         pars(54) = 0.; // mu
                         pars(55) = 0.; // Bmu
                         pars(56) = 0.; // xiS
                      }
  		    } else {
		      /// This is flavour violation with EXTPAR: mSUGRA BCs
		      /// with flavour violation
		      boundaryCondition = &flavourBcs;		    
		      if (pars.displayEnd() == 3) {
			double m0 = pars(1), m12 = pars(2), a0 = pars(3);
			double msq = m0 * m0;
			pars.setEnd(77);
			int i; for (i=1; i<=3; i++) pars(i) = m12;
			/// Fill in scalar mass squareds
			for (i=1; i<=5; i++) {
			  int num = (i-1) * 6 + 4;
			  pars(num)  = msq; 
			  pars(num+3)  = msq; 
			  pars(num+5)  = msq;
			}
			
			pars(62) = a0;
			pars(63) = msq; pars(64) = msq;
		      }
		      kw.setNumRpcBcs(65);
		    }
		  }
		  
		  if (!strcmp(modelIdent, "nonUniversal")) {
		    /// First, put parameters that depend not on
		    /// flavoured/unflavoured input
		    if (i == 0) { 
		      mgutGuess = d;
		      gaugeUnification = false;
		      // setting Minput=-1 should yield MSSM BCs at MSUSY
		      if (fabs(d + 1.0) < EPSTOL) {
			mgutGuess = 1.0e3;
			ewsbBCscale = true;
			QEWSB = 1.0;
			if (gaugeUnification) 
			  cout << "# Gauge unification ignored since pheno MSSM"
			       << " assumes BC set at QEWSB\n"; 
			gaugeUnification = false;
		      }
		    }
		    else if (i == 25) {
		      tanb = d;
		      if (!flavourViolation && pars.displayEnd() != 49) 
			pars.setEnd(49);
		      pars(i) = d;
		      r->setSetTbAtMX(true);
		    } 
		    else if (i == 23 || i == 26) {
		      r->useAlternativeEwsb(); 
		      if (i == 23) {
                         r->setMuCond(d); r->setSusyMu(d);
                         if (susy_model == NMSSM) {
                           if (pars.displayEnd() < 56) pars.setEnd(56);
                           nmssm_input.set(NMSSM_input::mu, d);
                           pars(54) = d;
                         }
                      }
		      if (i == 26) r->setMaCond(d); 
		    }
		    else if (!flavourViolation) {
		      if ((i > 0 && i <=  3) || (i >= 11 && i <= 13) || 
			  (i >= 21 && i <= 23) || (i == 26 || i == 25) 
			  || (i >= 31 && i <= 36) || 
			  (i >= 41 && i <= 49)) {
			if (pars.displayEnd() < 49) pars.setEnd(49);
			pars(i) = d;
                        if (susy_model == NMSSM) {
                           if (pars.displayEnd() < 56) pars.setEnd(56);
                           switch (i) {
                           case 21: nmssm_input.set(NMSSM_input::mHd2, d); break;
                           case 22: nmssm_input.set(NMSSM_input::mHu2, d); break;
                           case 23:
                             nmssm_input.set(NMSSM_input::mu, d);
                             pars(54) = d;
                             break;
                           }
                        }
  		      } else if ((61 <= i && i <= 70) || i == 24) {
                        switch (i) {
                        case 24:
                           nmssm_input.set(NMSSM_input::BmuOverCosBetaSinBeta, d);
                           pars(55) = d;
                           break;
                        case 63:
                           nmssm_input.set(NMSSM_input::Alambda, d);
                           pars(50) = d;
                           break;
                        case 64:
                           nmssm_input.set(NMSSM_input::Akappa, d);
                           pars(51) = d;
                           break;
                        case 67:
                           nmssm_input.set(NMSSM_input::xiS, d);
                           pars(56) = d;
                           break;
                        case 69:
                           nmssm_input.set(NMSSM_input::mPrimeS2, d);
                           // setting pars(52) = B' = mS'^2 / mu'
                           if (nmssm_input.is_set(NMSSM_input::muPrime)) {
                              const double muPrime = nmssm_input.get(NMSSM_input::muPrime);
                              if (!close(muPrime, 0.0, EPSTOL))
                                 pars(52) = d / muPrime;
                           }
                           break;
                        case 70:
                           nmssm_input.set(NMSSM_input::mS2, d);
                           pars(53) = d;
                           break;
                        }
                      } else {
                        cout << "WARNING: did not understand parameter " 
                             << i << " in non-flavoured EXTPAR inputs\n";
                      }
		    } else {
		      /// Have to translate the numbers from SLHA to your
		      /// convention with flavour violation
		      if ((i > 0 && i < 4)) pars(i) = d; 
		      else if (i == 31) pars(22) = sqr(d);
		      else if (i == 32) pars(25) = sqr(d);
		      else if (i == 33) pars(27) = sqr(d);
		      else if (i == 34) pars(28) = sqr(d);
		      else if (i == 35) pars(31) = sqr(d);
		      else if (i == 36) pars(33) = sqr(d);
		      else if (i == 41) pars(4) = sqr(d);
		      else if (i == 42) pars(7) = sqr(d);
		      else if (i == 43) pars(9) = sqr(d);
		      else if (i == 44) pars(10) = sqr(d);
		      else if (i == 45) pars(13) = sqr(d);
		      else if (i == 46) pars(15) = sqr(d);
		      else if (i == 47) pars(16) = sqr(d);
		      else if (i == 48) pars(19) = sqr(d);
		      else if (i == 49) pars(21) = sqr(d);
		      else if (i == 21) pars(63) = d;
		      else if (i == 22) pars(64) = d;
		      else if (i > 10 && i < 14) 
			cout << "WARNING: At,Ab,Atau are for SLHA1 only. "
			     << "Setting them to zero.\n"
			     << "Please use blocks TUIN, TDIN, TEIN for "
			     << "flavour violating SLHA input\n";
		      else {
			cout << "WARNING: did not understand parameter " 
			     << i << " in flavoured EXTPAR inputs\n";
		      }
		    }
		  }
		}
                else if (block == "QEXTPAR") {
                  int i; double d; kk >> i >> d;
                  if (susy_model == NMSSM) {
                     switch (i) {
                     case 61: // scale where to input lambda
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTlambda = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter lambda"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     case 62: // scale where to input kappa
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTkappa = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter kappa"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     case 65: // scale where to input <S>
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTsVev = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter <S>"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     case 66: // scale where to input xiF
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTxiF = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter xiF"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     case 68: // scale where to input mu'
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTmuPrime = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter mu'"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     default:
                        cout << "# WARNING: cannot use parameter " << i <<
                           " (set in QEXTPAR) as input at a different"
                           " scale (in the NMSSM) -- ignoring the scale choice\n";
                     }
                     continue;
                  }
                  cout << "# WARNING: cannot use parameter " << i <<
                     " (set in QEXTPAR) as input at a different"
                     " scale -- ignoring the scale choice\n";
                }
		else if (block == "VCKMIN") {
		  int i; double d; kk >> i >> d;
		  switch(i) {
		  case 1: lambdaW = d; break;
		  case 2: aCkm = d;   break;		  
		  case 3: rhobar = d; break;
		  case 4: etabar = d; break;
		  default:
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n";
		    break;
		  }
		}
		else if (block == "UMNSIN") {
		  int i; double d; kk >> i >> d;
		  switch(i) {
		  case 1: k.setThetaB12(asin(d)); break;
		  case 2: k.setThetaB23(asin(d)); break;		  
		  case 3: k.setThetaB13(asin(d)); break;
		  case 4: cout << "# Cannot yet do complex phases: ";
		    cout << "setting it to zero" << endl; break;
		  case 5: cout << "# Cannot yet do complex phases: ";
		    cout << "setting it to zero" << endl; break;
		  case 6: cout << "# Cannot yet do complex phases: ";
		    cout << "setting it to zero" << endl; break;
		  default:
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n";
		    break;
		  }
		}
		else if (block == "SMINPUTS") {
		  int i; double d; kk >> i >> d; 
		  switch (i) {
		  case 1: oneset.setAlpha(ALPHA, 1.0 / d); break;
		  case 2: GMU = d; break;
		  case 3: oneset.setAlpha(ALPHAS, d); break; 
		  case 4: oneset.setMu(d); m.setData(oneset); MZ = d; break;
		  case 5: oneset.setMass(mBottom, d); 
		    oneset.setMbMb(d); break;
		  case 6: oneset.setPoleMt(d); break;
		  case 7: oneset.setMass(mTau, d); 
		    oneset.setPoleMtau(d); break;
		  case 8: k.setMnuTau(d); break;
		  case 11: oneset.setMass(mElectron, d); k.setPoleMe(d); break;
		  case 12: k.setMnuMu(d); break;
		  case 13: oneset.setMass(mMuon, d); k.setPoleMmu(d); break;
		  case 14: k.setMnuTau(d); break;
		  case 21: oneset.setMass(mDown, d); k.setMd2GeV(d);
		    break;
		  case 22: oneset.setMass(mUp, d); k.setMu2GeV(d); break;
		  case 23: oneset.setMass(mStrange, d); k.setMs2GeV(d); break;
		  case 24: oneset.setMass(mCharm, d); k.setMcMc(d); break;
		  default: 
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		} 
		else if (block == "MSQ2IN") {
		  modelIdent = "nonUniversal";
		  int i, j; double d; kk >> i >> j >> d;
		  pars(positionOfSym(i, j) + 3) = d;
		}
		else if (block == "MSU2IN") {
		  modelIdent = "nonUniversal";
		  int i, j; double d; kk >> i >> j >> d;
		  pars(positionOfSym(i, j) + 9) = d;
		}
		else if (block == "MSD2IN") {
		  modelIdent = "nonUniversal";
		  int i, j; double d; kk >> i >> j >> d;
		  pars(positionOfSym(i, j) + 15) = d;
		}
		else if (block == "MSL2IN") {
		  modelIdent = "nonUniversal";
		  int i, j; double d; kk >> i >> j >> d;
		  pars(positionOfSym(i, j) + 21) = d;
	      }
		else if (block == "MSE2IN") {
		  modelIdent = "nonUniversal";
		  int i, j; double d; kk >> i >> j >> d;
		  pars(positionOfSym(i, j) + 27) = d;
		}
		else if (block == "TUIN") {
		  modelIdent = "nonUniversal";
		  int i, j; double d; kk >> i >> j >> d;
		  pars((i-1) * 3 + j + 33) = d;
		  slha2setTrilinear[(i-1) * 3 + j - 1] = true;
		}
		else if (block == "TDIN") {
		  modelIdent = "nonUniversal";
		  int i, j; double d; kk >> i >> j >> d;
		  pars((i-1) * 3 + j + 42) = d;
		  slha2setTrilinear[(i-1) * 3 + j + 8] = true;
		}
		else if (block == "TEIN") {
		  modelIdent = "nonUniversal";
		  int i, j; double d; kk >> i >> j >> d;
		  pars((i-1) * 3 + j + 51) = d;
		  slha2setTrilinear[(i-1) * 3 + j + 17] = true;
		}
		else if (block == "RVLAMLLEIN") {
		  int i,j,k; double d; kk >> i >> j >> k >> d;
		  if ((i > 0 && i <=  3) || (j > 0 && j <=  3) ||
		      (k > 0 && k <=  3)) {
		    kw.setLambda(LE, k, i, j, d);
		  }
		  else {
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << j << " " << k << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		}
		else if (block == "RVLAMLQDIN") {
		  int i,j,k; double d; kk >> i >> j >> k >> d;
		  if((i > 0 && i <=  3) || (j > 0 && j <=  3) ||
		     (k > 0 && k <=  3)) {
		    kw.setLambda(LD, k, i, j, d);
		  }
		  else {
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << j << " " << k << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		}
		else if (block == "RVLAMUDDIN") {
		  int i,j,k; double d; kk >> i >> j >> k >> d;
		  if((i > 0 && i <=  3) || (j > 0 && j <=  3) || 
		     (k > 0 && k <=  3)) {
		    kw.setLambda(LU, i, j, k, d);
		  }
		}
		else if (block == "RVTIN") {
		  int i,j,k; double d; kk >> i >> j >> k >> d;
		  if((i > 0 && i <=  3) || (j > 0 && j <=  3) || 
		     (k > 0 && k <=  3)) {
		    kw.setHr(LE, k, i, j, d);
		  }
		  else {
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << j << " " << k << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		}
		else if (block == "RVTPIN") {
		  int i,j,k; double d; kk >> i >> j >> k >> d;
		  if((i > 0 && i <=  3) || (j > 0 && j <=  3) || 
		     (k > 0 && k <=  3)) {
		    kw.setHr(LD, k, i, j, d);
		  }
		  else {
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << j << " " << k << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		}
		else if (block == "RVTPPIN") {
		  int i,j,k; double d; kk >> i >> j >> k >> d;
		  if((i > 0 && i <=  3) || (j > 0 && j <=  3) ||
		     (k > 0 && k <=  3)) {
		    kw.setHr(LU, i, j, k, d);
		  }
		  else {
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << j << " " << k << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		}
		else if (block == "RVKAPPAIN") {
		  int i; double d; kk >> i >> d;
		  if (i > 0 && i <=  3) {
		    kw.setKappa(i, d);
		}
		  else {
		    cout << "# WARNING: Don't understand data input " << i << d 
		       << " in block " << block << ": ignoring it\n"; break;
		  }
		}
		else if (block == "RVDIN") {
		  int i; double d; kk >> i >> d;
		  if (i > 0 && i <=  3) {
		    kw.setD(i, d);
		  }
		  else {
		    cout << "# WARNING: Don't understand data input " 
			 << i << d << " in block " << block 
			 << ": ignoring it\n"; break;
		  }
		}
		// input of sneutrino VEVs not supported yet
		else if (block == "RVSNUVEVIN") {
		  int i; double d; kk >> i >> d;
		  cout << "# WARNING: in block " << block 
		       << ": SOFTSUSY does not support setting of sneutrino VEVs"
		       << " yet : ignoring them\n"; break;
		}
		else if (block == "RVMLH1SQIN") {
		  int i; double d; kk >> i >> d;
		  if ((i > 0 && i <=  3)) {
		    kw.setMh1lSquared(i, d);
		  }
		  else {
		    cout << "# WARNING: Don't understand data input " 
			 << i << d << " in block " << block 
			 << ": ignoring it\n"; break;
		  }
		}
		else if (block == "SOFTSUSY") {
		  int i; double d; kk >> i >> d;
		  switch(i) {
		  case 1: TOLERANCE = d; break;
		  case 2: 
		    MIXING = int(d); 
		    //if (MIXING > 0) flavourViolation = true;
		    break;
		  case 3: PRINTOUT = int(d); break;
		  case 4: QEWSB = d; break;
		  case 5: INCLUDE_2_LOOP_SCALAR_CORRECTIONS = 
		      bool(int(d+EPSTOL)); break;
		  case 6: outputCharacteristics(int(d+EPSTOL)-1); break;  
		  case 7: {
		    int num = int(d+EPSTOL);
		    if (num != 1 && num!= 2) {
		      cout << "# WARNING: Can only set number of loops for"
			   << " higgs masses and REWSB to be 1 or 2 in "
			   << " BLOCK SOFTSUSY parameter 7, not " << num 
			   << ". Ignoring.\n";
		  } else {
		      numHiggsMassLoops = num;
		      numRewsbLoops = num;
		    }
		  }
		    break;
		  case 8: {
		  int num = int(d + EPSTOL);
		  if (num != 0 && num != 1) {
		    cout << "# Incorrect flag in BLOCK SOFTSUSY 8: " 
			 << num << ". Ignoring it.\n";
		  } else {
		    if (num == 1) susyRpvBCatMSUSY = true;
		    else susyRpvBCatMSUSY = false;
		  }
		  }
		    break;
		  case 9: {
		    int num = int(d + EPSTOL);		  
		    if (num != 0 && num != 1) {
		      cout << "# Incorrect flag in BLOCK SOFTSUSY 9: " 
			   << num << ". Ignoring it.\n";
		    } else {
		      if (num == 1) kw.setInvertedOutput();
		      else kw.setNormalOutput();		      
		  }
		  }
		    break;
		  case 10: {
		    int num = int(d + EPSTOL);
		    if (num != 0 && num != 1) {
		      cout << "# Incorrect flag in BLOCK SOFTSUSY 10: " 
			   << num << ". Ignoring it.\n";
		    } else if (num == 1) forceSlha1 = true;
		    else forceSlha1 = false;
		  }
		    break;
		  case 11: {
		    /// Set gravitino mass
		    r->setM32(d); 
		  }
		    break;
		  case 12: {
		    int num = int(d + EPSTOL);
		    if (num == 1) printRuledOutSpectra = true;
		    else if (num == 0) printRuledOutSpectra = false;
		    else cout << "# WARNING: Don't understand Block SOFTSUSY "
			      << "parameter 12 " << d << ". Ignoring it." 
			      << endl;
		  }
		    break;
		  case 13: {
		    int num = int(d + EPSTOL);
		    if (num == 1) mAFlag = true;		  
		  }
                    break;
                  case 15: {
                    int num = int(d + EPSTOL);
                    softsusy::NMSSMTools = num;
                  }
		    break;
                  case 16: {
                     int num = int(d + EPSTOL);
                     softsusy::MICROMEGAS = num;
                  }
                    break;
                  case 17: {
                     int num = int(d + EPSTOL);
                     softsusy::NMSDECAY = num;
                  }
		    break;
                  case 18: {
                    int num = int(d + EPSTOL);
                    if(num == 1) softsusy::SoftHiggsOut = true;
                  }
                     break;
#ifdef COMPILE_THREE_LOOP_RGE
		  case 19: {
                    int num = int(d + EPSTOL);
		    if (num == 1) USE_THREE_LOOP_RGE = true;
		    else if (num == 0) USE_THREE_LOOP_RGE = false;
		    else cout << "WARNING: incorrect setting for SOFTSUSY Block 19 (should be 0 or 1)\n";
		    break;			     
		  }
#endif
#ifdef COMPILE_TWO_LOOP_GAUGE_YUKAWA
		  case 20: {
                    int num = int(d + EPSTOL);
		    // AVB: can be set to just 1 Turn on all thresholds
		    //      can be set to 1 + 2 * ( flags for included thresholds)
		    //      to have a finer control over included thresholds
		    if (num % 2 == 1) {
		      USE_TWO_LOOP_GAUGE_YUKAWA = true;
		      if (num == 1) 
			r->included_thresholds = (ENABLE_TWO_LOOP_AS_AS_YUK |
						  ENABLE_TWO_LOOP_MT_AS | 
						  ENABLE_TWO_LOOP_MB_AS | 
						  ENABLE_TWO_LOOP_MB_YUK | 
						  ENABLE_TWO_LOOP_MTAU_YUK);
		      else 
			r->included_thresholds = ((num >> 1) & 
						  (ENABLE_TWO_LOOP_AS_AS_YUK | 
						   ENABLE_TWO_LOOP_MT_AS | 
						   ENABLE_TWO_LOOP_MB_AS | 
						   ENABLE_TWO_LOOP_MB_YUK | 
						   ENABLE_TWO_LOOP_MTAU_YUK));
				   
		    } else if (num == 0) USE_TWO_LOOP_GAUGE_YUKAWA = false;
		    else cout << "WARNING: incorrect setting for SOFTSUSY Block 20 (should be 0 or 1)\n";
		    break;
		  }
#endif
		  default:
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		}
		else {
		  cout << "# WARNING: don't recognise block " << block 
		       << ": ignoring all data in it" << endl;
		}
		// end if blocks
	      
	      } // end of data
	    } // end of no-comment
            
	} // end of file
          
        }
	else errorCall();
      }

      /// prepare CKM angles
      if (flavourViolation || RPVflag) k.setAngles(lambdaW, aCkm, rhobar, 
						   etabar);
      
      if (r->displayAltEwsb()) {
	if (strcmp(modelIdent, "splitgmsb")) {
	  //	boundaryCondition = &extendedSugraBcs2;
	  r->setSusyMu(pars(23)); 
	} else {
	  ostringstream ii;
	  ii << "Split GMSB BCs should not supported with alternative EWSB\n";
	  throw ii.str();
	  /// Split GMSB BCs: different
	  /*	r->setSusyMu(400.);
		r->setMuCond(400.);
		r->setMaCond(400.);*/
	}
	sgnMu = 0; // Flags different BCs
      }

      // set NMSSM boundary conditions
      if (susy_model == NMSSM) {
         softsusy::Z3 = nmssm_input.is_Z3_symmetric();

         if (flavourViolation) {
            string msg("# Error: flavour violation in the NMSSM is currenty"
                       " not supported\n");
            throw msg;
         }
         if (strcmp(modelIdent, "sugra") == 0) {
           // if we chose the soft Higgs masses as EWSB output, we
           // must not use msugra, because it would overwrite mHu2,
           // mHd2, mS2 at the GUT scale
           if (softsusy::SoftHiggsOut) {
             if (pars.size() != 6)
               pars.setEnd(6);
             pars(4) = nmssm_input.get(NMSSM_input::mu);
             pars(5) = nmssm_input.get(NMSSM_input::BmuOverCosBetaSinBeta);
             pars(6) = nmssm_input.get(NMSSM_input::xiS);
             nmssmBoundaryCondition = &NmssmSugraNoSoftHiggsMassBcs;
           } else {
             if (softsusy::Z3) {
               // Here we must use SemiMsugraBcs to avoid setting mS2
               // at the GUT scale
               if (pars.size() != 5)
                 pars.setEnd(5);
               pars(4) = pars(3); // sets Al to A0
               pars(5) = pars(3); // sets Ak to A0
               nmssmBoundaryCondition = &SemiMsugraBcs;
             } else {
               nmssmBoundaryCondition = &NmssmMsugraBcs;
             }
           }
         } else if (strcmp(modelIdent, "nonUniversal") == 0) {
           nmssmBoundaryCondition = &extendedNMSugraBcs;
           if (pars.size() != 56) {
              string msg("# Error: NMSSM non-minmal sugra boundary condition"
                         " chosen, but pars does not have 56 entries\n");
              throw msg;
           }
           if (nmssm_input.is_set(NMSSM_input::mu))
             pars(54) = nmssm_input.get(NMSSM_input::mu);
           if (nmssm_input.is_set(NMSSM_input::BmuOverCosBetaSinBeta))
             pars(55) = nmssm_input.get(NMSSM_input::BmuOverCosBetaSinBeta);
           if (nmssm_input.is_set(NMSSM_input::xiS))
             pars(56) = nmssm_input.get(NMSSM_input::xiS);
         } else {
            string msg("# Error: non-sugra boundary conditions for the NMSSM"
                       " are currently not supported\n");
            throw msg;
         }
      }

      if (RPVflag) {
	kw.rpvDisplay(pars);
	kw.setFlavourSoftsusy(k);
	r = &kw;
	
	if (boundaryCondition == &sugraBcs) 
	  boundaryCondition = &rpvSugraBcs;
	else if (boundaryCondition == &amsbBcs) 
	  boundaryCondition = &rpvAmsbBcs;
	else if (boundaryCondition == &gmsbBcs) 
	  boundaryCondition = &rpvGmsbBcs;
	else if (boundaryCondition == &extendedSugraBcs) {
	  boundaryCondition = &rpvExtendedSugraBcs;
	}
	else if (boundaryCondition == &extendedSugraBcs) {
	  boundaryCondition = &rpvExtendedSugraBcs;
	  kw.useAlternativeEwsb();
	}
	else {
	  ostringstream ii;
	  ii << "# ERROR: there is no RPV version for selected "
	     << " boundary condition. If you have" << endl
	     << "# *BOTH* Flavour violating *AND* RPV flags switched on, "
	     << " remove the flavour" << endl
	     << " violating one." << endl;
	  throw ii.str();
	}
      }
      // intput error checking  
      if (sgnMu != 1 && sgnMu != -1 && sgnMu != 0) {
	ostringstream ii;
	ii << "Incorrect input for sign(mu)=" << sgnMu <<endl;
	throw ii.str();
      }
      if (tanb < 1.0 || tanb > 5.0e2)  {
	ostringstream ii;
	ii << "Incorrect input for tan beta=" << tanb <<endl;
	throw ii.str();
      }
      
      oneset.toMz();

      cout << "susy_model=" << susy_model << endl;
      
    switch (susy_model) {
    case MSSM: {
      /// Switch off 3-loop RGEs etc
      double cs = 0.;
      double omega2=asin(2.), msqAv2 = 0.;  bool uni = gaugeUnification;
      USE_THREE_LOOP_RGE = false;   USE_TWO_LOOP_GAUGE_YUKAWA = false;
      MssmSoftsusy twoLoop;
      twoLoop.useAlternativeEwsb(); twoLoop.setMuCond(2500.); 
      twoLoop.setMaCond(1580.); 
      getMssmAndOmega(twoLoop, pars, tanb, sgnMu, oneset, mgutGuess, 
		      uni, omega2, msqAv2, boundaryCondition, ewsbBCscale, cs);
      //      cout << "twoLoop=" << twoLoop << " omega2=" << omega2 
      //   << " cs=" << cs << endl;

      double csOdd = 0., omegaOdd = 0., msqAvOdd = 0.;
      USE_THREE_LOOP_RGE = true;   USE_TWO_LOOP_GAUGE_YUKAWA = false;
      MssmSoftsusy oddLoop;
      oddLoop.useAlternativeEwsb(); oddLoop.setMuCond(2500.); 
      oddLoop.setMaCond(1580.); 
      getMssmAndOmega(oddLoop, pars, tanb, sgnMu, oneset, mgutGuess, 

		      uni, omegaOdd, msqAvOdd, boundaryCondition, ewsbBCscale, 
		      csOdd);
      //      cout << "oddLoop=" << oddLoop;

      /// Just 2-loop thresholds for strong coupling constant
      double omegaAs = asin(2.), msqAvAs = 0., csAs = 0.; mgutGuess = 2.5e3;
      USE_THREE_LOOP_RGE = false;
      MssmSoftsusy twoLoopAs; 
      twoLoopAs.useAlternativeEwsb(); twoLoopAs.setMuCond(2500.); 
      twoLoopAs.setMaCond(1580.); 
      twoLoopAs.included_thresholds |= ENABLE_TWO_LOOP_AS_AS_YUK;
      USE_TWO_LOOP_GAUGE_YUKAWA = true;
      getMssmAndOmega(twoLoopAs, pars, tanb, sgnMu, oneset, mgutGuess, 
		      uni, omegaAs, msqAvAs, boundaryCondition, ewsbBCscale, 
		      csAs); 
      //      cout << "twoLoopAs=" << twoLoopAs;

      /// Just 2-loop strong thresholds for mt
      USE_TWO_LOOP_GAUGE_YUKAWA = false;
      double omegaMt = asin(2.), msqAvMt = 0., csMt = 0.; mgutGuess = 2.5e3;
      MssmSoftsusy twoLoopMt; 
      twoLoopMt.useAlternativeEwsb(); twoLoopMt.setMuCond(2500.); 
      twoLoopMt.setMaCond(1580.); 
      twoLoopMt.included_thresholds |= ENABLE_TWO_LOOP_MT_AS;
      USE_TWO_LOOP_GAUGE_YUKAWA = true;
      getMssmAndOmega(twoLoopMt, pars, tanb, sgnMu, oneset, mgutGuess, 
		      uni, omegaMt, msqAvMt, boundaryCondition, ewsbBCscale,
		      csMt); 
      //    cout << twoLoopMt;    
    
      /// Just 2-loop for mb,mtau
      USE_TWO_LOOP_GAUGE_YUKAWA = false;
      double omegaMb = asin(2.), msqAvMb = 0., csMb = 0.; mgutGuess = 2.5e3;
      MssmSoftsusy twoLoopMb; 
      twoLoopMb.useAlternativeEwsb(); twoLoopMb.setMuCond(2500.); 
      twoLoopMb.setMaCond(1580.); 
      twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MB_AS;
      twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MB_YUK;
      twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MTAU_YUK;
      USE_TWO_LOOP_GAUGE_YUKAWA = true;
      getMssmAndOmega(twoLoopMb, pars, tanb, sgnMu, oneset, mgutGuess, 
		      uni, omegaMb, msqAvMb, boundaryCondition, ewsbBCscale, 
		      csMb); 
      //    cout << "twoLoopMb=" << twoLoopMb;    

      /// 3-loop etc ON
      double omega3 = asin(2.), msqAv3 = 0., cs3 = 0.; mgutGuess = 2.5e3;
      USE_THREE_LOOP_RGE = true; 
      USE_TWO_LOOP_GAUGE_YUKAWA = true; 
      MssmSoftsusy threeLoop;
      threeLoop.useAlternativeEwsb(); threeLoop.setMuCond(2500.); 
      threeLoop.setMaCond(1580.); 
      getMssmAndOmega(threeLoop, pars, tanb, sgnMu, oneset, mgutGuess, 
		      uni, omega3, msqAv3, boundaryCondition, ewsbBCscale, 
		      cs3); 
      
      //    cout << "threeLoop=" << threeLoop;
      
      writeTable(twoLoop, oddLoop, twoLoopAs, twoLoopMt, twoLoopMb, threeLoop, 
		 omega2, omegaOdd, omegaAs, omegaMt, omegaMb, omega3,
		 msqAv2, msqAvOdd, msqAvAs, msqAvMt, msqAvMb, msqAv3, 
		 cs, csOdd, csAs, csMt, csMb, cs3);
    }
      break;
    case NMSSM: {
      nmssm_input.check_setup();

      DoubleVector nmpars(nmssm_input.get_nmpars());
      nmssm.lowOrg(nmssmBoundaryCondition, mgutGuess, pars, nmpars, sgnMu,
                   tanb, oneset, gaugeUnification, ewsbBCscale);
      nmssm.lesHouchesAccordOutput(cout, modelIdent, pars, sgnMu, tanb, qMax,
                                   numPoints, ewsbBCscale);
      if (nmssm.displayProblem().test()) {
         cout << "# SOFTSUSY problem with NMSSM point: "
              << nmssm.displayProblem() << endl;
         return -1;
      }
      }
      break;
    default:
      cout << "# Error: unknown susy model " << susy_model
           << ", please check your MODSEL (entry 3) settings" << endl;
      break;
    }
  }
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }
  
  return 0;
}

