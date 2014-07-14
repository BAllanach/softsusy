
/** 
   Project:     SOFTSUSY 
   File:        scan.cpp
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

void printLineOut(double m0, double m12, double a0, double tanb, 
		  MssmSoftsusy twoLoop, MssmSoftsusy twoLoopAs, 
		  MssmSoftsusy & twoLoopMt, MssmSoftsusy & twoLoopMb, 
		  MssmSoftsusy & threeLoop, double omega2, double omegaAs, 
		  double omegaMt, double omegaMb, double omega3, 
		  double msqAv2, double msqAvAs, double msqAvMt, 
		  double msqAvMb, double msqAv3, double cs, double csAs, 
		  double csMt, double csMb, double cs3, double dA, 
		  double dAAs, double dAMt, double dAMb, double dA3, 
		  double dY, double dYAs, double dYMt, 
		  double dYMb, double dY3, double dYp, double dYpAs, 
		  double dYpMt, double dYpMb, double dYp3) {
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
	 << " " << threeLoop.displayGaugeCoupling(3)             ///< 134
	 << " " << cs                                             ///< 135
	 << " " << csAs                                           ///< 136
	 << " " << csMt                                           ///< 137
	 << " " << csMb                                           ///< 138
	 << " " << cs3                                            ///< 139
	 << " " << QEWSB                                          ///< 140
	 << " " << dA                                             ///< 141
	 << " " << dAAs                                           ///< 142
	 << " " << dAMt                                           ///< 143
	 << " " << dAMb                                           ///< 144
	 << " " << dA3                                            ///< 145
	 << " " << dY                                             ///< 146
	 << " " << dYAs                                           ///< 147
	 << " " << dYMt                                           ///< 148
	 << " " << dYMb                                           ///< 149
	 << " " << dY3                                            ///< 150
	 << " " << dYp                                            ///< 151
	 << " " << dYpAs                                          ///< 152
	 << " " << dYpMt                                          ///< 153
	 << " " << dYpMb                                          ///< 154
	 << " " << dYp3;                                          ///< 155
    if (twoLoop.displayProblem().test()) cout << "# 2-loop problem: " 
					      << twoLoop.displayProblem();
    if (threeLoop.displayProblem().test()) cout << "# 3-loop problem " 
						<< threeLoop.displayProblem(); 
    cout << endl;
}

double getCrossSection(MssmSoftsusy & r, char * fileName) {
  double msq = (r.displayPhys().mu(1, 1) + r.displayPhys().mu(2, 1) + 
	 r.displayPhys().md(1, 1) + r.displayPhys().mu(2, 1)) / 4;
  double mg = r.displayPhys().mGluino;
  double mt = r.displayDataSet().displayPoleMt();

  fstream fout("/home/bca20/code/prospino/input", ios::out); 
  fout << msq << " " << mg << " " << mt << " " << 1. << endl;
  fout.close();
  char buff[500];
  sprintf(buff, "cd /home/bca20/code/prospino; cat input | ./crosssec | grep -v NAN > output");
  int err = system(buff);
  double xs = 0.;
  if (!err) {
  fstream fin("/home/bca20/code/prospino/output", ios::in); 
  fin >> xs;
  } else  cout << "CROSS SECTION ERROR\n";
  
  cout << xs * 1.0e3;
  return xs * 1.0e3;
  /*  char buff[500];
  sprintf(buff, "cp lesHout ../../code/prospino2.1/prospino.in.les_houches; cd ../../code/prospino2.1/; ./prospino_2.run output > err");
  int err = system(buff);
  double gg = 0., ss = 0., sb=0., sg = 0.;
  if (!err) {
  fstream fin("/home/bca20/code/prospino2.1/prospino.dat", ios::in); 
  string o;
  fin >> o >> o >> o >> o >> o >> o >> o 
	>> o >> o >> o >> o >> sg >> o 
	>> o >> o >> o;
  fin >> o >> o >> o >> o >> o >> o >> o 
	>> o >> o >> o >> o >> ss >> o 
	>> o >> o >> o;
  fin >> o >> o >> o >> o >> o >> o >> o 
	>> o >> o >> o >> o >> sb >> o 
	>> o >> o >> o;
  fin >> o >> o >> o >> o >> o >> o >> o 
	>> o >> o >> o >> o >> gg >> o 
	>> o >> o >> o;
  }
  else cout << "CROSS SECTION ERROR\n";

  sprintf(buff, "cd ../../code/softsusy"); 
  err = system(buff);
  if (err) cout << "SHIT\n";

  return (gg + ss + sb + sg) * 1.0e3; ///< total in fb*/
}

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

void writeTable(MssmSoftsusy & twoLoop, MssmSoftsusy & twoLoopAs, 
		MssmSoftsusy & twoLoopMt, MssmSoftsusy & twoLoopMb, 
		MssmSoftsusy & threeLoop, double omega2, double omega2As, 
		double omega2Mt, double omega2Mb, double omega3, 
		double msqAv2, double msqAv2As, double msqAv2Mt, 
		double msqAv2Mb, double msqAv3, double cs, double csAs, 
		double csMt, double csMb, double cs3, double dAs, 
		double dAsAs, double dAsMt, double dAsMb, double dAs3, 
		double dY, double dYAs, double dYMt, 
		double dYMb, double dY3, double dYp, double dYpAs, 
		double dYpMt, double dYpMb, double dYp3) {
  cout << "\\begin{table}\n\\begin{center}\n\\begin{tabular}{|c|c|ccccccc|}"
       << "\\hline\nThreshold & RGEs  & $m_h$  & $m_{\\tilde g}$ & "
       << "$m_{{\\tilde q}}$ & $m_{\\chi_1^0}$  & $m_{\\chi_2^0}$ & "
       << "$m_{\\chi_3^0}$ & $m_{\\chi_4}^0$ \\\\ \\hline\nNone"
       << "               & 2 & ";
  if (omega2 != omega2) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A\\\\\n"); else   printf("%5.1f & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f &%4.0f\\\\\n",
	 twoLoop.displayPhys().mh0(1), 
	 twoLoop.displayPhys().mGluino, 
	 msqAv2, 
	 fabs(twoLoop.displayPhys().mneut(1)), 
	 fabs(twoLoop.displayPhys().mneut(2)), 
	 fabs(twoLoop.displayPhys().mneut(3)), 
	 fabs(twoLoop.displayPhys().mneut(4))
	 );
  cout << "$\\Delta \\alpha_s$  & 2 & ";
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
  cout << "$\\Delta m_t$      & 2 & ";
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
  cout << "$\\Delta m_b, m_\\tau$& 2 & ";
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
  cout << "$\\Delta$ All      & 3 & ";
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
       << "&& $m_{{\\tilde t}_L}$  & $m_{{\\tilde t}_R}$ &$m_{{\\tilde b}_L}$&"
       << "$m_{{\\tilde b}_R}$&$m_{{\\tilde \\tau}_L}$&$m_{{\\tilde \\tau}_R}$&"
       << "$m_{\\chi_1}^\\pm$ \\\\ \\hline\n"
       << "None             & 2 &";
  if (omega2 != omega2) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n"); else   printf("%4.0f & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f &%4.0f\\\\\n",
	 twoLoop.displayPhys().mu(1, 3), 
	 twoLoop.displayPhys().mu(2, 3), 
	 twoLoop.displayPhys().md(1, 3), 
	 twoLoop.displayPhys().md(2, 3), 
	 twoLoop.displayPhys().me(1, 3), 
	 twoLoop.displayPhys().me(2, 3), 
	 twoLoop.displayPhys().mch(1)
	 );
  cout << "$\\Delta \\alpha_s$  & 2 & ";
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
  cout << "$\\Delta m_t$      & 2 & ";
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
  cout << "$\\Delta m_b, m_\\tau$& 2 & ";
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
  cout << "$\\Delta$ All      & 3 & ";
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
  cout << "      &  & $g_3(M_{SUSY})$ & $Y_t(M_{SUSY})$ & "
       << " $Y_b(M_{SUSY})$ & $Y_\\tau(M_{SUSY})$  & $\\mu(M_{SUSY})$"
       << "    & $\\Omega_{CDM} h^2$ & $\\sigma_{SUSY}^{TOT}$\\\\ \\hline\n"
       << " None                   & 2 & ";

  if (omega2 != omega2) printf("N/A & N/A & N/A & N/A & N/A & N/A &\\\\\n"); else printf("%5.3f & %5.3f & %5.3f & %5.3f & %4.0f & %5.1f & %5.1f\\\\\n",
	 twoLoop.displayGaugeCoupling(3), 
	 twoLoop.displayYukawaElement(YU, 3, 3), 
	 twoLoop.displayYukawaElement(YD, 3, 3), 
	 twoLoop.displayYukawaElement(YE, 3, 3), 
	 twoLoop.displaySusyMu(),
											 omega2, cs
	 );

  cout << "$\\Delta \\alpha_s$  & 2 & ";
  if (omega2As != omega2As) printf("N/A & N/A & N/A & N/A & N/A & N/A &\\\\\n"); else printf("%+5.3f & %+5.3f & %+5.3f & %+5.3f & %+4.0f & %+5.1f & %+5.1f\\\\\n",
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
  cout << "$\\Delta m_t$      & 2 & ";
    if (omega2Mt != omega2Mt) printf("N/A & N/A & N/A & N/A & N/A & N/A& N/A\\\\\n"); else printf("%+5.3f & %+5.3f & %+5.3f & %+5.3f & %+4.0f & %+5.1f & %+5.1f\\\\\n",
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
    cout << "$\\Delta m_b, m_\\tau$& 2 & ";
    if (omega2Mb != omega2Mb) printf("N/A & N/A & N/A & N/A & N/A & N/A & N/A\\\\\n"); else printf("%+5.3f & %+5.3f & %+5.3f & %+5.3f & %+4.0f & %+5.3f & %+5.1f\\\\\n",
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
    cout << "$\\Delta$ All      & 3 & ";
    if (omega3 != omega3) printf("N/A & N/A & N/A & N/A & N/A & N/A & \\\\\n"); 
    else printf("%+5.3f & %+5.3f & %+5.3f & %+5.3f & %+4.0f & %+5.1f &%+5.1f\\\\\n",
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

    twoLoop.runto(twoLoop.displayMxBC());
    twoLoopAs.runto(twoLoop.displayMxBC());
    twoLoopMt.runto(twoLoop.displayMxBC());
    twoLoopMb.runto(twoLoop.displayMxBC());
    threeLoop.runto(twoLoop.displayMxBC());

  cout << "\n%\n\\hline"
       << "&& $M_{GUT}/10^{16}$& $1/\\alpha_{GUT}$&$\\Delta (\\alpha)$ & $\\Delta Y_{b\\tau}$ & $\\Delta Y_{tb}$&"
       << "&"
       << " \\\\ \\hline"
       << " None                   & 2 & ";
  if (omega2 != omega2) printf("N/A & N/A & N/A & N/A & N/A & N/A &\\\\\n"); 
  else printf("%5.3f & %5.3f & %5.3f & %5.3f & %5.3f& & \\\\\n", 
	      twoLoop.displayMxBC() * 1.0e-16,
	      (4 * PI) / sqr(twoLoop.displayGaugeCoupling(1)),
	      (sqr(twoLoop.displayGaugeCoupling(3)) - sqr(twoLoop.displayGaugeCoupling(1))) / (4 * PI), 
	      twoLoop.displayYukawaElement(YD, 3, 3) - twoLoop.displayYukawaElement(YE, 3, 3), 
	      twoLoop.displayYukawaElement(YU, 3, 3) - 	 twoLoop.displayYukawaElement(YD, 3, 3));

  cout << "$\\Delta \\alpha_s$  & 2 & ";
  if (omega2As != omega2As) printf("N/A & N/A & N/A & N/A & & &\\\\\n"); 
  else printf(" %5.3f &  %5.3f &  %5.3f & %5.3f &%5.3f & & \\\\\n", 
	      (twoLoopAs.displayMxBC() - twoLoop.displayMxBC()) * 1.0e-16,
	      (4 * PI) / sqr(twoLoopAs.displayGaugeCoupling(1)) - 
	      (4 * PI) / sqr(twoLoop.displayGaugeCoupling(1)),
	      (sqr(twoLoopAs.displayGaugeCoupling(3)) - 
	       sqr(twoLoopAs.displayGaugeCoupling(1))) / (4 * PI),
	      twoLoopAs.displayYukawaElement(YD, 3, 3) - 
	      twoLoopAs.displayYukawaElement(YE, 3, 3) , 
	      (twoLoopAs.displayYukawaElement(YU, 3, 3) - 	 
	       twoLoopAs.displayYukawaElement(YD, 3, 3)));

  cout << "$\\Delta m_t$  & 2 & ";
  if (omega2Mt != omega2Mt) printf("N/A & N/A & N/A & N/A &N/A  & &\\\\\n"); 
  else printf(" %5.3f &  %5.3f &  %5.3f & %5.3f &%5.3f & & \\\\\n" , 
	      (twoLoopMt.displayMxBC() - twoLoop.displayMxBC()) * 1.0e-16,
	      4 * PI / sqr(twoLoopMt.displayGaugeCoupling(1)) / - 
	      4 * PI / sqr(twoLoop.displayGaugeCoupling(1)),
	      (sqr(twoLoopMt.displayGaugeCoupling(3)) - 
	       sqr(twoLoopMt.displayGaugeCoupling(1))) / (4 * PI),
	      twoLoopMt.displayYukawaElement(YD, 3, 3) - 
	      twoLoopMt.displayYukawaElement(YE, 3, 3) , 
	      (twoLoopMt.displayYukawaElement(YU, 3, 3) - 	 
	       twoLoopMt.displayYukawaElement(YD, 3, 3)));

  cout << "$\\Delta m_b,m_\\tau$  & 2 & ";
  if (omega2Mb != omega2Mb) printf("N/A & N/A & N/A & N/A & & &\\\\\n"); 
  else printf(" %5.3f &  %5.3f &  %5.3f & %5.3f &%5.3f & & \\\\\n",
	      (twoLoopMb.displayMxBC() - twoLoop.displayMxBC()) * 1.0e-16,
	      4 * PI / sqr(twoLoopMb.displayGaugeCoupling(1)) - 
	      4 * PI / sqr(twoLoop.displayGaugeCoupling(1)),
	      (sqr(twoLoopMb.displayGaugeCoupling(3)) - 
	       sqr(twoLoopMb.displayGaugeCoupling(1))) / (4 * PI),
	      twoLoopMb.displayYukawaElement(YD, 3, 3) - 
	      twoLoopMb.displayYukawaElement(YE, 3, 3) , 
	      (twoLoopMb.displayYukawaElement(YU, 3, 3) - 	 
	       twoLoopMb.displayYukawaElement(YD, 3, 3)));
  
  cout << "$\\Delta$ All  & 3 & ";
  if (omega3 != omega3) printf("N/A & N/A & N/A & N/A & & &\\\\\n"); 
  else printf(" %5.3f &  %5.3f &  %5.3f & %5.3f &%5.3f & & \\\\\n", 
	      (threeLoop.displayMxBC() - twoLoop.displayMxBC()) * 1.0e-16,
	      4 * PI / sqr(threeLoop.displayGaugeCoupling(1)) - 
	      4 * PI / sqr(twoLoop.displayGaugeCoupling(1)),
	      (sqr(threeLoop.displayGaugeCoupling(3)) - 
	       sqr(threeLoop.displayGaugeCoupling(1))) / (4 * PI),
	      threeLoop.displayYukawaElement(YD, 3, 3) - 
	      threeLoop.displayYukawaElement(YE, 3, 3), 
	      (threeLoop.displayYukawaElement(YU, 3, 3) - 	 
	       threeLoop.displayYukawaElement(YD, 3, 3)));
  

  cout << "\\hline\n\\end{tabular}\n\\end{center}\n";
}

/// Returns the object along with omega. Oneset should already be fixed at MZ
void getCmssmAndOmega(MssmSoftsusy & r, DoubleVector & pars, const double tanb, 
			      const int sgnMu, const QedQcd & oneset, 
			      double mGutGuess, bool uni, double & omega, 
		      double & msqAv, double & cs, double & deltaAs, 
		      double & deltaY, double & deltaYp) {
  r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
  r.setData(oneset);

  /// Produces SLHA output file
  char fileName[500]; 
  sprintf(fileName,"lesHout");
  produceSLHAfile(r, fileName, sgnMu, tanb, pars);
  cs = getCrossSection(r, fileName);
  if (!r.displayProblem().test()) 
    omega = doDarkMatter(pars, tanb, sgnMu, fileName);

  msqAv = (r.displayPhys().mu(2, 1) + 
	   r.displayPhys().mu(1, 1) +
	   r.displayPhys().md(2, 1) + 
	   r.displayPhys().md(1, 1)) * 
    0.25;
  
  MssmSoftsusy c(r);
  c.runto(c.displayMxBC());
  deltaAs = sqr(c.displayGaugeCoupling(3)) / (4.0 * PI) - 
    sqr(c.displayGaugeCoupling(1)) / (4.0 * PI);
  deltaY  = c.displayYukawaElement(YD, 3, 3) - 
    c.displayYukawaElement(YE, 3, 3);
  deltaYp = c.displayYukawaElement(YU, 3, 3) -
    c.displayYukawaElement(YD, 3, 3);

  remove(fileName); 
  return;
}

int main(int argc, char *argv[]) {
  TOLERANCE = 1.0e-5;
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

  if (argc == 6) {
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
    double omega2=asin(2.), msqAv2 = 0., cs = 0., dAs = 0., dY = 0., dYp = 0.;  
    USE_THREE_LOOP_RGE = false;   USE_TWO_LOOP_THRESHOLD = false;
    MssmSoftsusy twoLoop;
    getCmssmAndOmega(twoLoop, pars, tanb, sgnMu, oneset, mGutGuess, 
		     uni, omega2, msqAv2, cs, dAs, dY, dYp);

    /// Just 2-loop thresholds for strong coupling constant
    double omegaAs = asin(2.), msqAvAs = 0., csAs = 0., dAsAs = 0., 
      dYAs = 0., dYpAs = 0.; mGutGuess = 2.e16;
    MssmSoftsusy twoLoopAs; 
    twoLoopAs.included_thresholds |= ENABLE_TWO_LOOP_AS_AS_YUK;
    USE_TWO_LOOP_THRESHOLD = true;
    getCmssmAndOmega(twoLoopAs, pars, tanb, sgnMu, oneset, mGutGuess, 
    		     uni, omegaAs, msqAvAs, csAs, dAsAs, dYAs, dYpAs); 

    /// Just 2-loop strong thresholds for mt
    USE_TWO_LOOP_THRESHOLD = false;
    double omegaMt = asin(2.), msqAvMt = 0., csMt = 0., dAsMt = 0., 
      dYMt = 0., dYpMt = 0.; mGutGuess = 2.e16;
    MssmSoftsusy twoLoopMt; 
    twoLoopMt.included_thresholds |= ENABLE_TWO_LOOP_MT_AS;
    USE_TWO_LOOP_THRESHOLD = true;
    getCmssmAndOmega(twoLoopMt, pars, tanb, sgnMu, oneset, mGutGuess, 
    		     uni, omegaMt, msqAvMt, csMt, dAsMt, dYMt, dYpMt); 

    /// Just 2-loop for mb,mtau
    USE_TWO_LOOP_THRESHOLD = false;
    double omegaMb = asin(2.), msqAvMb = 0., csMb = 0., dAsMb = 0., 
      dYMb = 0., dYpMb = 0.; mGutGuess = 2.e16;
    MssmSoftsusy twoLoopMb; 
    twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MB_AS;
    twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MB_YUK;
    twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MTAU_YUK;
    USE_TWO_LOOP_THRESHOLD = true;
    getCmssmAndOmega(twoLoopMb, pars, tanb, sgnMu, oneset, mGutGuess, 
    		     uni, omegaMb, msqAvMb, csMb, dAsMb,
		     dYMb, dYpMb); 

    /// 3-loop etc ON
    double omega3 = asin(2.), msqAv3 = 0., cs3 = 0., dAs3 = 0., 
      dY3 = 0., dYp3 = 0.; mGutGuess = 2.0e16;
    USE_THREE_LOOP_RGE = true;
    USE_TWO_LOOP_THRESHOLD = true;
    MssmSoftsusy threeLoop;
    getCmssmAndOmega(threeLoop, pars, tanb, sgnMu, oneset, mGutGuess, 
		     uni, omega3, msqAv3, cs3, dAs3, dY3, dYp3); 
    //    cout << threeLoop << cs3 << endl; 

    /*    printLineOut(m0, m12, a0, tanb, twoLoop, twoLoopAs, twoLoopMt, 
		 twoLoopMb, threeLoop, 
		 omega2, omegaAs, omegaMt, omegaMb, omega3,
		 msqAv2, msqAvAs, msqAvMt, msqAvMb, msqAv3, cs, csAs, csMt, 
		 csMb, cs3, dAs, dAsAs, dAsMt, dAsMb, dAs3, dY, dYAs, dYMt, 
		 dYMb, dY3, dYp, dYpAs, dYpMt, dYpMb, dYp3);*/
    writeTable(twoLoop, twoLoopAs, twoLoopMt, twoLoopMb, threeLoop, 
	       omega2, omegaAs, omegaMt, omegaMb, omega3,
	       msqAv2, msqAvAs, msqAvMt, msqAvMb, msqAv3, cs, csAs, csMt, 
	       csMb, cs3, dAs, dAsAs, dAsMt, dAsMb, dAs3, dY, dYAs, dYMt, 
	       dYMb, dY3, dYp, dYpAs, dYpMt, dYpMb, dYp3);
  } else 
    if (argc == 5) { /// Scan in m0
      m12 = atof(argv[1]);
	a0 = atof(argv[2]);
	tanb = atof(argv[3]);
	sgnMu = atoi(argv[4]);
    int numPoints = 10;
    double startM0 = 7100., endM0 = 7400.;
    int i; for (i=0; i<=numPoints; i++) {
      m0 = (endM0 - startM0) / double(numPoints) * double(i) + startM0;
    /// Preparation for calculation: set up object and input parameters
    DoubleVector pars(3); 
    pars(1) = m0; pars(2) = m12; pars(3) = a0;
    bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
    
    /// Switch off 3-loop RGEs etc
    double omega2=0., msqAv2 = 0., cs = 0., dAs = 0., 
      dY = 0., dYp = 0.;  
    USE_THREE_LOOP_RGE = true;   USE_TWO_LOOP_THRESHOLD = true;
    MssmSoftsusy twoLoop;
    getCmssmAndOmega(twoLoop, pars, tanb, sgnMu, oneset, mGutGuess, 
		     uni, omega2, msqAv2, cs, dAs, dY, dYp);
    cout << m0 << " " << m12 << " " << a0 << " " << tanb << " " << sgnMu 
	 << " " << omega2 << " " << twoLoop.displayPhys().mh0(1) << endl;      
    }} else if (argc == 1) { /// scan in tan beta
      QedQcd oneset;      ///< See "lowe.h" for default definitions parameters
      
      /// most important Standard Model inputs: 
      /// you may change these and recompile
      double alphasMZ = 0.1187, mtop = 173.2, mbmb = 4.18;
      oneset.setAlpha(ALPHAS, alphasMZ);
      oneset.setPoleMt(mtop);
      oneset.setMbMb(mbmb);
      oneset.toMz();      ///< Runs SM fermion masses to MZ
      
      m0 = 3000.; m12 = 2000.; a0 = -6000.; tanb = 30.; sgnMu = 1;
      double tStart = 2.0, tEnd = 60.0; int numPoints = 40;
      for (int i=0; i<=numPoints; i++) {
	tanb = (tEnd - tStart) * double(i) / double(numPoints) + tStart;
	/*m0 = (mEnd - mStart) * double(i) / double(numPoints) + mStart;
	m12 = m0;
	a0 = -2.0 * m0;*/

	/// Preparation for calculation: set up object and input parameters
	DoubleVector pars(3); 
	pars(1) = m0; pars(2) = m12; pars(3) = a0;
	bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
	
	/// Switch off 3-loop RGEs etc
	double omega2=asin(2.), msqAv2 = 0., cs = 0., dA = 0.,
	  dY = 0., dYp = 0.;  
	USE_THREE_LOOP_RGE = false;   USE_TWO_LOOP_THRESHOLD = false;
	MssmSoftsusy twoLoop;
	getCmssmAndOmega(twoLoop, pars, tanb, sgnMu, oneset, mGutGuess, 
			 uni, omega2, msqAv2, cs, dA, dY, dYp);

	/// Just 2-loop thresholds for strong coupling constant
	double omegaAs = asin(2.), msqAvAs = 0., csAs = 0., dAAs = 0., 
	  dYAs = 0, dYpAs = 0.; mGutGuess = 2.e16; 
	MssmSoftsusy twoLoopAs; 
	twoLoopAs.included_thresholds |= ENABLE_TWO_LOOP_AS_AS_YUK;
	USE_TWO_LOOP_THRESHOLD = true;
	getCmssmAndOmega(twoLoopAs, pars, tanb, sgnMu, oneset, mGutGuess, 
			 uni, omegaAs, msqAvAs, csAs, dAAs, dYAs, dYpAs); 
    
	/// Just 2-loop strong thresholds for mt
	USE_TWO_LOOP_THRESHOLD = false;
	double omegaMt = asin(2.), msqAvMt = 0., csMt = 0., dAMt = 0., 
	  dYMt = 0, dYpMt = 0.; mGutGuess = 2.e16;
	MssmSoftsusy twoLoopMt; 
	twoLoopMt.included_thresholds |= ENABLE_TWO_LOOP_MT_AS;
	USE_TWO_LOOP_THRESHOLD = true;
	getCmssmAndOmega(twoLoopMt, pars, tanb, sgnMu, oneset, mGutGuess, 
			 uni, omegaMt, msqAvMt, csMt, dAMt, 
	  dYMt, dYpMt); 
	
	/// Just 2-loop for mb,mtau
	USE_TWO_LOOP_THRESHOLD = false;
	double omegaMb = asin(2.), msqAvMb = 0., csMb = 0., dAMb = 0., 
	  dYMb = 0, dYpMb = 0.; mGutGuess = 2.e16;
	MssmSoftsusy twoLoopMb; 
	twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MB_AS;
	twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MB_YUK;
	twoLoopMb.included_thresholds |= ENABLE_TWO_LOOP_MTAU_YUK;
	USE_TWO_LOOP_THRESHOLD = true;
	getCmssmAndOmega(twoLoopMb, pars, tanb, sgnMu, oneset, mGutGuess, 
			 uni, omegaMb, msqAvMb, csMb, dAMb, 
	  dYMb, dYpMb); 
	
	/// 3-loop etc ON
	double omega3 = asin(2.), msqAv3 = 0., cs3 = 0., dA3 = 0., 
	  dY3 = 0., dYp3 = 0.; mGutGuess = 2.0e16;
	USE_THREE_LOOP_RGE = true;
	USE_TWO_LOOP_THRESHOLD = true;
	MssmSoftsusy threeLoop;
	getCmssmAndOmega(threeLoop, pars, tanb, sgnMu, oneset, mGutGuess, 
			 uni, omega3, msqAv3, cs3, dA3, 
	  dY3, dYp3); 
	
	printLineOut(m0, m12, a0, tanb, twoLoop, twoLoopAs, twoLoopMt, 
		     twoLoopMb, threeLoop, 
		     omega2, omegaAs, omegaMt, omegaMb, omega3,
		     msqAv2, msqAvAs, msqAvMt, msqAvMb, msqAv3, cs, csAs, csMt, 
		     csMb, cs3, dA, dAAs, dAMt, dAMb, dA3, dY, dYAs, dYMt, 
		     dYMb, dY3, dYp, dYpAs, dYpMt, dYpMb, dYp3);
      }
    }
    }
 
  

  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}

