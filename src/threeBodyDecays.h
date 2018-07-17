/** \file threeBodyDecays.h
   - Project:     SOFTSUSY 
   - Author:      Tom Cridge, Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Code calculates three-body decay modes and prints out an SLHA format
          file with them in. For R-parity conserving NMSSM/MSSM. 
	  See arXiv://1703.09717
*/

#ifndef THREEBODYDECAYS_H
#define THREEBODYDECAYS_H

#include "nmssmsoftsusy.h"
#include "softsusy.h"
#include "physpars.h"
#include "lowe.h"
#include "def.h"
#include "softpars.h"
#include "softsusy.h"
#include "flavoursoft.h"
#include "susy.h"
#include "particle.h"
#include "twoBodyDecays.h"
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include <complex>

using namespace std;

/// Approximate accuracy with which 3 body decays are calculated
const double accuracy = 0.01; 

/// Change mcpole to alter quark "pole" masses used in h -> qq with QCD
/// corrections 
const double mcpole = 1.40;
/// Change mspole to alter quark "pole" masses used in h -> qq with QCD
/// corrections 
const double mspole = 0.19; 
/// When mass splitting is below this scale, chargino decays involving hadrons
/// are calculated rather than to quarks
const double hadronicScale = 1.5;

void printRowPDG(ostream & cout, double x);
void OutputNoPWs(ostream & cout, Particle & P);
void OutputYesPWs(ostream & cout, Particle & P);

  ///Integral calculating functions for 1->3 decays
  double Zsfintegralsum(double m1, double m2, double msf, double mf, double mz, double min, double max, double Nsteps, int adaptive, double approx);
  double Aintegralsum(double m1, double m2, double mz, double mA, double mf, double min, double max, double Nsteps, int adaptive, double approx);
  double G3integralsum(double m1, double m2, double m3, double mt, double mb, double min, double max, double Nsteps, int adaptive, double approx);
  double G2integralsum(double m1, double m2, double m3, double mt, double mb, double min, double max, double Nsteps, int adaptive, double approx);
double chiprimeintegralsum(double m1, double m2, double m3, double m4, double mq, double min, double max, double Nsteps, int adaptive, double approx);
  double Yintegralsum(double m1, double m2, double m3, double m4, double mq, double min, double max, double Nsteps, int adaptive, double approx);
  double Xintegralsum(double m1, double m2, double m3, double m4, double mq, double min, double max, double Nsteps, int adaptive, double approx);
  double zetaintegralsum(double m1, double m2, double m3, double m4, double mq, double min, double max, double Nsteps, int adaptive, double approx);
  double chiintegralsum(double m1, double m2, double m3, double m4, double mq, double min, double max, double Nsteps, int adaptive, double approx);
  double rhointegralsum(double m1, double m2, double m3, double m4, double mq, double min, double max, double Nsteps, int adaptive, double approx);
  double xsiintegralsum(double m1, double m2, double m3, double m4, double mq, double min, double max, double Nsteps, int adaptive, double approx);
  double phitildaintegralsum(double m1, double m2, double m3, double m4, double mq, double min, double max, double Nsteps, int adaptive, double approx);
  double psitildaintegralsum(double m1, double m2, double m3, double m4, double mq, double min, double max, double Nsteps, int adaptive, double approx);

  double Zintegralsum(double m1, double m2, double mz, double mf, double min, double max, double Nsteps, int adaptive, double approx);
  double G1integralsum(double m1, double m2, double m3, double mq, double min, double max, double Nsteps, int adaptive, double approx);
  
  double G4integralsum(double m1, double m2, double m3, double m4, double mt, double mb, double min, double max, double Nsteps, int adaptive, double approx);
  double G5integralsum(double m1, double m2, double m3, double m4, double mt, double mb, double min, double max, double Nsteps, int adaptive, double approx);
  double G6integralsum(double m1, double m2, double m3, double m4, double mt, double mb, double min, double max, double Nsteps, int adaptive, double approx);
  double G7integralsum(double m1, double m2, double m3, double m4, double mt, double mb, double min, double max, double Nsteps, int adaptive, double approx);
  double G8integralsum(double m1, double m2, double m3, double m4, double mt, double mb, double min, double max, double Nsteps, int adaptive, double approx);
  
  double Jintegralsum(double m1, double m2, double msf, double mphi, double mf, double min, double max, double Nsteps, double approx, int adaptive, int AorhorH);
double ghHintegral(double m1, double m2, double mf, double mh, double mH, double g, double gp, double alpha, int neutralinoi, int neutralinoj, DoubleMatrix & mixNeut, double E); 
double hHintegral (double m1, double m2, double mf, double mh, double mH, double fromE, double toE, double stepE, double g, double gp, double alpha, DoubleMatrix & mixNeut, int neutralinoi, int neutralinoj);
										///Hintegral
																				//only done via usual method (not via a separate integrating function compute_areai) as can't pass a DoubleMatrix & in a function pointer



  ///1 to 3 decay functions calling dgauss
  double gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (double mgluino, double mneutralino, double msqL, double msqR, double mquark, double g, double gp, DoubleMatrix & mixNeut, double alphas, char uord, int neut, bool onetothree);
double gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogenlimit (double mgluino, double mneutralino, double msqL, double msqR, double mquark, double g, double gp, DoubleMatrix & mixNeut, double alphas, char uord, int neut, bool onetothree);
  double gluinoamplitudedecaydgaussneutralinottbar (double mgluino, double mst1, double mst2, double mneutralino, double mt, double mWboson, double g, double gp, double thetat, double beta, double alphas, DoubleMatrix & mixNeut, double runmt, int neutralino, bool onetothree, char torb);
  double gluinoamplitudedecaydgausschartbbar (double mgluino, double mst1, double mst2, double msb1, double msb2, double mtop, double mbottom, double mchar, double alphas, double thetat, double thetab, double MWBoson, double g, double gp, double gammaL, double gammaR, double beta, double runmt, double runmb, int chargino, bool onetothree);
  double neutralinoamplitudedecaydgaussneutralinoffbar (double mneutralinoi, double msf1, double msf2, double mZboson, double mhiggsl, double mhiggsH, double mhiggsA, double mneutralinoj, double mf, double alphas, double thetaq, double mWboson, double g, double gp, double alpha, double beta, double runmq, DoubleMatrix & mixNeut, int ineutralino, int jneutralino, bool onetothree, char uordornuorl);
  double neutralinoamplitudedecaycharginoffprimebar (double mneutralinoi, double msfp1, double msfp2, double msf1, double msf2, double mWboson, double mHP, double mcharginoj, double mfp, double mf, double thetaq, double thetaqp, double g, double gp, double alpha, double beta, double thetaL2, double thetaR2, double runmqp, double runmq, DoubleMatrix & mixNeut, int ineutralino, int jchargino, bool onetothree, char qorl, char norc); 


/// dgauss integrands for numerical integration of phase space in 1->3 decays
/// Make sure you have set useInDGauss before calling this function
double fdgauss(double x);
double gpsitildadgauss(double Et);
double gphitildadgauss(double Et);
double gxsidgauss (double Et);
double grhodgauss (double Et);
double gchidgauss (double Et);
double gzetadgauss (double Et);
double gXdgauss (double Et);
double gYdgauss (double Et);
double gchiprimedgauss (double Et);
double gG1dgauss (double Et);
double gG4dgauss (double Et);
double gG5dgauss(double Et);
double gG6dgauss(double Et);
double gG7dgauss(double Et);
double gG8dgauss(double Et);
double gG2dgauss(double Eb);
double gG3dgauss(double Eb);
double gZdgauss(double E);
double ghHdgauss (double E);
double gAdgauss (double E);
double gZsfdgauss (double s);
double gJdgauss(double s);

double gneutineutjffZ1dgauss(double s);
double gneutineutjffZ2dgauss(double s);
double gneutineutjffZ3dgauss(double s);
double gneutineutjffZ4dgauss(double s);

double gintegralhdgauss(double E);
double gintegralHdgauss(double E);
double gintegralh1dgauss(double E);
double gintegralh2dgauss(double E);
double gintegralh3dgauss(double E);
double gintegralh4dgauss(double E);  
double gintegralH1dgauss(double E);
double gintegralH2dgauss(double E);
  double gintegralH3dgauss(double E);
  double gintegralH4dgauss(double E);
  double gintegralhH1dgauss(double E);
  double gintegralhH2dgauss(double E);
  double gintegralhH3dgauss(double E);
  double gintegralhH4dgauss(double E);
  double gintegralA1dgauss(double E);
  double gintegralA2dgauss(double E);
  double gintegralA3dgauss(double E);
  double gintegralA4dgauss(double E);
  double gintegral1Zsfdgauss(double E);
  double gintegral2Zsfdgauss(double E);
  double gintegral3Zsfdgauss(double E);
  double gintegral4Zsfdgauss(double E);
  double gintegral5Zsfdgauss(double E);
  double gintegral6Zsfdgauss(double E);
  double gintegral7Zsfdgauss(double E);
  double gintegral8Zsfdgauss(double E);
  double gintegral1hsfdgauss(double E);
  double gintegral2hsfdgauss(double E);
  double gintegral3hsfdgauss(double E);
  double gintegral4hsfdgauss(double E);
  double gintegral5hsfdgauss(double E);
  double gintegral6hsfdgauss(double E);
  double gintegral7hsfdgauss(double E);
  double gintegral8hsfdgauss(double E);
  double gintegral1Hsfdgauss(double E);
  double gintegral2Hsfdgauss(double E);
  double gintegral3Hsfdgauss(double E);
  double gintegral4Hsfdgauss(double E);
  double gintegral5Hsfdgauss(double E);
  double gintegral6Hsfdgauss(double E);
  double gintegral7Hsfdgauss(double E);
  double gintegral8Hsfdgauss(double E);
  double gintegral1ZAdgauss(double E);
  double gintegral2ZAdgauss(double E);
  double gintegral3ZAdgauss(double E);
  double gintegral4ZAdgauss(double E);
  double gneutineutjffgA1dgauss(double E);
  double gneutineutjffgA2dgauss(double E);
  double gneutineutjffgA3dgauss(double E);
  double gneutineutjffgA4dgauss(double E);

  double gneuticharjffpW1dgauss(double E);
  double gneuticharjffpW2dgauss(double E);
  double gneuticharjffpHpm1dgauss(double E);
  double gneuticharjffpHpm2dgauss(double E);
  double gneuticharjffpHpm3dgauss(double E);
  double gneuticharjffpHpm4dgauss(double E);
  double gneuticharjffp1sf1sf2dgauss(double E);
  double gneuticharjffp2sf1sf2dgauss(double E);
  double gneuticharjffp3sf1sf2dgauss(double E);
  double gneuticharjffp4sf1sf2dgauss(double E);
  double gneuticharjffp1sfp1sf2dgauss(double E);
  double gneuticharjffp2sfp1sf2dgauss(double E);
  double gneuticharjffp3sfp1sf2dgauss(double E);
  double gneuticharjffp4sfp1sf2dgauss(double E);
  double gneuticharjffp5sfp1sf2dgauss(double E);
  double gneuticharjffp6sfp1sf2dgauss(double E);
  double gneuticharjffp7sfp1sf2dgauss(double E);
  double gneuticharjffp8sfp1sf2dgauss(double E);
  double gneuticharjffp1WHpmdgauss(double E);
  double gneuticharjffp2WHpmdgauss(double E);
  double gneuticharjffp3WHpmdgauss(double E);
  double gneuticharjffp4WHpmdgauss(double E);
  double gneuticharjffpW1Sfpdgauss(double E);
  double gneuticharjffpW2Sfpdgauss(double E);
  double gneuticharjffpW3Sfpdgauss(double E);
  double gneuticharjffpW4Sfpdgauss(double E);
  double gneuticharjffpW5Sfpdgauss(double E);
  double gneuticharjffpW6Sfpdgauss(double E);
  double gneuticharjffpW7Sfpdgauss(double E);
  double gneuticharjffpW8Sfpdgauss(double E);
  double gneuticharjffpHg1dgauss(double E);
  double gneuticharjffpHg2dgauss(double E);
  double gneuticharjffpHg3dgauss(double E);
  double gneuticharjffpHg4dgauss(double E);
  double gneuticharjffp1gsfpdgauss(double E);
  double gneuticharjffp2gsfpdgauss(double E);
  double gneuticharjffp3gsfpdgauss(double E);
  double gneuticharjffp4gsfpdgauss(double E);
  double gneuticharjffp5gsfpdgauss(double E);
  double gneuticharjffp6gsfpdgauss(double E);
  double gneuticharjffp7gsfpdgauss(double E);
  double gneuticharjffp8gsfpdgauss(double E);
  double gneuticharjffp1sfpsfpdgauss(double E);
  double gneuticharjffp21sfpsfpdgauss(double E);
  double gneuticharjffp31sfpsfpdgauss(double E);
  double gneuticharjffp41sfpsfpdgauss(double E);

  DoubleVector squarkmixcharginocouplings (double g, double theta, double beta, double gammaL, double gammaR, double runmt, double runmb, double mWboson, int torb);
  DoubleVector higgsphisamecharginocouplings(double alpha, double beta, double thetaL, double thetaR);
  DoubleVector higgsphidifcharginocouplings (double alpha, double beta, double thetaL, double thetaR);
  DoubleVector higgshsquarksamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mupq, double mdownq);
  DoubleVector higgshsquarkdiffhandcouplings (double mWboson, double g, double alpha, double beta, double mupq, double mdownq, double greekmu, double Aup, double Adown);
  DoubleVector higgsHsquarksamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mupq, double mdownq);
  DoubleVector higgsHsquarkdiffhandcouplings (double mWboson, double g, double alpha, double beta, double mupq, double mdownq, double greekmu, double Aup, double Adown);
  DoubleVector higgshsleptonsamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mel);
  DoubleVector higgshsleptondiffhandcouplings (double mWboson, double g, double alpha, double beta, double mel, double greekmu, double Ae);
  DoubleVector higgsHsleptonsamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mel) ;
  DoubleVector higgsHsleptondiffhandcouplings (double mWboson, double g, double alpha, double beta, double mel, double greekmu, double Ae);
  DoubleVector higgsHplussquarkcouplings (double mWboson, double g, double beta, double mup, double mdown, double greekmu, double Au, double Ad);
  DoubleVector squarkmixcharginocouplings (double g, double theta, double beta, double gammaL, double gammaR, double runmt, double runmb, double mWboson, double mch1, double mch2, int torb);
DoubleVector higgsphisamecharginocouplings (double alpha, double beta, double thetaL, double thetaR);
DoubleVector higgsphidifcharginocouplings (double alpha, double beta, double thetaL, double thetaR);
DoubleVector higgshsquarksamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mupq, double mdownq);
DoubleVector higgshsquarkdiffhandcouplings (double mWboson, double g, double alpha, double beta, double mupq, double mdownq, double greekmu, double Aup, double Adown);
DoubleVector higgsHsquarksamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mupq, double mdownq);
DoubleVector higgsHsquarkdiffhandcouplings (double mWboson, double g, double alpha, double beta, double mupq, double mdownq, double greekmu, double Aup, double Adown);
DoubleVector higgshsleptonsamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mel);
DoubleVector higgshsleptondiffhandcouplings (double mWboson, double g, double alpha, double beta, double mel, double greekmu, double Ae);
DoubleVector higgsHsleptonsamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mel);
DoubleVector higgsHsleptondiffhandcouplings (double mWboson, double g, double alpha, double beta, double mel, double greekmu, double Ae);
DoubleVector higgsHplussquarkcouplings (double mWboson, double g, double beta, double mup, double mdown, double greekmu, double Au, double Ad);
///Function Declarations used to calculate Partial Widths in decays.cpp
DoubleVector foftau(double mpart, double mcomp);  
DoubleVector goftau(double mpart, double mcomp);
double chToN2piInt(double qSq, const DoubleVector & v);
double charginoToNeutralino2pion(const MssmSoftsusy * m);
Complex fofqsq(double qSq);
Complex bw(double mSq, double gamma, double qSq);
  
#endif
