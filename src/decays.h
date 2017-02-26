/** \file decays.h
   - Project:     SOFTSUSY 
   - Author:      Tom Cridge, Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Code calculates decay modes and prints out an SLHA format file with 
          them in. For R-parity conserving NMSSM/MSSM. See arXiv://1702.?????
*/

#ifndef DECAYS_H
#define DECAYS_H

#include "nmssmsoftsusy.h"
#include "softsusy.h"
#include "physpars.h"
#include "lowe.h"
#include "softpars.h"
#include "softsusy.h"
#include "flavoursoft.h"
#include "susy.h"
#include "particle.h"
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <stdlib.h>
#include <vector>

using namespace std;

/// Approximate accuracy with which 3 body decays are calculated
const double accuracy = 0.01; 

/// Calculate Decays does all decay table calculations and outputs
void calculateDecays(ostream & out, MssmSoftsusy * r,
		     const NmssmSoftsusy & nmssm, 
		     bool nmssmIsIt);


///Dgauss numerical integration integrand
double fdgauss(double x);

/// Change mcpole to alter quark "pole" masses used in h -> qq with QCD
/// corrections 
const double mcpole = 1.40;
/// Change mspole to alter quark "pole" masses used in h -> qq with QCD
/// corrections 
const double mspole = 0.19; 

///Function Declarations used to calculate Partial Widths in decays.cpp
double gluinoamplitudedecay (double m1, double m2, double m3,
			     double alphastrong); 
double gluinoamplitudedecaymix (double m1, double m2, double m3,
				double alphastrong, double squarkmix,
				double theta);
double squarkamplitudedecaygluino (double m1, double m2, double m3,
				   double alphastrong);
double squarkamplitudedecaygluinomix (double m1, double m2, double m3,
				      double alphastrong, double squarkmix,
				      double theta);
double squarkamplitudedecaycharginoW1 (double m1, double m2, double m3,
				       double g, double gamma);
double squarkamplitudedecaycharginoW2 (double m1, double m2, double m3,
				       double g, double gamma);
double squark1amplitudedecaycharginoW1mix
(double m1, double m2, double m3, double g, double gammaL, double gammaR,
 double theta, double beta, double mWboson, double runmt, double runmb,
 double torb); 
double squark1amplitudedecaycharginoW2mix
(double m1, double m2, double m3, double g, double gammaL, double gammaR,
 double theta, double beta, double mWboson, double runmt, double runmb,
 double torb);
double squark2amplitudedecaycharginoW1mix
(double m1, double m2, double m3, double g, double gammaL, double gammaR,
 double theta, double beta, double mWboson, double mup, double mdown,
 double torb);
double squark2amplitudedecaycharginoW2mix
(double m1, double m2, double m3, double g, double gammaL, double gammaR,
 double theta, double beta, double mWboson, double mup, double mdown,
 double torb);
double squarkLamplitudedecayneutralino
(double m1, double m2, double m3, double g, double gprime,
 DoubleMatrix & mixNeut, int neutralino, int uord ); 
  double squarkRamplitudedecayneutralino (double m1, double m2, double m3, double g, double gprime, DoubleMatrix & mixNeut, int neutralino, int uord );
  
  double testrealsquark3amplitudedecayneutralino (double m1, double m2, double m3, double mWboson, double theta, double beta, DoubleMatrix & mixNeut, double g, double gp, double runmt, int squark , int oneortwo,  int neutralino);
  double squark3amplitudedecaysquark3Wboson (double m1, double m2, double m3, double g, double thetat, double thetab, int m1torb, int m1oneortwo, int m3torb, int m3oneortwo);
  double squark3amplitudedecaychargedHiggssquark3 (double m1, double m2, double m3, double g, double mWboson, double beta, double thetat, double thetab, double greekmu, double At, double Ab, double mt, double mb, int t1or2, int b1or2);
  double squark32amplitudedecayneutralHiggssquark3 (double m1, double m2, double m3, double g, double gp, double mWboson, double beta, double alpha, double thetat, double thetab, double greekmu, double At, double Ab, double mt, double mb, int torb, char phi);
  double squark32amplitudedecaysquark3Zboson (double m1, double m2, double m3, double g, double gp, double theta);
  double sleptonamplitudedecayleptonneutralino (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, char LorR, int neutralino);
  double sneutrinoamplitudedecayneutrinoneutralino (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int neutralino);
  double sleptonamplitudedecaychargino (double m1, double m2, double m3, double g, double theta, int chargino);
  double stauamplitudedecaytauneutralino (double m1, double m2, double m3, double g, double gp, double mWboson, DoubleMatrix & mixNeut, double theta, double beta, int oneortwo, int neutralino);
  double stausneutrinoamplitudedecaytauneutrinoneutralino (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int neutralino); 
  double stauamplitudedecaynutauchargino (double m1, double m2, double m3, double g, double mWboson, double theta, double thetaL, double beta, double mtau, int oneortwo, int chargino);
  double stausneutrinoamplitudedecaytauchargino (double m1, double m2, double m3, double g, double mWboson, double beta, double thetaL, double thetaR, int chargino);
  double stauamplitudedecaysnustauHminus (double m1, double m2, double m3, double g, double mWboson, double beta, double thetatau, double mtau, double mu, double Atau, int oneortwo);
  double stauamplitudedecaysnustauWboson (double m1, double m2, double m3, double g, double thetatau, int oneortwo);
  double stau2amplitudedecaystau1Zboson (double m1, double m2, double m3, double g, double gp, double thetatau);
  double stau2amplitudedecaystau1phi (double m1, double m2, double m3, double g, double gp, double thetatau, double beta, double alpha, double mWboson, double mtau, double greekmu, double Atau, char phi);
  double charginoamplitudedecayquarksquarkL (double m1, double m2, double m3, double g, double theta, int chargino);
  double charginoamplitudedecayquarksquarkmix (double m1, double m2, double m3, double g, double theta, double thetaL, double thetaR, double beta, double runmt, double runmb, double mWboson, int chargino, int upordowntypesquark, int oneortwo);
  double charginoamplitudedecayleptonsleptonL (double m1, double m2, double m3, double g, double thetaLorR, int chargino);
  double charginoamplitudedecaysnutautau (double m1, double m2, double m3, double g, double thetaL, double thetaR, double beta, double mWboson, int chargino);
  double charginoamplitudedecaystaunutau (double m1, double m2, double m3, double g, double thetaL, double thetaR, double thetatau, double beta, double mWboson, double mtau, int oneortwo, int chargino);
  double charginoamplitudedecayWbosonneutralino (double m1, double m2, double m3, double g, double thetaL, double thetaR, DoubleMatrix & mixNeut, int chargino, int neutralino);
  double charginoamplitudedecayHminusneutralino (double m1, double m2, double m3, double g, double gp, double thetaL, double thetaR, double beta, DoubleMatrix & mixNeut, int chargino, int neutralino);
  double chargino2amplitudedecaychargino1Zboson (double m1, double m2, double m3, double g, double gp, double thetaL, double thetaR);
  double chargino2amplitudedecaychargino1neutHiggs (double m1, double m2, double m3, double g, double gp, double thetaL, double thetaR, double beta, double alpha, char phi);
  double neutralinoamplitudedecayquarksquarkLorR (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int uordtype , char LorR, int neutralino);
  double neutralinoamplitudedecayleptonsleptonLorR (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, char LorR, int neutralino);
  double neutralinoamplitudedecayneutrinosneutrinoL (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int neutralino);
  double neutralinoamplitudedecaysquark3quarkmix (double m1, double m2, double m3, double mWboson, double theta, double beta, DoubleMatrix & mixNeut, double g, double gp, double runmq, int squark , int oneortwo,  int neutralino);
  double neutralinoamplitudedecaystautau (double m1, double m2, double m3, double mWboson, double theta, double beta, DoubleMatrix & mixNeut, double g, double gp, int oneortwo,  int neutralino);
  double neutralinoamplitudedecaycharginoWboson (double m1, double m2, double m3, double g, double thetaL, double thetaR, DoubleMatrix & mixNeut, int neutralino,  int chargino);
  double neutralinoamplitudedecaycharginoHplus (double m1, double m2, double m3, double g, double gp, double beta, double thetaL, double thetaR, DoubleMatrix & mixNeut, int neutralino,  int chargino);
  double neutralinoamplitudedecayneutralinoZboson (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int ineutralino,  int fneutralino);
  double neutralinoamplitudedecayneutralinoneutHiggs (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, double mixingangle, int ineutralino,  int fneutralino, char phi);

  double higgslorHamplitudedecayquarkantiquark (double m1, double m2, double g, double alpha, double beta, double mWboson, int uord, char lorH, DoubleMatrix & CPEMix, bool nmssmIsIt, bool QCD, double alphas);
  double higgsAamplitudedecayquarkantiquark (double m1, double m2, double g, double beta, double mWboson, int uord, bool QCD, double alphas); 
  double higgsAamplitudedecayquarkantiquarkNMSSM (double m1, double m2, double beta, DoubleMatrix & CPOMix, int uord, int higgs, bool QCD, double alphas);
  double higgsphiamplitudedecayneutralinoneutralino (double m1, double m2, double m3, double g, double tanthetaW, double mixingangle, DoubleMatrix & mixNeut, int ineutralino, int fneutralino, char phi);
  double higgsphiamplitudedecaysamechargino (double m1, double m2, double g, double thetaL, double thetaR, double alpha, double beta, int chargino, char phi); ///this function calls the function higgsphicharginocouplings to calculate the couplings for it
  double higgsphiamplitudedecaydifchargino (double m1, double m2, double m3, double g, double thetaL, double thetaR, double alpha, double beta, char phi);
  double higgshamplitudedecayAA (double m1, double m2, double g, double gp, double alpha, double beta, double mWboson);
  double higgsHamplitudedecayhh (double m1, double m2, double g, double gp, double alpha, double beta, double mWboson);
  double higgsHamplitudedecayAA (double m1, double m2, double g, double gp, double alpha, double beta, double mWboson);
  double higgsHamplitudedecayHplusHminus (double m1, double m2, double g, double gp, double alpha, double beta, double mWboson);
  double higgshamplitudedecayhiggsAZboson (double m1, double m2, double m3, double g, double gp, double alpha, double beta);
  double higgsHamplitudedecayhiggsAZboson (double m1, double m2, double m3, double g, double gp, double alpha, double beta);
  double higgsAamplitudedecayhiggshZboson (double m1, double m2, double m3, double g, double gp, double alpha, double beta);
  double higgsAamplitudedecayhiggsHZboson (double m1, double m2, double m3, double g, double gp, double alpha, double beta);
  double higgshamplitudedecay2squarksamehand (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mupq, double mdownq, int sq);
  double higgshamplitudedecay2squarkdiffhand (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson, double mupq, double mdownq, double greekmu, double Aup, double Adown, int sq);
  double higgsHamplitudedecay2squarksamehand (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mupq, double mdownq, int sq);
  double higgsHamplitudedecay2squarkdiffhand (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson, double mupq, double mdownq, double greekmu, double Aup, double Adown, int sq);
  double higgshamplitudedecay2sleptonsamehand (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mel, int sl);
  double higgshamplitudedecay2sleptondiffhand (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson, double mel, double greekmu, double Aelectron, int sl);
  double higgsHamplitudedecay2sleptonsamehand (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mel, int sl);
  double higgsHamplitudedecay2sleptondiffhand (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson, double mel, double greekmu, double Aelectron, int sl);
  double higgshamplitudedecaystop1stop1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgshamplitudedecaystop2stop2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgshamplitudedecaystop1stop2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgshamplitudedecaysbottom1sbottom1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgshamplitudedecaysbottom2sbottom2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgshamplitudedecaysbottom1sbottom2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgsHamplitudedecaystop1stop1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgsHamplitudedecaystop2stop2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgsHamplitudedecaystop1stop2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgsHamplitudedecaysbottom1sbottom1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgsHamplitudedecaysbottom2sbottom2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgsHamplitudedecaysbottom1sbottom2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta);
  double higgshamplitudedecaystau1stau1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta);
  double higgshamplitudedecaystau2stau2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta);
  double higgshamplitudedecaystau1stau2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta);
  double higgshtestamplitudedecaystau1stau1 (double m1, double m2, double thetatau, double g, double gp, double mWboson, double alpha, double beta, double mtau, double greekmu, double Atau);
  double higgsHamplitudedecaystau1stau1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta);
  double higgsHamplitudedecaystau2stau2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta);
  double higgsHamplitudedecaystau1stau2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta);
  double higgsAamplitudedecaysfermions (double m1, double m2, double m3, double g, double mWboson, double mf, double greekmu, double Asf, double beta, char uord);
  double higgsHplusamplitudedecayquarkantiquark (double m1, double m2, double m3, double g, double mWboson, double beta, DoubleMatrix & VCKM, int quark, int antiquark);
  double higgsHplusamplitudedecayneutralinochargino (double m1, double m2, double m3, double g, double gp, double beta, double thetaL, double thetaR, DoubleMatrix & mixNeut, int neutralino,  int chargino);
  double higgsHplusamplitudedecayneutralinocharginosusyhitway (double m1, double m2, double m3, double g, double gp, double beta, double thetaL, double thetaR, DoubleMatrix & mixNeut, int neutralino,  int chargino);
  double higgsHplusamplitudedecayWbosonhiggsh (double m1, double m2, double m3, double g, double alpha, double beta);
  DoubleVector higgsHplusamplitudedecaysquarksquark (double m1, double m2, double m3, double g, double beta, double mWboson, double mup, double mdown, double greekmu, double Aup, double Adown);
  DoubleVector higgsHplusamplitudedecaysquarksquarkmix (double m1, double m2, double m3, double g, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double thetatop, double thetabottom);
  double higgsesamplitudedecaygammagammatotal(double m1, double g, double gprime, double alphaEmrun, double mWboson, double polemw, double alpha, double beta, double mtop, double mbottom, double mcharm, double mtau, double mHpm, double mstop1, double mstop2, double msbottom1, double msbottom2, double mstau1, double mstau2, double mchargino1, double mchargino2, double thetaL, double thetaR, double thetat, double thetab, double thetatau, double greekmu, double Atop, double Abottom, double Atau, char higgstype);
  double higgsesamplitudedecaygluongluontotal(double m1, double g, double gs, double gprime, double mWboson, double alpha, double beta, double mtop, double mbottom, double mcharm, double mstop1, double mstop2, double msbottom1, double msbottom2, double thetat, double thetab, double greekmu, double Atop, double Abottom, double mstrange, double mscharmL, double mscharmR, double msstrangeL, double msstrangeR, double Acharm, double Astrange, double mup, double mdown, double msupL, double msupR, double msdownL, double msdownR, double Aup, double Adown, char higgstype, bool QCD);
  ///double higgsamplitudedecayVVstar (double m1, double mboson, double g, double gp, double beta, double alpha, char Vtype);
  DoubleVector higgshamplitudedecayVV(double m1, double mWboson, double mZboson, double g, double gp, double alpha, double beta, char Vtype, DoubleMatrix & CPEMix, bool nmssmIsIt);
  DoubleVector higgsHamplitudedecayVV(double m1, double mWboson, double mZboson, double g, double gp, double alpha, double beta, char Vtype, DoubleMatrix & CPEMix, bool nmssmIsIt);
  DoubleVector higgsH3amplitudedecayVVNMSSM(double m1, double mWboson, double mZboson, double g, double gp, double alpha, double beta, char Vtype, DoubleMatrix & CPEMix, bool nmssmIsIt);
  double higgsAamplitudedecayHpmWboson(double m1, double mWboson, double mHpm, double g, double thetaA, int pseudoscalar, bool nmssmIsIt);
  double stop2amplitudedecaystop1CPevenhiggsNMSSM (double mst2, double mst1, double mh, double mt , double thetat, DoubleMatrix & CPEMix, double beta, double mWboson, double g, double gp, double At, double mueff, double lam, int higgs);
  double stop2amplitudedecaystop1CPoddhiggsNMSSM (double mst2, double mst1, double ma, double mt , double thetat, DoubleMatrix & CPOMix, double beta, double mWboson, double g, double At, double mueff, double lam, int higgs);
  double sbottom2amplitudedecaysbottom1CPevenhiggsNMSSM (double msb2, double msb1, double mh, double mb , double thetab, DoubleMatrix & CPEMix, double beta, double mWboson, double g, double gp, double Ab, double mueff, double lam, int higgs);
  double sbottom2amplitudedecaysbottom1CPoddhiggsNMSSM (double msb2, double msb1, double ma, double mb , double thetab, DoubleMatrix & CPOMix, double beta, double mWboson, double g , double Ab, double mueff, double lam, int higgs);
  double stau2amplitudedecaystau1CPevenhiggsNMSSM (double mstau2, double mstau1, double mh, double mtau, double thetatau, DoubleMatrix & CPEMix, double beta, double mWboson, double g, double gp, double Atau, double mueff, double lam, int higgs);
  double stau2amplitudedecaystau1CPoddhiggsNMSSM (double mstau2, double mstau1, double ma, double mtau, double thetatau, DoubleMatrix & CPOMix, double beta, double mWboson, double g , double Atau, double mueff, double lam, int higgs);
  double chargino2amplitudedecaychargino1CPevenhiggsNMSSM (double mchar2, double mchar1, double mh, double g, double lam, double thetaL, double thetaR, DoubleMatrix & CPEMix, int higgs);
  double chargino2amplitudedecaychargino1CPoddhiggsNMSSM (double mchar2, double mchar1, double mA, double g, double lam, double thetaL, double thetaR, DoubleMatrix & CPOMix, int higgs);
  double neutralinoamplitudedecaycharginoWNMSSM (double mneut, double mchar, double mWboson, double g, double thetaL, double thetaR, DoubleMatrix & mixNeut, int neutralino, int chargino);
  double neutralinoamplitudedecayneutralinoZNMSSM (double mneuti, double mneutj, double mZboson, double g, double gp, DoubleMatrix & mixNeut, int neutralinoi, int neutralinoj);
  double neutralinoamplitudecharginoHpmNMSSM (double mneut, double mchar, double mHp, double g, double gp, double beta, double thetaL, double thetaR, double lam, DoubleMatrix & mixNeut, int neutralino, int chargino);
  double neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (double mneuti, double mneutj, double mhiggs, double g, double gp, double lam, double kappa, DoubleMatrix & mixNeut, DoubleMatrix & CPEMix, int neutralinoi, int neutralinoj, int higgs);
  double neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (double mneuti, double mneutj, double ma, double g, double gp, double lam, double kappa, DoubleMatrix & mixNeut, DoubleMatrix & CPOMix, int neuti, int neutj, int higgsa);
  double neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (double mneut, double msf, double mf, double g, double gp, DoubleMatrix & mixNeut, int neut, char type, char LorR);
  double neutralinoamplitudestoptopNMSSM (double mneut, double mst, double mt, double g, double gp, double thetat, double beta, double mWboson, DoubleMatrix & mixNeut, double runmt, int neut, int stop);
  double neutralinoamplitudesbottombottomNMSSM (double mneut, double msb, double mb, double g, double gp, double thetab, double beta, double mWboson, DoubleMatrix & mixNeut, double runmb, int neut, int sbottom);
  double neutralinoamplitudestautauNMSSM (double mneut, double mstau, double mtau, double g, double gp, double thetatau, double beta, double mWboson, DoubleMatrix & mixNeut, double runmtau, int neut, int stau);
  double neutralinoamplitudestauneutrinotauneutrinoNMSSM (double mneut, double mstaunu, double mtaunu, double g, double gp, DoubleMatrix & mixNeut, int neut);

  double squarkamplitudedecayquarkneutralinoNMSSM (double m1, double mq, double mneut, double g, double gp, DoubleMatrix & mixNeut, char uord, char LorR, int neut);
  double sleptonamplitudedecayleptonneutralinoNMSSM (double m1, double ml, double mneut, double g, double gp, DoubleMatrix & mixNeut, char uord, char LorR, int neut);
  double stopamplitudedecaytopneutralinoNMSSM (double m1, double mt, double mneut, double g, double gp, double thetat, DoubleMatrix & mixNeut, double runmt, double mWboson, double beta, int stop, int neut);
  double sbottomamplitudedecaybottomneutralinoNMSSM (double m1, double mb, double mneut, double g, double gp, double thetab, DoubleMatrix & mixNeut, double runmb, double mWboson, double beta, int sbottom, int neut);
  double stauamplitudedecaytauneutralinoNMSSM (double m1, double mtau, double mneut, double g, double gp, double thetatau, DoubleMatrix & mixNeut, double runmtau, double mWboson, double beta, int stau, int neut);
  double charginoiamplitudedecayneutralinojHpmNMSSM (double mchar, double mneut, double mHpm, double g, double gp, double thetaL, double thetaR, double beta, DoubleMatrix & mixNeut, double lam, int chargino, int neut);
  double charginoiamplitudedecayneutralinojWNMSSM (double mchar, double mneut, double mWboson, double g, double gp, double thetaL, double thetaR, DoubleMatrix & mixNeut, int chargino, int neut);
  double HpmamplitudecharginojneutralinoiNMSSM (double mHp, double mchar, double mneut, double g, double gp, double beta, double thetaL, double thetaR, double lam, DoubleMatrix & mixNeut, int neutralino, int chargino);
  double snutauamplitudedecaynutauneutralinoNMSSM (double m1, double mneut, double g, double gp, DoubleMatrix & mixNeut, int neutralino);

  double higgsesamplitudedecayZbosonphotontotal(double m1, double mZboson, double g, double gprime, double alphaEmrun, double polemw, double runmw, double alpha, double beta, double mtop, double mbottom, double mcharm, double mstrange, double mstop1, double mstop2, double msbottom1, double msbottom2, double mHplus, double thetat, double thetab, double greekmu, double Atop, double Abottom, char higgstype);
double gluinoamplitudedecaydgausscharginoqqpbarfirsttwogen (double mgluino, double mchargino, double mquark, double mquarkp, double msqL, double msqpL, double g, double thetaL, double thetaR, double alphas, int charg, bool onetothree);
/// BEN
double gluinoamplitudedecay1to3neutfirsttwogen (double m1, double m2, double m3, double m4, double m5, double g, double gp, DoubleMatrix & mixNeut, double alphas, char uord, int neut, int Nsteps, int adaptive, bool onetothree, double approx);
  double gluinoamplitudedecay1to3charfirsttwogen (double m1, double m2, double m3, double m4, double m5, double m6, double g, double thetaL, double thetaR, double alphas, int charg, int Nsteps, int adaptive, bool onetothree, double approx);
  double gluinoamplitudedecay1to3neutttbar (double m1, double m2, double m3, double m4, double m5, double mw, double g, double gp, double thetat, double beta,double alphas, DoubleMatrix & mixNeut, double runmq, int neutralino, int Nsteps, int adaptive, bool onetothree, double approx);
  double gluinoamplitudedecay1to3neutbbbar (double m1, double m2, double m3, double m4, double m5, double mw, double g, double gp, double thetab, double beta,double alphas, DoubleMatrix & mixNeut, double runmq, int neutralino, int Nsteps, int adaptive, bool onetothree, double approx);
  double gluinoamplitudedecaychartbbar (double m1, double m2, double m3, double m4, double m5, double m6, double m7, double m8, double alphas, double thetat, double thetab, double mw, double g, double gp, double gammaL, double gammaR, double beta, double runmt, double runmb, int chargino, int Nsteps, int adaptive, bool onetothree, double approx);
  double neutralinoamplitudedecayneutffbar (double m1, double m2, double mf, double msf1, double msf2, double mz, double mh, double mH, double mA, double runmf, double mw, double thetaf, double beta, double alpha, double g, double gp, DoubleMatrix & mixNeut, int neutralinoj, int neutralinoi, char qorl, char uord, int Nsteps, int adaptive, bool onetothree, double approx);
  double gluinoamplitudedecaygravitino (double m1, double mgrav, double MPlreduced, int gravonoff, int gluNLSP);
  double squarkamplitudedecaygravitino(double m1, double mgrav, double mquark, double MPlreduced, int gravonoff, int squNLSP);
  double neutralinoamplitudedecayphotongravitino(double m1, double mgrav, double MPlreduced, DoubleMatrix & mixNeut, double g, double gp, int neutralino, int gravonoff, int neutNLSP);
  double neutralinoamplitudedecayZgravitino(double m1, double mZ, double mgrav, double MPlreduced, DoubleMatrix & mixNeut, double g, double gp, double beta, int neutralino, int gravonoff, int neutNLSP);
  double neutralinoamplitudedecayphigravitino(double m1, double mphi, double mgrav, double MPlreduced, DoubleMatrix & mixNeut, double alpha, double beta, int neutralino, int gravonoff, char phi, int neutNLSP);

 ///NMSSM functions (where MSSM functions haven't been recycled)
  double higgsAamplitudedecaysfermionsNMSSM (double m1, double m2, double m3, double g, double mWboson, double mf, double Asf, double beta, double lam, double mueff, DoubleMatrix & CPOMix, char uord, int pseudoscalar);
  double higgsAamplitudedecaysamecharginoNMSSM (double m1, double m2, double g, double thetaL, double thetaR, double alpha, double lam, DoubleMatrix & CPOMix, int chargino, int pseudoscalar);
  double higgsAamplitudedecaydifcharginoNMSSM (double m1, double m2, double m3, double g, double thetaL, double thetaR, double alpha, double lam, DoubleMatrix & CPOMix, int pseudoscalar);
  double higgsAamplitudedecayneutralinoneutralinoNMSSM(double m1, double m2, double m3, double g, double tanthetaW, double lam, double kappa, DoubleMatrix & CPOMix, DoubleMatrix & mixNeut, int ineutralino, int fneutralino, int pseudoscalar);
  double higgsAamplitudedecayhiggshorHZbosonNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double thetaA, DoubleMatrix & CPEMix, int pseudoscalar, int higgs);
  double higgslHamplitudedecayquarkantiquarkNMSSM (double MSSMamplitude, DoubleMatrix & CPEMix, double alpha, int higgs, char uord); ///takes MSSM amplitude and applies relevant NMSSM prefactor to get NMSSM amplitude
  double higgshamplitudedecay2squarksamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sq);
  double higgshamplitudedecay2sleptonsamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sl);
  double higgsHamplitudedecay2squarksamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sq);
  double higgsHamplitudedecay2sleptonsamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sl);
  double higgsH3amplitudedecay2squarksamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sq);
  double higgsH3amplitudedecay2sleptonsamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sl);
  double higgshamplitudedecay2squarkdiffhandNMSSM (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson,double mq, double Aq, double mueff, double lam, DoubleMatrix & CPEMix, int sq, int higgs);
  double higgsphiamplitudedecaysamecharginoNMSSM (double m1, double m2, double g, double thetaL, double thetaR, double lam, DoubleMatrix & CPEMix, int chargino, int higgs);
  double higgsphiamplitudedecaydiffcharginoNMSSM (double m1, double m2, double m3, double g, double thetaL, double thetaR, double lam, DoubleMatrix & CPEMix, int higgs);
  double higgsAamplitudedecaygammagammaNMSSM (double m1, double g, double gprime, double alpha, double mWboson, DoubleMatrix & CPOMix, double beta, double mtop, double mbottom, double mcharm, double mtau, double mch1, double mch2, double thetaL, double thetaR, double lam, int higgs);
  double higgsAamplitudedecaygluongluonNMSSM (double m1, double g, double gs, double alphas, double mWboson, DoubleMatrix & CPOMix, double beta, double mtop, double mbottom, double mcharm, double lam, int higgs, bool QCD);
  double higgsAamplitudedecayZgammaNMSSM (double m1, double g, double gp, double alpha, double mWboson, double mZboson, DoubleMatrix & CPOMix, double beta, double mtop, double mbottom, double mcharm, double mch1, double mch2, double thetaL, double thetaR, double lam, int higgs);
  double higgsCPevenamplitudedecaygammagammaNMSSM(double m1, double mtop, double mbottom, double mcharm, double mtau, double mWboson, double mHpm, double mchar1, double mchar2, double mscharmL, double mscharmR, double mstop1, double mstop2, double msstrangeL, double msstrangeR, double msbottom1, double msbottom2, double msmuonL, double msmuonR, double mstau1, double mstau2, DoubleMatrix & CPEMix, double beta, double g, double gp, double alpha, double thetat, double thetab, double thetatau, double thetaL, double thetaR, double At, double Ab, double Atau, double mu, double mueff, double lam, double kappa, double Alambda, int higgs);
  double higgsCPevenamplitudedecaygluongluonNMSSM(double m1, double mtop, double mbottom, double mcharm, double mWboson, double mscharmL, double mscharmR, double mstop1, double mstop2, double msstrangeL, double msstrangeR, double msbottom1, double msbottom2, double msupL, double msupR, double msdownL, double msdownR, double runmt, double runmb, DoubleMatrix & CPEMix, double beta, double g, double gp, double gs, double alphas, double thetat, double thetab, double thetaL, double thetaR, double At, double Ab, double mu, double mueff, double lam, double kappa, double Alambda, int higgs, bool QCD);
  double higgshamplitudedecayZgammaNMSSM (double m1, double g, double gp, double alpha, double mWboson, double mZboson, double mHpm, DoubleMatrix & CPEMix, double beta, double mtop, double mbottom, double mcharm, double mch1, double mch2, double thetaL, double thetaR, double lam, double kappa, double Alambda, double greekmu, double mueff, int higgs);
  double higgshamplitudedecayneutineutjNMSSM (double m1, double mneuti, double mneutj, double g, double gp, DoubleMatrix & CPEMix, DoubleMatrix & mixNeut, double lam, double kappa, int neuti, int neutj, int higgs); 
  double higgsCPevenamplitudedecayAANMSSM(double m1, double mA1, double mA2, double mWboson, double runmt, double runmb, double g, double gp, double beta, DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, double lam, double kappa, double Alambda, double Akappa, double mueff, int higgs, int pseudoscalar1, int pseudoscalar2);
  double higgsCPevenamplitudedecaypseudoscalarZNMSSM (double m1, double mA, double mZboson, double g, double gp, double beta, DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, int higgs, int pseudoscalar);
  double higgsCPevenamplitudedecayHpHmNMSSM (double m1, double mHpm, double mWboson,  double g, double gp, double mtop, double mbottom, double beta, double lam, double mueff, double kappa, double Alambda, DoubleMatrix & CPEMix, int higgs);
  double higgsCPevenamplitudedecayhhorhHorHH(double m1, double mh1, double mh2, double g, double gp, double runmw, double beta, double lam, double Alambda, double kappa, double Akappa, double mueff, DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, int higgs1, int higgs2);
  double higgsCPevenamplitudedecayhhorhHorHHNMSSM(double m1, double mh1, double mh2, double g, double gp, double mWboson, double mtop, double mbottom, double beta, double lam, double Alambda, double kappa, double Akappa, double mueff, DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, int higgs1, int higgs2, int starthiggs);
  double higgsA2amplitudedecayA1CPevenNMSSM(double m1, double mA1, double mh, double mWboson, double runmt, double runmb, double g, double gp, double beta, DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, double lam, double kappa, double Alambda, double Akappa, double mueff, int higgs);
  double higgsCPevenamplitudedecayWHpmNMSSM (double m1, double mWboson, double mHpm, double beta, double g, DoubleMatrix & CPEMix, int higgs);
  double higgsCPevenamplitudedecaystopistopiNMSSM (double m1, double mstopi, double thetat, double runmt, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double At, double mueff, double lam, int stop, int higgs);
  double higgsCPevenamplitudedecaystopistopjNMSSM (double m1, double mstopi, double mstopj, double thetat, double runmt, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double At, double mueff, double lam, int higgs);
  double higgsCPevenamplitudedecaysbottomisbottomiNMSSM (double m1, double msbottomi, double thetab, double runmb, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double Ab, double mueff, double lam, int sbottom, int higgs);
  double higgsCPevenamplitudedecaysbottomisbottomjNMSSM (double m1, double msbottomi, double msbottomj, double thetab, double runmb, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double Ab, double mueff, double lam, int higgs);
  double higgsCPevenamplitudedecaystauistauiNMSSM (double m1, double mstaui, double thetatau, double runmtau, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double Atau, double mueff, double lam, int stau, int higgs);
  double higgsCPevenamplitudedecaystauistaujNMSSM (double m1, double mstaui, double mstauj, double thetatau, double runmtau, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double Atau, double mueff, double lam, int higgs);
  double higgsCPevenamplitudedecaysnusnuNMSSM (double m1, double msnu, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, int higgs);

  double higgsesamplitudedecaygammagammatotal(double m1, double g, double gprime, double alphaEmrun, double mWboson, double polemw, double alpha, double beta, double mtop, double mbottom, double mcharm, double mtau, double mHpm, double mstop1, double mstop2, double msbottom1, double msbottom2, double mstau1, double mstau2, double mchargino1, double mchargino2, double thetaL, double thetaR, double thetat, double thetab, double thetatau, double greekmu, double Atop, double Abottom, double Atau, char higgstype);
  
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


  DoubleVector foftau(double mpart, double mcomp);  
  
  
  void printRowPDG(ostream & cout, double x);
  void OutputNoPWs(ostream & cout, Particle & P, double BRTol);
  void OutputYesPWs(ostream & cout, Particle & P, double BRTol);

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
																										///only done via usual method (not via a separate integrating function compute_areai) as can't pass a DoubleMatrix & in a function pointer



  ///1 to 3 decay functions calling dgauss
  double gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (double mgluino, double mneutralino, double msqL, double msqR, double mquark, double g, double gp, DoubleMatrix & mixNeut, double alphas, char uord, int neut, bool onetothree);
  double gluinoamplitudedecaydgaussneutralinottbar (double mgluino, double mst1, double mst2, double mneutralino, double mt, double mWboson, double g, double gp, double thetat, double beta, double alphas, DoubleMatrix & mixNeut, double runmt, int neutralino, bool onetothree, char torb);
  // double gluinoamplitudedecaydgaussneutralinobbbarorig(double mgluino, double msb1, double msb2, double mneutralino, double mb, double mWboson, double g, double gp, double thetab, double beta, double alphas, DoubleMatrix & mixNeut, double runmb, int neutralino, bool onetothree);
  double gluinoamplitudedecaydgausschartbbar (double mgluino, double mst1, double mst2, double msb1, double msb2, double mtop, double mbottom, double mchar, double alphas, double thetat, double thetab, double MWBoson, double g, double gp, double gammaL, double gammaR, double beta, double runmt, double runmb, int chargino, bool onetothree);
  // double gluinoamplitudedecaydgausschartbbarorig (double mgluino, double mst1, double mst2, double msb1, double msb2, double mtop, double mbottom, double mchar, double alphas, double thetat, double thetab, double MWBoson, double g, double gp, double gammaL, double gammaR, double beta, double runmt, double runmb, int chargino, bool onetothree);
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


#endif
