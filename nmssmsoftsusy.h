/** \file nmssmsoftsusy.cpp
    Project: NMSSMSOFTSUSY 
    Author: Ben Allanach, Peter Athron, Lewis Tunstall, Alexander Voigt 
    Manual: TBW
    Webpage:  https://github.com/Expander/softsusy.git 
*/


#ifndef NMSSMSOFTSUSY_H
#define NMSSMSOFTSUSY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <def.h>
#include <utils.h>
#include <numerics.h>
#include <physpars.h> 
#include <lowe.h>
#include <nmssmsoftpars.h>
#include <softsusy.h>
#include "mssmUtils.h"
//#include <nmssm2loop.h>
using namespace softsusy;
using namespace std;
/* class NmssmSoftsusy;  */
/* std::istream & operator >>(std::istream &left, NmssmSoftsusy &s);\ */

/// Contains all supersymmetric NMSSM parameters, incorporating R_p NMSSM
class NmssmSoftsusy: public Softsusy<SoftParsNmssm> {
private:
  
   double tSOVSMs;  ///< New Nmssm DRbar tadpole(MSusy): incl 2 loops
   double tSOVSMs1loop; ///<New Nmssm DRbar tadpole(MSusy): excl 2 loops

public:
//  void (*boundaryCondition)(NmssmSoftsusy &, const DoubleVector &);
  /// Default constructor fills object with zeroes
  NmssmSoftsusy();
  /// Constructor sets SUSY parameters only from another object
  NmssmSoftsusy(const NmssmSusy &);
  /// Constructor copies another object
  NmssmSoftsusy(const NmssmSoftsusy &);
  /// Sets all parameters from s, sp, mu is the mu superpotential parameter, l
  /// is the number of loops used for RG evolution, t is the thresholds
  /// accuracy parameter, mg is the gravitino mass, hv is the Higgs VEV
  /// parameter.
  NmssmSoftsusy(const SoftParsNmssm & s, const sPhysical & sp, double mu, int l, int t, double hv);
  /// Set all data in the object equal to another
  const NmssmSoftsusy & operator=(const NmssmSoftsusy & s);
  
  double displayTadpoleSMs() const; ///< displays t_s/s tadpole
  double displayTadpoleSMs1loop() const; ///< displays t_2/v_s tadpole@1 loop
  //PA: obtains NMSSM H1-sfermion-sfermion couplings 
  //for 3rd generation sfermions
  void  H1SfSfCouplings(DoubleMatrix & lTS1Lr, DoubleMatrix & lBS1Lr, 
			DoubleMatrix  & lTauS1Lr, double gmzOcthW, 
			double mu,  double cosb, double v1) const;
//PA: obtains NMSSM H2-sfermion-sfermion couplings 
  //for 3rd generation sfermions
  void H2SfSfCouplings(DoubleMatrix & lTS1Lr, DoubleMatrix & lBS1Lr, 
		       DoubleMatrix  & lTauS1Lr, double gmzOcthW, 
		       double mu,  double sinb) const;
  //PA: obtains NMSSM S-sfermion-sfermion couplings 
  //for 3rd generation sfermion
  void SSfSfCouplings(DoubleMatrix & lTS3Lr, DoubleMatrix & lBS3Lr, 
		      DoubleMatrix & lTauS3Lr, double lam) const;
  //PA:routine to calculate sfermiom contributions to (16 \pi^2) ts / s 
  double doCalcTadSSfermions(DoubleMatrix lTS3Lr, DoubleMatrix lBS3Lr, 
			     DoubleMatrix lTauS3Lr, double q, double s) const;
  //PA: for loop corrections, helps adding Higgs corrections in a tidy way
  void assignHiggs(DoubleVector & higgsm, DoubleVector & higgsa, 
                     DoubleVector & higgsc) const;
  //PA: NMSSM routine to obtain Higgs loop parts of (16 \pi^2) t1/v1
//Includes goldstone bosons. 
  double doCalcTad1Higgs(double q, double costhDRbar, double g, double tanb) const;
  //PA: NMSSM routine to obtain Higgs loop parts of (16 \pi^2) t2/v2
  //Includes goldstone bosons.
  double doCalcTad2Higgs(double q, double costhDRbar, double g, double tanb) const;
  //PA: NMSSM routine to obtain Higgs loop parts of (16 \pi^2) ts/s
  //Includes goldstone bosons.
  double doCalcTadSHiggs(double q, double tanb) const;
  //PA: NMSSM routine to obtain Neutralino loop parts of (16 \pi^2) t1/v1
  double doCalcTad1Neutralinos(double q, double costhDRbar, double g, 
                               double tanb) const;
  //PA: NMSSM routine to obtain Neutralino loop parts of (16 \pi^2) t2/v2
  double doCalcTad2Neutralinos(double q, double costhDRbar, 
			     double g, double sinb) const;
  //PA: NMSSM routine to obtain Neutralino loop parts of (16 \pi^2) tS/s
  double doCalcTadSNeutralinos(double q,  double lam, double kap) const;
  //PA: NMSSM routine to obtain Chargino loop parts of (16 \pi^2) tS/s
  double doCalcTadSCharginos(double q,  double lam) const;
/// Does the calculation of one-loop pieces of \f$ t_1 / v_1 \f$ 
  virtual double doCalcTadpole1oneLoop(double mt, double sinthDRbar) const;
  /// Does the calculation of one-loop pieces of \f$ t_2 / v_2 \f$ 
  virtual double doCalcTadpole2oneLoop(double mt, double sinthDRbar) const;
  /// Does the calculation of one-loop pieces of \f$ t_s / s \f$ 
  virtual double doCalcTadpoleSoneLoop(double mt, double sinthDRbar) const;

/// Calculates and sets the one-loop pieces of \f$ t_1 / v_1 \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  virtual void calcTadpole1Ms1loop(double mt, double sinthDRbar);
  /// Calculates then sets the one-loop pieces of \f$ t_2 / v_2 \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  virtual void calcTadpole2Ms1loop(double mt, double sinthDRbar);
/// Calculates then sets the one-loop pieces of \f$ t_s / s \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  /// Calculates then sets the one-loop pieces of \f$ t_s / s \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  virtual void calcTadpoleSMs1loop(double mt, double sinthDRbar);
 



  /// Organises tree-level calculation of all sparticle masses and mixings
  virtual void calcDrBarPars();
  
  /// Returns tree-level up squark mass matrix in "mass" for NMSSM.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mtrun=DR bar top
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
void treeUpSquark(DoubleMatrix & mass, double mtrun, 
                    double pizztMS, double sinthDRbarMS, 
                    int family);

/// Returns tree-level down squark mass matrix in "mass" for NMSSM.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mbrun=DR bar bottom
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
  void treeDownSquark(DoubleMatrix & mass, double mbrun, double pizztMS, 
		double sinthDRbarMS, int family);
  /// Returns tree-level down squark mass matrix in "mass" for NMSSM.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mTrun=DR bar tau
  /// mass, family=generation of slepton, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w
  void treeChargedSlepton(DoubleMatrix & mass, double mTrun, double pizztMS, 
		double sinthDRbarMS, int family);
  /// LCT: new routine to allocate NMSSM chargino masses, 
  //Returns tree-level chargino mass matrix in the NMSSM 
  void calcDrBarCharginos(DoubleMatrix & mass, double beta, double mw);
  
/// LCT: new routine for NMSSM neutralino masses, 
  //Returns tree-level Neutralino mass matrix in the NMSSM 
  void calcDrBarNeutralinos(DoubleMatrix & mass, double beta, double mz, double mw, double sinthDRbar);
  //PA:  calls  calcDrBarCharginos and calcDrBarNeutralinos, 
  // performs diagonalisation and fills eg with appropriate values.
  void calcDrBarGauginos(double beta, double mw, double mz, double sinth, drBarPars & eg);

  void calcDrBarHiggs(double beta, double mz2, double mw2, double sinthDRbar, drBarPars & eg);
 
  //PA:: fixes The CP odd mixing matrix with the conventions 
  // Degrassi and Slavich arXiv:0907.4682
  void DegrassiSlavicMix(DoubleMatrix & P) const; 
  /// Calculates Higgs contribution to the transverse part of Z self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piZZTHiggs(double p, double q,	double thetaWDRbar) const;
  /// Calculates neutralino contrib. to the transverse part of Z self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piZZTNeutralinos(double p, double Q, double thetaWDRbar) const;
  /// Calculates transverse part of Z self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piZZT(double p, double Q, bool usePoleMt = false) const;
   /// Calculates Higgs contribution to the transverse part of W self-energy: 
  //for p=external momentum, q=renormalisation scale
  virtual double piWWTHiggs(double p, double q, double thetaWDRbar) const;
  /// Calculates gaugino contr. to transverse part of W self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piWWTgauginos(double p, double Q, double thetaWDRbar) const;
  /// Calculates transverse part of W self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piWWT(double p, double Q, bool usePoleMt = false) const;
  //PA: Obtains trilnear couplings of s1-higgs-higgs for use in loop functions
   void getS1HiggsTriCoup(DoubleMatrix & sss1, DoubleMatrix & pps1, DoubleMatrix & hphps1, double thetaWDRbar) const; 
  //PA: Obtains trilnear couplings of s2-higgs-higgs for use in loop functions
  void getS2HiggsTriCoup(DoubleMatrix & sss2, DoubleMatrix & pps2, DoubleMatrix & hphps2, double thetaWDRbar) const; 
  //PA: Obtains trilnear couplings of s2-higgs-higgs for use in loop functions
  void getS3HiggsTriCoup(DoubleMatrix & sss3, DoubleMatrix & pps3, DoubleMatrix & hphps3) const; 
  //PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s1Higgs(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s2Higgs(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis2s2Higgs(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s3Higgs(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis2s3Higgs(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis3s3Higgs(double p, double q) const;
  void getS1NeutralinoCoup(ComplexMatrix & aChi, ComplexMatrix & bChi) const;
  //PA:Calculates (16 Pi^2) times Neutralino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  void getS2NeutralinoCoup(ComplexMatrix & aChi, ComplexMatrix & bChi) const;
  //PA:Calculates (16 Pi^2) times Neutralino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  void getS3NeutralinoCoup(ComplexMatrix & aChi, ComplexMatrix & bChi) const;
  //PA:Calculates (16 Pi^2) times Neutralino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s1Neutralinos(double p, double q) const;
  //PA:Calculates (16 Pi^2) times Neutralino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s2Neutralinos(double p, double q) const;
  //PA:Calculates (16 Pi^2) times Neutralino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis2s2Neutralinos(double p, double q) const;
  //PA:Calculates (16 Pi^2) times Neutralino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s3Neutralinos(double p, double q) const;
   //PA:Calculates (16 Pi^2) times Neutralino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis2s3Neutralinos(double p, double q) const;
  //PA:Calculates (16 Pi^2) times Neutralino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis3s3Neutralinos(double p, double q) const;
  //PA:Calculates (16 Pi^2) times Chargino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s3Charginos(double p, double q) const; 
  double pis2s3Charginos(double p, double q) const;
  double pis3s3Charginos(double p, double q) const;
  
   //PA:Calculates (16 Pi^2) times sfermion contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s3Sfermions(double p, double q, DoubleMatrix ls1tt,  
                         DoubleMatrix ls1bb,  DoubleMatrix ls1tautau, 
                         DoubleMatrix ls3tt,  DoubleMatrix ls3bb,  
                         DoubleMatrix ls3tautau) const;
  double pis2s3Sfermions(double p, double q, DoubleMatrix ls2tt,
                         DoubleMatrix ls2bb,  DoubleMatrix ls2tautau, 
                         DoubleMatrix ls3tt,  DoubleMatrix ls3bb,  
                         DoubleMatrix ls3tautau) const;
  double pis3s3Sfermions(double p, double q, DoubleMatrix ls3tt,  
                         DoubleMatrix ls3bb, DoubleMatrix ls3tautau) const;
     
  /// Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis1s1(double p, double q) const;
  /// Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis1s2(double p, double q) const;
  /// Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis2s2(double p, double q) const;
  // Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis1s3(double p, double q) const;
  // Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis2s3(double p, double q) const;
  // Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis3s3(double p, double q) const;

 //PA: A print method used in development.  I find it useful and easier to read than couting the normal display function or calling printlong etc.    
  void printall();

};
inline NmssmSoftsusy::NmssmSoftsusy()
                     : Softsusy<SoftParsNmssm>(), tSOVSMs(0.0), tSOVSMs1loop(0.0)  {}


inline NmssmSoftsusy::NmssmSoftsusy(const NmssmSoftsusy & s)
		     : Softsusy<SoftParsNmssm>(s),
                     tSOVSMs(s.tSOVSMs), tSOVSMs1loop(s.tSOVSMs1loop)  {
  setPars(121);   
}


inline NmssmSoftsusy::NmssmSoftsusy(const NmssmSusy &s)
                     : Softsusy<SoftParsNmssm>(s),tSOVSMs(0.0), tSOVSMs1loop(0.0)  {
  setPars(121);
}


inline NmssmSoftsusy::NmssmSoftsusy
(const SoftParsNmssm & s, const sPhysical & sp, double mu, int l, int t,
 double hv): Softsusy<SoftParsNmssm>(s, sp, mu, l, t, hv),tSOVSMs(0.0), tSOVSMs1loop(0.0)  {
 setPars(121);
  
}


inline double NmssmSoftsusy::displayTadpoleSMs() const {
  return tSOVSMs; 
}

double NmssmSoftsusy::displayTadpoleSMs1loop() const {
  return tSOVSMs1loop; 
}

#endif
