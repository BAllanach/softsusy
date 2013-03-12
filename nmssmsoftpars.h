
/** \file softpars.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   \brief Soft SUSY breaking parameters

*/

#ifndef NMSSMSOFTPARS_H
#define NMSSMSOFTPARS_H

#include <cmath>
//#include <susy.h>
#include <def.h>
#include <linalg.h>
#include <utils.h>
#include <numerics.h>
#include <nmssmsusy.h>
#include <softpars.h>
using namespace softsusy;

/// PA: Number of parameters contained in the NMSSM RGEs (5 more soft paraemeters than the MSSM).
const static int numSoftParsNmssm = 83 + numNMssmPars;

/// Soft SUSY breaking parameters and beta functions.
class SoftParsNmssm: public SoftPars<nMssmSusy, nmsBrevity>
{
private:
  // NMSSM Trilinear soft terms not present in MSSM
  double alambda, akappa;
  // NMSSM soft breaking masses \f$ m_S^2 |S|^2 + \frac{1}{2}m{'2}_S S^2 + \zeta_S S + h.c. \f$ repectively.
  double mSsq, mSpsq, zeta_s;

public:
  using SoftPars<nMssmSusy, nmsBrevity>::displayGaugino;
  using SoftPars<nMssmSusy, nmsBrevity>::setM32;
  using SoftPars<nMssmSusy, nmsBrevity>::setM3Squared;
  using SoftPars<nMssmSusy, nmsBrevity>::displayM3Squared;
  using SoftPars<nMssmSusy, nmsBrevity>::setMh1Squared;
  using SoftPars<nMssmSusy, nmsBrevity>::displayMh1Squared;
  using SoftPars<nMssmSusy, nmsBrevity>::setMh2Squared;
  using SoftPars<nMssmSusy, nmsBrevity>::displayMh2Squared;
  using SoftPars<nMssmSusy, nmsBrevity>::yTildes;
  using SoftPars<nMssmSusy, nmsBrevity>::u1R_PQflip;
  using SoftPars<nMssmSusy, nmsBrevity>::universalGauginos;
  using SoftPars<nMssmSusy, nmsBrevity>::displayTrilinear;
  using SoftPars<nMssmSusy, nmsBrevity>::setTrilinearElement;
  using SoftPars<nMssmSusy, nmsBrevity>::setTrilinearMatrix;
  using SoftPars<nMssmSusy, nmsBrevity>::displaySoftMassSquared;
  using SoftPars<nMssmSusy, nmsBrevity>::displaySoftA;

  /// Default constructor fills object with zeroes
  SoftParsNmssm();
  /// Constructor fills SUSY conserving parts with another object, all
  /// SUSY breaking parameters set to zero
  SoftParsNmssm(const nMssmSusy &);
  /// Constructor sets all parameters equal to those in another object
  SoftParsNmssm(const SoftParsNmssm &);
  /// Sets all parameters equal to those in another object
  const SoftParsNmssm & operator=(const SoftParsNmssm & s);
  /// Constructor sets RPC SUSY parameters to s, gaugino masses to mG,
  /// trilinears to aU, aD, aE for au, ad, ae
  /// trilnears respectively,  \f$m_Q^2\f$=mQl, \f$m_U^2\f$=mUr,
  /// \f$m_D^2\f$=mDr, \f$m_L^2\f$=mLl, \f$m_E^2\f$=mEr, \f$ m_3^2\f$=m3sq,
  /// \f$m_{H_1}^2\f$=mH1sq, \f$m_{H_2}^2\f$=mH2sq, mu parameter, number of
  /// loops=l, and threshold parameter=t
  SoftParsNmssm(const nMssmSusy & s, const DoubleVector & mG, const
	       DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix
	       & aE, const double & alambda, const double & akappa,
               const DoubleMatrix & mQl, const DoubleMatrix & mUr, const
	       DoubleMatrix & mDr, const DoubleMatrix & mLl, const
	       DoubleMatrix & mEr, double m3sq, double mH1sq, double mH2sq,
               double mSsq, double mSpsq, double zeta_s,
	       double mGravitino, double mu, int l, int t);

  /// Returns whole object as a const
  const SoftParsNmssm & displaySoftPars() const;

  //PA: Return trilinear soft mass $a_\lambda S H_u H_d$
  double displayTrialambda() const;
  //PA: Return trillinear soft mass $a_\kappa S S S$
  double displayTriakappa() const;
  //PA: Return trilinear soft mass in "SUGRA style" $ \lambda A_\lambda S H_u H_d$
  double displaySoftAlambda() const;
  //PA: Return trillinear soft mass in "SUGRA style" $\kappa A_\kappa S S S$
  double displaySoftAkappa() const;

  double displayMsSquared() const; ///< Return \f$m_{S}^2\f$=mSsq
  double displayMspSquared() const; ///< Return \f$m_{S'}^2\f$=mSpsq
  double displayZeta_s() const; ///< Return zeta_s i.e \f$zeta_s S\f$
  /// Return contents of object in a vector: for RG evolution
  virtual const DoubleVector display() const;

  /// Sets whole thing equal to another object
  void setSoftPars(SoftParsNmssm const &);
  //PA: Set trilinear SUSY breaking parameter alambda
  void setTrialambda(double al);
  //PA: Set trilinear SUSY breaking parameter akappa
  void setTriakappa(double ak);
  void setMsSquared(double); //PA: <sets soft scalar mass \f$m_S^2 |S|^2\f$
  void setMspSquared(double); //PA: <sets soft scalar mass \f$m_S^2 S^2\f$
  void setZeta_s(double); //PA: <sets soft breaking \f$ zeta_s S\f$
  /// Sets total set of RGE parameters equal to elements of a vector
  void set(const DoubleVector &);
  //  void setSusy(const nMssmSusy &);

  /// Returns double vector containing numerical beta functions of parameters
  DoubleVector beta() const;
  /// Returns numerical beta functions of parameters
  SoftParsNmssm beta2() const;
  /// Returns derivatives of anomalous dimensions of fields with respect to
  /// renormalisation scale in MSSM for: RH leptons, LH leptons, LH quarks, RH
  /// up quarks, RH down quarks, H1 and H2 respectively
  void anomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
		      DoubleMatrix & gQQ, DoubleMatrix & gUU,
		      DoubleMatrix & gDD,
		      double & gH1H1, double & gH2H2, double & gSS) const;

  /// Reads in universal boundary conditions at the current scale:
  /// m0, M1/2, A0, B-parameter and mu
  void universal(double m0,  double m12,  double a0,  double mu,
		 double m3sq);
  /// Reads in universal boundary conditions at the current scale: m0, M1/2, A0
  void standardSugra(double m0,  double m12, double a0);
  /// Reads in universal boundary conditions at the current scale: m0, M1/2, A, mS, Al and Ak
  void standardsemiSugra(double m0,  double mS, double m12, double a0, double Al, double Ak);
/// Sets all flavour-diagonal SUSY breaking scalar masses to m0
  void universalScalars(double m0);
  /// Sets singet soft mass to mS and all other flavour-diagonal SUSY breaking scalar masses to m0
  void semiuniversalScalars(double m0, double mS);
  /// Sets all SUSY breaking trilinear couplings to a0
  void universalTrilinears(double a0);
 /// Sets singlet trilinears to ak  and al and all SUSY breaking trilinear couplings to a0
  void semiuniversalTrilinears(double a0, double al, double ak);
  /// Boundary conditions to be applied at messenger scale for Gauge mediated
  /// SUSY breaking (see hep-ph/9703211 for example), n5 is the number of
  /// 5-plets, mMess is the messenger scale and lambda is the GMSB scale
  void minimalGmsb(int n5, double lambda, double mMess, double cgrav);

  /// Reads in soft SUSY breaking parameters from a file
  void inputSoftParsOnly();
};

/// Formatted ouput of whole object
ostream & operator <<(ostream &left, const SoftParsNmssm &s);
/// Formatted input of whole object
istream & operator >>(istream &left, SoftParsNmssm &s);

inline SoftParsNmssm::SoftParsNmssm()
  : SoftPars<nMssmSusy, nmsBrevity>(), alambda(0.0), akappa(0.0),
  mSsq(0.0), mSpsq(0.0), zeta_s(0.0) {
  setPars(numSoftParsNmssm);
  setMu(0.0);
  setLoops(0);
  setThresholds(0);
}

inline SoftParsNmssm::SoftParsNmssm(const SoftParsNmssm & s)
  : SoftPars<nMssmSusy, nmsBrevity>(s),
  alambda(s.displayTrialambda()), akappa(displayTriakappa()),
  mSsq(s.displayMsSquared()), mSpsq(s.displayMspSquared()),
  zeta_s(s.displayZeta_s()) {

  setPars(numSoftParsNmssm);
  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
}

inline SoftParsNmssm::SoftParsNmssm(const nMssmSusy &s)
    : SoftPars<nMssmSusy, nmsBrevity>(s), alambda(0.0), akappa(0.0),
    mSsq(0.0), mSpsq(0.0), zeta_s(0.0) {
      setPars(numSoftParsNmssm);
      setMu(s.displayMu());
      setLoops(s.displayLoops());
      setThresholds(s.displayThresholds());
}

inline SoftParsNmssm::SoftParsNmssm
(const nMssmSusy & s, const DoubleVector & mG, const
 DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix & aE,
 const double & aL, const double & aK, const DoubleMatrix & mQl,
 const DoubleMatrix & mUr, const DoubleMatrix & mDr, const
 DoubleMatrix & mLl, const DoubleMatrix & mEr, double m3sqn, double mH1sq,
 double mH2sq, double mSsq, double mSpsq, double zeta_s,
 double mg, double mu, int l, int t)
   : SoftPars<nMssmSusy, nmsBrevity>(s, mG, aU, aD, aE, mQl, mUr, mDr, mLl, mEr, m3sqn, mH1sq, mH2sq, mg, mu, l, t)
   , alambda(aL), akappa(aK), mSsq(mSsq), mSpsq(mSpsq), zeta_s(zeta_s) {
      setPars(numSoftParsNmssm);
      setMu(mu);
      setLoops(l);
      setThresholds(t);
}

inline const SoftParsNmssm & SoftParsNmssm::displaySoftPars() const { return *this; }

inline double SoftParsNmssm::displayMsSquared() const { return mSsq; }

inline double SoftParsNmssm::displayMspSquared() const { return mSpsq; }

inline double SoftParsNmssm::displayZeta_s() const { return zeta_s; }

inline double SoftParsNmssm::displayTrialambda() const { return alambda; }

inline double SoftParsNmssm::displayTriakappa() const { return akappa; }

inline void SoftParsNmssm::setMsSquared(double f) { mSsq = f; }
inline void SoftParsNmssm::setMspSquared(double f) { mSpsq = f; }
inline void SoftParsNmssm::setZeta_s(double f) { zeta_s = f; }
inline void SoftParsNmssm::setSoftPars(SoftParsNmssm const & s) { *this = s; }

#endif
