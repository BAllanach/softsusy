
/** \file nmssmsoftpars.h
    - Project: Next-to-minimal SOFTSUSY 
    - Author: Ben Allanach, Peter Athron, Lewis Tunstall, Alexander Voigt,
    Anthony Williams
    - Manual: http://arxiv.org/abs/1311.7659
    - Webpage: http://hepforge.cedar.ac.uk/softsusy/
    
    \brief Soft SUSY breaking parameters for NMSSM
*/

#ifndef NMSSMSOFTPARS_H
#define NMSSMSOFTPARS_H

#include <cmath>
//#include <susy.h>
#include "def.h"
#include "linalg.h"
#include "utils.h"
#include "numerics.h"
#include "nmssmsusy.h"
#include "softpars.h"

namespace softsusy {
  
  /// PA: Number of parameters contained in the NMSSM RGEs (5 more soft 
  /// paremeters than the MSSM).
  const static int numSoftParsNmssm = 83 + numNMssmPars;

  /** start of RGE functions **/
  /// Returns beta functions to NMSSM soft parameters (only)
  class SoftParsNmssm; class NmssmSoftsusy;
  void addBetaSoftParsNmssm
  (nmsBrevity & a, const NmssmSusyPars & n, const MssmSusy & m, 
   const MssmSoftPars & msoft, const SoftParsNmssm & nsoft, DoubleVector & dmG, 
   double & dmH1sq, double & dmH2sq, double & dm3sq, double & dmSsq, 
   double & dmSpsq, double & dxiSsq, DoubleMatrix & dmq, DoubleMatrix & dmu, 
   DoubleMatrix & dmd, DoubleMatrix & dme, DoubleMatrix dml, 
   DoubleMatrix & dhu, DoubleMatrix & dhd, DoubleMatrix & dhe, double & dhlam, 
   double & dhkap);
  /** end of RGE functions **/
  
  /** start of SoftParsNmssm object **/
  /// Contains only the NMSSM soft parameters, NOT the MSSM ones
  class SoftParsNmssm {
  private:
    /// NMSSM Trilinear soft terms not present in MSSM
    double alambda, akappa;
    /// NMSSM soft breaking masses 
    /// \f$ m_S^2 |S|^2 + \frac{1}{2}m{'2}_S S^2 + \xiS S + h.c. \f$ 
    /// repectively.
    double mSsq, mSpsq, xiS;
  public:
    /// Default constructor fills object with zeroes
    SoftParsNmssm()
      : alambda(0.0), akappa(0.0), mSsq(0.0), mSpsq(0.0), xiS(0.0) {};
    SoftParsNmssm(const SoftParsNmssm & s)
    : alambda(s.displayTrialambda()), akappa(s.displayTriakappa()), 
      mSsq(s.displayMsSquared()), mSpsq(s.displayMspSquared()), 
      xiS(s.displayXiS()) {};
    SoftParsNmssm(double al, double ak, double ms, double msp, double x)
      : alambda(al), akappa(ak), mSsq(ms), mSpsq(msp), xiS(x) {}
    const SoftParsNmssm & operator=(const SoftParsNmssm & s);
    
    const SoftParsNmssm & displaySoftParsNmssm() const { return *this; }
    //PA: Return trilinear soft mass $a_\lambda S H_u H_d$
    double displayTrialambda() const { return alambda; };
    //PA: Return trillinear soft mass $a_\kappa S S S$
    double displayTriakappa() const { return akappa; };
    /// Return \f$m_{S}^2\f$=mSsq
    double displayMsSquared() const { return mSsq; }; 
    /// Return \f$m_{S'}^2\f$=mSpsq
    double displayMspSquared() const { return mSpsq; }; 
    /// Return xiS i.e \f$xiS S\f$
    double displayXiS() const { return xiS; }; 
    /// PA: Return trilinear soft mass in "SUGRA style" 
    /// \f$ \lambda A_\lambda S H_u H_d\f$
    double displaySoftAlambda(double lam) const;
    /// PA: Return trillinear soft mass in "SUGRA style" 
    /// \f$\kappa A_\kappa S S S\f$
    double displaySoftAkappa(double ka) const;
    /// Whole object output in a doublevector
    const DoubleVector display(const NmssmSusy & n, const MssmSoftPars & s) 
      const;
    /// Increments k and fills the DoubleVector with entries
    void display(DoubleVector & y, int & k) const;

    /// PA: Set trilinear SUSY breaking parameter alambda
    void setTrialambda(double al) { alambda = al; };
    /// PA: Set trilinear SUSY breaking parameter akappa
    void setTriakappa(double ak) { akappa = ak; };
    /// PA: Sets soft scalar mass \f$m_S^2 |S|^2\f$
    void setMsSquared(double f) { mSsq = f; }; 
    /// PA: Sets soft scalar mass \f$m_S^2 S^2\f$
    void setMspSquared(double f) { mSpsq = f; }; 
    /// PA: Sets soft breaking \f$ xiS S\f$
    void setXiS(double f) { xiS = f; }; 
    void setSoftParsNmssm(SoftParsNmssm const & s) { *this = s; }
    void set(const DoubleVector &, NmssmSusyPars & n);    

    /// Returns numerical beta functions of parameters and Brevity
    /*    SoftParsNmssm beta2(nmsBrevity&, const NmssmSusy & betaNmssmSusy,
	  const MssmSoftPars & betaMssmSoftPars) const;*/
    /// Returns double vector containing numerical beta functions of parameters
    //    DoubleVector beta() const;

    /// Returns derivatives of anomalous dimensions of fields with respect to
    /// renormalisation scale in MSSM for: RH leptons, LH leptons, LH quarks, RH
    /// up quarks, RH down quarks, H1 and H2 respectively
    void anomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
			DoubleMatrix & gQQ, DoubleMatrix & gUU,
			DoubleMatrix & gDD,
			double & gH1H1, double & gH2H2, double & gSS) const;
    
    /// adds NMSSM amsb contribution to soft masses.
    void addAmsb(double maux, const NmssmSusy & Nms, MssmSoftPars & m);
    /// Reads in universal boundary conditions at the current scale:
    /// m0, M1/2, A
    void standardSugra(double m0, double m12, double a0, 
		       const NmssmSusyPars & n, const MssmSusy & m,
		       const MssmSoftPars & msoft);
    /// Reads in universal boundary conditions at the current scale: 
    /// m0, M1/2, A, mS, Al and Ak
    void standardsemiSugra(double m0,  double m12, double a0, double Al, 
			   double Ak, const NmssmSusy & n, 
			   const MssmSoftPars & m, double mS=0.0);
    /// Sets all flavour-diagonal SUSY breaking scalar masses to m0
    void universalScalars(double m0);
    /// Sets singlet soft mass to mS and all other flavour-diagonal 
    /// SUSY breaking scalar masses to m0
    void semiuniversalScalars(double m0);
    /// Sets all NMSSM SUSY breaking trilinear couplings to a0
    void universalTrilinears(double a0, const NmssmSusy & n);
    /// Sets singlet trilinears to ak  and al and all SUSY breaking 
    /// trilinear couplings to a0
    void semiuniversalTrilinears(double a0, double al, double ak, 
				 const NmssmSusy & n);
    
    /// Reads in soft SUSY breaking parameters from a file
    void inputSoftParsOnly();
  };
  
  /// Soft SUSY breaking parameters and beta functions.
  /// NB It is NOT clear you need this: delete if possible
  class SoftParsNmssmRGE: public RGE, public SoftParsNmssm, 
			  public NmssmSusy, public MssmSoftPars {
  private:
  public:
    /// Default constructor fills object with zeroes
    SoftParsNmssmRGE();
    /// Constructor fills SUSY conserving parts with another object, all
    /// SUSY breaking parameters set to zero
    //  SoftParsNmssmRGE(const NmssmSusy &);
    /// Constructor sets all parameters equal to those in another object
    SoftParsNmssmRGE(const SoftParsNmssmRGE &);
    /// Constructor sets all parameters equal to those in the base class
    //  SoftParsNmssmRGE(const SoftParsNmssmRGE &);
    /// Sets all parameters equal to those in another object
    const SoftParsNmssmRGE & operator=(const SoftParsNmssmRGE & s);
    /// Constructor sets RPC SUSY parameters to s, gaugino masses to mG,
    /// trilinears to aU, aD, aE for au, ad, ae
    /// trilnears respectively,  \f$m_Q^2\f$=mQl, \f$m_U^2\f$=mUr,
    /// \f$m_D^2\f$=mDr, \f$m_L^2\f$=mLl, \f$m_E^2\f$=mEr, \f$ m_3^2\f$=m3sq,
    /// \f$m_{H_1}^2\f$=mH1sq, \f$m_{H_2}^2\f$=mH2sq, mu parameter, number of
    /// loops=l, and threshold parameter=t
    SoftParsNmssmRGE(const DoubleVector & mG, const
		     DoubleMatrix & aU, const DoubleMatrix & aD, 
		     const DoubleMatrix & aE, const double & alambda,
		     const double & akappa, const DoubleMatrix & mQl, 
		     const DoubleMatrix & mUr, const DoubleMatrix & mDr, 
		     const DoubleMatrix & mLl, const DoubleMatrix & mEr, 
		     double m3sq, double mH1sq, double mH2sq,
		     double mSsq, double mSpsq, double xiS,
		     double mGravitino, int l, int t);
    
    /// Returns whole object as a const
    const SoftParsNmssmRGE & displaySoftPars() const;
        
    /// Return contents of object in a vector: for RG evolution
    //  virtual const DoubleVector display() const;
    
    /// Sets whole thing equal to another object
    void setSoftPars(SoftParsNmssmRGE const &);
    /// Sets total set of RGE parameters equal to elements of a vector
    void set(const DoubleVector & y) { 
      NmssmSusyPars s;
      SoftParsNmssm::set(y, s); 
      setNmssmSusyPars(s);
      MssmSusy::set(y);
    }
    //  void setSusy(const NmssmSusy &);
    
    /// Returns double vector containing numerical beta functions of parameters
    DoubleVector beta() const;
  };
  
  /// Formatted ouput of whole object
  ostream & operator <<(ostream &left, const SoftParsNmssmRGE &s);
  ostream & operator <<(ostream &left, const SoftParsNmssm &s);
  /// Formatted input of whole object
  istream & operator >>(istream &left, SoftParsNmssmRGE &s);
  
 

  inline SoftParsNmssmRGE::SoftParsNmssmRGE()
    : SoftParsNmssm() {}
  
  inline SoftParsNmssmRGE::SoftParsNmssmRGE(const SoftParsNmssmRGE & s)
    : SoftParsNmssm(s.displaySoftParsNmssm()) {}
  

  /*  inline SoftParsNmssmRGE::SoftParsNmssmRGE(const NmssmSusy &s)
    : SoftParsNmssmRGE(s), alambda(0.0), akappa(0.0),
      mSsq(0.0), mSpsq(0.0), xiS(0.0) {
    setLoops(s.displayLoops());
    setThresholds(s.displayThresholds());
    }*/
  
  /*  inline SoftParsNmssm::SoftParsNmssm
  (const DoubleVector & mG, const DoubleMatrix & aU, const DoubleMatrix & aD, 
   const DoubleMatrix & aE, const double & aL, const double & aK, 
   const DoubleMatrix & mQl, const DoubleMatrix & mUr, 
   const DoubleMatrix & mDr, const DoubleMatrix & mLl, 
   const DoubleMatrix & mEr, double m3sqn, double mH1sq,
   double mH2sq, double mSsq, double mSpsq, double xiS,
   double mg, int l, int t)
    : MssmSoftPars(mG, aU, aD, aE, mQl, mUr, mDr, mLl, mEr, m3sqn, mH1sq, 
		   mH2sq, mg), alambda(aL), akappa(aK), mSsq(mSsq), 
		   mSpsq(mSpsq), xiS(xiS) { setNmssmApprox(l, t); }*/

  /*  inline SoftParsNmssmRGE::SoftParsNmssmRGE(const DoubleVector & mG, const
   DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix & aE,
   const double & aL, const double & aK, const DoubleMatrix & mQl,
   const DoubleMatrix & mUr, const DoubleMatrix & mDr, const
   DoubleMatrix & mLl, const DoubleMatrix & mEr, double m3sqn, double mH1sq,
   double mH2sq, double mSsq, double mSpsq, double xiS,
					    double mg, int l, int t)
    : SoftParsNmssm(mG, aU, aD, aE, aL, aK, mQl, mUr, mDr, mLl, mEr, 
		    m3sqn, mH1sq, mH2sq, mSsq, mSpsq, xiS, mg, l, t) {
		    }*/
  
    
} // namespace softsusy

#endif
