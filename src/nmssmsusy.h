/** \file nmssmsusy.h
    - Project:     Next-to-Minimal SOFTSUSY
    - Author: Ben Allanach, Peter Athron, Lewis Tunstall, Alexander Voigt,
    Anthony Williams
    - Manual: http://arxiv.org/abs/1311.7659
    - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
    
    \brief NmssmSusy NMSSM SUSY couplings and tan beta.
    NmssmSusy contains all NMSSM SUSY couplings and tan beta, as
    well as their beta functions
*/

#ifndef NMSSMSUSY_H
#define NMSSMSUSY_H

#include "susy.h"

namespace softsusy {  
  /// with the 5 new nmssm parameters we now have 38 in total
  const static int numNMssmPars = numSusyPars + 5;
  
  /// Contains data needed in beta function calculation to make it faster
  struct nmsBrevity : public sBrevity {
    double lsq, ksq, l4, k4;
    
    /// Constructor fills struct with zeroes by default
    nmsBrevity();
    /// Constructor sets struct to be equal to another
    nmsBrevity(const nmsBrevity&);
    /// Sets struct to be equal to another
    const nmsBrevity & operator=(const nmsBrevity &);
    /// Fill MSSM part of the struct
    const nmsBrevity & operator=(const sBrevity &);
    
    /// Calculates all the data given Yukawa matrices yu, yd and ye and gauge
    /// couplings g, lambda and kappa
    void calculate(const DoubleMatrix & yu, const DoubleMatrix & yd, const
		   DoubleMatrix & ye, const DoubleVector & g,
		   const double & lam, const double & kap);
  }; 
  
  inline nmsBrevity::nmsBrevity()
    : sBrevity()
    , lsq(0.0), ksq(0.0), l4(0.0), k4(0.0)
  {}
  
  inline nmsBrevity::nmsBrevity(const nmsBrevity &s)
    : sBrevity(s)
    , lsq(s.lsq), ksq(s.ksq), l4(s.l4), k4(s.k4)
  {}
  /** end of nmsBrevity **/

  /** Start of functions used in the RGEs **/
  class NmssmSusyPars; class NmssmSusy;
  /// Adds NMSSM one loop pieces onto the wave function renomalisation
  void addOneLpAnomNmssm(double & gH1H1, double & gH2H2, double & gSS, 
			 double lambda, double kappa);
  /// Adds NMSSM two loop pieces onto the wave function renomalisation. a
  /// should be calculated before it is called. 
  void addTwoLpAnomNmssm(DoubleMatrix & gEE, DoubleMatrix & gLL,
			 DoubleMatrix & gQQ, DoubleMatrix & gDD,
			 DoubleMatrix & gUU, double & gH1H1, double &
			 gH2H2, double & gSS, nmsBrevity & a);
  /// Organises the calculation of the NMSSM+MSSM pieces of the anomalous 
  /// dimensions
  void anomalousDimensionNmssmSusy(const MssmSusy & s, const NmssmSusyPars & n,
				   DoubleMatrix & gEE, 
				   DoubleMatrix & gLL,DoubleMatrix & gQQ, 
				   DoubleMatrix & gUU,
				   DoubleMatrix & gDD, DoubleVector & dg, 
				   double & gH1H1, double & gH2H2,
				   double & gSS, nmsBrevity & a);
  /// beta functions
  NmssmSusy betaNmssmSusy(nmsBrevity & a, const MssmSusy & s, 
			  const NmssmSusyPars & n);
  /// Outputs beta function coefficients for MSSM gauge coupling evolution in
  /// arguments.
  void nmsetBetas(DoubleMatrix &, DoubleVector  &, DoubleVector  &, DoubleVector
		  &, DoubleVector  &, DoubleVector  & );
  /// Outputs beta function coefficients for lambda.
  void setBetaLambda(DoubleVector&);
  /** end of RGE functions **/

  /// Contains NMSSM-only part of RPC SUSY parameters
  class NmssmSusyPars {
  private:
    /// new nmssm parameters, lambda, kappa appearing as superpotential
    /// terms, lambda S H_u H_d and \frac{1}{3} kappa S^3 and sVev is
    /// the vev of the singlet superfield defined by sVev = \sqrt{2}<S>.
    double lambda, kappa, sVev;
    /// new parameters appearing the general MSSM, but not the ususal
    /// form of the NMSSM which solves the mu-problem nor the MSSM.
    /// They appear in superpotential terms xiF S and \frac{1}{2} mupr
    /// S^2
    double xiF, mupr;
  public:
    NmssmSusyPars(); ///< Constructor fills object with zeroes by default
    NmssmSusyPars(const NmssmSusyPars & );
    NmssmSusyPars(double sv, double lambda, double kappa, double xiF, 
		  double mupr);
    const NmssmSusyPars & operator=(const NmssmSusyPars & s);
    virtual ~NmssmSusyPars() {};

    /// sets DRbar running singlet vev.
    void setSvev(double s) { sVev = s; };
    /// sets the \lambda S H_u H_d coupling
    void setLambda(double l) { lambda = l; };
    /// sets the \frac{1}{3} \kappa S^3 coupling
    void setKappa(double k) { kappa = k; };
    /// sets the \frac{1}{2} mupr S^2 coupling
    void setMupr(double m) { mupr = m; };
    /// sets the xiF S coupling
    void setXiF(double x) { xiF = x; };
    /// sets the whole object
    void setNmssmSusyPars(const NmssmSusyPars &s);
    void setNmssmSusyPars(const DoubleVector & y);
    void set(const DoubleVector & y, int & k);

    /// returns DRbar running Singlet Higgs vev
    double displaySvev() const { return sVev; };
    /// returns superpotential parameter lambda
    double displayLambda() const { return lambda; };
    /// returns superpotential parameter lambda
    double displayKappa() const { return kappa; };
    /// returns mupr superpotential parameter
    double displayMupr() const { return mupr; };
    /// returns xiF superpotential parameter
    double displayXiF() const { return xiF; };
    /// returns whole object
    const NmssmSusyPars & displayNmssmSusyPars() const { return *this; };
    /// whole object in a vector y, starting at index k
    void display(DoubleVector & y, int & k) const; 

    /// Outputs one-loop anomalous dimensions gii given matrix inputs.
    /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
    /// respectively. Note that we use the convention (for matrices in terms of
    /// gamma's): gamma^Li_Lj = M_ij for LH fields and
    /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
    /// conjugates of the RH fields). a should already be defined.
    void getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		      DoubleMatrix & gQQ, DoubleMatrix & gDD,
		      DoubleMatrix & gUU, double & gH1H1, double &
		      gH2H2, double & gSS, nmsBrevity & a, int loops) const;
    /// Outputs two-loop anomlous dimensions gii given matrix inputs.
    /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
    /// respectively. Note that we use the convention (for matrices in terms of
    /// gamma's): gamma^Li_Lj = M_ij for LH fields and
    /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
    /// conjugates of the RH fields). a should already be defined.
    void getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		      DoubleMatrix & gQQ, DoubleMatrix & gDD,
		      DoubleMatrix & gUU, double & gH1H1, double &
		      gH2H2, double & gSS, nmsBrevity & a, int loops) const;
    /// Outputs wave function renormalisation for SUSY parameters and gauge beta
    /// functions up to 2 loops. Also calculates and outputs a.
    /// IO parameters: RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1
    /// and H2 respectively.
    /// g^Li_Lj = m_{ij} for LH fields
    /// g^Ei_Ej = m_{ji} for RH fields
    void anomalousDimension
    (const MssmSusy & s, DoubleMatrix & gEE, DoubleMatrix & gLL,
     DoubleMatrix & gQQ, DoubleMatrix & gUU, DoubleMatrix & gDD, 
     DoubleVector & dg, double & gH1H1, double & gH2H2, double & gSS, 
     nmsBrevity & a, int loops)  const {
      anomalousDimensionNmssmSusy(s, displayNmssmSusyPars(), gEE, gLL, gQQ, gUU,
				  gDD, dg, gH1H1, gH2H2, gSS, a);
    };
  }; ///< end of NmssmSusyPars

  /// Contains all SUSY RPC-NMSSM parameters but is not an RGE object
  class NmssmSusy: public MssmSusy, public NmssmSusyPars {
  private:
    Approx nmssmSusyApprox; ///< Number of loops and thresholds
  public:
    NmssmSusy(); ///< Constructor fills object with zeroes by default
    /// Constructor sets object to be equal to another
    NmssmSusy(const NmssmSusy &);
    NmssmSusy(const MssmSusy & m, const NmssmSusyPars & nsp);
    /// PA: Constructor given Yukawa matrices u,d,e, gauge couplings v, mu
    /// parameter=m, tan beta=tb, lambda, kappa, mupr, xiF,
    // number of loops in
    /// RG evolution l and thresholds parameter t
    NmssmSusy(const DoubleMatrix & u,
	      const DoubleMatrix & d, 
	      const DoubleMatrix & e, const DoubleVector & v, double m, 
	      double tb,  double hv, int mix, int l, int t, double sv,
	      double lambda, double kappa, double xiF, 
	      double mupr);
    NmssmSusy(double lambda, double kappa, double sv, double xiF, 
	      double mupr);
    //    NmssmSusy(const NmssmSusyRGE & nms);
    
    inline const NmssmSusy & displayNmssmSusy() const { return *this; };
    
    virtual ~NmssmSusy(); ///< Default destructor
    
    /// sets object to be equal to another
    const NmssmSusy & operator=(const NmssmSusy & s);
    
    void setNmssmLoops(double l) { nmssmSusyApprox.setLoops(l); };
    void setNmssmApprox(int l, int t) { nmssmSusyApprox.setLoops(l); 
      nmssmSusyApprox.setThresholds(t); }; 
    void setNmssmApprox(const Approx & a) { nmssmSusyApprox = a; };
    void setNmssmSusy(const NmssmSusy & y) { 
      setNmssmApprox(y.displayNmssmSusyApprox()); 
      setMssmSusy(y.displayMssmSusy());
      setNmssmSusyPars(y.displayNmssmSusyPars());
    };    

    int displayNmssmLoops() const { return nmssmSusyApprox.displayLoops(); };
    int displayNmssmThresholds() const { 
      return nmssmSusyApprox.displayThresholds(); 
    };
    Approx displayNmssmSusyApprox() const { return nmssmSusyApprox; };
    /// Returns all parameters as elements of a vector
    const DoubleVector display() const;

    void set(const DoubleVector & y, int & k);
    void set(const DoubleVector & y);
    /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
    NmssmSusy beta(nmsBrevity &) const;
    /// Outputs one-loop anomlous dimensions gii given matrix inputs.
    /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
    /// respectively. Note that we use the convention (for matrices in terms of
    /// gamma's): gamma^Li_Lj = M_ij for LH fields and
    /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
    /// conjugates of the RH fields). a should already be defined.
    void getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		      DoubleMatrix & gQQ, DoubleMatrix & gDD,
		      DoubleMatrix & gUU, double & gH1H1, double &
		      gH2H2, double & gSS, nmsBrevity & a) const;
    /// Outputs two-loop anomlous dimensions gii given matrix inputs.
    /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
    /// respectively. Note that we use the convention (for matrices in terms of
    /// gamma's): gamma^Li_Lj = M_ij for LH fields and
    /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
    /// conjugates of the RH fields). a should already be defined.
    void getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		      DoubleMatrix & gQQ, DoubleMatrix & gDD,
		      DoubleMatrix & gUU, double & gH1H1, double &
		      gH2H2, double & gSS, nmsBrevity & a) const;
    /// Outputs wave function renormalisation for SUSY parameters and gauge beta
    /// functions up to 2 loops. Also calculates and outputs a.
    /// IO parameters: RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1
    /// and H2 respectively.
    /// g^Li_Lj = m_{ij} for LH fields
    /// g^Ei_Ej = m_{ji} for RH fields
    /// Returns the effective mu parameter
    double displayMuEff() const { return displaySusyMu() + 
	displayLambda() * displaySvev() / sqrt(2.0); };
  }; ///< end of NmssmSusy
  
  /// Contains all supersymmetric RPC-MSSM parameters and is an RGE
  class NmssmSusyRGE: public RGE, public NmssmSusy {
  private:
  public:
    NmssmSusyRGE(); ///< Constructor fills object with zeroes by default
    /// Constructor sets object to be equal to another
    NmssmSusyRGE(const NmssmSusyRGE &);
    /// PA: Constructor given Yukawa matrices u,d,e, gauge couplings v, mu
    /// parameter=m, tan beta=tb, lambda, kappa, mupr, xiF,
    //renormalisation scale MU, number of loops in
    /// RG evolution l and thresholds parameter t
    NmssmSusyRGE(const DoubleMatrix & u, const DoubleMatrix & d, const
		 DoubleMatrix & e, const DoubleVector & v, double m,
		 double tb, double MU, int l, int t, double h, int mix,
		 double s, double lambda, double kappa, double xiF,
		 double mupr);
    
    NmssmSusyRGE(const MssmSusy &m);
    
    virtual ~NmssmSusyRGE(); ///< Default destructor
    
    /// sets object to be equal to another
    const NmssmSusyRGE & operator=(const NmssmSusyRGE & s);
    /// set the MssmSusy part
    const NmssmSusyRGE & operator=(const NmssmSusy & s);
    /// sets object to be equal to another
    void setSusy(const NmssmSusyRGE &s);
    void set(const DoubleVector &y);
    
    /// Copies MSSM Yukawa matrices and gauge couplings from s only
    void setSomePars(const NmssmSusyRGE & s);

    /// Returns whole object as a const
    inline const NmssmSusyRGE & displaySusy() const { return *this; };
    /// Returns all parameters as elements of a vector
    const DoubleVector display() const;
    
    /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
    DoubleVector beta() const;
  }; ///< end of NmssmSusyRGE
  
  
  /// Formatted output
  ostream & operator <<(ostream &, const NmssmSusyRGE &);
  ostream & operator <<(ostream &, const NmssmSusy &);
  ostream & operator <<(ostream &, const NmssmSusyPars &);
  
  /// Formatted input
  istream & operator >>(istream &left, NmssmSusyRGE &s);
  
} ///< namespace softsusy

#endif
