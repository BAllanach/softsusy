
/** \file nmssmsusy.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Peter Athron, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief nMssmSusy contains all NMSSM SUSY couplings and tan beta, as
   well as their beta functions
*/

#ifndef NMSSMSUSY_H
#define NMSSMSUSY_H

#include <susy.h>
using namespace softsusy;

// with the 5 new nmssm parameters we now have 38 in total
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


/// Contains all supersymmetric RPC-MSSM parameters and RGEs
class nMssmSusy: public MssmSusy
{
private:
  /// new nmssm parameters, lambda, kappa appearing as superpotential
  /// terms, lambda S H_u H_d and \frac{1}{3} kappa S^3 and sVev is
  /// the vev of the singlet superfield defined by sVev = \sqrt{2}<S>.
  double lambda, kappa, sVev;
  /// new parameters appearing the general MSSM, but not the ususal
  /// form of the NMSSM which solves the mu-problem nor the MSSM.
  /// They appear in superpotential terms zeta S and \frac{1}{2} mu_s
  /// S^2
  double zeta, mu_s;

public:
  nMssmSusy(); ///< Constructor fills object with zeroes by default
  /// Constructor sets object to be equal to another
  nMssmSusy(const nMssmSusy &);
  /// PA: Constructor given Yukawa matrices u,d,e, gauge couplings v, mu
  /// parameter=m, tan beta=tb, lambda, kappa, mu_s, zeta,
  //renormalisation scale MU, number of loops in
  /// RG evolution l and thresholds parameter t
  nMssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
            DoubleMatrix & e, const DoubleVector & v, double m,
            double tb, double MU, int l, int t, double h, double s,
            double lambda, double kappa, double zeta, double mu_s);
  virtual ~nMssmSusy(); ///< Default destructor

  /// sets object to be equal to another
  const nMssmSusy & operator=(const nMssmSusy & s);
  /// set the MssmSusy part
  const nMssmSusy & operator=(const MssmSusy & s);
  /// sets object to be equal to another
  void setSusy(const nMssmSusy &s);

  /// sets DRbar running singlet vev.
  void setSvev(double s);
  /// Copies MSSM Yukawa matrices and gauge couplings from s only
  void setSomePars(const nMssmSusy & s);
  /// sets the \lambda S H_u H_d coupling
  void setLambda(double);
  /// sets the \frac{1}{3} \kappa S^3 coupling
  void setKappa(double);
  /// sets the \frac{1}{2} mu_s S^2 coupling
  void setMu_s(double);
   /// sets the zeta S coupling
  void setZeta(double);
  /// Sets all RGE parameters to elements of vector
  void set(const DoubleVector &);

  /// returns DRbar running Singlet Higgs vev
  double displaySvev() const;
  /// Returns whole object as a const
  inline const nMssmSusy & displaySusy() const;
  /// returns superpotential parameter lambda
  double displayLambda() const;
  /// returns superpotential parameter lambda
  double displayKappa() const;
  /// returns mu_s superpotential parameter
  double displayMu_s() const;
  /// returns zeta superpotential parameter
  double displayZeta() const;
  /// Returns all parameters as elements of a vector
  const DoubleVector display() const;

  /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
  DoubleVector beta() const;
  /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
  nMssmSusy beta(nmsBrevity &) const;
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
  void anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
			  DoubleMatrix & gQQ, DoubleMatrix & gUU,
			  DoubleMatrix & gDD, DoubleVector & dg,
			  double & gH1H1, double & gH2H2, double & gSS,
			  nmsBrevity & a) const;
};

/// Formatted output
ostream & operator <<(ostream &, const nMssmSusy &);

/// Formatted input
istream & operator >>(istream &left, nMssmSusy &s);

/// Outputs beta function coefficients for MSSM gauge coupling evolution in
/// arguments.
void nmsetBetas(DoubleMatrix &, DoubleVector  &, DoubleVector  &, DoubleVector
                &, DoubleVector  &, DoubleVector  & );

inline const nMssmSusy & nMssmSusy::displaySusy() const { return *this; }

inline double nMssmSusy::displaySvev() const { return sVev; }

inline void nMssmSusy::setSvev(double h) { sVev = h; }
inline void nMssmSusy::setMu_s(double f) { mu_s = f; }
inline void nMssmSusy::setZeta(double z) { zeta = z; }

inline void nMssmSusy::setLambda(double l) { lambda = l; }
inline void nMssmSusy::setKappa(double k) { kappa = k; }

inline double nMssmSusy::displayMu_s() const { return mu_s; }
inline double nMssmSusy::displayZeta() const { return zeta; }
inline double nMssmSusy::displayLambda() const { return lambda; }
inline double nMssmSusy::displayKappa() const { return kappa; }

#endif
