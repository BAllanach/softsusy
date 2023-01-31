
/** \file nmssmsusy.cpp
    - Project:     SOFTSUSY
    - Author:      Ben Allanach, Peter Athron, Alexander Voigt
    - Manual:      hep-ph/0104145
    - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
    - Description: All NMSSM SUSY respecting parameters
*/

#include "nmssmsusy.h"
#include <cassert>

namespace softsusy {

  /*** start of nmsBrevity **/  
  const nmsBrevity & nmsBrevity::operator=(const sBrevity &s) {
    if (this == &s) return *this;
    sBrevity::operator=(s);
    return *this;
  }
  
  const nmsBrevity & nmsBrevity::operator=(const nmsBrevity &s) {
    if (this == &s) return *this;
    sBrevity::operator=(s);
    lsq = s.lsq; ksq = s.ksq; l4 = s.l4; k4 = s.k4;
    return *this;
  }
  
  void nmsBrevity::calculate(const DoubleMatrix & yu, const DoubleMatrix & yd,
			     const DoubleMatrix & ye, const DoubleVector & g,
			     const double & lam, const double & kap) {
    sBrevity::calculate(yu, yd, ye, g);
    lsq = lam * lam; ksq = kap * kap;
    l4 = lsq * lsq; k4 = ksq * ksq;
  }
  /** end of nmsBrevity **/  

  /** start of RGE functions **/
  void anomalousDimensionNmssmSusy(const MssmSusy & s, const NmssmSusyPars & n,
				   DoubleMatrix & gEE, 
				   DoubleMatrix & gLL,DoubleMatrix & gQQ, 
				   DoubleMatrix & gUU, DoubleMatrix & gDD, 
				   DoubleVector & dg, double & gH1H1, 
				   double & gH2H2, double & gSS, 
				   nmsBrevity & a) {
    double kap = n.displayKappa(), lam = n.displayLambda();
    // Constants for gauge running
    static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3), clBeta(3);
    static DoubleMatrix babBeta(3, 3);
    if (bBeta(1) < 1.0e-5) // Constants not set yet
      nmsetBetas(babBeta, cuBeta, cdBeta, ceBeta, clBeta, bBeta);
    
    // nmsBrevity a contains all of the shortcutted matrices etc;
    a.calculate(s.displayYukawaMatrix(YU), s.displayYukawaMatrix(YD), 
		s.displayYukawaMatrix(YE), s.displayGauge(), lam, kap);
    
    // For calculational brevity
    double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT, &lsq = a.lsq;
    DoubleVector &gsq=a.gsq, &g3=a.g3;
    
    // 1 loop contributions:
    if (s.displayMssmLoops() > 0) {
      static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
      s.getOneLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);
      addOneLpAnomNmssm(gH1H1, gH2H2, gSS, lam, kap);
      dg = oneO16Pisq * g3 * bBeta;
    }
    
    if (s.displayMssmLoops() > 1) {
      s.getTwoLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);
      addTwoLpAnomNmssm(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, gSS, a);
      const static double twolp = 4.010149318236068e-5;
      dg = dg + g3 * (babBeta * gsq - cuBeta * uuT - cdBeta *
		  ddT - ceBeta * eeT - clBeta * lsq) * twolp;
      //PA: checked numerically 14/9/2012
    }
  }

  void addOneLpAnomNmssm(double & gH1H1, double & gH2H2, double & gSS, 
			 double lambda, double kappa) {
    const double ksq = sqr(kappa), lsq = sqr(lambda);
    static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
    
    gH1H1 += oneO16Pisq * lsq;
    gH2H2 += oneO16Pisq * lsq;
    gSS   += oneO16Pisq * (2.0 * lsq + 2.0 * ksq);
  }
    
  // adds two-loop anomalous dimension contribution to gii given matrix inputs
  // g^Li_Lj = m_{ij} for LH fields
  // g^Ei_Ej = m_{ji} for RH fields CHECKED: 23/5/02
  void addTwoLpAnomNmssm(DoubleMatrix & gEE, DoubleMatrix & gLL,
			 DoubleMatrix & gQQ, DoubleMatrix & gDD,
			 DoubleMatrix & gUU, double & gH1H1, double &
			 gH2H2, double & gSS, nmsBrevity & a) {
    // For calculational brevity
    DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2,
      &u2t=a.u2t, &d2t=a.d2t, &e2t=a.e2t;
    double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
    double &lsq = a.lsq, &ksq = a.ksq, &l4 = a.l4, &k4 = a.k4;
    DoubleVector &gsq = a.gsq;
    
    // Everything gets the (1/16pi^2)^2 factor at the bottom
    double h1h1, h2h2;
    
    // Two-loop pure gauge anom dimensions
    double ss = 1.2 * lsq * gsq(1) + 6.0 * lsq * gsq(2);
    
    // Two-loop pure Yukawa contributions
    const double s = (eeT + 3.0 * ddT);
    
    const DoubleMatrix ll(- lsq * e2);
    const DoubleMatrix ee(- 2.0 * lsq * e2t);
    const DoubleMatrix qq(- lsq * (u2 + d2));
    const DoubleMatrix dd(- 2.0 * lsq * d2t);
    const DoubleMatrix uu(- 2.0 * lsq * u2t);
    h1h1 = - (3.0 * l4 + 2.0 * lsq * ksq + 3.0 * lsq * uuT);
    h2h2 = - (3.0 * l4 + 2.0 * lsq * ksq + lsq * s);
    ss = ss - (4.0 * l4 + 8.0 * k4 + 8.0 * lsq * ksq + 2.0 * lsq * s
	       + 6.0 * lsq * uuT);
    
    const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2
    
    gLL = gLL + twolp * ll;
    gEE = gEE + twolp * ee;
    gQQ = gQQ + twolp * qq;
    gDD = gDD + twolp * dd;
    gUU = gUU + twolp * uu;
    gH1H1 = gH1H1 + twolp * h1h1;
    gH2H2 = gH2H2 + twolp * h2h2;
    gSS = gSS + twolp * ss;
  }

  // Outputs derivatives (DRbar scheme) in the form of ds. a contains the
  // matrices calculated that are handy for computation.
  // W=  LL Y^E H1 ER + QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1
  // is the superpotential. Consistent with Allanach, Dedes, Dreiner
  // hep-ph/9902251 and Barger, Berger and Ohmann hep-ph/9209232, 9311269
  // EXCEPT for the sign of smu, which is opposite. These equations are also
  // valid for W=  - LL Y^E H1 ER - QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1, the
  // New SOFTSUSY convention
  NmssmSusy betaNmssmSusy(nmsBrevity & a, const MssmSusy & s, 
			  const NmssmSusyPars & n) {
    // Wave function renormalisations: convention for g**(i, j) is that i is the
    // LOWER index and j the upper in our paper hep-ph/9902251
    static DoubleMatrix gEE(3, 3), gLL(3, 3), gQQ(3, 3), gDD(3, 3),
      gUU(3, 3);
    
    double gH1H1=0.0, gH2H2=0.0, gSS = 0.0;
    static DoubleVector dg(1,3);
    
    // keep this option in order to interface with extensions
    anomalousDimensionNmssmSusy(s, n, gEE, gLL, gQQ, gUU,
				gDD, dg, gH1H1, gH2H2, gSS, a); 
    
    // To keep this a const function
    const DoubleMatrix &u1 = s.displayYukawaMatrix(YU).display(),
      &d1 = s.displayYukawaMatrix(YD).display(),
      &e1 = s.displayYukawaMatrix(YE).display();
    
    // contain derivatives of up, down quarks and leptons
    static DoubleMatrix du(3, 3), dd(3, 3), de(3, 3);
    static double dl, dk, dz, dmupr;
    // mu parameter derivatives
    double dmu;
    
    // RGEs of SUSY parameters
    du = u1 * (gUU + gH2H2) + gQQ * u1;
    dd = d1 * (gDD + gH1H1) + gQQ * d1;
    de = e1 * (gEE + gH1H1) + gLL * e1;
    dl = n.displayLambda() * (gH1H1 + gH2H2 + gSS);
    dk = n.displayKappa() * (3.0 * gSS);

    dmu = s.displaySusyMu() * (gH1H1 + gH2H2);
    dmupr = 2.0 * n.displayMupr() * gSS;
    dz = n.displayXiF() * gSS;
    
    // Following is from hep-ph/9308335: scalar H anomalous dimensions (as
    // opposed to the chiral superfield one - see hep-ph/0111209).
    // Additional contribution from Feynman gauge running at two-loops of tan
    // beta: we need this to link up with BPMZ: hep-ph/0112251
    const double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
    const DoubleVector &gsq=a.gsq, &g4 = a.g4;
    const DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2, &d2t=a.d2t;
    const double t = (d2 * u2).trace();
    const double &lsq = a.lsq, &ksq = a.ksq, &l4 = a.l4, &k4 = a.k4;
    static const double oneLoop = 1.0 / (16.0 * sqr(PI));
    double sH1H1 = oneLoop * (3.0 * ddT + eeT + lsq);
    double sH2H2 = oneLoop * (3.0 * uuT + lsq);
    
    const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2
    if (s.displayMssmLoops() > 1) {
      const double g4terms = 1.035 * g4(1) + 0.45 * gsq(1) * gsq(2) + 
	5.875 * g4(2);
      sH1H1 = sH1H1 + twolp *
	(-(3.0 * (e2 * e2).trace() + 9.0 * (d2t * d2t).trace() + 3.0 * t) +
	 (16 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT
	 + g4terms
	 - 2 * ksq * lsq - 3 * l4 - 3 * lsq * uuT);
      sH2H2 = sH2H2 + twolp *
	(- (9.0 * (u2 * u2).trace() + 3.0 * t) +
	 (16 * gsq(3) + 0.8 * gsq(1)) * uuT
	 + g4terms
	 - 2 * ksq * lsq - 3 * l4 - lsq * (3 * ddT + eeT));
    }
    
    double cosb2 = sqr(cos(atan(s.displayTanb()))), sinb2 = 1.0 - cosb2;
    double feynman = 1.5 * gsq(2) + 0.3 * gsq(1);
    /// One-loop RGEs in Feynman gauge
    double dt = s.displayTanb() * (sH1H1 - sH2H2);
    double dHvev = s.displayHvev() *
      (cosb2 * (-sH1H1 + feynman * oneLoop) +
       sinb2 * (-sH2H2 + feynman * oneLoop));
    double dSvev = - n.displaySvev() * oneLoop * 2.0 * (lsq + ksq);
    
    if (s.displayMssmLoops() > 1) {
      /// Two-loop pieces
      dt = dt + s.displayTanb() * twolp * 
	(3.0 * ddT + eeT - 3.0 * uuT) * feynman;
      dHvev = dHvev - s.displayHvev() * twolp * (cosb2 * (3.0 * ddT + eeT) +
						 sinb2 * 3.0 * uuT + lsq) * 
	feynman	+ s.displayHvev() * twolp * 4.5 * g4(2);
      dSvev = dSvev + 
	n.displaySvev() * twolp * (8 * k4 + 8 * ksq * lsq + 2 * lsq * 
				 (-2 * feynman + 2 * lsq + 3 * ddT + eeT + 
				  3.0 * uuT));
    }
   
    // Contains all susy derivatives:
    NmssmSusy ds(du, dd, de, dg, dmu, dt, dHvev, s.displayMixing(),
		 s.displayMssmLoops(),
		 s.displayMssmApprox().displayThresholds(), 
		 dSvev, dl, dk, dz, dmupr);
    
    return ds;
  }
  
  void nmsetBetas(DoubleMatrix & babBeta, DoubleVector &cuBeta, DoubleVector
		  & cdBeta, DoubleVector & ceBeta, DoubleVector & clBeta,
		  DoubleVector & bBeta) {
    setBetas(babBeta, cuBeta, cdBeta, ceBeta, bBeta);
    setBetaLambda(clBeta);
  }
  
  void setBetaLambda(DoubleVector& clBeta) {
    clBeta(1) = 6.0 / 5.0;  clBeta(2) = 2.0; clBeta(3) = 0.0;
  }
  /** end of RGE functions **/
  
  /** start of NmssmSusyPars **/
  NmssmSusyPars::NmssmSusyPars()
    : lambda(0.), kappa(0.), sVev(0.), xiF(0.), mupr(0.) {}

  NmssmSusyPars::NmssmSusyPars(const NmssmSusyPars & s) {
    lambda = s.displayLambda();
    kappa  = s.displayKappa();
    sVev   = s.displaySvev();
    xiF    = s.displayXiF();
    mupr   = s.displayMupr();
  }

  NmssmSusyPars::NmssmSusyPars(double l, double k, double s, double x, double m)
    : lambda(l), kappa(k), sVev(s), xiF(x), mupr(m) {
  }

  const NmssmSusyPars & NmssmSusyPars::operator=(const NmssmSusyPars &s) {
    if (this == &s) return *this;
    lambda = s.displayLambda();
    kappa  = s.displayKappa();
    sVev   = s.displaySvev();
    xiF    = s.displayXiF();
    mupr   = s.displayMupr();
    return *this;
  }

  void NmssmSusyPars::setNmssmSusyPars(const NmssmSusyPars &s) {
    sVev   = s.displaySvev();
    lambda = s.displayLambda();
    kappa  = s.displayKappa();
    mupr   = s.displayMupr();
    xiF    = s.displayXiF();
  }

  void NmssmSusyPars::setNmssmSusyPars(const DoubleVector & y) {
    setSvev(y.display(34));
    setLambda(y.display(35));
    setKappa(y.display(36));
    setMupr(y.display(37));
    setXiF(y.display(38));
  }
  /** end of NmssmSusyPars **/

  /** start of NmssmSusy **/
  NmssmSusy::NmssmSusy(double s, double l, double k, double z, 
			       double m) 
    : NmssmSusyPars(l, k, s, z, m) {}
  
  const NmssmSusy & NmssmSusy::operator=(const NmssmSusy & s) {
    if (this == &s) return *this;
    setNmssmSusyPars(s.displayNmssmSusyPars());
    setSvev(s.displaySvev());
    setLambda(s.displayLambda());
    setKappa(s.displayKappa());
    setMupr(s.displayMupr());
    setXiF(s.displayXiF());
    nmssmSusyApprox = s.displayNmssmSusyApprox();
    return *this;
  }

  NmssmSusy::NmssmSusy()
    :  MssmSusy(), NmssmSusyPars(), nmssmSusyApprox() {}
  
  NmssmSusy::NmssmSusy(const MssmSusy & m, const NmssmSusyPars & nsp) 
    : MssmSusy(m), NmssmSusyPars(nsp), nmssmSusyApprox(m.displayMssmApprox()) {}

  NmssmSusy::NmssmSusy(const NmssmSusy & n)
    : MssmSusy(n.displayMssmSusy()), NmssmSusyPars(n.displayNmssmSusyPars()),
      nmssmSusyApprox(n.displayNmssmSusyApprox()) {}

  NmssmSusy::NmssmSusy(const DoubleMatrix & u,
		       const DoubleMatrix & d, const
		       DoubleMatrix & e, const DoubleVector & v, 
		       double m,
		       double tb,  double hv, int mix, int l, int t, double sv,
		       double lam, double kap, double z, 
		       double mup)
    : MssmSusy(u, d, e, v, m, tb, hv, mix),
      NmssmSusyPars(lam, kap, sv, z, mup) {
    setNmssmApprox(l, t);
  }

  //  NmssmSusy::NmssmSusy(const NmssmSusyRGE & nms)
  //    : MssmSusy(nms.displayMssmSusy()) {}

  // outputs one-loop anomlous dimensions gii given matrix inputs
  // Note that we use the convention (for matrices in terms of gamma's)
  // gamma^Li_Lj = M_ij for LH fields and
  // gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
  // conjugates of the RH fields): CHECKED 23/5/02
  void NmssmSusy::getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
			       DoubleMatrix & gQQ, DoubleMatrix & gDD,
			       DoubleMatrix & gUU, double & gH1H1, double &
			       gH2H2, double & gSS, nmsBrevity & a) const {
    MssmSusy::getOneLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);
    addOneLpAnomNmssm(gH1H1, gH2H2, gSS, displayLambda(), displayKappa());
  }
  
  
  // adds two-loop anomalous dimension contribution to gii given matrix inputs
  // g^Li_Lj = m_{ij} for LH fields
  // g^Ei_Ej = m_{ji} for RH fields CHECKED: 23/5/02
  void NmssmSusy::getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
			       DoubleMatrix & gQQ, DoubleMatrix & gDD,
			       DoubleMatrix & gUU, double & gH1H1, double &
			       gH2H2, double & gSS, nmsBrevity & a) const {
    MssmSusy::getTwoLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);
    addTwoLpAnomNmssm(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, gSS, a);
  }

  NmssmSusy NmssmSusy::beta(nmsBrevity & a) const {
    return betaNmssmSusy(a, displayMssmSusy(), displayNmssmSusyPars());
  }
  /** end of NmssmSusy **/

  /** start of NmssmSusyRGE **/
  NmssmSusyRGE::NmssmSusyRGE()
    : NmssmSusy() {
    setPars(numNMssmPars);
  }
  
  NmssmSusyRGE::NmssmSusyRGE(const NmssmSusyRGE &s)
    : NmssmSusy(s.displayNmssmSusy()) {
    setSusy(s.displayMssmSusy());
    setPars(numNMssmPars);
  }
  
  NmssmSusyRGE::NmssmSusyRGE(const MssmSusy &m)
    : NmssmSusy() {
    setSusy(m);
    setPars(numNMssmPars);
  }
  
  
  NmssmSusyRGE::NmssmSusyRGE(const DoubleMatrix & u,
			     const DoubleMatrix & d, const
			     DoubleMatrix & e, const DoubleVector & v, double m,
			     double tb, double MU, int l, int t, double hv, 
			     int mix, double sv,
			     double lam, double kap, double mpr, double z)
    : NmssmSusy(u, d, e, v, m, tb, hv, mix, l, t, sv, lam, kap, z, mpr) {
    setPars(numNMssmPars);
    setMu(MU);
    setNmssmApprox(l, t);
  }
  
  
  NmssmSusyRGE::~NmssmSusyRGE() {}
  NmssmSusy::~NmssmSusy() {}
  
  const NmssmSusyRGE & NmssmSusyRGE::operator=(const NmssmSusyRGE & s) {
    if (this == &s) return *this;
    MssmSusy::operator=(s.displayMssmSusy());
    NmssmSusy::operator=(s.displayNmssmSusy());
    
    return *this;
  }
  
  
  const NmssmSusyRGE & NmssmSusyRGE::operator=(const NmssmSusy & s) {
    setMssmSusy(s);
    NmssmSusy::operator=(s);

    return *this;
  }
  
  void NmssmSusyRGE::setSomePars(const NmssmSusyRGE & s) {
    MssmSusy::setSomePars(s);
  }
    
  void NmssmSusyPars::display(DoubleVector & y, int & k) const {
    y(k) = displaySvev(); k++;
    y(k) = displayLambda(); k++;
    y(k) = displayKappa(); k++;
    y(k) = displayMupr(); k++;
    y(k) = displayXiF(); k++;
  }

  void NmssmSusyPars::set(const DoubleVector & y, int & k) {
    setSvev(y.display(k)); k++;
    setLambda(y.display(k)); k++;
    setKappa(y.display(k)); k++;
    setMupr(y.display(k)); k++;
    setXiF(y.display(k)); k++;
  }

  void NmssmSusyRGE::set(const DoubleVector &y) { 
    int k = 1; NmssmSusy::set(y, k);
  }

  const DoubleVector NmssmSusy::display() const {
    DoubleVector y(MssmSusy::display());
    assert(y.displayStart() == 1 && y.displayEnd() == numSusyPars);
    y.setEnd(numNMssmPars);
    int k = 34; 
    NmssmSusyPars::display(y, k);
    
    return y;
  }

  const DoubleVector NmssmSusyRGE::display() const {
    return NmssmSusy::display();
  }
  
  void NmssmSusy::set(const DoubleVector & y, int & k) {
    assert(y.displayEnd() - y.displayStart() + 1 >= numNMssmPars);
    MssmSusy::set(y, k); k++;
    NmssmSusyPars::set(y, k);
  }
  
  void NmssmSusy::set(const DoubleVector & y) { int k =1;
    set(y, k); }

  ostream & operator <<(ostream &left, const NmssmSusyPars &s) {
    left << "singlet VEV: " << s.displaySvev()
	 << " lambda: " << s.displayLambda()
	 << " kappa: " << s.displayKappa() << endl
	 << " mupr: " << s.displayMupr()
	 << " xiF: " << s.displayXiF()
	 << '\n';
    return left;
  }

  ostream & operator <<(ostream &left, const NmssmSusy &s) {
    left << s.displayMssmSusy()
	 << s.displayNmssmSusyPars();
    return left;
  }

  ostream & operator <<(ostream &left, const NmssmSusyRGE &s) {
    left << "NMSSM SUSY parameters at Q: " << s.displayMu() << endl
	 << s.displayNmssmSusy();
    return left;
  }
  
  void NmssmSusyRGE::setSusy(const NmssmSusyRGE & s) {
    MssmSusy::setMssmSusy(s);
    setLambda(s.displayLambda());
    setKappa(s.displayKappa());
    setSvev(s.displaySvev());
    setMupr(s.displayMupr());
    setXiF(s.displayXiF());
  }
  
  istream & operator >>(istream &left, NmssmSusyRGE &s) {
    MssmSusy ms;
    left >> ms;
    s = ms;
    double mupr = 0.0, xiF = 0.0, sv = 0.0, lambda = 0.0, kappa = 0.0;
    s.setLambda(lambda);
    s.setKappa(kappa);
    s.setSvev(sv);
    s.setMupr(mupr);
    s.setXiF(xiF);
    return left;
  }  
  
  // Outputs derivatives vector y[n] for SUSY parameters: interfaces to
  // integration routines
  DoubleVector NmssmSusyRGE::beta() const {
    static nmsBrevity a;
    
    // calculate the derivatives
    static NmssmSusyRGE ds;
    
    ds = NmssmSusy::beta(a);

    return ds.display(); // convert to a long vector
  }
  /** end of NmssmSusyRGE **/
  
} // namespace softsusy
