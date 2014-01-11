/** \file susy.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: All SUSY respecting parameters

*/

#include "susy.h"

namespace softsusy {

#define HR "---------------------------------------------------------------\n"

const sBrevity & sBrevity::operator=(const sBrevity &s) {
  if (this == &s) return *this;
  dt = s.dt; ut = s.ut; et = s.et; 
  u2 = s.u2; d2 = s.d2; e2 = s.e2; 
  u2t = s.u2t; e2t = s.e2t; d2t = s.d2t; 
  gsq = s.gsq; g3 = s.g3; g4 = s.g4;
  uuT = s.uuT; ddT = s.ddT; eeT = s.eeT;
  d1 = s.d1; e1 = s.e1; u1 = s.u1;
  return *this;
}

void sBrevity::calculate(const DoubleMatrix & yu, const DoubleMatrix & yd,
			  const DoubleMatrix & ye, const DoubleVector & g) {
  static DoubleVector g1(1, 3);
  g1 = g.display();
  u1 = yu.display(); 
  d1 = yd.display(); 
  e1 = ye.display();
  dt = d1.transpose(); ut = u1.transpose(); et = e1.transpose();
  u2 = u1 * ut; d2 = d1 * dt; e2 = e1 * et; 
  u2t = ut * u1; d2t = dt * d1; e2t = et * e1;
  uuT = u2.trace(), ddT = d2.trace(), eeT = e2.trace();
  gsq = g1 * g1; g3 = gsq * g1; g4 = g3 * g1;
}


MssmSusy::MssmSusy()
  : u(3, 3), d(3, 3), e(3, 3), g(3), smu(0.0), tanb(0.0), hVev(0.0) {
    setPars(numSusyPars);
    setMu(0.0);
    setLoops(2);
    setThresholds(0);
}

MssmSusy::MssmSusy(const MssmSusy &s)
  : u(s.u), d(s.d), e(s.e), g(s.g), smu(s.smu), tanb(s.tanb), hVev(s.hVev) { 
    setPars(numSusyPars);
    setMu(s.displayMu()); 
    setLoops(s.displayLoops());
    setThresholds(s.displayThresholds());
}

MssmSusy::MssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
		     DoubleMatrix & e, const DoubleVector & v, double m,
		     double tb, double MU, int l, int t, double hv)
  : u(u), d(d), e(e), g(v), smu(m), tanb(tb), hVev(hv) { 
    setPars(numSusyPars);
    setMu(MU); 
    setLoops(l);
    setThresholds(t);
}

const MssmSusy & MssmSusy::operator=(const MssmSusy & s) {
  if (this == &s) return *this;
  u = s.u;
  d = s.d;
  e = s.e;
  smu = s.smu;
  tanb = s.tanb;
  g = s.g;
  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
  hVev = s.hVev;
  return *this;
}

void MssmSusy::setSomePars(const MssmSusy & s) {
  u = s.u;
  d = s.d;
  e = s.e;
  g = s.g;
}

void MssmSusy::setYukawaElement(yukawa k, int i, int j, double f) { 
  switch(k) {
  case YU: u(i, j) = f; break;
  case YD: d(i, j) = f; break;
  case YE: e(i, j) = f; break;
  default: 
    ostringstream ii;
    ii << "MssmSusy::set called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
}

void MssmSusy::setYukawaMatrix(yukawa k, const DoubleMatrix & m) { 
  switch(k) {
  case YU: u = m; break;
  case YD: d = m; break;
  case YE: e = m; break;
  default: 
    ostringstream ii;
    ii << "MssmSusy::set called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
}

double MssmSusy::displayYukawaElement(yukawa k, int i, int j) const {
  switch(k) {
  case YU: return u.display(i, j); break;
  case YD: return d.display(i, j); break;
  case YE: return e.display(i, j); break;
  default: 
    ostringstream ii;
    ii << "MssmSusy::display called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
  return 0.0;
}

const DoubleMatrix & MssmSusy::displayYukawaMatrix(yukawa k) const {
  switch(k) {
  case YU: return u; break;
  case YD: return d; break;
  case YE: return e; break;
  default: 
    ostringstream ii;    
    ii << "MssmSusy::display called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
}

const DoubleVector MssmSusy::display() const {
  DoubleVector y(numSusyPars);
  int i, j, k=0;
  for (i=1; i<=3; i++)    
    for (j=1; j<=3; j++) {
      k++;
      y(k) = u.display(i, j);
      y(k+9) = d.display(i, j);
      y(k+18) = e.display(i, j);
    }
  k=27;
  for (i=1; i<=3; i++) {
    k++;
    y(k) = g.display(i);
  }
  y(31) = smu;
  y(32) = tanb;
  y(33) = hVev;
  return y;
}

void MssmSusy::set(const DoubleVector & y) {
  int i, j, k=0;
  for (i=1; i<=3; i++)    
    for (j=1; j<=3; j++){
      k++;
      u(i, j) = y.display(k);
      d(i, j) = y.display(k+9);
      e(i, j) = y.display(k+18);
    }
  k=27;
  for (i=1; i<=3; i++) {
    k++;
    g(i) = y.display(k);
  }
  smu = y.display(31);
  tanb = y.display(32);
  hVev = y.display(33);
}

double MssmSusy::displayTanb() const { return tanb; }

ostream & operator <<(ostream &left, const MssmSusy &s) {
  left << "Supersymmetric parameters at Q: " << s.displayMu() << endl;
  left << " Y^U" << s.displayYukawaMatrix(YU) << " Y^D" <<
    s.displayYukawaMatrix(YD) << " Y^E" << s.displayYukawaMatrix(YE);
  left << "higgs VEV: " << s.displayHvev() 
       << " tan beta: " << s.displayTanb() << " smu: " << s.displaySusyMu() << 
    "\n";
  left << "g1: " << s.displayGaugeCoupling(1) << " g2: " <<
    s.displayGaugeCoupling(2) << " g3: " << 
    s.displayGaugeCoupling(3) << endl; 
  left << "thresholds: " << s.displayThresholds() 
       << " #loops: " << s.displayLoops() << '\n';
  return left;
}

void MssmSusy::setSusy(const MssmSusy & s) {
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
  setMu(s.displayMu());
  setYukawaMatrix(YU, s.displayYukawaMatrix(YU)); 
  setYukawaMatrix(YD, s.displayYukawaMatrix(YD)); 
  setYukawaMatrix(YE, s.displayYukawaMatrix(YE)); 
  setHvev(s.displayHvev());
  setTanb(s.displayTanb());
  setSusyMu(s.displaySusyMu());
  setAllGauge(s.displayGauge());
}

istream & operator >>(istream &left, MssmSusy &s) {
  string c;
  DoubleMatrix u(3, 3), d(3, 3), e(3, 3);
  double g1, g2, g3, smu, mu, tanb, hv;
  int loops, thresh;
  left >> c >> c >> c >> c >> mu;
  left >> c >> u >> c >> d >> c >> e >> c >> c >> hv;
  left >> c >> c >> tanb >> c >> smu;
  left >> c >> g1 >> c >> g2 >> c >> g3;
  left >> c >> thresh >> c >> loops;
  s.setYukawaMatrix(YU, u);
  s.setYukawaMatrix(YD, d);
  s.setYukawaMatrix(YE, e);
  s.setHvev(hv);
  s.setTanb(tanb);
  s.setGaugeCoupling(1, g1);
  s.setGaugeCoupling(2, g2);
  s.setGaugeCoupling(3, g3);
  s.setThresholds(thresh);
  s.setSusyMu(smu);
  s.setMu(mu);
  s.setLoops(loops);
  return left;
}



// Outputs derivatives (DRbar scheme) in the form of ds. a contains the
// matrices calculated that are handy for computation.
// W=  LL Y^E H1 ER + QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1
// is the superpotential. Consistent with Allanach, Dedes, Dreiner
// hep-ph/9902251 and Barger, Berger and Ohmann hep-ph/9209232, 9311269
// EXCEPT for the sign of smu, which is opposite. These equations are also
// valid for W=  - LL Y^E H1 ER - QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1, the
// new SOFTSUSY convention
MssmSusy MssmSusy::beta(sBrevity & a) const {
  // Wave function renormalisations: convention for g**(i, j) is that i is the
  // LOWER index and j the upper in our paper hep-ph/9902251
  static DoubleMatrix gEE(3, 3), gLL(3, 3), gQQ(3, 3), gDD(3, 3), 
    gUU(3, 3);
  
  double gH1H1=0.0, gH2H2=0.0;
  static DoubleVector dg(1,3);
  
  // keep this option in order to interface with RPVSUSY  
  anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg, gH1H1, gH2H2, a);
  
  // To keep this a const function
  const DoubleMatrix &u1 = u.display(), &d1 = d.display(), &e1 = e.display();
  
  // contain derivatives of up, down quarks and leptons
  static DoubleMatrix du(3, 3), dd(3, 3), de(3, 3); 
  // mu parameter derivatives
  double dmu;
  
  // RGEs of SUSY parameters
  du = u1 * (gUU + gH2H2) + gQQ * u1;
  dd = d1 * (gDD + gH1H1) + gQQ * d1;
  de = e1 * (gEE + gH1H1) + gLL * e1;
  
  dmu = smu * (gH1H1 + gH2H2);

  // Following is from hep-ph/9308335: scalar H anomalous dimensions (as
  // opposed to the chiral superfield one - see hep-ph/0111209).
  // Additional contribution from Feynman gauge running at two-loops of tan
  // beta: we need this to link up with BPMZ: hep-ph/0112251
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq, &g4 = a.g4;
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2, &d2t=a.d2t;
  double t = (d2 * u2).trace();
  static const double oneLoop = 1.0 / (16.0 * sqr(PI));
  double sH1H1 = oneLoop * (3.0 * ddT + eeT);
  double sH2H2 = oneLoop * 3.0 * uuT;

  const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2
  if (displayLoops() > 1) {
  const double g4terms = 1.035 * g4(1) + 0.45 * gsq(1) * gsq(2) + 5.875 * g4(2);
    sH1H1 = sH1H1 + twolp * 
      (-(3.0 * (e2 * e2).trace() + 9.0 * (d2t * d2t).trace() + 3.0 * t) + 
       (16 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT  + g4terms);
    sH2H2 = sH2H2 + twolp *
      (- (9.0 * (u2 * u2).trace() + 3.0 * t) +
       (16 * gsq(3) + 0.8 * gsq(1)) * uuT+ g4terms);
  }
  
  double cosb2 = sqr(cos(atan(tanb))), sinb2 = 1.0 - cosb2;
  double feynman = 1.5 * gsq(2) + 0.3 * gsq(1);
  /// One-loop RGEs in Feynman gauge
  double dt = displayTanb() * (sH1H1 - sH2H2);
  double dHvev = hVev * 
    (cosb2 * (-sH1H1 + feynman * oneLoop) + 
     sinb2 * (-sH2H2 + feynman * oneLoop)); 
  if (displayLoops() > 1) {
    /// Two-loop pieces
    dt = dt + displayTanb() * twolp * (3.0 * ddT + eeT - 3.0 * uuT) * feynman;
    dHvev = dHvev - hVev * twolp * (cosb2 * (3.0 * ddT + eeT) +
				    sinb2 * 3.0 * uuT) * feynman
                  + hVev * twolp * 4.5 * g4(2);
  }
  // Contains all susy derivatives:
  MssmSusy ds(du, dd, de, dg, dmu, dt, displayMu(), displayLoops(),
	       displayThresholds(), dHvev); 

  return ds;
}

void setBetas(DoubleMatrix & babBeta, DoubleVector &cuBeta, DoubleVector
	       & cdBeta, DoubleVector & ceBeta, DoubleVector & bBeta) {
  // 1 loop gauge beta fns
  bBeta(1) = 33.0 / 5.0; bBeta(2) = 1.0; bBeta(3) = -3.0; 
  
  // Extra sleptons included in vectorlike rep.s: 3 ER + 2 LL
  //#ifdef SLEPTONS
  //bBeta(1) = bBeta(1) + (3.0 * 1.2 + 0.6 * 2.0);
  //bBeta(2) = bBeta(2) + 2.0;
  //#endif
  
  // Next come the two loop MSSM constants for gauge beta fns
  babBeta(1, 1) = 199.0 / 25.0; babBeta(1, 2) = 27.0 / 5.0; 
  babBeta(1, 3) = 88.0 / 5.0; 
  babBeta(2, 1) = 9.0 / 5.0;    babBeta(2, 2) = 25.0;       
  babBeta(2, 3) = 24.0;
  babBeta(3, 1) = 11.0 / 5.0;   babBeta(3, 2) = 9.0;        
  babBeta(3, 3) = 14.0;
  cuBeta(1) = 26.0 / 5.0; cuBeta(2) = 6.0; cuBeta(3) = 4.0;
  cdBeta(1) = 14.0 / 5.0; cdBeta(2) = 6.0; cdBeta(3) = 4.0;
  ceBeta(1) = 18.0 / 5.0; ceBeta(2) = 2.0; ceBeta(3) = 0.0;
}

// outputs one-loop anomlous dimensions gii given matrix inputs
// Note that we use the convention (for matrices in terms of gamma's)
// gamma^Li_Lj = M_ij for LH fields and
// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
// conjugates of the RH fields): CHECKED 23/5/02
void MssmSusy::getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
				DoubleMatrix & gQQ, DoubleMatrix & gDD,
				DoubleMatrix & gUU, double & gH1H1, double &
				gH2H2, sBrevity & a) const {
  // For calculational brevity
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2, &u2t=a.u2t, &d2t=a.d2t,
    &e2t=a.e2t;      
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq;
  
  static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
  
  gEE = oneO16Pisq * (2.0 * e2t - 1.2 * gsq(1));
  gLL = oneO16Pisq * (e2 - (0.3 * gsq(1) + 1.5 * gsq(2)));
  gQQ = oneO16Pisq * (d2 + u2 - (gsq(1) / 30.0 + 1.5 * gsq(2) + 8 *
				    gsq(3) / 3.0));
  gUU = oneO16Pisq * (2.0 * u2t - (8 * gsq(1) / 15.0 + 8 * gsq(3) /
				      3.0)); 
  gDD = oneO16Pisq * (2.0 * d2t - 
			 (2 * gsq(1) / 15.0 + 8 * gsq(3) / 3.0));
  gH1H1 = oneO16Pisq * (3.0 * ddT + eeT - (0.3 * gsq(1) + 1.5 *
					      gsq(2)));
  gH2H2 = oneO16Pisq * (3.0 * uuT - (0.3 * gsq(1) + 1.5 * gsq(2)));
}

// adds two-loop anomalous dimension contribution to gii given matrix inputs
// g^Li_Lj = m_{ij} for LH fields
// g^Ei_Ej = m_{ji} for RH fields CHECKED: 23/5/02
void MssmSusy::getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
				DoubleMatrix & gQQ, DoubleMatrix & gDD,
				DoubleMatrix & gUU, double & gH1H1, double &
				gH2H2, sBrevity & a) const {
  // For calculational brevity
  DoubleMatrix &dt=a.dt, &ut=a.ut, &u2=a.u2, &d2=a.d2, &e2=a.e2,
    &u2t=a.u2t, &d2t=a.d2t, &e2t=a.e2t, &u1=a.u1, &d1=a.d1;      
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq, &g4=a.g4;
  
  // Everything gets the (1/16pi^2)^2 factor at the bottom
  DoubleMatrix ee(3, 3), ll(3, 3), qq(3, 3), dd(3, 3), uu(3, 3); 
  
  // Two-loop pure gauge anom dimensions
  double h1h1 = (3.75 * g4(2) + 2.07 * g4(1) + 0.9 * gsq(2) * gsq(1));
  double h2h2 = h1h1;
  ll = h1h1;
  ee = (234. * g4(1) / 25.0);
  qq = (-8.0 * g4(3) / 9.0 + 3.75 * g4(2) + 199.0 * g4(1) / 900.0 + 8.0 *
	gsq(3) * gsq(2) + 8 * gsq(3) * gsq(1) / 45.0 + 0.1 * gsq(1) *
	gsq(2));
  dd = (-8.0 * g4(3) / 9.0 + 202.0 / 225.0 * g4(1) + 32.0 / 45.0 *
	gsq(3) * gsq(1));
  uu = (-8.0 * g4(3) / 9.0 + 856.0 / 225.0 * g4(1) + 128.0 / 45.0 *
	gsq(3) * gsq(1));

  ll = ll + 1.2 * gsq(1) * e2;
  ee = ee + (6.0 * gsq(2) - 1.2 * gsq(1)) * e2t;
  qq = qq + 0.4 * gsq(1) * (d2 + 2.0 * u2);
  dd = dd + (6.0 * gsq(2) + 0.4 * gsq(1)) * d2t;
  uu = uu + (6.0 * gsq(2) - 0.4 * gsq(1)) * u2t;

  h1h1 = h1h1 + (16 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) *
    eeT; 
  h2h2 = h2h2 + (16 * gsq(3) + 0.8 * gsq(1)) * uuT;

  // Two-loop pure Yukawa contributions
  double s = (eeT + 3.0 * ddT), t = (d2 * u2).trace();

  ll = ll - (2.0 * e2 * e2 + s * e2);
  ee = ee - (2.0 * e2t * e2t + 2.0 * s * e2t);
  qq = qq - (2.0 * d2 * d2 + d2 * s + 2.0 * u2 * u2 + 3.0 * uuT * u2);
  dd = dd - (2.0 * d2t * d2t + 2.0 * (dt * u2 * d1) + 2 * s * d2t);
  uu = uu - (2.0 * u2t * u2t + 2.0 * (ut * d2 * u1) + 6.0 * uuT * u2t);
  h1h1 = h1h1 - (3.0 * (e2 * e2).trace() + 9.0 * (d2t * d2t).trace() +
		 3.0 * t);
  h2h2 = h2h2 - (9.0 * (u2 * u2).trace() + 3.0 * t);

  const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2
  
  gLL = gLL + twolp * ll;
  gEE = gEE + twolp * ee;
  gQQ = gQQ + twolp * qq;
  gDD = gDD + twolp * dd;
  gUU = gUU + twolp * uu;
  gH1H1 = gH1H1 + twolp * h1h1;
  gH2H2 = gH2H2 + twolp * h2h2;
}

// Outputs wave function renormalisation for SUSY parameters and gauge beta
// functions up to 2 loops. 
void MssmSusy::anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
				    DoubleMatrix & gQQ, DoubleMatrix & gUU,
				    DoubleMatrix & gDD, DoubleVector & dg, 
				    double & gH1H1, double & gH2H2, 
				    sBrevity & a)  const {
  // Constants for gauge running
  static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3);
  static DoubleMatrix babBeta(3, 3);
  if (bBeta(1) < 1.0e-5) // Constants not set yet
    setBetas(babBeta, cuBeta, cdBeta, ceBeta, bBeta);
  
  //  sBrevity a contains all of the shortcutted matrices etc;
  a.calculate(u.display(), d.display(), e.display(), g.display());
  
  // For calculational brevity
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq, &g3=a.g3;
  
  // 1 loop contributions: 
  if (displayLoops() > 0) {
    static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI)); 
    
    getOneLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);
    dg = oneO16Pisq * g3 * bBeta;  
  } 
  
  if (displayLoops() > 1) { 
    getTwoLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a); 
							      
    const static double twolp = 4.010149318236068e-5; 
    
    dg = dg + g3 * (babBeta * gsq - cuBeta * uuT - cdBeta *
		    ddT - ceBeta * eeT) * twolp;  
  }
  if (displayLoops() > 2) {
    DoubleVector & g4 = a.g4;
    DoubleMatrix & u2t = a.u2t, &d2t = a.d2t;
    const static double threelp = 2.53945672191370e-7;
    // from http://www.liv.ac.uk/~dij/betas/allgennb.log
    // only strong coupling constant + third family yukawa corrections are added
    getThreeLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a); 
    
    dg(3) = dg(3) +  g3(3) * ( 347.0/3.0 * g4(3) 
			       - 104.0/3.0 * gsq(3) * (uuT + ddT)
			       + 18.0 * (uuT * uuT + ddT * ddT) 
			       + 6.0 * (ddT * eeT)
			       + 12.0 * ( (u2t*u2t).trace() + 
					  (d2t*d2t).trace() ) 
			       + 8.0 * (u2t*d2t).trace())*threelp; 
  }
}

// Outputs derivatives vector y[n] for SUSY parameters: interfaces to
// integration routines
DoubleVector MssmSusy::beta() const {
  static sBrevity a;
  
  // calculate the derivatives
  static MssmSusy ds;
  
  ds = beta(a);

  return ds.display(); // convert to a long vector
}



// r should be valid AT mt
void MssmSusy::setDiagYukawas(const QedQcd & r, double vev) {

  double v1, v2; // Higgs VEVs

  v1 = vev * cos(atan(displayTanb()));
  v2 = vev * sin(atan(displayTanb()));

  DoubleMatrix u1(3, 3), d1(3, 3), e1(3, 3);
  
  double invv2, invv1; 
  invv2 = 1.0 / v2; invv1 = 1.0 / v1;
  u1(1, 1) = r.displayMass(mUp) * invv2;
  u1(2, 2) = r.displayMass(mCharm) * invv2;
  u1(3, 3) = r.displayMass(mTop) * invv2;
  
  d1(1, 1) = r.displayMass(mDown) * invv1;
  d1(2, 2) = r.displayMass(mStrange) * invv1;
  d1(3, 3) = r.displayMass(mBottom) * invv1;
  e1(1, 1) = r.displayMass(mElectron) * invv1;
  e1(2, 2) = r.displayMass(mMuon) * invv1;
  e1(3, 3) = r.displayMass(mTau) * invv1;
  
  setYukawaMatrix(YU, u1); 
  setYukawaMatrix(YD, d1);
  setYukawaMatrix(YE, e1);
}
  
// mix = 0 for all mixing in downs...at present this is the only possibility.
// Takes diagonal quark Yukawa matrices and mixes them up according to the CKM
// matrix assuming:
// mix=2, all mixing is in down sector
// mix=1, all mixing is in up sector
void MssmSusy::quarkMixing(const DoubleMatrix & CKM, int mix) {
  switch(mix) {
    case 1:       
      setYukawaMatrix(YU, CKM.transpose() * displayYukawaMatrix(YU) * CKM);
      break;
    case 2: 
      setYukawaMatrix(YD, CKM * displayYukawaMatrix(YD) * CKM.transpose()); 
      break;
     
    default:
    ostringstream ii;
    ii << "Error. MssmSusy::quarkMixing called with mix=" << mix;
    throw ii.str();
  }
}

void MssmSusy::getQuarkMixedYukawas(const QedQcd & r, const DoubleMatrix &
				    CKM, int mix, double vev) { 
  setDiagYukawas(r, vev);
  quarkMixing(CKM, mix);
}

// outputs object QedQcd & r valid at 1 GeV from SUSY data at mt, at
// present, only diagonal masses are handled. 
void MssmSusy::getMasses(QedQcd & r, double vev) const {
  double v1, v2;
  v1 = vev * cos(atan(displayTanb()));
  v2 = vev * sin(atan(displayTanb()));
  
  DoubleMatrix u1(displayYukawaMatrix(YU)), d1(displayYukawaMatrix(YD)),
    e1(displayYukawaMatrix(YE));
  r.setMass(mUp, v2 * u1(1, 1));
  r.setMass(mCharm, v2 * u1(2, 2));
  r.setMass(mTop, v2 * u1(3, 3));
  r.setMass(mDown, v1 * d1(1, 1));
  r.setMass(mStrange, v1 * d1(2, 2));
  r.setMass(mBottom, v1 * d1(3, 3));
  r.setMass(mElectron, v1 * e1(1, 1));
  r.setMass(mMuon, v1 * e1(2, 2));
  r.setMass(mTau, v1 * e1(3, 3));
}

#undef HR

// Rotates to quark mass basis, returning the mixing matrices defined as 
// yu_diag = vul yu vur^T  
// yd_diag = vdl yd vdr^T 
// All matrices should be 3 by 3
void MssmSusy::diagQuarkBasis(DoubleMatrix & vdl, DoubleMatrix & vdr, 
			DoubleMatrix & vul, DoubleMatrix & vur) const {
  DoubleMatrix u(3, 3), v(3, 3);
  DoubleVector ydDiag(3), yuDiag(3);
  displayYukawaMatrix(YU).diagonalise(u, v, yuDiag);
  vul = u.transpose(); vur = v.transpose();

  displayYukawaMatrix(YD).diagonalise(u, v, ydDiag);
  vdl = u.transpose(); vdr = v.transpose();
}

/// By Bednyakov, see arXiv:1009.5455
void MssmSusy::getThreeLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
				DoubleMatrix & gQQ, DoubleMatrix & gDD,
				DoubleMatrix & gUU, double & gH1H1, double &
				gH2H2, sBrevity & a) const {
  DoubleVector &gsq=a.gsq;
 
  /// powers of gauge couplings
  double a1 = gsq(1),  a2 = gsq(2),   a3 = gsq(3);   
  double a12 = a1*gsq(1),  a22 = a2*gsq(2),   a32 = a3*gsq(3);   
  double a13 = a12*gsq(1), a23 = a22*gsq(2),   a33 = a32*gsq(3);   

  const static double threelp = 2.53945672191370e-7; ///< 1/(16 pi^2)^3
  const static double kz = 7.21234141895757; ///< 6 Zeta(3)

  /// For calculational brevity
  /// NB Change notation to that of J&J  hep-ph/0408128 (Y->Y^T , etc)
  DoubleMatrix &u2t=a.u2, &d2t=a.d2, &e2t=a.e2;

  DoubleMatrix u2tu2t = u2t*u2t; 
  DoubleMatrix d2td2t = d2t*d2t; 
  DoubleMatrix e2te2t = e2t*e2t; 

  /// Everything gets the (1/16pi^2)^3 factor at the bottom
  DoubleMatrix ee(3, 3), ll(3, 3), qq(3, 3), dd(3, 3), uu(3, 3); 

  double ht = displayYukawaElement(YU, 3, 3), ht2 = sqr(ht), 
    ht4 = sqr(ht2), ht6 = ht2*ht4;
  double htau = displayYukawaElement(YE, 3, 3), htau2 = sqr(htau), 
    htau4 = sqr(htau2), htau6 = htau2*htau4;
  double hb = displayYukawaElement(YD, 3, 3), hb2 = sqr(hb), 
    hb4 = sqr(hb2), hb6 = hb2*hb4;
  
  const static double O45= .02222222222222222222 ;
  const static double O27= .03703703703703703703 ;
  const static double O6= .16666666666666666666 ;
  const static double O3375= .00029629629629629629 ;
  const static double O375= .00266666666666666666 ;
  const static double O1500= .00066666666666666666 ;
  const static double O75= .01333333333333333333 ;
  const static double O300= .00333333333333333333 ;
  const static double O9= .11111111111111111111 ;
  const static double O675= .00148148148148148148 ;
  const static double O150= .00666666666666666666 ;
  const static double O12= .08333333333333333333 ;
  const static double O450= .00222222222222222222 ;
  const static double O225= .00444444444444444444 ;
  const static double O180= .00555555555555555555 ;
  const static double O900= .00111111111111111111 ;
  const static double O60= .01666666666666666666 ;
  const static double O13500= .00007407407407407407 ;
  const static double O15= .06666666666666666666 ;
  const static double O3= .33333333333333333333 ;

  double dll3 = htau6 * (7. + kz) + hb2*ht2*htau2 * (+ 6) 
    + hb4*htau2 * (+ 9)
    + a3*hb2*htau2 * (- 32+ 8*kz)
    + a2*htau4 * (+ 3+ 3*kz)
    + a2*hb2*htau2 * (+ 18)
    + a22*htau2 * (- 18.75- 5.25*kz)
    + a22*ht2 * (- 22.5)
    + a22*hb2 * (- 22.5)
    + a22*a3 * (+ 90- 18*kz)
    + a23 * (+ 86.25+ 26.25*kz)
    + a1*htau4 * (+ 9- 1.8*kz)
    + a1*hb2*htau2 * (+ 8- 2*kz)
    + a1*a2*htau2 * (- 8.1+ 2.7*kz)
    + a1*a22 * (+ 2.25- 1.35*kz)
    + a12*htau2 * (- 30.15+ 0.27*kz)
    + a12*ht2 * (- 3.9)
    + a12*hb2 * (- 2.1)
    + a12*a3 * (+ 13.2- 2.64*kz)
    + a12*a2 * (+ 0.27- 0.81*kz)
    + a13 * (+ 18.39- 1.194*kz);

  double dtautau = 
    + htau6 * (+ 18+ 2*kz)
    + hb2*htau4 * (- 6)
    + hb2*ht2*htau2 * (+ 12)
    + hb4*htau2 * (+ 18)
    + a3*hb2*htau2 * (- 64+ 16*kz)
    + a2*htau4 * (+ 18- 6*kz)
    + a2*hb2*htau2 * (+ 27- 9*kz)
    + a22*htau2 * (- 43.5- 1.5*kz)
    + a1*htau4 * (+ 3.6+ 3.6*kz)
    + a1*hb2*htau2 * (+ 21.4+ 1.4*kz)
    + a1*a2*htau2 * (- 27+ 5.4*kz)
    + a12*htau2 * (- 40.86- 2.7*kz)
    + a12*ht2 * (- 15.6)
    + a12*hb2 * (- 8.4)
    + a12*a3 * (+ 52.8- 10.56*kz)
    + a12*a2 * (+ 16.2- 3.24*kz)
    + a13 * (+ 63.192- 4.776*kz);
  
  double dbb = 
    + hb2*htau4 * (+ 10)
    + hb2*ht2*htau2 * (- 2)
    + hb2*ht4 * (+ 18)
    + hb4*htau2 * (- 10)
    + hb4*ht2 * (+ 2)
    + hb6 * (+ 30+ 2*kz)
    + a3*hb2*htau2 * (+ 16)
    + a3*hb2*ht2 * (+ 64*O3)
    + a3*hb4 * (+ 16*O3+ 16*kz)
    + a32*ht2 * (- 80*O3)
    + a32*hb2 * (- 64*O3- 272*O9*kz)
    + a33 * (+ 2720*O27+ 160*O3*kz)
    + a2*hb2*htau2 * (+ 9- 3*kz)
    + a2*hb2*ht2 * (+ 9- 3*kz)
    + a2*hb4 * (+ 36- 12*kz)
    + a2*a3*hb2 * (- 88+ 16*kz)
    + a2*a32 * (+ 60- 12*kz)
    + a22*hb2 * (- 43.5- 1.5*kz)
    + a1*hb2*htau2 * (- 3+ kz)
    + a1*hb2*ht2 * (- 29*O15+ 0.6*kz)
    + a1*hb4 * (+ 20*O3- 0.8*kz)
    + a1*a3*hb2 * (- 4.8+ 16*O9*kz)
    + a1*a32 * (+ 452*O45- 44*O15*kz)
    + a1*a2*hb2 * (- 8.6+ 1.4*kz)
    + a12*htau2 * (- 1.2)
    + a12*ht2 * (- 26*O15)
    + a12*hb2 * (- 73*O6- 7*O450*kz)
    + a12*a3 * (+ 728*O225- 88*O75*kz)
    + a12*a2 * (+ 1.8- 0.36*kz)
    + a13 * (+ 5629*O675- 199*O375*kz);
  
  double dtt = 
    + ht6 * (+ 30+ 2*kz)
    + hb2*ht2*htau2 * (+ 4)
    + hb2*ht4 * (+ 2)
    + hb4*ht2 * (+ 18)
    + a3*ht4 * (+ 16*O3+ 16*kz)
    + a3*hb2*ht2 * (+ 64*O3)
    + a32*ht2 * (- 64*O3- 272*O9*kz)
    + a32*hb2 * (- 80*O3)
    + a33 * (+ 2720*O27+ 160*O3*kz)
    + a2*ht4 * (+ 36- 12*kz)
    + a2*hb2*ht2 * (+ 9- 3*kz)
  + a2*a3*ht2 * (- 88+ 16*kz)
    + a2*a32 * (+ 60- 12*kz)
    + a22*ht2 * (- 43.5- 1.5*kz)
    + a1*ht4 * (+ 20*O3+ 2.4*kz)
    + a1*hb2*ht2 * (+ 19*O15+ 0.6*kz)
    + a1*a3*ht2 * (- 8*O15- 112*O45*kz)
    + a1*a32 * (- 172*O45- 44*O15*kz)
    + a1*a2*ht2 * (- 13.4+ 2.6*kz)
    + a12*htau2 * (- 4.8)
    + a12*ht2 * (- 3437*O150- 247*O450*kz)
    + a12*hb2 * (- 56*O15)
    + a12*a3 * (+ 2144*O225- 352*O75*kz)
    + a12*a2 * (+ 7.2- 1.44*kz)
    + a13 * (+ 106868*O3375- 796*O375*kz);
  
  double dqq3 = 
    + ht6 * (+ 15+ kz)
    + hb2*htau4 * (+ 5)
    + hb2*ht4 * (+ 10)
    + hb4*htau2 * (- 4)
    + hb4*ht2 * (+ 10)
    + hb6 * (+ 15+ kz)
    + a3*ht4 * (+ 40*O3+ 8*kz)
    + a3*hb2*htau2 * (+ 8)
    + a3*hb4 * (+ 40*O3+ 8*kz)
    + a32*ht2 * (- 24- 136*O9*kz)
    + a32*hb2 * (- 24- 136*O9*kz)
    + a33 * (+ 2720*O27+ 160*O3*kz)
    + a2*ht4 * (+ 15+ 3*kz)
    + a2*hb2*htau2 * (+ 6)
    + a2*hb4 * (+ 15+ 3*kz)
    + a2*a3*ht2 * (- 4)
    + a2*a3*hb2 * (- 4)
    + a2*a32 * (+ 8- 12*kz)
    + a22*htau2 * (- 7.5)
    + a22*ht2 * (- 33.75- 5.25*kz)
    + a22*hb2 * (- 33.75- 5.25*kz)
    + a22*a3 * (+ 50- 18*kz)
    + a23 * (+ 86.25+ 26.25*kz)
    + a1*ht4 * (+ 17*O3- 1.8*kz)
    + a1*hb2*htau2 * (- 1.6+ 0.4*kz)
    + a1*hb4 * (+ 11*O3- kz)
    + a1*a3*ht2 * (- 13.6+ 64*O45*kz)
    + a1*a3*hb2 * (- 76*O15+ 64*O45*kz)
    + a1*a32 * (+ 608*O45- 44*O15*kz)
    + a1*a2*ht2 * (- 5.9+ 1.5*kz)
    + a1*a2*hb2 * (- 4.1+ 0.3*kz)
    + a1*a2*a3 * (- 1.6)
    + a1*a22 * (+ 6.25- 1.35*kz)
    + a12*htau2 * (- 0.3)
    + a12*ht2 * (- 12.99+ 143*O900*kz)
    + a12*hb2 * (- 1969*O300+ 7*O180*kz)
    + a12*a3 * (+ 194*O225- 22*O75*kz)
    + a12*a2 * (+ 0.11- 0.09*kz)
    + a13 * (+ 28457*O13500- 199*O1500*kz);

  double h2h2 = 
    + ht6 * (+ 57+ 3*kz)
    + hb2*ht2*htau2 * (+ 6)
    + hb4*ht2 * (+ 27)
    + a3*ht4 * (+ 72- 24*kz)
    + a3*hb2*ht2 * (+ 24- 8*kz)
    + a32*ht2 * (- 160*O3- 8*O3*kz)
    + a2*ht4 * (+ 9+ 9*kz)
    + a2*hb2*ht2 * (+ 18)
    + a2*a3*ht2 * (- 132+ 24*kz)
    + a22*htau2 * (- 7.5)
    + a22*ht2 * (- 56.25- 15.75*kz)
    + a22*hb2 * (- 22.5)
    + a22*a3 * (+ 90- 18*kz)
    + a23 * (+ 86.25+ 26.25*kz)
    + a1*ht4 * (+ 11.4- 0.6*kz)
    + a1*hb2*ht2 * (+ 1.2+ 0.2*kz)
    + a1*a3*ht2 * (- 124*O3+ 104*O15*kz)
    + a1*a2*ht2 * (- 5.7+ 2.1*kz)
    + a1*a22 * (+ 2.25- 1.35*kz)
    + a12*htau2 * (- 2.7)
    + a12*ht2 * (- 2123*O60- 13*O60*kz)
    + a12*hb2 * (- 2.1)
    + a12*a3 * (+ 13.2- 2.64*kz)
    + a12*a2 * (+ 0.27- 0.81*kz)
    + a13 * (+ 18.39- 1.194*kz);

  double h1h1 = 
    + htau6 * (+ 7+ kz)
    + hb2*htau4 * (+ 18)
    + hb2*ht4 * (+ 27)
    + hb4*htau2 * (+ 18)
    + hb6 * (+ 57+ 3*kz)
    + a3*hb2*ht2 * (+ 24- 8*kz)
    + a3*hb4 * (+ 72- 24*kz)
    + a32*hb2 * (- 160*O3- 8*O3*kz)
    + a2*htau4 * (+ 3+ 3*kz)
    + a2*hb2*ht2 * (+ 18)
    + a2*hb4 * (+ 9+ 9*kz)
    + a2*a3*hb2 * (- 132+ 24*kz)
    + a22*htau2 * (- 18.75- 5.25*kz)
    + a22*ht2 * (- 22.5)
    + a22*hb2 * (- 56.25- 15.75*kz)
    + a22*a3 * (+ 90- 18*kz)
    + a23 * (+ 86.25+ 26.25*kz)
    + a1*htau4 * (+ 9- 1.8*kz)
    + a1*hb2*ht2 * (- 2.4+ 1.4*kz)
    + a1*hb4 * (+ 3+ 1.8*kz)
    + a1*a3*hb2 * (- 284*O15+ 56*O15*kz)
    + a1*a2*htau2 * (- 8.1+ 2.7*kz)
    + a1*a2*hb2 * (- 0.3- 1.5*kz)
    + a1*a22 * (+ 2.25- 1.35*kz)
    + a12*htau2 * (- 30.15+ 0.27*kz)
    + a12*ht2 * (- 3.9)
    + a12*hb2 * (- 175*O12- 77*O300*kz)
    + a12*a3 * (+ 13.2- 2.64*kz)
    + a12*a2 * (+ 0.27- 0.81*kz)
    + a13 * (+ 18.39- 1.194*kz);
  
  gLL(3,3) = gLL(3,3) + threelp * dll3;
  gEE(3,3) = gEE(3,3) + threelp * dtautau;
  gQQ(3,3) = gQQ(3,3) + threelp * dqq3;
  gDD(3,3) = gDD(3,3) + threelp * dbb;
  gUU(3,3) = gUU(3,3) + threelp * dtt;
  gH1H1 = gH1H1 + threelp * h1h1;
  gH2H2 = gH2H2 + threelp * h2h2;
}

} // namespace softsusy
