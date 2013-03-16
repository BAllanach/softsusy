/** \file softpars.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: includes all MSSM parameters in the Lagrangian

*/

#include "nmssmsoftpars.h"

const SoftParsNmssm & SoftParsNmssm::operator=(const SoftParsNmssm & s) {
  if (this == &s) return *this;
  SoftPars<nMssmSusy, nmsBrevity>::operator=(s);
  alambda = s.alambda;
  akappa = s.akappa;
  mSsq = s.mSsq;
  mSpsq = s.mSpsq;
  zeta_s = s.zeta_s;
  return *this;
}

double SoftParsNmssm::displaySoftAlambda() const {
   if (fabs(displayLambda()) < 1.0e-100) {
    ostringstream ii;
    ii << "WARNING: asking for SoftParsNmssm::displaySoftAlambda() where Yukawa coupling is " <<
       fabs(displayLambda()) << endl;
    throw ii.str();
  } 
  double temp = displayTrialambda() / displayLambda();
   return temp;
}

double SoftParsNmssm::displaySoftAkappa() const {
   if (fabs(displayKappa()) < 1.0e-100) {
    ostringstream ii;
    ii << "WARNING: asking for SoftParsNmssm::displaySoftAkappa() where Yukawa coupling is " <<
       fabs(displayKappa()) << endl;
    throw ii.str();
  } 
   double temp = displayTriakappa() / displayKappa();
   return temp;
}

const DoubleVector SoftParsNmssm::display() const {
  DoubleVector y(nMssmSusy::display());
  y.setEnd(numSoftParsNmssm);
  int i, j, k=numNMssmPars;
  for (i=1; i<=3; i++) {
    k++;
    y(k) = displayGaugino(i);
  }
  
  for (i=1; i<=3; i++)    
    for (j=1; j<=3; j++) {
      k++;
      y(k) = displayTrilinear(UA, i, j);
      y(k+9) = displayTrilinear(DA, i, j);
      y(k+18) = displayTrilinear(EA, i, j);
      y(k+27) = displaySoftMassSquared(mQl, i, j);
      y(k+36) = displaySoftMassSquared(mUr, i, j);
      y(k+45) = displaySoftMassSquared(mDr, i, j);
      y(k+54) = displaySoftMassSquared(mLl, i, j);
      y(k+63) = displaySoftMassSquared(mEr, i, j);
    }
  y(k+64) = displayM3Squared();
  y(k+65) = displayMh1Squared();
  y(k+66) = displayMh2Squared();
  y(k+67) = mSsq;
  y(k+68) = mSpsq;
  y(k+69) = zeta_s;
  y(k+70) = alambda;
  y(k+71) = akappa;

  return y;
}

void SoftParsNmssm::setTrialambda(double al) {
   alambda = al;
}

void SoftParsNmssm::setTriakappa(double ak) {
   akappa = ak;
}

void SoftParsNmssm::set(const DoubleVector & y) {
  nMssmSusy::set(y);
  int i, j, k=numNMssmPars;
  for (i=1; i<=3; i++) {
    k++;
    setGauginoMass(i, y.display(k));
  }
  
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      k++;
      setTrilinearElement(UA, i, j, y.display(k));
      setTrilinearElement(DA, i, j, y.display(k+9));
      setTrilinearElement(EA, i, j, y.display(k+18));
      setSoftMassElement(mQl, i, j, y.display(k+27));
      setSoftMassElement(mUr, i, j, y.display(k+36));
      setSoftMassElement(mDr, i, j, y.display(k+45));
      setSoftMassElement(mLl, i, j, y.display(k+54));
      setSoftMassElement(mEr, i, j, y.display(k+63));
    }
  setM3Squared(y.display(k+64));
  setMh1Squared(y.display(k+65));
  setMh2Squared(y.display(k+66));
  mSsq = y.display(k+67);
  mSpsq = y.display(k+68);
  zeta_s = y.display(k+69);
  alambda =  y.display(k+70);
  akappa = y.display(k+71);
}


// Outputs derivatives (DRbar scheme) in the form of dsoft
// thresholds = 0 and NOTHING is decoupled.
// thresholds = 2 and SUSY params/gauge couplings are decoupled at sparticle
// thresholds.
// CHECKED: 24/05/02
SoftParsNmssm SoftParsNmssm::beta2() const {
 // Constants for gauge running
   static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3), clBeta(3);
  static DoubleMatrix babBeta(3, 3);
  
  if (bBeta(1) < 1.0e-5) // Constants not set yet
     nmsetBetas(babBeta, cuBeta, cdBeta, ceBeta, clBeta, bBeta);
  
  // For calculational brevity: 
  static nmsBrevity a;
  // convert to beta functions
  static nMssmSusy dsb;
  
  // calculate derivatives for full SUSY spectrum. Brevity calculations come
  // out encoded in a
  dsb = nMssmSusy::beta(a);
  SoftPars<nMssmSusy, nmsBrevity> base(SoftPars<nMssmSusy, nmsBrevity>::beta2());
  
  // To keep this a const function: TIME SAVINGS
  DoubleMatrix &u1=a.u1, &d1=a.d1, &e1=a.e1;
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2,
    &u2t=a.u2t, &d2t=a.d2t, &e2t=a.e2t, &dt=a.dt, &ut=a.ut, &et=a.et;      
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  double  &lsq = a.lsq, &ksq = a.ksq, &l4 = a.l4, &k4 =a.k4; 
  DoubleVector &gsq=a.gsq;
  const double &lam = displayLambda();
  const double &kap = displayKappa();
  // Best to make these const DoubleMatrix & type
  const DoubleMatrix &hu = displayTrilinear(UA), &hd = displayTrilinear(DA),
    &he = displayTrilinear(EA);
  const double &hlam = displayTrialambda();
 const double &hkap = displayTriakappa();
 const double &mu_s = displayMu_s();
 const double &zeta = displayZeta();  
 const double &smu = displaySusyMu();
  static DoubleMatrix hut(3, 3), hdt(3, 3), het(3, 3);
  static DoubleMatrix hu2(3, 3), hd2(3, 3), he2(3, 3);
  static DoubleMatrix hu2t(3, 3), hd2t(3, 3), he2t(3, 3);
  const DoubleMatrix &mq = displaySoftMassSquared(mQl),
    &mu = displaySoftMassSquared(mUr), &md = displaySoftMassSquared(mDr),
    &ml = displaySoftMassSquared(mLl), &me = displaySoftMassSquared(mEr); 
  static DoubleVector mG(displayGaugino()); 
  static DoubleVector msq(1, 3), gsqM(1, 3), gMsq(1, 3);
  

  hut = hu.transpose(); hdt = hd.transpose(); het = he.transpose();
  hu2 = hu * hut; hd2 = hd * hdt; he2 = he * het;
  double huT = hu2.trace(), hdT = hd2.trace(), heT = he2.trace();
  hu2t = hut * hu; hd2t = hdt * hd; he2t = het * he;
  double mqT = mq.trace(), muT = mu.trace(), mdT = md.trace(), meT =
    me.trace(), mlT = ml.trace();
  
  double huuT = (hu * ut).trace(), hddT = (hd * dt).trace(), 
    heeT = (he * et).trace(); 
  mG = displayGaugino(); msq = mG * mG; gsqM = gsq * mG; gMsq = gsq * msq;
 
  const double m3sq = displayM3Squared();
  const double mH1sq = displayMh1Squared();
  const double mH2sq = displayMh2Squared();
  
  // derivatives of soft parameters
  static DoubleVector dmG(1, 3);
  static double dmH1sq, dmH2sq, dm3sq;
  static double dmSsq=0, dmSpsq=0, dz_s=0;
  static DoubleMatrix dmq(3, 3), dmu(3, 3), dmd(3, 3), dme(3, 3), dml(3, 3);
  static DoubleMatrix dhu(3, 3), dhd(3, 3), dhe(3, 3);
  static double dhlam=0, dhkap=0;
  
//PA:  defined with an additional factor of \lambda^2  compared to literature to avoid division and ensure running can be tested with lambda = 0.
  //PA: for speed issues may wish to set differently in depending on value of mixing
  double Mlamsq = lsq * (mH2sq + mH1sq + mSsq) + sqr(hlam); 		   
  
  double Mkapsq = 3.0 * ksq * mSsq + sqr(hkap);			   
    // double Xu = (u2 * mq).trace() + (u2t * mu).trace() + mH2sq * uuT + huT; 
    // double Xd = mH1sq * ddT + (d2 * mq).trace() + (d2t * md).trace() + hdT;
    // double Xe = mH1sq * eeT + (e2 * ml).trace() + (e2t * me).trace() + heT;
 double Xu = (mq * u2).trace() + (mu * u2t).trace() + mH2sq * uuT + huT; 
    double Xd = mH1sq * ddT + (mq * d2).trace() + (md * d2t).trace() + hdT;
    double Xe = mH1sq * eeT + (ml * e2).trace() + (me * e2t).trace() + heT;


    //PA: these trace factors appear frequently, so defined here to simplify expressions 
    double Ytr = 3.0 * uuT + 3.0 * ddT + eeT;
    double aYtr = 3.0 * huuT + 3.0 * hddT + heeT;

  const double ONEO16Pisq = 1.0 / (16.0 * sqr(PI));
  if (displayLoops() > 0) {
    const static double sixteenO3 = 16.0 / 3.0, oneO15 = 1.0 /
      15.0;
    
    dmG = gsqM * 2.0 * bBeta;
    //PA: checked 14/9/2012
    double curlyS = mH2sq - mH1sq + mqT - mlT - 2.0 * muT + mdT + meT;
    
    dm3sq = base.displayM3Squared();

    const double dm3sqnm = + 4.0 * displaySusyMu() * lam * hlam + 6.0 * lsq * m3sq + 2.0 * lam * kap * mSpsq;
    dm3sq += dm3sqnm * ONEO16Pisq;
   
    dmH1sq = 2.0 * 
      (-0.3 * gsq(1) * curlyS - 0.6 * gMsq(1) - 3.0 * gMsq(2) +
       3.0 * (mH1sq * ddT + (d2 * mq).trace() + (d2t * md).trace() + hdT) +
       (mH1sq * eeT + (e2 * ml).trace() + (e2t * me).trace() + heT))
       + 2.0 * Mlamsq;
    

    dmH2sq = 2.0 * 
      (0.3 * gsq(1) * curlyS - 0.6 * gMsq(1) - 3.0 * gMsq(2) +
       3.0 * (mH2sq * uuT + (u2 * mq).trace() + (u2t * mu).trace() + huT))
       + 2.0 * Mlamsq;
  
    dmSsq = 4.0 * Mlamsq + 4.0 * Mkapsq;
   
    
    dmSpsq = 4.0 * lam * (lam * mSpsq + 2.0 * mu_s * hlam) 
       + 8.0 * kap * (kap * mSpsq + mu_s * hkap) + 8.0 * lam * kap * m3sq;
   
    //PA: checked numerically 14/9/2012

    dmq = 2.0 *
      (0.1 * gsq(1) * curlyS - oneO15 * gMsq(1) - 3.0 * gMsq(2) -
       sixteenO3 * gMsq(3) +
       0.5 * (u2 * mq + mq * u2 + 2.0 * (u1 * mu * ut + mH2sq * u2 + hu2))
       + 0.5 * (d2 * mq + mq * d2 + 2.0 * (d1 * md * dt + mH1sq * d2 +
					   hd2)));
    
    dml = 2.0 *
      ( -0.3 * gsq(1) * curlyS - 0.6 * gMsq(1) - 3.0 * gMsq(2) + 
	0.5 * (e2 * ml + ml * e2 + 2.0 * (e1 * me * et + mH1sq * e2 + he2)));
    
    dmu = 2.0 *
      (-0.4 * gsq(1) * curlyS - 16.0 * oneO15 * gMsq(1) - sixteenO3 *
       gMsq(3) + 
       (u2t * mu + mu * u2t + 2.0 * (ut * mq * u1 + mH2sq * u2t + hu2t)));
    
    dmd = 2.0 * (0.2 * gsq(1) * curlyS 
		 - 4.0 * oneO15 * gMsq(1) - sixteenO3 * gMsq(3) + 
		 (d2t * md + md * d2t + 2.0 * (dt * mq * d1 + mH1sq * d2t +
					       hd2t)));
    
    dme = 2.0 *
      (0.6 * gsq(1) * curlyS - 36.0 * oneO15 * gMsq(1) +
       e2t * me + me * e2t + 2.0 * (et * ml * e1 + mH1sq * e2t + he2t));
    
    dhu = (3.0 * uuT 
	   - sixteenO3 * gsq(3) - 3.0 * gsq(2) - 13.0 * oneO15 * gsq(1) + lsq)
      * hu
      + 2.0 * (3.0 * huuT + 13.0 * oneO15 * gsqM(1) + 3 * gsqM(2) + sixteenO3
	     * gsqM(3)) * u1 + 
      4.0 * hu * u2t + 5.0 * u2 * hu + 2.0 * hd * dt * u1 + d2 * hu 
       + 2.0 * lam * hlam * u1;
    
    dhd = (eeT + 3.0 * ddT -sixteenO3 * gsq(3) - 3.0 * gsq(2) - 
	   7.0 * oneO15 * gsq(1) + lsq) * hd 
      + 2.0 * (7.0 * oneO15 * gsqM(1) + 3.0 * gsqM(2) + sixteenO3
	     * gsqM(3) + heeT + 3.0 * hddT) * d1 +
      4.0 * hd * d2t + 5.0 * d2 * hd + 
      2.0 * hu * ut * d1 + u2 * hd + 2.0 * lam * hlam * d1;
    
    dhe = (eeT + 3.0 * ddT - 3.0 * gsq(2) - 1.8 * gsq(1) + lsq) * he
      + (3.6 * gsqM(1) + 6.0 * gsqM(2) + 2.0 * heeT + 6.0 * hddT) * e1 +
      4.0 * he * e2t + 5.0 * e2 * he + 2.0 * lam * hlam * e1;

    dhlam = (Ytr + 4.0 * lsq + 2.0 * ksq
             - 0.6 * gsq(1) - 3.0 * gsq(2)) * hlam
       + lam * ( 2.0 * aYtr  + 8.0 * hlam * lam
                 + 4.0 * hkap * kap + 1.2 * gsqM(1) + 6.0 * gsqM(2) );
   
    dhkap = 18.0 * hkap * ksq + 12.0 * lam * kap * hlam 
       + 6.0 * hkap * lsq;
  
    //tadpole
    dz_s = 2.0 * zeta_s * (lsq +  ksq) 
       + 2.0 * mu_s * (2.0 * lam * m3sq + kap * mSpsq)
       + 4.0 * zeta * (lam * hlam + kap * hkap)
       + 4.0 * lam * smu * (mH2sq + mH1sq) + 4.0 * kap * mu_s * mSsq
       + 4.0 * hlam * m3sq + 2.0 * hkap * mSpsq;

    // convert to proper derivatives: 
    dmG    = ONEO16Pisq * dmG;
    dmH1sq *= ONEO16Pisq;
    dmH2sq *= ONEO16Pisq;
    dmSsq  *= ONEO16Pisq;
    dmSpsq *= ONEO16Pisq;
    dmq    *= ONEO16Pisq;
    dml    *= ONEO16Pisq;
    dmu    *= ONEO16Pisq;
    dmd    *= ONEO16Pisq;
    dme    *= ONEO16Pisq;
    dhu    *= ONEO16Pisq;
    dhd    *= ONEO16Pisq;
    dhe    *= ONEO16Pisq;
    dhkap  *= ONEO16Pisq;
    dhlam  *= ONEO16Pisq;
    dz_s   *= ONEO16Pisq; 
  }
  
  // two-loop contributions. I got these from hep-ph/9311340. WIth respect to
  // their notation, the Yukawas and h's are TRANSPOSED. Gaugino masses are
  // identical, as are soft masses. B(mV) = mu(BBO) B(BBO)
  if (displayLoops() > 1) {
    static DoubleVector dmG2(1, 3);
    static DoubleMatrix dmq2(3, 3), dmu2(3, 3), dmd2(3, 3), dme2(3, 3), 
      dml2(3, 3);
    static DoubleMatrix dhu2(3, 3), dhd2(3, 3), dhe2(3, 3);
    static DoubleVector sigma(3);
    static double dhlam2, dhkap2; 
    const double oneO16Pif = sqr(ONEO16Pisq);
    
    DoubleVector &g4=a.g4;

    double mq3 = displaySoftMassSquared(mQl, 3, 3);
    double mu3 = displaySoftMassSquared(mUr, 3, 3);
    double md3 = displaySoftMassSquared(mDr, 3, 3);
    double ml3 = displaySoftMassSquared(mLl, 3, 3);
    double me3 = displaySoftMassSquared(mEr, 3, 3);
    double ht = displayYukawaElement(YU, 3, 3), ht2 = sqr(ht), ht4 = sqr(ht2);
    double htau = displayYukawaElement(YE, 3, 3), htau2 = sqr(htau), 
      htau4 = sqr(htau2);
    double hb = displayYukawaElement(YD, 3, 3), hb2 = sqr(hb), hb4 = sqr(hb2);
    double Ut = displayTrilinear(UA, 3, 3), Ut2 = sqr(Ut);
    double Ub = displayTrilinear(DA, 3, 3), Ub2 = sqr(Ub);
    double Utau = displayTrilinear(EA, 3, 3), Utau2 = sqr(Utau);
   
    double sP;
    if (MIXING < 1) {
      sP = ///< dominant 3rd family approximation
	-(3.0 * mH2sq + mq3 - 4.0 * mu3) * ht2 +  
	(3.0 * mH1sq - mq3 - 2.0 * md3) * hb2 + 
	(mH1sq + ml3 - 2.0 * me3) * htau2 +
	(1.5 * gsq(2) + 0.3 * gsq(1)) * (mH2sq - mH1sq - mlT) +
	(8.0 / 3.0 * gsq(3) + 1.5 * gsq(2) + gsq(1) / 30.0) * mqT -
	(16.0 / 3.0 * gsq(3) + 16.0 / 15.0 * gsq(1)) * muT +
	(8.0 / 3.0 * gsq(3) + 2.0 / 15.0 * gsq(1)) * mdT +       
         1.2 * gsq(1) * meT + lsq * (mH1sq - mH2sq);  //checked
    } else {
      sP = ///< fully mixed case -- slower
	-(3.0 * mH2sq * uuT + (mq * u2).trace()) 
	+ 4.0 * (u1 * mu * ut).trace() + 
	3.0 * mH1sq * ddT - (mq * d2).trace() - 2.0 * (d1 * md * dt).trace() + 
	mH1sq * eeT + (ml * e2).trace() -
	(2.0 * e1 * me * et).trace() +
	(1.5 * gsq(2) + 0.3 * gsq(1)) * (mH2sq - mH1sq - mlT) +
	(8.0 / 3.0 * gsq(3) + 1.5 * gsq(2) + gsq(1) / 30.0) * mqT -
	(16.0 / 3.0 * gsq(3) + 16.0 / 15.0 * gsq(1)) * muT +
	(8.0 / 3.0 * gsq(3) + 2.0 / 15.0 * gsq(1)) * mdT +       
         1.2 * gsq(1) * meT + lsq * (mH1sq - mH2sq);  //checked
      }


    sigma(1) = 0.2 * gsq(1) *
      (3.0 * (mH1sq + mH2sq) + mqT + 3.0 * mlT + 8.0 * muT + 2.0 * mdT + 6.0 *
       meT);
    sigma(2) = gsq(2) * (mH1sq + mH2sq + (3.0 * mqT + mlT));
    sigma(3) = gsq(3) * (2.0 * mqT + muT + mdT); // these 3 checked

    DoubleMatrix u4(u2 * u2), d4(d2 * d2), e4(e2 * e2), u4t(u2t * u2t), d4t(d2t * d2t);

    if (INCLUDE_2_LOOP_SCALAR_CORRECTIONS) {
    /* new dominant 3-family version */
      if (MIXING < 1) {
	dmq2 = -(2.0 * mq + 8.0 * mH2sq) * u4 - 
	  4.0 * u1 * mu * u2t * ut - 
	  4.0 * u2 * mq * u2 - 4.0 * u2 * u1 * mu * ut - 2.0 * u4 * mq -
	  (2.0 * mq + 8.0 * mH1sq) * d4 - 4.0 * d1 * md * d2t * dt - 
	  4.0 * d2 * mq * d2 - 4.0 * d2 * d1 * md * dt - 2.0 * d4 * mq -
	  ((mq + 4.0 * mH2sq) * u2 + 2.0 * u1 * mu * ut + u2 * mq) * uuT 
	  * 3.0 -
	  ((mq + 4.0 * mH1sq) * d2 + 2.0 * d1 * md * dt + d2 * mq) * 
	  (3.0 * ddT + eeT) -
	  6.0 * u2 * ((mq * u2).trace() + (u1 * mu * ut).trace()) -
	  d2 * ((6.0 * mq * d2).trace() + (6.0 * d1 * md * dt).trace() + 
		(2.0 * ml * e2).trace() + (2.0 * e1 * me * et).trace()) -
	  4.0 * (u2 * hu2 + hu2 * u2 + u1 * hu2t * ut + hu * u2t * hut) -
	  4.0 * (d2 * hd2 + hd2 * d2 + d1 * hd2t * dt + hd * d2t * hdt) -
	  hu2 * uuT * 6.0 - u2 * 6.0 * Ut2 - hu * ut * 6.0 * huuT -
	  u1 * hut * 6.0 * huuT - hd2 * (6.0 * ddT + 2.0 * eeT) - 
	  d2 * (6.0 * Ub2 + 2.0 * Utau2) -
	  hd * dt * (6.0 * hddT + 2.0 * heeT) - 
	  d1 * hdt * (6.0 * hddT + 2.0 * heeT) + 0.4 * gsq(1) * 
	  ((2.0 * mq +  4.0 * mH2sq) * u2 + 4.0 * u1 * mu * ut +
	   2.0 * u2 * mq + 4.0 * hu2 - 4.0 * mG(1) * hu * ut - 4.0
	   * mG(1) * u1 * hut + 8.0 * msq(1) * u2 +
	   (mq + 2.0 * mH1sq) * d2 + 2.0 * d1 * md * dt +
	   d2 * mq + 2.0 * hd2 - 2.0 * mG(1) * hd * dt - 
	   2.0 * mG(1) * d1 * hdt + 4.0 * msq(1) * d2) +
	  0.4 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 32.0 * gsq(3) *
	  gsq(2) * (msq(3) + msq(2) + mG(2) * mG(3)) +
	  32.0 / 45.0 * gsq(3) * gsq(1) * 
	  (msq(3) + msq(1) + mG(1) * mG(3)) + 33.0
	  * g4(2) * msq(2) + 0.4 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) *
						      mG(2)) +
	  199.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
	  3.0 * gsq(2) * sigma(2) + gsq(1) * sigma(1) / 15.0 // checked
         - lsq * (2.0 * ut * mu * u1 + mq * u2 + u2 * mq + 2.0 * mH2sq * u2
                     + 2.0 * hu2) - 2.0 * Mlamsq * u2  
           - 2.0 * lam * hlam * (u1 * hut + hu * ut)
            - lsq * (2.0 * dt * md * d1 + mq * d2 + d2 * mq + 2.0 * mH1sq * d2
                     + 2.0 * hd2) - 2.0 * Mlamsq * d2  
            - 2.0 * lam * hlam * (d1 * hdt + hd * dt);

	dml2 = -(2.0 * ml + 8.0 * mH1sq) * e2 * e2 - 4.0 * e1 * me * e2t * et -
	  4.0 * e2 * ml * e2 - 4.0 * e2 * e1 * me * et - 2.0 * e2 * e2 * ml -
	  ((ml + 4.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml) * 
	  (3.0 * ddT + eeT) -
	  e2 * (6.0 * mq * d2 + 6.0 * d1 * md * dt + 2.0 * ml * e2 + 
		2.0 * e1 * me * et).trace() -
	  4.0 * (e2 * he2 + he2 * e2 + e1 * he2t * et + he * e2t * het) -
	  he2 * (6.0 * ddT + 2.0 * eeT) - e2 * (6.0 * Ub2 + 2.0 * Utau2) -
	  he * et * (6.0 * hddT + 2.0 * heeT) - 
	  e1 * het * (6.0 * hddT + 2.0 * heeT) +
	  1.2 * gsq(2) * 
	  ((ml + 2.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml +
	   2.0 * he2 - 2.0 * mG(1) * he * et - 2.0 * mG(1) * e1 *
			  het + 4.0 * msq(1) * e2) -
	  1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 3.6 * gsq(2) * gsq(1) *
	  (msq(2) + msq(1) + mG(1) * mG(2)) + 621.0 / 25.0 * g4(1) * msq(1) +
	  3.0 * gsq(2) * sigma(2) + 0.6 * gsq(1) * sigma(1) //checked
           - lsq * (2.0 * et * me * e1 + ml * e2 + e2 * ml 
                          + 2.0 * mH1sq * e2 + 2.0 * he2) - 2.0 * Mlamsq * e2  
           - 2.0 * lam * hlam * (e1 * het + he * et);
	
dmu2 = -(2.0 * mu + 8.0 * mH2sq) * u4t - 4.0 * ut * mq * u2 * u1 
	  - 4.0 * u2t * mu * u2t - 4.0 * u2t * ut * mq * u1 - 
	  2.0 * u4t * mu - (2.0 * mu + 4.0 * mH2sq + 4.0 * mH1sq) 
	  * ut * d2 * u1 - 4.0 * ut * mq * d2 * u1 - 
	  4.0 * ut * d1 * md * dt * u1 -
	  4.0 * ut * d2 * mq * u1 - 2.0 * ut * d2 * u1 * mu -
	  ((mu + 4.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu) * 
	  6.0 * uuT -
	  12.0 * u2t * (mq * u2 + u1 * mu * ut).trace() -
	  4.0 * (hu2t * u2t + u2t * hu2t + hut * u2 * hu + ut * hu2 * u1) -
	  4.0 * (hut * hd * dt * u1 + ut * d1 * hdt * hu + hut * d2 * hu + 
	     ut * hd2 * u1) -
	  12.0 * (hu2t * uuT + u2t * Ut2 + hut * u1 * huuT +
		  ut * hu * huuT) +
	  (6.0 * gsq(2) - 0.4 * gsq(1)) * 
	  ((mu + 2.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu + 
	   2.0 * hu2t) +
	  12.0 * gsq(2) * (2.0 * msq(2) * u2t - mG(2) * hut * u1 - 
			   mG(2) * ut * hu) -
	  0.8 * gsq(1) * (2.0 * msq(1) * u2t - mG(1) * hut * u1 - 
			  mG(1) * ut * hu)-
	  1.6 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 
	  512.0 / 45.0 * gsq(1) * gsq(3) * (msq(3) + msq(1) + mG(1) * mG(3)) +
	  3424.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
	  16.0 / 15.0 * gsq(1) * sigma(1) // checked
         - 2.0 * lsq * (2.0 * ut * mq * u1 + mu * u2 + u2 * mu 
                        + 2.0 * mH2sq * u2 + 2.0 * hu2) - 4.0 * Mlamsq * u2  
           - 4.0 * lam * hlam * (ut * hu + hut * u1);
	dmd2 = -(2.0 * md + 8.0 * mH1sq) * d4t - 
	  4.0 * dt * mq * d1 * d2t -
	  4.0 * d2t * md * d2t - 4.0 * d2t * dt * mq * d1 - 
	  2.0 * d4t * md -
	  (2.0 * md + 4.0 * mH2sq + 4.0 * mH1sq) * dt * u2 * d1 -
	  4.0 * dt * mq * u2 * d1 - 4.0 * dt * u1 * mu * ut * d1 - 
	  4.0 * dt * u2 * mq * d1 - 2.0 * dt * u2 * d1 * md - 
	  ((md + 4.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md) * 
	  (6.0 * ddT + 2.0 * eeT) -
	  4.0 * d2t * (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2 
		       + e1 * me * et).trace() -
	  4.0 * (hd2t * d2t + d2t * hd2t + hdt * d2 * hd + dt * hd2 * d1) -
	  4.0 * (hdt * hu * ut * d1 + dt * u1 * hut * hd + hdt * u2 * hd + 
		 dt * hu2 * d1) -
	  4.0 * hd2t * (3.0 * ddT + eeT) - 4.0 * d2t * (3.0 * hdT + heT) -
	  4.0 * hdt * d1 * (3.0 * hddT + heeT) - 
	  4.0 * dt * hd * (3.0 * hddT + heeT) +
	  (6.0 * gsq(2) + 0.4 * gsq(1)) * 
	  ((md + 2.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md + 
	   2.0 * hd2t) +
	  12.0 * gsq(2) * (2.0 * msq(2) * d2t - mG(2) * hdt * d1 
			   - mG(2) * dt * hd) +
	  0.8 * gsq(1) * (2.0 * msq(1) * d2t - mG(1) * hdt * d1 - 
			  mG(1) * dt * hd)+
	  0.8 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 
	  128.0 / 45.0 * gsq(3) * gsq(1) * (msq(3) + msq(1) + mG(1) * mG(3)) +
	  808.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
	  4.0 / 15.0 * gsq(1) * sigma(1)// checked
           - 2.0 * lsq * (2.0 * dt * mq * d1 + md * d2 + d2 * md 
                          + 2.0 * mH1sq * d2 + 2.0 * hd2) - 4.0 * Mlamsq * d2  
           - 4.0 * lam * hlam * (dt * hd + hdt * d1); 

	dme2 = -(2.0 * me + 8.0 * mH1sq) * e2t * e2t - 
	  4.0 * et * ml * e2 * e1 -
	  4.0 * e2t * me * e2t - 4.0 * e2t * et * ml * e1 - 
	  2.0 * e2t * e2t * me -
	  ((me + 4.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me) *
	  (6.0 * ddT + 2.0 * eeT) - 4.0 * e2t * 
	  (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2 + 
	   e1 * me * et).trace() -
	  4.0 * (he2t * e2t + e2t * he2t + het * e2 * he + et * he2 * e1) -
	  4.0 * he2t * (3.0 * ddT + eeT) - 4.0 * e2t * (3.0 * Ub2 + Utau2) -
	  4.0 * het * e1 * (3.0 * hddT + heeT) - 4.0 * et * he * 
	  (3.0 * hddT + heeT) + (6.0 * gsq(2) - 1.2 * gsq(1)) * 
	  ((me + 2.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me + 
	   2.0 * he2t) +
	  12.0 * gsq(2) * (2.0 * msq(2) * e2t - mG(2) * het * e1 
			   - mG(2) * et * he) -
	  2.4 * gsq(1) * (2.0 * msq(1) * e2t - mG(1) * het * e1 - 
			  mG(1) * et * he)
	  + 2.4 * gsq(1) * sP + 2808. / 25.0 * g4(1) * msq(1) + 
	  2.4 * gsq(1) * sigma(1) // checked
        - 2.0 * lsq * (2.0 * et * ml * e1 + me * e2 + e2 * me + 2.0 * mH1sq * e2
                       + 2.0 * he2) - 4.0 * Mlamsq * e2  
           - 4.0 * lam * hlam * (et * he + het * e1); 
	dhu2 = 
	  (-3.0 * (3.0 * ht4 + ht2 * hb2) - d2 * (3.0 * ddT + eeT)
	   -15.0 * u2 * uuT - 6.0 * u4 - 2.0 * d4 -
	   4.0 * u2 * d2 + (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT +
	   12.0 * gsq(2) * u2 + 0.4 * gsq(1) * d2 - 16.0 / 9.0 * g4(3) +
	   8.0 * gsq(3) * gsq(2) + 136.0 / 45.0 * gsq(3) * gsq(1) + 
	   7.5 * g4(2) + 
	   gsq(2) * gsq(1) + 2743.0 / 450.0 * g4(1)) * hu +
	  (-6.0 * (6.0 * Ut * ht2 * ht + Ut * hb2 * ht + Ub * ht2 * hb) -
	   18.0 * u2 * huuT - d2 * (6.0 * hddT + 2.0 * heeT) - 
	   12.0 * hu * ut * uuT - hd * dt * (6.0 * ddT + 2.0 * eeT) - 
	   6.0 * hu * u2t * ut - 8.0 * u2 * hu * ut - 4.0 * hd * d2t * dt -
	   4.0 * d2 * hd * dt - 2.0 * hu * ut * d2 - 4.0 * u2 * hd * dt +
	   (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT + 
	   (6.0 * gsq(2) + 1.2 * gsq(1)) * hu * ut + 0.8 * gsq(1) * hd * dt -
	   (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
	   (12.0 * gsq(2) * mG(2) + 0.8 * gsq(1) * mG(1)) * u2 -
	   0.8 * gsq(1) * mG(1) * d2 + 64.0 / 9.0 * g4(3) * mG(3) - 
	   16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) - 
	   272.0 / 45.0 * gsq(3) * gsq(1) * (mG(1) + mG(3)) - 
	   30.0 * g4(2) * mG(2) - 2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 
	   5486.0 / 225.0 * g4(1) * mG(1)) * u1 // checked
           - lsq * hu * (3.0 * lsq + 2.0 * ksq + 3.0 * ddT + eeT)
           - lsq * (5.0 * u2 * hu + 4.0 * hu * u2t + d2 * hu + 2.0 * hd * dt * u1)
           - 2.0 * lam * hlam * u1 * (3.0 * lsq + 2.0 * ksq + 3.0 * ddT + eeT)
           - 2.0 * lsq * u1 * (3.0 * lam * hlam + 2.0 * kap * hkap 
                               + 3.0 * hddT + heeT)
           - 2.0 * lam * hlam * (3.0 * u2 * u1 + d2 * u1); 
        
	dhd2 =  // checked + debugged 19.05.04 by M R Ramage (thanks Mike!)
	  (-3.0 * (3.0 * hb4 + ht2 * hb2 + htau4) -
	   3.0 * u2 * uuT - 5.0 * d2 * (3.0 * ddT + eeT) - 6.0 * d4 -
	   2.0 * u4 - 4.0 * d2 * u2 + (16.0 * gsq(3) - 
					    0.4 * gsq(1)) * ddT +
	   1.2 * gsq(1) * eeT + 0.8 * gsq(1) * u2 + 
	   (12.0 * gsq(2) + 1.2 * gsq(1)) * d2 - 16.0 / 9.0 * g4(3) + 
	   8.0 * gsq(3) * gsq(2) + 8.0 / 9.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
	   gsq(2) * gsq(1) + 287.0 / 90.0 * g4(1)) * hd +
	  (-6.0 * (6.0 * Ub * hb2 * hb + Ut * hb2 * ht + Ub * ht2 * hb +
		   2.0 * Utau * htau2 * htau) - 6.0 * u2 * huuT -
	   6.0 * d2 * (3.0 * hddT + heeT) - 6.0 * hu * ut * uuT -
	   4.0 * hd * dt * (3.0 * ddT + eeT) - 6.0 * hd * d2t * dt -
	   8.0 * d2 * hd * dt - 4.0 * u2 * hu * ut - 4.0 * hu * u2t * ut -
	   4.0 * d2 * hu * ut - 2.0 * hd * dt * u2 + 
	   (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT +
	   1.6 * gsq(1) * hu * ut + (6.0 * gsq(2) + 1.2 * gsq(1)) * hd * dt -
	   (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT - 
	   2.4 * gsq(1) * mG(1) * eeT - 
	   (12.0 * gsq(2) * mG(2) + 1.6 * gsq(1) * mG(1)) * d2 -
	   1.6 * gsq(1) * mG(1) * u2 + 64.0 / 9.0 * g4(3) * mG(3) - 
	   16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) - 
	   16.0 / 9.0 * gsq(3) * gsq(1) * (mG(3) + mG(1)) - 
	   30.0 * g4(2) * mG(2) - 
	   2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 
	   574.0 / 45.0 * g4(1) * mG(1)) * d1
           - lsq * hd * (3.0 * lsq + 2.0 * ksq + 3.0 * uuT)
           - lsq * (5.0 * d2 * hd + 4.0 * hd * d2t + u2 * hd + 2.0 * hu * ut * d1)
           - 2.0 * lam * hlam * d1 * (3.0 * lsq + 2.0 * ksq + 3.0 * uuT)
           - 2.0 * lsq * d1 * (3.0 * lam * hlam + 2.0 * kap * hkap + 3.0 * huuT)
           - 2.0 * lam * hlam * (3.0 * d2 * d1 + u2 * d1); 

	dhe2 = 
	  (-3.0 * (3.0 * hb4 + ht2 * hb2 + htau4) -
	   5.0 * e2 * (3.0 * ddT + eeT) - 6.0 * e2 * e2 + 
	   (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT +
	   (12.0 * gsq(2) - 1.2 * gsq(1)) * eeT + 7.5 * g4(2) + 
	   1.8 * gsq(2) * gsq(1) + 13.5 * g4(1)) * he +
	  (-6.0 * (6.0 * Ub * hb2 * hb + Ut * hb2 * ht + Ub * ht2 * hb +
		   2.0 * Utau * htau2 * htau) - 
	   4.0 * he * et * (3.0 * ddT + eeT) - 6.0 * e2 * (3.0 * hddT + heeT) -
	   6.0 * he * e2t * et - 8.0 * e2 * he * et + 
	   (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT + 
	   (6.0 * gsq(2) + 1.2 * gsq(1)) * he * et -
	   (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT - 
	   2.4 * gsq(1) * mG(1) * eeT - 12.0 * gsq(2) * mG(2) * e2 -
	   30.0 * g4(2) * mG(2) - 3.6 * gsq(2) * gsq(1) * (mG(1) + mG(2)) -
	   54.0 * g4(1) * mG(1)) * e1 // checked
           - lsq * he * (3.0 * lsq + 2.0 * ksq + 3.0 * uuT)
           - lsq * (3.0 * e2 * he + 3.0 * he * e2t)
           - 2.0 * lam * hlam * e1 * (3.0 * lsq + 2.0 * ksq + 3.0 * uuT)
           - 2.0 * lsq * e1 * (3.0 * lam * hlam + 2.0 * kap * hkap + 3.0 * huuT)
           - 6.0 * lam * hlam *  e2 * e1; 
      

 //Mixing < 1  // PA: warning in  3rd family terms the trilinears are Ut etc and the Yukawas are ht, hb, htau 
        

  dhlam2 = - 50.0 * l4 * hlam - 36.0 * lam * (Ut * ht * ht2)
           - 36.0 * lam * (Ub * hb * hb2)  
           - 12.0 * lam * (Utau * htau * htau2) 
           - 9.0 * hlam * (ht2 * ht2)
           - 9.0 * hlam * (hb2 * hb2)
           - 3.0 * hlam * (htau2 * htau2)
           - 8.0 * k4 * hlam - 32.0 * lam * kap * ksq * hkap 
           - 12.0 * lsq * ksq * hlam - 6.0 * lsq *(hlam * Ytr + lam * aYtr) 
           - 24.0 * lsq * kap * (kap * hlam  +  lam * hkap) 
           - 12.0 * lam * (Ut * ht * hb2 + Ub * hb * ht2)
           - 3.0 * lsq * hlam * Ytr - 6.0 * hlam * (ht2 * hb2)
           + 2.4 * gsq(1) * lsq * (1.5 * hlam - mG(1)) 
           + 1.6 * gsq(1) * lam * (huuT - mG(1) * uuT) 
           - 0.8 * gsq(1) * lam * (hddT - mG(1) * ddT)
           + 2.4 * gsq(1) * lam * (heeT - mG(1) * eeT)
           + 0.4 * gsq(1) * hlam * (2.0 * uuT - ddT + 3.0 * eeT)
           + 12.0 * gsq(2) * lsq * hlam * (1.5 * hlam - mG(2))
           + 32.0 * gsq(3) * lam * (huuT - mG(3) * uuT) 
           + 32.0 * gsq(3) * lam * (hddT - mG(3) * ddT)
           + 16.0 * gsq(3) * hlam * (uuT + ddT) 
           + 0.02 * g4(1) * (207.0 * hlam - 828.0 * lam * mG(1))
           + 0.5 * g4(2) * (15.0 * hlam - 60.0 * lam * mG(2))
           + 1.8 * gsq(1) * gsq(2) * (hlam - 2.0* (mG(1) + mG(2)));
 

} else {
	// NB you should introduce speedy computation here. For example sum
	// traces, not add then trace. Also surely some of the products are
	// repeated, eg u2 * u2 = u4, d2 * d2 = d4, u1 * mu * ut, d1 * md *
	// dt, d2 * d1, d2 * mq
         dmq2 = -(2.0 * mq + 8.0 * mH2sq) * u4 - 4.0 * u1 * mu * u2t * ut 
            - 4.0 * u2 * mq * u2 - 4.0 * u2 * u1 * mu * ut - 2.0 * u4 * mq -
            (2.0 * mq + 8.0 * mH1sq) * d4 - 4.0 * d1 * md * d2t * dt - 
            4.0 * d2 * mq * d2 - 4.0 * d2 * d1 * md * dt - 2.0 * d4 * mq -
            ((mq + 4.0 * mH2sq) * u2 + 2.0 * u1 * mu * ut + u2 * mq) * uuT * 3.0 
            -((mq + 4.0 * mH1sq) * d2 + 2.0 * d1 * md * dt + d2 * mq) * 
            (3.0 * ddT + eeT) -
            6.0 * u2 * ((mq * u2).trace() + (u1 * mu * ut).trace()) -
            d2 * (6.0 * (mq * d2).trace() + 6.0 * (d1 * md * dt).trace() + 
                  2.0 * (ml * e2).trace() + 2.0 * (e1 * me * et).trace()) -
            4.0 * (u2 * hu2 + hu2 * u2 + u1 * hu2t * ut + hu * u2t * hut) -
            4.0 * (d2 * hd2 + hd2 * d2 + d1 * hd2t * dt + hd * d2t * hdt) -
            hu2 * uuT * 6.0 - u2 * 6.0 * huT - 
            hu * ut * 6.0 * huuT - u1 * hut * 6.0 * huuT - hd2 * 
            (6.0 * ddT + 2.0 * eeT) - 
            d2 * (6.0 * hdT + 2.0 * heT) -
            hd * dt * (6.0 * hddT + 2.0 * heeT) - 
            d1 * hdt * (6.0 * hddT + 2.0 * heeT) +
            0.4 * gsq(1) * 
            ((2.0 * mq +  4.0 * mH2sq) * u2 + 4.0 * u1 * mu * ut +
             2.0 * u2 * mq + 4.0 * hu2 - 4.0 * mG(1) * hu * ut - 4.0
             * mG(1) * u1 * hut + 8.0 * msq(1) * u2 +
             (mq + 2.0 * mH1sq) * d2 + 2.0 * d1 * md * dt +
             d2 * mq + 2.0 * hd2 - 2.0 * mG(1) * hd * dt - 
             2.0 * mG(1) * d1 * hdt + 4.0 * msq(1) * d2) +
            0.4 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 32.0 * gsq(3) *
            gsq(2) * (msq(3) + msq(2) + mG(2) * mG(3)) +
            32.0 / 45.0 * gsq(3) * gsq(1) * (msq(3) + msq(1) + 
					   mG(1) * mG(3)) + 33.0
            * g4(2) * msq(2) + 0.4 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) *
						      mG(2)) +
            199.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
            3.0 * gsq(2) * sigma(2) + gsq(1) * sigma(1) / 15.0 // checked
            - lsq * (2.0 * ut * mu * u1 + mq * u2 + u2 * mq + 2.0 * mH2sq * u2
                     + 2.0 * hu2) - 2.0 * Mlamsq * u2  
           - 2.0 * lam * hlam * (u1 * hut + hu * ut)
            - lsq * (2.0 * dt * md * d1 + mq * d2 + d2 * mq + 2.0 * mH1sq * d2
                     + 2.0 * hd2) - 2.0 * Mlamsq * d2  
            - 2.0 * lam * hlam * (d1 * hdt + hd * dt);
        

         //PA: checked 14/9/2012
         
         
         dml2 = -(2.0 * ml + 8.0 * mH1sq) * e2 * e2 - 4.0 * e1 * me * e2t * et -
            4.0 * e2 * ml * e2 - 4.0 * e2 * e1 * me * et - 2.0 * e2 * e2 * ml -
            ((ml + 4.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml) * 
            (3.0 * ddT + eeT) -
            e2 * (6.0 * mq * d2 + 6.0 * d1 * md * dt + 2.0 * ml * e2 + 
                  2.0 * e1 * me * et).trace() -
            4.0 * (e2 * he2 + he2 * e2 + e1 * he2t * et + he * e2t * het) -
            he2 * (6.0 * ddT + 2.0 * eeT) - 
            e2 * (6.0 * hdT + 2.0 * heT) -
            he * et * (6.0 * hddT + 2.0 * heeT) - 
            e1 * het * (6.0 * hddT + 2.0 * heeT) +
            1.2 * gsq(1) * ((ml + 2.0 * mH1sq) * e2 + 2.0 * e1 * me * et + 
                            e2 * ml +
                            2.0 * he2 - 2.0 * mG(1) * he * et - 
                            2.0 * mG(1) * e1 *
                            het + 4.0 * msq(1) * e2) -
            1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 3.6 * gsq(2) * gsq(1) *
            (msq(2) + msq(1) + mG(1) * mG(2)) + 621.0 / 25.0 * g4(1) * msq(1) +
            3.0 * gsq(2) * sigma(2) + 0.6 * gsq(1) * sigma(1) //checked
            - lsq * (2.0 * et * me * e1 + ml * e2 + e2 * ml + 2.0 * mH1sq * e2
                     + 2.0 * he2) - 2.0 * Mlamsq * e2  
            - 2.0 * lam * hlam * (e1 * het + he * et);	
         
         //PA: checked numerically 18/9/2012
         
         dmu2 = -(2.0 * mu + 8.0 * mH2sq) * u4t - 
            4.0 * ut * mq * u2 * u1 -
            4.0 * u2t * mu * u2t - 4.0 * u2t * ut * mq * u1 - 
            2.0 * u4t * mu - (2.0 * mu + 4.0 * mH2sq + 
                              4.0 * mH1sq) * ut * d2
            * u1 - 4.0 * ut * mq * d2 * u1 - 4.0 * ut * d1 * md * dt * u1 -
            4.0 * ut * d2 * mq * u1 - 2.0 * ut * d2 * u1 * mu -
            ((mu + 4.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + 
             u2t * mu) * 6.0 * uuT -
            12.0 * u2t * (mq * u2 + u1 * mu * ut).trace() -
            4.0 * (hu2t * u2t + u2t * hu2t + hut * u2 * hu + ut * hu2 * u1) -
            4.0 * (hut * hd * dt * u1 + ut * d1 * hdt * hu + hut * d2 * hu + 
                   ut * hd2 * u1) -
            12.0 * (hu2t * uuT + u2t * huT + hut * u1 * huuT +
                    ut * hu * huuT) +
            (6.0 * gsq(2) - 0.4 * gsq(1)) * 
            ((mu + 2.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu + 
             2.0 * hu2t) +
            12.0 * gsq(2) * (2.0 * msq(2) * u2t - mG(2) * hut * u1 - 
                             mG(2) * ut * hu) -
            0.8 * gsq(1) * (2.0 * msq(1) * u2t - mG(1) * hut * u1 - 
                            mG(1) * ut * hu)-
            1.6 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 
            512.0 / 45.0 * gsq(1) * gsq(3) * (msq(3) + msq(1) + mG(1) * mG(3)) +
            3424.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
            16.0 / 15.0 * gsq(1) * sigma(1) // checked
            - 2.0 * lsq * (2.0 * ut * mq * u1 + mu * u2t + u2t * mu 
                           + 2.0 * mH2sq * u2t + 2.0 * hu2t) - 4.0 * Mlamsq * u2t  
            - 4.0 * lam * hlam * (ut * hu + hut * u1);
         //PA: checked and debugged 18/9/2012 
        
         dmd2 = -(2.0 * md + 8.0 * mH1sq) * d4t - 
            4.0 * dt * mq * d1 * d2t -
            4.0 * d2t * md * d2t - 4.0 * d2t * dt * mq * d1 - 
            2.0 * d4t * md -
            (2.0 * md + 4.0 * mH2sq + 4.0 * mH1sq) * dt * u2 * d1 -
            4.0 * dt * mq * u2 * d1 - 4.0 * dt * u1 * mu * ut * d1 - 
            4.0 * dt * u2 * mq * d1 - 2.0 * dt * u2 * d1 * md - 
            ((md + 4.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md) * 
            (6.0 * ddT + 2.0 * eeT) -
            4.0 * d2t * (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2 
                         + e1 * me * et).trace() -
            4.0 * (hd2t * d2t + d2t * hd2t + hdt * d2 * hd + dt * hd2 * d1) -
            4.0 * (hdt * hu * ut * d1 + dt * u1 * hut * hd + hdt * u2 * hd + 
                   dt * hu2 * d1) -
            4.0 * hd2t * (3.0 * ddT + eeT) - 4.0 * d2t * (3.0 * hdT + heT) -
            4.0 * hdt * d1 * (3.0 * hddT + heeT) - 
            4.0 * dt * hd * (3.0 * hddT + heeT) +
            (6.0 * gsq(2) + 0.4 * gsq(1)) * 
            ((md + 2.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md + 
             2.0 * hd2t) +
            12.0 * gsq(2) * (2.0 * msq(2) * d2t - mG(2) * hdt * d1 
                             - mG(2) * dt * hd) +
            0.8 * gsq(1) * (2.0 * msq(1) * d2t - mG(1) * hdt * d1 - 
                            mG(1) * dt * hd)+
            0.8 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 
            128.0 / 45.0 * gsq(3) * gsq(1) * (msq(3) + msq(1) + mG(1) * mG(3)) +
            808.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
            4.0 / 15.0 * gsq(1) * sigma(1) // checked
            - 2.0 * lsq * (2.0 * dt * mq * d1 + md * d2t + d2t * md 
                           + 2.0 * mH1sq * d2t + 2.0 * hd2t) - 4.0 * Mlamsq * d2t  
            - 4.0 * lam * hlam * (dt * hd + hdt * d1); 
         //PA: checked numerically 18/9/2012
         
         dme2 = -(2.0 * me + 8.0 * mH1sq) * e2t * e2t - 
            4.0 * et * ml * e2 * e1 -
            4.0 * e2t * me * e2t - 4.0 * e2t * et * ml * e1 - 
            2.0 * e2t * e2t * me -
            ((me + 4.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me) *
            (6.0 * ddT + 2.0 * eeT) - 4.0 * e2t * 
            (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2 + 
             e1 * me * et).trace() -
            4.0 * (he2t * e2t + e2t * he2t + het * e2 * he + et * he2 * e1) -
            4.0 * he2t * (3.0 * ddT + eeT) - 4.0 * e2t * (3.0 * hdT + heT) -
            4.0 * het * e1 * (3.0 * hddT + heeT) - 4.0 * et * he * 
            (3.0 * hddT + heeT) + (6.0 * gsq(2) - 1.2 * gsq(1)) * 
            ((me + 2.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me + 
             2.0 * he2t) +
            12.0 * gsq(2) * (2.0 * msq(2) * e2t - mG(2) * het * e1 
                             - mG(2) * et * he) -
            2.4 * gsq(1) * (2.0 * msq(1) * e2t - mG(1) * het * e1 - 
                            mG(1) * et * he)
            + 2.4 * gsq(1) * sP + 2808. / 25.0 * g4(1) * msq(1) + 
            2.4 * gsq(1) * sigma(1) //checked
            - 2.0 * lsq * (2.0 * et * ml * e1 + me * e2t + e2t * me 
                           + 2.0 * mH1sq * e2t + 2.0 * he2t) - 4.0 * Mlamsq * e2t  
            - 4.0 * lam * hlam * (et * he + het * e1);
        //PA: checked 18/9/2012
 
    /* Old full 3-family version 
    dhd2 =  // checked + debugged 19.05.04 by M R Ramage (thanks Mike!)
      (-3.0 * (3.0 * d2t * d2t + ut * d2 * u1 + e2t * e2t).trace() -
       3.0 * u2 * uuT - 5.0 * d2 * (3.0 * ddT + eeT) - 6.0 * d2 * d2 -
       2.0 * u2 * u2 - 4.0 * d2 * u2 + (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT +
       1.2 * gsq(1) * eeT + 0.8 * gsq(1) * u2 + 
       (12.0 * gsq(2) + 1.2 * gsq(1)) * d2 - 16.0 / 9.0 * g4(3) + 
       8.0 * gsq(3) * gsq(2) + 8.0 / 9.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
       gsq(2) * gsq(1) + 287.0 / 90.0 * g4(1)) * hd +
      (-6.0 * (6.0 * hdt * d2 * d1 + hut * d2 * u1 + hdt * u2 * d1 +
	       2.0 * het * e2 * e1).trace() - 6.0 * u2 * huuT -
       6.0 * d2 * (3.0 * hddT + heeT) - 6.0 * hu * ut * uuT -
       4.0 * hd * dt * (3.0 * ddT + eeT) - 6.0 * hd * d2t * dt -
       8.0 * d2 * hd * dt - 4.0 * u2 * hu * ut - 4.0 * hu * u2t * ut -
       4.0 * d2 * hu * ut - 2.0 * hd * dt * u2 + 
       (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT +
       1.6 * gsq(1) * hu * ut + (6.0 * gsq(2) + 1.2 * gsq(1)) * hd * dt -
       (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT - 
       2.4 * gsq(1) * mG(1) * eeT - 
       (12.0 * gsq(2) * mG(2) + 1.6 * gsq(1) * mG(1)) * d2 -
       1.6 * gsq(1) * mG(1) * u2 + 64.0 / 9.0 * g4(3) * mG(3) - 
       16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) - 
       16.0 / 9.0 * gsq(3) * gsq(1) * (mG(3) + mG(1)) - 30.0 * g4(2) * mG(2) - 
       2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 574.0 / 45.0 * g4(1) * mG(1)) 
       * d1; 
    */      

	dhu2 = 
           (-3.0 * (3.0 * u4 + ut * d2 * u1).trace() - 
            d2 * (3.0 * ddT + eeT)
            -15.0 * u2 * uuT - 6.0 * u4 - 2.0 * d4 -
            4.0 * u2 * d2 + (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT +
            12.0 * gsq(2) * u2 + 0.4 * gsq(1) * d2 - 16.0 / 9.0 * g4(3) +
            8.0 * gsq(3) * gsq(2) + 136.0 / 45.0 * gsq(3) * gsq(1) + 
            7.5 * g4(2) + 
            gsq(2) * gsq(1) + 2743.0 / 450.0 * g4(1)) * hu +
           (-6.0 * (6.0 * hut * u2 * u1 + hut * d2 * u1 + 
                    hdt * u2 * d1).trace() -
            18.0 * u2 * huuT - d2 * (6.0 * hddT + 2.0 * heeT) - 
            12.0 * hu * ut * uuT - hd * dt * (6.0 * ddT + 2.0 * eeT) - 
            6.0 * hu * u2t * ut - 8.0 * u2 * hu * ut - 4.0 * hd * d2t * dt -
            4.0 * d2 * hd * dt - 2.0 * hu * ut * d2 - 4.0 * u2 * hd * dt +
            (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT + 
            (6.0 * gsq(2) + 1.2 * gsq(1)) * hu * ut + 0.8 * gsq(1) * hd * dt -
            (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
            (12.0 * gsq(2) * mG(2) + 0.8 * gsq(1) * mG(1)) * u2 -
            0.8 * gsq(1) * mG(1) * d2 + 64.0 / 9.0 * g4(3) * mG(3) - 
            16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) - 
            272.0 / 45.0 * gsq(3) * gsq(1) * (mG(1) + mG(3)) - 
            30.0 * g4(2) * mG(2) - 2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 
            5486.0 / 225.0 * g4(1) * mG(1)) * u1 // checked
           - lsq * hu * (3.0 * lsq + 2.0 * ksq + 3.0 * ddT + eeT)
           - lsq * (5.0 * u2 * hu + 4.0 * hu * u2t + d2 * hu + 2.0 * hd * dt * u1)
           - 2.0 * lam * hlam * u1 * (3.0 * lsq + 2.0 * ksq + 3.0 * ddT + eeT)
           - 2.0 * lsq * u1 * (3.0 * lam * hlam + 2.0 * kap * hkap 
                               + 3.0 * hddT + heeT)
           - 2.0 * lam * hlam * (3.0 * u2 * u1 + d2 * u1); 

	dhd2 =  // checked + debugged 19.05.04 by M R Ramage (thanks Mike!)
           (-3.0 * (3.0 * d4t + ut * d2 * u1 + e2t * e2t).trace() -
            3.0 * u2 * uuT - 5.0 * d2 * (3.0 * ddT + eeT) - 6.0 * d4 -
            2.0 * u4 - 4.0 * d2 * u2 + (16.0 * gsq(3) 
                                        - 0.4 * gsq(1)) * ddT +
            1.2 * gsq(1) * eeT + 0.8 * gsq(1) * u2 + 
            (12.0 * gsq(2) + 1.2 * gsq(1)) * d2 - 16.0 / 9.0 * g4(3) + 
            8.0 * gsq(3) * gsq(2) + 8.0 / 9.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
            gsq(2) * gsq(1) + 287.0 / 90.0 * g4(1)) * hd +
           (-6.0 * (6.0 * hdt * d2 * d1 + hut * d2 * u1 + hdt * u2 * d1 +
                    2.0 * het * e2 * e1).trace() - 6.0 * u2 * huuT -
            6.0 * d2 * (3.0 * hddT + heeT) - 6.0 * hu * ut * uuT -
            4.0 * hd * dt * (3.0 * ddT + eeT) - 6.0 * hd * d2t * dt -
            8.0 * d2 * hd * dt - 4.0 * u2 * hu * ut - 4.0 * hu * u2t * ut -
            4.0 * d2 * hu * ut - 2.0 * hd * dt * u2 + 
            (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT +
            1.6 * gsq(1) * hu * ut + (6.0 * gsq(2) + 1.2 * gsq(1)) * hd * dt -
            (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT - 
            2.4 * gsq(1) * mG(1) * eeT - 
            (12.0 * gsq(2) * mG(2) + 1.6 * gsq(1) * mG(1)) * d2 -
            1.6 * gsq(1) * mG(1) * u2 + 64.0 / 9.0 * g4(3) * mG(3) - 
            16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) - 
            16.0 / 9.0 * gsq(3) * gsq(1) * (mG(3) + mG(1)) - 
            30.0 * g4(2) * mG(2) - 
            2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 
            574.0 / 45.0 * g4(1) * mG(1)) 
           * d1
           - lsq * hd * (3.0 * lsq + 2.0 * ksq + 3.0 * uuT)
           - lsq * (5.0 * d2 * hd + 4.0 * hd * d2t + u2 * hd + 2.0 * hu * ut * d1)
           - 2.0 * lam * hlam * d1 * (3.0 * lsq + 2.0 * ksq + 3.0 * uuT)
           - 2.0 * lsq * d1 * (3.0 * lam * hlam + 2.0 * kap * hkap + 3.0 * huuT)
           - 2.0 * lam * hlam * (3.0 * d2 * d1 + u2 * d1); 

	dhe2 = 
           (-3.0 * (3.0 * d4t + ut * d2 * u1 + e2t * e2t).trace() -
            5.0 * e2 * (3.0 * ddT + eeT) - 6.0 * e2 * e2 + 
            (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT +
            (12.0 * gsq(2) - 1.2 * gsq(1)) * e2 + 7.5 * g4(2) + 
            1.8 * gsq(2) * gsq(1) + 13.5 * g4(1)) * he +
           (-6.0 * (6.0 * hdt * d2 * d1 + hut * d2 * u1 + hdt * u2 * d1 +
                    2.0 * het * e2 * e1).trace() - 
            4.0 * he * et * (3.0 * ddT + eeT) - 6.0 * e2 * (3.0 * hddT + heeT) -
            6.0 * he * e2t * et - 8.0 * e2 * he * et + 
            (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT + 
            (6.0 * gsq(2) + 1.2 * gsq(1)) * he * et -
            (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT - 
            2.4 * gsq(1) * mG(1) * eeT - 12.0 * gsq(2) * mG(2) * e2 -
            30.0 * g4(2) * mG(2) - 3.6 * gsq(2) * gsq(1) * (mG(1) + mG(2)) -
            54.0 * g4(1) * mG(1)) * e1;
	
  // checked
       dhe2 = dhe2 - lsq * he * (3.0 * lsq + 2.0 * ksq + 3.0 * uuT)
           - lsq * (5.0 * e2 * he + 4.0 * he * e2t)
           - 2.0 * lam * hlam * e1 * (3.0 * lsq + 2.0 * ksq + 3.0 * uuT)
           - 2.0 * lsq * e1 * (3.0 * lam * hlam + 2.0 * kap * hkap + 3.0 * huuT)
           - 6.0 * lam * hlam *  e2 * e1; 

	
        dhlam2 = - 50.0 * l4 * hlam - 36.0 * lam * (hu * ut * u2).trace()
           - 36.0 * lam * (hd * dt * d2).trace()  
           - 12.0 * lam * (he * et * e2).trace() 
          - hlam * (3.0 * (e4).trace() + 9.0 * (u4).trace() + 9.0 * (d4).trace())
           - 8.0 * k4 * hlam - 32.0 * lam * kap * ksq * hkap 
           - 12.0 * lsq * ksq * hlam - 6.0 * lsq * (hlam * Ytr + lam * aYtr) 
           - 24.0 * lsq * kap * (kap * hlam  +  lam * hkap) 
           - 12.0 * lam * (hu * ut * d2 + hd * dt * u2).trace()
           - 3.0 * lsq * hlam * Ytr - 6.0 * hlam * (u2 * d2).trace()
           + 2.4 * gsq(1) * lsq * (1.5 * hlam - lam * mG(1)) 
	  + 1.6 * gsq(1) * lam * (huuT - mG(1) * uuT) 
           - 0.8 * gsq(1) * lam * (hddT - mG(1) * ddT)
           + 2.4 * gsq(1) * lam * (heeT - mG(1) * eeT)
           + 0.4 * gsq(1) * hlam * (2.0 * uuT - ddT + 3.0 * eeT)
           + 12.0 * gsq(2) * lsq * (1.5 * hlam - lam * mG(2))
           + 32.0 * gsq(3) * lam * (huuT - mG(3) * uuT) 
           + 32.0 * gsq(3) * lam * (hddT - mG(3) * ddT)
           + 16.0 * gsq(3) * hlam * (uuT + ddT) 
           + 0.02 * g4(1) * (207.0 * hlam - 828.0 * lam * mG(1))
           + 0.5 * g4(2) * (15.0 * hlam - 60.0 * lam * mG(2))
           + 1.8 * gsq(1) * gsq(2) * (hlam - 2.0 * lam * (mG(1) + mG(2)));
        

      }
      dhkap2 = - 120.0 * k4 * hkap - 12.0 * l4 * hkap 
         - 48.0 * lsq * lam * kap * hlam
         - 48.0 * ksq * kap * lam * hlam - 48.0 * lsq * ksq * hkap
         - 24.0 * lsq * ksq * hkap - 12.0 * lam * kap * (hlam * Ytr 
                                                         +  lam * aYtr) 
         - 6.0 * lsq * hkap * Ytr 
         + 7.2 * gsq(1) * lam * (kap * hlam  + 0.5 * lam * hkap 
                                 - lam * kap * mG(1)) 
         + 36.0 * gsq(2) * lam * (kap * hlam  +  0.5 * lam * hkap 
                                  - lam * kap * mG(2)); 
      
    // add onto one-loop beta functions
    dmq2 *= oneO16Pif; dmq += dmq2;
    dml2 *= oneO16Pif; dml += dml2;
    dmu2 *= oneO16Pif; dmu += dmu2;
    dmd2 *= oneO16Pif; dmd += dmd2;
    dme2 *= oneO16Pif; dme += dme2;
    dhu2 *= oneO16Pif; dhu += dhu2;
    dhd2 *= oneO16Pif; dhd += dhd2;
    dhe2 *= oneO16Pif; dhe += dhe2;
    dhlam2 *=oneO16Pif; dhlam += dhlam2; 
    dhkap2 *=oneO16Pif; dhkap += dhkap2;

    }
    
    // Default is to include these 2-loop corrections anyhow because they can 
    // be so important: gauginos + higgs 
    dmG2 = 2.0 * gsq * 
       (babBeta * gsqM + mG * (babBeta * gsq) + 
        cuBeta * huuT - cuBeta * mG * uuT + 
        cdBeta * hddT - cdBeta * mG * ddT + 
        ceBeta * heeT - ceBeta * mG * eeT +
        clBeta * lam * hlam - clBeta * mG * lsq); // checked

    double dmH1sq2, dm3sq2, dmH2sq2;
    double dmSsq2, dmSpsq2, dz_s2;
    if (MIXING < 1) {
      /// The following are valid in the third-family approximation -- they are
      /// much faster, and should be a good approximation to the no-mixed case
      dm3sq2 = m3sq * 
         (-3.0 * (3.0 * ht4 + 3.0 * hb4 + 2.0 * hb2 * ht2 + htau4) +
          (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT + 
          (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 
          1.2 * gsq(1) * eeT + 7.5 * g4(2) + 1.8 * gsq(1) * gsq(2) + 
          207. / 50. * g4(1)) +
         displaySusyMu() * 
         (-12.0 * (3.0 * Ut * ht * ht2 + 3.0 * Ub * hb * hb2 + Ut * hb2 * ht +
                   Ub * ht2 * hb + Utau * htau2 * htau) +
          (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT + 
          (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT -
          (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
          (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
          2.4 * gsq(1) * mG(1) * eeT - 30.0 * g4(2) * mG(2) - 
          18.0 / 5.0 * gsq(1)
          * gsq(2) * (mG(1) + mG(2)) - 414. / 25.0 * g4(1) * mG(1)) // checked
         - 14.0 * l4 * m3sq - 32.0 * smu * lam * lsq * hlam
         - 15.0 * lsq * m3sq * uuT - 18.0 * lsq * smu * huuT 
         - 6.0 * smu * lam * hlam * uuT
         - 15.0 * lsq * m3sq * ddT - 18.0 * lsq * smu * hddT
         - 6.0 * smu * lam * hlam * ddT
         - 5.0 * lsq * m3sq * eeT - 6.0 * lsq * smu * heeT
         - 2.0 * smu * lam * hlam * eeT
         - 4.0 * lsq * ksq * m3sq - 8.0 * lam * ksq * smu * hlam
         - 8.0 * lsq * kap * smu * hkap - 8.0 * lam * kap * (ksq + lsq) * mSpsq
         - 8.0 * (lsq * kap * hlam + ksq * lam * hkap) * mu_s  
         + 2.4 * gsq(1) * lsq * (m3sq - smu * mG(1))
         + 12 * gsq(2) * lsq * (m3sq - smu * mG(2));

      //PA: checked 18/9/2012

      dmH2sq2 = -6.0 * 
         (6.0 * (mH2sq + mq3 + mu3) * ht4 + 
          hb2 * ht2 * (mH1sq + mH2sq + 2.0 * mq3 + mu3 + md3) + 
          ht2 * 12.0 * Ut2 + ht2 * Ub2 + hb2 * Ut2 + 2.0 * hb * ht * Ut * Ub) +
         (32.0 * gsq(3) + 1.6 * gsq(1)) * 
         ((mH2sq + mq3 + mu3) * ht2 + Ut2) +
         32.0 * gsq(3) * (2.0 * msq(3) * ht2 - 2.0 * mG(3) * ht * Ut) +
         1.6 * gsq(1) * (2.0 * msq(1) * ht2 - 2.0 * mG(1) * ht * Ut) +
         1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 
         18.0 / 5.0 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
         621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) +
	0.6 * gsq(1) * sigma(1); // checked
         - 12.0 * lsq * (Mlamsq + hlam * hlam) 
         - 6.0 * (lsq * Xd + ddT * Mlamsq + 2.0 * lam * hlam * hddT)
         - 2.0 * (lsq * Xe + eeT * Mlamsq + 2.0 * lam * hlam * heeT)
         - 4.0 * (ksq * Mlamsq + lsq * Mkapsq + 2.0 * lam * kap * hlam * hkap);
      
      
      dmH1sq2 = -6.0 *
         (6.0 * (mH1sq + mq3 + md3) * hb4 + 
          (mH1sq + mH2sq + 2.0 * mq3 + mu3 + md3) * ht2 * hb2 + 
          2.0 * (mH1sq + ml3 + me3) * htau4 +
          12.0 * Ub2 * hb2 + hb2 * Ut2 + ht2 * Ub2 + 2.0 * Ut * ht * Ub * hb +
          4.0 * htau2 * Utau2) +
         (32.0 * gsq(3) - 0.8 * gsq(1)) * ((mH1sq + mq3 + md3) * hb2 + Ub2) +
         32.0 * gsq(3) * (2.0 * msq(3) * ddT - 2.0 * mG(3) * hddT) -
         0.8 * gsq(1) * (2.0 * msq(1) * ddT - 2.0 * mG(1) * hddT) +
         2.4 * gsq(1) * (((mH1sq + ml3 + me3) * htau2 + Utau2) +
                         2.0 * msq(1) * eeT - 2.0 * mG(1) * heeT) 
         - 1.2 * gsq(1) * sP +
         33.0 * g4(2) * msq(2) + 
         3.6 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
         621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) + 
         0.6 * gsq(1) * sigma(1) // checked
         - 12.0 * lsq * (Mlamsq + hlam * hlam) 
         - 6.0 * (lsq * Xu + uuT * Mlamsq + 2.0 * lam * hlam * huuT)
         - 4.0 * (ksq * Mlamsq + lsq * Mkapsq + 2.0 * lam * kap * hlam * hkap);
      
      dmSsq2 = - 16.0 * lsq * (Mlamsq + sqr(hlam)) 
         - 32.0 * ksq * (Mkapsq + sqr(hkap))
         - 12.0 * (Mlamsq * uuT + lsq * Xu + 2.0 * lam * hlam * huuT)
         - 12.0 * (Mlamsq * ddT + lsq * Xd + 2.0 * lam * hlam * hddT)
         - 4.0 * (Mlamsq * eeT + lsq * Xe + 2.0 * lam * hlam * heeT)
         - 16.0 * (ksq * Mlamsq + lsq * Mkapsq + 2.0 * lam * kap * hlam * hkap)
         + 2.4 * gsq(1) * (Mlamsq - 2.0 * lam * hlam * mG(1) + 2.0 * lsq * msq(1)) 
         + 12.0 * gsq(2) * (Mlamsq - 2.0 * lam * hlam * mG(2) + 2.0 * lsq * msq(2));
      
      dmSpsq2 = - 8.0 * (l4 * mSpsq + 4.0 * lsq * lam * mu_s * hlam)
         - 16.0 * (2.0 * k4 * mSpsq + 5.0 * ksq * kap * mu_s * hkap)
         - 32.0 * ksq * lsq * mSpsq - 48.0 * lam * ksq * mu_s * hlam
         - 32.0 * lsq * kap * mu_s * hkap - 8.0 * lam * hlam * mu_s * Ytr  
         - 4.0 * lsq * ( mSpsq * Ytr + 2.0 * mu_s * aYtr)
         - 16.0 * lsq * kap * (lam * m3sq + smu * hlam)
         - 8.0 * lam * kap * (m3sq * Ytr + smu * aYtr)
         + 4.8 * lam * kap * gsq(1) * (m3sq - smu * mG(1))
         + 24.0 * lam * kap * gsq(2) * (m3sq - smu * mG(2))
         + 2.4 * gsq(1) * (lsq * mSpsq + 2.0 * mu_s * lam * hlam 
                           - 2.0 * mu_s * mG(1))
         + 12.0 * gsq(2) * (lsq * mSpsq + 2.0 * mu_s * lam * hlam 
                           - 2.0 * mu_s * mG(2));
         

   

         
                        
                        
                        
      
    } else {
      /// In the mixed case, we need to use the slower full 3-family
      /// expressions
      dmH2sq2 = -6.0 * 
	(6.0 * (mH2sq + mq) * u4 + 6.0 * u1 * mu * u2t * ut +
	 (mH1sq + mH2sq + mq) * u2 * d2 + u1 * mu * ut * d2 +
	 u2 * mq * d2 + u2 * d1 * md * dt + 6.0 * hu2 * u2 +
	 6.0 * hu * u2t * hut + hd2 * u2 + d2 * hu2 + 
	 hd * dt * u1 * hut + d1 * hdt * hu * ut).trace() +
	(32.0 * gsq(3) + 1.6 * gsq(1)) * 
	((mH2sq + mq) * u2 + u1 * mu * ut + hu2).trace() +

	32.0 * gsq(3) * (2.0 * msq(3) * uuT - 2.0 * mG(3) * huuT) +
	1.6 * gsq(1) * (2.0 * msq(1) * uuT - 2.0 * mG(1) * huuT) +
	1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 
	18.0 / 5.0 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
	621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) +
	0.6 * gsq(1) * sigma(1) // checked
         - 12.0 * lsq * (Mlamsq + hlam * hlam) 
         - 6.0 * (lsq * Xd + ddT * Mlamsq + 2.0 * lam * hlam * hddT)
         - 2.0 * (lsq * Xe + eeT * Mlamsq + 2.0 * lam * hlam * heeT)
         - 4.0 * (ksq * Mlamsq + lsq * Mkapsq + 2.0 * lam * kap * hlam * hkap); 
	 

     //PA: checked numerically 19/9/2012

      dmH1sq2 = -6.0 *
	(6.0 * (mH1sq + mq) * d4 + 6.0 * dt * md * d2 * d1 +
	 (mH1sq + mH2sq + mq) * u2 * d2 + u1 * mu * ut * d2 +
	 u2 * mq * d2 + u2 * d1 * md * dt + 2.0 * (mH1sq + ml) * e2 * e2 +
	 2.0 * e1 * me * e2t * et + 6.0 * hd2 * d2 + 6.0 * hd * d2t * hdt + 
	 hu2 * d2 + u2 * hd2 + hu * ut * d1 * hdt + u1 * hut * hd * dt +
	 2.0 * he2 * e2 + 2.0 * he * e2t * het).trace() +
	(32.0 * gsq(3) - 0.8 * gsq(1)) * ((mH1sq + mq) * d2 + d1 * md * dt +
					  hd2).trace() +
	32.0 * gsq(3) * (2.0 * msq(3) * ddT - 2.0 * mG(3) * hddT) -
	0.8 * gsq(1) * (2.0 * msq(1) * ddT - 2.0 * mG(1) * hddT) +
	2.4 * gsq(1) * (((mH1sq + ml) * e2 + e1 * me * et + he2).trace() +
			2.0 * msq(1) * eeT - 2.0 * mG(1) * heeT) - 
	1.2 * gsq(1) * sP +
	33.0 * g4(2) * msq(2) + 
	3.6 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
	621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) + 
	0.6 * gsq(1) *
	sigma(1) // checked+corrected 27/9/05
          - 12.0 * lsq * (Mlamsq + hlam * hlam) 
         - 6.0 * (lsq * Xu + uuT * Mlamsq + 2.0 * lam * hlam * huuT)
         - 4.0 * (ksq * Mlamsq + lsq * Mkapsq + 2.0 * lam * kap * hlam * hkap);
      
      //PA: checked numerically 19/9/2012

      dm3sq2 = m3sq * 
	(-3.0 * (3.0 * u4t + 3.0 * d4t + 2.0 * ut * d2
		 * u1 + e2t * e2t).trace() +
	 (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT + 
	 (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 
	 1.2 * gsq(1) * eeT + 7.5 * g4(2) + 1.8 * gsq(1) * gsq(2) + 
	 207. / 50. * g4(1)) +
	displaySusyMu() * 
	(-12.0 * (3.0 * hut * u2 * u1 + 3.0 * hdt * d2 * d1 + hut * d2 * u1 +
		  hdt * u2 * d1 + het * e2 * e1).trace() +
	 (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT + 
	 (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT -
	 (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
	 (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
	 2.4 * gsq(1) * mG(1) * eeT - 30.0 * g4(2) * mG(2) 
	 - 18.0 / 5.0 * gsq(1)
	 * gsq(2) * (mG(1) + mG(2)) - 414. / 25.0 * g4(1) * mG(1)); // checked
      
      //PA: new nmssm parts  
      dm3sq2 = dm3sq2 - 14.0 * l4 * m3sq - 32.0 * smu * lam * lsq * hlam
         - 15.0 * lsq * m3sq * uuT - 18.0 * lsq * smu * huuT 
         - 6.0 * smu * lam * hlam * uuT
         - 15.0 * lsq * m3sq * ddT - 18.0 * lsq * smu * hddT
         - 6.0 * smu * lam * hlam * ddT
         - 5.0 * lsq * m3sq * eeT - 6.0 * lsq * smu * heeT
         - 2.0 * smu * lam * hlam * eeT
         - 4.0 * lsq * ksq * m3sq - 8.0 * lam * ksq * smu * hlam
         - 8.0 * lsq * kap * smu * hkap - 8.0 * lam * kap * (ksq + lsq) * mSpsq
         - 8.0 * (lsq * kap * hlam + ksq * lam * hkap) * mu_s  
         + 2.4 * gsq(1) * lsq * (m3sq - smu * mG(1))
         + 12 * gsq(2) * lsq * (m3sq - smu * mG(2));
  

      dmSsq2 = - 16.0 * lsq * (Mlamsq + sqr(hlam)) 
         - 32.0 * ksq * (Mkapsq + sqr(hkap))
         - 12.0 * (Mlamsq * uuT + lsq * Xu + 2.0 * lam * hlam * huuT)
         - 12.0 * (Mlamsq * ddT + lsq * Xd + 2.0 * lam * hlam * hddT)
         - 4.0 * (Mlamsq * eeT + lsq * Xe + 2.0 * lam * hlam * heeT)
         - 16.0 * (ksq * Mlamsq + lsq * Mkapsq + 2.0 * lam * kap * hlam * hkap)
         + 2.4 * gsq(1) * (Mlamsq - 2.0 * lam * hlam * mG(1) 
                           + 2.0 * lsq * msq(1)) 
         + 12.0 * gsq(2) * (Mlamsq - 2.0 * lam * hlam * mG(2) 
                            + 2.0 * lsq * msq(2));
 
//PA: checked numerically 19/9/2012

    
      dmSpsq2 = - 8.0 * (l4 * mSpsq + 4.0 * lsq * lam * mu_s * hlam)
         - 16.0 * (2.0 * k4 * mSpsq + 5.0 * ksq * kap * mu_s * hkap)
         - 32.0 * ksq * lsq * mSpsq - 48.0 * lam * ksq * mu_s * hlam
         - 32.0 * lsq * kap * mu_s * hkap
         - 4.0 * lsq * (mSpsq * Ytr + 2.0 * mu_s * aYtr) 
         - 8.0 * lam * hlam * mu_s * Ytr
         - 16.0 * lsq * kap * (lam * m3sq + smu * hlam)
         - 8.0 * lam * kap * (m3sq * Ytr + smu * aYtr)
         + 4.8 * lam * kap * gsq(1) * (m3sq - smu * mG(1))
         + 24.0 * lam * kap * gsq(2) * (m3sq - smu * mG(2))
         + 2.4 * gsq(1) * (lsq * mSpsq + 2.0 * mu_s * lam * hlam 
                           - 2.0 * lsq * mu_s * mG(1))
         + 12.0 * gsq(2) * (lsq * mSpsq + 2.0 * mu_s * lam * hlam 
                           - 2.0 * lsq * mu_s * mG(2));

      //PA checked 19/9/2012
    }
      dz_s2 = - 4.0 * l4 * zeta_s - 16.0 * lam * lsq * hlam * zeta 
         - 8.0 * k4 * zeta_s - 32.0 * kap * ksq * hkap * zeta 
         - 2.0 * lsq * zeta * (zeta_s * Ytr + 2.0 * aYtr) 
	- 4.0 * lam * hlam * zeta * Ytr
	- 8.0 * lsq * ksq * zeta_s 
	- 16.0 * lam * kap * (kap * hlam + lam * hkap) * zeta
         - 4.0 * lam * (m3sq * (mu_s * Ytr + aYtr)
                        + smu * (3.0 * (Xu + Xd) + Xe) 
                        + smu * (mH1sq + mH2sq) * Ytr 
                        + smu * mu_s * aYtr)
         - 4.0 * hlam * (m3sq * Ytr + smu * aYtr)
         - 8.0 * lsq * (m3sq * (2.0 * hlam + lam * mu_s) + hlam * smu * mu_s) 
	- 8.0 * lam * (Mlamsq + sqr(hlam) + lsq * (mH1sq + mH2sq)) * smu
         - 8.0 * lsq * kap * mu_s * (mSpsq + 2.0 * mSsq) 
         - 8.0 * (kap * mu_s * Mlamsq + (kap * lam * hlam + lsq * hkap) * mSpsq
                  + lam * kap * hlam * sqr(mu_s) + lam * hlam * hkap * mu_s )
          - 8.0 * ksq * (mSpsq * (2.0 * hkap + kap * mu_s) + hkap * mu_s * mu_s)
	- 8.0 * kap * (mu_s * Mkapsq + sqr(hkap) * mu_s + 2.0 * ksq * mu_s * mSsq)
         + 1.2 * gsq(1) * (3.0 * m3sq - 2.0 * smu * mG(1)) * hlam
         + 1.2 * gsq(1) * lam * (3.0 * m3sq * (mu_s - mG(1)) + lam * zeta_s
                                + 2.0 * smu * (mH1sq + mH2sq - mu_s * mG(1) 
                                               + 2.0 * msq(1))
                                + 2.0 * zeta * (hlam - lam * mG(1)))
         + 3.0 * gsq(2) * (3.0 * m3sq - 2.0 * smu * mG(2)) * hlam
         + 3.0 * gsq(2) * lam * (3.0 * m3sq * (mu_s - mG(2)) + lam * zeta_s
                                + 2.0 * smu * (mH1sq + mH2sq - mu_s * mG(2) 
                                               + 2.0 * msq(2))
                                + 2.0 * zeta * (hlam - lam * mG(2)));

      
    dmG = dmG + dmG2 * oneO16Pif;
    dm3sq = dm3sq + dm3sq2 * oneO16Pif;
    dmH1sq = dmH1sq + dmH1sq2 * oneO16Pif;
    dmH2sq = dmH2sq + dmH2sq2 * oneO16Pif;
    dmSsq = dmSsq + dmSsq2 * oneO16Pif;
    dmSpsq = dmSpsq + dmSpsq2 * oneO16Pif;
    dz_s = dz_s + dz_s2 * oneO16Pif;
   
}  
  

  SoftParsNmssm dsoft(dsb, dmG, dhu, dhd, dhe, dhlam, dhkap, dmq, dmu,
                      dmd, dml, dme, dm3sq, dmH1sq, dmH2sq, dmSsq, dmSpsq, dz_s,
		     displayGravitino(), displayMu(),
		     displayLoops(), displayThresholds());
  return dsoft;
}

// Outputs derivatives vector y[109] for SUSY parameters: interfaces to
// integration routines
DoubleVector SoftParsNmssm::beta() const {
  // calculate the derivatives
  static SoftParsNmssm dsoft; dsoft = beta2();

  return dsoft.display(); // convert to a long vector
}

// Outputs derivatives of anomalous dimensions, from which the running can be
// derived. 
void SoftParsNmssm::anomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
				  DoubleMatrix & gQQ, DoubleMatrix & gUU,
				  DoubleMatrix & gDD, 
                                   double & gH1H1, double & gH2H2, double & gSS)  const {
  // Constants for gauge running
   static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3), clBeta(3);
  static DoubleMatrix babBeta(3, 3);
  
  if (bBeta(1) < 1.0e-5) // Constants not set yet
     nmsetBetas(babBeta, cuBeta, cdBeta, ceBeta, clBeta, bBeta);
  
  // For calculational brevity: 
  static nmsBrevity a;
  // convert to beta functions
  static nMssmSusy dsb;
  
  // calculate derivatives for full SUSY spectrum. Brevity calculations come
  // out encoded in a
  dsb = nMssmSusy::beta(a);
  
  static DoubleVector g1(3);
  g1 = displayGauge();
  
  // To keep this a const function: TIME SAVINGS
  DoubleMatrix &dt=a.dt, &ut=a.ut, &et=a.et;      
  DoubleVector &gsq=a.gsq; // &g3=a.g3, &g4=a.g4;
  
  static DoubleMatrix hu(3, 3), hd(3, 3), he(3, 3);
  static DoubleMatrix hut(3, 3), hdt(3, 3), het(3, 3);
  static DoubleVector gsqM(1, 3);
  
  hu = displayTrilinear(UA); hd = displayTrilinear(DA); 
  he = displayTrilinear(EA); 
  hut = hu.transpose(); hdt = hd.transpose(); het = he.transpose();
  
  double huuT = (hu * ut).trace(), hddT = (hd * dt).trace(), 
    heeT = (he * et).trace(); 
  gsqM = gsq * displayGaugino(); 
  
  if (displayLoops() > 0) { // CHECKED: agrees with our conventions
    const static double eightO3 = 8.0 / 3.0, oneO30 = 1.0 / 30.0;
    gLL = - (he * et + 0.3 * gsqM(1) + 1.5 * gsqM(2));
    gEE = - (2.0 * (et * he) + 1.2 * gsqM(1));
    gQQ = - (hd * dt + hu * ut + oneO30 * gsqM(1) + 1.5 * gsqM(2) + 
	     eightO3 * gsqM(3));
    gDD = - (2.0 * dt * hd + 4.0 * oneO30 * gsqM(1) + 
	     eightO3 * gsqM(3));
    gUU = - (2.0 * ut * hu + 16.0 * oneO30 * gsqM(1) +
	     eightO3 * gsqM(3));
    gH1H1 = - (3.0 * hddT + heeT + 0.3 * gsqM(1) + 1.5 * gsqM(2));
    gH2H2 = - (3.0 * huuT + 0.3 * gsqM(1) + 1.5 * gsqM(2));

    const static double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
    gH1H1 = gH1H1 * oneO16Pisq;
    gH2H2 = gH2H2 * oneO16Pisq;
    gEE  *= oneO16Pisq;    
    gLL  *= oneO16Pisq;
    gQQ  *= oneO16Pisq;
    gUU  *= oneO16Pisq;
    gDD  *= oneO16Pisq;
  }
}

// Reads in universal boundary conditions at the current scale:
// m0, M1/2, A0, B and sign of mu
void SoftParsNmssm::universal(double m0,  double m12,  double a0,  double mu,
			      double m3sq) {
  standardSugra(m0, m12, a0);  
  setSusyMu(mu);
  setM3Squared(m3sq);
}

//PA: for fully constrained models
void SoftParsNmssm::universalScalars(double m0) {
  // scalar masses
  DoubleMatrix ID(3, 3), mm0(3, 3);
  int i; for (i=1; i<=3; i++) ID(i, i) = 1.0;
  mm0 = ID * sqr(m0);
  setSoftMassMatrix(mQl, mm0); setSoftMassMatrix(mUr, mm0);
  setSoftMassMatrix(mDr, mm0); setSoftMassMatrix(mLl, mm0);
  setSoftMassMatrix(mEr, mm0);
  setMh1Squared(sqr(m0)); setMh2Squared(sqr(m0));
  setMsSquared(sqr(m0));
}
//PA: for semi constrained models
void SoftParsNmssm::semiuniversalScalars(double m0, double mS) {
  // scalar masses
  DoubleMatrix ID(3, 3), mm0(3, 3);
  int i; for (i=1; i<=3; i++) ID(i, i) = 1.0;
  mm0 = ID * sqr(m0);
  setSoftMassMatrix(mQl, mm0); setSoftMassMatrix(mUr, mm0);
  setSoftMassMatrix(mDr, mm0); setSoftMassMatrix(mLl, mm0);
  setSoftMassMatrix(mEr, mm0);
  setMh1Squared(sqr(m0)); setMh2Squared(sqr(m0));
  setMsSquared(sqr(mS));
}

//PA:: for fully constrained
void SoftParsNmssm::universalTrilinears(double a0)  {  
  // trilinears
  setTrilinearMatrix(UA, a0 * displayYukawaMatrix(YU)); 
  setTrilinearMatrix(DA, a0 * displayYukawaMatrix(YD));
  setTrilinearMatrix(EA, a0 * displayYukawaMatrix(YE));
  setTrialambda(a0 * displayLambda()); 
  setTriakappa(a0 * displayKappa());
}
 //PA:: for semi constrained models 
//(allows both Alambda and Akappa to be separately specified)
void SoftParsNmssm::semiuniversalTrilinears(double a0, double al, double ak)  {  
  // trilinears
  setTrilinearMatrix(UA, a0 * displayYukawaMatrix(YU)); 
  setTrilinearMatrix(DA, a0 * displayYukawaMatrix(YD));
  setTrilinearMatrix(EA, a0 * displayYukawaMatrix(YE));
   setTrialambda(al * displayLambda()); 
   setTriakappa(ak * displayKappa());
}
 


// Input m0, NOT m0 squared.
void SoftParsNmssm::standardSugra(double m0,  double m12, double a0) {
  /*  if (m0 < 0.0) {
    ostringstream ii;
    ii << "m0=" << m0 << " passed to universal boundary" <<
      "conditions illegally negative.";
    throw ii.str();
    } Deleted on request from A Pukhov */
  universalScalars(m0);
  universalGauginos(m12);
  universalTrilinears(a0);
}
//PA: for semi constrained models
void SoftParsNmssm::standardsemiSugra(double m0,  double mS, double m12, double a0, double Al, double Ak) {
  /*  if (m0 < 0.0) {
    ostringstream ii;
    ii << "m0=" << m0 << " passed to universal boundary" <<
      "conditions illegally negative.";
    throw ii.str();
    } Deleted on request from A Pukhov */
   semiuniversalScalars(m0, mS);
  universalGauginos(m12);
  semiuniversalTrilinears(a0, Al, Ak);
}


#define HR "---------------------------------------------------------------\n"

ostream & operator <<(ostream &left, const SoftParsNmssm &s) {
  left << "SUSY breaking MSSM parameters at Q: " << s.displayMu() << endl;
  left << " UA" << s.displayTrilinear(UA) 
       << " UD" << s.displayTrilinear(DA) 
       << " UE" << s.displayTrilinear(EA); 
  left << " alambda: " << s.displayTrialambda()
       << " akappa: " << s.displayTriakappa();
  left << " mQLsq" << s.displaySoftMassSquared(mQl) 
       << " mURsq" << s.displaySoftMassSquared(mUr) 
       << " mDRsq" << s.displaySoftMassSquared(mDr) 
       << " mLLsq" << s.displaySoftMassSquared(mLl) 
       << " mSEsq" << s.displaySoftMassSquared(mEr);
  left << "m3sq: " << s.displayM3Squared() 
       << " mH1sq: " << s.displayMh1Squared() 
       << " mH2sq: " << s.displayMh2Squared() 
       << " mSsq: "  << s.displayMsSquared()
       << " mSPsq: "  << s.displayMspSquared()
      << " zeta_s: "  << s.displayZeta_s()
       << '\n';
  left << "Gaugino masses" << s.displayGaugino();
  left << s.displaySusy();
  return left;
}

#undef HR

void SoftParsNmssm::inputSoftParsOnly() {
  SoftPars<nMssmSusy, nmsBrevity>::inputSoftParsOnly();
  char c[70];
  cin >> c >> mSsq >> c >> mSpsq >> c >> zeta_s;
  cin >> c >> alambda
      >> c >> akappa;
}

istream & operator >>(istream &left, SoftParsNmssm &s) {
  char c[70];

  left >> c >> c >> c >> c >> c >> c >> c;
  DoubleMatrix ua(3, 3), da(3, 3), ea(3, 3);
  left >> c >> ua
       >> c >> da
       >> c >> ea;
  s.setTrilinearMatrix(UA, ua);
  s.setTrilinearMatrix(DA, da);
  s.setTrilinearMatrix(EA, ea);
  double al, ak;
  left >> c >> al
       >> c >> ak;
  s.setTrialambda(al); 
  s.setTriakappa(ak);
  DoubleMatrix mqlsq(3, 3), mursq(3, 3), mdrsq(3, 3), mllsq(3, 3), mersq(3, 3);
  left >> c >> mqlsq
       >> c >> mursq
       >> c >> mdrsq
       >> c >> mllsq
       >> c >> mersq;
  s.setSoftMassMatrix(mQl, mqlsq); 
  s.setSoftMassMatrix(mUr, mursq); 
  s.setSoftMassMatrix(mDr, mdrsq); 
  s.setSoftMassMatrix(mLl, mllsq); 
  s.setSoftMassMatrix(mEr, mersq); 
  double m3sq, mh1sq, mh2sq;
  left >> c >> m3sq >> c >> mh1sq >> c >> mh2sq;
  s.setM3Squared(m3sq); s.setMh1Squared(mh1sq); s.setMh2Squared(mh2sq);
  double mSsq, mSpsq, zeta_s;
  left >> c >> mSsq >> c >> mSpsq >> c >> zeta_s;
  s.setMsSquared(mSsq); s.setMspSquared(mSpsq); s.setZeta_s(zeta_s);
  DoubleVector mg(3);
  left >> c >> mg; 
  s.setAllGauginos(mg);
  nMssmSusy ss;
  left >> ss;   s.setSusy(ss);
  return left;
}

// Boundary conditions to be applied at messenger scale for Gauge mediated
// SUSY breaking (see hep-ph/9703211 for example)
void SoftParsNmssm::minimalGmsb(int n5, double lambda, double mMess, 
			       double cgrav) {

// Modified thresholds by JEL 1-26-04 to accomodate numerical infinities

  const double epstol = 1.0e-4;
  double x = lambda / mMess;

  double f, g;

  if(fabs(x) < epstol) {
    g = 1.0 + x*x/6.0 + sqr(x*x)/15.0;
    f = 1.0 + x*x/36.0 + 11.0*sqr(x*x)/450.0;
  }
  else if(fabs(x-1.0) < 0.0001) {
    g  =  log(4.0);
    f  = -sqr(PI)/6.0 + log(4.0) + 0.5*sqr(log(4.0));
    g -=  0.0008132638905771205626;
    f -= -0.0049563838821509165200;
  }
  else {
    g = 1.0 / sqr(x) * 
      ((1.0 + x) * log(1.0 + x) + (1.0 - x) * log(1.0 - x));
    f = (1.0 + x) / sqr(x) * 
    (log(1.0 + x) - 2.0 * dilog(x / (1.0 + x)) + 0.5 * 
     dilog(2.0 * x / (1.0 + x))) + 
     (1.0 - x) / sqr(x) * (log(1.0 - x) - 2.0 * dilog(-x / (1.0 - x)) +
			 0.5 * dilog(-2.0 * x / (1.0 - x)));
  }

  double n5d = double(n5);

  double m1, m2, m3;
  m1 = n5d * sqr(displayGaugeCoupling(1)) / (16.0 * sqr(PI)) * lambda * g; 
  m2 = n5d * sqr(displayGaugeCoupling(2)) / (16.0 * sqr(PI)) * lambda * g; 
  m3 = n5d * sqr(displayGaugeCoupling(3)) / (16.0 * sqr(PI)) * lambda * g; 
  setGauginoMass(1, m1);   setGauginoMass(2, m2);   setGauginoMass(3, m3);

  setM32(2.37e-19 * lambda * mMess * cgrav);

  double g1f = sqr(sqr(displayGaugeCoupling(1)));
  double g2f = sqr(sqr(displayGaugeCoupling(2)));
  double g3f = sqr(sqr(displayGaugeCoupling(3)));

  double mursq, mdrsq, mersq, mqlsq, mllsq;
  mursq = 2.0 * f * sqr(lambda) * n5d * 
    (4.0 / 3.0 * g3f + 0.6 * 4.0 / 9.0 * g1f) 
    / sqr(16.0 * sqr(PI));
  mdrsq = 2.0 * f * sqr(lambda) * n5d * 
    (4.0 / 3.0 * g3f + 0.6 * 1.0 / 9.0 * g1f) 
    / sqr(16.0 * sqr(PI));
  mersq = 2.0 * f * sqr(lambda) * n5d * 
    (0.6 * g1f) 
    / sqr(16.0 * sqr(PI));
  mqlsq = 2.0 * f * sqr(lambda) * n5d * 
    (4.0 / 3.0 * g3f + 0.75 * g2f + 0.6 * g1f / 36.0) 
    / sqr(16.0 * sqr(PI));
  mllsq = 2.0 * f * sqr(lambda) * n5d * 
    (                  0.75 * g2f + 0.6 * 0.25 * g1f) 
    / sqr(16.0 * sqr(PI));

  // You need Higgs masses too!

  DoubleMatrix id(3, 3);
  id(1, 1) = 1.0; id(2, 2) = 1.0; id(3, 3) = 1.0;

  setSoftMassMatrix(mQl, mqlsq * id);
  setSoftMassMatrix(mUr, mursq * id);
  setSoftMassMatrix(mDr, mdrsq * id);
  setSoftMassMatrix(mLl, mllsq * id);  
  setMh1Squared(mllsq);
  setMh2Squared(mllsq);
  setSoftMassMatrix(mEr, mersq * id);

  universalTrilinears(0.0);
}

