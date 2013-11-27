
/** \file rgeroutines.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: odds and sods for my own machievellian needs

*/

#include <src/rgeroutines.h>

namespace softsusy {

double mlqspc(const MssmSoftsusy & r) {
  double chi = fabs(r.displayPhys().mneut.display(1)), 
    xi = fabs(r.displayPhys().mneut.display(2)), 
    l = (r.displayPhys().me.display(2, 1)),
    q  = ((r.displayPhys().mu.display(1, 1) +
	      r.displayPhys().md.display(1, 1)) * 0.5);
  return mlqspc(q, xi, l, chi);
}

double mlqspc(double msqL, double mchi20, double mer, double mchi10) {
  double chi = sqr(mchi10), 
    xi = sqr(mchi20), 
    l = sqr(mer),
    q  = sqr(msqL);

  return sqrt(fabs((q - xi) * (l - chi) / (2.0 * l - chi)));
}

double mlqMaxFar(const MssmSoftsusy & r) {

  double mchi20 = r.displayPhys().mneut.display(2), 
    mchi10 = r.displayPhys().mneut.display(1),
    msqL = 0.5 * (r.displayPhys().mu.display(1, 1) +
		  r.displayPhys().md.display(1, 1)), 
    mer = r.displayPhys().me.display(2, 1);

  return mlqMaxFar(msqL, mchi20, mer, mchi10);
}

double mlqMaxFar(double msqL, double mchi20, double mer, double mchi10) {

  double sqrtArg = (sqr(msqL) - sqr(mchi20)) * (sqr(mer) - sqr(mchi10)) 
    / (sqr(mer));

  double answer = sqrtArg;

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double mllqMax(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut.display(1), 
    mchi20 = r.displayPhys().mneut.display(2), 
    msqL = 0.5 * (r.displayPhys().mu.display(1, 1) + 
		  r.displayPhys().md.display(1, 1)), 
    mer = r.displayPhys().me.display(2, 1);

  double q = sqr(msqL), xi = sqr(mchi20), chi = sqr(mchi10), sl = sqr(mer);
  double m1 = (q - xi) * (xi - chi) / xi;
  double m2 = (q - sl) * (sl - chi) / sl;
  double m3 = (q * sl - xi * chi) * (xi - sl) / (xi * sl);

  if (sqr(sl) < q * chi && q * chi < sqr(xi) && sqr(xi) * chi < q * sqr(sl)) {

    return (msqL - mchi10);
  }

  double answer = maximum(m3, maximum(m1, m2));

  if (answer > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double mllqMax(double msqL, double mchi20, double mer, double mchi10) {

  double q = sqr(msqL), xi = sqr(mchi20), chi = sqr(mchi10), sl = sqr(mer);
  double m1 = (q - xi) * (xi - chi) / xi;
  double m2 = (q - sl) * (sl - chi) / sl;
  double m3 = (q * sl - xi * chi) * (xi - sl) / (xi * sl);

  if (sqr(sl) < q * chi && q * chi < sqr(xi) && sqr(xi) * chi < q * sqr(sl)) {
    return (msqL - mchi10);
  }

  double answer = maximum(m3, maximum(m1, m2));

  if (answer > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double mllbMin(const MssmSoftsusy & r) {

  double mchi10 = fabs(r.displayPhys().mneut.display(1)), 
    mchi20 = fabs(r.displayPhys().mneut.display(2)), 
    msqL = minimum(r.displayPhys().md.display(1, 3), r.displayPhys().md.display(2, 3)), 
    mer = r.displayPhys().me.display(2, 1); 
 
  double sqrtArg = (sqr(mchi10) * sqr(mchi10) + sqr(mer) * sqr(mer)) * 
	     sqr(sqr(mchi20) + sqr(mer)) + 2.0 * sqr(mchi10) * sqr(mer) *
	     (sqr(mchi20) * sqr(mchi20) - 6.0 * sqr(mchi20) * sqr(mer) + 
	      sqr(mer) * sqr(mer));

  double answer = 
    ( - sqr(mchi10) * sqr(mchi20) * sqr(mchi20) + 
      3.0 * sqr(mchi10) * sqr(mchi20) * sqr(mer) - 
      sqr(mchi20) * sqr(mchi20) * sqr(mer) - 
      sqr(mchi20) * sqr(mer) * sqr(mer) - 
      sqr(mchi10) * sqr(mchi20) * sqr(msqL) - 
      sqr(mchi10) * sqr(mer) * sqr(msqL) + 
      3.0 * sqr(mchi20) * sqr(mer) * sqr(msqL) - 
      sqr(mer) * sqr(mer) * sqr(msqL) + (sqr(mchi20) - sqr(msqL)) 
      * sqrt(fabs(sqrtArg)))
    / (4.0 * sqr(mchi20) * sqr(mer));

  if (sqrtArg > 0.0 && answer > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double chiSqSU3(const MssmSoftsusy & r, double & chiMllMax, 
		double & chiMllqMax, double & chiMllqMin, 
		double & chiMlqMin, double & chiMlqMax, 
		double & chiMtautauMax, double msl, double msq) { 

  // exptl data taken from ATLAS TDR: 1 fb^-1
  /*  double mllMaxExp = 99.7, mllMaxErr = 2.3;
  double mllqMaxExp = 517., mllqMaxErr = 37.5;
  double mllqMinExp = 265., mllqMinErr = 25.0;
  double mlqMinExp = 333., mlqMinErr = 15.0;
  double mlqMaxExp = 445., mlqMaxErr = 23.6;
  double mtautauMaxExp = 102., mtautauMaxErr = 19.2;*/

  /// Checking increase by factor of 3 in precision: what happens?
  double mllMaxExp = 99.7, mllMaxErr = 2.3 / 3.;
  double mllqMaxExp = 517., mllqMaxErr = 37.5 / 3.;
  double mllqMinExp = 265., mllqMinErr = 25.0 / 3.;
  double mlqMinExp = 333., mlqMinErr = 15.0 / 3.;
  double mlqMaxExp = 445., mlqMaxErr = 23.6 / 3.;
  double mtautauMaxExp = 102., mtautauMaxErr = 19.2 / 3.;

  // need absolute masses for kinematics
  double mchi10 = fabs(r.displayPhys().mneut.display(1));
  double mchi20 = fabs(r.displayPhys().mneut.display(2));

  // these need special treatment
  double mlqmax = maximum(mlqMax(msq, mchi20, msl, mchi10), 
			  mlqMaxFar(msq, mchi20, msl, mchi10));
  double mlqmin = minimum(mlqMax(msq, mchi20, msl, mchi10), 
			  mlqspc(msq, mchi20, msl, mchi10));

  // Warning - in general, you should check to see if they are available!
  chiMllMax = sqr((mllMaxExp - mllMax(msq, mchi20, msl, mchi10)) / mllMaxErr);
  chiMllqMax = sqr((mllqMaxExp - mllqMax(msq, mchi20, msl, mchi10)) 
		   / mllqMaxErr);
  chiMllqMin = sqr((mllqMinExp - mllqMin(msq, mchi20, msl, mchi10)) 
		   / mllqMinErr);
  chiMlqMin = sqr((mlqMinExp - mlqmin) / mlqMinErr);
  chiMlqMax = sqr((mlqMaxExp - mlqmax) / mlqMaxErr);
  chiMtautauMax = sqr((mtautauMaxExp - mtautauMax(r)) / mtautauMaxErr);

  /// DEBUG
  /*
  cout << "mll=" << mllMax(r) << " mllqmax=" << mllqMax(r) << " mllqmin" << mllqMin(r) << " mlqmin=" << mlqmin << " mlqmax=" << mlqmax << endl;
  
  cout << r << "chisq: mllmax=" << chiMllMax << " mllqMax=" << chiMllqMax 
       << " mllqmin=" << chiMllqMin << " mlqlow=" << chiMlqMin 
       << " mlqmax=" << chiMlqMax << endl;
  */

  double chisq = chiMllMax + chiMllqMax + chiMllqMin + chiMlqMin + 
    chiMlqMax + chiMtautauMax;

  return chisq;
}

double mttMax(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut.display(1), 
    mchi20 = r.displayPhys().mneut.display(2), 
    mstau1 = minimum(r.displayPhys().me.display(1, 3),  
		     r.displayPhys().me.display(2, 3));

  double sqArg1 = sqr(mchi20) - sqr(mstau1);
  double sqArg2 = 1.0 - sqr(mchi10) / sqr(mstau1);

  if (sqArg1 > 0.0 && sqArg2 > 0.0) return sqrt(sqArg1 * sqArg2);
  else 
    return -sqrt(fabs(sqArg1 * sqArg2));
}

// The following edge parameters are defined in Paige, Bachacou, Hinchcliffe
// PRD 62 015009 (2000). They return minus values if the chain doesn't exist
double mhqMax(const MssmSoftsusy & r) {

  double mh0 = r.displayPhys().mh0(1), mchi10 = r.displayPhys().mneut.display(1), 
    mchi20 = r.displayPhys().mneut.display(2), 
    msqL = 0.5 * (r.displayPhys().mu.display(1, 1) + r.displayPhys().md.display(1, 1));

  double sqrtArg = sqr(sqr(mchi20) - sqr(mh0) - sqr(mchi10))
    - 4.0 * sqr(mh0) * sqr(mchi10);

  double answer = sqr(mh0) + (sqr(msqL) - sqr(mchi20)) * 
      (sqr(mchi20) + sqr(mh0) - sqr(mchi10) + 
       sqrt(fabs(sqrtArg))) / (2.0 * sqr(mchi20));

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(answer);
}

double mllMax(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut.display(1), 
    mchi20 = r.displayPhys().mneut.display(2), 
    mer = r.displayPhys().me.display(2, 1),
    msqL = 0.;
  
  return mllMax(msqL, mchi20, mer, mchi10);
}


double mllMax(double /* msqL */, double mchi20, double mer, double mchi10) {
  double sqrtArg = (sqr(mchi20) - sqr(mer)) *
    (sqr(mer) - sqr(mchi10));

  double answer = fabs(sqrtArg) /  sqr(mer);

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(answer);
}

double mtautauMax(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut.display(1), 
    mchi20 = r.displayPhys().mneut.display(2), 
    mer = minimum(r.displayPhys().me.display(2, 3), 
		  r.displayPhys().me.display(1, 3));

  double sqrtArg = (sqr(mchi20) - sqr(mer)) *
    (sqr(mer) - sqr(mchi10));

  double answer = fabs(sqrtArg) /  sqr(mer);

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(answer);
}

double edgeRatio(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut.display(1), mchi20 = r.displayPhys().mneut.display(2),
    mer = r.displayPhys().me.display(2, 1);

  double sqrtArg = (sqr(mchi20) - sqr(mer)) / (sqr(mchi20) - sqr(mchi10));

  double answer = fabs(sqrtArg);

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(answer);
}

double mlqMax(const MssmSoftsusy & r) {

  double mchi10 = fabs(r.displayPhys().mneut.display(1)), 
    mchi20 = fabs(r.displayPhys().mneut.display(2)), 
    msqL = 0.5 * (r.displayPhys().mu.display(1, 1) + 
		  r.displayPhys().md.display(1, 1)), 
    mer = r.displayPhys().me.display(2, 1);

  return mlqMax(msqL, mchi20, mer, mchi10);
}

double mlqMax(double msqL, double mchi20, double mer, double /* mchi10 */) {
  double sqrtArg = (sqr(msqL) - sqr(mchi20)) * (sqr(mchi20) - sqr(mer)) 
    / (sqr(mchi20));

  double answer = sqrtArg;

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double mllqMin(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut.display(1), 
    mchi20 = r.displayPhys().mneut.display(2), 
    msqL = 0.5 * (r.displayPhys().mu.display(1, 1) + 
		  r.displayPhys().md.display(1, 1)), 
    mer = r.displayPhys().me.display(2, 1);

  return mllqMin(msqL, mchi20, mer, mchi10);
}

double mllqMin(double msqL, double mchi20, double mer, double mchi10) {

  double sqrtArg = (sqr(mchi10) * sqr(mchi10) + sqr(mer) * sqr(mer)) * 
	     sqr(sqr(mchi20) + sqr(mer)) + 2.0 * sqr(mchi10) * sqr(mer) *
	     (sqr(mchi20) * sqr(mchi20) - 6.0 * sqr(mchi20) * sqr(mer) + 
	      sqr(mer) * sqr(mer));

  double answer = 
    ( - sqr(mchi10) * sqr(mchi20) * sqr(mchi20) + 
      3.0 * sqr(mchi10) * sqr(mchi20) * sqr(mer) - 
      sqr(mchi20) * sqr(mchi20) * sqr(mer) - 
      sqr(mchi20) * sqr(mer) * sqr(mer) - 
      sqr(mchi10) * sqr(mchi20) * sqr(msqL) - 
      sqr(mchi10) * sqr(mer) * sqr(msqL) + 
      3.0 * sqr(mchi20) * sqr(mer) * sqr(msqL) - 
      sqr(mer) * sqr(mer) * sqr(msqL) + (sqr(mchi20) - sqr(msqL)) 
      * sqrt(fabs(sqrtArg)))
    / (4.0 * sqr(mchi20) * sqr(mer));

  if (sqrtArg > 0.0 && answer > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double mhqMin(const MssmSoftsusy & r) {

  double mh0 = r.displayPhys().mh0(1), mchi10 = r.displayPhys().mneut.display(1), 
    mchi20 = r.displayPhys().mneut.display(2), 
    msqL = 0.5 * (r.displayPhys().mu.display(1, 1) + 
		  r.displayPhys().md.display(1, 1));

  double sqrtArg = sqr(sqr(mchi20) - sqr(mh0) - sqr(mchi10)) - 
    4.0 * sqr(mchi10) * sqr(mh0);

  double answer = 0.5 / sqr(mchi20) * 
    (sqr(msqL) - sqr(mchi20)) * (sqr(mchi20) + sqr(mh0) - sqr(mchi10) -
     sqrt(fabs(sqrtArg)));

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

} // namespace softsusy
