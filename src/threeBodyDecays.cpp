/** \file threeBodyDecays.cpp
   Project:     SOFTSUSY 
   Author:      Tom Cridge, Ben Allanach
   Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "threeBodyDecays.h"

/*static double m1 = 0.,m2 = 0.,m3 = 0.,m4 = 0.,mq = 0.,m5 = 0.,m6 = 0.,
  m7 = 0., m8 = 0., MZboson = 0., MWboson = 0., mh = 0., mH = 0.,
  mA = 0., mphi = 0., g1 = 0., g2 = 0., alphamix = 0., betavac = 0.;*/

/// First, do hadronic decays
double charginoToNeutralino1pion(const MssmSoftsusy * m) {
  double mchi1 = fabs(m->displayPhys().mch(1)),
    mneut1 = fabs(m->displayPhys().mneut(1));
  if (mchi1 < mneut1 + mpiplus) return 0.;
  if (mchi1 - mneut1 - mpiplus > hadronicScale) return 0.;

  double OL11 = cos(m->displayPhys().thetaL),
    OR11 = cos(m->displayPhys().thetaR); 
  double kpi = sqrt(lambda(sqr(mchi1), sqr(mneut1), sqr(mpiplus))) * 0.5
    / (mchi1);
  double width = sqr(fpi) * sqr(GMU) * kpi / (4.0 * PI * sqr(mchi1)) *
    ( sqr(OL11 + OR11) * ( sqr(sqr(mchi1) - sqr(mneut1)) -
			   sqr(mpiplus) * sqr(mchi1 - mneut1) ) +
      sqr(OL11 - OR11) * ( sqr(sqr(mchi1) - sqr(mneut1)) -
			   sqr(mpiplus) * sqr(mchi1 + mneut1))
      );

  return width; 
}

