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

double gluinoamplitudedecay (double m1, double m2, double m3, double alphastrong) {
  double squareratio, amplitudeW;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {
      squareratio = 1 + pow((m2/m1),2) - pow((m3/m1),2);
      amplitudeW = 1./4*alphastrong*(1/(2*m1))*squareratio*
	sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
  }
  return amplitudeW;
}


double gluinoamplitudedecaymix (double m1, double m2, double m3, double alphastrong, double squarkmix, double theta) {
  double squareratio, amplitudeW=0, squareratiomix1, squareratiomix2;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {
    squareratio = 1 + pow((m2/m1),2) - pow((m3/m1),2);
    if (squarkmix ==1) { 
      squareratiomix1= squareratio - 2*sin(2*theta)*m2/m1;
      amplitudeW = (alphastrong*1/4)*squareratiomix1*(1/(2*m1))*sqrt(lambda(sqr(m1), sqr(m2), sqr(m3))); 
    }
    else if (squarkmix ==2) {
      squareratiomix2 = squareratio + 2*sin(2*theta)*m2/m1;
      amplitudeW = (alphastrong*1/4)*squareratiomix2*(1/(2*m1))*sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
    }
    else {
      throw ("problem: squarkmix must be 1 or 2 in gluinoamplitudedecaymix\n");
    }
  }
  return amplitudeW;
}


double squarkamplitudedecaygluino (double m1, double m2, double m3,
				   double alphastrong) {
  double squareratio, amplitudeW;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {  
    squareratio = 1 - pow((m2/m1),2) - pow((m3/m1),2);
    amplitudeW = 4./3*alphastrong*(1/(2*m1))*squareratio*
      sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
  }
  return amplitudeW;
}

double squarkamplitudedecaygluinomix (double m1, double m2, double m3,
				      double alphastrong, double squarkmix,
				      double theta) {
  double squareratiomix, amplitudeW=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {
    if (squarkmix == 1) {
      squareratiomix = 1- pow(m2/m1,2) - pow(m3/m1,2) +2*sin(2*theta)*m2*m3/(pow(m1,2));
      amplitudeW = 4./3*alphastrong*1/(2*m1)*squareratiomix*
	sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
    }
    else if (squarkmix == 2) {
      squareratiomix = 1- pow(m2/m1,2) - pow(m3/m1,2) -2*sin(2*theta)*m2*m3/(pow(m1,2));
      amplitudeW = 4./3*alphastrong*1/(2*m1)*squareratiomix*
	sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
    }
    else {
      throw ("problem: quarkmix must be 1 or 2 in squarkamplitudedecaygluinomix\n");
      errorflag = -1;
    }
  }
  return amplitudeW;
}

