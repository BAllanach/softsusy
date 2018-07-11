/** \file twoBodyDecays.cpp
   Project:     SOFTSUSY 
   Author:      Tom Cridge, Ben Allanach
   Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "twoBodyDecays.h"

const double GFosqrt2 = GMU / root2;

/// First, do hadronic decays
double charginoToNeutralino1pion(const MssmSoftsusy * m) {
  double mchi1 = fabs(m->displayPhys().mch(1));
  double mneut1 = fabs(m->displayPhys().mneut(1));
  if (mchi1 < mneut1 + mpiplus) return 0.;
  if (mchi1 - mneut1 - mpiplus > hadronicScale) return 0.;

  Complex OL11 = -1.0 / root2 * m->displayDrBarPars().nBpmz.display(1, 4) *
    m->displayDrBarPars().vBpmz(1, 2).conj() +
    m->displayDrBarPars().nBpmz.display(1, 2) *
    m->displayDrBarPars().vBpmz(1, 1).conj();
  Complex OR11 = +1.0 / root2 *
    m->displayDrBarPars().nBpmz.display(1, 3).conj() *
    m->displayDrBarPars().uBpmz(1, 2) +
    m->displayDrBarPars().nBpmz.display(1, 2).conj() *
    m->displayDrBarPars().uBpmz(1, 1);

  double kpi = sqrt(lambda(sqr(mchi1), sqr(mneut1), sqr(mpiplus))) * 0.5
    / (mchi1);
  double width = sqr(fpi) * sqr(GMU) * kpi / (4.0 * PI * sqr(mchi1)) *
    ( abs((OL11 + OR11) * (OL11 + OR11)) * ( sqr(sqr(mchi1) - sqr(mneut1)) -
				   sqr(mpiplus) * sqr(mchi1 - mneut1) ) +
      abs((OL11 - OR11) * (OL11 - OR11)) * ( sqr(sqr(mchi1) - sqr(mneut1)) -
				    sqr(mpiplus) * sqr(mchi1 + mneut1))
      );

  return width; 
}

/// First, do hadronic decays
/*double charginoToNeutralino21pion(const MssmSoftsusy * m) {
  double mchi1 = fabs(m->displayPhys().mch(1));
  double mneut2 = fabs(m->displayPhys().mneut(2));
  if (mchi1 < mneut2 + mpiplus) return 0.;
  if (mchi1 - mneut2 - mpiplus > hadronicScale) return 0.;

  Complex OL21 = -1.0 / root2 * m->displayDrBarPars().nBpmz.display(2, 4) *
    m->displayDrBarPars().vBpmz(1, 2).conj() +
    m->displayDrBarPars().nBpmz.display(2, 2) *
    m->displayDrBarPars().vBpmz(1, 1).conj();
  Complex OR21 = +1.0 / root2 *
    m->displayDrBarPars().nBpmz.display(2, 3).conj() *
    m->displayDrBarPars().uBpmz(1, 2) +
    m->displayDrBarPars().nBpmz.display(2, 2).conj() *
    m->displayDrBarPars().uBpmz(1, 1);

  double kpi = sqrt(lambda(sqr(mchi1), sqr(mneut2), sqr(mpiplus))) * 0.5
    / (mchi1);
  double width = sqr(fpi) * sqr(GMU) * kpi / (4.0 * PI * sqr(mchi1)) *
    ( abs((OL21 + OR21) * (OL21 + OR21)) * ( sqr(sqr(mchi1) - sqr(mneut2)) -
				   sqr(mpiplus) * sqr(mchi1 - mneut2) ) +
      abs((OL21 - OR21) * (OL21 - OR21)) * ( sqr(sqr(mchi1) - sqr(mneut2)) -
				    sqr(mpiplus) * sqr(mchi1 + mneut2))
      );

  return width; 
  }*/


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
    }
  }
  return amplitudeW;
}

double squarkamplitudedecaycharginoW1 (double m1, double m2, double m3, double g, double gamma) {
  double squareratio, amplitudeW;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {
    squareratio = 1 - pow(fabs(m3)/m1,2) - pow(m2/m1,2);
    amplitudeW = pow(g,2)*pow(sin(gamma),2)/(16*M_PI)*squareratio*
      sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)))/m1;
  }
  return amplitudeW;
}
		  
double squarkamplitudedecaycharginoW2 (double m1, double m2, double m3, double g, double gamma) {
  double squareratio, amplitudeW;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {
    squareratio = 1 - pow(fabs(m3)/m1,2) - pow(m2/m1,2);
    amplitudeW = pow(g,2)*pow(cos(gamma),2)/(16*M_PI)*squareratio*
      sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)))/m1;      
  }
  return amplitudeW;
}


double squark1amplitudedecaycharginoW1mix (double m1, double m2, double m3, double g, double gammaL, double gammaR,  double theta, double beta, double mWboson, double runmt, double runmb, double torb) /// the variable torb depends on if it is stop (torb=1) or sbottom (torb =2) decaying and changes AprimeuW1 to AprimedW1 accordingly
{
  double squareratio=0, angular1=0, angular2=0, amplitudeW=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {
    squareratio = 1.0 - pow(fabs(m3)/m1,2) - pow(m2/m1,2);
    angular1 = squarkmixcharginocouplings(g, theta, beta, gammaL, gammaR, runmt, runmb, mWboson, m3, 0, torb)(1);
    angular2 = squarkmixcharginocouplings(g, theta, beta, gammaL, gammaR, runmt, runmb, mWboson, m3, 0, torb)(2);
    amplitudeW = 1.0/(16*PI*m1)*(angular1*squareratio + m3*m2/(pow(m1,2))*angular2)*sqrt(lambda(sqr(m1), sqr(m2), sqr(m3))); //test
  }
  return amplitudeW;
}



double squark1amplitudedecaycharginoW2mix (double m1, double m2, double m3, double g, double gammaL, double gammaR,  double theta, double beta, double mWboson, double runmt, double runmb, double torb) /// the variable torb depends on if it is stop (torb=1) or sbottom (torb =2) decaying and changes AprimeuW2 to AprimedW2 accordingly
{
  double squareratio=0, angular1=0, angular2=0, amplitudeW=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {
    squareratio = 1- pow(fabs(m3)/m1,2) - pow(m2/m1,2);
    angular1 = squarkmixcharginocouplings (g, theta, beta, gammaL, gammaR, runmt, runmb, mWboson, 0, m3, torb)(3);
    angular2 = squarkmixcharginocouplings (g, theta, beta, gammaL, gammaR, runmt, runmb, mWboson, 0, m3, torb)(4);
		
    amplitudeW = 1.0/(16*PI*m1)*sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)))*(angular1*squareratio + m3*m2/(pow(m1,2))*angular2);
  }
  return amplitudeW;
}



double squark2amplitudedecaycharginoW1mix (double m1, double m2, double m3, double g, double gammaL, double gammaR,  double theta, double beta, double mWboson, double runmt, double runmb, double torb) /// the variable torb depends on if it is stop (torb=1) or sbottom (torb =2) decaying and changes AprimeuW1 to AprimedW1 accordingly
{
  double squareratio=0, angular1=0, angular2=0, amplitudeW=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {
    squareratio = 1- pow(fabs(m3)/m1,2) - pow(m2/m1,2);
    angular1 = squarkmixcharginocouplings (g, theta, beta, gammaL, gammaR, runmt, runmb, mWboson, m3, 0, torb)(5);
    angular2 = squarkmixcharginocouplings (g, theta, beta, gammaL, gammaR, runmt, runmb, mWboson, m3, 0, torb)(6);
    amplitudeW = 1./(16*PI*m1)*sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)))*(angular1*squareratio + m3*m2/(pow(m1,2))*angular2);
   
  }
  return amplitudeW;
}

double squark2amplitudedecaycharginoW2mix (double m1, double m2, double m3, double g, double gammaL, double gammaR,  double theta, double beta, double mWboson, double runmt, double runmb, double torb) /// the variable torb depends on if it is stop (torb=1) or sbottom (torb =2) decaying and changes AprimeuW2 to AprimedW2 accordingly
{
  double squareratio=0, angular1=0, angular2=0, amplitudeW=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else {
    squareratio = 1- pow(fabs(m3)/m1,2) - pow(m2/m1,2);
    angular1 = squarkmixcharginocouplings (g, theta, beta, gammaL, gammaR, runmt, runmb, mWboson, 0, m3, torb)(7);
    angular2 = squarkmixcharginocouplings (g, theta, beta, gammaL, gammaR, runmt, runmb, mWboson, 0, m3, torb)(8);
    
    amplitudeW = 1.0/(m1*16*PI)*(angular1*squareratio +
				 m3*m2/(pow(m1,2))*angular2) *
      sqrt(lambda(sqr(m1), sqr(m2), sqr(m3))); 
  }
  return amplitudeW;
}

///neutralino takes values of 1, 2, 3, or 4 and denotes which neutralino mass eigenstate we decay into, uord takes value 1 for up quarks (or charm quarks) and -1 for down quarks (or strange quarks) in order to change the sign of the g*neutralinomixingmatrix[2][element] in AqZ accordingly
double squarkLamplitudedecayneutralino (double m1, double m2, double m3,
					double g, double gprime,
					DoubleMatrix & mixNeut,
					int neutralino, int uord) {
  double squareratio, AqZ, lam, amplitudeW;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
	  amplitudeW = 0;
    	}
	else {
	  squareratio = 1 - pow(fabs(m3)/m1,2) - pow(m2/m1,2);
	  AqZ = 1/(root2)*(uord*g*mixNeut(neutralino,2) + gprime*mixNeut(neutralino,1)/3); 
	  lam = sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
	  amplitudeW = pow(AqZ,2)/(16*PI*m1)*squareratio*lam; 
	}
	return amplitudeW;
}

double squarkRamplitudedecayneutralino (double m1, double m2, double m3, double g, double gprime, DoubleMatrix & mixNeut, int neutralino, int uord ) ///neutralino takes values of 1, 2, 3, or 4 and denotes which neutralino mass eigenstate we decay into, uord takes value 1 for up quarks (or charm quarks) and -1 for down quarks (or strange quarks) in order to change the sign and magnitude of the coefficient BqZ accordingly
{
  double squareratio, BqZ, uordchanger=0, lam, amplitudeW;
	if (fabs(m1) < fabs(m2) +fabs(m3)) {
	  amplitudeW = 0;
    	}
	else {
	  squareratio = 1 - pow(m3/m1,2) - pow(m2/m1,2);
	  if ( uord == 1) { uordchanger =1;}
	  else if (uord == -1) { uordchanger =-0.5;}
	  else {
	    throw("problem: uord must be 1 or -1 in squarkRamplitudedecayneutralino");
	  }
	  BqZ = 1/(root2)*4./3*uordchanger*gprime*mixNeut(neutralino,1); ///Following changes in AqZ suggested by SUSYHIT
	  lam = sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
	  amplitudeW = pow(BqZ,2)/(16*PI*m1)*squareratio*lam; /// need to check this formula, unsure about BqZ and also if it should be mod squared or just the squared I have done here!
	  
	}
	return amplitudeW;
}


double squark3amplitudedecayneutralino (double m1, double m2, double m3, double mWboson, double theta, double beta, DoubleMatrix & mixNeut, double g, double gp, double runmq, int squark , int oneortwo,  int neutralino) {

  double amplitudeW, lam, fq=0, alphatilda=0, betatilda=0, a=0, b=0;
  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  
  else { 
    lam = sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
    if (squark == 1) /// we have stops
      {
	fq = g*runmq/(root2*mWboson*sin(beta));
	if (oneortwo == 1) /// we have stop1s
	  {
	    alphatilda = (cos(theta)*1/(root2)*(-g*mixNeut(neutralino,2) - gp/3*mixNeut(neutralino,1)) - fq*sin(theta)*mixNeut(neutralino,4)); /// this is just simpler (rearranged) form of atopr(1,i) from susyhit *g
	    betatilda = (4/(3*root2)*gp*mixNeut(neutralino,1)*sin(theta) - fq*mixNeut(neutralino,4)*cos(theta)); /// this is just simpler (rearranged) form of btopr(1,i) from susyhit * g
	  }
     	else if (oneortwo == 2) /// we have stop2s
	  {
	    alphatilda = (-sin(theta)*1/(root2)*(mixNeut(neutralino,1)*gp/3 + mixNeut(neutralino,2)*g) + cos(theta)*fq*mixNeut(neutralino,4)); /// again simplified form of atopr(2,i) * g from susyhit
	    betatilda = -4/(3*root2)*cos(theta)*mixNeut(neutralino,1)*gp - fq*sin(theta)*mixNeut(neutralino,4); /// again simplified form of btopr(2,i) * g from susyhit
	  }
	else {
	  throw("problem: stop must be stop1 or stop2\n"); 
	}
      }
    else if ( squark ==2) /// we have sbottoms
      {
	fq = g*runmq/(root2*mWboson*cos(beta));
	
	if ( oneortwo == 1) /// we have sbottom1s
	  {
	    alphatilda = (1/root2)*cos(theta)*(-mixNeut(neutralino,1)*gp/3 + mixNeut(neutralino,2)*g) - sin(theta)*mixNeut(neutralino,3)*fq; /// just simplified form of abot(1,i) * g from susyhit
	    betatilda = sin(theta)*2/(3*root2)*(-mixNeut(neutralino,1)*gp) - cos(theta)*fq*mixNeut(neutralino,3);
	  }
	else if (oneortwo == 2) /// we have sbottom2s
	  {
	    alphatilda = sin(theta)*(mixNeut(neutralino,1)*-gp*(1/(3*root2)) + mixNeut(neutralino,2)*g*1/(root2)) + cos(theta)*fq*mixNeut(neutralino,3); /// simplified abot(2,i)*g from susyhit
	    betatilda = cos(theta)*(2/(3*root2))*gp*(mixNeut(neutralino,1)) - sin(theta)*fq*mixNeut(neutralino,3); /// simplified bbot(2,i)*g from susyhit
	  }
	  
	else {
	  throw("problem: must be sbottom1s or sbottom2s\n");
	}
      }
      
    else {
      throw("problem: third generation squarks must be stops or sbottoms\n");
    }
          
    a = 0.5*(alphatilda + betatilda);
    b = 0.5*(betatilda - alphatilda);
    double squareplus = 1 - pow((m2+m3)/m1,2);
    double squareminus = 1 - pow((m2-m3)/m1,2);
                        
    amplitudeW = 1./(m1*8*PI)*lam*(pow(a,2)*squareplus + pow(b,2)*squareminus);
  }

  return amplitudeW;
  
}
      
double squark3amplitudedecaysquark3Wboson (double m1, double m2, double m3, double g, double thetat, double thetab, int m1torb, int m1oneortwo, int m3torb, int m3oneortwo) /// m1torb tells the function if the initial squark is a stop or sbottom, m1oneortwo tells it whether it's the lighter or heavier squark, similarly for the final state squark with m3torb and m3oneortwo
{
  double amplitudeW, lam, capthetai=0, capthetaf=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  else { 
    lam = sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
    
    if(m1torb == 1) /// we have an initial stop
      {
	if(m1oneortwo == 1) /// we have an initial stop1
	  {
	    capthetai = pow(cos(thetat),2);
	  }
	else if (m1oneortwo ==2) // we have an initial stop2
	  {
	    capthetai = pow(sin(thetat),2);
	  }
	else {
	  throw("problem: m1oneortwo must be 1 or 2 in squark3amplitudedecaysquark3Wboson!\n");
	}
      }
    else if(m1torb == 2) /// we have an initial sbottom
      {
	if(m1oneortwo == 1) /// we have an initial sbottom1
	  {
	    capthetai = pow(cos(thetab),2);
	  }
	else if (m1oneortwo == 2) // we have an initial sbottom2
	  {
	    capthetai = pow(sin(thetab),2);
	  }
	else {throw("problem: m1oneortwo must be 1 or 2 in squark3amplitudedecaysquark3Wboson!\n");
	}
      }
    else { throw("problem: m1torb must be 1 or 2 in squark3amplitudedecaysquark3Wboson!");
    }
      

    if(m3torb == 1) /// we have a final state stop
      {
	if(m3oneortwo == 1) // we have a final state stop1
	  {
	    capthetaf = pow(cos(thetat),2);
	  }
	else if (m3oneortwo ==2) /// we have a final state stop2
	  {
	    capthetaf = pow(sin(thetat),2);
	  }
	else { throw("problem: m3oneortwo must be 1 or 2 in squark3amplitudedecaysquark3Wboson!");
	}
      }
    else if(m3torb == 2) /// we have a final state sbottom
      {
	if(m3oneortwo == 1) /// we have a final state sbottom1
	  {
	    capthetaf = pow(cos(thetab),2);
	  }
	else if (m3oneortwo == 2) /// we have a final state sbottom2
	  {
	    capthetaf = pow(sin(thetab),2);
	  }
	else { throw("problem: m3oneortwo must be 1 or 2 in squark3amplitudedecaysquark3Wboson!");
	}
      }
    else { throw("problem: m3torb must be 1 or 2 in squark3amplitudedecaysquark3Wboson!");
    }
    amplitudeW = sqr(g)/(32*PI*m1*sqr(m1))/pow(m2,2)*pow(lam,3)*capthetai*capthetaf;
  }
  return amplitudeW;
}
      
double squark3amplitudedecaychargedHiggssquark3 (double m1, double m2, double m3, double g, double mWboson, double beta, double thetat, double thetab, double greekmu, double At, double Ab, double mt, double mb, int t1or2, int b1or2) {
  double amplitudeW, lam, A, A11, A12, A21, A22, combo1, combo2, combo3, combo4;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    lam = sqrt(lambda(sqr(m1), sqr(m2), sqr(m3)));
    combo1 = tan(beta) + 1/(tan(beta));
    combo2 = greekmu + At/(tan(beta));
    combo3 = greekmu + Ab*tan(beta);
    combo4 = (pow(mb,2)*tan(beta) + pow(mt,2)/(tan(beta))) - pow(mWboson,2)*sin(2*beta);
      
    A11 = g/(root2*mWboson)*(mt*mb*combo1*sin(thetat)*sin(thetab) + mt*combo2*sin(thetat)*cos(thetab) + mb*combo3*sin(thetab)*cos(thetat) + combo4*cos(thetat)*cos(thetab));
    A12 = g/(root2*mWboson)*(mt*mb*combo1*sin(thetat)*(-cos(thetab)) + mt*combo2*sin(thetat)*sin(thetab) + mb*combo3*(-cos(thetab))*cos(thetat) + combo4*cos(thetat)*sin(thetab));
    A21 = g/(root2*mWboson)*(mt*mb*combo1*(-cos(thetat))*sin(thetab) + mt*combo2*(-cos(thetat))*cos(thetab) + mb*combo3*sin(thetab)*sin(thetat) + combo4*sin(thetat)*cos(thetab));
    A22 = g/(root2*mWboson)*(mt*mb*combo1*(-cos(thetat))*(-cos(thetab)) + mt*combo2*-(cos(thetat))*sin(thetab) + mb*combo3*(-cos(thetab))*sin(thetat) + combo4*sin(thetat)*sin(thetab));
    // std::cout << "t1or2 = " << t1or2 << std::endl;
    if(t1or2 == 1) /// we have an initial stop1
      {
	if (b1or2 == 1) /// we have a final state sbottom1
	  {
	    A = A11;
	  }
	else if (b1or2 == 2) /// we have a final state sbottom2
	  {
	    A = A12;
	  }
	else {
	  throw("problem:b1or2 must be 1 or 2 in squark3amplitudedecaychargedHiggssquark3\n");
	}
      }
    else if(t1or2 == 2) /// we have an initial stop2
      {
	if (b1or2 == 1) /// we have a final state sbottom1
	  {
	    A = A21;
	  }
	else if (b1or2 == 2) /// we have a final state sbottom2
	  {
	    A = A22;
	  }
	else {
	  throw("problem:b1or2 must be 1 or 2 in squark3amplitudedecaychargedHiggssquark3\n");
	}
      }
    else {
      throw("problem:t1or2 must be 1 or 2 in squark3amplitudedecaychargedHiggssquark3\n");
    }
    
    amplitudeW = pow(A,2)*lam/(16*PI* m1 * sqr(m1));
  }
  return amplitudeW;
  
}


double squark32amplitudedecayneutralHiggssquark3 (double m1, double m2, double m3, double g, double gp, double mWboson, double beta, double alpha, double thetat, double thetab, double greekmu, double At, double Ab, double mt, double mb, int torb, char phi) {
  double amplitudeW=0, squareplus, squareminus, lambda, A, B, Ah, AH, AA, Bh, BH, BA, combo1, combo2, combo3, combo4, combo5, combo6, combo7;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in squark32amplitudedecayneutralHiggssquark3\n");
    }
    combo1 = 1 - 5./3 * pow(gp/g,2);
    combo2 = -1 +1./3 * pow(gp/g,2);
    combo3 = At*cos(alpha) + greekmu*sin(alpha);
    combo4 = -At*sin(alpha) + greekmu*cos(alpha);
    combo5 = At/(tan(beta)) + greekmu;
    combo6 = -Ab*sin(alpha) - greekmu*cos(alpha);
    combo7 = Ab*cos(alpha) - greekmu*sin(alpha);
      
    Ah = g*mWboson*sin(beta+alpha)/4 * combo1*sin(2*thetat) + g*mt*cos(2*thetat)/(2*mWboson*sin(beta))*combo3;
    AH = -g*mWboson*cos(beta+alpha)/4 * combo1*sin(2*thetat) - g*mt*cos(2*thetat)/(2*mWboson*sin(beta))*combo4; /// possible minus sign issue here, both terms are same sign in susyhit -> changed second sign to - ?? whereas Tata and Baer has a signs of two terms opposite - note as it's squared anyway only relative sign between the two terms matters
    AA = g*mt*combo5/(2*mWboson);

    Bh = g*mWboson*sin(beta+alpha)/4 * combo2*sin(2*thetab) + g*mb*cos(2*thetab)/(2*mWboson*cos(beta))*combo6;
    BH = -g*mWboson*cos(beta+alpha)/4 * combo2*sin(2*thetab) + g*mb*cos(2*thetab)/(2*mWboson*cos(beta))*combo7;
    BA = g*mb*(Ab*tan(beta)+greekmu)/(2*mWboson);

    if(phi == 'h') {
      A=Ah;
      B=Bh;
    }
      
    else if (phi == 'H') {
      A=AH;
      B=BH;
    }

    else if (phi == 'A') {
      A=AA;
      B=BA;
    }
      
    else {
      throw("problem: phi must be one of h, H, or A in squark32amplitudedecayneutraliHiggssquark3 \n");
    }
    
    if (torb == 1) ///stop2 decay
      {
	amplitudeW = pow(A,2)*lambda/(16*PI*m1);
      }
    
    else if (torb == 2) ///sbottom2 decay
      {
	amplitudeW = pow(B,2)*lambda/(16*PI*m1);
      }
    else {
      throw("problem:torb must be 1 or 2 in squark32amplitudedecayneutralHiggssquark3\n");
    }
  }
  return amplitudeW;
}

	     
double squark32amplitudedecaysquark3Zboson (double m1, double m2, double m3, double g, double gp, double theta) {
  double amplitudeW, squareplus, squareminus, lambda, angular, costhetaW;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in squark32amplitudedecaysquark3Zboson\n");
    }
    angular = pow(cos(theta)*sin(theta),2);
    costhetaW = g/pow((pow(g,2)+pow(gp,2)),0.5);
    
    amplitudeW = pow(g,2)*pow(m1,3)/(64*PI*pow(m3*costhetaW,2))*pow(lambda,3)*angular;
  }
  return amplitudeW;
}

double sleptonamplitudedecayleptonneutralino (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, char LorR, int neutralino) {
  double amplitudeW, squareplus, squareminus, squareratio, lambda, A, B, C;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+fabs(m3))/m1,2);
    squareminus = 1 - pow((m2-fabs(m3))/m1,2);
    squareratio = 1 - pow(m3/m1,2) - pow(m2/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in sleptonamplitudedecayleptonneutralino\n");
    }
    A = -1/(root2) * (g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1));
    B = -root2*gp*mixNeut(neutralino,1);
    if(LorR == 'L') {
      C=A;
    }
      
    else if (LorR == 'R') {
      C=B;
    }
    else {
      throw("problem: LorR must be L or R in sleptonamplitudedecayleptonneutralino\n");
    }
      
   
    amplitudeW = pow(C,2)*m1*fabs(squareratio)/(16*PI) * lambda;
  }
  return amplitudeW;
}


double sneutrinoamplitudedecayneutrinoneutralino (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int neutralino) {
  double amplitudeW, squareplus, squareminus, squareratio, lambda, A;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+fabs(m3))/m1,2);
    squareminus = 1 - pow((m2-fabs(m3))/m1,2);
    squareratio = 1 - pow(m3/m1,2) - pow(m2/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in sneutrinoamplitudedecayneutrinoneutralino\n");
    }
    A = 1/(root2)*(g*mixNeut(neutralino,2) - gp*mixNeut(neutralino,1));

    amplitudeW = pow(A,2)*m1*squareratio*lambda/(16*PI); /// note here as m2=0 then lambda and squareratio each reduce to (1-pow(m3/m1,2)) 
      }
  return amplitudeW;
}


double sleptonamplitudedecaychargino (double m1, double m2, double m3, double g, double theta, int chargino) ///for both sleptonL decays to charginos + neutrinos and for sneutrino decays to lepton + charginos - just change theta from thetaL in first case to thetaR in second case
{
  double amplitudeW, squareplus, squareminus, squareratio, lambda, trigtheta=0;  
  
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareratio = 1 - pow(m3/m1,2) - pow(m2/m1,2);
    squareplus = 1 - pow((m2+fabs(m3))/m1,2);
    squareminus = 1 - pow((m2-fabs(m3))/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in sleptonamplitudedecaychargino\n");
    }
    if (chargino == 1) {
      trigtheta = sin(theta);
    }
    else if (chargino ==2) {
      trigtheta = cos(theta);
    }
    else {
      throw("problem: chargino must be a 1 or 2 in sleptonamplitudedecaychargino\n");
    }	    
    amplitudeW = pow(g*trigtheta,2)*m1*squareratio*lambda/(16*PI); /// note in the slepton decay to neutrinos case m2 is zero so lambda reduces to squareratio - giving squareratio squared
      }
  return amplitudeW;
}

      
double stauamplitudedecaytauneutralino (double m1, double m2, double m3, double g, double gp, double mWboson, DoubleMatrix & mixNeut, double theta, double beta, int oneortwo, int neutralino) {
  double amplitudeW, squareplus, squareminus, lambda, ftau, alphatilda=0, betatilda=0, a, b;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+fabs(m3))/m1,2);
    squareminus = 1 - pow((m2-fabs(m3))/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stauamplitudedecaytauneutralino\n");
    }
    ftau = g*m2/(root2*mWboson*cos(beta));
    if (oneortwo == 1) /// we have a stau1 decaying
      {
	alphatilda = 1/(root2)*sin(theta)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1)) + ftau*mixNeut(neutralino,3)*cos(theta);
	betatilda =  root2*gp*mixNeut(neutralino,1)*-cos(theta) + ftau*mixNeut(neutralino,3)*sin(theta);
      }
    else if (oneortwo ==2) /// we have a stau2 decaying
      {
	alphatilda = 1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*-cos(theta) + ftau*mixNeut(neutralino,3)*sin(theta);
	betatilda = root2*gp*mixNeut(neutralino,1)*sin(theta) + ftau*mixNeut(neutralino,3)*cos(theta); 
      }
    else {
      throw("problem: oneortwo must be 1 or 2 in stauamplitudedecaytauneutralino\n");
    }
    a = 0.5*(alphatilda + betatilda);
    b = 0.5*(betatilda - alphatilda);
    amplitudeW = m1*lambda/(8*PI)*(pow(a,2)*squareplus + pow(b,2)*squareminus);
  }
  return amplitudeW;
}

double stausneutrinoamplitudedecaytauneutrinoneutralino (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int neutralino) {
  double amplitudeW, squareplus, squareminus, squareratio, lambda, A;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  
  else { 
    squareratio = 1 - pow(m3/m1,2) - pow(m2/m1,2);
    squareplus = 1 - pow((m2+fabs(m3))/m1,2);
    squareminus = 1 - pow((m2-fabs(m3))/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stausneutrinoamplitudedecaytauneutrinoneutralino\n");
    }
    A = 1/(root2)*(g*-mixNeut(neutralino,2) + gp*mixNeut(neutralino,1));
    amplitudeW = m1/(16*PI)*pow(A,2)*lambda*squareratio; /// note lambda and squareratio reduce to same thing here as m2 = 0
  }
  return amplitudeW;
}

double stauamplitudedecaynutauchargino (double m1, double m2, double m3, double g, double mWboson, double theta, double thetaL, double beta, double mtau, int oneortwo, int chargino) {
  double amplitudeW = 0, squareplus = 0, squareminus = 0, squareratio = 0, lambda = 0, ftau = 0, A=0, B=0, mixingpart = 0;  
  
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareratio = 1 - pow(m3/m1,2) - pow(m2/m1,2);
    squareplus = 1 - pow((m2+fabs(m3))/m1,2);
    squareminus = 1 - pow((m2-fabs(m3))/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stauamplitudedecaynutauchargino\n");
    }
    ftau = g*mtau/(root2*mWboson*cos(beta));
      
      if (chargino ==1)	{
	A = g*sin(thetaL);
	B = -ftau*cos(thetaL);
	
      }

      else if(chargino ==2) {
	A = g*cos(thetaL);
	B = ftau*sin(thetaL);
      }
            
      else {
	throw("problem: chargino must be a 1 or 2 in stauamplitudedecaynutauchargino\n");
      }
      
      if(oneortwo == 2) { 
	mixingpart = -A*sin(theta) + B*cos(theta);
      }
      else if (oneortwo == 1) {
	mixingpart = A*cos(theta) + B*sin(theta);
      }
      else {
	throw("problem: oneortwo must be 1 or 2 in stauamplitudedecaynutauchargino\n");
      }
            
      amplitudeW = pow(mixingpart,2)*m1*lambda*squareratio/(16*PI); ///again lambda and squareratio reduce to same thing as m2 = 0
    }

  return amplitudeW;
}


double stausneutrinoamplitudedecaytauchargino (double m1, double m2, double m3, double g, double mWboson, double beta, double thetaL, double thetaR, int chargino) {
  double amplitudeW, squareplus, squareminus, squareratio, lambda, ftau, A, B, sumsquare;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareratio = 1 - pow(m3/m1,2) - pow(m2/m1,2);
    squareplus = 1 - pow((m2+fabs(m3))/m1,2);
    squareminus = 1 - pow((m2-fabs(m3))/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stausneutrinoamplitudedecaytauchargino\n");
    }
    ftau = g*m2/(root2*mWboson*cos(beta)); ///as m2 is mtau which is what we need to calculate ftau
    
    if(chargino==1) {
      A = g*sin(thetaR);
      B = -ftau*cos(thetaL);
    }
    else if (chargino==2) {
      A = g*cos(thetaR);
      B = ftau*sin(thetaL);
    }
    else {
      throw("problem: chargino must be a 1 or 2 in stausneutrinoamplitudedecaytauchargino\n");
    }   
    sumsquare = pow(A,2) + pow(B,2);
    amplitudeW = m1*lambda/(16*PI)*(sumsquare*squareratio + 4*m2*m3/(pow(m1,2))*B*A);
  }
  return amplitudeW;
}


double stauamplitudedecaysnustauHminus (double m1, double m2, double m3, double g, double mWboson, double beta, double thetatau, double mtau, double greekmu, double Atau, int oneortwo) ///Also does decay mode snustau to Hplus and stau1/2 just with m1, m2 and m3 permuted.

{
  double amplitudeW, squareplus, squareminus, lambda, A, combo1, combo2;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stauamplitudedecaysnustauHminus\n");
    }
    
    combo1 = pow(mtau,2)*tan(beta) - pow(mWboson,2)*sin(2*beta);
    combo2 = mtau*(greekmu+Atau*tan(beta));

      if (oneortwo == 1) {
	A = g/(root2*mWboson)*(combo1*sin(thetatau) - combo2*cos(thetatau));
      }
      else if (oneortwo == 2) {
	A = g/(root2*mWboson)*(-combo1*cos(thetatau) - sin(thetatau)*combo2);
      }
      else {
	throw("problem: oneortwo must be 1 or 2 in stauamplitudedecaysnustauHminus\n");
      }
      
      amplitudeW = pow(A,2)*lambda/(16*PI*m1);
  }
  return amplitudeW;
}


double stauamplitudedecaysnustauWboson (double m1, double m2, double m3, double g, double thetatau, int oneortwo) ///m3 must be mw here
  
{
  double amplitudeW, squareplus, squareminus, lambda, mixangle=0;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stauamplitudedecaysnustauWboson\n");
    }
    
    if (oneortwo == 1) {
      mixangle = sin(thetatau);
    }
      else if (oneortwo == 2) {
	mixangle = cos(thetatau);
      }
      else {
	throw("problem oneortwo must be 1 or 2 in stauamplitudedecaysnustauWboson\n");
      }
    
    amplitudeW = pow(g*mixangle,2)*pow(m1,3)/(32*PI*pow(m3,2))*pow(lambda,3);
  }
  return amplitudeW;
}

double stau2amplitudedecaystau1Zboson (double m1, double m2, double m3, double g, double gp, double thetatau) ///m3 must be mz here
  
{
  double amplitudeW, squareplus, squareminus, lambda, mixangle, costhetaW;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
    }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stau2amplitudedecaystau1Zboson\n");
    }
    mixangle = pow(cos(thetatau)*sin(thetatau),2);
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    
    amplitudeW = pow(g,2)*mixangle*pow(m1,3)*pow(lambda,3)/(64*PI*pow(m3*costhetaW,2));
  }
  return amplitudeW;
}

double stau2amplitudedecaystau1phi (double m1, double m2, double m3, double g, double gp, double thetatau, double beta, double alpha, double mWboson, double mtau, double greekmu, double Atau, char phi) {
  double amplitudeW, squareplus, squareminus, lambda, combo1, combo2, combo3, combo4, combo5, combo6, Acoeff;  

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stau2amplitudedecaystau1phi\n");
    }
    combo1 = -1 +3*pow(gp/g,2); /// is just -1 + 3*pow(tanthetaW,2)
    combo2 = g*mWboson/4;
    combo3 = g*mtau/(2*mWboson);
    combo4 = greekmu*cos(alpha) + Atau*sin(alpha);
    combo5 = -greekmu*sin(alpha) + Atau*cos(alpha);
    combo6 = greekmu + Atau*tan(beta);

    if (phi == 'h') {
      Acoeff = -combo2*sin(beta+alpha)*sin(2*thetatau)*combo1 + combo3/(cos(beta))*cos(2*thetatau)*combo4;
    }
    else if (phi == 'H') {
      Acoeff = combo2*cos(beta+alpha)*sin(2*thetatau)*combo1 - combo3/(cos(beta))*cos(2*thetatau)*combo5;
    }
    else if (phi == 'A') {
      Acoeff = combo3*combo6;
    }
    else {
      throw("problem: phi can only be an h, H or A instau2amplitudedecaystau1phi\n");
    }    
    amplitudeW = pow(Acoeff,2)*lambda/(16*PI*m1);
  }
  return amplitudeW;
}


double charginoamplitudedecayquarksquarkL (double m1, double m2, double m3, double g, double theta, int chargino) ///quark mass is m2

{
  double amplitudeW, squareplus, squareminus, alteredsquareratio, lambda, scriptA;

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    alteredsquareratio = 1 - pow(m3/m1,2) + pow(m2/m1,2); /// altered as now + quarkmass/initialmass ^2 not minus
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in charginoamplitudedecayquarksquarkL\n");
    }
    
    if (chargino == 1) /// chargino1 decaying
      {
	scriptA = g*sin(theta); /// theta is thetaL if sdownL and thetaR if supL
      }
    else if (chargino == 2) ///chargino2 decaying
      {
	scriptA = g*cos(theta); /// theta is thetaL if sdownL and thetaR if supL
      }
    else {
      throw("problem: chargino must be a 1 or 2 in charginoamplitudedecayquarksquarkL\n");
    }
    amplitudeW = 3*fabs(m1)*lambda/(32*PI)*(pow(scriptA,2)*alteredsquareratio);
  }
  return amplitudeW;
}

double charginoamplitudedecayquarksquarkmix (double m1, double m2, double m3, double g, double theta, double thetaL, double thetaR, double beta, double runmt, double runmb, double mWboson, int chargino, int upordowntypesquark, int oneortwo) ///quark mass is m2

{
  double amplitudeW, squareplus, squareminus, alteredsquareratio, lambda, combo1=0, combo2=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    alteredsquareratio = 1 - pow(m3/m1,2) + pow(m2/m1,2); ///altered as now + quarkmass/initialmass ^2 not minus
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in charginoamplitudedecayquarksquarkmix\n");
    }
    
    int torb = upordowntypesquark;

    if (chargino == 1) {
      if (oneortwo == 1) {
    	combo1 = squarkmixcharginocouplings (g, theta, beta, thetaL, thetaR, runmt, runmb, mWboson, m3, 0, torb)(1);
    	combo2 = -squarkmixcharginocouplings (g, theta, beta, thetaL, thetaR, runmt, runmb, mWboson, m3, 0, torb)(2);
      }
      else if (oneortwo == 2) {
    	combo1 = squarkmixcharginocouplings (g, theta, beta, thetaL, thetaR, runmt, runmb, mWboson, m3, 0, torb)(5);
    	combo2 = -squarkmixcharginocouplings (g, theta, beta, thetaL, thetaR, runmt, runmb, mWboson, m3, 0, torb)(6);
      }
      else {
	throw("problem: oneortwo must be a 1 or 2 in charginoamplitudedecayquarksquarkmix\n");
      }
    }
    else if (chargino == 2) {
      if (oneortwo == 1) {
    	combo1 = squarkmixcharginocouplings (g, theta, beta, thetaL, thetaR, runmt, runmb, mWboson, 0, m3, torb)(3);
    	combo2 = -squarkmixcharginocouplings (g, theta, beta, thetaL, thetaR, runmt, runmb, mWboson, 0, m3, torb)(4);
      }
      else if (oneortwo == 2) {
    	combo1 = squarkmixcharginocouplings (g, theta, beta, thetaL, thetaR, runmt, runmb, mWboson, 0, m3, torb)(7);
    	combo2 = -squarkmixcharginocouplings (g, theta, beta, thetaL, thetaR, runmt, runmb, mWboson, 0, m3, torb)(8);
      }
      else { 
	throw("problem: oneortwo must be a 1 or 2 in charginoamplitudedecayquarksquarkmix\n");
      }
    }
    else {
      throw("problem: chargino must be a 1 or 2 in charginoamplitudedecayquarksquarkmix\n");
    }

    amplitudeW = 3*fabs(m1)*lambda/(32*PI)*(combo1*alteredsquareratio + combo2*m2/m1);
  }
  return amplitudeW;
}


double charginoamplitudedecayleptonsleptonL (double m1, double m2, double m3, double g, double thetaLorR, int chargino) ///lepton mass is m2, use thetaR for decays to sneutrinos and leptons, use thetaL for decays to sleptons neutrinos

{
  double amplitudeW, squareplus, squareminus, alteredsquareratio, lambda, A=0;

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    alteredsquareratio = 1 - pow(m3/m1,2) + pow(m2/m1,2); /// altered as now + quarkmass/initialmass ^2 not minus
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in charginoamplitudedecayleptonsleptonL\n");
    }
      
    if (chargino == 1) { 
      A = -g*sin(thetaLorR);
    }
    else if (chargino == 2) {
      A = -g*cos(thetaLorR);
    }
    else {
      throw("problem: chargino must be a 1 or 2 in charginoamplitudedecayleptonsleptonL\n");
    }

    amplitudeW = fabs(m1)*lambda/(32*PI)*pow(A,2)*alteredsquareratio;
  }
  return amplitudeW;
}


double charginoamplitudedecaysnutautau (double m1, double m2, double m3, double g, double thetaL, double thetaR, double beta, double mWboson, int chargino) ///m2 must be tau mass

{
  double amplitudeW, squareplus, squareminus, alteredsquareratio, lambda, ftau, A=0, Bprimeprime=0;

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    alteredsquareratio = 1 - pow(m3/m1,2) + pow(m2/m1,2); /// altered as now + quarkmass/initialmass ^2 not minus
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);	 
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in charginoamplitudedecaysnutautau\n");
    }
    ftau = g*m2/(root2*mWboson*cos(beta));
    
    if (chargino == 1) {
      A = g*sin(thetaR);
      Bprimeprime = - ftau*cos(thetaL);
    }
    else if (chargino == 2) {
      A = g*cos(thetaR);
      Bprimeprime = ftau*sin(thetaL);
    }
    else { 
      throw("problem: chargino must be a 1 or 2 in charginoamplitudedecaysnutautau\n");
    }

    amplitudeW = fabs(m1)*lambda/(32*PI)*((pow(A,2)+pow(Bprimeprime,2))*alteredsquareratio + 4*A*Bprimeprime*m2/fabs(m1));
  }
  return amplitudeW;
}
      

double charginoamplitudedecaystaunutau (double m1, double m2, double m3, double g, double thetaL, double thetaR, double thetatau, double beta, double mWboson, double mtau, int oneortwo, int chargino) ///m2 must be nutau mass (i.e. 0)

{
  double amplitudeW, squareplus, squareminus, alteredsquareratio, lambda, ftau, A=0, Bprimeprime=0, scriptA=0;

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  
  else { 
    alteredsquareratio = 1 - pow(m3/m1,2) + pow(m2/m1,2); /// altered as now + quarkmass/initialmass ^2 not minus
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);	 /// reduces to 1 - pow(m3/m1,2) as m2 = 0
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in charginoamplitudedecaystaunutau\n");
    }
    ftau = g*mtau/(root2*mWboson*cos(beta));
      
      if (chargino == 1) {
	A = -g*sin(thetaL);
	Bprimeprime = -ftau*cos(thetaL);
	if (oneortwo == 1) /// decaying to a stau1
	  {
	    scriptA = A*sin(thetatau) + Bprimeprime*cos(thetatau);
	  }
	else if (oneortwo == 2) /// decaying to a stau2
	  {
	    scriptA = -A*cos(thetatau) + Bprimeprime*sin(thetatau);
	  }
	else {
	  throw("problem: oneortwo must be a 1 or 2 in charginoamplitudedecaystaunutau\n");
	}
      }
      else if (chargino == 2) {
	A = -g*cos(thetaL);
	Bprimeprime = ftau*sin(thetaL);
	
	if (oneortwo == 1) /// decaying to a stau1
	  {
	    scriptA = A*sin(thetatau) + Bprimeprime*cos(thetatau);
	  }
	else if (oneortwo == 2) /// decaying to a stau2
	  {
	    scriptA= -A*cos(thetatau) + Bprimeprime*sin(thetatau);
	  }
	else {
	  throw("problem: oneortwo must be a 1 or 2 in charginoamplitudedecaystaunutau\n");
	}
      }
      else {
	throw("problem: chargino must be a 1 or 2 in charginoamplitudedecaystaunutau\n");
      }      
      amplitudeW = pow(scriptA,2)*fabs(m1)*lambda*alteredsquareratio/(32*PI);
  }
  return amplitudeW;
}



double charginoamplitudedecayWbosonneutralino (double m1, double m2, double m3, double g, double thetaL, double thetaR, DoubleMatrix & mixNeut, int chargino, int neutralino) ///m2 must be Wboson mass

{
  double amplitudeW, squareplus, squareminus, squarecombo1, squarecombo2, lambda, X=0, Y=0;

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in charginoamplitudedecayWbosonneutralino\n");
    }
    squarecombo1 = pow(m1,2) + pow(m3,2) - pow(m2,2);
    squarecombo2 = (pow(pow(m1,2) - pow(m3,2),2) - pow(m2,4))/pow(m2,2);
      
    if (chargino == 1) {
      X = 0.5*((cos(thetaR)*mixNeut(neutralino,4)/(root2) - sin(thetaR)*mixNeut(neutralino,2)) - cos(thetaL)*mixNeut(neutralino,3)/(root2) - sin(thetaL)*mixNeut(neutralino,2));
      Y = 0.5*(-(cos(thetaR)*mixNeut(neutralino,4)/(root2) - sin(thetaR)*mixNeut(neutralino,2)) - cos(thetaL)*mixNeut(neutralino,3)/(root2) - sin(thetaL)*mixNeut(neutralino,2));
    }
      
    else if (chargino == 2) {
      X = 0.5*((-sin(thetaR)*mixNeut(neutralino,4)/(root2) - cos(thetaR)*mixNeut(neutralino,2)) + sin(thetaL)*mixNeut(neutralino,3)/(root2) - cos(thetaL)*mixNeut(neutralino,2));
      Y = 0.5*((sin(thetaR)*mixNeut(neutralino,4)/(root2) + cos(thetaR)*mixNeut(neutralino,2)) + sin(thetaL)*mixNeut(neutralino,3)/(root2) - cos(thetaL)*mixNeut(neutralino,2));
    }
    else {
      throw("problem: chargino must be a 1 or 2 in charginoamplitudedecayWbosonneutralino\n");
    }
    amplitudeW = pow(g,2)/(16*PI*fabs(m1))*lambda*((pow(X,2)+pow(Y,2))*(squarecombo1+squarecombo2) - 6*(pow(X,2) - pow(Y,2))*(m1)*(m3));
  }
  return amplitudeW;
}


double charginoamplitudedecayHminusneutralino (double m1, double m2, double m3, double g, double gp, double thetaL, double thetaR, double beta, DoubleMatrix & mixNeut, int chargino, int neutralino) ///m2 must be Hminus mass

{
  double amplitudeW, squareplus, squareminus, lambda, squarecombo1, A1=0, A2=0, A3=0, A4=0, a=0, b=0;

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+fabs(m3))/m1,2);
    squareminus = 1 - pow((m2-fabs(m3))/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in charginoamplitudedecayHminusneutralino\n");
    }
    squarecombo1 = pow(m1,2) + pow(m3,2) - pow(m2,2);
    
    A1 = 1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*sin(thetaR) - g*mixNeut(neutralino,4)*cos(thetaR);
    A2 = -1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*cos(thetaR) - g*mixNeut(neutralino,4)*sin(thetaR);
    A3 = 1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*sin(thetaL) + g*mixNeut(neutralino,3)*cos(thetaL);
    A4 = -1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*cos(thetaL) + g*mixNeut(neutralino,3)*sin(thetaL);
    
    if (chargino == 1) {
      a = 0.5*(-cos(beta)*A2 + sin(beta)*A4);
      b = 0.5*(-cos(beta)*A2 - sin(beta)*A4);
    }
      
    else if (chargino == 2) {
      a = 0.5*(-cos(beta)*A1 + sin(beta)*A3);
      b = 0.5*(-cos(beta)*A1 - sin(beta)*A3);
    }
    else { 
      throw("problem: chargino must be a 1 or 2 in charginoamplitudedecayHminusneutralino\n");
    }
    amplitudeW = lambda/(16*PI*fabs(m1))*((pow(a,2)+pow(b,2))*squarecombo1 + 2*(pow(a,2) - pow(b,2))*m1*m3);
  }
  return amplitudeW;
}



double chargino2amplitudedecaychargino1Zboson (double m1, double m2, double m3, double g, double gp, double thetaL, double thetaR) ///m2 must be Zboson mass

{
  double amplitudeW, squareplus, squareminus, lambda, squarecombo1, squarecombo2, matelem=0, x=0, y=0, e;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+fabs(m3))/fabs(m1),2);
    squareminus = 1 - pow((m2-fabs(m3))/fabs(m1),2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in chargino2amplitudedecaychargino1Zboson\n");
    }
    squarecombo1 = pow(m1,2) + pow(m3,2) - pow(m2,2);
    squarecombo2 = (pow(pow(m1,2) - pow(m3,2),2) - pow(m2,4))/pow(m2,2); 
    e = g*gp/(pow(pow(g,2)+pow(gp,2),0.5));
      
    x = 0.5*(sin(thetaL)*cos(thetaL) - sin(thetaR)*cos(thetaR));
    y = 0.5*(sin(thetaL)*cos(thetaL) + sin(thetaR)*cos(thetaR));
    matelem = (pow(x,2) + pow(y,2))*(squarecombo1 + squarecombo2) + 6*(pow(x,2) - pow(y,2))*m1*m3;
    
    amplitudeW = pow(e,2)*lambda/(64*PI*fabs(m1))*pow(g/gp + gp/g,2)*matelem;
  }
  return amplitudeW;
}


double chargino2amplitudedecaychargino1neutHiggs (double m1, double m2, double m3, double g, double gp, double thetaL, double thetaR, double beta, double alpha, char phi) ///m2 must be neutral Higgs mass
  
{
  double amplitudeW, squareplus, squareminus, lambda, squarecombo1, matelem=0, S=0, P=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+fabs(m3))/fabs(m1),2);
    squareminus = 1 - pow((m2-fabs(m3))/fabs(m1),2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in chargino2amplitudedecaychargino1neutHiggs\n");
    }
    squarecombo1 = pow(m1,2) - pow(m2,2) + pow(m3,2);
    
    if(phi == 'h') {
      S= 0.5*(-sin(thetaR)*sin(thetaL)*sin(alpha) - cos(thetaL)*cos(thetaR)*cos(alpha) + sin(thetaL)*sin(thetaR)*cos(alpha) + cos(thetaL)*cos(thetaR)*sin(alpha));
      P = 0.5*(sin(thetaR)*sin(thetaL)*sin(alpha) + cos(thetaL)*cos(thetaR)*cos(alpha) + sin(thetaL)*sin(thetaR)*cos(alpha) + cos(thetaL)*cos(thetaR)*sin(alpha));
      ///  throw("DECAY to h\n"); 
	}
    
    else if(phi == 'H') {
      S= 0.5*(sin(thetaR)*sin(thetaL)*cos(alpha) - cos(thetaL)*cos(thetaR)*sin(alpha) + sin(thetaL)*sin(thetaR)*sin(alpha) - cos(thetaL)*cos(thetaR)*cos(alpha));
      P = 0.5*(-sin(thetaR)*sin(thetaL)*cos(alpha) + cos(thetaL)*cos(thetaR)*sin(alpha) + sin(thetaL)*sin(thetaR)*sin(alpha) - cos(thetaL)*cos(thetaR)*cos(alpha));
      ///  throw("DECAY to H\n"); 
    }
      
    else if (phi == 'A') {
      S = 0.5*(sin(thetaR)*sin(thetaL)*sin(beta) - cos(thetaL)*cos(thetaR)*cos(beta) - sin(thetaL)*sin(thetaR)*cos(beta) + cos(thetaL)*cos(thetaR)*sin(beta));
      P = 0.5*(-sin(thetaR)*sin(thetaL)*sin(beta) + cos(thetaL)*cos(thetaR)*cos(beta) - sin(thetaL)*sin(thetaR)*cos(beta) + cos(thetaL)*cos(thetaR)*sin(beta));
      /// throw("DECAY to A\n"); 
    }
    
    else {
      throw("problem: phi must be one of h, H, or A in chargino2amplitudedecaychargino1neutHiggs\n");
    }    
    matelem = (pow(S,2) + pow(P,2))*squarecombo1 + 2*(pow(S,2)-pow(P,2))*m1*m3;
    
    amplitudeW = pow(g,2)*lambda/(32*PI*fabs(m1))*matelem;
  }
  return amplitudeW;
}


double neutralinoamplitudedecayquarksquarkLorR (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int uordtype , char LorR, int neutralino) ///m2 must be quark mass

{
  double amplitudeW, squareplus, squareminus, lambda, alteredsquareratio, A=0, B=0, C=0;

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecayquarksquarkLorR\n");
    }
    alteredsquareratio = 1 + pow(m2/m1,2) - pow(m3/m1,2);
    
    if (uordtype == 1) /// up type squark (i.e. u squark or c squark)
      {
	A = 1/(root2)*(-g*mixNeut(neutralino,2) - gp/3*mixNeut(neutralino,1));
	B = -4/(3*root2)*gp*mixNeut(neutralino,1);
      }
    
    else if (uordtype == 2) /// down type squark (i.e. d squark or s squark)
      {
	A = 1/(root2)*(g*mixNeut(neutralino,2) - gp/3*mixNeut(neutralino,1));
	B = 2/(3*root2)*gp*mixNeut(neutralino,1);
      }
    else {
      throw("problem: uordtype must be a 1 or 2 in neutralinoamplitudedecayquarksquarLorR\n");
    }
    if(LorR == 'L') /// so get wino and zino couplings as LH squark
      {
	C = A;
      }
    
    else if (LorR == 'R') /// so get only Zino couplings as RH
      {
	C = B;
      }
    else {
      throw("problem: LorR must be a L or R in neutralinoamplitudedecayquarksquarLorR\n");
    }    
    amplitudeW = 3*pow(C,2)*fabs(m1)*lambda/(32*PI)*alteredsquareratio;
  }
  return amplitudeW;
}

	  

double neutralinoamplitudedecayleptonsleptonLorR (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, char LorR, int neutralino) ///m2 must be lepton mass
  
{
  double amplitudeW, squareplus, squareminus, lambda, alteredsquareratio, A=0, B=0, C=0;
    
  
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecayleptonsleptonLorR\n");
    }
    alteredsquareratio = 1 + pow(m2/m1,2) - pow(m3/m1,2);
    
    A = -1/(root2) * (g*-mixNeut(neutralino,2) + gp*-mixNeut(neutralino,1));
    B = root2*gp*mixNeut(neutralino,1);
    
    if(LorR == 'L' ) /// so get wino and zino couplings as LH slepton
      {
	C = A;
      }
    
    else if (LorR == 'R') /// so get only Zino couplings as RH
      {
	C = B;
      }    
    else {
      throw("problem: LorR must be a L or R in neutralinoamplitudedecayleptonsleptonLorR\n");
    }  
    amplitudeW = pow(C,2)*fabs(m1)*lambda*alteredsquareratio/(32*PI);
  }
  return amplitudeW;
}



double neutralinoamplitudedecayneutrinosneutrinoL (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int neutralino) ///m2 must be neutrino mass (i.e. 0)

{
  double amplitudeW, squareplus, squareminus, lambda, alteredsquareratio, A=0;
    

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecayneutrinosneutrinoL\n");
    }
    alteredsquareratio = 1 + pow(m2/m1,2) - pow(m3/m1,2);
    
    A = 1/(root2) * (g*-mixNeut(neutralino,2) + gp*mixNeut(neutralino,1));
    
    amplitudeW = pow(A,2)*fabs(m1)*lambda*alteredsquareratio/(32*PI);
  }
  return amplitudeW;
}




double neutralinoamplitudedecaysquark3quarkmix (double m1, double m2, double m3, double mWboson, double theta, double beta, DoubleMatrix & mixNeut, double g, double gp, double runmq, int squark , int oneortwo,  int neutralino) /// m2 must be quark mass
{
  double amplitudeW, squareplus, squareminus, masscombo1, masscombo2, lambda, fq=0, alphatilda=0, betatilda=0, a=0, b=0;

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecaysquark3quarkmix\n");
    }
    masscombo1 = pow(1+m2/(m1),2) - pow(m3/m1,2);
    masscombo2 = pow(1-m2/(m1),2) - pow(m3/m1,2);      
    
    
    if (squark == 1) /// we have stops
      {
	fq = g*runmq/(root2*mWboson*sin(beta));
	if (oneortwo == 1) /// we have stop1s
	  {
	    alphatilda = (cos(theta)*1/(root2)*(-g*mixNeut(neutralino,2) - gp/3*mixNeut(neutralino,1)) - fq*sin(theta)*mixNeut(neutralino,4)); /// this is just simpler (rearranged) form of atopr(1,i) from susyhit *g
	    betatilda = (4/(3*root2)*gp*mixNeut(neutralino,1)*sin(theta) - fq*mixNeut(neutralino,4)*cos(theta)); /// this is just simpler (rearranged) form of btopr(1,i) from susyhit * g
	  }
	
	else if (oneortwo == 2) /// we have stop2s
	  {
	    alphatilda = (sin(theta)*1/(root2)*(-mixNeut(neutralino,1)*gp/3 - mixNeut(neutralino,2)*g) + cos(theta)*fq*mixNeut(neutralino,4)); /// again simplified form of atopr(2,i) * g from susyhit
	    betatilda = 4/(3*root2)*cos(theta)*-mixNeut(neutralino,1)*gp - fq*sin(theta)*mixNeut(neutralino,4); /// again simplified form of btopr(2,i) * g from susyhit
	  }
	
	else {
	  throw("problem: oneortwo must be 1 or 2 in neutralinoamplitudedecaysquark3quarkmix\n");
	}
      }
    
    else if ( squark ==2) /// we have sbottoms
      {
	fq = g*runmq/(root2*mWboson*cos(beta));
	
	if ( oneortwo == 1) /// we have sbottom1s
	  {
	    alphatilda = (1/root2)*cos(theta)*(-mixNeut(neutralino,1)*gp/3 + mixNeut(neutralino,2)*g) - sin(theta)*mixNeut(neutralino,3)*fq; /// just simplified form of abot(1,i) * g from susyhit
	    betatilda = sin(theta)*2/(3*root2)*(-mixNeut(neutralino,1)*gp) - cos(theta)*fq*mixNeut(neutralino,3);
	  }
	
	else if (oneortwo == 2) /// we have sbottom2s
	  {
	    alphatilda = sin(theta)*(mixNeut(neutralino,1)*-gp*(1/(3*root2)) + mixNeut(neutralino,2)*g*1/(root2)) + cos(theta)*fq*mixNeut(neutralino,3); /// simplified abot(2,i)*g from susyhit
	    betatilda = cos(theta)*(2/(3*root2))*gp*(mixNeut(neutralino,1)) - sin(theta)*fq*mixNeut(neutralino,3); /// simplified bbot(2,i)*g from susyhit
	  }
	else {
	  throw("problem: oneortwo must be 1 or 2 in neutralinoamplitudedecaysquark3quarkmix\n");
	}
	}
      
    else {
      throw("problem: squark must be 1 or 2 in neutralinoamplitudedecaysquark3quarkmix\n");
    }    
    a = 0.5*(alphatilda + betatilda);
    b = 0.5*(betatilda - alphatilda);
    amplitudeW = 3*lambda*fabs(m1)/(16*PI)*(pow(a,2)*masscombo1 + pow(b,2)*masscombo2);
    
  }
  return amplitudeW;
}



double neutralinoamplitudedecaystautau (double m1, double m2, double m3, double mWboson, double theta, double beta, DoubleMatrix & mixNeut, double g, double gp, int oneortwo,  int neutralino) /// m2 must be tau (i.e. lepton) mass
{

  double amplitudeW, squareplus, squareminus, lambda, ftau=0, alphatilda=0, betatilda=0, a=0, b=0, factor1=0, factor2=0;

  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecaystautau\n");
    }
    
    ftau = g*m2/(root2*mWboson*cos(beta));
    if (oneortwo == 1) /// we are decaying into a stau1
      {
	alphatilda = 1/(root2)*sin(theta)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1)) + ftau*mixNeut(neutralino,3)*cos(theta);
	betatilda =  root2*gp*mixNeut(neutralino,1)*-cos(theta) + ftau*mixNeut(neutralino,3)*sin(theta);
      }
    else if (oneortwo ==2) /// we are decaying into a stau2
      {
	alphatilda = 1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*-cos(theta) + ftau*mixNeut(neutralino,3)*sin(theta);
	betatilda = -root2*gp*mixNeut(neutralino,1)*sin(theta) - ftau*mixNeut(neutralino,3)*cos(theta);
      }
    else {
       throw("problem: oneortwo must be 1 or 2 in neutralinoamplitudedecaystautau\n"); 
    }

    a = 0.5*(alphatilda + betatilda);
    b = 0.5*(betatilda - alphatilda);
    factor1 = (pow(m1+m2,2)-pow(m3,2))/pow(m1,2);
    factor2 = (pow(m1-m2,2)-pow(m3,2))/pow(m1,2);

    amplitudeW = fabs(m1)*lambda/(16*PI)*(pow(a,2)*factor1 + pow(b,2)*factor2);
  }
  return amplitudeW;
}



double neutralinoamplitudedecaycharginoWboson (double m1, double m2, double m3, double g, double thetaL, double thetaR, DoubleMatrix & mixNeut, int neutralino,  int chargino) /// m2 must be W boson mass (expect to be pole mass 80.4...), neutralino is i in T&B whilst chargino is j
{

  double amplitudeW, squareplus, squareminus, squarecombo1=0, squarecombo2=0, lambda, X=0, Y=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+fabs(m3))/m1,2);
    squareminus = 1 - pow((m2-fabs(m3))/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecaycharginoWboson\n");
    }    
  
    squarecombo1 = pow(m1,2) + pow(m3,2) - pow(m2,2);
    squarecombo2 = (pow(pow(m1,2) - pow(m3,2),2) - pow(m2,4))/pow(m2,2);
    
    if (chargino == 1) {
      X = 0.5*((cos(thetaR)*mixNeut(neutralino,4)/(root2) - sin(thetaR)*mixNeut(neutralino,2)) - cos(thetaL)*mixNeut(neutralino,3)/(root2) - sin(thetaL)*mixNeut(neutralino,2));
      Y = 0.5*(-(cos(thetaR)*mixNeut(neutralino,4)/(root2) - sin(thetaR)*mixNeut(neutralino,2)) - cos(thetaL)*mixNeut(neutralino,3)/(root2) - sin(thetaL)*mixNeut(neutralino,2));
	}
      
    else if (chargino == 2) {
      X = 0.5*((-sin(thetaR)*mixNeut(neutralino,4)/(root2) - cos(thetaR)*mixNeut(neutralino,2)) + sin(thetaL)*mixNeut(neutralino,3)/(root2) - cos(thetaL)*mixNeut(neutralino,2));
      Y = 0.5*((sin(thetaR)*mixNeut(neutralino,4)/(root2) + cos(thetaR)*mixNeut(neutralino,2)) + sin(thetaL)*mixNeut(neutralino,3)/(root2) - cos(thetaL)*mixNeut(neutralino,2));
    }
    else {
      throw("problem: chargino must be 1 or 2 in neutralinoamplitudedecaycharginoWboson\n");
    }
    amplitudeW = pow(g,2)/(16*PI*fabs(m1))*lambda*((pow(X,2)+pow(Y,2))*(squarecombo1+squarecombo2) - 6*(pow(X,2) - pow(Y,2))*m1*m3);
      
  }
  return amplitudeW;
}


double neutralinoamplitudedecaycharginoHplus (double m1, double m2, double m3, double g, double gp, double beta, double thetaL, double thetaR, DoubleMatrix & mixNeut, int neutralino,  int chargino) /// m2 must be Hplus mass, neutralino is i in T&B whilst chargino is j
{
  
  double amplitudeW, squareplus, squareminus, squarecombo1=0, lambda, A1=0, A2=0, A3=0, A4=0,a=0, b=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
     amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecaycharginoHplus\n");
    }    
    
    squarecombo1 = pow(m1,2) + pow(m3,2) - pow(m2,2);

    A1 = 1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*sin(thetaR) - g*mixNeut(neutralino,4)*cos(thetaR);
    A2 = -1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*cos(thetaR) - g*mixNeut(neutralino,4)*sin(thetaR);
    A3 = 1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*sin(thetaL) + g*mixNeut(neutralino,3)*cos(thetaL);
    A4 = -1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*cos(thetaL) + g*mixNeut(neutralino,3)*sin(thetaL);
          
    if (chargino == 1) {
      a = 0.5*(-cos(beta)*A2 + sin(beta)*A4);
      b = 0.5*(-cos(beta)*A2 - sin(beta)*A4);
    }
     
    else if (chargino == 2) {
      a = 0.5*(-cos(beta)*A1 + sin(beta)*A3);
      b = 0.5*(-cos(beta)*A1 - sin(beta)*A3);
    }
    else {
      throw("problem: chargino must be 1 or 2 in neutralinoamplitudedecaycharginoHplus\n");
    }
    amplitudeW = lambda/(16*PI*fabs(m1))*((pow(a,2)+pow(b,2))*squarecombo1 + 2*(pow(a,2) - pow(b,2))*m1*m3);
  }
  return amplitudeW;
}


double neutralinoamplitudedecayneutralinoZboson (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, int ineutralino,  int fneutralino) /// m2 must be Z mass, ineutralino is i in T&B whilst fneutralino is j
{
  
  double amplitudeW, squareplus, squareminus, squarecombo1=0, squarecombo2=0, lambda, Wij=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecayneutralinoZboson\n");
    }    
    
    squarecombo1 = pow(m1,2) + pow(m3,2) - pow(m2,2);
    squarecombo2 = (pow(pow(m1,2)-pow(m3,2),2) -pow(m2,4))/(pow(m2,2));
    
    Wij = 0.25*pow(pow(g,2)+pow(gp,2),0.5)*(mixNeut(ineutralino,4)*mixNeut(fneutralino,4) - mixNeut(ineutralino,3)*mixNeut(fneutralino,3));

    amplitudeW = pow(Wij,2)*lambda/(4*PI*fabs(m1))*(squarecombo1 + squarecombo2 +6*m1*m3);
    }
  return amplitudeW;
}

      
double neutralinoamplitudedecayneutralinoneutHiggs (double m1, double m2, double m3, double g, double gp, DoubleMatrix & mixNeut, double mixingangle, int ineutralino,  int fneutralino, char phi) /// m2 must be phi mass, ineutralino is i in T&B whilst fneutralino is j
{

  double amplitudeW=0, squareplus, squareminus, squarecombo=0, lambda, Xij=0, Xji=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    squarecombo = pow(m1,2) + pow(m3,2) - pow(m2,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecayneutralinoneutHiggs\n");
    }    
        
    if(phi == 'h') /// here mixingangle is alpha
      {
	Xij = -0.5*(mixNeut(ineutralino,3)*-sin(mixingangle) - mixNeut(ineutralino,4)*cos(mixingangle))*(-g*mixNeut(fneutralino,2) + gp*mixNeut(fneutralino,1));
	Xji = -0.5*(mixNeut(fneutralino,3)*-sin(mixingangle) - mixNeut(fneutralino,4)*cos(mixingangle))*(-g*mixNeut(ineutralino,2) + gp*mixNeut(ineutralino,1));	
	amplitudeW = pow(Xij+Xji,2)*lambda/(16*PI*fabs(m1))*(squarecombo + 2*m1*m3);
      }
    
    else if (phi == 'H') /// here mixingangle is alpha
      {
	Xij = -0.5*(mixNeut(ineutralino,3)*cos(mixingangle) - mixNeut(ineutralino,4)*sin(mixingangle))*(-g*mixNeut(fneutralino,2) + gp*mixNeut(fneutralino,1));
	Xji = -0.5*(mixNeut(fneutralino,3)*cos(mixingangle) - mixNeut(fneutralino,4)*sin(mixingangle))*(-g*mixNeut(ineutralino,2) + gp*mixNeut(ineutralino,1));
	amplitudeW = pow(Xij+Xji,2)*lambda/(16*PI*fabs(m1))*(squarecombo + 2*m1*m3);
      }
    
    else if (phi == 'A') /// here mixingangle is beta
      {
	Xij = 0.5*(mixNeut(ineutralino,3)*sin(mixingangle) - mixNeut(ineutralino,4)*cos(mixingangle))*(-g*mixNeut(fneutralino,2) + gp*mixNeut(fneutralino,1));
	Xji = 0.5*(mixNeut(fneutralino,3)*sin(mixingangle) - mixNeut(fneutralino,4)*cos(mixingangle))*(-g*mixNeut(ineutralino,2) + gp*mixNeut(ineutralino,1));
	amplitudeW = pow(Xij+Xji,2)*lambda/(16*PI*fabs(m1))*(squarecombo - 2*m1*m3);
      }
    else { 
      throw("problem: phi must be h or H or A in neutralinoamplitudedecayneutralinoneutHiggs\n");
    }
  }
  return amplitudeW;
}
double higgslorHamplitudedecayquarkantiquark (double m1, double m2, double g, double alpha, double beta, double mWboson, int uord, char lorH, DoubleMatrix & CPEMix, bool NMSSMmodel, bool QCD, double alphas) /// uord indicates if it's an up type quark (1) or down type quark (0) to choose which trig functions to act on the mixing angles (which parts of the mixing matrices), and int lorH tells the program if it's a light higgs (l) or heavy higgs (H) decaying in order to change the trig functions acting on the mixing angles; note can also use this formula for decays of higgs to leptons as these are in this sense down-type
{

  double amplitudeW=0, squarecombo1=0, angular=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }

  else { 
    squarecombo1 = 1 - 4*pow(m2/m1,2);
    if (lorH == 'l') {
	if (uord == 1) { 
	  if (NMSSMmodel == false) {
	    angular = pow(cos(alpha)/sin(beta),2);
	  }
	  else if (NMSSMmodel == true) {
	    angular = pow(CPEMix(1,1)/sin(beta),2);
	  }
	}
	else if (uord == 0) {
	  if (NMSSMmodel == false) {
	    angular = pow(-sin(alpha)/cos(beta),2);
	  }
	  else if (NMSSMmodel == true) {
	    angular = pow(CPEMix(1,2)/cos(beta),2);
	  }
	}
	else {
	  throw("problem: uord must be 1 or 0 in higgslorHamplitudedecayquarkantiquark\n");
	}
      }
    else if (lorH == 'H') {
	if (uord == 1) {
	  if (NMSSMmodel == false) {
	    angular = pow(-sin(alpha)/sin(beta),2);
	  }
	  else if (NMSSMmodel ==true) {
	    angular = pow(CPEMix(2,1)/sin(beta),2);
	  }
	}
	else if (uord == 0) {
	  if (NMSSMmodel == false) {
	    angular = pow(cos(alpha)/cos(beta),2);
	  }
	  else if (NMSSMmodel == true) {
	    angular = pow(CPEMix(2,2)/cos(beta),2);
	  }
	}
	else {
	  throw("problem: uord must be 1 or 0 in higgslorHamplitudedecayquarkantiquark\n");
	}
    }

    else if (NMSSMmodel ==true && lorH == 'M') ///M for max i.e. H3 heaviest higgs in NMSSM
      {
	if (uord == 1) /// u
	  {
	    angular = pow(CPEMix(3,1)/sin(beta),2);
	  }
	else if (uord == 0) /// d
	  {
	    angular = pow(CPEMix(3,2)/cos(beta),2);
	  }
	else {
	  throw("problem: uord must be 1 or 0 in higgslorHamplitudedecayquarkantiquark\n");
	}
      }
    else {
      throw("problem: lorH must be h or H or M (i.e. H3), where M is only for NMSSM, in higgslorHamplitudedecayquarkantiquark\n");
    }

    amplitudeW = GFosqrt2*3*m1/(4*PI)*pow(m2,2)*angular*pow(squarecombo1,1.5);

    if (QCD == true) {
      double higgsCPevenamplitudedecayqqbarQCDcorrections (double amplitude, double alphas, double x);
      double x = 0;
      x = m2/m1;
      amplitudeW = higgsCPevenamplitudedecayqqbarQCDcorrections(amplitudeW, alphas, x);
    }
    else {
      amplitudeW = amplitudeW;
    }
  }
  return amplitudeW;
}   


double higgsCPevenamplitudedecayqqbarQCDcorrections (double amplitude, double alphas, double x) {
  double beta = 0, A = 0, amplitudeW = 0, corrfactor = 0;
  beta = pow(1-4*pow(x,2),0.5);
  A = (1+pow(beta,2))*(4*dilog((1-beta)/(1+beta)) + 2*dilog((beta-1)/(beta+1)) - 3*log((1+beta)/(1-beta))*log(2/(1+beta)) - 2*log((1+beta)/(1-beta))*log(beta)) - 3*beta*log(4/(1-pow(beta,2))) - 4*beta*log(beta);
  corrfactor = 4*alphas/(3*PI)*(A/beta + (3+34*pow(beta,2)-13*pow(beta,4))/(16*pow(beta,3))*log((1+beta)/(1-beta)) + 3/(8*pow(beta,2))*(7*pow(beta,2)-1));
  amplitudeW = amplitude*(1+corrfactor);
  return amplitudeW;
}




double higgsAamplitudedecayquarkantiquark (double m1, double m2, double g, double beta, double mWboson, int uord, bool QCD, double alphas) /// uord indicates if it's an up type quark (1) or down type quark (0) to choose which trig functions to act on the mixing angles (which parts of the mixing matrices), note can also use this formula for decays of higgs to leptons as these are in this sense down-type
{

  double amplitudeW=0, squarecombo1=0, angular=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }

  else { 
    squarecombo1 = 1 - 4*pow(m2/m1,2);
    
    if (uord == 1) { 
      angular = pow(1/(tan(beta)),2);
    }
    else if (uord == 0) {
      angular = pow(tan(beta),2);
    }
    else {
      throw("problem: uord must be 1 or 0 in higgsAamplitudedecayquarkantiquark\n");
    }
    amplitudeW = 3*GFosqrt2/(4*PI)*angular*(pow(m2,2))*m1*pow(squarecombo1,0.5);

    if (QCD == true) {
      double higgsCPoddamplitudedecayqqbarQCDcorrections (double amplitude, double alphas, double x);
      double x = 0;
      x = m2/m1;
      amplitudeW = higgsCPoddamplitudedecayqqbarQCDcorrections(amplitudeW, alphas, x);
    }
    else {
      amplitudeW = amplitudeW;
    }
  }
  
  return amplitudeW;
}

double higgsCPoddamplitudedecayqqbarQCDcorrections (double amplitude, double alphas, double x) {
  double beta = 0, A = 0, amplitudeW = 0, corrfactor = 0;
  beta = pow(1-4*pow(x,2),0.5);
  A = (1+pow(beta,2))*(4*dilog((1-beta)/(1+beta)) + 2*dilog((beta-1)/(beta+1)) - 3*log((1+beta)/(1-beta))*log(2/(1+beta)) - 2*log((1+beta)/(1-beta))*log(beta)) - 3*beta*log(4/(1-pow(beta,2))) - 4*beta*log(beta);
  corrfactor = 4*alphas/(3*PI)*(A/beta + (19+2*pow(beta,2)+3*pow(beta,4))/(16*beta)*log((1+beta)/(1-beta)) + 3./8*(7-pow(beta,2)));
  amplitudeW = amplitude*(1+corrfactor);
  return amplitudeW;
}

double higgsAamplitudedecayquarkantiquarkNMSSM (double m1, double m2, double beta, DoubleMatrix & CPOMix, int uord, int higgs, bool QCD, double alphas) /// uord indicates if it's an up type quark (1) or down type quark (0) to choose which trig functions to act on the mixing angles (which parts of the mixing matrices), note can also use this formula for decays of higgs to leptons as these are in this sense down-type
{

  double amplitudeW=0, squarecombo1=0, angular=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }

  else { 
    squarecombo1 = 1 - 4*pow(m2/m1,2);
    
    if (uord == 1) { 
      angular = pow(CPOMix(higgs,1)/(sin(beta)),2);
    }
    else if (uord == 0) {
      angular = pow(CPOMix(higgs,2)/cos(beta),2);
    }
    else {
      throw("problem: uord must be 1 or 0 in higgsAamplitudedecayquarkantiquarkNMSSM\n");
    }
    amplitudeW = 3*GFosqrt2/(4*PI)*angular*(pow(m2,2))*m1*pow(squarecombo1,0.5);

    if (QCD == true) {
      double higgsCPoddamplitudedecayqqbarQCDcorrections (double amplitude, double alphas, double x);
      double x = 0;
      x = m2/m1;
      amplitudeW = higgsCPoddamplitudedecayqqbarQCDcorrections(amplitudeW, alphas, x);
    }
    else {
      amplitudeW = amplitudeW;
    }
  }
  
  return amplitudeW;
}


double higgsphiamplitudedecayneutralinoneutralino (double m1, double m2, double m3, double g, double tanthetaW, double mixingangle, DoubleMatrix & mixNeut, int ineutralino, int fneutralino, char phi) /// phi tells it whether a "h", "H" or "A" is decaying 
{

  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, Xij=0, Xji=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);        
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsphiamplitudedecayneutralinoneutralino\n");
    }    
     
    if(phi == 'h') /// here mixingangle is alpha
      {
	Xij = -0.5*(mixNeut(ineutralino,3)*-sin(mixingangle) - mixNeut(ineutralino,4)*cos(mixingangle))*(-mixNeut(fneutralino,2) + tanthetaW*mixNeut(fneutralino,1));
	Xji = -0.5*(mixNeut(fneutralino,3)*-sin(mixingangle) - mixNeut(fneutralino,4)*cos(mixingangle))*(-mixNeut(ineutralino,2) + tanthetaW*mixNeut(ineutralino,1)); ///these are actually Xij/g and Xji/g really
	if (ineutralino == fneutralino) {
	  amplitudeW = 0.5*pow(g,2)*fabs(m1)/(8*PI)*(pow(Xij+Xji,2))*squareplus*lambda;
	}
	else 
	  amplitudeW = pow(g,2)*fabs(m1)/(8*PI)*(pow(Xij+Xji,2))*squareplus*lambda;
      }
    
    else if (phi == 'H') /// here mixingangle is alpha
      {
	Xij = -0.5*(mixNeut(ineutralino,3)*cos(mixingangle) - mixNeut(ineutralino,4)*sin(mixingangle))*(-mixNeut(fneutralino,2) + tanthetaW*mixNeut(fneutralino,1));
	Xji = -0.5*(mixNeut(fneutralino,3)*cos(mixingangle) - mixNeut(fneutralino,4)*sin(mixingangle))*(-mixNeut(ineutralino,2) + tanthetaW*mixNeut(ineutralino,1));
	if (ineutralino == fneutralino) {
	  amplitudeW = 0.5*pow(g,2)*fabs(m1)/(8*PI)*(pow(Xij+Xji,2))*squareplus*lambda;
	}
	else 
	  amplitudeW = pow(g,2)*fabs(m1)/(8*PI)*(pow(Xij+Xji,2))*squareplus*lambda;
      }
    
    else if (phi == 'A') /// here mixingangle is beta
      {
	Xij = 0.5*(mixNeut(ineutralino,3)*sin(mixingangle) - mixNeut(ineutralino,4)*cos(mixingangle))*(-mixNeut(fneutralino,2) + tanthetaW*mixNeut(fneutralino,1));
	Xji = 0.5*(mixNeut(fneutralino,3)*sin(mixingangle) - mixNeut(fneutralino,4)*cos(mixingangle))*(-mixNeut(ineutralino,2) + tanthetaW*mixNeut(ineutralino,1));
	if (ineutralino == fneutralino) {
	  amplitudeW = 0.5*pow(g,2)*fabs(m1)/(8*PI)*(pow(Xij+Xji,2))*squareminus*lambda;
	}
	else {
	  amplitudeW = pow(g,2)*fabs(m1)/(8*PI)*(pow(Xij+Xji,2))*squareminus*lambda;
	}
      }
    else {
      throw("problem: phi must be h or H or A in higgsphiamplitudedecayneutralinoneutralino\n");
    }
  }

  return amplitudeW;
}
     


double higgsAamplitudedecayneutralinoneutralinoNMSSM (double m1, double m2, double m3, double g, double tanthetaW, double lam, double kappa, DoubleMatrix & CPOMix, DoubleMatrix & mixNeut, int ineutralino, int jneutralino, int pseudoscalar)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, delta = 0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsAamplitudedecayneutralinoneutralinoNMSSM\n");
    }    

    coupling = lam/(root2)*(CPOMix(pseudoscalar,1)*(mixNeut(ineutralino,3)*mixNeut(jneutralino,5) + mixNeut(ineutralino,5)*mixNeut(jneutralino,3)) + CPOMix(pseudoscalar,2)*(mixNeut(ineutralino,4)*mixNeut(jneutralino,5) + mixNeut(ineutralino,5)*mixNeut(jneutralino,4)) + CPOMix(pseudoscalar,3)*(mixNeut(ineutralino,3)*mixNeut(jneutralino,4)+mixNeut(jneutralino,3)*mixNeut(ineutralino,4))) - root2*kappa*CPOMix(pseudoscalar,3)*mixNeut(ineutralino,5)*mixNeut(jneutralino,5) - tanthetaW*g/2*(-CPOMix(pseudoscalar,1)*(mixNeut(ineutralino,1)*mixNeut(jneutralino,4) + mixNeut(ineutralino,4)*mixNeut(jneutralino,1)) + CPOMix(pseudoscalar,2)*(mixNeut(ineutralino,1)*mixNeut(jneutralino,3) + mixNeut(ineutralino,3)*mixNeut(jneutralino,1))) - g/2*(CPOMix(pseudoscalar,1)*(mixNeut(ineutralino,2)*mixNeut(jneutralino,4) + mixNeut(ineutralino,4)*mixNeut(jneutralino,2)) - CPOMix(pseudoscalar,2)*(mixNeut(ineutralino,2)*mixNeut(jneutralino,3) + mixNeut(ineutralino,3)*mixNeut(jneutralino,2)));

    if (ineutralino == jneutralino) { delta = 1;}
    else if (ineutralino != jneutralino) { delta = 2;}

    amplitudeW = m1/(16*PI)*squareminus*lambda*pow(coupling,2)*delta;

  }
  return amplitudeW;
}


double higgsphiamplitudedecaysamechargino (double m1, double m2, double g, double thetaL, double thetaR, double alpha, double beta, int chargino, char phi) /// phi tells it whether a "h", "H" or "A" is decaying
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, S=0, index=0;
  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m2)/m1,2);
    squareminus = 1 - pow((m2-m2)/m1,2); /// so ofcourse squarminus is just 1, just keeping link with other amplitudes
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsphiamplitudedecaysamechargino\n");
    }    
    
    DoubleVector SWcoupling(6);
     for (int i=1; i<=6; i++) {
     SWcoupling(i) = 0;
   }
     
     SWcoupling = higgsphisamecharginocouplings(alpha, beta, thetaL, thetaR);
     
     if (phi == 'h' && chargino == 1) {
       S = SWcoupling(1);
       index = 3;
     }
     else if (phi == 'h' && chargino == 2) {
       S = SWcoupling(2);
       index = 3;
     }
     else if (phi == 'H' && chargino == 1) {
       S = SWcoupling(3);
       index = 3;
     }
     else if (phi == 'H' && chargino == 2) {
       S = SWcoupling(4);
       index = 3;
     }
     else if (phi == 'A' && chargino == 1) {
       S = SWcoupling(5);
       index = 1;
     }
     else if (phi == 'A' && chargino == 2) {
       S = SWcoupling(6);
       index = 1;
     }
     else {
       throw("problem: phi must be h or H or A and chargino must be 1 or 2 in higgsphiamplitudedecaysamechargino\n");
     }    
    amplitudeW = pow(g,2)/(4*PI)*pow(S,2)*fabs(m1)*pow(lambda,index);
  }  
  return amplitudeW; 
}


double higgsphiamplitudedecaysamecharginoNMSSM (double m1, double m2, double g, double thetaL, double thetaR, double lam, DoubleMatrix & CPEMix, int chargino, int higgs)
{

  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m2)/m1,2);
    squareminus = 1 - pow((m2-m2)/m1,2); /// so ofcourse squarminus is just 1, just keeping link with other amplitudes
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsphiamplitudedecaysamecharginoNMSSM\n");
    }    

    if (chargino == 1) {
      coupling = (lam/(root2)*CPEMix(higgs,3)*cos(thetaL)*cos(thetaR) + g/(root2)*(CPEMix(higgs,1)*sin(thetaL)*cos(thetaR) + CPEMix(higgs,2)*cos(thetaL)*sin(thetaR)));
    }
    else if (chargino == 2) {
      coupling = (lam/(root2)*CPEMix(higgs,3)*sin(thetaL)*sin(thetaR) - g/(root2)*(CPEMix(higgs,1)*cos(thetaL)*sin(thetaR) + CPEMix(higgs,2)*sin(thetaL)*cos(thetaR)));
    }
    else {
      throw("problem: phi must be h or H or A in higgsphiamplitudedecaysamecharginoNMSSM\n");
    }
    
    amplitudeW = m1/(8*PI)*pow(lambda,3)*pow(coupling,2);

  }
  return amplitudeW;
}


double higgsAamplitudedecaysamecharginoNMSSM (double m1, double m2, double g, double thetaL, double thetaR, double alpha, double lam, DoubleMatrix & CPOMix, int chargino, int pseudoscalar)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, S=0;

  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m2)/m1,2);
    squareminus = 1 - pow((m2-m2)/m1,2); /// so ofcourse squarminus is just 1, just keeping link with other amplitudes
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsAamplitudedecaysamecharginoNMSSM\n");
    }    

    if (chargino == 1) {
      S = (lam/(root2)*CPOMix(pseudoscalar,3)*cos(thetaL)*cos(thetaR) - g/(root2)*(CPOMix(pseudoscalar,1)*sin(thetaL)*cos(thetaR) + CPOMix(pseudoscalar,2)*cos(thetaL)*sin(thetaR)));
    }
    else if (chargino == 2) {
      S = (lam/(root2)*CPOMix(pseudoscalar,3)*sin(thetaL)*sin(thetaR) + g/(root2)*(CPOMix(pseudoscalar,1)*cos(thetaL)*sin(thetaR) + CPOMix(pseudoscalar,2)*sin(thetaL)*cos(thetaR))); 
    }
    else { 
      throw("problem: chargino must be 1 or 2 in higgsAamplitudedecaysamecharginoNMSSM\n");
    }
    amplitudeW = m1/(8*PI)*lambda*pow(S,2);

  }
  return amplitudeW;
}

double higgsphiamplitudedecaydiffcharginoNMSSM (double m1, double m2, double m3, double g, double thetaL, double thetaR, double lam, DoubleMatrix & CPEMix, int higgs)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling1=0, coupling2=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsphiamplitudedecaydiffcharginoNMSSM\n");
    }  

    coupling1 = (lam/(root2)*CPEMix(higgs,3)*cos(thetaL)*sin(thetaR) + g/(root2)*(CPEMix(higgs,1)*sin(thetaL)*sin(thetaR) - CPEMix(higgs,2)*cos(thetaL)*cos(thetaR)));
    coupling2 = (lam/(root2)*CPEMix(higgs,3)*sin(thetaL)*cos(thetaR) - g/(root2)*(CPEMix(higgs,1)*cos(thetaL)*cos(thetaR) - CPEMix(higgs,2)*sin(thetaL)*sin(thetaR)));

    amplitudeW = m1/(16*PI)*lambda*((pow(coupling1,2)+pow(coupling2,2))*0.5*(squareplus+squareminus) + coupling1*coupling2*(squareminus-squareplus));
  }
  return amplitudeW;
}
 

double higgsphiamplitudedecaydifchargino (double m1, double m2, double m3, double g, double thetaL, double thetaR, double alpha, double beta, char phi) /// phi tells it whether a "h", "H" or "A" is decaying
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, S=0, P=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2); 
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsphiamplitudedecaydifchargino\n");
    }  
    
    DoubleVector SPWcoupling(6);
     for (int i=1; i<=6; i++) {
     SPWcoupling(i) = 0;
   }
     
     SPWcoupling = higgsphidifcharginocouplings(alpha, beta, thetaL, thetaR);
     
     if (phi == 'h') {
       S = SPWcoupling(1);
       P = SPWcoupling(2);
     }
     else if (phi == 'H') {
       S = SPWcoupling(3);
       P = SPWcoupling(4);
     }
     else if (phi == 'A') {
       S = SPWcoupling(5);
       P = SPWcoupling(6);
     }
     else {
       throw("problem: phi must be h or H or A in higgsphiamplitudedecaydifchargino\n");
     }     
     amplitudeW = pow(g,2)/(16*PI)*fabs(m1)*lambda*(pow(S,2)*squareplus + pow(P,2)*squareminus);
  }  
  
  return amplitudeW; 
};

double higgsAamplitudedecaydifcharginoNMSSM (double m1, double m2, double m3, double g, double thetaL, double thetaR, double alpha, double lam, DoubleMatrix & CPOMix, int pseudoscalar)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, C1 = 0, C2 = 0;
    
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2); 
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsAamplitudedecaydifcharginoNMSSM\n");
    }  

    C1 = (lam/root2*CPOMix(pseudoscalar,3)*cos(thetaL)*sin(thetaR) - pow(2,-0.5)*g*(CPOMix(pseudoscalar,1)*sin(thetaL)*sin(thetaR) - CPOMix(pseudoscalar,2)*cos(thetaL)*cos(thetaR)));

    C2 = (lam/root2*CPOMix(pseudoscalar,3)*sin(thetaL)*cos(thetaR) - pow(2,-0.5)*g*(CPOMix(pseudoscalar,1)*-cos(thetaL)*cos(thetaR) + CPOMix(pseudoscalar,2)*sin(thetaR)*sin(thetaL)));
    
    amplitudeW = lambda/(16*PI)*m1*((pow(C1,2) + pow(C2,2))*0.5*(squareplus+squareminus) + 4*C1*C2*0.25*(squareminus-squareplus));

  }
  return amplitudeW;
}


double higgshamplitudedecayAA (double m1, double m2, double g, double gp, double alpha, double beta, double mWboson) ///calculates partial width for h->AA
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, costhetaW=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m2)/m1,2);
    squareminus = 1 - pow((m2-m2)/m1,2); ///obviously 1, just for comparison and consistency with other decay modes
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecayAA\n");
    }  
    
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    coupling = g*mWboson/(4*pow(costhetaW,2))*sin(beta+alpha)*cos(2*beta);
    
    amplitudeW = pow(coupling,2)/(8*PI*fabs(m1))*lambda;
  }
  return amplitudeW;
}


double higgsHamplitudedecayhh (double m1, double m2, double g, double gp, double alpha, double beta, double mWboson) ///calculates partial width for h->hh
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, costhetaW=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m2)/m1,2);
    squareminus = 1 - pow((m2-m2)/m1,2); ///obviously 1, just for comparison and consistency with other decay modes
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecayhh\n");
    }  
    
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    coupling = g*mWboson/(4*pow(costhetaW,2))*(cos(2*alpha)*cos(beta+alpha) - 2*sin(2*alpha)*sin(beta+alpha));
    
    amplitudeW = pow(coupling,2)/(8*PI*fabs(m1))*lambda;
  }
  return amplitudeW;
}


double higgsHamplitudedecayAA (double m1, double m2, double g, double gp, double alpha, double beta, double mWboson) ///calculates partial width for H->AA
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, costhetaW=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m2)/m1,2);
    squareminus = 1 - pow((m2-m2)/m1,2); ///obviously 1, just for comparison and consistency with other decay modes
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecayAA\n");
    }  
    
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    coupling = -g*mWboson/(4*pow(costhetaW,2))*(cos(beta+alpha)*cos(2*beta));
    
    amplitudeW = pow(coupling,2)/(8*PI*fabs(m1))*lambda;
  }
  return amplitudeW;
}


double higgsHamplitudedecayHplusHminus (double m1, double m2, double g, double gp, double alpha, double beta, double mWboson) ///calculates partial width for H->H+H-
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, costhetaW=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m2)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m2)/m1,2);
    squareminus = 1 - pow((m2-m2)/m1,2); ///obviously 1, just for comparison and consistency with other decay modes
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecayHplusHminus\n");
    }  
    
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    coupling = g*mWboson*(cos(beta-alpha) - cos(beta+alpha)*cos(2*beta)/(2*pow(costhetaW,2)));
    
    amplitudeW = pow(coupling,2)/(16*PI*fabs(m1))*lambda;
  }
  return amplitudeW;
}


double higgshamplitudedecayhiggsAZboson (double m1, double m2, double m3, double g, double gp, double alpha, double beta) ///calculates partial width for h->AZ, m2 must be Zboson mass
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, costhetaW=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecayhiggsAZboson\n");
    }  
    
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    coupling =(g/(costhetaW))*(cos(beta-alpha));
    
    amplitudeW = pow(coupling,2)*fabs(pow(m1,3))/(64*PI*pow(m2,2))*pow(lambda,3);
  }
  return amplitudeW;
}


double higgsHamplitudedecayhiggsAZboson (double m1, double m2, double m3, double g, double gp, double alpha, double beta) ///calculates partial width for H->AZ, m2 must be Zboson mass
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, costhetaW=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecayhiggsAZboson\n");
    }  
    
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    coupling =(g/(costhetaW))*(sin(beta-alpha));
    
    amplitudeW = pow(coupling,2)*fabs(pow(m1,3))/(64*PI*pow(m2,2))*pow(lambda,3);
  }
  return amplitudeW;
}



double higgsAamplitudedecayhiggshZboson (double m1, double m2, double m3, double g, double gp, double alpha, double beta) ///calculates partial width for A->hZ, m2 must be Zboson mass
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, costhetaW=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsAamplitudedecayhiggshZboson\n");
    }  
    
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    coupling =(g/(costhetaW))*(cos(beta-alpha));
    
    amplitudeW = pow(coupling,2)*fabs(pow(m1,3))/(64*PI*pow(m2,2))*pow(lambda,3);
  }
  return amplitudeW;
}




double higgsAamplitudedecayhiggsHZboson (double m1, double m2, double m3, double g, double gp, double alpha, double beta) ///calculates partial width for A->HZ, m2 must be Zboson mass
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, costhetaW=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsAamplitudedecayhiggsHZboson\n");
    }  
    
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    coupling =(g/(costhetaW))*(sin(beta-alpha));
    
    amplitudeW = pow(coupling,2)*fabs(pow(m1,3))/(64*PI*pow(m2,2))*pow(lambda,3);
  }
  return amplitudeW;
}


double higgsAamplitudedecayhiggshorHZbosonNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double thetaA, DoubleMatrix & CPEMix, int pseudoscalar, int higgs) ///m2 must be Z mass
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsAamplitudedecayhiggshorHZbosonNMSSM\n");
    }  
    
    if (pseudoscalar == 1 && higgs == 1) {
      coupling = (CPEMix(1,1)*cos(beta) - CPEMix(1,2)*sin(beta))*cos(thetaA);
    }
    else if (pseudoscalar == 2 && higgs == 1) {
      coupling = (CPEMix(1,1)*cos(beta) - CPEMix(1,2)*sin(beta))*sin(thetaA);
    }
    else if (pseudoscalar == 1 && higgs == 2) {
      coupling = (CPEMix(2,1)*cos(beta) - CPEMix(2,2)*sin(beta))*cos(thetaA);
    }
    else if (pseudoscalar == 2 && higgs == 2) {
      coupling = (CPEMix(2,1)*cos(beta) - CPEMix(2,2)*sin(beta))*sin(thetaA);
    }
    else if (pseudoscalar == 1 && higgs == 3) {
      coupling = (CPEMix(3,1)*cos(beta) - CPEMix(3,2)*sin(beta))*cos(thetaA);
    }
    else if (pseudoscalar == 2 && higgs == 3) {
      coupling = (CPEMix(3,1)*cos(beta) - CPEMix(3,2)*sin(beta))*sin(thetaA);
    }
    else { 
      throw("problem: higgs must be 1 or 2 or 3 and pseudoscalar must be 1 or 2 in higgsAamplitudedecayhiggshorHZbosonNMSSM\n");
    }

    amplitudeW = GFosqrt2*pow(m1,3)/(8*PI)*pow(lambda,3)*pow(coupling,2);
  }
  return amplitudeW;
}



double higgshamplitudedecay2squarksamehand (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mupq, double mdownq, int sq) ///calculates partial width for h->squark squark with no mixing and squarks of same handedness, therefore for first two generations, int sq tells the function which squarks it's decaying into - uL, dL, uR, dR for sq = 1,2,3,4 respectively
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecay2squarksamehand\n");
    }  
     DoubleVector hsqsqcoupling(4);
     for (int i=1; i<=4; i++) {
     hsqsqcoupling(i) = 0;
     }

     hsqsqcoupling = higgshsquarksamehandcouplings(mWboson, g, gp, alpha, beta, mupq, mdownq);
     if (sq == 1) {
      coupling = hsqsqcoupling(1);
    }
    else if (sq == 2) {
      coupling = hsqsqcoupling(2);
    }
    else if (sq == 3) {
      coupling = hsqsqcoupling(3);
    }
    else if (sq == 4) {
      coupling = hsqsqcoupling(4);
    }    
    else {
      throw("problem: sq must be 1, 2, 3, or 4 in higgshamplitudedecay2squarksamehand\n");
    }
    
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgshamplitudedecay2squarksamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sq)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
    
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecay2squarksamehandNMSSM\n");
    }  

    if (sq == 1) ///uLuL
      {
	coupling = g*(mWboson*(0.5-pow(gp/g,2)/6)*(sin(beta)*CPEMix(1,1) - cos(beta)*CPEMix(1,2)) - pow(mq,2)*CPEMix(1,1)/(mWboson*sin(beta)));
      }
    else if (sq == 2) ///dLdL
      {
	coupling = g*(mWboson*(-0.5-pow(gp/g,2)/6)*(sin(beta)*CPEMix(1,1) - cos(beta)*CPEMix(1,2)) - pow(mq,2)*CPEMix(1,2)/(mWboson*cos(beta)));
      }
    else if (sq == 3) ///uRuR
      {
	coupling = g*(2*mWboson/3*pow(gp/g,2)*(sin(beta)*CPEMix(1,1) - cos(beta)*CPEMix(1,2)) - pow(mq,2)*CPEMix(1,1)/(mWboson*sin(beta)));
      }
    else if (sq == 4) ///dRdR
      {
	coupling = g*(-mWboson/3*pow(gp/g,2)*(sin(beta)*CPEMix(1,1) - cos(beta)*CPEMix(1,2)) - pow(mq,2)*CPEMix(1,2)/(mWboson*cos(beta)));
      }
    else {
      throw("problem: sq must be 1, 2, 3, or 4 in higgshamplitudedecay2squarksamehandNMSSM\n");
    }    
    amplitudeW = 3*lambda*pow(coupling,2)/(16*PI*m1);
  }
  return amplitudeW;
}

double higgshamplitudedecay2sleptonsamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sl)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
    
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecay2sleptonsamehandNMSSM\n");
    }  

    if (sl == 1) ///snuLsnuL - drop terms proportional to ml^2 here
      {
	coupling = g*(mWboson*(0.5+pow(gp/g,2)/2)*(sin(beta)*CPEMix(1,1) - cos(beta)*CPEMix(1,2)));
      }
    else if (sl == 2) ///slLslL
      {
	coupling = g*(mWboson*(0.5-pow(gp/g,2)/2)*(sin(beta)*CPEMix(1,1) - cos(beta)*CPEMix(1,2)));
      }
    else if (sl == 3) ///slRslR
      {
	coupling = g*(mWboson*pow(gp/g,2)*(sin(beta)*CPEMix(1,1) - cos(beta)*CPEMix(1,2)));
      }
    else {
      throw("problem: sl must be 1, 2, or 3 in higgshamplitudedecay2sleptonsamehandNMSSM\n");
    }    
    amplitudeW = lambda*pow(coupling,2)/(16*PI*m1);

  }
  return amplitudeW;
}


double higgsHamplitudedecay2sleptonsamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sl)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
    
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);

    if (squareplus*squareminus < 0) {
      throw("lambda will give nan in higgsHamplitudedecay2sleptonsamehandNMSSM");
    }  

    if (sl == 1) ///snuLsnuL - drop terms proportional to ml^2 here
      {
	coupling = g*(mWboson*(0.5+pow(gp/g,2)/2)*(sin(beta)*CPEMix(2,1) - cos(beta)*CPEMix(2,2)));
      }
    else if (sl == 2) ///slLslL
      {
	coupling = g*(mWboson*(0.5-pow(gp/g,2)/2)*(sin(beta)*CPEMix(2,1) - cos(beta)*CPEMix(2,2)));
      }
    else if (sl == 3) ///slRslR
      {
	coupling = g*(mWboson*pow(gp/g,2)*(sin(beta)*CPEMix(2,1) - cos(beta)*CPEMix(2,2)));
      }
    else {
      throw("problem: sl must be 1, 2, or 3 in higgsHamplitudedecay2sleptonsamehandNMSSM\n");
    }
    
    amplitudeW = lambda*pow(coupling,2)/(16*PI*m1);
  }
  return amplitudeW;
}

double higgsH3amplitudedecay2sleptonsamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sl)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
    
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsH3amplitudedecay2sleptonsamehandNMSSM\n");
    }  

    if (sl == 1) ///snuLsnuL - drop terms proportional to ml^2 here
      {
	coupling = g*(mWboson*(0.5+pow(gp/g,2)/2)*(sin(beta)*CPEMix(3,1) - cos(beta)*CPEMix(3,2)));
      }
    else if (sl == 2) ///slLslL
      {
	coupling = g*(mWboson*(0.5-pow(gp/g,2)/2)*(sin(beta)*CPEMix(3,1) - cos(beta)*CPEMix(3,2)));
      }
    else if (sl == 3) ///slRslR
      {
	coupling = g*(mWboson*pow(gp/g,2)*(sin(beta)*CPEMix(3,1) - cos(beta)*CPEMix(3,2)));
      }
    else {
      throw("problem: sl must be 1, 2, or 3 in higgsH3amplitudedecay2sleptonsamehandNMSSM\n");
    }
    
    amplitudeW = lambda*pow(coupling,2)/(16*PI*m1); 
  }
  return amplitudeW;
}

double higgshamplitudedecay2squarkdiffhand (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson, double mupq, double mdownq, double greekmu, double Aup, double Adown, int sq) ///calculates partial width for h->squark squark with no mixing and squarks of different handedness, therefore for first two generations, int sq tells the function which squarks it's decaying into - uL, dL, uR, dR for sq = 1,2,3,4 respectively
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecay2squarkdiffhand\n");
    }  

     DoubleVector hsqsqcoupling(2);
     for (int i=1; i<=2; i++) {
     hsqsqcoupling(i) = 0;
     }

     hsqsqcoupling = higgshsquarkdiffhandcouplings(mWboson, g, alpha, beta, mupq, mdownq, greekmu, Aup, Adown);
   
    if (sq == 1) {
      coupling = hsqsqcoupling(1);
    }
    else if (sq == 2) {
      coupling = hsqsqcoupling(2);
    }
    else {
      throw("problem: sq must be 1, or 2 in higgshamplitudedecay2squarkdiffhand\n");
    }
        
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgshamplitudedecay2squarkdiffhandNMSSM (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson,double mq, double Aq, double mueff, double lam, DoubleMatrix & CPEMix, int sq, int higgs)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
    
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecay2squarkdiffhandNMSSM\n");
    }  

    if (sq == 1) ///uLuR
      {
	coupling = mq*g/(2*mWboson*sin(beta))*(Aq*CPEMix(higgs,1) - mueff*CPEMix(higgs,2) - lam*root2*mWboson*cos(beta)/g*CPEMix(higgs,3));
      }
    if (sq == 2) ///dLdR
      {
	coupling = mq*g/(2*mWboson*cos(beta))*(Aq*CPEMix(higgs,2) - mueff*CPEMix(higgs,1) - lam*root2*mWboson*sin(beta)/g*CPEMix(higgs,3));
      }
    else {
      throw("problem: sq must be 1, or 2 in higgshamplitudedecay2squarkdiffhandNMSSM\n");
    }
    amplitudeW = 3/(16*PI*m1)*lambda*pow(coupling,2);
  }
  return amplitudeW;
}



double higgsHamplitudedecay2squarksamehand (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mupq, double mdownq, int sq) ///calculates partial width for H->squark squark with no mixing and squarks of same handedness, therefore for first two generations, int sq tells the function which squarks it's decaying into - uL, dL, uR, dR for sq = 1,2,3,4 respectively
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecay2squarksamehand\n");
    }  

     DoubleVector Hsqsqcoupling(4);
     for (int i=1; i<=4; i++) {
     Hsqsqcoupling(i) = 0;
     }

     Hsqsqcoupling = higgsHsquarksamehandcouplings(mWboson, g, gp, alpha, beta, mupq, mdownq);

    if (sq == 1) {
      coupling = Hsqsqcoupling(1);
    }
    else if (sq == 2) {
      coupling = Hsqsqcoupling(2);
    }
    else if (sq == 3) {
      coupling = Hsqsqcoupling(3);
    }
    else if (sq == 4) {
      coupling = Hsqsqcoupling(4);
    }    
    else {
      throw("problem: sq must be 1, 2, 3 or 4 in higgsHamplitudedecay2squarksamehand\n");
    }
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgsHamplitudedecay2squarksamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sq)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
    
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecay2squarksamehandNMSSM\n");
    }  

    if (sq == 1) ///uLuL
      {
	coupling = g*(mWboson*(0.5-pow(gp/g,2)/6)*(sin(beta)*CPEMix(2,1) - cos(beta)*CPEMix(2,2)) - pow(mq,2)*CPEMix(2,1)/(mWboson*sin(beta)));
      }
    else if (sq == 2) ///dLdL
      {
	coupling = g*(mWboson*(-0.5-pow(gp/g,2)/6)*(sin(beta)*CPEMix(2,1) - cos(beta)*CPEMix(2,2)) - pow(mq,2)*CPEMix(2,2)/(mWboson*cos(beta)));
      }
    else if (sq == 3) ///uRuR
      {
	coupling = g*(2*mWboson/3*pow(gp/g,2)*(sin(beta)*CPEMix(2,1) - cos(beta)*CPEMix(2,2)) - pow(mq,2)*CPEMix(2,1)/(mWboson*sin(beta)));
      }
    else if (sq == 4) ///dRdR
      {
	coupling = g*(-mWboson/3*pow(gp/g,2)*(sin(beta)*CPEMix(2,1) - cos(beta)*CPEMix(2,2)) - pow(mq,2)*CPEMix(2,2)/(mWboson*cos(beta)));
      }
    else {
      throw("problem: sq must be 1, 2, 3 or 4 in higgsHamplitudedecay2squarksamehandNMSSM\n");
    }
    
    amplitudeW = 3*lambda*pow(coupling,2)/(16*PI*m1); 
  }
  return amplitudeW;
}

double higgsH3amplitudedecay2squarksamehandNMSSM (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson,double mq, DoubleMatrix & CPEMix, int sq)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
    
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsH3amplitudedecay2squarksamehandNMSSM\n");
    }  
    
    if (sq == 1) ///uLuL
      {
	coupling = g*(mWboson*(0.5-pow(gp/g,2)/6)*(sin(beta)*CPEMix(3,1) - cos(beta)*CPEMix(3,2)) - pow(mq,2)*CPEMix(3,1)/(mWboson*sin(beta)));
      }
    else if (sq == 2) ///dLdL
      {
	coupling = g*(mWboson*(-0.5-pow(gp/g,2)/6)*(sin(beta)*CPEMix(3,1) - cos(beta)*CPEMix(3,2)) - pow(mq,2)*CPEMix(3,2)/(mWboson*cos(beta)));
      }
    else if (sq == 3) ///uRuR
      {
	coupling = g*(2*mWboson/3*pow(gp/g,2)*(sin(beta)*CPEMix(3,1) - cos(beta)*CPEMix(3,2)) - pow(mq,2)*CPEMix(3,1)/(mWboson*sin(beta)));
      }
    else if (sq == 4) ///dRdR
      {
	coupling = g*(-mWboson/3*pow(gp/g,2)*(sin(beta)*CPEMix(3,1) - cos(beta)*CPEMix(3,2)) - pow(mq,2)*CPEMix(3,2)/(mWboson*cos(beta)));
      }
    else {
      throw("problem: sq must be 1, 2, 3 or 4 in higgsH3amplitudedecay2squarksamehandNMSSM\n");
    }    
    amplitudeW = 3*lambda*pow(coupling,2)/(16*PI*m1); 
  }
  return amplitudeW;
}


double higgsHamplitudedecay2squarkdiffhand (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson, double mupq, double mdownq, double greekmu, double Aup, double Adown, int sq) ///calculates partial width for H->squark squark with no mixing and squarks of different handedness, therefore for first two generations, int sq tells the function which squarks it's decaying into - uL, dL, uR, dR for sq = 1,2,3,4 respectively
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecay2squarkdiffhand\n");
    }  

     DoubleVector Hsqsqcoupling(2);
     for (int i=1; i<=2; i++) {
     Hsqsqcoupling(i) = 0;
     }

     Hsqsqcoupling = higgsHsquarkdiffhandcouplings(mWboson, g, alpha, beta, mupq, mdownq, greekmu, Aup, Adown);
   
    if (sq == 1) {
      coupling = Hsqsqcoupling(1);
    }
    else if (sq == 2) {
      coupling = Hsqsqcoupling(2);
    }
    else { 
      throw("problem: sq must be 1 or 2 in higgsHamplitudedecay2squarkdiffhand\n");
    }        
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgshamplitudedecay2sleptonsamehand (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mel, int sl) ///calculates partial width for h->slepton slepton with no mixing and sleptons of same handedness, therefore for first two generations, int sl tells the function which sleptons it's decaying into - nuL, eL, eR for sl = 1,2,3 respectively
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecay2sleptonsamehand\n");
    }  

     DoubleVector hslslcoupling(3);
     for (int i=1; i<=3; i++) {
     hslslcoupling(i) = 0;
     }

     hslslcoupling = higgshsleptonsamehandcouplings (mWboson, g, gp, alpha, beta, mel);

    if (sl == 1) {
      coupling = hslslcoupling(1);
    }
    else if (sl == 2) {
      coupling = hslslcoupling(2);
    }
    else if (sl == 3) {
      coupling = hslslcoupling(3);
    }
    else {
      throw("problem: sl must be 1, 2, or 3 in higgshamplitudedecay2sleptonsamehand\n");
    }        
    amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgshamplitudedecay2sleptondiffhand (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson, double mel, double greekmu, double Aelectron, int sl) ///calculates partial width for h->slepton slepton with no mixing and sleptons of different handedness, therefore for first two generations, int sl tells the function which sleptons it's decaying into - nuL, eL, eR for sl = 1,2,3 respectively
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecay2sleptondiffhand\n");
    }  

     DoubleVector hslslcoupling(1);
     for (int i=1; i<=1; i++) {
     hslslcoupling(i) = 0;
     }

     hslslcoupling = higgshsleptondiffhandcouplings (mWboson, g, alpha, beta, mel, greekmu, Aelectron);
   
    if (sl == 1) {
      coupling = hslslcoupling(1);
    }
    else {
      throw("problem: sl must be 1 in higgshamplitudedecay2sleptondiffhand\n");
    } 
    amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgsHamplitudedecay2sleptonsamehand (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mel, int sl) ///calculates partial width for H->slepton slepton with no mixing and sleptons of same handedness, therefore for first two generations, int sl tells the function which sleptons it's decaying into - nuL, eL, eR for sl = 1,2,3 respectively
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecay2sleptonsamehand\n");
    }  

     DoubleVector Hslslcoupling(3);
     for (int i=1; i<=3; i++) {
     Hslslcoupling(i) = 0;
     }

     Hslslcoupling = higgsHsleptonsamehandcouplings (mWboson, g, gp, alpha, beta, mel);

    if (sl == 1) {
      coupling = Hslslcoupling(1);
    }
    else if (sl == 2) {
      coupling = Hslslcoupling(2);
    }
    else if (sl == 3) {
      coupling = Hslslcoupling(3);
    }
    else {
      throw("problem: sl must be 1, 2, or 3 in higgsHamplitudedecay2sleptonsamehand\n");
    }   
    amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgsHamplitudedecay2sleptondiffhand (double m1, double m2, double m3, double g, double alpha, double beta, double mWboson, double mel, double greekmu, double Aelectron, int sl) ///calculates partial width for H->slepton slepton with no mixing and sleptons of different handedness, therefore for first two generations, int sl tells the function which sleptons it's decaying into - nuL, eL, eR for sl = 1,2,3 respectively
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecay2sleptondiffhand\n");
    }  

     DoubleVector Hslslcoupling(1);
     for (int i=1; i<=1; i++) {
     Hslslcoupling(i) = 0;
     }

     Hslslcoupling = higgsHsleptondiffhandcouplings (mWboson, g, alpha, beta, mel, greekmu, Aelectron);
   
    if (sl == 1) {
      coupling = Hslslcoupling(1);
    }
    else {
      throw("problem: sl must be 1 in higgsHamplitudedecay2sleptondiffhand\n");
    }        
    amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgshamplitudedecaystop1stop1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for h->stop1 antistop1
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecaystop1stop1\n");
    }  

     DoubleVector hst1st1samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     hst1st1samehandcoupling(i) = 0;
     }

     hst1st1samehandcoupling = higgshsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = hst1st1samehandcoupling(1);
     couplingRR = hst1st1samehandcoupling(3);

     DoubleVector hst1st1diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     hst1st1diffhandcoupling(i) = 0;
     }

     hst1st1diffhandcoupling = higgshsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = hst1st1diffhandcoupling(1);

     coupling = couplingLL*pow(cos(theta),2) + couplingRR*pow(sin(theta),2) - 2*couplingLR*cos(theta)*sin(theta);
        
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgshamplitudedecaystop2stop2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for h->stop2 antistop2
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecaystop2stop2\n");
    }  

     DoubleVector hst2st2samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     hst2st2samehandcoupling(i) = 0;
     }

     hst2st2samehandcoupling = higgshsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = hst2st2samehandcoupling(1);
     couplingRR = hst2st2samehandcoupling(3);

     DoubleVector hst2st2diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     hst2st2diffhandcoupling(i) = 0;
     }

     hst2st2diffhandcoupling = higgshsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = hst2st2diffhandcoupling(1);

     coupling = couplingLL*pow(sin(theta),2) + couplingRR*pow(cos(theta),2) + 2*couplingLR*cos(theta)*sin(theta);
       
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgshamplitudedecaystop1stop2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for h->stop1 antistop2 or stop2 antistop1 as they are the same
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecay2stop1stop2\n");
    }  

     DoubleVector hst1st2samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     hst1st2samehandcoupling(i) = 0;
     }

     hst1st2samehandcoupling = higgshsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = hst1st2samehandcoupling(1);
     couplingRR = hst1st2samehandcoupling(3);

     DoubleVector hst1st2diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     hst1st2diffhandcoupling(i) = 0;
     }

     hst1st2diffhandcoupling = higgshsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = hst1st2diffhandcoupling(1);

     coupling = couplingLL*sin(theta)*cos(theta) - couplingRR*cos(theta)*sin(theta) + couplingLR*cos(2*theta);
       
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgshamplitudedecaysbottom1sbottom1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for h->sbottom1 antisbottom1
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecaysbottom1sbottom1\n");
    }  

     DoubleVector hsb1sb1samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     hsb1sb1samehandcoupling(i) = 0;
     }

     hsb1sb1samehandcoupling = higgshsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = hsb1sb1samehandcoupling(2);
     couplingRR = hsb1sb1samehandcoupling(4);

     DoubleVector hsb1sb1diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     hsb1sb1diffhandcoupling(i) = 0;
     }

     hsb1sb1diffhandcoupling = higgshsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = hsb1sb1diffhandcoupling(2);

     coupling = couplingLL*pow(cos(theta),2) + couplingRR*pow(sin(theta),2) - 2*couplingLR*cos(theta)*sin(theta);
        
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgshamplitudedecaysbottom2sbottom2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for h->sbottom2 antisbottom2
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecaysbottom2sbottom2\n");
    }  

     DoubleVector hsb2sb2samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     hsb2sb2samehandcoupling(i) = 0;
     }

     hsb2sb2samehandcoupling = higgshsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = hsb2sb2samehandcoupling(2);
     couplingRR = hsb2sb2samehandcoupling(4);

     DoubleVector hsb2sb2diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     hsb2sb2diffhandcoupling(i) = 0;
     }

     hsb2sb2diffhandcoupling = higgshsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = hsb2sb2diffhandcoupling(2);

     coupling = couplingLL*pow(sin(theta),2) + couplingRR*pow(cos(theta),2) + 2*couplingLR*cos(theta)*sin(theta);
       
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgshamplitudedecaysbottom1sbottom2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for h->sbottom1 antisbottom2 or sbottom2 antisbottom1 as they are the same
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecaysbottom1sbottom2\n");
    }  

     DoubleVector hsb1sb2samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     hsb1sb2samehandcoupling(i) = 0;
     }

     hsb1sb2samehandcoupling = higgshsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = hsb1sb2samehandcoupling(2);
     couplingRR = hsb1sb2samehandcoupling(4);

     DoubleVector hsb1sb2diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     hsb1sb2diffhandcoupling(i) = 0;
     }

     hsb1sb2diffhandcoupling = higgshsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = hsb1sb2diffhandcoupling(2);

     coupling = couplingLL*sin(theta)*cos(theta) - couplingRR*cos(theta)*sin(theta) + couplingLR*cos(2*theta);
       
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgsHamplitudedecaystop1stop1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for H->stop1 antistop1
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecaystop1stop1\n");
    }  

     DoubleVector Hst1st1samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     Hst1st1samehandcoupling(i) = 0;
     }

     Hst1st1samehandcoupling = higgsHsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = Hst1st1samehandcoupling(1);
     couplingRR = Hst1st1samehandcoupling(3);

     DoubleVector Hst1st1diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     Hst1st1diffhandcoupling(i) = 0;
     }

     Hst1st1diffhandcoupling = higgsHsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = Hst1st1diffhandcoupling(1);

     coupling = couplingLL*pow(cos(theta),2) + couplingRR*pow(sin(theta),2) - 2*couplingLR*cos(theta)*sin(theta);
        
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgsHamplitudedecaystop2stop2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for H->stop2 antistop2
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecaystop2stop2\n");
    }  

     DoubleVector Hst2st2samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     Hst2st2samehandcoupling(i) = 0;
     }

     Hst2st2samehandcoupling = higgsHsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = Hst2st2samehandcoupling(1);
     couplingRR = Hst2st2samehandcoupling(3);

     DoubleVector Hst2st2diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     Hst2st2diffhandcoupling(i) = 0;
     }

     Hst2st2diffhandcoupling = higgsHsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = Hst2st2diffhandcoupling(1);

     coupling = couplingLL*pow(sin(theta),2) + couplingRR*pow(cos(theta),2) + 2*couplingLR*cos(theta)*sin(theta);
       
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgsHamplitudedecaystop1stop2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for H->stop1 antistop2 or stop2 antistop1 as they are the same
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecaystop1stop2\n");
    }  

     DoubleVector Hst1st2samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     Hst1st2samehandcoupling(i) = 0;
     }

     Hst1st2samehandcoupling = higgsHsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = Hst1st2samehandcoupling(1);
     couplingRR = Hst1st2samehandcoupling(3);

     DoubleVector Hst1st2diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     Hst1st2diffhandcoupling(i) = 0;
     }

     Hst1st2diffhandcoupling = higgsHsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = Hst1st2diffhandcoupling(1);

     coupling = couplingLL*sin(theta)*cos(theta) - couplingRR*cos(theta)*sin(theta) + couplingLR*cos(2*theta);
       
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgsHamplitudedecaysbottom1sbottom1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for H->sbottom1 antisbottom1
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecaysbottom1sbottom1\n");
    }  

     DoubleVector Hsb1sb1samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     Hsb1sb1samehandcoupling(i) = 0;
     }

     Hsb1sb1samehandcoupling = higgsHsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = Hsb1sb1samehandcoupling(2);
     couplingRR = Hsb1sb1samehandcoupling(4);
     
     DoubleVector Hsb1sb1diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     Hsb1sb1diffhandcoupling(i) = 0;
     }

     Hsb1sb1diffhandcoupling = higgsHsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = Hsb1sb1diffhandcoupling(2);

     coupling = couplingLL*pow(cos(theta),2) + couplingRR*pow(sin(theta),2) - 2*couplingLR*cos(theta)*sin(theta);
        
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgsHamplitudedecaysbottom2sbottom2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for H->sbottom2 antisbottom2
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecaysbottom2sbottom2\n");
    }  

     DoubleVector Hsb2sb2samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     Hsb2sb2samehandcoupling(i) = 0;
     }

     Hsb2sb2samehandcoupling = higgsHsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = Hsb2sb2samehandcoupling(2);
     couplingRR = Hsb2sb2samehandcoupling(4);

     DoubleVector Hsb2sb2diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     Hsb2sb2diffhandcoupling(i) = 0;
     }

     Hsb2sb2diffhandcoupling = higgsHsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = Hsb2sb2diffhandcoupling(2);

     coupling = couplingLL*pow(sin(theta),2) + couplingRR*pow(cos(theta),2) + 2*couplingLR*cos(theta)*sin(theta);
       
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgsHamplitudedecaysbottom1sbottom2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double theta) ///calculates partial width for H->sbottom1 antisbottom2 or sbottom2 antisbottom1 as they are the same
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecaysbottom1sbottom2\n");
    }  

     DoubleVector Hsb1sb2samehandcoupling(4);
     for (int i=1; i<=4; i++) {
     Hsb1sb2samehandcoupling(i) = 0;
     }

     Hsb1sb2samehandcoupling = higgsHsquarksamehandcouplings (mWboson, g, gp, alpha, beta, mtop, mbottom);
     couplingLL = Hsb1sb2samehandcoupling(2);
     couplingRR = Hsb1sb2samehandcoupling(4);

     DoubleVector Hsb1sb2diffhandcoupling(2);
     for (int i=1; i<=2; i++) {
     Hsb1sb2diffhandcoupling(i) = 0;
     }

     Hsb1sb2diffhandcoupling = higgsHsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
     couplingLR = Hsb1sb2diffhandcoupling(2);

     coupling = couplingLL*sin(theta)*cos(theta) - couplingRR*cos(theta)*sin(theta) + couplingLR*cos(2*theta);
       
    amplitudeW = 3*pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgshamplitudedecaystau1stau1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta) ///calculates partial width for h->stau1 antistau1
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecaystau1stau1\n");
    }  

     DoubleVector hstau1stau1samehandcoupling(3);
     for (int i=1; i<=3; i++) {
     hstau1stau1samehandcoupling(i) = 0;
     }

     hstau1stau1samehandcoupling = higgshsleptonsamehandcouplings (mWboson, g, gp, alpha, beta, mtau);
     couplingLL = hstau1stau1samehandcoupling(2);
     couplingRR = hstau1stau1samehandcoupling(3);

     DoubleVector hstau1stau1diffhandcoupling(1);
     for (int i=1; i<=1; i++) {
     hstau1stau1diffhandcoupling(i) = 0;
     }

     hstau1stau1diffhandcoupling = higgshsleptondiffhandcouplings (mWboson, g, alpha, beta, mtau, greekmu, Atau);
     couplingLR = hstau1stau1diffhandcoupling(1);

     coupling = couplingLL*pow(sin(theta),2) + couplingRR*pow(cos(theta),2) + 2*couplingLR*cos(theta)*sin(theta);

     amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgshamplitudedecaystau2stau2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta) ///calculates partial width for h->stau2 antistau2
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecaystau2stau2\n");
    }  

     DoubleVector hstau2stau2samehandcoupling(3);
     for (int i=1; i<=3; i++) {
     hstau2stau2samehandcoupling(i) = 0;
     }

     hstau2stau2samehandcoupling = higgshsleptonsamehandcouplings (mWboson, g, gp, alpha, beta, mtau);
     couplingLL = hstau2stau2samehandcoupling(2);
     couplingRR = hstau2stau2samehandcoupling(3);

     DoubleVector hstau2stau2diffhandcoupling(1);
     for (int i=1; i<=1; i++) {
     hstau2stau2diffhandcoupling(i) = 0;
     }

     hstau2stau2diffhandcoupling = higgshsleptondiffhandcouplings (mWboson, g, alpha, beta, mtau, greekmu, Atau);
     couplingLR = hstau2stau2diffhandcoupling(1);

     coupling = couplingLL*pow(cos(theta),2) + couplingRR*pow(sin(theta),2) - 2*couplingLR*cos(theta)*sin(theta);

     amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}

double higgshamplitudedecaystau1stau2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta) ///calculates partial width for h->stau1 antistau2 or h->stau2 antistau1 as these widths have the same value
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecaystau1stau2\n");
    }  

     DoubleVector hstau1stau2samehandcoupling(3);
     for (int i=1; i<=3; i++) {
     hstau1stau2samehandcoupling(i) = 0;
     }

     hstau1stau2samehandcoupling = higgshsleptonsamehandcouplings (mWboson, g, gp, alpha, beta, mtau);
     couplingLL = hstau1stau2samehandcoupling(2);
     couplingRR = hstau1stau2samehandcoupling(3);

     DoubleVector hstau1stau2diffhandcoupling(1);
     for (int i=1; i<=1; i++) {
     hstau1stau2diffhandcoupling(i) = 0;
     }

     hstau1stau2diffhandcoupling = higgshsleptondiffhandcouplings (mWboson, g, alpha, beta, mtau, greekmu, Atau);
     couplingLR = hstau1stau2diffhandcoupling(1);

     coupling = -couplingLL*cos(theta)*sin(theta) + couplingRR*sin(theta)*cos(theta) + couplingLR*cos(2*theta);
     amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgsHamplitudedecaystau1stau1 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta) ///calculates partial width for H->stau1 antistau1
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecaystau1stau1\n");
    }  

     DoubleVector Hstau1stau1samehandcoupling(3);
     for (int i=1; i<=3; i++) {
     Hstau1stau1samehandcoupling(i) = 0;
     }

     Hstau1stau1samehandcoupling = higgsHsleptonsamehandcouplings (mWboson, g, gp, alpha, beta, mtau);
     couplingLL = Hstau1stau1samehandcoupling(2);
     couplingRR = Hstau1stau1samehandcoupling(3);

     DoubleVector Hstau1stau1diffhandcoupling(1);
     for (int i=1; i<=1; i++) {
     Hstau1stau1diffhandcoupling(i) = 0;
     }

     Hstau1stau1diffhandcoupling = higgsHsleptondiffhandcouplings (mWboson, g, alpha, beta, mtau, greekmu, Atau);
     couplingLR = Hstau1stau1diffhandcoupling(1);

     coupling = couplingLL*pow(sin(theta),2) + couplingRR*pow(cos(theta),2) + 2*couplingLR*cos(theta)*sin(theta);

     amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgsHamplitudedecaystau2stau2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta) ///calculates partial width for H->stau2 antistau2
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecaystau2stau2\n");
    } 

     DoubleVector Hstau2stau2samehandcoupling(3);
     for (int i=1; i<=3; i++) {
     Hstau2stau2samehandcoupling(i) = 0;
     }

     Hstau2stau2samehandcoupling = higgsHsleptonsamehandcouplings (mWboson, g, gp, alpha, beta, mtau);
     couplingLL = Hstau2stau2samehandcoupling(2);
     couplingRR = Hstau2stau2samehandcoupling(3);

     DoubleVector Hstau2stau2diffhandcoupling(1);
     for (int i=1; i<=1; i++) {
     Hstau2stau2diffhandcoupling(i) = 0;
     }

     Hstau2stau2diffhandcoupling = higgsHsleptondiffhandcouplings (mWboson, g, alpha, beta, mtau, greekmu, Atau);
     couplingLR = Hstau2stau2diffhandcoupling(1);

     coupling = couplingLL*pow(cos(theta),2) + couplingRR*pow(sin(theta),2) - 2*couplingLR*cos(theta)*sin(theta);

     amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgsHamplitudedecaystau1stau2 (double m1, double m2, double m3, double g, double gp, double alpha, double beta, double mWboson, double mtau, double greekmu, double Atau, double theta) ///calculates partial width for H->stau1 antistau2 or H->stau2 antistau1 as these widths have the same value
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, couplingLL=0, couplingRR=0, couplingLR = 0;
  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
  }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHamplitudedecaystau1stau2\n");
    } 

     DoubleVector Hstau1stau2samehandcoupling(3);
     for (int i=1; i<=3; i++) {
     Hstau1stau2samehandcoupling(i) = 0;
     }

     Hstau1stau2samehandcoupling = higgsHsleptonsamehandcouplings (mWboson, g, gp, alpha, beta, mtau);
     couplingLL = Hstau1stau2samehandcoupling(2);
     couplingRR = Hstau1stau2samehandcoupling(3);

     DoubleVector Hstau1stau2diffhandcoupling(1);
     for (int i=1; i<=1; i++) {
     Hstau1stau2diffhandcoupling(i) = 0;
     }

     Hstau1stau2diffhandcoupling = higgsHsleptondiffhandcouplings (mWboson, g, alpha, beta, mtau, greekmu, Atau);
     couplingLR = Hstau1stau2diffhandcoupling(1);

     coupling = -couplingLL*cos(theta)*sin(theta) + couplingRR*sin(theta)*cos(theta) + couplingLR*cos(2*theta);

    amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgsAamplitudedecaysfermions (double m1, double m2, double m3, double g, double mWboson, double mf, double greekmu, double Asf, double beta, char uord) ///calculates partial width for A->sfermion1 sfermion2, these sfermions must be of the same type, note can't have decays to sfermion1 and sfermion 1 or to sfermion2 and sfermion2 by CP conservation, also note there's no dependence on sfermion mixing angle in these decays, mf is the mass of the corresponding fermion. char uord tells it whether we have up type (so Asf*cot(beta)) or down type (so Asf*tan(beta)). Finally note that the extra factor of three from colour for squarks is added above by multiplying the answer this function gives by 3
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0;  

  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
      }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsAamplitudedecaysfermions\n");
    } 
    
    if(uord == 'u') {
      coupling = g*mf/(2*mWboson)*(greekmu + Asf/(tan(beta)));
    }
    else if (uord == 'd') {
      coupling = g*mf/(2*mWboson)*(greekmu + Asf*tan(beta));
    }
    else {
      throw("problem: uord must be u or d in higgsAamplitudedecaysfermions\n");
    }    
    amplitudeW = pow(coupling,2)/(16*PI*m1)*lambda;
  }
  return amplitudeW;
}


double higgsAamplitudedecaysfermionsNMSSM (double m1, double m2, double m3, double g, double mWboson, double mf, double Asf, double beta, double lam, double mueff, DoubleMatrix & CPOMix, char uord, int pseudoscalar) ///mf here should be runmf and mWboson should be runmw
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, coupling=0, fq = 0;  

  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
      }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsAamplitudedecaysfermionsNMSSM\n");
    } 

    if(uord == 'u') {
      fq = g*mf/(root2*mWboson*sin(beta));
      coupling = fq/root2*(Asf*CPOMix(pseudoscalar,1) + mueff*CPOMix(pseudoscalar,2) + lam*root2*mWboson*cos(beta)/g*CPOMix(pseudoscalar,3));
    }
    else if (uord == 'd') {
      fq = g*mf/(root2*mWboson*cos(beta));
      coupling = fq/root2*(mueff*CPOMix(pseudoscalar,1) + Asf*CPOMix(pseudoscalar,2) + lam*root2*mWboson*sin(beta)/g*CPOMix(pseudoscalar,3));
    } 
    else {
      throw("problem: uord must be u or d in higgsAamplitudedecaysfermionsNMSSM\n");
    }    
    amplitudeW = 1/(16*PI*m1)*lambda*pow(coupling,2);
  }

  return amplitudeW;
}


double higgsHplusamplitudedecayquarkantiquark (double m1, double m2, double m3, double g, double mWboson, double beta, DoubleMatrix & VCKM, int quark, int antiquark) ///calculates partial width for H+ ->quark1 antiquark2, this is the same as H- -> antiquark1 quark2, note the up type quark (quark1) is m2, the down type quark (quark2) is m3, the int quark is 1,2,3 (u,c,t) and int antiquark is 1,2,3 (d,s,b) then this selects VCKM matrix element VCKM(quark, antiquark)
{
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, massbetacombination=0, CKM=0;  

  if (fabs(m1) < fabs(m2) +fabs(m3)){
    amplitudeW = 0;
      }
  
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHplusamplitudedecayquarkantiquark\n");
    } 
    CKM = VCKM(quark,antiquark);

    massbetacombination = (pow(m3*tan(beta),2) + pow(m2/(tan(beta)),2))*(pow(m1,2) - pow(m2,2) - pow(m3,2)) - 4*pow(m2*m3,2);
    
    amplitudeW = 3*GFosqrt2*pow(CKM,2)/(4*PI*m1)*massbetacombination*lambda;

  }
  return amplitudeW;
}



double higgsHplusamplitudedecayneutralinochargino (double m1, double m2, double m3, double g, double gp, double beta, double thetaL, double thetaR, DoubleMatrix & mixNeut, int neutralino,  int chargino) /// Calculates the partial width for decays of the charged H+ higgs to a neutralino and a chragino Wtilda+ where neutralino is j in T&B whilst chargino is i
{
  
  double amplitudeW, squareplus, squareminus, lambda, A1=0, A2=0, A3=0, A4=0,a=0, b=0;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHplusamplitudedecayneutralinochargino\n");
    } 

    A1 = 1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*sin(thetaR) - g*mixNeut(neutralino,4)*cos(thetaR);
    A2 = -1/(root2)*(-g*mixNeut(neutralino,2) - gp*mixNeut(neutralino,1))*-cos(thetaR) - g*mixNeut(neutralino,4)*sin(thetaR);
    A3 = 1/(root2)*(g*mixNeut(neutralino,2) + gp*mixNeut(neutralino,1))*sin(thetaL) + g*mixNeut(neutralino,3)*cos(thetaL);
    A4 = -1/(root2)*(-g*mixNeut(neutralino,2) - gp*mixNeut(neutralino,1))*-cos(thetaL) + g*mixNeut(neutralino,3)*sin(thetaL);

    
    if (chargino == 1) {
      a = 0.5*(-cos(beta)*A2 + sin(beta)*A4);
      b = 0.5*(-cos(beta)*A2 - sin(beta)*A4);
    }
     
    else if (chargino == 2) {
      a = 0.5*(-cos(beta)*A1 + sin(beta)*A3);
      b = 0.5*(-cos(beta)*A1 - sin(beta)*A3);
    }
    else {
      throw("problem: chargino must be 1 or 2 in higgsHplusamplitudedecayneutralinochargino\n");
    } 
    amplitudeW = 1/(8*PI*fabs(m1))*lambda*((pow(a,2)+pow(b,2))*(pow(m1,2)-pow(m2,2)-pow(m3,2)) - 2*(pow(a,2)-pow(b,2))*m2*m3);
  }
  return amplitudeW;
}



double higgsHplusamplitudedecayWbosonhiggsh (double m1, double m2, double m3, double g, double alpha, double beta) /// Calculates the partial width for a charged Higgs boson H+ to decay to a Wboson and a neutral light scalar higgs (i.e. the SM-like higgs) h, m2 must be the W boson mass
{
  
  double amplitudeW, squareplus, squareminus, lambda;
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    amplitudeW = 0;
  }

  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHplusamplitudedecayWbosonhiggsh\n");
    } 
    
    // amplitudeW = pow(g,2)*pow(cos(beta-alpha),2)*pow(m1,3)/(64*PI*pow(m2,2))*pow(lambda,3);
    amplitudeW = GFosqrt2*pow(cos(beta-alpha),2)*pow(m1,3)/(8*PI)*pow(lambda,3);
  }
  return amplitudeW;
}


///calculates partial width for
///Hplus->squarki antisquarkj where i,j are each L/R so no mixing, i.e. first
///two generations of squarks 
DoubleVector higgsHplusamplitudedecaysquarksquark (double m1, double m2, double m3, double g, double beta, double mWboson, double mup, double mdown, double greekmu, double Aup, double Adown) {
  double squareplus=0, squareminus=0, lambda=0;
  DoubleVector amplitudeW(4);
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    for(int i=1; i<=4; i++) {
      amplitudeW(i) = 0;
    }
  }  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHplusamplitudedecaysquarksquark\n");
    } 

     DoubleVector Hplussquarksquarkcoupling(4);
     for (int i=1; i<=4; i++) {
     Hplussquarksquarkcoupling(i) = 0;
     }

     Hplussquarksquarkcoupling = higgsHplussquarkcouplings (mWboson, g, beta, mup, mdown, greekmu, Aup, Adown);
     
     amplitudeW(1) = pow(Hplussquarksquarkcoupling(1),2)*3/(16*PI*m1)*lambda;
     amplitudeW(2) = pow(Hplussquarksquarkcoupling(2),2)*3/(16*PI*m1)*lambda;
     amplitudeW(3) = pow(Hplussquarksquarkcoupling(3),2)*3/(16*PI*m1)*lambda;
     amplitudeW(4) = pow(Hplussquarksquarkcoupling(4),2)*3/(16*PI*m1)*lambda;

  }
  return amplitudeW;
}


DoubleVector higgsHplusamplitudedecaysquarksquarkmix (double m1, double m2, double m3, double g, double beta, double mWboson, double mtop, double mbottom, double greekmu, double Atop, double Abottom, double thetatop, double thetabottom) ///calculates partial width for Hplus->squarki antisquarkj where i,j are each 1/2 so mixing included, i.e. third generations of squarks
{
  double squareplus=0, squareminus=0, lambda=0, coupling11=0, coupling22=0, coupling12=0, coupling21=0; 
  DoubleVector amplitudeW(4);
  for (int i=1; i<=4; i++) {
    amplitudeW(i) = 0;
  }    
  
  if (fabs(m1) < fabs(m2) +fabs(m3)) {
    for(int i=1; i<=4; i++) {
      amplitudeW(i) = 0;
    }
  }
  else { 
    squareplus = 1 - pow((m2+m3)/m1,2);
    squareminus = 1 - pow((m2-m3)/m1,2);
    lambda = pow(squareplus*squareminus,0.5);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsHplusamplitudedecaysquarksquarkmix\n");
    } 

     DoubleVector Hplussquarksquarkcoupling(4);
     for (int i=1; i<=4; i++) {
     Hplussquarksquarkcoupling(i) = 0;
     }

     Hplussquarksquarkcoupling = higgsHplussquarkcouplings (mWboson, g, beta, mtop, mbottom, greekmu, Atop, Abottom);
     
     coupling11 = Hplussquarksquarkcoupling(1)*cos(thetatop)*cos(thetabottom) + Hplussquarksquarkcoupling(2)*sin(thetatop)*sin(thetabottom) - Hplussquarksquarkcoupling(3)*cos(thetatop)*sin(thetabottom) - Hplussquarksquarkcoupling(4)*sin(thetatop)*cos(thetabottom);
     coupling22 = Hplussquarksquarkcoupling(1)*sin(thetatop)*sin(thetabottom) + Hplussquarksquarkcoupling(2)*cos(thetatop)*cos(thetabottom) + Hplussquarksquarkcoupling(3)*sin(thetatop)*cos(thetabottom) + Hplussquarksquarkcoupling(4)*cos(thetatop)*sin(thetabottom);
     coupling12 = Hplussquarksquarkcoupling(1)*cos(thetatop)*sin(thetabottom) - Hplussquarksquarkcoupling(2)*sin(thetatop)*cos(thetabottom) + Hplussquarksquarkcoupling(3)*cos(thetatop)*cos(thetabottom) - Hplussquarksquarkcoupling(4)*sin(thetatop)*sin(thetabottom);
     coupling21 = Hplussquarksquarkcoupling(1)*sin(thetatop)*cos(thetabottom) - Hplussquarksquarkcoupling(2)*cos(thetatop)*sin(thetabottom) - Hplussquarksquarkcoupling(3)*sin(thetatop)*sin(thetabottom) + Hplussquarksquarkcoupling(4)*cos(thetatop)*cos(thetabottom);

     amplitudeW(1) = pow(coupling11,2)*3/(16*PI*m1)*lambda;
     amplitudeW(2) = pow(coupling22,2)*3/(16*PI*m1)*lambda;
     amplitudeW(3) = pow(coupling12,2)*3/(16*PI*m1)*lambda;
     amplitudeW(4) = pow(coupling21,2)*3/(16*PI*m1)*lambda;

  }
  return amplitudeW;
}

double higgsesamplitudedecaygammagammatotal(double m1, double g, double gprime, double alphaEmrun, double mWboson, double polemw, double alpha, double beta, double mtop, double mbottom, double mcharm, double mtau, double mHpm, double mstop1, double mstop2, double msbottom1, double msbottom2, double mstau1, double mstau2, double mchargino1, double mchargino2, double thetaL, double thetaR, double thetat, double thetab, double thetatau, double greekmu, double Atop, double Abottom, double Atau, char higgstype) /// function that calculates the partial width for h->gamma gamma at 1-loop (forbidden at tree-level)
{
  double prefactor=0, Itr=0, Iti=0, Ist1r=0, Ist1i=0, Ist2r=0, Ist2i=0, Ibr=0, Ibi=0, Isb1r=0, Isb1i=0, Isb2r=0, Isb2i=0, Icr=0, Ici=0, Itaur=0, Itaui=0, Istau1r=0, Istau1i=0, Istau2r=0, Istau2i=0, IWr=0, IWi=0, IHpmr=0, IHpmi=0, Ichar1r=0, Ichar1i=0, Ichar2r=0, Ichar2i=0, matelemmodsquare=0, amplitudeW=0;

  Itr = higgsmatrixelementgammagammaviatops (m1, mtop, alpha, beta, higgstype)(1);
  Iti = higgsmatrixelementgammagammaviatops (m1, mtop, alpha, beta, higgstype)(2);
  Ist1r = -higgsmatrixelementgammagammaviastops (m1, mstop1, mstop2, mtop, mbottom, mWboson, thetat, g, gprime, alpha, beta, greekmu, Atop, Abottom, higgstype)(1);
  Ist1i = -higgsmatrixelementgammagammaviastops (m1, mstop1, mstop2, mtop, mbottom, mWboson, thetat, g, gprime, alpha, beta, greekmu, Atop, Abottom, higgstype)(2);
  Ist2r = -higgsmatrixelementgammagammaviastops (m1, mstop1, mstop2, mtop, mbottom, mWboson, thetat, g, gprime, alpha, beta, greekmu, Atop, Abottom, higgstype)(3);
  Ist2i = -higgsmatrixelementgammagammaviastops (m1, mstop1, mstop2, mtop, mbottom, mWboson, thetat, g, gprime, alpha, beta, greekmu, Atop, Abottom, higgstype)(4);
  Ibr = higgsmatrixelementgammagammaviabottoms (m1, mbottom, alpha, beta, higgstype)(1);
  Ibi = higgsmatrixelementgammagammaviabottoms (m1, mbottom, alpha, beta, higgstype)(2);
  Isb1r = -higgsmatrixelementgammagammaviasbottoms (m1, msbottom1, msbottom2, mbottom, mtop, mWboson, thetab, g, gprime, alpha, beta, Atop, Abottom, greekmu, higgstype)(1);
  Isb1i = -higgsmatrixelementgammagammaviasbottoms (m1, msbottom1, msbottom2, mbottom, mtop, mWboson, thetab, g, gprime, alpha, beta, Atop, Abottom, greekmu, higgstype)(2);
  Isb2r = -higgsmatrixelementgammagammaviasbottoms (m1, msbottom1, msbottom2, mbottom, mtop, mWboson, thetab, g, gprime, alpha, beta, Atop, Abottom, greekmu, higgstype)(3);
  Isb2i = -higgsmatrixelementgammagammaviasbottoms (m1, msbottom1, msbottom2, mbottom, mtop, mWboson, thetab, g, gprime, alpha, beta, Atop, Abottom, greekmu, higgstype)(4);
  Icr = higgsmatrixelementgammagammaviacharms (m1, mcharm, alpha, beta, higgstype)(1);
  Ici = higgsmatrixelementgammagammaviacharms (m1, mcharm, alpha, beta, higgstype)(2);
  Itaur = higgsmatrixelementgammagammaviataus (m1, mtau, alpha, beta, higgstype)(1);
  Itaui = higgsmatrixelementgammagammaviataus (m1, mtau, alpha, beta, higgstype)(2);
  Istau1r = -higgsmatrixelementgammagammaviastaus (m1, mstau1, mstau2, mtau, mWboson, thetatau, g, gprime, alpha, beta, greekmu, Atau, higgstype)(1);
  Istau1i = -higgsmatrixelementgammagammaviastaus (m1, mstau1, mstau2, mtau, mWboson, thetatau, g, gprime, alpha, beta, greekmu, Atau, higgstype)(2);
  Istau2r = -higgsmatrixelementgammagammaviastaus (m1, mstau1, mstau2, mtau, mWboson, thetatau, g, gprime, alpha, beta, greekmu, Atau, higgstype)(3);
  Istau2i = -higgsmatrixelementgammagammaviastaus (m1, mstau1, mstau2, mtau, mWboson, thetatau, g, gprime, alpha, beta, greekmu, Atau, higgstype)(4);
  IWr = higgsmatrixelementgammagammaviaWbosons (m1, polemw, alpha, beta, g, gprime, higgstype)(1);
  IWi = higgsmatrixelementgammagammaviaWbosons (m1, polemw, alpha, beta, g, gprime, higgstype)(2);
  IHpmr = higgsmatrixelementgammagammaviaHpms (m1, mHpm, mWboson, alpha, beta, g, gprime, higgstype)(1);
  IHpmi = higgsmatrixelementgammagammaviaHpms (m1, mHpm, mWboson, alpha, beta, g, gprime, higgstype)(2);
  Ichar1r = higgsmatrixelementgammagammaviachargino1s (m1, mchargino1, mWboson, alpha, beta, thetaL, thetaR, higgstype)(1);
  Ichar1i = higgsmatrixelementgammagammaviachargino1s (m1, mchargino1, mWboson, alpha, beta, thetaL, thetaR, higgstype)(2);
  Ichar2r = higgsmatrixelementgammagammaviachargino2s (m1, mchargino2, mWboson, alpha, beta, thetaL, thetaR, higgstype)(1);
  Ichar2i = higgsmatrixelementgammagammaviachargino2s (m1, mchargino2, mWboson, alpha, beta, thetaL, thetaR, higgstype)(2);

  DoubleVector matelemsum(2);
  matelemsum(1) = Itr + Ist1r + Ist2r + Ibr + Isb1r + Isb2r + Icr + Itaur + Istau1r + Istau2r + IWr + IHpmr + Ichar1r + Ichar2r;
  matelemsum(2) = Iti + Ist1i + Ist2i + Ibi + Isb1i + Isb2i + Ici + Itaui + Istau1i + Istau2i + IWi + IHpmi + Ichar1i + Ichar2i;
  
  prefactor = GFosqrt2*pow(alphaEmrun,2)*pow(m1,3)/(128*pow(PI,3));
  matelemmodsquare = pow(matelemsum(1),2) + pow(matelemsum(2),2);
  if (matelemmodsquare != matelemmodsquare) {
    throw("nan in matelemmodsquare in higgsesamplitudedecaygammagammatotal\n");
  }
  amplitudeW = prefactor*matelemmodsquare;
     
  return amplitudeW;
}

///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for tops in the loop.
DoubleVector higgsmatrixelementgammagammaviatops (double m1, double mtop, double alpha, double beta, char higgstype) {
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f(3), F1over2(2), It(2);
  double Rt=0;
  for (int i=1; i<=3; i++) {
    f(i) = 0;
  }

  for (int j=1; j<=2; j++) {
    F1over2(j) = 0;
    It(j) = 0;
  }
    
  f = foftau(mtop, m1);
  if (higgstype == 'h')
    {
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
      Rt = cos(alpha)/sin(beta);
    }
  else if (higgstype == 'H')
    {
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
      Rt = sin(alpha)/sin(beta);
    }
  else if (higgstype == 'A')
    {
      F1over2(1) = -2*f(3)*f(1);
      F1over2(2) = -2*f(3)*f(2);
      Rt = 1/(tan(beta));
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviatops\n");
  }
   
  It(1) = 3*4./9*Rt*F1over2(1);
  It(2) = 3*4./9*Rt*F1over2(2);
    return It;
}

DoubleVector higgsmatrixelementgammagammaviastops (double m1, double mstop1, double mstop2, double mtop, double mbottom, double mWboson, double thetat, double g, double gprime, double alpha, double beta,double greekmu, double Atop, double Abottom, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for stop1s and for stop2s in the loop.
{
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f1(3), f2(3), F01(2), F02(2), Ist1(2), Ist2(2), Istopsboth(4);
  double Rst1=0, Rst2=0, RstL=0, RstR=0, RstLR=0, RstL1=0, RstL2=0, RstR1=0, RstR2=0, RstLR1=0, RstLR2=0;

  for (int i=1; i<=3;i++) {
    f1(i) = 0;
    f2(i) = 0;
  }
  for (int j=1; j<=2;j++) {
    F01(j) = 0;
    F02(j) = 0;
    Ist1(j) = 0;
    Ist2(j) = 0;
  }
  for (int k=1; k<=4; k++) {
    Istopsboth(k) = 0;
  }

  DoubleVector ststsamehandcoupling(4);
  for (int i=1; i<=4; i++) {
  ststsamehandcoupling(i) = 0;
  }

  DoubleVector ststdiffhandcoupling(2);
  for (int i=1; i<=2; i++) {
  ststdiffhandcoupling(i) = 0;
  }
  
  f1 = foftau(mstop1, m1);
  f2 = foftau(mstop2, m1);
  if (higgstype == 'h')
    {
      ststsamehandcoupling = higgshsquarksamehandcouplings (mWboson, g, gprime, alpha, beta, mtop, mbottom);
      RstL = ststsamehandcoupling(1);
      RstL1 = RstL*mWboson/(g*pow(mstop1,2));
      RstL2 = RstL*mWboson/(g*pow(mstop2,2));
      RstR = ststsamehandcoupling(3);
      RstR1 = RstR*mWboson/(g*pow(mstop1,2));
      RstR2 = RstR*mWboson/(g*pow(mstop2,2));
      ststdiffhandcoupling = higgshsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
      RstLR = ststdiffhandcoupling(1);
      RstLR1 = RstLR*mWboson/(g*pow(mstop1,2));
      RstLR2 = RstLR*mWboson/(g*pow(mstop2,2));
      Rst1 = RstL1*pow(cos(thetat),2) + RstR1*pow(sin(thetat),2) - 2*RstLR1*cos(thetat)*sin(thetat);
      Rst2 = RstL2*pow(sin(thetat),2) + RstR2*pow(cos(thetat),2) + 2*RstLR2*cos(thetat)*sin(thetat);
    }
  else if (higgstype == 'H')
    {
      ststsamehandcoupling = higgsHsquarksamehandcouplings (mWboson, g, gprime, alpha, beta, mtop, mbottom);
      RstL = ststsamehandcoupling(1);
      RstL1 = RstL*mWboson/(g*pow(mstop1,2));
      RstL2 = RstL*mWboson/(g*pow(mstop2,2));
      RstR = ststsamehandcoupling(3);
      RstR1 = RstR*mWboson/(g*pow(mstop1,2));
      RstR2 = RstR*mWboson/(g*pow(mstop2,2));
      ststdiffhandcoupling = higgsHsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
      RstLR = ststdiffhandcoupling(1);
      RstLR1 = RstLR*mWboson/(g*pow(mstop1,2));
      RstLR2 = RstLR*mWboson/(g*pow(mstop2,2));
      Rst1 = RstL1*pow(cos(thetat),2) + RstR1*pow(sin(thetat),2) - 2*RstLR1*cos(thetat)*sin(thetat);
      Rst2 = RstL2*pow(sin(thetat),2) + RstR2*pow(cos(thetat),2) + 2*RstLR2*cos(thetat)*sin(thetat);
    }
  else if (higgstype == 'A')
    {
      RstL=0;
      RstR=0;
      Rst1 = cos(thetat)*RstL + sin(thetat)*RstR;
      Rst2 = -sin(thetat)*RstL + cos(thetat)*RstR;
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviastops\n");
  }
  
  F01(1) = f1(3)*(1-f1(3)*f1(1));
  F01(2) = f1(3)*(-f1(3)*f1(2));
  F02(1) = f2(3)*(1-f2(3)*f2(1));
  F02(2) = f2(3)*(-f2(3)*f2(2));

  Ist1(1) = 3*4./9*Rst1*F01(1);
  Ist1(2) = 3*4./9*Rst1*F01(2);
  Ist2(1) = 3*4./9*Rst2*F02(1);
  Ist2(2) = 3*4./9*Rst2*F02(2);
  Istopsboth(1) = Ist1(1);
  Istopsboth(2) = Ist1(2);
  Istopsboth(3) = Ist2(1);
  Istopsboth(4) = Ist2(2);
  return Istopsboth;
}  



DoubleVector higgsmatrixelementgammagammaviabottoms (double m1, double mbottom, double alpha, double beta, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for bottoms in the loop.
{
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f(3), F1over2(2), Ib(2);
  double Rb=0;
  for (int i=1; i<=3;i++) {
    f(i) = 0;
  }
  for (int j=1; j<=2;j++) {
    F1over2(j) = 0;
    Ib(j) = 0;
  }
    
  f = foftau(mbottom, m1);
  if (higgstype == 'h')
    {
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
      Rb = -sin(alpha)/cos(beta);
    }
  else if (higgstype == 'H')
    {
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
      Rb = cos(alpha)/cos(beta);
    }
  else if (higgstype == 'A')
    {
      F1over2(1) = -2*f(3)*f(1);
      F1over2(2) = -2*f(3)*f(2);
      Rb = tan(beta);
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviabottoms\n");
  }   
  Ib(1) = 3*1./9*Rb*F1over2(1);
  Ib(2) = 3*1./9*Rb*F1over2(2);
  return Ib;
}  


 DoubleVector higgsmatrixelementgammagammaviasbottoms (double m1, double msbottom1, double msbottom2, double mbottom, double mtop, double mWboson, double thetab, double g, double gprime, double alpha, double beta, double Atop, double Abottom, double greekmu, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for sbottom1s and for sbottom2s in the loop.
{
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f1(3), f2(3), F01(2), F02(2), Isb1(2), Isb2(2), Isbottomsboth(4);
  double Rsb1=0, Rsb2=0, RsbL=0, RsbR=0, RsbLR=0, RsbL1=0, RsbL2=0, RsbR1=0, RsbR2=0, RsbLR1=0, RsbLR2=0;
  for (int i=1; i<=3;i++) {
    f1(i) = 0;
    f2(i) = 0;
  }
  for (int j=1; j<=2;j++) {
    F01(j) = 0;
    F02(j) = 0;
    Isb1(j) = 0;
    Isb2(j) = 0;
  }
  for (int k=1; k<=4; k++) {
    Isbottomsboth(k) = 0;
  }

  DoubleVector sbsbsamehandcoupling(4);
  for (int i=1; i<=4; i++) {
  sbsbsamehandcoupling(i) = 0;
  }

  DoubleVector sbsbdiffhandcoupling(2);
  for (int i=1; i<=2; i++) {
  sbsbdiffhandcoupling(i) = 0;
  }
  

  f1 = foftau(msbottom1, m1);
  f2 = foftau(msbottom2, m1);

  if (higgstype == 'h')
    {
      sbsbsamehandcoupling = higgshsquarksamehandcouplings (mWboson, g, gprime, alpha, beta, mtop, mbottom);
      RsbL = sbsbsamehandcoupling(2);
      RsbL1 = RsbL*mWboson/(g*pow(msbottom1,2));
      RsbL2 = RsbL*mWboson/(g*pow(msbottom2,2));
      RsbR = sbsbsamehandcoupling(4);
      RsbR1 = RsbR*mWboson/(g*pow(msbottom1,2));
      RsbR2 = RsbR*mWboson/(g*pow(msbottom2,2));
      sbsbdiffhandcoupling = higgshsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
      RsbLR = sbsbdiffhandcoupling(2);
      RsbLR1 = RsbLR*mWboson/(g*pow(msbottom1,2));
      RsbLR2 = RsbLR*mWboson/(g*pow(msbottom2,2));
      Rsb1 = RsbL1*pow(cos(thetab),2) + RsbR1*pow(sin(thetab),2) - 2*RsbLR1*cos(thetab)*sin(thetab);
      Rsb2 = RsbL2*pow(sin(thetab),2) + RsbR2*pow(cos(thetab),2) + 2*RsbLR2*cos(thetab)*sin(thetab);
    }
  else if (higgstype == 'H')
    {
      sbsbsamehandcoupling = higgsHsquarksamehandcouplings (mWboson, g, gprime, alpha, beta, mtop, mbottom);
      RsbL = sbsbsamehandcoupling(2);
      RsbL1 = RsbL*mWboson/(g*pow(msbottom1,2));
      RsbL2 = RsbL*mWboson/(g*pow(msbottom2,2));
      RsbR = sbsbsamehandcoupling(4);
      RsbR1 = RsbR*mWboson/(g*pow(msbottom1,2));
      RsbR2 = RsbR*mWboson/(g*pow(msbottom2,2));
      sbsbdiffhandcoupling = higgsHsquarkdiffhandcouplings (mWboson, g, alpha, beta, mtop, mbottom, greekmu, Atop, Abottom);
      RsbLR = sbsbdiffhandcoupling(2);
      RsbLR1 = RsbLR*mWboson/(g*pow(msbottom1,2));
      RsbLR2 = RsbLR*mWboson/(g*pow(msbottom2,2));
      Rsb1 = RsbL1*pow(cos(thetab),2) + RsbR1*pow(sin(thetab),2) - 2*RsbLR1*cos(thetab)*sin(thetab);
      Rsb2 = RsbL2*pow(sin(thetab),2) + RsbR2*pow(cos(thetab),2) + 2*RsbLR2*cos(thetab)*sin(thetab);
    }
  else if (higgstype == 'A')
    {
      RsbL=0;
      RsbR=0;
      Rsb1 = cos(thetab)*RsbL + sin(thetab)*RsbR;
      Rsb2 = -sin(thetab)*RsbL + cos(thetab)*RsbR;
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviasbottoms\n");
  }  
  F01(1) = f1(3)*(1-f1(3)*f1(1));
  F01(2) = f1(3)*(-f1(3)*f1(2));
  F02(1) = f2(3)*(1-f2(3)*f2(1));
  F02(2) = f2(3)*(-f2(3)*f2(2));

  Isb1(1) = 3*1./9*Rsb1*F01(1);
  Isb1(2) = 3*1./9*Rsb1*F01(2);
  Isb2(1) = 3*1./9*Rsb2*F02(1);
  Isb2(2) = 3*1./9*Rsb2*F02(2);
  Isbottomsboth(1) = Isb1(1);
  Isbottomsboth(2) = Isb1(2);
  Isbottomsboth(3) = Isb2(1);
  Isbottomsboth(4) = Isb2(2);
  return Isbottomsboth;
}  

DoubleVector higgsmatrixelementgammagammaviacharms (double m1, double mcharm, double alpha, double beta, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for charms in the loop.
{
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f(3), F1over2(2), Ic(2);
  double Rc=0;
  for (int i=1; i<=3; i++) {
    f(i) = 0;
  }

  for (int j=1; j<=2; j++) {
    F1over2(j) = 0;
    Ic(j) = 0;
  }
    
  f = foftau(mcharm, m1);
  if (higgstype == 'h')
    {
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
      Rc = cos(alpha)/sin(beta);
    }
  else if (higgstype == 'H')
    {
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
      Rc = sin(alpha)/sin(beta);
    }
  else if (higgstype == 'A')
    {
      F1over2(1) = -2*f(3)*f(1);
      F1over2(2) = -2*f(3)*f(2);
      Rc = 1/(tan(beta));
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviacharms\n");
  }
  Ic(1) = 3*4./9*Rc*F1over2(1);
  Ic(2) = 3*4./9*Rc*F1over2(2);
  return Ic;
}

DoubleVector higgsmatrixelementgammagammaviataus (double m1, double mtau, double alpha, double beta, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for taus in the loop.
{
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f(3), F1over2(2), Itau(2);
  double Rtau=0;
  for (int i=1; i<=3;i++) {
    f(i) = 0;
  }
  for (int j=1; j<=2;j++) {
    F1over2(j) = 0;
    Itau(j) = 0;
  }
    
  f = foftau(mtau, m1);
  if (higgstype == 'h')
    {
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
      Rtau = -sin(alpha)/cos(beta);
    }
  else if (higgstype == 'H')
    {
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
      Rtau = cos(alpha)/cos(beta);
    }
  else if (higgstype == 'A')
    {
      F1over2(1) = -2*f(3)*f(1);
      F1over2(2) = -2*f(3)*f(2);
      Rtau = tan(beta);
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviataus\n");
  }  
  Itau(1) = 1*Rtau*F1over2(1);
  Itau(2) = 1*Rtau*F1over2(2);
  return Itau;
}  

DoubleVector higgsmatrixelementgammagammaviastaus (double m1, double mstau1, double mstau2, double mtau, double mWboson, double thetatau, double g, double gprime, double alpha, double beta, double greekmu, double Atau, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for stau1s and for stau2s in the loop.
{
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f1(3), f2(3), F01(2), F02(2), Istau1(2), Istau2(2), Istausboth(4);
  double Rstau1=0, Rstau2=0, RstauL=0, RstauR=0, RstauLR=0, RstauL1=0, RstauL2=0, RstauR1=0, RstauR2=0, RstauLR1=0, RstauLR2=0;
  
  for (int i=1; i<=3;i++) {
    f1(i) = 0;
    f2(i) = 0;
  }
  for (int j=1; j<=2;j++) {
    F01(j) = 0;
    F02(j) = 0;
    Istau1(j) = 0;
    Istau2(j) = 0;
  }
  for (int k=1; k<=4; k++) {
    Istausboth(k) = 0;
  }
  
  DoubleVector staustausamehandcoupling(3);
  for (int i=1; i<=3; i++) {
  staustausamehandcoupling(i) = 0;
  }

  DoubleVector staustaudiffhandcoupling(1);
  for (int i=1; i<=1; i++) {
  staustaudiffhandcoupling(i) = 0;
  }

  f1 = foftau(mstau1, m1);
  f2 = foftau(mstau2, m1);
  if (higgstype == 'h')
    {
      staustausamehandcoupling = higgshsleptonsamehandcouplings (mWboson, g, gprime, alpha, beta, mtau);
      RstauL = staustausamehandcoupling(2);
      RstauL1 = RstauL*mWboson/(g*pow(mstau1,2));
      RstauL2 = RstauL*mWboson/(g*pow(mstau2,2));
      RstauR = staustausamehandcoupling(3);
      RstauR1 = RstauR*mWboson/(g*pow(mstau1,2));
      RstauR2 = RstauR*mWboson/(g*pow(mstau2,2));
      staustaudiffhandcoupling = higgshsleptondiffhandcouplings (mWboson, g, alpha, beta, mtau, greekmu, Atau);
      RstauLR = staustaudiffhandcoupling(1);
      RstauLR1 = RstauLR*mWboson/(g*pow(mstau1,2));
      RstauLR2 = RstauLR*mWboson/(g*pow(mstau2,2));
      Rstau1 = RstauL1*pow(sin(thetatau),2) + RstauR1*pow(cos(thetatau),2) + 2*RstauLR1*cos(thetatau)*sin(thetatau);
      Rstau2 = RstauL2*pow(cos(thetatau),2) + RstauR2*pow(sin(thetatau),2) - 2*RstauLR2*cos(thetatau)*sin(thetatau);
    }
  else if (higgstype == 'H')
    {
      staustausamehandcoupling = higgsHsleptonsamehandcouplings (mWboson, g, gprime, alpha, beta, mtau);
      RstauL = staustausamehandcoupling(2);
      RstauL1 = RstauL*mWboson/(g*pow(mstau1,2));
      RstauL2 = RstauL*mWboson/(g*pow(mstau2,2));
      RstauR = staustausamehandcoupling(3);
      RstauR1 = RstauR*mWboson/(g*pow(mstau1,2));
      RstauR2 = RstauR*mWboson/(g*pow(mstau2,2));
      staustaudiffhandcoupling = higgsHsleptondiffhandcouplings (mWboson, g, alpha, beta, mtau, greekmu, Atau);
      RstauLR = staustaudiffhandcoupling(1);
      RstauLR1 = RstauLR*mWboson/(g*pow(mstau1,2));
      RstauLR2 = RstauLR*mWboson/(g*pow(mstau2,2));
      Rstau1 = RstauL1*pow(sin(thetatau),2) + RstauR1*pow(cos(thetatau),2) + 2*RstauLR1*cos(thetatau)*sin(thetatau);
      Rstau2 = RstauL2*pow(cos(thetatau),2) + RstauR2*pow(sin(thetatau),2) - 2*RstauLR2*cos(thetatau)*sin(thetatau);
    }
  else if (higgstype == 'A')
    {
      RstauL=0;
      RstauR=0;
      Rstau1 = sin(thetatau)*RstauL - cos(thetatau)*RstauR;
      Rstau2 = cos(thetatau)*RstauL + sin(thetatau)*RstauR;
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviastaus\n");
  }
  
  F01(1) = f1(3)*(1-f1(3)*f1(1));
  F01(2) = f1(3)*(-f1(3)*f1(2));
  F02(1) = f2(3)*(1-f2(3)*f2(1));
  F02(2) = f2(3)*(-f2(3)*f2(2));

  Istau1(1) = 1*Rstau1*F01(1);
  Istau1(2) = 1*Rstau1*F01(2);
  Istau2(1) = 1*Rstau2*F02(1);
  Istau2(2) = 1*Rstau2*F02(2);
  Istausboth(1) = Istau1(1);
  Istausboth(2) = Istau1(2);
  Istausboth(3) = Istau2(1);
  Istausboth(4) = Istau2(2);
  return Istausboth;
} 

DoubleVector higgsmatrixelementgammagammaviaWbosons (double m1, double mWboson, double alpha, double beta, double g, double gprime, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for W bosons in the loop.
{
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f(3), F1(2), IW(2);
  double RW=0;
  for (int i=1; i<=3;i++) {
    f(i) = 0;
  }
  for (int j=1; j<=2;j++) {
    F1(j) = 0;
    IW(j) = 0;
  }
    
  f = foftau(mWboson, m1);
  if (higgstype == 'h')
    {
      RW = sin(beta-alpha);
    }
  else if (higgstype == 'H')
    {
      RW = cos(beta-alpha);
    }
  else if (higgstype == 'A')
    {
      RW = 0;
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviaWbosons\n");
  }
  F1(1) = 2 + 3*f(3) + 3*f(3)*(2-f(3))*f(1);
  F1(2) = 3*f(3)*(2-f(3))*f(2);
  
  IW(1) = RW*F1(1);
  IW(2) = RW*F1(2);
  return IW;
}  

DoubleVector higgsmatrixelementgammagammaviaHpms (double m1, double mHpm, double mWboson, double alpha, double beta, double g, double gprime, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for H+-s in the loop.
{
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f(3), F0(2), IHpm(2);
  double RHpm=0, costhetaW=0;
  for (int i=1; i<=3;i++) {
    f(i) = 0;
  }
  for (int j=1; j<=2;j++) {
    F0(j) = 0;
    IHpm(j) = 0;
  }
    
  f = foftau(mHpm, m1);

  costhetaW = pow(pow(g,2)/(pow(g,2)+pow(gprime,2)),0.5);
  if (higgstype == 'h')
    {
      RHpm = sin(beta-alpha) + cos(2*beta)*sin(beta+alpha)/(2*pow(costhetaW,2));
    }
  else if (higgstype == 'H')
    {
      RHpm = cos(beta-alpha) - cos(2*beta)*cos(beta+alpha)/(2*pow(costhetaW,2));
    }
  else if (higgstype == 'A')
    {
      RHpm = 0;
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviaHpms\n");
  }

  F0(1) = f(3)*(1-f(3)*f(1));
  F0(2) = f(3)*(f(3)*f(2));
  
  IHpm(1) = RHpm*F0(1)*pow(mWboson/mHpm,2);
  IHpm(2) = RHpm*F0(2)*pow(mWboson/mHpm,2);
  return IHpm;

}


DoubleVector higgsmatrixelementgammagammaviachargino1s (double m1, double mchargino1, double mWboson, double alpha, double beta, double thetaL, double thetaR, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for chargino1s in the loop.
{
  DoubleVector higgsphisamecharginocouplings (double alpha, double beta, double thetaL, double thetaR);
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f(3), F1over2(2), Ichargino1(2), Charcouplings(6);
  double Rchargino1=0;
  for (int i=1; i<=3;i++) {
    f(i) = 0;
  }
  for (int j=1; j<=2;j++) {
    F1over2(j) = 0;
    Ichargino1(j) = 0;
  }
  for (int k=1; k<=6;k++) {
    Charcouplings(k) = 0;
  }
  double thetaL2 = -thetaL + PI/2;
  double thetaR2 = -thetaR + PI/2;
    
  f = foftau(mchargino1, m1);
  Charcouplings = higgsphisamecharginocouplings (alpha, beta, thetaL2, thetaR2);

  if (higgstype == 'h')
    {
      Rchargino1 = root2*((cos(thetaL2)*sin(thetaR2))*-sin(alpha) + (sin(thetaL2)*cos(thetaR2))*cos(alpha));
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));

    }
  else if (higgstype == 'H')
    {
      Rchargino1 = root2*((cos(thetaL2)*sin(thetaR2))*cos(alpha) + (sin(thetaL2)*cos(thetaR2))*sin(alpha));
      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
    }
  else if (higgstype == 'A')
    {
      Rchargino1 = -root2*(-(cos(thetaL2)*sin(thetaR2))*sin(beta) - (sin(thetaL2)*cos(thetaR2))*cos(beta));
      F1over2(1) = -2*f(3)*f(1);
      F1over2(2) = -2*f(3)*f(2);
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviachargino1s\n");
  }
  
  Ichargino1(1) = Rchargino1*F1over2(1)*mWboson/mchargino1;
  Ichargino1(2) = Rchargino1*F1over2(2)*mWboson/mchargino1;

  return Ichargino1;

}

DoubleVector higgsmatrixelementgammagammaviachargino2s (double m1, double mchargino2, double mWboson, double alpha, double beta, double thetaL, double thetaR, char higgstype) ///function that calculates the part of the matrix element that differs depending on the loop particles, in this case it calculates it for chargino2s in the loop.
{
  DoubleVector foftau(double mpart, double mcomp);
  DoubleVector f(3), F1over2(2), Ichargino2(2), Charcouplings(6);
  double Rchargino2=0;
  for (int i=1; i<=3;i++) {
    f(i) = 0;
  }
  for (int j=1; j<=2;j++) {
    F1over2(j) = 0;
    Ichargino2(j) = 0;
  }
  for (int k=1; k<=6;k++) {
    Charcouplings(k) = 0;
  }
  double thetaL2 = -thetaL + PI/2;
  double thetaR2 = -thetaR + PI/2;
    
  f = foftau(mchargino2, m1);

  Charcouplings = higgsphisamecharginocouplings (alpha, beta, thetaL2, thetaR2);


  if (higgstype == 'h')
    {
      Rchargino2 = root2*(-sin(thetaL2)*cos(thetaR2)*-sin(alpha) + (-cos(thetaL2)*sin(thetaR2))*cos(alpha));

      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));

    }
  else if (higgstype == 'H')
    {
      Rchargino2 = root2*(-sin(thetaL2)*cos(thetaR2)*cos(alpha) + (-cos(thetaL2)*sin(thetaR2))*sin(alpha));

      F1over2(1) = -2*f(3)*(1+(1-f(3))*f(1));
      F1over2(2) = -2*f(3)*((1-f(3))*f(2));
    }
  else if (higgstype == 'A')
    {
      Rchargino2 = -root2*(sin(thetaL2)*cos(thetaR2)*sin(beta) + cos(thetaL2)*sin(thetaR2)*cos(beta));
      F1over2(1) = -2*f(3)*f(1);
      F1over2(2) = -2*f(3)*f(2);
    }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementgammagammaviachargino2s\n");
  }

  Ichargino2(1) = Rchargino2*F1over2(1)*mWboson/mchargino2;

  Ichargino2(2) = Rchargino2*F1over2(2)*mWboson/mchargino2;
  return Ichargino2;

}



double higgsesamplitudedecaygluongluontotal(double m1, double g, double gs, double gprime, double mWboson, double alpha, double beta, double mtop, double mbottom, double mcharm, double mstop1, double mstop2, double msbottom1, double msbottom2, double thetat, double thetab, double greekmu, double Atop, double Abottom, double mstrange, double mscharmL, double mscharmR, double msstrangeL, double msstrangeR, double Acharm, double Astrange, double mup, double mdown, double msupL, double msupR, double msdownL, double msdownR, double Aup, double Adown, char higgstype, bool QCD) /// function that calculates the partial width for h->gluon gluon at 1-loop (forbidden at tree-level)
{
  double prefactor=0, Itr=0, Iti=0, Ist1r=0, Ist1i=0, Ist2r=0, Ist2i=0, Ibr=0, Ibi=0, Isb1r=0, Isb1i=0, Isb2r=0, Isb2i=0, Icr=0, Ici=0, IssLr = 0, IssLi = 0, IssRr = 0, IssRi = 0, IscLr = 0, IscLi = 0, IscRr = 0, IscRi = 0, IsdLr = 0, IsdLi = 0, IsdRr = 0, IsdRi = 0, IsuLr = 0, IsuLi = 0, IsuRr = 0, IsuRi = 0, matelemmodsquare=0, amplitudeW=0;
  /// Now as the only difference is in the prefactor I can use the same functions for each of the loop contributions as in the higgs to gamma gamma case.
  DoubleVector higgsmatrixelementgammagammaviatops (double m1, double mtop, double alpha, double beta, char higgstype);
  DoubleVector higgsmatrixelementgammagammaviastops (double m1, double mstop1, double mstop2, double mtop, double mbottom, double mWboson, double thetat, double g, double gprime, double alpha, double beta, double greekmu, double Atop, double Abottom, char higgstype);
  DoubleVector higgsmatrixelementgammagammaviabottoms (double m1, double mbottom, double alpha, double beta, char higgstype);
  DoubleVector higgsmatrixelementgammagammaviasbottoms (double m1, double msbottom1, double msbottom2, double mbottom, double mtop, double mWboson, double thetab, double g, double gprime, double alpha, double beta, double Atop, double Abottom, double greekmu, char higgstype);
  DoubleVector higgsmatrixelementgammagammaviacharms (double m1, double mcharm, double alpha, double beta, char higgstype);


  Itr = 1* higgsmatrixelementgammagammaviatops (m1, mtop, alpha, beta, higgstype)(1);
  Iti = 1* higgsmatrixelementgammagammaviatops (m1, mtop, alpha, beta, higgstype)(2);
  Ist1r = -higgsmatrixelementgammagammaviastops (m1, mstop1, mstop2, mtop, mbottom, mWboson, thetat, g, gprime, alpha, beta, greekmu, Atop, Abottom, higgstype)(1);

  Ist1i = -higgsmatrixelementgammagammaviastops (m1, mstop1, mstop2, mtop, mbottom, mWboson, thetat, g, gprime, alpha, beta, greekmu, Atop, Abottom, higgstype)(2);
  Ist2r = -higgsmatrixelementgammagammaviastops (m1, mstop1, mstop2, mtop, mbottom, mWboson, thetat, g, gprime, alpha, beta, greekmu, Atop, Abottom, higgstype)(3);
  Ist2i = -higgsmatrixelementgammagammaviastops (m1, mstop1, mstop2, mtop, mbottom, mWboson, thetat, g, gprime, alpha, beta, greekmu, Atop, Abottom, higgstype)(4);
  Ibr = 4*higgsmatrixelementgammagammaviabottoms (m1, mbottom, alpha, beta, higgstype)(1);
  Ibi = 4*higgsmatrixelementgammagammaviabottoms (m1, mbottom, alpha, beta, higgstype)(2);

  Isb1r = -4*higgsmatrixelementgammagammaviasbottoms (m1, msbottom1, msbottom2, mbottom, mtop, mWboson, thetab, g, gprime, alpha, beta, Atop, Abottom, greekmu, higgstype)(1);
  Isb1i = -4*higgsmatrixelementgammagammaviasbottoms (m1, msbottom1, msbottom2, mbottom, mtop, mWboson, thetab, g, gprime, alpha, beta, Atop, Abottom, greekmu, higgstype)(2);
  Isb2r = -4*higgsmatrixelementgammagammaviasbottoms (m1, msbottom1, msbottom2, mbottom, mtop, mWboson, thetab, g, gprime, alpha, beta, Atop, Abottom, greekmu, higgstype)(3);
  Isb2i = -4*higgsmatrixelementgammagammaviasbottoms (m1, msbottom1, msbottom2, mbottom, mtop, mWboson, thetab, g, gprime, alpha, beta, Atop, Abottom, greekmu, higgstype)(4);
  Icr = 1*higgsmatrixelementgammagammaviacharms (m1, mcharm, alpha, beta, higgstype)(1);
  Ici = 1*higgsmatrixelementgammagammaviacharms (m1, mcharm, alpha, beta, higgstype)(2);

  IscLr = -higgsmatrixelementgammagammaviastops (m1, mscharmL, mscharmR, mcharm, mstrange, mWboson, 0, g, gprime, alpha, beta, greekmu, Acharm, Astrange, higgstype)(1);
  IscLi = -higgsmatrixelementgammagammaviastops (m1, mscharmL, mscharmR, mcharm, mstrange, mWboson, 0, g, gprime, alpha, beta, greekmu, Acharm, Astrange, higgstype)(2);
  IscRr = -higgsmatrixelementgammagammaviastops (m1, mscharmL, mscharmR, mcharm, mstrange, mWboson, 0, g, gprime, alpha, beta, greekmu, Acharm, Astrange, higgstype)(3);
  IscRi = -higgsmatrixelementgammagammaviastops (m1, mscharmL, mscharmR, mcharm, mstrange, mWboson, 0, g, gprime, alpha, beta, greekmu, Acharm, Astrange, higgstype)(4);
  IssLr = -4*higgsmatrixelementgammagammaviasbottoms (m1, msstrangeL, msstrangeR, mstrange, mcharm, mWboson, 0, g, gprime, alpha, beta, Acharm, Astrange, greekmu, higgstype)(1);
  IssLi = -4*higgsmatrixelementgammagammaviasbottoms (m1, msstrangeL, msstrangeR, mstrange, mcharm, mWboson, 0, g, gprime, alpha, beta, Acharm, Astrange, greekmu, higgstype)(2);
  IssRr = -4*higgsmatrixelementgammagammaviasbottoms (m1, msstrangeL, msstrangeR, mstrange, mcharm, mWboson, 0, g, gprime, alpha, beta, Acharm, Astrange, greekmu, higgstype)(3);
  IssRi = -4*higgsmatrixelementgammagammaviasbottoms (m1, msstrangeL, msstrangeR, mstrange, mcharm, mWboson, 0, g, gprime, alpha, beta, Acharm, Astrange, greekmu, higgstype)(4);
  IsuLr = -higgsmatrixelementgammagammaviastops (m1, msupL, msupR, mup, mdown, mWboson, 0, g, gprime, alpha, beta, greekmu, Aup, Adown, higgstype)(1);
  IsuLi = -higgsmatrixelementgammagammaviastops (m1, msupL, msupR, mup, mdown, mWboson, 0, g, gprime, alpha, beta, greekmu, Aup, Adown, higgstype)(2);
  IsuRr = -higgsmatrixelementgammagammaviastops (m1, msupL, msupR, mup, mdown, mWboson, 0, g, gprime, alpha, beta, greekmu, Aup, Adown, higgstype)(3);
  IsuRi = -higgsmatrixelementgammagammaviastops (m1, msupL, msupR, mup, mdown, mWboson, 0, g, gprime, alpha, beta, greekmu, Aup, Adown, higgstype)(4);
  IsdLr = -4*higgsmatrixelementgammagammaviasbottoms (m1, msdownL, msdownR, mdown, mup, mWboson, 0, g, gprime, alpha, beta, Aup, Adown, greekmu, higgstype)(1);
  IsdLi = -4*higgsmatrixelementgammagammaviasbottoms (m1, msdownL, msdownR, mdown, mup, mWboson, 0, g, gprime, alpha, beta, Aup, Adown, greekmu, higgstype)(2);
  IsdRr = -4*higgsmatrixelementgammagammaviasbottoms (m1, msdownL, msdownR, mdown, mup, mWboson, 0, g, gprime, alpha, beta, Aup, Adown, greekmu, higgstype)(3);
  IsdRi = -4*higgsmatrixelementgammagammaviasbottoms (m1, msdownL, msdownR, mdown, mup, mWboson, 0, g, gprime, alpha, beta, Aup, Adown, greekmu, higgstype)(4);

  DoubleVector matelemsum(2);
  matelemsum(1) = Itr + Ist1r + Ist2r + Ibr + Isb1r + Isb2r + Icr + IscLr + IscRr + IssLr + IssRr + IsuLr + IsuRr + IsdLr + IsdRr;
  matelemsum(2) = Iti + Ist1i + Ist2i + Ibi + Isb1i + Isb2i + Ici + IscLi + IscRi + IssLi + IssRi + IsuLi + IsuRi + IsdLi + IsdRi;

  prefactor = pow(gs,4)*GFosqrt2/(128*pow(PI,5)*16)*pow(m1,3)*9./8;

  matelemmodsquare = pow(matelemsum(1),2) + pow(matelemsum(2),2);

  amplitudeW = prefactor*matelemmodsquare;

  double SMTOTRE = 0, SMTOTIM = 0, SQTOTRE = 0, SQTOTIM = 0;
  SQTOTRE = Ist1r + Ist2r + Isb1r + Isb2r + IscLr + IscRr + IssLr + IssRr + IsuLr + IsuRr + IsdLr + IsdRr;
  SQTOTIM = Ist1i + Ist2i + Isb1i + Isb2i + IscLi + IscRi + IssLi + IssRi + IsuLi + IsuRi + IsdLi + IsdRi;
  SMTOTRE = matelemsum(1) - SQTOTRE; SMTOTIM = matelemsum(2) - SQTOTIM;

  if (QCD == true) {
    int NF = 0;
    if (higgstype == 'h') { NF = 5;}
    else if (higgstype == 'H' || higgstype == 'A') { NF = 6;}
    else {
      throw("Problem: higgstype must be 'h', 'H' or 'A' in higgsesamplitudedecaygluongluontotal as in MSSM!\n");
    }
    amplitudeW = hggQCDcorrections(amplitudeW, pow(gs,2)/(4*PI), NF, higgstype, prefactor, SMTOTRE, SMTOTIM, SQTOTRE, SQTOTIM)(1);

  }
  else if (QCD == false) {
    amplitudeW = amplitudeW;
  }
      
  return amplitudeW;

}

///Function for QCD corrections to h->gg
DoubleVector hggQCDcorrections(double amplitudeW, double alphas, int Nf, char higgs, double prefactor, double SMtotr, double SMtoti, double sqtotr, double sqtoti) {
  double amplitudeafterFQCD = 0, amplitudefromSQCD = 0, amplitude = 0, hggstargcc = 0, hggstargbb = 0;
  double fermionQCDcorrections(double amplitudeW, double alphasnow, double alphasprev, int Nf, char higgs);
  double susyQCDcorrections(double prefactor, double alphas, double totr, double toti, double sqtotr, double sqtoti);
  DoubleVector Returns(3);
  for (int i = 1; i<=3; i++) {
    Returns(i) = 0;
  }

  amplitudeafterFQCD = fermionQCDcorrections(amplitudeW,alphas,alphas,Nf,higgs);

  amplitudefromSQCD = susyQCDcorrections(prefactor, alphas, SMtotr, SMtoti, sqtotr, sqtoti);

  amplitude = amplitudeafterFQCD + amplitudefromSQCD;
  
  Returns(1) = amplitude; Returns(2) = hggstargcc; Returns(3) = hggstargbb;
  return Returns;
}

///Function for the fermionic QCD corrections to e.g. h->gg
double fermionQCDcorrections(double amplitudeW, double alphasnow, double alphasprev, int Nf, char higgs) ///alphasprev required so can multiply by ratio accounting for change in alphas effect in prefactor of hgg in gluongluontotal function above
{
  double amplitude = 0;
  if (higgs == 'h' || higgs == 'H') {
    amplitude = amplitudeW*(1+alphasnow/PI*(95./4-7./6*Nf))*pow(alphasnow/alphasprev,2);
  }
  else if (higgs == 'A') {
    amplitude = amplitudeW*(1+alphasnow/PI*(97./4-7./6*Nf))*pow(alphasnow/alphasprev,2);
  }
  else {
    throw("problem: higgstype must be 'h' or 'H' or 'A' in fermionQCDcorrections as must be in MSSM\n");
  }
  return amplitude;
}

///Function for the susy QCD corrections to h->gg, i.e. the corrections to the squark loops
double susyQCDcorrections(double prefactor, double alphas, double SMtotr, double SMtoti, double sqtotr, double sqtoti) ///
{
  double amplitude = 0, factor = 0;
  factor = 17./6*alphas/PI; ///additional correction at NLO to squark loops
  amplitude = prefactor*(SMtotr*sqtotr+pow(sqtotr,2)+SMtoti*sqtoti+pow(sqtoti,2))*factor; ///Note the (SMr*sqr+pow(sqr,2)+SMi*sqi+pow(sqi,2)) is just Re[c.c. of total matrix element x matrix element part from squarks] this is the interference of the total matrix element with the squark parts of the matrix element which get extra alphas susy qcd corrections. Note amplitude here must be added to the amplitude inc FQCD corrections to get the total FQCD and SQCD corrected amplitude.

  return amplitude;
}

double higgsamplitudedecayVVstar (double m1, double mboson, double g, double gp, double beta, double alpha, char Vtype, DoubleMatrix & CPEMix, bool nmssmIsIt, int higgs) ///Function that calculates higgs to VV* to Vff'bar. Formula derived as in Marciano and Keung, massless fermion limit and zero W width considered.
{
  double prefactor=0, epsilon=0, integrals=0, a=0, b=0, c=0, coupling = 0, amplitudeW=0;
  double sin2thetaW = pow(gp,2)/(pow(g,2)+pow(gp,2));
  double cos2thetaW = pow(g,2)/(pow(g,2)+pow(gp,2));

  if (higgs == 1) {
    if (nmssmIsIt == false) { coupling = sin(beta-alpha);}
    else if (nmssmIsIt == true) { coupling = CPEMix(1,1)*sin(beta) + CPEMix(1,2)*cos(beta);}
  }
  else if (higgs == 2) {
     if (nmssmIsIt == false) { coupling = cos(beta-alpha);}
    else if (nmssmIsIt == true) { coupling = CPEMix(2,1)*sin(beta) + CPEMix(2,2)*cos(beta);}
  }
  else if (higgs == 2) {
    if (nmssmIsIt == false) { coupling = 0;} ///Only get H3 in NMSSM
    else if (nmssmIsIt == true) { coupling = CPEMix(3,1)*sin(beta) + CPEMix(3,2)*cos(beta);}
  }
  else {
    throw("problem: higgs must be 1 or 2 in higgsamplituddecayVVstar\n");
  }

  if(m1 < 2*mboson && m1 > mboson) {
    if (Vtype == 'W') {
      prefactor = 3*pow(g,4)*m1/(512*pow(PI,3))*pow(coupling,2);
    }
    else if (Vtype == 'Z') {
      prefactor = pow(g,4)*m1/(2048*pow(PI,3)*pow(cos2thetaW,2))*(7 - 40*sin2thetaW/3 + 160*pow(sin2thetaW,2)/9)*pow(coupling,2);
    }
    else {
      throw("problem: Vtype must be W or Z in higgsamplitudedecayVVstar\n");
    }
  epsilon = mboson/m1;
  a = 3*(1-8*pow(epsilon,2)+20*pow(epsilon,4))/(pow(4*pow(epsilon,2)-1,0.5))*acos((3*pow(epsilon,2)-1)/(2*pow(epsilon,3)));
  b = (1-pow(epsilon,2))*(47./2*pow(epsilon,2)-13./2+1/(pow(epsilon,2)));
  c = 3*(1-6*pow(epsilon,2)+4*pow(epsilon,4))*log(epsilon);
  integrals = a - b - c;
  amplitudeW = prefactor*integrals;
  }
  else if(m1> 2*mboson || m1< mboson) {
     throw("problem: in higgsamplitudedecayVVstar m1 either greater than 2*mboson so both on-shell, or m1 less than mboson so both off-shell, this formula is for one on-shell vector boson and one off-shell vector boson (of the same type of course)\n"); 
  }
  return amplitudeW;
}

DoubleVector higgshamplitudedecayVV(double m1, double mWboson, double mZboson, double g, double gp, double alpha, double beta, char Vtype, DoubleMatrix & CPEMix, bool nmssmIsIt) ///Function that calculates the light Higgs decays to two vector bosons, both on-shell or one off-shell
{
  DoubleVector Returns(2);
  for (int i=1; i<=2; i++) {
    Returns(i) = 0;
  }
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, prefactor=0, m2=0, m3=0, massratio=0, coupling = 0;
  if (Vtype == 'W')
    {
      if (fabs(m1) < fabs(mWboson)) {
	amplitudeW = 0; ///I don't calculate cases where both are off-shell
      }
      else if (fabs(m1) < 2*fabs(mWboson)) {
	amplitudeW = higgsamplitudedecayVVstar (m1, mWboson, g, gp, beta, alpha, 'W', CPEMix, nmssmIsIt, 1);
	  Returns(2) = 1;
	}
      else if (fabs(m1) >= 2*fabs(mWboson)) {
	m2 = mWboson, m3 = mWboson;
	squareplus = 1 - pow(m3/m1+m2/m1,2);
	squareminus = 1 - pow(m3/m1-m2/m1,2);
	lambda = pow(squareplus*squareminus,0.5);
	if (squareplus*squareminus < 0) {
	  throw("lambda will give nan in higgshamplitudedecayVV for W\n");
	} 
	// prefactor = pow(g,2)*pow(m1,3)/(64*PI*pow(mWboson,2));
	prefactor = GFosqrt2*pow(m1,3)/(8*PI);
	massratio = 4*pow(mWboson/m1,2);
	if (nmssmIsIt == false) { coupling = sin(beta-alpha);}
	else if (nmssmIsIt == true) { coupling = CPEMix(1,1)*sin(beta) + CPEMix(1,2)*cos(beta);}
	amplitudeW = lambda*prefactor*(1-massratio + 3./4*pow(massratio,2))*pow(coupling,2);
	Returns(2) = 2;
      }

    }
    else if (Vtype == 'Z') {
      if (fabs(m1) < fabs(mZboson)) {
	amplitudeW = 0; ///I don't calculate cases where both are off-shell
      }
      else if (fabs(m1) < 2*fabs(mZboson)) {
	amplitudeW = higgsamplitudedecayVVstar (m1, mZboson, g, gp, beta, alpha, 'Z', CPEMix, nmssmIsIt, 1);
	Returns(2) = 1;
      }
      else if (fabs(m1) >= 2*fabs(mZboson)) {
	m2 = mZboson, m3 = mZboson;
	squareplus = 1 - pow(m3/m1+m2/m1,2);
	squareminus = 1 - pow(m3/m1-m2/m1,2);
	lambda = pow(squareplus*squareminus,0.5);
	if (squareplus*squareminus < 0) {
	  throw("lambda will give nan in higgshamplitudedecayVV for Z\n");
	} 
	// prefactor = pow(g,2)*pow(m1,3)/(128*PI*pow(mWboson,2));	
	prefactor = GFosqrt2*pow(m1,3)/(16*PI);
	massratio = 4*pow(mZboson/m1,2);
	if (nmssmIsIt == false) { coupling = sin(beta-alpha);}
	else if (nmssmIsIt == true) { coupling = CPEMix(1,1)*sin(beta) + CPEMix(1,2)*cos(beta);}
	amplitudeW = lambda*prefactor*(1-massratio + 3./4*pow(massratio,2))*pow(coupling,2);
	Returns(2) = 2;
      }
    }
    else {
      throw("problem: Vtype must be W or Z in higgshamplitudedecayVV\n");
    }
  Returns(1) = amplitudeW;
  return Returns;
}
 

DoubleVector higgsHamplitudedecayVV(double m1, double mWboson, double mZboson, double g, double gp, double alpha, double beta, char Vtype, DoubleMatrix & CPEMix, bool nmssmIsIt) ///Function that calculates the Heavy Higgs decays to two vector bosons, assuming both are on-shell.
{
  DoubleVector Returns(2);
  for (int i=1; i<=2; i++) {
    Returns(i) = 0;
  }
  double higgsamplitudedecayVVstar (double m1, double mboson, double g, double gp, double beta, double alpha, char Vtype , DoubleMatrix & CPEMix, bool nmssmIsIt, int higgs);
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, prefactor=0, m2=0, m3=0, massratio=0, coupling = 0;
  if (Vtype == 'W')
    {
      if (fabs(m1) < fabs(mWboson)) {
	amplitudeW = 0;
      }
      else if (fabs(m1) < 2*fabs(mWboson)) {
	amplitudeW = higgsamplitudedecayVVstar (m1, mWboson, g, gp, beta, alpha, 'W', CPEMix, nmssmIsIt, 2);
	  Returns(2) = 1;
	}
      else if (fabs(m1) >= 2*fabs(mWboson)) {
	m2 = mWboson, m3 = mWboson;
	squareplus = 1 - pow(m3/m1+m2/m1,2);
	squareminus = 1 - pow(m3/m1-m2/m1,2);
	lambda = pow(squareplus*squareminus,0.5);
	if (squareplus*squareminus < 0) {
	  throw("lambda will give nan in higgsHamplitudedecayVV for W\n");
	} 
	// prefactor = pow(g,2)*pow(m1,3)/(64*PI*pow(mWboson,2));
	prefactor = GFosqrt2*pow(m1,3)/(8*PI);
	massratio = 4*pow(mWboson/m1,2);
	if (nmssmIsIt == false) { coupling = cos(beta-alpha);}
	else if (nmssmIsIt == true) { coupling = CPEMix(2,1)*sin(beta) + CPEMix(2,2)*cos(beta);}
	amplitudeW = lambda*prefactor*(1-massratio + 3./4*pow(massratio,2))*pow(coupling,2);
	Returns(2) = 2;
      }

    }
    else if (Vtype == 'Z') {
      if (fabs(m1) < fabs(mZboson)) {
	amplitudeW = 0;
      }
      else if (fabs(m1) < 2*fabs(mZboson)) {
	amplitudeW = higgsamplitudedecayVVstar (m1, mZboson, g, gp, beta, alpha, 'Z', CPEMix, nmssmIsIt, 2);
	Returns(2) = 1;
      }
      else if (fabs(m1) >= 2*fabs(mZboson)) {
	m2 = mZboson, m3 = mZboson;
	squareplus = 1 - pow(m3/m1+m2/m1,2);
	squareminus = 1 - pow(m3/m1-m2/m1,2);
	lambda = pow(squareplus*squareminus,0.5);
	if (squareplus*squareminus < 0) {
	  throw("lambda will give nan in higgsHamplitudedecayVV for Z\n");
	} 
	// prefactor = pow(g,2)*pow(m1,3)/(128*PI*pow(mWboson,2));	
	prefactor = GFosqrt2*pow(m1,3)/(16*PI);
	massratio = 4*pow(mZboson/m1,2);
	if (nmssmIsIt == false) { coupling = cos(beta-alpha);}
	else if (nmssmIsIt == true) { coupling = CPEMix(2,1)*sin(beta) + CPEMix(2,2)*cos(beta);}
	amplitudeW = lambda*prefactor*(1-massratio + 3./4*pow(massratio,2))*pow(coupling,2);
	Returns(2) = 2;
      }
    }
    else {
      throw("problem: Vtype must be W or Z in higgsHamplitudedecayVV\n");
    }
  Returns(1) = amplitudeW;
   return Returns;
}

DoubleVector higgsH3amplitudedecayVVNMSSM(double m1, double mWboson, double mZboson, double g, double gp, double alpha, double beta, char Vtype, DoubleMatrix & CPEMix, bool nmssmIsIt) ///Function that calculates the Heavy Higgs decays to two vector bosons, assuming both are on-shell.
{
  DoubleVector Returns(2);
  for (int i=1; i<=2; i++) {
    Returns(i) = 0;
  }
  double higgsamplitudedecayVVstar (double m1, double mboson, double g, double gp, double beta, double alpha, char Vtype , DoubleMatrix & CPEMix, bool nmssmIsIt, int higgs);
  double amplitudeW=0, squareplus=0, squareminus=0, lambda=0, prefactor=0, m2=0, m3=0, massratio=0, coupling = 0;
  if (Vtype == 'W')
    {
      if (fabs(m1) < fabs(mWboson)) {
	amplitudeW = 0;
      }
      else if (fabs(m1) < 2*fabs(mWboson)) {
	amplitudeW = higgsamplitudedecayVVstar (m1, mWboson, g, gp, beta, alpha, 'W', CPEMix, nmssmIsIt, 2);
	  Returns(2) = 1;
	}
      else if (fabs(m1) >= 2*fabs(mWboson)) {
	m2 = mWboson, m3 = mWboson;
	squareplus = 1 - pow(m3/m1+m2/m1,2);
	squareminus = 1 - pow(m3/m1-m2/m1,2);
	lambda = pow(squareplus*squareminus,0.5);
	if (squareplus*squareminus < 0) {
	  throw("lambda will give nan in higgsH3amplitudedecayVVNMSSM for W\n");
	} 
	prefactor = GFosqrt2*pow(m1,3)/(8*PI);
	massratio = 4*pow(mWboson/m1,2);
	if (nmssmIsIt == false) { coupling = 0;} ///Only get H3 in NMSSM
	else if (nmssmIsIt == true) { coupling = CPEMix(3,1)*sin(beta) + CPEMix(3,2)*cos(beta);}
	amplitudeW = lambda*prefactor*(1-massratio + 3./4*pow(massratio,2))*pow(coupling,2);
	Returns(2) = 2;
      }

    }
    else if (Vtype == 'Z') {
      if (fabs(m1) < fabs(mZboson)) {
	amplitudeW = 0;
      }
      else if (fabs(m1) < 2*fabs(mZboson)) {
	amplitudeW = higgsamplitudedecayVVstar (m1, mZboson, g, gp, beta, alpha, 'Z', CPEMix, nmssmIsIt, 2);
	Returns(2) = 1;
      }
      else if (fabs(m1) >= 2*fabs(mZboson)) {
	m2 = mZboson, m3 = mZboson;
	squareplus = 1 - pow(m3/m1+m2/m1,2);
	squareminus = 1 - pow(m3/m1-m2/m1,2);
	lambda = pow(squareplus*squareminus,0.5);
	if (squareplus*squareminus < 0) {
	  throw("lambda will give nan in higgsH3amplitudedecayVVNMSSM for Z\n");
	} 
	prefactor = GFosqrt2*pow(m1,3)/(16*PI);	
	massratio = 4*pow(mZboson/m1,2);
	if (nmssmIsIt == false) { coupling = 0;} ///Only get H3 in NMSSM
	else if (nmssmIsIt == true) { coupling = CPEMix(3,1)*sin(beta) + CPEMix(3,2)*cos(beta);}
	amplitudeW = lambda*prefactor*(1-massratio + 3./4*pow(massratio,2))*pow(coupling,2);
	Returns(2) = 2;
      }
    }
    else {
      throw("problem: Vtype must be W or Z in higgsH3amplitudedecayVV\n");
    }
  Returns(1) = amplitudeW;

  return Returns;
}


double higgsesamplitudedecayZbosonphotontotal(double m1, double mZboson, double g, double gprime, double alphaEmrun, double polemw, double runmw, double alpha, double beta, double mtop, double mbottom, double mcharm, double mstrange, double mstop1, double mstop2, double msbottom1, double msbottom2, double mHplus, double thetat, double thetab, double greekmu, double Atop, double Abottom, char higgstype) /// function that calculates the partial width for h->Z gamma at 1-loop (forbidden at tree-level)
{
  double prefactor=0, Itr=0, Iti=0, Ibr=0, Ibi=0, Icr=0, Ici=0, Isr=0, Isi=0, IWr=0, IWi=0, IHpr=0, IHpi=0, matelemmodsquare=0, amplitudeW=0;
 
  if (m1 < mZboson) {
    amplitudeW = 0;
  }
  else {
 
     DoubleVector higgsmatrixelementZgammaviafermions (double m1, double mferm, double mZboson, double mWboson, double alpha, double beta, double g, double gprime, double Qferm, double I3ferm, char fermtype, char higgstype);
     DoubleVector higgsmatrixelementZgammaviaWbosons (double m1, double mWboson, double mZboson, double alpha, double beta, double g, double gprime, char higgstype);
     DoubleVector higgsmatrixelementZgammaviaHplus (double m1, double mWboson, double mZboson, double mHplus, double alpha, double beta, double g, double gprime, char higgstype);

     Itr = 1* higgsmatrixelementZgammaviafermions (m1, mtop, mZboson, polemw, alpha, beta, g, gprime, 2./3, 1./2, 'u', higgstype)(1);
     Iti = 1* higgsmatrixelementZgammaviafermions (m1, mtop, mZboson, polemw, alpha, beta, g, gprime, 2./3, 1./2, 'u', higgstype)(2);

     Ibr = 1* higgsmatrixelementZgammaviafermions (m1, mbottom, mZboson, polemw, alpha, beta, g, gprime, -1./3, -1./2, 'd', higgstype)(1);
     Ibi = 1* higgsmatrixelementZgammaviafermions (m1, mbottom, mZboson, polemw, alpha, beta, g, gprime, -1./3, -1./2, 'd', higgstype)(2);

     Icr = 1* higgsmatrixelementZgammaviafermions (m1, mcharm, mZboson, polemw, alpha, beta, g, gprime, 2./3, 1./2, 'u', higgstype)(1);
     Ici = 1* higgsmatrixelementZgammaviafermions (m1, mcharm, mZboson, polemw, alpha, beta, g, gprime, 2./3, 1./2, 'u', higgstype)(2);
     Isr = 1* higgsmatrixelementZgammaviafermions (m1, mstrange, mZboson, polemw, alpha, beta, g, gprime, -1./3, -1./2, 'd', higgstype)(1);
     Isi = 1* higgsmatrixelementZgammaviafermions (m1, mstrange, mZboson, polemw, alpha, beta, g, gprime, -1./3, -1./2, 'd', higgstype)(2);
     IWr = higgsmatrixelementZgammaviaWbosons (m1, polemw, mZboson, alpha, beta, g, gprime, higgstype)(1);
     IWi = higgsmatrixelementZgammaviaWbosons (m1, polemw, mZboson, alpha, beta, g, gprime, higgstype)(2);
     IHpr = higgsmatrixelementZgammaviaHplus (m1, polemw, mZboson, mHplus, alpha, beta, g, gprime, higgstype)(1);
     IHpi = higgsmatrixelementZgammaviaHplus (m1, polemw, mZboson, mHplus, alpha, beta, g, gprime, higgstype)(2);  
        
     DoubleVector matelemsum(2);
     matelemsum(1) = Itr + Ibr + Icr + Isr + IWr + IHpr;
     matelemsum(2) = Iti + Ibi + Ici + Isi + IWi + IHpi;
     prefactor = pow(m1,3)/(64*pow(PI,3))*pow(1-pow(mZboson/m1,2),3)*GFosqrt2*pow(alphaEmrun,2);
         
     matelemmodsquare = pow(matelemsum(1),2) + pow(matelemsum(2),2);
     amplitudeW = prefactor*matelemmodsquare;
     
  }
  return amplitudeW;
}

DoubleVector higgsmatrixelementZgammaviafermions (double m1, double mferm, double mZboson, double mWboson, double alpha, double beta, double g, double gprime, double Qferm, double I3ferm, char fermtype, char higgstype) ///Calculates contribution of fermion loop to higgs to Z gamma decay
{
  double sinthetaW=0, costhetaW=0, Rhf=0, factor1=0, Integral1r=0, Integral1i=0, Integral2r=0, Integral2i=0;
  DoubleVector amplitude(2);
  DoubleVector ftau(3), flambda(3), gtau(3), glambda(3);
  for (int i=1; i<=3;i++) {
    ftau(i) = 0;
    flambda(i) = 0;
    gtau(i) = 0;
    glambda(i) = 0;
  }
  for (int i=1; i<=2;i++) {
    amplitude(i) = 0;
  }
  costhetaW = mWboson/mZboson;
  sinthetaW = pow(1 - pow(costhetaW,2),0.5);
  factor1 = -2*Qferm*(I3ferm-2*Qferm*pow(sinthetaW,2));
  if (higgstype == 'h') {
    if (fermtype == 'u') {
      Rhf = cos(alpha)/(sin(beta));
    }
    else if (fermtype == 'd') {
      Rhf = -sin(alpha)/cos(beta);
      }
    else {
      throw("problem: fermtype must be u or d in higgsmatrixelementZgammaviafermions\n");
    }
  }
  else if (higgstype == 'H') {
    if (fermtype == 'u') {
      Rhf = sin(alpha)/(sin(beta));
    }
    else if (fermtype == 'd') {
      Rhf = cos(alpha)/(cos(beta));
    }
    else {
      throw("problem: fermtype must be u or d in higgsmatrixelementZgammaviafermions\n");
    }
  }
  else if (higgstype == 'A') {
    if (fermtype == 'u') {
      Rhf = 1/tan(beta);
    }
    else if (fermtype == 'd') {
      Rhf = tan(beta);
    }
    else {
      throw("problem: fermtype must be u or d in higgsmatrixelementZgammaviafermions\n");
    }
  }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementZgammaviafermions\n");
  }
  ftau = foftau(mferm, m1);
  flambda = foftau(mferm, mZboson);
  gtau = goftau(mferm, m1);
  glambda = goftau (mferm, mZboson); 
  if (higgstype == 'h' || higgstype == 'H') {
    Integral1r = ftau(3)*flambda(3)/(2*(ftau(3)-flambda(3))) + pow(ftau(3)*flambda(3),2)/(2*pow((ftau(3)-flambda(3)),2))*(ftau(1)-flambda(1)) + pow(ftau(3),2)*flambda(3)/(pow(ftau(3)-flambda(3),2))*(gtau(1)-glambda(1));
    Integral1i = pow(ftau(3)*flambda(3),2)/(2*pow((ftau(3)-flambda(3)),2))*(ftau(2)-flambda(2)) + pow(ftau(3),2)*flambda(3)/(pow(ftau(3)-flambda(3),2))*(gtau(2)-glambda(2));
  }
  else if (higgstype == 'A') {
    Integral1r = 0;
    Integral1i = 0;
  }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementZgammaviafermions\n");
  }
  Integral2r = -ftau(3)*flambda(3)/(2*(ftau(3)-flambda(3)))*(ftau(1)-flambda(1));
  Integral2i = -ftau(3)*flambda(3)/(2*(ftau(3)-flambda(3)))*(ftau(2)-flambda(2));  
  
  amplitude(1) = 3*Rhf*factor1/(sinthetaW*costhetaW)*(Integral1r-Integral2r);
  amplitude(2) = 3*Rhf*factor1/(sinthetaW*costhetaW)*(Integral1i-Integral2i);
  return amplitude;
}

DoubleVector higgsmatrixelementZgammaviaWbosons (double m1, double mWboson, double mZboson, double alpha, double beta, double g, double gprime, char higgstype) ///Calculates contribution of W boson loop to higgs to Z gamma decay
{
  double tanthetaW=0, RhW=0, factor1=0, factor2=0, Integral1r=0, Integral1i=0, Integral2r=0, Integral2i=0;
  DoubleVector amplitude(2);

  DoubleVector ftau(3), flambda(3), gtau(3), glambda(3);
  for (int i=1; i<=3;i++) {
    ftau(i) = 0;
    flambda(i) = 0;
    gtau(i) = 0;
    glambda(i) = 0;
  }
  for (int i=1; i<=2;i++) {
    amplitude(i) = 0;
  }
  tanthetaW = pow((1-pow((mWboson/mZboson),2))/(pow(mWboson/mZboson,2)),0.5);
  if (higgstype == 'h') {
    RhW = sin(beta-alpha);
    }
  else if (higgstype == 'H') {
    RhW = cos(beta-alpha);
  }
  else if (higgstype == 'A') {
    RhW = 0;
  }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementZgammaviaWbosons\n");
  }
  ftau = foftau(mWboson, m1);
  flambda = foftau(mWboson, mZboson);
  gtau = goftau(mWboson, m1);
  glambda = goftau (mWboson, mZboson); 
  
  factor1 = 4*(3-pow(tanthetaW,2));
  factor2 = (1+2/ftau(3))*pow(tanthetaW,2) - (5+2/ftau(3));

  Integral1r = ftau(3)*flambda(3)/(2*(ftau(3)-flambda(3))) + pow(ftau(3)*flambda(3),2)/(2*pow((ftau(3)-flambda(3)),2))*(ftau(1)-flambda(1)) + pow(ftau(3),2)*flambda(3)/(pow(ftau(3)-flambda(3),2))*(gtau(1)-glambda(1));
  Integral1i =pow(ftau(3)*flambda(3),2)/(2*pow((ftau(3)-flambda(3)),2))*(ftau(2)-flambda(2)) + pow(ftau(3),2)*flambda(3)/(pow(ftau(3)-flambda(3),2))*(gtau(2)-glambda(2));
  Integral2r = -ftau(3)*flambda(3)/(2*(ftau(3)-flambda(3)))*(ftau(1)-flambda(1));
  Integral2i = -ftau(3)*flambda(3)/(2*(ftau(3)-flambda(3)))*(ftau(2)-flambda(2));
  amplitude(1) = -RhW*1/tanthetaW*(factor1*Integral2r + factor2*Integral1r);
  amplitude(2) = -RhW*1/tanthetaW*(factor1*Integral2i + factor2*Integral1i);

  return amplitude;
}

DoubleVector higgsmatrixelementZgammaviaHplus (double m1, double mWboson, double mZboson, double mHplus, double alpha, double beta, double g, double gprime, char higgstype) ///Calculates contribution of charged higgs boson loop to higgs to Z gamma decay
{
  double sinthetaW=0, costhetaW=0, RhHplus=0, factor1=0, Integral1r=0, Integral1i=0;
  DoubleVector amplitude(2);
  DoubleVector ftau(3), flambda(3), gtau(3), glambda(3);
  for (int i=1; i<=3;i++) {
    ftau(i) = 0;
    flambda(i) = 0;
    gtau(i) = 0;
    glambda(i) = 0;
  }
  for (int i=1; i<=2;i++) {
    amplitude(i) = 0;
  }
  costhetaW = mWboson/mZboson;
  sinthetaW = pow(1 - pow(costhetaW,2),0.5);
  if (higgstype == 'h') {
    RhHplus = sin(beta-alpha) + cos(2*beta)*sin(beta+alpha)/(2*pow(costhetaW,2));
    }
  else if (higgstype == 'H') {
    RhHplus = -(cos(beta-alpha) - cos(2*beta)*cos(beta+alpha)/(2*pow(costhetaW,2)));
  }
  else if (higgstype == 'A') {
    RhHplus =0;
  }
  else {
    throw("problem: higgstype must be h or H or A in higgsmatrixelementZgammaviaHplus\n");
  }
  ftau = foftau(mHplus, m1);
  flambda = foftau(mHplus, mZboson);
  gtau = goftau(mHplus, m1);
  glambda = goftau (mHplus, mZboson); 
  
  factor1 = 1 - 2*pow(sinthetaW,2);

  Integral1r = ftau(3)*flambda(3)/(2*(ftau(3)-flambda(3))) + pow(ftau(3)*flambda(3),2)/(2*pow((ftau(3)-flambda(3)),2))*(ftau(1)-flambda(1)) + pow(ftau(3),2)*flambda(3)/(pow(ftau(3)-flambda(3),2))*(gtau(1)-glambda(1));
  Integral1i = pow(ftau(3)*flambda(3),2)/(2*pow((ftau(3)-flambda(3)),2))*(ftau(2)-flambda(2)) + pow(ftau(3),2)*flambda(3)/(pow(ftau(3)-flambda(3),2))*(gtau(2)-glambda(2));
  
  amplitude(1) = RhHplus*factor1/(costhetaW*sinthetaW)*Integral1r*pow(mWboson/mHplus,2);
  amplitude(2) = RhHplus*factor1/(costhetaW*sinthetaW)*Integral1i*pow(mWboson/mHplus,2);
  
  return amplitude;
}

	  
/// Decays to gravitinos

double gluinoamplitudedecaygravitino (double m1, double mgrav, double MPlreduced, int gravonoff, int gluNLSP) /// Function that calculates the decays to gravitinos of gluinos if int gravonoff == 1, otherwise decays to gravitinos are off
{
  double amplitudeW = 0;
  if (gluNLSP == 0) { amplitudeW = 0;}
  else if (gluNLSP == 1) {
    if(gravonoff == 0) { amplitudeW=0;}
    else if (gravonoff == 1 && fabs(m1) > mgrav) { ///Comment this line, last line and line 3 below this out if always want decays to gravitinos output
      amplitudeW =  pow(fabs(m1),5)/(48*PI*pow(mgrav*MPlreduced,2));
    }
  }
  return amplitudeW;
}

double squarkamplitudedecaygravitino(double m1, double mgrav, double mquark, double MPlreduced, int gravonoff, int squNLSP) /// Function that calculates the decays to gravitinos of squark if int gravonoff == 1, otherwise decays to gravitinos are off
{
  double amplitudeW = 0;
  if (squNLSP == 0) { amplitudeW = 0;}
  else if (squNLSP == 1) {
    if(gravonoff == 0) { amplitudeW=0;}
    else if (gravonoff == 1 && fabs(m1) > mgrav) { ///Comment this line, last line and line 3 below this out if always want decays to gravitinos output
      amplitudeW =  pow(pow(m1,2)-pow(mquark,2),4)/(48*PI*pow(m1,3)*pow(MPlreduced*mgrav,2));
    }
  }
  return amplitudeW;
}

double neutralinoamplitudedecayphotongravitino(double m1, double mgrav, double MPlreduced, DoubleMatrix & mixNeut, double g, double gp, int neutralino, int gravonoff, int neutNLSP)  /// Function that calculates the decays to gravitino and photon of a neutralino if int gravonoff == 1, otherwise decays to gravitinos are off
{
  double amplitudeW = 0, coupling = 0, costhetaW = 0, sinthetaW = 0;
  if (neutNLSP == 0) { amplitudeW =0;}
  else if (neutNLSP == 1) {
    if(gravonoff == 0) { amplitudeW=0;}
    else if (gravonoff == 1 && fabs(m1) > mgrav) { ///Comment this line, last line and line 3 below this out if always want decays to gravitinos output
      costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
      sinthetaW = gp/(pow(pow(g,2)+pow(gp,2),0.5));
      coupling = mixNeut(neutralino,1)*costhetaW + mixNeut(neutralino,2)*sinthetaW;
      amplitudeW =  pow(coupling,2)*pow(fabs(m1),5)/(48*PI*pow(MPlreduced*mgrav,2));
    }
  }
  return amplitudeW;
}

double neutralinoamplitudedecayZgravitino(double m1, double mZ, double mgrav, double MPlreduced, DoubleMatrix & mixNeut, double g, double gp, double beta, int neutralino, int gravonoff, int neutNLSP)  /// Function that calculates the decays to gravitino and Z boson of a neutralino if int gravonoff == 1, otherwise decays to gravitinos are off
{
  double amplitudeW = 0, coupling = 0, costhetaW = 0, sinthetaW = 0;
  if (neutNLSP == 0 || fabs(m1) < mZ + mgrav) { amplitudeW = 0;}
  else if (neutNLSP == 1) {
    if(gravonoff == 0) { amplitudeW=0;}
    else if (gravonoff == 1) { ///Comment this line, last line and line 3 below this out if always want decays to gravitinos output
      costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));
      sinthetaW = gp/(pow(pow(g,2)+pow(gp,2),0.5));
      coupling = 2*pow((mixNeut(neutralino,1)*sinthetaW - mixNeut(neutralino,2)*costhetaW),2)+ pow((mixNeut(neutralino,4)*sin(beta) - mixNeut(neutralino,3)*cos(beta)),2);

      amplitudeW =  coupling*pow(pow(m1,2)-pow(mZ,2),4)/(96*PI*pow(MPlreduced*mgrav,2)*pow(fabs(m1),3));
    }
  }
  return amplitudeW;
}

double neutralinoamplitudedecayphigravitino(double m1, double mphi, double mgrav, double MPlreduced, DoubleMatrix & mixNeut, double alpha, double beta, int neutralino, int gravonoff, char phi, int neutNLSP)  /// Function that calculates the decays to gravitino and Z boson of a neutralino if int gravonoff == 1, otherwise decays to gravitinos are off
{
  double amplitudeW = 0, coupling = 0;
  if (neutNLSP == 0 || fabs(m1) < mphi) { amplitudeW = 0;}
  else if (neutNLSP == 1) {
    if(gravonoff == 0) { amplitudeW=0;}
    else if (gravonoff == 1 && fabs(m1) > mgrav) { ///Comment this line, last line and line 3 below this out if always want decays to gravitinos output
    if (phi =='h') {
      coupling = 1/(pow(6,0.5)*MPlreduced*mgrav)*(mixNeut(neutralino,4)*cos(alpha) - mixNeut(neutralino,3)*sin(alpha));
    }
    else if (phi =='H') {
      coupling = 1/(pow(6,0.5)*MPlreduced*mgrav)*(mixNeut(neutralino,4)*sin(alpha)+mixNeut(neutralino,3)*cos(alpha));
    }
    else if (phi == 'A') {
      coupling = 1/(pow(6,0.5)*MPlreduced*mgrav)*(mixNeut(neutralino,4)*cos(beta)+mixNeut(neutralino,3)*sin(beta));
    }
    else {  throw("problem: phi must be h or H or A in neutralinoamplitudedecayphigravitino\n");  }
    amplitudeW =  pow(coupling,2)/(16*PI*pow(fabs(m1),3))*pow((pow(m1,2)-pow(mphi,2)),4);
    }
  }
  return amplitudeW;
}

void slhaDecays(ostream & fout, vector<Particle> & decayTable,
		bool outputPartialWidths) {
  vector<Particle>::iterator i;
  for (i=decayTable.begin(); i<decayTable.end(); i++)
    if (!outputPartialWidths) OutputNoPWs(fout, *i);
    else OutputYesPWs(fout, *i);
}
