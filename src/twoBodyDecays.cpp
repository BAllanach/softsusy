/** \file threeBodyDecays.cpp
   Project:     SOFTSUSY 
   Author:      Tom Cridge, Ben Allanach
   Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "twoBodyDecays.h"

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
  DoubleVector squarkmixcharginocouplings (double g, double theta, double beta, double gammaL, double gammaR, double runmt, double runmb, double mWboson, double mch1, double mch2, int torb);
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
  DoubleVector squarkmixcharginocouplings (double g, double theta, double beta, double gammaL, double gammaR, double runmt, double runmb, double mWboson, double mch1, double mch2, int torb);
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
  DoubleVector squarkmixcharginocouplings (double g, double theta, double beta, double gammaL, double gammaR, double runmt, double runmb, double mWboson, double mch1, double mch2, int torb);

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
  DoubleVector squarkmixcharginocouplings (double g, double theta, double beta, double gammaL, double gammaR, double runmt, double runmb, double mWboson, double mch1, double mch2, int torb);
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
	  AqZ = 1/(pow(2,0.5))*(uord*g*mixNeut(neutralino,2) + gprime*mixNeut(neutralino,1)/3); 
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
	  BqZ = 1/(pow(2,0.5))*4./3*uordchanger*gprime*mixNeut(neutralino,1); ///Following changes in AqZ suggested by SUSYHIT
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
	fq = g*runmq/(pow(2,0.5)*mWboson*sin(beta));
	if (oneortwo == 1) /// we have stop1s
	  {
	    alphatilda = (cos(theta)*1/(pow(2,0.5))*(-g*mixNeut(neutralino,2) - gp/3*mixNeut(neutralino,1)) - fq*sin(theta)*mixNeut(neutralino,4)); /// this is just simpler (rearranged) form of atopr(1,i) from susyhit *g
	    betatilda = (4/(3*pow(2,0.5))*gp*mixNeut(neutralino,1)*sin(theta) - fq*mixNeut(neutralino,4)*cos(theta)); /// this is just simpler (rearranged) form of btopr(1,i) from susyhit * g
	  }
     	else if (oneortwo == 2) /// we have stop2s
	  {
	    alphatilda = (-sin(theta)*1/(pow(2,0.5))*(mixNeut(neutralino,1)*gp/3 + mixNeut(neutralino,2)*g) + cos(theta)*fq*mixNeut(neutralino,4)); /// again simplified form of atopr(2,i) * g from susyhit
	    betatilda = -4/(3*pow(2,0.5))*cos(theta)*mixNeut(neutralino,1)*gp - fq*sin(theta)*mixNeut(neutralino,4); /// again simplified form of btopr(2,i) * g from susyhit
	  }
	else {
	  throw("problem: stop must be stop1 or stop2\n"); 
	}
      }
    else if ( squark ==2) /// we have sbottoms
      {
	fq = g*runmq/(pow(2,0.5)*mWboson*cos(beta));
	
	if ( oneortwo == 1) /// we have sbottom1s
	  {
	    alphatilda = (1/pow(2,0.5))*cos(theta)*(-mixNeut(neutralino,1)*gp/3 + mixNeut(neutralino,2)*g) - sin(theta)*mixNeut(neutralino,3)*fq; /// just simplified form of abot(1,i) * g from susyhit
	    betatilda = sin(theta)*2/(3*pow(2,0.5))*(-mixNeut(neutralino,1)*gp) - cos(theta)*fq*mixNeut(neutralino,3);
	  }
	else if (oneortwo == 2) /// we have sbottom2s
	  {
	    alphatilda = sin(theta)*(mixNeut(neutralino,1)*-gp*(1/(3*pow(2,0.5))) + mixNeut(neutralino,2)*g*1/(pow(2,0.5))) + cos(theta)*fq*mixNeut(neutralino,3); /// simplified abot(2,i)*g from susyhit
	    betatilda = cos(theta)*(2/(3*pow(2,0.5)))*gp*(mixNeut(neutralino,1)) - sin(theta)*fq*mixNeut(neutralino,3); /// simplified bbot(2,i)*g from susyhit
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
      
    A11 = g/(pow(2,0.5)*mWboson)*(mt*mb*combo1*sin(thetat)*sin(thetab) + mt*combo2*sin(thetat)*cos(thetab) + mb*combo3*sin(thetab)*cos(thetat) + combo4*cos(thetat)*cos(thetab));
    A12 = g/(pow(2,0.5)*mWboson)*(mt*mb*combo1*sin(thetat)*(-cos(thetab)) + mt*combo2*sin(thetat)*sin(thetab) + mb*combo3*(-cos(thetab))*cos(thetat) + combo4*cos(thetat)*sin(thetab));
    A21 = g/(pow(2,0.5)*mWboson)*(mt*mb*combo1*(-cos(thetat))*sin(thetab) + mt*combo2*(-cos(thetat))*cos(thetab) + mb*combo3*sin(thetab)*sin(thetat) + combo4*sin(thetat)*cos(thetab));
    A22 = g/(pow(2,0.5)*mWboson)*(mt*mb*combo1*(-cos(thetat))*(-cos(thetab)) + mt*combo2*-(cos(thetat))*sin(thetab) + mb*combo3*(-cos(thetab))*sin(thetat) + combo4*sin(thetat)*sin(thetab));

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
    if(t1or2 == 2) /// we have an initial stop2
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
    
    amplitudeW = pow(A,2)*lam/(16*PI*m1 * sqr(m1));
  }
  return amplitudeW;
  
}


