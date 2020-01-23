/** \file threeBodyDecays.cpp
   Project:     SOFTSUSY 
   Author:      Tom Cridge, Ben Allanach
   Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   Webpage:     http://softsusy.hepforge.org/

*/

#include "threeBodyDecays.h"

using namespace std;

const double GFosqrt2 = GMU / root2;

///Integrand functions for 1->3 decays
double Zfunc(double m1, double mq, double m, double Etbarmax, double Etbarmin) ///required in many of the 1->3 integrals
{
  double Z = 0;
  Z = (sqr(m1)+sqr(mq)-2*fabs(m1)*Etbarmax - sqr(m))/(sqr(m1)+sqr(mq)-2*fabs(m1)*Etbarmin-sqr(m));
    if (Z <=0) {
      std::cout << "May have nan issue if do log(Z) as Z<=0! See Zfunc used in 1->3 decays\n" << endl;
  }
  return Z;
}

DoubleVector Etbarmaxmin(double m1, double m2, double massq, double Et) ///required for many of the 1->3 integrals, this gives E max and min of t bar for given Et therefore different to overall limits of integration on Et which are Etmax, Etmin
{
  DoubleVector Etbarsupremum(2);
  double pt=0, zet=0, A=0, B=0;
  pt = sqrt(sqr(Et)-sqr(massq));
  zet = 2.0*sqr(massq)+sqr(m1)-sqr(m2)-2*fabs(m1)*Et;
  A = sqr(m1)+sqr(massq)-2*fabs(m1)*Et;
  B = (sqr(pt*zet)-4*sqr(pt*massq)*A);
  // std::cout << "Et = " << Et << std::endl;
  // std::cout << "massq = " << massq << std::endl;
  // std::cout << "pt = " << pt << std::endl;
  // std::cout << "zet = " << zet << std::endl;
  // std::cout << "pt*zet = " << pt*zet << std::endl;
  // std::cout << "pt*massq = " << pt*massq << std::endl;
  //  if (printOutNow) cout << "DEBUG: pt=" << pt << " zet=" << zet << " A=" << A << " B=" << B << endl;
  if ( B < 0) {
    B = 0; /// B may become very very small and negative at the tails of this squareroot, this causes problems as we get sqrt(-ve) therefore set B to 0 here as it's very small anyway so has negligible effect on the overall answer.
  }
 
  Etbarsupremum(1) = (zet*(fabs(m1)-Et)+sqrt(B))/(2*A); ///Etbarmax
  Etbarsupremum(2) = (zet*(fabs(m1)-Et)-sqrt(B))/(2*A); ///Etbarmin
  // std::cout << "Etbarsupremum(1) = " << Etbarsupremum(1) << std::endl;
  // std::cout << "Etbarsupremum(2) = " << Etbarsupremum(2) << std::endl;
  // std::cout << "sqrt(B) = " << sqrt(B)<< std::endl;
  if (Etbarsupremum(1) != Etbarsupremum(1) || Etbarsupremum(2) != Etbarsupremum(2)) {
    std::cout << "problem: Etbar gives nan! See Etbarmaxmin used in 1->3 decays\n" << endl;
  }
    return Etbarsupremum;
}


DoubleVector Ebbarmaxmin (double mass1, double mass2, double mass3, double mass4, double Et) ///function required for gluino 3 body decays to chargino and q q'bar, mass1 is gluino mass, mass2 is top mass, mass3 is bottom mass, mass4 is chargino mass
{
  DoubleVector Ebbar(2);
  for (int i=1; i<=2; i++) { Ebbar(i)=0;}
  double pt = 0, squareplus = 0, squareminus = 0, lambda = 0, A = 0;
  pt = sqrt(sqr(Et)-sqr(mass2));
  A = sqr(mass1)+sqr(mass2)-2*Et*fabs(mass1);
  squareplus = A - sqr(mass3+mass4);
  squareminus = A - sqr(mass3-mass4);

  if (squareplus < 0 && fabs(squareplus) < 1e-8) {
    squareplus = 0; ///Avoid issues of finite precision meaning you get a very small negative squareplus at Amin (smin) rather than 0, this can cause issues when you take the squareroot of lambda, giving a nan
  }

  lambda = sqrt(squareplus*squareminus);
  if (squareplus*squareminus < 0) {
    std::cout << "problem: lambda will give nan in Ebbarmaxmin used in 1->3 decays\n"<< endl;
  } 
  Ebbar(1) = ((sqr(mass1)+sqr(mass2)-2*fabs(mass1)*Et-pow(mass3,2)-pow(mass4,2))*(fabs(mass1)-Et) + pt*lambda)/(2*(sqr(mass1)+sqr(mass2)-2*Et*fabs(mass1))); ///I have -pow(mass3,2) here whereas T&B have + pow(mass3,2), I have changed the sign to agree with SPheno, note however this actually makes negligible difference due to the smallness of the b mass
  Ebbar(2) = ((sqr(mass1)+sqr(mass2)-2*fabs(mass1)*Et-pow(mass3,2)-pow(mass4,2))*(fabs(mass1)-Et) - pt*lambda)/(2*(sqr(mass1)+sqr(mass2)-2*Et*fabs(mass1))); ///I have -pow(mass3,2) here whereas T&B have + pow(mass3,2), I have changed the sign to agree with SPheno, note however this actually makes negligible difference due to the smallness of the b mass

  return Ebbar;
}

 double Xfunc (double mass1, double mass2, double mass3, double mass4, double mass5, double Et) ///required for gluino 3 body decays to chargino and q q'bar, mass1 is gluino mass, mass2 is top mass, mass3 is bottom mass, mass4 is chargino mass, mass5 is the sbottom mass
{
  double X=0;
  DoubleVector Ebbarmaxmin (double mass1, double mass2, double mass3, double mass4, double Et);
  DoubleVector Ebbar(2);
  for (int i=1; i<=2; i++) { Ebbar(i)=0;}
  Ebbar = Ebbarmaxmin (mass1, mass2, mass3, mass4, Et);
  X = (pow(mass5,2)+2*Ebbar(1)*fabs(mass1)-sqr(mass1))/(pow(mass5,2)+2*Ebbar(2)*fabs(mass1)-sqr(mass1));
  return X;
}


///Functions for 1->3 decays via dgauss:
double gluinoamplitudedecaydgausscharginoqqpbarfirsttwogen (double mgluino, double mchargino, double mquark, double mquarkp, double msqL, double msqpL, double g, double thetaL, double thetaR, double alphas, int charg, bool onetothree)
{
  double amplitudeW=0, Au=0, Ad=0, psiu=0, psid=0, phiIu=0, phiId=0, upper=0, upperd=0;
  upper = (sqr(mgluino)+sqr(mquark)-sqr(mquarkp)-2*mquarkp*fabs(mchargino)-pow(mchargino,2))/(2*mgluino);
  upperd = (sqr(mgluino)+sqr(mquarkp)-sqr(mquark)-2*mquark*fabs(mchargino)-pow(mchargino,2))/(2*mgluino);
  // std::cout << "mgluino = " <<mgluino << std::endl;
  // std::cout << "mchargino = " << mchargino <<std::endl;
  // std::cout << "mquark = " << mquark << std::endl;
  // std::cout << "mquarkp = " << mquarkp << std::endl;
  if (onetothree == false) {
    amplitudeW = 0;
  }
  else if (onetothree == true) {
    
    if(mgluino < fabs(mchargino)+mquark+mquarkp) {
      amplitudeW = 0;
    }
    else if (mgluino > msqL + mquark || mgluino > msqpL + mquarkp) {
	amplitudeW = 0; /// 1->3 decay not relevant here if one to two decay open
      }

    else {
      if (charg == 1) /// chargino1 decaying
      {
	Ad = g*sin(thetaR); 
	Au = g*sin(thetaL);
      }
    else if (charg == 2) ///chargino2 decaying
      {
	Ad = g*cos(thetaR);
	Au = g*cos(thetaL);
      }
    else {
      throw("problem: charg must be 1 or 2 in gluinoamplitudedecaydgausscharginoqqpbarfirsttwogen\n");
    }
      // std::cout << "mquark = " << mquark << std::endl;
      // std::cout << "upper = " << upper << std::endl;
      m1 = mgluino, m2 = msqL, m3 = msqL, m4 = mchargino, mq = mquark;
      psiu = dgauss(gpsitildadgauss,mquark,upper,accuracy)*1/(sqr(PI)*m1);
      m1 = mgluino, m2 = msqpL, m3 = msqpL, m4 = mchargino, mq = mquarkp;
      psid = dgauss(gpsitildadgauss,mquarkp,upperd,accuracy)*1/(sqr(PI)*m1);
      m1 = mgluino, m2 = msqL, m3 = msqpL, m4 = mchargino, mq = mquarkp; ///Note mq isn't important really for these modes as they are first two gen quarks
      // std::cout << "psiu = " << psiu << std::endl;
      // std::cout << "psid = " << psid << std::endl;
      // std::cout << "mq = " << mq << std::endl;
      // std::cout << "mquark = " << mquark << std::endl;
      // std::cout << "upper = " << upper << std::endl;
      // std::cout << "mquarkp = " << mquarkp << std::endl;
      // std::cout << "upperd = " << upperd << std::endl;
      // // std::cout << "mq = " << mq << std::endl;
      // std::cout << "upper = " << upper << std::endl;
	
      phiIu = dgauss(gphitildadgauss,mquark,upper,accuracy)*1/(sqr(PI)*m1);
      // std::cout << "phiIu = " << phiIu << std::endl;
      phiId = dgauss(gphitildadgauss,mquarkp,upperd,accuracy)*1/(sqr(PI)*m1);
      // std::cout << "phiId = " << phiId << std::endl;

      
      amplitudeW = alphas/(16*sqr(PI))*(pow(Ad,2)*psiu + pow(Au,2)*psid + Au*Ad*(phiIu+phiId));
      // std::cout << "amplitudeW = " << amplitudeW << std::endl;
    }
  }
  return amplitudeW;
}

///Function that runs alphas to desired scale using 1-loop renormalisation group equations, inc. MSSM gluinos and squarks (hence 2/3 not 2/7)
double alphasrun (double mu, double mu0, double alphasmu0) {
  double alphasmu = (2*PI/3)/(2*PI/(3*alphasmu0) + log(mu/mu0)); ///MSSM
  // double alphasmu = (2*PI/7)/(2*PI/(7*alphasmu0) + log(mu/mu0)); ///SM
  return alphasmu;
}

///Function that mimics hdecay's ALPHAS running with the link to Lambda QCD
double alphasrunlambdaQCD (double mu, double LAMBDA, double Nf) {
  double alphasmu2 = 24*PI/((33-2*Nf)*log(mu/LAMBDA))*(1 - ((12*(153-19*Nf)/sqr(33-2*Nf))*log(2*log(mu/LAMBDA))/log(sqr(mu/LAMBDA))));
  return alphasmu2;
}

///Function that runs alpha to desired scale using 1-loop renormalisation group equations, inc. MSSM particles hence 1 rather than -19/6 of SM
double alpharun (double mu, double mu0, double alphamu0) {
  double alphamu = (2*PI)/(2*PI/(alphamu0) - log(mu/mu0)); ///MSSM
  // double alphamu = (12*PI)/19/(12*PI/(19*alphamu0) + log(mu/mu0)); ///SM
  return alphamu;
}


// Outputs a space before if greater than zero, a minus otherwise, also outputs spaces after depending on no. of digits in PDG code
// Useful for outputting negative numbers in rows, Tom Cridge added to output PDG codes with no. of space depending on length of PDG (usually 1 or 2 for SM particles whereas 7 for SUSY particles -  therefore output 5 extra spaces for SM)
void printRowPDG(ostream & fout, double x) {

  /// make it return a character when you've worked out the equivalent of printf

  double underflow = 1.0e-120;
  if (fabs(x) < underflow) x = 0.0; /// Traps -0.0
  if (x >= 0.0) fout << " " << x;
  else fout << x;
  if (fabs(x)<10) fout << "      ";
  else if (fabs(x)<100) fout << "     ";
  else if (fabs(x)<1000) fout << "    ";
}


void OutputNoPWs(ostream & fout, Particle & P) ///Outputs the decay table into the leshouchesOutput file with no PWs (partial widths) given, just branching ratios
 {
   fout << left << setw(6) << "#" << setw(12) << "PDG" << setw(18) << "Width" << endl;
   fout << "DECAY " << setw(12) << fixed << setprecision(0) << P.PDG << setw(12) << scientific << setprecision(8) <<  P.total_width << "   " << "# " << P.name << " decays" << endl;
   if (1-P.three_width/P.total_width > minBR) {
     fout << left << setw(6) << "#" << setw(18) << "BR" << setw(6) << "NDA" << setw(12) << left << "PDG1" << setw(12) << " PDG2" << setw(18) << "Comments" << endl;
     for (int k=0; k<P.No_of_Decays; k++) {
       if( P.Array_Decays[k][2] != 0 && P.Array_Decays[k][5] > minBR && P.Array_Decays[k][2] > 0 && P.Array_Decays[k][3] == 2) {
	 fout << left << setw(6) << " " << setw(18) << scientific << setprecision(8) << P.Array_Decays[k][5] << setprecision(0) << setw(6) << fixed << P.Array_Decays[k][3];  printRowPDG(fout, P.Array_Decays[k][0]); fout << "   "; printRowPDG(fout, P.Array_Decays[k][1]); fout << "   " << left << setprecision(0) << setw(15) << P.Array_Comments[k] << endl;
       }
     }
   }
   if (P.three_width/P.total_width > minBR) {
     fout << left << setw(6) << "#" << setw(18) << "BR" << setw(8) << "NDA" << setw(12) << left << " PDG1" << setw(12) << " PDG2" << setw(12) << " PDG3 " << setw(18) << "Comments" << endl;
     for (int k=0; k<P.No_of_Decays; k++) {
       if( P.Array_Decays[k][2] != 0 && P.Array_Decays[k][5] > minBR && P.Array_Decays[k][2] > 0 && P.Array_Decays[k][3] == 3) {
	 fout << left << setw(6) << " " << setw(18) << scientific << setprecision(8) << P.Array_Decays[k][5] << setprecision(0) << setw(6) << fixed << P.Array_Decays[k][3] << setw(2) << " ";  printRowPDG(fout, P.Array_Decays[k][0]); fout << "    "; printRowPDG(fout,P.Array_Decays[k][1]); fout << "    "; printRowPDG(fout,P.Array_Decays[k][4]); fout << "   " << left << setprecision(0) << setw(25) << P.Array_Comments[k] << endl;
       }
     }
   }
   fout << "#" << endl;
 }


void OutputYesPWs(ostream & fout, Particle & P) ///Outputs the decay table into the leshouchesOutput file with PWs (partial widths) given after the comments column so as not to affect SLHA form
 {
   fout << left << setw(6) << "#" << setw(12) << "PDG" << setw(18) << "Width" << endl;
   fout << "DECAY " << setw(12) << fixed << setprecision(0) << P.PDG << setw(12) << scientific << setprecision(8) <<  P.total_width << "   " << "# " << P.name << " decays" << endl;

   if (1-P.three_width/P.total_width > minBR) {
     fout << left << setw(6) << "# " << setw(20) << "BR " << setw(6) << "NDA " << setw(12) << "PDG1 " << setw(11) << "PDG2" << setw(18) << "Comments" << setw(18) << "PW" << endl;
     for (int k=0; k<P.No_of_Decays; k++) {
       if( P.Array_Decays[k][2] != 0 && P.Array_Decays[k][5] > minBR && P.Array_Decays[k][2] > 0 && P.Array_Decays[k][3] == 2) {
     	 fout << left << setw(6) << " " << setprecision(8) << P.Array_Decays[k][5] << setprecision(0) << setw(6) << fixed << " " << P.Array_Decays[k][3] << setw(4) << " ";  printRowPDG(fout, P.Array_Decays[k][0]); fout << "    "; printRowPDG(fout,P.Array_Decays[k][1]); fout << "    "; fout << left << setprecision(0) << setw(30) << P.Array_Comments[k] << "    " << scientific << setprecision(8) << setw(18) << P.Array_Decays[k][2] << endl;
       }
     }
     fout << "#" << endl; 
   }
   if (P.three_width/P.total_width > minBR) {
     fout << left << setw(6) << "# " << setw(20) << "BR" << setw(6) << "NDA" << setw(12) << "PDG1 " << setw(11) << "PDG2" << setw(12) << "PDG3" << setw(18) << "Comments" << setw(28) << "PW" << endl;
     for (int k=0; k<P.No_of_Decays; k++) {
       if( P.Array_Decays[k][2] != 0 && P.Array_Decays[k][5] > minBR && P.Array_Decays[k][2] > 0 && P.Array_Decays[k][3] == 3) {
     	 fout << left << setw(6) << " " << setprecision(8) << P.Array_Decays[k][5] << setprecision(0) << setw(8) << fixed << " " << P.Array_Decays[k][3] << setw(2) << " ";  printRowPDG(fout, P.Array_Decays[k][0]); fout << "    "; printRowPDG(fout,P.Array_Decays[k][1]); fout << "    "; printRowPDG(fout,P.Array_Decays[k][4]); fout << "   " << left << setprecision(0) << setw(38) << P.Array_Comments[k] << "    "<< setw(18) << scientific << setprecision(8) << P.Array_Decays[k][2] << endl;
       }
     }
     fout << "#" << endl; 
   }
 }



double gpsitildadgauss(double Et) {
  double gpsitildadgauss = 0, pt = 0, squareplus = 0, squareminus = 0, lambda = 0, A = 0;
  A = sqr(m1)+sqr(mq)-2*fabs(m1)*Et;
  pt = pow(sqr(Et)-sqr(mq),0.5);
  // std::cout << "Et = " << Et << std::endl;
  // std::cout << "mq = " << mq <<std::endl;
  // std::cout << "pt = " << pt << std::endl;
  squareplus = A - pow(fabs(m4)+mq,2);
  squareminus = A - pow(fabs(m4)-mq,2);
  if (squareplus <0) ///< this can happen erronesouly at very end of range due to finite precision used, squareplus should actually then be very very small and +ve
    { squareplus = 0;} ///< set to zero to avoid nan problem in lambda, note that given squareplus very very very small anyway here this should not affect the accuracy of the integral
  lambda = sqrt(squareplus*squareminus);
  if (lambda != lambda) {
    throw("problem: nan in lambda in gpsitildadgauss used in 1->3 decays\n");
  }
  gpsitildadgauss = sqr(PI)*fabs(m1)*pt*Et*lambda/(A)*(sqr(m1)-sqr(m4)-2*fabs(m1)*Et)/((A-sqr(m2))*(A-sqr(m3)));


  return gpsitildadgauss;
}


double gphitildadgauss(double Et) {
  double gphitildadgauss = 0, A = 0, Z=0;
  DoubleVector Etbarmaxmin (double m1, double m2, double massq, double Et);
  double Zfunc(double m1, double mq, double m, double Etbarmax, double Etbarmin);
  DoubleVector Etbar(2);
  // std::cout << "gphitildadgauss: " << std::endl;
  // std::cout << "mq = " << mq << std::endl;
  // std::cout << "Et = " << Et << std::endl;
  for (int i=1; i<=2; i++) { Etbar(i) = 0;}
  Etbar = Etbarmaxmin(m1, m4, mq, Et);
  Z = Zfunc(m1, mq, m3, Etbar(1), Etbar(2));
  A = sqr(m1)+sqr(mq)-2*fabs(m1)*Et;
  gphitildadgauss = 0.5*sqr(PI)*fabs(m1)*fabs(m4)/(A-sqr(m2))*(-(Etbar(1)-Etbar(2)) - (sqr(m4)-sqr(mq)+2*Et*fabs(m1)-sqr(m3))/(2*fabs(m1))*log(Z));
  return gphitildadgauss;
}

double gphitildadgausslimit(double Et) {
  //This is gphitilda with the logZ expanded in the limit mqL (=m2=m3) >> rest of masses or equally in the limit of very compressed spectra where (Etbarmax - Etbarmin) << rest of masses (i.e. pt is small cf gluino mass). phiL/R should be positive definite as the integrand for gluino 3 body decays to neutralinos and a quark antiquark pair is proportional to it, yet in very compressed spectra you can find fine cancellations between the -(Etbar(1)-Etbar(2)) and the - (sqr(m4)-sqr(mq)+2*Et*fabs(m1)-sqr(m3))/(2*fabs(m1))*log(Z) pieces in the usual expression in gphitildadgauss. These fine cancellations cause random positive and negative values (due to finite precision), potentially leaving phiL/R overall negative. In this case we therefore call this form here, which is expanded and rearranged to explicitly remove the fine cancellations (the function is always positive) and so allow better evaluation (no numerical precision issues).
  double gphitildadgauss = 0, A = 0, Z=0, pt = 0;
  DoubleVector Etbarmaxmin (double m1, double m2, double massq, double Et);
  double Zfunc(double m1, double mq, double m, double Etbarmax, double Etbarmin);
  DoubleVector Etbar(2);
  // std::cout << "gphitildadgauss: " << std::endl;
  // std::cout << "mq = " << mq << std::endl;
  // std::cout << "Et = " << Et << std::endl;
  for (int i=1; i<=2; i++) { Etbar(i) = 0;}
  Etbar = Etbarmaxmin(m1, m4, mq, Et);
  Z = Zfunc(m1, mq, m3, Etbar(1), Etbar(2));
  A = sqr(m1)+sqr(mq)-2*fabs(m1)*Et;
  // std::cout << "m1 = " << m1 << " m2 = " << m2 << " m4 = " << m4 << " mq = " << mq << std::endl;
  pt = sqrt(sqr(Et)-sqr(mq));
  // std::cout << "Et = " << Et << " pt = " << pt << std::endl;
  if (pt/m1 > 0.1 && m2/m1 < 5) { //The choice of the exact limit is somewhat arbitrary and may be varied by the user here to attempt to avoid any discontinuities arising upon switching to the limiting form in the else statement
    gphitildadgauss = 0.5*sqr(PI)*fabs(m1)*fabs(m4)/(A-sqr(m2))*(-(Etbar(1)-Etbar(2)) - (sqr(m4)-sqr(mq)+2*Et*fabs(m1)-sqr(m3))/(2*fabs(m1))*log(Z));
    // std::cout << "normal" << std::endl;
  }
  else {
    gphitildadgauss = 0.5*sqr(PI)*fabs(m1)*fabs(m4)/(A-sqr(m2))*(Etbar(1)-Etbar(2))*(sqr(m4)-2*sqr(mq)-sqr(m1)+2*Et*fabs(m1)+2*fabs(m1)*Etbar(2))/(sqr(m1)+sqr(mq)-2*fabs(m1)*Etbar(2)-sqr(m2));
    // std::cout << "limit" << std::endl;
  }
  // std::cout << "gphitildadgauss = " << gphitildadgauss << std::endl;
  return gphitildadgauss;
}



double gxsidgauss (double Et)
{
  double gxsidgauss = 0, Z=0, A=0;
  double Zfunc(double m1, double mq, double m, double Etbarmax, double Etbarmin);
  DoubleVector Etbarmaxmin (double m1, double m2, double massq, double Et);
  DoubleVector Etbar(2);
  //fout << "m1 = " << m1 << " m2 = " << m2 << " m3 = " << m3 << " m4 = " << m4 << " mq = " << mq << endl;
  for (int i=1; i<=2; i++) { Etbar(i) = 0;}
  Etbar = Etbarmaxmin(m1, m4, mq, Et);
  Z = Zfunc(m1,mq,m3,Etbar(1),Etbar(2));
  A = sqr(m1)+sqr(mq)-2*fabs(m1)*Et;
  gxsidgauss = 0.5*sqr(PI)/(A-sqr(m2))*((Etbar(1)-Etbar(2))-(sqr(m1)-sqr(mq)-2*fabs(m1)*Et+sqr(m3))/(2*m1)*log(Z));
  return gxsidgauss;
}


double grhodgauss (double Et)
{
  double grhodgauss = 0, Z=0;
  double Zfunc(double m1, double mq, double m, double Etbarmax, double Etbarmin);
  DoubleVector Etbar(2);
  //fout << "m1 = " << m1 << " m2 = " << m2 << " m3 = " << m3 << " m4 = " << m4 << " mq = " << mq << endl;
  for (int i=1; i<=2; i++) { Etbar(i) = 0;}
  DoubleVector Etbarmaxmin (double m1, double m2, double massq, double Et);
  Etbar = Etbarmaxmin(m1, m4, mq, Et);
  Z = Zfunc(m1,mq,m3,Etbar(1),Etbar(2));
  grhodgauss = -sqr(PI)/(2*fabs(m1))*1/(sqr(m1)+sqr(mq)-2*fabs(m1)*Et-sqr(m2))*log(Z);
  return grhodgauss;
}


double gchidgauss (double Et) {
  double gchidgauss = 0, pt=0, A=0, squareplus=0, squareminus=0, lambda=0;
  A = sqr(m1)+sqr(mq)-2*fabs(m1)*Et;
  squareplus = A - pow((fabs(m4) + mq),2);
  squareminus = A - pow((fabs(m4) - mq),2);
  if (squareplus < 0 && fabs(squareplus) < 1e-8) squareplus = 0.;
  lambda = sqrt(squareplus*squareminus);
  if (lambda != lambda) {
    throw("problem: nan in lambda in gchidgauss used in 1->3 decays\n");
  }
  pt = pow(sqr(Et) - sqr(mq),0.5);
  if (pt != pt) {
    throw("problem: nan in pt in gchidgauss used in 1->3 decays\n");
  }
  gchidgauss = sqr(PI)*fabs(m1)*pt*Et*lambda/A*1/((sqr(m1)+sqr(mq)-2*fabs(m1)*Et-sqr(m2))*(sqr(m1)+sqr(mq)-2*fabs(m1)*Et-sqr(m3)));
  return gchidgauss;
}


double gzetadgauss (double Et)
{
  double gzetadgauss = 0, A=0;
  DoubleVector Etbar(2);
  //  DoubleVector Etbarmaxmin(double m1, double m2, double massq, double Et);
  Etbar = Etbarmaxmin(m1, m4, mq, Et);
  A = sqr(m1)+sqr(mq)-2.0*fabs(m1)*Et;

  gzetadgauss = sqr(PI)*(Etbar(1)-Etbar(2))/((A-sqr(m2))*(A-sqr(m3)));
  //  if (printOutNow) cout << "Et=" << Et << " Etbarmaxmin=" << Etbar << " DEBUG: gzetadgauss=" << gzetadgauss << endl;
  return gzetadgauss;
}


double gXdgauss (double Et) {
  double gXdgauss = 0, pt=0, A=0, B=0, squareplus=0, squareminus=0, lambda=0;
  A = sqr(m1)+sqr(mq)-2*fabs(m1)*Et;
  B = sqr(m1)-sqr(m4)-2*fabs(m1)*Et;
  pt = pow(sqr(Et) - sqr(mq),0.5);
  squareplus = A - pow((fabs(m4) + mq),2);
  squareminus = A - pow((fabs(m4) - mq),2);
  if (squareplus < 0 && fabs(squareplus) < 1e-8) squareplus = 0.;
  if (squareplus*squareminus < 0) 
    throw("problem: lambda will give nan in gXdgauss used in 1->3 decays\n");
  
  lambda = sqrt(squareplus * squareminus);
  gXdgauss = 0.5*sqr(PI)*pt*B/A*lambda*1/((sqr(m1)+sqr(mq)-2*fabs(m1)*Et-sqr(m2))*(sqr(m1)+sqr(mq)-2*fabs(m1)*Et-sqr(m3)));
  return gXdgauss;
}


double gYdgauss (double Et)
{
  double gYdgauss = 0, A=0, Z=0;
  double Zfunc(double m1, double mq, double m, double Etbarmax, double Etbarmin);
  DoubleVector Etbar(2);
  for (int i=1; i<=2; i++) { Etbar(i) = 0;}

  Etbar = Etbarmaxmin(m1, m4, mq, Et);
  Z = Zfunc(m1,mq,m3,Etbar(1),Etbar(2));
  A = sqr(m1)+sqr(mq)-2*fabs(m1)*Et;
  gYdgauss = 0.5*sqr(PI)*1/(A-sqr(m2))*((Etbar(1)-Etbar(2))*A + 1/(2*fabs(m1))*(sqr(m1)*sqr(m4)-sqr(m1)*sqr(m3)+pow(mq,4)+2*fabs(m1)*Et*sqr(m3)-sqr(m3)*sqr(mq))*log(Z));
  return gYdgauss;
}


double gchiprimedgauss (double Et)
{
  double gchiprimedgauss = 0, A=0, Z=0;
  DoubleVector Etbar = Etbarmaxmin(m1, m4, mq, Et);
  double Zfunc(double m1, double mq, double m, double Etbarmax, double Etbarmin);
  Z = Zfunc(m1,mq,m2,Etbar(1),Etbar(2));
  A = sqr(m1)+sqr(mq)-2*fabs(m1)*Et;
  gchiprimedgauss = -0.5*sqr(PI)*Et/(A-sqr(m3))*log(Z);
  return gchiprimedgauss;
}



double gG1dgauss(double Et) ///m1 = mgluino, m2 = mstopi, m6 = mtop, m8 = mcharginoj
{
  double gG1dgauss = 0, A=0, pt=0;
  pt = pow(sqr(Et)-pow(m6,2),0.5);
  A = sqr(m1)+pow(m6,2)-2*fabs(m1)*Et;
  gG1dgauss = fabs(m1)*pt*Et*pow((A-pow(m8,2)),2)/(pow(A-sqr(m2),2)*A);
  return gG1dgauss;
}


double gG4dgauss(double Et) ///m1 = mgluino, m2 = mstopi, m4 = msbottomi, m6 = mtop, m7 = mbottom , m8 = mcharginoj
{
  double gG4dgauss = 0, A=0, X=0, Ebbarmax=0, Ebbarmin=0;
  Ebbarmax = Ebbarmaxmin (m1, m6, m7, m8, Et)(1);
  Ebbarmin = Ebbarmaxmin (m1, m6, m7, m8, Et)(2);
  X = Xfunc(m1, m6, m7, m8, m4, Et);
  A = sqr(m1)+pow(m6,2)-2*fabs(m1)*Et;
  gG4dgauss = (m1)*(m8)*1/(A-sqr(m2))*((Ebbarmax-Ebbarmin) - (sqr(m4)+pow(m6,2)-2*Et*fabs(m1)-pow(m8,2))/(2*fabs(m1))*log(X));
  return gG4dgauss;
}


double gG5dgauss(double Et) ///m1 = mgluino, m2 = mstopi, m4 = msbottomi, m6 = mtop, m7 = mbottom , m8 = mcharginoj
{
  double gG5dgauss = 0, A=0, X=0;
  X = Xfunc(m1, m6, m7, m8, m4, Et);
  A = sqr(m1)+pow(m6,2)-2*fabs(m1)*Et;
  
  gG5dgauss = (fabs(m1)/m1)*m6/2*(A-pow(m8,2))/(A-sqr(m2))*log(X);
  return gG5dgauss;
}


double gG6dgauss(double Et) ///m1 = mgluino, m2 = mstopi, m4 = msbottomi, m6 = mtop, m7 = mbottom , m8 = mcharginoj
{
  double gG6dgauss = 0, A=0, X=0, Ebbarmax=0, Ebbarmin=0;
  Ebbarmax = Ebbarmaxmin (m1, m6, m7, m8, Et)(1);
  Ebbarmin = Ebbarmaxmin (m1, m6, m7, m8, Et)(2);
  X = Xfunc(m1, m6, m7, m8, m4, Et);
  A = sqr(m1)+pow(m6,2)-2*fabs(m1)*Et;
  gG6dgauss = 0.5/(A-sqr(m2))*((fabs(m1)*(A-pow(m8,2))-(sqr(m4)-sqr(m1))/fabs(m1)*-A)*log(X) + 2*(-A)*(Ebbarmax-Ebbarmin));
  return gG6dgauss;
}


double gG7dgauss(double Et)
{
  double gG7dgauss = 0, A=0, X=0, Ebbarmax=0, Ebbarmin=0;

  Ebbarmax = Ebbarmaxmin (m1, m6, m7, m8, Et)(1);
  Ebbarmin = Ebbarmaxmin (m1, m6, m7, m8, Et)(2);
  X = Xfunc(m1, m6, m7, m8, m4, Et);
  A = sqr(m1)+pow(m6,2)-2*fabs(m1)*Et;
  gG7dgauss = 0.5*(m8)*m6*1/(A-sqr(m2))*(2*(Ebbarmax-Ebbarmin)-(sqr(m4)-sqr(m1))/fabs(m1)*log(X));
  return gG7dgauss;
}


double gG8dgauss(double Et)
{
  double gG8dgauss = 0, A=0, Ebbarmax=0, Ebbarmin=0;
  Ebbarmax = Ebbarmaxmin (m1, m6, m7, m8, Et)(1);
  Ebbarmin = Ebbarmaxmin (m1, m6, m7, m8, Et)(2);
  A = sqr(m1)+pow(m6,2)-2*fabs(m1)*Et;
  gG8dgauss = ((m1))*m6*(A-pow(m8,2))*(Ebbarmax-Ebbarmin)/((A-sqr(m2))*(A-sqr(m3)));
  return gG8dgauss;
}


double gG2dgauss(double Ebbar)
{ 
  double gG2dgauss = 0, A=0, squareplus=0, squareminus=0, lambda=0;
  A = sqr(m1)+pow(m7,2)-2*fabs(m1)*Ebbar;
  squareplus = A - pow((m8)+m6,2);
  squareminus = A - pow((m8)-m6,2);
  lambda = sqrt(squareplus*squareminus);
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gG2dgauss used in 1->3 decays" );
  } 
  gG2dgauss = fabs(m1)*pow(Ebbar,2)*lambda*(A-pow(m6,2)-pow(m8,2))/(pow(A-sqr(m4),2)*A);
  return gG2dgauss;
}


double gG3dgauss(double Ebbar)
{
  double gG3dgauss = 0, A=0, squareplus=0, squareminus=0, lambda=0;
  A = sqr(m1)+pow(m7,2)-2*fabs(m1)*Ebbar;
  squareplus = A - pow((m8)+m6,2);
  squareminus = A - pow((m8)-m6,2);
  lambda = sqrt(squareplus*squareminus);
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gG3dgauss used in 1->3 decays" );
  } 

  gG3dgauss = pow(Ebbar,2)*lambda*4*fabs(m1)*fabs(m8)*fabs(m6)/(pow(A-sqr(m4),2)*A);
  return gG3dgauss;
}


double gZdgauss(double E) ///m1 = mZi, m4 = mZj, m2 = mstop1, m3 = mstop2, mq = mt, MZboson = mz
{
  double gZdgauss = 0, Bf=0, D=0, F=0, G=0, H=0;
  Bf = pow(1-4*sqr(mq)/(sqr(m1)+sqr(m4)-2*E*fabs(m1)),0.5);
  D = sqr(m1)+sqr(m4)-pow(MZboson,2)-2*E*fabs(m1);
  F = sqr(m1)+sqr(m4)-2*fabs(m1)*fabs(m4);
  G = pow(E,2)+sqr(m4) + Bf/3*(pow(E,2)-sqr(m4));
  H = fabs(m4)*(fabs(m1)/m1)*(sqr(m1)+sqr(m4)-2*sqr(mq));
  gZdgauss = Bf*pow(pow(E,2)-sqr(m4),0.5)/(pow(D,2))*(E*F-fabs(m1)*G+H);

  return gZdgauss;
}


double ghHdgauss (double E) ///m1 = mZi, m4 = mZj, m2 = mstop1, m3 = mstop2, mq = mt, MZboson = mz
{
   double ghHdgauss = 0, Bf=0, D=0, F=0, G=0, Xijh=0, XijH=0, Xjih=0, XjiH=0, H=0;
   Bf = pow(1-4*sqr(mq)/(sqr(m1)+sqr(m4)-2*E*m1),0.5);
   D = sqr(m1)+sqr(m4)-2*E*m1-2*sqr(mq);
   F = sqr(m1)+sqr(m4)-2*m1*E-pow(mh,2);
   G = sqr(m1)+sqr(m4)-2*m1*E-pow(mH,2);
   H = pow(E,2) - sqr(m4);
   
   Xjih = -0.5*(fabs(m1)/m1)*(fabs(m4)/m4)*(NeutMIX(neutralinoj,3)*-sin(alphamix)-NeutMIX(neutralinoj,4)*cos(alphamix))*(g1*-NeutMIX(neutralinoi,2)-g2*-NeutMIX(neutralinoi,1));
   XjiH = -0.5*(fabs(m1)/m1)*(fabs(m4)/m4)*(NeutMIX(neutralinoj,3)*cos(alphamix)-NeutMIX(neutralinoj,4)*sin(alphamix))*(g1*-NeutMIX(neutralinoi,2)-g2*-NeutMIX(neutralinoi,1));
   Xijh = -0.5*(fabs(m1)/m1)*(fabs(m4)/m4)*(NeutMIX(neutralinoi,3)*-sin(alphamix)-NeutMIX(neutralinoi,4)*cos(alphamix))*(g1*-NeutMIX(neutralinoj,2)-g2*-NeutMIX(neutralinoj,1));
   XijH = -0.5*(fabs(m1)/m1)*(fabs(m4)/m4)*(NeutMIX(neutralinoi,3)*cos(alphamix)-NeutMIX(neutralinoi,4)*sin(alphamix))*(g1*-NeutMIX(neutralinoj,2)-g2*-NeutMIX(neutralinoj,1));
   ghHdgauss = fabs(m1)*Bf*pow(H,0.5)*D*(E+fabs(m4)*(fabs(m1)/m1))*pow(((-sin(alphamix)*(Xjih+Xijh))/F + cos(alphamix)*(XjiH+XijH)/G),2);
   
   return ghHdgauss;
}


double gAdgauss(double E) ///m1 = mZi, m4 = mZj, m2 = mstop1, m3 = mstop2, mq = mt
{
  double gAdgauss = 0, Bf=0, D=0, F=0;
  Bf = pow(1-4*sqr(mq)/(sqr(m1)+sqr(m4)-2*E*m1),0.5);
  D = sqr(m1)+sqr(m4)-2*E*m1-2*sqr(mq);
  F = sqr(m1)+sqr(m4)-2*m1*E-pow(mA,2);
  gAdgauss = fabs(m1)*Bf*pow(pow(E,2)-sqr(m4),0.5)*D*(E-fabs(m4)*(fabs(m1)/m1))/pow(F,2);
  return gAdgauss;
}

double gZsfdgauss(double s) ///m2 = msfi where i is the sfermion contribution interference with Z you are considering
{
  double gZsfintegral = 0, EQ=0, Q=0, Qprime=0, musquared=0, B=0, C=0, D=0, F=0, G=0, H=0;
  EQ = (s + sqr(m1)-sqr(m4))/(2*m1);
  Q = pow(pow(EQ,2)-s,0.5);
  Qprime = Q*pow(1-4*sqr(mq)/s,0.5);
  musquared = s + sqr(m2) - sqr(m4)- sqr(mq);
  B = m1*EQ+sqr(m2)-sqr(m1)-s-sqr(mq);
  C = sqr(m2)-sqr(m4)-sqr(mq);
  D = sqr(m2)-sqr(m1)-sqr(mq);
  F = fabs(m1)*fabs(m4)*(s-2*sqr(mq));
  G = m1*(EQ+Qprime)-musquared;
  H = m1*(EQ-Qprime)-musquared;
  gZsfintegral = sqr(PI)/2*1/fabs(m1)*1/(s-pow(MZboson,2))*(-0.5*Qprime*B-1/(4*m1)*(C*D + F)*log(G/H));
  
  return gZsfintegral;
}

double gJdgauss (double s)
{
  double gJdgauss = 0, EQ=0, Q=0, Qprime=0, musquared=0, B=0, G=0, H=0;
  EQ = (s + sqr(m1)-sqr(m4))/(2*m1);
  Q = pow(pow(EQ,2)-s,0.5);
  Qprime = Q*pow(1-4*sqr(mq)/s,0.5);
  musquared = s + sqr(m2) - sqr(m4)- sqr(mq);
  B = s*sqr(m2)-sqr(mq)*(sqr(m1)+sqr(m4))+pow(-1,AorhorH)*fabs(m1)*fabs(m4)*(s-2*sqr(mq));
  G = m1*(EQ+Qprime)-musquared;
  H = m1*(EQ-Qprime)-musquared;
  gJdgauss = 1/(s-pow(mphi,2))*(0.5*s*Qprime + B/(4*m1)*log(G/H));
  return gJdgauss;
}	

double gneutineutjffZ1dgauss(double s) ///m1 = mneuti, m4 = mneutj, mq = mf, MZboson = mZboson
{
  double Z1dgauss = 0, lambda1 = 0, lambda2 = 0, squareplus1 = 0, squareminus1 = 0;
  squareplus1 = s - pow(fabs(m1)+fabs(m4),2);
  squareminus1 = s - pow(fabs(m1)-fabs(m4),2);
  if (squareminus1 > 0 && fabs(squareminus1) < 1e-9) { 
    squareminus1 = 0; ///Set to 0 to avoid numerical precision causing a very small but positive squareminus1 at smax when it should be exactly 0, this can cause problems as it would make lambda1 the sqrt of a negative number (as squareplus1 is negative), avoid this by setting to 0
  }
  if (squareplus1*squareminus1 < 0) {
    throw( "problem: lambda1 will give nan in gneutineutjffZ1dgauss used in 1->3 decays" );
  } 
  lambda1 = pow(squareplus1*squareminus1,0.5);
  lambda2 = pow(s*(s-4*sqr(mq)),0.5);    
  Z1dgauss = 1/(3*pow(s,2)*pow(s-pow(MZboson,2),2))*(-2*pow(s,4) + (sqr(m1) + sqr(m4) + 2*sqr(mq))*pow(s,3) + (pow(sqr(m1)-sqr(m4),2) - 2*(sqr(m1)+sqr(m4))*2*sqr(mq))*pow(s,2) + 2*sqr(mq)*pow(sqr(m1)-sqr(m4),2)*s)*1/s*lambda1*lambda2;
  return Z1dgauss;
}

double gneutineutjffZ2dgauss(double s)
{
  double Z2dgauss = 0, lambda1 = 0, lambda2 = 0, squareplus1 = 0, squareminus1 = 0;
  squareplus1 = s - pow(fabs(m1)+fabs(m4),2);
  squareminus1 = s - pow(fabs(m1)-fabs(m4),2);
  if (squareminus1 > 0 && fabs(squareminus1) < 1e-9) { 
    squareminus1 = 0; ///Set to 0 to avoid numerical precision causing a very small but positive squareminus1 at smax when it should be exactly 0, this can cause problems as it would make lambda1 the sqrt of a negative number (as squareplus1 is negative), avoid this by setting to 0
  }
  if (squareplus1*squareminus1 < 0) {
    throw( "problem: lambda1 will give nan in gneutineutjffZ2dgauss used in 1->3 decays" );
  } 
  lambda1 = pow(squareplus1*squareminus1,0.5);
  lambda2 = pow(s*(s-4*sqr(mq)),0.5);    
  Z2dgauss = 1/(s*pow(s-pow(MZboson,2),2))*lambda1*lambda2*(s-2*sqr(mq));
  return Z2dgauss;
}

double gneutineutjffZ3dgauss(double s)
{
  double Z3dgauss = 0, lambda1 = 0, lambda2 = 0, squareplus1 = 0, squareminus1 = 0;
  squareplus1 = s - pow(fabs(m1)+fabs(m4),2);
  squareminus1 = s - pow(fabs(m1)-fabs(m4),2);
  if (squareminus1 > 0 && fabs(squareminus1) < 1e-9) { 
    squareminus1 = 0; ///Set to 0 to avoid numerical precision causing a very small but positive squareminus1 at smax when it should be exactly 0, this can cause problems as it would make lambda1 the sqrt of a negative number (as squareplus1 is negative), avoid this by setting to 0
  }
  if (squareplus1*squareminus1 < 0) {
    throw( "problem: lambda1 will give nan in gneutineutjffZ3dgauss used in 1->3 decays" );
  } 
  lambda1 = pow(squareplus1*squareminus1,0.5);
  lambda2 = pow(s*(s-4*sqr(mq)),0.5);    
  Z3dgauss = 1/(s*pow(s-pow(MZboson,2),2))*lambda1*lambda2*(-s+sqr(m1)+sqr(m4));
  return Z3dgauss;
}

double gneutineutjffZ4dgauss(double s)
{
  double Z4dgauss = 0, lambda1 = 0, lambda2 = 0, squareplus1 = 0, squareminus1 = 0;
  squareplus1 = s - pow(fabs(m1)+fabs(m4),2);
  squareminus1 = s - pow(fabs(m1)-fabs(m4),2);
  if (squareminus1 > 0 && fabs(squareminus1) < 1e-9) { 
    squareminus1 = 0; ///Set to 0 to avoid numerical precision causing a very small but positive squareminus1 at smax when it should be exactly 0, this can cause problems as it would make lambda1 the sqrt of a negative number (as squareplus1 is negative), avoid this by setting to 0
  }
  if (squareplus1*squareminus1 < 0) {
    throw( "problem: lambda1 will give nan in gneutineutjffZ4dgauss used in 1->3 decays" );
  } 
  lambda1 = pow(squareplus1*squareminus1,0.5);
  lambda2 = pow(s*(s-4*sqr(mq)),0.5);    
  Z4dgauss = 1/(s*pow(s-pow(MZboson,2),2))*lambda1*lambda2;
  return Z4dgauss;
}


double gintegralhdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl
{
  double gintegralhdgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*m1*E;
  gintegralhdgauss = (pow(pow(E,2)-sqr(m4),0.5)*pow(s-4*sqr(mq),0.5))/(pow(s,0.5)*pow(s-pow(mh,2),2))*(E-fabs(m4))*(s-4*sqr(mq));

  return gintegralhdgauss;
}

double gintegralHdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mH = mhiggsH
{
  double gintegralHdgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*m1*E;
  gintegralHdgauss = (pow(pow(E,2)-sqr(m4),0.5)*pow(s-4*sqr(mq),0.5))/(pow(s,0.5)*pow(s-pow(mH,2),2))*(E-fabs(m4))*(s-4*sqr(mq));
   
  return gintegralHdgauss;
}


double gintegralh1dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl
{
  double gintegralh1dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralh1dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)/(pow(s,0.5)*pow(s-pow(mh,2),2));
  return gintegralh1dgauss;
}

double gintegralh2dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl
{
  double gintegralh2dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralh2dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)*(s-2*sqr(mq))/(pow(s,0.5)*pow(s-pow(mh,2),2));

  return gintegralh2dgauss;
}


double gintegralh3dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl
{
  double gintegralh3dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralh3dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)*2*fabs(m1)*E/(pow(s,0.5)*pow(s-pow(mh,2),2));
  return gintegralh3dgauss;
}


double gintegralh4dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl
{
  double gintegralh4dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralh4dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)*2*fabs(m1)*E*(s-2*sqr(mq))/(pow(s,0.5)*pow(s-pow(mh,2),2));
  return gintegralh4dgauss;
}

double gintegralH1dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mH = mhiggsH
{
  double gintegralH1dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralH1dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)/(pow(s,0.5)*pow(s-pow(mH,2),2));
  return gintegralH1dgauss;
}

double gintegralH2dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mH = mhiggsH
{
  double gintegralH2dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralH2dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)*(s-2*sqr(mq))/(pow(s,0.5)*pow(s-pow(mH,2),2));
  return gintegralH2dgauss;
}


double gintegralH3dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mH = mhiggsH
{
  double gintegralH3dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralH3dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)*2*fabs(m1)*E/(pow(s,0.5)*pow(s-pow(mH,2),2));
  return gintegralH3dgauss;
}


double gintegralH4dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl
{
  double gintegralH4dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralH4dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)*2*fabs(m1)*E*(s-2*sqr(mq))/(pow(s,0.5)*pow(s-pow(mH,2),2));
  return gintegralH4dgauss;
}


double gintegralhH1dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl, mH = mhiggsH
{
  double gintegralhH1dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralhH1dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)/(pow(s,0.5)*(s-pow(mh,2))*(s-pow(mH,2)));
  return gintegralhH1dgauss;
}

double gintegralhH2dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl, mH = mhiggsH
{
  double gintegralhH2dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralhH2dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)*(s-2*sqr(mq))/(pow(s,0.5)*(s-pow(mh,2))*(s-pow(mH,2)));
  return gintegralhH2dgauss;
}


double gintegralhH3dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl, mH = mhiggsH
{
  double gintegralhH3dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralhH3dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)*2*fabs(m1)*E/(pow(s,0.5)*(s-pow(mh,2))*(s-pow(mH,2)));
  return gintegralhH3dgauss;
}


double gintegralhH4dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mh = mhiggsl, mH = mhiggsH
{
  double gintegralhH4dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralhH4dgauss = 2*fabs(m1)*pow(s-4*sqr(mq),0.5)*pow(pow(E,2)-sqr(m4),0.5)*2*fabs(m1)*2*fabs(m1)*E*(s-2*sqr(mq))/(pow(s,0.5)*(s-pow(mh,2))*(s-pow(mH,2)));
  return gintegralhH4dgauss;
}

double gintegralA1dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mA = mhiggsA
{
  double gintegralA1dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralA1dgauss = 4*sqr(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(s-4*sqr(mq),0.5)/(pow(s,0.5)*pow(s-pow(mA,2),2));
  return gintegralA1dgauss;
}


double gintegralA2dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mA = mhiggsA
{
  double gintegralA2dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralA2dgauss = 4*sqr(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(s-4*sqr(mq),0.5)*(s-2*sqr(mq))/(pow(s,0.5)*pow(s-pow(mA,2),2));
  return gintegralA2dgauss;
}


double gintegralA3dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mA = mhiggsA
{
  double gintegralA3dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralA3dgauss = 4*sqr(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(s-4*sqr(mq),0.5)*2*fabs(m1)*E/(pow(s,0.5)*pow(s-pow(mA,2),2));
  return gintegralA3dgauss;
}


double gintegralA4dgauss (double E) ///m1 = mZi, m4 = mZj, mq = mt, MZboson = mz, mA = mhiggsA
{
  double gintegralA4dgauss = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegralA4dgauss = 4*sqr(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(s-4*sqr(mq),0.5)*2*fabs(m1)*E*(s-2*sqr(mq))/(pow(s,0.5)*pow(s-pow(mA,2),2));
  return gintegralA4dgauss;
}

double gintegral1Zsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, MZboson = mz, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral1Zsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral1Zsfdgauss used in 1->3 decays" );
  }

  gintegral1Zsfdgauss = 1/(s-pow(MZboson,2))*(-2*fabs(m1)*pow(1-4*sqr(mq)/s,0.5)*pow(pow(E,2)-sqr(m4),0.5) - (sqr(m2)-sqr(mq)+sqr(m4) - 2*fabs(m1)*E)*log(logarg));
  
  return gintegral1Zsfdgauss;
}

double gintegral2Zsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, MZboson = mz, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral2Zsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral2Zsfdgauss used in 1->3 decays" );
  }  
  gintegral2Zsfdgauss = 1/(s-pow(MZboson,2))*(2*fabs(m1)*pow(1-4*sqr(mq)/s,0.5)*pow(pow(E,2)-sqr(m4),0.5) + (sqr(m2)+sqr(m1)- 2*fabs(m1)*E - sqr(mq))*log(logarg));
   
  return gintegral2Zsfdgauss;
}

double gintegral3Zsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, MZboson = mz, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral3Zsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral3Zsfdgauss used in 1->3 decays" );
  }

  gintegral3Zsfdgauss = 1/(s-pow(MZboson,2))*((sqr(m1) + 2*sqr(mq) + sqr(m4) -1.5*sqr(m2) - 0.5*(sqr(mq) + fabs(m1)*E + fabs(m1)*pow(1-4*sqr(mq)/s,0.5)*pow(pow(E,2)-sqr(m4),0.5)))*(sqr(mq) + fabs(m1)*E + fabs(m1)*pow(1-4*sqr(mq)/s,0.5)*pow(pow(E,2)-sqr(m4),0.5)-sqr(m2)) - (sqr(m1) + 2*sqr(mq) + sqr(m4) -1.5*sqr(m2) - 0.5*(sqr(mq) + fabs(m1)*E - fabs(m1)*pow(1-4*sqr(mq)/s,0.5)*pow(pow(E,2)-sqr(m4),0.5)))*(sqr(mq) + fabs(m1)*E - fabs(m1)*pow(1-4*sqr(mq)/s,0.5)*pow(pow(E,2)-sqr(m4),0.5)-sqr(m2)) + (sqr(m1) + sqr(mq) - sqr(m2))*(sqr(m2)-sqr(mq)-sqr(m4))*log(logarg));
   
  return gintegral3Zsfdgauss;
}

double gintegral4Zsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, MZboson = mz, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral4Zsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral4Zsfdgauss used in 1->3 decays" );
  }
  
  gintegral4Zsfdgauss = 1/(s-pow(MZboson,2))*(2*fabs(m1)*pow(1-4*sqr(mq)/s,0.5)*pow(pow(E,2)-sqr(m4),0.5) + (sqr(m2)-sqr(mq)-sqr(m4))*log(logarg));
  return gintegral4Zsfdgauss;
}    

double gintegral5Zsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, MZboson = mz, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral5Zsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral5Zsfdgauss used in 1->3 decays" );
  }
  
  gintegral5Zsfdgauss = -1/(s-pow(MZboson,2))*(2*fabs(m1)*pow(1-4*sqr(mq)/s,0.5)*pow(pow(E,2)-sqr(m4),0.5) + (sqr(m2)-sqr(mq)-sqr(m1))*log(logarg));
  return gintegral5Zsfdgauss;
}

double gintegral6Zsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, MZboson = mz, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral6Zsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral6Zsfdgauss used in 1->3 decays" );
  }
  
  gintegral6Zsfdgauss = 1/(s-pow(MZboson,2))*(s- 2*sqr(mq))*log(logarg);
  return gintegral6Zsfdgauss;
}

double gintegral7Zsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, MZboson = mz, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral7Zsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral7Zsfdgauss used in 1->3 decays" );
  }
  
  gintegral7Zsfdgauss = 1/(s-pow(MZboson,2))*(2*fabs(m1)*E)*log(logarg);
  return gintegral7Zsfdgauss;
}

double gintegral8Zsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, MZboson = mz, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral8Zsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral8Zsfdgauss used in 1->3 decays" );
  }

  gintegral8Zsfdgauss = 1/(s-pow(MZboson,2))*log(logarg);
  return gintegral8Zsfdgauss;
}



double gintegral1hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral1hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral1hsfdgauss used in 1->3 decays" ); 
  }

  gintegral1hsfdgauss = 2/(s-pow(mh,2))*(2*s*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2)*s - sqr(mq)*(sqr(m1)+sqr(m4)))*log(logarg));
  
  return gintegral1hsfdgauss;
}


double gintegral2hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral2hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral2hsfdgauss used in 1->3 decays" ); 
  }

  gintegral2hsfdgauss = -1/(s-pow(mh,2))*(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2) + sqr(m4) -2*fabs(m1)*E - sqr(mq))*log(logarg));
     
  return gintegral2hsfdgauss;
}



double gintegral3hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral3hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral3hsfdgauss used in 1->3 decays" ); 
  }

  gintegral3hsfdgauss = 1/(s-pow(mh,2))*(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2) + sqr(m1) -2*fabs(m1)*E - sqr(mq))*log(logarg));
     
  return gintegral3hsfdgauss;
}



double gintegral4hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral4hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral4hsfdgauss used in 1->3 decays" ); 
  }

  gintegral4hsfdgauss = 1/(s-pow(mh,2))*(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2) - sqr(mq) - sqr(m4))*log(logarg));
     
  return gintegral4hsfdgauss;
}


double gintegral5hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral5hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral5hsfdgauss used in 1->3 decays" ); 
  }

  gintegral5hsfdgauss = -1/(s-pow(mh,2))*(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2) - sqr(mq) - sqr(m1))*log(logarg));
     
  return gintegral5hsfdgauss;
}


double gintegral6hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral6hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral6hsfdgauss used in 1->3 decays" ); 
  }

  gintegral6hsfdgauss = 1/(s-pow(mh,2))*(s-2*sqr(mq))*log(logarg);
     
  return gintegral6hsfdgauss;
}


double gintegral7hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral7hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral7hsfdgauss used in 1->3 decays" ); 
  }

  gintegral7hsfdgauss = 1/(s-pow(mh,2))*(2*fabs(m1)*E)*log(logarg);
     
  return gintegral7hsfdgauss;
}


double gintegral8hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral8hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral8hsfdgauss used in 1->3 decays" ); 
  }

  gintegral8hsfdgauss = 1/(s-pow(mh,2))*log(logarg);
     
  return gintegral8hsfdgauss;
}


double gintegral1Hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral1Hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral1Hsfdgauss used in 1->3 decays" ); 
  }

  gintegral1Hsfdgauss = 2/(s-pow(mH,2))*(2*s*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2)*s - sqr(mq)*(sqr(m1)+sqr(m4)))*log(logarg));
   
   return gintegral1Hsfdgauss;
}


double gintegral2Hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral2Hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral2Hsfdgauss used in 1->3 decays" ); 
  }

  gintegral2Hsfdgauss = -1/(s-pow(mH,2))*(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2) + sqr(m4) -2*fabs(m1)*E - sqr(mq))*log(logarg));
     
  return gintegral2Hsfdgauss;
}



double gintegral3Hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral3Hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral3Hsfdgauss used in 1->3 decays" ); 
  }

  gintegral3Hsfdgauss = 1/(s-pow(mH,2))*(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2) + sqr(m1) -2*fabs(m1)*E - sqr(mq))*log(logarg));
     
  return gintegral3Hsfdgauss;
}



double gintegral4Hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral4Hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral4Hsfdgauss used in 1->3 decay" );
  }

  gintegral4Hsfdgauss = 1/(s-pow(mH,2))*(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2) - sqr(mq) - sqr(m4))*log(logarg));
     
  return gintegral4Hsfdgauss;
}


double gintegral5Hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral5Hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral5Hsfdgauss used in 1->3 decays" ); 
  }

  gintegral5Hsfdgauss = -1/(s-pow(mH,2))*(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5) + (sqr(m2) - sqr(mq) - sqr(m1))*log(logarg));
     
  return gintegral5Hsfdgauss;
}


double gintegral6Hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral6Hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral6Hsfdgauss used in 1->3 decays" ); 
  }

  gintegral6Hsfdgauss = 1/(s-pow(mH,2))*(s-2*sqr(mq))*log(logarg);
     
  return gintegral6Hsfdgauss;
}


double gintegral7Hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral7Hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral7Hsfdgauss used in 1->3 decays" ); 
  }

  gintegral7Hsfdgauss = 1/(s-pow(mH,2))*(2*fabs(m1)*E)*log(logarg);
     
  return gintegral7Hsfdgauss;
}


double gintegral8Hsfdgauss (double E) ///m1 = mZi, m4 = mZj, mq = mf, mhiggsl = mh, m2 = msfi where i is the sfermion index you're considering for that interference (i = 1,2)
{
  double gintegral8Hsfdgauss = 0, logarg=0, EQ = 0, Qprime = 0, musquared = 0, s=0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  musquared = s + sqr(m2) - sqr(m4) - sqr(mq);
  EQ = (s+sqr(m1) - sqr(m4))/(2*fabs(m1));
  Qprime = pow(pow(EQ,2) - s,0.5)*pow(1 - 4*sqr(mq)/s,0.5);
  logarg = (fabs(m1)*(EQ + Qprime)-musquared)/(fabs(m1)*(EQ - Qprime) - musquared); ///argument of the log
  if (logarg < 0) {
    throw( "Problem: will get nan as logarg < 0 in gintegral8Hsfdgauss used in 1->3 decays" ); 
  }

  gintegral8Hsfdgauss = 1/(s-pow(mH,2))*log(logarg);
     
  return gintegral8Hsfdgauss;
}


double gintegral1ZAdgauss(double E) /// m1 = mZi, m4 = mZj, mq = mf, mA = mhiggsA, m2 = msfi, MZboson = mz
{
  double gintegral1ZAdgauss = 0, s = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegral1ZAdgauss = 1/((s-pow(MZboson,2))*(s-pow(mA,2)))*2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5)*(sqr(m1) - fabs(m1)*E);
  if (pow(E,2)-sqr(m4)<0) {
    throw( "problem: pow(E,2)-sqr(m4)< 0 so sqrt gives nan in gintegral1ZAdgauss used in 1->3 decays" );
  }
  if (1-4*sqr(mq)/s<0) {
    throw( "problem: 1-4*sqr(mq)/s< 0 so sqrt gives nan in gintegral1ZAdgauss used in 1->3 decays" );
  }
  return gintegral1ZAdgauss;
}

double gintegral2ZAdgauss(double E) /// m1 = mZi, m4 = mZj, mq = mf, mA = mhiggsA, m2 = msfi, MZboson = mz
{
  double gintegral2ZAdgauss = 0, s = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegral2ZAdgauss = -1/((s-pow(MZboson,2))*(s-pow(mA,2)))*2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5)*(sqr(m4) - fabs(m1)*E);
  if (pow(E,2)-sqr(m4)<0) {
    throw( "problem: pow(E,2)-sqr(m4)< 0 so sqrt gives nan in gintegral2ZAdgauss used in 1->3 decays" );
  }
  if (1-4*sqr(mq)/s<0) {
    throw( "problem: 1-4*sqr(mq)/s< 0 so sqrt gives nan in gintegral2ZAdgauss used in 1->3 decays" );
  }

  return gintegral2ZAdgauss;
}

double gintegral3ZAdgauss(double E) /// m1 = mZi, m4 = mZj, mq = mf, mA = mhiggsA, m2 = msfi, MZboson = mz
{
  double gintegral3ZAdgauss = 0, s = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegral3ZAdgauss = 1/((s-pow(MZboson,2))*(s-pow(mA,2)))*2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5)*(sqr(m1)-fabs(m1)*E);
  if (pow(E,2)-sqr(m4)<0) {
    throw( "problem: pow(E,2)-sqr(m4)< 0 so sqrt gives nan in gintegral3ZAdgauss used in 1->3 decays" );
  }
  if (1-4*sqr(mq)/s<0) {
    throw( "problem: 1-4*sqr(mq)/s< 0 so sqrt gives nan in gintegral3ZAdgauss used in 1->3 decays" );
  }

  return gintegral3ZAdgauss;
}

double gintegral4ZAdgauss(double E) /// m1 = mZi, m4 = mZj, mq = mf, mA = mhiggsA, m2 = msfi, MZboson = mz
{
  double gintegral4ZAdgauss = 0, s = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gintegral4ZAdgauss = -1/((s-pow(MZboson,2))*(s-pow(mA,2)))*2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5)*(sqr(m4) - fabs(m1)*E);
  if (pow(E,2)-sqr(m4)<0) {
    throw( "problem: pow(E,2)-sqr(m4)< 0 so sqrt gives nan in gintegral4ZAdgauss used in 1->3 decays" );
  }
  if (1-4*sqr(mq)/s<0) {
    throw( "problem: 1-4*sqr(mq)/s< 0 so sqrt gives nan in gintegral4ZAdgauss used in 1->3 decays" );
  }

  return gintegral4ZAdgauss;
}

double gneutineutjffgA1dgauss(double E) /// m1 = mZi, m4 = mZj, mq = mf, mA = mhiggsA, MZboson = mz
{
  double gneutineutjffgA1dgauss = 0, s = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gneutineutjffgA1dgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5)/((s-pow(MZboson,2))*(s-pow(mA,2)));
  if (pow(E,2)-sqr(m4)<0) {
    throw( "problem: pow(E,2)-sqr(m4)< 0 so sqrt gives nan in gneutineutjffgA1dgauss used in 1->3 decays" );
  }
  if (1-4*sqr(mq)/s<0) {
    throw( "problem: 1-4*sqr(mq)/s< 0 so sqrt gives nan in gneutineutjffgA1dgauss used in 1->3 decays" );
  }  

  return gneutineutjffgA1dgauss;
}

double gneutineutjffgA2dgauss(double E) /// m1 = mZi, m4 = mZj, mq = mf, mA = mhiggsA, MZboson = mz
{
  double gneutineutjffgA2dgauss = 0, s = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gneutineutjffgA2dgauss = 2*fabs(m1)*(s-2*sqr(mq))*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5)/((s-pow(MZboson,2))*(s-pow(mA,2)));
  if (pow(E,2)-sqr(m4)<0) {
    throw( "problem: pow(E,2)-sqr(m4)< 0 so sqrt gives nan in gneutineutjffgA2dgauss used in 1->3 decays" );
  }
  if (1-4*sqr(mq)/s<0) {
    throw( "problem: 1-4*sqr(mq)/s< 0 so sqrt gives nan in gneutineutjffgA2dgauss used in 1->3 decays" );
  }  
  
  return gneutineutjffgA2dgauss;
}

double gneutineutjffgA3dgauss(double E) /// m1 = mZi, m4 = mZj, mq = mf, mA = mhiggsA, MZboson = mz
{
  double gneutineutjffgA3dgauss = 0, s = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gneutineutjffgA3dgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5)*2*fabs(m1)*E/((s-pow(MZboson,2))*(s-pow(mA,2)));
  if (pow(E,2)-sqr(m4)<0) {
    throw( "problem: pow(E,2)-sqr(m4)< 0 so sqrt gives nan in gneutineutjffgA3dgauss used in 1->3 decays" );
  }
  if (1-4*sqr(mq)/s<0) {
    throw( "problem: 1-4*sqr(mq)/s< 0 so sqrt gives nan in gneutineutjffgA3dgauss used in 1->3 decays" );
  }  
  
  return gneutineutjffgA3dgauss;
}

double gneutineutjffgA4dgauss(double E) /// m1 = mZi, m4 = mZj, mq = mf, mA = mhiggsA, MZboson = mz
{
  double gneutineutjffgA4dgauss = 0, s = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  gneutineutjffgA4dgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*pow(1-4*sqr(mq)/s,0.5)*(s-2*sqr(mq))*2*fabs(m1)*E/((s-pow(MZboson,2))*(s-pow(mA,2)));
  if (pow(E,2)-sqr(m4)<0) {
    throw( "problem: pow(E,2)-sqr(m4)< 0 so sqrt gives nan in gneutineutjffgA4dgauss used in 1->3 decays" );
  }
  if (1-4*sqr(mq)/s<0) {
    throw( "problem: 1-4*sqr(mq)/s< 0 so sqrt gives nan in gneutineutjffgA4dgauss used in 1->3 decays" );
  } 
  
  
  return gneutineutjffgA4dgauss;
}


double gneuticharjffpW1dgauss(double E) ///m1 = mZi, m2 = mWj, m3 = mf, m4 = mfp, MWboson = mw
{
  double gneuticharjffpW1dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW1dgauss used in 1->3 decays" );
  } 

  gneuticharjffpW1dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)*(-2*pow(s,4) + (sqr(m1) + sqr(m2) + sqr(m3) + sqr(m4))*pow(s,3) + (pow(sqr(m1)-sqr(m2),2) + pow(sqr(m3)-sqr(m4),2) - 2*(sqr(m1)+sqr(m2))*(sqr(m3)+sqr(m4)))*pow(s,2) + ((sqr(m1)+sqr(m2))*pow(sqr(m3)-sqr(m4),2) + (sqr(m3)+sqr(m4))*pow(sqr(m1)-sqr(m2),2))*s - 2*pow(sqr(m1)-sqr(m2),2)*pow(sqr(m3)-sqr(m4),2))*1/(3*pow(s,2))*1/(pow(s-pow(MWboson,2),2));
  return gneuticharjffpW1dgauss;
}

double gneuticharjffpW2dgauss(double E) ///m1 = mZi, m2 = mWj, m3 = mf, m4 = mfp, MWboson = mw
{
  double gneuticharjffpW2dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW2dgauss used in 1->3 decays" );
  }   
  gneuticharjffpW2dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)*(s-sqr(m3)-sqr(m4))/(pow(s-pow(MWboson,2),2));

  return gneuticharjffpW2dgauss;
}


double gneuticharjffpHpm1dgauss(double E) ///m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mHpm; //sometimes mass order differs as called in different cases
{
  double gneuticharjffpHpm1dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+fabs(m4),2);
  squareminus = s - pow(fabs(m3)-fabs(m4),2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpHpm1dgauss used in 1->3 decays" );
  } 

  gneuticharjffpHpm1dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)/(pow(s-pow(m5,2),2));
  return gneuticharjffpHpm1dgauss;
}

double gneuticharjffpHpm2dgauss(double E) ///m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mHpm; //sometimes mass order differs as called in different cases
{
  double gneuticharjffpHpm2dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+fabs(m4),2);
  squareminus = s - pow(fabs(m3)-fabs(m4),2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpHpm2dgauss used in 1->3 decays" );
  } 

  gneuticharjffpHpm2dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)*(s-sqr(m3)-sqr(m4))/(pow(s-pow(m5,2),2));
  return gneuticharjffpHpm2dgauss;
}

double gneuticharjffpHpm3dgauss(double E) ///m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mHpm; //sometimes mass order differs as called in different cases
{
  double gneuticharjffpHpm3dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+fabs(m4),2); 
  squareminus = s - pow(fabs(m3)-fabs(m4),2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpHpm3dgauss used in 1->3 decays" );
  } 

  gneuticharjffpHpm3dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)*2*fabs(m1)*E/(pow(s-pow(m5,2),2));
  return gneuticharjffpHpm3dgauss;
}

double gneuticharjffpHpm4dgauss(double E) ///m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mHpm; //sometimes mass order differs as called in different cases
{
  double gneuticharjffpHpm4dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+fabs(m4),2);
  squareminus = s - pow(fabs(m3)-fabs(m4),2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
      cout.precision(10);
      throw( "problem: lambda will give nan in gneuticharjffpHpm4dgauss used in 1->3 decays" );
  } 

  gneuticharjffpHpm4dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)*2*fabs(m1)*E*(s-sqr(m3)-sqr(m4))/(pow(s-pow(m5,2),2));

  return gneuticharjffpHpm4dgauss;
}


double gneuticharjffp1sf1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msf1, m6 = msf2
{
  double gneuticharjffp1sf1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m4,2); //inserted fabs
  squareminus = s - pow(fabs(m3)-m4,2); //inserted fabs
  
  if (squareplus < 0 && fabs(squareplus/s) < 1e-1) {squareplus = 0;} ///Avoid small negative values of squareplus at upper boundary where theoretically s = 0 exactly.
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp1sf1sf2dgauss used in 1->3 decays" );
  } 

  gneuticharjffp1sf1sf2dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)/((s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp1sf1sf2dgauss;
}

double gneuticharjffp2sf1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msf1, m6 = msf2
{
  double gneuticharjffp2sf1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m4,2); //inserted fabs
  squareminus = s - pow(fabs(m3)-m4,2); //inserted fabs

  if (squareplus < 0 && fabs(squareplus/s) < 1e-1) {squareplus = 0;} ///Avoid small negative values of squareplus at upper boundary where theoretically s = 0 exactly.
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp2sf1sf2dgauss used in 1->3 decays" );
  } 

  gneuticharjffp2sf1sf2dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)*(s-sqr(m3)-sqr(m4))/((s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp2sf1sf2dgauss;
}
 
double gneuticharjffp3sf1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msf1, m6 = msf2
{
  double gneuticharjffp3sf1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m4,2); //inserted fabs
  squareminus = s - pow(fabs(m3)-m4,2); //inserted fabs

  if (squareplus < 0 && fabs(squareplus/s) < 1e-1) {squareplus = 0;} ///Avoid small negative values of squareplus at upper boundary where theoretically s = 0 exactly.
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp3sf1sf2dgauss used in 1->3 decays" );
  } 

  gneuticharjffp3sf1sf2dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)*2*fabs(m1)*E/((s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp3sf1sf2dgauss;
} 

double gneuticharjffp4sf1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msf1, m6 = msf2
{
  double gneuticharjffp4sf1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m4,2); //inserted fabs
  squareminus = s - pow(fabs(m3)-m4,2); //inserted fabs

  if (squareplus < 0 && fabs(squareplus/s) < 1e-1) {squareplus = 0;} ///Avoid small negative values of squareplus at upper boundary where theoretically s = 0 exactly.
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp4sf1sf2dgauss used in 1->3 decays" );
  } 

  gneuticharjffp4sf1sf2dgauss = 2*fabs(m1)/s*sqrt(squareplus*squareminus)*pow(pow(E,2)-sqr(m2),0.5)*2*fabs(m1)*E*(s-sqr(m3)-sqr(m4))/((s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp4sf1sf2dgauss;
}
 

double gneuticharjffp1sfp1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf2
{
  double gneuticharjffp1sfp1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0, numerator = 0, denominator = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m2,2);
  squareminus = s - pow(fabs(m3)-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp1sfp1sf2dgauss used in 1->3 decays" );
  } 

  numerator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s + 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  denominator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s - 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  Z = numerator/denominator;
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp1sfp1sf2dgauss used in 1->3 decays" );
  }

  gneuticharjffp1sfp1sf2dgauss = (2*(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus) + (pow(m6,2)*s - sqr(m1)*sqr(m3) - sqr(m2)*sqr(m4))*log(Z)))/(s-pow(m5,2));
  return gneuticharjffp1sfp1sf2dgauss;
}

double gneuticharjffp2sfp1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf2
{
  double gneuticharjffp2sfp1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0, numerator = 0, denominator = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m2,2); //inserted fabs
  squareminus = s - pow(fabs(m3)-m2,2); //inserted fabs
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp2sfp1sf2dgauss used in 1->3 decays" );
  } 
  numerator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s + 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  denominator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s - 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  Z = numerator/denominator;
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp2sfp1sf2dgauss used in 1->3 decays" );
  }
  
  gneuticharjffp2sfp1sf2dgauss = -(2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s + (pow(m6,2)-2*fabs(m1)*E + sqr(m4)-sqr(m3))*log(Z))/(s-pow(m5,2));
  return gneuticharjffp2sfp1sf2dgauss;
}

double gneuticharjffp3sfp1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf2
{
  double gneuticharjffp3sfp1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0, numerator = 0, denominator = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m2,2);
  squareminus = s - pow(fabs(m3)-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp3sfp1sf2dgauss used in 1->3 decays" );
  } 
  numerator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s + 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  denominator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s - 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  Z = numerator/denominator;
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp3sfp1sf2dgauss used in 1->3 decays" );
  }
  gneuticharjffp3sfp1sf2dgauss = (2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s + (pow(m6,2)+sqr(m1)-2*fabs(m1)*E-sqr(m2))*log(Z))/(s-pow(m5,2));
  return gneuticharjffp3sfp1sf2dgauss;
}

double gneuticharjffp4sfp1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf2
{
  double gneuticharjffp4sfp1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0, numerator = 0, denominator = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m2,2); //inserted fabs
  squareminus = s - pow(fabs(m3)- m2,2); //inserted fabs
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp4sfp1sf2dgauss used in 1->3 decays" );
  } 
  numerator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s + 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  denominator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s - 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  Z = numerator/denominator;
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp4sfp1sf2dgauss used in 1->3 decays" );
  }
  gneuticharjffp4sfp1sf2dgauss = (2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s + (pow(m6,2)-sqr(m3)-sqr(m4))*log(Z))/(s-pow(m5,2));
  return gneuticharjffp4sfp1sf2dgauss;
}

double gneuticharjffp5sfp1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf2
{
  double gneuticharjffp5sfp1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0, numerator = 0, denominator = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m2,2); //inserted fabs
  squareminus = s - pow(fabs(m3)-m2,2); //inserted fabs
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp5sfp1sf2dgauss used in 1->3 decays" );
  } 
  numerator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s + 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  denominator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s - 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  Z = numerator/denominator;
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp5sfp1sf2dgauss used in 1->3 decays" );
  }
  gneuticharjffp5sfp1sf2dgauss = (-2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - (pow(m6,2)-sqr(m1)-sqr(m2))*log(Z))/(s-pow(m5,2));
  return gneuticharjffp5sfp1sf2dgauss;
}

double gneuticharjffp6sfp1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf2
{
  double gneuticharjffp6sfp1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0, numerator = 0, denominator = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m2,2); //inserted fabs
  squareminus = s - pow(fabs(m3)-m2,2); //inserted fabs
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp6sfp1sf2dgauss used in 1->3 decays" );
  } 
  numerator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s + 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  denominator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s - 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  Z = numerator/denominator;
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp6sfp1sf2dgauss used in 1->3 decays" );
  }
  gneuticharjffp6sfp1sf2dgauss = log(Z)*(s-sqr(m2)-sqr(m3))/(s-pow(m5,2));
  return gneuticharjffp6sfp1sf2dgauss;
}

double gneuticharjffp7sfp1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf2
{
  double gneuticharjffp7sfp1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0, numerator = 0, denominator = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m2,2); //inserted fabs
  squareminus = s - pow(fabs(m3)-m2,2);  //inserted fabs
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
 if (squareplus*squareminus < 0) {
     cout.precision(10);
     // throw( "squareplus = " << squareplus << std::endl);
     // throw( "squareminus = " << squareminus << std::endl);
     throw( "problem: lambda will give nan in gneuticharjffp7sfp1sf2dgauss used in 1->3 decays" );
  } 

  numerator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s + 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  denominator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s - 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  Z = numerator/denominator;
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp7sfp1sf2dgauss used in 1->3 decays" );
  }
  gneuticharjffp7sfp1sf2dgauss = log(Z)*2*fabs(m1)*E/(s-pow(m5,2));
  return gneuticharjffp7sfp1sf2dgauss;
}

double gneuticharjffp8sfp1sf2dgauss(double E) ///m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf2
{
  double gneuticharjffp8sfp1sf2dgauss = 0, s = 0, squareplus = 0, squareminus = 0, numerator = 0, denominator = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(fabs(m3)+m2,2); //inserted fabs
  squareminus = s - pow(fabs(m3)-m2,2); //inserted fabs
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp8sfp1sf2dgauss used in 1->3 decays" );
  } 
  numerator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s + 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  denominator = 0.5*(sqr(m2)+sqr(m3)+2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s - 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s - 2*pow(m6,2));
  Z = numerator/denominator;
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp8sfp1sf2dgauss used in 1->3 decays" );
  }
  gneuticharjffp8sfp1sf2dgauss = log(Z)/(s-pow(m5,2));
  return gneuticharjffp8sfp1sf2dgauss;
}


double gneuticharjffp1WHpmdgauss(double E) /// m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHpm
{
  double gneuticharjffp1WHpmdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0;
  s= sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp1WHpmdgauss used in 1->3 decays" );
  } 
  A = 2*fabs(m1)*E + sqr(m3) + sqr(m4) - (sqr(m1)-sqr(m2))*(sqr(m3)-sqr(m4))/s;
  B = 2*fabs(m1)/s*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus);
  gneuticharjffp1WHpmdgauss = (-0.5*(A*B) + (sqr(m1)+sqr(m3))*B)/((s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp1WHpmdgauss;
}

double gneuticharjffp2WHpmdgauss(double E) /// m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHpm
{
  double gneuticharjffp2WHpmdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0;
  s= sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp2WHpmdgauss used in 1->3 decays" );
  } 
  A = 2*fabs(m1)*E + sqr(m3) + sqr(m4) - (sqr(m1)-sqr(m2))*(sqr(m3)-sqr(m4))/s;
  B = 2*fabs(m1)/s*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus);
  
  gneuticharjffp2WHpmdgauss = (0.5*(A*B) - (sqr(m2)+sqr(m4))*B)/((s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp2WHpmdgauss;
}

double gneuticharjffp3WHpmdgauss(double E) /// m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHpm
{
  double gneuticharjffp3WHpmdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0;
  s= sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp3WHpmdgauss used in 1->3 decays" );
  } 
  A = 2*fabs(m1)*E + sqr(m3) + sqr(m4) - (sqr(m1)-sqr(m2))*(sqr(m3)-sqr(m4))/s;
  B = 2*fabs(m1)/s*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus);
  gneuticharjffp3WHpmdgauss = (0.5*(A*B)+(sqr(m1)-2*fabs(m1)*E-sqr(m3))*B)/((s-pow(m5,2))*(s-pow(m6,2)));

  return gneuticharjffp3WHpmdgauss;
}

double gneuticharjffp4WHpmdgauss(double E) /// m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHpm
{
  double gneuticharjffp4WHpmdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0;
  s= sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp4WHpmdgauss used in 1->3 decays");
  } 
  A = 2*fabs(m1)*E + sqr(m3) + sqr(m4) - (sqr(m1)-sqr(m2))*(sqr(m3)-sqr(m4))/s;
  B = 2*fabs(m1)/s*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus);
  
  gneuticharjffp4WHpmdgauss = (-0.5*(A*B)-(sqr(m2)-2*fabs(m1)*E-sqr(m4))*B)/((s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp4WHpmdgauss;
}

double gneuticharjffpW1Sfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp
{
  double gneuticharjffpW1Sfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m2+m3,2);
  squareminus = s - pow(m2-m3,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW1Sfpdgauss used in 1->3 decays" );
  } 
  A = sqr(m2) + sqr(m3) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (A + B - 2*pow(m6,2))/(A - B - 2*pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffpW1Sfpdgauss used in 1->3 decays");
  }
  gneuticharjffpW1Sfpdgauss = (-B - (pow(m6,2)+sqr(m4)-2*fabs(m1)*E-sqr(m3))*log(Z))/(s-pow(m5,2));
  return gneuticharjffpW1Sfpdgauss;
}

double gneuticharjffpW2Sfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp
{
  double gneuticharjffpW2Sfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m2+m3,2);
  squareminus = s - pow(m2-m3,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW2Sfpdgauss used in 1->3 decays");
  } 
  A = sqr(m2) + sqr(m3) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (A + B - 2*pow(m6,2))/(A - B - 2*pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffpW2Sfpdgauss used in 1->3 decays" );
  }
  gneuticharjffpW2Sfpdgauss = (B + (pow(m6,2)+sqr(m1)-2*fabs(m1)*E-sqr(m2))*log(Z))/(s-pow(m5,2));
  return gneuticharjffpW2Sfpdgauss;
}

double gneuticharjffpW3Sfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp
{
  double gneuticharjffpW3Sfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m2+m3,2);
  squareminus = s - pow(m2-m3,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW3Sfpdgauss used in 1->3 decays");
  } 
  A = sqr(m2) + sqr(m3) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (A + B - 2*pow(m6,2))/(A - B - 2*pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffpW3Sfpdgauss used in 1->3 decays" );
  }
  gneuticharjffpW3Sfpdgauss = ((sqr(m1)+sqr(m3)+sqr(m2)+sqr(m4)-1.5*pow(m6,2)-0.25*(A+B))*(0.5*(A+B)-pow(m6,2)) - (sqr(m1)+sqr(m3)+sqr(m2)+sqr(m4)-1.5*pow(m6,2)-0.25*(A-B))*(0.5*(A-B)-pow(m6,2)) + (sqr(m1)+sqr(m2)-pow(m6,2))*(pow(m6,2)-sqr(m3)-sqr(m4))*log(Z))/(s-pow(m5,2));
  return gneuticharjffpW3Sfpdgauss;
}


double gneuticharjffpW4Sfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp
{
  double gneuticharjffpW4Sfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m2+m3,2);
  squareminus = s - pow(m2-m3,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW4Sfpdgauss used in 1->3 decays");
  } 
  A = sqr(m2) + sqr(m3) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (A + B - 2*pow(m6,2))/(A - B - 2*pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffpW4Sfpdgauss used in 1->3 decays" );
  }
  gneuticharjffpW4Sfpdgauss = (B + (pow(m6,2)-sqr(m3)-sqr(m4))*log(Z))/(s-pow(m5,2));
  return gneuticharjffpW4Sfpdgauss;
}

double gneuticharjffpW5Sfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp
{
  double gneuticharjffpW5Sfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m2+m3,2);
  squareminus = s - pow(m2-m3,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW5Sfpdgauss used in 1->3 decays");
  } 
  A = sqr(m2) + sqr(m3) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (A + B - 2*pow(m6,2))/(A - B - 2*pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffpW5Sfpdgauss used in 1->3 decays" );
  }  
  gneuticharjffpW5Sfpdgauss = (-B - (pow(m6,2)-sqr(m1)-sqr(m2))*log(Z))/(s-pow(m5,2));
  return gneuticharjffpW5Sfpdgauss;
}

double gneuticharjffpW6Sfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp
{
  double gneuticharjffpW6Sfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m2+m3,2);
  squareminus = s - pow(m2-m3,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW6Sfpdgauss used in 1->3 decays");
  } 
  A = sqr(m2) + sqr(m3) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (A + B - 2*pow(m6,2))/(A - B - 2*pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffpW6Sfpdgauss used in 1->3 decays" );
  }
  gneuticharjffpW6Sfpdgauss = (s-sqr(m2)-sqr(m3))*log(Z)/(s-pow(m5,2));
  return gneuticharjffpW6Sfpdgauss;
}

double gneuticharjffpW7Sfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp
{
  double gneuticharjffpW7Sfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m2+m3,2);
  squareminus = s - pow(m2-m3,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW7Sfpdgauss used in 1->3 decays");
  } 
  A = sqr(m2) + sqr(m3) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (A + B - 2*pow(m6,2))/(A - B - 2*pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffpW7Sfpdgauss used in 1->3 decays" );
  }  
  gneuticharjffpW7Sfpdgauss = 2*fabs(m1)*E*log(Z)/(s-pow(m5,2));
  return gneuticharjffpW7Sfpdgauss;
}

double gneuticharjffpW8Sfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp
{
  double gneuticharjffpW8Sfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m2+m3,2);
  squareminus = s - pow(m2-m3,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpW8Sfpdgauss used in 1->3 decays");
  } 
  A = sqr(m2) + sqr(m3) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (A + B - 2*pow(m6,2))/(A - B - 2*pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffpW8Sfpdgauss used in 1->3 decays" );
  }
  gneuticharjffpW8Sfpdgauss = log(Z)/(s-pow(m5,2));
  return gneuticharjffpW8Sfpdgauss;
}

double gneuticharjffpHg1dgauss(double E) ///m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHP
{
  double gneuticharjffpHg1dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpHg1dgauss used in 1->3 decays" );
  } 
  gneuticharjffpHg1dgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus)/(s*(s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffpHg1dgauss;
}

double gneuticharjffpHg2dgauss(double E) ///m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHP
{
  double gneuticharjffpHg2dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpHg2dgauss used in 1->3 decays" );
  }   
  gneuticharjffpHg2dgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus)*(s-sqr(m3)-sqr(m4))/(s*(s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffpHg2dgauss;
}

double gneuticharjffpHg3dgauss(double E) ///m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHP
{
  double gneuticharjffpHg3dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpHg3dgauss used in 1->3 decays" );
  } 
  gneuticharjffpHg3dgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus)*2*fabs(m1)*E/(s*(s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffpHg3dgauss;
}

double gneuticharjffpHg4dgauss(double E) ///m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHP
{
  double gneuticharjffpHg4dgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m4,2);
  squareminus = s - pow(m3-m4,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffpHg4dgauss used in 1->3 decays" );
  } 
  gneuticharjffpHg4dgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus)*2*fabs(m1)*E*(s-sqr(m3)-sqr(m4))/(s*(s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffpHg4dgauss;
}


double gneuticharjffp1gsfpdgauss(double E) ///m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1
{
  double gneuticharjffp1gsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m2,2);
  squareminus = s - pow(m3-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp1gsfpdgauss used in 1->3 decays" );
  } 
  A = sqr(m3) + sqr(m2) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (0.5*(A+B)-pow(m6,2))/(0.5*(A-B)-pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp1gsfpdgauss used in 1->3 decays" );
  }
  gneuticharjffp1gsfpdgauss = 2*(s*B + (pow(m6,2)*s - sqr(m1)*sqr(m3) - sqr(m2)*sqr(m4))*log(Z))/(s-pow(m5,2));
  return gneuticharjffp1gsfpdgauss;
}

double gneuticharjffp2gsfpdgauss(double E) ///m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1
{
  double gneuticharjffp2gsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m2,2);
  squareminus = s - pow(m3-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp2gsfpdgauss used in 1->3 decays" );
  } 
  A = sqr(m3) + sqr(m2) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (0.5*(A+B)-pow(m6,2))/(0.5*(A-B)-pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp2gsfpdgauss used in 1->3 decays" );
  }
  gneuticharjffp2gsfpdgauss = (-B - (pow(m6,2) + sqr(m4) - 2*fabs(m1)*E - sqr(m3))*log(Z))/(s-pow(m5,2));
  return gneuticharjffp2gsfpdgauss;
}

double gneuticharjffp3gsfpdgauss(double E) ///m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1
{
  double gneuticharjffp3gsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m2,2);
  squareminus = s - pow(m3-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp3gsfpdgauss used in 1->3 decays" );
  } 
  A = sqr(m3) + sqr(m2) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (0.5*(A+B)-pow(m6,2))/(0.5*(A-B)-pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp3gsfpdgauss used in 1->3 decays" );
  }
  gneuticharjffp3gsfpdgauss = (B + (pow(m6,2)+sqr(m1)-2*fabs(m1)*E - sqr(m2))*log(Z))/(s-pow(m5,2));
  return gneuticharjffp3gsfpdgauss;
}

double gneuticharjffp4gsfpdgauss(double E) ///m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1
{
  double gneuticharjffp4gsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m2,2);
  squareminus = s - pow(m3-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp4gsfpdgauss used in 1->3 decays" );
  } 
  A = sqr(m3) + sqr(m2) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (0.5*(A+B)-pow(m6,2))/(0.5*(A-B)-pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp4gsfpdgauss used in 1->3 decays" );
  }
  gneuticharjffp4gsfpdgauss = (B + (pow(m6,2)-sqr(m3)-sqr(m4))*log(Z))/(s-pow(m5,2));
  return gneuticharjffp4gsfpdgauss;
}

double gneuticharjffp5gsfpdgauss(double E) ///m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1
{
  double gneuticharjffp5gsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m2,2);
  squareminus = s - pow(m3-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp5gsfpdgauss used in 1->3 decays" );
  } 
  A = sqr(m3) + sqr(m2) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (0.5*(A+B)-pow(m6,2))/(0.5*(A-B)-pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp5gsfpdgauss used in 1->3 decays" );
  }
  gneuticharjffp5gsfpdgauss = (-B-(pow(m6,2)-sqr(m1)-sqr(m2))*log(Z))/(s-pow(m5,2));
  return gneuticharjffp5gsfpdgauss;
}

double gneuticharjffp6gsfpdgauss(double E) ///m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1
{
  double gneuticharjffp6gsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m2,2);
  squareminus = s - pow(m3-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp6gsfpdgauss used in 1->3 decays" );
  } 
  A = sqr(m3) + sqr(m2) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (0.5*(A+B)-pow(m6,2))/(0.5*(A-B)-pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp6gsfpdgauss used in 1->3 decays" );
  }
  gneuticharjffp6gsfpdgauss = (s-sqr(m2)-sqr(m3))*log(Z)/(s-pow(m5,2));
  return gneuticharjffp6gsfpdgauss;
}

double gneuticharjffp7gsfpdgauss(double E) ///m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1
{
  double gneuticharjffp7gsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m2,2);
  squareminus = s - pow(m3-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp7gsfpdgauss used in 1->3 decays" );
  } 
  A = sqr(m3) + sqr(m2) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (0.5*(A+B)-pow(m6,2))/(0.5*(A-B)-pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp5gsfpdgauss used in 1->3 decays" );
  }
  gneuticharjffp7gsfpdgauss = 2*fabs(m1)*E*log(Z)/(s-pow(m5,2));
  return gneuticharjffp7gsfpdgauss;
}

double gneuticharjffp8gsfpdgauss(double E) ///m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1
{
  double gneuticharjffp8gsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0, A = 0, B = 0, Z = 0;
  s = sqr(m1) + sqr(m4) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+m2,2);
  squareminus = s - pow(m3-m2,2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
    throw( "problem: lambda will give nan in gneuticharjffp8gsfpdgauss used in 1->3 decays" );
  } 
  A = sqr(m3) + sqr(m2) + 2*fabs(m1)*E + (sqr(m1)-sqr(m4))*(sqr(m3)-sqr(m2))/s;
  B = 2*fabs(m1)*pow(pow(E,2)-sqr(m4),0.5)*sqrt(squareplus*squareminus)/s;
  Z = (0.5*(A+B)-pow(m6,2))/(0.5*(A-B)-pow(m6,2));
  if (Z < 0) {
    throw( "problem: Z<0 will give log(Z) as nan in gneuticharjffp8gsfpdgauss used in 1->3 decays" );
  }
  gneuticharjffp8gsfpdgauss = log(Z)/(s-pow(m5,2));
  return gneuticharjffp8gsfpdgauss;
}

double gneuticharjffp1sfpsfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = -mf, m4 = charginoj, m5 = msfp1, m6 = msfp2
{
  double gneuticharjffp1sfpsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+fabs(m4),2);
  squareminus = s - pow(m3-fabs(m4),2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
      cout.precision(10);
      throw( "problem: lambda will give nan in gneuticharjffp1sfpsfpdgauss used in 1->3 decays" );
  } 
  gneuticharjffp1sfpsfpdgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus)/(s*(s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp1sfpsfpdgauss;
}

double gneuticharjffp2sfpsfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = -mf, m4 = charginoj, m5 = msfp1, m6 = msfp2
{
  double gneuticharjffp2sfpsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+fabs(m4),2);
  squareminus = s - pow(m3-fabs(m4),2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
      cout.precision(10);
      throw( "problem: lambda will give nan in gneuticharjffp2sfpsfpdgauss used in 1->3 decays" );
  } 
  gneuticharjffp2sfpsfpdgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus)*(s-sqr(m3)-sqr(m4))/(s*(s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp2sfpsfpdgauss;
}

double gneuticharjffp3sfpsfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = -mf, m4 = charginoj, m5 = msfp1, m6 = msfp2
{
  double gneuticharjffp3sfpsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+fabs(m4),2);
  squareminus = s - pow(m3-fabs(m4),2);
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
      cout.precision(10);
      throw( "problem: lambda will give nan in gneuticharjffp3sfpsfpdgauss used in 1->3 decays" );
  } 
  gneuticharjffp3sfpsfpdgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus)*2*fabs(m1)*E/(s*(s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp3sfpsfpdgauss;
}

double gneuticharjffp4sfpsfpdgauss(double E) /// m1 = mneutralinoi, m2 = mfp, m3 = -mf, m4 = charginoj, m5 = msfp1, m6 = msfp2
{
  double gneuticharjffp4sfpsfpdgauss = 0, s = 0, squareplus = 0, squareminus = 0;
  s = sqr(m1) + sqr(m2) - 2*fabs(m1)*E;
  squareplus = s - pow(m3+fabs(m4),2); //inserted fabs
  squareminus = s - pow(m3-fabs(m4),2); //inserted fabs
  if (fabs(squareplus) < 1e-6) {
      squareplus = 0.0; //Catch nan due to numerical precision
  }
  if (squareplus*squareminus < 0) {
      cout.precision(10);
      throw( "problem: lambda will give nan in gneuticharjffp4sfpsfpdgauss used in 1->3 decays " );
  } 
  
  gneuticharjffp4sfpsfpdgauss = 2*fabs(m1)*pow(pow(E,2)-sqr(m2),0.5)*sqrt(squareplus*squareminus)*2*fabs(m1)*E*(s-sqr(m3)-sqr(m4))/(s*(s-pow(m5,2))*(s-pow(m6,2)));
  return gneuticharjffp4sfpsfpdgauss;
}


double gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogen (double mgluino, double mneutralino, double msqL, double msqR, double mquark, double g, double gp, DoubleMatrix & mixNeut, double alphas, char uord, int neut, bool onetothree)/// m1 is mgluino, m2 is neutralinoi mass, m3 is sqL mass, m4 is sqR mass, m5 is quark mass but assumed zero in calculation here, just used to check allowed for now; char uord tells us if the quark is u type 'u' or d type 'd', int neut tells us which neutralino it is
{
  double amplitudeW=0, phiL=0, phiR=0, psiL=0, psiR=0, A = 0, B = 0, from = 0, upper = 0;
  int i = neut;
  from = mquark;
  upper = (sqr(mgluino)-2*mquark*fabs(mneutralino)-pow(mneutralino,2))/(2*mgluino);
  if(onetothree == false)
    {
      amplitudeW = 0;
    }
  else if (onetothree == true)
    {  
      if(mgluino < fabs(mneutralino) + 2*mquark) {
      	amplitudeW =0; 
      }
      else if (mgluino > msqL + mquark || mgluino > msqR + mquark) {
      	amplitudeW = 0; /// 1->3 decay not relevant here if one to two decay open
      }
      else {
	if ( uord == 'u') {
	  A = 1/(root2)*(g*-mixNeut(i,2)+gp/3*-mixNeut(i,1));
	  B = 4/(3*root2)*gp*-mixNeut(i,1);
	}
	else if (uord == 'd') {
	  A = 1/(root2)*(-g*-mixNeut(i,2) + gp/3*-mixNeut(i,1));
	  B = -2/(3*root2)*gp*-mixNeut(i,1);
	}
	else {
	  throw("problem: uord must be u or d in gluinoamplitudedecaydgaussneutralinoqqbarfirsttwogen");
	}

	m1 = mgluino; mq = mquark; m4 = mneutralino; m2 = msqL; m3 = msqL;
	
	psiL = dgauss(gpsitildadgauss,from,upper,accuracy)*1/(sqr(PI)*m1);
	phiL = dgauss(gphitildadgauss,from,upper,accuracy)*1/(sqr(PI)*m1);
	m1 = mgluino; mq = mquark; m4 = mneutralino; m2 = msqR; m3 = msqR;
	
	psiR = dgauss(gpsitildadgauss,from,upper,accuracy)*1/(sqr(PI)*m1);
	phiR = dgauss(gphitildadgauss,from,upper,accuracy)*1/(sqr(PI)*m1);

	if ((mneutralino > 0 && mgluino > 0) || (mneutralino < 0 && mgluino < 0)) {
	  amplitudeW = alphas/(8*sqr(PI))*(pow(A,2)*(psiL + phiL) + pow(B,2)*(psiR + phiR));
	}
	else if ((mneutralino < 0 || mgluino < 0) && mneutralino*mgluino<0) {
	  amplitudeW = alphas/(8*sqr(PI))*(pow(A,2)*(psiL - phiL) + pow(B,2)*(psiR - phiR));
	}

      }
    }
  return amplitudeW;
}

double gluinoamplitudedecaydgaussneutralinoqqpbarfirsttwogenlimit (double mgluino, double mneutralino, double msqL, double msqR, double mquark, double g, double gp, DoubleMatrix & mixNeut, double alphas, char uord, int neut, bool onetothree)/// m1 is mgluino, m2 is neutralinoi mass, m3 is sqL mass, m4 is sqR mass, m5 is quark mass but assumed zero in calculation here, just used to check allowed for now; char uord tells us if the quark is u type 'u' or d type 'd', int neut tells us which neutralino it is
{
  double amplitudeW=0, phiL=0, phiR=0, psiL=0, psiR=0, A = 0, B = 0, from = 0, upper = 0;
  int i = neut;
  from = mquark;
  upper = (sqr(mgluino)-2*mquark*fabs(mneutralino)-pow(mneutralino,2))/(2*mgluino);
  // std::cout << "from = " << from << " upper = " << upper << std::endl;
  if(onetothree == false)
    {
      amplitudeW = 0;
    }
  else if (onetothree == true)
    {  
      if(mgluino < fabs(mneutralino) + 2*mquark) {
      	amplitudeW =0; 
      }
      else if (mgluino > msqL + mquark || mgluino > msqR + mquark) {
      	amplitudeW = 0; /// 1->3 decay not relevant here if one to two decay open
      }
      else {
	if ( uord == 'u') {
	  A = 1/(root2)*(g*-mixNeut(i,2)+gp/3*-mixNeut(i,1));
	  B = 4/(3*root2)*gp*-mixNeut(i,1);
	}
	else if (uord == 'd') {
	  A = 1/(root2)*(-g*-mixNeut(i,2) + gp/3*-mixNeut(i,1));
	  B = -2/(3*root2)*gp*-mixNeut(i,1);
	}
	else {
	  throw("problem: uord must be u or d in gluinoamplitudedecaydgaussneutralinoqqbarfirsttwogen");
	}

	m1 = mgluino; mq = mquark; m4 = mneutralino; m2 = msqL; m3 = msqL;
	
	psiL = dgauss(gpsitildadgauss,from,upper,accuracy)*1/(sqr(PI)*m1);

	phiL = dgauss(gphitildadgausslimit,from,upper,accuracy)*1/(sqr(PI)*m1);
	
	m1 = mgluino; mq = mquark; m4 = mneutralino; m2 = msqR; m3 = msqR;
	
	psiR = dgauss(gpsitildadgauss,from,upper,accuracy)*1/(sqr(PI)*m1);
	phiR = dgauss(gphitildadgausslimit,from,upper,accuracy)*1/(sqr(PI)*m1);
	
	// std::cout << "A = " << A << std::endl;
	// std::cout << "B = " << B << std::endl;
	// std::cout << "psiL = " << psiL << std::endl;
	// std::cout << "psiR = " << psiR << std::endl;
	// std::cout << "phiL = " << phiL << std::endl;
	// std::cout << "phiR = " << phiR << std::endl;
	
	if ((mneutralino > 0 && mgluino > 0) || (mneutralino < 0 && mgluino < 0)) {
	  amplitudeW = alphas/(8*sqr(PI))*(pow(A,2)*(psiL + phiL) + pow(B,2)*(psiR + phiR));
	}
	else if ((mneutralino < 0 || mgluino < 0) && mneutralino*mgluino<0) {
	  amplitudeW = alphas/(8*sqr(PI))*(pow(A,2)*(psiL - phiL) + pow(B,2)*(psiR - phiR));
	}

      }
    }
  return amplitudeW;
}


double gluinoamplitudedecaydgaussneutralinottbar (double mgluino, double mst1, double mst2, double mneutralino, double mt, double mWboson, double g, double gp, double thetat, double beta, double alphas, DoubleMatrix & mixNeut, double runmt, int neutralino, bool onetothree, char torb) ///calculates PW for gluino -> neutralino + q qbar pair where q are t 
{
  double amplitudeW=0, Gammast1=0, Gammast2=0, Gammast1st2=0;
  double AtZ=0, BtZ=0, ft=0, from=0, upper=0;

  from = mt;
  upper = (sqr(mgluino)-2*mt*fabs(mneutralino)-pow(mneutralino,2))/(2*mgluino);
  double gs = pow(alphas*4*PI,0.5);

  if (fabs(mgluino) < fabs(mneutralino) + 2*mt || fabs(mgluino)> mst1 + mt || fabs(mgluino) > mst2 + mt || onetothree == false) { amplitudeW = 0;}
  else {
    if (torb == 't') {
      ft = g*runmt/(root2*mWboson*sin(beta));
    }
    
    else if (torb == 'b') {
      ft = g*runmt/(root2*mWboson*cos(beta));
    }
    else {
      throw("problem: torb be t or b in gluinoamplitudedecaydgaussneutralinottbar");
    }

    double psiLLst1=0, chiLLst1=0, phiLLst1=0, rhoLLst1=0, xsiLLst1=0;
    
    m1 = mgluino, mq = mt, m2 = mst1, m3 = mst1, m4 = mneutralino;

    psiLLst1 = dgauss(gpsitildadgauss,from,upper,accuracy);
    phiLLst1 = dgauss(gphitildadgauss,from,upper,accuracy);
    chiLLst1 = dgauss(gchidgauss,from,upper,accuracy);
    rhoLLst1 = dgauss(grhodgauss,from,upper,accuracy);
    xsiLLst1 = dgauss(gxsidgauss,from,upper,accuracy);
 
 
    double zetaL1R1st1=0, XL1R1st1=0;
    zetaL1R1st1 = dgauss(gzetadgauss,from,upper,accuracy);
    XL1R1st1 = dgauss(gXdgauss,from,upper,accuracy);
   
    double YL1R2st1=0;
    YL1R2st1 = dgauss(gYdgauss,from,upper,accuracy);
        
    double psiLLst2=0, chiLLst2=0, phiLLst2=0, rhoLLst2=0, xsiLLst2=0;

    m1 = mgluino, mq = mt, m2 = mst2, m3 = mst2, m4 = mneutralino;
    psiLLst2 = dgauss(gpsitildadgauss,from,upper,accuracy);
    phiLLst2 = dgauss(gphitildadgauss,from,upper,accuracy);
    chiLLst2 = dgauss(gchidgauss,from,upper,accuracy); 
    rhoLLst2 = dgauss(grhodgauss,from,upper,accuracy);
    xsiLLst2 = dgauss(gxsidgauss,from,upper,accuracy);
 
    double zetaL1R1st2=0, XL1R1st2=0;
    zetaL1R1st2 = dgauss(gzetadgauss,from,upper,accuracy);
    XL1R1st2 = dgauss(gXdgauss,from,upper,accuracy);
    
   
    double YL1R2st2=0;
    YL1R2st2 = dgauss(gYdgauss,from,upper,accuracy);

 
    double phitildaLLst1st2=0, rhotildaLLst1st2=0, xsiLLst1st2=0;
    m1 = mgluino, mq = mt, m2 = mst1, m3 = mst2, m4 = mneutralino;
    phitildaLLst1st2 = dgauss(gphitildadgauss,from,upper,accuracy);
    rhotildaLLst1st2 = dgauss(grhodgauss,from,upper,accuracy);
    xsiLLst1st2 = dgauss(gxsidgauss,from,upper,accuracy);

 
    double zetaLRst1st2=0, XLRst1st2=0, YLRst1st2=0, chiprimeLRst1st2=0;
    zetaLRst1st2 = dgauss(gzetadgauss,from,upper,accuracy);
    XLRst1st2 = dgauss(gXdgauss,from,upper,accuracy);
    YLRst1st2 = dgauss(gYdgauss,from,upper,accuracy);
    chiprimeLRst1st2 = dgauss(gchiprimedgauss,from,upper,accuracy);
    
 
 
    Complex ast1alpha1(0.0,0.0), ast1beta1(0.0,0.0), ast2alpha1(0.0,0.0), ast2beta1(0.0,0.0), aAtZ(0.0,0.0), aBtZ(0.0,0.0), aft(0.0,0.0);
    ast1alpha1 = Complex(1.0,2.0);
    double pm = 0;
    
    if (torb == 't') {

      AtZ = g/(root2)*(-mixNeut(neutralino,2)) + gp/(3*root2)*(-mixNeut(neutralino,1));
      BtZ = (4./3)*gp/(root2)*(-mixNeut(neutralino,1));

      if (mneutralino >=0) {
	ast1alpha1 = Complex(AtZ*cos(thetat) - ft*mixNeut(neutralino,4)*sin(thetat),0.0);
	ast1beta1 = Complex(ft*mixNeut(neutralino,4)*cos(thetat) + BtZ*sin(thetat),0.0);
	ast2alpha1 = Complex(AtZ*sin(thetat)+ft*mixNeut(neutralino,4)*cos(thetat),0.0);
	ast2beta1 = Complex(ft*mixNeut(neutralino,4)*sin(thetat)-BtZ*cos(thetat),0.0);
	pm = 1;
      } 
      
      if (mneutralino < 0)
	{
	  ast1alpha1 = Complex(0,-(AtZ*cos(thetat) - ft*mixNeut(neutralino,4)*sin(thetat)));
	  ast1beta1 = Complex(0,-(ft*mixNeut(neutralino,4)*cos(thetat) + BtZ*sin(thetat)));
	  ast2alpha1 = Complex(0,-(AtZ*sin(thetat)+ft*mixNeut(neutralino,4)*cos(thetat)));
	  ast2beta1 = Complex(0,-(ft*mixNeut(neutralino,4)*sin(thetat)-BtZ*cos(thetat)));
	  pm = -1;
	}
      
    }

    else if (torb == 'b') {

      AtZ = g/(root2)*(mixNeut(neutralino,2)) + gp/(3*root2)*(-mixNeut(neutralino,1));
      BtZ = (2./3)*gp/(root2)*(mixNeut(neutralino,1));

      if (mneutralino >=0) {
	ast1alpha1 = Complex(AtZ*cos(thetat) - ft*mixNeut(neutralino,3)*sin(thetat),0.0);
	ast1beta1 = Complex(ft*mixNeut(neutralino,3)*cos(thetat) + BtZ*sin(thetat),0.0);
	ast2alpha1 = Complex(AtZ*sin(thetat)+ft*mixNeut(neutralino,3)*cos(thetat),0.0);
	ast2beta1 = Complex(ft*mixNeut(neutralino,3)*sin(thetat)-BtZ*cos(thetat),0.0);
	pm = 1;
      } 
      
      if (mneutralino < 0)
	{
	  ast1alpha1 = Complex(0,-(AtZ*cos(thetat) - ft*mixNeut(neutralino,3)*sin(thetat)));
	  ast1beta1 = Complex(0,-(ft*mixNeut(neutralino,3)*cos(thetat) + BtZ*sin(thetat)));
	  ast2alpha1 = Complex(0,-(AtZ*sin(thetat)+ft*mixNeut(neutralino,3)*cos(thetat)));
	  ast2beta1 = Complex(0,-(ft*mixNeut(neutralino,3)*sin(thetat)-BtZ*cos(thetat)));
	  pm = -1;
	}
    }
    else {
      throw("problem: torb be t or b in gluinoamplitudedecaydgaussneutralinottbar");
    }    
   
      ///Note the effect of the complex couplings and the pm factors for negative masses cancel out as the couplings always appear in pairs so multiplying them gives an extra minus sign when they are purely imaginary, but this extra minus sign is cancelled out by the extra minus sign in the pm factor. Therefore really the additional minus signs come from the fact I've used the neutralino mass itself throughout my calculation (rather than its absolute value) which therefore naturally introduces additional minus signs.
    
    Complex aGammast1 = Complex(0.0,0.0), aextraGammast1 = Complex(0.0,0.0);

    aGammast1 = (ast1alpha1*ast1alpha1+ast1beta1*ast1beta1)*psiLLst1*pm + 4*pm*mt*mneutralino*chiLLst1*ast1alpha1*ast1beta1 - 4*sin(thetat)*cos(thetat)*(ast1alpha1*ast1alpha1 + ast1beta1*ast1beta1)*mgluino*mt*XL1R1st1*pm - 8*pm*sin(thetat)*cos(thetat)*(ast1alpha1*ast1beta1)*mgluino*mt*mt*mneutralino*zetaL1R1st1 - 2*pm*sin(thetat)*cos(thetat)*ast1alpha1*ast1beta1*YL1R2st1 + pm*(ast1alpha1*ast1alpha1*pow(cos(thetat),2) + ast1beta1*ast1beta1*pow(sin(thetat),2))*phiLLst1 - pm*2*mt*mt*sin(thetat)*cos(thetat)*ast1alpha1*ast1beta1*xsiLLst1 + pm*mgluino*mt*ast1alpha1*ast1beta1*xsiLLst1 - pm*mgluino*mt*ast1alpha1*ast1beta1*pow(mneutralino,2)*rhoLLst1 + pm*mgluino*mt*mt*mneutralino*(pow(sin(thetat),2)*ast1alpha1*ast1alpha1 + pow(cos(thetat),2)*ast1beta1*ast1beta1)*rhoLLst1;
    
    Gammast1 = aGammast1.real();

    double extraGammast1 = 0;
    aextraGammast1 = -pm*mneutralino*mt*sin(thetat)*cos(thetat)*(ast1alpha1*ast1alpha1 + ast1beta1*ast1beta1)*sqr(mgluino)*rhoLLst1 + pm*mneutralino*mt*sin(thetat)*cos(thetat)*(ast1alpha1*ast1alpha1 + ast1beta1*ast1beta1)*xsiLLst1;
    extraGammast1 = aextraGammast1.real();

    Complex aGammast2 = Complex(0.0,0.0), aextraGammast2 = Complex(0.0,0.0);
    aGammast2 = pm*(ast2alpha1*ast2alpha1 + ast2beta1*ast2beta1)*psiLLst2 + 4*pm*mt*mneutralino*chiLLst2*ast2alpha1*ast2beta1 + 4*mgluino*mt*sin(thetat)*cos(thetat)*(ast2alpha1*ast2alpha1 + ast2beta1*ast2beta1)*XL1R1st2*pm + 8*pm*sin(thetat)*cos(thetat)*ast2alpha1*ast2beta1*mgluino*pow(mt,2)*mneutralino*zetaL1R1st2 + pm*2*sin(thetat)*cos(thetat)*ast2alpha1*ast2beta1*YL1R2st2 + pm*(pow(sin(thetat),2)*ast2alpha1*ast2alpha1 + pow(cos(thetat),2)*ast2beta1*ast2beta1)*phiLLst2 + pm*xsiLLst2*(2*pow(mt,2)*sin(thetat)*cos(thetat)*ast2alpha1*ast2beta1 + mgluino*mt*ast2alpha1*ast2beta1) + pm*rhoLLst2*(-mgluino*mt*pow(mneutralino,2)*ast2alpha1*ast2beta1) + pm*rhoLLst2*mgluino*mneutralino*pow(mt,2)*(pow(cos(thetat),2)*ast2alpha1*ast2alpha1 + pow(sin(thetat),2)*ast2beta1*ast2beta1);
    Gammast2 = aGammast2.real();

    double extraGammast2 = 0;
    // extraGammast2 = sqr(PI)/(8*pow(mgluino*gs,2))*(8*pow(mgluino*gs/PI,2)*(rhoLLst2*sqr(mgluino)*mneutralino*mt*sin(thetat)*cos(thetat)*(pow(st2alpha1,2)+pow(st2beta1,2)) + mneutralino*mt*sin(thetat)*cos(thetat)*(pow(st2alpha1,2) + pow(st2beta1,2))*xsiLLst2)); /// Note no chiprime term whereas Spheno has a chiprime term -> missing term relative to T&B
    aextraGammast2 = pm*(xsiLLst2*-mneutralino*mt*sin(thetat)*cos(thetat)*(ast2alpha1*ast2alpha1 + ast2beta1*ast2beta1) + rhoLLst2*sqr(mgluino)*mneutralino*mt*sin(thetat)*cos(thetat)*(ast2alpha1*ast2alpha1 + ast2beta1*ast2beta1));
    extraGammast2 = aextraGammast2.real();

    Complex aGammast1st2 = Complex(0.0,0.0), aextraGammast1st2 = Complex(0.0,0.0);
    aGammast1st2 = sqr(PI)/(8*pow(mgluino*gs,2))*(32*pow(mgluino*gs/PI,2)*pm*mgluino*mt*(pow(cos(thetat),2)-pow(sin(thetat),2))*(ast1alpha1*ast2alpha1 + ast1beta1*ast2beta1)*XLRst1st2 + 32*pm*pow(mgluino*gs/PI,2)*mgluino*pow(mt,2)*mneutralino*(ast1alpha1*ast2beta1 + ast1beta1*ast2alpha1)*(pow(cos(thetat),2)-pow(sin(thetat),2))*zetaLRst1st2 + pm*16*pow(mgluino*gs/PI,2)*(ast1beta1*ast2alpha1*pow(cos(thetat),2) - pow(sin(thetat),2)*ast2beta1*ast1alpha1)*YLRst1st2 + 16*pm*pow(mgluino*gs/PI,2)*sin(thetat)*cos(thetat)*(ast1alpha1*ast2alpha1 - ast1beta1*ast2beta1)*phitildaLLst1st2 + 16*pm*pow(mgluino*gs/PI,2)*mt*mneutralino*(ast1alpha1*ast2alpha1 - ast1beta1*ast2beta1)*chiprimeLRst1st2 + pm*16*pow(mgluino*gs/PI,2)*pow(mt,2)*(pow(cos(thetat),2)*ast1alpha1*ast2beta1 - pow(sin(thetat),2)*ast1beta1*ast2alpha1)*xsiLLst1st2 + 16*pm*pow(mgluino*gs/PI,2)*mgluino*pow(mt,2)*mneutralino*(ast1beta1*ast2beta1-ast1alpha1*ast2alpha1)*sin(thetat)*cos(thetat)*rhotildaLLst1st2);

    Gammast1st2 =  aGammast1st2.real(); ///Gammast1st2 terms as from T&B but with sign change for X term and also the factor of -(cos^2(thetat)-sin^2(thetat))(st1alpha1*st2alpha1+st1beta1*st2beta1) T&B has in the chiprime term has been changed to (cos^2(thetat)+sin^2(thetat))(st1alpha1*st2alpha1 - st1beta1*st2beta1) = 1*(st1alpha1*st2alpha1 - st1beta1*st2beta1), in order in both cases to agree with SPHENO
    double extraGammast1st2 = 0; 
    aextraGammast1st2 = (sqr(PI)/(8*pow(mgluino*gs,2))*(-32*pm*pow(mgluino*gs/PI,2)*sin(thetat)*cos(thetat)*mgluino*mt*(ast1alpha1*ast2beta1-ast1beta1*ast2alpha1)*chiprimeLRst1st2 -16*pm*pow(mgluino*gs/PI,2)*sin(thetat)*cos(thetat)*mgluino*mt*(ast1alpha1*ast2beta1 - ast1beta1*ast2alpha1)*xsiLLst1st2 + 16*pm*pow(mgluino*gs/PI,2)*mt*mneutralino*(pow(sin(thetat),2)*ast1alpha1*ast2alpha1 - pow(cos(thetat),2)*ast1beta1*ast2beta1)*xsiLLst1st2 + 16*pm*pow(mgluino*gs/PI,2)*pow(mgluino,3)*mt*(ast1alpha1*ast2beta1-ast1beta1*ast2alpha1)*sin(thetat)*cos(thetat)*rhotildaLLst1st2 - 16*pm*pow(mgluino*gs/PI,2)*sqr(mgluino)*mt*mneutralino*(pow(sin(thetat),2)*ast1alpha1*ast2alpha1 - pow(cos(thetat),2)*ast1beta1*ast2beta1)*rhotildaLLst1st2)); ///Extra terms SPHENO has in Gammast1st2 not present in T&B

    extraGammast1st2 = aextraGammast1st2.real();

    amplitudeW = (alphas)/(8*pow(PI,4)*mgluino)*(Gammast1+Gammast2+Gammast1st2+extraGammast1+extraGammast2+extraGammast1st2);
  }

  return amplitudeW;
}




double gluinoamplitudedecaydgausschartbbar (double mgluino, double mst1, double mst2, double msb1, double msb2, double mtop, double mbottom, double mchar, double alphas, double thetat, double thetab, double MWboson, double g, double gp, double gammaL, double gammaR, double beta, double runmt, double runmb, int chargino, bool onetothree)
{
  double Gammast1 = 0, Gammast2 = 0, Gammast1st2 = 0 , Gammasb1 =0, Gammasb2 = 0, Gammast1sb1 = 0, Gammast1sb2 = 0, Gammast2sb1 = 0, Gammast2sb2 = 0, from = 0, upper = 0, fromb = 0, upperb = 0, amplitudeW = 0, sumsquarest1 = 0, sumsquarest2 = 0, sumsquaresb1 = 0, sumsquaresb2 = 0, alphasb1ch = 0, alphasb2ch = 0, alphast1ch = 0, alphast2ch = 0, betasb1ch = 0, betasb2ch = 0, betast1ch = 0, betast2ch = 0;

  if (mgluino > mbottom + msb1 || mgluino > mbottom + msb2 || mgluino > mtop + mst1 || mgluino > mtop + mst2 || mgluino < mtop + mbottom + mchar || onetothree == false) {amplitudeW = 0;}
  else {
    from = mtop;
    upper = (sqr(mgluino)+pow(mtop,2)-pow(fabs(mchar)+mbottom,2))/(2*mgluino);
    fromb = mbottom;
    upperb = (sqr(mgluino)-pow(mtop+fabs(mchar),2))/(2*mgluino);
  
    DoubleVector couplingst(16);
    DoubleVector couplingsb(16);

    for (int i=1; i<=16; i++) {
      couplingst(i) = 0;
      couplingsb(i) = 0;
    }

    // couplingst = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,1);
    // couplingsb = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,2);

    if (chargino == 1) {
      sumsquarest1 = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,1)(1);
      sumsquarest2 = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,1)(5);
      sumsquaresb1 = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,2)(1);
      sumsquaresb2 = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,2)(5);
      alphast1ch = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,1)(9);
      alphast2ch = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,1)(13);
      betast1ch = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,1)(11);
      betast2ch = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,1)(15);
      alphasb1ch = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,2)(9);
      alphasb2ch = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,2)(13);
      betasb1ch = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,2)(11);
      betasb2ch = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,mchar,0,2)(15);
    }

    else if (chargino == 2) {
      sumsquarest1 = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,1)(3);
      sumsquarest2 = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,1)(7);
      sumsquaresb1 = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,2)(3);
      sumsquaresb2 = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,2)(7);
      alphast1ch = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,1)(10);
      alphast2ch = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,1)(14);
      betast1ch = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,1)(12);
      betast2ch = squarkmixcharginocouplings(g,thetat,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,1)(16);
      alphasb1ch = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,2)(10);
      alphasb2ch = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,2)(14);
      betasb1ch = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,2)(12);
      betasb2ch = squarkmixcharginocouplings(g,thetab,beta,gammaL,gammaR,runmt,runmb,MWboson,0,mchar,2)(16);
    }
    else {
      throw("problem: chargino must be 1 or 2 in gluinoamplitudedecaydgausschartbbar");
    }
    double G1st1 = 0, G1st2 = 0, G2sb1 = 0, G2sb2 = 0, G3sb1 = 0, G3sb2 = 0, G4st1sb1 = 0, G4st1sb2 = 0, G4st2sb1 = 0, G4st2sb2 = 0, G5st1sb1 = 0, G5st1sb2 = 0, G5st2sb1 = 0, G5st2sb2 = 0, G6st1sb1 = 0, G6st1sb2 = 0, G6st2sb1 = 0,G6st2sb2 = 0, G7st1sb1 = 0, G7st1sb2 = 0, G7st2sb1 = 0, G7st2sb2 = 0, G8st1 = 0, G8st2 = 0, G8st1st2 = 0;

    m1 = mgluino, m2 = mst1, m3 = mst1, m6 = mtop, m8 = mchar;
    G1st1 = dgauss(gG1dgauss,from,upper,accuracy);
    G8st1 = dgauss(gG8dgauss,from,upper,accuracy);
    
    Gammast1 = sumsquarest1*(G1st1 - sin(2*thetat)*G8st1); ///Different to T&B, they have + sin(2*thetat)*G8, I follow SPheno
    
    m1 = mgluino, m2 = mst2, m3 = mst2, m6 = mtop, m8 = mchar;
    G1st2 = dgauss(gG1dgauss,from,upper,accuracy);
    G8st2 = dgauss(gG8dgauss,from,upper,accuracy);

    Gammast2 = sumsquarest2*(G1st2 + sin(2*thetat)*G8st2); ///Different to T&B, they have - sin(2*thetat)*G8, I follow SPheno
   
    m1 = mgluino, m2 = mst1, m3 = mst2, m6 = mtop, m8 = mchar;
    G8st1st2 = dgauss(gG8dgauss,from,upper,accuracy);

    Gammast1st2 = 2*(alphast1ch*alphast2ch+betast1ch*betast2ch)*cos(2*thetat)*G8st1st2; ///Global minus sign difference cf T&B, I follow SPheno    
    
    m1 = mgluino, m4 = msb1, m6 = mtop, m7 = mbottom, m8 = mchar;
    G2sb1 = dgauss(gG2dgauss,fromb,upperb,accuracy);
    G3sb1 = dgauss(gG3dgauss,fromb,upperb,accuracy);

    Gammasb1 = sumsquaresb1*G2sb1 + alphasb1ch*betasb1ch*G3sb1; ///Different to T&B, they have - alphasb1ch*betasb1ch*G3sb1, I follow SPheno
    
    m1 = mgluino, m4 = msb2, m6 = mtop, m7 = mbottom, m8 = mchar;
    G2sb2 = dgauss(gG2dgauss,fromb,upperb,accuracy);
    G3sb2 = dgauss(gG3dgauss,fromb,upperb,accuracy);

    Gammasb2 = sumsquaresb2*G2sb2 + alphasb2ch*betasb2ch*G3sb2; ///Different to T&B, they have - alphasb2ch*betasb2ch*G3sb2, I follow SPheno
    ///Taken mb -> 0 limit in squared matrix element (but not in phase space) so no sb1sb2 interference term
    
    m1 = mgluino, m2 = mst1, m4 = msb1, m6 = mtop, m7 = mbottom , m8 = mchar;
    G4st1sb1 = dgauss(gG4dgauss,from,upper,accuracy);
    G5st1sb1 = dgauss(gG5dgauss,from,upper,accuracy);
    G6st1sb1 = dgauss(gG6dgauss,from,upper,accuracy);
    G7st1sb1 = dgauss(gG7dgauss,from,upper,accuracy);
    
    Gammast1sb1 = (cos(thetat)*sin(thetab)*alphasb1ch*betast1ch+sin(thetat)*cos(thetab)*betasb1ch*alphast1ch)*G6st1sb1 - (cos(thetat)*cos(thetab)*alphasb1ch*alphast1ch + sin(thetat)*sin(thetab)*betasb1ch*betast1ch)*G4st1sb1 - (cos(thetat)*cos(thetab)*betasb1ch*alphast1ch + sin(thetat)*sin(thetab)*alphasb1ch*betast1ch)*G5st1sb1 + (cos(thetat)*sin(thetab)*betasb1ch*betast1ch + sin(thetat)*cos(thetab)*alphasb1ch*alphast1ch)*G7st1sb1; ///The sign in front of the G7st1sb1 and G5st1sb1 terms has been changed from that given in Baer and Tata in order to agree with SPheno.
    
    m1 = mgluino, m2 = mst1, m4 = msb2, m6 = mtop, m7 = mbottom, m8 = mchar;
    G4st1sb2 = dgauss(gG4dgauss,from,upper,accuracy);
    G5st1sb2 = dgauss(gG5dgauss,from,upper,accuracy);
    G6st1sb2 = dgauss(gG6dgauss,from,upper,accuracy);
    G7st1sb2 = dgauss(gG7dgauss,from,upper,accuracy);

    Gammast1sb2 = (cos(thetat)*-cos(thetab)*alphasb2ch*betast1ch+sin(thetat)*sin(thetab)*betasb2ch*alphast1ch)*G6st1sb2 - (cos(thetat)*sin(thetab)*alphasb2ch*alphast1ch - sin(thetat)*cos(thetab)*betasb2ch*betast1ch)*G4st1sb2 - (cos(thetat)*sin(thetab)*betasb2ch*alphast1ch - sin(thetat)*cos(thetab)*alphasb2ch*betast1ch)*G5st1sb2 + (cos(thetat)*-cos(thetab)*betasb2ch*betast1ch + sin(thetat)*sin(thetab)*alphasb2ch*alphast1ch)*G7st1sb2; ///The sign in front of the G7st1sb1 and G5st1sb1 terms has been changed from that given in Baer and Tata in order to agree with SPheno.

    
    m1 = mgluino, m2 = mst2, m4 = msb1, m6 = mtop, m7 = mbottom, m8 = mchar;
    G4st2sb1 = dgauss(gG4dgauss,from,upper,accuracy);
    G5st2sb1 = dgauss(gG5dgauss,from,upper,accuracy);
    G6st2sb1 = dgauss(gG6dgauss,from,upper,accuracy);
    G7st2sb1 = dgauss(gG7dgauss,from,upper,accuracy);

    Gammast2sb1 = (sin(thetat)*sin(thetab)*alphasb1ch*betast2ch-cos(thetat)*cos(thetab)*betasb1ch*alphast2ch)*G6st2sb1 - (sin(thetat)*cos(thetab)*alphasb1ch*alphast2ch - cos(thetat)*sin(thetab)*betasb1ch*betast2ch)*G4st2sb1 - (sin(thetat)*cos(thetab)*betasb1ch*alphast2ch - cos(thetat)*sin(thetab)*alphasb1ch*betast2ch)*G5st2sb1 + (sin(thetat)*sin(thetab)*betasb1ch*betast2ch - cos(thetat)*cos(thetab)*alphasb1ch*alphast2ch)*G7st2sb1; ///The sign in front of the G7st1sb1 and G5st1sb1 terms has been changed from that given in Baer and Tata in order to agree with SPheno.

    m1 = mgluino, m2 = mst2, m4 = msb2, m6 = mtop, m7 = mbottom, m8 = mchar;
    G4st2sb2 = dgauss(gG4dgauss,from,upper,accuracy);
    G5st2sb2 = dgauss(gG5dgauss,from,upper,accuracy);
    G6st2sb2 = dgauss(gG6dgauss,from,upper,accuracy);
    G7st2sb2 = dgauss(gG7dgauss,from,upper,accuracy);
    
    Gammast2sb2 = (sin(thetat)*-cos(thetab)*alphasb2ch*betast2ch-cos(thetat)*sin(thetab)*betasb2ch*alphast2ch)*G6st2sb2 - (sin(thetat)*sin(thetab)*alphasb2ch*alphast2ch + cos(thetat)*cos(thetab)*betasb2ch*betast2ch)*G4st2sb2 - (sin(thetat)*sin(thetab)*betasb2ch*alphast2ch + cos(thetat)*cos(thetab)*alphasb2ch*betast2ch)*G5st2sb2 + (sin(thetat)*-cos(thetab)*betasb2ch*betast2ch - cos(thetat)*sin(thetab)*alphasb2ch*alphast2ch)*G7st2sb2; ///The sign in front of the G7st1sb1 and G5st1sb1 terms has been changed from that given in Baer and Tata in order to agree with SPheno.

    amplitudeW = alphas/(16*sqr(PI)*mgluino)*(Gammast1 + Gammast2 + Gammast1st2 + Gammasb1 + Gammasb2 + Gammast1sb1 + Gammast2sb1 + Gammast1sb2 + Gammast2sb2);
  }
  
  return amplitudeW;
}


double neutralinoamplitudedecaydgaussneutralinoffbar (double mneutralinoi, double msf1, double msf2, double mZboson, double mhiggsl, double mhiggsH, double mhiggsA, double mneutralinoj, double mf, double alphas, double thetaq, double mWboson, double g, double gp, double alpha, double beta, double runmq, DoubleMatrix & mixNeut, int ineutralino, int jneutralino, bool onetothree, char uordornuorl)
{
  double GammaZ = 0, Gammahsf1 = 0, Gammahsf2 = 0, GammaHsf1 = 0, GammaHsf2 = 0, GammaAsf1 = 0, GammaAsf2 = 0, GammaZsf1 = 0, GammaZsf2 = 0, amplitudeW = 0;

  if (fabs(mneutralinoi) > mf + msf1 || fabs(mneutralinoi) > mf + msf2 || fabs(mneutralinoi) > fabs(mneutralinoj) + mhiggsl || fabs(mneutralinoi) > fabs(mneutralinoj) + mhiggsH || fabs(mneutralinoi) > fabs(mneutralinoj) + mhiggsA || fabs(mneutralinoi) > fabs(mneutralinoj) + mZboson || fabs(mneutralinoi) < fabs(mneutralinoj) + mf + mf || onetothree == false) {
    amplitudeW = 0;
  }
  else {

    double from = 0, to = 0, fromz = 0, toz = 0, fq = 0, AZi = 0, BZi = 0, sf1alpha1Zi = 0, sf1beta1Zi = 0, sf2alpha1Zi = 0, sf2beta1Zi = 0, AZj = 0, BZj = 0, sf1alpha1Zj = 0, sf1beta1Zj = 0, sf2alpha1Zj = 0, sf2beta1Zj = 0, alphaf = 0, betaf = 0, XijA = 0, XjiA = 0, Wij = 0, Xijh = 0, Xjih = 0, XijH = 0, XjiH = 0, Nc = 0, trigofalphah = 0, trigofalphaH = 0, Aq = 0 , goldstoneffcoup = 0;
        
    from = mf;
    to = (pow(mneutralinoi,2) - 2*mf*fabs(mneutralinoj) - pow(mneutralinoj,2))/(2*fabs(mneutralinoi));

    double ri = 0, rj = 0;
    if (mneutralinoi >= 0) { ri = 1;}
    else if (mneutralinoi < 0) { ri = -1;} ///correction factor for negative masses
    if (mneutralinoj >= 0) { rj = 1;}
    else if (mneutralinoj < 0) { rj = -1;} ///correction factor for negative masses
    
    if (uordornuorl == 'u') {
      fq = g*runmq/(root2*mWboson*sin(beta));
      AZi = g/(root2)*(-mixNeut(ineutralino,2)) + gp/(3*root2)*(-mixNeut(ineutralino,1));
      BZi = (4./3)*gp/(root2)*(-mixNeut(ineutralino,1));
      sf1alpha1Zi = AZi*cos(thetaq) - fq*mixNeut(ineutralino,4)*sin(thetaq);
      sf1beta1Zi = fq*mixNeut(ineutralino,4)*cos(thetaq) + BZi*sin(thetaq);
      sf2alpha1Zi = (AZi*sin(thetaq)+fq*mixNeut(ineutralino,4)*cos(thetaq));
      sf2beta1Zi = fq*mixNeut(ineutralino,4)*sin(thetaq)-BZi*cos(thetaq);
      AZj = g/(root2)*(-mixNeut(jneutralino,2)) + gp/(3*root2)*(-mixNeut(jneutralino,1));
      BZj = (4./3)*gp/(root2)*(-mixNeut(jneutralino,1));
      sf1alpha1Zj = AZj*cos(thetaq) - fq*mixNeut(jneutralino,4)*sin(thetaq);
      sf1beta1Zj = fq*mixNeut(jneutralino,4)*cos(thetaq) + BZj*sin(thetaq);
      sf2alpha1Zj = (AZj*sin(thetaq)+fq*mixNeut(jneutralino,4)*cos(thetaq));
      sf2beta1Zj = fq*mixNeut(jneutralino,4)*sin(thetaq)-BZj*cos(thetaq);
      alphaf =-5*gp/(g*12) + 0.25*(g/gp);
      betaf = -0.25*(gp/g + g/gp);
      Nc = 3;
      trigofalphah = cos(alpha);
      trigofalphaH = sin(alpha);
      Aq = g*runmq/(mWboson*tan(beta));
      goldstoneffcoup = -fq*sin(beta)/root2;
    }
    
    else if (uordornuorl == 'd') {
      fq = g*runmq/(root2*mWboson*cos(beta));
      AZi = g/(root2)*(mixNeut(ineutralino,2)) + gp/(3*root2)*(-mixNeut(ineutralino,1));
      BZi = (2./3)*gp/(root2)*(mixNeut(ineutralino,1));
      sf1alpha1Zi = AZi*cos(thetaq) - fq*mixNeut(ineutralino,3)*sin(thetaq);
      sf1beta1Zi = fq*mixNeut(ineutralino,3)*cos(thetaq) + BZi*sin(thetaq);
      sf2alpha1Zi = (AZi*sin(thetaq)+fq*mixNeut(ineutralino,3)*cos(thetaq));
      sf2beta1Zi = fq*mixNeut(ineutralino,3)*sin(thetaq)-BZi*cos(thetaq);
      AZj = g/(root2)*(mixNeut(jneutralino,2)) + gp/(3*root2)*(-mixNeut(jneutralino,1));
      BZj = (2./3)*gp/(root2)*(mixNeut(jneutralino,1));
      sf1alpha1Zj = AZj*cos(thetaq) - fq*mixNeut(jneutralino,3)*sin(thetaq);
      sf1beta1Zj = fq*mixNeut(jneutralino,3)*cos(thetaq) + BZj*sin(thetaq);
      sf2alpha1Zj = (AZj*sin(thetaq)+fq*mixNeut(jneutralino,3)*cos(thetaq));
      sf2beta1Zj = fq*mixNeut(jneutralino,3)*sin(thetaq)-BZj*cos(thetaq);
      alphaf = gp/(g*12) - 0.25*(g/gp);
      betaf = 0.25*(gp/g + g/gp);
      Nc = 3;
      trigofalphah = -sin(alpha);
      trigofalphaH = cos(alpha);
      Aq = g*runmq*tan(beta)/(mWboson);
      goldstoneffcoup = fq*cos(beta)/sqrt(2);
    }

    else if (uordornuorl == 'n') {
      fq = 0;
      AZi = g/(root2)*(-mixNeut(ineutralino,2)) + gp/(root2)*(mixNeut(ineutralino,1));
      BZi = 0;
      sf1alpha1Zi = AZi*cos(thetaq) - fq*mixNeut(ineutralino,4)*sin(thetaq);
      sf1beta1Zi = fq*mixNeut(ineutralino,4)*cos(thetaq) + BZi*sin(thetaq);
      sf2alpha1Zi = (AZi*sin(thetaq)+fq*mixNeut(ineutralino,4)*cos(thetaq));
      sf2beta1Zi = fq*mixNeut(ineutralino,4)*sin(thetaq)-BZi*cos(thetaq);
      AZj = g/(root2)*(-mixNeut(jneutralino,2)) + gp/(root2)*(mixNeut(jneutralino,1));
      BZj = 0;
      sf1alpha1Zj = AZj*cos(thetaq) - fq*mixNeut(jneutralino,4)*sin(thetaq);
      sf1beta1Zj = fq*mixNeut(jneutralino,4)*cos(thetaq) + BZj*sin(thetaq);
      sf2alpha1Zj = (AZj*sin(thetaq)+fq*mixNeut(jneutralino,4)*cos(thetaq));
      sf2beta1Zj = fq*mixNeut(jneutralino,4)*sin(thetaq)-BZj*cos(thetaq);
      alphaf = 0.25*(gp/g + g/gp);
      betaf = -0.25*(gp/g + g/gp);
      Nc = 1;
      trigofalphah = cos(alpha);
      trigofalphaH = sin(alpha);
      Aq = g*runmq/(tan(beta)*mWboson);
      goldstoneffcoup = -fq*sin(beta)/sqrt(2);
    }

    else if (uordornuorl == 'l') {
      fq = g*runmq/(root2*mWboson*cos(beta));
      AZi = g/(root2)*(mixNeut(ineutralino,2)) + gp/(root2)*(mixNeut(ineutralino,1));
      BZi = root2*gp*mixNeut(ineutralino,1);
      sf1alpha1Zi = AZi*cos(thetaq) - fq*mixNeut(ineutralino,3)*sin(thetaq);
      sf1beta1Zi = fq*mixNeut(ineutralino,3)*cos(thetaq) + BZi*sin(thetaq);
      sf2alpha1Zi = (AZi*sin(thetaq)+fq*mixNeut(ineutralino,3)*cos(thetaq));
      sf2beta1Zi = fq*mixNeut(ineutralino,3)*sin(thetaq)-BZi*cos(thetaq);
      AZj = g/(root2)*(mixNeut(jneutralino,2)) + gp/(root2)*(mixNeut(jneutralino,1));
      BZj = root2*gp*mixNeut(jneutralino,1);
      sf1alpha1Zj = AZj*cos(thetaq) - fq*mixNeut(jneutralino,3)*sin(thetaq);
      sf1beta1Zj = fq*mixNeut(jneutralino,3)*cos(thetaq) + BZj*sin(thetaq);
      sf2alpha1Zj = (AZj*sin(thetaq)+fq*mixNeut(jneutralino,3)*cos(thetaq));
      sf2beta1Zj = fq*mixNeut(jneutralino,3)*sin(thetaq)-BZj*cos(thetaq);

      alphaf = 0.75*gp/g - 0.25*g/gp;
      betaf = 0.25*(gp/g + g/gp);
      Nc = 1;
      trigofalphah = -sin(alpha);
      trigofalphaH = cos(alpha);
      Aq = g*runmq*tan(beta)/(mWboson);
      goldstoneffcoup = fq*cos(beta)/sqrt(2);
    }
    else {
      throw("problem: uordornuorl must be u or d or n or l in neutralinoamplitudedecaydgaussneutralinoffbar");
    }

    double YZisf1sf1Zj = 0, YZisf2sf2Zj = 0, psitildaZisf1sf2Zj = 0, phitildaZisf1sf2Zj = 0, YZisf1sf2Zj = 0;

    m1 = mneutralinoi, m2 = msf1, m3 = msf1, m4 = mneutralinoj, mq = mf;
    YZisf1sf1Zj = dgauss(gYdgauss,from,to,accuracy);

    m1 = mneutralinoi, m2 = msf2, m3 = msf2, m4 = mneutralinoj, mq = mf;
    YZisf2sf2Zj = dgauss(gYdgauss,from,to,accuracy);

    m1 = mneutralinoi, m2 = msf1, m3 = msf2, m4 = mneutralinoj, mq = mf;
    psitildaZisf1sf2Zj = dgauss(gpsitildadgauss,from,to,accuracy);
    phitildaZisf1sf2Zj = dgauss(gphitildadgauss,from,to,accuracy);
    YZisf1sf2Zj = dgauss(gYdgauss,from,to,accuracy);

    fromz = fabs(mneutralinoj);
    toz = (pow(mneutralinoi,2) + pow(mneutralinoj,2) - 4*pow(mf,2))/(2*fabs(mneutralinoi));
    
    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, MZboson = mZboson, mh = mhiggsl, mH = mhiggsH, mA = mhiggsA, g1 = g, g2 = gp, NeutMIX = mixNeut, neutralinoi = ineutralino, neutralinoj = jneutralino, alphamix = alpha,mq = mf;

    Wij = 0.25*pow(pow(g,2)+pow(gp,2),0.5)*(mixNeut(ineutralino,4)*mixNeut(jneutralino,4) - mixNeut(ineutralino,3)*mixNeut(jneutralino,3));
    double sinthetaW = 0;
    sinthetaW = gp/(pow(pow(g,2)+pow(gp,2),0.5));

    double intZ1 = 0, intZ2 = 0, intZ3 = 0, intZ4 = 0, sminz = 0, smaxz = 0;
    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, MZboson = mZboson;
    sminz = 4*pow(mf,2);
    smaxz = pow(fabs(mneutralinoi)-fabs(mneutralinoj),2);
    intZ1 = dgauss(gneutineutjffZ1dgauss,sminz,smaxz,accuracy);
    intZ2 = dgauss(gneutineutjffZ2dgauss,sminz,smaxz,accuracy);
    intZ3 = dgauss(gneutineutjffZ3dgauss,sminz,smaxz,accuracy);
    intZ4 = dgauss(gneutineutjffZ4dgauss,sminz,smaxz,accuracy);

    GammaZ = 64*pow(g*sinthetaW,2)*pow(Wij,2)*(-4*fabs(mneutralinoi)*fabs(mneutralinoj)*pow(mf,2)*(pow(alphaf,2)-pow(betaf,2))*intZ4*-rj*ri + pow(mf,2)*(pow(alphaf,2)-pow(betaf,2))*intZ3 - fabs(mneutralinoi)*fabs(mneutralinoj)*(pow(alphaf,2) + pow(betaf,2))*intZ2*-rj*ri + 0.5*(pow(alphaf,2) + pow(betaf,2))*intZ1);

    ///First sf t sf u msf1-msf1, msf1-msf2, msf2-msf2 components:
     m1 = mneutralinoi, m2 = msf1, m3 = msf1, m4 = mneutralinoj, mq = mf;
     double xsiZisf1sf1Zj = 0, rhotildaZisf1sf1Zj = 0, chiprimeZisf1sf1Zj = 0, phitildaZisf1sf1Zj = 0, psitildaZisf1sf1Zj = 0, zetaZisf1sf1Zj = 0, XZisf1sf1Zj = 0, chitildaZisf1sf1Zj = 0;
    phitildaZisf1sf1Zj = dgauss(gphitildadgauss,from,to,accuracy);
    psitildaZisf1sf1Zj = dgauss(gpsitildadgauss,from,to,accuracy);
    YZisf1sf1Zj = dgauss(gYdgauss,from,to,accuracy);
    xsiZisf1sf1Zj = dgauss(gxsidgauss,from,to,accuracy);
    rhotildaZisf1sf1Zj = dgauss(grhodgauss,from,to,accuracy);
    chiprimeZisf1sf1Zj = dgauss(gchiprimedgauss,from,to,accuracy);
    /// the problem is in here - DEBUG
    zetaZisf1sf1Zj = dgauss(gzetadgauss,from,to,accuracy);
    XZisf1sf1Zj = dgauss(gXdgauss,from,to,accuracy);
    chitildaZisf1sf1Zj = dgauss(gchidgauss,from,to,accuracy);

    m1 = mneutralinoi, m2 = msf2, m3 = msf2, m4 = mneutralinoj, mq = mf;
    double xsiZisf2sf2Zj = 0, rhotildaZisf2sf2Zj = 0, chiprimeZisf2sf2Zj = 0, phitildaZisf2sf2Zj = 0, psitildaZisf2sf2Zj = 0, zetaZisf2sf2Zj = 0, XZisf2sf2Zj = 0, chitildaZisf2sf2Zj = 0;
    phitildaZisf2sf2Zj = dgauss(gphitildadgauss,from,to,accuracy);
    psitildaZisf2sf2Zj = dgauss(gpsitildadgauss,from,to,accuracy);
    YZisf2sf2Zj = dgauss(gYdgauss,from,to,accuracy);
    xsiZisf2sf2Zj = dgauss(gxsidgauss,from,to,accuracy);
    rhotildaZisf2sf2Zj = dgauss(grhodgauss,from,to,accuracy);
    chiprimeZisf2sf2Zj = dgauss(gchiprimedgauss,from,to,accuracy);
    zetaZisf2sf2Zj = dgauss(gzetadgauss,from,to,accuracy);
    XZisf2sf2Zj = dgauss(gXdgauss,from,to,accuracy);
    chitildaZisf2sf2Zj = dgauss(gchidgauss,from,to,accuracy);

    m1 = mneutralinoi, m2 = msf1, m3 = msf2, m4 = mneutralinoj, mq = mf;

    double xsiZisf1sf2Zj = 0, rhotildaZisf1sf2Zj = 0, chiprimeZisf1sf2Zj = 0, zetaZisf1sf2Zj = 0, XZisf1sf2Zj = 0, chitildaZisf1sf2Zj = 0;
    xsiZisf1sf2Zj = dgauss(gxsidgauss,from,to,accuracy);
    rhotildaZisf1sf2Zj = dgauss(grhodgauss,from,to,accuracy);
    chiprimeZisf1sf2Zj = dgauss(gchiprimedgauss,from,to,accuracy);
    zetaZisf1sf2Zj = dgauss(gzetadgauss,from,to,accuracy);
    XZisf1sf2Zj = dgauss(gXdgauss,from,to,accuracy);
    chitildaZisf1sf2Zj = dgauss(gchidgauss,from,to,accuracy);

    double Gammasftsfumsf1msf2 = 0, Gammasftsfumsf1msf1 = 0, Gammasftsfumsf2msf2 = 0;
    Gammasftsfumsf1msf1 = -2*(8*(sf1alpha1Zi*sf1beta1Zi*sf1beta1Zj*sf1alpha1Zj + sf1beta1Zi*sf1alpha1Zi*sf1alpha1Zj*sf1beta1Zj)*pow(mneutralinoi/PI,2)*YZisf1sf1Zj*ri -(sf1alpha1Zi*sf1alpha1Zi*sf1alpha1Zj*sf1alpha1Zj + sf1beta1Zi*sf1beta1Zi*sf1beta1Zj*sf1beta1Zj)*8*pow(mneutralinoi/PI,2)*phitildaZisf1sf1Zj*ri*rj + (sf1beta1Zi*sf1alpha1Zi*sf1beta1Zj*sf1alpha1Zj + sf1alpha1Zi*sf1beta1Zi*sf1alpha1Zj*sf1beta1Zj)*pow(mf,2)*8*pow(mneutralinoi/PI,2)*xsiZisf1sf1Zj*ri -((sf1alpha1Zi*sf1alpha1Zi*sf1beta1Zj*sf1alpha1Zj + sf1beta1Zi*sf1beta1Zi*sf1alpha1Zj*sf1beta1Zj)*fabs(mneutralinoi)*mf)*(8*pow(mneutralinoi/PI,2)*xsiZisf1sf1Zj - 4*pow(mneutralinoi/PI,2)*(pow(mneutralinoi,2)+pow(mneutralinoj,2))*rhotildaZisf1sf1Zj*ri + 8*pow(mneutralinoi/PI,2)*chiprimeZisf1sf1Zj) + fabs(mneutralinoj)*rj*mf*(sf1beta1Zi*sf1alpha1Zi*sf1alpha1Zj*sf1alpha1Zj + sf1alpha1Zi*sf1beta1Zi*sf1beta1Zj*sf1beta1Zj)*(-8*pow(mneutralinoi/PI,2)*xsiZisf1sf1Zj + 8*pow(mneutralinoi,4)/sqr(PI)*rhotildaZisf1sf1Zj*ri - 8*pow(mneutralinoi/PI,2)*chiprimeZisf1sf1Zj) - (sf1beta1Zi*sf1beta1Zi*sf1beta1Zj*sf1alpha1Zj + sf1alpha1Zi*sf1alpha1Zi*sf1alpha1Zj*sf1beta1Zj)*fabs(mneutralinoi)*ri*mf*(4*pow(mneutralinoi/PI,2)*(pow(mneutralinoi,2)-pow(mneutralinoj,2))*rhotildaZisf1sf1Zj - 8*pow(mneutralinoi/PI,2)*chiprimeZisf1sf1Zj) + (sf1alpha1Zi*sf1beta1Zi*sf1alpha1Zj*sf1alpha1Zj + sf1beta1Zi*sf1alpha1Zi*sf1beta1Zj*sf1beta1Zj)*mf*fabs(mneutralinoj)*rj*ri*8*pow(mneutralinoi/PI,2)*chiprimeZisf1sf1Zj -2*(sf1beta1Zi*sf1beta1Zi*sf1alpha1Zj*sf1alpha1Zj  + sf1alpha1Zi*sf1alpha1Zi*sf1beta1Zj*sf1beta1Zj)*pow(mf,2)*ri*rj*fabs(mneutralinoi)*fabs(mneutralinoj)*4*pow(mneutralinoi/PI,2)*rhotildaZisf1sf1Zj);

    Gammasftsfumsf2msf2 = -2*(8*(sf2alpha1Zi*sf2beta1Zi*sf2beta1Zj*sf2alpha1Zj + sf2beta1Zi*sf2alpha1Zi*sf2alpha1Zj*sf2beta1Zj)*pow(mneutralinoi/PI,2)*YZisf2sf2Zj*ri -(sf2alpha1Zi*sf2alpha1Zi*sf2alpha1Zj*sf2alpha1Zj + sf2beta1Zi*sf2beta1Zi*sf2beta1Zj*sf2beta1Zj)*8*pow(mneutralinoi/PI,2)*phitildaZisf2sf2Zj*ri*rj + (sf2beta1Zi*sf2alpha1Zi*sf2beta1Zj*sf2alpha1Zj + sf2alpha1Zi*sf2beta1Zi*sf2alpha1Zj*sf2beta1Zj)*pow(mf,2)*8*pow(mneutralinoi/PI,2)*xsiZisf2sf2Zj*ri -((sf2alpha1Zi*sf2alpha1Zi*sf2beta1Zj*sf2alpha1Zj + sf2beta1Zi*sf2beta1Zi*sf2alpha1Zj*sf2beta1Zj)*fabs(mneutralinoi)*mf)*(8*pow(mneutralinoi/PI,2)*xsiZisf2sf2Zj - 4*pow(mneutralinoi/PI,2)*(pow(mneutralinoi,2)+pow(mneutralinoj,2))*rhotildaZisf2sf2Zj*ri + 8*pow(mneutralinoi/PI,2)*chiprimeZisf2sf2Zj) + fabs(mneutralinoj)*rj*mf*(sf2beta1Zi*sf2alpha1Zi*sf2alpha1Zj*sf2alpha1Zj + sf2alpha1Zi*sf2beta1Zi*sf2beta1Zj*sf2beta1Zj)*(-8*pow(mneutralinoi/PI,2)*xsiZisf2sf2Zj + 8*pow(mneutralinoi,4)/sqr(PI)*rhotildaZisf2sf2Zj*ri - 8*pow(mneutralinoi/PI,2)*chiprimeZisf2sf2Zj) - (sf2beta1Zi*sf2beta1Zi*sf2beta1Zj*sf2alpha1Zj + sf2alpha1Zi*sf2alpha1Zi*sf2alpha1Zj*sf2beta1Zj)*fabs(mneutralinoi)*ri*mf*(4*pow(mneutralinoi/PI,2)*(pow(mneutralinoi,2)-pow(mneutralinoj,2))*rhotildaZisf2sf2Zj - 8*pow(mneutralinoi/PI,2)*chiprimeZisf2sf2Zj) + (sf2alpha1Zi*sf2beta1Zi*sf2alpha1Zj*sf2alpha1Zj + sf2beta1Zi*sf2alpha1Zi*sf2beta1Zj*sf2beta1Zj)*mf*fabs(mneutralinoj)*rj*ri*8*pow(mneutralinoi/PI,2)*chiprimeZisf2sf2Zj -2*(sf2beta1Zi*sf2beta1Zi*sf2alpha1Zj*sf2alpha1Zj  + sf2alpha1Zi*sf2alpha1Zi*sf2beta1Zj*sf2beta1Zj)*pow(mf,2)*fabs(mneutralinoi)*fabs(mneutralinoj)*ri*rj*4*pow(mneutralinoi/PI,2)*rhotildaZisf2sf2Zj);

    Gammasftsfumsf1msf2 = -2*(8*(sf1alpha1Zi*sf2beta1Zi*sf1beta1Zj*sf2alpha1Zj + sf1beta1Zi*sf2alpha1Zi*sf1alpha1Zj*sf2beta1Zj)*pow(mneutralinoi/PI,2)*YZisf1sf2Zj -(sf1alpha1Zi*sf2alpha1Zi*sf1alpha1Zj*sf2alpha1Zj + sf1beta1Zi*sf2beta1Zi*sf1beta1Zj*sf2beta1Zj)*8*pow(mneutralinoi/PI,2)*phitildaZisf1sf2Zj*ri*rj + (sf1beta1Zi*sf2alpha1Zi*sf1beta1Zj*sf2alpha1Zj + sf1alpha1Zi*sf2beta1Zi*sf1alpha1Zj*sf2beta1Zj)*pow(mf,2)*8*pow(mneutralinoi/PI,2)*xsiZisf1sf2Zj*ri -((sf1alpha1Zi*sf2alpha1Zi*sf1beta1Zj*sf2alpha1Zj + sf1beta1Zi*sf2beta1Zi*sf1alpha1Zj*sf2beta1Zj)*fabs(mneutralinoi)*mf)*(8*pow(mneutralinoi/PI,2)*xsiZisf1sf2Zj - 4*pow(mneutralinoi/PI,2)*(pow(mneutralinoi,2)+pow(mneutralinoj,2))*rhotildaZisf1sf2Zj*ri + 8*pow(mneutralinoi/PI,2)*chiprimeZisf1sf2Zj) + fabs(mneutralinoj)*rj*mf*(sf1beta1Zi*sf2alpha1Zi*sf1alpha1Zj*sf2alpha1Zj + sf1alpha1Zi*sf2beta1Zi*sf1beta1Zj*sf2beta1Zj)*(-8*pow(mneutralinoi/PI,2)*xsiZisf1sf2Zj + 8*pow(mneutralinoi,4)/sqr(PI)*rhotildaZisf1sf2Zj*ri - 8*pow(mneutralinoi/PI,2)*chiprimeZisf1sf2Zj) - (sf1beta1Zi*sf2beta1Zi*sf1beta1Zj*sf2alpha1Zj + sf1alpha1Zi*sf2alpha1Zi*sf1alpha1Zj*sf2beta1Zj)*fabs(mneutralinoi)*ri*mf*(4*pow(mneutralinoi/PI,2)*(pow(mneutralinoi,2)-pow(mneutralinoj,2))*rhotildaZisf1sf2Zj - 8*pow(mneutralinoi/PI,2)*chiprimeZisf1sf2Zj) + (sf1alpha1Zi*sf2beta1Zi*sf1alpha1Zj*sf2alpha1Zj + sf1beta1Zi*sf2alpha1Zi*sf1beta1Zj*sf2beta1Zj)*mf*fabs(mneutralinoj)*rj*8*pow(mneutralinoi/PI,2)*chiprimeZisf1sf2Zj -2*(sf1beta1Zi*sf2beta1Zi*sf1alpha1Zj*sf2alpha1Zj  + sf1alpha1Zi*sf2alpha1Zi*sf1beta1Zj*sf2beta1Zj)*pow(mf,2)*fabs(mneutralinoi)*fabs(mneutralinoj)*ri*rj*4*pow(mneutralinoi/PI,2)*rhotildaZisf1sf2Zj);

    ///sf sf diag contribution 2x as t and u same
    double Gammasf1sf1diag = 0, Gammasf2sf2diag = 0;
    Gammasf1sf1diag = (pow(sf1alpha1Zi,2)+pow(sf1beta1Zi,2))*(pow(sf1alpha1Zj,2)+pow(sf1beta1Zj,2))*8*pow(mneutralinoi/PI,2)*psitildaZisf1sf1Zj + rj*4*(pow(sf1alpha1Zi,2)+pow(sf1beta1Zi,2))*sf1alpha1Zj*sf1beta1Zj*mf*fabs(mneutralinoj)*8*pow(mneutralinoi/PI,2)*chitildaZisf1sf1Zj - ri*4*(pow(sf1alpha1Zj,2)+pow(sf1beta1Zj,2))*sf1alpha1Zi*sf1beta1Zi*fabs(mneutralinoi)*mf*8*pow(mneutralinoi/PI,2)*XZisf1sf1Zj - ri*rj*16*sf1alpha1Zi*sf1beta1Zi*sf1alpha1Zj*sf1beta1Zj*fabs(mneutralinoi)*pow(mf,2)*fabs(mneutralinoj)*4*pow(mneutralinoi/PI,2)*zetaZisf1sf1Zj;
    Gammasf2sf2diag = (pow(sf2alpha1Zi,2)+pow(sf2beta1Zi,2))*(pow(sf2alpha1Zj,2)+pow(sf2beta1Zj,2))*8*pow(mneutralinoi/PI,2)*psitildaZisf2sf2Zj + rj*4*(pow(sf2alpha1Zi,2)+pow(sf2beta1Zi,2))*sf2alpha1Zj*sf2beta1Zj*mf*fabs(mneutralinoj)*8*pow(mneutralinoi/PI,2)*chitildaZisf2sf2Zj - 4*ri*(pow(sf2alpha1Zj,2)+pow(sf2beta1Zj,2))*sf2alpha1Zi*sf2beta1Zi*fabs(mneutralinoi)*mf*8*pow(mneutralinoi/PI,2)*XZisf2sf2Zj - ri*rj*16*sf2alpha1Zi*sf2beta1Zi*sf2alpha1Zj*sf2beta1Zj*fabs(mneutralinoi)*pow(mf,2)*fabs(mneutralinoj)*4*pow(mneutralinoi/PI,2)*zetaZisf2sf2Zj;

    ///sf sf non diag 2x as t and u same
    double Gammasfsfnondiag = 0;
    Gammasfsfnondiag = 16*pow(mneutralinoi/PI,2)*((sf1beta1Zi*sf2beta1Zi + sf1alpha1Zi*sf2alpha1Zi)*(sf1alpha1Zj*sf2alpha1Zj + sf2beta1Zj*sf1beta1Zj)*psitildaZisf1sf2Zj + rj*2*(sf1beta1Zi*sf2beta1Zi + sf1alpha1Zi*sf2alpha1Zi)*(sf2alpha1Zj*sf1beta1Zj + sf1alpha1Zj*sf2beta1Zj)*mf*fabs(mneutralinoj)*chitildaZisf1sf2Zj - ri*2*(sf1alpha1Zi*sf2beta1Zi + sf2alpha1Zi*sf1beta1Zi)*(sf1alpha1Zj*sf2alpha1Zj + sf2beta1Zj*sf1beta1Zj)*fabs(mneutralinoi)*mf*XZisf1sf2Zj - ri*rj*2*(sf1alpha1Zi*sf2beta1Zi + sf2alpha1Zi*sf1beta1Zi)*(sf2alpha1Zj*sf1beta1Zj + sf1alpha1Zj*sf2beta1Zj)*pow(mf,2)*fabs(mneutralinoi)*fabs(mneutralinoj)*zetaZisf1sf2Zj);
   
    double Gammasftot = 0;
    Gammasftot = 2*Gammasf1sf1diag + 2*Gammasf2sf2diag + 2*Gammasfsfnondiag +Gammasftsfumsf1msf1 + Gammasftsfumsf2msf2 + 2*Gammasftsfumsf1msf2;

    double Gammah = 0, GammaH = 0, GammahHinterf = 0, integralh1 = 0, integralh2 = 0, integralh3 = 0, integralh4 =0, integralH1 = 0, integralH2 = 0, integralH3 = 0, integralH4 = 0, integralhH1 = 0, integralhH2 = 0, integralhH3 = 0, integralhH4 = 0;
    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, MZboson = mZboson, mh = mhiggsl, mH = mhiggsH, mA = mhiggsA, g1 = g, g2 = gp, NeutMIX = mixNeut, neutralinoi = ineutralino, neutralinoj = jneutralino, alphamix = alpha,mq = mf;

    integralh1 = dgauss(gintegralh1dgauss,fromz,toz,accuracy);
    integralh2 = dgauss(gintegralh2dgauss,fromz,toz,accuracy);
    integralh3 = dgauss(gintegralh3dgauss,fromz,toz,accuracy);
    integralh4 = dgauss(gintegralh4dgauss,fromz,toz,accuracy);

    integralH1 = dgauss(gintegralH1dgauss,fromz,toz,accuracy);
    integralH2 = dgauss(gintegralH2dgauss,fromz,toz,accuracy);
    integralH3 = dgauss(gintegralH3dgauss,fromz,toz,accuracy);
    integralH4 = dgauss(gintegralH4dgauss,fromz,toz,accuracy);

    integralhH1 = dgauss(gintegralhH1dgauss,fromz,toz,accuracy);
    integralhH2 = dgauss(gintegralhH2dgauss,fromz,toz,accuracy);
    integralhH3 = dgauss(gintegralhH3dgauss,fromz,toz,accuracy);
    integralhH4 = dgauss(gintegralhH4dgauss,fromz,toz,accuracy);

    Gammah = (2*pow(Xijh + Xjih,2)*pow(fq*trigofalphah,2)*(integralh4 -2*sqr(mq)*integralh3 +rj*2*mneutralinoi*fabs(mneutralinoj)*integralh2 - rj*4*mneutralinoi*fabs(mneutralinoj)*sqr(mq)*integralh1));
    
    GammaH = (2*pow(XijH+ XjiH,2)*pow(fq*trigofalphaH,2)*(integralH4 - 2*sqr(mq)*integralH3 +rj*2*mneutralinoi*fabs(mneutralinoj)*integralH2 - rj*4*mneutralinoi*fabs(mneutralinoj)*sqr(mq)*integralH1));

    GammahHinterf = 2*2*(Xijh + Xjih)*(XijH + XjiH)*pow(fq,2)*trigofalphah*trigofalphaH*(integralhH4 - 2*sqr(mq)*integralhH3 + 2*ri*rj*fabs(mneutralinoi)*fabs(mneutralinoj)*integralhH2 - 4*ri*rj*fabs(mneutralinoi)*fabs(mneutralinoj)*sqr(mq)*integralhH1);

    ///Pseudoscalar A contribution:

    double GammaA = 0;
    double integralA1 = 0, integralA2 = 0, integralA3 = 0, integralA4 = 0;
    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, MZboson = mZboson, mh = mhiggsl, mH = mhiggsH, mA = mhiggsA, g1 = g, g2 = gp, NeutMIX = mixNeut, neutralinoi = ineutralino, neutralinoj = jneutralino, alphamix = alpha,mq = mf;

    integralA1 = dgauss(gintegralA1dgauss,fromz,toz,accuracy);
    integralA2 = dgauss(gintegralA2dgauss,fromz,toz,accuracy);
    integralA3 = dgauss(gintegralA3dgauss,fromz,toz,accuracy);
    integralA4 = dgauss(gintegralA4dgauss,fromz,toz,accuracy);

    GammaA = pow(XijA + XjiA,2)*pow(Aq,2)*(integralA4 + 2*sqr(mq)*integralA3 + 2*fabs(mneutralinoi)*fabs(mneutralinoj)*integralA2*-rj*ri + 4*sqr(mq)*fabs(mneutralinoi)*fabs(mneutralinoj)*integralA1*-rj*ri);

    double Gammagoldstone = 0; ///Have an additional pseudoscalar goldstone with mass of Z due to being in Feynman gauge, to add in the longitudinal polarisation of the Z
    double integralgoldstone1 = 0, integralgoldstone2 = 0, integralgoldstone3 = 0, integralgoldstone4 = 0;
    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, MZboson = mZboson, mh = mhiggsl, mH = mhiggsH, mA = mhiggsA, g1 = g, g2 = gp, NeutMIX = mixNeut, neutralinoi = ineutralino, neutralinoj = jneutralino, alphamix = alpha,mq = mf;
    mA = mZboson; ///set so can use same integrals as for pseudoscalar A boson but now with Z mass

    double goldstoneneutineutjcoup = 0;
    goldstoneneutineutjcoup = -0.5*((gp*mixNeut(neutralinoi,1)-g*mixNeut(neutralinoi,2))*(mixNeut(neutralinoj,3)*-cos(beta) - mixNeut(neutralinoj,4)*sin(beta)) + (gp*mixNeut(neutralinoj,1)-g*mixNeut(neutralinoj,2))*(mixNeut(neutralinoi,3)*-cos(beta) - mixNeut(neutralinoi,4)*sin(beta)));

    integralgoldstone1 = dgauss(gintegralA1dgauss,fromz,toz,accuracy);
    integralgoldstone2 = dgauss(gintegralA2dgauss,fromz,toz,accuracy);
    integralgoldstone3 = dgauss(gintegralA3dgauss,fromz,toz,accuracy);
    integralgoldstone4 = dgauss(gintegralA4dgauss,fromz,toz,accuracy);

    Gammagoldstone = pow(goldstoneneutineutjcoup,2)*pow(goldstoneffcoup,2)*4*(integralgoldstone4 + 2*sqr(mq)*integralgoldstone3 + 2*fabs(mneutralinoi)*fabs(mneutralinoj)*integralgoldstone2*-rj*ri + 4*sqr(mq)*fabs(mneutralinoi)*fabs(mneutralinoj)*integralgoldstone1*-rj*ri);

    double coupcombo1Zsf1 = 0, coupcombo2Zsf1 = 0, coupcombo3Zsf1 = 0, coupcombo4Zsf1 = 0, coupcombo5Zsf1 = 0, coupcombo6Zsf1 = 0, coupcombo7Zsf1 = 0, coupcombo8Zsf1 = 0; 
    coupcombo1Zsf1 = -2*2*Wij*g*sinthetaW*(-sf1alpha1Zi*(alphaf-betaf)*sf1beta1Zj + sf1beta1Zi*(alphaf+betaf)*sf1alpha1Zj)*mf*fabs(mneutralinoi);
    coupcombo2Zsf1 = 2*2*Wij*g*sinthetaW*(-sf1alpha1Zi*(alphaf+betaf)*sf1beta1Zj + sf1beta1Zi*(alphaf-betaf)*sf1alpha1Zj)*mf*fabs(mneutralinoj)*-ri*rj;
    coupcombo3Zsf1 = -2*2*Wij*g*sinthetaW*(sf1beta1Zi*(alphaf+betaf)*sf1beta1Zj - sf1alpha1Zi*(alphaf - betaf)*sf1alpha1Zj)*ri;
    coupcombo4Zsf1 = -4*2*Wij*g*sinthetaW*(-sf1alpha1Zi*(alphaf+betaf)*sf1beta1Zj + sf1beta1Zi*(alphaf-betaf)*sf1alpha1Zj)*fabs(mneutralinoi)*mf;
    coupcombo5Zsf1 = 4*2*Wij*g*sinthetaW*(-sf1alpha1Zi*(alphaf-betaf)*sf1beta1Zj + sf1beta1Zi*(alphaf+betaf)*sf1alpha1Zj)*mf*fabs(mneutralinoj)*-ri*rj;
    coupcombo6Zsf1 = 2*2*Wij*g*sinthetaW*(sf1beta1Zi*(alphaf+betaf)*sf1beta1Zj - sf1alpha1Zi*(alphaf-betaf)*sf1alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj)*-rj;
    coupcombo7Zsf1 = -2*2*Wij*g*sinthetaW*(sf1beta1Zi*(alphaf-betaf)*sf1beta1Zj - sf1alpha1Zi*(alphaf+betaf)*sf1alpha1Zj)*pow(mf,2)*ri;
    coupcombo8Zsf1 = 8*2*Wij*g*sinthetaW*(sf1beta1Zi*(alphaf-betaf)*sf1beta1Zj - sf1alpha1Zi*(alphaf+betaf)*sf1alpha1Zj)*pow(mf,2)*fabs(mneutralinoi)*fabs(mneutralinoj)*-rj;

    double integral1Zsf1 = 0, integral2Zsf1 = 0, integral3Zsf1 = 0, integral4Zsf1 = 0, integral5Zsf1 = 0, integral6Zsf1 = 0, integral7Zsf1 = 0, integral8Zsf1 = 0, integral1Zsf2 = 0, integral2Zsf2 = 0, integral3Zsf2 = 0, integral4Zsf2 = 0, integral5Zsf2 = 0, integral6Zsf2 = 0, integral7Zsf2 = 0, integral8Zsf2 = 0;

    double fromzsfE = 0, tozsfE = 0;
    fromzsfE = fabs(mneutralinoj);
    tozsfE = (pow(mneutralinoi,2) + pow(mneutralinoj,2) - 4*pow(mf,2))/(2*fabs(mneutralinoi));
    
    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, MZboson = mZboson, m2 = msf1;
    integral1Zsf1 = 2*fabs(m1)*dgauss(gintegral1Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2Zsf1 = 2*fabs(m1)*dgauss(gintegral2Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3Zsf1 = 2*fabs(m1)*dgauss(gintegral3Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4Zsf1 = 2*fabs(m1)*dgauss(gintegral4Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5Zsf1 = 2*fabs(m1)*dgauss(gintegral5Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6Zsf1 = 2*fabs(m1)*dgauss(gintegral6Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7Zsf1 = 2*fabs(m1)*dgauss(gintegral7Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8Zsf1 = 2*fabs(m1)*dgauss(gintegral8Zsfdgauss,fromzsfE,tozsfE,accuracy);

    GammaZsf1 = (coupcombo1Zsf1*integral1Zsf1 + coupcombo2Zsf1*integral2Zsf1 + coupcombo3Zsf1*integral3Zsf1 + coupcombo4Zsf1*integral4Zsf1 + coupcombo5Zsf1*integral5Zsf1 + coupcombo6Zsf1*integral6Zsf1 + coupcombo7Zsf1*integral7Zsf1 + coupcombo8Zsf1*integral8Zsf1)*ri;

    double coupcombo1Zsf2 = 0, coupcombo2Zsf2 = 0, coupcombo3Zsf2 = 0, coupcombo4Zsf2 = 0, coupcombo5Zsf2 = 0, coupcombo6Zsf2 = 0, coupcombo7Zsf2 = 0, coupcombo8Zsf2 = 0;
    coupcombo1Zsf2 = 2*2*Wij*g*sinthetaW*(-sf2alpha1Zi*(alphaf-betaf)*sf2beta1Zj + sf2beta1Zi*(alphaf+betaf)*sf2alpha1Zj)*mf*fabs(mneutralinoi);
    coupcombo2Zsf2 = -2*2*Wij*g*sinthetaW*(-sf2alpha1Zi*(alphaf+betaf)*sf2beta1Zj + sf2beta1Zi*(alphaf-betaf)*sf2alpha1Zj)*mf*fabs(mneutralinoj)*-ri*rj;
    coupcombo3Zsf2 = 2*2*Wij*g*sinthetaW*(sf2beta1Zi*(alphaf+betaf)*sf2beta1Zj - sf2alpha1Zi*(alphaf - betaf)*sf2alpha1Zj)*ri;
    coupcombo4Zsf2 = 4*2*Wij*g*sinthetaW*(-sf2alpha1Zi*(alphaf+betaf)*sf2beta1Zj + sf2beta1Zi*(alphaf-betaf)*sf2alpha1Zj)*fabs(mneutralinoi)*mf;
    coupcombo5Zsf2 = -4*2*Wij*g*sinthetaW*(-sf2alpha1Zi*(alphaf-betaf)*sf2beta1Zj + sf2beta1Zi*(alphaf+betaf)*sf2alpha1Zj)*mf*fabs(mneutralinoj)*-ri*rj;
    coupcombo6Zsf2 =  -2*2*Wij*g*sinthetaW*(sf2beta1Zi*(alphaf+betaf)*sf2beta1Zj - sf2alpha1Zi*(alphaf-betaf)*sf2alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj)*-rj;
    coupcombo7Zsf2 = 2*2*Wij*g*sinthetaW*(sf2beta1Zi*(alphaf-betaf)*sf2beta1Zj - sf2alpha1Zi*(alphaf+betaf)*sf2alpha1Zj)*pow(mf,2)*ri;
    coupcombo8Zsf2 = -8*2*Wij*g*sinthetaW*(sf2beta1Zi*(alphaf-betaf)*sf2beta1Zj - sf2alpha1Zi*(alphaf+betaf)*sf2alpha1Zj)*pow(mf,2)*fabs(mneutralinoi)*fabs(mneutralinoj)*-rj;

    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, MZboson = mZboson, m2 = msf2;
    integral1Zsf2 = 2*fabs(m1)*dgauss(gintegral1Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2Zsf2 = 2*fabs(m1)*dgauss(gintegral2Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3Zsf2 = 2*fabs(m1)*dgauss(gintegral3Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4Zsf2 = 2*fabs(m1)*dgauss(gintegral4Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5Zsf2 = 2*fabs(m1)*dgauss(gintegral5Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6Zsf2 = 2*fabs(m1)*dgauss(gintegral6Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7Zsf2 = 2*fabs(m1)*dgauss(gintegral7Zsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8Zsf2 = 2*fabs(m1)*dgauss(gintegral8Zsfdgauss,fromzsfE,tozsfE,accuracy);

    GammaZsf2 = (coupcombo1Zsf2*integral1Zsf2 + coupcombo2Zsf2*integral2Zsf2 + coupcombo3Zsf2*integral3Zsf2 + coupcombo4Zsf2*integral4Zsf2 + coupcombo5Zsf2*integral5Zsf2 + coupcombo6Zsf2*integral6Zsf2 + coupcombo7Zsf2*integral7Zsf2 + coupcombo8Zsf2*integral8Zsf2)*ri;

    double coupcombo1hsf1 = 0, coupcombo2hsf1 = 0, coupcombo3hsf1 = 0, coupcombo4hsf1 = 0, coupcombo5hsf1 = 0, coupcombo6hsf1 = 0, coupcombo7hsf1 = 0, coupcombo8hsf1 = 0, coupcombo1Hsf1 = 0, coupcombo2Hsf1 = 0, coupcombo3Hsf1 = 0, coupcombo4Hsf1 = 0, coupcombo5Hsf1 = 0, coupcombo6Hsf1 = 0, coupcombo7Hsf1 = 0, coupcombo8Hsf1 = 0, coupcombo1hsf2 = 0, coupcombo2hsf2 = 0, coupcombo3hsf2 = 0, coupcombo4hsf2 = 0, coupcombo5hsf2 = 0, coupcombo6hsf2 = 0, coupcombo7hsf2 = 0, coupcombo8hsf2 = 0, coupcombo1Hsf2 = 0, coupcombo2Hsf2 = 0, coupcombo3Hsf2 = 0, coupcombo4Hsf2 = 0, coupcombo5Hsf2 = 0, coupcombo6Hsf2 = 0, coupcombo7Hsf2 = 0, coupcombo8Hsf2 = 0;

    coupcombo1hsf1 = 0.5*(Xijh+Xjih)*fq*trigofalphah/(root2)*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*-ri*rj;
    coupcombo2hsf1 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*mq*-rj;
    coupcombo3hsf1 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*mq*fabs(mneutralinoj)*ri;
    coupcombo4hsf1 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(-sf1beta1Zi*sf1beta1Zj - sf1alpha1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*mq*-rj;
    coupcombo5hsf1 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*-mq*fabs(mneutralinoj)*ri;
    coupcombo6hsf1 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(-sf1alpha1Zi*sf1beta1Zj - sf1beta1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj);
    coupcombo7hsf1 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*-sqr(mq)*-ri*rj;
    coupcombo8hsf1 = 2*(Xijh+Xjih)*fq*trigofalphah/(root2)*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*sqr(mq)*fabs(mneutralinoj);

    coupcombo1Hsf1 = 0.5*(XijH+XjiH)*fq*trigofalphaH/(root2)*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*-ri*rj;
    coupcombo2Hsf1 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*mq*-rj;
    coupcombo3Hsf1 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*mq*fabs(mneutralinoj)*ri;
    coupcombo4Hsf1 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(-sf1beta1Zi*sf1beta1Zj - sf1alpha1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*mq*-rj;
    coupcombo5Hsf1 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*-mq*fabs(mneutralinoj)*ri;
    coupcombo6Hsf1 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(-sf1alpha1Zi*sf1beta1Zj - sf1beta1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj);
    coupcombo7Hsf1 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*-sqr(mq)*-ri*rj;
    coupcombo8Hsf1 = 2*(XijH+XjiH)*fq*trigofalphaH/(root2)*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*sqr(mq)*fabs(mneutralinoj);

    coupcombo1hsf2 = 0.5*(Xijh+Xjih)*fq*trigofalphah/(root2)*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*-ri*rj;
    coupcombo2hsf2 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(sf2beta1Zi*sf2beta1Zj + sf2alpha1Zi*sf2alpha1Zj)*fabs(mneutralinoi)*mq*-rj;
    coupcombo3hsf2 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(sf2beta1Zi*sf2beta1Zj + sf2alpha1Zi*sf2alpha1Zj)*mq*fabs(mneutralinoj)*ri;
    coupcombo4hsf2 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(-sf2beta1Zi*sf2beta1Zj - sf2alpha1Zi*sf2alpha1Zj)*fabs(mneutralinoi)*mq*-rj;
    coupcombo5hsf2 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(sf2beta1Zi*sf2beta1Zj + sf2alpha1Zi*sf2alpha1Zj)*-mq*fabs(mneutralinoj)*ri;
    coupcombo6hsf2 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(-sf2alpha1Zi*sf2beta1Zj - sf2beta1Zi*sf2alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj);
    coupcombo7hsf2 = (Xijh+Xjih)*fq*trigofalphah/(root2)*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*-sqr(mq)*-ri*rj;
    coupcombo8hsf2 = 2*(Xijh+Xjih)*fq*trigofalphah/(root2)*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*fabs(mneutralinoi)*sqr(mq)*fabs(mneutralinoj);

    coupcombo1Hsf2 = 0.5*(XijH+XjiH)*fq*trigofalphaH/(root2)*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*-ri*rj;
    coupcombo2Hsf2 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(sf2beta1Zi*sf2beta1Zj + sf2alpha1Zi*sf2alpha1Zj)*fabs(mneutralinoi)*mq*-rj;
    coupcombo3Hsf2 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(sf2beta1Zi*sf2beta1Zj + sf2alpha1Zi*sf2alpha1Zj)*mq*fabs(mneutralinoj)*ri;
    coupcombo4Hsf2 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(-sf2beta1Zi*sf2beta1Zj - sf2alpha1Zi*sf2alpha1Zj)*fabs(mneutralinoi)*mq*-rj;
    coupcombo5Hsf2 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(sf2beta1Zi*sf2beta1Zj + sf2alpha1Zi*sf2alpha1Zj)*-mq*fabs(mneutralinoj)*ri;
    coupcombo6Hsf2 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(-sf2alpha1Zi*sf2beta1Zj - sf2beta1Zi*sf2alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj);
    coupcombo7Hsf2 = (XijH+XjiH)*fq*trigofalphaH/(root2)*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*-sqr(mq)*-ri*rj;
    coupcombo8Hsf2 = 2*(XijH+XjiH)*fq*trigofalphaH/(root2)*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*fabs(mneutralinoi)*sqr(mq)*fabs(mneutralinoj);
    
    double integral1hsf1 = 0, integral2hsf1 = 0, integral3hsf1 = 0, integral4hsf1 = 0, integral5hsf1 = 0, integral6hsf1 = 0, integral7hsf1 = 0, integral8hsf1 = 0, integral1hsf2 = 0, integral2hsf2 = 0, integral3hsf2 = 0, integral4hsf2 = 0, integral5hsf2 = 0, integral6hsf2 = 0, integral7hsf2 = 0, integral8hsf2 = 0, integral1Hsf1 = 0, integral2Hsf1 = 0, integral3Hsf1 = 0, integral4Hsf1 = 0, integral5Hsf1 = 0, integral6Hsf1 = 0, integral7Hsf1 = 0, integral8Hsf1 = 0, integral1Hsf2 = 0, integral2Hsf2 = 0, integral3Hsf2 = 0, integral4Hsf2 = 0, integral5Hsf2 = 0, integral6Hsf2 = 0, integral7Hsf2 = 0, integral8Hsf2 = 0;

    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mh = mhiggsl, mH = mhiggsH, mA = mhiggsA, m2 = msf1;
    integral1hsf1 = 2*fabs(m1)*dgauss(gintegral1hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2hsf1 = 2*fabs(m1)*dgauss(gintegral2hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3hsf1 = 2*fabs(m1)*dgauss(gintegral3hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4hsf1 = 2*fabs(m1)*dgauss(gintegral4hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5hsf1 = 2*fabs(m1)*dgauss(gintegral5hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6hsf1 = 2*fabs(m1)*dgauss(gintegral6hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7hsf1 = 2*fabs(m1)*dgauss(gintegral7hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8hsf1 = 2*fabs(m1)*dgauss(gintegral8hsfdgauss,fromzsfE,tozsfE,accuracy);

    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mh = mhiggsl, mH = mhiggsH, mA = mhiggsA, m2 = msf1;
    integral1Hsf1 = 2*fabs(m1)*dgauss(gintegral1Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2Hsf1 = 2*fabs(m1)*dgauss(gintegral2Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3Hsf1 = 2*fabs(m1)*dgauss(gintegral3Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4Hsf1 = 2*fabs(m1)*dgauss(gintegral4Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5Hsf1 = 2*fabs(m1)*dgauss(gintegral5Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6Hsf1 = 2*fabs(m1)*dgauss(gintegral6Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7Hsf1 = 2*fabs(m1)*dgauss(gintegral7Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8Hsf1 = 2*fabs(m1)*dgauss(gintegral8Hsfdgauss,fromzsfE,tozsfE,accuracy);

    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mh = mhiggsl, mH = mhiggsH, mA = mhiggsA, m2 = msf2;
    integral1hsf2 = 2*fabs(m1)*dgauss(gintegral1hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2hsf2 = 2*fabs(m1)*dgauss(gintegral2hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3hsf2 = 2*fabs(m1)*dgauss(gintegral3hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4hsf2 = 2*fabs(m1)*dgauss(gintegral4hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5hsf2 = 2*fabs(m1)*dgauss(gintegral5hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6hsf2 = 2*fabs(m1)*dgauss(gintegral6hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7hsf2 = 2*fabs(m1)*dgauss(gintegral7hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8hsf2 = 2*fabs(m1)*dgauss(gintegral8hsfdgauss,fromzsfE,tozsfE,accuracy);

    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mh = mhiggsl, mH = mhiggsH, mA = mhiggsA, m2 = msf2;
    integral1Hsf2 = 2*fabs(m1)*dgauss(gintegral1Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2Hsf2 = 2*fabs(m1)*dgauss(gintegral2Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3Hsf2 = 2*fabs(m1)*dgauss(gintegral3Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4Hsf2 = 2*fabs(m1)*dgauss(gintegral4Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5Hsf2 = 2*fabs(m1)*dgauss(gintegral5Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6Hsf2 = 2*fabs(m1)*dgauss(gintegral6Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7Hsf2 = 2*fabs(m1)*dgauss(gintegral7Hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8Hsf2 = 2*fabs(m1)*dgauss(gintegral8Hsfdgauss,fromzsfE,tozsfE,accuracy);

    Gammahsf1 = coupcombo1hsf1*integral1hsf1+coupcombo2hsf1*integral2hsf1+coupcombo3hsf1*integral3hsf1+coupcombo4hsf1*integral4hsf1+coupcombo5hsf1*integral5hsf1+coupcombo6hsf1*integral6hsf1+coupcombo7hsf1*integral7hsf1+coupcombo8hsf1*integral8hsf1;
    GammaHsf1 = coupcombo1Hsf1*integral1Hsf1+coupcombo2Hsf1*integral2Hsf1+coupcombo3Hsf1*integral3Hsf1+coupcombo4Hsf1*integral4Hsf1+coupcombo5Hsf1*integral5Hsf1+coupcombo6Hsf1*integral6Hsf1+coupcombo7Hsf1*integral7Hsf1+coupcombo8Hsf1*integral8Hsf1;
    Gammahsf2 = coupcombo1hsf2*integral1hsf2+coupcombo2hsf2*integral2hsf2+coupcombo3hsf2*integral3hsf2+coupcombo4hsf2*integral4hsf2+coupcombo5hsf2*integral5hsf2+coupcombo6hsf2*integral6hsf2+coupcombo7hsf2*integral7hsf2+coupcombo8hsf2*integral8hsf2;
    GammaHsf2 = coupcombo1Hsf2*integral1Hsf2+coupcombo2Hsf2*integral2Hsf2+coupcombo3Hsf2*integral3Hsf2+coupcombo4Hsf2*integral4Hsf2+coupcombo5Hsf2*integral5Hsf2+coupcombo6Hsf2*integral6Hsf2+coupcombo7Hsf2*integral7Hsf2+coupcombo8Hsf2*integral8Hsf2;
    
    double coupcombo1Asf1 = 0, coupcombo2Asf1 = 0, coupcombo3Asf1 = 0, coupcombo4Asf1 = 0, coupcombo5Asf1 = 0, coupcombo6Asf1 = 0, coupcombo7Asf1 = 0, coupcombo8Asf1 = 0, coupcombo1Asf2 = 0, coupcombo2Asf2 = 0, coupcombo3Asf2 = 0, coupcombo4Asf2 = 0, coupcombo5Asf2 = 0, coupcombo6Asf2 = 0, coupcombo7Asf2 = 0, coupcombo8Asf2 = 0;  
    
    coupcombo1Asf1 = 0.5*(XijA + XjiA)*Aq/2*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*ri;
    coupcombo2Asf1 = -(XijA + XjiA)*Aq/2*fabs(mneutralinoi)*mq*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj);
    coupcombo3Asf1 = (XijA + XjiA)*Aq/2*(-sf1beta1Zi*sf1beta1Zj - sf1alpha1Zi*sf1alpha1Zj)*mq*fabs(mneutralinoj)*-ri*rj;
    coupcombo4Asf1 = (XijA + XjiA)*Aq/2*(-sf1beta1Zi*sf1beta1Zj - sf1alpha1Zi*sf1alpha1Zj)*mq*fabs(mneutralinoi);
    coupcombo5Asf1 = -(XijA + XjiA)*Aq/2*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*mq*fabs(mneutralinoj)*-ri*rj;
    coupcombo6Asf1 = (XijA + XjiA)*Aq/2*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj)*-rj;
    coupcombo7Asf1 = (XijA + XjiA)*Aq/2*sqr(mq)*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*ri;
    coupcombo8Asf1 = (XijA + XjiA)*Aq*fabs(mneutralinoi)*sqr(mq)*fabs(mneutralinoj)*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*-rj;

    coupcombo1Asf2 = 0.5*(XijA + XjiA)*Aq/2*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*ri;
    coupcombo2Asf2 = -(XijA + XjiA)*Aq/2*fabs(mneutralinoi)*mq*(sf2beta1Zi*sf2beta1Zj + sf2alpha1Zi*sf2alpha1Zj);
    coupcombo3Asf2 = (XijA + XjiA)*Aq/2*(-sf2beta1Zi*sf2beta1Zj - sf2alpha1Zi*sf2alpha1Zj)*mq*fabs(mneutralinoj)*-ri*rj;
    coupcombo4Asf2 = (XijA + XjiA)*Aq/2*(-sf2beta1Zi*sf2beta1Zj - sf2alpha1Zi*sf2alpha1Zj)*mq*fabs(mneutralinoi);
    coupcombo5Asf2 = -(XijA + XjiA)*Aq/2*(sf2beta1Zi*sf2beta1Zj + sf2alpha1Zi*sf2alpha1Zj)*mq*fabs(mneutralinoj)*-ri*rj;
    coupcombo6Asf2 = (XijA + XjiA)*Aq/2*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj)*-rj;
    coupcombo7Asf2 = (XijA + XjiA)*Aq/2*sqr(mq)*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*ri;
    coupcombo8Asf2 = (XijA + XjiA)*Aq*fabs(mneutralinoi)*sqr(mq)*fabs(mneutralinoj)*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj)*-rj;

    double integral1Asf1 = 0, integral2Asf1 = 0, integral3Asf1 = 0, integral4Asf1 = 0, integral5Asf1 = 0, integral6Asf1 = 0, integral7Asf1 = 0, integral8Asf1 = 0, integral1Asf2 = 0, integral2Asf2 = 0, integral3Asf2 = 0, integral4Asf2 = 0, integral5Asf2 = 0, integral6Asf2 = 0, integral7Asf2 = 0, integral8Asf2 = 0; 

    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mh = mhiggsA, m2 = msf1; ///here set mh = mA as using the integrands written for h with mh -> mA
    integral1Asf1 = 2*fabs(m1)*dgauss(gintegral1hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2Asf1 = 2*fabs(m1)*dgauss(gintegral2hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3Asf1 = 2*fabs(m1)*dgauss(gintegral3hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4Asf1 = 2*fabs(m1)*dgauss(gintegral4hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5Asf1 = 2*fabs(m1)*dgauss(gintegral5hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6Asf1 = 2*fabs(m1)*dgauss(gintegral6hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7Asf1 = 2*fabs(m1)*dgauss(gintegral7hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8Asf1 = 2*fabs(m1)*dgauss(gintegral8hsfdgauss,fromzsfE,tozsfE,accuracy);

    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mh = mhiggsA, m2 = msf2; ///here set mh = mA as using the integrands written for h with mh -> mA
    integral1Asf2 = 2*fabs(m1)*dgauss(gintegral1hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2Asf2 = 2*fabs(m1)*dgauss(gintegral2hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3Asf2 = 2*fabs(m1)*dgauss(gintegral3hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4Asf2 = 2*fabs(m1)*dgauss(gintegral4hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5Asf2 = 2*fabs(m1)*dgauss(gintegral5hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6Asf2 = 2*fabs(m1)*dgauss(gintegral6hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7Asf2 = 2*fabs(m1)*dgauss(gintegral7hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8Asf2 = 2*fabs(m1)*dgauss(gintegral8hsfdgauss,fromzsfE,tozsfE,accuracy);

    GammaAsf1 = coupcombo1Asf1*integral1Asf1 + coupcombo2Asf1*integral2Asf1 + coupcombo3Asf1*integral3Asf1 + coupcombo4Asf1*integral4Asf1+ coupcombo5Asf1*integral5Asf1 + coupcombo6Asf1*integral6Asf1 + coupcombo7Asf1*integral7Asf1 + coupcombo8Asf1*integral8Asf1;

    GammaAsf2 = coupcombo1Asf2*integral1Asf2 + coupcombo2Asf2*integral2Asf2 + coupcombo3Asf2*integral3Asf2 + coupcombo4Asf2*integral4Asf2+ coupcombo5Asf2*integral5Asf2 + coupcombo6Asf2*integral6Asf2 + coupcombo7Asf2*integral7Asf2 + coupcombo8Asf2*integral8Asf2;

    ///Goldstone - sf interference:
    double Gammagsf1 = 0, Gammagsf2 = 0;

    double coupcombo1gsf1 = 0, coupcombo2gsf1 = 0, coupcombo3gsf1 = 0, coupcombo4gsf1 = 0, coupcombo5gsf1 = 0, coupcombo6gsf1 = 0, coupcombo7gsf1 = 0, coupcombo8gsf1 = 0; 
    coupcombo1gsf1 = 0.5*goldstoneneutineutjcoup*goldstoneffcoup*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj);
    coupcombo2gsf1 = -goldstoneneutineutjcoup*goldstoneffcoup*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*mf*ri;
    coupcombo3gsf1 = goldstoneneutineutjcoup*goldstoneffcoup*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*mf*fabs(mneutralinoj)*rj;
    coupcombo4gsf1 = -goldstoneneutineutjcoup*goldstoneffcoup*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*mf*ri;
    coupcombo5gsf1 = goldstoneneutineutjcoup*goldstoneffcoup*(sf1beta1Zi*sf1beta1Zj + sf1alpha1Zi*sf1alpha1Zj)*fabs(mneutralinoj)*mf*rj;
    coupcombo6gsf1 = goldstoneneutineutjcoup*goldstoneffcoup*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj)*-rj*ri;
    coupcombo7gsf1 = goldstoneneutineutjcoup*goldstoneffcoup*(sf1alpha1Zi*sf1beta1Zj + sf1beta1Zi*sf1alpha1Zj)*pow(mf,2);
    coupcombo8gsf1 = 2*goldstoneneutineutjcoup*goldstoneffcoup*(-sf1alpha1Zi*sf1beta1Zj - sf1beta1Zi*sf1alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj)*pow(mf,2)*rj*ri;

    double integral1gsf1 = 0, integral2gsf1 = 0, integral3gsf1 = 0, integral4gsf1 = 0, integral5gsf1 = 0, integral6gsf1 = 0, integral7gsf1 = 0, integral8gsf1 = 0, integral1gsf2 = 0, integral2gsf2 = 0, integral3gsf2 = 0, integral4gsf2 = 0, integral5gsf2 = 0, integral6gsf2 = 0, integral7gsf2 = 0, integral8gsf2 = 0; 

    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mh = mZboson, m2 = msf1; ///here set mh = mZ as using the integrands written for h with mh -> mA -> mZ
    integral1gsf1 = 2*fabs(m1)*dgauss(gintegral1hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2gsf1 = 2*fabs(m1)*dgauss(gintegral2hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3gsf1 = 2*fabs(m1)*dgauss(gintegral3hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4gsf1 = 2*fabs(m1)*dgauss(gintegral4hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5gsf1 = 2*fabs(m1)*dgauss(gintegral5hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6gsf1 = 2*fabs(m1)*dgauss(gintegral6hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7gsf1 = 2*fabs(m1)*dgauss(gintegral7hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8gsf1 = 2*fabs(m1)*dgauss(gintegral8hsfdgauss,fromzsfE,tozsfE,accuracy);

    Gammagsf1 = (coupcombo1gsf1*integral1gsf1 + coupcombo2gsf1*integral2gsf1 + coupcombo3gsf1*integral3gsf1 + coupcombo4gsf1*integral4gsf1 + coupcombo5gsf1*integral5gsf1 + coupcombo6gsf1*integral6gsf1 + coupcombo7gsf1*integral7gsf1 + coupcombo8gsf1*integral8gsf1);

    double coupcombo1gsf2 = 0, coupcombo2gsf2 = 0, coupcombo3gsf2 = 0, coupcombo4gsf2 = 0, coupcombo5gsf2 = 0, coupcombo6gsf2 = 0, coupcombo7gsf2 = 0, coupcombo8gsf2 = 0; 
    coupcombo1gsf2 = 0.5*goldstoneneutineutjcoup*goldstoneffcoup*(sf2alpha1Zi*sf2beta1Zj + sf2beta1Zi*sf2alpha1Zj);
    coupcombo2gsf2 = -goldstoneneutineutjcoup*goldstoneffcoup*(-sf2beta1Zi*-sf2beta1Zj + -sf2alpha1Zi*-sf2alpha1Zj)*fabs(mneutralinoi)*-mf*-ri;
    coupcombo3gsf2 = -goldstoneneutineutjcoup*goldstoneffcoup*(-sf2beta1Zi*sf2beta1Zj + sf2alpha1Zi*-sf2alpha1Zj)*mf*fabs(mneutralinoj)*rj;
    coupcombo4gsf2 = -goldstoneneutineutjcoup*goldstoneffcoup*(-sf2beta1Zi*-sf2beta1Zj - sf2alpha1Zi*-sf2alpha1Zj)*fabs(mneutralinoi)*mf*ri;
    coupcombo5gsf2 = goldstoneneutineutjcoup*goldstoneffcoup*(-sf2beta1Zi*-sf2beta1Zj - sf2alpha1Zi*-sf2alpha1Zj)*fabs(mneutralinoj)*-mf*-rj;
    coupcombo6gsf2 = goldstoneneutineutjcoup*goldstoneffcoup*(-sf2alpha1Zi*-sf2beta1Zj - sf2beta1Zi*-sf2alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj)*-rj*ri;
    coupcombo7gsf2 = goldstoneneutineutjcoup*goldstoneffcoup*(-sf2alpha1Zi*-sf2beta1Zj - sf2beta1Zi*-sf2alpha1Zj)*pow(mf,2);
    coupcombo8gsf2 = 2*goldstoneneutineutjcoup*goldstoneffcoup*(sf2alpha1Zi*-sf2beta1Zj + sf2beta1Zi*-sf2alpha1Zj)*fabs(mneutralinoi)*fabs(mneutralinoj)*-pow(mf,2)*-rj*ri;

    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mh = mZboson, m2 = msf2; ///here set mh = mZ as using the integrands written for h with mh -> mA -> mZ
    integral1gsf2 = 2*fabs(m1)*dgauss(gintegral1hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral2gsf2 = 2*fabs(m1)*dgauss(gintegral2hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral3gsf2 = 2*fabs(m1)*dgauss(gintegral3hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral4gsf2 = 2*fabs(m1)*dgauss(gintegral4hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral5gsf2 = 2*fabs(m1)*dgauss(gintegral5hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral6gsf2 = 2*fabs(m1)*dgauss(gintegral6hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral7gsf2 = 2*fabs(m1)*dgauss(gintegral7hsfdgauss,fromzsfE,tozsfE,accuracy);
    integral8gsf2 = 2*fabs(m1)*dgauss(gintegral8hsfdgauss,fromzsfE,tozsfE,accuracy);

    Gammagsf2 = (coupcombo1gsf2*integral1gsf2 + coupcombo2gsf2*integral2gsf2 + coupcombo3gsf2*integral3gsf2 + coupcombo4gsf2*integral4gsf2 + coupcombo5gsf2*integral5gsf2 + coupcombo6gsf2*integral6gsf2 + coupcombo7gsf2*integral7gsf2 + coupcombo8gsf2*integral8gsf2);

    double coupcombo1ZA = 0, coupcombo2ZA = 0, coupcombo3ZA = 0, coupcombo4ZA = 0;
    coupcombo1ZA = 4*Wij*(XijA+XjiA)*Aq*g*sinthetaW*betaf*fabs(mneutralinoj)*mq*-ri*rj;
    coupcombo2ZA = 4*Wij*(XijA+XjiA)*Aq*g*sinthetaW*betaf*fabs(mneutralinoi)*mq;
    coupcombo3ZA = 4*Wij*(XijA+XjiA)*Aq*g*sinthetaW*betaf*mq*fabs(mneutralinoj)*-ri*rj;
    coupcombo4ZA = 4*Wij*(XijA+XjiA)*Aq*g*sinthetaW*betaf*fabs(mneutralinoi)*mq;
    
    double integral1ZA = 0, integral2ZA = 0, integral3ZA = 0, integral4ZA = 0;
    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mA = mhiggsA, m2 = msf1, MZboson = mZboson;
    integral1ZA = 2*fabs(m1)*dgauss(gintegral1ZAdgauss,fromzsfE,tozsfE,accuracy);
    integral2ZA = 2*fabs(m1)*dgauss(gintegral2ZAdgauss,fromzsfE,tozsfE,accuracy);
    integral3ZA = 2*fabs(m1)*dgauss(gintegral3ZAdgauss,fromzsfE,tozsfE,accuracy);
    integral4ZA = 2*fabs(m1)*dgauss(gintegral4ZAdgauss,fromzsfE,tozsfE,accuracy);

    double GammaZA = 0;
    
    GammaZA = coupcombo1ZA*integral1ZA + coupcombo2ZA*integral2ZA + coupcombo3ZA*integral3ZA + coupcombo4ZA*integral4ZA;

    double coupcombo1Zg = 0, coupcombo2Zg = 0, coupcombo3Zg = 0, coupcombo4Zg = 0;
    coupcombo1Zg = 8*Wij*goldstoneffcoup*goldstoneneutineutjcoup*g*sinthetaW*betaf*fabs(mneutralinoj)*mq*-rj;
    coupcombo2Zg = 8*Wij*goldstoneffcoup*goldstoneneutineutjcoup*g*sinthetaW*betaf*fabs(mneutralinoi)*mq*ri;
    coupcombo3Zg = 8*Wij*goldstoneffcoup*goldstoneneutineutjcoup*g*sinthetaW*betaf*fabs(mneutralinoj)*mq*-rj;
    coupcombo4Zg = 8*Wij*goldstoneffcoup*goldstoneneutineutjcoup*g*sinthetaW*betaf*fabs(mneutralinoi)*mq*ri;
    
    double integral1Zg = 0, integral2Zg = 0, integral3Zg = 0, integral4Zg = 0;
    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, mA = mZboson, m2 = msf1, MZboson = mZboson;
    integral1Zg = 2*fabs(m1)*dgauss(gintegral1ZAdgauss,fromzsfE,tozsfE,accuracy);
    integral2Zg = 2*fabs(m1)*dgauss(gintegral2ZAdgauss,fromzsfE,tozsfE,accuracy);
    integral3Zg = 2*fabs(m1)*dgauss(gintegral3ZAdgauss,fromzsfE,tozsfE,accuracy);
    integral4Zg = 2*fabs(m1)*dgauss(gintegral4ZAdgauss,fromzsfE,tozsfE,accuracy);

    double GammaZg = 0;
    
    GammaZg = coupcombo1Zg*integral1Zg + coupcombo2Zg*integral2Zg + coupcombo3Zg*integral3Zg + coupcombo4Zg*integral4Zg;

    ///Goldstone-A interference
    double GammagA = 0;

    ///Coupling combos:
    double coupcombogA1 = 0, coupcombogA2 = 0, coupcombogA3 = 0, coupcombogA4 = 0;
    coupcombogA1 = -2*goldstoneneutineutjcoup*(XijA+XjiA)*ri;
    coupcombogA2 = 2*goldstoneneutineutjcoup*(XijA+XjiA)*rj;
    coupcombogA3 = -goldstoneffcoup*Aq;
    coupcombogA4 = goldstoneffcoup*Aq;

    double integralgA1 = 0, integralgA2 = 0, integralgA3 = 0, integralgA4 = 0;
    m1 = mneutralinoi, m4 = mneutralinoj, mq = mf, MZboson = mZboson, mA = mhiggsA;
    integralgA1 = 2*fabs(m1)*dgauss(gneutineutjffgA1dgauss,fromz,toz,accuracy);
    integralgA2 = 2*fabs(m1)*dgauss(gneutineutjffgA2dgauss,fromz,toz,accuracy);
    integralgA3 = 2*fabs(m1)*dgauss(gneutineutjffgA3dgauss,fromz,toz,accuracy);
    integralgA4 = 2*fabs(m1)*dgauss(gneutineutjffgA4dgauss,fromz,toz,accuracy);
    
    GammagA = coupcombogA1*coupcombogA3*integralgA4 - 2*coupcombogA1*coupcombogA4*sqr(mq)*integralgA3 + 2*coupcombogA2*coupcombogA3*fabs(mneutralinoi)*fabs(mneutralinoj)*integralgA2 -4*coupcombogA2*coupcombogA4*sqr(mq)*fabs(mneutralinoi)*fabs(mneutralinoj)*integralgA1;

    amplitudeW = Nc/(512*pow(PI,3)*pow(fabs(mneutralinoi),3))*(GammaZ + Gammah + GammaH + GammaA + GammahHinterf + Gammasftot - 4*Gammahsf1 - 4*Gammahsf2 - 4*GammaHsf1 - 4*GammaHsf2 - 4*GammaAsf1 - 4*GammaAsf2 + 4*GammaZsf1 - 4*GammaZsf2 - 4*GammaZA + 2*GammagA - 4*GammaZg + Gammagoldstone - 4*Gammagsf1 -4*Gammagsf2);    
}
  
  return amplitudeW;
}

double neutralinoamplitudedecaycharginoffprimebar (double mneutralinoi, double msfp1, double msfp2, double msf1, double msf2, double mWboson, double mHP, double mcharginoj, double mfp, double mf, double thetaqp, double thetaq, double g, double gp, double alpha, double beta, double thetaL2, double thetaR2, double runmqp, double runmq, DoubleMatrix & mixNeut, int ineutralino, int jchargino, bool onetothree, char qorl, char norc) {
  
  double amplitudeW = 0;

  if (fabs(mneutralinoi) > msf1 + mf || fabs(mneutralinoi) > msf2 + mf || fabs(mneutralinoi) > msfp1 + mfp || fabs(mneutralinoi) > msfp2 + mfp || fabs(mneutralinoi) > fabs(mcharginoj) + mHP ||fabs(mneutralinoi) > fabs(mcharginoj) + mWboson || fabs(mneutralinoi) < fabs(mcharginoj) + mf + mfp || onetothree == false){
    amplitudeW = 0;
  }

  else{
    double GammaW = 0, GammaHpm = 0, Gammagoldstone = 0, Gammasfp1 = 0, Gammasfp2 = 0, Gammasf1 = 0, Gammasf2 = 0, Gammasfp1sf1 = 0, Gammasfp1sf2 = 0, Gammasfp2sf1 = 0, Gammasfp2sf2 = 0, GammaWHpm = 0, GammaWgoldstone = 0, GammaWSfp1 = 0, GammaWSfp2 = 0, GammaWSf1 = 0, GammaWSf2 = 0, GammaHgoldstone = 0, Gammagsfp1 = 0, Gammagsfp2 = 0, GammaHpmsfp1 = 0, GammaHpmsfp2 = 0, Gammagsf1 = 0, Gammagsf2 = 0, GammaHpmsf1 = 0, GammaHpmsf2 = 0, Gammasfpsfp = 0;

    double charneutWcoupL = 0, charneutWcoupR = 0, coupHpmcharneutL = 0, coupHpmcharneutR = 0, coupHpm1charneutL = 0, coupHpm1charneutR = 0, coupHpm2charneutL = 0, coupHpm2charneutR = 0, fd = 0, fu = 0, coupHpm1ffpd = 0, coupHpm1ffpu = 0, coupHpm2ffpd = 0, coupHpm2ffpu = 0;
    fu = g*runmqp/(root2*sin(beta)*mWboson); ///Just usual yukawa for up type
    fd = g*runmq/(root2*cos(beta)*mWboson); ///Just usual yukawa for down type
    double ri = 0, rj = 0, rc = 0;
    if (mneutralinoi >= 0) { ri = 1;}
    else if (mneutralinoi < 0) { ri = -1;} ///correction factor for negative masses
    if (mcharginoj >= 0) { rj = 1;}
    else if (mcharginoj < 0) { rj = -1;} ///correction factor for negative masses
    if (norc == 'n') { rc = 1;}
    else if (norc == 'c') { rc = -1;} ///correction factor for if it's chargino -> neutralino fp fbar rather than neutralino -> chargino fpbar f

    if ( jchargino == 1) {
      charneutWcoupL = g*(sin(thetaL2)*mixNeut(ineutralino,2) + cos(thetaL2)*mixNeut(ineutralino,3)/(root2));
      charneutWcoupR = g*(sin(thetaR2)*mixNeut(ineutralino,2) - cos(thetaR2)*mixNeut(ineutralino,4)/(root2));
      coupHpmcharneutL = (g*sin(thetaR2)*mixNeut(ineutralino,4) + cos(thetaR2)/(root2)*(gp*mixNeut(ineutralino,1)+g*mixNeut(ineutralino,2)));
      coupHpmcharneutR = (g*sin(thetaL2)*mixNeut(ineutralino,3) - cos(thetaL2)/(root2)*(gp*mixNeut(ineutralino,1)+g*mixNeut(ineutralino,2)));
    }

    else if ( jchargino == 2) {
      charneutWcoupL = g*(cos(thetaL2)*mixNeut(ineutralino,2) - sin(thetaL2)*mixNeut(ineutralino,3)/(root2));
      charneutWcoupR = g*(cos(thetaR2)*mixNeut(ineutralino,2) + sin(thetaR2)*mixNeut(ineutralino,4)/(root2));
      coupHpmcharneutL = (g*cos(thetaR2)*mixNeut(ineutralino,4) - sin(thetaR2)/(root2)*(gp*mixNeut(ineutralino,1)+g*mixNeut(ineutralino,2)));
      coupHpmcharneutR = (g*cos(thetaL2)*mixNeut(ineutralino,3) + sin(thetaL2)/(root2)*(gp*mixNeut(ineutralino,1)+g*mixNeut(ineutralino,2)));
    }
    else {
      throw("problem: jchargino must be 1 or 2 in neutralinoamplitudedecaycharginoffprimebar");
    }

    double AZiu = 0, BZiu = 0, sf1alpha1Ziu = 0, sf1beta1Ziu = 0, sf2alpha1Ziu = 0, sf2beta1Ziu = 0;
    double AZid = 0, BZid = 0, sf1alpha1Zid = 0, sf1beta1Zid = 0, sf2alpha1Zid = 0, sf2beta1Zid = 0;
    double Nc = 0;
    
    if (qorl == 'q') {
      AZiu = g/(root2)*(-mixNeut(ineutralino,2)) + gp/(3*root2)*(-mixNeut(ineutralino,1));
      BZiu = (4./3)*gp/(root2)*(-mixNeut(ineutralino,1));
      sf1alpha1Ziu = AZiu*cos(thetaqp)*rj*-rc*ri - fu*mixNeut(ineutralino,4)*sin(thetaqp);
      sf1beta1Ziu = rc*fu*mixNeut(ineutralino,4)*cos(thetaqp) - BZiu*sin(thetaqp);
      sf2alpha1Ziu = -ri*fu*mixNeut(ineutralino,4)*cos(thetaqp) + rc*AZiu*sin(thetaqp);
      sf2beta1Ziu = -BZiu*cos(thetaqp)*rj*-rc*ri + fu*mixNeut(ineutralino,4)*sin(thetaqp);

      AZid = g/(root2)*(mixNeut(ineutralino,2)) + gp/(3*root2)*(-mixNeut(ineutralino,1));
      BZid = (2./3)*gp/(root2)*(mixNeut(ineutralino,1));
      sf1alpha1Zid = AZid*cos(thetaq)*rj*-rc*ri - fd*mixNeut(ineutralino,4)*sin(thetaq);
      sf1beta1Zid = fd*mixNeut(ineutralino,3)*cos(thetaq) - ri*BZid*sin(thetaq);
      sf2alpha1Zid = fd*mixNeut(ineutralino,3)*cos(thetaq)*rc - rc*ri*AZid*sin(thetaq);
      sf2beta1Zid = -BZid*cos(thetaq)*rj*-rc*ri + fd*mixNeut(ineutralino,4)*sin(thetaq);
            
      Nc = 3;
    }
    else if (qorl == 'l') {
      AZiu = g/(root2)*(-mixNeut(ineutralino,2)) + gp/(root2)*(mixNeut(ineutralino,1));
      BZiu = 0;
      sf1alpha1Ziu = AZiu*cos(thetaqp)*rj*-rc*ri - fu*mixNeut(ineutralino,4)*sin(thetaqp);
      sf1beta1Ziu = -rc*fu*mixNeut(ineutralino,4)*cos(thetaqp) - BZiu*sin(thetaqp);
      sf2alpha1Ziu = 0; ///No snuR exists so only one sfp contribution
      sf2beta1Ziu = 0; ///No snuR exists so only one sfp contribution
         
      AZid = g/(root2)*(mixNeut(ineutralino,2)) + gp/(root2)*(mixNeut(ineutralino,1));
      BZid = root2*gp*mixNeut(ineutralino,1);
      sf1alpha1Zid = AZid*cos(thetaq)*rj*-rc*ri - fd*mixNeut(ineutralino,3)*sin(thetaq);
      sf1beta1Zid = fd*mixNeut(ineutralino,3)*cos(thetaq) - ri*BZid*sin(thetaq);
      sf2alpha1Zid = fd*mixNeut(ineutralino,3)*cos(thetaq)*rc - rc*ri*AZid*sin(thetaq);
      sf2beta1Zid = -BZid*cos(thetaq)*rj*-rc*ri + fd*mixNeut(ineutralino,3)*sin(thetaq);

      Nc = 1;
    }
    else {
      throw("problem: qorl must be q or l in neutralinoamplitudedecaycharginoffprimebar");
    }

    ///GammaW
    double intW1 = 0, intW2 = 0, from = 0, to = 0;
    from = fabs(mcharginoj);
    to = 1/(2*fabs(mneutralinoi))*(pow(mneutralinoi,2)+pow(mcharginoj,2)-pow(mf,2)-pow(mfp,2)-2*mf*mfp);
    m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, MWboson = mWboson;
    intW1 = 2*fabs(m1)*dgauss(gneuticharjffpW1dgauss,from,to,accuracy);
    intW2 = 2*fabs(m1)*dgauss(gneuticharjffpW2dgauss,from,to,accuracy);
    GammaW = -8*charneutWcoupL*charneutWcoupR*pow(g,2)/2*fabs(mneutralinoi)*fabs(mcharginoj)*intW2*ri*rj + 2*(pow(charneutWcoupL,2)+pow(charneutWcoupR,2))*pow(g,2)/2*intW1;

    double coupcombo1Hpm1 = 0, coupcombo2Hpm1 = 0, coupcombo3Hpm1 = 0, coupcombo4Hpm1 = 0, int1Wpm = 0, int2Wpm = 0, int3Wpm = 0, int4Wpm = 0;
    double coupcombo1Hpm2 = 0, coupcombo2Hpm2 = 0, coupcombo3Hpm2 = 0, coupcombo4Hpm2 = 0, int1Hpm = 0, int2Hpm = 0, int3Hpm = 0, int4Hpm = 0;
    ///Hpm1 contribution (W+ goldstone):
    coupHpm1charneutL = coupHpmcharneutL*sin(beta);
    coupHpm1charneutR = coupHpmcharneutR*-cos(beta);

    coupHpm1ffpu = fu*sin(beta);
    coupHpm1ffpd = fd*-cos(beta);

    coupcombo1Hpm1 = pow(coupHpm1charneutL,2) + pow(coupHpm1charneutR,2);
    coupcombo2Hpm1 = coupHpm1charneutL*coupHpm1charneutR*ri;
    coupcombo3Hpm1 = pow(coupHpm1ffpu,2) + pow(coupHpm1ffpd,2);
    coupcombo4Hpm1 = coupHpm1ffpu*coupHpm1ffpd;

    m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson;
    int1Wpm = 2*fabs(m1)*dgauss(gneuticharjffpHpm1dgauss,from,to,accuracy);
    int2Wpm = 2*fabs(m1)*dgauss(gneuticharjffpHpm2dgauss,from,to,accuracy);
    int3Wpm = 2*fabs(m1)*dgauss(gneuticharjffpHpm3dgauss,from,to,accuracy);
    int4Wpm = 2*fabs(m1)*dgauss(gneuticharjffpHpm4dgauss,from,to,accuracy);   

    ///Hpm2 contribution (Actual Hpm):
    coupHpm2charneutL = coupHpmcharneutL*cos(beta);
    coupHpm2charneutR = coupHpmcharneutR*sin(beta);
 
    coupHpm2ffpu = fu*cos(beta);
    coupHpm2ffpd = fd*sin(beta);

    coupcombo1Hpm2 = pow(coupHpm2charneutL,2) + pow(coupHpm2charneutR,2);
    coupcombo2Hpm2 = coupHpm2charneutL*coupHpm2charneutR*ri;
    coupcombo3Hpm2 = pow(coupHpm2ffpu,2) + pow(coupHpm2ffpd,2);
    coupcombo4Hpm2 = coupHpm2ffpu*coupHpm2ffpd;

    m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mHP;
    int1Hpm = 2*fabs(m1)*dgauss(gneuticharjffpHpm1dgauss,from,to,accuracy);
    int2Hpm = 2*fabs(m1)*dgauss(gneuticharjffpHpm2dgauss,from,to,accuracy);
    int3Hpm = 2*fabs(m1)*dgauss(gneuticharjffpHpm3dgauss,from,to,accuracy);
    int4Hpm = 2*fabs(m1)*dgauss(gneuticharjffpHpm4dgauss,from,to,accuracy);  
  
    Gammagoldstone = coupcombo1Hpm1*coupcombo3Hpm1*int4Wpm - 4*coupcombo1Hpm1*coupcombo4Hpm1*int3Wpm*mf*mfp + 4*coupcombo2Hpm1*coupcombo3Hpm1*int2Wpm*fabs(mneutralinoi)*fabs(mcharginoj)*rj - 16*coupcombo2Hpm1*coupcombo4Hpm1*int1Wpm*mf*mfp*fabs(mneutralinoi)*fabs(mcharginoj)*rj;

    GammaHpm = coupcombo1Hpm2*coupcombo3Hpm2*int4Hpm - 4*coupcombo1Hpm2*coupcombo4Hpm2*int3Hpm*mf*mfp + 4*coupcombo2Hpm2*coupcombo3Hpm2*int2Hpm*fabs(mneutralinoi)*fabs(mcharginoj)*rj - 16*coupcombo2Hpm2*coupcombo4Hpm2*int1Hpm*mf*mfp*fabs(mneutralinoi)*fabs(mcharginoj)*rj;
 
    ///Sfp Sfp diagonal (Remember fp is u-type fermion)

    double alphasfp1char = 0, betasfp1char = 0, alphasf1char = 0, betasf1char = 0, alphasfp2char = 0, betasfp2char = 0, alphasf2char = 0, betasf2char = 0;

    if (jchargino == 1) {
      alphasfp1char = -g*sin(thetaR2)*cos(thetaqp) + fu*cos(thetaR2)*sin(thetaqp);
      betasfp1char = -fd*cos(thetaL2)*cos(thetaqp)*rc;
      alphasf1char = -g*sin(thetaL2)*cos(thetaq) + fd*cos(thetaL2)*rc*sin(thetaq);
      betasf1char = -fu*cos(thetaR2)*cos(thetaq);
      alphasfp2char = rc*g*sin(thetaR2)*sin(thetaqp) + fu*cos(thetaR2)*-cos(thetaqp);
      betasfp2char = -fd*cos(thetaL2)*sin(thetaqp)*rc;
      alphasf2char = -fd*cos(thetaL2)*cos(thetaq) + g*sin(thetaL2)*sin(thetaq);
      betasf2char = -fu*cos(thetaR2)*sin(thetaq);
    }
    else if (jchargino == 2) {
      alphasfp1char = -g*cos(thetaR2)*cos(thetaqp) - fu*sin(thetaR2)*sin(thetaqp);
      betasfp1char = fd*sin(thetaL2)*cos(thetaqp)*rc;
      alphasf1char = -g*cos(thetaL2)*cos(thetaq) - fd*sin(thetaL2)*rc*sin(thetaq);
      betasf1char = fu*sin(thetaR2)*cos(thetaq);
      alphasfp2char = rc*g*cos(thetaR2)*sin(thetaqp) - fu*sin(thetaR2)*-cos(thetaq);
      betasfp2char = fd*sin(thetaL2)*sin(thetaqp)*rc;
      alphasf2char = fd*sin(thetaL2)*cos(thetaq) + g*cos(thetaL2)*sin(thetaq);
      betasf2char = fu*sin(thetaR2)*sin(thetaq);
    }
    else {
      throw("problem: jchargino must be 1 or 2 in neutralinoamplitudedecaycharginoffprimebar");
    }

    double coupcombo1sfp1 = 0, coupcombo2sfp1 = 0, coupcombo3sfp1 = 0, coupcombo4sfp1 = 0;
    coupcombo1sfp1 = pow(sf1alpha1Ziu,2) + pow(sf1beta1Ziu,2);
    coupcombo2sfp1 = -sf1alpha1Ziu*sf1beta1Ziu*ri;
    coupcombo3sfp1 = pow(alphasfp1char,2) + pow(betasfp1char,2);
    coupcombo4sfp1 = -alphasfp1char*betasfp1char;

    double Eupper = 0;
    Eupper = 1/(2*fabs(mneutralinoi))*(pow(mneutralinoi,2) + pow(mfp,2) - pow(mf,2) - pow(mcharginoj,2) -2*mf*fabs(mcharginoj));
    double int1sfp1 = 0, int2sfp1 = 0, int3sfp1 = 0, int4sfp1 = 0;
    m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = msfp1;
    int1sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpHpm1dgauss,mfp,Eupper,accuracy);
    int2sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpHpm2dgauss,mfp,Eupper,accuracy);
    int3sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpHpm3dgauss,mfp,Eupper,accuracy);
    int4sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpHpm4dgauss,mfp,Eupper,accuracy);

    Gammasfp1 = coupcombo1sfp1*coupcombo3sfp1*int4sfp1 + 4*coupcombo1sfp1*coupcombo4sfp1*-mf*fabs(mcharginoj)*int3sfp1 + 4*coupcombo2sfp1*coupcombo3sfp1*fabs(mneutralinoi)*mfp*int2sfp1 + 16*coupcombo2sfp1*coupcombo4sfp1*fabs(mneutralinoi)*mfp*-mf*fabs(mcharginoj)*int1sfp1;

    double coupcombo1sfp2 = 0, coupcombo2sfp2 = 0, coupcombo3sfp2 = 0, coupcombo4sfp2 = 0;
    coupcombo1sfp2 = pow(sf2alpha1Ziu,2) + pow(sf2beta1Ziu,2);
    coupcombo2sfp2 = sf2alpha1Ziu*sf2beta1Ziu;
    coupcombo3sfp2 = pow(alphasfp2char,2) + pow(betasfp2char,2);
    coupcombo4sfp2 = alphasfp2char*betasfp2char;

    double int1sfp2 = 0, int2sfp2 = 0, int3sfp2 = 0, int4sfp2 = 0;
    m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = msfp2;
    int1sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpHpm1dgauss,mfp,Eupper,accuracy);
    int2sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpHpm2dgauss,mfp,Eupper,accuracy);
    int3sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpHpm3dgauss,mfp,Eupper,accuracy);
    int4sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpHpm4dgauss,mfp,Eupper,accuracy);

    Gammasfp2 = coupcombo1sfp2*coupcombo3sfp2*int4sfp2 + 4*coupcombo1sfp2*coupcombo4sfp2*-mf*fabs(mcharginoj)*int3sfp2*rc + 4*coupcombo2sfp2*coupcombo3sfp2*fabs(mneutralinoi)*mfp*int2sfp2*rc + 16*coupcombo2sfp2*coupcombo4sfp2*fabs(mneutralinoi)*mfp*-mf*fabs(mcharginoj)*int1sfp2;

    ///Sf Sf diagonal
    double coupcombo1sf1 = 0, coupcombo2sf1 = 0, coupcombo3sf1 = 0, coupcombo4sf1 = 0;
    coupcombo1sf1 = pow(sf1alpha1Zid,2) + pow(sf1beta1Zid,2);
    coupcombo2sf1 = -ri*sf1alpha1Zid*sf1beta1Zid;
    coupcombo3sf1 = pow(alphasf1char,2) + pow(betasf1char,2);
    coupcombo4sf1 = -alphasf1char*betasf1char*rj;

    double Eupper2 = 0;
    Eupper2 = 1/(2*fabs(mneutralinoi))*(pow(mneutralinoi,2) + pow(mf,2) - pow(mfp,2) - pow(mcharginoj,2) -2*mfp*fabs(mcharginoj));
    double int1sf1 = 0, int2sf1 = 0, int3sf1 = 0, int4sf1 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msf1;
    int1sf1 = 2*fabs(m1)*dgauss(gneuticharjffpHpm1dgauss,mf,Eupper2,accuracy);
    int2sf1 = 2*fabs(m1)*dgauss(gneuticharjffpHpm2dgauss,mf,Eupper2,accuracy);
    int3sf1 = 2*fabs(m1)*dgauss(gneuticharjffpHpm3dgauss,mf,Eupper2,accuracy);
    int4sf1 = 2*fabs(m1)*dgauss(gneuticharjffpHpm4dgauss,mf,Eupper2,accuracy);

    Gammasf1 = coupcombo1sf1*coupcombo3sf1*int4sf1 + 4*coupcombo1sf1*coupcombo4sf1*-mfp*fabs(mcharginoj)*int3sf1*rc*rj + 4*coupcombo2sf1*coupcombo3sf1*fabs(mneutralinoi)*mf*int2sf1*rc + 16*coupcombo2sf1*coupcombo4sf1*fabs(mneutralinoi)*mf*-mfp*fabs(mcharginoj)*rj*int1sf1;

    double coupcombo1sf2 = 0, coupcombo2sf2 = 0, coupcombo3sf2 = 0, coupcombo4sf2 = 0;
    coupcombo1sf2 = pow(sf2alpha1Zid,2) + pow(sf2beta1Zid,2);
    coupcombo2sf2 = sf2alpha1Zid*sf2beta1Zid;
    coupcombo3sf2 = pow(alphasf2char,2) + pow(betasf2char,2);
    coupcombo4sf2 = alphasf2char*betasf2char*rc;

    double int1sf2 = 0, int2sf2 = 0, int3sf2 = 0, int4sf2 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msf2;
    int1sf2 = 2*fabs(m1)*dgauss(gneuticharjffpHpm1dgauss,mf,Eupper2,accuracy);
    int2sf2 = 2*fabs(m1)*dgauss(gneuticharjffpHpm2dgauss,mf,Eupper2,accuracy);
    int3sf2 = 2*fabs(m1)*dgauss(gneuticharjffpHpm3dgauss,mf,Eupper2,accuracy);
    int4sf2 = 2*fabs(m1)*dgauss(gneuticharjffpHpm4dgauss,mf,Eupper2,accuracy);

    Gammasf2 = coupcombo1sf2*coupcombo3sf2*int4sf2 + 4*coupcombo1sf2*coupcombo4sf2*-mfp*fabs(mcharginoj)*int3sf2*rc + 4*coupcombo2sf2*coupcombo3sf2*fabs(mneutralinoi)*mf*int2sf2*rc + 16*coupcombo2sf2*coupcombo4sf2*fabs(mneutralinoi)*mf*-mfp*fabs(mcharginoj)*int1sf2;

    m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msf1, m6 = msf2;
    ///Sfp t sf u interference
    ///Sfp 1 sf 2
    
    double coupcombo1sfp1sf2 = 0, coupcombo2sfp1sf2 = 0, coupcombo3sfp1sf2 = 0, coupcombo4sfp1sf2 = 0, coupcombo5sfp1sf2 = 0, coupcombo6sfp1sf2 = 0, coupcombo7sfp1sf2 = 0, coupcombo8sfp1sf2 = 0;
    coupcombo1sfp1sf2 = -0.5*ri*(-sf1alpha1Ziu*sf2beta1Zid*betasfp1char*alphasf2char + sf1beta1Ziu*sf2alpha1Zid*alphasfp1char*betasf2char);
    coupcombo2sfp1sf2 = fabs(mneutralinoi)*fabs(mcharginoj)*(ri*sf1alpha1Ziu*sf2alpha1Zid*alphasfp1char*alphasf2char + sf1beta1Ziu*sf2beta1Zid*betasfp1char*betasf2char);
    coupcombo3sfp1sf2 = ri*mf*mfp*(sf1beta1Ziu*sf2alpha1Zid*betasfp1char*alphasf2char - rc*sf1alpha1Ziu*sf2beta1Zid*alphasfp1char*betasf2char);
    coupcombo8sfp1sf2 = ri*2*fabs(mneutralinoi)*fabs(mcharginoj)*mf*mfp*(rc*sf1beta1Ziu*sf2beta1Zid*alphasfp1char*alphasf2char + ri*sf1alpha1Ziu*sf2alpha1Zid*betasfp1char*betasf2char);
    if (norc == 'n') {
    coupcombo4sfp1sf2 = ri*fabs(mneutralinoi)*mf*(sf1alpha1Ziu*sf2alpha1Zid*betasfp1char*alphasf2char - sf1beta1Ziu*sf2beta1Zid*alphasfp1char*betasf2char);
    coupcombo5sfp1sf2 = mfp*fabs(mcharginoj)*(ri*sf1beta1Ziu*sf2alpha1Zid*alphasfp1char*alphasf2char + sf1alpha1Ziu*sf2beta1Zid*betasfp1char*betasf2char)*rj;
    coupcombo6sfp1sf2 = -ri*fabs(mneutralinoi)*mfp*(-sf1beta1Ziu*sf2beta1Zid*betasfp1char*alphasf2char + sf1alpha1Ziu*sf2alpha1Zid*alphasfp1char*betasf2char);
    coupcombo7sfp1sf2 = fabs(mcharginoj)*mf*(ri*sf1alpha1Ziu*sf2beta1Zid*alphasfp1char*alphasf2char + sf1beta1Ziu*sf2alpha1Zid*betasfp1char*betasf2char)*rj;
    }
    else if (norc == 'c') {
      coupcombo4sfp1sf2 = -(fabs(mneutralinoi)*mf*(-rj*alphasfp1char*alphasf2char*sf1beta1Ziu*sf2alpha1Zid + betasfp1char*betasf2char*sf1alpha1Ziu*sf2beta1Zid));
      coupcombo5sfp1sf2 = -(fabs(mcharginoj)*mfp*(alphasf2char*betasfp1char*sf1alpha1Ziu*sf2alpha1Zid - rj*alphasfp1char*betasf2char*sf1beta1Ziu*sf2beta1Zid));
      coupcombo6sfp1sf2 = -(-fabs(mneutralinoi)*mfp*(-betasfp1char*betasf2char*sf1beta1Ziu*sf2alpha1Zid - alphasfp1char*alphasf2char*sf1alpha1Ziu*sf2beta1Zid));
      coupcombo7sfp1sf2 = -(-fabs(mcharginoj)*mf*(-alphasfp1char*betasf2char*sf1alpha1Ziu*sf2alpha1Zid + rj*betasfp1char*alphasf2char*sf1beta1Ziu*sf2beta1Zid));
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1sfp1sf2 = 0, int2sfp1sf2 = 0, int3sfp1sf2 = 0, int4sfp1sf2 = 0, int5sfp1sf2 = 0, int6sfp1sf2 = 0, int7sfp1sf2 = 0, int8sfp1sf2 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf2;
    int1sfp1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp1sfp1sf2dgauss,mfp,Eupper,accuracy);
    int2sfp1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp2sfp1sf2dgauss,mfp,Eupper,accuracy);
    int3sfp1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp3sfp1sf2dgauss,mfp,Eupper,accuracy);
    int4sfp1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp4sfp1sf2dgauss,mfp,Eupper,accuracy);
    int5sfp1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp5sfp1sf2dgauss,mfp,Eupper,accuracy);
    int6sfp1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp6sfp1sf2dgauss,mfp,Eupper,accuracy);
    int7sfp1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp7sfp1sf2dgauss,mfp,Eupper,accuracy);
    int8sfp1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp8sfp1sf2dgauss,mfp,Eupper,accuracy);

    Gammasfp1sf2 = coupcombo1sfp1sf2*int1sfp1sf2 + coupcombo2sfp1sf2*int2sfp1sf2 + coupcombo3sfp1sf2*int3sfp1sf2 + coupcombo4sfp1sf2*int4sfp1sf2 + coupcombo5sfp1sf2*int5sfp1sf2 + coupcombo6sfp1sf2*int6sfp1sf2 + coupcombo7sfp1sf2*int7sfp1sf2 + coupcombo8sfp1sf2*int8sfp1sf2;

    ///Sfp 1 Sf 1 interference

    double coupcombo1sfp1sf1 = 0, coupcombo2sfp1sf1 = 0, coupcombo3sfp1sf1 = 0, coupcombo4sfp1sf1 = 0, coupcombo5sfp1sf1 = 0, coupcombo6sfp1sf1 = 0, coupcombo7sfp1sf1 = 0, coupcombo8sfp1sf1 = 0;
    coupcombo1sfp1sf1 = -0.5*(ri*sf1alpha1Ziu*sf1beta1Zid*betasfp1char*alphasf1char + sf1beta1Ziu*sf1alpha1Zid*alphasfp1char*betasf1char)*rj;
    coupcombo2sfp1sf1 = -fabs(mneutralinoi)*fabs(mcharginoj)*(sf1alpha1Ziu*sf1alpha1Zid*alphasfp1char*alphasf1char + ri*sf1beta1Ziu*sf1beta1Zid*betasfp1char*betasf1char)*ri*rj;
    coupcombo3sfp1sf1 = -rc*mf*mfp*(sf1beta1Ziu*sf1alpha1Zid*betasfp1char*alphasf1char + ri*sf1alpha1Ziu*sf1beta1Zid*alphasfp1char*betasf1char)*rj;
    coupcombo8sfp1sf1 = -2*rc*fabs(mneutralinoi)*fabs(mcharginoj)*mf*mfp*(ri*sf1beta1Ziu*sf1beta1Zid*alphasfp1char*alphasf1char + sf1alpha1Ziu*sf1alpha1Zid*betasfp1char*betasf1char)*ri*rj;
    if (norc == 'n') {
    coupcombo4sfp1sf1 = fabs(mneutralinoi)*mf*(-sf1alpha1Ziu*sf1alpha1Zid*betasfp1char*alphasf1char - ri*sf1beta1Ziu*sf1beta1Zid*alphasfp1char*betasf1char)*ri;
    coupcombo5sfp1sf1 = -mfp*fabs(mcharginoj)*(sf1beta1Ziu*sf1alpha1Zid*alphasfp1char*alphasf1char + ri*sf1alpha1Ziu*sf1beta1Zid*betasfp1char*betasf1char)*rj;
    coupcombo6sfp1sf1 = fabs(mneutralinoi)*mfp*(-ri*sf1beta1Ziu*sf1beta1Zid*betasfp1char*alphasf1char - sf1alpha1Ziu*sf1alpha1Zid*alphasfp1char*betasf1char)*ri;
    coupcombo7sfp1sf1 = -fabs(mcharginoj)*mf*(ri*sf1alpha1Ziu*sf1beta1Zid*alphasfp1char*alphasf1char + sf1beta1Ziu*sf1alpha1Zid*betasfp1char*betasf1char)*rj;
    }
    else if (norc == 'c') {
      coupcombo4sfp1sf1 = fabs(mneutralinoi)*mf*(alphasfp1char*alphasf1char*sf1beta1Ziu*sf1alpha1Zid + betasfp1char*betasf1char*sf1alpha1Ziu*sf1beta1Zid)*rj;
      coupcombo5sfp1sf1 = -fabs(mcharginoj)*mfp*(alphasf1char*betasfp1char*sf1alpha1Ziu*sf1alpha1Zid + alphasfp1char*betasf1char*sf1beta1Ziu*sf1beta1Zid)*rj;
      coupcombo6sfp1sf1 = fabs(mneutralinoi)*mfp*(-betasfp1char*betasf1char*sf1beta1Ziu*sf1alpha1Zid - alphasfp1char*alphasf1char*sf1alpha1Ziu*sf1beta1Zid)*rj;
      coupcombo7sfp1sf1 = -fabs(mcharginoj)*mf*(-alphasfp1char*betasf1char*sf1alpha1Ziu*sf1alpha1Zid - betasfp1char*alphasf1char*sf1beta1Ziu*sf1beta1Zid)*rj;
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }
    
    double int1sfp1sf1 = 0, int2sfp1sf1 = 0, int3sfp1sf1 = 0, int4sfp1sf1 = 0, int5sfp1sf1 = 0, int6sfp1sf1 = 0, int7sfp1sf1 = 0, int8sfp1sf1 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp1, m6 = msf1;
    int1sfp1sf1 = 2*fabs(m1)*dgauss(gneuticharjffp1sfp1sf2dgauss,mfp,Eupper,accuracy);
    int2sfp1sf1 = 2*fabs(m1)*dgauss(gneuticharjffp2sfp1sf2dgauss,mfp,Eupper,accuracy);
    int3sfp1sf1 = 2*fabs(m1)*dgauss(gneuticharjffp3sfp1sf2dgauss,mfp,Eupper,accuracy);
    int4sfp1sf1 = 2*fabs(m1)*dgauss(gneuticharjffp4sfp1sf2dgauss,mfp,Eupper,accuracy);
    int5sfp1sf1 = 2*fabs(m1)*dgauss(gneuticharjffp5sfp1sf2dgauss,mfp,Eupper,accuracy);
    int6sfp1sf1 = 2*fabs(m1)*dgauss(gneuticharjffp6sfp1sf2dgauss,mfp,Eupper,accuracy);
    int7sfp1sf1 = 2*fabs(m1)*dgauss(gneuticharjffp7sfp1sf2dgauss,mfp,Eupper,accuracy);
    int8sfp1sf1 = 2*fabs(m1)*dgauss(gneuticharjffp8sfp1sf2dgauss,mfp,Eupper,accuracy);

    Gammasfp1sf1 = coupcombo1sfp1sf1*int1sfp1sf1 + coupcombo2sfp1sf1*int2sfp1sf1 + coupcombo3sfp1sf1*int3sfp1sf1 + coupcombo4sfp1sf1*int4sfp1sf1 + coupcombo5sfp1sf1*int5sfp1sf1 + coupcombo6sfp1sf1*int6sfp1sf1 + coupcombo7sfp1sf1*int7sfp1sf1 + coupcombo8sfp1sf1*int8sfp1sf1;

    ///Sfp 2 Sf 2 interference
    
    double coupcombo1sfp2sf2 = 0, coupcombo2sfp2sf2 = 0, coupcombo3sfp2sf2 = 0, coupcombo4sfp2sf2 = 0, coupcombo5sfp2sf2 = 0, coupcombo6sfp2sf2 = 0, coupcombo7sfp2sf2 = 0, coupcombo8sfp2sf2 = 0;
    coupcombo8sfp2sf2 = -2*fabs(mneutralinoi)*fabs(mcharginoj)*mf*mfp*(sf2beta1Ziu*sf2beta1Zid*alphasfp2char*alphasf2char - sf2alpha1Ziu*sf2alpha1Zid*betasfp2char*betasf2char)*ri*rj;

    if (norc == 'n') {
      coupcombo1sfp2sf2 = 0.5*(-sf2alpha1Ziu*sf2beta1Zid*betasfp2char*alphasf2char + ri*sf2beta1Ziu*sf2alpha1Zid*alphasfp2char*betasf2char)*rj;
      coupcombo2sfp2sf2 = fabs(mneutralinoi)*fabs(mcharginoj)*(sf2alpha1Ziu*sf2alpha1Zid*alphasfp2char*alphasf2char - sf2beta1Ziu*sf2beta1Zid*betasfp2char*betasf2char)*rj;
      coupcombo3sfp2sf2 = -mf*mfp*(-ri*sf2beta1Ziu*sf2alpha1Zid*betasfp2char*alphasf2char + sf2alpha1Ziu*sf2beta1Zid*alphasfp2char*betasf2char);
      coupcombo4sfp2sf2 = -fabs(mneutralinoi)*mf*(sf2alpha1Ziu*sf2alpha1Zid*betasfp2char*alphasf2char - ri*sf2beta1Ziu*sf2beta1Zid*alphasfp2char*betasf2char)*rj;
      coupcombo5sfp2sf2 = -mfp*fabs(mcharginoj)*(sf2beta1Ziu*sf2alpha1Zid*alphasfp2char*alphasf2char + ri*rj*sf2alpha1Ziu*sf2beta1Zid*betasfp2char*betasf2char)*rj*ri;
      coupcombo6sfp2sf2 = fabs(mneutralinoi)*mfp*(ri*sf2beta1Ziu*sf2beta1Zid*betasfp2char*alphasf2char - sf2alpha1Ziu*sf2alpha1Zid*alphasfp2char*betasf2char);
      coupcombo7sfp2sf2 = fabs(mcharginoj)*mf*(sf2alpha1Ziu*sf2beta1Zid*alphasfp2char*alphasf2char - sf2beta1Ziu*sf2alpha1Zid*betasfp2char*betasf2char)*rj;
 
    }
    else if (norc == 'c') {
      coupcombo1sfp2sf2 = -0.5*(-rj*sf2alpha1Ziu*sf2beta1Zid*betasfp2char*alphasf2char - sf2beta1Ziu*sf2alpha1Zid*alphasfp2char*betasf2char);
      coupcombo2sfp2sf2 = -rj*fabs(mneutralinoi)*fabs(mcharginoj)*(sf2alpha1Ziu*sf2alpha1Zid*alphasfp2char*alphasf2char - rj*sf2beta1Ziu*sf2beta1Zid*betasfp2char*betasf2char)*ri;
      coupcombo3sfp2sf2 = mf*mfp*(-rj*sf2beta1Ziu*sf2alpha1Zid*betasfp2char*alphasf2char - sf2alpha1Ziu*sf2beta1Zid*alphasfp2char*betasf2char);
      coupcombo4sfp2sf2 = rj*fabs(mneutralinoi)*mf*(alphasfp2char*alphasf2char*sf2beta1Ziu*sf2alpha1Zid + betasfp2char*betasf2char*sf2alpha1Ziu*sf2beta1Zid);
      coupcombo5sfp2sf2 = rj*fabs(mcharginoj)*mfp*(alphasf2char*betasfp2char*sf2alpha1Ziu*sf2alpha1Zid - alphasfp2char*betasf2char*sf2beta1Ziu*sf2beta1Zid)*rj;
      coupcombo6sfp2sf2 = -rj*fabs(mneutralinoi)*mfp*(-betasfp2char*betasf2char*sf2beta1Ziu*sf2alpha1Zid - alphasfp2char*alphasf2char*sf2alpha1Ziu*sf2beta1Zid);
      coupcombo7sfp2sf2 = fabs(mcharginoj)*mf*(-alphasfp2char*betasf2char*sf2alpha1Ziu*sf2alpha1Zid + betasfp2char*alphasf2char*sf2beta1Ziu*sf2beta1Zid);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }
    
    double int1sfp2sf2 = 0, int2sfp2sf2 = 0, int3sfp2sf2 = 0, int4sfp2sf2 = 0, int5sfp2sf2 = 0, int6sfp2sf2 = 0, int7sfp2sf2 = 0, int8sfp2sf2 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp2, m6 = msf2;
    int1sfp2sf2 = 2*fabs(m1)*dgauss(gneuticharjffp1sfp1sf2dgauss,mfp,Eupper,accuracy);
    int2sfp2sf2 = 2*fabs(m1)*dgauss(gneuticharjffp2sfp1sf2dgauss,mfp,Eupper,accuracy);
    int3sfp2sf2 = 2*fabs(m1)*dgauss(gneuticharjffp3sfp1sf2dgauss,mfp,Eupper,accuracy);
    int4sfp2sf2 = 2*fabs(m1)*dgauss(gneuticharjffp4sfp1sf2dgauss,mfp,Eupper,accuracy);
    int5sfp2sf2 = 2*fabs(m1)*dgauss(gneuticharjffp5sfp1sf2dgauss,mfp,Eupper,accuracy);
    int6sfp2sf2 = 2*fabs(m1)*dgauss(gneuticharjffp6sfp1sf2dgauss,mfp,Eupper,accuracy);
    int7sfp2sf2 = 2*fabs(m1)*dgauss(gneuticharjffp7sfp1sf2dgauss,mfp,Eupper,accuracy);
    int8sfp2sf2 = 2*fabs(m1)*dgauss(gneuticharjffp8sfp1sf2dgauss,mfp,Eupper,accuracy);

    Gammasfp2sf2 = coupcombo1sfp2sf2*int1sfp2sf2 + coupcombo2sfp2sf2*int2sfp2sf2 + coupcombo3sfp2sf2*int3sfp2sf2 + coupcombo4sfp2sf2*int4sfp2sf2 + coupcombo5sfp2sf2*int5sfp2sf2 + coupcombo6sfp2sf2*int6sfp2sf2 + coupcombo7sfp2sf2*int7sfp2sf2 + coupcombo8sfp2sf2*int8sfp2sf2;

    ///Sfp 2 Sf 1 interference
    
    double coupcombo1sfp2sf1 = 0, coupcombo2sfp2sf1 = 0, coupcombo3sfp2sf1 = 0, coupcombo4sfp2sf1 = 0, coupcombo5sfp2sf1 = 0, coupcombo6sfp2sf1 = 0, coupcombo7sfp2sf1 = 0, coupcombo8sfp2sf1 = 0;
    coupcombo8sfp2sf1 = ri*rc*2*fabs(mneutralinoi)*fabs(mcharginoj)*mf*mfp*(sf2beta1Ziu*sf1beta1Zid*alphasfp2char*alphasf1char - sf2alpha1Ziu*sf1alpha1Zid*betasfp2char*betasf1char);
    if (norc == 'n') {
      coupcombo1sfp2sf1 = ri*0.5*(-sf2alpha1Ziu*sf1beta1Zid*betasfp2char*alphasf1char + sf2beta1Ziu*sf1alpha1Zid*alphasfp2char*betasf1char);
      coupcombo2sfp2sf1 = -fabs(mneutralinoi)*fabs(mcharginoj)*(sf2alpha1Ziu*sf1alpha1Zid*alphasfp2char*alphasf1char - sf2beta1Ziu*sf1beta1Zid*betasfp2char*betasf1char);
      coupcombo3sfp2sf1 = ri*mf*mfp*(-sf2beta1Ziu*sf1alpha1Zid*betasfp2char*alphasf1char + sf2alpha1Ziu*sf1beta1Zid*alphasfp2char*betasf1char);
      coupcombo4sfp2sf1 = fabs(mneutralinoi)*mf*(sf2alpha1Ziu*sf1alpha1Zid*betasfp2char*alphasf1char + ri*sf2beta1Ziu*sf1beta1Zid*alphasfp2char*betasf1char);
      coupcombo5sfp2sf1 = ri*mfp*fabs(mcharginoj)*(sf2beta1Ziu*sf1alpha1Zid*alphasfp2char*alphasf1char - sf2alpha1Ziu*sf1beta1Zid*betasfp2char*betasf1char)*rj;
      coupcombo6sfp2sf1 = -ri*fabs(mneutralinoi)*mfp*(sf2beta1Ziu*sf1beta1Zid*betasfp2char*alphasf1char + ri*sf2alpha1Ziu*sf1alpha1Zid*alphasfp2char*betasf1char);
      coupcombo7sfp2sf1 = -ri*fabs(mcharginoj)*mf*(ri*sf2alpha1Ziu*sf1beta1Zid*alphasfp2char*alphasf1char + sf2beta1Ziu*sf1alpha1Zid*betasfp2char*betasf1char)*rj;
    }
    else if (norc == 'c') {
      coupcombo1sfp2sf1 = -0.5*(-sf2alpha1Ziu*sf1beta1Zid*betasfp2char*alphasf1char - sf2beta1Ziu*sf1alpha1Zid*alphasfp2char*betasf1char);
      coupcombo2sfp2sf1 = fabs(mneutralinoi)*fabs(mcharginoj)*(sf2alpha1Ziu*sf1alpha1Zid*alphasfp2char*alphasf1char + sf2beta1Ziu*sf1beta1Zid*betasfp2char*betasf1char)*ri;
      coupcombo3sfp2sf1 = mf*mfp*(-sf2beta1Ziu*sf1alpha1Zid*betasfp2char*alphasf1char - sf2alpha1Ziu*sf1beta1Zid*alphasfp2char*betasf1char);
      coupcombo4sfp2sf1 = -fabs(mneutralinoi)*mf*(alphasfp2char*alphasf1char*sf2beta1Ziu*sf1alpha1Zid + betasfp2char*betasf1char*sf2alpha1Ziu*sf1beta1Zid);
      coupcombo5sfp2sf1 = fabs(mcharginoj)*mfp*(alphasf1char*betasfp2char*sf2alpha1Ziu*sf1alpha1Zid + alphasfp2char*betasf1char*sf2beta1Ziu*sf1beta1Zid);
      coupcombo6sfp2sf1 = -fabs(mneutralinoi)*mfp*(-betasfp2char*betasf1char*sf2beta1Ziu*sf1alpha1Zid - alphasfp2char*alphasf1char*sf2alpha1Ziu*sf1beta1Zid);
      coupcombo7sfp2sf1 = fabs(mcharginoj)*mf*(-alphasfp2char*betasf1char*sf2alpha1Ziu*sf1alpha1Zid - betasfp2char*alphasf1char*sf2beta1Ziu*sf1beta1Zid);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1sfp2sf1 = 0, int2sfp2sf1 = 0, int3sfp2sf1 = 0, int4sfp2sf1 = 0, int5sfp2sf1 = 0, int6sfp2sf1 = 0, int7sfp2sf1 = 0, int8sfp2sf1 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msfp2, m6 = msf1;
    int1sfp2sf1 = 2*fabs(m1)*dgauss(gneuticharjffp1sfp1sf2dgauss,mfp,Eupper,accuracy);
    int2sfp2sf1 = 2*fabs(m1)*dgauss(gneuticharjffp2sfp1sf2dgauss,mfp,Eupper,accuracy);
    int3sfp2sf1 = 2*fabs(m1)*dgauss(gneuticharjffp3sfp1sf2dgauss,mfp,Eupper,accuracy);
    int4sfp2sf1 = 2*fabs(m1)*dgauss(gneuticharjffp4sfp1sf2dgauss,mfp,Eupper,accuracy);
    int5sfp2sf1 = 2*fabs(m1)*dgauss(gneuticharjffp5sfp1sf2dgauss,mfp,Eupper,accuracy);
    int6sfp2sf1 = 2*fabs(m1)*dgauss(gneuticharjffp6sfp1sf2dgauss,mfp,Eupper,accuracy);
    int7sfp2sf1 = 2*fabs(m1)*dgauss(gneuticharjffp7sfp1sf2dgauss,mfp,Eupper,accuracy);
    int8sfp2sf1 = 2*fabs(m1)*dgauss(gneuticharjffp8sfp1sf2dgauss,mfp,Eupper,accuracy);

    Gammasfp2sf1 = coupcombo1sfp2sf1*int1sfp2sf1 + coupcombo2sfp2sf1*int2sfp2sf1 + coupcombo3sfp2sf1*int3sfp2sf1 + coupcombo4sfp2sf1*int4sfp2sf1 + coupcombo5sfp2sf1*int5sfp2sf1 + coupcombo6sfp2sf1*int6sfp2sf1 + coupcombo7sfp2sf1*int7sfp2sf1 + coupcombo8sfp2sf1*int8sfp2sf1;

    ///W-Hpm interference
    double coupcombo1WHpm = 0, coupcombo2WHpm = 0, coupcombo3WHpm = 0, coupcombo4WHpm = 0;
    coupcombo1WHpm = (charneutWcoupR*coupHpm2charneutR + charneutWcoupL*coupHpm2charneutL)*-g/(root2)*coupHpm2ffpu*fabs(mcharginoj)*mfp*rc;
    coupcombo2WHpm = (charneutWcoupL*coupHpm2charneutR + charneutWcoupR*coupHpm2charneutL)*g/(root2)*-coupHpm2ffpd*fabs(mneutralinoi)*-mf*ri*rc*rj;
    coupcombo3WHpm = (charneutWcoupR*coupHpm2charneutR + charneutWcoupL*coupHpm2charneutL)*g/(root2)*-coupHpm2ffpd*fabs(mcharginoj)*-mf*rc;
    coupcombo4WHpm = (charneutWcoupL*coupHpm2charneutR + charneutWcoupR*coupHpm2charneutL)*-g/(root2)*coupHpm2ffpu*fabs(mneutralinoi)*mfp*ri*rc*rj;

    double int1WHpm = 0, int2WHpm = 0, int3WHpm = 0, int4WHpm = 0;
    double Eupper3 = 0;
    Eupper3 = 1/(2*fabs(mneutralinoi))*(pow(mneutralinoi,2) + pow(mcharginoj,2) - pow(mf,2) - pow(mfp,2) - 2*mf*mfp);
    m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHP;
    int1WHpm = 2*fabs(m1)*dgauss(gneuticharjffp1WHpmdgauss,fabs(mcharginoj), Eupper3, accuracy);
    int2WHpm = 2*fabs(m1)*dgauss(gneuticharjffp2WHpmdgauss,fabs(mcharginoj), Eupper3, accuracy);
    int3WHpm = 2*fabs(m1)*dgauss(gneuticharjffp3WHpmdgauss,fabs(mcharginoj), Eupper3, accuracy);
    int4WHpm = 2*fabs(m1)*dgauss(gneuticharjffp4WHpmdgauss,fabs(mcharginoj), Eupper3, accuracy);

    GammaWHpm = coupcombo1WHpm*int1WHpm*rc*rj + coupcombo2WHpm*int2WHpm*rc*rj + coupcombo3WHpm*int3WHpm + coupcombo4WHpm*int4WHpm;

    ///W-goldstone interference
    double coupcombo1Wg = 0, coupcombo2Wg = 0, coupcombo3Wg = 0, coupcombo4Wg = 0;
    coupcombo1Wg = (charneutWcoupR*coupHpm1charneutR + charneutWcoupL*coupHpm1charneutL)*-g/(root2)*coupHpm1ffpu*fabs(mcharginoj)*mfp*rc;
    coupcombo2Wg = (charneutWcoupL*coupHpm1charneutR + charneutWcoupR*coupHpm1charneutL)*g/(root2)*-coupHpm1ffpd*fabs(mneutralinoi)*-mf*ri*rc*rj;
    coupcombo3Wg = (charneutWcoupR*coupHpm1charneutR + charneutWcoupL*coupHpm1charneutL)*g/(root2)*-coupHpm1ffpd*fabs(mcharginoj)*-mf*rc;
    coupcombo4Wg = (charneutWcoupL*coupHpm1charneutR + charneutWcoupR*coupHpm1charneutL)*-g/(root2)*coupHpm1ffpu*fabs(mneutralinoi)*mfp*ri*rc*rj;

    double int1Wg = 0, int2Wg = 0, int3Wg = 0, int4Wg = 0;
    Eupper3 = 1/(2*fabs(mneutralinoi))*(pow(mneutralinoi,2) + pow(mcharginoj,2) - pow(mf,2) - pow(mfp,2) - 2*mf*mfp);
    m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mWboson;
    int1Wg = 2*fabs(m1)*dgauss(gneuticharjffp1WHpmdgauss,fabs(mcharginoj), Eupper3, accuracy);
    int2Wg = 2*fabs(m1)*dgauss(gneuticharjffp2WHpmdgauss,fabs(mcharginoj), Eupper3, accuracy);
    int3Wg = 2*fabs(m1)*dgauss(gneuticharjffp3WHpmdgauss,fabs(mcharginoj), Eupper3, accuracy);
    int4Wg = 2*fabs(m1)*dgauss(gneuticharjffp4WHpmdgauss,fabs(mcharginoj), Eupper3, accuracy);

    GammaWgoldstone = coupcombo1Wg*int1Wg + coupcombo2Wg*int2Wg + coupcombo3Wg*int3Wg + coupcombo4Wg*int4Wg;
    
    ///W Sfp 1 interference
    double coupcombo1Wsfp1 = 0, coupcombo2Wsfp1 = 0, coupcombo3Wsfp1 = 0, coupcombo4Wsfp1 = 0, coupcombo5Wsfp1 = 0, coupcombo6Wsfp1 = 0, coupcombo7Wsfp1 = 0, coupcombo8Wsfp1 = 0;
    coupcombo7Wsfp1 = -2*charneutWcoupL*sf1beta1Ziu*g/(root2)*betasfp1char*mfp*mf;
    coupcombo8Wsfp1 = 8*charneutWcoupR*sf1beta1Ziu*g/(root2)*betasfp1char*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj)*ri*rj;
    coupcombo1Wsfp1 = 2*charneutWcoupL*sf1alpha1Ziu*-g/(root2)*betasfp1char*fabs(mneutralinoi)*mf*rj;
    coupcombo2Wsfp1 = rc*2*charneutWcoupL*-sf1beta1Ziu*g/(root2)*alphasfp1char*mfp*fabs(mcharginoj)*rj;
    coupcombo3Wsfp1 = rc*2*charneutWcoupR*sf1alpha1Ziu*g/(root2)*alphasfp1char*ri*rj;
    coupcombo4Wsfp1 = rc*ri*4*charneutWcoupR*sf1beta1Ziu*g/(root2)*alphasfp1char*fabs(mneutralinoi)*mfp;
    coupcombo5Wsfp1 = 4*charneutWcoupR*sf1alpha1Ziu*g/(root2)*betasfp1char*mf*fabs(mcharginoj)*ri;
    coupcombo6Wsfp1 = -rc*2*charneutWcoupL*sf1alpha1Ziu*g/(root2)*alphasfp1char*fabs(mneutralinoi)*fabs(mcharginoj);

    double intW1Sfp1 = 0, intW2Sfp1 = 0, intW3Sfp1 = 0, intW4Sfp1 = 0, intW5Sfp1 = 0, intW6Sfp1 = 0, intW7Sfp1 = 0, intW8Sfp1 = 0;
    m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1;
    intW1Sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpW1Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW2Sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpW2Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW3Sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpW3Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW4Sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpW4Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW5Sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpW5Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW6Sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpW6Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW7Sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpW7Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW8Sfp1 = 2*fabs(m1)*dgauss(gneuticharjffpW8Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);

    GammaWSfp1 = coupcombo1Wsfp1*intW1Sfp1 + coupcombo2Wsfp1*intW2Sfp1 + coupcombo3Wsfp1*intW3Sfp1 + coupcombo4Wsfp1*intW4Sfp1 + coupcombo5Wsfp1*intW5Sfp1 + coupcombo6Wsfp1*intW6Sfp1 + coupcombo7Wsfp1*intW7Sfp1 + coupcombo8Wsfp1*intW8Sfp1;

    /// W Sfp 2 interference
    double coupcombo1Wsfp2 = 0, coupcombo2Wsfp2 = 0, coupcombo3Wsfp2 = 0, coupcombo4Wsfp2 = 0, coupcombo5Wsfp2 = 0, coupcombo6Wsfp2 = 0, coupcombo7Wsfp2 = 0, coupcombo8Wsfp2 = 0;
    coupcombo1Wsfp2 = rc*2*charneutWcoupL*sf2alpha1Ziu*-g/(root2)*betasfp2char*fabs(mneutralinoi)*-mf*rj;
    coupcombo2Wsfp2 = -ri*2*charneutWcoupL*-sf2beta1Ziu*g/(root2)*alphasfp2char*mfp*fabs(mcharginoj)*rj;
    coupcombo3Wsfp2 = 2*charneutWcoupR*sf2alpha1Ziu*g/(root2)*alphasfp2char*ri;
    coupcombo4Wsfp2 = -4*charneutWcoupR*sf2beta1Ziu*g/(root2)*alphasfp2char*fabs(mneutralinoi)*mfp;
    coupcombo5Wsfp2 = -rc*4*charneutWcoupR*sf2alpha1Ziu*g/(root2)*betasfp2char*mf*fabs(mcharginoj)*ri;
    coupcombo6Wsfp2 = -2*charneutWcoupL*sf2alpha1Ziu*g/(root2)*alphasfp2char*fabs(mneutralinoi)*fabs(mcharginoj)*rj;
    coupcombo7Wsfp2 = 2*charneutWcoupL*sf2beta1Ziu*g/(root2)*betasfp2char*mfp*mf*rj*-ri*rc;
    coupcombo8Wsfp2 = 8*charneutWcoupR*sf2beta1Ziu*g/(root2)*betasfp2char*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj)*rc;

    double intW1Sfp2 = 0, intW2Sfp2 = 0, intW3Sfp2 = 0, intW4Sfp2 = 0, intW5Sfp2 = 0, intW6Sfp2 = 0, intW7Sfp2 = 0, intW8Sfp2 = 0;
    m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp2;
    intW1Sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpW1Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW2Sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpW2Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW3Sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpW3Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW4Sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpW4Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW5Sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpW5Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW6Sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpW6Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW7Sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpW7Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW8Sfp2 = 2*fabs(m1)*dgauss(gneuticharjffpW8Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);

    GammaWSfp2 = coupcombo1Wsfp2*intW1Sfp2 + coupcombo2Wsfp2*intW2Sfp2 + coupcombo3Wsfp2*intW3Sfp2 + coupcombo4Wsfp2*intW4Sfp2 + coupcombo5Wsfp2*intW5Sfp2 + coupcombo6Wsfp2*intW6Sfp2 + coupcombo7Wsfp2*intW7Sfp2 + coupcombo8Wsfp2*intW8Sfp2;

    ///W Sf1 interference
    double coupcombo1Wsf1 = 0, coupcombo2Wsf1 = 0, coupcombo3Wsf1 = 0, coupcombo4Wsf1 = 0, coupcombo5Wsf1 = 0, coupcombo6Wsf1 = 0, coupcombo7Wsf1 = 0, coupcombo8Wsf1 = 0;

    coupcombo6Wsf1 = 2*charneutWcoupR*sf1alpha1Zid*g/(root2)*alphasf1char*fabs(mneutralinoi)*fabs(mcharginoj)*-rc;
    coupcombo7Wsf1 = -2*charneutWcoupR*sf1beta1Zid*g/(root2)*betasf1char*mf*mfp;
    coupcombo8Wsf1 = -8*charneutWcoupL*sf1beta1Zid*g/(root2)*betasf1char*fabs(mneutralinoi)*mf*mfp*fabs(mcharginoj)*rj*-ri;
    if (norc == 'n') {
      coupcombo1Wsf1 = -2*charneutWcoupR*sf1alpha1Zid*-g/(root2)*betasf1char*fabs(mneutralinoi)*-mfp;
      coupcombo2Wsf1 = 2*charneutWcoupR*-sf1beta1Zid*g/(root2)*alphasf1char*mf*fabs(mcharginoj)*rj;
      coupcombo3Wsf1 = 2*charneutWcoupR*sf1alpha1Zid*g/(root2)*alphasf1char;
      coupcombo4Wsf1 = ri*4*charneutWcoupL*sf1beta1Zid*g/(root2)*alphasf1char*fabs(mneutralinoi)*mf;
      coupcombo5Wsf1 = ri*4*charneutWcoupL*sf1alpha1Zid*g/(root2)*betasf1char*mfp*fabs(mcharginoj)*rj;
    }
    else if (norc == 'c') {
      coupcombo1Wsf1 = 2*charneutWcoupR*alphasf1char*-g/root2*sf1beta1Zid*fabs(mneutralinoi)*-mfp*ri;
      coupcombo2Wsf1 = 2*charneutWcoupR*betasf1char*-g/root2*sf1alpha1Zid*mf*fabs(mcharginoj);
      coupcombo3Wsf1 = 2*charneutWcoupL*alphasf1char*-g/root2*sf1alpha1Zid*ri*rj;
      coupcombo4Wsf1 = -4*charneutWcoupL*betasf1char*-g/root2*sf1alpha1Zid*fabs(mneutralinoi)*mf*ri*rj;
      coupcombo5Wsf1 = -4*charneutWcoupL*alphasf1char*g/root2*sf1beta1Zid*mfp*fabs(mcharginoj)*rj;
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double intW1Sf1 = 0, intW2Sf1 = 0, intW3Sf1 = 0, intW4Sf1 = 0, intW5Sf1 = 0, intW6Sf1 = 0, intW7Sf1 = 0, intW8Sf1 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mfp, m4 = mcharginoj, m5 = mWboson, m6 = msf1;
    intW1Sf1 = 2*fabs(m1)*dgauss(gneuticharjffpW1Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW2Sf1 = 2*fabs(m1)*dgauss(gneuticharjffpW2Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW3Sf1 = 2*fabs(m1)*dgauss(gneuticharjffpW3Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW4Sf1 = 2*fabs(m1)*dgauss(gneuticharjffpW4Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW5Sf1 = 2*fabs(m1)*dgauss(gneuticharjffpW5Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW6Sf1 = 2*fabs(m1)*dgauss(gneuticharjffpW6Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW7Sf1 = 2*fabs(m1)*dgauss(gneuticharjffpW7Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW8Sf1 = 2*fabs(m1)*dgauss(gneuticharjffpW8Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);

    GammaWSf1 = coupcombo1Wsf1*intW1Sf1 + coupcombo2Wsf1*intW2Sf1 + coupcombo3Wsf1*intW3Sf1 + coupcombo4Wsf1*intW4Sf1 + coupcombo5Wsf1*intW5Sf1 + coupcombo6Wsf1*intW6Sf1 + coupcombo7Wsf1*intW7Sf1 + coupcombo8Wsf1*intW8Sf1;

    ///W Sf2 interference
    double coupcombo1Wsf2 = 0, coupcombo2Wsf2 = 0, coupcombo3Wsf2 = 0, coupcombo4Wsf2 = 0, coupcombo5Wsf2 = 0, coupcombo6Wsf2 = 0, coupcombo7Wsf2 = 0, coupcombo8Wsf2 = 0;
    coupcombo6Wsf2 = 2*charneutWcoupR*sf2alpha1Zid*g/(root2)*alphasf2char*fabs(mneutralinoi)*fabs(mcharginoj)*rj*rc;
    coupcombo7Wsf2 = -2*charneutWcoupR*sf2beta1Zid*g/(root2)*betasf2char*mf*mfp*rc;
    coupcombo8Wsf2 = 8*charneutWcoupL*sf2beta1Zid*g/(root2)*betasf2char*fabs(mneutralinoi)*mf*mfp*fabs(mcharginoj)*rj*rc*ri;
    if (norc == 'n') {
      coupcombo1Wsf2 = -2*charneutWcoupR*sf2alpha1Zid*-g/(root2)*betasf2char*fabs(mneutralinoi)*-mfp;
      coupcombo2Wsf2 = -2*charneutWcoupR*-sf2beta1Zid*g/(root2)*alphasf2char*mf*fabs(mcharginoj)*rj;
      coupcombo3Wsf2 = -2*charneutWcoupR*sf2alpha1Zid*g/(root2)*alphasf2char;
      coupcombo4Wsf2 = -4*charneutWcoupL*sf2beta1Zid*g/(root2)*alphasf2char*fabs(mneutralinoi)*mf*ri;
      coupcombo5Wsf2 = 4*charneutWcoupL*sf2alpha1Zid*g/(root2)*betasf2char*mfp*fabs(mcharginoj)*rj*ri;
    }
    else if (norc == 'c') {
      coupcombo1Wsf2 = 2*charneutWcoupR*alphasf2char*-g/root2*sf2beta1Zid*fabs(mneutralinoi)*mfp*ri*rj*rc;
      coupcombo2Wsf2 = 2*charneutWcoupR*betasf2char*g/root2*sf2alpha1Zid*mf*fabs(mcharginoj)*rc;
      coupcombo3Wsf2 = -2*charneutWcoupL*alphasf2char*g/root2*sf2alpha1Zid*ri*rc;
      coupcombo4Wsf2 = -4*charneutWcoupL*betasf2char*g/root2*sf2alpha1Zid*fabs(mneutralinoi)*mf*ri*rj*rc;
      coupcombo5Wsf2 = 4*charneutWcoupL*alphasf2char*g/root2*sf2beta1Zid*mfp*fabs(mcharginoj)*rc;
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double intW1Sf2 = 0, intW2Sf2 = 0, intW3Sf2 = 0, intW4Sf2 = 0, intW5Sf2 = 0, intW6Sf2 = 0, intW7Sf2 = 0, intW8Sf2 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mfp, m4 = mcharginoj, m5 = mWboson, m6 = msf2;
    intW1Sf2 = 2*fabs(m1)*dgauss(gneuticharjffpW1Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW2Sf2 = 2*fabs(m1)*dgauss(gneuticharjffpW2Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW3Sf2 = 2*fabs(m1)*dgauss(gneuticharjffpW3Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW4Sf2 = 2*fabs(m1)*dgauss(gneuticharjffpW4Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW5Sf2 = 2*fabs(m1)*dgauss(gneuticharjffpW5Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW6Sf2 = 2*fabs(m1)*dgauss(gneuticharjffpW6Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW7Sf2 = 2*fabs(m1)*dgauss(gneuticharjffpW7Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    intW8Sf2 = 2*fabs(m1)*dgauss(gneuticharjffpW8Sfpdgauss,fabs(mcharginoj),Eupper3,accuracy);

    GammaWSf2 = coupcombo1Wsf2*intW1Sf2 + coupcombo2Wsf2*intW2Sf2 + coupcombo3Wsf2*intW3Sf2 + coupcombo4Wsf2*intW4Sf2 + coupcombo5Wsf2*intW5Sf2 + coupcombo6Wsf2*intW6Sf2 + coupcombo7Wsf2*intW7Sf2 + coupcombo8Wsf2*intW8Sf2;

    ///H+ goldstone interference
    double coupcombo1Hg = 0, coupcombo2Hg = 0, coupcombo3Hg = 0, coupcombo4Hg = 0;
    coupcombo1Hg = coupHpm1charneutL*coupHpm2charneutL + coupHpm1charneutR*coupHpm2charneutR;
    coupcombo2Hg = (coupHpm1charneutR*coupHpm2charneutL + coupHpm1charneutL*coupHpm2charneutR)*ri*rj;
    coupcombo3Hg = coupHpm1ffpu*coupHpm2ffpu + coupHpm1ffpd*coupHpm2ffpd;
    coupcombo4Hg = coupHpm1ffpd*coupHpm2ffpu + coupHpm1ffpu*coupHpm2ffpd;

    double int1Hg = 0, int2Hg = 0, int3Hg = 0, int4Hg = 0;
    m1 = mneutralinoi, m2 = mcharginoj, m3 = mf, m4 = mfp, m5 = mWboson, m6 = mHP;
    int1Hg = 2*fabs(m1)*dgauss(gneuticharjffpHg1dgauss,fabs(mcharginoj),Eupper3,accuracy);
    int2Hg = 2*fabs(m1)*dgauss(gneuticharjffpHg2dgauss,fabs(mcharginoj),Eupper3,accuracy);
    int3Hg = 2*fabs(m1)*dgauss(gneuticharjffpHg3dgauss,fabs(mcharginoj),Eupper3,accuracy);
    int4Hg = 2*fabs(m1)*dgauss(gneuticharjffpHg4dgauss,fabs(mcharginoj),Eupper3,accuracy);

    GammaHgoldstone = coupcombo1Hg*coupcombo3Hg*int4Hg - 2*coupcombo1Hg*coupcombo4Hg*mf*mfp*int3Hg + 2*coupcombo2Hg*coupcombo3Hg*fabs(mneutralinoi)*fabs(mcharginoj)*int2Hg - 4*coupcombo2Hg*coupcombo4Hg*fabs(mneutralinoi)*fabs(mcharginoj)*mf*mfp*int1Hg*rj;
    
    ///goldstone - sfp1 interference
    double coupcombo1gsfp1 = 0, coupcombo2gsfp1 = 0, coupcombo3gsfp1 = 0, coupcombo4gsfp1 = 0, coupcombo5gsfp1 = 0, coupcombo6gsfp1 = 0, coupcombo7gsfp1 = 0, coupcombo8gsfp1 = 0; 
    if (norc == 'n') {
      coupcombo1gsfp1 = 0.5*(coupHpm1charneutR*sf1alpha1Ziu*coupHpm1ffpd*betasfp1char + coupHpm1charneutL*sf1beta1Ziu*coupHpm1ffpu*alphasfp1char);
      coupcombo2gsfp1 = (coupHpm1charneutR*sf1beta1Ziu*coupHpm1ffpu*betasfp1char + coupHpm1charneutL*sf1alpha1Ziu*coupHpm1ffpd*alphasfp1char)*fabs(mneutralinoi)*-mf;
      coupcombo3gsfp1 = (coupHpm1charneutL*-sf1beta1Ziu*coupHpm1ffpd*betasfp1char - coupHpm1charneutR*sf1alpha1Ziu*coupHpm1ffpu*alphasfp1char)*mfp*fabs(mcharginoj)*-ri;
      coupcombo4gsfp1 = (coupHpm1charneutR*-sf1beta1Ziu*coupHpm1ffpd*betasfp1char + coupHpm1charneutL*sf1alpha1Ziu*coupHpm1ffpu*-alphasfp1char)*fabs(mneutralinoi)*-mfp;
      coupcombo5gsfp1 = -(-coupHpm1charneutL*sf1beta1Ziu*coupHpm1ffpu*betasfp1char + ri*coupHpm1charneutR*sf1alpha1Ziu*coupHpm1ffpd*alphasfp1char)*mf*fabs(mcharginoj);
      coupcombo6gsfp1 = ri*(coupHpm1charneutL*sf1alpha1Ziu*coupHpm1ffpd*betasfp1char + coupHpm1charneutR*sf1beta1Ziu*coupHpm1ffpu*alphasfp1char)*fabs(mneutralinoi)*fabs(mcharginoj);
      coupcombo7gsfp1 = -(coupHpm1charneutR*sf1alpha1Ziu*coupHpm1ffpu*betasfp1char + coupHpm1charneutL*sf1beta1Ziu*coupHpm1ffpd*alphasfp1char)*mfp*mf;
      coupcombo8gsfp1 = 2*(coupHpm1charneutL*sf1alpha1Ziu*coupHpm1ffpu*betasfp1char - ri*coupHpm1charneutR*sf1beta1Ziu*coupHpm1ffpd*alphasfp1char)*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj);
    }
    else if (norc == 'c') {
      coupcombo1gsfp1 = -0.5*(coupHpm1charneutL*betasfp1char*coupHpm1ffpu*sf1alpha1Ziu + rj*coupHpm1charneutR*alphasfp1char*coupHpm1ffpd*sf1beta1Ziu);
      coupcombo2gsfp1 = (coupHpm1charneutL*alphasfp1char*coupHpm1ffpd*sf1alpha1Ziu + rj*coupHpm1charneutR*betasfp1char*coupHpm1ffpu*sf1beta1Ziu)*fabs(mneutralinoi)*mfp;
      coupcombo3gsfp1 = -(coupHpm1charneutR*alphasfp1char*coupHpm1ffpu*sf1alpha1Ziu + rj*coupHpm1charneutL*betasfp1char*coupHpm1ffpd*sf1beta1Ziu)*mf*fabs(mcharginoj)*rj;
      coupcombo4gsfp1 = -(coupHpm1charneutL*alphasfp1char*coupHpm1ffpu*sf1alpha1Ziu + rj*coupHpm1charneutR*betasfp1char*coupHpm1ffpd*sf1beta1Ziu)*fabs(mneutralinoi)*mf;
      coupcombo5gsfp1 = (rj*coupHpm1charneutR*alphasfp1char*coupHpm1ffpd*sf1alpha1Ziu + coupHpm1charneutL*betasfp1char*coupHpm1ffpu*sf1beta1Ziu)*mfp*fabs(mcharginoj);
      coupcombo6gsfp1 = -(rj*coupHpm1charneutR*betasfp1char*coupHpm1ffpu*sf1alpha1Ziu + coupHpm1charneutL*alphasfp1char*coupHpm1ffpd*sf1beta1Ziu)*fabs(mcharginoj)*fabs(mneutralinoi);
      coupcombo7gsfp1 = (rj*coupHpm1charneutL*betasfp1char*coupHpm1ffpd*sf1alpha1Ziu + coupHpm1charneutR*alphasfp1char*coupHpm1ffpu*sf1beta1Ziu)*mf*mfp*rj;
      coupcombo8gsfp1 = 2*(rj*coupHpm1charneutR*betasfp1char*coupHpm1ffpd*sf1alpha1Ziu + coupHpm1charneutL*alphasfp1char*coupHpm1ffpu*sf1beta1Ziu)*mf*mfp*fabs(mcharginoj)*fabs(mneutralinoi);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1gsfp1 = 0, int2gsfp1 = 0, int3gsfp1 = 0, int4gsfp1 = 0, int5gsfp1 = 0, int6gsfp1 = 0, int7gsfp1 = 0, int8gsfp1 = 0;
    m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp1;
    if (norc == 'n') {
      int1gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else if (norc == 'c') { ///swap which integral goes with which couplingcombo
      int1gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8gsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    Gammagsfp1 = coupcombo1gsfp1*int1gsfp1 + coupcombo2gsfp1*int2gsfp1 + coupcombo3gsfp1*int3gsfp1 + coupcombo4gsfp1*int4gsfp1 + coupcombo5gsfp1*int5gsfp1 + coupcombo6gsfp1*int6gsfp1 + coupcombo7gsfp1*int7gsfp1 + coupcombo8gsfp1*int8gsfp1;

    ///goldstone - sfp2 interference
    double coupcombo1gsfp2 = 0, coupcombo2gsfp2 = 0, coupcombo3gsfp2 = 0, coupcombo4gsfp2 = 0, coupcombo5gsfp2 = 0, coupcombo6gsfp2 = 0, coupcombo7gsfp2 = 0, coupcombo8gsfp2 = 0; 
    if (norc == 'n') {
      coupcombo1gsfp2 = -0.5*(ri*coupHpm1charneutR*sf2alpha1Ziu*coupHpm1ffpd*-betasfp2char + coupHpm1charneutL*sf2beta1Ziu*coupHpm1ffpu*alphasfp2char);
      coupcombo2gsfp2 = -ri*(coupHpm1charneutR*sf2beta1Ziu*coupHpm1ffpu*-betasfp2char + coupHpm1charneutL*sf2alpha1Ziu*coupHpm1ffpd*alphasfp2char)*fabs(mneutralinoi)*mf;
      coupcombo3gsfp2 = -(coupHpm1charneutL*-sf2beta1Ziu*coupHpm1ffpd*-betasfp2char - coupHpm1charneutR*sf2alpha1Ziu*coupHpm1ffpu*alphasfp2char)*mfp*fabs(mcharginoj)*rj;
      coupcombo4gsfp2 = -(coupHpm1charneutR*-sf2beta1Ziu*coupHpm1ffpd*-betasfp2char + ri*coupHpm1charneutL*sf2alpha1Ziu*coupHpm1ffpu*-alphasfp2char)*fabs(mneutralinoi)*mfp;
      coupcombo5gsfp2 = -(coupHpm1charneutL*sf2beta1Ziu*coupHpm1ffpu*-betasfp2char - coupHpm1charneutR*sf2alpha1Ziu*coupHpm1ffpd*alphasfp2char)*mf*fabs(mcharginoj)*rj;
      coupcombo6gsfp2 = -(coupHpm1charneutL*sf2alpha1Ziu*coupHpm1ffpd*-betasfp2char + ri*coupHpm1charneutR*sf2beta1Ziu*coupHpm1ffpu*alphasfp2char)*fabs(mneutralinoi)*fabs(mcharginoj);
      coupcombo7gsfp2 = (ri*coupHpm1charneutR*sf2alpha1Ziu*coupHpm1ffpu*-betasfp2char + coupHpm1charneutL*sf2beta1Ziu*coupHpm1ffpd*alphasfp2char)*mfp*mf;
      coupcombo8gsfp2 = 2*(coupHpm1charneutL*sf2alpha1Ziu*coupHpm1ffpu*-betasfp2char + ri*coupHpm1charneutR*sf2beta1Ziu*coupHpm1ffpd*alphasfp2char)*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj);
    }
    else if (norc == 'c') {
      coupcombo1gsfp2 = 0.5*(coupHpm1charneutL*betasfp2char*coupHpm1ffpu*sf2alpha1Ziu*rj + coupHpm1charneutR*alphasfp2char*coupHpm1ffpd*sf2beta1Ziu);
      coupcombo2gsfp2 = -(coupHpm1charneutL*alphasfp2char*coupHpm1ffpd*sf2alpha1Ziu + rj*coupHpm1charneutR*betasfp2char*coupHpm1ffpu*sf2beta1Ziu)*fabs(mneutralinoi)*mfp*rj;
      coupcombo3gsfp2 = (coupHpm1charneutR*alphasfp2char*coupHpm1ffpu*sf2alpha1Ziu + rj*coupHpm1charneutL*betasfp2char*coupHpm1ffpd*sf2beta1Ziu)*mf*fabs(mcharginoj);
      coupcombo4gsfp2 = (coupHpm1charneutL*alphasfp2char*coupHpm1ffpu*sf2alpha1Ziu + rj*coupHpm1charneutR*betasfp2char*coupHpm1ffpd*sf2beta1Ziu)*fabs(mneutralinoi)*mf*rj;
      coupcombo5gsfp2 = -(coupHpm1charneutR*alphasfp2char*coupHpm1ffpd*sf2alpha1Ziu + rj*coupHpm1charneutL*betasfp2char*coupHpm1ffpu*sf2beta1Ziu)*mfp*fabs(mcharginoj);
      coupcombo6gsfp2 = (rj*coupHpm1charneutR*betasfp2char*coupHpm1ffpu*sf2alpha1Ziu + coupHpm1charneutL*alphasfp2char*coupHpm1ffpd*sf2beta1Ziu)*fabs(mcharginoj)*fabs(mneutralinoi)*rj;
      coupcombo7gsfp2 = -(rj*coupHpm1charneutL*betasfp2char*coupHpm1ffpd*sf2alpha1Ziu + coupHpm1charneutR*alphasfp2char*coupHpm1ffpu*sf2beta1Ziu)*mf*mfp;
      coupcombo8gsfp2 = -2*(rj*coupHpm1charneutR*betasfp2char*coupHpm1ffpd*sf2alpha1Ziu + coupHpm1charneutL*alphasfp2char*coupHpm1ffpu*sf2beta1Ziu)*mf*mfp*fabs(mcharginoj)*fabs(mneutralinoi)*rj;
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1gsfp2 = 0, int2gsfp2 = 0, int3gsfp2 = 0, int4gsfp2 = 0, int5gsfp2 = 0, int6gsfp2 = 0, int7gsfp2 = 0, int8gsfp2 = 0;
    m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mWboson, m6 = msfp2;
    if (norc == 'n') {
      int1gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else if (norc == 'c') {
      int1gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8gsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    Gammagsfp2 = coupcombo1gsfp2*int1gsfp2 + coupcombo2gsfp2*int2gsfp2 + coupcombo3gsfp2*int3gsfp2 + coupcombo4gsfp2*int4gsfp2 + coupcombo5gsfp2*int5gsfp2 + coupcombo6gsfp2*int6gsfp2 + coupcombo7gsfp2*int7gsfp2 + coupcombo8gsfp2*int8gsfp2;

    ///H+ - sfp1 interference
    double coupcombo1Hpmsfp1 = 0, coupcombo2Hpmsfp1 = 0, coupcombo3Hpmsfp1 = 0, coupcombo4Hpmsfp1 = 0, coupcombo5Hpmsfp1 = 0, coupcombo6Hpmsfp1 = 0, coupcombo7Hpmsfp1 = 0, coupcombo8Hpmsfp1 = 0; 
    if (norc == 'n') {
      coupcombo1Hpmsfp1 = 0.5*(coupHpm2charneutR*sf1alpha1Ziu*coupHpm2ffpd*betasfp1char - ri*coupHpm2charneutL*sf1beta1Ziu*coupHpm2ffpu*alphasfp1char)*-ri;
      coupcombo2Hpmsfp1 = (coupHpm2charneutR*sf1beta1Ziu*coupHpm2ffpu*betasfp1char + coupHpm2charneutL*sf1alpha1Ziu*coupHpm2ffpd*alphasfp1char)*fabs(mneutralinoi)*-mf;
      coupcombo3Hpmsfp1 = (coupHpm2charneutL*-sf1beta1Ziu*coupHpm2ffpd*betasfp1char - coupHpm2charneutR*sf1alpha1Ziu*coupHpm2ffpu*alphasfp1char)*mfp*fabs(mcharginoj)*-ri;
      coupcombo4Hpmsfp1 = (-ri*coupHpm2charneutR*-sf1beta1Ziu*coupHpm2ffpd*betasfp1char + coupHpm2charneutL*sf1alpha1Ziu*coupHpm2ffpu*-alphasfp1char)*fabs(mneutralinoi)*-mfp;
      coupcombo5Hpmsfp1 = (coupHpm2charneutL*sf1beta1Ziu*coupHpm2ffpu*betasfp1char + coupHpm2charneutR*sf1alpha1Ziu*coupHpm2ffpd*alphasfp1char)*mf*fabs(mcharginoj)*-ri;
      coupcombo6Hpmsfp1 = ri*(coupHpm2charneutL*sf1alpha1Ziu*coupHpm2ffpd*betasfp1char + coupHpm2charneutR*sf1beta1Ziu*coupHpm2ffpu*alphasfp1char)*fabs(mneutralinoi)*fabs(mcharginoj);
      coupcombo7Hpmsfp1 = -(coupHpm2charneutR*sf1alpha1Ziu*coupHpm2ffpu*betasfp1char + coupHpm2charneutL*sf1beta1Ziu*coupHpm2ffpd*alphasfp1char)*mfp*mf;
      coupcombo8Hpmsfp1 = -2*ri*(coupHpm2charneutL*sf1alpha1Ziu*coupHpm2ffpu*betasfp1char + coupHpm2charneutR*sf1beta1Ziu*coupHpm2ffpd*alphasfp1char)*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj);
    }
    else if (norc == 'c') {
      coupcombo1Hpmsfp1 = rj*0.5*(coupHpm2charneutL*betasfp1char*coupHpm2ffpu*sf1alpha1Ziu + rj*coupHpm2charneutR*alphasfp1char*coupHpm2ffpd*sf1beta1Ziu);
      coupcombo2Hpmsfp1 = (coupHpm2charneutL*alphasfp1char*coupHpm2ffpd*sf1alpha1Ziu + rj*coupHpm2charneutR*betasfp1char*coupHpm2ffpu*sf1beta1Ziu)*fabs(mneutralinoi)*mfp;
      coupcombo3Hpmsfp1 = -(rj*coupHpm2charneutR*alphasfp1char*coupHpm2ffpu*sf1alpha1Ziu + coupHpm2charneutL*betasfp1char*coupHpm2ffpd*sf1beta1Ziu)*mf*fabs(mcharginoj);
      coupcombo4Hpmsfp1 = -(coupHpm2charneutL*alphasfp1char*coupHpm2ffpu*sf1alpha1Ziu + rj*coupHpm2charneutR*betasfp1char*coupHpm2ffpd*sf1beta1Ziu)*fabs(mneutralinoi)*mf;
      coupcombo5Hpmsfp1 = (rj*coupHpm2charneutR*alphasfp1char*coupHpm2ffpd*sf1alpha1Ziu + coupHpm2charneutL*betasfp1char*coupHpm2ffpu*sf1beta1Ziu)*mfp*fabs(mcharginoj);
      coupcombo6Hpmsfp1 = (coupHpm2charneutR*betasfp1char*coupHpm2ffpu*sf1alpha1Ziu + rj*coupHpm2charneutL*alphasfp1char*coupHpm2ffpd*sf1beta1Ziu)*fabs(mcharginoj)*fabs(mneutralinoi);
      coupcombo7Hpmsfp1 = -(rj*coupHpm2charneutL*betasfp1char*coupHpm2ffpd*sf1alpha1Ziu + coupHpm2charneutR*alphasfp1char*coupHpm2ffpu*sf1beta1Ziu)*mf*mfp;
      coupcombo8Hpmsfp1 = -2*(coupHpm2charneutR*betasfp1char*coupHpm2ffpd*sf1alpha1Ziu + rj*coupHpm2charneutL*alphasfp1char*coupHpm2ffpu*sf1beta1Ziu)*mf*mfp*fabs(mcharginoj)*fabs(mneutralinoi);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1Hpmsfp1 = 0, int2Hpmsfp1 = 0, int3Hpmsfp1 = 0, int4Hpmsfp1 = 0, int5Hpmsfp1 = 0, int6Hpmsfp1 = 0, int7Hpmsfp1 = 0, int8Hpmsfp1 = 0;
    m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mHP, m6 = msfp1;
    if (norc == 'n') {
      int1Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else if (norc == 'c') {
      int1Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8Hpmsfp1 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    GammaHpmsfp1 = coupcombo1Hpmsfp1*int1Hpmsfp1 + coupcombo2Hpmsfp1*int2Hpmsfp1 + coupcombo3Hpmsfp1*int3Hpmsfp1 + coupcombo4Hpmsfp1*int4Hpmsfp1 + coupcombo5Hpmsfp1*int5Hpmsfp1 + coupcombo6Hpmsfp1*int6Hpmsfp1 + coupcombo7Hpmsfp1*int7Hpmsfp1 + coupcombo8Hpmsfp1*int8Hpmsfp1;

    ///H+ - sfp2 interference
    double coupcombo1Hpmsfp2 = 0, coupcombo2Hpmsfp2 = 0, coupcombo3Hpmsfp2 = 0, coupcombo4Hpmsfp2 = 0, coupcombo5Hpmsfp2 = 0, coupcombo6Hpmsfp2 = 0, coupcombo7Hpmsfp2 = 0, coupcombo8Hpmsfp2 = 0; 
    if (norc == 'n') {
      coupcombo1Hpmsfp2 = -0.5*(-ri*coupHpm2charneutR*sf2alpha1Ziu*coupHpm2ffpd*-betasfp2char + coupHpm2charneutL*sf2beta1Ziu*coupHpm2ffpu*alphasfp2char);
      coupcombo2Hpmsfp2 = ri*(coupHpm2charneutR*sf2beta1Ziu*coupHpm2ffpu*-betasfp2char + coupHpm2charneutL*sf2alpha1Ziu*coupHpm2ffpd*alphasfp2char)*fabs(mneutralinoi)*-mf;
      coupcombo3Hpmsfp2 = -ri*(coupHpm2charneutL*-sf2beta1Ziu*coupHpm2ffpd*-betasfp2char - ri*coupHpm2charneutR*sf2alpha1Ziu*coupHpm2ffpu*alphasfp2char)*mfp*fabs(mcharginoj);
      coupcombo4Hpmsfp2 = -(ri*coupHpm2charneutR*-sf2beta1Ziu*coupHpm2ffpd*-betasfp2char - coupHpm2charneutL*sf2alpha1Ziu*coupHpm2ffpu*-alphasfp2char)*fabs(mneutralinoi)*mfp;
      coupcombo5Hpmsfp2 = -(coupHpm2charneutL*sf2beta1Ziu*coupHpm2ffpu*-betasfp2char + coupHpm2charneutR*sf2alpha1Ziu*coupHpm2ffpd*alphasfp2char)*mf*fabs(mcharginoj)*rj;
      coupcombo6Hpmsfp2 = -(coupHpm2charneutL*sf2alpha1Ziu*coupHpm2ffpd*-betasfp2char + ri*coupHpm2charneutR*sf2beta1Ziu*coupHpm2ffpu*alphasfp2char)*fabs(mneutralinoi)*fabs(mcharginoj);
      coupcombo7Hpmsfp2 = (ri*coupHpm2charneutR*sf2alpha1Ziu*coupHpm2ffpu*-betasfp2char + coupHpm2charneutL*sf2beta1Ziu*coupHpm2ffpd*alphasfp2char)*mfp*mf;
      coupcombo8Hpmsfp2 = 2*ri*(coupHpm2charneutL*sf2alpha1Ziu*coupHpm2ffpu*-betasfp2char + coupHpm2charneutR*sf2beta1Ziu*coupHpm2ffpd*alphasfp2char)*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj);
    }
    else if (norc == 'c') {
      coupcombo1Hpmsfp2 = 0.5*(rj*coupHpm2charneutL*betasfp2char*coupHpm2ffpu*sf2alpha1Ziu + coupHpm2charneutR*alphasfp2char*coupHpm2ffpd*sf2beta1Ziu);
      coupcombo2Hpmsfp2 = -(rj*coupHpm2charneutL*alphasfp2char*coupHpm2ffpd*sf2alpha1Ziu + coupHpm2charneutR*betasfp2char*coupHpm2ffpu*sf2beta1Ziu)*fabs(mneutralinoi)*mfp;
      coupcombo3Hpmsfp2 = (coupHpm2charneutR*alphasfp2char*coupHpm2ffpu*sf2alpha1Ziu + rj*coupHpm2charneutL*betasfp2char*coupHpm2ffpd*sf2beta1Ziu)*mf*fabs(mcharginoj);
      coupcombo4Hpmsfp2 = (rj*coupHpm2charneutL*alphasfp2char*coupHpm2ffpu*sf2alpha1Ziu + coupHpm2charneutR*betasfp2char*coupHpm2ffpd*sf2beta1Ziu)*fabs(mneutralinoi)*mf;
      coupcombo5Hpmsfp2 = -(coupHpm2charneutR*alphasfp2char*coupHpm2ffpd*sf2alpha1Ziu + rj*coupHpm2charneutL*betasfp2char*coupHpm2ffpu*sf2beta1Ziu)*mfp*fabs(mcharginoj);
      coupcombo6Hpmsfp2 = (coupHpm2charneutR*betasfp2char*coupHpm2ffpu*sf2alpha1Ziu + rj*coupHpm2charneutL*alphasfp2char*coupHpm2ffpd*sf2beta1Ziu)*fabs(mcharginoj)*fabs(mneutralinoi);
      coupcombo7Hpmsfp2 = -(rj*coupHpm2charneutL*betasfp2char*coupHpm2ffpd*sf2alpha1Ziu + coupHpm2charneutR*alphasfp2char*coupHpm2ffpu*sf2beta1Ziu)*mf*mfp;
      coupcombo8Hpmsfp2 = -2*(coupHpm2charneutR*betasfp2char*coupHpm2ffpd*sf2alpha1Ziu + rj*coupHpm2charneutL*alphasfp2char*coupHpm2ffpu*sf2beta1Ziu)*mf*mfp*fabs(mcharginoj)*fabs(mneutralinoi);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1Hpmsfp2 = 0, int2Hpmsfp2 = 0, int3Hpmsfp2 = 0, int4Hpmsfp2 = 0, int5Hpmsfp2 = 0, int6Hpmsfp2 = 0, int7Hpmsfp2 = 0, int8Hpmsfp2 = 0;
    m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = mHP, m6 = msfp2;
    if (norc == 'n') {
      int1Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else if (norc == 'c') {
      int1Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8Hpmsfp2 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    GammaHpmsfp2 = coupcombo1Hpmsfp2*int1Hpmsfp2 + coupcombo2Hpmsfp2*int2Hpmsfp2 + coupcombo3Hpmsfp2*int3Hpmsfp2 + coupcombo4Hpmsfp2*int4Hpmsfp2 + coupcombo5Hpmsfp2*int5Hpmsfp2 + coupcombo6Hpmsfp2*int6Hpmsfp2 + coupcombo7Hpmsfp2*int7Hpmsfp2 + coupcombo8Hpmsfp2*int8Hpmsfp2;

    ///goldstone - sf1 interference
    double coupcombo1gsf1 = 0, coupcombo2gsf1 = 0, coupcombo3gsf1 = 0, coupcombo4gsf1 = 0, coupcombo5gsf1 = 0, coupcombo6gsf1 = 0, coupcombo7gsf1 = 0, coupcombo8gsf1 = 0; 
    if (norc == 'n') {
      coupcombo1gsf1 = -ri*0.5*(-ri*coupHpm1charneutR*-sf1beta1Zid*coupHpm1ffpd*-alphasf1char + coupHpm1charneutL*-sf1alpha1Zid*coupHpm1ffpu*-betasf1char);
      coupcombo2gsf1 = (coupHpm1charneutR*-sf1alpha1Zid*coupHpm1ffpu*-alphasf1char + ri*coupHpm1charneutL*-sf1beta1Zid*coupHpm1ffpd*-betasf1char)*fabs(mneutralinoi)*-mfp;
      coupcombo3gsf1 = (coupHpm1charneutL*sf1alpha1Zid*coupHpm1ffpd*-alphasf1char - ri*coupHpm1charneutR*-sf1beta1Zid*coupHpm1ffpu*-betasf1char)*mf*fabs(mcharginoj)*-ri;
      coupcombo4gsf1 = (coupHpm1charneutR*sf1alpha1Zid*coupHpm1ffpd*-alphasf1char + ri*coupHpm1charneutL*sf1beta1Zid*coupHpm1ffpu*betasf1char)*fabs(mneutralinoi)*-mf;
      coupcombo5gsf1 = (coupHpm1charneutL*-sf1alpha1Zid*coupHpm1ffpu*-alphasf1char + ri*coupHpm1charneutR*-sf1beta1Zid*coupHpm1ffpd*-betasf1char)*mfp*fabs(mcharginoj)*-ri;
      coupcombo6gsf1 = (-coupHpm1charneutL*-sf1beta1Zid*coupHpm1ffpd*-alphasf1char + ri*coupHpm1charneutR*-sf1alpha1Zid*coupHpm1ffpu*-betasf1char)*fabs(mneutralinoi)*fabs(mcharginoj)*-ri;
      coupcombo7gsf1 = (-ri*coupHpm1charneutR*-sf1beta1Zid*coupHpm1ffpu*-alphasf1char + coupHpm1charneutL*-sf1alpha1Zid*coupHpm1ffpd*-betasf1char)*mfp*mf*ri;
      coupcombo8gsf1 = -2*(ri*coupHpm1charneutL*-sf1beta1Zid*coupHpm1ffpu*-alphasf1char + coupHpm1charneutR*-sf1alpha1Zid*coupHpm1ffpd*-betasf1char)*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj);
    }
    else if (norc == 'c') {
      coupcombo1gsf1 = 0.5*(rj*coupHpm1charneutL*betasf1char*coupHpm1ffpu*sf1alpha1Zid + coupHpm1charneutR*alphasf1char*coupHpm1ffpd*sf1beta1Zid);
      coupcombo2gsf1 = (rj*coupHpm1charneutL*alphasf1char*coupHpm1ffpd*sf1alpha1Zid + coupHpm1charneutR*betasf1char*coupHpm1ffpu*sf1beta1Zid)*fabs(mneutralinoi)*mf;
      coupcombo3gsf1 = -(coupHpm1charneutR*alphasf1char*coupHpm1ffpu*sf1alpha1Zid + rj*coupHpm1charneutL*betasf1char*coupHpm1ffpd*sf1beta1Zid)*mfp*fabs(mcharginoj);
      coupcombo4gsf1 = -(rj*coupHpm1charneutL*alphasf1char*coupHpm1ffpu*sf1alpha1Zid + coupHpm1charneutR*betasf1char*coupHpm1ffpd*sf1beta1Zid)*fabs(mneutralinoi)*mfp;
      coupcombo5gsf1 = (coupHpm1charneutR*alphasf1char*coupHpm1ffpd*sf1alpha1Zid + rj*coupHpm1charneutL*betasf1char*coupHpm1ffpu*sf1beta1Zid)*mf*fabs(mcharginoj);
      coupcombo6gsf1 = (coupHpm1charneutR*betasf1char*coupHpm1ffpu*sf1alpha1Zid + rj*coupHpm1charneutL*alphasf1char*coupHpm1ffpd*sf1beta1Zid)*fabs(mcharginoj)*fabs(mneutralinoi);
      coupcombo7gsf1 = -(rj*coupHpm1charneutL*betasf1char*coupHpm1ffpd*sf1alpha1Zid + coupHpm1charneutR*alphasf1char*coupHpm1ffpu*sf1beta1Zid)*mf*mfp;
      coupcombo8gsf1 = -2*(coupHpm1charneutR*betasf1char*coupHpm1ffpd*sf1alpha1Zid + rj*coupHpm1charneutL*alphasf1char*coupHpm1ffpu*sf1beta1Zid)*mf*mfp*fabs(mcharginoj)*fabs(mneutralinoi);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1gsf1 = 0, int2gsf1 = 0, int3gsf1 = 0, int4gsf1 = 0, int5gsf1 = 0, int6gsf1 = 0, int7gsf1 = 0, int8gsf1 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mfp, m4 = mcharginoj, m5 = mWboson, m6 = msf1;
    if (norc == 'n') {
      int1gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else if (norc == 'c') {
      int1gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8gsf1 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    Gammagsf1 = coupcombo1gsf1*int1gsf1 + coupcombo2gsf1*int2gsf1 + coupcombo3gsf1*int3gsf1 + coupcombo4gsf1*int4gsf1 + coupcombo5gsf1*int5gsf1 + coupcombo6gsf1*int6gsf1 + coupcombo7gsf1*int7gsf1 + coupcombo8gsf1*int8gsf1;

    ///goldstone - sf2 interference
    double coupcombo1gsf2 = 0, coupcombo2gsf2 = 0, coupcombo3gsf2 = 0, coupcombo4gsf2 = 0, coupcombo5gsf2 = 0, coupcombo6gsf2 = 0, coupcombo7gsf2 = 0, coupcombo8gsf2 = 0; 
    if (norc == 'n') {
      coupcombo1gsf2 = -0.5*(coupHpm1charneutR*-sf2beta1Zid*coupHpm1ffpd*-alphasf2char + ri*coupHpm1charneutL*-sf2alpha1Zid*coupHpm1ffpu*-betasf2char);
      coupcombo2gsf2 = -(coupHpm1charneutR*-sf2alpha1Zid*coupHpm1ffpu*-alphasf2char + ri*coupHpm1charneutL*-sf2beta1Zid*coupHpm1ffpd*-betasf2char)*fabs(mneutralinoi)*-mfp;
      coupcombo3gsf2 = (ri*coupHpm1charneutL*sf2alpha1Zid*coupHpm1ffpd*-alphasf2char - coupHpm1charneutR*-sf2beta1Zid*coupHpm1ffpu*-betasf2char)*mf*fabs(mcharginoj);
      coupcombo4gsf2 = -(coupHpm1charneutR*sf2alpha1Zid*coupHpm1ffpd*-alphasf2char + ri*coupHpm1charneutL*-sf2beta1Zid*coupHpm1ffpu*betasf2char)*fabs(mneutralinoi)*-mf;
      coupcombo5gsf2 = (ri*coupHpm1charneutL*-sf2alpha1Zid*coupHpm1ffpu*-alphasf2char + coupHpm1charneutR*-sf2beta1Zid*coupHpm1ffpd*-betasf2char)*mfp*fabs(mcharginoj);
      coupcombo6gsf2 = -(ri*coupHpm1charneutL*-sf2beta1Zid*coupHpm1ffpd*-alphasf2char + coupHpm1charneutR*-sf2alpha1Zid*coupHpm1ffpu*-betasf2char)*fabs(mneutralinoi)*fabs(mcharginoj)*rj;
      coupcombo7gsf2 = (coupHpm1charneutR*-sf2beta1Zid*coupHpm1ffpu*-alphasf2char + ri*coupHpm1charneutL*-sf2alpha1Zid*coupHpm1ffpd*-betasf2char)*mfp*mf;
      coupcombo8gsf2 = 2*(ri*coupHpm1charneutL*-sf2beta1Zid*coupHpm1ffpu*-alphasf2char + coupHpm1charneutR*-sf2alpha1Zid*coupHpm1ffpd*-betasf2char)*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj)*rj;
    }
    else if (norc == 'c') {
      coupcombo1gsf2 = -0.5*rj*(coupHpm1charneutL*betasf2char*coupHpm1ffpu*sf2alpha1Zid - coupHpm1charneutR*alphasf2char*coupHpm1ffpd*sf2beta1Zid);
      coupcombo2gsf2 = -(coupHpm1charneutL*alphasf2char*coupHpm1ffpd*sf2alpha1Zid - rj*coupHpm1charneutR*betasf2char*coupHpm1ffpu*sf2beta1Zid)*fabs(mneutralinoi)*mf;
      coupcombo3gsf2 = -(-rj*coupHpm1charneutR*alphasf2char*coupHpm1ffpu*sf2alpha1Zid + coupHpm1charneutL*betasf2char*coupHpm1ffpd*sf2beta1Zid)*mfp*fabs(mcharginoj);
      coupcombo4gsf2 = (coupHpm1charneutL*alphasf2char*coupHpm1ffpu*sf2alpha1Zid + rj*coupHpm1charneutR*betasf2char*coupHpm1ffpd*sf2beta1Zid)*fabs(mneutralinoi)*mfp;
      coupcombo5gsf2 = -(rj*coupHpm1charneutR*alphasf2char*coupHpm1ffpd*sf2alpha1Zid + coupHpm1charneutL*betasf2char*coupHpm1ffpu*sf2beta1Zid)*mf*fabs(mcharginoj);
      coupcombo6gsf2 = (coupHpm1charneutR*betasf2char*coupHpm1ffpu*sf2alpha1Zid + coupHpm1charneutL*alphasf2char*coupHpm1ffpd*sf2beta1Zid)*fabs(mcharginoj)*fabs(mneutralinoi);
      coupcombo7gsf2 = -rj*(coupHpm1charneutL*betasf2char*coupHpm1ffpd*sf2alpha1Zid + coupHpm1charneutR*alphasf2char*coupHpm1ffpu*sf2beta1Zid)*mf*mfp;
      coupcombo8gsf2 = -2*(rj*coupHpm1charneutR*betasf2char*coupHpm1ffpd*sf2alpha1Zid + coupHpm1charneutL*alphasf2char*coupHpm1ffpu*sf2beta1Zid)*mf*mfp*fabs(mcharginoj)*fabs(mneutralinoi);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1gsf2 = 0, int2gsf2 = 0, int3gsf2 = 0, int4gsf2 = 0, int5gsf2 = 0, int6gsf2 = 0, int7gsf2 = 0, int8gsf2 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mfp, m4 = mcharginoj, m5 = mWboson, m6 = msf2;
    if (norc == 'n') {
      int1gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else if (norc == 'c') {
      int1gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8gsf2 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    Gammagsf2 = coupcombo1gsf2*int1gsf2 + coupcombo2gsf2*int2gsf2 + coupcombo3gsf2*int3gsf2 + coupcombo4gsf2*int4gsf2 + coupcombo5gsf2*int5gsf2 + coupcombo6gsf2*int6gsf2 + coupcombo7gsf2*int7gsf2 + coupcombo8gsf2*int8gsf2;

    //Hpm - sf1 interference
    double coupcombo1Hpmsf1 = 0, coupcombo2Hpmsf1 = 0, coupcombo3Hpmsf1 = 0, coupcombo4Hpmsf1 = 0, coupcombo5Hpmsf1 = 0, coupcombo6Hpmsf1 = 0, coupcombo7Hpmsf1 = 0, coupcombo8Hpmsf1 = 0; 
    if (norc == 'n') {
      coupcombo1Hpmsf1 = 0.5*(coupHpm2charneutR*sf1beta1Zid*coupHpm2ffpd*alphasf1char + coupHpm2charneutL*sf1alpha1Zid*coupHpm2ffpu*betasf1char);
      coupcombo2Hpmsf1 = -(coupHpm2charneutR*sf1alpha1Zid*coupHpm2ffpu*alphasf1char + coupHpm2charneutL*sf1beta1Zid*coupHpm2ffpd*betasf1char)*fabs(mneutralinoi)*mfp;
      coupcombo3Hpmsf1 = (-coupHpm2charneutL*sf1alpha1Zid*coupHpm2ffpd*alphasf1char - coupHpm2charneutR*sf1beta1Zid*coupHpm2ffpu*betasf1char)*mf*fabs(mcharginoj)*-ri;
      coupcombo4Hpmsf1 = (-coupHpm2charneutR*sf1alpha1Zid*coupHpm2ffpd*alphasf1char - coupHpm2charneutL*sf1beta1Zid*coupHpm2ffpu*betasf1char)*fabs(mneutralinoi)*-mf;
      coupcombo5Hpmsf1 = (coupHpm2charneutL*-sf1alpha1Zid*coupHpm2ffpu*-alphasf1char - ri*coupHpm2charneutR*-sf1beta1Zid*coupHpm2ffpd*-betasf1char)*mfp*fabs(mcharginoj)*-ri;
      coupcombo6Hpmsf1 = -(-ri*coupHpm2charneutL*-sf1beta1Zid*coupHpm2ffpd*-alphasf1char + coupHpm2charneutR*-sf1alpha1Zid*coupHpm2ffpu*-betasf1char)*fabs(mneutralinoi)*fabs(mcharginoj);
      coupcombo7Hpmsf1 = (-coupHpm2charneutR*-sf1beta1Zid*coupHpm2ffpu*-alphasf1char + ri*coupHpm2charneutL*-sf1alpha1Zid*coupHpm2ffpd*-betasf1char)*mfp*mf;
      coupcombo8Hpmsf1 = 2*(-ri*coupHpm2charneutL*-sf1beta1Zid*coupHpm2ffpu*-alphasf1char + coupHpm2charneutR*-sf1alpha1Zid*coupHpm2ffpd*-betasf1char)*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj);
    }
    else if (norc == 'c') {
      coupcombo1Hpmsf1 = 0.5*(rj*coupHpm2charneutL*betasf1char*coupHpm2ffpu*sf1alpha1Zid + coupHpm2charneutR*alphasf1char*coupHpm2ffpd*sf1beta1Zid);
      coupcombo2Hpmsf1 = (rj*coupHpm2charneutL*alphasf1char*coupHpm2ffpd*sf1alpha1Zid + coupHpm2charneutR*betasf1char*coupHpm2ffpu*sf1beta1Zid)*fabs(mneutralinoi)*mf;
      coupcombo3Hpmsf1 = -(coupHpm2charneutR*alphasf1char*coupHpm2ffpu*sf1alpha1Zid + rj*coupHpm2charneutL*betasf1char*coupHpm2ffpd*sf1beta1Zid)*mfp*fabs(mcharginoj);
      coupcombo4Hpmsf1 = -(rj*coupHpm2charneutL*alphasf1char*coupHpm2ffpu*sf1alpha1Zid + coupHpm2charneutR*betasf1char*coupHpm2ffpd*sf1beta1Zid)*fabs(mneutralinoi)*mfp;
      coupcombo5Hpmsf1 = (coupHpm2charneutR*alphasf1char*coupHpm2ffpd*sf1alpha1Zid + rj*coupHpm2charneutL*betasf1char*coupHpm2ffpu*sf1beta1Zid)*mf*fabs(mcharginoj);
      coupcombo6Hpmsf1 = (coupHpm2charneutR*betasf1char*coupHpm2ffpu*sf1alpha1Zid + rj*coupHpm2charneutL*alphasf1char*coupHpm2ffpd*sf1beta1Zid)*fabs(mcharginoj)*fabs(mneutralinoi);
      coupcombo7Hpmsf1 = -(rj*coupHpm2charneutL*betasf1char*coupHpm2ffpd*sf1alpha1Zid + coupHpm2charneutR*alphasf1char*coupHpm2ffpu*sf1beta1Zid)*mf*mfp;
      coupcombo8Hpmsf1 = -2*(coupHpm2charneutR*betasf1char*coupHpm2ffpd*sf1alpha1Zid + rj*coupHpm2charneutL*alphasf1char*coupHpm2ffpu*sf1beta1Zid)*mf*mfp*fabs(mcharginoj)*fabs(mneutralinoi);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1Hpmsf1 = 0, int2Hpmsf1 = 0, int3Hpmsf1 = 0, int4Hpmsf1 = 0, int5Hpmsf1 = 0, int6Hpmsf1 = 0, int7Hpmsf1 = 0, int8Hpmsf1 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mfp, m4 = mcharginoj, m5 = mHP, m6 = msf1;
    if( norc == 'n') {
      int1Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else if (norc == 'c') {
      int1Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8Hpmsf1 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }    
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    GammaHpmsf1 = coupcombo1Hpmsf1*int1Hpmsf1 + coupcombo2Hpmsf1*int2Hpmsf1 + coupcombo3Hpmsf1*int3Hpmsf1 + coupcombo4Hpmsf1*int4Hpmsf1 + coupcombo5Hpmsf1*int5Hpmsf1 + coupcombo6Hpmsf1*int6Hpmsf1 + coupcombo7Hpmsf1*int7Hpmsf1 + coupcombo8Hpmsf1*int8Hpmsf1;

    ///Hpm - sf2 interference
    double coupcombo1Hpmsf2 = 0, coupcombo2Hpmsf2 = 0, coupcombo3Hpmsf2 = 0, coupcombo4Hpmsf2 = 0, coupcombo5Hpmsf2 = 0, coupcombo6Hpmsf2 = 0, coupcombo7Hpmsf2 = 0, coupcombo8Hpmsf2 = 0; 
    if (norc == 'n') {
      coupcombo1Hpmsf2 = -0.5*(coupHpm2charneutR*-sf2beta1Zid*coupHpm2ffpd*-alphasf2char + coupHpm2charneutL*-sf2alpha1Zid*coupHpm2ffpu*-betasf2char);
      coupcombo2Hpmsf2 = (coupHpm2charneutR*-sf2alpha1Zid*coupHpm2ffpu*-alphasf2char + coupHpm2charneutL*-sf2beta1Zid*coupHpm2ffpd*-betasf2char)*fabs(mneutralinoi)*mfp;
      coupcombo3Hpmsf2 = (coupHpm2charneutL*sf2alpha1Zid*coupHpm2ffpd*-alphasf2char - coupHpm2charneutR*-sf2beta1Zid*coupHpm2ffpu*-betasf2char)*mf*fabs(mcharginoj)*ri;
      coupcombo4Hpmsf2 = (coupHpm2charneutR*sf2alpha1Zid*coupHpm2ffpd*-alphasf2char + coupHpm2charneutL*-sf2beta1Zid*coupHpm2ffpu*betasf2char)*fabs(mneutralinoi)*mf;
      coupcombo5Hpmsf2 = (ri*coupHpm2charneutL*-sf2alpha1Zid*coupHpm2ffpu*-alphasf2char + coupHpm2charneutR*-sf2beta1Zid*coupHpm2ffpd*-betasf2char)*mfp*fabs(mcharginoj);
      coupcombo6Hpmsf2 = -(coupHpm2charneutL*-sf2beta1Zid*coupHpm2ffpd*-alphasf2char + coupHpm2charneutR*-sf2alpha1Zid*coupHpm2ffpu*-betasf2char)*fabs(mneutralinoi)*fabs(mcharginoj)*rj*ri;
      coupcombo7Hpmsf2 = (coupHpm2charneutR*-sf2beta1Zid*coupHpm2ffpu*-alphasf2char + coupHpm2charneutL*-sf2alpha1Zid*coupHpm2ffpd*-betasf2char)*mfp*mf;
      coupcombo8Hpmsf2 = 2*(ri*coupHpm2charneutL*-sf2beta1Zid*coupHpm2ffpu*-alphasf2char + coupHpm2charneutR*-sf2alpha1Zid*coupHpm2ffpd*-betasf2char)*fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj)*rj;
    }
    else if (norc == 'c') {
      coupcombo1Hpmsf2 = 0.5*(coupHpm2charneutL*betasf2char*coupHpm2ffpu*sf2alpha1Zid + rj*coupHpm2charneutR*alphasf2char*coupHpm2ffpd*sf2beta1Zid);
      coupcombo2Hpmsf2 = -(coupHpm2charneutL*alphasf2char*coupHpm2ffpd*sf2alpha1Zid - rj*coupHpm2charneutR*betasf2char*coupHpm2ffpu*sf2beta1Zid)*fabs(mneutralinoi)*mf;
      coupcombo3Hpmsf2 = -(-rj*coupHpm2charneutR*alphasf2char*coupHpm2ffpu*sf2alpha1Zid + coupHpm2charneutL*betasf2char*coupHpm2ffpd*sf2beta1Zid)*mfp*fabs(mcharginoj);
      coupcombo4Hpmsf2 = (coupHpm2charneutL*alphasf2char*coupHpm2ffpu*sf2alpha1Zid - rj*coupHpm2charneutR*betasf2char*coupHpm2ffpd*sf2beta1Zid)*fabs(mneutralinoi)*mfp;
      coupcombo5Hpmsf2 = -(rj*coupHpm2charneutR*alphasf2char*coupHpm2ffpd*sf2alpha1Zid + coupHpm2charneutL*betasf2char*coupHpm2ffpu*sf2beta1Zid)*mf*fabs(mcharginoj);
      coupcombo6Hpmsf2 = (coupHpm2charneutR*betasf2char*coupHpm2ffpu*sf2alpha1Zid + coupHpm2charneutL*alphasf2char*coupHpm2ffpd*sf2beta1Zid)*fabs(mcharginoj)*fabs(mneutralinoi);
      coupcombo7Hpmsf2 = -rj*(coupHpm2charneutL*betasf2char*coupHpm2ffpd*sf2alpha1Zid + coupHpm2charneutR*alphasf2char*coupHpm2ffpu*sf2beta1Zid)*mf*mfp;
      coupcombo8Hpmsf2 = -2*(-rj*coupHpm2charneutR*betasf2char*coupHpm2ffpd*sf2alpha1Zid + coupHpm2charneutL*alphasf2char*coupHpm2ffpu*sf2beta1Zid)*mf*mfp*fabs(mcharginoj)*fabs(mneutralinoi);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1Hpmsf2 = 0, int2Hpmsf2 = 0, int3Hpmsf2 = 0, int4Hpmsf2 = 0, int5Hpmsf2 = 0, int6Hpmsf2 = 0, int7Hpmsf2 = 0, int8Hpmsf2 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mfp, m4 = mcharginoj, m5 = mHP, m6 = msf2;
    if(norc == 'n') {
      int1Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }
    else if( norc == 'c') {
      int1Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp1gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int4Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp2gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int5Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp3gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int2Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp4gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int3Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp5gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int6Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp6gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int7Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp7gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
      int8Hpmsf2 = 2*fabs(m1)*dgauss(gneuticharjffp8gsfpdgauss,fabs(mcharginoj),Eupper3,accuracy);
    }  
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    GammaHpmsf2 = coupcombo1Hpmsf2*int1Hpmsf2 + coupcombo2Hpmsf2*int2Hpmsf2 + coupcombo3Hpmsf2*int3Hpmsf2 + coupcombo4Hpmsf2*int4Hpmsf2 + coupcombo5Hpmsf2*int5Hpmsf2 + coupcombo6Hpmsf2*int6Hpmsf2 + coupcombo7Hpmsf2*int7Hpmsf2 + coupcombo8Hpmsf2*int8Hpmsf2;

    ///Sfp Sfp interference
    double coupcombo1sfpsfp = 0, coupcombo2sfpsfp = 0, coupcombo3sfpsfp = 0, coupcombo4sfpsfp = 0;
    if (norc == 'n') {
      coupcombo1sfpsfp = (sf1beta1Ziu*sf2beta1Ziu + sf1alpha1Ziu*sf2alpha1Ziu)*ri;
      coupcombo2sfpsfp = (sf1alpha1Ziu*sf2beta1Ziu + sf1beta1Ziu*sf2alpha1Ziu)*ri;
      coupcombo3sfpsfp = (-alphasfp1char*alphasfp2char + betasfp1char*betasfp2char)*ri;
      coupcombo4sfpsfp = (betasfp1char*alphasfp2char - alphasfp1char*betasfp2char)*rj*ri;
    }
    else if (norc == 'c') {
      coupcombo1sfpsfp = -betasfp1char*betasfp2char + alphasfp1char*alphasfp2char*rj;
      coupcombo2sfpsfp = rj*alphasfp1char*betasfp2char - alphasfp2char*betasfp1char;
      coupcombo3sfpsfp = (-sf1alpha1Ziu*sf2alpha1Ziu + rj*sf1beta1Ziu*sf2beta1Ziu);
      coupcombo4sfpsfp = (-sf1alpha1Ziu*sf2beta1Ziu + rj*sf2alpha1Ziu*sf1beta1Ziu);
    }
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double int1sfpsfp = 0, int2sfpsfp = 0, int3sfpsfp = 0, int4sfpsfp = 0;
 
    m1 = mneutralinoi, m2 = mfp, m3 = mf, m4 = mcharginoj, m5 = msfp1, m6 = msfp2;
    if (norc == 'n') {
      int1sfpsfp = 2*fabs(m1)*dgauss(gneuticharjffp1sfpsfpdgauss,mfp,Eupper,accuracy);
      int2sfpsfp = 2*fabs(m1)*dgauss(gneuticharjffp2sfpsfpdgauss,mfp,Eupper,accuracy);
      int3sfpsfp = 2*fabs(m1)*dgauss(gneuticharjffp3sfpsfpdgauss,mfp,Eupper,accuracy);
      int4sfpsfp = 2*fabs(m1)*dgauss(gneuticharjffp4sfpsfpdgauss,mfp,Eupper,accuracy);
    }
    else if (norc == 'c') {
      int1sfpsfp = 2*fabs(m1)*dgauss(gneuticharjffp1sfpsfpdgauss,mfp,Eupper,accuracy);
      int3sfpsfp = 2*fabs(m1)*dgauss(gneuticharjffp2sfpsfpdgauss,mfp,Eupper,accuracy);
      int2sfpsfp = 2*fabs(m1)*dgauss(gneuticharjffp3sfpsfpdgauss,mfp,Eupper,accuracy);
      int4sfpsfp = 2*fabs(m1)*dgauss(gneuticharjffp4sfpsfpdgauss,mfp,Eupper,accuracy);
    }    
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    if (norc == 'n') {
      Gammasfpsfp = 4*coupcombo2sfpsfp*coupcombo4sfpsfp*-fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj)*int1sfpsfp + 2*coupcombo2sfpsfp*coupcombo3sfpsfp*fabs(mneutralinoi)*mfp*int2sfpsfp + 2*coupcombo1sfpsfp*coupcombo4sfpsfp*-mf*fabs(mcharginoj)*int3sfpsfp + coupcombo1sfpsfp*coupcombo3sfpsfp*int4sfpsfp;
    }
    else if (norc == 'c') {
      Gammasfpsfp = 4*coupcombo2sfpsfp*coupcombo4sfpsfp*-fabs(mneutralinoi)*mfp*mf*fabs(mcharginoj)*int1sfpsfp + 2*coupcombo2sfpsfp*coupcombo3sfpsfp*fabs(mneutralinoi)*mf*int2sfpsfp + 2*coupcombo1sfpsfp*coupcombo4sfpsfp*-mfp*fabs(mcharginoj)*int3sfpsfp + coupcombo1sfpsfp*coupcombo3sfpsfp*int4sfpsfp;
    } 
    else {
      throw("problem: norc must be n or c for neut or chargino respectively as decaying particle in neutralinoamplitudedecaycharginoffprimebar");
    }

    double Gammasf1sf2 = 0;
    ///Sf - Sf interference
    double coupcombo1sf1sf2 = 0, coupcombo2sf1sf2 = 0, coupcombo3sf1sf2 = 0, coupcombo4sf1sf2 = 0;
    coupcombo3sf1sf2 = -(rc*sf1alpha1Zid*sf2alpha1Zid + sf1beta1Zid*sf2beta1Zid)*rj;
    coupcombo4sf1sf2 = (ri*sf1beta1Zid*sf2alpha1Zid - sf1alpha1Zid*sf2beta1Zid)*rj*rc;
    coupcombo1sf1sf2 = -(betasf1char*-betasf2char + alphasf1char*alphasf2char)*ri;
    coupcombo2sf1sf2 = rc*(-alphasf1char*-betasf2char - betasf1char*alphasf2char);

    double int1sf1sf2 = 0, int2sf1sf2 = 0, int3sf1sf2 = 0, int4sf1sf2 = 0;
    m1 = mneutralinoi, m2 = mf, m3 = mcharginoj, m4 = mfp, m5 = msf1, m6 = msf2;
    int1sf1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp1sf1sf2dgauss,mf,Eupper2,accuracy);
    int2sf1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp2sf1sf2dgauss,mf,Eupper2,accuracy);
    int3sf1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp3sf1sf2dgauss,mf,Eupper2,accuracy);
    int4sf1sf2 = 2*fabs(m1)*dgauss(gneuticharjffp4sf1sf2dgauss,mf,Eupper2,accuracy);
    
      Gammasf1sf2 = (-rc*ri*coupcombo1sf1sf2*coupcombo3sf1sf2*int4sf1sf2 + 2*ri*coupcombo1sf1sf2*coupcombo4sf1sf2*mf*-fabs(mcharginoj)*int2sf1sf2 -rc*2*coupcombo2sf1sf2*coupcombo3sf1sf2*fabs(mneutralinoi)*mfp*int3sf1sf2 + 4*coupcombo2sf1sf2*coupcombo4sf1sf2*fabs(mneutralinoi)*mf*-fabs(mcharginoj)*mfp*int1sf1sf2);
   

    amplitudeW = Nc/(512*pow(PI,3)*pow(fabs(mneutralinoi),3))*(GammaW + Gammasf1 + Gammasf2 + Gammasfp1 + Gammasfp2 -2*Gammasfp1sf1 - 2*Gammasfp1sf2 - 2*Gammasfp2sf1 - 2*Gammasfp2sf2 + 2*GammaWHpm + 2*GammaWgoldstone+ GammaHpm + Gammagoldstone -2*GammaWSfp1 - 2*GammaWSfp2 - 2*GammaWSf1 - 2*GammaWSf2 + 2*GammaHgoldstone - 2*Gammagsfp1 - 2*Gammagsfp2 - 2*Gammagsf1 - 2*Gammagsf2 - 2*GammaHpmsfp1 - 2*GammaHpmsfp2 - 2*GammaHpmsf1 - 2*GammaHpmsf2 + 2*Gammasfpsfp - 2*Gammasf1sf2);
    // fout << "amplitudeW = " << amplitudeW << endl;
    // fout << "GammaW = " << GammaW << " Gammasf1 = " << Gammasf1 << " Gammasf2 = " << Gammasf2 << " Gammasfp1 = " << Gammasfp1 << " Gammasfp2 = " << Gammasfp2 << endl;
    // fout << "Gammasfp1sf1 = " << Gammasfp1sf1 << " Gammasfp1sf2 = " << Gammasfp1sf2 <<" Gammasfp2sf1 = " << Gammasfp2sf1 <<" Gammasfp2sf2 = " << Gammasfp2sf2 << endl;
    // fout << "GammaWHpm = " << GammaWHpm << " GammaWgoldstone = " << GammaWgoldstone << " GammaHpm = " << GammaHpm << " Gammagoldstone = " << Gammagoldstone << endl;
    // fout << "GammaWSfp1 = " << GammaWSfp1 << " GammaWSfp2 = " << GammaWSfp2 << " GammaWSf1 = " << GammaWSf1 << " GammaWSf2 = " << GammaWSf2 << endl;
    // fout << "GammaHgoldstone = " << GammaHgoldstone << " Gammagsfp1 = " << Gammagsfp1 << " Gammagsfp2 = " << Gammagsfp2 << " Gammagsf1 = " << Gammagsf1 << " Gammagsf2 = " << Gammagsf2 << endl;
    // fout << "GammaHpmsfp1 = " << GammaHpmsfp1 << " GammaHpmsfp2 = " << GammaHpmsfp2 << " GammaHpmsf1 = " << GammaHpmsf1 << " GammaHpmsf2 = " << GammaHpmsf2 << endl;
    // fout << "Gammasfpsfp = " << Gammasfpsfp << " Gammasf1sf2 = " << Gammasf1sf2 << endl;
  }

  return amplitudeW;
}
    



double higgsAamplitudedecaygammagammaNMSSM (double m1, double g, double gprime, double alpha, double mWboson, DoubleMatrix & CPOMix, double beta, double mtop, double mbottom, double mcharm, double mtau, double mch1, double mch2, double thetaL, double thetaR, double lam, int higgs)
{
  double prefactor=0, Itr=0, Iti=0, Ibr = 0, Ibi = 0, Icr = 0, Ici = 0, Itaur = 0, Itaui = 0, Ichar1r = 0, Ichar1i = 0, Ichar2r = 0, Ichar2i = 0, couplingt = 0, couplingb = 0, couplingc = 0, couplingtau = 0, couplingch1 = 0, couplingch2 = 0, kintr = 0, kinti = 0, kinbr =0, kinbi = 0, kincr = 0, kinci = 0, kintaur = 0, kintaui = 0, kinch1r = 0, kinch1i = 0, kinch2r = 0, kinch2i = 0, matelemmodsquare=0, amplitudeW=0;

  DoubleVector tfoftau(3), bfoftau(3), cfoftau(3), taufoftau(3), ch1foftau(3), ch2foftau(3);
  ///Initialise these components
  for (int i = 1; i <= 3; i++) {
    tfoftau(i) = 0; bfoftau(i) = 0; cfoftau(i) = 0; taufoftau(i) = 0; ch1foftau(i) = 0; ch2foftau(i) = 0;
  }
  tfoftau = foftau(mtop, m1); bfoftau = foftau(mbottom, m1); cfoftau = foftau(mcharm, m1); taufoftau = foftau(mtau, m1); ch1foftau = foftau(mch1, m1); ch2foftau = foftau(mch2, m1);

  kintr = tfoftau(3)*(tfoftau(1)); kinbr = bfoftau(3)*(bfoftau(1)); kincr = cfoftau(3)*(cfoftau(1)); kintaur = taufoftau(3)*(taufoftau(1));
  kinti = tfoftau(3)*(tfoftau(2)); kinbi = bfoftau(3)*(bfoftau(2)); kinci = cfoftau(3)*(cfoftau(2)); kintaui = taufoftau(3)*(taufoftau(2));
  
  kinch1r = ch1foftau(3)*(ch1foftau(1)); kinch2r = ch2foftau(3)*(ch2foftau(1));
  kinch1i = ch1foftau(3)*(ch1foftau(2)); kinch2i = ch2foftau(3)*(ch2foftau(2));

  couplingt = 4./3*CPOMix(higgs,1)/(sin(beta)); couplingc = 4./3*CPOMix(higgs,1)/(sin(beta)); couplingb = 1./3*CPOMix(higgs,2)/(cos(beta)); couplingtau = CPOMix(higgs,2)/(cos(beta));
  
  couplingch1 = 2*mWboson/(g*mch1)*(lam/(root2)*CPOMix(higgs,3)*cos(thetaL)*cos(thetaR) - g/(root2)*(CPOMix(higgs,1)*sin(thetaL)*cos(thetaR) + CPOMix(higgs,2)*cos(thetaL)*sin(thetaR)));
  couplingch2 = 2*mWboson/(g*mch2)*(lam/(root2)*CPOMix(higgs,3)*sin(thetaL)*sin(thetaR) + g/(root2)*(CPOMix(higgs,1)*cos(thetaL)*sin(thetaR) + CPOMix(higgs,2)*sin(thetaL)*cos(thetaR)));

  Itr = couplingt*kintr; Icr = couplingc*kincr; Ibr = couplingb*kinbr; Itaur = couplingtau*kintaur; Ichar1r = couplingch1*kinch1r; Ichar2r = couplingch2*kinch2r;
  Iti = couplingt*kinti; Ici = couplingc*kinci; Ibi = couplingb*kinbi; Itaui = couplingtau*kintaui; Ichar1i = couplingch1*kinch1i; Ichar2i = couplingch2*kinch2i;

  DoubleVector matelemsum(2);
  for (int i = 1; i <= 2; i++) {
    matelemsum(i) = 0;
  }
  matelemsum(1) = Itr + Ibr + Icr + Itaur + Ichar1r + Ichar2r;
  matelemsum(2) = Iti + Ibi + Ici + Itaui + Ichar1i + Ichar2i;

  prefactor = (GMU*pow(alpha,2))/(root2*32*pow(PI,3))*pow(m1,3);

  matelemmodsquare = pow(matelemsum(1),2) + pow(matelemsum(2),2);
  amplitudeW = prefactor*matelemmodsquare;
     
  return amplitudeW;

}


double higgsAamplitudedecaygluongluonNMSSM (double m1, double g, double gs, double alphas, double mWboson, DoubleMatrix & CPOMix, double beta, double mtop, double mbottom, double mcharm, double lam, int higgs, bool QCD)
{
  double prefactor=0, Itr=0, Iti=0, Ibr = 0, Ibi = 0, Icr = 0, Ici = 0, couplingt = 0, couplingb = 0, couplingc = 0, kintr = 0, kinti = 0, kinbr =0, kinbi = 0, kincr = 0, kinci = 0, matelemmodsquare=0, amplitudeW=0;


  DoubleVector tfoftau(3), bfoftau(3), cfoftau(3);
  ///Initialise these components
  for (int i = 1; i <= 3; i++) {
    tfoftau(i) = 0; bfoftau(i) = 0; cfoftau(i) = 0;
  }

  tfoftau = foftau(mtop, m1); bfoftau = foftau(mbottom, m1); cfoftau = foftau(mcharm, m1);

  kintr = tfoftau(3)*(tfoftau(1)); kinbr = bfoftau(3)*(bfoftau(1)); kincr = cfoftau(3)*(cfoftau(1)); 
  kinti = tfoftau(3)*(tfoftau(2)); kinbi = bfoftau(3)*(bfoftau(2)); kinci = cfoftau(3)*(cfoftau(2));

  couplingt = CPOMix(higgs,1)/(sin(beta)); couplingc = CPOMix(higgs,1)/(sin(beta)); couplingb = CPOMix(higgs,2)/(cos(beta));
  
  Itr = couplingt*kintr; Icr = couplingc*kincr; Ibr = couplingb*kinbr;
  Iti = couplingt*kinti; Ici = couplingc*kinci; Ibi = couplingb*kinbi;

  DoubleVector matelemsum(2);
  for (int i = 1; i <= 2; i++) {
    matelemsum(i) = 0;
  }
  matelemsum(1) = Itr + Ibr + Icr;
  matelemsum(2) = Iti + Ibi + Ici;

  prefactor = GMU*pow(alphas,2)/(root2*16*pow(PI,3))*pow(m1,3);

  matelemmodsquare = pow(matelemsum(1),2) + pow(matelemsum(2),2);
  amplitudeW = prefactor*matelemmodsquare;
  
  double SMTOTRE = 0, SMTOTIM = 0, SQTOTRE = 0, SQTOTIM = 0;
  SMTOTRE = Itr + Ibr + Icr;
  SQTOTRE = 0; ///No squark loop contributions as CP odd higgs
  SMTOTIM = Iti + Ibi + Ici;
  SQTOTIM = 0; ///No squark loop contributions as CP odd higgs

  if (QCD == true) {
    int NF = 0;
    NF = 6;

    amplitudeW = hggQCDcorrections(amplitudeW, alphas, NF, 'A', prefactor, SMTOTRE, SMTOTIM, SQTOTRE, SQTOTIM)(1); ///Pass hggQCDcorrections 'A' as it just matters whether it is CPeven (so 95/4 in FQCD) or CPodd (so 97/4 in FQCD)
  }
  else if (QCD == false) {
    amplitudeW = amplitudeW;
  }
     
  return amplitudeW;
}


double higgsAamplitudedecayZgammaNMSSM (double m1, double g, double gp, double alpha, double mWboson, double mZboson, DoubleMatrix & CPOMix, double beta, double mtop, double mbottom, double mcharm, double mch1, double mch2, double thetaL, double thetaR, double lam, int higgs)
{
  double prefactor=0, Itr=0, Iti=0, Ibr = 0, Ibi = 0, Icr = 0, Ici = 0, Ichar1r = 0, Ichar1i = 0, Ichar2r = 0, Ichar2i = 0, couplingt = 0, couplingb = 0, couplingc = 0, couplingch1 =0, couplingch2 = 0, kintr = 0, kinti = 0, kinbr =0, kinbi = 0, kincr = 0, kinci = 0, kinch1r = 0, kinch1i = 0, kinch2r = 0, kinch2i = 0, matelemmodsquare=0, amplitudeW=0;
  double sin2thW = 0, sinthW = 0, costhW = 0;

  DoubleVector tfoftau(3), bfoftau(3), cfoftau(3), ch1foftau(3), ch2foftau(3), tZfoftau(3), bZfoftau(3), cZfoftau(3), ch1Zfoftau(3), ch2Zfoftau(3);
  ///Initialise these components
  for (int i = 1; i <= 3; i++) {
    tfoftau(i) = 0; bfoftau(i) = 0; cfoftau(i) = 0; ch1foftau(i) = 0; ch2foftau(i) = 0;
    tZfoftau(i) = 0; bZfoftau(i) = 0; cZfoftau(i) = 0; ch1Zfoftau(i) = 0; ch2Zfoftau(i) = 0;
  }
  tfoftau = foftau(mtop, m1); bfoftau = foftau(mbottom, m1); cfoftau = foftau(mcharm, m1); ch1foftau = foftau(mch1, m1); ch2foftau = foftau(mch2, m1);
  tZfoftau = foftau(mtop, mZboson); bZfoftau = foftau(mbottom, mZboson); cZfoftau = foftau(mcharm, mZboson); ch1Zfoftau = foftau(mch1, mZboson); ch2Zfoftau = foftau(mch2, mZboson);

  sin2thW = pow(gp,2)/(pow(g,2)+pow(gp,2)); sinthW = pow(sin2thW,0.5); costhW = g/(pow(pow(g,2)+pow(gp,2),0.5));

  kintr = tfoftau(3)*tZfoftau(3)/(2*(tfoftau(3)-tZfoftau(3)))*(tfoftau(1)-tZfoftau(1)); kinbr = bfoftau(3)*bZfoftau(3)/(2*(bfoftau(3)-bZfoftau(3)))*(bfoftau(1)-bZfoftau(1)); kincr = cfoftau(3)*cZfoftau(3)/(2*(cfoftau(3)-cZfoftau(3)))*(cfoftau(1)-cZfoftau(1)); kinch1r = ch1foftau(3)*ch1Zfoftau(3)/(2*(ch1foftau(3)-ch1Zfoftau(3)))*(ch1foftau(1)-ch1Zfoftau(1)); kinch2r = ch2foftau(3)*ch2Zfoftau(3)/(2*(ch2foftau(3)-ch2Zfoftau(3)))*(ch2foftau(1)-ch2Zfoftau(1));
  kinti = tfoftau(3)*tZfoftau(3)/(2*(tfoftau(3)-tZfoftau(3)))*(tfoftau(2)-tZfoftau(2)); kinbi = bfoftau(3)*bZfoftau(3)/(2*(bfoftau(3)-bZfoftau(3)))*(bfoftau(2)-bZfoftau(2)); kinci = cfoftau(3)*cZfoftau(3)/(2*(cfoftau(3)-cZfoftau(3)))*(cfoftau(2)-cZfoftau(2)); kinch1i = ch1foftau(3)*ch1Zfoftau(3)/(2*(ch1foftau(3)-ch1Zfoftau(3)))*(ch1foftau(2)-ch1Zfoftau(2)); kinch2i = ch2foftau(3)*ch2Zfoftau(3)/(2*(ch2foftau(3)-ch2Zfoftau(3)))*(ch2foftau(2)-ch2Zfoftau(2));
 
  couplingt = -2*(1-8*sin2thW/3)/(sinthW*costhW)*CPOMix(higgs,1)/(sin(beta)); couplingc = -2*(1-8*sin2thW/3)/(sinthW*costhW)*CPOMix(higgs,1)/(sin(beta)); couplingb =(-1+4*sin2thW/3)/(sinthW*costhW)*CPOMix(higgs,2)/(cos(beta));

  
  couplingch1 = 4*mWboson/(mch1*g*sinthW*costhW)*(-pow(sin(thetaR),2) - 0.5*pow(cos(thetaR),2) - pow(sin(thetaL),2) - 0.5*pow(cos(thetaL),2) + 2*sin2thW)*(lam/(root2)*CPOMix(higgs,3)*cos(thetaL)*cos(thetaR) - g/(root2)*(CPOMix(higgs,1)*sin(thetaL)*cos(thetaR) + CPOMix(higgs,2)*cos(thetaL)*sin(thetaR)));

  couplingch2 = 4*mWboson/(mch2*g*sinthW*costhW)*(-pow(cos(thetaR),2) - 0.5*pow(sin(thetaR),2) - pow(cos(thetaL),2) - 0.5*pow(sin(thetaL),2) + 2*sin2thW)*(lam/(root2)*CPOMix(higgs,3)*sin(thetaL)*sin(thetaR) + g/(root2)*(CPOMix(higgs,1)*cos(thetaL)*sin(thetaR) + CPOMix(higgs,2)*sin(thetaL)*cos(thetaR)));

  
  Itr = couplingt*kintr; Icr = couplingc*kincr; Ibr = couplingb*kinbr; Ichar1r = couplingch1*kinch1r; Ichar2r = couplingch2*kinch2r; 
  Iti = couplingt*kinti; Ici = couplingc*kinci; Ibi = couplingb*kinbi; Ichar1i = couplingch1*kinch1i; Ichar2i = couplingch2*kinch2i;

  DoubleVector matelemsum(2);
  for (int i = 1; i <= 2; i++) {
    matelemsum(i) = 0;
  }
  matelemsum(1) = Itr + Ibr + Icr + Ichar1r + Ichar2r;
  matelemsum(2) = Iti + Ibi + Ici + Ichar1i + Ichar2i;

  // prefactor = pow(g,2)*pow(m1,3)*pow(alpha,2)/(512*pow(PI,3)*pow(mWboson,2))*pow((1-pow(mZboson/m1,2)),3);
  prefactor = GMU*pow(m1,3)*pow(alpha,2)/(root2*64*pow(PI,3))*pow((1-pow(mZboson/m1,2)),3);

  matelemmodsquare = pow(matelemsum(1),2) + pow(matelemsum(2),2);
  amplitudeW = prefactor*matelemmodsquare;
     
  return amplitudeW;
}


double higgsCPevenamplitudedecaygammagammaNMSSM(double m1, double mtop, double mbottom, double mcharm, double mtau, double mWboson, double mHpm, double mchar1, double mchar2, double mscharmL, double mscharmR, double mstop1, double mstop2, double msstrangeL, double msstrangeR, double msbottom1, double msbottom2, double msmuonL, double msmuonR, double mstau1, double mstau2, DoubleMatrix & CPEMix, double beta, double g, double gp, double alpha, double thetat, double thetab, double thetatau, double thetaL, double thetaR, double At, double Ab, double Atau, double mu, double mueff, double lam, double kappa, double Alambda, int higgs)
{
  double amplitudeW = 0, prefactor = 0, couplingt = 0, couplingb = 0, couplingc = 0, couplingtau = 0, couplingW = 0, couplingHpm = 0, couplingch1 = 0, couplingch2 = 0, couplingscL = 0, couplingscR = 0, couplingst1 = 0, couplingst2 = 0, couplingssL = 0, couplingssR = 0, couplingsb1 = 0, couplingsb2 = 0, couplingsmuL = 0, couplingsmuR = 0, couplingstau1 = 0, couplingstau2 = 0;
  double ft = 0, fb = 0, ftau = 0;
  DoubleVector tfoftau(3), bfoftau(3), cfoftau(3), taufoftau(3), Wfoftau(3), Hpmfoftau(3), ch1foftau(3), ch2foftau(3), scLfoftau(3), scRfoftau(3), ssLfoftau(3), ssRfoftau(3), st1foftau(3), st2foftau(3), sb1foftau(3), sb2foftau(3), stau1foftau(3), stau2foftau(3), smuLfoftau(3), smuRfoftau(3);
  ///Initialise
  for (int i = 1; i<=3; i++) {
    tfoftau(i) = 0; bfoftau(i) = 0; cfoftau(i) = 0, taufoftau(i) = 0, Wfoftau(i) = 0, Hpmfoftau(i) = 0, ch1foftau(i) = 0, ch2foftau(i) = 0, scLfoftau(i) = 0, scRfoftau(i) = 0, ssLfoftau(i) = 0, ssRfoftau(i) = 0, st1foftau(i) = 0, st2foftau(i) = 0, sb1foftau(i) = 0, sb2foftau(i) = 0, stau1foftau(i) = 0, stau2foftau(i) = 0, smuLfoftau(i) = 0, smuRfoftau(i) = 0;
  }
  double kintr = 0, kinbr = 0, kincr = 0, kintaur = 0, kinWr = 0, kinHpmr = 0, kinch1r = 0, kinch2r = 0, kinscLr = 0, kinscRr = 0, kinssLr = 0, kinssRr = 0, kinst1r = 0, kinst2r = 0, kinsb1r = 0, kinsb2r = 0, kinsmuLr = 0, kinsmuRr = 0, kinstau1r = 0, kinstau2r = 0;
  double kinti = 0, kinbi = 0, kinci = 0, kintaui = 0, kinWi = 0, kinHpmi = 0, kinch1i = 0, kinch2i = 0, kinscLi = 0, kinscRi = 0, kinssLi = 0, kinssRi = 0, kinst1i = 0, kinst2i = 0, kinsb1i = 0, kinsb2i = 0, kinsmuLi = 0, kinsmuRi = 0, kinstau1i = 0, kinstau2i = 0;
  double Itr = 0, Icr = 0, Ibr = 0, Itaur = 0, Ichar1r = 0, Ichar2r = 0, IWr = 0, IHpmr = 0, Iti = 0, Ici = 0, Ibi = 0, Itaui = 0, Ichar1i = 0, Ichar2i = 0, IWi = 0, IHpmi = 0, IscLr = 0, IscRr = 0, IssLr = 0, IssRr = 0, IsmuLr = 0, IsmuRr = 0, IscLi = 0, IscRi = 0, IssLi = 0, IssRi = 0, IsmuLi = 0, IsmuRi = 0, Ist1r = 0, Ist2r = 0, Isb1r = 0, Isb2r = 0, Istau1r = 0, Istau2r = 0, Ist1i = 0, Ist2i = 0, Isb1i = 0, Isb2i = 0, Istau1i = 0, Istau2i = 0;

  tfoftau = foftau(mtop, m1); bfoftau = foftau(mbottom, m1); cfoftau = foftau(mcharm, m1); taufoftau = foftau(mtau, m1); Wfoftau = foftau(mWboson, m1); Hpmfoftau = foftau(mHpm, m1); ch1foftau = foftau(mchar1, m1); ch2foftau = foftau(mchar2, m1); scLfoftau = foftau(mscharmL, m1); scRfoftau = foftau(mscharmR, m1); ssLfoftau = foftau(msstrangeL, m1); ssRfoftau = foftau(msstrangeR, m1); st1foftau = foftau(mstop1, m1); st2foftau = foftau(mstop2, m1); sb1foftau = foftau(msbottom1, m1); sb2foftau = foftau(msbottom2, m1); stau1foftau = foftau(mstau1, m1); stau2foftau = foftau(mstau2, m1); smuLfoftau = foftau(msmuonL, m1); smuRfoftau = foftau(msmuonR, m1);

  couplingt = 4./3*CPEMix(higgs,1)/sin(beta); couplingc = 4./3*CPEMix(higgs,1)/sin(beta); couplingb = CPEMix(higgs,2)/(3*cos(beta)); couplingtau = CPEMix(higgs,2)/cos(beta);
  couplingW = CPEMix(higgs,1)*sin(beta) + CPEMix(higgs,2)*cos(beta);

  couplingch1 = (lam/(root2)*CPEMix(higgs,3)*cos(thetaL)*cos(thetaR) + g/(root2)*(CPEMix(higgs,1)*sin(thetaL)*cos(thetaR) + CPEMix(higgs,2)*cos(thetaL)*sin(thetaR)))*1/(pow(GFosqrt2*2,0.5)*mchar1);
  couplingch2 = (lam/(root2)*CPEMix(higgs,3)*sin(thetaL)*sin(thetaR) - g/(root2)*(CPEMix(higgs,1)*cos(thetaL)*sin(thetaR) + CPEMix(higgs,2)*sin(thetaL)*cos(thetaR)))*1/(pow(GFosqrt2*2,0.5)*mchar2);

  couplingHpm = (lam*mueff/(root2)*(2*CPEMix(higgs,3)*pow(cos(beta),2) + 2*CPEMix(higgs,3)*pow(sin(beta),2)) - pow(lam,2)*mWboson*sin(beta)/g*2*CPEMix(higgs,2)*cos(beta)*sin(beta) - pow(lam,2)*mWboson*sin(beta)/g*2*CPEMix(higgs,1)*cos(beta)*sin(beta) + mueff*kappa*root2*2*CPEMix(higgs,3)*cos(beta)*sin(beta) + lam*Alambda/root2*2*CPEMix(higgs,3)*cos(beta)*sin(beta) + pow(gp,2)/4*mWboson/g*(sin(beta)*(2*CPEMix(higgs,1)*pow(cos(beta),2) - 2*CPEMix(higgs,1)*pow(sin(beta),2)) + cos(beta)*(2*CPEMix(higgs,2)*pow(sin(beta),2) - 2*CPEMix(higgs,2)*pow(cos(beta),2))) + g/4*mWboson*(sin(beta)*(2*CPEMix(higgs,1)*pow(cos(beta),2) + 2*CPEMix(higgs,1)*pow(sin(beta),2) + 2*2*CPEMix(higgs,2)*sin(beta)*cos(beta)) + cos(beta)*(2*CPEMix(higgs,2)*pow(cos(beta),2) + 2*CPEMix(higgs,2)*pow(sin(beta),2) + 2*2*CPEMix(higgs,1)*sin(beta)*cos(beta))) + lam/root2*0)/(2*pow(mHpm,2)*pow(2*GFosqrt2,0.5));

  couplingscL = 4./3*2*mWboson/(g*pow(mscharmL,2))*(pow(gp,2)/12+pow(g,2)/4)*2*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));
  couplingscR = 4./3*2*mWboson/(g*pow(mscharmR,2))*(pow(gp,2)/6)*2*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));
  couplingssL = 1./3*2*mWboson/(g*pow(msstrangeL,2))*2*(pow(gp,2)/12 + pow(g,2)/4)*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));
  couplingssR = 1./3*2*mWboson/(g*pow(msstrangeR,2))*(pow(gp,2)/6)*2*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));
  couplingsmuL = 2*mWboson/(g*pow(msmuonL,2))*2*mWboson/g*(-pow(gp,2)/4+pow(g,2)/4)*(sin(beta)*CPEMix(higgs,1) - CPEMix(higgs,2));
  couplingsmuR = 2*mWboson/(g*pow(msmuonR,2))*2*mWboson/g*(pow(gp,2)/2)*(sin(beta)*CPEMix(higgs,1) - CPEMix(higgs,2));

  ft = g*mtop/(root2*mWboson*sin(beta)); fb = g*mbottom/(root2*mWboson*cos(beta)); ftau = g*mtau/(root2*mWboson*cos(beta));

  couplingst1 = 1/(2*pow(2*GFosqrt2,0.5)*pow(mstop1,2))*(pow(cos(thetat),2)*root2*(pow(ft,2)*root2*mWboson*sin(beta)/g*CPEMix(higgs,1) + (pow(gp,2)/12 - pow(g,2)/4)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + pow(sin(thetat),2)*root2*(pow(ft,2)*root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - (pow(gp,2)/3)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + 2*sin(thetat)*cos(thetat)*ft/(root2)*(At*CPEMix(higgs,1) - mueff*CPEMix(higgs,2)-lam*root2*mWboson*cos(beta)/g*CPEMix(higgs,3)));

  couplingst2 = 1/(2*pow(2*GFosqrt2,0.5)*pow(mstop2,2))*(pow(sin(thetat),2)*root2*(pow(ft,2)*root2*mWboson*sin(beta)/g*CPEMix(higgs,1) + (pow(gp,2)/12 - pow(g,2)/4)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + pow(cos(thetat),2)*root2*(pow(ft,2)*root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - (pow(gp,2)/3)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) - 2*sin(thetat)*cos(thetat)*ft/(root2)*(At*CPEMix(higgs,1) - mueff*CPEMix(higgs,2)-lam*root2*mWboson*cos(beta)/g*CPEMix(higgs,3)));

  couplingsb1 = 1/(2*pow(2*GFosqrt2,0.5)*pow(msbottom1,2))*(pow(cos(thetab),2)*root2*(pow(fb,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + (pow(gp,2)/12 + pow(g,2)/4)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + pow(sin(thetab),2)*root2*(pow(fb,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + pow(gp,2)/6*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + 2*sin(thetab)*cos(thetab)*fb/(root2)*(-mueff*CPEMix(higgs,1)+Ab*CPEMix(higgs,2)-lam*root2*mWboson*sin(beta)/g*CPEMix(higgs,3)));

  couplingsb2 = 1/(2*pow(2*GFosqrt2,0.5)*pow(msbottom2,2))*(pow(sin(thetab),2)*root2*(pow(fb,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + (pow(gp,2)/12 + pow(g,2)/4)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + pow(cos(thetab),2)*root2*(pow(fb,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + pow(gp,2)/6*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) - 2*sin(thetab)*cos(thetab)*fb/(root2)*(-mueff*CPEMix(higgs,1)+Ab*CPEMix(higgs,2)-lam*root2*mWboson*sin(beta)/g*CPEMix(higgs,3)));			
  
  couplingstau1 = 1/(2*pow(2*GFosqrt2,0.5)*pow(mstau1,2))*(pow(sin(thetatau),2)*root2*(pow(ftau,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + (-pow(gp,2)/4+pow(g,2)/4)*root2*mWboson/g*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2))) + pow(cos(thetatau),2)*root2*(pow(ftau,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + pow(gp,2)/2*(root2*mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2))) - 2*sin(thetatau)*cos(thetatau)*ftau/(root2)*(-mueff*CPEMix(higgs,1) + Atau*CPEMix(higgs,2) - lam*root2*mWboson*sin(beta)/g*CPEMix(higgs,3)));

  couplingstau2 = 1/(2*pow(2*GFosqrt2,0.5)*pow(mstau2,2))*(pow(cos(thetatau),2)*root2*(pow(ftau,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + (-pow(gp,2)/4+pow(g,2)/4)*root2*mWboson/g*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2))) + pow(sin(thetatau),2)*root2*(pow(ftau,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + pow(gp,2)/2*(root2*mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2))) + 2*sin(thetatau)*cos(thetatau)*ftau/(root2)*(-mueff*CPEMix(higgs,1) + Atau*CPEMix(higgs,2) - lam*root2*mWboson*sin(beta)/g*CPEMix(higgs,3)));

  kintr = 2*tfoftau(3)*(1 + (1-tfoftau(3))*tfoftau(1)); kincr = 2*cfoftau(3)*(1 + (1-cfoftau(3))*cfoftau(1)); kinbr = 2*bfoftau(3)*(1 + (1-bfoftau(3))*bfoftau(1)); kintaur = 2*taufoftau(3)*(1 + (1-taufoftau(3))*taufoftau(1));
  kinti = 2*tfoftau(3)*((1-tfoftau(3))*tfoftau(2)); kinci = 2*cfoftau(3)*((1-cfoftau(3))*cfoftau(2)); kinbi = 2*bfoftau(3)*((1-bfoftau(3))*bfoftau(2)); kintaui = 2*taufoftau(3)*((1-taufoftau(3))*taufoftau(2));
  
  kinWr = -(2+3*Wfoftau(3) +3*Wfoftau(3)*(2-Wfoftau(3))*Wfoftau(1)); kinHpmr = Hpmfoftau(3)*(Hpmfoftau(3)*Hpmfoftau(1) -1); kinch1r = 2*ch1foftau(3)*(1+ (1-ch1foftau(3))*ch1foftau(1)); kinch2r = 2*ch2foftau(3)*(1+ (1-ch2foftau(3))*ch2foftau(1));
  kinWi = -(3*Wfoftau(3)*(2-Wfoftau(3))*Wfoftau(2)); kinHpmi = Hpmfoftau(3)*(Hpmfoftau(3)*Hpmfoftau(2)); kinch1i = 2*ch1foftau(3)*((1-ch1foftau(3))*ch1foftau(2)); kinch2i = 2*ch2foftau(3)*((1-ch2foftau(3))*ch2foftau(2));

  kinscLr = scLfoftau(3)*(scLfoftau(3)*scLfoftau(1)-1); kinscRr = scRfoftau(3)*(scRfoftau(3)*scRfoftau(1)-1); kinssLr = ssLfoftau(3)*(ssLfoftau(3)*ssLfoftau(1)-1); kinssRr = ssRfoftau(3)*(ssRfoftau(3)*ssRfoftau(1)-1); kinsmuLr = smuLfoftau(3)*(smuLfoftau(3)*smuLfoftau(1)-1); kinsmuRr = smuRfoftau(3)*(smuRfoftau(3)*smuRfoftau(1)-1);
  kinscLi = scLfoftau(3)*(scLfoftau(3)*scLfoftau(2)); kinscRi = scRfoftau(3)*(scRfoftau(3)*scRfoftau(2)); kinssLi = ssLfoftau(3)*(ssLfoftau(3)*ssLfoftau(2)); kinssRi = ssRfoftau(3)*(ssRfoftau(3)*ssRfoftau(2));  kinsmuLi = smuLfoftau(3)*(smuLfoftau(3)*smuLfoftau(2)); kinsmuRi = smuRfoftau(3)*(smuRfoftau(3)*smuRfoftau(2));

  kinst1r = st1foftau(3)*(st1foftau(3)*st1foftau(1)-1); kinst2r = st2foftau(3)*(st2foftau(3)*st2foftau(1)-1); kinsb1r = sb1foftau(3)*(sb1foftau(3)*sb1foftau(1)-1); kinsb2r = sb2foftau(3)*(sb2foftau(3)*sb2foftau(1)-1); kinstau1r = stau1foftau(3)*(stau1foftau(3)*stau1foftau(1)-1); kinstau2r = stau2foftau(3)*(stau2foftau(3)*stau2foftau(1)-1);
  kinst1i = st1foftau(3)*(st1foftau(3)*st1foftau(2)); kinst2i = st2foftau(3)*(st2foftau(3)*st2foftau(2)); kinsb1i = sb1foftau(3)*(sb1foftau(3)*sb1foftau(2)); kinsb2i = sb2foftau(3)*(sb2foftau(3)*sb2foftau(2)); kinstau1i = stau1foftau(3)*(stau1foftau(3)*stau1foftau(2)); kinstau2i = stau2foftau(3)*(stau2foftau(3)*stau2foftau(2));

  Itr = couplingt*kintr; Icr = couplingc*kincr; Ibr = couplingb*kinbr; Itaur = couplingtau*kintaur; Ichar1r = couplingch1*kinch1r; Ichar2r = couplingch2*kinch2r; IWr = couplingW*kinWr; IHpmr = couplingHpm*kinHpmr;
  Iti = couplingt*kinti; Ici = couplingc*kinci; Ibi = couplingb*kinbi; Itaui = couplingtau*kintaui; Ichar1i = couplingch1*kinch1i; Ichar2i = couplingch2*kinch2i; IWi = couplingW*kinWi; IHpmi = couplingHpm*kinHpmi;

  IscLr = couplingscL*kinscLr; IscRr = couplingscR*kinscRr; IssLr = couplingssL*kinssLr; IssRr = couplingssR*kinssRr; IsmuLr = couplingsmuL*kinsmuLr; IsmuRr = couplingsmuR*kinsmuRr;
  IscLi = couplingscL*kinscLi; IscRi = couplingscR*kinscRi; IssLi = couplingssL*kinssLi; IssRi = couplingssR*kinssRi; IsmuLi = couplingsmuL*kinsmuLi; IsmuRi = couplingsmuR*kinsmuRi;

  Ist1r = couplingst1*kinst1r; Ist2r = couplingst2*kinst2r; Isb1r = couplingsb1*kinsb1r; Isb2r = couplingsb2*kinsb2r; Istau1r = couplingstau1*kinstau1r; Istau2r = couplingstau2*kinstau2r;
  Ist1i = couplingst1*kinst1i; Ist2i = couplingst2*kinst2i; Isb1i = couplingsb1*kinsb1i; Isb2i = couplingsb2*kinsb2i; Istau1i = couplingstau1*kinstau1i; Istau2i = couplingstau2*kinstau2i;

  DoubleVector matelemsum(2); double matelemmodsquare = 0;
  for (int i = 1; i <= 2; i++) {
    matelemsum(i) = 0;
  }
  matelemsum(1) = Itr + Ibr + Icr + Itaur + Ichar1r + Ichar2r + IWr + IHpmr + IscLr + IscRr + IssLr + IssRr + IsmuLr + IsmuRr + Ist1r + Ist2r + Isb1r + Isb2r + Istau1r + Istau2r;
  matelemsum(2) = Iti + Ibi + Ici + Itaui + Ichar1i + Ichar2i + IWi + IHpmi + IscLi + IscRi + IssLi + IssRi + IsmuLi + IsmuRi + Ist1i + Ist2i + Isb1i + Isb2i + Istau1i + Istau2i;

  prefactor = GFosqrt2/(4*PI)*pow(m1,3)/2*pow(alpha/PI,2)/16;

  matelemmodsquare = pow(matelemsum(1),2) + pow(matelemsum(2),2);

  amplitudeW = prefactor*matelemmodsquare;

  return amplitudeW;
}



double higgsCPevenamplitudedecaygluongluonNMSSM(double m1, double mtop, double mbottom, double mcharm, double mWboson, double mscharmL, double mscharmR, double mstop1, double mstop2, double msstrangeL, double msstrangeR, double msbottom1, double msbottom2, double msupL, double msupR, double msdownL, double msdownR, double runmt, double runmb, DoubleMatrix & CPEMix, double beta, double g, double gp, double gs, double alphas, double thetat, double thetab, double thetaL, double thetaR, double At, double Ab, double mu, double mueff, double lam, double kappa, double Alambda, int higgs, bool QCD)
{
  double amplitudeW = 0, prefactor = 0, couplingt = 0, couplingb = 0, couplingc = 0, couplingscL = 0, couplingscR = 0, couplingst1 = 0, couplingst2 = 0, couplingssL = 0, couplingssR = 0, couplingsb1 = 0, couplingsb2 = 0, couplingsuL = 0, couplingsuR = 0, couplingsdL = 0, couplingsdR = 0;
  double ft = 0, fb = 0;
  DoubleVector tfoftau(3), bfoftau(3), cfoftau(3), taufoftau(3), Wfoftau(3), Hpmfoftau(3), ch1foftau(3), ch2foftau(3), scLfoftau(3), scRfoftau(3), ssLfoftau(3), ssRfoftau(3), st1foftau(3), st2foftau(3), sb1foftau(3), sb2foftau(3), stau1foftau(3), stau2foftau(3), suLfoftau(3), suRfoftau(3), sdLfoftau(3), sdRfoftau(3);

  ///Initialise
  for (int i = 1; i<=3; i++) {
    tfoftau(i) = 0; bfoftau(i) = 0; cfoftau(i) = 0, scLfoftau(i) = 0, scRfoftau(i) = 0, ssLfoftau(i) = 0, ssRfoftau(i) = 0, st1foftau(i) = 0, st2foftau(i) = 0, sb1foftau(i) = 0, sb2foftau(i) = 0, suLfoftau(i) = 0, suRfoftau(i) = 0, sdLfoftau(i) = 0, sdRfoftau(i) = 0;
  }
  double kintr = 0, kinbr = 0, kincr = 0, kinscLr = 0, kinscRr = 0, kinssLr = 0, kinssRr = 0, kinst1r = 0, kinst2r = 0, kinsb1r = 0, kinsb2r = 0, kinsuLr = 0, kinsuRr = 0, kinsdLr = 0, kinsdRr = 0;
  double kinti = 0, kinbi = 0, kinci = 0, kinscLi = 0, kinscRi = 0, kinssLi = 0, kinssRi = 0, kinst1i = 0, kinst2i = 0, kinsb1i = 0, kinsb2i = 0, kinsuLi = 0, kinsuRi = 0, kinsdLi = 0, kinsdRi = 0;
  double Itr = 0, Icr = 0, Ibr = 0, Iti = 0, Ici = 0, Ibi = 0, IscLr = 0, IscRr = 0, IssLr = 0, IssRr = 0, IscLi = 0, IscRi = 0, IssLi = 0, IssRi = 0, Ist1r = 0, Ist2r = 0, Isb1r = 0, Isb2r = 0, Ist1i = 0, Ist2i = 0, Isb1i = 0, Isb2i = 0, IsuLr = 0, IsuRr = 0, IsuLi = 0, IsuRi = 0, IsdLr = 0, IsdRr = 0, IsdLi = 0, IsdRi = 0;

  tfoftau = foftau(mtop, m1); bfoftau = foftau(mbottom, m1); cfoftau = foftau(mcharm, m1); scLfoftau = foftau(mscharmL, m1); scRfoftau = foftau(mscharmR, m1); ssLfoftau = foftau(msstrangeL, m1); ssRfoftau = foftau(msstrangeR, m1); st1foftau = foftau(mstop1, m1); st2foftau = foftau(mstop2, m1); sb1foftau = foftau(msbottom1, m1); sb2foftau = foftau(msbottom2, m1); suLfoftau = foftau(msupL, m1); suRfoftau = foftau(msupR, m1); sdLfoftau = foftau(msdownL, m1); sdRfoftau = foftau(msdownR, m1);

  couplingt = CPEMix(higgs,1)/sin(beta); couplingc = CPEMix(higgs,1)/sin(beta); couplingb = CPEMix(higgs,2)/(cos(beta));

  couplingscL = 2*mWboson/(g*pow(mscharmL,2))*(pow(gp,2)/12-pow(g,2)/4)*2*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));
  couplingscR = 2*mWboson/(g*pow(mscharmR,2))*(-pow(gp,2)/3)*2*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));

  couplingssL = 2*mWboson/(g*pow(msstrangeL,2))*2*(pow(gp,2)/12 + pow(g,2)/4)*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));
  couplingssR = 2*mWboson/(g*pow(msstrangeR,2))*(pow(gp,2)/6)*2*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));

  couplingsuL = 2*mWboson/(g*pow(msupL,2))*(pow(gp,2)/12-pow(g,2)/4)*2*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));
  couplingsuR = 2*mWboson/(g*pow(msupR,2))*(-pow(gp,2)/3)*2*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));
  couplingsdL = 2*mWboson/(g*pow(msdownL,2))*2*(pow(gp,2)/12 + pow(g,2)/4)*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));
  couplingsdR = 2*mWboson/(g*pow(msdownR,2))*(pow(gp,2)/6)*2*(mWboson/g)*(sin(beta)*CPEMix(higgs,1) - cos(beta)*CPEMix(higgs,2));

  ft = g*runmt/(root2*mWboson*sin(beta)); fb = g*runmb/(root2*mWboson*cos(beta));

  couplingst1 = mWboson/(g*pow(mstop1,2))*(pow(cos(thetat),2)*root2*(pow(ft,2)*root2*mWboson*sin(beta)/g*CPEMix(higgs,1) + (pow(gp,2)/12 - pow(g,2)/4)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + pow(sin(thetat),2)*root2*(pow(ft,2)*root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - (pow(gp,2)/3)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + 2*sin(thetat)*cos(thetat)*ft/(root2)*(At*CPEMix(higgs,1) - mueff*CPEMix(higgs,2)-lam*root2*mWboson*cos(beta)/g*CPEMix(higgs,3)));

  couplingst2 = mWboson/(g*pow(mstop2,2))*(pow(sin(thetat),2)*root2*(pow(ft,2)*root2*mWboson*sin(beta)/g*CPEMix(higgs,1) + (pow(gp,2)/12 - pow(g,2)/4)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + pow(cos(thetat),2)*root2*(pow(ft,2)*root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - (pow(gp,2)/3)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) - 2*sin(thetat)*cos(thetat)*ft/(root2)*(At*CPEMix(higgs,1) - mueff*CPEMix(higgs,2)-lam*root2*mWboson*cos(beta)/g*CPEMix(higgs,3)));

  couplingsb1 = mWboson/(g*pow(msbottom1,2))*(pow(cos(thetab),2)*root2*(pow(fb,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + (pow(gp,2)/12 + pow(g,2)/4)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + pow(sin(thetab),2)*root2*(pow(fb,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + pow(gp,2)/6*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + 2*sin(thetab)*cos(thetab)*fb/(root2)*(-mueff*CPEMix(higgs,1)+Ab*CPEMix(higgs,2)-lam*root2*mWboson*sin(beta)/g*CPEMix(higgs,3)));

  couplingsb2 = mWboson/(g*pow(msbottom2,2))*(pow(sin(thetab),2)*root2*(pow(fb,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + (pow(gp,2)/12 + pow(g,2)/4)*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) + pow(cos(thetab),2)*root2*(pow(fb,2)*root2*mWboson*cos(beta)/g*CPEMix(higgs,2) + pow(gp,2)/6*(root2*mWboson*sin(beta)/g*CPEMix(higgs,1) - root2*mWboson*cos(beta)/g*CPEMix(higgs,2))) - 2*sin(thetab)*cos(thetab)*fb/(root2)*(-mueff*CPEMix(higgs,1)+Ab*CPEMix(higgs,2)-lam*root2*mWboson*sin(beta)/g*CPEMix(higgs,3)));	
  
  kintr = 2*tfoftau(3)*(1 + (1-tfoftau(3))*tfoftau(1)); kincr = 2*cfoftau(3)*(1 + (1-cfoftau(3))*cfoftau(1)); kinbr = 2*bfoftau(3)*(1 + (1-bfoftau(3))*bfoftau(1));
  kinti = 2*tfoftau(3)*((1-tfoftau(3))*tfoftau(2)); kinci = 2*cfoftau(3)*((1-cfoftau(3))*cfoftau(2)); kinbi = 2*bfoftau(3)*((1-bfoftau(3))*bfoftau(2));
  
  kinscLr = scLfoftau(3)*(scLfoftau(3)*scLfoftau(1)-1); kinscRr = scRfoftau(3)*(scRfoftau(3)*scRfoftau(1)-1); kinssLr = ssLfoftau(3)*(ssLfoftau(3)*ssLfoftau(1)-1); kinssRr = ssRfoftau(3)*(ssRfoftau(3)*ssRfoftau(1)-1);
  kinscLi = scLfoftau(3)*(scLfoftau(3)*scLfoftau(2)); kinscRi = scRfoftau(3)*(scRfoftau(3)*scRfoftau(2)); kinssLi = ssLfoftau(3)*(ssLfoftau(3)*ssLfoftau(2)); kinssRi = ssRfoftau(3)*(ssRfoftau(3)*ssRfoftau(2));

  kinsuLr = suLfoftau(3)*(suLfoftau(3)*suLfoftau(1)-1); kinsuRr = suRfoftau(3)*(suRfoftau(3)*suRfoftau(1)-1); kinsdLr = sdLfoftau(3)*(sdLfoftau(3)*sdLfoftau(1)-1); kinsdRr = sdRfoftau(3)*(sdRfoftau(3)*sdRfoftau(1)-1);
  kinsuLi = suLfoftau(3)*(suLfoftau(3)*suLfoftau(2)); kinsuRi = suRfoftau(3)*(suRfoftau(3)*suRfoftau(2)); kinsdLi = sdLfoftau(3)*(sdLfoftau(3)*sdLfoftau(2)); kinsdRi = sdRfoftau(3)*(sdRfoftau(3)*sdRfoftau(2));

  kinst1r = st1foftau(3)*(st1foftau(3)*st1foftau(1)-1); kinst2r = st2foftau(3)*(st2foftau(3)*st2foftau(1)-1); kinsb1r = sb1foftau(3)*(sb1foftau(3)*sb1foftau(1)-1); kinsb2r = sb2foftau(3)*(sb2foftau(3)*sb2foftau(1)-1);
  kinst1i = st1foftau(3)*(st1foftau(3)*st1foftau(2)); kinst2i = st2foftau(3)*(st2foftau(3)*st2foftau(2)); kinsb1i = sb1foftau(3)*(sb1foftau(3)*sb1foftau(2)); kinsb2i = sb2foftau(3)*(sb2foftau(3)*sb2foftau(2));

  Itr = couplingt*kintr; Icr = couplingc*kincr; Ibr = couplingb*kinbr;
  Iti = couplingt*kinti; Ici = couplingc*kinci; Ibi = couplingb*kinbi;

  IscLr = couplingscL*kinscLr; IscRr = couplingscR*kinscRr; IssLr = couplingssL*kinssLr; IssRr = couplingssR*kinssRr;
  IscLi = couplingscL*kinscLi; IscRi = couplingscR*kinscRi; IssLi = couplingssL*kinssLi; IssRi = couplingssR*kinssRi;

  IsuLr = couplingsuL*kinsuLr; IsuRr = couplingsuR*kinsuRr; IsdLr = couplingsdL*kinsdLr; IsdRr = couplingsdR*kinsdRr;
  IsuLi = couplingsuL*kinsuLi; IsuRi = couplingsuR*kinsuRi; IsdLi = couplingsdL*kinsdLi; IsdRi = couplingsdR*kinsdRi;

  Ist1r = couplingst1*kinst1r; Ist2r = couplingst2*kinst2r; Isb1r = couplingsb1*kinsb1r; Isb2r = couplingsb2*kinsb2r;
  Ist1i = couplingst1*kinst1i; Ist2i = couplingst2*kinst2i; Isb1i = couplingsb1*kinsb1i; Isb2i = couplingsb2*kinsb2i;

  DoubleVector matelemsum(2); double matelemmodsquare = 0;
  for (int i = 1; i <= 2; i++) {
    matelemsum(i) = 0;
  }
  matelemsum(1) = Itr + Ibr + Icr + IscLr + IscRr + IssLr + IssRr + IsuLr + IsuRr + IsdLr + IsdRr + Ist1r + Ist2r + Isb1r + Isb2r;
  matelemsum(2) = Iti + Ibi + Ici + IscLi + IscRi + IssLi + IssRi + IsuLi + IsuRi + IsdLi + IsdRi + Ist1i + Ist2i + Isb1i + Isb2i;

  prefactor =  GFosqrt2*1/(4*PI)*pow(m1,3)/2*pow(alphas/PI,2)/8;

  matelemmodsquare = pow(matelemsum(1),2) + pow(matelemsum(2),2);
  
  amplitudeW = prefactor*matelemmodsquare;

  
  double SMTOTRE = 0, SMTOTIM = 0, SQTOTRE = 0, SQTOTIM = 0;
  SMTOTRE = Itr + Ibr + Icr;
  SQTOTRE = IscLr + IscRr + IssLr + IssRr + IsuLr + IsuRr + IsdLr + IsdRr + Ist1r + Ist2r + Isb1r + Isb2r;
  SMTOTIM = Iti + Ibi + Ici;
  SQTOTIM = IscLi + IscRi + IssLi + IssRi + IsuLi + IsuRi + IsdLi + IsdRi + Ist1i + Ist2i + Isb1i + Isb2i;

  if (QCD == true) {
    int NF = 0;
    if (m1 < mtop) { NF = 5;}
    else { NF = 6;}
    // std::cout << "NF = " << NF << std::endl;
    if (higgs == 1 || higgs == 2 || higgs == 3) {}
    else {
      throw("Problem - higgs must be 1, 2, 3 i.e 'h', 'H' or 'H3' in NMSSM!\n");
    }

    amplitudeW = hggQCDcorrections(amplitudeW, alphas, NF, 'h', prefactor, SMTOTRE, SMTOTIM, SQTOTRE, SQTOTIM)(1); ///Pass hggQCDcorrections 'h' as it just matters whether it is CPeven (so 95/4 in FQCD) or CPodd (so 97/4 in FQCD)
   
  }
  else if (QCD == false) {
    amplitudeW = amplitudeW;
  }
     
  return amplitudeW;
}


double higgshamplitudedecayZgammaNMSSM (double m1, double g, double gp, double alpha, double mWboson, double mZboson, double mHpm, DoubleMatrix & CPEMix, double beta, double mtop, double mbottom, double mcharm, double mch1, double mch2, double thetaL, double thetaR, double lam, double kappa, double Alambda, double greekmu, double mueff, int higgs)
{
  double prefactor=0, Itr=0, Iti=0, Ibr = 0, Ibi = 0, Icr = 0, Ici = 0, Ichar1r = 0, Ichar1i = 0, Ichar2r = 0, Ichar2i = 0, IWr = 0, IWi = 0, IHpmr = 0, IHpmi = 0, couplingt = 0, couplingb = 0, couplingc = 0, couplingch1 =0, couplingch2 = 0, couplingW = 0, couplingHpm = 0, kintr = 0, kinti = 0, kinbr =0, kinbi = 0, kincr = 0, kinci = 0, kinch1r = 0, kinch1i = 0, kinch2r = 0, kinch2i = 0, kinWr = 0, kinWi = 0, kinHpmr = 0, kinHpmi = 0, matelemmodsquare=0, amplitudeW=0;
  double sin2thW = 0, sinthW = 0, costhW = 0;

  if (m1 < mZboson) { amplitudeW = 0;}
  else {
    
    DoubleVector tfoftau(3), bfoftau(3), cfoftau(3), ch1foftau(3), ch2foftau(3), Wfoftau(3), Hpmfoftau(3), tZfoftau(3), bZfoftau(3), cZfoftau(3), ch1Zfoftau(3), ch2Zfoftau(3), WZfoftau(3), HpmZfoftau(3);
    DoubleVector tgoftau(3), bgoftau(3), cgoftau(3), ch1goftau(3), ch2goftau(3), Wgoftau(3), Hpmgoftau(3), tZgoftau(3), bZgoftau(3), cZgoftau(3), ch1Zgoftau(3), ch2Zgoftau(3), WZgoftau(3), HpmZgoftau(3);
    ///Initialise these components
    for (int i = 1; i <= 3; i++) {
      tfoftau(i) = 0; bfoftau(i) = 0; cfoftau(i) = 0; ch1foftau(i) = 0; ch2foftau(i) = 0; Wfoftau(i) = 0; Hpmfoftau(i) = 0; tZfoftau(i) = 0; bZfoftau(i) = 0; cZfoftau(i) = 0; ch1Zfoftau(i) = 0; ch2Zfoftau(i) = 0; WZfoftau(i) = 0; HpmZfoftau(i) = 0;
      tgoftau(i) = 0; bgoftau(i) = 0; cgoftau(i) = 0; ch1goftau(i) = 0; ch2goftau(i) = 0; Wgoftau(i) = 0; Hpmgoftau(i) = 0; tZgoftau(i) = 0; bZgoftau(i) = 0; cZgoftau(i) = 0; ch1Zgoftau(i) = 0; ch2Zgoftau(i) = 0; WZgoftau(i) = 0; HpmZgoftau(i) = 0;
    }
    tfoftau = foftau(mtop, m1); bfoftau = foftau(mbottom, m1); cfoftau = foftau(mcharm, m1); ch1foftau = foftau(mch1, m1); ch2foftau = foftau(mch2, m1); Wfoftau = foftau(mWboson, m1); Hpmfoftau = foftau(mHpm, m1);
    tZfoftau = foftau(mtop, mZboson); bZfoftau = foftau(mbottom, mZboson); cZfoftau = foftau(mcharm, mZboson); ch1Zfoftau = foftau(mch1, mZboson); ch2Zfoftau = foftau(mch2, mZboson); WZfoftau = foftau(mWboson, mZboson); HpmZfoftau = foftau(mHpm, mZboson);
    
    tgoftau = goftau(mtop, m1); bgoftau = goftau(mbottom, m1); cgoftau = goftau(mcharm, m1); ch1goftau = goftau(mch1, m1); ch2goftau = goftau(mch2, m1); Wgoftau = goftau(mWboson, m1); Hpmgoftau = goftau(mHpm, m1);
    tZgoftau = goftau(mtop, mZboson); bZgoftau = goftau(mbottom, mZboson); cZgoftau = goftau(mcharm, mZboson); ch1Zgoftau = goftau(mch1, mZboson); ch2Zgoftau = goftau(mch2, mZboson); WZgoftau = goftau(mWboson, mZboson); HpmZgoftau = goftau(mHpm, mZboson);

    sin2thW = pow(gp,2)/(pow(g,2)+pow(gp,2)); sinthW = pow(sin2thW,0.5); costhW = g/(pow(pow(g,2)+pow(gp,2),0.5));
    
    couplingt = -2*(1-8*sin2thW/3)/(sinthW*costhW)*CPEMix(higgs,1)/(sin(beta)); couplingc = -2*(1-8*sin2thW/3)/(sinthW*costhW)*CPEMix(higgs,1)/(sin(beta)); couplingb =(-1+4*sin2thW/3)/(sinthW*costhW)*CPEMix(higgs,2)/(cos(beta));

    couplingch1 = 4*mWboson/(mch1*g*sinthW*costhW)*(lam/(root2)*CPEMix(higgs,3)*cos(thetaL)*cos(thetaR) + g/(root2)*(CPEMix(higgs,1)*sin(thetaL)*cos(thetaR) + CPEMix(higgs,2)*cos(thetaL)*sin(thetaR)))*(-pow(sin(thetaR),2) - 0.5*pow(cos(thetaR),2) + 2*pow(sinthW,2) - pow(sin(thetaL),2) - 0.5*pow(cos(thetaL),2));
    couplingch2 = 4*mWboson/(mch2*g*sinthW*costhW)*(lam/(root2)*CPEMix(higgs,3)*sin(thetaL)*sin(thetaR) - g/(root2)*(CPEMix(higgs,1)*cos(thetaL)*sin(thetaR) + CPEMix(higgs,2)*sin(thetaL)*cos(thetaR)))*(-pow(cos(thetaR),2) - 0.5*pow(sin(thetaR),2) + 2*pow(sinthW,2) - pow(cos(thetaL),2) - 0.5*pow(sin(thetaL),2));

    couplingHpm = (1-2*sin2thW)/(sinthW*costhW*2*pow(mHpm,2))*1/(pow(GFosqrt2*2,0.5))*(lam*mueff/(root2)*(2*CPEMix(higgs,3)*pow(cos(beta),2) + 2*CPEMix(higgs,3)*pow(sin(beta),2)) - pow(lam,2)*mWboson*sin(beta)/g*2*CPEMix(higgs,2)*cos(beta)*sin(beta) - pow(lam,2)*mWboson*sin(beta)/g*2*CPEMix(higgs,1)*cos(beta)*sin(beta) + mueff*kappa*root2*2*CPEMix(higgs,3)*cos(beta)*sin(beta) + lam*Alambda/root2*2*CPEMix(higgs,3)*cos(beta)*sin(beta) + pow(gp,2)/4*mWboson/g*(sin(beta)*(2*CPEMix(higgs,1)*pow(cos(beta),2) - 2*CPEMix(higgs,1)*pow(sin(beta),2)) + cos(beta)*(2*CPEMix(higgs,2)*pow(sin(beta),2) - 2*CPEMix(higgs,2)*pow(cos(beta),2))) + g/4*mWboson*(sin(beta)*(2*CPEMix(higgs,1)*pow(cos(beta),2) + 2*CPEMix(higgs,1)*pow(sin(beta),2) + 2*2*CPEMix(higgs,2)*sin(beta)*cos(beta)) + cos(beta)*(2*CPEMix(higgs,2)*pow(cos(beta),2) + 2*CPEMix(higgs,2)*pow(sin(beta),2) + 2*2*CPEMix(higgs,1)*sin(beta)*cos(beta))) + lam/root2*0); ///ignored corrections for now - no top or bottom loop corrections ///MUP in NMSSMTools is a fine-tuning parameter - set 0

    couplingW = -g/gp*(CPEMix(higgs,1)*sin(beta) + CPEMix(higgs,2)*cos(beta));

    kintr = tfoftau(3)*tZfoftau(3)/(2*(tfoftau(3)-tZfoftau(3))) + pow(tfoftau(3)*tZfoftau(3),2)/(2*(pow(tfoftau(3)-tZfoftau(3),2)))*(tfoftau(1) - tZfoftau(1)) + pow(tfoftau(3),2)*tZfoftau(3)/(pow(tfoftau(3)-tZfoftau(3),2))*(tgoftau(1)-tZgoftau(1)) + tfoftau(3)*tZfoftau(3)/(2*(tfoftau(3)-tZfoftau(3)))*(tfoftau(1)-tZfoftau(1));
    kinti = pow(tfoftau(3)*tZfoftau(3),2)/(2*(pow(tfoftau(3)-tZfoftau(3),2)))*(tfoftau(2) - tZfoftau(2)) + pow(tfoftau(3),2)*tZfoftau(3)/(pow(tfoftau(3)-tZfoftau(3),2))*(tgoftau(2)-tZgoftau(2)) + tfoftau(3)*tZfoftau(3)/(2*(tfoftau(3)-tZfoftau(3)))*(tfoftau(2)-tZfoftau(2));
    kinbr = bfoftau(3)*bZfoftau(3)/(2*(bfoftau(3)-bZfoftau(3))) + pow(bfoftau(3)*bZfoftau(3),2)/(2*(pow(bfoftau(3)-bZfoftau(3),2)))*(bfoftau(1) - bZfoftau(1)) + pow(bfoftau(3),2)*bZfoftau(3)/(pow(bfoftau(3)-bZfoftau(3),2))*(bgoftau(1)-bZgoftau(1)) + bfoftau(3)*bZfoftau(3)/(2*(bfoftau(3)-bZfoftau(3)))*(bfoftau(1)-bZfoftau(1));
    kinbi = pow(bfoftau(3)*bZfoftau(3),2)/(2*(pow(bfoftau(3)-bZfoftau(3),2)))*(bfoftau(2) - bZfoftau(2)) + pow(bfoftau(3),2)*bZfoftau(3)/(pow(bfoftau(3)-bZfoftau(3),2))*(bgoftau(2)-bZgoftau(2)) + bfoftau(3)*bZfoftau(3)/(2*(bfoftau(3)-bZfoftau(3)))*(bfoftau(2)-bZfoftau(2));
    kincr = cfoftau(3)*cZfoftau(3)/(2*(cfoftau(3)-cZfoftau(3))) + pow(cfoftau(3)*cZfoftau(3),2)/(2*(pow(cfoftau(3)-cZfoftau(3),2)))*(cfoftau(1) - cZfoftau(1)) + pow(cfoftau(3),2)*cZfoftau(3)/(pow(cfoftau(3)-cZfoftau(3),2))*(cgoftau(1)-cZgoftau(1)) + cfoftau(3)*cZfoftau(3)/(2*(cfoftau(3)-cZfoftau(3)))*(cfoftau(1)-cZfoftau(1));
    kinci = pow(cfoftau(3)*cZfoftau(3),2)/(2*(pow(cfoftau(3)-cZfoftau(3),2)))*(cfoftau(2) - cZfoftau(2)) + pow(cfoftau(3),2)*cZfoftau(3)/(pow(cfoftau(3)-cZfoftau(3),2))*(cgoftau(2)-cZgoftau(2)) + cfoftau(3)*cZfoftau(3)/(2*(cfoftau(3)-cZfoftau(3)))*(cfoftau(2)-cZfoftau(2));

    kinWr = 4*(3-pow(gp/g,2))*-Wfoftau(3)*WZfoftau(3)/(2*(Wfoftau(3)-WZfoftau(3)))*(Wfoftau(1) - WZfoftau(1)) + ((1 + 2/Wfoftau(3))*pow(gp/g,2) - (5+ 2/Wfoftau(3)))*(Wfoftau(3)*WZfoftau(3)/(2*(Wfoftau(3)-WZfoftau(3))) + pow(Wfoftau(3)*WZfoftau(3),2)/(2*pow(Wfoftau(3)-WZfoftau(3),2))*(Wfoftau(1)-WZfoftau(1)) + pow(Wfoftau(3),2)*WZfoftau(3)/(pow(Wfoftau(3)-WZfoftau(3),2))*(Wgoftau(1)-WZgoftau(1)));
    kinWi = 4*(3-pow(gp/g,2))*-Wfoftau(3)*WZfoftau(3)/(2*(Wfoftau(3)-WZfoftau(3)))*(Wfoftau(2) - WZfoftau(2)) + ((1 + 2/Wfoftau(3))*pow(gp/g,2) - (5+ 2/Wfoftau(3)))*(pow(Wfoftau(3)*WZfoftau(3),2)/(2*pow(Wfoftau(3)-WZfoftau(3),2))*(Wfoftau(2)-WZfoftau(2)) + pow(Wfoftau(3),2)*WZfoftau(3)/(pow(Wfoftau(3)-WZfoftau(3),2))*(Wgoftau(2)-WZgoftau(2)));

    kinHpmr = Hpmfoftau(3)*HpmZfoftau(3)/(2*(Hpmfoftau(3)-HpmZfoftau(3))) + pow(Hpmfoftau(3)*HpmZfoftau(3),2)/(2*pow(Hpmfoftau(3)-HpmZfoftau(3),2))*(Hpmfoftau(1)-HpmZfoftau(1)) + pow(Hpmfoftau(3),2)*HpmZfoftau(3)/(pow(Hpmfoftau(3)-HpmZfoftau(3),2))*(Hpmgoftau(1)-HpmZgoftau(1));
    kinHpmi = pow(Hpmfoftau(3)*HpmZfoftau(3),2)/(2*pow(Hpmfoftau(3)-HpmZfoftau(3),2))*(Hpmfoftau(2)-HpmZfoftau(2)) + pow(Hpmfoftau(3),2)*HpmZfoftau(3)/(pow(Hpmfoftau(3)-HpmZfoftau(3),2))*(Hpmgoftau(2)-HpmZgoftau(2));
   
    kinch1r = ch1foftau(3)*ch1Zfoftau(3)/(2*(ch1foftau(3)-ch1Zfoftau(3))) + pow(ch1foftau(3)*ch1Zfoftau(3),2)/(2*pow(ch1foftau(3)-ch1Zfoftau(3),2))*(ch1foftau(1) - ch1Zfoftau(1)) + pow(ch1foftau(3),2)*ch1Zfoftau(3)/pow(ch1foftau(3) - ch1Zfoftau(3),2)*(ch1goftau(1) - ch1Zgoftau(1)) + ch1foftau(3)*ch1Zfoftau(3)/(2*(ch1foftau(3)-ch1Zfoftau(3)))*(ch1foftau(1)-ch1Zfoftau(1));
    kinch1i = pow(ch1foftau(3)*ch1Zfoftau(3),2)/(2*pow(ch1foftau(3)-ch1Zfoftau(3),2))*(ch1foftau(2) - ch1Zfoftau(2)) + pow(ch1foftau(3),2)*ch1Zfoftau(3)/pow(ch1foftau(3) - ch1Zfoftau(3),2)*(ch1goftau(2) - ch1Zgoftau(2)) + ch1foftau(3)*ch1Zfoftau(3)/(2*(ch1foftau(3)-ch1Zfoftau(3)))*(ch1foftau(2)-ch1Zfoftau(2));
    kinch2r = ch2foftau(3)*ch2Zfoftau(3)/(2*(ch2foftau(3)-ch2Zfoftau(3))) + pow(ch2foftau(3)*ch2Zfoftau(3),2)/(2*pow(ch2foftau(3)-ch2Zfoftau(3),2))*(ch2foftau(1) - ch2Zfoftau(1)) + pow(ch2foftau(3),2)*ch2Zfoftau(3)/pow(ch2foftau(3) - ch2Zfoftau(3),2)*(ch2goftau(1) - ch2Zgoftau(1)) + ch2foftau(3)*ch2Zfoftau(3)/(2*(ch2foftau(3)-ch2Zfoftau(3)))*(ch2foftau(1)-ch2Zfoftau(1));
    kinch2i = pow(ch2foftau(3)*ch2Zfoftau(3),2)/(2*pow(ch2foftau(3)-ch2Zfoftau(3),2))*(ch2foftau(2) - ch2Zfoftau(2)) + pow(ch2foftau(3),2)*ch2Zfoftau(3)/pow(ch2foftau(3) - ch2Zfoftau(3),2)*(ch2goftau(2) - ch2Zgoftau(2)) + ch2foftau(3)*ch2Zfoftau(3)/(2*(ch2foftau(3)-ch2Zfoftau(3)))*(ch2foftau(2)-ch2Zfoftau(2));
       
    Itr = couplingt*kintr; Icr = couplingc*kincr; Ibr = couplingb*kinbr; Ichar1r = couplingch1*kinch1r; Ichar2r = couplingch2*kinch2r; 
    Iti = couplingt*kinti; Ici = couplingc*kinci; Ibi = couplingb*kinbi; Ichar1i = couplingch1*kinch1i; Ichar2i = couplingch2*kinch2i; 
    IWr = couplingW*kinWr; IWi = couplingW*kinWi; IHpmr = couplingHpm*kinHpmr; IHpmi = couplingHpm*kinHpmi;

    DoubleVector matelemsum(2);
    for (int i = 1; i <= 2; i++) {
      matelemsum(i) = 0;
    }
    matelemsum(1) = Itr + Ibr + Icr + Ichar1r + Ichar2r + IWr + IHpmr;
    matelemsum(2) = Iti + Ibi + Ici + Ichar1i + Ichar2i + IWi + IHpmi;
    
    prefactor = GFosqrt2*pow(m1,3)*pow(alpha,2)/(64*pow(PI,3))*pow((1-pow(mZboson/m1,2)),3);

    matelemmodsquare = pow(matelemsum(1),2) + pow(matelemsum(2),2);
    amplitudeW = prefactor*matelemmodsquare;

  }  
  return amplitudeW;

}

 

double higgshamplitudedecayneutineutjNMSSM (double m1, double mneuti, double mneutj, double g, double gp, DoubleMatrix & CPEMix, DoubleMatrix & mixNeut, double lam, double kappa, int neuti, int neutj, int higgs) 
{
  double amplitudeW = 0;
  double coupling = 0, squareplus = 0, squareminus = 0, lambda = 0, factor = 0;
  if (m1 < fabs(mneuti)+fabs(mneutj)) { amplitudeW = 0;}
  else {
    if(neuti != neutj) { factor = 2;} ///accounts for decay being neuti neutj or neutj neuti when i!=j
    else if (neuti == neutj) { factor = 1;} ///for same neutralino they are QM indistinguishable so don't get factor of 2
    squareplus = 1 - pow(mneutj/m1+mneuti/m1,2);
    squareminus = 1 - pow(mneutj/m1-mneuti/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgshamplitudedecayneutineutjNMSSM\n");
    } 
    
    coupling = lam/(root2)*(CPEMix(higgs,1)*(mixNeut(neuti,3)*mixNeut(neutj,5) + mixNeut(neuti,5)*mixNeut(neutj,3)) + CPEMix(higgs,2)*(mixNeut(neuti,4)*mixNeut(neutj,5) + mixNeut(neuti,5)*mixNeut(neutj,4)) + CPEMix(higgs,3)*(mixNeut(neuti,4)*mixNeut(neutj,3)+mixNeut(neuti,3)*mixNeut(neutj,4))) - root2*kappa*CPEMix(higgs,3)*mixNeut(neuti,5)*mixNeut(neutj,5) + gp/2*(-CPEMix(higgs,1)*(mixNeut(neuti,1)*mixNeut(neutj,4) + mixNeut(neuti,4)*mixNeut(neutj,1)) + CPEMix(higgs,2)*(mixNeut(neuti,1)*mixNeut(neutj,3) + mixNeut(neuti,3)*mixNeut(neutj,1))) + g/2*(CPEMix(higgs,1)*(mixNeut(neuti,2)*mixNeut(neutj,4) + mixNeut(neuti,4)*mixNeut(neutj,2)) - CPEMix(higgs,2)*(mixNeut(neuti,2)*mixNeut(neutj,3) + mixNeut(neuti,3)*mixNeut(neutj,2)));

    amplitudeW = factor*squareplus*sqr(m1)*lambda*pow(coupling,2)/(16*PI*fabs(m1));
  }
  return amplitudeW;
}



double higgsAamplitudedecayHpmWboson(double m1, double mWboson, double mHpm, double g, double thetaA, int pseudoscalar, bool nmssmIsIt) ///Does A/A2 -> H+W- (note multiply by 2 if want to include H-W+ as well as H+W- in this)
{
  double amplitudeW = 0;
  if (m1 < mWboson + mHpm) {amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, coupling = 0;
    squareplus = 1 - pow(mHpm/m1 + mWboson/m1,2);
    squareminus = 1 - pow(mHpm/m1 - mWboson/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsAamplitudedecayHpmWboson\n");
    } 
    
    if (nmssmIsIt == false) { coupling = 1;}
    else if (nmssmIsIt == true) { 
      if( pseudoscalar == 1) {
	coupling = cos(thetaA);
      }
      else if (pseudoscalar == 2) {
	coupling = sin(thetaA);
      }
      else{
	throw("problem: pseudoscalar must be 1 or 2 in higgsAamplitudedecayHpmWboson");
      }
    }

    amplitudeW = GFosqrt2*pow(fabs(m1),3)*pow(coupling,2)/(8*PI)*pow(lambda,3);
  }
  return amplitudeW;
}


double CPEhCPOACPOACoupling  (DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, int higgs, int pseudoscalar1, int pseudoscalar2, int j, int k, int l)
{
  double coupling = 0;
  coupling = CPEMix(higgs,j)*(CPOMix(pseudoscalar1,k)*CPOMix(pseudoscalar2,l) + CPOMix(pseudoscalar1,l)*CPOMix(pseudoscalar2,k));
  return coupling;
}


double higgsCPevenamplitudedecayAANMSSM(double m1, double mA1, double mA2, double mWboson, double runmt, double runmb, double g, double gp, double beta, DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, double lam, double kappa, double Alambda, double Akappa, double mueff, int higgs, int pseudoscalar1, int pseudoscalar2) 
{
  double amplitudeW = 0;
  if (fabs(m1) < fabs(mA1) + fabs(mA2)) {
    amplitudeW = 0;
  }
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, hvev1 = 0, hvev2 = 0, prefactor = 0, coupling = 0;
    squareplus = 1 - pow(mA1/m1 + mA2/m1,2);
    squareminus = 1 - pow(mA1/m1 - mA2/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecayAANMSSM\n");
    } 
    
    prefactor = 1/(32*PI*m1);
    hvev1 = root2*mWboson*sin(beta)/g;
    hvev2 = root2*mWboson*cos(beta)/g;

    coupling = (pow(g,2)+pow(gp,2))/(4*root2)*(hvev1*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,1,1) - CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,2,2)) + hvev2*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,2,2) - CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,1,1))) + lam*Alambda/(root2)*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,2,3)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,1,3)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,1,2)) - kappa*Akappa/(root2)*CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,3,3) + pow(lam,2)/root2*(hvev1*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,2,2)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,3,3)) + hvev2*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,1,1)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,3,3)) + mueff/lam*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,1,1)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,2,2))) + pow(kappa,2)*root2*mueff/lam*CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,3,3) + lam*kappa/(root2)*(hvev1*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,3,3)-2*CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,2,3)) + hvev2*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,3,3)-2*CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,1,3)) + 2*mueff/lam*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,1,2)-CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,2,3) - CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,1,3)));

    amplitudeW = prefactor*lambda*pow(coupling,2); ///Multiplied by 2 later if A A2 to account for indistinguishability

  }
  return amplitudeW;
}
    

double higgsCPevenamplitudedecaypseudoscalarZNMSSM (double m1, double mA, double mZboson, double g, double gp, double beta, DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, int higgs, int pseudoscalar) 
{
  double amplitudeW = 0;
  if (m1 < mA + mZboson) { amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, coupling = 0;
    squareplus = 1 - pow(mA/m1+mZboson/m1,2);
    squareminus = 1 - pow(mA/m1-mZboson/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecaypseudoscalarZNMSSM\n");
    } 
    coupling = (CPEMix(higgs,1)*cos(beta) - CPEMix(higgs,2)*sin(beta))*(CPOMix(pseudoscalar,1)*cos(beta) - CPOMix(pseudoscalar,2)*sin(beta));
    amplitudeW = (pow(g,2)+pow(gp,2))*pow(m1,3)/(64*PI*pow(mZboson,2))*pow(coupling,2)*pow(lambda,3);

  }
  return amplitudeW;
}


double higgsCPevenamplitudedecayHpHmNMSSM (double m1, double mHpm, double mWboson,  double g, double gp, double mtop, double mbottom, double beta, double lam, double mueff, double kappa, double Alambda, DoubleMatrix & CPEMix, int higgs)
{
  double amplitudeW = 0;
  if (m1 < 2*mHpm) {amplitudeW = 0;}
  else {
    double squareplus = 0, lambda = 0, couplingHpm = 0;
    squareplus = 1 - pow(2*mHpm/m1,2);
    lambda = pow(squareplus,0.5);
    if (squareplus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecayHpHmNMSSM\n");
    } 
    couplingHpm = (lam*mueff/(root2)*(2*CPEMix(higgs,3)*pow(cos(beta),2) + 2*CPEMix(higgs,3)*pow(sin(beta),2)) - pow(lam,2)*mWboson*sin(beta)/g*2*CPEMix(higgs,2)*cos(beta)*sin(beta) - pow(lam,2)*mWboson*sin(beta)/g*2*CPEMix(higgs,1)*cos(beta)*sin(beta) + mueff*kappa*root2*2*CPEMix(higgs,3)*cos(beta)*sin(beta) + lam*Alambda/root2*2*CPEMix(higgs,3)*cos(beta)*sin(beta) + pow(gp,2)/4*mWboson/g*(sin(beta)*(2*CPEMix(higgs,1)*pow(cos(beta),2) - 2*CPEMix(higgs,1)*pow(sin(beta),2)) + cos(beta)*(2*CPEMix(higgs,2)*pow(sin(beta),2) - 2*CPEMix(higgs,2)*pow(cos(beta),2))) + g/4*mWboson*(sin(beta)*(2*CPEMix(higgs,1)*pow(cos(beta),2) + 2*CPEMix(higgs,1)*pow(sin(beta),2) + 2*2*CPEMix(higgs,2)*sin(beta)*cos(beta)) + cos(beta)*(2*CPEMix(higgs,2)*pow(cos(beta),2) + 2*CPEMix(higgs,2)*pow(sin(beta),2) + 2*2*CPEMix(higgs,1)*sin(beta)*cos(beta))) + lam/root2*0);

    amplitudeW = lambda*pow(couplingHpm,2)/(16*PI*m1);
  }
  return amplitudeW;
}

double hHH3Couplings(DoubleMatrix & CPEMix, int higgsno, int higgs1, int higgs2, int x, int y, int z) {
  double coupling = 0;
  coupling = CPEMix(higgsno,x)*CPEMix(higgs1,y)*CPEMix(higgs2,z) + CPEMix(higgsno,x)*CPEMix(higgs2,y)*CPEMix(higgs1,z) + CPEMix(higgs1,x)*CPEMix(higgsno,y)*CPEMix(higgs2,z) + CPEMix(higgs1,x)*CPEMix(higgsno,z)*CPEMix(higgs2,y) + CPEMix(higgs2,x)*CPEMix(higgsno,y)*CPEMix(higgs1,z) + CPEMix(higgs2,x)*CPEMix(higgs1,y)*CPEMix(higgsno,z);

  return coupling;
}

double higgsCPevenamplitudedecayhhorhHorHHNMSSM(double m1, double mh1, double mh2, double g, double gp, double mWboson, double mtop, double mbottom, double beta, double lam, double Alambda, double kappa, double Akappa, double mueff, DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, int higgs1, int higgs2, int starthiggs)
{
  double amplitudeW = 0;
  if(m1 < mh1 + mh2) {amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, coupling = 0, deltah1h2 = 0, hvev1 = 0, hvev2 = 0;
    squareplus = 1 - pow(mh1/m1+mh2/m1,2);
    squareminus = 1 - pow(mh1/m1-mh2/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecayhhorhHorHHNMSSM\n");
    } 
    if (higgs1 == higgs2) { deltah1h2 = 1;}
    else {deltah1h2 = 2;}
    prefactor = deltah1h2/(32*PI*m1);
    hvev1 = root2*mWboson*sin(beta)/g;
    hvev2 = root2*mWboson*cos(beta)/g;
    
    coupling = (pow(g,2)+pow(gp,2))/(4*root2)*(hvev1*(hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,1,1,1)-hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,1,2,2))+ hvev2*(hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,2,2,2) - hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,2,1,1))) - lam*Alambda/root2*hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,1,2,3) + kappa*Akappa/(3*root2)*hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,3,3,3) + pow(lam,2)/root2*(hvev1*(hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,1,2,2) + hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,1,3,3)) + hvev2*(hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,2,1,1) + hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,2,3,3)) + mueff/lam*(hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,3,1,1)+hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,3,2,2))) + pow(kappa,2)*root2*mueff/lam*hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,3,3,3) - lam*kappa/(root2)*(hvev1*hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,3,2,3) + hvev2*hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,3,1,3) + 2*mueff/lam*hHH3Couplings(CPEMix,starthiggs,higgs1,higgs2,1,2,3));

    amplitudeW = prefactor*lambda*pow(coupling,2);
  }
  return amplitudeW;
}



double higgsA2amplitudedecayA1CPevenNMSSM(double m1, double mA1, double mh, double mWboson, double runmt, double runmb, double g, double gp, double beta, DoubleMatrix & CPEMix, DoubleMatrix & CPOMix, double lam, double kappa, double Alambda, double Akappa, double mueff, int higgs) 
{
  double amplitudeW = 0;
  if (fabs(m1) < fabs(mA1) + fabs(mh)) {
    amplitudeW = 0;
  }
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, hvev1 = 0, hvev2 = 0, prefactor = 0, coupling = 0;
    squareplus = 1 - pow(mA1/m1 + mh/m1,2);
    squareminus = 1 - pow(mA1/m1 - mh/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsA2amplitudedecayA1CPevenNMSSM\n");
    } 
    
    prefactor = 1/(16*PI*m1);
    hvev1 = root2*mWboson*sin(beta)/g;
    hvev2 = root2*mWboson*cos(beta)/g;
    int pseudoscalar1 = 1, pseudoscalar2 = 2;

    coupling = (pow(g,2)+pow(gp,2))/(4*root2)*(hvev1*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,1,1) - CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,2,2)) + hvev2*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,2,2) - CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,1,1))) + lam*Alambda/(root2)*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,2,3)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,1,3)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,1,2)) - kappa*Akappa/(root2)*CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,3,3) + pow(lam,2)/root2*(hvev1*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,2,2)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,3,3)) + hvev2*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,1,1)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,3,3)) + mueff/lam*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,1,1)+CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,2,2))) + pow(kappa,2)*root2*mueff/lam*CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,3,3) + lam*kappa/(root2)*(hvev1*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,3,3)-2*CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,2,3)) + hvev2*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,3,3)-2*CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,1,3)) + 2*mueff/lam*(CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,3,1,2)-CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,1,2,3) - CPEhCPOACPOACoupling(CPEMix,CPOMix,higgs,pseudoscalar1,pseudoscalar2,2,1,3)));

    amplitudeW = prefactor*lambda*pow(coupling,2);
  }
  return amplitudeW;
}


double higgsCPevenamplitudedecayWHpmNMSSM (double m1, double mWboson, double mHpm, double beta, double g, DoubleMatrix & CPEMix, int higgs)
{
  double amplitudeW = 0;
  if (m1 < mWboson + mHpm) {amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0;
    squareplus = 1 - pow(mWboson/m1+mHpm/m1,2);
    squareminus = 1 - pow(mWboson/m1-mHpm/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecayWHpmNMSSM\n");
    } 
    
    amplitudeW = GFosqrt2*pow(m1,3)/(8*PI)*pow(CPEMix(higgs,1)*cos(beta)-CPEMix(higgs,2)*sin(beta),2)*pow(lambda,3);
  }
  return amplitudeW;
}


double higgsCPevenamplitudedecaystopistopiNMSSM (double m1, double mstopi, double thetat, double runmt, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double At, double mueff, double lam, int stop, int higgs) 
{
  double amplitudeW = 0;
  if (m1 < 2*mstopi) { amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, hvev1 = 0, hvev2 = 0, huq = 0, coupling = 0, angle1 = 0, angle2 = 0;
    squareplus = 1 - 4*pow(mstopi/m1,2);
    squareminus = 1;
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecaystopistopiNMSSM\n");
    } 
    hvev1 = (root2*mWboson*sin(beta))/g;
    hvev2 = (root2*mWboson*cos(beta))/g;
    huq = runmt/hvev1;
    if (stop == 1) {
      angle1 = cos(thetat); angle2 = sin(thetat);
    }
    else if (stop == 2) {
      angle1 = -sin(thetat); angle2 = cos(thetat);
    }
    else{
      throw("problem: stop must be 1 or 2 in higgsCPevenamplitudedecaystopistopiNMSSM\n");
    }
    coupling = pow(angle1,2)*root2*(pow(huq,2)*hvev1*CPEMix(higgs,1) + (pow(gp,2)/12 - pow(g,2)/4)*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2))) + pow(angle2,2)*root2*(pow(huq,2)*hvev1*CPEMix(higgs,1) - pow(gp,2)/3*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2))) + 2*angle1*angle2*huq/(root2)*(At*CPEMix(higgs,1) - mueff*CPEMix(higgs,2) - lam*hvev2*CPEMix(higgs,3));

    amplitudeW = 3/(16*PI*m1)*lambda*pow(coupling,2);
  }
  return amplitudeW;
}


double higgsCPevenamplitudedecaystopistopjNMSSM (double m1, double mstopi, double mstopj, double thetat, double runmt, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double At, double mueff, double lam, int higgs) 
{
  double amplitudeW = 0;
  if (m1 < mstopi+mstopj) { amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, hvev1 = 0, hvev2 = 0, huq = 0, coupling = 0;
    squareplus = 1 - pow(mstopi/m1+mstopj/m1,2);
    squareminus = 1- pow(mstopi/m1-mstopj/m1,2); 
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecaystopistopjNMSSM\n");
    } 
    hvev1 = (root2*mWboson*sin(beta))/g;
    hvev2 = (root2*mWboson*cos(beta))/g;
    huq = runmt/hvev1;

    coupling = cos(thetat)*sin(thetat)*(root2*(pow(huq,2)*hvev1*CPEMix(higgs,1) - pow(gp,2)/3*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2))) - root2*(pow(huq,2)*hvev1*CPEMix(higgs,1) + (pow(gp,2)/12 - pow(g,2)/4)*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2)))) + (pow(cos(thetat),2)-pow(sin(thetat),2))*huq/(root2)*(At*CPEMix(higgs,1) - mueff*CPEMix(higgs,2) - lam*hvev2*CPEMix(higgs,3));

    amplitudeW = 3/(16*PI*m1)*lambda*pow(coupling,2);
  }
  return amplitudeW;
}


double higgsCPevenamplitudedecaysbottomisbottomiNMSSM (double m1, double msbottomi, double thetab, double runmb, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double Ab, double mueff, double lam, int sbottom, int higgs) 
{
  double amplitudeW = 0;
  if (m1 < 2*msbottomi) { amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, hvev1 = 0, hvev2 = 0, hdq = 0, coupling = 0, angle1 = 0, angle2 = 0;
    squareplus = 1 - 4*pow(msbottomi/m1,2);
    squareminus = 1;
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecaysbottomisbottomiNMSSM\n");
    } 
    hvev1 = (root2*mWboson*sin(beta))/g;
    hvev2 = (root2*mWboson*cos(beta))/g;
    hdq = runmb/hvev2;
    if (sbottom == 1) {
      angle1 = cos(thetab); angle2 = sin(thetab);
    }
    else if (sbottom == 2) {
      angle1 = -sin(thetab); angle2 = cos(thetab);
    }
    else{
      throw("problem: sbottom must be 1 or 2 in higgsCPevenamplitudedecaysbottomisbottomiNMSSM");
    }

    coupling = pow(angle1,2)*root2*(pow(hdq,2)*hvev2*CPEMix(higgs,2) + (pow(gp,2)/12 + pow(g,2)/4)*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2))) + pow(angle2,2)*root2*(pow(hdq,2)*hvev2*CPEMix(higgs,2) + pow(gp,2)/6*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2))) + 2*angle1*angle2*hdq/(root2)*(Ab*CPEMix(higgs,2) - mueff*CPEMix(higgs,1) - lam*hvev1*CPEMix(higgs,3));

    amplitudeW = 3/(16*PI*m1)*lambda*pow(coupling,2);
  }
  return amplitudeW;
}


double higgsCPevenamplitudedecaysbottomisbottomjNMSSM (double m1, double msbottomi, double msbottomj, double thetab, double runmb, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double Ab, double mueff, double lam, int higgs) 
{
  double amplitudeW = 0;
  if (m1 < msbottomi+msbottomj) { amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, hvev1 = 0, hvev2 = 0, hdq = 0, coupling = 0;
    squareplus = 1 - pow(msbottomi/m1+msbottomj/m1,2);
    squareminus = 1- pow(msbottomi/m1-msbottomj/m1,2); 
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecaysbottomisbottomjNMSSM\n");
    } 
    hvev1 = (root2*mWboson*sin(beta))/g;
    hvev2 = (root2*mWboson*cos(beta))/g;
    hdq = runmb/hvev2;

    coupling = cos(thetab)*sin(thetab)*(root2*(pow(hdq,2)*hvev2*CPEMix(higgs,2) + pow(gp,2)/6*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2))) - root2*(pow(hdq,2)*hvev2*CPEMix(higgs,2) + (pow(gp,2)/12 + pow(g,2)/4)*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2)))) + (pow(cos(thetab),2)-pow(sin(thetab),2))*hdq/(root2)*(Ab*CPEMix(higgs,2) - mueff*CPEMix(higgs,1) - lam*hvev1*CPEMix(higgs,3));

    amplitudeW = 3/(16*PI*m1)*lambda*pow(coupling,2);
  }
  return amplitudeW;
}


double higgsCPevenamplitudedecaystauistauiNMSSM (double m1, double mstaui, double thetatau, double runmtau, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double Atau, double mueff, double lam, int stau, int higgs) 
{
  double amplitudeW = 0;
  if (m1 < 2*mstaui) { amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, hvev1 = 0, hvev2 = 0, hlq = 0, coupling = 0, angle1 = 0, angle2 = 0;
    squareplus = 1 - 4*pow(mstaui/m1,2);
    squareminus = 1;
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecaystauistauiNMSSM\n");
    } 
    hvev1 = (root2*mWboson*sin(beta))/g;
    hvev2 = (root2*mWboson*cos(beta))/g;
    hlq = runmtau/hvev2;
    if (stau == 1) {
      angle1 = sin(thetatau); angle2 = cos(thetatau);
    }
    else if (stau == 2) {
      angle1 = -cos(thetatau); angle2 = sin(thetatau);
    }
    else{
      throw("problem: stau must be 1 or 2 in higgsCPevenamplitudedecaystauistauiNMSSM");
    }

    coupling = pow(angle1,2)*root2*(pow(hlq,2)*hvev2*CPEMix(higgs,2) + (-pow(gp,2)/4 + pow(g,2)/4)*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2))) + pow(angle2,2)*root2*(pow(hlq,2)*hvev2*CPEMix(higgs,2) + pow(gp,2)/2*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2))) - 2*angle1*angle2*hlq/(root2)*(Atau*CPEMix(higgs,2) - mueff*CPEMix(higgs,1) - lam*hvev1*CPEMix(higgs,3));

    amplitudeW = 1/(16*PI*m1)*lambda*pow(coupling,2);
  }
  return amplitudeW;
}


double higgsCPevenamplitudedecaystauistaujNMSSM (double m1, double mstaui, double mstauj, double thetatau, double runmtau, double g, double gp, double mWboson, double beta, DoubleMatrix & CPEMix, double Atau, double mueff, double lam, int higgs) 
{
  double amplitudeW = 0;
  if (m1 < mstaui+mstauj) { amplitudeW = 0;}
  else {
    double squareplus = 0, squareminus = 0, lambda = 0, hvev1 = 0, hvev2 = 0, hlq = 0, coupling = 0;
    squareplus = 1 - pow(mstaui/m1+mstauj/m1,2);
    squareminus = 1- pow(mstaui/m1-mstauj/m1,2); 
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in higgsCPevenamplitudedecaystauistaujNMSSM\n");
    } 
    hvev1 = (root2*mWboson*sin(beta))/g;
    hvev2 = (root2*mWboson*cos(beta))/g;
    hlq = runmtau/hvev2;

    coupling = cos(thetatau)*sin(thetatau)*(root2*(pow(hlq,2)*hvev2*CPEMix(higgs,2) + pow(gp,2)/2*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2))) - root2*(pow(hlq,2)*hvev2*CPEMix(higgs,2) + (-pow(gp,2)/4 + pow(g,2)/4)*(hvev1*CPEMix(higgs,1) - hvev2*CPEMix(higgs,2)))) - (pow(sin(thetatau),2)-pow(cos(thetatau),2))*hlq/(root2)*(Atau*CPEMix(higgs,2) - mueff*CPEMix(higgs,1) - lam*hvev1*CPEMix(higgs,3));

    amplitudeW = 1/(16*PI*m1)*lambda*pow(coupling,2);
  }
  return amplitudeW;
}


double stop2amplitudedecaystop1CPevenhiggsNMSSM (double mst2, double mst1, double mh, double mt , double thetat, DoubleMatrix & CPEMix, double beta, double mWboson, double g, double gp, double At, double mueff, double lam, int higgs) 
{
  double amplitudeW = 0, coupling = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, ft = 0, hvev1 = 0, hvev2 = 0, CL = 0, CR = 0, CLR = 0;
  if(mst2 < mst1 + mh) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(16*PI*mst2);
    squareplus = 1 - pow(mst1/mst2 + mh/mst2,2);
    squareminus = 1 - pow(mst1/mst2 - mh/mst2,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stop2amplitudedecaystop1CPevenhiggsNMSSM\n");
    } 
    ft = g*mt/(root2*mWboson*sin(beta));
    hvev1 = root2*mWboson*sin(beta)/g;
    hvev2 = root2*mWboson*cos(beta)/g;
    CL = -root2*(pow(ft,2)*hvev1*CPEMix(higgs,1) + (pow(gp,2)/12-pow(g,2)/4)*(hvev1*CPEMix(higgs,1)-hvev2*CPEMix(higgs,2)));
    CR = -root2*(pow(ft,2)*hvev1*CPEMix(higgs,1) - pow(gp,2)/3*(hvev1*CPEMix(higgs,1)-hvev2*CPEMix(higgs,2)));
    CLR = -ft/root2*(At*CPEMix(higgs,1)-mueff*CPEMix(higgs,2)-lam*hvev2*CPEMix(higgs,3));
    coupling = cos(thetat)*sin(thetat)*(CR-CL) + (pow(cos(thetat),2)-pow(sin(thetat),2))*CLR;
    amplitudeW = prefactor*lambda*pow(coupling,2);
  }
  return amplitudeW;
}
    

double stop2amplitudedecaystop1CPoddhiggsNMSSM (double mst2, double mst1, double ma, double mt , double thetat, DoubleMatrix & CPOMix, double beta, double mWboson, double g , double At, double mueff, double lam, int higgs) 
{
  double amplitudeW = 0, coupling = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, ft = 0, hvev2 = 0, ALR = 0;
  if(mst2 < mst1 + ma) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(16*PI*mst2);
    squareplus = 1 - pow(mst1/mst2 + ma/mst2,2);
    squareminus = 1 - pow(mst1/mst2 - ma/mst2,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stop2amplitudedecaystop1CPoddhiggsNMSSM\n");
    } 
    ft = g*mt/(root2*mWboson*sin(beta));
    hvev2 = root2*mWboson*cos(beta)/g;
    ALR = -(ft/root2)*(At*CPOMix(higgs,1) + mueff*CPOMix(higgs,2) + lam*hvev2*CPOMix(higgs,3));
    coupling = (pow(cos(thetat),2)-pow(sin(thetat),2))*ALR;
    amplitudeW = prefactor*lambda*pow(coupling,2);

  }
  return amplitudeW;
}


double sbottom2amplitudedecaysbottom1CPevenhiggsNMSSM (double msb2, double msb1, double mh, double mb , double thetab, DoubleMatrix & CPEMix, double beta, double mWboson, double g, double gp, double Ab, double mueff, double lam, int higgs) 
{
  double amplitudeW = 0, coupling = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, fb = 0, hvev1 = 0, hvev2 = 0, CL = 0, CR = 0, CLR = 0;
  if(msb2 < msb1 + mh) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(16*PI*msb2);
    squareplus = 1 - pow(msb1/msb2 + mh/msb2,2);
    squareminus = 1 - pow(msb1/msb2 - mh/msb2,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in sbottom2amplitudedecaysbottom1CPevenhiggsNMSSM\n");
    } 
    fb = g*mb/(root2*mWboson*cos(beta));
    hvev1 = root2*mWboson*sin(beta)/g;
    hvev2 = root2*mWboson*cos(beta)/g;
    CL = -root2*(pow(fb,2)*hvev2*CPEMix(higgs,2) + (pow(gp,2)/12+pow(g,2)/4)*(hvev1*CPEMix(higgs,1)-hvev2*CPEMix(higgs,2)));
    CR = -root2*(pow(fb,2)*hvev2*CPEMix(higgs,2) + pow(gp,2)/6*(hvev1*CPEMix(higgs,1)-hvev2*CPEMix(higgs,2)));
    CLR = -fb/root2*(Ab*CPEMix(higgs,2)-mueff*CPEMix(higgs,1)-lam*hvev1*CPEMix(higgs,3));
    coupling = cos(thetab)*sin(thetab)*(CR-CL) + (pow(cos(thetab),2)-pow(sin(thetab),2))*CLR;
    amplitudeW = prefactor*lambda*pow(coupling,2);

  }
  return amplitudeW;
}


double sbottom2amplitudedecaysbottom1CPoddhiggsNMSSM (double msb2, double msb1, double ma, double mb , double thetab, DoubleMatrix & CPOMix, double beta, double mWboson, double g , double Ab, double mueff, double lam, int higgs) 
{
  double amplitudeW = 0, coupling = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, fb = 0, hvev1 = 0, ALR = 0;
  if(msb2 < msb1 + ma) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(16*PI*msb2);
    squareplus = 1 - pow(msb1/msb2 + ma/msb2,2);
    squareminus = 1 - pow(msb1/msb2 - ma/msb2,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in sbottom2amplitudedecaysbottom1CPoddhiggsNMSSM\n");
    } 
    fb = g*mb/(root2*mWboson*cos(beta));
    hvev1 = root2*mWboson*sin(beta)/g;
    ALR = -(fb/root2)*(Ab*CPOMix(higgs,2) + mueff*CPOMix(higgs,1) + lam*hvev1*CPOMix(higgs,3));
    coupling = (pow(cos(thetab),2)-pow(sin(thetab),2))*ALR;
    amplitudeW = prefactor*lambda*pow(coupling,2);

  }
  return amplitudeW;
}


double stau2amplitudedecaystau1CPevenhiggsNMSSM (double mstau2, double mstau1, double mh, double mtau, double thetatau, DoubleMatrix & CPEMix, double beta, double mWboson, double g, double gp, double Atau, double mueff, double lam, int higgs) 
{
  double amplitudeW = 0, coupling = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, ftau = 0, hvev1 = 0, hvev2 = 0, CL = 0, CR = 0, CLR = 0;
  if(mstau2 < mstau1 + mh) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(16*PI*mstau2);
    squareplus = 1 - pow(mstau1/mstau2 + mh/mstau2,2);
    squareminus = 1 - pow(mstau1/mstau2 - mh/mstau2,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stau2amplitudedecaystau1CPevenhiggsNMSSM\n");
    } 
    ftau = g*mtau/(root2*mWboson*cos(beta));
    hvev1 = root2*mWboson*sin(beta)/g;
    hvev2 = root2*mWboson*cos(beta)/g;
    CL = -root2*(pow(ftau,2)*hvev2*CPEMix(higgs,2) + (-pow(gp,2)/4+pow(g,2)/4)*(hvev1*CPEMix(higgs,1)-hvev2*CPEMix(higgs,2)));
    CR = -root2*(pow(ftau,2)*hvev2*CPEMix(higgs,2) + pow(gp,2)/2*(hvev1*CPEMix(higgs,1)-hvev2*CPEMix(higgs,2)));
    CLR = -ftau/root2*(Atau*CPEMix(higgs,2)-mueff*CPEMix(higgs,1)-lam*hvev1*CPEMix(higgs,3));
    coupling = cos(thetatau)*sin(thetatau)*(CR-CL) + (pow(cos(thetatau),2)-pow(sin(thetatau),2))*CLR;
    amplitudeW = prefactor*lambda*pow(coupling,2);
  }
  return amplitudeW;
}


double stau2amplitudedecaystau1CPoddhiggsNMSSM
(double mstau2, double mstau1, double ma, double mtau, double thetatau,
 DoubleMatrix & CPOMix, double beta, double mWboson, double g , double Atau,
 double mueff, double lam, int higgs) {
  double amplitudeW = 0, coupling = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, ftau = 0, hvev1 = 0, ALR = 0;
  if(mstau2 < mstau1 + ma) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(16*PI*mstau2);
    squareplus = 1 - pow(mstau1/mstau2 + ma/mstau2,2);
    squareminus = 1 - pow(mstau1/mstau2 - ma/mstau2,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stau2amplitudedecaystau1CPoddhiggsNMSSM\n");
    } 
    ftau = g*mtau/(root2*mWboson*cos(beta));
    hvev1 = root2*mWboson*sin(beta)/g;
    ALR = -(ftau/root2)*(Atau*CPOMix(higgs,2) + mueff*CPOMix(higgs,1) + lam*hvev1*CPOMix(higgs,3));
    coupling = (pow(cos(thetatau),2)-pow(sin(thetatau),2))*ALR;
    amplitudeW = prefactor*lambda*pow(coupling,2);

  }
  return amplitudeW;
}


double chargino2amplitudedecaychargino1CPevenhiggsNMSSM (double mchar2, double mchar1, double mh, double g, double lam, double thetaL, double thetaR, DoubleMatrix & CPEMix, int higgs) 
{
  double amplitudeW = 0, lambda = 0, squareplus = 0, squareminus = 0, coupling1 = 0, coupling2 = 0, coupling = 0, prefactor = 0;
  if (fabs(mchar2) < fabs(mchar1) + mh) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(32*PI*fabs(mchar2));
    coupling1 = (lam/(root2)*CPEMix(higgs,3)*cos(thetaL)*sin(thetaR) + g/(root2)*(CPEMix(higgs,1)*sin(thetaL)*sin(thetaR) - CPEMix(higgs,2)*cos(thetaL)*cos(thetaR)));
    coupling2 = (lam/(root2)*CPEMix(higgs,3)*sin(thetaL)*cos(thetaR) - g/(root2)*(CPEMix(higgs,1)*cos(thetaL)*cos(thetaR) - CPEMix(higgs,2)*sin(thetaL)*sin(thetaR)));
    squareplus = 1 - pow(mchar1/mchar2 + mh/mchar2,2);
    squareminus = 1 - pow(mchar1/mchar2 - mh/mchar2,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in chargino2amplitudedecaychargino1CPevenhiggsNMSSM\n");
    } 
    coupling = (pow(coupling1,2)+pow(coupling2,2))*(pow(mchar1,2)+pow(mchar2,2)-pow(mh,2)) + 4*coupling1*coupling2*mchar1*mchar2;
    amplitudeW = prefactor*lambda*coupling;

  }
  return amplitudeW;
}

 
double chargino2amplitudedecaychargino1CPoddhiggsNMSSM (double mchar2, double mchar1, double mA, double g, double lam, double thetaL, double thetaR, DoubleMatrix & CPOMix, int higgs) 
{
  double amplitudeW = 0, lambda = 0, squareplus = 0, squareminus = 0, C1 = 0, C2 = 0, coupling = 0, prefactor = 0;
  if (fabs(mchar2) < fabs(mchar1) + mA) {
    amplitudeW = 0;
  }
  else {
    int pseudoscalar = higgs;
    prefactor = 1/(32*PI*fabs(mchar2));
    C1 = (lam/root2*CPOMix(pseudoscalar,3)*cos(thetaL)*sin(thetaR) - pow(2,-0.5)*g*(CPOMix(pseudoscalar,1)*sin(thetaL)*sin(thetaR) - CPOMix(pseudoscalar,2)*cos(thetaL)*cos(thetaR)));
    C2 = (lam/root2*CPOMix(pseudoscalar,3)*sin(thetaL)*cos(thetaR) - pow(2,-0.5)*g*(CPOMix(pseudoscalar,1)*-cos(thetaL)*cos(thetaR) + CPOMix(pseudoscalar,2)*sin(thetaR)*sin(thetaL)));
    squareplus = 1 - pow(mchar1/mchar2 + mA/mchar2,2);
    squareminus = 1 - pow(mchar1/mchar2 - mA/mchar2,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in chargino2amplitudedecaychargino1CPevoddhiggsNMSSM\n");
    } 
    coupling = (pow(C1,2)+pow(C2,2))*(pow(mchar1,2)+pow(mchar2,2)-pow(mA,2)) + 4*C1*C2*mchar1*mchar2;
    amplitudeW = prefactor*lambda*coupling;
  }
  return amplitudeW;
} 


double neutralinoamplitudedecaycharginoWNMSSM (double mneut, double mchar, double mWboson, double g, double thetaL, double thetaR, DoubleMatrix & mixNeut, int neutralino, int chargino) 
{
  double amplitudeW = 0, squareplus = 0, squareminus = 0, lambda = 0, coupleL = 0, coupleR = 0, coupling = 0, prefactor = 0, V2 = 0, V1 = 0, U2 = 0, U1 = 0;
  if (fabs(mneut) < fabs(mchar) + mWboson) {
    amplitudeW = 0;
  }
  else {
    prefactor = pow(g,2)/(32*PI*fabs(mneut));
    
    if (chargino == 1) {
      V2 = cos(thetaR); V1 = sin(thetaR); U2 = cos(thetaL); U1 = sin(thetaL);
    }
    else if (chargino == 2) {
      V2 = sin(thetaR); V1 = -cos(thetaR); U2 = sin(thetaL); U1 = -cos(thetaL);
    }
    else{
      throw("problem: chargino must be 1 or 2 in neutralinoamplitudedecaycharginoWNMSSM");
    }

    coupleL = -1/root2*mixNeut(neutralino,4)*V2 + mixNeut(neutralino,2)*V1;
    coupleR = 1/root2*mixNeut(neutralino,3)*U2 + mixNeut(neutralino,2)*U1;

    squareplus = 1 - pow(fabs(mchar)/fabs(mneut) + mWboson/fabs(mneut),2);
    squareminus = 1 - pow(fabs(mchar)/fabs(mneut) - mWboson/fabs(mneut),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecaycharginoWNMSSM\n");
    } 

    coupling = -12*mneut*mchar*coupleL*coupleR + (pow(coupleL,2)+pow(coupleR,2))*((pow(mchar,2)+pow(mneut,2)-pow(mWboson,2)) + (pow(mneut,2)+pow(mWboson,2)-pow(mchar,2))*(pow(mneut,2)-pow(mchar,2)-pow(mWboson,2))/pow(mWboson,2));

    amplitudeW = prefactor*lambda*coupling;
  }

  return amplitudeW;
}
    

double neutralinoamplitudedecayneutralinoZNMSSM (double mneuti, double mneutj, double mZboson, double g, double gp, DoubleMatrix & mixNeut, int neutralinoi, int neutralinoj)
{
  double amplitudeW = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, costhetaW = 0, coupleZL = 0, coupleZR = 0, coupling = 0;
  if (fabs(mneuti) < fabs(mneutj) + mZboson) { 
    amplitudeW = 0;
  }
  else {
    prefactor = pow(g,2)/(32*PI*fabs(mneuti));
    squareplus = 1 - pow(fabs(mneutj)/fabs(mneuti) + mZboson/fabs(mneuti),2);
    squareminus = 1 - pow(fabs(mneutj)/fabs(mneuti) - mZboson/fabs(mneuti),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecayneutralinoZNMSSM\n");
    } 

    costhetaW = g/pow(pow(g,2)+pow(gp,2),0.5);

    coupleZR = 1/(2*costhetaW)*(mixNeut(neutralinoi,3)*mixNeut(neutralinoj,3) - mixNeut(neutralinoi,4)*mixNeut(neutralinoj,4));
    coupleZL = -coupleZR;

    coupling = -12*mneuti*mneutj*coupleZL*coupleZR + (pow(coupleZL,2)+pow(coupleZR,2))*((pow(mneuti,2)+pow(mneutj,2)-pow(mZboson,2)) + (pow(mneuti,2)-pow(mneutj,2)+pow(mZboson,2))*(pow(mneuti,2)-pow(mneutj,2)-pow(mZboson,2))/pow(mZboson,2));

    amplitudeW = prefactor*lambda*coupling;
  }
  return amplitudeW;
}


double neutralinoamplitudecharginoHpmNMSSM (double mneut, double mchar, double mHp, double g, double gp, double beta, double thetaL, double thetaR, double lam, DoubleMatrix & mixNeut, int neutralino, int chargino)
{
  double amplitudeW = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, coupHpZiWjL = 0, coupHpZiWjR = 0, coupling = 0, V1 = 0, V2 = 0, U1 = 0, U2 = 0;
  if (fabs(mneut) < fabs(mchar) + mHp) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(32*PI*fabs(mneut));
    squareplus = 1 - pow(fabs(mchar)/fabs(mneut) + mHp/fabs(mneut),2);
    squareminus = 1 - pow(fabs(mchar)/fabs(mneut) - mHp/fabs(mneut),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecaycharginoHpmNMSSM\n");
    } 
    
    if (chargino == 1) {
      V2 = cos(thetaR); V1 = sin(thetaR); U2 = cos(thetaL); U1 = sin(thetaL);
    }
    else if (chargino == 2) {
      V2 = sin(thetaR); V1 = -cos(thetaR); U2 = sin(thetaL); U1 = -cos(thetaL);
    }
    else{
      throw("problem: chargino must be 1 or 2 in neutralinoamplitudecharginoHpmNMSSM\n");
    }

    coupHpZiWjL = lam*cos(beta)*mixNeut(neutralino,5)*U2 - sin(beta)/(root2)*(gp*mixNeut(neutralino,1)+g*mixNeut(neutralino,2))*U2 + sin(beta)*g*mixNeut(neutralino,3)*U1;
    coupHpZiWjR = lam*sin(beta)*mixNeut(neutralino,5)*V2 + cos(beta)/(root2)*(gp*mixNeut(neutralino,1)+g*mixNeut(neutralino,2))*V2 + cos(beta)*g*mixNeut(neutralino,4)*V1;

    coupling = (pow(coupHpZiWjL,2)+pow(coupHpZiWjR,2))*(pow(mchar,2)+pow(mneut,2)-pow(mHp,2)) + 4*coupHpZiWjL*coupHpZiWjR*mneut*mchar;

    amplitudeW = prefactor*coupling*lambda;
  }
  return amplitudeW;
}


double neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM (double mneuti, double mneutj, double mhiggs, double g, double gp, double lam, double kappa, DoubleMatrix & mixNeut, DoubleMatrix & CPEMix, int neutralinoi, int neutralinoj, int higgs) 
{
  double amplitudeW = 0, prefactor = 0, squareplus = 0, squareminus = 0, lambda = 0, coupling = 0, coupHZiZj = 0;
  if (fabs(mneuti) < fabs(mneutj) + mhiggs) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(8*PI*fabs(mneuti));
    squareplus = 1 - pow(fabs(mneutj)/fabs(mneuti) + mhiggs/fabs(mneuti),2);
    squareminus = 1 - pow(fabs(mneutj)/fabs(mneuti) - mhiggs/fabs(mneuti),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecayneutralinoCPevenhiggsNMSSM\n");
    } 

    coupHZiZj = 0.5*(lam/(root2)*(CPEMix(higgs,1)*(mixNeut(neutralinoj,3)*mixNeut(neutralinoi,5) + mixNeut(neutralinoj,5)*mixNeut(neutralinoi,3)) + CPEMix(higgs,2)*(mixNeut(neutralinoj,4)*mixNeut(neutralinoi,5)+mixNeut(neutralinoi,4)*mixNeut(neutralinoj,5)) + CPEMix(higgs,3)*(mixNeut(neutralinoj,4)*mixNeut(neutralinoi,3) + mixNeut(neutralinoj,3)*mixNeut(neutralinoi,4))) - root2*kappa*CPEMix(higgs,3)*mixNeut(neutralinoj,5)*mixNeut(neutralinoi,5) + gp/2*(-CPEMix(higgs,1)*(mixNeut(neutralinoj,1)*mixNeut(neutralinoi,4)+mixNeut(neutralinoj,4)*mixNeut(neutralinoi,1)) + CPEMix(higgs,2)*(mixNeut(neutralinoj,1)*mixNeut(neutralinoi,3)+mixNeut(neutralinoi,1)*mixNeut(neutralinoj,3))) + g/2*(CPEMix(higgs,1)*(mixNeut(neutralinoj,2)*mixNeut(neutralinoi,4)+mixNeut(neutralinoj,4)*mixNeut(neutralinoi,2)) - CPEMix(higgs,2)*(mixNeut(neutralinoj,2)*mixNeut(neutralinoi,3)+mixNeut(neutralinoj,3)*mixNeut(neutralinoi,2))));

    coupling = 2*pow(coupHZiZj,2)*(pow(mneuti,2)+pow(mneutj,2)-pow(mhiggs,2)) + 4*pow(coupHZiZj,2)*mneuti*mneutj;

    amplitudeW = prefactor*coupling*lambda;
  }
  return amplitudeW;
}
    

double neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM (double mneuti, double mneutj, double ma, double g, double gp, double lam, double kappa, DoubleMatrix & mixNeut, DoubleMatrix & CPOMix, int neuti, int neutj, int higgsa) 
{
  double amplitudeW = 0, lambda = 0, squareplus = 0, squareminus = 0, prefactor = 0, coupling = 0, coupAZiZj = 0;
  if (fabs(mneuti) < fabs(mneutj) + ma) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(4*PI*fabs(mneuti));
    squareplus = 1 - pow(fabs(mneutj)/fabs(mneuti) + ma/fabs(mneuti),2);
    squareminus = 1 - pow(fabs(mneutj)/fabs(mneuti) - ma/fabs(mneuti),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecayneutralinoCPoddhiggsNMSSM\n");
    } 

    coupAZiZj = -0.5*(lam/root2*(CPOMix(higgsa,1)*(mixNeut(neutj,3)*mixNeut(neuti,5) + mixNeut(neuti,3)*mixNeut(neutj,5)) + CPOMix(higgsa,2)*(mixNeut(neutj,4)*mixNeut(neuti,5) + mixNeut(neuti,4)*mixNeut(neutj,5)) + CPOMix(higgsa,3)*(mixNeut(neutj,4)*mixNeut(neuti,3) + mixNeut(neutj,3)*mixNeut(neuti,4))) - root2*kappa*CPOMix(higgsa,3)*mixNeut(neuti,5)*mixNeut(neutj,5) - gp/2*(-CPOMix(higgsa,1)*(mixNeut(neutj,1)*mixNeut(neuti,4) + mixNeut(neuti,1)*mixNeut(neutj,4)) + CPOMix(higgsa,2)*(mixNeut(neutj,1)*mixNeut(neuti,3) + mixNeut(neuti,1)*mixNeut(neutj,3))) - g/2*(CPOMix(higgsa,1)*(mixNeut(neutj,2)*mixNeut(neuti,4) + mixNeut(neutj,4)*mixNeut(neuti,2)) - CPOMix(higgsa,2)*(mixNeut(neutj,2)*mixNeut(neuti,3)+mixNeut(neutj,3)*mixNeut(neuti,2))));

    coupling = pow(coupAZiZj,2)*(pow(mneuti,2)+pow(mneutj,2) - pow(ma,2)) - 2*pow(coupAZiZj,2)*mneuti*mneutj;

    amplitudeW = coupling*prefactor*lambda;
  }
  return amplitudeW;
}
    

double neutralinoamplitudedecaysfermionfermionfirst2genNMSSM (double mneut, double msf, double mf, double g, double gp, DoubleMatrix & mixNeut, int neut, char type, char LorR) ///type indicates type of fermion, 'u' for up type quark, 'd' for down type quark, 'l' for charged lepton, 'n' for neutrino, note ignoring yukawa of fermion here as first 2 gen
{
  double amplitudeW = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, costhetaW = 0, sinthetaW = 0, coup1 = 0, coup2 = 0, coupling = 0, q = 0, N = 0;
  DoubleVector c(5);
  for (int i=1;i<=5;i++) { c(i) = 0;}
  if (fabs(mneut) < mf + msf) {
    amplitudeW = 0;
  }
  else {
    squareplus = 1 - pow(msf/fabs(mneut) + mf/fabs(mneut),2);
    squareminus = 1 - pow(msf/fabs(mneut) - mf/fabs(mneut),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudedecaysfermionfermionfirst2genNMSSM\n");
    } 
    costhetaW = g/pow(pow(gp,2)+pow(g,2),0.5);
    sinthetaW = gp/pow(pow(gp,2)+pow(g,2),0.5);

    c(1) = mixNeut(neut,1)*costhetaW + mixNeut(neut,2)*sinthetaW;
    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    c(3) = mixNeut(neut,3);
    c(4) = mixNeut(neut,4);
    c(5) = mixNeut(neut,5);

    if (type == 'u') {
      N = 3;
      q = 2./3;
      coup1 = -root2*(q*c(1)*sinthetaW + (0.5-q*pow(sinthetaW,2))*c(2)/costhetaW);
      coup2 = -q*root2*sinthetaW*(c(2)*gp/g - c(1));
    }
    else if (type == 'd') {
      N = 3;
      q = -1./3;
      coup1 = root2*(c(1)*sinthetaW*-q + (0.5+q*pow(sinthetaW,2))*c(2)/costhetaW);
      coup2 = -q*root2*sinthetaW*(c(2)*gp/g - c(1));
    }
    else if (type == 'l') {
      N = 1;
      q = -1;
      coup1 = root2*(c(1)*sinthetaW*-q + (0.5-pow(sinthetaW,2))*c(2)/costhetaW);
      coup2 = -q*root2*sinthetaW*(c(2)*gp/g-c(1));
    }
    else if (type == 'n') {
      N = 1;
      q = 0;
      coup1 = -c(2)/(root2*costhetaW);
      coup2 = 0;
    }
    else{
      throw("problem: type must be u or d or l or n in neutralinoamplitudedecaysfermionsfermionfirst2genNMSSM\n");
    }
    
    if (LorR == 'L') {
      coupling = coup1;
    }
    else if (LorR == 'R') {
      coupling = coup2;
    }
    else{
      throw("problem: LorR must be L or R in neutralinoamplitudedecaysfermionsfermionfirst2genNMSSM\n");
    }
    prefactor = N*pow(g,2)/(32*PI*fabs(mneut));
    amplitudeW = prefactor*lambda*(pow(coupling,2))*(pow(mneut,2)-pow(msf,2)+pow(mf,2));
  }
  return amplitudeW;
}
      
  
double neutralinoamplitudestoptopNMSSM (double mneut, double mst, double mt, double g, double gp, double thetat, double beta, double mWboson, DoubleMatrix & mixNeut, double runmt, int neut, int stop) 
{
  double amplitudeW = 0, squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, ft = 0, coup1 = 0, coup2 = 0, sinthetaW = 0, costhetaW = 0;
  DoubleVector c(5);
  for (int i=1;i<=5;i++) { c(i) = 0;}
  if (fabs(mneut) < mst + mt) {
    amplitudeW = 0;
  }
  else {
    prefactor = 3*pow(g,2)/(32*PI*fabs(mneut));
    squareplus = 1 - pow(mst/fabs(mneut) + mt/fabs(mneut),2);
    squareminus = 1 - pow(mst/fabs(mneut) - mt/fabs(mneut),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudestoptopNMSSM\n");
    } 
    costhetaW = g/pow(pow(gp,2)+pow(g,2),0.5);
    sinthetaW = gp/pow(pow(gp,2)+pow(g,2),0.5);

    c(1) = mixNeut(neut,1)*costhetaW + mixNeut(neut,2)*sinthetaW;
    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    c(3) = mixNeut(neut,3);
    c(4) = mixNeut(neut,4);
    c(5) = mixNeut(neut,5);
    ft = g*runmt/(root2*mWboson*sin(beta));

    if ( stop == 1) {
      coup1 = cos(thetat)*root2*(-2./3*c(1)*sinthetaW + (-0.5 + 2./3*pow(sinthetaW,2))*c(2)/costhetaW) - sin(thetat)*ft/g*mixNeut(neut,4);
      coup2 = -2./3*sin(thetat)*root2*sinthetaW*(c(2)*gp/g-c(1)) - cos(thetat)*ft/g*mixNeut(neut,4);
    }
    else if (stop == 2) {
      coup1 = -sin(thetat)*root2*(-2./3*c(1)*sinthetaW + (-0.5+2./3*pow(sinthetaW,2))*c(2)/costhetaW) - cos(thetat)*ft/g*mixNeut(neut,4);
      coup2 = -2./3*cos(thetat)*root2*sinthetaW*(c(2)*gp/g - c(1)) + sin(thetat)*ft/g*mixNeut(neut,4);
    }
    else{
      throw("problem: stop must be 1 or 2 in neutralinoamplitudestoptopNMSSM");
    }
    amplitudeW = prefactor*lambda*((pow(coup1,2)+pow(coup2,2))*(pow(mneut,2)-pow(mst,2)+pow(mt,2)) + 4*mt*mneut*coup1*coup2);
  }
  return amplitudeW;
}

double neutralinoamplitudesbottombottomNMSSM (double mneut, double msb, double mb, double g, double gp, double thetab, double beta, double mWboson, DoubleMatrix & mixNeut, double runmb, int neut, int sbottom) 
{
  double amplitudeW = 0, squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, fb = 0, coup1 = 0, coup2 = 0, sinthetaW = 0, costhetaW = 0;
  DoubleVector c(5);
  for (int i=1;i<=5;i++) { c(i) = 0;}
  if (fabs(mneut) < msb + mb) {
    amplitudeW = 0;
  }
  else {
    prefactor = 3*pow(g,2)/(32*PI*fabs(mneut));
    squareplus = 1 - pow(msb/fabs(mneut) + mb/fabs(mneut),2);
    squareminus = 1 - pow(msb/fabs(mneut) - mb/fabs(mneut),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudesbottombottomNMSSM\n");
    } 
    costhetaW = g/pow(pow(gp,2)+pow(g,2),0.5);
    sinthetaW = gp/pow(pow(gp,2)+pow(g,2),0.5);

    c(1) = mixNeut(neut,1)*costhetaW + mixNeut(neut,2)*sinthetaW;
    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    c(3) = mixNeut(neut,3);
    c(4) = mixNeut(neut,4);
    c(5) = mixNeut(neut,5);
    fb = g*runmb/(root2*mWboson*cos(beta));

    if ( sbottom == 1) {
      coup1 = cos(thetab)*root2*(1./3*c(1)*sinthetaW + (0.5 - 1./3*pow(sinthetaW,2))*c(2)/costhetaW) - sin(thetab)*fb/g*mixNeut(neut,3);
      coup2 = 1./3*sin(thetab)*root2*sinthetaW*(c(2)*gp/g-c(1)) - cos(thetab)*fb/g*mixNeut(neut,3);
    }
    else if (sbottom == 2) {
      coup1 = -sin(thetab)*root2*(1./3*c(1)*sinthetaW + (0.5-1./3*pow(sinthetaW,2))*c(2)/costhetaW) - cos(thetab)*fb/g*mixNeut(neut,3);
      coup2 = 1./3*cos(thetab)*root2*sinthetaW*(c(2)*gp/g - c(1)) + sin(thetab)*fb/g*mixNeut(neut,3);
    }
    else{
      throw("problem: sbottom must be 1 or 2 in neutralinoamplitudesbottombottomNMSSM\n");
    }
    amplitudeW = prefactor*lambda*((pow(coup1,2)+pow(coup2,2))*(pow(mneut,2)-pow(msb,2)+pow(mb,2)) + 4*mb*mneut*coup1*coup2);
  }
  return amplitudeW;
}


double neutralinoamplitudestautauNMSSM (double mneut, double mstau, double mtau, double g, double gp, double thetatau, double beta, double mWboson, DoubleMatrix & mixNeut, double runmtau, int neut, int stau) 
{
  double amplitudeW = 0, squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, ftau = 0, coup1 = 0, coup2 = 0, sinthetaW = 0, costhetaW = 0;
  DoubleVector c(5);
  for (int i=1;i<=5;i++) { c(i) = 0;}
  if (fabs(mneut) < mstau + mtau) {
    amplitudeW = 0;
  }
  else {
    prefactor = pow(g,2)/(32*PI*fabs(mneut));
    squareplus = 1 - pow(mstau/fabs(mneut) + mtau/fabs(mneut),2);
    squareminus = 1 - pow(mstau/fabs(mneut) - mtau/fabs(mneut),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudestautauNMSSM\n");
    } 
    costhetaW = g/pow(pow(gp,2)+pow(g,2),0.5);
    sinthetaW = gp/pow(pow(gp,2)+pow(g,2),0.5);

    c(1) = mixNeut(neut,1)*costhetaW + mixNeut(neut,2)*sinthetaW;
    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    c(3) = mixNeut(neut,3);
    c(4) = mixNeut(neut,4);
    c(5) = mixNeut(neut,5);
    ftau = g*runmtau/(root2*mWboson*cos(beta));

    if ( stau == 1) {
      coup1 = cos(thetatau)*root2*(c(1)*sinthetaW + (0.5 - pow(sinthetaW,2))*c(2)/costhetaW) - sin(thetatau)*ftau/g*mixNeut(neut,3);
      coup2 = sin(thetatau)*root2*sinthetaW*(c(2)*gp/g-c(1)) - cos(thetatau)*ftau/g*mixNeut(neut,3);
    }
    else if (stau == 2) {
      coup1 = -sin(thetatau)*root2*(c(1)*sinthetaW + (0.5-pow(sinthetaW,2))*c(2)/costhetaW) - cos(thetatau)*ftau/g*mixNeut(neut,3);
      coup2 = cos(thetatau)*root2*sinthetaW*(c(2)*gp/g - c(1)) + sin(thetatau)*ftau/g*mixNeut(neut,3);
    }
    else{
      throw("problem: stau must be 1 or 2 in neutralinoamplitudestautauNMSSM");
    }
    amplitudeW = prefactor*lambda*((pow(coup1,2)+pow(coup2,2))*(pow(mneut,2)-pow(mstau,2)+pow(mtau,2)) + 4*mtau*mneut*coup1*coup2);
  }
  return amplitudeW;
}


double neutralinoamplitudestauneutrinotauneutrinoNMSSM (double mneut, double mstaunu, double mtaunu, double g, double gp, DoubleMatrix & mixNeut, int neut) 
{
  double amplitudeW = 0, squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, coup1 = 0, costhetaW = 0, sinthetaW = 0;
  DoubleVector c(5);
  for (int i=1;i<=5;i++) { c(i) = 0;}
  if (fabs(mneut) < mstaunu + mtaunu) {
    amplitudeW = 0;
  }
  else {
    prefactor = pow(g,2)/(32*PI*fabs(mneut));
    squareplus = 1 - pow(mstaunu/fabs(mneut) + mtaunu/fabs(mneut),2);
    squareminus = 1 - pow(mstaunu/fabs(mneut) - mtaunu/fabs(mneut),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in neutralinoamplitudestauneutrinotauneutrinoNMSSM\n");
    } 
    costhetaW = g/pow(pow(gp,2)+pow(g,2),0.5);
    sinthetaW = gp/pow(pow(gp,2)+pow(g,2),0.5);

    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    
    coup1 = -c(2)/(root2*costhetaW);
    
    amplitudeW = prefactor*lambda*pow(coup1,2)*(pow(mneut,2)-pow(mstaunu,2));
  }
  return amplitudeW;
}


double squarkamplitudedecayquarkneutralinoNMSSM (double m1, double mq, double mneut, double g, double gp, DoubleMatrix & mixNeut, char uord, char LorR, int neut) 
{
  double amplitudeW = 0;
  if (m1 < mq + fabs(mneut)) { amplitudeW = 0;}
  else {
    double coupling = 0, squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, sinthetaW = 0, costhetaW = 0;
    DoubleVector c(5);
    for (int i=1;i<=5;i++) { c(i) = 0;}

    sinthetaW = gp/(pow(pow(g,2)+pow(gp,2),0.5));
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));

    squareplus = 1 - pow(mq/m1 + fabs(mneut)/m1,2);
    squareminus = 1 - pow(mq/m1 - fabs(mneut)/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in squarkamplitudedecayquarkneutralinoNMSSM\n");
    } 

    c(1) = mixNeut(neut,1)*costhetaW + mixNeut(neut,2)*sinthetaW;
    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    c(3) = mixNeut(neut,3);
    c(4) = mixNeut(neut,4);
    c(5) = mixNeut(neut,5);

    prefactor = pow(g,2)/(16*PI*m1);
    
    if (uord == 'u' && LorR == 'L') {
      coupling = -root2*(2./3*c(1)*sinthetaW + (0.5 - 2./3*pow(sinthetaW,2))*c(2)/costhetaW);
    }
    else if (uord == 'u' && LorR == 'R') {
      coupling = -2./3*root2*sinthetaW*(c(2)*gp/g-c(1));
    }
    else if (uord == 'd' && LorR == 'L') {
      coupling = root2*(c(1)*1./3*sinthetaW + (0.5-1./3*pow(sinthetaW,2))*c(2)/costhetaW);
    }
    else if (uord == 'd' && LorR == 'R') {
      coupling = root2*sinthetaW*(c(2)*gp/g - c(1))*1./3;
    }
    else{
      throw("problem: uord must be u or d and LorR must be L or R in squarkamplitudedecayquarkneutralinoNMSSM\n");
    }

    amplitudeW = prefactor*pow(coupling,2)*(sqr(m1)-pow(mneut,2) - sqr(mq))*lambda;
    
  }
  return amplitudeW;
}

    
double sleptonamplitudedecayleptonneutralinoNMSSM (double m1, double ml, double mneut, double g, double gp, DoubleMatrix & mixNeut, char uord, char LorR, int neut) 
{
  double amplitudeW = 0;
  if (m1 < ml + fabs(mneut)) { amplitudeW = 0;}
  else {
    double coupling = 0, squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, sinthetaW = 0, costhetaW = 0;
    DoubleVector c(5);
    for (int i=1;i<=5;i++) { c(i) = 0;}

    sinthetaW = gp/(pow(pow(g,2)+pow(gp,2),0.5));
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));

    squareplus = 1 - pow(ml/m1 + fabs(mneut)/m1,2);
    squareminus = 1 - pow(ml/m1 - fabs(mneut)/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in sleptonamplitudedecayleptonneutralinoNMSSM\n");
    } 

    c(1) = mixNeut(neut,1)*costhetaW + mixNeut(neut,2)*sinthetaW;
    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    c(3) = mixNeut(neut,3);
    c(4) = mixNeut(neut,4);
    c(5) = mixNeut(neut,5);

    prefactor = pow(g,2)/(16*PI*m1);
    
    if (uord == 'u' && LorR == 'L') {
      coupling = -c(2)/(root2*costhetaW);
    }
    else if (uord == 'u' && LorR == 'R') {
      coupling = 0; ///no sneutrinoR
    }
    else if (uord == 'd' && LorR == 'L') {
      coupling = root2*(c(1)*sinthetaW + (0.5-pow(sinthetaW,2))*c(2)/costhetaW);
    }
    else if (uord == 'd' && LorR == 'R') {
      coupling = root2*sinthetaW*(c(2)*gp/g - c(1));
    }
    else{
      throw("problem: uord must be u or d and LorR must be L or R in sleptonamplitudedecayleptonneutralinoNMSSM\n");
    }

    amplitudeW = prefactor*pow(coupling,2)*(sqr(m1)-pow(mneut,2) - pow(ml,2))*lambda;
    
  }
  return amplitudeW;
}

double stopamplitudedecaytopneutralinoNMSSM (double m1, double mt, double mneut, double g, double gp, double thetat, DoubleMatrix & mixNeut, double runmt, double mWboson, double beta, int stop, int neut)
{
  double amplitudeW = 0;
  if (m1 < mt + fabs(mneut)) { amplitudeW = 0;}
  else {
    double coupling1 = 0, coupling2 = 0, squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, sinthetaW = 0, costhetaW = 0, ft = 0;
    DoubleVector c(5);
    for (int i=1;i<=5;i++) { c(i) = 0;}

    sinthetaW = gp/(pow(pow(g,2)+pow(gp,2),0.5));
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));

    squareplus = 1 - pow(mt/m1 + fabs(mneut)/m1,2);
    squareminus = 1 - pow(mt/m1 - fabs(mneut)/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stopamplitudedecaytopneutralinoNMSSM\n");
    } 
    ft = g*runmt/(root2*mWboson*sin(beta));

    c(1) = mixNeut(neut,1)*costhetaW + mixNeut(neut,2)*sinthetaW;
    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    c(3) = mixNeut(neut,3);
    c(4) = mixNeut(neut,4);
    c(5) = mixNeut(neut,5);

    prefactor = pow(g,2)/(16*PI*m1);
    
    if (stop == 1) {
      coupling1 = cos(thetat)*root2*(-2./3*c(1)*sinthetaW + (-0.5+2./3*pow(sinthetaW,2))*c(2)/costhetaW) - sin(thetat)*ft/g*mixNeut(neut,4);
      coupling2 = -2./3*sin(thetat)*root2*sinthetaW*(c(2)*gp/g - c(1)) - cos(thetat)*ft/g*mixNeut(neut,4);
    }
    else if (stop == 2) {
      coupling1 = -sin(thetat)*root2*(-2./3*c(1)*sinthetaW + (-0.5+2./3*pow(sinthetaW,2))*c(2)/costhetaW) - cos(thetat)*ft/g*mixNeut(neut,4);
      coupling2 = -2./3*cos(thetat)*root2*sinthetaW*(c(2)*gp/g - c(1)) + sin(thetat)*ft/g*mixNeut(neut,4);
    }
    else{
      throw("problem: stop must be 1 or 2 in stopamplitudedecaytopneutralinoNMSSM\n");
    }
    
    amplitudeW = prefactor*((pow(coupling1,2)+pow(coupling2,2))*(sqr(m1) - pow(mt,2) - pow(mneut,2)) - 4*coupling1*coupling2*mt*mneut)*lambda;

  }
  return amplitudeW;
}


double sbottomamplitudedecaybottomneutralinoNMSSM (double m1, double mb, double mneut, double g, double gp, double thetab, DoubleMatrix & mixNeut, double runmb, double mWboson, double beta, int sbottom, int neut)
{
  double amplitudeW = 0;
  if (m1 < mb + fabs(mneut)) { amplitudeW = 0;}
  else {
    double coupling1 = 0, coupling2 = 0, squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, sinthetaW = 0, costhetaW = 0, fb = 0;
    DoubleVector c(5);
    for (int i=1;i<=5;i++) { c(i) = 0;}

    sinthetaW = gp/(pow(pow(g,2)+pow(gp,2),0.5));
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));

    squareplus = 1 - pow(mb/m1 + fabs(mneut)/m1,2);
    squareminus = 1 - pow(mb/m1 - fabs(mneut)/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in sbottomamplitudedecaybottomneutralinoNMSSM\n");
    } 
    fb = g*runmb/(root2*mWboson*cos(beta));

    c(1) = mixNeut(neut,1)*costhetaW + mixNeut(neut,2)*sinthetaW;
    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    c(3) = mixNeut(neut,3);
    c(4) = mixNeut(neut,4);
    c(5) = mixNeut(neut,5);

    prefactor = pow(g,2)/(16*PI*m1);
    
    if (sbottom == 1) {
      coupling1 = cos(thetab)*root2*(1./3*c(1)*sinthetaW + (0.5-1./3*pow(sinthetaW,2))*c(2)/costhetaW) - sin(thetab)*fb/g*mixNeut(neut,3);
      coupling2 = 1./3*sin(thetab)*root2*sinthetaW*(c(2)*gp/g - c(1)) - cos(thetab)*fb/g*mixNeut(neut,3);
    }
    else if (sbottom == 2) {
      coupling1 = -sin(thetab)*root2*(1./3*c(1)*sinthetaW + (0.5-1./3*pow(sinthetaW,2))*c(2)/costhetaW) - cos(thetab)*fb/g*mixNeut(neut,3);
      coupling2 = 1./3*cos(thetab)*root2*sinthetaW*(c(2)*gp/g - c(1)) + sin(thetab)*fb/g*mixNeut(neut,3);
    }
    else{
      throw("problem: sbottom must be 1 or 2 in sbottomamplitudedecaybottomneutralinoNMSSM\n");
    }
    
    amplitudeW = prefactor*((pow(coupling1,2)+pow(coupling2,2))*(sqr(m1) - pow(mb,2) - pow(mneut,2)) - 4*coupling1*coupling2*mb*mneut)*lambda;
  }
  return amplitudeW;
}


double stauamplitudedecaytauneutralinoNMSSM (double m1, double mtau, double mneut, double g, double gp, double thetatau, DoubleMatrix & mixNeut, double runmtau, double mWboson, double beta, int stau, int neut)
{
  double amplitudeW = 0;
  if (m1 < mtau + fabs(mneut)) { amplitudeW = 0;}
  else {
    double coupling1 = 0, coupling2 = 0, squareplus = 0, squareminus = 0, lambda = 0, prefactor = 0, sinthetaW = 0, costhetaW = 0, ftau = 0;
    DoubleVector c(5);
    for (int i=1;i<=5;i++) { c(i) = 0;}

    sinthetaW = gp/(pow(pow(g,2)+pow(gp,2),0.5));
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));

    squareplus = 1 - pow(mtau/m1 + fabs(mneut)/m1,2);
    squareminus = 1 - pow(mtau/m1 - fabs(mneut)/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in stauamplitudedecaytauneutralinoNMSSM\n");
    } 
    ftau = g*runmtau/(root2*mWboson*cos(beta));

    c(1) = mixNeut(neut,1)*costhetaW + mixNeut(neut,2)*sinthetaW;
    c(2) = -mixNeut(neut,1)*sinthetaW + mixNeut(neut,2)*costhetaW;
    c(3) = mixNeut(neut,3);
    c(4) = mixNeut(neut,4);
    c(5) = mixNeut(neut,5);

    prefactor = pow(g,2)/(16*PI*m1);
    
    if (stau == 1) {
      coupling1 = cos(thetatau)*root2*(c(1)*sinthetaW + (0.5-pow(sinthetaW,2))*c(2)/costhetaW) - sin(thetatau)*ftau/g*mixNeut(neut,3);
      coupling2 = sin(thetatau)*root2*sinthetaW*(c(2)*gp/g - c(1)) - cos(thetatau)*ftau/g*mixNeut(neut,3);
    }
    else if (stau == 2) {
      coupling1 = -sin(thetatau)*root2*(c(1)*sinthetaW + (0.5-pow(sinthetaW,2))*c(2)/costhetaW) - cos(thetatau)*ftau/g*mixNeut(neut,3);
      coupling2 = cos(thetatau)*root2*sinthetaW*(c(2)*gp/g - c(1)) + sin(thetatau)*ftau/g*mixNeut(neut,3);
    }
    else{
      throw("problem: stau must be 1 or 2 in stauamplitudedecaytauneutralinoNMSSM\n");
    }
    
    amplitudeW = prefactor*((pow(coupling1,2)+pow(coupling2,2))*(sqr(m1) - pow(mtau,2) - pow(mneut,2)) - 4*coupling1*coupling2*mtau*mneut)*lambda;

  }
  return amplitudeW;
}


double charginoiamplitudedecayneutralinojHpmNMSSM (double mchar, double mneut, double mHpm, double g, double gp, double thetaL, double thetaR, double beta, DoubleMatrix & mixNeut, double lam, int chargino, int neut)
{
  double amplitudeW = 0;
  if (fabs(mchar) < fabs(mneut) + mHpm) { amplitudeW = 0;}
  else {
    double prefactor = 0, coupling1 = 0, coupling2 = 0, squareplus = 0, squareminus = 0, lambda = 0;
    squareplus = 1 - pow(mHpm/fabs(mchar) + fabs(mneut)/fabs(mchar),2);
    squareminus = 1 - pow(mHpm/fabs(mchar) - fabs(mneut)/fabs(mchar),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in charginoiamplitudedecayneutralinojHpmNMSSM\n");
    } 

    prefactor = pow(g,2)/(32*PI*fabs(mchar));
    
    if (chargino == 1) {
      coupling1 = 1/g*(lam*sin(beta)*mixNeut(neut,5)*cos(thetaR)+ cos(beta)/root2*(gp*mixNeut(neut,1) + g*mixNeut(neut,2))*cos(thetaR) + g*cos(beta)*mixNeut(neut,4)*sin(thetaR));
      coupling2 = 1/g*(lam*cos(beta)*mixNeut(neut,5)*cos(thetaL) - sin(beta)/root2*(gp*mixNeut(neut,1) + g*mixNeut(neut,2))*cos(thetaL) + g*sin(beta)*mixNeut(neut,3)*sin(thetaL));
    }
    else if (chargino == 2) {
      coupling1 = 1/g*(lam*sin(beta)*mixNeut(neut,5)*sin(thetaR)+ cos(beta)/root2*(gp*mixNeut(neut,1) + g*mixNeut(neut,2))*sin(thetaR) - g*cos(beta)*mixNeut(neut,4)*cos(thetaR));
      coupling2 = 1/g*(lam*cos(beta)*mixNeut(neut,5)*sin(thetaL) - sin(beta)/root2*(gp*mixNeut(neut,1) + g*mixNeut(neut,2))*sin(thetaL) - g*sin(beta)*mixNeut(neut,3)*cos(thetaL));
    }
    else{
      throw("problem: chargino must be 1 or 2 in charginoiamplitudedecayneutralinojHpmNMSSM\n");
    }

    amplitudeW = prefactor*lambda*((pow(coupling1,2)+pow(coupling2,2))*(pow(mneut,2)+pow(mchar,2)-pow(mHpm,2)) + 4*coupling1*coupling2*mneut*mchar);

  }
  return amplitudeW;
}


double charginoiamplitudedecayneutralinojWNMSSM (double mchar, double mneut, double mWboson, double g, double gp, double thetaL, double thetaR, DoubleMatrix & mixNeut, int chargino, int neut) 
{
  double amplitudeW = 0, squareplus = 0, squareminus = 0, lambda = 0, coupleL = 0, coupleR = 0, coupling = 0, prefactor = 0, V2 = 0, V1 = 0, U2 = 0, U1 = 0;
  if (fabs(mchar) < fabs(mneut) + mWboson) {
    amplitudeW = 0;
  }
  else {
    prefactor = pow(g,2)/(32*PI*fabs(mchar));
    
    if (chargino == 1) {
      V2 = cos(thetaR); V1 = sin(thetaR); U2 = cos(thetaL); U1 = sin(thetaL);
    }
    else if (chargino == 2) {
      V2 = sin(thetaR); V1 = -cos(thetaR); U2 = sin(thetaL); U1 = -cos(thetaL);
    }
    else{
      throw("problem: chargino must be 1 or 2 in charginoiamplitudedecayneutralinojWNMSSM\n");
    }

    coupleL = -1/root2*mixNeut(neut,4)*V2 + mixNeut(neut,2)*V1;
    coupleR = 1/root2*mixNeut(neut,3)*U2 + mixNeut(neut,2)*U1;

    squareplus = 1 - pow(fabs(mneut)/fabs(mchar) + mWboson/fabs(mchar),2);
    squareminus = 1 - pow(fabs(mneut)/fabs(mchar) - mWboson/fabs(mchar),2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in charginoiamplitudedecayneutralinojWNMSSM\n");
    } 

    coupling = -12*mneut*mchar*coupleL*coupleR + (pow(coupleL,2)+pow(coupleR,2))*((pow(mchar,2)+pow(mneut,2)-pow(mWboson,2)) + (pow(mchar,2)+pow(mWboson,2)-pow(mneut,2))*(pow(mchar,2)-pow(mneut,2)-pow(mWboson,2))/pow(mWboson,2));

    amplitudeW = prefactor*lambda*coupling;
  }
  return amplitudeW;
}


double HpmamplitudecharginojneutralinoiNMSSM (double mHp, double mchar, double mneut, double g, double gp, double beta, double thetaL, double thetaR, double lam, DoubleMatrix & mixNeut, int neutralino, int chargino)
{
  double amplitudeW = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, coupHpZiWjL = 0, coupHpZiWjR = 0, coupling = 0, V1 = 0, V2 = 0, U1 = 0, U2 = 0;

  if (mHp < fabs(mchar) + fabs(mneut)) {
    amplitudeW = 0;
  }
  else {
    prefactor = 1/(16*PI*mHp);
    squareplus = 1 - pow(fabs(mchar)/mHp + fabs(mneut)/mHp,2);
    squareminus = 1 - pow(fabs(mchar)/mHp - fabs(mneut)/mHp,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in HpmamplitudecharginojneutralinoiNMSSM\n");
    } 
    
    if (chargino == 1) {
      V2 = cos(thetaR); V1 = sin(thetaR); U2 = cos(thetaL); U1 = sin(thetaL);
    }
    else if (chargino == 2) {
      V2 = sin(thetaR); V1 = -cos(thetaR); U2 = sin(thetaL); U1 = -cos(thetaL);
    }
    else{
      throw("problem: chargino must be 1 or 2 in HpmamplitudecharginojneutralinoiNMSSM\n");
    }

    coupHpZiWjL = lam*cos(beta)*mixNeut(neutralino,5)*U2 - sin(beta)/(root2)*(gp*mixNeut(neutralino,1)+g*mixNeut(neutralino,2))*U2 + sin(beta)*g*mixNeut(neutralino,3)*U1;
    coupHpZiWjR = lam*sin(beta)*mixNeut(neutralino,5)*V2 + cos(beta)/(root2)*(gp*mixNeut(neutralino,1)+g*mixNeut(neutralino,2))*V2 + cos(beta)*g*mixNeut(neutralino,4)*V1;

    coupling = (pow(coupHpZiWjL,2)+pow(coupHpZiWjR,2))*(pow(mHp,2)-pow(mneut,2)-pow(mchar,2)) - 4*coupHpZiWjL*coupHpZiWjR*mneut*mchar;

    amplitudeW = prefactor*coupling*lambda;

  }
  return amplitudeW;
}

double snutauamplitudedecaynutauneutralinoNMSSM (double m1, double mneut, double g, double gp, DoubleMatrix & mixNeut, int neutralino)
{
  double amplitudeW = 0, prefactor = 0, lambda = 0, squareplus = 0, squareminus = 0, coup = 0, costhetaW = 0, sinthetaW = 0, c2 = 0; ///c2 is c(2) from above formulae
  if (m1 < fabs(mneut)) { amplitudeW = 0;}
  else {
    prefactor = pow(g,2)/(16*PI*m1);
    squareplus = 1 - pow(fabs(mneut)/m1,2);
    squareminus = 1 - pow(fabs(mneut)/m1,2);
    lambda = sqrt(squareplus*squareminus);
    if (squareplus*squareminus < 0) {
      throw ("problem: lambda will give nan in snutauamplitudedecaynutauneutralinoNMSSM\n");
    } 

    sinthetaW = gp/(pow(pow(g,2)+pow(gp,2),0.5));
    costhetaW = g/(pow(pow(g,2)+pow(gp,2),0.5));

    c2 = -mixNeut(neutralino,1)*sinthetaW + mixNeut(neutralino,2)*costhetaW;

    coup = -c2/(root2*costhetaW);

    amplitudeW = prefactor*lambda*(sqr(m1)-pow(mneut,2))*pow(coup,2);

  }
  return amplitudeW;
}


///Functions that calculate the couplings:
DoubleVector squarkmixcharginocouplings (double g, double theta, double beta, double gammaL, double gammaR, double runmt, double runmb, double mWboson, double mch1, double mch2, int torb)
{
  DoubleVector couplings(16);
  for (int i=1; i<=16; i++) { couplings(i) = 0;}

  double fu=0, fd=0, AprimeuW1=0, AprimedW1=0, BW1=0, BprimeW1=0, sq1AprimeW1=0, sq1ch1B1=0, sq1ch1combo=0, sq1ch1angular1=0, sq1ch1angular2=0, sq1ch1B2=0, AprimeuW2=0, AprimedW2=0, BW2=0, BprimeW2=0, sq1AprimeW2=0, sq1ch2B1=0, sq1ch2B2=0, sq1ch2combo=0, sq1ch2angular1=0, sq1ch2angular2=0, sq2AprimeW1=0, sq2ch1B1=0, sq2ch1B2=0, sq2ch1combo=0, sq2ch1angular1=0, sq2ch1angular2=0, sq2AprimeW2=0, sq2ch2B1=0, sq2ch2B2=0, sq2ch2combo=0, sq2ch2angular1=0, sq2ch2angular2=0;

  fu = g*runmt/(root2*mWboson*sin(beta));
  fd = g*runmb/(root2*mWboson*cos(beta));
  AprimeuW1 = -g*sin(gammaL);
  AprimedW1 = -g*sin(gammaR);
  BW1 = -fu*cos(gammaR);
  BprimeW1 = -fd*cos(gammaL);
  
    if (torb == 1 ) {
      sq1AprimeW1 = AprimedW1;
      sq1ch1B1 = BW1;
      sq1ch1B2 = BprimeW1;
    }
    else if (torb == 2) {
      sq1AprimeW1 = AprimeuW1;
      sq1ch1B1 = BprimeW1;
      sq1ch1B2 = BW1;
    }
    else {
      throw("problem: torb must be 1 or 2 in squarkmixcharginocouplings\n");
    }
    AprimeuW2 = -g*cos(gammaL);
    AprimedW2 = -g*cos(gammaR);
    BW2 = fu*sin(gammaR);
    BprimeW2 = fd*sin(gammaL);
    if (torb == 1 ) {
      sq1AprimeW2 = AprimedW2;
      sq1ch2B1 = BW2;
      sq1ch2B2 = BprimeW2;
    }
    else if (torb == 2) {
      sq1AprimeW2 = AprimeuW2;
      sq1ch2B1 = BprimeW2;
      sq1ch2B2 = BW2;
    }
    else {
      throw("problem: torb must be 1 or 2 in squarkmixcharginocouplings\n");
    }
    
    if (torb == 1 ) {
      sq2AprimeW1 = AprimedW1;
      sq2ch1B1 = BW1;
      sq2ch1B2 = BprimeW1;
    }
    else if (torb == 2) {
      sq2AprimeW1 = AprimeuW1;
      sq2ch1B1 = BprimeW1;
      sq2ch1B2 = BW1;
    }
    else {
      throw("problem: torb must be 1 or 2 in squarkmixcharginocouplings\n");
    }

    if (torb == 1 ) {
      sq2AprimeW2 = AprimedW2;
      sq2ch2B1 = BW2;
      sq2ch2B2 = BprimeW2;
    }
    else if (torb == 2) {
      sq2AprimeW2 = AprimeuW2;
      sq2ch2B1 = BprimeW2;
      sq2ch2B2 = BW2;
    }
    else {
      throw("problem: torb must be 1 or 2 in squarkmixcharginocouplings\n");
    }
		
    sq1ch1combo = sq1AprimeW1*cos(theta)-sq1ch1B1*sin(theta);
    sq1ch1angular1 = pow(sq1ch1combo,2) + pow(sq1ch1B2*cos(theta),2);
    sq1ch1angular2 = 4*sq1ch1combo*sq1ch1B2*cos(theta);

    sq1ch2combo = sq1AprimeW2*cos(theta)-sq1ch2B1*sin(theta);
    sq1ch2angular1 = pow(sq1ch2combo,2) + pow(sq1ch2B2*cos(theta),2);
    sq1ch2angular2 = 4*sq1ch2combo*sq1ch2B2*cos(theta);

    sq2ch1combo = sq2AprimeW1*sin(theta)+sq2ch1B1*cos(theta);
    sq2ch1angular1 = pow(sq2ch1combo,2) + pow(sq2ch1B2*sin(theta),2);
    sq2ch1angular2 = 4*sq2ch1combo*sq2ch1B2*sin(theta);

    sq2ch2combo = sq2AprimeW2*sin(theta) + sq2ch2B1*cos(theta);
    sq2ch2angular1 = pow(sq2ch2combo,2) + pow(sq2ch2B2*sin(theta),2);
    sq2ch2angular2 = 4*sq2ch2combo*sq2ch2B2*sin(theta);

    couplings(1) = sq1ch1angular1; couplings(2) = sq1ch1angular2; couplings(3) = sq1ch2angular1; couplings(4) = sq1ch2angular2; couplings(5) = sq2ch1angular1; couplings(6) = sq2ch1angular2; couplings(7) = sq2ch2angular1; couplings(8) = sq2ch2angular2; couplings(9) = sq1ch1combo; couplings(10) = sq1ch2combo; couplings(11) = sq1ch1B2*cos(theta); couplings(12) = sq1ch2B2*cos(theta); couplings(13) = sq2ch1combo; couplings(14) = sq2ch2combo; couplings(15) = sq2ch1B2*sin(theta); couplings(16) = sq2ch2B2*sin(theta);
    return couplings;
}

DoubleVector higgsphisamecharginocouplings (double alpha, double beta, double thetaL, double thetaR) /// calculates the couplings of a neutral higgs (h,H,A) to Wtildai Wtildai - i.e. to two charginos of the same type (mass)
{
  DoubleVector Scoupling(6);
  for (int i=1; i<=6; i++) {
    Scoupling(i) = 0;
  }
  DoubleVector Sih(2), SiH(2), SiA(2);
  for (int i=1; i<=2; i++) {
    Scoupling(i) = 0;
  }

  Sih(1) = 0.5*(-sin(alpha)*sin(thetaR)*cos(thetaL) + cos(alpha)*sin(thetaL)*cos(thetaR));
  Sih(2) = 0.5*(sin(alpha)*cos(thetaR)*sin(thetaL) - cos(alpha)*cos(thetaL)*sin(thetaR));
  SiH(1) = 0.5*(cos(alpha)*sin(thetaR)*cos(thetaL) + sin(alpha)*sin(thetaL)*cos(thetaR));
  SiH(2) = -0.5*(cos(alpha)*cos(thetaR)*sin(thetaL) + sin(alpha)*cos(thetaL)*sin(thetaR));
  SiA(1) = 0.5*(sin(thetaR)*cos(thetaL)*sin(beta) + sin(thetaL)*cos(thetaR)*cos(beta));
  SiA(2) = -0.5*(cos(thetaR)*sin(thetaL)*sin(beta) + cos(thetaL)*sin(thetaR)*cos(beta));
 
  Scoupling(1) = Sih(1); Scoupling(2) = Sih(2); Scoupling(3) = SiH(1); Scoupling(4) = SiH(2); Scoupling(5) = SiA(1); Scoupling(6) = SiA(2);
  return Scoupling;
}

DoubleVector higgsphidifcharginocouplings (double alpha, double beta, double thetaL, double thetaR) /// calculates the couplings of a netural higgs (h,H,A) to Wtildai Wtildaj - i.e. to two charginos of different type (mass)
{
  DoubleVector SPcoupling(6);
  for (int i=1; i<=6; i++) {
    SPcoupling(i) = 0;
  }
  double Sh=0, Ph=0, SH=0, PH=0, SA=0, PA=0;  

  Sh = 0.5*(sin(thetaR)*sin(thetaL)*sin(alpha) + cos(thetaL)*cos(thetaR)*cos(alpha) - sin(thetaL)*sin(thetaR)*cos(alpha) - cos(thetaL)*cos(thetaR)*sin(alpha ));
  Ph = 0.5*(-sin(thetaR)*sin(thetaL)*sin(alpha) - cos(thetaL)*cos(thetaR)*cos(alpha) - sin(thetaL)*sin(thetaR)*cos(alpha) - cos(thetaL)*cos(thetaR)*sin(alpha ));
  SH = 0.5*(-sin(thetaR)*sin(thetaL)*cos(alpha) + cos(thetaL)*cos(thetaR)*sin(alpha) - sin(thetaL)*sin(thetaR)*sin(alpha) + cos(thetaL)*cos(thetaR)*cos(alpha));
  PH = 0.5*(sin(thetaR)*sin(thetaL)*cos(alpha) - cos(thetaL)*cos(thetaR)*sin(alpha) - sin(thetaL)*sin(thetaR)*sin(alpha) + cos(thetaL)*cos(thetaR)*cos(alpha));
  SA = 0.5*(-sin(thetaR)*sin(thetaL)*sin(beta) + cos(thetaL)*cos(thetaR)*cos(beta) + sin(thetaL)*sin(thetaR)*cos(beta) - cos(thetaL)*cos(thetaR)*sin(beta));
  PA = 0.5*(sin(thetaR)*sin(thetaL)*sin(beta) - cos(thetaL)*cos(thetaR)*cos(beta) + sin(thetaL)*sin(thetaR)*cos(beta) - cos(thetaL)*cos(thetaR)*sin(beta));
 
  SPcoupling(1) = Sh; SPcoupling(2) = Ph; SPcoupling(3) = SH; SPcoupling(4) = PH; SPcoupling(5) = SA; SPcoupling(6) = PA;
  return SPcoupling;
}
 
DoubleVector higgshsquarksamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mupq, double mdownq) /// calculates the couplings of light scalar higgs h to two squarks of same handedness
{
  DoubleVector hsqsqcoupling(4);
  for (int i=1; i<=4; i++) {
    hsqsqcoupling(i) = 0;
  }
  double huLuL=0, huRuR=0, hdLdL=0, hdRdR=0;  
  huLuL = g*(mWboson*(0.5 - pow(gp/g,2)/6)*sin(alpha+beta) - pow(mupq,2)*cos(alpha)/(mWboson*sin(beta)));
  hdLdL = g*(mWboson*(-0.5 - pow(gp/g,2)/6)*sin(alpha+beta) + pow(mdownq,2)*sin(alpha)/(mWboson*cos(beta)));
  huRuR = g*(2*mWboson/3*pow(gp/g,2)*sin(alpha+beta) - pow(mupq,2)*cos(alpha)/(mWboson*sin(beta)));
  hdRdR = g*(-mWboson/3*pow(gp/g,2)*sin(alpha+beta) + pow(mdownq,2)*sin(alpha)/(mWboson*cos(beta)));
  hsqsqcoupling(1) = huLuL; hsqsqcoupling(2) = hdLdL; hsqsqcoupling(3) = huRuR; hsqsqcoupling(4) = hdRdR; 
  return hsqsqcoupling;
}    

DoubleVector higgshsquarkdiffhandcouplings (double mWboson, double g, double alpha, double beta, double mupq, double mdownq, double greekmu, double Aup, double Adown) /// calculates the couplings of light scalar higgs h to two squarks of different handedness
{
  DoubleVector hsqsqcoupling(2);
  for (int i=1; i<=2; i++) {
    hsqsqcoupling(i) = 0;
  }
  double huLuR=0, hdLdR=0;
  huLuR = g*mupq/(2*mWboson*sin(beta))*(greekmu*sin(alpha) + Aup*cos(alpha));
  hdLdR = g*mdownq/(2*mWboson*cos(beta))*(-greekmu*cos(alpha) - Adown*sin(alpha));

  hsqsqcoupling(1) = huLuR; hsqsqcoupling(2) = hdLdR; 
  return hsqsqcoupling;
}    

DoubleVector higgsHsquarksamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mupq, double mdownq) /// calculates the couplings of heavy scalar higgs H to two squarks of same handedness
{
  DoubleVector Hsqsqcoupling(4);
  for (int i=1; i<=4; i++) {
    Hsqsqcoupling(i) = 0;
  }
  double HuLuL=0, HuRuR=0, HdLdL=0, HdRdR=0;  
  HuLuL = g*(-mWboson*(0.5 - pow(gp/g,2)/6)*cos(alpha+beta) - pow(mupq,2)*sin(alpha)/(mWboson*sin(beta)));
  HdLdL = g*(mWboson*(0.5 + pow(gp/g,2)/6)*cos(alpha+beta) - pow(mdownq,2)*cos(alpha)/(mWboson*cos(beta)));
  HuRuR = g*(-2*mWboson/3*pow(gp/g,2)*cos(alpha+beta) - pow(mupq,2)*sin(alpha)/(mWboson*sin(beta)));
  HdRdR = g*(mWboson/3*pow(gp/g,2)*cos(alpha+beta) - pow(mdownq,2)*cos(alpha)/(mWboson*cos(beta)));
 
  Hsqsqcoupling(1) = HuLuL; Hsqsqcoupling(2) = HdLdL; Hsqsqcoupling(3) = HuRuR; Hsqsqcoupling(4) = HdRdR; 
  return Hsqsqcoupling;
}    

DoubleVector higgsHsquarkdiffhandcouplings (double mWboson, double g, double alpha, double beta, double mupq, double mdownq, double greekmu, double Aup, double Adown) /// calculates the couplings of heavy scalar higgs H to two squarks of different handedness
{
  DoubleVector Hsqsqcoupling(2);
  for (int i=1; i<=2; i++) {
    Hsqsqcoupling(2) = 0;
  }
  double HuLuR=0, HdLdR=0;
  HuLuR = g*mupq/(2*mWboson*sin(beta))*(-greekmu*cos(alpha) + Aup*sin(alpha));
  HdLdR = g*mdownq/(2*mWboson*cos(beta))*(-greekmu*sin(alpha) + Adown*cos(alpha));
  
  Hsqsqcoupling(1) = HuLuR; Hsqsqcoupling(2) = HdLdR; 
  return Hsqsqcoupling;
}    

DoubleVector higgshsleptonsamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mel) /// calculates the couplings of light scalar higgs h to two sleptons of same handedness
{
  DoubleVector hslslcoupling(3);
  for (int i=1; i<=3; i++) {
    hslslcoupling(3) = 0;
  }
  double hnuLnuL=0, heLeL=0, heReR=0;  
  hnuLnuL = g*(mWboson*(0.5 + pow(gp/g,2)/2)*sin(alpha+beta));
  heLeL = g*(mWboson*(-0.5 + pow(gp/g,2)/2)*sin(alpha+beta) + pow(mel,2)*sin(alpha)/(mWboson*cos(beta)));
  heReR = g*(-mWboson*pow(gp/g,2)*sin(alpha+beta) + pow(mel,2)*sin(alpha)/(mWboson*cos(beta)));
 
  hslslcoupling(1) = hnuLnuL; hslslcoupling(2) = heLeL; hslslcoupling(3) = heReR; 
  return hslslcoupling;
}    
  
DoubleVector higgshsleptondiffhandcouplings (double mWboson, double g, double alpha, double beta, double mel, double greekmu, double Ae) /// calculates the couplings of light scalar higgs h to two sleptons of different handedness
{
  DoubleVector hslslcoupling(1);
  for (int i=1; i<=1; i++) {
    hslslcoupling(1) = 0;
  }
  double heLeR=0;
  heLeR = g*mel/(2*mWboson*cos(beta))*(-greekmu*cos(alpha) - Ae*sin(alpha));
  
  hslslcoupling(1) = heLeR;
  return hslslcoupling;
}    

DoubleVector higgsHsleptonsamehandcouplings (double mWboson, double g, double gp, double alpha, double beta, double mel) /// calculates the couplings of heavy scalar higgs H to two sleptons of same handedness
{
  DoubleVector Hslslcoupling(3);
  for (int i=1; i<=3; i++) {
    Hslslcoupling(3) = 0;
  }
  double HnuLnuL=0, HeLeL=0, HeReR=0;  
  HnuLnuL = g*(-mWboson*(0.5 + pow(gp/g,2)/2)*cos(alpha+beta));
  HeLeL = g*(mWboson*(0.5 - pow(gp/g,2)/2)*cos(alpha+beta) - pow(mel,2)*cos(alpha)/(mWboson*cos(beta)));
  HeReR = g*(mWboson*pow(gp/g,2)*cos(alpha+beta) - pow(mel,2)*cos(alpha)/(mWboson*cos(beta)));
 
  Hslslcoupling(1) = HnuLnuL; Hslslcoupling(2) = HeLeL; Hslslcoupling(3) = HeReR; 
  return Hslslcoupling;
}    

DoubleVector higgsHsleptondiffhandcouplings (double mWboson, double g, double alpha, double beta, double mel, double greekmu, double Ae) /// calculates the couplings of heavy scalar higgs H to two sleptons of different handedness
{
  DoubleVector Hslslcoupling(1);
  for (int i=1; i<=1; i++) {
    Hslslcoupling(i) = 0;
  }
  double HeLeR=0;
  HeLeR = g*mel/(2*mWboson*cos(beta))*(-greekmu*sin(alpha) + Ae*cos(alpha));
  
  Hslslcoupling(1) = HeLeR;
  return Hslslcoupling;
}  

DoubleVector higgsHplussquarkcouplings (double mWboson, double g, double beta, double mup, double mdown, double greekmu, double Au, double Ad) /// calculates the couplings of charged scalar higgs H+ to two squarks
{
  DoubleVector Hplussqsqcoupling(4);
  for (int i=1; i<=4; i++) {
    Hplussqsqcoupling(i) = 0;
  }
  double HplusuLdL=0, HplusuRdR=0, HplusuLdR=0, HplusuRdL=0; 
  HplusuLdL = g/(root2)*(-mWboson*sin(2*beta) + (pow(mdown,2)*tan(beta)+pow(mup,2)/(tan(beta)))/mWboson);
  HplusuRdR = (g*mup*mdown*(tan(beta) + 1/(tan(beta))))/(root2*mWboson);
  HplusuLdR = -g*mdown/(root2*mWboson)*(Ad*tan(beta) + greekmu);
  HplusuRdL = -g*mup/(root2*mWboson)*(Au/(tan(beta)) + greekmu);
  
  Hplussqsqcoupling(1) = HplusuLdL;
  Hplussqsqcoupling(2) = HplusuRdR;
  Hplussqsqcoupling(3) = HplusuLdR;
  Hplussqsqcoupling(4) = HplusuRdL;
  
  return Hplussqsqcoupling;
}      

/// Function to calculate the gluino decay amplitudes
DoubleVector foftau(double mpart, double mcomp) ///f(tau) function for use in h->gamma gamma or Z gamma
{
  double tau=0, fr=0, fi=0, etap=0, etam=0;
  DoubleVector f(3);
  for (int i=1; i<=3; i++) {
    f(i) =0;
  }
  tau = 4*pow(mpart/mcomp,2);
  if (tau >=1) {
    fr = pow((asin(1/pow(tau,0.5))),2);
    fi = 0;
  }
  else if (tau<1) {
    etap = 1 + pow(1-tau,0.5);
    etam = 1 - pow(1-tau,0.5);
    if (etap/etam < 0) {
      throw("problem: log will give nan as etap/etam < 0\n");
    }
    if (etam == 0) {
      throw("problem: will get inf as etam = 0 so etap/etam = inf\n");
    }
    fr = -0.25*pow(log(etap/etam),2)+0.25*sqr(PI);
    fi = 0.5*PI*log(etap/etam);
  }
  f(1)=fr;
  f(2)=fi;
  f(3)=tau;
  return f;
}

DoubleVector goftau(double mpart, double mcomp) ///g(tau) function for use in h->Z gamma
{
  double tau=0, gr=0, gi=0, etap=0, etam=0;
  DoubleVector g(3);
  for (int i=1; i<=3; i++) {
    g(i) =0;
  }
  tau = 4*pow(mpart/mcomp,2);
  if (tau >=1) {
    gr = pow(tau-1,0.5)*((asin(1/pow(tau,0.5))));
    gi = 0;
  }
  else if (tau<1) {
    etap = 1 + pow(1-tau,0.5);
    etam = 1 - pow(1-tau,0.5);
    if (etap/etam < 0) {
      throw("problem: log will give nan as etap/etam < 0\n");
    }
    if (etam == 0) {
      throw("problem: will get inf as etam = 0 so etap/etam = inf\n");
    }
    gr = 0.5*pow(1-tau,0.5)*(log(etap/etam));
    gi = 0.5*pow(1-tau,0.5)*-PI;
  }
  g(1)=gr;
  g(2)=gi;
  g(3)=tau;
  return g;
}

Complex fofqsq(double qSq) {
  const double beta = -0.145, mRhoSq = sqr(0.773), gammaRho = 0.145,
    mRhoPrimeSq = sqr(1.37), gammaRhoPrime = 0.51;
  return (bw(mRhoSq, gammaRho, qSq) + beta *
	  bw(mRhoPrimeSq, gammaRhoPrime, qSq)) / (1.0 + beta);
}

Complex bw(double mSq, double gamma, double qSq) {
  double gammaQsq = gamma;
  double p = sqrt(qSq - 4.0 * sqr(mpiplus)) * 0.5;
  double prho = sqrt(mSq - 4.0 * sqr(mpiplus)) * 0.5;
  gammaQsq *= mSq * p * sqr(p) /
  (prho * sqr(prho) * qSq);
  
  return mSq / (mSq - qSq - Complex(0., 1.) * sqrt(qSq) * gammaQsq);
}

double chToN2piInt(double qSq, const DoubleVector & v) {
  Complex OL11(v.display(1), v.display(2)), OR11(v.display(3), v.display(4));
  double mch = v.display(5), mn = v.display(6);
  double integrand = sqr(fofqsq(qSq).mod()) *
    sqrt(1.0 - 4.0 * sqr(mpiplus) / qSq) * (1.0 - 4.0 * sqr(mpiplus) / qSq) *
    sqrt(lambda(sqr(mch), sqr(mn), qSq)) *
    ( abs(OL11 * OL11 + OR11 * OR11)  *
	   (qSq * (sqr(mch) + sqr(mn) - 2.0 * qSq) +
	    sqr(sqr(mch) - sqr(mn))) -
      12.0 * abs(OL11 * OR11) * qSq * mch * mn);
  return integrand;
}

double charginoToNeutralino2pion(const MssmSoftsusy * m) {
  double mchi1 = fabs(m->displayPhys().mch(1)),
    mneut1 = fabs(m->displayPhys().mneut(1));
  if (mchi1 < mneut1 + 2.0 * mpiplus) return 0.;
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

  DoubleVector v(6);
  v(1) = OL11.real(); v(2) = OL11.imag();
  v(3) = OR11.real(); v(4) = OR11.imag();
  v(5) = mchi1; v(6) = mneut1;
  
  double qSqstart = 4.0 * mpiplus * mpiplus;
  double qSqend = sqr(mchi1 - mneut1);
  double preFactor =sqr(GMU) / (192.0 * sqr(PI) * PI * mchi1 * sqr(mchi1));
  double tol = TOLERANCE;
  double integral = dgauss(chToN2piInt, v, qSqstart, qSqend, tol);

  return integral * preFactor;
}
