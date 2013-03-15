
#include "sarah_nmssm_susy.h"
#include "mathematica_wrappers.hpp"

#define CLASSNAME NMSSMSusyPars

NMSSMSusyPars::NMSSMSusyPars()
   : RGE()
   , Yu(3,3), Yd(3,3), Ye(3,3), Lambdax(0), Kappa(0), g1(0), g2(0), g3(0), vd(0
   ), vu(0), vS(0)

{
   setPars(numberOfParameters);
}

NMSSMSusyPars::NMSSMSusyPars(
   double scale_, double loops_, double thresholds_
   , DoubleMatrix Yu_, DoubleMatrix Yd_, DoubleMatrix Ye_, double Lambdax_,
   double Kappa_, double g1_, double g2_, double g3_, double vd_, double vu_,
   double vS_

)
   : RGE()
   , Yu(Yu_), Yd(Yd_), Ye(Ye_), Lambdax(Lambdax_), Kappa(Kappa_), g1(g1_), g2(
   g2_), g3(g3_), vd(vd_), vu(vu_), vS(vS_)

{
   setPars(numberOfParameters);
   setMu(scale_);
   setLoops(loops_);
   setThresholds(thresholds_);
}

DoubleVector NMSSMSusyPars::beta() const
{
   return calcBeta().display();
}

NMSSMSusyPars NMSSMSusyPars::calcBeta() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;

   const double traceYuAdjYu = trace(Yu*Adj(Yu));
   const double traceYdAdjYd = trace(Yd*Adj(Yd));
   const double traceYeAdjYe = trace(Ye*Adj(Ye));
   const double traceYdAdjYuYuAdjYd = trace(Yd*Adj(Yu)*Yu*Adj(Yd));
   const double traceYuAdjYuYuAdjYu = trace(Yu*Adj(Yu)*Yu*Adj(Yu));
   const double traceYdAdjYdYdAdjYd = trace(Yd*Adj(Yd)*Yd*Adj(Yd));
   const double traceYeAdjYeYeAdjYe = trace(Ye*Adj(Ye)*Ye*Adj(Ye));

   DoubleMatrix beta_Yu(3,3);
   DoubleMatrix beta_Yd(3,3);
   DoubleMatrix beta_Ye(3,3);
   double beta_Lambdax(0);
   double beta_Kappa(0);
   double beta_g1(0);
   double beta_g2(0);
   double beta_g3(0);
   double beta_vd(0);
   double beta_vu(0);
   double beta_vS(0);

   beta_Yu = oneOver16PiSqr*(Yu*((-13*Power(g1,2))/15. - 3*Power(g2,2) -
      (16*Power(g3,2))/3. + 3*traceYuAdjYu + Lambdax*Conj(Lambdax)) + Yu*Adj(Yd
      )*Yd + 3*(Yu*Adj(Yu)*Yu));
   beta_Yd = oneOver16PiSqr*(Yd*((-7*Power(g1,2))/15. - 3*Power(g2,2) - (
      16*Power(g3,2))/3. + 3*traceYdAdjYd + traceYeAdjYe + Lambdax*Conj(Lambdax
      )) + 3*(Yd*Adj(Yd)*Yd) + Yd*Adj(Yu)*Yu);
   beta_Ye = oneOver16PiSqr*(Ye*((-9*Power(g1,2))/5. - 3*Power(g2,2) + 3*
      traceYdAdjYd + traceYeAdjYe + Lambdax*Conj(Lambdax)) + 3*(Ye*Adj(Ye)*Ye))
      ;
   beta_Lambdax = oneOver16PiSqr*((-3*Power(g1,2)*Lambdax)/5. - 3*Power(
      g2,2)*Lambdax + 3*Lambdax*traceYdAdjYd + Lambdax*traceYeAdjYe + 3*Lambdax
      *traceYuAdjYu + 2*Kappa*Lambdax*Conj(Kappa) + 4*Power(Lambdax,2)*Conj(
      Lambdax));
   beta_Kappa = 6*Kappa*oneOver16PiSqr*(Kappa*Conj(Kappa) + Lambdax*Conj(
      Lambdax));
   beta_g1 = (33*Power(g1,3)*oneOver16PiSqr)/5.;
   beta_g2 = Power(g2,3)*oneOver16PiSqr;
   beta_g3 = -3*Power(g3,3)*oneOver16PiSqr;
   beta_vd = (oneOver16PiSqr*vd*(3*Power(g1,2) + 15*Power(g2,2) - 30*
      traceYdAdjYd - 10*traceYeAdjYe - 10*Lambdax*Conj(Lambdax)))/10.;
   beta_vu = (oneOver16PiSqr*vu*(3*(Power(g1,2) + 5*Power(g2,2) - 10*
      traceYuAdjYu) - 10*Lambdax*Conj(Lambdax)))/10.;
   beta_vS = -2*oneOver16PiSqr*vS*(Kappa*Conj(Kappa) + Lambdax*Conj(
      Lambdax));

   if (displayLoops() > 1) {
      beta_Yu = beta_Yu + twoLoop*(Yu*((2743*Power(g1,4))/450. + Power
         (g1,2)*Power(g2,2) + (15*Power(g2,4))/2. + (136*Power(g1,2)*Power(g3,2
         ))/45. + 8*Power(g2,2)*Power(g3,2) - (16*Power(g3,4))/9. - 3*
         traceYdAdjYuYuAdjYd + (4*Power(g1,2)*traceYuAdjYu)/5. + 16*Power(g3,2)
         *traceYuAdjYu - 9*traceYuAdjYuYuAdjYu - Lambdax*(3*traceYdAdjYd +
         traceYeAdjYe)*Conj(Lambdax) - 2*Kappa*Lambdax*Conj(Kappa)*Conj(Lambdax
         ) - 3*Power(Lambdax,2)*Power(Conj(Lambdax),2)) + ((2*Power(g1,2))/5. -
         3*traceYdAdjYd - traceYeAdjYe - Lambdax*Conj(Lambdax))*(Yu*Adj(Yd)*Yd
         ) + (2*Power(g1,2)*(Yu*Adj(Yu)*Yu))/5. + 6*Power(g2,2)*(Yu*Adj(Yu)*Yu)
         - 9*traceYuAdjYu*(Yu*Adj(Yu)*Yu) - 3*Lambdax*Conj(Lambdax)*(Yu*Adj(Yu
         )*Yu) - 2*(Yu*Adj(Yd)*Yd*Adj(Yd)*Yd) - 2*(Yu*Adj(Yd)*Yd*Adj(Yu)*Yu) -
         4*(Yu*Adj(Yu)*Yu*Adj(Yu)*Yu));
      beta_Yd = beta_Yd + twoLoop*(Yd*((287*Power(g1,4))/90. + Power(
         g1,2)*Power(g2,2) + (15*Power(g2,4))/2. + (8*Power(g1,2)*Power(g3,2))
         /9. + 8*Power(g2,2)*Power(g3,2) - (16*Power(g3,4))/9. - (2*Power(g1,2)
         *traceYdAdjYd)/5. + 16*Power(g3,2)*traceYdAdjYd - 9*
         traceYdAdjYdYdAdjYd - 3*traceYdAdjYuYuAdjYd + (6*Power(g1,2)*
         traceYeAdjYe)/5. - 3*traceYeAdjYeYeAdjYe - 3*Lambdax*traceYuAdjYu*Conj
         (Lambdax) - 2*Kappa*Lambdax*Conj(Kappa)*Conj(Lambdax) - 3*Power(
         Lambdax,2)*Power(Conj(Lambdax),2)) + ((4*Power(g1,2))/5. + 6*Power(g2,
         2) - 9*traceYdAdjYd - 3*traceYeAdjYe - 3*Lambdax*Conj(Lambdax))*(Yd*
         Adj(Yd)*Yd) + (4*Power(g1,2)*(Yd*Adj(Yu)*Yu))/5. - 3*traceYuAdjYu*(Yd*
         Adj(Yu)*Yu) - Lambdax*Conj(Lambdax)*(Yd*Adj(Yu)*Yu) - 4*(Yd*Adj(Yd)*Yd
         *Adj(Yd)*Yd) - 2*(Yd*Adj(Yu)*Yu*Adj(Yd)*Yd) - 2*(Yd*Adj(Yu)*Yu*Adj(Yu)
         *Yu));
      beta_Ye = beta_Ye + twoLoop*(Ye*((27*Power(g1,4))/2. + (9*Power(
         g1,2)*Power(g2,2))/5. + (15*Power(g2,4))/2. - (2*Power(g1,2)*
         traceYdAdjYd)/5. + 16*Power(g3,2)*traceYdAdjYd - 9*traceYdAdjYdYdAdjYd
         - 3*traceYdAdjYuYuAdjYd + (6*Power(g1,2)*traceYeAdjYe)/5. - 3*
         traceYeAdjYeYeAdjYe - 3*Lambdax*traceYuAdjYu*Conj(Lambdax) - 2*Kappa*
         Lambdax*Conj(Kappa)*Conj(Lambdax) - 3*Power(Lambdax,2)*Power(Conj(
         Lambdax),2)) + (6*Power(g2,2) - 9*traceYdAdjYd - 3*traceYeAdjYe - 3*
         Lambdax*Conj(Lambdax))*(Ye*Adj(Ye)*Ye) - 4*(Ye*Adj(Ye)*Ye*Adj(Ye)*Ye))
         ;
      beta_Lambdax = beta_Lambdax + -(Lambdax*twoLoop*(-207*Power(g1,4
         ) - 90*Power(g1,2)*Power(g2,2) - 375*Power(g2,4) + 20*Power(g1,2)*
         traceYdAdjYd - 800*Power(g3,2)*traceYdAdjYd + 450*traceYdAdjYdYdAdjYd
         + 300*traceYdAdjYuYuAdjYd - 60*Power(g1,2)*traceYeAdjYe + 150*
         traceYeAdjYeYeAdjYe - 40*Power(g1,2)*traceYuAdjYu - 800*Power(g3,2)*
         traceYuAdjYu + 450*traceYuAdjYuYuAdjYu + 400*Power(Kappa,2)*Power(Conj
         (Kappa),2) - 30*Lambdax*(2*Power(g1,2) + 10*Power(g2,2) - 15*
         traceYdAdjYd - 5*traceYeAdjYe - 15*traceYuAdjYu)*Conj(Lambdax) + 600*
         Kappa*Lambdax*Conj(Kappa)*Conj(Lambdax) + 500*Power(Lambdax,2)*Power(
         Conj(Lambdax),2)))/50.;
      beta_Kappa = beta_Kappa + (-6*Kappa*twoLoop*(20*Power(Kappa,2)*
         Power(Conj(Kappa),2) + 20*Kappa*Lambdax*Conj(Kappa)*Conj(Lambdax) +
         Lambdax*Conj(Lambdax)*(-3*Power(g1,2) - 15*Power(g2,2) + 15*
         traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu + 10*Lambdax*Conj(
         Lambdax))))/5.;
      beta_g1 = beta_g1 + (Power(g1,3)*twoLoop*(199*Power(g1,2) + 135*
         Power(g2,2) + 440*Power(g3,2) - 70*traceYdAdjYd - 90*traceYeAdjYe -
         130*traceYuAdjYu - 30*Lambdax*Conj(Lambdax)))/25.;
      beta_g2 = beta_g2 + (Power(g2,3)*twoLoop*(9*Power(g1,2) + 125*
         Power(g2,2) + 120*Power(g3,2) - 30*traceYdAdjYd - 10*traceYeAdjYe - 30
         *traceYuAdjYu - 10*Lambdax*Conj(Lambdax)))/5.;
      beta_g3 = beta_g3 + (Power(g3,3)*(11*Power(g1,2) + 45*Power(g2,2
         ) + 70*Power(g3,2) - 20*traceYdAdjYd - 20*traceYuAdjYu)*twoLoop)/5.;
      beta_vd = beta_vd + twoLoop*(-((-40*(Power(g1,2) - 40*Power(g3,2
         ))*traceYdAdjYd + 3*(69*Power(g1,4) + 30*Power(g1,2)*Power(g2,2) + 125
         *Power(g2,4) - 300*traceYdAdjYdYdAdjYd - 100*traceYdAdjYuYuAdjYd + 40*
         Power(g1,2)*traceYeAdjYe - 100*traceYeAdjYeYeAdjYe))*vd)/100. + 3*
         Lambdax*traceYuAdjYu*vd*Conj(Lambdax) + 2*Kappa*Lambdax*vd*Conj(Kappa)
         *Conj(Lambdax) + 3*Power(Lambdax,2)*vd*Power(Conj(Lambdax),2));
      beta_vu = beta_vu + twoLoop*(-((207*Power(g1,4) + 90*Power(g1,2)
         *Power(g2,2) + 375*Power(g2,4) - 300*traceYdAdjYuYuAdjYd + 80*(Power(
         g1,2) + 20*Power(g3,2))*traceYuAdjYu - 900*traceYuAdjYuYuAdjYu)*vu)
         /100. + Lambdax*(3*traceYdAdjYd + traceYeAdjYe)*vu*Conj(Lambdax) + 2*
         Kappa*Lambdax*vu*Conj(Kappa)*Conj(Lambdax) + 3*Power(Lambdax,2)*vu*
         Power(Conj(Lambdax),2));
      beta_vS = beta_vS + twoLoop*(8*Power(Kappa,2)*vS*Power(Conj(
         Kappa),2) + 8*Kappa*Lambdax*vS*Conj(Kappa)*Conj(Lambdax) + (2*Lambdax*
         vS*Conj(Lambdax)*(-3*Power(g1,2) - 15*Power(g2,2) + 15*traceYdAdjYd +
         5*traceYeAdjYe + 15*traceYuAdjYu + 10*Lambdax*Conj(Lambdax)))/5.);

   }


   return NMSSMSusyPars(displayMu(), displayLoops(), displayThresholds(),
                    beta_Yu, beta_Yd, beta_Ye, beta_Lambdax, beta_Kappa, beta_g1, beta_g2, beta_g3, beta_vd, beta_vu, beta_vS);
}

DoubleMatrix CLASSNAME::get_SqSq() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   DoubleMatrix anomDim(3,3);

   anomDim = oneOver16PiSqr*(Adj(Yd)*Yd + Adj(Yu)*Yu - ((Power(g1,2) + 45
      *Power(g2,2) + 80*Power(g3,2))*UNITMATRIX(3))/30.);

   if (displayLoops() > 1) {
      anomDim = anomDim + twoLoop*((4*Power(g1,2)*(Adj(Yu)*Yu))/5. -
         Lambdax*Conj(Lambdax)*(Adj(Yu)*Yu) - 2*(Adj(Yd)*Yd*Adj(Yd)*Yd) - 2*(
         Adj(Yu)*Yu*Adj(Yu)*Yu) + Adj(Yd)*Yd*((2*Power(g1,2))/5. - Lambdax*Conj
         (Lambdax) - 3*trace(Yd*Adj(Yd)) - trace(Ye*Adj(Ye))) - 3*(Adj(Yu)*Yu)*
         trace(Yu*Adj(Yu)) + ((199*Power(g1,4))/900. + (15*Power(g2,4))/4. + 8*
         Power(g2,2)*Power(g3,2) - (8*Power(g3,4))/9. + (Power(g1,2)*(9*Power(
         g2,2) + 16*Power(g3,2)))/90.)*UNITMATRIX(3));
   }

   return anomDim;
}

DoubleMatrix CLASSNAME::get_SlSl() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   DoubleMatrix anomDim(3,3);

   anomDim = oneOver16PiSqr*(Adj(Ye)*Ye - (3*(Power(g1,2) + 5*Power(g2,2)
      )*UNITMATRIX(3))/10.);

   if (displayLoops() > 1) {
      anomDim = anomDim + twoLoop*(-2*(Adj(Ye)*Ye*Adj(Ye)*Ye) + Adj(Ye
         )*Ye*((6*Power(g1,2))/5. - Lambdax*Conj(Lambdax) - 3*trace(Yd*Adj(Yd))
         - trace(Ye*Adj(Ye))) + (3*(69*Power(g1,4) + 30*Power(g1,2)*Power(g2,2
         ) + 125*Power(g2,4))*UNITMATRIX(3))/100.);
   }

   return anomDim;
}

double CLASSNAME::get_SHdSHd() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim(0);

   anomDim = oneOver16PiSqr*((-3*Power(g1,2))/10. - (3*Power(g2,2))/2. +
      Lambdax*Conj(Lambdax) + 3*trace(Yd*Adj(Yd)) + trace(Ye*Adj(Ye)));

   if (displayLoops() > 1) {
      anomDim = anomDim + twoLoop*((207*Power(g1,4))/100. + (9*Power(
         g1,2)*Power(g2,2))/10. + (15*Power(g2,4))/4. - 2*Kappa*Lambdax*Conj(
         Kappa)*Conj(Lambdax) - 3*Power(Lambdax,2)*Power(Conj(Lambdax),2) - (2*
         Power(g1,2)*trace(Yd*Adj(Yd)))/5. + 16*Power(g3,2)*trace(Yd*Adj(Yd)) +
         (6*Power(g1,2)*trace(Ye*Adj(Ye)))/5. - 3*Lambdax*Conj(Lambdax)*trace(
         Yu*Adj(Yu)) - 9*trace(Yd*Adj(Yd)*Yd*Adj(Yd)) - 3*trace(Yd*Adj(Yu)*Yu*
         Adj(Yd)) - 3*trace(Ye*Adj(Ye)*Ye*Adj(Ye)));
   }

   return anomDim;
}

double CLASSNAME::get_SHuSHu() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim(0);

   anomDim = oneOver16PiSqr*(Lambdax*Conj(Lambdax) - (3*(Power(g1,2) + 5*
      Power(g2,2) - 10*trace(Yu*Adj(Yu))))/10.);

   if (displayLoops() > 1) {
      anomDim = anomDim + twoLoop*((207*Power(g1,4))/100. + (9*Power(
         g1,2)*Power(g2,2))/10. + (15*Power(g2,4))/4. - 2*Kappa*Lambdax*Conj(
         Kappa)*Conj(Lambdax) - 3*Power(Lambdax,2)*Power(Conj(Lambdax),2) -
         Lambdax*Conj(Lambdax)*(3*trace(Yd*Adj(Yd)) + trace(Ye*Adj(Ye))) + (4*
         Power(g1,2)*trace(Yu*Adj(Yu)))/5. + 16*Power(g3,2)*trace(Yu*Adj(Yu)) -
         3*trace(Yd*Adj(Yu)*Yu*Adj(Yd)) - 9*trace(Yu*Adj(Yu)*Yu*Adj(Yu)));
   }

   return anomDim;
}

DoubleMatrix CLASSNAME::get_SdRSdR() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   DoubleMatrix anomDim(3,3);

   anomDim = oneOver16PiSqr*(2*(Conj(Yd)*Tp(Yd)) - (2*(Power(g1,2) + 20*
      Power(g3,2))*UNITMATRIX(3))/15.);

   if (displayLoops() > 1) {
      anomDim = anomDim + twoLoop*(-2*(Conj(Yd)*Tp(Yd)*Conj(Yd)*Tp(Yd)
         + Conj(Yd)*Tp(Yu)*Conj(Yu)*Tp(Yd)) + Conj(Yd)*Tp(Yd)*((2*Power(g1,2))
         /5. + 6*Power(g2,2) - 2*Lambdax*Conj(Lambdax) - 6*trace(Yd*Adj(Yd)) -
         2*trace(Ye*Adj(Ye))) + (2*(101*Power(g1,4) + 80*Power(g1,2)*Power(g3,2
         ) - 100*Power(g3,4))*UNITMATRIX(3))/225.);
   }

   return anomDim;
}

DoubleMatrix CLASSNAME::get_SuRSuR() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   DoubleMatrix anomDim(3,3);

   anomDim = oneOver16PiSqr*(2*(Conj(Yu)*Tp(Yu)) - (8*(Power(g1,2) + 5*
      Power(g3,2))*UNITMATRIX(3))/15.);

   if (displayLoops() > 1) {
      anomDim = anomDim + twoLoop*((-2*(5*(Conj(Yu)*Tp(Yd)*Conj(Yd)*Tp
         (Yu) + Conj(Yu)*Tp(Yu)*Conj(Yu)*Tp(Yu)) + Conj(Yu)*Tp(Yu)*(Power(g1,2)
         - 15*Power(g2,2) + 5*Lambdax*Conj(Lambdax) + 15*trace(Yu*Adj(Yu)))))
         /5. + (8*(107*Power(g1,4) + 80*Power(g1,2)*Power(g3,2) - 25*Power(g3,4
         ))*UNITMATRIX(3))/225.);
   }

   return anomDim;
}

DoubleMatrix CLASSNAME::get_SeRSeR() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   DoubleMatrix anomDim(3,3);

   anomDim = oneOver16PiSqr*(2*(Conj(Ye)*Tp(Ye)) - (6*Power(g1,2)*
      UNITMATRIX(3))/5.);

   if (displayLoops() > 1) {
      anomDim = anomDim + twoLoop*(-2*(Conj(Ye)*Tp(Ye)*Conj(Ye)*Tp(Ye)
         ) + Conj(Ye)*Tp(Ye)*((-6*Power(g1,2))/5. + 6*Power(g2,2) - 2*Lambdax*
         Conj(Lambdax) - 6*trace(Yd*Adj(Yd)) - 2*trace(Ye*Adj(Ye))) + (234*
         Power(g1,4)*UNITMATRIX(3))/25.);
   }

   return anomDim;
}

double CLASSNAME::get_SsRSsR() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim(0);

   anomDim = 2*oneOver16PiSqr*(Kappa*Conj(Kappa) + Lambdax*Conj(Lambdax))
      ;

   if (displayLoops() > 1) {
      anomDim = anomDim + twoLoop*(-8*Power(Kappa,2)*Power(Conj(Kappa)
         ,2) - 8*Kappa*Lambdax*Conj(Kappa)*Conj(Lambdax) - (2*Lambdax*Conj(
         Lambdax)*(-3*Power(g1,2) - 15*Power(g2,2) + 10*Lambdax*Conj(Lambdax) +
         15*trace(Yd*Adj(Yd)) + 5*trace(Ye*Adj(Ye)) + 15*trace(Yu*Adj(Yu))))
         /5.);
   }

   return anomDim;
}


const DoubleVector NMSSMSusyPars::display() const
{
   DoubleVector pars(numberOfParameters);

   pars(1) = Yu(1,1);
   pars(2) = Yu(1,2);
   pars(3) = Yu(1,3);
   pars(4) = Yu(2,1);
   pars(5) = Yu(2,2);
   pars(6) = Yu(2,3);
   pars(7) = Yu(3,1);
   pars(8) = Yu(3,2);
   pars(9) = Yu(3,3);
   pars(10) = Yd(1,1);
   pars(11) = Yd(1,2);
   pars(12) = Yd(1,3);
   pars(13) = Yd(2,1);
   pars(14) = Yd(2,2);
   pars(15) = Yd(2,3);
   pars(16) = Yd(3,1);
   pars(17) = Yd(3,2);
   pars(18) = Yd(3,3);
   pars(19) = Ye(1,1);
   pars(20) = Ye(1,2);
   pars(21) = Ye(1,3);
   pars(22) = Ye(2,1);
   pars(23) = Ye(2,2);
   pars(24) = Ye(2,3);
   pars(25) = Ye(3,1);
   pars(26) = Ye(3,2);
   pars(27) = Ye(3,3);
   pars(28) = Lambdax;
   pars(29) = Kappa;
   pars(30) = g1;
   pars(31) = g2;
   pars(32) = g3;
   pars(33) = vd;
   pars(34) = vu;
   pars(35) = vS;


   return pars;
}

void NMSSMSusyPars::set(const DoubleVector& v)
{
   Yu(1,1) = v(1);
   Yu(1,2) = v(2);
   Yu(1,3) = v(3);
   Yu(2,1) = v(4);
   Yu(2,2) = v(5);
   Yu(2,3) = v(6);
   Yu(3,1) = v(7);
   Yu(3,2) = v(8);
   Yu(3,3) = v(9);
   Yd(1,1) = v(10);
   Yd(1,2) = v(11);
   Yd(1,3) = v(12);
   Yd(2,1) = v(13);
   Yd(2,2) = v(14);
   Yd(2,3) = v(15);
   Yd(3,1) = v(16);
   Yd(3,2) = v(17);
   Yd(3,3) = v(18);
   Ye(1,1) = v(19);
   Ye(1,2) = v(20);
   Ye(1,3) = v(21);
   Ye(2,1) = v(22);
   Ye(2,2) = v(23);
   Ye(2,3) = v(24);
   Ye(3,1) = v(25);
   Ye(3,2) = v(26);
   Ye(3,3) = v(27);
   Lambdax = v(28);
   Kappa = v(29);
   g1 = v(30);
   g2 = v(31);
   g3 = v(32);
   vd = v(33);
   vu = v(34);
   vS = v(35);

}
