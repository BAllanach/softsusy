
#include "sarah_nmssm_softpars.h"
#include "mathematica_wrappers.hpp"

NMSSMSoftPars::NMSSMSoftPars()
   : NMSSMSusyPars()
   , TYu(3,3), TYd(3,3), TYe(3,3), TLambdax(0), TKappa(0), mq2(3,3), ml2(3,3),
   mHd2(0), mHu2(0), md2(3,3), mu2(3,3), me2(3,3), ms2(0), MassB(0), MassWB(0),
   MassG(0)

{
   setPars(numberOfParameters);
}

NMSSMSoftPars::NMSSMSoftPars(
   const NMSSMSusyPars& susyModel
   , DoubleMatrix TYu_, DoubleMatrix TYd_, DoubleMatrix TYe_, double TLambdax_,
   double TKappa_, DoubleMatrix mq2_, DoubleMatrix ml2_, double mHd2_, double
   mHu2_, DoubleMatrix md2_, DoubleMatrix mu2_, DoubleMatrix me2_, double ms2_,
   double MassB_, double MassWB_, double MassG_

)
   : NMSSMSusyPars(susyModel)
   , TYu(TYu_), TYd(TYd_), TYe(TYe_), TLambdax(TLambdax_), TKappa(TKappa_), mq2
   (mq2_), ml2(ml2_), mHd2(mHd2_), mHu2(mHu2_), md2(md2_), mu2(mu2_), me2(me2_)
   , ms2(ms2_), MassB(MassB_), MassWB(MassWB_), MassG(MassG_)

{
   setPars(numberOfParameters);
}

DoubleVector NMSSMSoftPars::beta() const
{
   return calcBeta().display();
}

NMSSMSoftPars NMSSMSoftPars::calcBeta() const
{
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;

   const double Tr11 = Sqrt(0.6)*g1*(-mHd2 + mHu2 - Conj(trace(ml2)) +
      Conj(trace(mq2)) + trace(md2) + trace(me2) - 2*trace(mu2));
   const double Tr2U111 = (Power(g1,2)*(3*mHd2 + 3*mHu2 + 3*Conj(trace(
      ml2)) + Conj(trace(mq2)) + 2*trace(md2) + 6*trace(me2) + 8*trace(mu2)))
      /10.;
   const double Tr31 = (g1*(-9*Power(g1,2)*mHd2 - 45*Power(g2,2)*mHd2 + 9
      *Power(g1,2)*mHu2 + 45*Power(g2,2)*mHu2 + 30*Lambdax*(mHd2 - mHu2)*Conj(
      Lambdax) - 9*(Power(g1,2) + 5*Power(g2,2))*Conj(trace(ml2)) + Power(g1,2)
      *Conj(trace(mq2)) + 45*Power(g2,2)*Conj(trace(mq2)) + 80*Power(g3,2)*Conj
      (trace(mq2)) + 4*Power(g1,2)*trace(md2) + 80*Power(g3,2)*trace(md2) + 36*
      Power(g1,2)*trace(me2) - 32*Power(g1,2)*trace(mu2) - 160*Power(g3,2)*
      trace(mu2) + 90*mHd2*trace(Yd*Adj(Yd)) + 30*mHd2*trace(Ye*Adj(Ye)) - 90*
      mHu2*trace(Yu*Adj(Yu)) - 60*trace(Yd*Adj(Yd)*Conj(md2)) - 30*trace(Yd*
      Conj(mq2)*Adj(Yd)) - 60*trace(Ye*Adj(Ye)*Conj(me2)) + 30*trace(Ye*Conj(
      ml2)*Adj(Ye)) + 120*trace(Yu*Adj(Yu)*Conj(mu2)) - 30*trace(Yu*Conj(mq2)*
      Adj(Yu))))/(20.*Sqrt(15));
   const double Tr22 = (mHd2 + mHu2 + Conj(trace(ml2)) + 3*Conj(trace(mq2
      )))/2.;
   const double Tr23 = (2*Conj(trace(mq2)) + trace(md2) + trace(mu2))/2.;

   const double traceYuAdjYu = trace(Yu*Adj(Yu));
   const double traceAdjYuTYu = trace(Adj(Yu)*TYu);
   const double traceYdAdjYd = trace(Yd*Adj(Yd));
   const double traceYeAdjYe = trace(Ye*Adj(Ye));
   const double traceYdAdjYuYuAdjYd = trace(Yd*Adj(Yu)*Yu*Adj(Yd));
   const double traceYuAdjYuYuAdjYu = trace(Yu*Adj(Yu)*Yu*Adj(Yu));
   const double traceAdjYdTYd = trace(Adj(Yd)*TYd);
   const double traceAdjYeTYe = trace(Adj(Ye)*TYe);
   const double traceYdAdjYuTYuAdjYd = trace(Yd*Adj(Yu)*TYu*Adj(Yd));
   const double traceYuAdjYdTYdAdjYu = trace(Yu*Adj(Yd)*TYd*Adj(Yu));
   const double traceYuAdjYuTYuAdjYu = trace(Yu*Adj(Yu)*TYu*Adj(Yu));
   const double traceYdAdjYdYdAdjYd = trace(Yd*Adj(Yd)*Yd*Adj(Yd));
   const double traceYeAdjYeYeAdjYe = trace(Ye*Adj(Ye)*Ye*Adj(Ye));
   const double traceYdAdjYdTYdAdjYd = trace(Yd*Adj(Yd)*TYd*Adj(Yd));
   const double traceYeAdjYeTYeAdjYe = trace(Ye*Adj(Ye)*TYe*Adj(Ye));
   const double traceconjTYdTpYd = trace(Conj(TYd)*Tp(Yd));
   const double traceconjTYeTpYe = trace(Conj(TYe)*Tp(Ye));
   const double traceconjTYdTpTYd = trace(Conj(TYd)*Tp(TYd));
   const double traceconjTYeTpTYe = trace(Conj(TYe)*Tp(TYe));
   const double tracemd2YdAdjYd = trace(md2*Yd*Adj(Yd));
   const double traceme2YeAdjYe = trace(me2*Ye*Adj(Ye));
   const double traceml2AdjYeYe = trace(ml2*Adj(Ye)*Ye);
   const double tracemq2AdjYdYd = trace(mq2*Adj(Yd)*Yd);
   const double traceconjTYuTpYu = trace(Conj(TYu)*Tp(Yu));
   const double traceconjTYuTpTYu = trace(Conj(TYu)*Tp(TYu));
   const double tracemq2AdjYuYu = trace(mq2*Adj(Yu)*Yu);
   const double tracemu2YuAdjYu = trace(mu2*Yu*Adj(Yu));
   const double traceYdAdjTYuTYuAdjYd = trace(Yd*Adj(TYu)*TYu*Adj(Yd));
   const double traceYdAdjYuTYuAdjTYd = trace(Yd*Adj(Yu)*TYu*Adj(TYd));
   const double traceYuAdjTYdTYdAdjYu = trace(Yu*Adj(TYd)*TYd*Adj(Yu));
   const double traceYuAdjYdTYdAdjTYu = trace(Yu*Adj(Yd)*TYd*Adj(TYu));
   const double tracemd2YdAdjYuYuAdjYd = trace(md2*Yd*Adj(Yu)*Yu*Adj(Yd))
      ;
   const double tracemq2AdjYdYdAdjYuYu = trace(mq2*Adj(Yd)*Yd*Adj(Yu)*Yu)
      ;
   const double tracemq2AdjYuYuAdjYdYd = trace(mq2*Adj(Yu)*Yu*Adj(Yd)*Yd)
      ;
   const double tracemu2YuAdjYdYdAdjYu = trace(mu2*Yu*Adj(Yd)*Yd*Adj(Yu))
      ;

   DoubleMatrix beta_TYu(3,3);
   DoubleMatrix beta_TYd(3,3);
   DoubleMatrix beta_TYe(3,3);
   double beta_TLambdax(0);
   double beta_TKappa(0);
   DoubleMatrix beta_mq2(3,3);
   DoubleMatrix beta_ml2(3,3);
   double beta_mHd2(0);
   double beta_mHu2(0);
   DoubleMatrix beta_md2(3,3);
   DoubleMatrix beta_mu2(3,3);
   DoubleMatrix beta_me2(3,3);
   double beta_ms2(0);
   double beta_MassB(0);
   double beta_MassWB(0);
   double beta_MassG(0);

   beta_TYu = oneOver16PiSqr*((-13*Power(g1,2)*TYu)/15. - 3*Power(g2,2)*
      TYu - (16*Power(g3,2)*TYu)/3. + 3*traceYuAdjYu*TYu + Lambdax*TYu*Conj(
      Lambdax) + Yu*((26*Power(g1,2)*MassB)/15. + (32*Power(g3,2)*MassG)/3. + 6
      *Power(g2,2)*MassWB + 6*traceAdjYuTYu + 2*TLambdax*Conj(Lambdax)) + TYu*
      Adj(Yd)*Yd + 5*(TYu*Adj(Yu)*Yu) + 2*(Yu*Adj(Yd)*TYd) + 4*(Yu*Adj(Yu)*TYu)
      );
   beta_TYd = oneOver16PiSqr*((-7*Power(g1,2)*TYd)/15. - 3*Power(g2,2)*
      TYd - (16*Power(g3,2)*TYd)/3. + 3*traceYdAdjYd*TYd + traceYeAdjYe*TYd +
      Lambdax*TYd*Conj(Lambdax) + Yd*((14*Power(g1,2)*MassB)/15. + (32*Power(g3
      ,2)*MassG)/3. + 6*Power(g2,2)*MassWB + 6*traceAdjYdTYd + 2*traceAdjYeTYe
      + 2*TLambdax*Conj(Lambdax)) + 5*(TYd*Adj(Yd)*Yd) + TYd*Adj(Yu)*Yu + 4*(Yd
      *Adj(Yd)*TYd) + 2*(Yd*Adj(Yu)*TYu));
   beta_TYe = oneOver16PiSqr*((-9*Power(g1,2)*TYe)/5. - 3*Power(g2,2)*TYe
      + 3*traceYdAdjYd*TYe + traceYeAdjYe*TYe + Lambdax*TYe*Conj(Lambdax) + Ye
      *((18*Power(g1,2)*MassB)/5. + 6*Power(g2,2)*MassWB + 6*traceAdjYdTYd + 2*
      traceAdjYeTYe + 2*TLambdax*Conj(Lambdax)) + 5*(TYe*Adj(Ye)*Ye) + 4*(Ye*
      Adj(Ye)*TYe));
   beta_TLambdax = oneOver16PiSqr*((6*Power(g1,2)*Lambdax*MassB)/5. + 6*
      Power(g2,2)*Lambdax*MassWB + 6*Lambdax*traceAdjYdTYd + 2*Lambdax*
      traceAdjYeTYe + 6*Lambdax*traceAdjYuTYu + 2*(2*Lambdax*TKappa + Kappa*
      TLambdax)*Conj(Kappa) + TLambdax*((-3*Power(g1,2))/5. - 3*Power(g2,2) + 3
      *traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + 12*Lambdax*Conj(Lambdax))
      );
   beta_TKappa = 6*oneOver16PiSqr*(3*Kappa*TKappa*Conj(Kappa) + (Lambdax*
      TKappa + 2*Kappa*TLambdax)*Conj(Lambdax));
   beta_mq2 = oneOver16PiSqr*(2*(Adj(TYd)*TYd) + 2*(Adj(TYu)*TYu) + 2*
      mHd2*(Adj(Yd)*Yd) + 2*mHu2*(Adj(Yu)*Yu) + mq2*Adj(Yd)*Yd + mq2*Adj(Yu)*Yu
      + 2*(Adj(Yd)*md2*Yd) + Adj(Yd)*Yd*mq2 + 2*(Adj(Yu)*mu2*Yu) + Adj(Yu)*Yu*
      mq2 + (g1*Tr11*UNITMATRIX(3))/Sqrt(15) - (2*Power(g1,2)*MassB*Conj(MassB)
      *UNITMATRIX(3))/15. - (32*Power(g3,2)*MassG*Conj(MassG)*UNITMATRIX(3))/3.
      - 6*Power(g2,2)*MassWB*Conj(MassWB)*UNITMATRIX(3));
   beta_ml2 = oneOver16PiSqr*(2*(Adj(TYe)*TYe) + 2*mHd2*(Adj(Ye)*Ye) +
      ml2*Adj(Ye)*Ye + 2*(Adj(Ye)*me2*Ye) + Adj(Ye)*Ye*ml2 - Sqrt(0.6)*g1*Tr11*
      UNITMATRIX(3) - (6*Power(g1,2)*MassB*Conj(MassB)*UNITMATRIX(3))/5. - 6*
      Power(g2,2)*MassWB*Conj(MassWB)*UNITMATRIX(3));
   beta_mHd2 = oneOver16PiSqr*(-(Sqrt(0.6)*g1*Tr11) + 6*traceconjTYdTpTYd
      + 2*traceconjTYeTpTYe + 6*tracemd2YdAdjYd + 2*traceme2YeAdjYe + 2*
      traceml2AdjYeYe + 6*tracemq2AdjYdYd + 6*mHd2*traceYdAdjYd + 2*mHd2*
      traceYeAdjYe + 2*Lambdax*mHd2*Conj(Lambdax) + 2*Lambdax*mHu2*Conj(Lambdax
      ) + 2*Lambdax*ms2*Conj(Lambdax) - (6*Power(g1,2)*MassB*Conj(MassB))/5. -
      6*Power(g2,2)*MassWB*Conj(MassWB) + 2*TLambdax*Conj(TLambdax));
   beta_mHu2 = oneOver16PiSqr*(Sqrt(0.6)*g1*Tr11 + 6*traceconjTYuTpTYu +
      6*tracemq2AdjYuYu + 6*tracemu2YuAdjYu + 6*mHu2*traceYuAdjYu + 2*Lambdax*
      mHd2*Conj(Lambdax) + 2*Lambdax*mHu2*Conj(Lambdax) + 2*Lambdax*ms2*Conj(
      Lambdax) - (6*Power(g1,2)*MassB*Conj(MassB))/5. - 6*Power(g2,2)*MassWB*
      Conj(MassWB) + 2*TLambdax*Conj(TLambdax));
   beta_md2 = oneOver16PiSqr*(4*(TYd*Adj(TYd)) + 4*mHd2*(Yd*Adj(Yd)) + 2*
      (md2*Yd*Adj(Yd)) + 4*(Yd*mq2*Adj(Yd)) + 2*(Yd*Adj(Yd)*md2) + (2*g1*Tr11*
      UNITMATRIX(3))/Sqrt(15) - (8*Power(g1,2)*MassB*Conj(MassB)*UNITMATRIX(3))
      /15. - (32*Power(g3,2)*MassG*Conj(MassG)*UNITMATRIX(3))/3.);
   beta_mu2 = oneOver16PiSqr*(4*(TYu*Adj(TYu)) + 4*mHu2*(Yu*Adj(Yu)) + 2*
      (mu2*Yu*Adj(Yu)) + 4*(Yu*mq2*Adj(Yu)) + 2*(Yu*Adj(Yu)*mu2) - (4*g1*Tr11*
      UNITMATRIX(3))/Sqrt(15) - (32*Power(g1,2)*MassB*Conj(MassB)*UNITMATRIX(3)
      )/15. - (32*Power(g3,2)*MassG*Conj(MassG)*UNITMATRIX(3))/3.);
   beta_me2 = oneOver16PiSqr*(2*(2*(TYe*Adj(TYe)) + 2*mHd2*(Ye*Adj(Ye)) +
      me2*Ye*Adj(Ye) + 2*(Ye*ml2*Adj(Ye)) + Ye*Adj(Ye)*me2) + 2*Sqrt(0.6)*g1*
      Tr11*UNITMATRIX(3) - (24*Power(g1,2)*MassB*Conj(MassB)*UNITMATRIX(3))/5.)
      ;
   beta_ms2 = 4*oneOver16PiSqr*(3*Kappa*ms2*Conj(Kappa) + Lambdax*(mHd2 +
      mHu2 + ms2)*Conj(Lambdax) + TKappa*Conj(TKappa) + TLambdax*Conj(TLambdax
      ));
   beta_MassB = (66*Power(g1,2)*MassB*oneOver16PiSqr)/5.;
   beta_MassWB = 2*Power(g2,2)*MassWB*oneOver16PiSqr;
   beta_MassG = -6*Power(g3,2)*MassG*oneOver16PiSqr;

   if (displayLoops() > 1) {
      beta_TYu = beta_TYu + twoLoop*((2743*Power(g1,4)*TYu)/450. +
         Power(g1,2)*Power(g2,2)*TYu + (15*Power(g2,4)*TYu)/2. + (136*Power(g1,
         2)*Power(g3,2)*TYu)/45. + 8*Power(g2,2)*Power(g3,2)*TYu - (16*Power(g3
         ,4)*TYu)/9. - 3*traceYdAdjYuYuAdjYd*TYu + (4*Power(g1,2)*traceYuAdjYu*
         TYu)/5. + 16*Power(g3,2)*traceYuAdjYu*TYu - 9*traceYuAdjYuYuAdjYu*TYu
         - 3*Lambdax*traceYdAdjYd*TYu*Conj(Lambdax) - Lambdax*traceYeAdjYe*TYu*
         Conj(Lambdax) - 2*Kappa*Lambdax*TYu*Conj(Kappa)*Conj(Lambdax) - 3*
         Power(Lambdax,2)*TYu*Power(Conj(Lambdax),2) - (2*Yu*(2743*Power(g1,4)*
         MassB + 225*Power(g1,2)*Power(g2,2)*MassB + 680*Power(g1,2)*Power(g3,2
         )*MassB + 680*Power(g1,2)*Power(g3,2)*MassG + 1800*Power(g2,2)*Power(
         g3,2)*MassG - 800*Power(g3,4)*MassG + 225*Power(g1,2)*Power(g2,2)*
         MassWB + 3375*Power(g2,4)*MassWB + 1800*Power(g2,2)*Power(g3,2)*MassWB
         - 180*Power(g1,2)*traceAdjYuTYu - 3600*Power(g3,2)*traceAdjYuTYu +
         675*traceYdAdjYuTYuAdjYd + 675*traceYuAdjYdTYdAdjYu + 180*Power(g1,2)*
         MassB*traceYuAdjYu + 3600*Power(g3,2)*MassG*traceYuAdjYu + 4050*
         traceYuAdjYuTYuAdjYu + 225*(Lambdax*(3*traceAdjYdTYd + traceAdjYeTYe)
         + TLambdax*(3*traceYdAdjYd + traceYeAdjYe))*Conj(Lambdax) + 450*(
         Lambdax*TKappa + Kappa*TLambdax)*Conj(Kappa)*Conj(Lambdax) + 1350*
         Lambdax*TLambdax*Power(Conj(Lambdax),2)))/225. + (2*Power(g1,2)*(TYu*
         Adj(Yd)*Yd))/5. - 3*traceYdAdjYd*(TYu*Adj(Yd)*Yd) - traceYeAdjYe*(TYu*
         Adj(Yd)*Yd) - Lambdax*Conj(Lambdax)*(TYu*Adj(Yd)*Yd) + 12*Power(g2,2)*
         (TYu*Adj(Yu)*Yu) - 15*traceYuAdjYu*(TYu*Adj(Yu)*Yu) - 5*Lambdax*Conj(
         Lambdax)*(TYu*Adj(Yu)*Yu) + (4*Power(g1,2)*(Yu*Adj(Yd)*TYd))/5. - 6*
         traceYdAdjYd*(Yu*Adj(Yd)*TYd) - 2*traceYeAdjYe*(Yu*Adj(Yd)*TYd) - 2*
         Lambdax*Conj(Lambdax)*(Yu*Adj(Yd)*TYd) - (2*(2*Power(g1,2)*MassB + 15*
         traceAdjYdTYd + 5*traceAdjYeTYe + 5*TLambdax*Conj(Lambdax))*(Yu*Adj(Yd
         )*Yd))/5. + (6*Power(g1,2)*(Yu*Adj(Yu)*TYu))/5. + 6*Power(g2,2)*(Yu*
         Adj(Yu)*TYu) - 12*traceYuAdjYu*(Yu*Adj(Yu)*TYu) - 4*Lambdax*Conj(
         Lambdax)*(Yu*Adj(Yu)*TYu) - (4*Power(g1,2)*MassB*(Yu*Adj(Yu)*Yu))/5. -
         12*Power(g2,2)*MassWB*(Yu*Adj(Yu)*Yu) - 18*traceAdjYuTYu*(Yu*Adj(Yu)*
         Yu) - 6*TLambdax*Conj(Lambdax)*(Yu*Adj(Yu)*Yu) - 2*(TYu*Adj(Yd)*Yd*Adj
         (Yd)*Yd) - 4*(TYu*Adj(Yd)*Yd*Adj(Yu)*Yu) - 6*(TYu*Adj(Yu)*Yu*Adj(Yu)*
         Yu) - 4*(Yu*Adj(Yd)*TYd*Adj(Yd)*Yd) - 4*(Yu*Adj(Yd)*TYd*Adj(Yu)*Yu) -
         4*(Yu*Adj(Yd)*Yd*Adj(Yd)*TYd) - 2*(Yu*Adj(Yd)*Yd*Adj(Yu)*TYu) - 8*(Yu*
         Adj(Yu)*TYu*Adj(Yu)*Yu) - 6*(Yu*Adj(Yu)*Yu*Adj(Yu)*TYu));
      beta_TYd = beta_TYd + twoLoop*((287*Power(g1,4)*TYd)/90. + Power
         (g1,2)*Power(g2,2)*TYd + (15*Power(g2,4)*TYd)/2. + (8*Power(g1,2)*
         Power(g3,2)*TYd)/9. + 8*Power(g2,2)*Power(g3,2)*TYd - (16*Power(g3,4)*
         TYd)/9. - (2*Power(g1,2)*traceYdAdjYd*TYd)/5. + 16*Power(g3,2)*
         traceYdAdjYd*TYd - 9*traceYdAdjYdYdAdjYd*TYd - 3*traceYdAdjYuYuAdjYd*
         TYd + (6*Power(g1,2)*traceYeAdjYe*TYd)/5. - 3*traceYeAdjYeYeAdjYe*TYd
         - 3*Lambdax*traceYuAdjYu*TYd*Conj(Lambdax) - 2*Kappa*Lambdax*TYd*Conj(
         Kappa)*Conj(Lambdax) - 3*Power(Lambdax,2)*TYd*Power(Conj(Lambdax),2) -
         (2*Yd*(287*Power(g1,4)*MassB + 45*Power(g1,2)*Power(g2,2)*MassB + 40*
         Power(g1,2)*Power(g3,2)*MassB + 40*Power(g1,2)*Power(g3,2)*MassG + 360
         *Power(g2,2)*Power(g3,2)*MassG - 160*Power(g3,4)*MassG + 45*Power(g1,2
         )*Power(g2,2)*MassWB + 675*Power(g2,4)*MassWB + 360*Power(g2,2)*Power(
         g3,2)*MassWB + 18*Power(g1,2)*traceAdjYdTYd - 720*Power(g3,2)*
         traceAdjYdTYd - 54*Power(g1,2)*traceAdjYeTYe - 18*Power(g1,2)*MassB*
         traceYdAdjYd + 720*Power(g3,2)*MassG*traceYdAdjYd + 810*
         traceYdAdjYdTYdAdjYd + 135*traceYdAdjYuTYuAdjYd + 54*Power(g1,2)*MassB
         *traceYeAdjYe + 270*traceYeAdjYeTYeAdjYe + 135*traceYuAdjYdTYdAdjYu +
         135*(Lambdax*traceAdjYuTYu + TLambdax*traceYuAdjYu)*Conj(Lambdax) + 90
         *(Lambdax*TKappa + Kappa*TLambdax)*Conj(Kappa)*Conj(Lambdax) + 270*
         Lambdax*TLambdax*Power(Conj(Lambdax),2)))/45. + (6*Power(g1,2)*(TYd*
         Adj(Yd)*Yd))/5. + 12*Power(g2,2)*(TYd*Adj(Yd)*Yd) - 15*traceYdAdjYd*(
         TYd*Adj(Yd)*Yd) - 5*traceYeAdjYe*(TYd*Adj(Yd)*Yd) - 5*Lambdax*Conj(
         Lambdax)*(TYd*Adj(Yd)*Yd) + (4*Power(g1,2)*(TYd*Adj(Yu)*Yu))/5. - 3*
         traceYuAdjYu*(TYd*Adj(Yu)*Yu) - Lambdax*Conj(Lambdax)*(TYd*Adj(Yu)*Yu)
         + (6*Power(g1,2)*(Yd*Adj(Yd)*TYd))/5. + 6*Power(g2,2)*(Yd*Adj(Yd)*TYd
         ) - 12*traceYdAdjYd*(Yd*Adj(Yd)*TYd) - 4*traceYeAdjYe*(Yd*Adj(Yd)*TYd)
         - 4*Lambdax*Conj(Lambdax)*(Yd*Adj(Yd)*TYd) - (2*(4*Power(g1,2)*MassB
         + 30*Power(g2,2)*MassWB + 45*traceAdjYdTYd + 15*traceAdjYeTYe + 15*
         TLambdax*Conj(Lambdax))*(Yd*Adj(Yd)*Yd))/5. + (8*Power(g1,2)*(Yd*Adj(
         Yu)*TYu))/5. - 6*traceYuAdjYu*(Yd*Adj(Yu)*TYu) - 2*Lambdax*Conj(
         Lambdax)*(Yd*Adj(Yu)*TYu) - (8*Power(g1,2)*MassB*(Yd*Adj(Yu)*Yu))/5. -
         6*traceAdjYuTYu*(Yd*Adj(Yu)*Yu) - 2*TLambdax*Conj(Lambdax)*(Yd*Adj(Yu
         )*Yu) - 6*(TYd*Adj(Yd)*Yd*Adj(Yd)*Yd) - 4*(TYd*Adj(Yu)*Yu*Adj(Yd)*Yd)
         - 2*(TYd*Adj(Yu)*Yu*Adj(Yu)*Yu) - 8*(Yd*Adj(Yd)*TYd*Adj(Yd)*Yd) - 6*(
         Yd*Adj(Yd)*Yd*Adj(Yd)*TYd) - 4*(Yd*Adj(Yu)*TYu*Adj(Yd)*Yd) - 4*(Yd*Adj
         (Yu)*TYu*Adj(Yu)*Yu) - 2*(Yd*Adj(Yu)*Yu*Adj(Yd)*TYd) - 4*(Yd*Adj(Yu)*
         Yu*Adj(Yu)*TYu));
      beta_TYe = beta_TYe + twoLoop*((27*Power(g1,4)*TYe)/2. + (9*
         Power(g1,2)*Power(g2,2)*TYe)/5. + (15*Power(g2,4)*TYe)/2. - (2*Power(
         g1,2)*traceYdAdjYd*TYe)/5. + 16*Power(g3,2)*traceYdAdjYd*TYe - 9*
         traceYdAdjYdYdAdjYd*TYe - 3*traceYdAdjYuYuAdjYd*TYe + (6*Power(g1,2)*
         traceYeAdjYe*TYe)/5. - 3*traceYeAdjYeYeAdjYe*TYe - 3*Lambdax*
         traceYuAdjYu*TYe*Conj(Lambdax) - 2*Kappa*Lambdax*TYe*Conj(Kappa)*Conj(
         Lambdax) - 3*Power(Lambdax,2)*TYe*Power(Conj(Lambdax),2) - (2*Ye*(135*
         Power(g1,4)*MassB + 9*Power(g1,2)*Power(g2,2)*MassB + 9*Power(g1,2)*
         Power(g2,2)*MassWB + 75*Power(g2,4)*MassWB + 2*Power(g1,2)*
         traceAdjYdTYd - 80*Power(g3,2)*traceAdjYdTYd - 6*Power(g1,2)*
         traceAdjYeTYe - 2*Power(g1,2)*MassB*traceYdAdjYd + 80*Power(g3,2)*
         MassG*traceYdAdjYd + 90*traceYdAdjYdTYdAdjYd + 15*traceYdAdjYuTYuAdjYd
         + 6*Power(g1,2)*MassB*traceYeAdjYe + 30*traceYeAdjYeTYeAdjYe + 15*
         traceYuAdjYdTYdAdjYu + 15*(Lambdax*traceAdjYuTYu + TLambdax*
         traceYuAdjYu)*Conj(Lambdax) + 10*(Lambdax*TKappa + Kappa*TLambdax)*
         Conj(Kappa)*Conj(Lambdax) + 30*Lambdax*TLambdax*Power(Conj(Lambdax),2)
         ))/5. - (6*Power(g1,2)*(TYe*Adj(Ye)*Ye))/5. + 12*Power(g2,2)*(TYe*Adj(
         Ye)*Ye) - 15*traceYdAdjYd*(TYe*Adj(Ye)*Ye) - 5*traceYeAdjYe*(TYe*Adj(
         Ye)*Ye) - 5*Lambdax*Conj(Lambdax)*(TYe*Adj(Ye)*Ye) + (6*Power(g1,2)*(
         Ye*Adj(Ye)*TYe))/5. + 6*Power(g2,2)*(Ye*Adj(Ye)*TYe) - 12*traceYdAdjYd
         *(Ye*Adj(Ye)*TYe) - 4*traceYeAdjYe*(Ye*Adj(Ye)*TYe) - 4*Lambdax*Conj(
         Lambdax)*(Ye*Adj(Ye)*TYe) - 6*(2*Power(g2,2)*MassWB + 3*traceAdjYdTYd
         + traceAdjYeTYe + TLambdax*Conj(Lambdax))*(Ye*Adj(Ye)*Ye) - 6*(TYe*Adj
         (Ye)*Ye*Adj(Ye)*Ye) - 8*(Ye*Adj(Ye)*TYe*Adj(Ye)*Ye) - 6*(Ye*Adj(Ye)*Ye
         *Adj(Ye)*TYe));
      beta_TLambdax = beta_TLambdax + twoLoop*((-414*Power(g1,4)*
         Lambdax*MassB)/25. - (18*Power(g1,2)*Power(g2,2)*Lambdax*MassB)/5. - (
         18*Power(g1,2)*Power(g2,2)*Lambdax*MassWB)/5. - 30*Power(g2,4)*Lambdax
         *MassWB + (207*Power(g1,4)*TLambdax)/50. + (9*Power(g1,2)*Power(g2,2)*
         TLambdax)/5. + (15*Power(g2,4)*TLambdax)/2. - (4*Power(g1,2)*Lambdax*
         traceAdjYdTYd)/5. + 32*Power(g3,2)*Lambdax*traceAdjYdTYd + (12*Power(
         g1,2)*Lambdax*traceAdjYeTYe)/5. + (8*Power(g1,2)*Lambdax*traceAdjYuTYu
         )/5. + 32*Power(g3,2)*Lambdax*traceAdjYuTYu + (4*Power(g1,2)*Lambdax*
         MassB*traceYdAdjYd)/5. - 32*Power(g3,2)*Lambdax*MassG*traceYdAdjYd - (
         2*Power(g1,2)*TLambdax*traceYdAdjYd)/5. + 16*Power(g3,2)*TLambdax*
         traceYdAdjYd - 36*Lambdax*traceYdAdjYdTYdAdjYd - 9*TLambdax*
         traceYdAdjYdYdAdjYd - 12*Lambdax*traceYdAdjYuTYuAdjYd - 6*TLambdax*
         traceYdAdjYuYuAdjYd - (12*Power(g1,2)*Lambdax*MassB*traceYeAdjYe)/5. +
         (6*Power(g1,2)*TLambdax*traceYeAdjYe)/5. - 12*Lambdax*
         traceYeAdjYeTYeAdjYe - 3*TLambdax*traceYeAdjYeYeAdjYe - 12*Lambdax*
         traceYuAdjYdTYdAdjYu - (8*Power(g1,2)*Lambdax*MassB*traceYuAdjYu)/5. -
         32*Power(g3,2)*Lambdax*MassG*traceYuAdjYu + (4*Power(g1,2)*TLambdax*
         traceYuAdjYu)/5. + 16*Power(g3,2)*TLambdax*traceYuAdjYu - 36*Lambdax*
         traceYuAdjYuTYuAdjYu - 9*TLambdax*traceYuAdjYuYuAdjYu - 8*Kappa*(4*
         Lambdax*TKappa + Kappa*TLambdax)*Power(Conj(Kappa),2) - (3*Lambdax*(2*
         Lambdax*(2*Power(g1,2)*MassB + 10*Power(g2,2)*MassWB + 15*
         traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu) + TLambdax*(-6*
         Power(g1,2) - 30*Power(g2,2) + 45*traceYdAdjYd + 15*traceYeAdjYe + 45*
         traceYuAdjYu) + 20*(2*Lambdax*TKappa + 3*Kappa*TLambdax)*Conj(Kappa))*
         Conj(Lambdax))/5. - 50*Power(Lambdax,2)*TLambdax*Power(Conj(Lambdax),2
         ));
      beta_TKappa = beta_TKappa + (-6*twoLoop*(100*Power(Kappa,2)*
         TKappa*Power(Conj(Kappa),2) + (Lambdax*TKappa*(-3*Power(g1,2) - 15*
         Power(g2,2) + 15*traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu + 60*
         Kappa*Conj(Kappa)) + 2*Kappa*(Lambdax*(3*Power(g1,2)*MassB + 15*Power(
         g2,2)*MassWB + 15*traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu)
         + TLambdax*(-3*Power(g1,2) - 15*Power(g2,2) + 15*traceYdAdjYd + 5*
         traceYeAdjYe + 15*traceYuAdjYu + 20*Kappa*Conj(Kappa))))*Conj(Lambdax)
         + 10*Lambdax*(Lambdax*TKappa + 4*Kappa*TLambdax)*Power(Conj(Lambdax),
         2)))/5.;
      beta_mq2 = beta_mq2 + twoLoop*((4*Power(g1,2)*(Adj(TYd)*TYd))/5.
         - 6*traceYdAdjYd*(Adj(TYd)*TYd) - 2*traceYeAdjYe*(Adj(TYd)*TYd) - 2*
         Lambdax*Conj(Lambdax)*(Adj(TYd)*TYd) - (4*Power(g1,2)*MassB*(Adj(TYd)*
         Yd))/5. - 6*traceAdjYdTYd*(Adj(TYd)*Yd) - 2*traceAdjYeTYe*(Adj(TYd)*Yd
         ) - 2*TLambdax*Conj(Lambdax)*(Adj(TYd)*Yd) + (8*Power(g1,2)*(Adj(TYu)*
         TYu))/5. - 6*traceYuAdjYu*(Adj(TYu)*TYu) - 2*Lambdax*Conj(Lambdax)*(
         Adj(TYu)*TYu) - (8*Power(g1,2)*MassB*(Adj(TYu)*Yu))/5. - 6*
         traceAdjYuTYu*(Adj(TYu)*Yu) - 2*TLambdax*Conj(Lambdax)*(Adj(TYu)*Yu) -
         6*traceconjTYdTpYd*(Adj(Yd)*TYd) - 2*traceconjTYeTpYe*(Adj(Yd)*TYd) -
         2*Lambdax*Conj(TLambdax)*(Adj(Yd)*TYd) + (4*Power(g1,2)*mHd2*(Adj(Yd)
         *Yd))/5. - 6*traceconjTYdTpTYd*(Adj(Yd)*Yd) - 2*traceconjTYeTpTYe*(Adj
         (Yd)*Yd) - 6*tracemd2YdAdjYd*(Adj(Yd)*Yd) - 2*traceme2YeAdjYe*(Adj(Yd)
         *Yd) - 2*traceml2AdjYeYe*(Adj(Yd)*Yd) - 6*tracemq2AdjYdYd*(Adj(Yd)*Yd)
         - 12*mHd2*traceYdAdjYd*(Adj(Yd)*Yd) - 4*mHd2*traceYeAdjYe*(Adj(Yd)*Yd
         ) - 4*Lambdax*mHd2*Conj(Lambdax)*(Adj(Yd)*Yd) - 2*Lambdax*mHu2*Conj(
         Lambdax)*(Adj(Yd)*Yd) - 2*Lambdax*ms2*Conj(Lambdax)*(Adj(Yd)*Yd) - 2*
         TLambdax*Conj(TLambdax)*(Adj(Yd)*Yd) - 6*traceconjTYuTpYu*(Adj(Yu)*TYu
         ) - 2*Lambdax*Conj(TLambdax)*(Adj(Yu)*TYu) + (8*Power(g1,2)*mHu2*(Adj(
         Yu)*Yu))/5. - 6*traceconjTYuTpTYu*(Adj(Yu)*Yu) - 6*tracemq2AdjYuYu*(
         Adj(Yu)*Yu) - 6*tracemu2YuAdjYu*(Adj(Yu)*Yu) - 12*mHu2*traceYuAdjYu*(
         Adj(Yu)*Yu) - 2*Lambdax*mHd2*Conj(Lambdax)*(Adj(Yu)*Yu) - 4*Lambdax*
         mHu2*Conj(Lambdax)*(Adj(Yu)*Yu) - 2*Lambdax*ms2*Conj(Lambdax)*(Adj(Yu)
         *Yu) - 2*TLambdax*Conj(TLambdax)*(Adj(Yu)*Yu) + (2*Power(g1,2)*(mq2*
         Adj(Yd)*Yd))/5. - 3*traceYdAdjYd*(mq2*Adj(Yd)*Yd) - traceYeAdjYe*(mq2*
         Adj(Yd)*Yd) - Lambdax*Conj(Lambdax)*(mq2*Adj(Yd)*Yd) + (4*Power(g1,2)*
         (mq2*Adj(Yu)*Yu))/5. - 3*traceYuAdjYu*(mq2*Adj(Yu)*Yu) - Lambdax*Conj(
         Lambdax)*(mq2*Adj(Yu)*Yu) + (4*Power(g1,2)*(Adj(Yd)*md2*Yd))/5. - 6*
         traceYdAdjYd*(Adj(Yd)*md2*Yd) - 2*traceYeAdjYe*(Adj(Yd)*md2*Yd) - 2*
         Lambdax*Conj(Lambdax)*(Adj(Yd)*md2*Yd) + (2*Power(g1,2)*(Adj(Yd)*Yd*
         mq2))/5. - 3*traceYdAdjYd*(Adj(Yd)*Yd*mq2) - traceYeAdjYe*(Adj(Yd)*Yd*
         mq2) - Lambdax*Conj(Lambdax)*(Adj(Yd)*Yd*mq2) + (8*Power(g1,2)*(Adj(Yu
         )*mu2*Yu))/5. - 6*traceYuAdjYu*(Adj(Yu)*mu2*Yu) - 2*Lambdax*Conj(
         Lambdax)*(Adj(Yu)*mu2*Yu) + (4*Power(g1,2)*(Adj(Yu)*Yu*mq2))/5. - 3*
         traceYuAdjYu*(Adj(Yu)*Yu*mq2) - Lambdax*Conj(Lambdax)*(Adj(Yu)*Yu*mq2)
         - 4*(Adj(TYd)*TYd*Adj(Yd)*Yd) - 4*(Adj(TYd)*Yd*Adj(Yd)*TYd) - 4*(Adj(
         TYu)*TYu*Adj(Yu)*Yu) - 4*(Adj(TYu)*Yu*Adj(Yu)*TYu) - 4*(Adj(Yd)*TYd*
         Adj(TYd)*Yd) - 4*(Adj(Yd)*Yd*Adj(TYd)*TYd) - 8*mHd2*(Adj(Yd)*Yd*Adj(Yd
         )*Yd) - 4*(Adj(Yu)*TYu*Adj(TYu)*Yu) - 4*(Adj(Yu)*Yu*Adj(TYu)*TYu) - 8*
         mHu2*(Adj(Yu)*Yu*Adj(Yu)*Yu) - 2*(mq2*Adj(Yd)*Yd*Adj(Yd)*Yd) - 2*(mq2*
         Adj(Yu)*Yu*Adj(Yu)*Yu) - 4*(Adj(Yd)*md2*Yd*Adj(Yd)*Yd) - 4*(Adj(Yd)*Yd
         *mq2*Adj(Yd)*Yd) - 4*(Adj(Yd)*Yd*Adj(Yd)*md2*Yd) - 2*(Adj(Yd)*Yd*Adj(
         Yd)*Yd*mq2) - 4*(Adj(Yu)*mu2*Yu*Adj(Yu)*Yu) - 4*(Adj(Yu)*Yu*mq2*Adj(Yu
         )*Yu) - 4*(Adj(Yu)*Yu*Adj(Yu)*mu2*Yu) - 2*(Adj(Yu)*Yu*Adj(Yu)*Yu*mq2)
         + 6*Power(g2,4)*Tr22*UNITMATRIX(3) + (32*Power(g3,4)*Tr23*UNITMATRIX(3
         ))/3. + (2*Power(g1,2)*Tr2U111*UNITMATRIX(3))/15. + (4*g1*Tr31*
         UNITMATRIX(3))/Sqrt(15) + (16*Power(g3,2)*(Power(g1,2)*(MassB + 2*
         MassG) + 15*(-8*Power(g3,2)*MassG + 3*Power(g2,2)*(2*MassG + MassWB)))
         *Conj(MassG)*UNITMATRIX(3))/45. + (Power(g1,2)*Power(g2,2)*MassB*Conj(
         MassWB)*UNITMATRIX(3))/5. + 16*Power(g2,2)*Power(g3,2)*MassG*Conj(
         MassWB)*UNITMATRIX(3) + (2*Power(g1,2)*Power(g2,2)*MassWB*Conj(MassWB)
         *UNITMATRIX(3))/5. + 33*Power(g2,4)*MassWB*Conj(MassWB)*UNITMATRIX(3)
         + 32*Power(g2,2)*Power(g3,2)*MassWB*Conj(MassWB)*UNITMATRIX(3) + (
         Power(g1,2)*Conj(MassB)*(180*(-(Adj(Yd)*TYd) + 2*MassB*(Adj(Yd)*Yd) -
         2*(Adj(Yu)*TYu) + 4*MassB*(Adj(Yu)*Yu)) + (597*Power(g1,2)*MassB + 5*(
         16*Power(g3,2)*(2*MassB + MassG) + 9*Power(g2,2)*(2*MassB + MassWB)))*
         UNITMATRIX(3)))/225.);
      beta_ml2 = beta_ml2 + twoLoop*((12*Power(g1,2)*(Adj(TYe)*TYe))
         /5. - 6*traceYdAdjYd*(Adj(TYe)*TYe) - 2*traceYeAdjYe*(Adj(TYe)*TYe) -
         2*Lambdax*Conj(Lambdax)*(Adj(TYe)*TYe) - (12*Power(g1,2)*MassB*(Adj(
         TYe)*Ye))/5. - 6*traceAdjYdTYd*(Adj(TYe)*Ye) - 2*traceAdjYeTYe*(Adj(
         TYe)*Ye) - 2*TLambdax*Conj(Lambdax)*(Adj(TYe)*Ye) - 6*traceconjTYdTpYd
         *(Adj(Ye)*TYe) - 2*traceconjTYeTpYe*(Adj(Ye)*TYe) - 2*Lambdax*Conj(
         TLambdax)*(Adj(Ye)*TYe) + (12*Power(g1,2)*mHd2*(Adj(Ye)*Ye))/5. - 6*
         traceconjTYdTpTYd*(Adj(Ye)*Ye) - 2*traceconjTYeTpTYe*(Adj(Ye)*Ye) - 6*
         tracemd2YdAdjYd*(Adj(Ye)*Ye) - 2*traceme2YeAdjYe*(Adj(Ye)*Ye) - 2*
         traceml2AdjYeYe*(Adj(Ye)*Ye) - 6*tracemq2AdjYdYd*(Adj(Ye)*Ye) - 12*
         mHd2*traceYdAdjYd*(Adj(Ye)*Ye) - 4*mHd2*traceYeAdjYe*(Adj(Ye)*Ye) - 4*
         Lambdax*mHd2*Conj(Lambdax)*(Adj(Ye)*Ye) - 2*Lambdax*mHu2*Conj(Lambdax)
         *(Adj(Ye)*Ye) - 2*Lambdax*ms2*Conj(Lambdax)*(Adj(Ye)*Ye) - 2*TLambdax*
         Conj(TLambdax)*(Adj(Ye)*Ye) + (6*Power(g1,2)*(ml2*Adj(Ye)*Ye))/5. - 3*
         traceYdAdjYd*(ml2*Adj(Ye)*Ye) - traceYeAdjYe*(ml2*Adj(Ye)*Ye) -
         Lambdax*Conj(Lambdax)*(ml2*Adj(Ye)*Ye) + (12*Power(g1,2)*(Adj(Ye)*me2*
         Ye))/5. - 6*traceYdAdjYd*(Adj(Ye)*me2*Ye) - 2*traceYeAdjYe*(Adj(Ye)*
         me2*Ye) - 2*Lambdax*Conj(Lambdax)*(Adj(Ye)*me2*Ye) + (6*Power(g1,2)*(
         Adj(Ye)*Ye*ml2))/5. - 3*traceYdAdjYd*(Adj(Ye)*Ye*ml2) - traceYeAdjYe*(
         Adj(Ye)*Ye*ml2) - Lambdax*Conj(Lambdax)*(Adj(Ye)*Ye*ml2) - 4*(Adj(TYe)
         *TYe*Adj(Ye)*Ye) - 4*(Adj(TYe)*Ye*Adj(Ye)*TYe) - 4*(Adj(Ye)*TYe*Adj(
         TYe)*Ye) - 4*(Adj(Ye)*Ye*Adj(TYe)*TYe) - 8*mHd2*(Adj(Ye)*Ye*Adj(Ye)*Ye
         ) - 2*(ml2*Adj(Ye)*Ye*Adj(Ye)*Ye) - 4*(Adj(Ye)*me2*Ye*Adj(Ye)*Ye) - 4*
         (Adj(Ye)*Ye*ml2*Adj(Ye)*Ye) - 4*(Adj(Ye)*Ye*Adj(Ye)*me2*Ye) - 2*(Adj(
         Ye)*Ye*Adj(Ye)*Ye*ml2) + 6*Power(g2,4)*Tr22*UNITMATRIX(3) + (6*Power(
         g1,2)*Tr2U111*UNITMATRIX(3))/5. - 4*Sqrt(0.6)*g1*Tr31*UNITMATRIX(3) +
         (3*Power(g2,2)*(55*Power(g2,2)*MassWB + 3*Power(g1,2)*(MassB + 2*
         MassWB))*Conj(MassWB)*UNITMATRIX(3))/5. + (3*Power(g1,2)*Conj(MassB)*(
         -20*(Adj(Ye)*TYe) + 40*MassB*(Adj(Ye)*Ye) + 3*(69*Power(g1,2)*MassB +
         5*Power(g2,2)*(2*MassB + MassWB))*UNITMATRIX(3)))/25.);
      beta_mHd2 = beta_mHd2 + (twoLoop*(Power(g1,2)*(621*Power(g1,2)*
         MassB + 90*Power(g2,2)*MassB + 45*Power(g2,2)*MassWB + 20*
         traceAdjYdTYd - 60*traceAdjYeTYe - 40*MassB*traceYdAdjYd + 120*MassB*
         traceYeAdjYe)*Conj(MassB) + 5*(3*Power(g2,2)*(55*Power(g2,2)*MassWB +
         3*Power(g1,2)*(MassB + 2*MassWB))*Conj(MassWB) - 2*(-15*Power(g2,4)*
         Tr22 - 3*Power(g1,2)*Tr2U111 + 2*Sqrt(15)*g1*Tr31 + 2*Power(g1,2)*
         traceconjTYdTpTYd - 80*Power(g3,2)*traceconjTYdTpTYd - 2*Power(g1,2)*
         MassB*traceconjTYdTpYd + 80*Power(g3,2)*MassG*traceconjTYdTpYd - 6*
         Power(g1,2)*traceconjTYeTpTYe + 6*Power(g1,2)*MassB*traceconjTYeTpYe +
         2*Power(g1,2)*tracemd2YdAdjYd - 80*Power(g3,2)*tracemd2YdAdjYd + 15*
         tracemd2YdAdjYuYuAdjYd - 6*Power(g1,2)*traceme2YeAdjYe - 6*Power(g1,2)
         *traceml2AdjYeYe + 2*Power(g1,2)*tracemq2AdjYdYd - 80*Power(g3,2)*
         tracemq2AdjYdYd + 15*tracemq2AdjYdYdAdjYuYu + 15*
         tracemq2AdjYuYuAdjYdYd + 15*tracemu2YuAdjYdYdAdjYu + 15*
         traceYdAdjTYuTYuAdjYd + 2*Power(g1,2)*mHd2*traceYdAdjYd - 80*Power(g3,
         2)*mHd2*traceYdAdjYd + 90*mHd2*traceYdAdjYdYdAdjYd + 15*
         traceYdAdjYuTYuAdjTYd + 15*mHd2*traceYdAdjYuYuAdjYd + 15*mHu2*
         traceYdAdjYuYuAdjYd - 6*Power(g1,2)*mHd2*traceYeAdjYe + 30*mHd2*
         traceYeAdjYeYeAdjYe + 15*traceYuAdjTYdTYdAdjYu + 15*
         traceYuAdjYdTYdAdjTYu + 30*Power(Lambdax,2)*(mHd2 + mHu2 + ms2)*Power(
         Conj(Lambdax),2) + 80*Power(g3,2)*traceAdjYdTYd*Conj(MassG) - 160*
         Power(g3,2)*MassG*traceYdAdjYd*Conj(MassG) + 15*Lambdax*traceAdjYuTYu*
         Conj(TLambdax) + 15*TLambdax*traceYuAdjYu*Conj(TLambdax) + 10*Conj(
         Kappa)*(Kappa*Lambdax*(mHd2 + mHu2 + 4*ms2)*Conj(Lambdax) + (Lambdax*
         TKappa + Kappa*TLambdax)*Conj(TLambdax)) + 5*Conj(Lambdax)*(2*(Lambdax
         *TKappa + Kappa*TLambdax)*Conj(TKappa) + 3*(Lambdax*traceconjTYuTpTYu
         + TLambdax*traceconjTYuTpYu + Lambdax*tracemq2AdjYuYu + Lambdax*
         tracemu2YuAdjYu + Lambdax*(mHd2 + 2*mHu2 + ms2)*traceYuAdjYu + 4*
         Lambdax*TLambdax*Conj(TLambdax))) + 90*trace(Yd*Adj(TYd)*TYd*Adj(Yd))
         + 90*trace(Yd*Adj(Yd)*TYd*Adj(TYd)) + 30*trace(Ye*Adj(TYe)*TYe*Adj(Ye)
         ) + 30*trace(Ye*Adj(Ye)*TYe*Adj(TYe)) + 90*trace(md2*Yd*Adj(Yd)*Yd*Adj
         (Yd)) + 30*trace(me2*Ye*Adj(Ye)*Ye*Adj(Ye)) + 30*trace(ml2*Adj(Ye)*Ye*
         Adj(Ye)*Ye) + 90*trace(mq2*Adj(Yd)*Yd*Adj(Yd)*Yd)))))/25.;
      beta_mHu2 = beta_mHu2 + (twoLoop*(Power(g1,2)*(621*Power(g1,2)*
         MassB + 90*Power(g2,2)*MassB + 45*Power(g2,2)*MassWB - 40*
         traceAdjYuTYu + 80*MassB*traceYuAdjYu)*Conj(MassB) + 5*(3*Power(g2,2)*
         (55*Power(g2,2)*MassWB + 3*Power(g1,2)*(MassB + 2*MassWB))*Conj(MassWB
         ) - 2*(-15*Power(g2,4)*Tr22 - 3*Power(g1,2)*Tr2U111 - 2*Sqrt(15)*g1*
         Tr31 - 4*Power(g1,2)*traceconjTYuTpTYu - 80*Power(g3,2)*
         traceconjTYuTpTYu + 4*Power(g1,2)*MassB*traceconjTYuTpYu + 80*Power(g3
         ,2)*MassG*traceconjTYuTpYu + 15*tracemd2YdAdjYuYuAdjYd + 15*
         tracemq2AdjYdYdAdjYuYu - 4*Power(g1,2)*tracemq2AdjYuYu - 80*Power(g3,2
         )*tracemq2AdjYuYu + 15*tracemq2AdjYuYuAdjYdYd + 15*
         tracemu2YuAdjYdYdAdjYu - 4*Power(g1,2)*tracemu2YuAdjYu - 80*Power(g3,2
         )*tracemu2YuAdjYu + 15*traceYdAdjTYuTYuAdjYd + 15*
         traceYdAdjYuTYuAdjTYd + 15*mHd2*traceYdAdjYuYuAdjYd + 15*mHu2*
         traceYdAdjYuYuAdjYd + 15*traceYuAdjTYdTYdAdjYu + 15*
         traceYuAdjYdTYdAdjTYu - 4*Power(g1,2)*mHu2*traceYuAdjYu - 80*Power(g3,
         2)*mHu2*traceYuAdjYu + 90*mHu2*traceYuAdjYuYuAdjYu + 30*Power(Lambdax,
         2)*(mHd2 + mHu2 + ms2)*Power(Conj(Lambdax),2) + 80*Power(g3,2)*
         traceAdjYuTYu*Conj(MassG) - 160*Power(g3,2)*MassG*traceYuAdjYu*Conj(
         MassG) + 15*Lambdax*traceAdjYdTYd*Conj(TLambdax) + 5*Lambdax*
         traceAdjYeTYe*Conj(TLambdax) + 15*TLambdax*traceYdAdjYd*Conj(TLambdax)
         + 5*TLambdax*traceYeAdjYe*Conj(TLambdax) + 5*Conj(Lambdax)*(3*Lambdax
         *traceconjTYdTpTYd + 3*TLambdax*traceconjTYdTpYd + Lambdax*
         traceconjTYeTpTYe + TLambdax*traceconjTYeTpYe + 3*Lambdax*
         tracemd2YdAdjYd + Lambdax*traceme2YeAdjYe + Lambdax*traceml2AdjYeYe +
         3*Lambdax*tracemq2AdjYdYd + 6*Lambdax*mHd2*traceYdAdjYd + 3*Lambdax*
         mHu2*traceYdAdjYd + 3*Lambdax*ms2*traceYdAdjYd + 2*Lambdax*mHd2*
         traceYeAdjYe + Lambdax*mHu2*traceYeAdjYe + Lambdax*ms2*traceYeAdjYe +
         2*(Lambdax*TKappa + Kappa*TLambdax)*Conj(TKappa) + 12*Lambdax*TLambdax
         *Conj(TLambdax)) + 10*Conj(Kappa)*(Kappa*Lambdax*(mHd2 + mHu2 + 4*ms2)
         *Conj(Lambdax) + (Lambdax*TKappa + Kappa*TLambdax)*Conj(TLambdax)) +
         90*trace(Yu*Adj(TYu)*TYu*Adj(Yu)) + 90*trace(Yu*Adj(Yu)*TYu*Adj(TYu))
         + 90*trace(mq2*Adj(Yu)*Yu*Adj(Yu)*Yu) + 90*trace(mu2*Yu*Adj(Yu)*Yu*Adj
         (Yu))))))/25.;
      beta_md2 = beta_md2 + twoLoop*((4*Power(g1,2)*(TYd*Adj(TYd)))/5.
         + 12*Power(g2,2)*(TYd*Adj(TYd)) - 12*traceYdAdjYd*(TYd*Adj(TYd)) - 4*
         traceYeAdjYe*(TYd*Adj(TYd)) - 4*Lambdax*Conj(Lambdax)*(TYd*Adj(TYd)) -
         12*traceconjTYdTpYd*(TYd*Adj(Yd)) - 4*traceconjTYeTpYe*(TYd*Adj(Yd))
         - 12*Power(g2,2)*Conj(MassWB)*(TYd*Adj(Yd)) - 4*Lambdax*Conj(TLambdax)
         *(TYd*Adj(Yd)) - (4*Power(g1,2)*MassB*(Yd*Adj(TYd)))/5. - 12*Power(g2,
         2)*MassWB*(Yd*Adj(TYd)) - 12*traceAdjYdTYd*(Yd*Adj(TYd)) - 4*
         traceAdjYeTYe*(Yd*Adj(TYd)) - 4*TLambdax*Conj(Lambdax)*(Yd*Adj(TYd)) +
         (4*Power(g1,2)*mHd2*(Yd*Adj(Yd)))/5. + 12*Power(g2,2)*mHd2*(Yd*Adj(Yd
         )) - 12*traceconjTYdTpTYd*(Yd*Adj(Yd)) - 4*traceconjTYeTpTYe*(Yd*Adj(
         Yd)) - 12*tracemd2YdAdjYd*(Yd*Adj(Yd)) - 4*traceme2YeAdjYe*(Yd*Adj(Yd)
         ) - 4*traceml2AdjYeYe*(Yd*Adj(Yd)) - 12*tracemq2AdjYdYd*(Yd*Adj(Yd)) -
         24*mHd2*traceYdAdjYd*(Yd*Adj(Yd)) - 8*mHd2*traceYeAdjYe*(Yd*Adj(Yd))
         - 8*Lambdax*mHd2*Conj(Lambdax)*(Yd*Adj(Yd)) - 4*Lambdax*mHu2*Conj(
         Lambdax)*(Yd*Adj(Yd)) - 4*Lambdax*ms2*Conj(Lambdax)*(Yd*Adj(Yd)) + 24*
         Power(g2,2)*MassWB*Conj(MassWB)*(Yd*Adj(Yd)) - 4*TLambdax*Conj(
         TLambdax)*(Yd*Adj(Yd)) + (2*Power(g1,2)*(md2*Yd*Adj(Yd)))/5. + 6*Power
         (g2,2)*(md2*Yd*Adj(Yd)) - 6*traceYdAdjYd*(md2*Yd*Adj(Yd)) - 2*
         traceYeAdjYe*(md2*Yd*Adj(Yd)) - 2*Lambdax*Conj(Lambdax)*(md2*Yd*Adj(Yd
         )) + (4*Power(g1,2)*(Yd*mq2*Adj(Yd)))/5. + 12*Power(g2,2)*(Yd*mq2*Adj(
         Yd)) - 12*traceYdAdjYd*(Yd*mq2*Adj(Yd)) - 4*traceYeAdjYe*(Yd*mq2*Adj(
         Yd)) - 4*Lambdax*Conj(Lambdax)*(Yd*mq2*Adj(Yd)) + (2*Power(g1,2)*(Yd*
         Adj(Yd)*md2))/5. + 6*Power(g2,2)*(Yd*Adj(Yd)*md2) - 6*traceYdAdjYd*(Yd
         *Adj(Yd)*md2) - 2*traceYeAdjYe*(Yd*Adj(Yd)*md2) - 2*Lambdax*Conj(
         Lambdax)*(Yd*Adj(Yd)*md2) - 4*(TYd*Adj(TYd)*Yd*Adj(Yd)) - 4*(TYd*Adj(
         TYu)*Yu*Adj(Yd)) - 4*(TYd*Adj(Yd)*Yd*Adj(TYd)) - 4*(TYd*Adj(Yu)*Yu*Adj
         (TYd)) - 4*(Yd*Adj(TYd)*TYd*Adj(Yd)) - 4*(Yd*Adj(TYu)*TYu*Adj(Yd)) - 4
         *(Yd*Adj(Yd)*TYd*Adj(TYd)) - 8*mHd2*(Yd*Adj(Yd)*Yd*Adj(Yd)) - 4*(Yd*
         Adj(Yu)*TYu*Adj(TYd)) - 4*mHd2*(Yd*Adj(Yu)*Yu*Adj(Yd)) - 4*mHu2*(Yd*
         Adj(Yu)*Yu*Adj(Yd)) - 2*(md2*Yd*Adj(Yd)*Yd*Adj(Yd)) - 2*(md2*Yd*Adj(Yu
         )*Yu*Adj(Yd)) - 4*(Yd*mq2*Adj(Yd)*Yd*Adj(Yd)) - 4*(Yd*mq2*Adj(Yu)*Yu*
         Adj(Yd)) - 4*(Yd*Adj(Yd)*md2*Yd*Adj(Yd)) - 4*(Yd*Adj(Yd)*Yd*mq2*Adj(Yd
         )) - 2*(Yd*Adj(Yd)*Yd*Adj(Yd)*md2) - 4*(Yd*Adj(Yu)*mu2*Yu*Adj(Yd)) - 4
         *(Yd*Adj(Yu)*Yu*mq2*Adj(Yd)) - 2*(Yd*Adj(Yu)*Yu*Adj(Yd)*md2) + (32*
         Power(g3,4)*Tr23*UNITMATRIX(3))/3. + (8*Power(g1,2)*Tr2U111*UNITMATRIX
         (3))/15. + (8*g1*Tr31*UNITMATRIX(3))/Sqrt(15) + (64*Power(g3,2)*(-30*
         Power(g3,2)*MassG + Power(g1,2)*(MassB + 2*MassG))*Conj(MassG)*
         UNITMATRIX(3))/45. + (4*Power(g1,2)*Conj(MassB)*(-45*(TYd*Adj(Yd)) +
         90*MassB*(Yd*Adj(Yd)) + 2*(303*Power(g1,2)*MassB + 40*Power(g3,2)*(2*
         MassB + MassG))*UNITMATRIX(3)))/225.);
      beta_mu2 = beta_mu2 + twoLoop*((-4*Power(g1,2)*(TYu*Adj(TYu)))
         /5. + 12*Power(g2,2)*(TYu*Adj(TYu)) - 12*traceYuAdjYu*(TYu*Adj(TYu)) -
         4*Lambdax*Conj(Lambdax)*(TYu*Adj(TYu)) - 12*traceconjTYuTpYu*(TYu*Adj
         (Yu)) - 12*Power(g2,2)*Conj(MassWB)*(TYu*Adj(Yu)) - 4*Lambdax*Conj(
         TLambdax)*(TYu*Adj(Yu)) + (4*Power(g1,2)*MassB*(Yu*Adj(TYu)))/5. - 12*
         Power(g2,2)*MassWB*(Yu*Adj(TYu)) - 12*traceAdjYuTYu*(Yu*Adj(TYu)) - 4*
         TLambdax*Conj(Lambdax)*(Yu*Adj(TYu)) - (4*Power(g1,2)*mHu2*(Yu*Adj(Yu)
         ))/5. + 12*Power(g2,2)*mHu2*(Yu*Adj(Yu)) - 12*traceconjTYuTpTYu*(Yu*
         Adj(Yu)) - 12*tracemq2AdjYuYu*(Yu*Adj(Yu)) - 12*tracemu2YuAdjYu*(Yu*
         Adj(Yu)) - 24*mHu2*traceYuAdjYu*(Yu*Adj(Yu)) - 4*Lambdax*mHd2*Conj(
         Lambdax)*(Yu*Adj(Yu)) - 8*Lambdax*mHu2*Conj(Lambdax)*(Yu*Adj(Yu)) - 4*
         Lambdax*ms2*Conj(Lambdax)*(Yu*Adj(Yu)) + 24*Power(g2,2)*MassWB*Conj(
         MassWB)*(Yu*Adj(Yu)) - 4*TLambdax*Conj(TLambdax)*(Yu*Adj(Yu)) - (2*
         Power(g1,2)*(mu2*Yu*Adj(Yu)))/5. + 6*Power(g2,2)*(mu2*Yu*Adj(Yu)) - 6*
         traceYuAdjYu*(mu2*Yu*Adj(Yu)) - 2*Lambdax*Conj(Lambdax)*(mu2*Yu*Adj(Yu
         )) - (4*Power(g1,2)*(Yu*mq2*Adj(Yu)))/5. + 12*Power(g2,2)*(Yu*mq2*Adj(
         Yu)) - 12*traceYuAdjYu*(Yu*mq2*Adj(Yu)) - 4*Lambdax*Conj(Lambdax)*(Yu*
         mq2*Adj(Yu)) - (2*Power(g1,2)*(Yu*Adj(Yu)*mu2))/5. + 6*Power(g2,2)*(Yu
         *Adj(Yu)*mu2) - 6*traceYuAdjYu*(Yu*Adj(Yu)*mu2) - 2*Lambdax*Conj(
         Lambdax)*(Yu*Adj(Yu)*mu2) - 4*(TYu*Adj(TYd)*Yd*Adj(Yu)) - 4*(TYu*Adj(
         TYu)*Yu*Adj(Yu)) - 4*(TYu*Adj(Yd)*Yd*Adj(TYu)) - 4*(TYu*Adj(Yu)*Yu*Adj
         (TYu)) - 4*(Yu*Adj(TYd)*TYd*Adj(Yu)) - 4*(Yu*Adj(TYu)*TYu*Adj(Yu)) - 4
         *(Yu*Adj(Yd)*TYd*Adj(TYu)) - 4*mHd2*(Yu*Adj(Yd)*Yd*Adj(Yu)) - 4*mHu2*(
         Yu*Adj(Yd)*Yd*Adj(Yu)) - 4*(Yu*Adj(Yu)*TYu*Adj(TYu)) - 8*mHu2*(Yu*Adj(
         Yu)*Yu*Adj(Yu)) - 2*(mu2*Yu*Adj(Yd)*Yd*Adj(Yu)) - 2*(mu2*Yu*Adj(Yu)*Yu
         *Adj(Yu)) - 4*(Yu*mq2*Adj(Yd)*Yd*Adj(Yu)) - 4*(Yu*mq2*Adj(Yu)*Yu*Adj(
         Yu)) - 4*(Yu*Adj(Yd)*md2*Yd*Adj(Yu)) - 4*(Yu*Adj(Yd)*Yd*mq2*Adj(Yu)) -
         2*(Yu*Adj(Yd)*Yd*Adj(Yu)*mu2) - 4*(Yu*Adj(Yu)*mu2*Yu*Adj(Yu)) - 4*(Yu
         *Adj(Yu)*Yu*mq2*Adj(Yu)) - 2*(Yu*Adj(Yu)*Yu*Adj(Yu)*mu2) + (32*Power(
         g3,4)*Tr23*UNITMATRIX(3))/3. + (32*Power(g1,2)*Tr2U111*UNITMATRIX(3))
         /15. - (16*g1*Tr31*UNITMATRIX(3))/Sqrt(15) - (128*Power(g3,2)*(15*
         Power(g3,2)*MassG - 2*Power(g1,2)*(MassB + 2*MassG))*Conj(MassG)*
         UNITMATRIX(3))/45. + (4*Power(g1,2)*Conj(MassB)*(45*(TYu*Adj(Yu) - 2*
         MassB*(Yu*Adj(Yu))) + 8*(321*Power(g1,2)*MassB + 40*Power(g3,2)*(2*
         MassB + MassG))*UNITMATRIX(3)))/225.);
      beta_me2 = beta_me2 + (2*twoLoop*(-5*(6*Power(g1,2)*(TYe*Adj(TYe
         )) - 30*Power(g2,2)*(TYe*Adj(TYe)) + 30*traceYdAdjYd*(TYe*Adj(TYe)) +
         10*traceYeAdjYe*(TYe*Adj(TYe)) + 10*Lambdax*Conj(Lambdax)*(TYe*Adj(TYe
         )) + 30*traceconjTYdTpYd*(TYe*Adj(Ye)) + 10*traceconjTYeTpYe*(TYe*Adj(
         Ye)) + 30*Power(g2,2)*Conj(MassWB)*(TYe*Adj(Ye)) + 10*Lambdax*Conj(
         TLambdax)*(TYe*Adj(Ye)) + (-6*Power(g1,2)*MassB + 30*Power(g2,2)*
         MassWB + 30*traceAdjYdTYd + 10*traceAdjYeTYe + 10*TLambdax*Conj(
         Lambdax))*(Ye*Adj(TYe)) + 2*(3*Power(g1,2)*mHd2 - 15*Power(g2,2)*mHd2
         + 15*traceconjTYdTpTYd + 5*traceconjTYeTpTYe + 15*tracemd2YdAdjYd + 5*
         traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*tracemq2AdjYdYd + 30*mHd2*
         traceYdAdjYd + 10*mHd2*traceYeAdjYe + 5*Lambdax*(2*mHd2 + mHu2 + ms2)*
         Conj(Lambdax) - 30*Power(g2,2)*MassWB*Conj(MassWB) + 5*TLambdax*Conj(
         TLambdax))*(Ye*Adj(Ye)) + 3*Power(g1,2)*(me2*Ye*Adj(Ye)) - 15*Power(g2
         ,2)*(me2*Ye*Adj(Ye)) + 15*traceYdAdjYd*(me2*Ye*Adj(Ye)) + 5*
         traceYeAdjYe*(me2*Ye*Adj(Ye)) + 5*Lambdax*Conj(Lambdax)*(me2*Ye*Adj(Ye
         )) + 6*Power(g1,2)*(Ye*ml2*Adj(Ye)) - 30*Power(g2,2)*(Ye*ml2*Adj(Ye))
         + 30*traceYdAdjYd*(Ye*ml2*Adj(Ye)) + 10*traceYeAdjYe*(Ye*ml2*Adj(Ye))
         + 10*Lambdax*Conj(Lambdax)*(Ye*ml2*Adj(Ye)) + 3*Power(g1,2)*(Ye*Adj(Ye
         )*me2) - 15*Power(g2,2)*(Ye*Adj(Ye)*me2) + 15*traceYdAdjYd*(Ye*Adj(Ye)
         *me2) + 5*traceYeAdjYe*(Ye*Adj(Ye)*me2) + 5*Lambdax*Conj(Lambdax)*(Ye*
         Adj(Ye)*me2) + 10*(TYe*Adj(TYe)*Ye*Adj(Ye)) + 10*(TYe*Adj(Ye)*Ye*Adj(
         TYe)) + 10*(Ye*Adj(TYe)*TYe*Adj(Ye)) + 10*(Ye*Adj(Ye)*TYe*Adj(TYe)) +
         20*mHd2*(Ye*Adj(Ye)*Ye*Adj(Ye)) + 5*(me2*Ye*Adj(Ye)*Ye*Adj(Ye)) + 10*(
         Ye*ml2*Adj(Ye)*Ye*Adj(Ye)) + 10*(Ye*Adj(Ye)*me2*Ye*Adj(Ye)) + 10*(Ye*
         Adj(Ye)*Ye*ml2*Adj(Ye)) + 5*(Ye*Adj(Ye)*Ye*Adj(Ye)*me2)) + 20*g1*(3*g1
         *Tr2U111 + Sqrt(15)*Tr31)*UNITMATRIX(3) + 6*Power(g1,2)*Conj(MassB)*(5
         *(TYe*Adj(Ye) - 2*MassB*(Ye*Adj(Ye))) + 234*Power(g1,2)*MassB*
         UNITMATRIX(3))))/25.;
      beta_ms2 = beta_ms2 + (-4*twoLoop*(120*Power(Kappa,2)*ms2*Power(
         Conj(Kappa),2) + 20*Power(Lambdax,2)*(mHd2 + mHu2 + ms2)*Power(Conj(
         Lambdax),2) + (Lambdax*(15*traceAdjYdTYd + 5*traceAdjYeTYe + 3*(Power(
         g1,2)*MassB + 5*Power(g2,2)*MassWB + 5*traceAdjYuTYu)) + TLambdax*(15*
         traceYdAdjYd + 5*traceYeAdjYe - 3*(Power(g1,2) + 5*Power(g2,2) - 5*
         traceYuAdjYu)))*Conj(TLambdax) + Conj(Lambdax)*(-3*Power(g1,2)*Lambdax
         *mHd2 - 15*Power(g2,2)*Lambdax*mHd2 - 3*Power(g1,2)*Lambdax*mHu2 - 15*
         Power(g2,2)*Lambdax*mHu2 - 3*Power(g1,2)*Lambdax*ms2 - 15*Power(g2,2)*
         Lambdax*ms2 + 15*Lambdax*traceconjTYdTpTYd + 15*TLambdax*
         traceconjTYdTpYd + 5*Lambdax*traceconjTYeTpTYe + 5*TLambdax*
         traceconjTYeTpYe + 15*Lambdax*traceconjTYuTpTYu + 15*TLambdax*
         traceconjTYuTpYu + 15*Lambdax*tracemd2YdAdjYd + 5*Lambdax*
         traceme2YeAdjYe + 5*Lambdax*traceml2AdjYeYe + 15*Lambdax*
         tracemq2AdjYdYd + 15*Lambdax*tracemq2AdjYuYu + 15*Lambdax*
         tracemu2YuAdjYu + 30*Lambdax*mHd2*traceYdAdjYd + 15*Lambdax*mHu2*
         traceYdAdjYd + 15*Lambdax*ms2*traceYdAdjYd + 10*Lambdax*mHd2*
         traceYeAdjYe + 5*Lambdax*mHu2*traceYeAdjYe + 5*Lambdax*ms2*
         traceYeAdjYe + 15*Lambdax*mHd2*traceYuAdjYu + 30*Lambdax*mHu2*
         traceYuAdjYu + 15*Lambdax*ms2*traceYuAdjYu + 3*Power(g1,2)*(-2*Lambdax
         *MassB + TLambdax)*Conj(MassB) + 15*Power(g2,2)*(-2*Lambdax*MassWB +
         TLambdax)*Conj(MassWB) + 20*Lambdax*TKappa*Conj(TKappa) + 20*Kappa*
         TLambdax*Conj(TKappa) + 40*Lambdax*TLambdax*Conj(TLambdax)) + 20*Conj(
         Kappa)*(Kappa*Lambdax*(mHd2 + mHu2 + 4*ms2)*Conj(Lambdax) + 4*Kappa*
         TKappa*Conj(TKappa) + (Lambdax*TKappa + Kappa*TLambdax)*Conj(TLambdax)
         )))/5.;
      beta_MassB = beta_MassB + (2*Power(g1,2)*twoLoop*(398*Power(g1,2
         )*MassB + 135*Power(g2,2)*MassB + 440*Power(g3,2)*MassB + 440*Power(g3
         ,2)*MassG + 135*Power(g2,2)*MassWB + 70*traceAdjYdTYd + 90*
         traceAdjYeTYe + 130*traceAdjYuTYu - 70*MassB*traceYdAdjYd - 90*MassB*
         traceYeAdjYe - 130*MassB*traceYuAdjYu - 30*(Lambdax*MassB - TLambdax)*
         Conj(Lambdax)))/25.;
      beta_MassWB = beta_MassWB + (2*Power(g2,2)*twoLoop*(9*Power(g1,2
         )*MassB + 120*Power(g3,2)*MassG + 9*Power(g1,2)*MassWB + 250*Power(g2,
         2)*MassWB + 120*Power(g3,2)*MassWB + 30*traceAdjYdTYd + 10*
         traceAdjYeTYe + 30*traceAdjYuTYu - 30*MassWB*traceYdAdjYd - 10*MassWB*
         traceYeAdjYe - 30*MassWB*traceYuAdjYu - 10*(Lambdax*MassWB - TLambdax)
         *Conj(Lambdax)))/5.;
      beta_MassG = beta_MassG + (2*Power(g3,2)*(11*Power(g1,2)*MassB +
         11*Power(g1,2)*MassG + 45*Power(g2,2)*MassG + 140*Power(g3,2)*MassG +
         45*Power(g2,2)*MassWB + 20*traceAdjYdTYd + 20*traceAdjYuTYu - 20*
         MassG*traceYdAdjYd - 20*MassG*traceYuAdjYu)*twoLoop)/5.;

   }


   NMSSMSusyPars susyModelBeta(NMSSMSusyPars::calcBeta());

   return NMSSMSoftPars(susyModelBeta, beta_TYu, beta_TYd, beta_TYe, beta_TLambdax, beta_TKappa, beta_mq2, beta_ml2, beta_mHd2, beta_mHu2, beta_md2, beta_mu2, beta_me2, beta_ms2, beta_MassB, beta_MassWB, beta_MassG);
}

const DoubleVector NMSSMSoftPars::display() const
{
   DoubleVector pars(NMSSMSusyPars::display());
   pars.setEnd(numberOfParameters);

   pars(36) = TYu(1,1);
   pars(37) = TYu(1,2);
   pars(38) = TYu(1,3);
   pars(39) = TYu(2,1);
   pars(40) = TYu(2,2);
   pars(41) = TYu(2,3);
   pars(42) = TYu(3,1);
   pars(43) = TYu(3,2);
   pars(44) = TYu(3,3);
   pars(45) = TYd(1,1);
   pars(46) = TYd(1,2);
   pars(47) = TYd(1,3);
   pars(48) = TYd(2,1);
   pars(49) = TYd(2,2);
   pars(50) = TYd(2,3);
   pars(51) = TYd(3,1);
   pars(52) = TYd(3,2);
   pars(53) = TYd(3,3);
   pars(54) = TYe(1,1);
   pars(55) = TYe(1,2);
   pars(56) = TYe(1,3);
   pars(57) = TYe(2,1);
   pars(58) = TYe(2,2);
   pars(59) = TYe(2,3);
   pars(60) = TYe(3,1);
   pars(61) = TYe(3,2);
   pars(62) = TYe(3,3);
   pars(63) = TLambdax;
   pars(64) = TKappa;
   pars(65) = mq2(1,1);
   pars(66) = mq2(1,2);
   pars(67) = mq2(1,3);
   pars(68) = mq2(2,1);
   pars(69) = mq2(2,2);
   pars(70) = mq2(2,3);
   pars(71) = mq2(3,1);
   pars(72) = mq2(3,2);
   pars(73) = mq2(3,3);
   pars(74) = ml2(1,1);
   pars(75) = ml2(1,2);
   pars(76) = ml2(1,3);
   pars(77) = ml2(2,1);
   pars(78) = ml2(2,2);
   pars(79) = ml2(2,3);
   pars(80) = ml2(3,1);
   pars(81) = ml2(3,2);
   pars(82) = ml2(3,3);
   pars(83) = mHd2;
   pars(84) = mHu2;
   pars(85) = md2(1,1);
   pars(86) = md2(1,2);
   pars(87) = md2(1,3);
   pars(88) = md2(2,1);
   pars(89) = md2(2,2);
   pars(90) = md2(2,3);
   pars(91) = md2(3,1);
   pars(92) = md2(3,2);
   pars(93) = md2(3,3);
   pars(94) = mu2(1,1);
   pars(95) = mu2(1,2);
   pars(96) = mu2(1,3);
   pars(97) = mu2(2,1);
   pars(98) = mu2(2,2);
   pars(99) = mu2(2,3);
   pars(100) = mu2(3,1);
   pars(101) = mu2(3,2);
   pars(102) = mu2(3,3);
   pars(103) = me2(1,1);
   pars(104) = me2(1,2);
   pars(105) = me2(1,3);
   pars(106) = me2(2,1);
   pars(107) = me2(2,2);
   pars(108) = me2(2,3);
   pars(109) = me2(3,1);
   pars(110) = me2(3,2);
   pars(111) = me2(3,3);
   pars(112) = ms2;
   pars(113) = MassB;
   pars(114) = MassWB;
   pars(115) = MassG;


   return pars;
}

void NMSSMSoftPars::set(const DoubleVector& v)
{
   NMSSMSusyPars::set(v);

   TYu(1,1) = v(36);
   TYu(1,2) = v(37);
   TYu(1,3) = v(38);
   TYu(2,1) = v(39);
   TYu(2,2) = v(40);
   TYu(2,3) = v(41);
   TYu(3,1) = v(42);
   TYu(3,2) = v(43);
   TYu(3,3) = v(44);
   TYd(1,1) = v(45);
   TYd(1,2) = v(46);
   TYd(1,3) = v(47);
   TYd(2,1) = v(48);
   TYd(2,2) = v(49);
   TYd(2,3) = v(50);
   TYd(3,1) = v(51);
   TYd(3,2) = v(52);
   TYd(3,3) = v(53);
   TYe(1,1) = v(54);
   TYe(1,2) = v(55);
   TYe(1,3) = v(56);
   TYe(2,1) = v(57);
   TYe(2,2) = v(58);
   TYe(2,3) = v(59);
   TYe(3,1) = v(60);
   TYe(3,2) = v(61);
   TYe(3,3) = v(62);
   TLambdax = v(63);
   TKappa = v(64);
   mq2(1,1) = v(65);
   mq2(1,2) = v(66);
   mq2(1,3) = v(67);
   mq2(2,1) = v(68);
   mq2(2,2) = v(69);
   mq2(2,3) = v(70);
   mq2(3,1) = v(71);
   mq2(3,2) = v(72);
   mq2(3,3) = v(73);
   ml2(1,1) = v(74);
   ml2(1,2) = v(75);
   ml2(1,3) = v(76);
   ml2(2,1) = v(77);
   ml2(2,2) = v(78);
   ml2(2,3) = v(79);
   ml2(3,1) = v(80);
   ml2(3,2) = v(81);
   ml2(3,3) = v(82);
   mHd2 = v(83);
   mHu2 = v(84);
   md2(1,1) = v(85);
   md2(1,2) = v(86);
   md2(1,3) = v(87);
   md2(2,1) = v(88);
   md2(2,2) = v(89);
   md2(2,3) = v(90);
   md2(3,1) = v(91);
   md2(3,2) = v(92);
   md2(3,3) = v(93);
   mu2(1,1) = v(94);
   mu2(1,2) = v(95);
   mu2(1,3) = v(96);
   mu2(2,1) = v(97);
   mu2(2,2) = v(98);
   mu2(2,3) = v(99);
   mu2(3,1) = v(100);
   mu2(3,2) = v(101);
   mu2(3,3) = v(102);
   me2(1,1) = v(103);
   me2(1,2) = v(104);
   me2(1,3) = v(105);
   me2(2,1) = v(106);
   me2(2,2) = v(107);
   me2(2,3) = v(108);
   me2(3,1) = v(109);
   me2(3,2) = v(110);
   me2(3,3) = v(111);
   ms2 = v(112);
   MassB = v(113);
   MassWB = v(114);
   MassG = v(115);

}
