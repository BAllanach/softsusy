#include <twoloopmtmb.h>

/* Two-loop O(\alpha_s^2) MSSM corrections to the pole masses of heavy quarks
* A.Bednyakov, A.Onishchenko, V.Velizhanin, O.Veretin
* *
* This file contains the two-loop MSSM corrections to the relation
* between pole and running masses of t-quark (defined in Eq.(57))
* *
* zt2, zt3 - Zeta-functions
* mmsb1 - first  sbottom mass square
* mmsb2 - second sbottom mass square
* mmst1 - first  stop mass square
* mmst2 - second stop mass square
* mt - b-quark mass
* mmt - b-quark mass square
* mgl - gluino mass
* mmgl - gluino mass square
* csb, cs2b, cs4b - cosine mixing angle for sbottom sector
* snb, sn2b, sn4b - sine mixing angle for sbottom sector
* cst, cs2t, cs4t - cosine mixing angle for stop sector
* snt, sn2t, sn4t - sine mixing angle for stop sectort
* mmsusy - two first genration squarks mass square
* mmu - scale square
* *
* functions fin(mm1,mm2) defined as
* if mm1>mm2
*  fin(mm1,mm2)=(-7/2 - (7*mm2)/(2*mm1) +
*   dilog(mm2/mm1) - (mm2*dilog(mm2/mm1))/mm1 -
*   (3*mm2*Log(mm1))/mm1 - Log(mm1)*Log(mm1 - mm2) +
*   (mm2*Log(mm1)*Log(mm1 - mm2))/mm1 + (3*mm2*Log(mm2))/mm1 -
*   Log(mm1)*Log(mm2) + (2*mm2*Log(mm1)*Log(mm2))/mm1 +
*   Log(mm1 - mm2)*Log(mm2) - (mm2*Log(mm1 - mm2)*Log(mm2))/mm1 - PI^2/4 +
*   (mm2*PI^2)/(12*mm1) + Log(mm1)^2 - (3*mm2*Log(mm1)^2)/(2*mm1) -
*   (mm2*Log(mm2)**2)/(2*mm1));
* if mm1<mm2
*  fin(mm1,mm2)=(-7/2 - (7*mm2)/(2*mm1) -
*   dilog(mm1/mm2) + (mm2*dilog(mm1/mm2))/mm1 -
*   (3*mm2*Log(mm1))/mm1 + (3*mm2*Log(mm2))/mm1 +
*   (mm2*Log(mm1)*Log(mm2))/mm1 - Log(mm1)*Log(-mm1 + mm2) +
*   (mm2*Log(mm1)*Log(-mm1 + mm2))/mm1 + Log(mm2)*Log(-mm1 + mm2) -
*   (mm2*Log(mm2)*Log(-mm1 + mm2))/mm1 + PI**2/12 -
*   (mm2*PI**2)/(4*mm1) + Log(mm1)**2/2 - (mm2*Log(mm1)**2)/mm1 -
*   Log(mm2)**2/2);
* if mm1=mm2
*    fin(mm1,mm1)=-7-PI**2/6;
*
* where dilog(x)=dilog(x,2)=Spence(1-x) dilogarifm function */

double MssmSoftsusy::twoLpMt() {
  const double zt2 = sqr(PI) / 6.;
  double mmsb1 = sqr(displayDrBarPars().md(1, 3));
  double mmsb2 = sqr(displayDrBarPars().md(2, 3));
  double mmst1 = sqr(displayDrBarPars().mu(1, 3));
  double mmst2 = sqr(displayDrBarPars().mu(2, 3));
  double mgl = displayGaugino(3);
  double mmgl = sqr(mgl);
  double mt = displayDrBarPars().mt;
  double mmt = sqr(mt);
  double mb = displayDrBarPars().mb;
  double mmb = sqr(mb);
  double csb = cos(displayDrBarPars().thetab), 
    cs2b = cos(displayDrBarPars().thetab * 2.), 
    cs4b = cos(4 * displayDrBarPars().thetab);
  double snb = sin(displayDrBarPars().thetab), 
    sn2b = sin(displayDrBarPars().thetab * 2.), 
    sn4b = sin(4 * displayDrBarPars().thetab);
  double cst = cos(displayDrBarPars().thetat), 
    cs2t = cos(displayDrBarPars().thetat * 2.), 
    cs4t = cos(4 * displayDrBarPars().thetat);
  double snt = sin(displayDrBarPars().thetat), 
    sn2t = sin(displayDrBarPars().thetat * 2.), 
    sn4t = sin(4 * displayDrBarPars().thetat);
  double mmu = sqr(displayMu());
  /// average of first 2 generations squark mass
  double msq = 0.125 * (displayDrBarPars().mu(1, 1) + 
			displayDrBarPars().mu(2, 1) + 
			displayDrBarPars().md(1, 1) + 
			displayDrBarPars().md(2, 1) + 		       
			displayDrBarPars().mu(1, 2) + 
			displayDrBarPars().mu(2, 2) + 
			displayDrBarPars().md(1, 2) + 
			displayDrBarPars().md(2, 2));
  double mmsusy = sqr(msq);

  double lnMglSq = log(mmgl);
  double lnMsbSq = log(mmsb1);
  double lnMsb2Sq = log(mmsb2);
  double lnMst1Sq = log(mmst1);
  double lnMst2Sq = log(mmst2);
  double resmt =

       + sqr(cs2t) * (
          - 640/9
          - 128/9*zt2
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,1)*sn2t * (
          + 32/3*mmsusy/mt*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,1) * (
          + 16/3*mmst1
          - 8*mmsusy
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          + 32/3*mmst1*mmsusy/mt*mgl
          - 32/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,2) * (
          - 56/3*mmst1*mmsusy
          + 56/3*sqr(mmst1)
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,3) * (
          - 32/3*sqr(mmst1)*mmsusy
          + 32/3*pow(mmst1,3)
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,1)*sn2t * (
          - 32/3*mmsusy/mt*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,1) * (
          + 16/3*mmst2
          - 8*mmsusy
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          - 32/3*mmst2*mmsusy/mt*mgl
          + 32/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,2) * (
          - 56/3*mmst2*mmsusy
          + 56/3*sqr(mmst2)
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,3) * (
          - 32/3*sqr(mmst2)*mmsusy
          + 32/3*pow(mmst2,3)
          )

       + fin(mmgl,mmsusy) * (
          - 16/3
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 4/3*mmsb1/mt*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,1) * (
          - 1/3*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb1*mmst1/mt*mgl
          + 4/3*sqr(mmsb1)/mt*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,2) * (
          + mmsb1*sqr(mmst1
          - mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb1*sqr(mmst1)
          - 4/3*sqr(mmsb1)*mmst1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 4/3*mmsb1/mt*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,1) * (
          - 1/3*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb1*mmst2/mt*mgl
          - 4/3*sqr(mmsb1)/mt*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,2) * (
          + mmsb1*sqr(mmst2
          - mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb1*sqr(mmst2)
          - 4/3*sqr(mmsb1)*mmst2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 4/3*mmsb2/mt*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,1) * (
          - 1/3*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb2*mmst1/mt*mgl
          + 4/3*sqr(mmsb2)/mt*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,2) * (
          + mmsb2*sqr(mmst1
          - mmsb2)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb2*sqr(mmst1)
          - 4/3*sqr(mmsb2)*mmst1
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 4/3*mmsb2/mt*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,1) * (
          - 1/3*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb2*mmst2/mt*mgl
          - 4/3*sqr(mmsb2)/mt*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,2) * (
          + mmsb2*sqr(mmst2
          - mmsb2)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb2*sqr(mmst2)
          - 4/3*sqr(mmsb2)*mmst2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 88/9*mmst1/mt*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 5/3*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 154/9*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 34/9*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1) * (
          + 22/9*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,2) * (
          + 12*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,3) * (
          + 16/3*pow(mmst1,3)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 16/9*mmst1/mt*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 5/3*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 26/9*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 34/9*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1) * (
          + 16/9*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmst1*mmst2/mt*mgl
          - 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 26/9*mmst1*mmst2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,2) * (
          - 17/9*mmst1*sqr(mmst2
          - mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,3) * (
          + 4/3*mmst1*sqr(mmst2)
          - 4/3*sqr(mmst1)*mmst2
          )

       + fin(mmst1,mmgl)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 128/9*mmst1
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,1) * (
          - 2/3*mmst1
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb1*mmst1/mt*mgl
          + 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,2) * (
          + mmsb1*mmst1
          - 7/3*sqr(mmst1)
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb1*sqr(mmst1)
          - 4/3*pow(mmst1,3)
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,1) * (
          - 2/3*mmst1
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb2*mmst1/mt*mgl
          + 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,2) * (
          + mmsb2*mmst1
          - 7/3*sqr(mmst1)
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb2*sqr(mmst1)
          - 4/3*pow(mmst1,3)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1)*sn2t * (
          - 4/9*mmst1/mt*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 11/9*mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 8/9*sqr(mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1) * (
          + mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmst1*mmst2/mt*mgl
          + 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 26/9*sqr(mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,2) * (
          + mmst1*mmst2
          + 5/9*sqr(mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,3) * (
          + 4/3*sqr(mmst1)*mmst2
          - 4/3*pow(mmst1,3)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1)*sn2t * (
          + 4/9*mmst1/mt*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 11/9*mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 8/9*sqr(mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1) * (
          + 1/9*mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,2)*sn2t * (
          - 4/3*mmst1*mmst2/mt*mgl
          + 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 26/9*mmst1*mmst2
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,2) * (
          + 5/9*mmst1*sqr(mmst2
          + mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,3) * (
          - 4/3*mmst1*sqr(mmst2)
          + 4/3*sqr(mmst1)*mmst2
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,1) * (
          - 16/3*mmst1
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          - 32/3*mmst1*mmsusy/mt*mgl
          + 32/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,2) * (
          + 8*mmst1*mmsusy
          - 56/3*sqr(mmst1)
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,3) * (
          + 32/3*sqr(mmst1)*mmsusy
          - 32/3*pow(mmst1,3)
          )

       + fin(mmst2,mmgl)*sqr(cs2t) * (
          - 128/9
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 16/9*mmst2/mt*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 26/9*mmst1
          + 11/9*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 26/9*sqr(mmst1)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 34/9*sqr(mmst1)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1) * (
          - 34/9*mmst1
          - 2*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmst1*mmst2/mt*mgl
          + 4/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 26/9*mmst1*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,2) * (
          - 17/9*mmst1*sqr(mmst2
          - mmst2)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,3) * (
          - 4/3*mmst1*sqr(mmst2)
          + 4/3*sqr(mmst1)*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 88/9*mmst2/mt*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 154/9*mmst1
          - 139/9*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 154/9*sqr(mmst1)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 34/9*sqr(mmst1)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1) * (
          + 34/9*mmst1
          + 56/9*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,2) * (
          + 12*sqr(mmst2)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,3) * (
          + 16/3*pow(mmst2,3)
          )

       + fin(mmst2,mmgl)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 128/9*mmst1
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,1) * (
          - 2/3*mmst2
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb1*mmst2/mt*mgl
          - 4/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,2) * (
          + mmsb1*mmst2
          - 7/3*sqr(mmst2)
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb1*sqr(mmst2)
          - 4/3*pow(mmst2,3)
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,1) * (
          - 2/3*mmst2
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb2*mmst2/mt*mgl
          - 4/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,2) * (
          + mmsb2*mmst2
          - 7/3*sqr(mmst2)
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb2*sqr(mmst2)
          - 4/3*pow(mmst2,3)
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,1) * (
          - 16/3*mmst2
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          + 32/3*mmst2*mmsusy/mt*mgl
          - 32/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,2) * (
          + 8*mmst2*mmsusy
          - 56/3*sqr(mmst2)
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,3) * (
          + 32/3*sqr(mmst2)*mmsusy
          - 32/3*pow(mmst2,3)
          )

       + lnMglSq*sqr(cs2t) * (
          + 128/3
          )

       + sqr(lnMglSq)*sqr(cs2t) * (
          - 64/9
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*sn2t * (
          + 2/3*mmsb1/mt*mgl
          + 2/3*mmsb2/mt*mgl
          + 10*mmst1/mt*mgl
          + 2/3*mmst2/mt*mgl
          + 16/3*mmsusy/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 53/3*mmst1
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 223/9*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 43/9*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t * (
          - 16/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          + 32/9*pow(mmst1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          - 16/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1) * (
          - 1/2*mmsb1
          - 1/2*mmsb2
          - 73/6*mmst1
          - 1/2*mmst2
          + 4/3*mmsusy
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2)*sn2t * (
          + 2/3*mmsb1*mmst1/mt*mgl
          + 2/3*mmsb2*mmst1/mt*mgl
          + 2/3*mmst1*mmst2/mt*mgl
          - 16*mmst1*mmsusy/mt*mgl
          - 122/9*sqr(mmst1)/mt*mgl
          + 32/3*sqr(mmsusy)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 53/6*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2) * (
          - 7/6*mmsb1*mmst1
          - 7/6*mmsb2*mmst1
          - 7/6*mmst1*mmst2
          + 52/3*mmst1*mmsusy
          + 151/9*sqr(mmst1)
          - 8*sqr(mmsusy)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,3)*sn2t * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,3) * (
          - 2/3*mmsb1*sqr(mmst1)
          - 2/3*mmsb2*sqr(mmst1)
          - 32/3*mmst1*sqr(mmsusy)
          - 2/3*sqr(mmst1)*mmst2
          + 16*sqr(mmst1)*mmsusy
          + 46/3*pow(mmst1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,4) * (
          + 4/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmsb1/mt*mgl
          - 2/3*mmsb2/mt*mgl
          + 2/9*mmst1/mt*mgl
          - 98/9*mmst2/mt*mgl
          - 16/3*mmsusy/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 223/9*mmst1
          - 64/9*mmst2
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t * (
          - 8/3*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 223/9*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 43/9*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          - 32/9*pow(mmst1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          + 16/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1) * (
          - 1/2*mmsb1
          - 1/2*mmsb2
          + 109/18*mmst1
          - 101/18*mmst2
          + 4/3*mmsusy
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*sn2t * (
          - 2/3*mmsb1*mmst2/mt*mgl
          - 2/3*mmsb2*mmst2/mt*mgl
          - 14/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          + 16*mmst2*mmsusy/mt*mgl
          + 38/3*sqr(mmst2)/mt*mgl
          - 32/3*sqr(mmsusy)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 53/6*sqr(mmst2)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          - 32/9*pow(mmst1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2) * (
          - 7/6*mmsb1*mmst2
          - 7/6*mmsb2*mmst2
          + 11/18*mmst1*mmst2
          + 8/3*sqr(mmst1)
          + 52/3*mmst2*mmsusy
          + 53/3*sqr(mmst2)
          - 8*sqr(mmsusy)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,3)*sn2t * (
          + 8/9*pow(mmst2,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,3) * (
          - 2/3*mmsb1*sqr(mmst2)
          - 2/3*mmsb2*sqr(mmst2)
          - 2/3*mmst1*sqr(mmst2)
          - 32/3*mmst2*sqr(mmsusy)
          + 16*sqr(mmst2)*mmsusy
          + 46/3*pow(mmst2,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,4) * (
          + 4/9*pow(mmst2,4)
          )

       + sqr(lnMglSq) * (
          - 166/9
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,1)*sn2t * (
          - 4/3*mmsb1/mt*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,1) * (
          + mmsb1
          - 4/3*mmst1
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb1*mmst1/mt*mgl
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,2) * (
          + 7/3*mmsb1*mmst1
          - 14/3*sqr(mmst1)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,3) * (
          + 4/3*mmsb1*sqr(mmst1)
          - 8/3*pow(mmst1,3)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,1)*sn2t * (
          + 4/3*mmsb1/mt*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,1) * (
          + mmsb1
          - 4/3*mmst2
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb1*mmst2/mt*mgl
          - 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,2) * (
          + 7/3*mmsb1*mmst2
          - 14/3*sqr(mmst2)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,3) * (
          + 4/3*mmsb1*sqr(mmst2)
          - 8/3*pow(mmst2,3)
          )

       + lnMglSq*lnMsbSq * (
          + 4/3
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,1)*sn2t * (
          - 4/3*mmsb2/mt*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,1) * (
          + mmsb2
          - 4/3*mmst1
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb2*mmst1/mt*mgl
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,2) * (
          + 7/3*mmsb2*mmst1
          - 14/3*sqr(mmst1)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,3) * (
          + 4/3*mmsb2*sqr(mmst1)
          - 8/3*pow(mmst1,3)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,1)*sn2t * (
          + 4/3*mmsb2/mt*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,1) * (
          + mmsb2
          - 4/3*mmst2
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb2*mmst2/mt*mgl
          - 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,2) * (
          + 7/3*mmsb2*mmst2
          - 14/3*sqr(mmst2)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,3) * (
          + 4/3*mmsb2*sqr(mmst2)
          - 8/3*pow(mmst2,3)
          )

       + lnMglSq*lnMsb2Sq * (
          + 4/3
          )

       + lnMglSq*lnMst1Sq*sn2t * (
          - 208/9/mt*mgl
          )

       + lnMglSq*lnMst1Sq*sqr(cs2t) * (
          - 64/9
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*sn2t * (
          - 28/3*mmst1/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 121/9*mmst1
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 287/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 43/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          - 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          + 16/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1) * (
          - 34/3*mmst1
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*sn2t * (
          + 172/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 11/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2) * (
          - 130/3*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,3)*sn2t * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          + 16/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,3) * (
          - 212/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,4) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*sn2t * (
          + 20/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 5/3*mmst1
          - 22/9*mmst2
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          - 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 31/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 43/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          - 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          - 16/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1) * (
          - 34/9*mmst1
          + 2/9*mmst2
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*sn2t * (
          + 4*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          - 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          - 52/9*sqr(mmst2)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2) * (
          + 37/9*mmst1*mmst2
          - 8/3*sqr(mmst1)
          + 10/9*sqr(mmst2)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          - 16/9*mmst1*sqr(mmst2)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,3) * (
          + 28/9*mmst1*sqr(mmst2)
          - 8/3*pow(mmst2,3)
          )

       + lnMglSq*lnMst1Sq*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 256/9*mmgl
          + 256/9*mmst1
          )

       + lnMglSq*lnMst1Sq * (
          + 52/9
          )

       + lnMglSq*lnMst2Sq*sn2t * (
          + 208/9/mt*mgl
          )

       + lnMglSq*lnMst2Sq*sqr(cs2t) * (
          + 64/3
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*sn2t * (
          - 28/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 16/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 53/9*mmst1
          - 16/9*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          - 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 31/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 43/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          - 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          + 16/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1) * (
          + 61/9*mmst1
          + 25/9*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*sn2t * (
          - 28/9*mmst1*mmst2/mt*mgl
          + 32/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          - 52/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2) * (
          + 53/9*mmst1*mmst2
          + 2*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          - 16/9*sqr(mmst1)*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,3) * (
          + 28/9*sqr(mmst1)*mmst2
          - 8/3*pow(mmst1,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*sn2t * (
          - 8/9*mmst1/mt*mgl
          + 92/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 16/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 287/9*mmst1
          + 166/9*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 287/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 43/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          - 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          - 16/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1) * (
          - 59/9*mmst1
          - 161/9*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          - 164/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 11/9*sqr(mmst2)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2) * (
          - 16/9*mmst1*mmst2
          - 8/3*sqr(mmst1)
          - 398/9*sqr(mmst2)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,3)*sn2t * (
          - 16/9*pow(mmst2,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          + 16/9*pow(mmst2,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,3) * (
          - 212/9*pow(mmst2,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,4) * (
          - 8/9*pow(mmst2,4)
          )

       + lnMglSq*lnMst2Sq*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 256/9*mmgl
          - 256/9*mmst1
          )

       + lnMglSq*lnMst2Sq * (
          + 52/9
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst1,1)*sn2t * (
          - 32/3*mmsusy/mt*mgl
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst1,1) * (
          - 8/3*mmsusy
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          + 32*mmst1*mmsusy/mt*mgl
          - 64/3*sqr(mmsusy)/mt*mgl
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst1,2) * (
          - 104/3*mmst1*mmsusy
          + 16*sqr(mmsusy)
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst1,3) * (
          + 64/3*mmst1*sqr(mmsusy)
          - 32*sqr(mmst1)*mmsusy
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst2,1)*sn2t * (
          + 32/3*mmsusy/mt*mgl
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst2,1) * (
          - 8/3*mmsusy
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          - 32*mmst2*mmsusy/mt*mgl
          + 64/3*sqr(mmsusy)/mt*mgl
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst2,2) * (
          - 104/3*mmst2*mmsusy
          + 16*sqr(mmsusy)
          )

       + lnMglSq*Log(mmsusy)*den(mmgl - mmst2,3) * (
          + 64/3*mmst2*sqr(mmsusy)
          - 32*sqr(mmst2)*mmsusy
          )

       + lnMglSq*Log(mmt)*den(mmgl - mmst1,1)*sn2t * (
          + 16/3*mmst1/mt*mgl
          )

       + lnMglSq*Log(mmt)*den(mmgl - mmst1,1) * (
          - 16/3*mmst1
          )

       + lnMglSq*Log(mmt)*den(mmgl - mmst1,2) * (
          - 8/3*sqr(mmst1)
          )

       + lnMglSq*Log(mmt)*den(mmgl - mmst2,1)*sn2t * (
          - 16/3*mmst2/mt*mgl
          )

       + lnMglSq*Log(mmt)*den(mmgl - mmst2,1) * (
          - 16/3*mmst2
          )

       + lnMglSq*Log(mmt)*den(mmgl - mmst2,2) * (
          - 8/3*sqr(mmst2)
          )

       + lnMglSq*Log(mmt) * (
          - 40/3
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,1)*sn2t * (
          - 16*mmst1/mt*mgl
          + 16/9*mmst2/mt*mgl
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 16*mmst1
          + 16/9*mmst2
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 128/9*sqr(mmst1)
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,1) * (
          + 332/9*mmst1
          - 16/9*mmst2
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,2)*sn2t * (
          + 16/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          - 32/3*sqr(mmst1)
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,2) * (
          - 32/9*mmst1*mmst2
          + 178/9*sqr(mmst1)
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,3)*sqr(cs2t) * (
          + 16/9*sqr(mmst1)*mmst2
          - 16/9*pow(mmst1,3)
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst1,3) * (
          - 16/9*sqr(mmst1)*mmst2
          + 8/9*pow(mmst1,3)
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,1)*sn2t * (
          - 16/9*mmst1/mt*mgl
          + 16*mmst2/mt*mgl
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 16*mmst1
          - 16/9*mmst2
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 128/9*sqr(mmst1)
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,1) * (
          - 16/9*mmst1
          + 332/9*mmst2
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,2)*sn2t * (
          - 16/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          - 32/3*sqr(mmst2)
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,2) * (
          - 32/9*mmst1*mmst2
          + 178/9*sqr(mmst2)
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,3)*sqr(cs2t) * (
          + 16/9*mmst1*sqr(mmst2)
          - 16/9*pow(mmst2,3)
          )

       + lnMglSq*Log(mmu)*den(mmgl - mmst2,3) * (
          - 16/9*mmst1*sqr(mmst2)
          + 8/9*pow(mmst2,3)
          )

       + lnMglSq*Log(mmu) * (
          + 36
          )

       + lnMglSq*den(mmgl - mmst1,1)*sn2t * (
          - 8/3*mmsb1/mt*mgl
          - 8/3*mmsb2/mt*mgl
          - 176/3*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          + 128/3*mmsusy/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 428/9*mmst1
          + 16/9*mmst2
          )

       + lnMglSq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 796/9*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 76/3*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst1,1) * (
          + 2*mmsb1
          + 2*mmsb2
          + 290/9*mmst1
          + 2/9*mmst2
          - 16*mmsusy
          )

       + lnMglSq*den(mmgl - mmst1,2)*sn2t * (
          - 8/3*mmsb1*mmst1/mt*mgl
          - 8/3*mmsb2*mmst1/mt*mgl
          - 8/9*mmst1*mmst2/mt*mgl
          - 64/3*mmst1*mmsusy/mt*mgl
          + 56/9*sqr(mmst1)/mt*mgl
          + 32*sqr(mmsusy)/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          - 238/9*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1) * (
          - 8/9*pow(mmst1,3)
          )

       + lnMglSq*den(mmgl - mmst1,2) * (
          + 14/3*mmsb1*mmst1
          + 14/3*mmsb2*mmst1
          + 10/9*mmst1*mmst2
          + 16/3*mmst1*mmsusy
          - 290/9*sqr(mmst1)
          - 24*sqr(mmsusy)
          )

       + lnMglSq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          + 16/9*sqr(mmst1)*mmst2
          - 16/9*pow(mmst1,3)
          )

       + lnMglSq*den(mmgl - mmst1,3) * (
          + 8/3*mmsb1*sqr(mmst1)
          + 8/3*mmsb2*sqr(mmst1)
          - 32*mmst1*sqr(mmsusy)
          + 8/9*sqr(mmst1)*mmst2
          + 64/3*sqr(mmst1)*mmsusy
          - 200/9*pow(mmst1,3)
          )

       + lnMglSq*den(mmgl - mmst2,1)*sn2t * (
          + 8/3*mmsb1/mt*mgl
          + 8/3*mmsb2/mt*mgl
          + 8/9*mmst1/mt*mgl
          + 176/3*mmst2/mt*mgl
          - 128/3*mmsusy/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 812/9*mmst1
          + 368/9*mmst2
          )

       + lnMglSq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 796/9*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 76/3*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst2,1) * (
          + 2*mmsb1
          + 2*mmsb2
          - 226/9*mmst1
          + 62/9*mmst2
          - 16*mmsusy
          )

       + lnMglSq*den(mmgl - mmst2,2)*sn2t * (
          + 8/3*mmsb1*mmst2/mt*mgl
          + 8/3*mmsb2*mmst2/mt*mgl
          + 8/9*mmst1*mmst2/mt*mgl
          + 64/3*mmst2*mmsusy/mt*mgl
          - 56/9*sqr(mmst2)/mt*mgl
          - 32*sqr(mmsusy)/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          - 238/9*sqr(mmst2)
          )

       + lnMglSq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          + 8/9*pow(mmst1,3)
          )

       + lnMglSq*den(mmgl - mmst2,2) * (
          + 14/3*mmsb1*mmst2
          + 14/3*mmsb2*mmst2
          + 2/9*mmst1*mmst2
          - 8/9*sqr(mmst1)
          + 16/3*mmst2*mmsusy
          - 298/9*sqr(mmst2)
          - 24*sqr(mmsusy)
          )

       + lnMglSq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          + 16/9*mmst1*sqr(mmst2)
          - 16/9*pow(mmst2,3)
          )

       + lnMglSq*den(mmgl - mmst2,3) * (
          + 8/3*mmsb1*sqr(mmst2)
          + 8/3*mmsb2*sqr(mmst2)
          + 8/9*mmst1*sqr(mmst2)
          - 32*mmst2*sqr(mmsusy)
          + 64/3*sqr(mmst2)*mmsusy
          - 200/9*pow(mmst2,3)
          )

       + lnMglSq * (
          + 232/3
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,1)*sn2t * (
          + 2/3*mmsb1/mt*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,1) * (
          - 7/6*mmsb1
          + 2/3*mmst1
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmsb1*mmst1/mt*mgl
          - 4/3*sqr(mmsb1)/mt*mgl
          - 4/3*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,2) * (
          - 4*mmsb1*sqr(mmst1
          + mmsb1)
          + 7/3*sqr(mmst1)
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,3) * (
          - 8/3*mmsb1*sqr(mmst1)
          + 4/3*sqr(mmsb1)*mmst1
          + 4/3*pow(mmst1,3)
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmsb1/mt*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,1) * (
          - 7/6*mmsb1
          + 2/3*mmst2
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,2)*sn2t * (
          - 8/3*mmsb1*mmst2/mt*mgl
          + 4/3*sqr(mmsb1)/mt*mgl
          + 4/3*sqr(mmst2)/mt*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,2) * (
          - 4*mmsb1*sqr(mmst2
          + mmsb1)
          + 7/3*sqr(mmst2)
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,3) * (
          - 8/3*mmsb1*sqr(mmst2)
          + 4/3*sqr(mmsb1)*mmst2
          + 4/3*pow(mmst2,3)
          )

       + sqr(lnMsbSq) * (
          - 1/3
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmst1,1) * (
          + 4/3*mmsb1
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmst1,2)*sn2t * (
          - 4*mmsb1*mmst1/mt*mgl
          + 8/3*sqr(mmsb1)/mt*mgl
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmst1,2) * (
          + 17/3*mmsb1*mmst1
          - 2*sqr(mmsb1)
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmst1,3) * (
          + 4*mmsb1*sqr(mmst1)
          - 8/3*sqr(mmsb1)*mmst1
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmst2,1) * (
          + 4/3*mmsb1
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          + 4*mmsb1*mmst2/mt*mgl
          - 8/3*sqr(mmsb1)/mt*mgl
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmst2,2) * (
          + 17/3*mmsb1*mmst2
          - 2*sqr(mmsb1)
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmst2,3) * (
          + 4*mmsb1*sqr(mmst2)
          - 8/3*sqr(mmsb1)*mmst2
          )

       + lnMsbSq*Log(mmt) * (
          - 2/3
          )

       + lnMsbSq*den(mmgl - mmst1,1)*sn2t * (
          + 8/3*mmsb1/mt*mgl
          )

       + lnMsbSq*den(mmgl - mmst1,1) * (
          + 2*mmst1
          )

       + lnMsbSq*den(mmgl - mmst1,2)*sn2t * (
          + 4*sqr(mmsb1)/mt*mgl
          - 4*sqr(mmst1)/mt*mgl
          )

       + lnMsbSq*den(mmgl - mmst1,2) * (
          + 4/3*mmsb1*mmst1
          - 3*sqr(mmsb1)
          + 7*sqr(mmst1)
          )

       + lnMsbSq*den(mmgl - mmst1,3) * (
          - 4*sqr(mmsb1)*mmst1
          + 4*pow(mmst1,3)
          )

       + lnMsbSq*den(mmgl - mmst2,1)*sn2t * (
          - 8/3*mmsb1/mt*mgl
          )

       + lnMsbSq*den(mmgl - mmst2,1) * (
          + 2*mmst2
          )

       + lnMsbSq*den(mmgl - mmst2,2)*sn2t * (
          - 4*sqr(mmsb1)/mt*mgl
          + 4*sqr(mmst2)/mt*mgl
          )

       + lnMsbSq*den(mmgl - mmst2,2) * (
          + 4/3*mmsb1*mmst2
          - 3*sqr(mmsb1)
          + 7*sqr(mmst2)
          )

       + lnMsbSq*den(mmgl - mmst2,3) * (
          - 4*sqr(mmsb1)*mmst2
          + 4*pow(mmst2,3)
          )

       + lnMsbSq * (
          + 1/9
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,1)*sn2t * (
          + 2/3*mmsb2/mt*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,1) * (
          - 7/6*mmsb2
          + 2/3*mmst1
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmsb2*mmst1/mt*mgl
          - 4/3*sqr(mmsb2)/mt*mgl
          - 4/3*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,2) * (
          - 4*mmsb2*sqr(mmst1
          + mmsb2)
          + 7/3*sqr(mmst1)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,3) * (
          - 8/3*mmsb2*sqr(mmst1)
          + 4/3*sqr(mmsb2)*mmst1
          + 4/3*pow(mmst1,3)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmsb2/mt*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,1) * (
          - 7/6*mmsb2
          + 2/3*mmst2
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,2)*sn2t * (
          - 8/3*mmsb2*mmst2/mt*mgl
          + 4/3*sqr(mmsb2)/mt*mgl
          + 4/3*sqr(mmst2)/mt*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,2) * (
          - 4*mmsb2*sqr(mmst2
          + mmsb2)
          + 7/3*sqr(mmst2)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,3) * (
          - 8/3*mmsb2*sqr(mmst2)
          + 4/3*sqr(mmsb2)*mmst2
          + 4/3*pow(mmst2,3)
          )

       + sqr(lnMsb2Sq) * (
          - 1/3
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmst1,1) * (
          + 4/3*mmsb2
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmst1,2)*sn2t * (
          - 4*mmsb2*mmst1/mt*mgl
          + 8/3*sqr(mmsb2)/mt*mgl
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmst1,2) * (
          + 17/3*mmsb2*mmst1
          - 2*sqr(mmsb2)
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmst1,3) * (
          + 4*mmsb2*sqr(mmst1)
          - 8/3*sqr(mmsb2)*mmst1
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmst2,1) * (
          + 4/3*mmsb2
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          + 4*mmsb2*mmst2/mt*mgl
          - 8/3*sqr(mmsb2)/mt*mgl
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmst2,2) * (
          + 17/3*mmsb2*mmst2
          - 2*sqr(mmsb2)
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmst2,3) * (
          + 4*mmsb2*sqr(mmst2)
          - 8/3*sqr(mmsb2)*mmst2
          )

       + lnMsb2Sq*Log(mmt) * (
          - 2/3
          )

       + lnMsb2Sq*den(mmgl - mmst1,1)*sn2t * (
          + 8/3*mmsb2/mt*mgl
          )

       + lnMsb2Sq*den(mmgl - mmst1,1) * (
          + 2*mmst1
          )

       + lnMsb2Sq*den(mmgl - mmst1,2)*sn2t * (
          + 4*sqr(mmsb2)/mt*mgl
          - 4*sqr(mmst1)/mt*mgl
          )

       + lnMsb2Sq*den(mmgl - mmst1,2) * (
          + 4/3*mmsb2*mmst1
          - 3*sqr(mmsb2)
          + 7*sqr(mmst1)
          )

       + lnMsb2Sq*den(mmgl - mmst1,3) * (
          - 4*sqr(mmsb2)*mmst1
          + 4*pow(mmst1,3)
          )

       + lnMsb2Sq*den(mmgl - mmst2,1)*sn2t * (
          - 8/3*mmsb2/mt*mgl
          )

       + lnMsb2Sq*den(mmgl - mmst2,1) * (
          + 2*mmst2
          )

       + lnMsb2Sq*den(mmgl - mmst2,2)*sn2t * (
          - 4*sqr(mmsb2)/mt*mgl
          + 4*sqr(mmst2)/mt*mgl
          )

       + lnMsb2Sq*den(mmgl - mmst2,2) * (
          + 4/3*mmsb2*mmst2
          - 3*sqr(mmsb2)
          + 7*sqr(mmst2)
          )

       + lnMsb2Sq*den(mmgl - mmst2,3) * (
          - 4*sqr(mmsb2)*mmst2
          + 4*pow(mmst2,3)
          )

       + lnMsb2Sq * (
          + 1/9
          )

       + lnMst1Sq*sn2t * (
          + 280/9/mt*mgl
          )

       + lnMst1Sq*sqr(cs2t) * (
          + 64/9
          )

       + sqr(lnMst1Sq)*sn2t * (
          + 8/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*sn2t * (
          - 2/9*mmst1/mt*mgl
          - 4/9*mmst2/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 14/9*mmst1
          - 11/9*mmst2
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*cs2t * (
          - 16/9*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 77/9*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 13/9*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1) * (
          - 2/3*mmsb1
          - 2/3*mmsb2
          + 419/18*mmst1
          + mmst2
          - 16/3*mmsusy
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,2)*sn2t * (
          + 2*mmsb1*mmst1/mt*mgl
          - 4/3*sqr(mmsb1)/mt*mgl
          + 2*mmsb2*mmst1/mt*mgl
          - 4/3*sqr(mmsb2)/mt*mgl
          + 2*mmst1*mmst2/mt*mgl
          + 16*mmst1*mmsusy/mt*mgl
          - 86/9*sqr(mmst1)/mt*mgl
          - 4/3*sqr(mmst2)/mt*mgl
          - 32/3*sqr(mmsusy)/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 26/9*mmst1*mmst2
          - 85/18*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,2) * (
          - 17/6*mmsb1*sqr(mmst1
          + mmsb1)
          - 17/6*mmsb2*sqr(mmst1
          + mmsb2)
          + 1/18*mmst1*mmst2
          - 68/3*mmst1*mmsusy
          + 92/3*sqr(sqr(mmst1)
          + mmst2)
          + 8*sqr(mmsusy)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,3)*sn2t * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,3)*sqr(cs2t) * (
          - 16/9*pow(mmst1,3)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,3) * (
          - 2*mmsb1*sqr(mmst1)
          + 4/3*sqr(mmsb1)*mmst1
          - 2*mmsb2*sqr(mmst1)
          + 4/3*sqr(mmsb2)*mmst1
          + 4/3*mmst1*sqr(mmst2)
          + 32/3*mmst1*sqr(mmsusy)
          - 2*sqr(mmst1)*mmst2
          - 16*sqr(mmst1)*mmsusy
          + 110/9*pow(mmst1,3)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,4) * (
          + 4/9*pow(mmst1,4)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmst1/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 13/9*mmst1
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 13/9*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 13/9*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1) * (
          + 17/18*mmst1
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,2) * (
          - 2/3*mmst1*mmst2
          )

       + sqr(lnMst1Sq)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 128/9*mmgl
          - 64/9*mmst1
          )

       + sqr(lnMst1Sq) * (
          + 41/9
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*sn2t * (
          - 8/9*mmst1/mt*mgl
          + 8/3*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 16/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 5/9*mmst1
          + 38/9*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 5/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 17/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          - 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          - 16/9*pow(mmst1,4)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1) * (
          - 11/3*mmst1
          - 34/9*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*sn2t * (
          - 20/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          + 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 28/3*mmst1*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,4)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2) * (
          - 11/3*mmst1*mmst2
          - 8/9*sqr(mmst1)
          - 2*sqr(mmst2)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          + 16/9*sqr(mmst1)*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,3) * (
          - 8/3*mmst1*sqr(mmst2)
          + 20/9*sqr(mmst1)*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*sn2t * (
          - 8/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 11/9*mmst1
          + 22/9*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 5/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 17/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          - 32/9*pow(mmst1,3)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          + 16/9*pow(mmst1,4)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1) * (
          + 17/9*mmst1
          - 2/9*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          - 4*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          + 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          + 52/9*sqr(mmst2)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          - 32/9*pow(mmst1,3)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,4)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2) * (
          - 25/9*mmst1*mmst2
          + 8/3*sqr(mmst1)
          - 10/9*sqr(mmst2)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          + 16/9*mmst1*sqr(mmst2)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,3) * (
          - 28/9*mmst1*sqr(mmst2)
          + 8/3*pow(mmst2,3)
          )

       + lnMst1Sq*Log(mmsusy)*den(mmgl - mmst1,1) * (
          + 32/3*mmsusy
          )

       + lnMst1Sq*Log(mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          - 32*mmst1*mmsusy/mt*mgl
          + 64/3*sqr(mmsusy)/mt*mgl
          )

       + lnMst1Sq*Log(mmsusy)*den(mmgl - mmst1,2) * (
          + 136/3*mmst1*mmsusy
          - 16*sqr(mmsusy)
          )

       + lnMst1Sq*Log(mmsusy)*den(mmgl - mmst1,3) * (
          - 64/3*mmst1*sqr(mmsusy)
          + 32*sqr(mmst1)*mmsusy
          )

       + lnMst1Sq*Log(mmt)*den(mmgl - mmst1,1)*sn2t * (
          - 16/3*mmst1/mt*mgl
          )

       + lnMst1Sq*Log(mmt)*den(mmgl - mmst1,1) * (
          + 16/3*mmst1
          )

       + lnMst1Sq*Log(mmt)*den(mmgl - mmst1,2) * (
          + 8/3*sqr(mmst1)
          )

       + lnMst1Sq*Log(mmt) * (
          - 2/3
          )

       + lnMst1Sq*Log(mmu)*sn2t * (
          + 64/9/mt*mgl
          )

       + lnMst1Sq*Log(mmu)*sqr(cs2t) * (
          + 64/9
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,1)*sn2t * (
          + 16*mmst1/mt*mgl
          - 16/9*mmst2/mt*mgl
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 16*mmst1
          - 16/9*mmst2
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 128/9*sqr(mmst1)
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,1) * (
          - 332/9*mmst1
          + 16/9*mmst2
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,2)*sn2t * (
          - 16/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          + 32/3*sqr(mmst1)
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,2) * (
          + 32/9*mmst1*mmst2
          - 178/9*sqr(mmst1)
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,3)*sqr(cs2t) * (
          - 16/9*sqr(mmst1)*mmst2
          + 16/9*pow(mmst1,3)
          )

       + lnMst1Sq*Log(mmu)*den(mmgl - mmst1,3) * (
          + 16/9*sqr(mmst1)*mmst2
          - 8/9*pow(mmst1,3)
          )

       + lnMst1Sq*Log(mmu)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 128/9*mmst1
          )

       + lnMst1Sq*Log(mmu) * (
          - 128/9
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*sn2t * (
          + 172/3*mmst1/mt*mgl
          - 28/9*mmst2/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 16/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 229/9*mmst1
          - 49/9*mmst2
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t * (
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 718/9*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 34/3*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst1,1) * (
          - 2*mmsb1
          - 2*mmsb2
          - 13/3*mmst1
          + 43/9*mmst2
          - 16*mmsusy
          )

       + lnMst1Sq*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmsb1*mmst1/mt*mgl
          - 4*sqr(mmsb1)/mt*mgl
          + 8/3*mmsb2*mmst1/mt*mgl
          - 4*sqr(mmsb2)/mt*mgl
          + 8/9*mmst1*mmst2/mt*mgl
          + 64/3*mmst1*mmsusy/mt*mgl
          + 52/9*sqr(mmst1)/mt*mgl
          - 4*sqr(mmst2)/mt*mgl
          - 32*sqr(mmsusy)/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 110/9*mmst1*mmst2
          + 16*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1) * (
          + 8/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst1,2) * (
          - 6*mmsb1*mmst1
          + 3*sqr(mmsb1)
          - 6*mmsb2*mmst1
          + 3*sqr(mmsb2)
          + 56/9*mmst1*mmst2
          - 48*mmst1*mmsusy
          + 187/9*sqr(mmst1)
          + 3*sqr(mmst2)
          + 24*sqr(mmsusy)
          )

       + lnMst1Sq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          - 16/9*sqr(mmst1)*mmst2
          + 16/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst1,3) * (
          - 8/3*mmsb1*sqr(mmst1)
          + 4*sqr(mmsb1)*mmst1
          - 8/3*mmsb2*sqr(mmst1)
          + 4*sqr(mmsb2)*mmst1
          + 4*mmst1*sqr(mmst2)
          + 32*mmst1*sqr(mmsusy)
          - 8/9*sqr(mmst1)*mmst2
          - 64/3*sqr(mmst1)*mmsusy
          + 92/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*sn2t * (
          - 16/3*mmst1/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 6*mmst1
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 26/3*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 34/3*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst2,1) * (
          + 52/9*mmst1
          )

       + lnMst1Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 16/9*mmst1*mmst2
          )

       + lnMst1Sq*den(mmgl - mmst2,2) * (
          - 40/9*mmst1*mmst2
          )

       + lnMst1Sq*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 128/3*mmgl
          - 640/9*mmst1
          )

       + lnMst1Sq * (
          + 5/9
          )

       + lnMst2Sq*sn2t * (
          - 280/9/mt*mgl
          )

       + lnMst2Sq*sqr(cs2t) * (
          - 64
          )

       + sqr(lnMst2Sq)*sn2t * (
          - 8/mt*mgl
          )

       + sqr(lnMst2Sq)*sqr(cs2t) * (
          - 64/9
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1)*sn2t * (
          + 4/9*mmst1/mt*mgl
          + 2/9*mmst2/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 8/3*mmst1
          - 11/9*mmst2
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 13/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 13/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1) * (
          - 14/9*mmst1
          + 1/2*mmst2
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmst1*mmst2/mt*mgl
          - 4/3*sqr(mmst1)/mt*mgl
          - 4/3*sqr(mmst2)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 26/9*mmst1*mmst2
          + 26/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,2) * (
          - 10/9*mmst1*mmst2
          - 5/9*sqr(sqr(mmst1)
          + mmst2)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,3) * (
          + 4/3*mmst1*sqr(mmst2)
          - 8/3*sqr(mmst1)*mmst2
          + 4/3*pow(mmst1,3)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*sn2t * (
          + 2/3*mmst2/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 16/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 77/9*mmst1
          - 34/3*mmst2
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*cs2t * (
          - 16/9*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 77/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 13/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1) * (
          - 2/3*mmsb1
          - 2/3*mmsb2
          + 13/9*mmst1
          + 149/6*mmst2
          - 16/3*mmsusy
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,2)*sn2t * (
          - 2*mmsb1*mmst2/mt*mgl
          + 4/3*sqr(mmsb1)/mt*mgl
          - 2*mmsb2*mmst2/mt*mgl
          + 4/3*sqr(mmsb2)/mt*mgl
          + 2/3*mmst1*mmst2/mt*mgl
          - 16*mmst2*mmsusy/mt*mgl
          + 74/9*sqr(mmst2)/mt*mgl
          + 32/3*sqr(mmsusy)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*sqr(mmst2)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 137/18*sqr(mmst2)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,2) * (
          - 17/6*mmsb1*sqr(mmst2
          + mmsb1)
          - 17/6*mmsb2*sqr(mmst2
          + mmsb2)
          + 1/2*mmst1*mmst2
          - 68/3*mmst2*mmsusy
          + 281/9*sqr(mmst2)
          + 8*sqr(mmsusy)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,3)*sn2t * (
          + 8/9*pow(mmst2,3)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,3)*sqr(cs2t) * (
          - 16/9*pow(mmst2,3)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,3) * (
          - 2*mmsb1*sqr(mmst2)
          + 4/3*sqr(mmsb1)*mmst2
          - 2*mmsb2*sqr(mmst2)
          + 4/3*sqr(mmsb2)*mmst2
          + 2/3*mmst1*sqr(mmst2)
          + 32/3*mmst2*sqr(mmsusy)
          - 16*sqr(mmst2)*mmsusy
          + 98/9*pow(mmst2,3)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,4) * (
          + 4/9*pow(mmst2,4)
          )

       + sqr(lnMst2Sq)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 128/9*mmgl
          + 64/9*mmst1
          )

       + sqr(lnMst2Sq) * (
          + 41/9
          )

       + lnMst2Sq*Log(mmsusy)*den(mmgl - mmst2,1) * (
          + 32/3*mmsusy
          )

       + lnMst2Sq*Log(mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          + 32*mmst2*mmsusy/mt*mgl
          - 64/3*sqr(mmsusy)/mt*mgl
          )

       + lnMst2Sq*Log(mmsusy)*den(mmgl - mmst2,2) * (
          + 136/3*mmst2*mmsusy
          - 16*sqr(mmsusy)
          )

       + lnMst2Sq*Log(mmsusy)*den(mmgl - mmst2,3) * (
          - 64/3*mmst2*sqr(mmsusy)
          + 32*sqr(mmst2)*mmsusy
          )

       + lnMst2Sq*Log(mmt)*den(mmgl - mmst2,1)*sn2t * (
          + 16/3*mmst2/mt*mgl
          )

       + lnMst2Sq*Log(mmt)*den(mmgl - mmst2,1) * (
          + 16/3*mmst2
          )

       + lnMst2Sq*Log(mmt)*den(mmgl - mmst2,2) * (
          + 8/3*sqr(mmst2)
          )

       + lnMst2Sq*Log(mmt) * (
          - 2/3
          )

       + lnMst2Sq*Log(mmu)*sn2t * (
          - 64/9/mt*mgl
          )

       + lnMst2Sq*Log(mmu)*sqr(cs2t) * (
          - 64/9
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,1)*sn2t * (
          + 16/9*mmst1/mt*mgl
          - 16*mmst2/mt*mgl
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 16*mmst1
          + 16/9*mmst2
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 128/9*sqr(mmst1)
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,1) * (
          + 16/9*mmst1
          - 332/9*mmst2
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,2)*sn2t * (
          + 16/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          + 32/3*sqr(mmst2)
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,2) * (
          + 32/9*mmst1*mmst2
          - 178/9*sqr(mmst2)
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,3)*sqr(cs2t) * (
          - 16/9*mmst1*sqr(mmst2)
          + 16/9*pow(mmst2,3)
          )

       + lnMst2Sq*Log(mmu)*den(mmgl - mmst2,3) * (
          + 16/9*mmst1*sqr(mmst2)
          - 8/9*pow(mmst2,3)
          )

       + lnMst2Sq*Log(mmu)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 128/9*mmst1
          )

       + lnMst2Sq*Log(mmu) * (
          - 128/9
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*sn2t * (
          + 4/9*mmst1/mt*mgl
          + 52/9*mmst2/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst2/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 37/3*mmst1
          + 19/3*mmst2
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 26/3*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 14*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,3)
          )

       + lnMst2Sq*den(mmgl - mmst1,1) * (
          - 137/9*mmst1
          - 23/3*mmst2
          )

       + lnMst2Sq*den(mmgl - mmst1,2)*sn2t * (
          - 4*sqr(mmst1)/mt*mgl
          + 4*sqr(mmst2)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 94/9*mmst1*mmst2
          + 26/3*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst1,2) * (
          - 82/9*mmst1*mmst2
          - 5/3*sqr(mmst1)
          - 3*sqr(mmst2)
          )

       + lnMst2Sq*den(mmgl - mmst1,3) * (
          - 4*mmst1*sqr(mmst2)
          + 4*pow(mmst1,3)
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*sn2t * (
          + 8/3*mmst1/mt*mgl
          - 520/9*mmst2/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          - 16/9*mmst2/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 734/9*mmst1
          - 152/3*mmst2
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t * (
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 718/9*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 14*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,3)
          )

       + lnMst2Sq*den(mmgl - mmst2,1) * (
          - 2*mmsb1
          - 2*mmsb2
          + 50/3*mmst1
          + 52/9*mmst2
          - 16*mmsusy
          )

       + lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          - 8/3*mmsb1*mmst2/mt*mgl
          + 4*sqr(mmsb1)/mt*mgl
          - 8/3*mmsb2*mmst2/mt*mgl
          + 4*sqr(mmsb2)/mt*mgl
          - 8/9*mmst1*mmst2/mt*mgl
          - 64/3*mmst2*mmsusy/mt*mgl
          - 16/9*sqr(mmst2)/mt*mgl
          + 32*sqr(mmsusy)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          + 74/3*sqr(mmst2)
          )

       + lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          - 8/9*pow(mmst1,3)
          )

       + lnMst2Sq*den(mmgl - mmst2,2) * (
          - 6*mmsb1*mmst2
          + 3*sqr(mmsb1)
          - 6*mmsb2*mmst2
          + 3*sqr(mmsb2)
          + 22/9*mmst1*mmst2
          + 8/9*sqr(mmst1)
          - 48*mmst2*mmsusy
          + 20*sqr(mmst2)
          + 24*sqr(mmsusy)
          )

       + lnMst2Sq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          - 16/9*mmst1*sqr(mmst2)
          + 16/9*pow(mmst2,3)
          )

       + lnMst2Sq*den(mmgl - mmst2,3) * (
          - 8/3*mmsb1*sqr(mmst2)
          + 4*sqr(mmsb1)*mmst2
          - 8/3*mmsb2*sqr(mmst2)
          + 4*sqr(mmsb2)*mmst2
          - 8/9*mmst1*sqr(mmst2)
          + 32*mmst2*sqr(mmsusy)
          - 64/3*sqr(mmst2)*mmsusy
          + 128/9*pow(mmst2,3)
          )

       + lnMst2Sq*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 128/3*mmgl
          + 640/9*mmst1
          )

       + lnMst2Sq * (
          + 5/9
          )

       + sqr(Log(mmsusy))*den(mmgl - mmst1,1)*sn2t * (
          + 16/3*mmsusy/mt*mgl
          )

       + sqr(Log(mmsusy))*den(mmgl - mmst1,1) * (
          - 4*mmsusy
          )

       + sqr(Log(mmsusy))*den(mmgl - mmst1,2) * (
          - 16/3*mmst1*mmsusy
          )

       + sqr(Log(mmsusy))*den(mmgl - mmst2,1)*sn2t * (
          - 16/3*mmsusy/mt*mgl
          )

       + sqr(Log(mmsusy))*den(mmgl - mmst2,1) * (
          - 4*mmsusy
          )

       + sqr(Log(mmsusy))*den(mmgl - mmst2,2) * (
          - 16/3*mmst2*mmsusy
          )

       + sqr(Log(mmsusy)) * (
          + 8/3
          )

       + Log(mmsusy)*Log(mmt) * (
          - 16/3
          )

       + Log(mmsusy)*den(mmgl - mmst1,1)*sn2t * (
          - 128/3*mmsusy/mt*mgl
          )

       + Log(mmsusy)*den(mmgl - mmst1,1) * (
          + 32*mmsusy
          )

       + Log(mmsusy)*den(mmgl - mmst1,2) * (
          + 128/3*mmst1*mmsusy
          )

       + Log(mmsusy)*den(mmgl - mmst2,1)*sn2t * (
          + 128/3*mmsusy/mt*mgl
          )

       + Log(mmsusy)*den(mmgl - mmst2,1) * (
          + 32*mmsusy
          )

       + Log(mmsusy)*den(mmgl - mmst2,2) * (
          + 128/3*mmst2*mmsusy
          )

       + Log(mmsusy) * (
          + 152/9
          )

       + Log(mmt)*Log(mmu) * (
          + 64/3
          )

       + Log(mmt)*den(mmgl - mmst1,1) * (
          + 8/3*mmst1
          )

       + Log(mmt)*den(mmgl - mmst2,1) * (
          + 8/3*mmst2
          )

       + Log(mmt) * (
          + 8
          )

       + Log(mmu)*sqr(cs2t) * (
          + 128/9
          )

       + sqr(Log(mmu)) * (
          - 130/9
          )

       + Log(mmu)*den(mmgl - mmst1,1)*sn2t * (
          + 8/9*mmst1/mt*mgl
          - 16/9*mmst2/mt*mgl
          )

       + Log(mmu)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + Log(mmu)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 88/9*mmst1
          - 8/3*mmst2
          )

       + Log(mmu)*den(mmgl - mmst1,1) * (
          - 58/3*mmst1
          + 8/3*mmst2
          )

       + Log(mmu)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 16/9*mmst1*mmst2
          + 16/9*sqr(mmst1)
          )

       + Log(mmu)*den(mmgl - mmst1,2) * (
          + 16/9*mmst1*mmst2
          - 8/9*sqr(mmst1)
          )

       + Log(mmu)*den(mmgl - mmst2,1)*sn2t * (
          + 16/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + Log(mmu)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + Log(mmu)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 8/3*mmst1
          + 88/9*mmst2
          )

       + Log(mmu)*den(mmgl - mmst2,1) * (
          + 8/3*mmst1
          - 58/3*mmst2
          )

       + Log(mmu)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 16/9*mmst1*mmst2
          + 16/9*sqr(mmst2)
          )

       + Log(mmu)*den(mmgl - mmst2,2) * (
          + 16/9*mmst1*mmst2
          - 8/9*sqr(mmst2)
          )

       + Log(mmu) * (
          - 932/9
          )

       + den(mmgl - mmst1,1)*sn2t * (
          + 4/3*mmsb1/mt*mgl*zt2
          + 8*mmsb1/mt*mgl
          + 4/3*mmsb2/mt*mgl*zt2
          + 8*mmsb2/mt*mgl
          + 88/9*mmst1/mt*mgl*zt2
          + 676/9*mmst1/mt*mgl
          + 4/3*mmst2/mt*mgl*zt2
          + 56/9*mmst2/mt*mgl
          + 32/3*mmsusy/mt*mgl*zt2
          + 64*mmsusy/mt*mgl
          )

       + den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 41/9*mmst1*zt2
          + 439/9*mmst1
          - 8/3*mmst2
          )

       + den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 20*sqr(mmst1)*zt2
          - 140*sqr(mmst1)
          )

       + den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 20/3*sqr(mmst1)*zt2
          + 428/9*sqr(mmst1)
          )

       + den(mmgl - mmst1,1) * (
          - mmsb1*zt2
          - 6*mmsb1
          - mmsb2*zt2
          - 6*mmsb2
          + 50/9*mmst1*zt2
          + 61/9*mmst1
          - mmst2*zt2
          - 10/3*mmst2
          - 8*mmsusy*zt2
          - 48*mmsusy
          )

       + den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 16/9*mmst1*mmst2
          + 16/9*sqr(mmst1)
          )

       + den(mmgl - mmst1,2) * (
          - 4/3*mmsb1*mmst1*zt2
          - 8*mmsb1*mmst1
          - 4/3*mmsb2*mmst1*zt2
          - 8*mmsb2*mmst1
          - 4/3*mmst1*mmst2*zt2
          - 56/9*mmst1*mmst2
          - 32/3*mmst1*mmsusy*zt2
          - 64*mmst1*mmsusy
          + 44/3*sqr(mmst1)*zt2
          + 868/9*sqr(mmst1)
          )

       + den(mmgl - mmst1,3) * (
          + 16/3*pow(mmst1,3)*zt2
          + 112/3*pow(mmst1,3)
          )

       + den(mmgl - mmst2,1)*sn2t * (
          - 4/3*mmsb1/mt*mgl*zt2
          - 8*mmsb1/mt*mgl
          - 4/3*mmsb2/mt*mgl*zt2
          - 8*mmsb2/mt*mgl
          - 4/3*mmst1/mt*mgl*zt2
          - 56/9*mmst1/mt*mgl
          - 88/9*mmst2/mt*mgl*zt2
          - 676/9*mmst2/mt*mgl
          - 32/3*mmsusy/mt*mgl*zt2
          - 64*mmsusy/mt*mgl
          )

       + den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 20*mmst1*zt2
          - 428/3*mmst1
          - 139/9*mmst2*zt2
          - 821/9*mmst2
          )

       + den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 20*sqr(mmst1)*zt2
          + 140*sqr(mmst1)
          )

       + den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 20/3*sqr(mmst1)*zt2
          - 428/9*sqr(mmst1)
          )

       + den(mmgl - mmst2,1) * (
          - mmsb1*zt2
          - 6*mmsb1
          - mmsb2*zt2
          - 6*mmsb2
          + 17/3*mmst1*zt2
          + 398/9*mmst1
          + 110/9*mmst2*zt2
          + 163/3*mmst2
          - 8*mmsusy*zt2
          - 48*mmsusy
          )

       + den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 16/9*mmst1*mmst2
          + 16/9*sqr(mmst2)
          )

       + den(mmgl - mmst2,2) * (
          - 4/3*mmsb1*mmst2*zt2
          - 8*mmsb1*mmst2
          - 4/3*mmsb2*mmst2*zt2
          - 8*mmsb2*mmst2
          - 4/3*mmst1*mmst2*zt2
          - 56/9*mmst1*mmst2
          - 32/3*mmst2*mmsusy*zt2
          - 64*mmst2*mmsusy
          + 44/3*sqr(mmst2)*zt2
          + 868/9*sqr(mmst2)
          )

       + den(mmgl - mmst2,3) * (
          + 16/3*pow(mmst2,3)*zt2
          + 112/3*pow(mmst2,3)
          )

       - 77/2
          + 8/9*zt2
         ;  

  return resmt;
}
