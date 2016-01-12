/* 
   Computes all MSSM beta functions at 1-loop or 2-loop order.

   All formulas from hep-ph/9311340 (as implemented in a Mathematica
   readable file by SPM).

   Three loop is available in some Jack and Jones paper, maybe to be
   implemented later.

   Also, betaLambda is not yet implemented. (Not really necessary.)
*/

#include "supermodel.h"

#define frac1o3 0.33333333333333333333333L
#define frac2o3 0.66666666666666666666667L
#define frac4o3 1.33333333333333333333333L
#define frac5o3 1.66666666666666666666667L
#define frac8o3 2.66666666666666666666667L 
#define frac10o3 3.33333333333333333333333L 
#define frac11o3 3.66666666666666666666667L
#define frac14o3 4.66666666666666666666667L
#define frac16o3 5.33333333333333333333333L
#define frac22o3 7.33333333333333333333333L 
#define frac26o3 8.66666666666666666666667L
#define frac28o3 9.33333333333333333333333L
#define frac32o3 10.6666666666666666666667L
#define frac52o3 17.3333333333333333333333L
#define frac88o3 29.3333333333333333333333L
#define frac128o3 42.6666666666666666666667L 
#define frac176o3 58.6666666666666666666667L 
#define frac1o9 0.111111111111111111111111L
#define frac2o9 0.222222222222222222222222L
#define frac4o9 0.444444444444444444444444L
#define frac7o9 0.777777777777777777777778L
#define frac8o9 0.888888888888888888888889L
#define frac13o9 1.44444444444444444444444L  
#define frac14o9 1.55555555555555555555556L
#define frac16o9 1.77777777777777777777778L
#define frac26o9 2.88888888888888888888889L
#define frac32o9 3.55555555555555555555556L
#define frac64o9 7.11111111111111111111111L 
#define frac199o9 22.1111111111111111111111L
#define frac398o9 44.2222222222222222222222L
#define frac1o18 0.0555555555555555555555556L
#define frac32o27 1.18518518518518518518519L
#define frac40o27 1.48148148148148148148148L
#define frac80o27 2.96296296296296296296296L 
#define frac128o27 4.74074074074074074074074L 
#define frac136o27 5.03703703703703703703704L
#define frac199o27 7.37037037037037037037037L
#define frac272o27 10.0740740740740740740741L
#define frac512o27 18.9629629629629629629630L 
#define frac808o27 29.9259259259259259259259L
#define frac3424o27 126.814814814814814814815L 
#define frac2870o81 35.4320987654320987654321L
#define frac5486o81 67.7283950617283950617284L 
#define frac1435o162 8.85802469135802469135802L
#define frac2743o162 16.9320987654320987654321L

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Full MSSM beta functions. */

int SUMO_Betas (int loop_order)
{ 
  TSIL_REAL temp;
  TSIL_REAL Sscalars, Sprimescalars, sigsca1, sigsca2, sigsca3;
  TSIL_REAL g3sq = g3 * g3;
  TSIL_REAL gsq = g * g;
  TSIL_REAL gpsq = gp * gp;
  TSIL_REAL ytopsq = ytop * ytop;
  TSIL_REAL ybotsq = ybot * ybot;
  TSIL_REAL ytausq = ytau * ytau;
  TSIL_COMPLEX M3g3sq = (M3) * g3sq;
  TSIL_COMPLEX M2gsq = (M2) * gsq;
  TSIL_COMPLEX M1gpsq = (M1) * gpsq;
  TSIL_COMPLEX atyt = atop * ytop;
  TSIL_COMPLEX abyb = abot * ybot;
  TSIL_COMPLEX atauytau = atau * ytau;

  TSIL_REAL trm2Q, trm2L, trm2u, trm2d, trm2e;  
  TSIL_REAL atsq = TSIL_CABS(atop * atop);
  TSIL_REAL absq = TSIL_CABS(abot * abot);
  TSIL_REAL atausq = TSIL_CABS(atau * atau);
  TSIL_REAL M3sqg3sq;
  TSIL_REAL M2sqgsq;
  TSIL_REAL M1sqgpsq;
  TSIL_REAL M1M2gpsqgsq;
  TSIL_REAL M2M3gsqg3sq; 
  TSIL_REAL M1M3gpsqg3sq; 
  TSIL_REAL comboT3, comboT2, comboT1;
  TSIL_REAL comboB3, comboB2, comboB1;
  TSIL_REAL comboTAU2, comboTAU1;
  TSIL_REAL comboAAbt, comboAAbtau;
  TSIL_REAL termT, termB, termTAU;

  beta_g3 = -3.0L * SUMO_oneloopfactor;
  beta_g = SUMO_oneloopfactor;
  beta_gp = 11.0L * SUMO_oneloopfactor;
  
  beta_M3 = -6.0L * SUMO_oneloopfactor * M3g3sq;
  beta_M2 = 2.0L * SUMO_oneloopfactor * M2gsq;
  beta_M1 = 22.0L * SUMO_oneloopfactor * M1gpsq;

  beta_ytop = SUMO_oneloopfactor * (
               6.0L * ytopsq + ybotsq - frac16o3 * g3sq 
               - 3.0L * gsq - frac13o9 * gpsq);

  beta_atop = (atop) * (beta_ytop) +
               SUMO_oneloopfactor * (ytop) * 
               (frac26o9 * M1gpsq + 6.0L * M2gsq + frac32o3 * M3g3sq 
                  + 2.0L * abyb + 12.0L * atyt);

  beta_ybot = SUMO_oneloopfactor * (
               6.0L * ybotsq + ytopsq + ytausq - frac16o3 * g3sq 
               - 3.0L * gsq - frac7o9 * gpsq);

  beta_abot = (abot) * (beta_ybot) +
               SUMO_oneloopfactor * (ybot) * 
               (frac14o9 * M1gpsq + 6.0L * M2gsq + frac32o3*M3g3sq 
                  + 2.0L * (atyt + atauytau) + 12.0L * abyb);

  beta_ytau = SUMO_oneloopfactor * (
               4.0L * ytausq + 3.0L * (ybotsq - gsq - gpsq)); 

  beta_atau = (atau) * (beta_ytau) +
               SUMO_oneloopfactor * (ytau) * 
               (6.0L *  (M1gpsq + M2gsq + abyb) + 8.0L * atauytau);
  
  beta_mu = SUMO_oneloopfactor * (
             3.0L * (ytopsq + ybotsq  - gsq) - gpsq + ytausq);

  beta_b = (b) * (beta_mu) +
            SUMO_oneloopfactor * (mu) * 
            (6.0L * (abyb + atyt + M2gsq) + 2.0L * (atauytau + M1gpsq));

  temp = 0.75L * gsq + 0.25L * gpsq;

  beta_vu = SUMO_oneloopfactor * (temp - 3.0L * ytopsq);

  beta_vd = SUMO_oneloopfactor * (temp - 3.0L * ybotsq - ytausq);

  M3sqg3sq = TSIL_CABS(M3g3sq * M3);
  M2sqgsq = TSIL_CABS(M2gsq * M2);
  M1sqgpsq = TSIL_CABS(M1gpsq * M1);

  trm2Q = m2Q[0] + m2Q[1] + m2Q[2];
  trm2L = m2L[0] + m2L[1] + m2L[2];
  trm2u = m2u[0] + m2u[1] + m2u[2];
  trm2d = m2d[0] + m2d[1] + m2d[2];
  trm2e = m2e[0] + m2e[1] + m2e[2];

  Sscalars = gpsq * (m2Hu - m2Hd + trm2e - trm2L
                     + trm2Q + trm2d - 2.0L * trm2u);

  beta_m2Hu = SUMO_oneloopfactor * (
               -6.0L * M2sqgsq - 2.0L * M1sqgpsq + Sscalars);

  beta_m2Hd = beta_m2L12 = beta_m2L3 = 
         beta_m2Hu - 2.0L * SUMO_oneloopfactor * Sscalars;

  beta_m2Q12 = beta_m2Q3 = SUMO_oneloopfactor * (
                -frac32o3 * M3sqg3sq - 6.0L * M2sqgsq 
                - frac2o9 * M1sqgpsq + frac1o3 * Sscalars);

  beta_m2u12 = beta_m2u3 = SUMO_oneloopfactor * (
                -frac32o3 * M3sqg3sq - frac32o9 * M1sqgpsq 
                - frac4o3 * Sscalars);

  beta_m2d12 = beta_m2d3 = SUMO_oneloopfactor * (
                -frac32o3 * M3sqg3sq - frac8o9 * M1sqgpsq 
                + frac2o3 * Sscalars);

  beta_m2e12 = beta_m2e3 = SUMO_oneloopfactor * (
                -8.0L * M1sqgpsq + 2.0L * Sscalars);

  termT = ytopsq * (m2Q[2] + m2u[2] + m2Hu) + atsq;
  termB = ybotsq * (m2Q[2] + m2d[2] + m2Hd) + absq;
  termTAU = ytausq * (m2L[2] + m2e[2] + m2Hd) + atausq;

  temp = SUMO_oneloopfactor * 2.0L * termT;

  beta_m2Hu += 3.0L * temp;
  beta_m2u3 += 2.0L * temp;
  beta_m2Q3 += temp;

  temp = SUMO_oneloopfactor * 2.0L * termB;

  beta_m2Hd += 3.0L * temp;
  beta_m2d3 += 2.0L * temp;
  beta_m2Q3 += temp;

  temp = SUMO_oneloopfactor * 2.0L * termTAU;

  beta_m2Hd += temp;
  beta_m2e3 += 2.0L * temp;
  beta_m2L3 += temp;
                
  if (2 == loop_order) {
    temp = SUMO_twoloopfactor * (
                   14.0L * g3sq + 9.0L * gsq + frac11o3 * gpsq
                   -4.0L * (ybotsq + ytopsq));

    beta_g3 += temp;

    beta_M3 += 2.0L * temp * M3g3sq 
                    + SUMO_twoloopfactor * g3sq * 
                  (28.0L * M3g3sq + 18.0L * M2gsq + frac22o3 * M1gpsq  
                   + 8.0L * (abyb + atyt));

    temp = SUMO_twoloopfactor * (
                  25.0L * gsq + 24.0L * g3sq + 3.0L * gpsq
                  -6.0L * (ybotsq + ytopsq) - 2.0L * ytausq);

    beta_g += temp;

    beta_M2 += 2.0L * temp * M2gsq 
                    + SUMO_twoloopfactor * gsq * (
                  50.L * M2gsq + 6.0L * M1gpsq  + 48.0L * M3g3sq 
                  + 12.0L * (abyb + atyt) + 4.0L * atauytau);

    temp = SUMO_twoloopfactor * (
                  frac199o9 * gpsq + 9.0L * gsq + frac88o3 * g3sq 
                  - frac14o3 * ybotsq - frac26o3 * ytopsq - 6.0L * ytausq);  

    beta_gp += temp;

    beta_M1 += 2.0L * temp * M1gpsq
                + SUMO_twoloopfactor * gpsq * (
                  frac398o9 * M1gpsq + 18.0L * M2gsq + frac176o3 * M3g3sq
                  + frac28o3 * abyb + 12.0L * atauytau + frac52o3 * atyt); 

    temp = SUMO_twoloopfactor * ( 
                  g3sq * (-frac16o9 * g3sq + frac136o27 * gpsq 
                  + 8.0L * gsq) + 7.5L * gsq * gsq 
                  + gpsq * (frac2743o162 * gpsq + frac5o3 * gsq ) 
                  + ytopsq * (16.0L * g3sq + 2.0L * gpsq + 6.0L * gsq 
                  - 5.0L * ybotsq - 22.0L * ytopsq) 
                  + ybotsq * (frac2o3 * gpsq - 5.0L * ybotsq - ytausq)); 

    beta_ytop += temp;

    beta_atop += (atop) * temp 
                  + SUMO_twoloopfactor * (ytop) * (
                  + atyt * (12.0L * gsq + 32.0L * g3sq + 4.0L * gpsq 
                      - 10.0L * ybotsq - 88.0L * ytopsq)
                  + abyb * (frac4o3 * gpsq - 2.0L * ytausq
                      - 10.0L * ytopsq - 20.L * ybotsq)
                  - 2.0L * atauytau * ybotsq
                  + M3g3sq * (frac64o9 * g3sq - frac272o27 * gpsq 
                      - 16.0L * gsq - 32.0L * ytopsq)
                  + M2gsq * (-16.0L * g3sq - frac10o3 * gpsq 
                      - 30.0L * gsq - 12.0L * ytopsq)
                  + M1gpsq * (-frac272o27 * g3sq - frac5486o81 * gpsq 
                     - frac10o3 * gsq - 4.0L * ytopsq -frac4o3 * ybotsq));
     
    temp = SUMO_twoloopfactor * (
                    ybotsq * (16.0L * g3sq + frac2o3 * gpsq  
                    + 6.0L * gsq - 22.0L * ybotsq - 3.0L * ytausq) 
                    + ytopsq * (frac4o3 * gpsq - 5.0L * (ybotsq + ytopsq))
                    + ytausq * (2.0L * gpsq - 3.0L * ytausq) 
                    + g3sq * (-frac16o9 * g3sq + frac40o27 * gpsq
                    + 8.0L * gsq) + 7.5L * gsq * gsq
                    + gpsq * (frac1435o162 * gpsq + frac5o3 * gsq)); 

    beta_ybot += temp;

    beta_abot += (abot) * temp 
                  + SUMO_twoloopfactor * (ybot) * (
                    + atyt * (frac8o3 * gpsq - 10.0L * ybotsq 
                         - 20.0L * ytopsq)
                    + abyb * (12.0L * gsq + 32.0L * g3sq + frac4o3 * gpsq 
                         - 88.0L * ybotsq - 10.0L * ytopsq - 6.0L * ytausq)
                    + atauytau * (4.0L * gpsq  - 6.0L * ybotsq  
                         - 12.0L * ytausq)
                    + M3g3sq * (frac64o9 * g3sq - frac80o27 * gpsq 
                        - 16.0L * gsq - 32.0L * ybotsq)
                    + M2gsq * (-16.0L * g3sq - frac10o3 * gpsq 
                        - 30.0L * gsq - 12.0L * ybotsq)
                    + M1gpsq * (-frac80o27 * g3sq - frac2870o81 * gpsq 
                        - frac10o3 * gsq - frac4o3 * ybotsq 
                        - frac8o3 * ytopsq - 4.0L * ytausq));

    temp = SUMO_twoloopfactor * ( 
                    ytausq * (2.0L * gpsq + 6.0L * gsq
                       - 9.0L * ybotsq - 10.0L * ytausq) 
                  + ybotsq * (16.0L * g3sq - frac2o3 * gpsq 
                       - 3.0L * ytopsq - 9.0L * ybotsq) 
                  + gpsq * (37.5L * gpsq + 3.0L * gsq) + 7.5L * gsq * gsq);

    beta_ytau += temp;

    beta_atau += (atau) * temp 
                  + SUMO_twoloopfactor * (ytau) * (
                 - atyt * 6.0L * ybotsq
                 + atauytau * (12.0L * gsq  + 4.0L * gpsq  
                   - 18.0L * ybotsq  - 40.0L * ytausq)
                 + abyb * (32.0L * g3sq - frac4o3 * gpsq - 6.0L * ytopsq
                   -36.0L * ybotsq - 18.0L * ytausq)
                 - 32.0L * M3g3sq * ybotsq
                 - M2gsq * (6.0L * gpsq + 30.0L * gsq + 12.0L * ytausq)
                 - M1gpsq * (150.0L * gpsq + 6.0L * gsq 
                     - frac4o3 * ybotsq + 4.0L * ytausq));
                    
    temp = SUMO_twoloopfactor * (
                  ytopsq * (16.0L * g3sq + frac4o3 * gpsq 
                   - 6.0L * ybotsq - 9.0L *ytopsq)
                + ybotsq * (16.0L * g3sq - frac2o3 * gpsq 
                   - 9.0L * ybotsq) 
                + ytausq * (2.0L * gpsq - 3.0L * ytausq) 
                + gpsq * (11.5L *gpsq  + 3.0L * gsq) + 7.5L * gsq * gsq);

    beta_mu += temp;

    beta_b += (b) * temp 
                 + SUMO_twoloopfactor * (mu) * (
                 atyt * (32.0L * g3sq + frac8o3 * gpsq 
                    - 12.0L * ybotsq - 36.0L * ytopsq)
                 + abyb * (32.0L * g3sq - frac4o3 * gpsq 
                    - 36.0L * ybotsq - 12.0L * ytopsq)
                 + atauytau * (4.0L * gpsq - 12.0L * ytausq)
                 - M3g3sq * 32.0L * (ytopsq + ybotsq)
                 - M2gsq * (30.0L * gsq + 6.0L * gpsq)
                 + M1gpsq * (-6.0L * gsq - 46.0L * gpsq + frac4o3 * ybotsq 
                   - frac8o3 * ytopsq - 4.0L * ytausq));

    temp = -gsq * (3.0L * gsq + 0.75L * gpsq) - 2.875L * gpsq * gpsq;

    beta_vu += SUMO_twoloopfactor * (
                ytopsq * (-16.0L * g3sq - frac4o3 * gpsq 
                          + 9.0L * ytopsq + 3.0L * ybotsq) + temp);

    beta_vd += SUMO_twoloopfactor * (
                ybotsq * (-16.0L * g3sq + frac2o3 * gpsq 
                          + 9.0L * ybotsq + 3.0L * ytopsq)
                + ytausq * (3.0L * ytausq - 2.0L * gpsq) + temp);

    Sprimescalars =  gpsq * (
      frac8o3 * g3sq * (trm2d + trm2Q - 2.0L * trm2u) 
      + 1.5L * gsq * (trm2Q - trm2L + m2Hu - m2Hd) + 
      + gpsq * (frac2o9 * trm2d + frac1o18 * trm2Q - frac16o9 * trm2u
                + 2.0L * trm2e + 0.5L * (m2Hu - m2Hd - trm2L))    
      + ytopsq * (4.0L * m2u[2] - 3.0L * m2Hu - m2Q[2])
      + ybotsq * (3.0L * m2Hd - m2Q[2] - 2.0L * m2d[2])
      + ytausq * (m2Hd + m2L[2] - 2.0L * m2e[2]));

    sigsca1 = gpsq * gpsq * (m2Hu + m2Hd + trm2L + 2.0L * trm2e
                  + frac1o3 * trm2Q + frac2o3 * trm2d + frac8o3 * trm2u);
    sigsca2 = gsq * gsq * (trm2L + m2Hu + m2Hd + 3.0L * trm2Q);
    sigsca3 = g3sq * g3sq * (trm2d + trm2u + 2.0L * trm2Q);

    M1M2gpsqgsq = TSIL_CREAL(M1gpsq * SUMO_CONJ (M2gsq));
    M2M3gsqg3sq = TSIL_CREAL(M2gsq * SUMO_CONJ (M3g3sq)); 
    M1M3gpsqg3sq = TSIL_CREAL(M1gpsq * SUMO_CONJ (M3g3sq)); 

    temp = SUMO_twoloopfactor * (
           M3sqg3sq * (32.0L * gsq - frac128o3 * g3sq + frac32o27 * gpsq)
           + M2sqgsq * (33.0L * gsq + 32.0L * g3sq + frac2o3 * gpsq)
           + M1sqgpsq * (frac2o3*gsq + frac32o27*g3sq + frac199o27 * gpsq)
           + frac2o3 * M1M2gpsqgsq
           + frac32o27 * M1M3gpsqg3sq
           + 32.0L * M2M3gsqg3sq
           + frac1o9 * sigsca1 + 3.0L * sigsca2 + frac16o3 * sigsca3
           + frac2o3 * Sprimescalars);

    beta_m2Q12 += temp;
    beta_m2Q3 += temp;

    temp = SUMO_twoloopfactor * (
           M3sqg3sq * (-frac128o3 * g3sq + frac512o27 * gpsq)
           + M1sqgpsq * (frac512o27*g3sq + frac3424o27 * gpsq)
           + frac512o27 * M1M3gpsqg3sq
           + frac16o9 * sigsca1 + frac16o3 * sigsca3
           - frac8o3 * Sprimescalars);

    beta_m2u12 += temp;
    beta_m2u3 += temp;

    temp = SUMO_twoloopfactor * (
           M3sqg3sq * (-frac128o3 * g3sq + frac128o27 * gpsq)
           + M1sqgpsq * (frac128o27*g3sq + frac808o27 * gpsq)
           + frac128o27 * M1M3gpsqg3sq
           + frac4o9 * sigsca1 + frac16o3 * sigsca3
           + frac4o3 * Sprimescalars);

    beta_m2d12 += temp;
    beta_m2d3 += temp;

    temp = SUMO_twoloopfactor * (
           M2sqgsq * (33.0L * gsq + 6.0L * gpsq)
           + M1sqgpsq * (6.0L * gsq + 69.0L * gpsq)
           + 6.0L * M1M2gpsqgsq
           + sigsca1 + 3.0L * sigsca2 -2.0L * Sprimescalars);

    beta_m2L12 += temp;
    beta_m2L3 += temp;
    beta_m2Hd += temp;
    beta_m2Hu += temp + 4.0L * SUMO_twoloopfactor * Sprimescalars;

    temp = SUMO_twoloopfactor * (
           M1sqgpsq * 312.0L * gpsq + 4.0L * (sigsca1 + Sprimescalars));
    
    beta_m2e12 += temp;
    beta_m2e3 += temp;

    comboT3 = 64.0L * (ytopsq * M3sqg3sq 
             - TSIL_CREAL(atyt * SUMO_CONJ(M3g3sq)));

    comboT2 = 24.0L * (ytopsq * M2sqgsq 
             - TSIL_CREAL(atyt * SUMO_CONJ(M2gsq)));

    comboT1 = frac16o3 * (ytopsq * M1sqgpsq 
             - TSIL_CREAL(atyt * SUMO_CONJ(M1gpsq)));

    comboB3 = 64.0L * (ybotsq * M3sqg3sq 
             - TSIL_CREAL(abyb * SUMO_CONJ(M3g3sq)));

    comboB2 = 24.0L * (ybotsq * M2sqgsq 
             - TSIL_CREAL(abyb * SUMO_CONJ(M2gsq)));

    comboB1 = frac8o3 * (ybotsq * M1sqgpsq 
             - TSIL_CREAL(abyb * SUMO_CONJ(M1gpsq)));

    comboTAU2 = 24.0L * (ytausq * M2sqgsq
                - TSIL_CREAL(atauytau * SUMO_CONJ(M2gsq)));

    comboTAU1 = 8.0L * (ytausq * M1sqgpsq
                - TSIL_CREAL(atauytau * SUMO_CONJ(M1gpsq)));

    comboAAbt = -8.0L * TSIL_CREAL(atyt * SUMO_CONJ(abyb));

    comboAAbtau = -4.0L * TSIL_CREAL(abyb * SUMO_CONJ(atauytau));

    beta_m2Q3 += SUMO_twoloopfactor * (
       comboT1 + comboB1 + comboAAbtau
       - atsq * 20.0L * ytopsq
       + termT * (frac8o3 * gpsq - 20.0L * ytopsq) 
       + termB * (frac4o3 * gpsq - 20.0L * ybotsq - 2.0L * ytausq) 
       - (20.0L * absq + termTAU * 2.0L) * ybotsq);

    beta_m2u3 += SUMO_twoloopfactor * (
       comboT2 - 0.5L * comboT1 + comboAAbt
       + termT * (12.0L * gsq - frac4o3 * gpsq - 4.0L * ybotsq 
                  -32.0L * ytopsq)
       - (atsq * 32.0L + termB * 4.0L) * ytopsq);

    beta_m2d3 += SUMO_twoloopfactor * (
       comboB1 + comboB2 + comboAAbt + 2.0L * comboAAbtau
       + termB * (12.0L * gsq + frac4o3 * gpsq 
             - 4.0L * ytopsq - 4.0L * ytausq - 32.0L * ybotsq)  
       - (32.0L * absq + 4.0L * (termT + termTAU)) * ybotsq );
       
    beta_m2L3 += SUMO_twoloopfactor * (
       comboTAU1 + 3.0L * comboAAbtau
       - (6.0L * termB + 12.0L * atausq ) * ytausq
       + termTAU * (4.0L * gpsq -6.0L * ybotsq - 12.0L * ytausq) );
              
    beta_m2e3 += SUMO_twoloopfactor * (
       comboTAU2 - comboTAU1 + 6.0L * comboAAbtau
       - (16.0L * atausq + 12.0L * termB) * ytausq 
       + termTAU * (12.0L * gsq -4.0L * gpsq - 12.0L * ybotsq 
                   - 16.0L * ytausq) );

    beta_m2Hu += SUMO_twoloopfactor * (
      comboT1 + comboT3 + 1.5L * comboAAbt
       - (36.0L * atsq + 6.0L * termB) * ytopsq 
       + termT * (32.0L * g3sq + frac8o3 * gpsq - 6.0L * ybotsq
                  - 36.0L * ytopsq));

    beta_m2Hd += SUMO_twoloopfactor * (
       comboB3 - comboB1 + comboTAU1 + 1.5L * comboAAbt
       + termB * (32.0L * g3sq - 36.0L * ybotsq - 6.0L * ytopsq 
                  - frac4o3 * gpsq) 
       + termTAU * (4.0L * gpsq - 12.0L * ytausq)
       - (6.0L * termT +  36.0L * absq) * ybotsq
       - 12.0L * atausq * ytausq );       
  } 

  beta_g3 *= g3 * g3sq;
  beta_g *= g * gsq;
  beta_gp *= gp * gpsq;

  beta_ytop *= ytop;
  beta_ybot *= ybot;
  beta_ytau *= ytau;

  beta_mu *= mu;

  beta_vu *= vu;
  beta_vd *= vd;

  return 0;
}  
