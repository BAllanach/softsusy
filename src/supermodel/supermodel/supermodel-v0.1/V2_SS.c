#include "supermodel.h"

SUMO_COMPLEX SUMO_V2_SS (void)
{ 
  int i,j;
  TSIL_COMPLEX tempj;
  TSIL_COMPLEX temp[21];
  TSIL_REAL qq = Q2;  
  TSIL_COMPLEX V2_00 = 0.0L;
  TSIL_COMPLEX V2_pm = 0.0L;
  TSIL_COMPLEX V2_0p = 0.0L;
  TSIL_COMPLEX V2_sf0g2 = 0.0L;
  TSIL_COMPLEX V2_sfpg2 = 0.0L;
  TSIL_COMPLEX V2_sfsfg2= 0.0L;
  TSIL_COMPLEX V2_sfsfgp2= 0.0L;
  TSIL_COMPLEX V2_sf0y2 = 0.0L;
  TSIL_COMPLEX V2_sfpy2 = 0.0L;
  TSIL_COMPLEX V2_sfsfy2= 0.0L;

  TSIL_REAL h0 = m2_phi0[0];
  TSIL_REAL H0 = m2_phi0[1];
  TSIL_REAL G0 = m2_phi0[2];
  TSIL_REAL A0 = m2_phi0[3];
  TSIL_REAL Gp = m2_phip[0];
  TSIL_REAL Hp = m2_phip[1];

  TSIL_REAL m2_sf[21] = {
    m2_snu[0], m2_snu[1], m2_snu[2],
    m2_seL[0], m2_seL[1],
    m2_suL[0], m2_suL[1],
    m2_sdL[0], m2_sdL[1],
    m2_stau[0], m2_stau[1],
    m2_stop[0], m2_stop[1],
    m2_sbot[0], m2_sbot[1],
    m2_suR[0], m2_suR[1],
    m2_sdR[0], m2_sdR[1],
    m2_seR[0], m2_seR[1]};

  TSIL_REAL x_sf[21] = {
    0.5L, 0.5L, 0.5L, -0.5L, -0.5L, 0.5L, 0.5L, -0.5L, -0.5L,
    -0.5L * Lstau[0] * Lstauc[0], 
    -0.5L * Lstau[1] * Lstauc[1],
     0.5L * Lstop[0] * Lstopc[0],  
     0.5L * Lstop[1] * Lstopc[1],
    -0.5L * Lsbot[0] * Lsbotc[0], 
    -0.5L * Lsbot[1] * Lsbotc[1],
    0.L, 0.L, 0.L, 0.L, 0.L, 0.L};

  TSIL_REAL xp_sf[21] = {
   -0.5L, -0.5L, -0.5L, -0.5L, -0.5L, 
    1.0L/6.0L, 1.0L/6.0L, 1.0L/6.0L, 1.0L/6.0L, 
    -0.5L * Lstau[0] * Lstauc[0] + Rstau[0] * Rstauc[0],
    -0.5L * Lstau[1] * Lstauc[1] + Rstau[1] * Rstauc[1],
(Lstop[0] * Lstopc[0] - 4.0L * Rstop[0] * Rstopc[0])/6.0L,
(Lstop[1] * Lstopc[1] - 4.0L * Rstop[1] * Rstopc[1])/6.0L,
(Lsbot[0] * Lsbotc[0] + 2.0L * Rsbot[0] * Rsbotc[0])/6.0L,
(Lsbot[1] * Lsbotc[1] + 2.0L * Rsbot[1] * Rsbotc[1])/6.0L,
    -2.0L/3.0L, -2.0L/3.0L, 1.0L/3.0L, 1.0L/3.0L, 1.0L, 1.0L};

  TSIL_REAL n_sf[21] = {1.L, 1.L, 1.L, 1.L, 1.L, 
                        3.L, 3.L, 3.L, 3.L, 
                        1.L, 1.L, 3.L, 3.L, 3.L, 3.L,  
                        3.L, 3.L, 3.L, 3.L, 1.L, 1.L};

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.15) */

  V2_00 += 3.0L * c2alpha * c2alpha *
             (SUMO_Fvac_SS (h0, h0, qq)
            + SUMO_Fvac_SS (H0, H0, qq));

  V2_00 += 3.0L * c2beta0 * c2beta0 *
             (SUMO_Fvac_SS (G0, G0, qq)
            + SUMO_Fvac_SS (A0, A0, qq));

  V2_00 += (4.0L - 6.0L * c2alpha * c2alpha) *
           SUMO_Fvac_SS (h0, H0, qq);

  V2_00 += (4.0L - 6.0L * c2beta0 * c2beta0) *
           SUMO_Fvac_SS (G0, A0, qq);

  V2_00 += 2.0L * c2alpha * c2beta0 * (
             SUMO_Fvac_SS (h0, A0, qq)
           + SUMO_Fvac_SS (H0, G0, qq)
           - SUMO_Fvac_SS (h0, G0, qq)
           - SUMO_Fvac_SS (H0, A0, qq));

  V2_00 *= (g2plusgp2)/32.0L;

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.16) */

  V2_pm += c2betapm * c2betapm *
            (SUMO_Fvac_SS (Gp, Gp, qq)
           + SUMO_Fvac_SS (Hp, Hp, qq));

  V2_pm += (1.0L - 2.0L * c2betapm * c2betapm) *
           SUMO_Fvac_SS (Gp, Hp, qq);

  V2_pm *= (g2plusgp2)/4.0L;

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.17) */

  V2_0p +=  (g2 * (1.0L - s2alpha * s2betapm)
           + gp2 * c2alpha * c2betapm) *
     (SUMO_Fvac_SS (h0, Hp, qq) + SUMO_Fvac_SS (H0, Gp, qq));

  V2_0p +=  (g2 * (1.0L + s2alpha * s2betapm)
           - gp2 * c2alpha * c2betapm) *
     (SUMO_Fvac_SS (h0, Gp, qq) + SUMO_Fvac_SS (H0, Hp, qq));

  V2_0p +=  (g2 * (1.0L - s2beta0 * s2betapm)
           + gp2 * c2beta0 * c2betapm) *
     (SUMO_Fvac_SS (A0, Hp, qq) + SUMO_Fvac_SS (G0, Gp, qq));

  V2_0p += (g2 * (1.0L + s2beta0 * s2betapm)
           - gp2 * c2beta0 * c2betapm) *
     (SUMO_Fvac_SS (A0, Gp, qq) + SUMO_Fvac_SS (G0, Hp, qq));

  V2_0p *= 0.125L;

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.18) */

  for (j=0; j<4; j++) {
    tempj = (kd[j] * kdc[j] - ku[j] * kuc[j]);
    for (i=0; i<21; i++) {
      V2_sf0g2 += n_sf[i] * (x_sf[i] * g2 - xp_sf[i] * gp2) *
                  tempj *
                  SUMO_Fvac_SS (m2_sf[i], m2_phi0[j], qq);
  }}

  V2_sf0g2 *= 0.25L;

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.19) */

  for (j=0; j<2; j++) {
    tempj = (kup[j] * kup[j] - kdp[j] * kdp[j]);
    for (i=0; i<21; i++) {      
      V2_sfpg2 += n_sf[i] * (x_sf[i] * g2 + xp_sf[i] * gp2) *
               tempj *
               SUMO_Fvac_SS (m2_sf[i], m2_phip[j], qq);
  }}

  V2_sfpg2 *= 0.5L;

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.20) */

  /* In following, make use of F_SS (x,y) = A(x) A(y) */

  for (i=0; i<15; i++) {
    temp[i] = n_sf[i] * x_sf[i] * TSIL_A (m2_sf[i], qq);
  }

  for (i=1; i<15; i++) {
    for (j=0; j<i; j++) {
      V2_sfsfg2 += temp[i] * temp[j];
  }}

  for (i=0; i<15; i++) {
    V2_sfsfg2 += (0.5L + 0.5L/n_sf[i]) * temp[i] * temp[i];
  }

  V2_sfsfg2 += 0.75L * Lstop[0] * Lstopc[0]
               * Lstop[1] * Lstopc[1] *
               SUMO_Fvac_SS (m2_stop[0], m2_stop[1], qq);

  V2_sfsfg2 += 0.75L * Lsbot[0] * Lsbotc[0]
               * Lsbot[1] * Lsbotc[1] *
               SUMO_Fvac_SS (m2_sbot[0], m2_sbot[1], qq);

  V2_sfsfg2 += 0.25L * Lstau[0] * Lstauc[0]
               * Lstau[1] * Lstauc[1] *
               SUMO_Fvac_SS (m2_stau[0], m2_stau[1], qq);

  for (i=0; i<2; i++) {
    V2_sfsfg2 += 
      1.5L * SUMO_Fvac_SS (m2_suL[i], m2_sdL[i], qq);

    V2_sfsfg2 += 
      0.5L * SUMO_Fvac_SS (m2_seL[i], m2_snu[i], qq);

    V2_sfsfg2 += 0.5L * Lstau[i] * Lstauc[i] *
      SUMO_Fvac_SS (m2_stau[i], m2_snu[2], qq);

    for (j=0; j<2; j++) {
      V2_sfsfg2 += 1.5L * Lstop[i] * Lstopc[i]
               * Lsbot[j] * Lsbotc[j] *
               SUMO_Fvac_SS (m2_stop[i], m2_sbot[j], qq);
    }
  }

  V2_sfsfg2 *= g2;

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.21) */

  for (i=0; i<21; i++) {
    temp[i] = n_sf[i] * xp_sf[i] * TSIL_A (m2_sf[i], qq);
  }

  for (i=1; i<21; i++) {
    for (j=0; j<i; j++) {
      V2_sfsfgp2 += temp[i] * temp[j];
  }}

  for (i=0; i<21; i++) {
    V2_sfsfgp2 += (0.5L + 0.5L/n_sf[i]) * temp[i] * temp[i];
  }

  V2_sfsfgp2 += (1.0L/12.0L) * 
    (Lstop[0] * Lstopc[1] - 4.0L * Rstop[0] * Rstopc[1]) *
    (Lstopc[0] * Lstop[1] - 4.0L * Rstopc[0] * Rstop[1]) *
               SUMO_Fvac_SS (m2_stop[0], m2_stop[1], qq);

  V2_sfsfgp2 += (1.0L/12.0L) * 
    (Lsbot[0] * Lsbotc[1] +2.0L * Rsbot[0] * Rsbotc[1]) *
    (Lsbotc[0] * Lsbot[1] +2.0L * Rsbotc[0] * Rsbot[1]) *
               SUMO_Fvac_SS (m2_sbot[0], m2_sbot[1], qq);

  V2_sfsfgp2 += 0.25L * 
    (Lstau[0] * Lstauc[1] -2.0L * Rstau[0] * Rstauc[1]) *
    (Lstauc[0] * Lstau[1] -2.0L * Rstauc[0] * Rstau[1]) *
               SUMO_Fvac_SS (m2_stau[0], m2_stau[1], qq);

  V2_sfsfgp2 *= gp2;

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.22) */

  for (i=0; i<2; i++) {
    for (j=0; j<4; j++) {
      V2_sf0y2 += 
        1.5L * ytop2 * ku[j] * kuc[j] *
        SUMO_Fvac_SS (m2_stop[i], m2_phi0[j], qq);

      V2_sf0y2 +=
        1.5L * ybot2 * kd[j] * kdc[j] *
        SUMO_Fvac_SS (m2_sbot[i], m2_phi0[j], qq);

      V2_sf0y2 +=
        0.5L * ytau2 * kd[j] * kdc[j] *
        SUMO_Fvac_SS (m2_stau[i], m2_phi0[j], qq);
  }}
  
  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.23) */

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      V2_sfpy2 += 
        3.0L * (ytop2 * kup[j] * kup[j] *
        Rstop[i] * Rstopc[i] +
        ybot2 * kdp[j] * kdp[j] *
        Lstop[i] * Lstopc[i]) *
        SUMO_Fvac_SS (m2_stop[i], m2_phip[j], qq);

      V2_sfpy2 +=
        3.0L * (ytop2 * kup[j] * kup[j] *
        Lsbot[i] * Lsbotc[i] +
        ybot2 * kdp[j] * kdp[j] *
        Rsbot[i] * Rsbotc[i]) *
        SUMO_Fvac_SS (m2_sbot[i], m2_phip[j], qq);

      V2_sfpy2 += ytau2 * kdp[j] * kdp[j] *
        Rstau[i] * Rstauc[i] *
        SUMO_Fvac_SS (m2_stau[i], m2_phip[j], qq);
  }}

  for (j=0; j<2; j++) {
      V2_sfpy2 += ytau2 * kdp[j] * kdp[j] *
        SUMO_Fvac_SS (m2_snu[2], m2_phip[j], qq);
  }

  /* ------------------------------------------------------------------- */
  /* hep-ph/0206136 eq. (3.24) */

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {

      V2_sfsfy2 += 3.0L * (ytop2 *
    (Lstop[i] * Lstopc[i] * Rstop[j] * Rstopc[j]
    + 3.0L * Lstop[i] * Rstopc[i] * Rstop[j] * Lstopc[j]))
    * SUMO_Fvac_SS (m2_stop[i], m2_stop[j], qq);

      V2_sfsfy2 += 3.0L * (ybot2 *
    (Lsbot[i] * Lsbotc[i] * Rsbot[j] * Rsbotc[j]
    + 3.0L * Lsbot[i] * Rsbotc[i] * Rsbot[j] * Lsbotc[j]))
    * SUMO_Fvac_SS (m2_sbot[i], m2_sbot[j], qq);

      V2_sfsfy2 += 3.0L * (ytop2 *
    Rstop[i] * Rstopc[i] * Lsbot[j] * Lsbotc[j]
                   + ybot2 *
    Lstop[i] * Lstopc[i] * Rsbot[j] * Rsbotc[j]) *
    SUMO_Fvac_SS (m2_stop[i], m2_sbot[j], qq);

      V2_sfsfy2 += 6.0L * (ytau) * (ybot) * TSIL_CREAL(
    Lsbot[i] * Rsbotc[i] * Rstau[j] * Lstauc[j]) *
    SUMO_Fvac_SS (m2_sbot[i], m2_stau[j], qq);

      V2_sfsfy2 += 0.5L * ytau2
    * (Lstau[i] * Rstau[j] + Rstau[i] * Lstau[j])
    * (Lstauc[i] * Rstauc[j] + Rstauc[i] * Lstauc[j])
    * SUMO_Fvac_SS (m2_stau[i], m2_stau[j], qq);

    }}

  for (i=0; i<2; i++) {
    V2_sfsfy2 += ytau2 * Rstau[i] * Rstauc[i] *
                SUMO_Fvac_SS (m2_snu[2], m2_stau[i], qq);
  }

  /* ------------------------------------------------------------------- */

  return (SUMO_twoloopfactor * (V2_00 + V2_pm + V2_0p
				     + V2_sf0g2
				     + V2_sfpg2
				     + V2_sfsfg2
				     + V2_sfsfgp2
				     + V2_sf0y2
				     + V2_sfpy2
				     + V2_sfsfy2
				      ));
}
