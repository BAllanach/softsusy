/* Default Standard Model input parameters, from the 2012 Review of
   Particle Properties. All units approriate powers of GeV. */

#ifndef _SUMO_SMDATA_H_
#define _SUMO_SMDATA_H_

/* Inverse EM coupling constant (MSBAR, 5 flavors, scale M_Z) */
#define ALPHA_MZ_INV_DEFAULT 127.944

/* Inverse fine-structure constant */
#define ALPHA_0_INV_DEFAULT 137.035999074

/* Fermi constant */
#define G_FERMI_DEFAULT     1.166378e-5

/* SU(3) coupling constant (MSBAR, 5 flavors, Q=M_Z).
   World average given in S. Bethke, arXiv:0908.1135. */
#define ALPHA_S_DEFAULT     0.1184

/* Z boson pole mass */
#define M_Z_DEFAULT         91.1876

/* W boson pole mass */
#define M_W_DEFAULT         80.385

/* Bottom quark MSbar mass (scale M_B) */
#define MBMB_DEFAULT       4.18

/* Bottom quark pole mass, from PDGlive */
#define M_BOT_DEFAULT       4.78

/* Top quark pole mass */      
#define M_TOP_DEFAULT       173.5

/* Tau lepton pole mass */
#define M_TAU_DEFAULT       1.77682

#endif /* sumo_smdata.h */
