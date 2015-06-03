/* Neutral Higgs self energy and pole masses */

#include "supermodel.h"
#include "self_scalar.h"

int SUMO_phi0_pole_masses (int which, int approx_oneloop, int approx_twoloop);

/* ___________________________________________________________________ */
/* The approximation options for neutral Higgs scalars are below.
   These can be mixed and matched, with a few limitations, to wit:
   for H0, G0, and A0, all options involving the effective potential 
   approximations are invalid. Also, approx_oneloop=0 and 
   approx_oneloop=4 are each invalid unless approx_twoloop=0.

One-loop approximation options:
  approx_oneloop = 0 is no 1-loop contribution
  approx_oneloop = 1 is 1-loop effective potential approximation
  approx_oneloop = 2 is full 1-loop self-energy
  approx_oneloop = 3 is full 1-loop self-energy with tree masses
                     expanded around pole masses
  approx_oneloop = 4 is full 1-loop self-energy with tree masses 
                     replaced by pole masses

Two-loop approximation options:
  approx_twoloop = 0 is no 2-loop contribution
  approx_twoloop = 1 is 2-loop effective potential approximation
  approx_twoloop = 2 is 2-loop QCD self-energy
  approx_twoloop = 3 is full known 2-loop self-energy
  approx_twoloop = 4 is 2-loop QCD self-energy with tree masses 
                     replaced by pole masses, 
  approx_twoloop = 5 is full known 2-loop self-energy with tree 
                     masses replaced by pole masses, 
  approx_twoloop = 6 is 2-loop QCD self-energy, plus 2-loop 
                     effective potential for the rest
  approx_twoloop = 7 is full known 2-loop self-energy, plus 2-loop 
                     effective potential for the rest
  approx_twoloop = 8 is 2-loop QCD self-energy, with tree masses
                     replaced by pole masses, plus 2-loop 
                     effective potential for the rest
  approx_twoloop = 9 is full known 2-loop self-energy, with tree
                     masses replaced by pole masses, plus 2-loop 
                     effective potential for the rest

For h0, the 32 valid combos for (approx_oneloop, approx_twoloop) are:
  tree-level = (0,0)
  one-loop = (i,0) for i=1,2,3,4
  two-loop = (i,j) for i=1,2,3 and j=1,2,3,4,5,6,7,8,9

For H0, A0, and G0 the 12 valid combos are:
  tree-level = (0,0)
  one-loop = (i,0) for i=2,3,4
  two-loop = (i,j) for i=2,3 and j = 2,3,4,5

The most likely to be used for h0 are probably: 
  (2,1) = used for solid line, Figure 1, 0211366,
          and dot-dashed line, Figure 3, 0405022,
          and dashed line, Figure 4 0405022 (with
          imaginary part for tau loops commented out
          of one-loop self-energy).
  (2,7) = used for solid line, Figure 3, 0405022,
          and solid line, Figure 4 0405022 (with
          imaginary part for tau loops commented out 
          of one-loop self-energy).
  (3,1) = fast 2-loop approximation
  (2,9) = nearly state-of-the-art but slow
  (3,9) = state-of-the-art but slow
*/

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
/* For now the only allowed options are:
   (0,0): tree
   (2,0): 1-loop self energy
   (2,1): 1-loop SE plus 2-loop EPA (fast)
   (2,7): 1-loop SE plus known 2-loop SE plus EPA for rest (slow)

   These are obtained for loop_order = 0, 1, 1.5 and 2, respectively.
*/

int SUMO_h0Pole (float loop_order)
{
  int approx_oneloop, approx_twoloop;
  char funcname[] = "SUMO_h0Pole";

  if (loop_order < 0.9)
    approx_oneloop = approx_twoloop = 0;
  else if (loop_order < 1.1) {
    approx_oneloop = 2; approx_twoloop = 0; }
  else if (loop_order < 1.6) {
    approx_oneloop = 2; approx_twoloop = 1; }
  else if (loop_order < 2.1) {
    approx_oneloop = 2; approx_twoloop = 7; }
  else
    SUMO_Error (funcname, "Invalid loop order specified.", 22);

  return SUMO_phi0_pole_masses (0, approx_oneloop, approx_twoloop);
}

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
/* For now the only allowed options are:
   (0,0): tree
   (2,0): 1-loop self energy
   (2,3): 1-loop SE plus known 2-loop SE

   These are obtained for loop_order = 0, 1, 2, respectively.
*/

int SUMO_H0Pole (int loop_order)
{
  int approx_oneloop = 2, approx_twoloop = 3;
  char funcname[] = "SUMO_H0Pole";

  if (loop_order < 0 || loop_order > 2)
    SUMO_Error (funcname, "Invalid loop order specified.", 22);

  if (loop_order == 1) 
    approx_twoloop = 0;
  else if (loop_order == 0)
    approx_oneloop = approx_twoloop = 0;

  return SUMO_phi0_pole_masses (1, approx_oneloop, approx_twoloop);
}

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
/* For now the only allowed options are:
   (0,0): tree
   (2,0): 1-loop self energy
   (2,3): 1-loop SE plus known 2-loop SE

   These are obtained for loop_order = 0, 1, 2, respectively.
*/

int SUMO_G0Pole (int loop_order)
{
  int approx_oneloop = 2, approx_twoloop = 3;
  char funcname[] = "SUMO_G0Pole";

  if (loop_order < 0 || loop_order > 2)
    SUMO_Error (funcname, "Invalid loop order specified.", 22);

  if (loop_order == 1) 
    approx_twoloop = 0;
  else if (loop_order == 0)
    approx_oneloop = approx_twoloop = 0;

  return SUMO_phi0_pole_masses (2, approx_oneloop, approx_twoloop);
}

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
/* For now the only allowed options are:
   (0,0): tree
   (2,0): 1-loop self energy
   (2,3): 1-loop SE plus known 2-loop SE

   These are obtained for loop_order = 0, 1, 2, respectively.
*/

int SUMO_A0Pole (int loop_order)
{
  int approx_oneloop = 2, approx_twoloop = 3;
  char funcname[] = "SUMO_A0Pole";

  if (loop_order < 0 || loop_order > 2)
    SUMO_Error (funcname, "Invalid loop order specified.", 22);

  if (loop_order == 1) 
    approx_twoloop = 0;
  else if (loop_order == 0)
    approx_oneloop = approx_twoloop = 0;

  return SUMO_phi0_pole_masses (3, approx_oneloop, approx_twoloop);
}

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
int SUMO_phi0_pole_masses (int which, int approx_oneloop, int approx_twoloop)
{
  char funcname[] = "SUMO_phi0_pole_masses";
  SUMO_COMPLEX se[3][2][2];
  SUMO_COMPLEX se_con_oneloop[2][2];
  SUMO_COMPLEX se_con_twoloop[2][2];
  SUMO_COMPLEX setot[2][2];
  SUMO_EIGEN_COMP a;
  int i,j,k,ii,jj;
  int offset = 0;

  if ((approx_oneloop < 0) || (approx_oneloop > 4))
    SUMO_Error(funcname, "Invalid one-loop approximation option specified.", 3);

  if ((approx_twoloop < 0) || (approx_twoloop > 9))
    SUMO_Error(funcname, "Invalid two-loop approximation option specified.", 3);

  if ((0 == approx_oneloop) && (approx_twoloop > 0))
    SUMO_Error(funcname, "Invalid combination: approx_oneloop=0 and approx_twoloop>0.", 3);

  if ((4 == approx_oneloop) && (approx_twoloop > 0))
    SUMO_Error(funcname, "Invalid combination: approx_oneloop=4 and approx_twoloop>0.", 3);

  if (which > 1) offset = 2;

  /* Tree-level masses: */
  se[0][0][0] = m2_phi0[offset];
  se[0][1][1] = m2_phi0[1+offset];
  se[0][0][1] = se[0][1][0] = 0;

  /* This could be restored if we want to leave EPA as a 1-loop option: */
  /* if (1 == approx_oneloop) { */
  /*   SUMO_V_derivs_numerical (1, 0.2L, 0.2L); */
  /*   se_con_oneloop[0][0] = 0.5L * (calpha * calpha * d2Vdvu2[1] */
  /*                       + salpha * salpha * d2Vdvd2[1]) */
  /*                 - calpha * salpha * d2Vdvuvd[1]; */
  /*   se_con_oneloop[1][1] = 0.5L * (calpha * calpha * d2Vdvd2[1] */
  /*                       + salpha * salpha * d2Vdvu2[1]) */
  /*                 + calpha * salpha * d2Vdvuvd[1]; */
  /*   se_con_oneloop[0][1] = 0.5L * (calpha * salpha * (d2Vdvu2[1] - d2Vdvd2[1]) + */
  /*                 (calpha * calpha - salpha * salpha) * d2Vdvuvd[1]); */
  /* } else { */
  /*   se_con_oneloop[0][0] = 0; */
  /*   se_con_oneloop[1][1] = 0; */
  /*   se_con_oneloop[0][1] = 0; */
  /* } */

  se_con_oneloop[0][0] = 0;
  se_con_oneloop[1][1] = 0;
  se_con_oneloop[0][1] = 0;

  if ((1 == approx_twoloop) || (approx_twoloop >= 6)) {
    SUMO_V_derivs_numerical (2, 0.2L, 0.2L);    
    se_con_twoloop[0][0] = 0.5L * (calpha * calpha * d2Vdvu2[2] 
                        + salpha * salpha * d2Vdvd2[2]) 
                  - calpha * salpha * d2Vdvuvd[2];
    se_con_twoloop[1][1] = 0.5L * (calpha * calpha * d2Vdvd2[2] 
                        + salpha * salpha * d2Vdvu2[2])
                  + calpha * salpha * d2Vdvuvd[2]; 
    se_con_twoloop[0][1] = 0.5L * (calpha * salpha * (d2Vdvu2[2] - d2Vdvd2[2]) +
                  (calpha * calpha - salpha * salpha) * d2Vdvuvd[2]);
  } else {
    se_con_twoloop[0][0] = 0;
    se_con_twoloop[1][1] = 0;
    se_con_twoloop[0][1] = 0;
  }

  /* Don't iterate if not necessary. SHOULD NOT BE NECESSARY FOR THESE OPTIONS */
  /* if ((2 == which) || ((approx_oneloop < 2) && (approx_twoloop < 2))) maxiters = 1;   */

  if (2 == which) M2_phi0[which] = 0;

  for (i=0; i<2; i++) {
    for (j=i; j<2; j++) {

      se[1][i][j] = se_con_oneloop[i][j];
      se[2][i][j] = se_con_twoloop[i][j];

      ii = i+offset;
      jj = j+offset;

      if (2 == approx_oneloop) 
	se[1][i][j] = SUMO_oneloopfactor * pi1_phi0 (ii, jj, m2_phi0[which]);

      if ((3 == approx_twoloop) || (7 == approx_twoloop))
	se[2][i][j] += SUMO_twoloopfactor * pi2QCD_phi0 (ii, jj, m2_phi0[which])
	             + SUMO_twoloopfactor * pi2nonQCD_phi0 (ii, jj, m2_phi0[which]);

      if (7 == approx_twoloop)
	se[2][i][j] += - SUMO_twoloopfactor * pi2QCD_phi0 (ii, jj, 0)
	               - SUMO_twoloopfactor * pi2nonQCD_phi0 (ii, jj, 0);
    }}

  se[1][1][0] = se[1][0][1];
  se[2][1][0] = se[2][0][1];

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      setot[i][j] = se[0][i][j] + se[1][i][j] + se[2][i][j];
    }}

  SUMO_DiagonalizeComp (setot[0], 2, &a);
  CM2_phi0[which] = a.eigenvalues[which-offset];

  M2_phi0[which] = SUMO_CREAL(CM2_phi0[which]);

  if (2 == which) {
    if (SUMO_FABS(M2_phi0[which]) < TSIL_TOL)
      M2_phi0[which] = TSIL_TOL;
  } else if (M2_phi0[which] < 0.0001)
    M2_phi0[which] = 0.0001;

  /* M_phi0 = SUMO_SGNSQRT(M2_phi0[which]); */
  /* Gamma_phi0 = -SUMO_CIMAG(CM2_phi0[which])/M_phi0; */

  return 0;
}
