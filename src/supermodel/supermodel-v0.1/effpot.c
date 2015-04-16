/* Compute effective potential through two loops.  Also included is
   the wrapper function needed by the minimization routines. */

#include "supermodel.h"

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Wrapper function that is called by the minimization routine. */

TSIL_REAL potfunc (TSIL_REAL v[2])
{
  v[0] = TSIL_FABS(v[0]);
  v[1] = TSIL_FABS(v[1]);

  vu = v[0];
  vd = v[1];

  return SUMO_CREAL
    (SUMO_Veff (v[0], v[1], Veff_loop_order));
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Minimizes Veff as a function of (vu,vd) and updates Veff as well as
   (vu,vd) in the model struct. */

void SUMO_Minimize_Veff (int loop_order)
{
  char funcname[] = "SUMO_Minimize_Veff";

  int ndim = 2;
  int i, iter;
  TSIL_REAL q[2], **xi, fret, bar[2][2];
  TSIL_REAL ftol = TSIL_TOL;

  /* For simplex method minimization */
  TSIL_REAL **p, gaak[3][2];
  TSIL_REAL y[3];
  int nfunc;
  TSIL_REAL lam = 27.0;
  /* End of declarations for simplex method */

  Veff_loop_order = loop_order;

  /* Starting point for Powell search */
  q[0] = 175.0*sinbeta;
  q[1] = 175.0*cosbeta;
  
  xi = (TSIL_REAL **) calloc (ndim, sizeof(TSIL_REAL *));
  for (i=0; i<2; i++)
    xi[i] = bar[i];

  /* Starting directions are just coordinate axes */
  xi[0][0] = 1.0;
  xi[0][1] = 0.0;
  xi[1][0] = 0.0;
  xi[1][1] = 1.0;
  
  PowellMin (q, xi, ndim, ftol, &iter, &fret, potfunc);

  printf("Minimizing effective potential (Powell)...\n");
  printf("Minimum at (vu,vd) = (%Lf, %Lf)\tFunc = %Le\n", q[0], q[1], fret);
  printf("We did %d iterations\n", iter);

  vu = TSIL_FABS(q[0]);
  vd = TSIL_FABS(q[1]);
  tanbeta = vu/vd;
  SUMO_Evaluate_Veff (loop_order);
  SUMO_Pole_masses_from_tree ();
  
  /* ===== Test with simplex method ===== */
/*   p = (TSIL_REAL **) calloc (3, sizeof(TSIL_REAL *)); */

/*   for (i=0; i<3; i++) */
/*     p[i] = gaak[i]; */

/*   /\* Set up initial simplex *\/ */
/*   p[0][0] = 100.0; */
/*   p[0][1] = 301.0; */

/*   p[1][0] = p[0][0]; */
/*   p[1][1] = p[0][1] + lam; */

/*   p[2][0] = p[0][0] + lam; */
/*   p[2][1] = p[0][1]; */
  
/*   for (i=0; i<3; i++) */
/*     y[i] = potfunc (p[i]); */

/*   ftol = TSIL_SQRT(ftol); */
/*   ftol = 1.e-6; */
/*   SimplexMin (p, y, ndim, ftol, potfunc, &nfunc); */

/*   printf("\nSimplex Method (2D) Results:\n"); */
/*   for (i=0; i<3; i++) */
/*     printf("Vertex %d: (%Lf, %Lf)\tFunc = %Le\n", */
/* 	   i, p[i][0], p[i][1], y[i]); */
/*   printf("We did %d function evaluations\n", nfunc); */

/*   model->vu = TSIL_FABS(p[0][0]); */
/*   model->vd = TSIL_FABS(p[0][1]); */
/*   model->tanbeta = (model->vu)/(model->vd); */
/*   SUMO_Evaluate_Veff (model, loop_order); */

  /* ===== End of test with simplex method ===== */

}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Evaluates tree, one- and two-loop contributions to the effective
   potential, setting the values in the model struct. */

int SUMO_Evaluate_Veff (int loop_order)
{
  char funcname[] = "SUMO_Evaluate_Veff";
  SUMO_COMPLEX res;

  if (loop_order < 0 || loop_order > 2)
    SUMO_Error(funcname, "Invalid loop order.", 5);

  Veff_loop_order = loop_order;
  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  res = SUMO_V_tree ();
  if (loop_order > 0) res += SUMO_V_oneloop ();
  if (loop_order > 1) res += SUMO_V_twoloop ();

  /* Put this back eventually?  */
  /* SUMO_V_field_independent (model);  */

  Veff = res;

  return 0;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Returns Veff (vu,vd). Leaves original values of vu,vd in model 
   unchanged.  */

SUMO_COMPLEX SUMO_Veff (SUMO_REAL vu_in, 
                        SUMO_REAL vd_in,
                        int loop_order)
{
  SUMO_REAL vu_orig, vd_orig;
  SUMO_COMPLEX res;
  char funcname[] = "SUMO_Veff";

  if (loop_order < 0 || loop_order > 2)
    SUMO_Error(funcname, "Invalid loop order specified.", 5);

  /* Save original values */
  vu_orig = vu;
  vd_orig = vd;
  vu = vu_in;
  vd = vd_in;

  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

/*   model->is_updated = 0; */
/*   model->are_tree_masses_updated = 0; */
/*   model->are_tree_couplings_updated = 0; */
  
  res = SUMO_V_tree ();
  if (loop_order > 0) res += SUMO_V_oneloop ();
  if (loop_order > 1) res += SUMO_V_twoloop ();

  vu = vu_orig;
  vd = vd_orig;

  /* Why not just do this, so the model is left as it (perhaps)
     was? */
  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

/*   model->is_updated = 0; */
/*   model->are_tree_masses_updated = 0; */
/*   model->are_tree_couplings_updated = 0; */

  return res;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------- */
/* Tree-level effective potential. */

SUMO_COMPLEX SUMO_V_tree () 
{
  SUMO_COMPLEX res;

  if (NO == is_updated) SUMO_Update ();
  
  /* Eventually this should be fixed in SUMO_Update, but first need to
     figure out phase rotations of other parameters. */
  if (TSIL_CABS(1.0L - b/TSIL_CABS(b)) > TSIL_TOL)
    printf("Error: b must be real, positive to evaluate V.");

  res = -2.L * TSIL_CABS(b) * vuvd
        + (mu2 + m2Hu) * vu2
        + (mu2 + m2Hd) * vd2
    + g2plusgp2 * (vu2 - vd2) * (vu2 - vd2)/8.0L;

  return res;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* One-loop contribution to effective potential. */

SUMO_COMPLEX SUMO_V_oneloop () 
{
  int i;
  TSIL_COMPLEX res = 0.0;
  char funcname[] = "SUMO_V_oneloop";

  if (NO == are_tree_masses_updated) SUMO_Tree_Masses ();

  /* DGR appears that the gluino is missing? ADDED 5/10 */
  res += -16.0L * SUMO_OneLoopVacuum (m2_gluino, Q2);

  res +=
    3.0L * SUMO_OneLoopVacuum (m2_Z, Q2)
    + 6.0L * SUMO_OneLoopVacuum (m2_W, Q2)
    - 12.0L * SUMO_OneLoopVacuum (m2_top, Q2)
    - 12.0L * SUMO_OneLoopVacuum (m2_bot, Q2)
    - 4.0L * SUMO_OneLoopVacuum (m2_tau, Q2)
    ;

  for (i=0; i<4; i++)
    res += -2.0L * SUMO_OneLoopVacuum (m2_Neut[i], Q2);

  for (i=0; i<2; i++)
    res += -4.0L * SUMO_OneLoopVacuum (m2_Char[i], Q2);

  for (i=0; i<4; i++)
    res += SUMO_OneLoopVacuum (m2_phi0[i], Q2);

  for (i=0; i<2; i++)
    res += 2.0L * SUMO_OneLoopVacuum (m2_phip[i], Q2);

  for (i=0; i<2; i++) {
    res += 6.0L * SUMO_OneLoopVacuum (m2_stop[i], Q2);
    res += 6.0L * SUMO_OneLoopVacuum (m2_sbot[i], Q2);
    res += 2.0L * SUMO_OneLoopVacuum (m2_stau[i], Q2);
  }

  for (i=0; i<2; i++) {
    res += 6.0L * SUMO_OneLoopVacuum (m2_suL[i], Q2);
    res += 6.0L * SUMO_OneLoopVacuum (m2_sdL[i], Q2);
    res += 6.0L * SUMO_OneLoopVacuum (m2_suR[i], Q2);
    res += 6.0L * SUMO_OneLoopVacuum (m2_sdR[i], Q2);
    res += 2.0L * SUMO_OneLoopVacuum (m2_seL[i], Q2);
    res += 2.0L * SUMO_OneLoopVacuum (m2_seR[i], Q2);
  }

  for (i=0; i<3; i++)
    res += 2.0L * SUMO_OneLoopVacuum (m2_snu[i], Q2);

  res *= SUMO_oneloopfactor;

  return res;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Two-loop contributions to the effective potential. */
/* Ref: hep-ph/0206136 */

SUMO_COMPLEX SUMO_V_twoloop () 
{
  SUMO_COMPLEX res = 0.0;

  if (NO == are_tree_couplings_updated) SUMO_Tree_Couplings ();

  res += SUMO_V2_strong ();
  res += SUMO_V2_SSS ();
  res += SUMO_V2_SS ();
  res += SUMO_V2_FFS ();
  res += SUMO_V2_SSV ();
  res += SUMO_V2_VS ();
  res += SUMO_V2_VVS ();
  res += SUMO_V2_FFV ();
  res += SUMO_V2_gauge ();
  
  return res;
}
