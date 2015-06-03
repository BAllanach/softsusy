/* ---------------------------------------------------------------- */
/* Implements a 6-stage, 5th-order Cash-Karp-Fehlberg Runge-Kutta step
   for renormalization group running.  The purpose of this (as opposed
   to classical 4-stage 4th-order Runge-Kutta) is to enable efficient
   adaptive step size control.

   The criterion for deciding whether the step size is small enough is
   presently:

   |(5th order change) - (4th order embedded change)|
   < (delta_err_fac * |5th order change| + |value|) 
     * val_err_fac * TSIL_EPSILON

   with delta_err_fac and val_err_fac #defined below.

   Otherwise, the step is rejected, unless force_step == 1.

   Provided it finds the step size to be small enough, OR the argument
   force_step==1, then RG_rk6 will return 1 and increment the RG scale
   from Q to Q exp(*dt).

   The dependent variables and their derivatives in *foo are also
   updated in preparation for the next step.

   If the step size needs to be decreased, then RG_rk6 will return 0,
   and not increment the independent variable and not change the
   dependent variables or their derivatives.

   RG_rk6 also replaces *dt with an estimate of the optimal step size,
   regardless of whether that is going to be a retry of this step if
   it failed, or the next step if this one passed. But, it never tries
   to increase or decrease the step size by more than a factor of
   2. And, the calling function can and will reject the suggested step
   size if it gets to be too small or too large. */
/* ---------------------------------------------------------------- */

#include "supermodel.h"

#define SafetyFactor 0.9L
#define delta_err_fac 2000.L
#define val_err_fac 5000.L  /* These can be adjusted. */

/* Here are the Butcher coefficients for Cash-Karp-Fehlberg: */

/* Confusingly, these are Numerical Recipes a_i */
#define ButchCKFc2 0.2L
#define ButchCKFc3 0.3L
#define ButchCKFc4 0.6L
#define ButchCKFc5 1.0L
#define ButchCKFc6 0.875L

/* Confusingly, these are Numerical Recipes b_ij */
#define ButchCKFa21 0.2L
#define ButchCKFa31 0.075L
#define ButchCKFa32 0.225L
#define ButchCKFa41 0.3L
#define ButchCKFa42 -0.9L
#define ButchCKFa43 1.2L
#define ButchCKFa51 -0.203703703703703703703703703704L
#define ButchCKFa52 2.5L
#define ButchCKFa53 -2.59259259259259259259259259259L
#define ButchCKFa54 1.29629629629629629629629629630L
#define ButchCKFa61 0.0294958043981481481481481481481L
#define ButchCKFa62 0.341796875L
#define ButchCKFa63 0.0415943287037037037037037037037L
#define ButchCKFa64 0.400345413773148148148148148148L
#define ButchCKFa65 0.061767578125L

/* Confusingly, these are Numerical Recipes c_i */
#define ButchCKFb1 0.0978835978835978835978835978836L
#define ButchCKFb2 0.0L
#define ButchCKFb3 0.402576489533011272141706924316L 
#define ButchCKFb4 0.210437710437710437710437710438L
#define ButchCKFb5 0.0L
#define ButchCKFb6 0.289102202145680406549971767363L

/* These are Numerical Recipes c_i - c_i^* */
#define ButchCKFe1 -0.00429377480158730158730158730159L
#define ButchCKFe2 0.0L 
#define ButchCKFe3 0.0186685860938578329882677708765L 
#define ButchCKFe4 -0.0341550268308080808080808080808L 
#define ButchCKFe5 -0.0193219866071428571428571428571L 
#define ButchCKFe6 0.0391022021456804065499717673631L 

#define N_Rparams 25
#define N_Cparams 8
#define N_Rparams_gYv 8

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

int RG_rk6 (TSIL_REAL *dt,
	    int loop_order,
	    int force_step)
{
  int i;

  TSIL_REAL *Rparams[N_Rparams] = {
      &gp, &g, &g3,
      &ytop, &ybot, &ytau,
      &m2Q[0], &m2Q[1], &m2Q[2],
      &m2L[0], &m2L[1], &m2L[2],
      &m2u[0], &m2u[1], &m2u[2],
      &m2d[0], &m2d[1], &m2d[2],
      &m2e[0], &m2e[1], &m2e[2],
      &m2Hu, &m2Hd,
      &vu, &vd};

  TSIL_COMPLEX *Cparams[N_Cparams] = { 
      &mu,
      &M1, &M2, &M3,
      &atop, &abot, &atau,
      &b};

  TSIL_REAL *Rparams_beta[N_Rparams] = {
      &beta_gp, &beta_g, &beta_g3,
      &beta_ytop, &beta_ybot, &beta_ytau,
      &beta_m2Q12, &beta_m2Q12, &beta_m2Q3,
      &beta_m2L12, &beta_m2L12, &beta_m2L3,
      &beta_m2u12, &beta_m2u12, &beta_m2u3,
      &beta_m2d12, &beta_m2d12, &beta_m2d3,
      &beta_m2e12, &beta_m2e12, &beta_m2e3,
      &beta_m2Hu, &beta_m2Hd,
      &beta_vu, &beta_vd};

  TSIL_COMPLEX *Cparams_beta[N_Cparams] = { 
      &beta_mu,
      &beta_M1, &beta_M2, &beta_M3,
      &beta_atop, &beta_abot, &beta_atau,
      &beta_b};

  TSIL_REAL    Rparams_delta[N_Rparams];
  TSIL_COMPLEX Cparams_delta[N_Cparams];

  TSIL_REAL    Rparams_err[N_Rparams];
  TSIL_REAL    Cparams_err[N_Cparams];

  TSIL_REAL    k1_Rparams[N_Rparams];
  TSIL_COMPLEX k1_Cparams[N_Cparams];

  TSIL_REAL    k2_Rparams[N_Rparams];
  TSIL_COMPLEX k2_Cparams[N_Cparams];

  TSIL_REAL    k3_Rparams[N_Rparams];
  TSIL_COMPLEX k3_Cparams[N_Cparams];

  TSIL_REAL    k4_Rparams[N_Rparams];
  TSIL_COMPLEX k4_Cparams[N_Cparams];

  TSIL_REAL    k5_Rparams[N_Rparams];
  TSIL_COMPLEX k5_Cparams[N_Cparams];

  TSIL_REAL    k6_Rparams[N_Rparams];
  TSIL_COMPLEX k6_Cparams[N_Cparams];

  TSIL_REAL next_step_size;
  TSIL_REAL temp, maxerr; 
  TSIL_REAL step_rescale = 0.19L; /* This can only increase. */
  int status = 0; 

  /* Set the starting values */

  TSIL_REAL Rparams_start[25] = {
      gp, g, g3,
      ytop, ybot, ytau,
      m2Q[0], m2Q[1], m2Q[2],
      m2L[0], m2L[1], m2L[2],
      m2u[0], m2u[1], m2u[2],
      m2d[0], m2d[1], m2d[2],
      m2e[0], m2e[1], m2e[2],
      m2Hu, m2Hd,
      vu, vd};

  TSIL_COMPLEX Cparams_start[8] = { 
      mu,
      M1, M2, M3,
      atop, abot, atau,
      b};

  /* Find k1 values using existing values for beta functions: */

  for (i=0; i < N_Rparams; i++)
    k1_Rparams[i] = (*dt) * (*(Rparams_beta[i]));

  for (i=0; i < N_Cparams; i++)
    k1_Cparams[i] = (*dt) * (*(Cparams_beta[i]));

  /* BEGIN STAGE 2. */

  /* Adjust data values */
  for (i=0; i < N_Rparams; i++)
    *(Rparams[i]) = Rparams_start[i] + ButchCKFa21 * k1_Rparams[i];

  for (i=0; i < N_Cparams; i++)
    *(Cparams[i]) = Cparams_start[i] * ButchCKFa21 * k1_Cparams[i];

  /* Find k2 values: */

  SUMO_Betas (loop_order);

  for (i=0; i < N_Rparams; i++)
    k2_Rparams[i] = (*dt) * (*(Rparams_beta[i]));

  for (i=0; i < N_Cparams; i++)
    k2_Cparams[i] = (*dt) * (*(Cparams_beta[i]));

  /* BEGIN STAGE 3. */

  /* Adjust data values */
  for (i=0; i < N_Rparams; i++)
    *(Rparams[i]) = Rparams_start[i] + ButchCKFa31 * k1_Rparams[i]
                                     + ButchCKFa32 * k2_Rparams[i];

  for (i=0; i < N_Cparams; i++)
    *(Cparams[i]) = Cparams_start[i] + ButchCKFa31 * k1_Cparams[i]
                                     + ButchCKFa32 * k2_Cparams[i];

  /* Find k3 values: */

  SUMO_Betas (loop_order);

  for (i=0; i < N_Rparams; i++)
    k3_Rparams[i] = (*dt) * (*(Rparams_beta[i]));

  for (i=0; i < N_Cparams; i++)
    k3_Cparams[i] = (*dt) * (*(Cparams_beta[i]));

  /* BEGIN STAGE 4. */

  for (i=0; i < N_Rparams; i++)
    *(Rparams[i]) = Rparams_start[i] + ButchCKFa41 * k1_Rparams[i]
                                     + ButchCKFa42 * k2_Rparams[i]
                                     + ButchCKFa43 * k3_Rparams[i];

  for (i=0; i < N_Cparams; i++)
    *(Cparams[i]) = Cparams_start[i] + ButchCKFa41 * k1_Cparams[i]
                                     + ButchCKFa42 * k2_Cparams[i]
                                     + ButchCKFa43 * k3_Cparams[i];

  /* Find k4 values: */

  SUMO_Betas (loop_order);

  for (i=0; i < N_Rparams; i++)
    k4_Rparams[i] = (*dt) * (*(Rparams_beta[i]));

  for (i=0; i < N_Cparams; i++)
    k4_Cparams[i] = (*dt) * (*(Cparams_beta[i]));

  /* BEGIN STAGE 5. */

  for (i=0; i < N_Rparams; i++)
    *(Rparams[i]) = Rparams_start[i] + ButchCKFa51 * k1_Rparams[i]
                                     + ButchCKFa52 * k2_Rparams[i]
                                     + ButchCKFa53 * k3_Rparams[i]
                                     + ButchCKFa54 * k4_Rparams[i];

  for (i=0; i < N_Cparams; i++)
    *(Cparams[i]) = Cparams_start[i] + ButchCKFa51 * k1_Cparams[i]
                                     + ButchCKFa52 * k2_Cparams[i]
                                     + ButchCKFa53 * k3_Cparams[i]
                                     + ButchCKFa54 * k4_Cparams[i];

  /* Find k5 values: */

  SUMO_Betas (loop_order);

  for (i=0; i < N_Rparams; i++)
    k5_Rparams[i] = (*dt) * (*(Rparams_beta[i]));

  for (i=0; i < N_Cparams; i++)
    k5_Cparams[i] = (*dt) * (*(Cparams_beta[i]));

  /* BEGIN STAGE 6. */

  for (i=0; i < N_Rparams; i++)
    *(Rparams[i]) = Rparams_start[i] + ButchCKFa61 * k1_Rparams[i]
                                     + ButchCKFa62 * k2_Rparams[i]
                                     + ButchCKFa63 * k3_Rparams[i]
                                     + ButchCKFa64 * k4_Rparams[i]
                                     + ButchCKFa65 * k5_Rparams[i];

  for (i=0; i < N_Cparams; i++)
    *(Cparams[i]) = Cparams_start[i] + ButchCKFa61 * k1_Cparams[i]
                                     + ButchCKFa62 * k2_Cparams[i]
                                     + ButchCKFa63 * k3_Cparams[i]
                                     + ButchCKFa64 * k4_Cparams[i]
                                     + ButchCKFa65 * k5_Cparams[i];

  /* Find k6 values: */

  SUMO_Betas (loop_order);

  for (i=0; i < N_Rparams; i++)
    k6_Rparams[i] = (*dt) * (*(Rparams_beta[i]));

  for (i=0; i < N_Cparams; i++)
    k6_Cparams[i] = (*dt) * (*(Cparams_beta[i]));

  /* DONE WITH STAGES. */

  /* Compute 5th-order result changes in data values */

  for (i=0; i < N_Rparams; i++)
    Rparams_delta[i] =   ButchCKFb1 * k1_Rparams[i]
                       + ButchCKFb2 * k2_Rparams[i]
                       + ButchCKFb3 * k3_Rparams[i]
                       + ButchCKFb4 * k4_Rparams[i]
                       + ButchCKFb5 * k5_Rparams[i]
                       + ButchCKFb6 * k6_Rparams[i];

  for (i=0; i < N_Cparams; i++)
    Cparams_delta[i] =   ButchCKFb1 * k1_Cparams[i]
                       + ButchCKFb2 * k2_Cparams[i]
                       + ButchCKFb3 * k3_Cparams[i]
                       + ButchCKFb4 * k4_Cparams[i]
                       + ButchCKFb5 * k5_Cparams[i]
                       + ButchCKFb6 * k6_Cparams[i];

  /* Compute estimated error of each dependent variable, and keep
     track of the maximum estimated step rescaling found. */

  for (i=0; i < N_Rparams; i++)
    Rparams_err[i] = TSIL_FABS(ButchCKFe1 * k1_Rparams[i]
                             + ButchCKFe2 * k2_Rparams[i]
                             + ButchCKFe3 * k3_Rparams[i]
                             + ButchCKFe4 * k4_Rparams[i]
                             + ButchCKFe5 * k5_Rparams[i]
                             + ButchCKFe6 * k6_Rparams[i]);

  for (i=0; i < N_Cparams; i++)
    Cparams_err[i] = TSIL_CABS(ButchCKFe1 * k1_Cparams[i]
                             + ButchCKFe2 * k2_Cparams[i]
                             + ButchCKFe3 * k3_Cparams[i]
                             + ButchCKFe4 * k4_Cparams[i]
                             + ButchCKFe5 * k5_Cparams[i]
                             + ButchCKFe6 * k6_Cparams[i]);

  maxerr = 0.0L;

  for (i=0; i < N_Rparams; i++){
    temp = Rparams_err[i]/(
           delta_err_fac * TSIL_FABS(Rparams_delta[i]) 
           + TSIL_FABS(*(Rparams[i])) + TSIL_EPSILON);

    if (temp > maxerr) maxerr = temp;	
  }

  for (i=0; i < N_Cparams; i++){
    temp = Cparams_err[i]/(
           delta_err_fac * TSIL_CABS(Cparams_delta[i]) 
           + TSIL_CABS(*(Cparams[i])) + TSIL_EPSILON);

    if (temp > maxerr) maxerr = temp;	
  }

  maxerr = maxerr/(val_err_fac * TSIL_EPSILON);

  /* Set data back to original values in preparation for possible update. */

  for (i=0; i < N_Rparams; i++)
    *(Rparams[i]) = Rparams_start[i];

  for (i=0; i < N_Cparams; i++)
    *(Cparams[i]) = Cparams_start[i];

  /* Now, if the error was acceptable, OR the step is being forced, we
     increment the independent variable and the data values and the
     derivatives, and get set to report status = success. */

  if ((maxerr < 1.0L) || (1 == force_step)) 
    {
      status = 1;

      /* Update parameter values, set up beta functions for next
	 step. */

      Q *= TSIL_EXP(*dt);

      for (i=0; i < N_Rparams; i++)
	*(Rparams[i]) += Rparams_delta[i];

      for (i=0; i < N_Cparams; i++)
	*(Cparams[i]) += Cparams_delta[i];
    }
  else 
    {
      status = 0;
    }

  SUMO_Betas (loop_order); 

  /* Predict the appropriate next step size. */
  step_rescale = TSIL_SQRT(TSIL_SQRT(1.0L/maxerr));

  next_step_size = (*dt) * (SafetyFactor * step_rescale);

  /* Don't let the next step size be bigger or smaller than the
     present size by more than a factor of 4. */

  if (TSIL_FABS(next_step_size) > 4.0L*TSIL_FABS(*dt))
    next_step_size = 4.0L * (*dt);

  if (TSIL_FABS(next_step_size) < 0.25L*TSIL_FABS(*dt))
    next_step_size = 0.25L * (*dt);

  /* Recommend the new step size to the calling function. */
  *dt = next_step_size;

  return status;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

