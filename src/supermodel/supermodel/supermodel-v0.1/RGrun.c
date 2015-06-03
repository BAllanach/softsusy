/* RG running of the basic paramters. */

#include "supermodel.h"

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Runs all MSSM parameters to Q_final at the specified loop order. */

int SUMO_RGrun (TSIL_REAL  Q_final, 
		int        loop_order)
{
  TSIL_REAL t_final, t, dt;
  int min_steps, max_steps;
  int force_step;
  int goodsteps, badsteps; 
  int rk6status; /* 1 for success or forced; 0 for need retry. */
  TSIL_REAL Q_init = Q;

  t_final = TSIL_LOG(Q_final/Q_init);

  if (TSIL_FABS(t_final) < TSIL_TOL) return 0;
  min_steps = 1 + ((int) (5.0L * TSIL_FABS(t_final)));
  max_steps = 1 + ((int) (125.0L * TSIL_FABS(t_final)));

  dt = t_final/(1 + ((int) (124.0L * TSIL_FABS(t_final))));
  t = 0.0L;

  SUMO_Betas (loop_order);

  goodsteps = badsteps = 0;   

  while ( TSIL_FABS(dt) < 0.5 * TSIL_FABS(t_final - t) ) 
    {
      if ( TSIL_FABS(dt) < TSIL_FABS(t_final/max_steps) ) 
	{
	  force_step = 1;
	  dt = t_final/((TSIL_REAL) max_steps);
	}
      else force_step = 0;
      
      if ( TSIL_FABS(dt) > TSIL_FABS(t_final)/min_steps )
	dt = (t_final)/min_steps;

      rk6status = RG_rk6 (&dt, loop_order, force_step);
      t = TSIL_LOG ((Q)/Q_init);
      
      if (1 == rk6status) {        
	goodsteps += (1-force_step);
	badsteps += force_step;
      }
    }

  /*
    printf("goodsteps = %d\n",goodsteps);
    printf("badsteps  = %d\n",badsteps);
    printf("Don't be surprised if there are a lot of badsteps... \n");
  */

  /* The remaining distance is less than 2.0 times the step size.  So,
     take exactly two more steps. */

  dt = 0.5L*(t_final - t);
  RG_rk6 (&dt, loop_order, 1);
  t = TSIL_LOG ((Q)/Q_init);

  /* Arrange final step to land exactly on t_final, and force it. */
  dt = t_final - t;
  RG_rk6 (&dt, loop_order, 1);

/*   is_updated = 0; */
/*   are_tree_masses_updated = 0; */
/*   are_tree_couplings_updated = 0; */
  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();
  SUMO_Pole_masses_from_tree ();

  /* Return a status code eventually */
  return 0;
}
