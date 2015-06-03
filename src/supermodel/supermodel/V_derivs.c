/* Routines for computing the 1st and 2nd derivatives of the tree,
   1-loop, and 2-loop effective potentials. Uses a 5x5 grid (except
   the 4 corner points, which aren't used, so actually 21 evaluations
   of the potential are done), with step sizes in vu and vd as inputs.
*/

#include "supermodel.h"

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Derivatives of the tree-level potential are trivial, so do them
   analytically. */

int SUMO_Vtree_derivs (void)
{
  TSIL_REAL mmu, mmd;

  if (0 == is_updated) SUMO_Update();

  mmu = m2Hu + mu2;
  mmd = m2Hd + mu2;

  dVdvu[0] = 2.0L * (mmd * vd - b * vu) 
           + 0.5L * vd * (vd * vd - vu * vu) * g2plusgp2;

  dVdvd[0] = 2.0L * (mmu * vu - b * vd) 
           + 0.5L * vu * (vd * vd - vu * vu) * g2plusgp2;

  d2Vdvd2[0] = 2.0L * mmd + g2plusgp2 * (1.5L * vd * vd - 0.5L * vu * vu);

  d2Vdvu2[0] = 2.0L * mmu + g2plusgp2 * (1.5L * vu * vu - 0.5L * vd * vd);

  d2Vdvuvd[0] = -2.0L * b - g2plusgp2 * vu * vd; 

  return 0;
}
/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

int SUMO_V_derivs_numerical (int loop_order, 
			     TSIL_REAL delta_vd, 
			     TSIL_REAL delta_vu)
{
  int i,j,k;
  TSIL_REAL Vgrid[5][5];
  int ijlist[21][2] = {         {-2,-1}, {-2,0}, {-2,1},
                       {-1,-2}, {-1,-1}, {-1,0}, {-1,1}, {-1,2},
                       {0, -2},  {0,-1},          {0,1},  {0,2},
                       {1, -2},  {1,-1},  {1,0},  {1,1},  {1,2},
                                 {2,-1},  {2,0},  {2,1},
                       {0,0}}; 
  TSIL_REAL vu_0 = vu;
  TSIL_REAL vd_0 = vd;
  char funcname[] = "SUMO_V_derivs_numerical";

  Veff_loop_order = loop_order;

  for (k=0; k<21; k++) {
    i = ijlist[k][0];
    j = ijlist[k][1];
    vd = vd_0 + delta_vd * ((TSIL_REAL) i);
    vu = vu_0 + delta_vu * ((TSIL_REAL) j);

    is_updated = 0;
    are_tree_masses_updated = 0;
    are_tree_couplings_updated = 0;
    
    if (2 == loop_order) Vgrid[i+2][j+2] = TSIL_CREAL(SUMO_V_twoloop ());
    else if (1 == loop_order) Vgrid[i+2][j+2] = TSIL_CREAL(SUMO_V_oneloop ()); 
    else if (0 == loop_order) Vgrid[i+2][j+2] = TSIL_CREAL(SUMO_V_tree ());    
    else TSIL_Error(funcname,
                    "Loop order must be 0, 1, or 2.", 
		    2);
  }

  dVdvd[loop_order] = (8.0L * (Vgrid[3][2] - Vgrid[1][2])
                 + (Vgrid[0][2] - Vgrid[4][2]))/(12.0L * delta_vd);

  dVdvu[loop_order] = (8.0L * (Vgrid[2][3] - Vgrid[2][1])
                 + (Vgrid[2][0] - Vgrid[2][4]))/(12.0L * delta_vu);

  d2Vdvd2[loop_order] = (16.0L * (Vgrid[3][2] + Vgrid[1][2])
                     - Vgrid[4][2] - Vgrid[0][2]
             - 30.0L * Vgrid[2][2])/(12.0L * delta_vd * delta_vd);

  d2Vdvu2[loop_order] = (16.0L * (Vgrid[2][3] + Vgrid[2][1])
                     - Vgrid[2][4] - Vgrid[2][0]              
             - 30.0L * Vgrid[2][2])/(12.0L * delta_vu * delta_vu);

  d2Vdvuvd[loop_order] = (10.0L * (Vgrid[3][3] + Vgrid[1][1]
                      - Vgrid[3][1] - Vgrid[1][3])
              + Vgrid[4][1] + Vgrid[0][3] + Vgrid[3][0] + Vgrid[1][4]
              - Vgrid[4][3] - Vgrid[0][1] - Vgrid[3][4] - Vgrid[1][0]
              )/(24.0L * delta_vu * delta_vd);

  vu = vu_0;
  vd = vd_0;
  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  return 0;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

int SUMO_V_firstderivs_numerical (int loop_order, 
				  TSIL_REAL delta_vd, 
				  TSIL_REAL delta_vu)
{
  int i,j,k;
  TSIL_REAL Vgrid[12];  
  int ijlist[12][2] = {{-3,0},{-2,0},{-1,0},{1,0},{2,0},{3,0},
                       {0,-3},{0,-2},{0,-1},{0,1},{0,2},{0,3}};
  TSIL_REAL vd_0 = vd;
  TSIL_REAL vu_0 = vu;
  char funcname[] = "SUMO_V_firstderivs_numerical";

  for (k=0; k<12; k++) {
    i = ijlist[k][0];
    j = ijlist[k][1];
    vd = vd_0 + delta_vd * ((TSIL_REAL) i);
    vu = vu_0 + delta_vu * ((TSIL_REAL) j);
    is_updated = 0;
    are_tree_masses_updated = 0;
    are_tree_couplings_updated = 0;

    if (2 == loop_order) Vgrid[k] = TSIL_CREAL(SUMO_V_twoloop ());     
    else if (1 == loop_order) Vgrid[k] = TSIL_CREAL(SUMO_V_oneloop ()); 
    else if (0 == loop_order) Vgrid[k] = TSIL_CREAL(SUMO_V_tree ());    
    else TSIL_Error(funcname,
                    "Loop order must be 0, 1, or 2.", 2);
  }

  dVdvd[loop_order] = (45.0L * (Vgrid[3] - Vgrid[2]) 
		       + 9.0L * (Vgrid[1] - Vgrid[4])
		       + Vgrid[5] - Vgrid[0])/(60.0L * delta_vd);

  dVdvu[loop_order] = (45.0L * (Vgrid[9] - Vgrid[8]) 
		       + 9.0L * (Vgrid[7] - Vgrid[10])
		       + Vgrid[11] - Vgrid[6])/(60.0L * delta_vu);

  vd = vd_0;
  vu = vu_0;
  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  return 0;
}
