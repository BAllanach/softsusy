/* Miscellaneous functions */

#include "supermodel.h"

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Basic initialization of the model struct. This should be called
   first thing on any such struct, before doing any calculation. In
   normal operation, it is called by SUMO_Read_SLHA. */

#include "sumo_smdata.h"

int SUMO_Initialize ()
{
  /* The most generic values here: */
  /* particleContent = MSSM; */
  /* modelType       = GENERAL; */

  is_updated = 0;
  are_tree_masses_updated = 0;
  are_tree_couplings_updated = 0;

  /* Set the SM data from defaults. The input SLHA file will typically
     override most or all of these except M_bot_phys and
     M2_bot_phys. */
  alpha_MZ_MS  = 1.0L/ALPHA_MZ_INV_DEFAULT;
  alpha0       = 1.0L/ALPHA_0_INV_DEFAULT;
  G_Fermi      = G_FERMI_DEFAULT;
  alphaS_MZ_MS = ALPHA_S_DEFAULT;
  mb_mb_MS     = MBMB_DEFAULT;
  M_Z_phys     = M_Z_DEFAULT;
  M_W_phys     = M_W_DEFAULT;
  M_top_phys   = M_TOP_DEFAULT;
  M_bot_phys   = M_BOT_DEFAULT; /* Pole mass! Not mb(mb). */
  M_tau_phys   = M_TAU_DEFAULT;
  M2_Z_phys = M_Z_DEFAULT*M_Z_DEFAULT;
  M2_W_phys = M_W_DEFAULT*M_W_DEFAULT;
  M2_top_phys = M_TOP_DEFAULT*M_TOP_DEFAULT;
  M2_bot_phys = M_BOT_DEFAULT*M_BOT_DEFAULT;
  M2_tau_phys = M_TAU_DEFAULT*M_TAU_DEFAULT;

  return 0;
}

/* ------------------------------------------------------------------ */ 
/* ------------------------------------------------------------------ */ 
/* This function first does a phase rotation on mu and b to make b
   real and positive. Then it computes a bunch of useful stuff. This
   function should always be called before computing the effective
   potential!  When in doubt, call it anyway, it's fast and doesn't
   hurt anything. */

int SUMO_Update ()
{
  TSIL_COMPLEX bphase;

  if (TSIL_CABS(b) > TSIL_TOL)
    bphase = (b)/TSIL_CABS(b);
  else
    bphase = 1.0;

  b = SUMO_CABS(b);
  mu = mu/bphase;

  Q2 = Q*Q;

  g2 = g*g;
  gp2 = gp*gp;
  g2plusgp2 = g2 + gp2;

  if (g2plusgp2 < TSIL_TOL) e2 = 0.0;
  else e2 = g2*gp2/g2plusgp2;

  e = SUMO_SQRT(e2);

  vu2 = vu*vu;
  vd2 = vd*vd;
  v2 = vu2 + vd2;
  vuvd = vu*vd;

  if (vd > TSIL_TOL) tanbeta = vu/vd; else tanbeta = 666;

  if (v2 > TSIL_TOL) {
    sinbeta = vu/TSIL_SQRT(v2);
    cosbeta = vd/TSIL_SQRT(v2);
  } else{
    sinbeta = 0;
    cosbeta = 0;
  }

  ytop2 = ytop*ytop;
  ybot2 = ybot*ybot;
  ytau2 = ytau*ytau;

  muc = SUMO_CONJ(mu);
  mu2 = TSIL_CREAL(mu*muc);
  atopc = SUMO_CONJ(atop);
  abotc = SUMO_CONJ(abot);
  atauc = SUMO_CONJ(atau);

  is_updated = 1;
  return 0;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Gives the signed square root of a possibly negative real number    */

TSIL_REAL SUMO_SGNSQRT (TSIL_REAL x)
{
  if (x < 0.0L) 
    return (-TSIL_SQRT(-x));
  else 
    return (TSIL_SQRT(x));
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

char *SUMO_Name (void) 
{
  return SUMO_NAME;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

char *SUMO_Version (void) 
{
  return SUMO_VERSION;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

/* void SUMO_Note (char *func, char *msg) */
/* { */
/*   fprintf (stderr, "NOTE (%s): %s\n", func, msg); */
/* } */

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SUMO_Warn (char *func, char *msg)
{
  TSIL_Warn (func, msg);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
#include <execinfo.h>

void SUMO_Error (char *func, char *msg, int code)
{
  void *callstack[128];
  int i, frames = backtrace(callstack, 128);
  char **strs = backtrace_symbols(callstack, frames);

  fprintf(stdout, "%s: %s Exiting...\n", func, msg);
  fprintf(stderr, "%s: %s Exiting...\n", func, msg);
  for (i=0; i<frames; i++)
    fprintf(stdout, "%s\n", strs[i]);

  free(strs);
  exit(code);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

TSIL_REAL SUMO_AbsSq (TSIL_COMPLEX z)
{
  return (TSIL_REAL) (z * TSIL_CONJ(z));
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Returns TRUE (1) is the numbers are equal and FALSE (0) if they are
   not equal.  */

int SUMO_FPCompare (TSIL_REAL x, TSIL_REAL y)
{
  TSIL_REAL tmp;
  TSIL_REAL absx, absy;

  absx = TSIL_FABS(x); absy = TSIL_FABS(y);

  /* First check for 0 = 0? */
  if (absx < 1000*TSIL_TOL) {
    if (absy < 1000*TSIL_TOL) return TRUE;
    else return FALSE;
  }

  /* Make x the one with the larger abs value: */
  if (absx < absy) {
    tmp = y; y = x; x = tmp;
  }
  if (TSIL_FABS(x-y) < absx*TSIL_TOL) return TRUE;
  else return FALSE;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Sets all pole masses equal to the corresponding tree masses. Should
   only be called once, just to have seed values. Note that M_W_phys,
   M_Z_phys, M_top_phys, M_bot_phys, and M_tau_phys and the
   corresponding squared masses are Standard Model inputs, so are
   already set either by default or read in from the SLHA input file.
   This was formerly SUMO_Tree_to_Pole.
*/

int SUMO_Pole_masses_from_tree () 
{
  int i;

  M2_Z = m2_Z;
  M2_W = m2_W;

  M2_top = m2_top;
  M_top = m_top;

  M2_bot = m2_bot;
  M_bot = m_bot;

  M2_tau = m2_tau;
  M_tau = m_tau;

  M2_gluino = m2_gluino;
  M_gluino = m_gluino;
  
  for (i=0; i<4; i++) {
    M2_Neut[i] = m2_Neut[i];
    M_Neut[i] = m_Neut[i];
  }

  for (i=0; i<2; i++) {
    M2_stop[i] = m2_stop[i];
    M2_sbot[i] = m2_sbot[i];
    M2_stau[i] = m2_stau[i];
    M2_seL[i] = m2_seL[i];
    M2_seR[i] = m2_seR[i];
    M2_suL[i] = m2_suL[i];
    M2_suR[i] = m2_suR[i];
    M2_sdL[i] = m2_sdL[i];
    M2_sdR[i] = m2_sdR[i];
    M2_Char[i] = m2_Char[i];
    M_Char[i] = m_Char[i];
  }

  for (i=0; i<3; i++) {
    M2_snu[i] = m2_snu[i];
  }

  M2_phi0[0] = m2_phi0[0];
  M2_phi0[1] = m2_phi0[1];
  M2_phi0[2] = 0.0L; /* The true all-orders pole mass of G^0. */
  M2_phi0[3] = m2_phi0[3];

  M2_phip[0] = 0.0L; /* The true all-orders pole mass of G^\pm. */
  M2_phip[1] = m2_phip[1];

  return 0;
}

/* ------------------------------------------------------------------ */
/* Sets all gluino, neutralino, chargino, Higgs, and sfermion tree
   masses equal to the corresponding physical masses. Used for
   calculating self-energy functions with physical masses substituted
   for tree masses. The t,b,tau,Z,W masses come from physical inputs
   rather than calculation, and the G^0 = phi0[2] and G^+ = phip[1]
   masses are set to 0 (the Landau gauge value).  Not to be confused
   with the old function SUMO_Tree_to_Pole(), which has now been
   renamed SUMO_Pole_masses_from_tree(), above.
*/

int SUMO_Tree_masses_from_pole () 
{
  int i;

  m_top = M_top_phys;
  m_bot = M_bot_phys;
  m_tau = M_tau_phys;
  m2_top = M_top_phys * M_top_phys;
  m2_bot = M_bot_phys * M_bot_phys;
  m2_tau = M_tau_phys * M_tau_phys;
  m2_Z = M_Z_phys * M_Z_phys;
  m2_W = M_W_phys * M_W_phys;

  m2_gluino = M2_gluino;
  m_gluino = M_gluino;
  
  for (i=0; i<4; i++) {
    m2_Neut[i] = M2_Neut[i];
    m_Neut[i] = M_Neut[i];
  }

  for (i=0; i<2; i++) {
    m2_stop[i] = M2_stop[i];
    m2_sbot[i] = M2_sbot[i];
    m2_stau[i] = M2_stau[i];
    m2_seL[i] = M2_seL[i];
    m2_seR[i] = M2_seR[i];
    m2_suL[i] = M2_suL[i];
    m2_suR[i] = M2_suR[i];
    m2_sdL[i] = M2_sdL[i];
    m2_sdR[i] = M2_sdR[i];
    m2_Char[i] = M2_Char[i];
    m_Char[i] = M_Char[i];
  }

  for (i=0; i<3; i++) {
    m2_snu[i] = M2_snu[i];
  }

  m2_phi0[0] = M2_phi0[0];
  m2_phi0[1] = M2_phi0[1];
  m2_phi0[2] = 0.0L; /* The true all-orders pole mass of G^0. */
  m2_phi0[3] = M2_phi0[3];

  m2_phip[0] = 0.0L; /* The true all-orders pole mass of G^\pm. */
  m2_phip[1] = M2_phip[1];

  return 0;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* A'(x), not in TSIL API.                                            */

TSIL_COMPLEX SUMO_Ap (TSIL_REAL x, TSIL_REAL qq) 
{
  if (TSIL_FABS(x) < TSIL_TOL) return 0.0;
  if (x > 0) return (TSIL_LOG(x/qq));  
  return (TSIL_LOG(-x/qq) + I*PI);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* For experimentation at this stage: sets parameter values to match
   those in hep-ph/0206136 section 6 */

int SUMO_SetTestModel (void) 
{
  Q = 640.L;

  /* For testing, these are the 1-loop values including strong
     contribs only: */
/*   vu = 172.610120L; */
/*   vd = 17.415145L; */

  /* For testing, these are the 2-loop values including strong
     contribs only: */
/*   vu = 177.040132L; */
/*   vd = 17.765925L; */

  vu = 172.L;
  vd = 17.2L;
  gp = 0.36L;
  g = 0.65L;
  g3 = 1.06L;
  ytop = 0.90L;
  ybot = 0.13L;
  ytau = 0.10L;
  M1 = 150.L;
  M2 = 280.L;
  M3 = 800.L;
  atop = -600.L;
  abot = -150.L;
  atau = -40.L;
  m2Q[1] = m2Q[0] = 780.L*780.L;
  m2u[1] = m2u[0] = 740.L*740.L;
  m2d[1] = m2d[0] = 735.L*735.L;
  m2L[1] = m2L[0] = 280.L*280.L;
  m2e[1] = m2e[0] = 200.L*200.L;
  m2Q[2] = 700.L*700.L;
  m2u[2] = 580.L*580.L;
  m2d[2] = 725.L*725.L;
  m2L[2] = 270.L*270.L;
  m2e[2] = 195.L*195.L;
  m2Hu = -500.L*500.L;
  m2Hd = 270.L*270.L;
  mu = 504.1811202L;
  b = 33937.10367865087L;
  LambdaVacuum = 0.L;

  tanbeta = vu/vd;

  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  /* Implement status code eventually... */
  return 0;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* For experimentation at this stage. Parameters from susana input
   file:

   600_2500_-5000_Softsusy.slha 

*/

int SUMO_SetTestModel2 (void) 
{
  Q = 1192.9812400;
  gp = 0.3618556; g = 0.6377560; g3 = 1.0443837;
  ytop = 0.8405106; ybot = 0.3727013; ytau = 0.3042641;
  vu = 172.0511758; vd = 5.8975940;
  tanbeta = 29.1731131;
  mu = 1997.7883300; b = 274833.6174415;
  m2Hu = -4000758.9500000; m2Hd = 2282924.3500000;
  M1 = 260.3053070; M2 = 481.0945230; M3 = 1326.7861500;
  atop = -2324.4248712; abot = -1913.4992102; atau = -1348.9331542;
  m2Q[2] = pow(1778.4864000,2); m2u[2] = pow(818.4605890,2);
  m2d[2] = pow(2325.6185000,2); 
  m2L[2] = pow(2293.6290200,2); m2e[2] = pow(2020.6285100,2);
  m2Q[1] = pow(2720.9045100,2); m2u[1] = pow(2708.9807700,2);
  m2d[1] = pow(2707.6696400,2);
  m2L[1] = pow(2518.9490200,2); m2e[1] = pow(2503.5243300,2);
  m2Q[0] = pow(2721.0360800,2); m2u[0] = pow(2709.0064100,2);
  m2d[0] = pow(2707.9107700,2);
  m2L[0] = pow(2519.6904000,2); m2e[0] = pow(2505.0246300,2);

  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  return 0;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* For experimentation at this stage. Parameters from susana input
   file:

   600_2500_-5000_SuSpect.slha 

*/

int SUMO_SetTestModel3 (void) 
{
  Q = 1185.2958500;
  gp = 0.3618578; g = 0.6372790; g3 = 1.0407687;
  ytop = 0.8418495; ybot = 0.3855710; ytau = 0.3024594;
  vu = 172.3771469; vd = 5.9095630;
  tanbeta = 29.1691869;
  mu = 1991.0474400; b = 269580.4222518;
  m2Hu = -3973340.2400000; m2Hd = 2119124.6100000;
  M1 = 262.5403450; M2 = 482.4322400; M3 = 1318.0386700;
  atop = -2296.7082491; abot = -1957.5889206; atau = -1334.6813211;
  m2Q[2] = pow(1758.2342000,2); m2u[2] = pow(817.2576430,2);
  m2d[2] = pow(2295.1127400,2);
  m2L[2] = pow(2297.6521100,2); m2e[2] = pow(2029.3275800,2);
  m2Q[1] = pow(2717.5845000,2); m2u[1] = pow(2705.5160100,2);
  m2d[1] = pow(2704.4123800,2);
  m2L[1] = pow(2519.6745100,2); m2e[1] = pow(2505.0860100,2);
  m2Q[0] = pow(2717.5845000,2); m2u[0] = pow(2705.5160100,2);
  m2d[0] = pow(2704.4123800,2);
  m2L[0] = pow(2519.6745100,2); m2e[0] = pow(2505.0860100,2);

  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  return 0;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* For experimentation at this stage. Parameters from Softsusy program
   main.cpp, via printout in higherorder.c.

   MSUGRA model with:
   m12 = 500, a0 = 0, tanbeta = 10, m0 = 125
*/

int SUMO_SetTestModel4 (void) 
{
  gp = 3.625155e-01;
  g = 6.430183e-01;
  g3 = 1.065820e+00;
  ytop = 8.354106e-01;
  ybot = 1.315183e-01;
  ytau = 1.006065e-01;
  M1 = 2.092662e+02;
  M2 = 3.881985e+02;
  M3 = 1.117739e+03;
  atop = -7.680666e+02;
  abot = -1.734235e+02;
  atau = -2.996008e+01;
  m2Hu = -3.659024e+05;
  m2Hd = 1.096987e+05;
  mu = 6.145860e+02;
  b = 5.268885e+04;
  Q = 8.844737e+02;
  m2Q[0] = 1.046745e+06;
  m2u[0] = 9.710204e+05;
  m2d[0] = 9.620489e+05;
  m2L[0] = 1.253124e+05;
  m2e[0] = 4.925768e+04;
  m2Q[1] = 1.046740e+06;
  m2u[1] = 9.710151e+05;
  m2d[1] = 9.620436e+05;
  m2L[1] = 1.253100e+05;
  m2e[1] = 4.925288e+04;
  m2Q[2] = 8.960003e+05;
  m2u[2] = 6.720852e+05;
  m2d[2] = 9.527687e+05;
  m2L[2] = 1.246002e+05;
  m2e[2] = 4.780542e+04;

  /* These are the Softsusy values... */
  vu = 1.718402e+02;
  vd = 1.773515e+01;

  SUMO_Update ();
  SUMO_Tree_Masses ();
  SUMO_Tree_Couplings ();

  /* SUMO_Minimize_Veff (2); */

  return 0;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Works the same as TSIL_GetFunction if interp == NO. Threshold
   interpolation performed for V functions when interp == YES. */

TSIL_COMPLEX SUMO_GetFunction (TSIL_DATA *foo, const char *which, 
			       int interp)
{
  TSIL_REAL arg1, arg2, delta, snew;
  TSIL_COMPLEX Vplus, Vminus;
  TSIL_DATA gaak;

  /* This is cut and pasted from tsil_names.h: */
  const char *vname[4][2] = {{"Vzxyv","Vzxvy"},
			     {"Vuyxv","Vuyvx"},
			     {"Vxzuv","Vxzvu"},
			     {"Vyuzv","Vyuvz"}};

  /* If no interp requested, or not a V function, just return the
     usual thing: */
  if (interp == NO || strncmp (which, "V", 1) != 0)
    return TSIL_GetFunction (foo, which);

  /* Check for a threshold case: */
  if (   !strcmp(which, vname[0][0]) || !strcmp(which, vname[0][1])
      || !strcmp(which, vname[2][0]) || !strcmp(which, vname[2][1])) {
    arg1 = foo->z; arg2 = foo->x;
  }
  else if (   !strcmp(which, vname[1][0]) || !strcmp(which, vname[1][1])
           || !strcmp(which, vname[3][0]) || !strcmp(which, vname[3][1])) {
    arg1 = foo->u; arg2 = foo->y;
  }
  else {
    printf("This can never happen!!!\n"); exit(234);
  }

  delta = foo->s/TSIL_POW(TSIL_SQRT(arg1)+TSIL_SQRT(arg2),2) - 1.0L;

  if (TSIL_FABS(delta) > THRESH_TOL)
    return TSIL_GetFunction (foo, which);

  /* If we get here we interpolate: */
  TSIL_SetParameters (&gaak, foo->x, foo->y, foo->z, foo->u, foo->v, foo->qq);
  snew = (1.0L + THRESH_TOL)*(foo->s);

  TSIL_Evaluate (&gaak, snew);
  Vplus = TSIL_GetFunction (&gaak, which);

  snew = (1.0L - THRESH_TOL)*(foo->s);
  TSIL_Evaluate (&gaak, snew);
  Vminus = TSIL_GetFunction (&gaak, which);

  return 0.5L*(1.0L + delta/THRESH_TOL)*Vplus + 
         0.5L*(1.0L - delta/THRESH_TOL)*Vminus;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Works the same as TSIL_GetFunctionR if interp == NO. Threshold
   interpolation performed for V functions when interp == YES. */

/* Takes a TSIL_RESULT * as first arg, rather than TSIL_DATA *. */

TSIL_COMPLEX SUMO_GetFunctionR (TSIL_RESULT *foo, const char *which, 
				int interp)
{
  TSIL_REAL arg1, arg2, delta, snew;
  TSIL_COMPLEX Vplus, Vminus;
  TSIL_DATA gaak;

  /* This is cut and pasted from tsil_names.h: */
  const char *vname[4][2] = {{"Vzxyv","Vzxvy"},
			     {"Vuyxv","Vuyvx"},
			     {"Vxzuv","Vxzvu"},
			     {"Vyuzv","Vyuvz"}};

  /* If no interp requested, or not a V function, just return the
     usual thing: */
  if (interp == NO || strncmp (which, "V", 1) != 0)
    return TSIL_GetFunctionR (foo, which);

  /* Check for a threshold case: */
  if (   !strcmp(which, vname[0][0]) || !strcmp(which, vname[0][1])
      || !strcmp(which, vname[2][0]) || !strcmp(which, vname[2][1])) {
    arg1 = foo->z; arg2 = foo->x;
  }
  else if (   !strcmp(which, vname[1][0]) || !strcmp(which, vname[1][1])
           || !strcmp(which, vname[3][0]) || !strcmp(which, vname[3][1])) {
    arg1 = foo->u; arg2 = foo->y;
  }
  else {
    printf("This can never happen!!!\n"); exit(234);
  }

  delta = foo->s/TSIL_POW(TSIL_SQRT(arg1)+TSIL_SQRT(arg2),2) - 1.0L;

  if (TSIL_FABS(delta) > THRESH_TOL)
    return TSIL_GetFunctionR (foo, which);

  /* If we get here we interpolate: */
  TSIL_SetParameters (&gaak, foo->x, foo->y, foo->z, foo->u, foo->v, foo->qq);
  snew = (1.0L + THRESH_TOL)*(foo->s);
  TSIL_Evaluate (&gaak, snew);
  Vplus = TSIL_GetFunction (&gaak, which);

  snew = (1.0L - THRESH_TOL)*(foo->s);
  TSIL_Evaluate (&gaak, snew);
  Vminus = TSIL_GetFunction (&gaak, which);

  return 0.5L*(1.0L + delta/THRESH_TOL)*Vplus + 
         0.5L*(1.0L - delta/THRESH_TOL)*Vminus;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

TSIL_COMPLEX SUMO_Bp (TSIL_REAL x,
		      TSIL_REAL y,
		      TSIL_COMPLEX s,
		      TSIL_REAL qq,
		      int interp)
{
  TSIL_COMPLEX snew, Bpplus, Bpminus;
  TSIL_REAL delta;

  if (interp == NO) return TSIL_Bp (x,y,s,qq);

  delta = TSIL_CABS(s)/TSIL_POW(TSIL_SQRT(x)+TSIL_SQRT(y), 2) - 1.0L;

  if (TSIL_FABS(delta) > THRESH_TOL)
    return TSIL_Bp (x,y,s,qq);

  /* If we get here we interpolate: */
  snew = (1.0L + THRESH_TOL)*s;
  Bpplus = TSIL_Bp (x,y,snew,qq);

  snew = (1.0L - THRESH_TOL)*s;
  Bpminus = TSIL_Bp (x,y,snew,qq);

  return 0.5L*(1.0L + delta/THRESH_TOL)*Bpplus + 
         0.5L*(1.0L - delta/THRESH_TOL)*Bpminus;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

TSIL_COMPLEX SUMO_dBds (TSIL_REAL x,
			TSIL_REAL y,
			TSIL_COMPLEX s,
			TSIL_REAL qq,
			int interp)
{
  TSIL_COMPLEX snew, dBdsplus, dBdsminus;
  TSIL_REAL delta;

  if (interp == NO) return TSIL_dBds (x,y,s,qq);

  delta = TSIL_CABS(s)/TSIL_POW(TSIL_SQRT(x)+TSIL_SQRT(y), 2) - 1.0L;

  if (TSIL_FABS(delta) > THRESH_TOL)
    return TSIL_dBds (x,y,s,qq);

  /* If we get here we interpolate: */
  snew = (1.0L + THRESH_TOL)*s;
  dBdsplus = TSIL_dBds (x,y,snew,qq);

  snew = (1.0L - THRESH_TOL)*s;
  dBdsminus = TSIL_dBds (x,y,snew,qq);

  return 0.5L*(1.0L + delta/THRESH_TOL)*dBdsplus + 
         0.5L*(1.0L - delta/THRESH_TOL)*dBdsminus;
}
