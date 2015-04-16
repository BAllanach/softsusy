/* Miscellaneous functions */

#include "supermodel.h"

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Basic initialization of the model struct. This should be called
   first thing on any such struct, before doing any calculation. In
   normal operation, it is called by SUMO_Read_SLHA. */

#include "sumo_smdata.h"

/* int SUMO_Initialize ()  */
/* { */
  /* The most generic values here: */
  /* particleContent = MSSM; */
  /* modelType       = GENERAL; */

  /* is_updated = 0; */
  /* are_tree_masses_updated = 0; */
  /* are_tree_couplings_updated = 0; */

  /* Set the SM data from defaults. The input SLHA file will typically
     override most or all of these except M_bot_phys and
     M2_bot_phys. */
/*   alpha_MZ_MS  = 1.0L/ALPHA_MZ_INV_DEFAULT; */
/*   alpha0       = 1.0L/ALPHA_0_INV_DEFAULT; */
/*   G_Fermi      = G_FERMI_DEFAULT; */
/*   alphaS_MZ_MS = ALPHA_S_DEFAULT; */
/*   mb_mb_MS     = MBMB_DEFAULT; */
/*   M_Z_phys     = M_Z_DEFAULT; */
/*   M_W_phys     = M_W_DEFAULT; */
/*   M_top_phys   = M_TOP_DEFAULT; */
/*   M_bot_phys   = M_BOT_DEFAULT; /\* Pole mass! Not mb(mb). *\/ */
/*   M_tau_phys   = M_TAU_DEFAULT; */
/*   M2_Z_phys = M_Z_DEFAULT*M_Z_DEFAULT; */
/*   M2_W_phys = M_W_DEFAULT*M_W_DEFAULT; */
/*   M2_top_phys = M_TOP_DEFAULT*M_TOP_DEFAULT; */
/*   M2_bot_phys = M_BOT_DEFAULT*M_BOT_DEFAULT; */
/*   M2_tau_phys = M_TAU_DEFAULT*M_TAU_DEFAULT; */

/*   return 0; */
/* } */

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

/* int SUMO_Delta (int i, int j) */
/* { */
/*   if (i == j) return 1; */
/*   else return 0; */
/* } */

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

TSIL_REAL SUMO_AbsSq (TSIL_COMPLEX z)
{
  return (TSIL_REAL) (z * TSIL_CONJ(z));
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
/* Sets all gluino, neutralino, chargino, Higgs, and sfermion 
   tree masses equal to the corresponding physical masses. Used for 
   calculating self-energy functions with physical masses substituted
   for tree masses. The t,b,tau,Z,W masses come from physical inputs
   rather than calculation, and the G^0 = phi0[2] and G^+ = phip[1] 
   masses are set to 0 (the Landau gauge value).
   Not to be confused with the old function SUMO_Tree_to_Pole(), which
   has now been renamed SUMO_Pole_masses_from_tree(), above.
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
