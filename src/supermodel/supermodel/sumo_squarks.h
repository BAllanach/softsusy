/* Self energy functions for squarks. */

#ifndef SUMO_SQUARKS_H
#define SUMO_SQUARKS_H

#include "../tsil/tsil.h"

/* Functions defined in se1_squark.c: */
TSIL_COMPLEX pi1_stop (int, int, TSIL_REAL);
TSIL_COMPLEX pi1_sbot (int, int, TSIL_REAL);
TSIL_COMPLEX pi1_sdL (int, TSIL_REAL); 
TSIL_COMPLEX pi1_sdR (int, TSIL_REAL); 
TSIL_COMPLEX pi1_suL (int, TSIL_REAL); 
TSIL_COMPLEX pi1_suR (int, TSIL_REAL);
TSIL_COMPLEX Pi1_stop (int, int, TSIL_REAL);
TSIL_COMPLEX Pi1_sbot (int, int, TSIL_REAL);
TSIL_COMPLEX Pi1_sdL (int, TSIL_REAL); 
TSIL_COMPLEX Pi1_sdR (int, TSIL_REAL); 
TSIL_COMPLEX Pi1_suL (int, TSIL_REAL); 
TSIL_COMPLEX Pi1_suR (int, TSIL_REAL);

/* Functions defined in se2_stop.c: */
TSIL_COMPLEX pi20tilde_stop (int, TSIL_REAL);
TSIL_COMPLEX pi21tilde_stop (int, TSIL_REAL);
TSIL_COMPLEX pi22tilde_stop (int, TSIL_REAL);
/* TSIL_COMPLEX Pi20tilde_stop (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi21tilde_stop (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi22tilde_stop (int, TSIL_REAL); */

/* Functions defined in se2_sbot.c: */
TSIL_COMPLEX pi20tilde_sbot (int, TSIL_REAL);
TSIL_COMPLEX pi21tilde_sbot (int, TSIL_REAL);
TSIL_COMPLEX pi22tilde_sbot (int, TSIL_REAL);
/* TSIL_COMPLEX Pi20tilde_sbot (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi21tilde_sbot (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi22tilde_sbot (int, TSIL_REAL); */

/* Self energy functions defined in se2_suL.c: */
TSIL_COMPLEX pi20tilde_suL (int, TSIL_REAL);
TSIL_COMPLEX pi21tilde_suL (int, TSIL_REAL);
TSIL_COMPLEX pi22tilde_suL (int, TSIL_REAL);
/* TSIL_COMPLEX Pi20tilde_suL (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi21tilde_suL (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi22tilde_suL (int, TSIL_REAL); */

/* Self energy functions defined in se2_suR.c: */
TSIL_COMPLEX pi20tilde_suR (int, TSIL_REAL);
TSIL_COMPLEX pi21tilde_suR (int, TSIL_REAL);
TSIL_COMPLEX pi22tilde_suR (int, TSIL_REAL);
/* TSIL_COMPLEX Pi20tilde_suR (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi21tilde_suR (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi22tilde_suR (int, TSIL_REAL); */

/* Self energy functions defined in se2_sdL.c: */
TSIL_COMPLEX pi20tilde_sdL (int, TSIL_REAL);
TSIL_COMPLEX pi21tilde_sdL (int, TSIL_REAL);
TSIL_COMPLEX pi22tilde_sdL (int, TSIL_REAL);
/* TSIL_COMPLEX Pi20tilde_sdL (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi21tilde_sdL (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi22tilde_sdL (int, TSIL_REAL); */

/* Self energy functions defined in se2_sdR.c: */
TSIL_COMPLEX pi20tilde_sdR (int, TSIL_REAL);
TSIL_COMPLEX pi21tilde_sdR (int, TSIL_REAL);
TSIL_COMPLEX pi22tilde_sdR (int, TSIL_REAL);
/* TSIL_COMPLEX Pi20tilde_sdR (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi21tilde_sdR (int, TSIL_REAL); */
/* TSIL_COMPLEX Pi22tilde_sdR (int, TSIL_REAL); */

#endif
