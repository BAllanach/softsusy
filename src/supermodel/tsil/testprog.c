/* Test Suite Driver Program */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

long double strtold (const char *, char **);

#include "internal.h"
#include "tsil_testparams.h"

#define Doperms 0
/* Set to 1 (0) if you do (don't) want permutations of the masses not
   involving v to be tested. */

#define Doextraperms 0
/* Set to 1 (0) if you do (don't) want permutations involving v to be
   tested. */

FILE *fp;

enum { FAIL, WARN, PASS };

/* **************************************************************** */

void TSIL_SkipLine (FILE * foo)
{
  int s;

  while ((s = getc (foo)) != '\n')
    ;

  return;
}

/* **************************************************************** */

void TSIL_Compare (const char   *name,
		   TSIL_COMPLEX actual,
		   TSIL_COMPLEX computed,
		   TSIL_REAL    allow_pass,
		   TSIL_REAL    allow_warn,
		   int          *result)
{
  TSIL_REAL a_re, a_im, c_re, c_im, magnitude, err;
  int foo;

  a_re = TSIL_CREAL (actual);
  a_im = TSIL_CIMAG (actual);
  c_re = TSIL_CREAL (computed);
  c_im = TSIL_CIMAG (computed);
  magnitude = TSIL_CABS (actual) + TSIL_TOL;

  /* DGR */
  if (TSIL_IsInfinite (actual))
    {
      if (TSIL_IsInfinite (computed))
	foo = PASS * PASS;
      else
	foo = FAIL;
    }
  else
    {
      /* Check Real part */
      err = TSIL_FABS (a_re - c_re) / magnitude;

      if (err < allow_pass)
	foo = PASS;
      else if (err < allow_warn)
	foo = WARN;
      else {
/* 	printf("\nFailure in re part: err = %Le\n", (long double) err); */
	foo = FAIL;
      }
    
      /* Check Imaginary part */
      err = TSIL_FABS (a_im - c_im) / magnitude;

      if (err < allow_pass)
	foo *= PASS;
      else if (err < allow_warn)
	foo *= WARN;
      else {
/* 	printf("\nFailure in im part: err = %Le\n", (long double) err); */
	foo *= FAIL;
      }
    }

  if (foo == 4)
    *result = PASS;
  else if (foo == 1 || foo == 2)
    {
      *result = WARN;
      printf ("\nWARN\n");
      printf ("Expected for %s: ", name);
      TSIL_cprintfM (actual);
      printf ("\n");
      printf ("Obtained for %s: ", name);
      TSIL_cprintfM (computed);
      printf ("\n");
    }
  else if (foo == 0)
    {
      *result = FAIL;
      printf ("\nFAIL\n");
      printf ("Expected for %s: ", name);
      TSIL_cprintfM (actual);
      printf ("\n");
      printf ("Obtained for %s: ", name);
      TSIL_cprintfM (computed);
      printf ("\n");
    }
  else
    printf ("NOPE! Can't EVER get here in TSIL_Compare!!!\n");

  return;
}

/* **************************************************************** */

TSIL_REAL TSIL_GetReal ()
{
  TSIL_REAL val;

  while (fgetc (fp) != '=')
    ;

#if defined(TSIL_SIZE_DOUBLE)
  fscanf (fp, "%lf;", &val);
#else
  fscanf (fp, "%Lf;", &val);
#endif
  return val;
}

/* **************************************************************** */

TSIL_COMPLEX TSIL_GetComplex ()
{
  TSIL_COMPLEX val;
  TSIL_REAL    im;
  char         s[50];
  const char   inf[] = "ComplexInfinity;";

  while (fgetc (fp) != '=')
    ;

  /* This automatically skips leading white space: */
  fscanf (fp, "%s", s);

  if (strcmp (s, inf) == 0)
    val = TSIL_Infinity;
  else
    {
      val = (TSIL_COMPLEX) strtold (s, (char **) NULL);
#if defined(TSIL_SIZE_DOUBLE)
      fscanf (fp, " + %lf I;", &im);
#else
      fscanf (fp, " + %Lf I;", &im);
#endif
      val += I * im;
    }

  return val;
}

/* ******************************************************************* */

void TSIL_SwapR (TSIL_REAL * px, TSIL_REAL * py)
{
  TSIL_REAL temp;

  temp = *px; *px = *py; *py = temp;
}

/* ******************************************************************* */

void TSIL_SwapC (TSIL_COMPLEX * px, TSIL_COMPLEX * py)
{
  TSIL_COMPLEX temp;

  temp = *px; *px = *py; *py = temp;
}

/* ******************************************************************* */

void TSIL_PermuteResults (TSIL_DATA * foo, int p)
{
  int i;

  if (1 == p)
    {
      TSIL_SwapR (&(foo->x), &(foo->y));
      TSIL_SwapR (&(foo->z), &(foo->u));
      TSIL_SwapC (&(foo->B[xz].value), &(foo->B[yu].value));
      TSIL_SwapC (&(foo->S[vyz].value), &(foo->S[uxv].value));
      TSIL_SwapC (&(foo->T[vyz].value), &(foo->T[vxu].value));
      TSIL_SwapC (&(foo->T[uxv].value), &(foo->T[zyv].value));
      TSIL_SwapC (&(foo->T[xuv].value), &(foo->T[yzv].value));
      TSIL_SwapC (&(foo->Tbar[vyz].value), &(foo->Tbar[vxu].value));
      TSIL_SwapC (&(foo->Tbar[uxv].value), &(foo->Tbar[zyv].value));
      TSIL_SwapC (&(foo->Tbar[xuv].value), &(foo->Tbar[yzv].value));
      TSIL_SwapC (&(foo->U[zxyv].value), &(foo->U[uyxv].value));
      TSIL_SwapC (&(foo->U[xzuv].value), &(foo->U[yuzv].value));
      TSIL_SwapC (&(foo->V[zxyv].value), &(foo->V[uyxv].value));
      TSIL_SwapC (&(foo->V[xzuv].value), &(foo->V[yuzv].value));

      for (i = 0; i < 3; i++)
	{
	  TSIL_SwapC (&(foo->S[vyz].bold[i]), &(foo->S[uxv].bold[i]));
	  TSIL_SwapC (&(foo->T[vyz].bold[i]), &(foo->T[vxu].bold[i]));
	  TSIL_SwapC (&(foo->T[uxv].bold[i]), &(foo->T[zyv].bold[i]));
	  TSIL_SwapC (&(foo->T[xuv].bold[i]), &(foo->T[yzv].bold[i]));
	  TSIL_SwapC (&(foo->U[zxyv].bold[i]), &(foo->U[uyxv].bold[i]));
	  TSIL_SwapC (&(foo->U[xzuv].bold[i]), &(foo->U[yuzv].bold[i]));
	  TSIL_SwapC (&(foo->V[zxyv].bold[i]), &(foo->V[uyxv].bold[i]));
	  TSIL_SwapC (&(foo->V[xzuv].bold[i]), &(foo->V[yuzv].bold[i]));
	}
    }

  if (2 == p)
    {
      TSIL_SwapR (&(foo->x), &(foo->z));
      TSIL_SwapR (&(foo->y), &(foo->u));
      TSIL_SwapC (&(foo->S[vyz].value), &(foo->S[uxv].value));
      TSIL_SwapC (&(foo->T[vyz].value), &(foo->T[vxu].value));
      TSIL_SwapC (&(foo->T[uxv].value), &(foo->T[yzv].value));
      TSIL_SwapC (&(foo->T[xuv].value), &(foo->T[zyv].value));
      TSIL_SwapC (&(foo->Tbar[vyz].value), &(foo->Tbar[vxu].value));
      TSIL_SwapC (&(foo->Tbar[uxv].value), &(foo->Tbar[yzv].value));
      TSIL_SwapC (&(foo->Tbar[xuv].value), &(foo->Tbar[zyv].value));
      TSIL_SwapC (&(foo->U[zxyv].value), &(foo->U[xzuv].value));
      TSIL_SwapC (&(foo->U[uyxv].value), &(foo->U[yuzv].value));
      TSIL_SwapC (&(foo->V[zxyv].value), &(foo->V[xzuv].value));
      TSIL_SwapC (&(foo->V[uyxv].value), &(foo->V[yuzv].value));

      for (i = 0; i < 3; i++)
	{
	  TSIL_SwapC (&(foo->S[vyz].bold[i]), &(foo->S[uxv].bold[i]));
	  TSIL_SwapC (&(foo->T[vyz].bold[i]), &(foo->T[vxu].bold[i]));
	  TSIL_SwapC (&(foo->T[uxv].bold[i]), &(foo->T[yzv].bold[i]));
	  TSIL_SwapC (&(foo->T[xuv].bold[i]), &(foo->T[zyv].bold[i]));
	  TSIL_SwapC (&(foo->U[zxyv].bold[i]), &(foo->U[xzuv].bold[i]));
	  TSIL_SwapC (&(foo->U[uyxv].bold[i]), &(foo->U[yuzv].bold[i]));
	  TSIL_SwapC (&(foo->V[zxyv].bold[i]), &(foo->V[xzuv].bold[i]));
	  TSIL_SwapC (&(foo->V[uyxv].bold[i]), &(foo->V[yuzv].bold[i]));
	}
    }

  if (3 == p)
    {
      TSIL_SwapR (&(foo->x), &(foo->u));
      TSIL_SwapR (&(foo->y), &(foo->z));
      TSIL_SwapC (&(foo->B[xz].value), &(foo->B[yu].value));
      TSIL_SwapC (&(foo->T[uxv].value), &(foo->T[xuv].value));
      TSIL_SwapC (&(foo->T[zyv].value), &(foo->T[yzv].value));
      TSIL_SwapC (&(foo->Tbar[uxv].value), &(foo->Tbar[xuv].value));
      TSIL_SwapC (&(foo->Tbar[zyv].value), &(foo->Tbar[yzv].value));
      TSIL_SwapC (&(foo->U[zxyv].value), &(foo->U[yuzv].value));
      TSIL_SwapC (&(foo->U[xzuv].value), &(foo->U[uyxv].value));
      TSIL_SwapC (&(foo->V[zxyv].value), &(foo->V[yuzv].value));
      TSIL_SwapC (&(foo->V[xzuv].value), &(foo->V[uyxv].value));
      for (i = 0; i < 3; i++)
	{
	  TSIL_SwapC (&(foo->T[uxv].bold[i]), &(foo->T[xuv].bold[i]));
	  TSIL_SwapC (&(foo->T[zyv].bold[i]), &(foo->T[yzv].bold[i]));
	  TSIL_SwapC (&(foo->U[zxyv].bold[i]), &(foo->U[yuzv].bold[i]));
	  TSIL_SwapC (&(foo->U[xzuv].bold[i]), &(foo->U[uyxv].bold[i]));
	  TSIL_SwapC (&(foo->V[zxyv].bold[i]), &(foo->V[yuzv].bold[i]));
	  TSIL_SwapC (&(foo->V[xzuv].bold[i]), &(foo->V[uyxv].bold[i]));
	}
    }

  return;
}

/* ******************************************************************* */

int main (int argc, char *argv[])
{
#include "tsil_names.h"

  TSIL_REAL x, y, z, u, v, s, qq;
  TSIL_COMPLEX val;
  TSIL_DATA result;
  int i, j, k, p, foo, tally[3], nwarn, nfail, npass;
  char c;
  TSIL_COMPLEX Mvalue, Uvalue[4], Vvalue[4], Svalue[2], Tvalue[6], Bvalue[2];
  TSIL_COMPLEX Tbarvalue[6], UUvalue[4][3], VVvalue[4][3], SSvalue[2][3];
  TSIL_COMPLEX TTvalue[6][3];
  TSIL_REAL xx, yy, zz, uu, vv;

  /* gcc complains unless the following are initialized; this may be a
     gcc bug? Fortunately, it doesn't hurt to initialize them. */
  int permU = 0;
  int permS = 0;
  int permT11 = 0;
  int permT12 = 0;
  int permT21 = 0;
  int permT22 = 0;
  int permT31 = 0;
  int permT32 = 0;

  if (argc == 1)
    TSIL_Error ("main", "Must supply test data filename(s)...", 2);

  nwarn = nfail = npass = 0;
  
  printf("===== TSIL TEST SUITE =====\n");
#if defined (TSIL_TEST_STU)
  printf("** Testing STU Evaluation only!\n");
#elif defined (TSIL_TEST_ST)
  printf("** Testing ST Evaluation only!\n");
#endif

  /* We don't need no stinkin' TSIL_Warnings here. */
  fclose (stderr);

  /* Loop over input files */
  for (i = 1; i < argc; i++)
    {
      if ((fp = fopen(argv[i], "r")) == NULL) {
	TSIL_Warn ("Test program", "Invalid file name");
	continue;
      }

      printf ("\nTest %d: ", i);
      printf ("%s\n", argv[i]);
      fflush (stdout);

      /* Skip any lines starting with '(' (comments): */
      while ((c = fgetc (fp)) == '(')
	TSIL_SkipLine (fp);

      /* Put back the last character after we find a non-comment line */
      ungetc ((int) c, fp);

      /* Read in parameters */
      xx = x = TSIL_GetReal ();
      yy = y = TSIL_GetReal ();
      zz = z = TSIL_GetReal ();
      uu = u = TSIL_GetReal ();
      vv = v = TSIL_GetReal ();
      s  = TSIL_GetReal ();
      qq = TSIL_GetReal ();

      /* Calculate everything... */
#if defined(TSIL_TEST_STU)
      TSIL_SetParametersSTU (&result, x, z, u, v, qq);
#elif defined(TSIL_TEST_ST)
      TSIL_SetParametersST (&result, x, u, v, qq);
#else
      TSIL_SetParameters (&result, x, y, z, u, v, qq);
#endif
      TSIL_Evaluate (&result, s);

      for (j = 0; j < 3; j++)
	tally[j] = 0;


      /* ...and test results: */

#if defined(TSIL_TEST_STU)
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* U */
      for (j = 2; j < 4; j+=2)
	{
	  val = TSIL_GetComplex ();
	  Uvalue[j] = result.U[j].value;
	  TSIL_Compare (uname[j][0], val, Uvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Uvalue[2], TSIL_GetFunction(&result,"Uxzuv"));

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* T */
      for (j = 1; j < 6; j+=2)
	{
	  val = TSIL_GetComplex ();
	  Tvalue[j] = result.T[j].value;
	  TSIL_Compare (tname[j][0], val, Tvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	  val = TSIL_GetComplex ();
	}
      TSIL_CheckConsistent (Tvalue[1], TSIL_GetFunction(&result,"Tuxv"));
      TSIL_CheckConsistent (Tvalue[3], TSIL_GetFunction(&result,"Txuv"));
      TSIL_CheckConsistent (Tvalue[5], TSIL_GetFunction(&result,"Tvxu"));

      /* S */
      for (j = 1; j < 2; j+=2)
	{
	  val = TSIL_GetComplex ();
	  Svalue[j] = result.S[j].value;
	  TSIL_Compare (sname[j][0], val, Svalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Svalue[1], TSIL_GetFunction(&result,"Suxv"));

      /* B */
      for (j = 0; j < 2; j+=2)
	{
	  val = TSIL_GetComplex ();
	  Bvalue[j] = result.B[j].value;
	  TSIL_Compare (bname[j][0], val, Bvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Bvalue[0], TSIL_GetFunction(&result,"Bxz"));

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* V */
      for (j = 2; j < 4; j+=2)
	{
	  val = TSIL_GetComplex ();
	  Vvalue[j] = result.V[j].value;
	  TSIL_Compare (vname[j][0], val, Vvalue[j], TSIL_PASS_V, TSIL_WARN_V, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Vvalue[2], TSIL_GetFunction(&result,"Vxzuv"));
      val = TSIL_GetComplex ();

      /* Tbar */
      val = TSIL_GetComplex ();
      for (j = 1; j < 6; j+=2)
	{
	  val = TSIL_GetComplex ();
	  Tbarvalue[j] = result.Tbar[j].value;
	  TSIL_Compare (tbarname[j][0], val, Tbarvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	  val = TSIL_GetComplex ();
	}

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* UU */
      for (j = 2; j < 4; j+=2)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      UUvalue[j][k] = result.U[j].bold[k];
	      TSIL_Compare (uuname[j][k], val, UUvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
	      tally[foo]++;
	    }
	}

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* VV */
      for (j = 2; j < 4; j+=2)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      VVvalue[j][k] = result.V[j].bold[k];
	      TSIL_Compare (vvname[j][k], val, VVvalue[j][k], TSIL_PASS_V, TSIL_WARN_V, &foo);
	      tally[foo]++;
	    }
	}
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* TT */
      for (j = 1; j < 6; j+=2)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      TTvalue[j][k] = result.T[j].bold[k];
	      TSIL_Compare (ttname[j][k], val, TTvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
	      tally[foo]++;
	    }
	  val = TSIL_GetComplex ();
	  val = TSIL_GetComplex ();
	  val = TSIL_GetComplex ();
	}

      /* SS */
      for (j = 1; j < 2; j+=2)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      SSvalue[j][k] = result.S[j].bold[k];
	      TSIL_Compare (ssname[j][k], val, SSvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
	      tally[foo]++;
	    }
	}
#elif defined(TSIL_TEST_ST)
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* U */
      for (j = 2; j < 4; j+=2)
	{
	  val = TSIL_GetComplex ();
	}

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* T */
      for (j = 1; j < 6; j+=2)
	{
	  val = TSIL_GetComplex ();
	  Tvalue[j] = result.T[j].value;
	  TSIL_Compare (tname[j][0], val, Tvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	  val = TSIL_GetComplex ();
	}
      TSIL_CheckConsistent (Tvalue[1], TSIL_GetFunction(&result,"Tuxv"));
      TSIL_CheckConsistent (Tvalue[3], TSIL_GetFunction(&result,"Txuv"));
      TSIL_CheckConsistent (Tvalue[5], TSIL_GetFunction(&result,"Tvxu"));

      /* S */
      for (j = 1; j < 2; j+=2)
	{
	  val = TSIL_GetComplex ();
	  Svalue[j] = result.S[j].value;
	  TSIL_Compare (sname[j][0], val, Svalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Svalue[1], TSIL_GetFunction(&result,"Suxv"));

      /* B */
      for (j = 0; j < 2; j+=2)
	{
	  val = TSIL_GetComplex ();
	}

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* V */
      for (j = 2; j < 4; j+=2)
	{
	  val = TSIL_GetComplex ();
	}
      val = TSIL_GetComplex ();

      /* Tbar */
      val = TSIL_GetComplex ();
      for (j = 1; j < 6; j+=2)
	{
	  val = TSIL_GetComplex ();
	  Tbarvalue[j] = result.Tbar[j].value;
	  TSIL_Compare (tbarname[j][0], val, Tbarvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	  val = TSIL_GetComplex ();
	}

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* UU */
      for (j = 2; j < 4; j+=2)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	    }
	}

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* VV */
      for (j = 2; j < 4; j+=2)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	    }
	}
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();

      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      val = TSIL_GetComplex ();
      /* TT */
      for (j = 1; j < 6; j+=2)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      TTvalue[j][k] = result.T[j].bold[k];
	      TSIL_Compare (ttname[j][k], val, TTvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
	      tally[foo]++;
	    }
	  val = TSIL_GetComplex ();
	  val = TSIL_GetComplex ();
	  val = TSIL_GetComplex ();
	}

      /* SS */
      for (j = 1; j < 2; j+=2)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      SSvalue[j][k] = result.S[j].bold[k];
	      TSIL_Compare (ssname[j][k], val, SSvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
	      tally[foo]++;
	    }
	}
#else
      /* M */
      val = TSIL_GetComplex ();
      Mvalue = result.M.value;
      TSIL_CheckConsistent (Mvalue, TSIL_GetFunction(&result,"M"));
      TSIL_Compare ("M", val, Mvalue, TSIL_PASS, TSIL_WARN, &foo);
      tally[foo]++;

      /* U */
      for (j = 0; j < 4; j++)
	{
	  val = TSIL_GetComplex ();
	  Uvalue[j] = result.U[j].value;
	  TSIL_Compare (uname[j][0], val, Uvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Uvalue[0], TSIL_GetFunction(&result,"Uzxyv"));
      TSIL_CheckConsistent (Uvalue[1], TSIL_GetFunction(&result,"Uuyxv"));
      TSIL_CheckConsistent (Uvalue[2], TSIL_GetFunction(&result,"Uxzuv"));
      TSIL_CheckConsistent (Uvalue[3], TSIL_GetFunction(&result,"Uyuzv"));

      /* T */
      for (j = 0; j < 6; j++)
	{
	  val = TSIL_GetComplex ();
	  Tvalue[j] = result.T[j].value;
	  TSIL_Compare (tname[j][0], val, Tvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Tvalue[0], TSIL_GetFunction(&result,"Tvyz"));
      TSIL_CheckConsistent (Tvalue[1], TSIL_GetFunction(&result,"Tuxv"));
      TSIL_CheckConsistent (Tvalue[2], TSIL_GetFunction(&result,"Tyzv"));
      TSIL_CheckConsistent (Tvalue[3], TSIL_GetFunction(&result,"Txuv"));
      TSIL_CheckConsistent (Tvalue[4], TSIL_GetFunction(&result,"Tzyv"));
      TSIL_CheckConsistent (Tvalue[5], TSIL_GetFunction(&result,"Tvxu"));

      /* S */
      for (j = 0; j < 2; j++)
	{
	  val = TSIL_GetComplex ();
	  Svalue[j] = result.S[j].value;
	  TSIL_Compare (sname[j][0], val, Svalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Svalue[0], TSIL_GetFunction(&result,"Svyz"));
      TSIL_CheckConsistent (Svalue[1], TSIL_GetFunction(&result,"Suxv"));

      /* B */
      for (j = 0; j < 2; j++)
	{
	  val = TSIL_GetComplex ();
	  Bvalue[j] = result.B[j].value;
	  TSIL_Compare (bname[j][0], val, Bvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Bvalue[0], TSIL_GetFunction(&result,"Bxz"));
      TSIL_CheckConsistent (Bvalue[1], TSIL_GetFunction(&result,"Byu"));

      /* V */
      for (j = 0; j < 4; j++)
	{
	  val = TSIL_GetComplex ();
	  Vvalue[j] = result.V[j].value;
	  TSIL_Compare (vname[j][0], val, Vvalue[j], TSIL_PASS_V, TSIL_WARN_V, &foo);
	  tally[foo]++;
	}
      TSIL_CheckConsistent (Vvalue[0], TSIL_GetFunction(&result,"Vzxyv"));
      TSIL_CheckConsistent (Vvalue[1], TSIL_GetFunction(&result,"Vuyxv"));
      TSIL_CheckConsistent (Vvalue[2], TSIL_GetFunction(&result,"Vxzuv"));
      TSIL_CheckConsistent (Vvalue[3], TSIL_GetFunction(&result,"Vyuzv"));

      /* Tbar */
      for (j = 0; j < 6; j++)
	{
	  val = TSIL_GetComplex ();
	  Tbarvalue[j] = result.Tbar[j].value;
	  TSIL_Compare (tbarname[j][0], val, Tbarvalue[j], TSIL_PASS, TSIL_WARN, &foo);
	  tally[foo]++;
	}

      /* UU */
      for (j = 0; j < 4; j++)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      UUvalue[j][k] = result.U[j].bold[k];
	      TSIL_Compare (uuname[j][k], val, UUvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
	      tally[foo]++;
	    }
	}

      /* VV */
      for (j = 0; j < 4; j++)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      VVvalue[j][k] = result.V[j].bold[k];
	      TSIL_Compare (vvname[j][k], val, VVvalue[j][k], TSIL_PASS_V, TSIL_WARN_V, &foo);
	      tally[foo]++;
	    }
	}

      /* TT */
      for (j = 0; j < 6; j++)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      TTvalue[j][k] = result.T[j].bold[k];
	      TSIL_Compare (ttname[j][k], val, TTvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
	      tally[foo]++;
	    }
	}

      /* SS */
      for (j = 0; j < 2; j++)
	{
	  for (k = 0; k < 3; k++)
	    {
	      val = TSIL_GetComplex ();
	      SSvalue[j][k] = result.S[j].bold[k];
	      TSIL_Compare (ssname[j][k], val, SSvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
	      tally[foo]++;
	    }
	}
#endif

      if (tally[FAIL] > 0)
	{
	  nfail++;
	  printf ("FAILED FOR INPUT PARAMETERS:\n");
	  printf ("x = %.6lf;\n", (double) x);
	  printf ("y = %.6lf;\n", (double) y);
	  printf ("z = %.6lf;\n", (double) z);
	  printf ("u = %.6lf;\n", (double) u);
	  printf ("v = %.6lf;\n", (double) v);
	  printf ("s = %.6lf;\n", (double) s);
	  printf ("qq = %.6lf;\n", (double) qq);
	}
      else if (tally[WARN] > 0)
        {
	  printf ("PASS WITH WARNING\n");
	  nwarn++;
        }
      else
	{
	  printf ("PASS\n");
	  npass++;
	}

      if (1 == Doperms)
	{
	  for (p = 0; p < 4; p++)
	    {
	      for (j = 0; j < 3; j++)
		tally[j] = 0;

	      x = xx;
	      y = yy;
	      z = zz;
	      u = uu;
	      v = vv;

	      if (1 == p)
		{
		  TSIL_SwapR (&x, &y);
		  TSIL_SwapR (&z, &u);
		}

	      if (2 == p)
		{
		  TSIL_SwapR (&x, &z);
		  TSIL_SwapR (&y, &u);
		}

	      if (3 == p)
		{
		  TSIL_SwapR (&x, &u);
		  TSIL_SwapR (&y, &z);
		}

	      /* Calculate everything... */
	      TSIL_SetParameters (&result, x, y, z, u, v, qq);
	      TSIL_Evaluate (&result, s);

	      TSIL_PermuteResults (&result, p);

	      printf ("\nTest %d, permutation %d (", i, p);
	      printf ("%s", argv[i]);
	      printf ("): ");

	      /* ...and check that permuted results agree as they
                 should: */

	      /* M */
	      val = result.M.value;
	      TSIL_Compare ("M", val, Mvalue, TSIL_PASS, TSIL_WARN, &foo);
	      tally[foo]++;

	      /* U */
	      for (j = 0; j < 4; j++)
		{
		  val = result.U[j].value;
		  TSIL_Compare (uname[j][0], val, Uvalue[j], TSIL_PASS, TSIL_WARN, &foo);
		  tally[foo]++;
		}

	      /* T */
	      for (j = 0; j < 6; j++)
		{
		  val = result.T[j].value;
		  TSIL_Compare (tname[j][0], val, Tvalue[j], TSIL_PASS, TSIL_WARN, &foo);
		  tally[foo]++;
		}

	      /* S */
	      for (j = 0; j < 2; j++)
		{
		  val = result.S[j].value;
		  TSIL_Compare (sname[j][0], val, Svalue[j], TSIL_PASS, TSIL_WARN, &foo);
		  tally[foo]++;
		}

	      /* B */
	      for (j = 0; j < 2; j++)
		{
		  val = result.B[j].value;
		  TSIL_Compare (bname[j][0], val, Bvalue[j], TSIL_PASS, TSIL_WARN, &foo);
		  tally[foo]++;
		}

	      /* V */
	      for (j = 0; j < 4; j++)
		{
		  val = result.V[j].value;
		  TSIL_Compare (vname[j][0], val, Vvalue[j], TSIL_PASS_V, TSIL_WARN_V, &foo);
		  tally[foo]++;
		}

	      /* Tbar */
	      for (j = 0; j < 6; j++)
		{
		  val = result.Tbar[j].value;
		  TSIL_Compare (tbarname[j][0], val, Tbarvalue[j], TSIL_PASS, TSIL_WARN, &foo);
		  tally[foo]++;
		}

	      /* UU */
	      for (j = 0; j < 4; j++)
		{
		  for (k = 0; k < 3; k++)
		    {
		      val = result.U[j].bold[k];
		      TSIL_Compare (uuname[j][k], val, UUvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		}

	      /* VV */
	      for (j = 0; j < 4; j++)
		{
		  for (k = 0; k < 3; k++)
		    {
		      val = result.V[j].bold[k];
		      TSIL_Compare (vvname[j][k], val, VVvalue[j][k], TSIL_PASS_V, TSIL_WARN_V, &foo);
		      tally[foo]++;
		    }
		}

	      /* TT */
	      for (j = 0; j < 6; j++)
		{
		  for (k = 0; k < 3; k++)
		    {
		      val = result.T[j].bold[k];
		      TSIL_Compare (ttname[j][k], val, TTvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		}

	      /* SS */
	      for (j = 0; j < 2; j++)
		{
		  for (k = 0; k < 3; k++)
		    {
		      val = result.S[j].bold[k];
		      TSIL_Compare (ssname[j][k], val, SSvalue[j][k], TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		}

	      if (tally[FAIL] > 0)
		{
		  nfail++;
		  printf ("\nFAILED FOR INPUT PARAMETERS:\n");
		  printf ("x = %.3lf;\n", (double) x);
		  printf ("y = %.3lf;\n", (double) y);
		  printf ("z = %.3lf;\n", (double) z);
		  printf ("u = %.3lf;\n", (double) u);
		  printf ("v = %.3lf;\n", (double) v);
		  printf ("s = %.3lf;\n", (double) s);
		  printf ("qq = %.3lf;\n", (double) qq);
		}
	      else if (tally[WARN] > 0)
                {
		  nwarn++;
		  printf ("\nPASS WITH WARNING\n");
                }
	      else
		{
		  printf ("\nPASS\n");
		  npass++;
		}
	    }
	}

      if (1 == Doextraperms)
	{
	  for (p = 4; p < 8; p++)
	    {
	      x = xx;
	      y = yy;
	      z = zz;
	      u = uu;
	      v = vv;

	      if (p == 4)
		{
		  TSIL_SwapR (&x, &v);
		  permU = 1;
		  permT11 = 1;
		  permT12 = 1;
		  permT21 = 3;
		  permT22 = 5;
		  permT31 = 5;
		  permT32 = 3;
		  permS = 1;
		}
	      if (p == 5)
		{
		  TSIL_SwapR (&y, &v);
		  permU = 0;
		  permT11 = 0;
		  permT12 = 2;
		  permT21 = 2;
		  permT22 = 0;
		  permT31 = 4;
		  permT32 = 4;
		  permS = 0;
		}
	      if (p == 6)
		{
		  TSIL_SwapR (&z, &v);
		  permU = 3;
		  permT11 = 0;
		  permT12 = 4;
		  permT21 = 4;
		  permT22 = 0;
		  permT31 = 2;
		  permT32 = 2;
		  permS = 0;
		}
	      if (p == 7)
		{
		  TSIL_SwapR (&u, &v);
		  permU = 2;
		  permT11 = 1;
		  permT12 = 5;
		  permT21 = 5;
		  permT22 = 1;
		  permT31 = 3;
		  permT32 = 3;
		  permS = 1;
		}

	      for (j = 0; j < 3; j++)
		tally[j] = 0;

	      /* Calculate everything... */
	      TSIL_SetParameters (&result, x, y, z, u, v, qq);
	      TSIL_Evaluate (&result, s);

	      printf ("\nTest %d, permutation %d (", i, p);
	      printf ("%s", argv[i]);
	      printf ("): ");

	      /* ...and check results: */

	      /* U */
	      for (j = 0; j < 4; j++)
		{
		  val = Uvalue[j];
		  if (j == permU)
		    {
		      TSIL_Compare (uname[j][0], val, result.U[j].value, TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		}

	      /* T */
	      for (j = 0; j < 6; j++)
		{
		  val = Tvalue[j];
		  if (j == permT11)
		    {
		      TSIL_Compare (tname[j][0], val, result.T[permT12].value,
                                 TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		  if (j == permT21)
		    {
		      TSIL_Compare (tname[j][0], val, result.T[permT22].value,
				 TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		  if (j == permT31)
		    {
		      TSIL_Compare (tname[j][0], val, result.T[permT32].value,
				 TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		}

	      /* S */
	      for (j = 0; j < 2; j++)
		{
		  val = Svalue[j];
		  if (j == permS)
		    {
		      TSIL_Compare (sname[j][0], val, result.S[j].value, 
                                    TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		}

	      /* V */
	      for (j = 0; j < 4; j++)
		{
		  val = Vvalue[j];
		  if (j == permU)
		    {
		      TSIL_Compare (vname[j][0], val, result.V[j].value, 
                                    TSIL_PASS_V, TSIL_WARN_V, &foo);
		      tally[foo]++;
		    }
		}

	      /* Tbar */
	      for (j = 0; j < 6; j++)
		{
		  val = Tbarvalue[j];
		  if (j == permT11)
		    {
		      TSIL_Compare (tbarname[j][0], val, result.Tbar[permT12].value,
				 TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		  if (j == permT21)
		    {
		      TSIL_Compare (tbarname[j][0], val, result.Tbar[permT22].value,
				 TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		  if (j == permT31)
		    {
		      TSIL_Compare (tbarname[j][0], val, result.Tbar[permT32].value,
				 TSIL_PASS, TSIL_WARN, &foo);
		      tally[foo]++;
		    }
		}

	      /* UU */
	      for (j = 0; j < 4; j++)
		{
		  for (k = 0; k < 3; k++)
		    {
		      val = UUvalue[j][k];
		      if (j == permU)
			{
			  TSIL_Compare (uuname[j][k], val, result.U[j].bold[k],
				     TSIL_PASS, TSIL_WARN, &foo);
			  tally[foo]++;
			}
		    }
		}

	      /* VV */
	      for (j = 0; j < 4; j++)
		{
		  for (k = 0; k < 3; k++)
		    {
		      val = VVvalue[j][k];
		      if (j == permU)
			{
			  TSIL_Compare (vvname[j][k], val, result.V[j].bold[k],
				     TSIL_PASS_V, TSIL_WARN_V, &foo);
			  tally[foo]++;
			}
		    }
		}

	      /* TT */
	      for (j = 0; j < 6; j++)
		{
		  for (k = 0; k < 3; k++)
		    {
		      val = TTvalue[j][k];
		      if (j == permT11)
			{
			  TSIL_Compare (ttname[j][k], val,
				     result.T[permT12].bold[k], 
                                      TSIL_PASS, TSIL_WARN, &foo);
			  tally[foo]++;
			}
		      if (j == permT21)
			{
			  TSIL_Compare (ttname[j][k], val,
				     result.T[permT22].bold[k], TSIL_PASS, TSIL_WARN, &foo);
			  tally[foo]++;
			}
		      if (j == permT31)
			{
			  TSIL_Compare (ttname[j][k], val,
				     result.T[permT32].bold[k], TSIL_PASS, TSIL_WARN, &foo);
			  tally[foo]++;
			}
		    }
		}

	      /* SS */
	      for (j = 0; j < 2; j++)
		{
		  for (k = 0; k < 3; k++)
		    {
		      val = SSvalue[j][k];
		      if (j == permS)
			{
			  TSIL_Compare (ssname[j][k], val, result.S[j].bold[k],
				     TSIL_PASS, TSIL_WARN, &foo);
			  tally[foo]++;
			}
		    }
		}

	      if (tally[FAIL] > 0)
		{
		  nfail++;
		  printf ("\nFAILED FOR INPUT PARAMETERS:\n");
		  printf ("x = %.16lf;\n", (double) x);
		  printf ("y = %.16lf;\n", (double) y);
		  printf ("z = %.16lf;\n", (double) z);
		  printf ("u = %.16lf;\n", (double) u);
		  printf ("v = %.16lf;\n", (double) v);
		  printf ("s = %.16lf;\n", (double) s);
		  printf ("qq = %.16lf;\n", (double) qq);
		}
	      else if (tally[WARN] > 0)
                {
         	  nwarn++;
		  printf ("\nPASS\n");
                }
	      else
		{
		  printf ("\nPASS\n");
		  npass++;
		}
	    }
	}

      fclose (fp);

      printf ("\n===== Done with input file ");
      printf ("%s", argv[i]);
      printf (" =====\n");
    }

  printf ("\n== FINAL RESULTS ==\n");
  printf ("Total input files: %d\n", (argc - 1));
  printf ("Total tests performed: %d\n",
	  (1 + Doperms * 3 + Doextraperms * 4) * (argc - 1));
  printf ("Pass: %d\n", npass);
  printf ("Warn: %d\n", nwarn);
  printf ("Fail: %d\n", nfail);

  return 0;
}
