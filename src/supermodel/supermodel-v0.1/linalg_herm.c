/*
  This file contains routines to diagonalize Hermitian matrices,
  returning the eigenvalues and orthonormal eigenvectors in a
  consistent order.

  The eigenpair struct is organized so that eigenvectors[j][i]
  contains the jth element of the ith eigenvector. In other words,
  since indices are in the standard order (row, column), the ith
  column is the ith eigenvector of the input matrix. The eigenvalues
  are in the same order as the eigenvectors.

  Diagonalization is accomplished either (1) by using the explicit
  solution (in the 2x2 case) or; (2) by first constructing a
  corresponding real, symmetric matrix of dimension (2 x size). This
  is then reduced to tri-diagonal form and diagonaized by a QL
  algorithm with implicit shifts. Eigenvalues and -vectors of the
  original matrix are produced in pairs. Routines are based on those
  found in Numerical Recipes.
*/

#include "supermodel.h"

void tred2(TSIL_REAL **, int, TSIL_REAL [], TSIL_REAL []);
void tqli(TSIL_REAL [], TSIL_REAL [], int, TSIL_REAL **);

/* Globals in the current scope */
SUMO_REAL **rm;
int sizex2;


/* ------------------------------------------------------------------ */

int SUMO_AreEigenvectorsTheSame (int i, int j)
{
  int k, size = sizex2/2;
  SUMO_COMPLEX check = 0.0 + 0.0*I;

  /* Test orthogonality */
  for (k=0; k<size; k++)
    check += (rm[k][i] - I*rm[k+size][i])*(rm[k][j] + I*rm[k+size][j]);

  if (TSIL_CABS(check) < 1.e-9)
    return NO;
  else
    return YES;
}

/* ------------------------------------------------------------------ */

void SUMO_SwapEigenvectors (int i, int j)
{
  int k;
  SUMO_REAL tmp;

  for (k=0; k<sizex2; k++) {
    tmp = rm[k][i];
    rm[k][i] = rm[k][j];
    rm[k][j] = tmp;
  }
}

/* ------------------------------------------------------------------ */

int SUMO_DiagonalizeHerm (SUMO_COMPLEX *mat,
			  int size,
			  SUMO_EIGEN_HERM *res)
{
  int i, j, k;
  SUMO_COMPLEX tempmax, m00, m01, m10, m11, tmp;
  SUMO_REAL *d, *e, p, rtmp;
  int cur, next, nextidentical;

/*   static int foo=0; */
/*   printf("count = %d\n", ++foo); */

  /* If 2x2, just use explicit results... */
  if (size == 2) {

/*     printf("In DiagonalizeHerm, computing explicitly (2x2)...\n"); */

    m00 = *mat;
    m01 = *(mat + 1);
    m10 = *(mat + 2);
    m11 = *(mat + 3);

    /* Fill up the result struct: */
    res->size = size;

    if (SUMO_CABS(m10) > TSIL_TOL) {
      tmp = SUMO_CSQRT(m00*m00 + 4.0L*m01*m10 - 2.0L*m00*m11 + m11*m11);
      res->eigenvalues[0] = SUMO_CREAL(0.5L*(m00 + m11 - tmp));
      res->eigenvalues[1] = SUMO_CREAL(0.5L*(m00 + m11 + tmp));
      res->eigenvectors[0][0] = 0.5L*(m00 - m11 - tmp)/m10;
      res->eigenvectors[1][0] = 1.0L;
      res->eigenvectors[0][1] = 0.5L*(m00 - m11 + tmp)/m10;
      res->eigenvectors[1][1] = 1.0L;

      /* Normalize eigenvectors */
      for (j=0; j<size; j++) {
	tmp = 0.0;
	for (i=0; i<size; i++)
	  tmp += res->eigenvectors[i][j]*SUMO_CONJ(res->eigenvectors[i][j]);
	tmp = SUMO_CSQRT(tmp);
	for (i=0; i<size; i++)
	  res->eigenvectors[i][j] /= tmp;
      }
    }
    else {
      res->eigenvalues[0] = SUMO_CREAL(m00);
      res->eigenvalues[1] = SUMO_CREAL(m11);
      res->eigenvectors[0][0] = 1.0L;
      res->eigenvectors[1][0] = 0.0L;
      res->eigenvectors[0][1] = 0.0L;
      res->eigenvectors[1][1] = 1.0L;
    }

    /* Sort eigenvalues in ascending order if necessary */
    if (res->eigenvalues[0] > res->eigenvalues[1]) {
      rtmp = res->eigenvalues[0];
      res->eigenvalues[0] = res->eigenvalues[1];
      res->eigenvalues[1] = rtmp;
      for (i=0; i<size; i++) {
	tmp = res->eigenvectors[i][0];
	res->eigenvectors[i][0] = res->eigenvectors[i][1];
	res->eigenvectors[i][1] = tmp;
      }
    }
  }
  else {
    /* Matrices bigger than 2x2 are handled numerically... */

/*     printf("In DiagonalizeHerm, diagonalizing numerically (QL)...\n"); */

    /* Size of corresponding real symmetric matrix: */
    sizex2 = 2*size;

    /* Allocate local arrays: */
    rm = (SUMO_REAL **) malloc (sizex2 * sizeof(SUMO_REAL *));
    rm[0] = (SUMO_REAL *) malloc (sizex2 * sizex2 * sizeof(SUMO_REAL));
    for (i=1; i<sizex2; i++) {
      rm[i] = rm[0] + sizex2*i;
    }
    d = (SUMO_REAL *) malloc (sizex2 * sizeof(SUMO_REAL));
    e = (SUMO_REAL *) malloc (sizex2 * sizeof(SUMO_REAL));

    /* Build the corresponding real symmetric matrix: */

    /* Upper left block */
    for (i=0; i<size; i++)
      for (j=0; j<size; j++)
	rm[i][j] = SUMO_CREAL(mat[i*size + j]);

    /* Lower left block */
    for (i=size; i<sizex2; i++)
      for (j=0; j<size; j++)
	rm[i][j] = SUMO_CIMAG(mat[(i-size)*size + j]);

    /* Upper right block */
    for (i=0; i<size; i++)
      for (j=size; j<sizex2; j++)
	rm[i][j] = rm[j][i];

    /* Lower right block */
    for (i=size; i<sizex2; i++)
      for (j=size; j<sizex2; j++)
	rm[i][j] = rm[i-size][j-size];


  /* DGR testing */
/*   printf("In DiagonalizeHerm...\n"); */
/*   printf("\nInput matrix:"); */
/*   for (i=0; i<size; i++) { */
/*     printf("\n"); */
/*       for (j=0; j<size; j++) */
/* 	printf("% Lf + % Lf I\t\t", SUMO_CREAL(mat[i*size+j]), SUMO_CIMAG(mat[i*size+j])); */
/*   } */

/*   printf("\n\nEquivalent real matrix:\n"); */
/*   for (i=0; i<sizex2; i++) { */
/*     printf("\n"); */
/*     for (j=0; j<sizex2; j++) */
/*       printf("% Lf\t", rm[i][j]); */
/*   } */
/*   printf("\n"); */

    /* Diagonalize it: */
    tred2 (rm, sizex2, d, e);
    tqli (d, e, sizex2, rm);

/*   printf("\nRaw eigenvalues:\n\n"); */
/*   for (i=0; i<sizex2; i++)  */
/*     printf("% Lf\n", d[i]); */

/*   printf("\nRaw eigenvectors:\n"); */
/*   for (i=0; i<sizex2; i++) { */
/*     printf("\n"); */
/*     for (j=0; j<sizex2; j++) */
/*       printf("% Lf\t", rm[i][j]); */
/*   } */
/*   printf("\n"); */

    /* Extract eigenvalues and vectors: */

    /* First, sort in ascending order */
    for (i=0; i<sizex2-1; i++) {
      p = d[k=i];
      for (j=i+1; j<sizex2; j++)
	if (d[j] < p) p = d[k=j];
      if (k != i) {
	d[k] = d[i];
	d[i] = p;
	for (j=0; j<sizex2; j++) {
	  p = rm[j][i];
	  rm[j][i] = rm[j][k];
	  rm[j][k] = p;
	}
      }
    }

/*   printf("\nSorted eigenvalues:\n\n"); */
/*   for (i=0; i<sizex2; i++)  */
/*     printf("% Lf\n", d[i]); */

/*   printf("\nSorted eigenvectors:\n"); */
/*   for (i=0; i<sizex2; i++) { */
/*     printf("\n"); */
/*     for (j=0; j<sizex2; j++) */
/*       printf("% Lf\t", rm[i][j]); */
/*   } */
/*   printf("\n"); */


    /* Next, insure that corresponding eigenvectors are adjacent: */
    cur = 0;
    while (cur < sizex2 - 2) { /* Don't want to do the last pair */
      next = cur + 1;
      if (SUMO_FABS(d[cur] - d[next]) > 1.e-6)
	SUMO_Error ("linalg", "This can't happen!", 222);

      for (nextidentical = next; nextidentical < sizex2; nextidentical++) {
	if (SUMO_AreEigenvectorsTheSame (cur, nextidentical) == NO)
	  ;
	else {
	  if (nextidentical != next)
	    SUMO_SwapEigenvectors (next, nextidentical);
	  break;
	}
      }
      cur += 2;
    }

/*   printf("\nRe-Sorted eigenvalues:\n\n"); */
/*   for (i=0; i<sizex2; i++)  */
/*     printf("% Lf\n", d[i]); */

/*   printf("\nRe-Sorted eigenvectors:\n"); */
/*   for (i=0; i<sizex2; i++) { */
/*     printf("\n"); */
/*     for (j=0; j<sizex2; j++) */
/*       printf("% Lf\t", rm[i][j]); */
/*   } */
/*   printf("\n"); */

/*   exit(0); */

    /* Fill up the result struct: */
    res->size = size;
    for (i=0; i<size; i++) {
      j = 2*i;
      res->eigenvalues[i] = d[j];
      for (k=0; k<size; k++)
	res->eigenvectors[k][i] = rm[k][j] + I*rm[k+size][j];
    }
    free ((void *) e);
    free ((void *) d);
    free ((void *) rm[0]);
    free ((void *) rm);
  }

  /* Finally, make the largest entry in each eigenvector real and
     positive, just to establish a standard presentation. */

  for (i = 0; i < size; i++){
    tempmax = res->eigenvectors[0][i];
    for (j = 1; j < size; j++){
      if (SUMO_CABS(res->eigenvectors[j][i]) > SUMO_CABS(tempmax)) 
	tempmax = res->eigenvectors[j][i];
    }
    for (j = 0; j < size; j++){
      res->eigenvectors[j][i] *= SUMO_CABS(tempmax)/tempmax;
    }
  }
  return 0;
}
