/*
  Reduction of a real symmetric matrix to tridiagonal form using the
  Householder technique (based on the routine from Numerical Recipes).
  Modified to restore standard C array indexing and pointer arguments.
*/

#include "supermodel.h"

void tred2 (SUMO_REAL **a, int n, SUMO_REAL *d, SUMO_REAL *e)
{
  int l,k,j,i;
  SUMO_REAL scale,hh,h,g,f;
  
  for (i=n-1; i>=1; i--) {
    l = i-1;
    h = scale = 0.0;
    if (l > 0) {
      for (k=0; k<=l; k++)
	scale += SUMO_FABS(a[i][k]);
      if (scale == 0.0)
	e[i] = a[i][l];
      else {
	for (k=0; k<=l; k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f = a[i][l];
	g = (f >= 0.0 ? -SUMO_SQRT(h) : SUMO_SQRT(h));
	e[i] = scale*g;
	h -= f*g;
	a[i][l] = f-g;
	f = 0.0;
	for (j=0; j<=l; j++) {
	  a[j][i] = a[i][j]/h;
	  g = 0.0;
	  for (k=0; k<=j; k++)
	    g += a[j][k]*a[i][k];
	  for (k=j+1; k<=l; k++)
	    g += a[k][j]*a[i][k];
	  e[j] = g/h;
	  f += e[j]*a[i][j];
	}
	hh = f/(h+h);
	for (j=0; j<=l; j++) {
	  f = a[i][j];
	  e[j] = g = e[j]-hh*f;
	  for (k=0; k<=j; k++)
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i] = a[i][l];
    d[i] = h;
  }

  d[0] = 0.0;
  e[0] = 0.0;

  for (i=0; i<n; i++) {
    l = i-1;
    if (d[i]) {
      for (j=0; j<=l; j++) {
	g = 0.0;
	for (k=0; k<=l; k++)
	  g += a[i][k]*a[k][j];
	for (k=0; k<=l; k++)
	  a[k][j] -= g*a[k][i];
      }
    }
    d[i] = a[i][i];
    a[i][i] = 1.0;
    for (j=0; j<=l; j++) a[j][i] = a[i][j] = 0.0;
  }
}
