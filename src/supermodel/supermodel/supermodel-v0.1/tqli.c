/*
  Diagonalization of a real tridiagonal matrix using the QL algorithm
  with implicit shifts (based on a routine from Numerical Recipes).
  Modified to include macros for SIGN and pythag, and to restore
  standard C indexing.
*/

#include "supermodel.h"

#define SIGN(a,b) ((b)<0 ? -SUMO_FABS(a) : SUMO_FABS(a))
#define pythag(a,b) (SUMO_SQRT(a*a+b*b))

void tqli (SUMO_REAL *d, SUMO_REAL *e, int n, SUMO_REAL **z)
{
  int m,l,iter,i,k;
  SUMO_REAL s,r,p,g,f,dd,c,b;

  for (i=1; i<n; i++)
    e[i-1] = e[i];
  e[n-1] = 0.0;
  for (l=0; l<n; l++) {
    iter = 0;
    do {
      for (m=l; m<n-1; m++) {
	dd = SUMO_FABS(d[m]) + SUMO_FABS(d[m+1]);
	if ((SUMO_REAL)(SUMO_FABS(e[m]) + dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) 
	  SUMO_Error("tqli", "Too many iterations.", 42);
	g = (d[l+1] - d[l])/(2.0*e[l]);
	r = pythag(g,1.0);
	g = d[m] - d[l] + e[l]/(g + SIGN(r,g));
	s = c = 1.0;
	p = 0.0;
	for (i=m-1;i>=l;i--) {
	  f = s*e[i];
	  b = c*e[i];
	  e[i+1] = (r = pythag(f,g));     
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m] = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = d[i+1] - p;
	  r = (d[i] - g)*s + 2.0*c*b;
	  d[i+1] = g + (p = s*r);
	  g = c*r - b;
	  for (k=0; k<n; k++) {
	    f = z[k][i+1];
	    z[k][i+1] = s*z[k][i] + c*f;
	    z[k][i] = c*z[k][i] - s*f;
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    } while (m != l);
  }
}
