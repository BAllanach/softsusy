/*
  Finds eigenvalues and right eigenvectors of a general complex
  matrix. Algorithm is that given by Eberlein in Wilkinson and
  Reinsch, _Handbook for Automatic Computation v.2_ (Springer,
  1971). If the matrix is normal (hermitian or skew-hermitian) then
  the algorithm reduces to the Jacobi method.  In the 2x2 case the
  explicit analytic formulae are used.

  **AT PRESENT ONLY EIGENVALUES ARE COMPUTED FOR THIS CASE**

  NOTE: The matrix mat *must* be defined in the calling program so
  that it occupies a contiguous block of memory!! The declaration here
  (TSIL_COMPLEX *mat) is designed so that the routine will work with
  any matrix dimension, with statically declared matrices in the
  calling program.  This requires that the matrix be "flattened" in
  the present routine.
*/

#include "supermodel.h"

/* ------------------------------------------------------------------ */

int SUMO_DiagonalizeComp (SUMO_COMPLEX *mat,
			  int size,
			  SUMO_EIGEN_COMP *res)
{
  TSIL_COMPLEX m00, m01, m10, m11, tmp;
  TSIL_REAL a[size][size], z[size][size];
/*   TSIL_REAL t[size][size], u[size][size]; */
  TSIL_REAL en[size];
  TSIL_REAL eps, tau, tem, tep, max, hj, hr, hi, g;
  TSIL_REAL te, tee, br, bi, er, ei, dr, di, root;
  TSIL_REAL root1, root2, ca, sa, cb, sb, cx, sx, sh, ch,eta;
  TSIL_REAL tse, nd, nc, s, c, cot2x, d, de, sig, cotx, cos2a;
  TSIL_REAL sin2a, tanh, c1r, c2r, s1r, s2r, c1i, c2i, s1i, s2i;
  TSIL_REAL isw, b, e, zim, zik, aik, aim;
/*   TSIL_REAL tim, tik, uim, uik; */
  TSIL_REAL ami, aki, zki, zmi;
  int i, j, k, m, n;
  int it, mark;

/*   static int foo=0; */
/*   printf("count = %d\n", ++foo); */

  /* If 2x2, just use explicit results... */
  if (size == 2) {

    m00 = *mat;
    m01 = *(mat + 1);
    m10 = *(mat + 2);
    m11 = *(mat + 3);

    /* Fill up the result struct: */
    res->size = size;

    if (SUMO_CABS(m10) > TSIL_TOL) {
      tmp = SUMO_CSQRT(m00*m00 + 4.0L*m01*m10 - 2.0L*m00*m11 + m11*m11);
      res->eigenvalues[0] = 0.5L*(m00 + m11 - tmp);
      res->eigenvalues[1] = 0.5L*(m00 + m11 + tmp);
/*       res->eigenvectors[0][0] = 0.5L*(m00 - m11 - tmp)/m10; */
/*       res->eigenvectors[1][0] = 1.0L; */
/*       res->eigenvectors[0][1] = 0.5L*(m00 - m11 + tmp)/m10; */
/*       res->eigenvectors[1][1] = 1.0L; */

      /* Normalize eigenvectors (?) */
/*       for (j=0; j<size; j++) { */
/* 	tmp = 0.0; */
/* 	for (i=0; i<size; i++) */
/* 	  tmp += res->eigenvectors[i][j]*SUMO_CONJ(res->eigenvectors[i][j]); */
/* 	tmp = SUMO_CSQRT(tmp); */
/* 	for (i=0; i<size; i++) */
/* 	  res->eigenvectors[i][j] /= tmp; */
/*       } */
    }
    else {
      if (TSIL_CABS(m00 - m11) > TSIL_TOL) {
	res->eigenvalues[0] = m00;
	res->eigenvalues[1] = m11;
/* 	res->eigenvectors[0][0] = 1.0L; */
/* 	res->eigenvectors[1][0] = 0.0L; */
/* 	res->eigenvectors[0][1] = m01/(m11 - m00); */
/* 	res->eigenvectors[1][1] = 1.0L; */

	/* Normalize last eigenvector (?) */
/* 	tmp = 0.0; */
/* 	for (i=0; i<size; i++) */
/* 	  tmp += res->eigenvectors[i][1]*SUMO_CONJ(res->eigenvectors[i][1]); */
/* 	tmp = SUMO_CSQRT(tmp); */
/* 	for (i=0; i<size; i++) */
/* 	  res->eigenvectors[i][1] /= tmp; */
      }
      else {
	/* This is a degenerate case; hopefully it never happens! */
	res->eigenvalues[0] = m00;
	res->eigenvalues[1] = m00;
/* 	res->eigenvectors[0][0] = 1.0L; */
/* 	res->eigenvectors[1][0] = 0.0L; */
/* 	res->eigenvectors[0][1] = 1.0L; */
/* 	res->eigenvectors[1][1] = 0.0L; */
      }
    }
    /* Sort eigenvalues in ascending order (of real part) if necessary */
    if (SUMO_CREAL(res->eigenvalues[0]) > SUMO_CREAL(res->eigenvalues[1])) {
      tmp = res->eigenvalues[0];
      res->eigenvalues[0] = res->eigenvalues[1];
      res->eigenvalues[1] = tmp;
/*       for (i=0; i<size; i++) { */
/* 	tmp = res->eigenvectors[i][0]; */
/* 	res->eigenvectors[i][0] = res->eigenvectors[i][1]; */
/* 	res->eigenvectors[i][1] = tmp; */
/*       } */
    }
  }
  else {   /* Bigger than 2x2 handled numerically */

    n = size;
    it = 0;
    mark = FALSE;
    eps = 1.e-14;

    /* Basic setup */
    for (i=0; i<n; i++) {
/*       t[i][i] = 1.0; */
/*       u[i][i] = 0.0; */
/*       for (j=i+1; j<n; j++) */
/* 	t[i][j] = t[j][i] = u[i][j] = u[j][i] = 0.0; */
      for (j=0; j<n; j++) {
	a[i][j] = SUMO_CREAL(mat[i*n + j]);
	z[i][j] = SUMO_CIMAG(mat[i*n + j]);
      }
    }
    for (it=0; it<35; it++) {
      if (mark == TRUE) goto done;
      tau = 0.0;
      for (k=0; k<n; k++) {
	tem = 0.0;
	for (i=0; i<n; i++)
	  if (i != k) tem += SUMO_FABS(a[i][k]) + SUMO_FABS(z[i][k]);
	tau += tem;
	en[k] = tem + SUMO_FABS(a[k][k]) + SUMO_FABS(z[k][k]);
      }
      for (k=0; k<n-1; k++) {
	max = en[k];
	i = k;
	for (j=k+1; j<n; j++) {
	  if (en[j] > max) {
	    max = en[j];
	    i = j;
	  }
	}
	if (i != k) {
	  en[i] = en[k];
	  for (j=0; j<n; j++) {
	    tep = a[k][j]; a[k][j] = a[i][j]; a[i][j] = tep;
	    tep = z[k][j]; z[k][j] = z[i][j]; z[i][j] = tep;
	  }
	  for (j=0; j<n; j++) {
	    tep = a[j][k]; a[j][k] = a[j][i]; a[j][i] = tep;
	    tep = z[j][k]; z[j][k] = z[j][i]; z[j][i] = tep;
/* 	    tep = t[j][k]; t[j][k] = t[j][i]; t[j][i] = tep; */
/* 	    tep = u[j][k]; u[j][k] = u[j][i]; u[j][i] = tep; */
	  }
	}
      }
      if (tau < 100.0*eps) goto done;

      /* Now we begin a sweep of the matrix */
      mark = TRUE;
      for (k=0; k<n-1; k++) {
	for (m=k+1; m<n; m++) {
	  hj = hr = hi = g = 0.0;
	  for (i=0; i<n; i++) {
	    if (i!=k && i!=m) {
	      hr += a[k][i]*a[m][i] + z[k][i]*z[m][i] -
		a[i][k]*a[i][m] - z[i][k]*z[i][m];
	      hi += z[k][i]*a[m][i] - a[k][i]*z[m][i] -
		a[i][k]*z[i][m] + z[i][k]*a[i][m];
	      te = a[i][k]*a[i][k] + z[i][k]*z[i][k] +
		a[m][i]*a[m][i] + z[m][i]*z[m][i];
	      tee = a[i][m]*a[i][m] + z[i][m]*z[i][m] +
		a[k][i]*a[k][i] + z[k][i]*z[k][i];
	      g += te + tee;
	      hj += -te + tee;
	    }
	  }
	  br = a[k][m] + a[m][k]; bi = z[k][m] + z[m][k];
	  er = a[k][m] - a[m][k]; ei = z[k][m] - z[m][k];
	  dr = a[k][k] - a[m][m]; di = z[k][k] - z[m][m];
	  te = br*br + ei*ei + dr*dr;
	  tee = bi*bi + er*er + di*di;
	  if (te >= tee) {
	    isw = 1.0; c = br; s = ei; d = dr; de = di;
	    root2 = SUMO_SQRT(te);
	  }
	  else {
	    isw = -1.0; c = bi; s = -er; d = di; de = dr;
	    root2 = SUMO_SQRT(tee);
	  }
	  root1 = SUMO_SQRT(s*s + c*c);
	  if (d >= 0.0)
	    sig = 1.0;
	  else
	    sig = -1.0;
	  sa = 0.0;
	  if (c >= 0.0)
	    ca = 1.0;
	  else
	    ca = -1.0;
	  if (root1 < eps) {
	    sx = sa = 0.0;
	    cx = ca = 1.0;
	    if (isw > 0.0) {
	      e = er;
	      b = bi;
	    }
	    else {
	      e = ei;
	      b = br;
	    }
	    nd = d*d + de*de;
	    goto enter1;
	  }
	  if (SUMO_FABS(s) > eps) {
	    ca = c/root1;
	    sa = s/root1;
	  }
	  cot2x = d/root1;
	  cotx = cot2x + (sig*SUMO_SQRT(1.0L + cot2x*cot2x));
	  sx = sig/SUMO_SQRT(1.0 + cotx*cotx);
	  cx = sx*cotx;
	
	  /* Find rotated elements */
	  eta = (er*br + bi*ei)/root1;
	  tse = (br*bi - er*ei)/root1;
	  te = sig*(-root1*de + tse*d)/root2;
	  tee = (d*de + root1*tse)/root2;
	  nd = root2*root2 + tee*tee;
	  tee = hj*cx*sx;
	  cos2a = ca*ca - sa*sa;
	  sin2a = 2.0L*ca*sa;
	  tem = hr*cos2a + hi*sin2a;
	  tep = hi*cos2a - hr*sin2a;
	  hr = cx*cx*hr - sx*sx*tem - ca*tee;
	  hi = cx*cx*hi + sx*sx*tep - sa*tee;
	  b = isw*te*ca + eta*sa;
	  e = ca*eta - isw*te*sa;
	
	enter1:
	  s = hr - sig*root2*e;
	  c = hi - sig*root2*b;
	  root = SUMO_SQRT(c*c + s*s);
	  if (root < eps) {
	    cb = ch = 1.0;
	    sb = sh = 0.0;
	    goto trans;
	  }
	  cb = -c/root; 
	  sb = s/root;
	  tee = cb*b - e*sb;
	  nc = tee*tee;
	  tanh = root/(g + 2.0L*(nc + nd));
	  ch = 1.0L/SUMO_SQRT(1.0L - tanh*tanh);
	  sh = ch*tanh;
	
	  /* prepare for transformation */
	trans:
	  tem = sx*sh*(sa*cb - sb*ca);
	  c1r = cx*ch - tem;
	  c2r = cx*ch + tem;
	  c1i = c2i = -sx*sh*(ca*cb + sa*sb);
	  tep = sx*ch*ca;
	  tem = cx*sh*sb;
	  s1r = tep - tem;
	  s2r = -tep - tem;
	  tep = sx*ch*sa;
	  tem = cx*sh*cb;
	  s1i = tep + tem;
	  s2i = tep - tem;

	  /* Decide whether to make transformation */
	  tem = SUMO_SQRT(s1r*s1r + s1i*s1i);
	  tep = SUMO_SQRT(s2r*s2r + s2i*s2i);
	  
	  if (tem > eps || tep > eps) {
	    mark = FALSE;
	    /* Transformation on the left... */
	    for (i=0; i<n; i++) {
	      aki = a[k][i]; ami = a[m][i];
	      zki = z[k][i]; zmi = z[m][i];
	      a[k][i] = c1r*aki - c1i*zki + s1r*ami - s1i*zmi;
	      z[k][i] = c1r*zki + c1i*aki + s1r*zmi + s1i*ami;
	      a[m][i] = s2r*aki - s2i*zki + c2r*ami - c2i*zmi;
	      z[m][i] = s2r*zki + s2i*aki + c2r*zmi + c2i*ami;
	    }
	    /* Transformation on the right */
	    for (i=0; i<n; i++) {
	      aik = a[i][k]; aim = a[i][m];
	      zik = z[i][k]; zim = z[i][m];
	      a[i][k] = c2r*aik - c2i*zik - s2r*aim + s2i*zim;
	      z[i][k] = c2r*zik + c2i*aik - s2r*zim - s2i*aim;
	      a[i][m] = -s1r*aik + s1i*zik + c1r*aim - c1i*zim;
	      z[i][m] = -s1r*zik - s1i*aik + c1r*zim + c1i*aim;
	      
/* 	      tik = t[i][k]; tim = t[i][m]; */
/* 	      uik = u[i][k]; uim = u[i][m]; */
/* 	      t[i][k] = c2r*tik - c2i*uik - s2r*tim + s2i*uim; */
/* 	      u[i][k] = c2r*uik + c2i*tik - s2r*uim - s2i*tim; */
	    
/* 	      t[i][m] = -s1r*tik + s1i*uik + c1r*tim - c1i*uim; */
/* 	      u[i][m] = -s1r*uik - s1i*tik + c1r*uim + c1i*tim; */
	    }
	  }
	}
      }
    }
    /* If we get here we had too many iterations... */
    /*   SUMO_Error("DiagComplex", "Too many iterations!", 333); */
    SUMO_Warn("DiagComplex", "Too many iterations!");

  done:
    /*   printf("(%d iterations)\n", it+1); */
    res->size = n;
    for (i=n-1; i>=0; i--) {
      res->eigenvalues[n-1-i] = a[i][i] + I*z[i][i];
/*       for (j=0; j<n; j++) */
/* 	res->eigenvectors[j][n-1-i] = t[j][i] + I*u[j][i]; */
    }
  }

  /* Now make the largest entry in each eigenvector real and positive, 
     just to establish a standard presentation. */

/*   for (i=0; i<n; i++){ */
/*     tep = res->eigenvectors[0][i]; */
/*     for (j=1; j<n; j++){ */
/*       if (TSIL_CABS(res->eigenvectors[j][i]) > TSIL_CABS(tep)) */
/* 	tep = res->eigenvectors[j][i]; */
/*     } */
/*     for (j=0; j<n; j++){ */
/*       res->eigenvectors[j][i] *= TSIL_CABS(tep)/tep; */
/*     } */
/*   } */
  return 0;
}
