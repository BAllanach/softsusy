#ifndef _SUMO_LINALG_H_
#define _SUMO_LINALG_H_


/* -------------------------------------------------------------------- */

struct SUMO_eigen_2x2 
{
  TSIL_REAL eigenvalues[2];
  TSIL_COMPLEX eigenvectors[2][2];
};

typedef struct SUMO_eigen_2x2 SUMO_EIGEN_2x2;

/* -------------------------------------------------------------------- */

struct SUMO_eigen_4x4 
{
  TSIL_REAL eigenvalues[4];
  TSIL_COMPLEX eigenvectors[4][4];
};

typedef struct SUMO_eigen_4x4 SUMO_EIGEN_4x4;

/* -------------------------------------------------------------------- */

struct SUMO_eigen_6x6 
{
  TSIL_REAL eigenvalues[6];
  TSIL_COMPLEX eigenvectors[6][6];
};

typedef struct SUMO_eigen_6x6 SUMO_EIGEN_6x6;

/* -------------------------------------------------------------------- */
/* Unified struct, slightly wasteful of memory but... */

struct SUMO_eigen_result
{
  int size;
  TSIL_COMPLEX eigenvalues[6];
  TSIL_COMPLEX eigenvectors[6][6];
};

typedef struct SUMO_eigen_result SUMO_EIGEN_RESULT;

/* -------------------------------------------------------------------- */
/* Unified struct, slightly wasteful of memory but... */

struct SUMO_eigen_herm
{
  int          size;
  TSIL_REAL    eigenvalues[4];
  TSIL_COMPLEX eigenvectors[4][4];
};

typedef struct SUMO_eigen_herm SUMO_EIGEN_HERM;

/* -------------------------------------------------------------------- */
/* Unified struct, slightly wasteful of memory but... */

struct SUMO_eigen_comp
{
  int          size;
  TSIL_COMPLEX eigenvalues[4];
  TSIL_COMPLEX eigenvectors[4][4];
};

typedef struct SUMO_eigen_comp SUMO_EIGEN_COMP;




#endif /* sumo_linalg.h */
