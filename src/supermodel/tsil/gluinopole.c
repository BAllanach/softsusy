/*   === gluinopole.c ===
 *
 *   This version is dated October 6, 2005.
 *
 *   Implements sections III, IV, and V.A of hep-ph/0509115
 *
 *   To use this, you first need to install TSIL, see:
 *   http://www.niu.edu/spmartin/TSIL
 *   or
 *   http://faculty.otterbein.edu/DRobertson/tsil
 *
 *   To compile this program: 
 *
 *   gcc -o gp gluinopole.c -L<dir> -ltsil -lm -DTSIL_SIZE_LONG
 *
 *   <dir> is the directory where libtsil.a is, and tsil.h needs to
 *   be copied to the current directory. 
 *
 *   No command line arguments, the input data is contained in main().
 *   Run as: ./gp
 *
 *   For authorship and citation information, see:
 *   http://www.niu.edu/spmartin/gluinopole
 */

#include <stdio.h>
#include "tsil.h"  /* TSIL header file */

/* Facts about QCD: */
#define CG 3.L
#define Cq (4.L/3.L)
#define Iq 0.5L

/* Other definitions */
#define PI 3.1415926535897932385L
#define PI2 9.8696044010893586188L
#define ln2 0.69314718055994530942L
#define kappa (1.L/(16.L*PI2))
#define Zeta3 1.2020569031595942854L

/* The following struct type should eventually be replaced by a grander 
struct containing the values of ALL parameters of the MSSM. */

struct MSSM_data {
   TSIL_REAL qq; /* Renormalization scale squared */
   TSIL_REAL g3; /* QCD gauge coupling */
   TSIL_REAL gluino; /* Gluino running squared mass */
   TSIL_REAL quark[6]; /* Quark squared masses in the order d,u,s,c,b,t */
   TSIL_REAL squark[6][2]; /* Squark squared masses */
   TSIL_COMPLEX Rsq[6][2]; /* Squark mixing matrices */
   TSIL_COMPLEX Lsq[6][2]; /* Squark mixing matrices */
};

typedef struct MSSM_data MSSM_DATA;

int bFS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL s, TSIL_REAL qq,
         TSIL_COMPLEX *BFS, TSIL_COMPLEX *BfS,
         TSIL_COMPLEX *BpFS, TSIL_COMPLEX *BpfS);

int bFV (TSIL_REAL x, TSIL_REAL v, TSIL_REAL s, TSIL_REAL qq, TSIL_REAL xi,
         TSIL_COMPLEX *BFV, TSIL_COMPLEX *BfV);

int mSFFSF (TSIL_DATA *input_data, 
            TSIL_COMPLEX *MSFFSF, TSIL_COMPLEX *MSFFSf, TSIL_COMPLEX *MSFfSF, 
            TSIL_COMPLEX *MSFfSf, TSIL_COMPLEX *MSffSF, TSIL_COMPLEX *MSffSf); 

int mSSFFS (TSIL_DATA *input_data, 
            TSIL_COMPLEX *MSSFFS, TSIL_COMPLEX *MSSFfS, TSIL_COMPLEX *MSSffS); 

int yFSSS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u,
           TSIL_REAL s, TSIL_REAL qq,
           TSIL_COMPLEX *YFSSS, TSIL_COMPLEX *YfSSS);

int vFSSSS (TSIL_DATA *input_data1, TSIL_DATA *input_data2,
            TSIL_COMPLEX *VFSSSS, TSIL_COMPLEX *VfSSSS);

int vSFFFS (TSIL_DATA *input_data1, TSIL_DATA *input_data2,
            TSIL_COMPLEX *VSFFFS, TSIL_COMPLEX *VSFFfS, TSIL_COMPLEX *VSFfFS,
            TSIL_COMPLEX *VSFffS, TSIL_COMPLEX *VSffFS, TSIL_COMPLEX *VSfffS);

int vFSSFF (TSIL_DATA *input_data1, TSIL_DATA *input_data2,
            TSIL_COMPLEX *VFSSFF, TSIL_COMPLEX *VFSSff, 
            TSIL_COMPLEX *VfSSFF, TSIL_COMPLEX *VfSSff);

int f1f4 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq, 
          TSIL_COMPLEX *resultf1, TSIL_COMPLEX *resultf4);

int f2f3f5f6 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq, 
              TSIL_COMPLEX *resultf2, TSIL_COMPLEX *resultf3,
              TSIL_COMPLEX *resultf5, TSIL_COMPLEX *resultf6);

TSIL_REAL F1 (TSIL_REAL x, TSIL_REAL qq);
TSIL_REAL F2 (TSIL_REAL x, TSIL_REAL qq);
TSIL_REAL F3 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL qq);
TSIL_REAL F4 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL qq);

TSIL_REAL f (TSIL_REAL r);

int gluinopole (MSSM_DATA *mssm_data,
                TSIL_REAL *Mpole1, TSIL_REAL *Gammapole1,
                TSIL_REAL *Mpole2, TSIL_REAL *Gammapole2);

/* ------------------------------------------------------------------- */

int main ()
{
  MSSM_DATA mssm_data;
  int i,j;
  TSIL_REAL Mpole1, Gammapole1, Mpole2, Gammapole2;
  TSIL_REAL rescalefactor;

  fclose(stderr);  /* Suppress irrelevant warning messages. */

  mssm_data.g3 = TSIL_SQRT(0.095L*4.L*PI);
  mssm_data.gluino = TSIL_POW(500.L, 2);
  mssm_data.qq = mssm_data.gluino;

  mssm_data.quark[0] = 0.L; 
  mssm_data.quark[1] = 0.L; 
  mssm_data.quark[2] = 0.L;
  mssm_data.quark[3] = 0.L;
  mssm_data.quark[4] = TSIL_POW(3.L, 2);  /* b quark running squared mass */
  mssm_data.quark[5] = TSIL_POW(165.L, 2);  /* t quark running squared mass */

  /* top squark masses and mixings */
  mssm_data.squark[5][0] = TSIL_POW(300.L, 2); 
  mssm_data.squark[5][1] = TSIL_POW(700.L, 2); 
  mssm_data.Rsq[5][0] = 0.8L; 
  mssm_data.Lsq[5][0] = 0.6L;
  mssm_data.Rsq[5][1] = -0.6L; 
  mssm_data.Lsq[5][1] = 0.8L;
 
  /* u,d,s,c,b squark masses and mixings */
  for (i=0;i<5;i++) { 
    mssm_data.squark[i][0] = TSIL_POW(600.L, 2); 
    mssm_data.squark[i][1] = TSIL_POW(700.L, 2); 
    mssm_data.Rsq[i][0] = 1.L; 
    mssm_data.Lsq[i][0] = 0.L;
    mssm_data.Rsq[i][1] = 0.L; 
    mssm_data.Lsq[i][1] = 1.L;
  }

  /* Following is useful for checks. */
  /* rescalefactor = 1.0L/TSIL_POW(500.L, 2); */
  rescalefactor = 1.0L;
  mssm_data.gluino *= rescalefactor;
  mssm_data.qq *= rescalefactor;
  for (i=0;i<6;i++) { 
    mssm_data.quark[i] *= rescalefactor; 
    mssm_data.squark[i][0] *= rescalefactor; 
    mssm_data.squark[i][1] *= rescalefactor; 
  }

  gluinopole (&mssm_data, &Mpole1, &Gammapole1, &Mpole2, &Gammapole2);

  printf("Tree-level gluino mass:     %.16lf\n", 
         (double) TSIL_SQRT(mssm_data.gluino)); 
  for (i=0; i<6; i++){
   for (j=0; j<2; j++){
    if ((i==0) && (j==0))
      printf("Tree-level squark masses:   %.16lf\n", 
          (double) TSIL_SQRT(mssm_data.squark[i][j])); 
    else 
      printf("                            %.16lf\n", 
          (double) TSIL_SQRT(mssm_data.squark[i][j])); 
  }}
  printf("One-loop gluino pole mass:  %.16lf", (double) Mpole1);
  printf("    Width =  %.16lf\n", (double) Gammapole1);
  printf("Two-loop gluino pole mass:  %.16lf", (double) Mpole2);
  printf("    Width =  %.16lf\n", (double) Gammapole2);
  printf("Delta    gluino pole mass:  %.16lf\n", (double) (Mpole2 - Mpole1));

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (5.1)-(5.3) of hep-ph/0509115 */
int gluinopole (MSSM_DATA *mssm_data,
                TSIL_REAL *mpole1, TSIL_REAL *gammapole1,
                TSIL_REAL *mpole2, TSIL_REAL *gammapole2)
{
  TSIL_REAL    x, y, z, g3, qq, quark[6], squark[6][2];
  TSIL_REAL    mgluino, mquark[6];
  TSIL_COMPLEX Rsq[6][2], Lsq[6][2], Rsqc[6][2], Lsqc[6][2], reLRc[6][2];
  TSIL_DATA    int_data[6][2][2];
  TSIL_COMPLEX deltas1, s1, deltas2, s2;
  TSIL_COMPLEX pi1Vtilde, pi1Stilde, pi1leftover;
  TSIL_COMPLEX pi2tildeMSFFSF, pi2tildeVSFFFS, pi2tildeVFSSFF;
  TSIL_COMPLEX pi2tildeYFSSS, pi2tilde1vector, pi2tilde2vector;
  TSIL_COMPLEX BFS[6][2], BfS[6][2], BpFS[6][2], BpfS[6][2];
  TSIL_COMPLEX Lj, Ljc, Rj, Rjc, Lk, Lkc, Rk, Rkc, Ln, Lnc, Rn, Rnc;
  TSIL_COMPLEX MSFFSF, MSFFSf, MSFfSF, MSFfSf, MSffSF, MSffSf;
  TSIL_COMPLEX VSFFFS, VSFFfS, VSFfFS, VSFffS, VSffFS, VSfffS;
  TSIL_COMPLEX VFSSFF, VFSSff, VfSSFF, VfSSff;
  TSIL_COMPLEX YFSSS, YfSSS;
  TSIL_COMPLEX f1result, f2result, f3result, f4result, f5result, f6result;
  TSIL_COMPLEX BFV, BfV;
  int i,j,k,n;

  /* Read in the parameters from the mssm_data struct. */
  x = mssm_data->gluino;
  g3 = mssm_data->g3;
  qq = mssm_data->qq;
  mgluino = TSIL_SQRT(x);

  for (i=0;i<6;i++) {
    quark[i] = mssm_data->quark[i];
    mquark[i] = TSIL_SQRT (quark[i]);
    for (j=0;j<2;j++) {
      squark[i][j] = mssm_data->squark[i][j];
      Rsq[i][j] = mssm_data->Rsq[i][j];
      Lsq[i][j] = mssm_data->Lsq[i][j];
      Rsqc[i][j] = TSIL_CREAL(Rsq[i][j]) - I*TSIL_CIMAG(Rsq[i][j]);
      Lsqc[i][j] = TSIL_CREAL(Lsq[i][j]) - I*TSIL_CIMAG(Lsq[i][j]);
      reLRc[i][j] = TSIL_CREAL(Lsq[i][j]*Rsqc[i][j]);
  }}

  /* One-loop part starts here. */

  pi1Vtilde = CG * x * (10.L - 6.L * TSIL_LOG(x/qq));

  /* Alternate version on next two lines, just for fun. 
  bFV (x, 0.L, x, qq, 0.L, &BFV, &BfV);
  pi1Vtilde = 2.L * CG * (BFV + x* BfV);
  */

  for (i=0; i<6; i++)
  {
   z = quark[i];
   for (j=0; j<2; j++)
   {
     y = squark[i][j];
     bFS (z,y,x,qq, 
          &(BFS[i][j]), &(BfS[i][j]), &(BpFS[i][j]), &(BpfS[i][j]));
  }}

  pi1Stilde = 0.L;

  for (i=0; i<6; i++){
   for (j=0; j<2; j++){
     pi1Stilde += 4.L * Iq * (BFS[i][j]  
        -2.L * reLRc[i][j] * mgluino * mquark[i] * BfS[i][j]);
  }}

  pi1leftover = 0.L;

  for (i=0; i<6; i++){
   for (j=0; j<2; j++){
    for (n=0; n<6; n++){
     for (k=0; k<2; k++){

      pi1leftover += 16.L * Iq * Iq * 
       (BFS[i][j]-2.L * reLRc[i][j] * mgluino * mquark[i] * BfS[i][j])*
       (BpFS[n][k]-2.L * reLRc[n][k] * mgluino * mquark[n] * BpfS[n][k]);

      pi1leftover += -4.L * Iq * Iq * BFS[i][j] * BFS[n][k]/x;

      pi1leftover += 16.L * Iq * Iq * reLRc[i][j] * reLRc[n][k] * 
                     mquark[i] * mquark[n] * BfS[i][j] * BfS[n][k];
  }}}}

  /* Two-loop part starts here. */

  for (i=0; i<6; i++){
   for (j=0; j<2; j++){
    for (k=0; k<2; k++){
      TSIL_SetParameters (&(int_data[i][j][k]), 
        squark[i][j], quark[i], quark[i], squark[i][k], x, qq);
      TSIL_Evaluate (&int_data[i][j][k], x);
  }}}

  pi2tildeMSFFSF = 0.L;
  pi2tildeVSFFFS = 0.L;
  pi2tildeVFSSFF = 0.L;
  pi2tildeYFSSS  = 0.L;

  for (i=0; i<6; i++){
   for (j=0; j<2; j++){
    for (k=0; k<2; k++){

      Lj = Lsq[i][j];
      Ljc = Lsqc[i][j];
      Rj = Rsq[i][j];
      Rjc = Rsqc[i][j];
      Lk = Lsq[i][k];
      Lkc = Lsqc[i][k];
      Rk = Rsq[i][k];
      Rkc = Rsqc[i][k];

      mSFFSF (&int_data[i][j][k], 
              &MSFFSF, &MSFFSf, &MSFfSF, &MSFfSf, &MSffSF, &MSffSf);

      pi2tildeMSFFSF += 8.L * Iq * (2.L * Cq - CG) * (
        Lj*Rjc*Lkc*Rk*MSFFSF 
        - reLRc[i][k] * mgluino * mquark[i] * MSFfSF
        + 0.5L*(Lj*Ljc*Lk*Lkc + Rj*Rjc*Rk*Rkc) * x * MSFFSf
        + Lj*Ljc*Rk*Rkc * quark[i] * MSffSF
        - reLRc[i][j] * mgluino * mquark[i] * MSFfSf
        + TSIL_CREAL(Lj*Rjc*Lk*Rkc) * x * quark[i] * MSffSf );

      vSFFFS (&int_data[i][j][k], &int_data[i][j][k], 
              &VSFFFS, &VSFFfS, &VSFfFS, &VSFffS, &VSffFS, &VSfffS);

      pi2tildeVSFFFS += 8.L * Iq * Cq * (
        (Lj*Ljc*Lk*Lkc + Rj*Rjc*Rk*Rkc) * VSFFFS
        - 2.L * reLRc[i][j] * mgluino * mquark[i] * VSFfFS
        + (Lj*Rjc*Lkc*Rk + Ljc*Rj*Lk*Rkc) * x * VSFFfS 
        + (Lj*Ljc*Rk*Rkc + Rj*Rjc*Lk*Lkc) * quark[i] * VSffFS
        - 2.L * reLRc[i][k] * mgluino * mquark[i]  * VSFffS
        + (Lj*Rjc*Lk*Rkc + Ljc*Rj*Lkc*Rk) * x * quark[i] * VSfffS );

      vFSSFF (&int_data[i][j][k], &int_data[i][k][j], 
              &VFSSFF, &VFSSff, &VfSSFF, &VfSSff);

      pi2tildeVFSSFF += 16.L * Iq * Cq * (
        -reLRc[i][j] * mgluino * mquark[i]  * VFSSff 
        +TSIL_CREAL(Lj*Rjc*Lk*Rkc + Lj*Ljc*Rk*Rkc) * x * quark[i] * VfSSff );

      if (j==k) {
        pi2tildeVFSSFF += 8.L * Iq * Cq * (VFSSFF
          -2.L * reLRc[i][j] * mgluino * mquark[i]  * VfSSFF);
      }

      for (n=0; n<2; n++){
        Ln = Lsq[i][n];
        Lnc = Lsqc[i][n];
        Rn = Rsq[i][n];
        Rnc = Rsqc[i][n];

        yFSSS(quark[i], squark[i][j], squark[i][k], squark[i][n], x, qq,
             &YFSSS, &YfSSS);

        if (j==k){
         pi2tildeYFSSS += 4.L*Iq*Cq*(Lj*Lnc - Rj*Rnc)*(Ljc*Ln - Rjc*Rn)*YFSSS;
        }
 
        pi2tildeYFSSS += -8.L*Iq*Cq* TSIL_CREAL(Ljc*Rk* 
          (Lj*Lnc - Rj*Rnc)*(Lkc*Ln - Rkc*Rn))*mgluino*mquark[i]*YfSSS;
      }
  }}}    

  pi2tilde1vector = 0.L;

  for (i=0; i<6; i++){
   for (j=0; j<2; j++){
     f1f4(x,quark[i],squark[i][j],qq, &f1result, &f4result);
     pi2tilde1vector += (4.L * Cq - 2.L * CG) * Iq * (f1result 
       - 2.L * reLRc[i][j] * mgluino * mquark[i] * f4result);
  }}

  for (i=0; i<6; i++){
   for (j=0; j<2; j++){
     f2f3f5f6(x,quark[i],squark[i][j],qq, 
              &f2result, &f3result, &f5result, &f6result);
     /* f2 and f5 are not actually used for anything in this example. */
     pi2tilde1vector += 4.L * CG * Iq * (f3result 
       - 2.L * reLRc[i][j] * mgluino * mquark[i] * f6result);
  }}

  pi2tilde2vector = CG * CG * (F1(x,qq) + F2(x,qq) + F3(x,x,qq));
  
  for (i=0;i<6;i++) {
    pi2tilde2vector += CG * Iq * 2.L * F3(x, quark[i], qq);
    for (j=0;j<2;j++) {
      pi2tilde2vector += CG * Iq * F4(x, squark[i][j], qq);
  }}

  /*
  printf("pi1leftover "); TSIL_cprintf(pi1leftover); printf("\n");
  printf("pi2tildeMSFFSF "); TSIL_cprintf(pi2tildeMSFFSF); printf("\n");
  printf("pi2tildeVSFFFS "); TSIL_cprintf(pi2tildeVSFFFS); printf("\n");
  printf("pi2tildeVFSSFF "); TSIL_cprintf(pi2tildeVFSSFF); printf("\n");
  printf("pi2tildeYFSSS "); TSIL_cprintf(pi2tildeYFSSS); printf("\n");
  printf("pi2tilde1vector "); TSIL_cprintf(pi2tilde1vector); printf("\n");
  printf("pi2tilde2vector "); TSIL_cprintf(pi2tilde2vector); printf("\n");
  */

  /* Assemble results to obtain the gluino pole mass */

  deltas1 = kappa*g3*g3*(pi1Vtilde + pi1Stilde);

  deltas2 = kappa*kappa*g3*g3*g3*g3*(pi1leftover + pi2tildeMSFFSF 
                + pi2tildeVSFFFS + pi2tildeVFSSFF + pi2tildeYFSSS
                + pi2tilde1vector + pi2tilde2vector);

  s1 = x + deltas1;
  s2 = s1 + deltas2;

  *mpole1 = TSIL_SQRT(TSIL_CREAL(s1));
  *gammapole1 = -TSIL_CIMAG(s1)/(*mpole1);

  *mpole2 = TSIL_SQRT(TSIL_CREAL(s2));
  *gammapole2 = -TSIL_CIMAG(s2)/(*mpole2);

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (3.3), (3.4), (3.7), (3.8) of hep-ph/0509115 */
int bFS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL s, TSIL_REAL qq,
         TSIL_COMPLEX *BFS, TSIL_COMPLEX *BfS,  
         TSIL_COMPLEX *BpFS, TSIL_COMPLEX *BpfS)
{
  *BFS = ((y-x-s)*TSIL_B(x,y,s,qq) -TSIL_A(x,qq) +TSIL_A(y,qq))/2.L;

  if (((x/(y+s)) < TSIL_TOL) && (TSIL_FABS((y-s)/(y+s)) < TSIL_TOL))
  {
    *BpFS = 0.5L * TSIL_LOG(y/qq) - 1.L;
  }
  else
  {  
    *BpFS = ((y-x-s)*TSIL_dBds(x,y,s,qq) - TSIL_B(x,y,s,qq))/2.L;
  }

  if ((x/(y+s)) > TSIL_TOL)
  {
    *BfS = -TSIL_B(x,y,s,qq);
    *BpfS = -TSIL_dBds(x,y,s,qq);
  }
  else 
  {
    *BfS = 0.L;
    *BpfS = 0.L;    
  }

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (3.5), (3.6) of hep-ph/0509115 */
int bFV (TSIL_REAL x, TSIL_REAL v, TSIL_REAL s, TSIL_REAL qq, TSIL_REAL xi,
         TSIL_COMPLEX *BFV, TSIL_COMPLEX *BfV)
{
  if ((v/(x+s)) > TSIL_TOL)
  {
    *BFV = (v-s-x)*TSIL_B(x, v, s, qq) + TSIL_A(v,qq) - TSIL_A(x,qq)
           + ( (v*(s+x) - (x-s)*(x-s))*TSIL_B(x, v, s, qq) 
               -(xi*v*(s+x) - (x-s)*(x-s))*TSIL_B(x, xi*v, s, qq) 
             + (x-s)*TSIL_A(v,qq) - (x-s)*TSIL_A(xi*v,qq) )/(2.L*v);
  }
  else
  {
    *BFV = -s + xi * (s - TSIL_A(x,qq) - (s+x)*TSIL_B(0, x, s, qq));
  }

  if ((x/(v+s)) > TSIL_TOL)
  {
    *BfV = 3.L*TSIL_B(x, v, s, qq) + xi * TSIL_B(x, xi*v, s, qq);
  }
  else
  {
    *BfV = 0.L;
  }

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.6)-(4.11) of hep-ph/0509115 */
int mSFFSF (TSIL_DATA *input_data, 
            TSIL_COMPLEX *MSFFSF, TSIL_COMPLEX *MSFFSf, TSIL_COMPLEX *MSFfSF, 
            TSIL_COMPLEX *MSFfSf, TSIL_COMPLEX *MSffSF, TSIL_COMPLEX *MSffSf) 
{
  TSIL_COMPLEX M, Uxzuv, Uyuzv, Uzxyv, Uuyxv, Suvx, Svyz, Buy, Bxz;
  TSIL_REAL x, y, z, u, v, s;

  x = input_data->x;
  y = input_data->y;
  z = input_data->z;
  u = input_data->u;
  v = input_data->v;
  s = input_data->s;

  M = TSIL_GetFunction(input_data, "M");
  Uxzuv = TSIL_GetFunction(input_data, "Uxzuv");
  Uyuzv = TSIL_GetFunction(input_data, "Uyuzv");
  Uzxyv = TSIL_GetFunction(input_data, "Uzxyv");
  Uuyxv = TSIL_GetFunction(input_data, "Uuyxv");
  Suvx = TSIL_GetFunction(input_data, "Suxv");
  Svyz = TSIL_GetFunction(input_data, "Svyz");
  Bxz = TSIL_GetFunction(input_data, "Bxz");
  Buy = TSIL_GetFunction(input_data, "Byu");
  
  *MSFFSF = (s*Buy*Bxz - (s*v - u*x + y*z)*M - Suvx + Svyz + y*Uuyxv 
            + z*Uxzuv - u*Uyuzv - x*Uzxyv)/2.L;

  *MSFFSf = (Buy*Bxz + (-s + u - v + x)*M - Uyuzv - Uzxyv)/2.L;

  *MSFfSF = (Buy*Bxz - (v - x + y)*M + Uxzuv - Uyuzv)/2.L;

  *MSFfSf = ((-s + u - y)*M + Uxzuv - Uzxyv)/2.L;
 
  *MSffSF = ((u + x - y - z)*M + Uuyxv + Uxzuv - Uyuzv - Uzxyv)/2.L;

  *MSffSf = -M;

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.12)-(4.14) of hep-ph/0509115 */
int mSSFFS (TSIL_DATA *input_data, 
            TSIL_COMPLEX *MSSFFS, TSIL_COMPLEX *MSSFfS, TSIL_COMPLEX *MSSffS)
{
  TSIL_COMPLEX M, Uyuzv, Uzxyv, Uuyxv, Buy, Bxz;
  TSIL_REAL x, y, z, u, v, s;

  x = input_data->x;
  y = input_data->y;
  z = input_data->z;
  u = input_data->u;
  v = input_data->v;
  s = input_data->s;

  M = TSIL_GetFunction(input_data, "M");
  Uyuzv = TSIL_GetFunction(input_data, "Uyuzv");
  Uzxyv = TSIL_GetFunction(input_data, "Uzxyv");
  Uuyxv = TSIL_GetFunction(input_data, "Uuyxv");
  Bxz = TSIL_GetFunction(input_data, "Bxz");
  Buy = TSIL_GetFunction(input_data, "Byu");

  *MSSFFS = 0.5L*((v-z-u)*M + Uzxyv + Uuyxv - Bxz*Buy);
  *MSSFfS = 0.5L*((x-z-s)*M - Uyuzv + Uuyxv);
  *MSSffS = -M;

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.15)-(4.18) of hep-ph/0509115 */
int yFSSS (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u,
           TSIL_REAL s, TSIL_REAL qq,
           TSIL_COMPLEX *YFSSS, TSIL_COMPLEX *YfSSS)
{
  TSIL_COMPLEX Bxy, Bxz, Bxyp, Ay, Az, Au;

  Bxy = TSIL_B (x,y,s,qq);
  Ay = TSIL_A (y,qq);
  Au = TSIL_A (u,qq);

  if(TSIL_FABS((y-z)/(y+z)) < TSIL_TOL) 
  {
    Bxyp = TSIL_Bp (y,x,s,qq);

    if (((x/(s+y)) < TSIL_TOL) && (TSIL_FABS((y-s)/(y+s)) < TSIL_TOL))
    {
      *YFSSS = Au;
    }
    else
    {        
      *YFSSS = Au*(1.L + Ay/y + Bxy + (y - x - s)*Bxyp)/2.L;
    }

    if (x/(s+y) > TSIL_TOL) 
    {
      *YfSSS = -Au * Bxyp;
    }
    else *YfSSS = 0;
  }  
  else 
  {
    Bxz = TSIL_B (x,z,s,qq);
    Az = TSIL_A (z,qq);
    *YFSSS = Au*(Ay - Az + (y-x-s)*Bxy - (z-x-s)*Bxz)/(2.L*(y - z));
   if (x/(s+y) > TSIL_TOL) {*YfSSS = Au*(Bxz - Bxy)/(y - z);}
     else *YfSSS = 0;
  }

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.19)-(4.22) of hep-ph/0509115 */
int vFSSSS (TSIL_DATA *input_data1, TSIL_DATA *input_data2,
            TSIL_COMPLEX *VFSSSS, TSIL_COMPLEX *VfSSSS)
{
  TSIL_COMPLEX Uzxyv, Vzxyv, Ipxyv, Uzuyv, Ixyv, Iuyv;
  TSIL_REAL x, y, z, u, v, s, qq;
  TSIL_REAL mustbe0;

  qq = input_data1->qq;
  s = input_data1->s;
  x = input_data1->x;
  u = input_data2->x;
  y = input_data1->y;
  v = input_data1->v;
  z = input_data1->z;

  mustbe0 = TSIL_FABS (input_data1->qq - input_data2->qq) 
            + TSIL_FABS (input_data1->z - input_data2->z) 
            + TSIL_FABS (input_data1->y - input_data2->y) 
            + TSIL_FABS (input_data1->v - input_data2->v) 
            + TSIL_FABS (input_data1->s - input_data2->s);
  
  if (mustbe0 > TSIL_TOL) {
    TSIL_Error ("vFSSSS", "Wrong data arguments, exiting.", 1);
  }

  Uzxyv = TSIL_GetFunction(input_data1, "Uzxyv");
  Ixyv = TSIL_I2 (x,y,v,qq);

  if (TSIL_FABS((x-u)/(x+u)) < TSIL_TOL)
  {    
    Vzxyv = TSIL_GetFunction(input_data1, "Vzxyv");
    Ipxyv = TSIL_I2p (x,y,v,qq);

    if (((z/(x+s)) < TSIL_TOL) && (TSIL_FABS((x-s)/(x+s)) < TSIL_TOL))
    {
      *VFSSSS = 0.5L*(-Uzxyv - Ipxyv);
    }
    else
    {
      *VFSSSS = 0.5L*((x-z-s)*Vzxyv - Uzxyv - Ipxyv);
    }

    if ((z/(x+s)) > TSIL_TOL)
    {
      *VfSSSS = -Vzxyv;
    }
    else 
    {
      *VfSSSS = 0.L;
    }
  }
  else
  {
    Uzuyv = TSIL_GetFunction(input_data2, "Uzxyv");
    Iuyv = TSIL_I2 (u,y,v,qq);

    *VFSSSS = 0.5L*(Ixyv - Iuyv+ Uzxyv*(x-z-s) + Uzuyv*(s-u+z))/(u-x);

    *VfSSSS = (Uzuyv - Uzxyv)/(u - x);
  }

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.23)-(4.34) of hep-ph/0509115 */
int vSFFFS (TSIL_DATA *input_data1, TSIL_DATA *input_data2, 
            TSIL_COMPLEX *VSFFFS, TSIL_COMPLEX *VSFFfS, TSIL_COMPLEX *VSFfFS, 
            TSIL_COMPLEX *VSFffS, TSIL_COMPLEX *VSffFS, TSIL_COMPLEX *VSfffS) 
{
  TSIL_COMPLEX Uxzuv, Vxzuv, Izuv, Ipzuv, Bxz, Bpzx;
  TSIL_COMPLEX Uxyuv, Iyuv, Bxy; 
  TSIL_COMPLEX Ax, Ay, Az, Au, Av;
  TSIL_COMPLEX Suxv, uTuxv, xTxvu, vTvux, umvsmxVbarx0uv;
  TSIL_REAL x, y, z, u, v, s, qq;
  TSIL_REAL mustbe0;

  qq = input_data1->qq;
  x = input_data1->x;
  z = input_data1->z;
  y = input_data2->z;
  u = input_data1->u;
  v = input_data1->v;
  s = input_data1->s;

  mustbe0 = TSIL_FABS (input_data1->qq - input_data2->qq) 
            + TSIL_FABS (input_data1->x - input_data2->x) 
            + TSIL_FABS (input_data1->u - input_data2->u) 
            + TSIL_FABS (input_data1->v - input_data2->v) 
            + TSIL_FABS (input_data1->s - input_data2->s);
  
  if (mustbe0 > TSIL_TOL) {
    TSIL_Error ("vSFFFS", "Wrong data arguments, exiting.", 1);
  }

  Uxzuv = TSIL_GetFunction(input_data1, "Uxzuv");
  if (u > TSIL_TOL) uTuxv = u*TSIL_GetFunction(input_data1, "Tuxv");
    else uTuxv = 0.L;
  if (x > TSIL_TOL) xTxvu = x*TSIL_GetFunction(input_data1, "Txuv");
    else xTxvu = 0.L;
  if (v > TSIL_TOL) vTvux = v*TSIL_GetFunction(input_data1, "Tvxu");
    else vTvux = 0.L;
  Suxv = TSIL_GetFunction(input_data1, "Suxv");
  Izuv = TSIL_I2 (z,u,v,qq);
  Bxz = TSIL_B (x,z,s,qq);
  Ax = TSIL_A (x,qq);
  Az = TSIL_A (z,qq);
  Ay = TSIL_A (y,qq);
  Au = TSIL_A (u,qq);
  Av = TSIL_A (v,qq);

  if (TSIL_FABS ((z - y)/(s + z + y)) < TSIL_TOL)
  {
    if ((z/s) > TSIL_TOL) {

      Vxzuv = TSIL_GetFunction(input_data1, "Vxzuv");
      Ipzuv = TSIL_I2p (z,u,v,qq);
      Bpzx = TSIL_Bp (z,x,s,qq);
 
      *VSFFFS = ((Av - Au)*(Bpzx*(s - x + z) + 1.L + Bxz + Az/z) +
                - Suxv + (Ipzuv - Vxzuv*(s - x + z))*(v + z - u) + 
                + Izuv + Uxzuv*(s - u + v - x + 2.L*z))/4.L;

      *VSFFfS = Uxzuv - Vxzuv*z;

      *VSFfFS = ((Av - Au) * Bpzx + Uxzuv + (u - v - z) * Vxzuv)/2.L;

      *VSFffS = (Ipzuv + Uxzuv + Vxzuv*(x - z - s))/2.L;

      *VSffFS = ((Au - Av)*Bxz*(s - x) + 2.L*xTxvu*(v - u) + 
                 + uTuxv*(s + v - u - x) - vTvux*(s + u - v - x) + 
           + (v - u)*(2.L*Suxv - Izuv - Ax + x + u + v -s/4.L) +
           + Av*(u - v + z) - Au*(v - u + z) + z*Ipzuv*(v + z - u) + 
           - ((Au - Av)*z*Bpzx +  Vxzuv*z*(v + z - u))*(s - x + z) +
           + Uxzuv*((s - x)*(u - v) + z*z))/(4.L*z*z);

      *VSfffS = -Vxzuv;
    }
    else 
    {
      if (TSIL_FABS((u-v)/(u+v)) < TSIL_TOL) 
      {
        *VSFFFS = ((s-x) * Uxzuv - Suxv + Izuv)/4.L;
      }
      else if (TSIL_FABS((s-x)/(s+x)) < TSIL_TOL)
      {
        *VSFFFS = ((Au - u)*(Av - v) + u*v + (u+v)*Izuv/2.L)/(2.L*(u-v))
                  + (2.L*(u-v+x)*xTxvu/x -3.L*vTvux + 5.L*uTuxv
                  + Izuv -3.L*Au + Av - Ax + 3.L*u - v
                  + 2.L*(Av - Au)*Ax/x + 0.75L*x)/8.L;
      }
      else 
      {
        umvsmxVbarx0uv = (Av - Au)*Bxz + (s-x)*(v*(Au*Bxz + uTuxv) 
          - u*(Av*Bxz + vTvux))/((u-v)*(u-v)) 
          + (u - v)*x*(Au + Av + Ax + Izuv - 2.L*Suxv - uTuxv - vTvux 
          - u - v - 0.75L*x - 2.L*xTxvu)/((s-x)*(s-x)) 
          + (x*(vTvux - uTuxv) + (Av - Au)*(Ax - x + 2.L*Bxz*x) 
          + (v-u)*(xTxvu - x/4.L))/(s - x) + ((u+v)*(xTxvu + Suxv 
          + (x-s) * Bxz/2.L - Ax/2.L + x + u + v - 5.L*s/8.L) + Au*Av 
          -2.L*(u*Av + v*Au) + u*vTvux + v*uTuxv)/(u - v);

        *VSFFFS = (umvsmxVbarx0uv + (s-x+v-u) * Uxzuv
                 - Suxv + Izuv + 2.L*(Av-Au)*(s-Ax-x*Bxz)/(s-x))/4.L;
      }
      
      *VSFFfS = Uxzuv;
      *VSFfFS = 0.L;
      *VSFffS = 0.L;
      *VSffFS = 0.L;
      *VSfffS = 0.L;
    }
  }
  else
  {
    Uxyuv = TSIL_GetFunction(input_data2, "Uxzuv");
    Iyuv = TSIL_I2 (y,u,v,qq);
    Bxy = TSIL_B (x,y,s,qq);

    *VSFFFS = ((Au - Av)*(Bxy*(s - x + y) - Bxz*(s - x + z) + Ay - Az)  
              + (Uxzuv*(s - x + z) + Izuv)*(v + z - u) 
              - (Uxyuv*(s - x + y) + Iyuv)*(v + y - u) 
              + Suxv*(y - z))/(4.L*(z-y)); 

   *VSFFfS = (z*Uxzuv - y*Uxyuv)/(z-y);

   *VSFfFS = ((Au - Av)*(Bxy - Bxz) +(u-v-y)*Uxyuv -(u-v-z)*Uxzuv)/(2.L*(z-y));

   *VSFffS = (Izuv - Iyuv + Uxzuv*(s-x+z) - Uxyuv*(s-x+y))/(2.L*(z-y));

   if ((y > TSIL_TOL) && (z > TSIL_TOL)) {
     *VSffFS = ((Au - Av)*(Ay + Bxy*(s - x + y)) 
             + (u - v - y)*(Iyuv + Uxyuv*(s - x + y)))/(4.L*y*(z - y)) 
             + ((Au - Av)*(Az + Bxz*(s - x + z)) 
             + (u - v - z)*(Izuv + Uxzuv*(s - x + z)))/(4.L*z*(y - z))  
             + (2.L*(v-u)*xTxvu + (-s - u + v + x)*vTvux 
             + (s - u + v - x)*uTuxv 
             + (u - v)*(Au + Av + Ax + s/4.L - 2.L*Suxv - u - v - x) 
             )/(4.L*y*z);
   } 
   else *VSffFS = 0.L;

   *VSfffS = (Uxzuv - Uxyuv)/(z-y);
  }

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.35)-(4.42) of hep-ph/0509115 */
int vFSSFF (TSIL_DATA *input_data1, TSIL_DATA *input_data2,
            TSIL_COMPLEX *VFSSFF, TSIL_COMPLEX *VFSSff, 
            TSIL_COMPLEX *VfSSFF, TSIL_COMPLEX *VfSSff)
{
  TSIL_COMPLEX Uzxyv, Svyz, Bxz, Bpxz, Ax, Ay, Az, Au, Av;
  TSIL_COMPLEX Vzxyv, Ipxyv, Uzuyv, Ixyv, Iuyv, Buz;
  TSIL_REAL x, y, z, u, v, s, qq;
  TSIL_REAL mustbe0;

  qq = input_data1->qq;
  s = input_data1->s;
  x = input_data1->x;
  u = input_data2->x;
  y = input_data1->y;
  v = input_data1->v;
  z = input_data1->z;

  mustbe0 = TSIL_FABS (input_data1->qq - input_data2->qq) 
            + TSIL_FABS (input_data1->z - input_data2->z) 
            + TSIL_FABS (input_data1->y - input_data2->y) 
            + TSIL_FABS (input_data1->v - input_data2->v) 
            + TSIL_FABS (input_data1->s - input_data2->s);
  
  if (mustbe0 > TSIL_TOL) {
    TSIL_Error ("vFSSFF", "Wrong data arguments, exiting.", 1);
  }

  Uzxyv = TSIL_GetFunction(input_data1, "Uzxyv");
  Svyz = TSIL_GetFunction(input_data1, "Svyz");
  Ixyv = TSIL_I2 (x,y,v,qq);
  Bxz = TSIL_B (x,z,s,qq);
  Ax = TSIL_A (x,qq);
  Ay = TSIL_A (y,qq);
  Av = TSIL_A (v,qq);

  if (TSIL_FABS((x-u)/(x+u)) < TSIL_TOL)
  {    
    Vzxyv = TSIL_GetFunction(input_data1, "Vzxyv");
    Ipxyv = TSIL_I2p (x,y,v,qq);
    Bpxz = TSIL_Bp (x,z,s,qq);

    if (((z/(x+s)) < TSIL_TOL) && (TSIL_FABS((x-s)/(x+s)) < TSIL_TOL))
    {
      *VFSSFF = (Svyz - Ixyv + Uzxyv*(v - x + y) 
                -2.L*(Av + Ay) + (v - x + y)*Ipxyv)/2.L;

      *VFSSff = Uzxyv + Ipxyv;
    }
    else
    {
      *VFSSFF = (Svyz - Ixyv + Uzxyv*(s + v - 2.L*x + y + z) 
                + (Av + Ay)*((s - x + z)*Bpxz - 1.L - Bxz - Ax/x)
                + (v - x + y)*((s - x + z)*Vzxyv + Ipxyv))/2.L;

      *VFSSff = Uzxyv + Ipxyv + (s - x + z)*Vzxyv;
    }


    if ((z/(x+s)) > TSIL_TOL)
    {
      *VfSSFF = Uzxyv + (Av + Ay)*Bpxz + (v - x + y)*Vzxyv;
      *VfSSff = 2.L*Vzxyv;
    }
    else 
    {
      *VfSSFF = 0.L;
      *VfSSff = 0.L;
    }
  }
  else
  {
    Uzuyv = TSIL_GetFunction(input_data2, "Uzxyv");
    Iuyv = TSIL_I2 (u,y,v,qq);
    Buz = TSIL_B (u,z,s,qq);
    Au = TSIL_A (u,qq);

    *VFSSFF = 0.5L*Svyz + (Iuyv*(u - v - y) + Ixyv*(v - x + y) 
              + (Av + Ay)*(Au - Ax + Bxz*(s - x + z) - Buz*(s - u + z))
              + Uzuyv*(-u + v + y)*(s - u + z) - Uzxyv*(v - x + y)*(s - x + z)
              )/(2.L*(x - u));

    *VFSSff = (Ixyv -Iuyv+ Uzxyv*(x - z - s) + Uzuyv*(s - u + z))/(x - u);

    *VfSSFF = ((Av + Ay)*(Bxz - Buz) + Uzxyv*(x - y - v) 
              + Uzuyv*(v + y - u))/(x - u);

    *VfSSff = 2.L*(Uzuyv - Uzxyv)/(x - u);
  }

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.59), (4.65), (4.70) of hep-ph/0509115 */
int f1f4 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq, 
          TSIL_COMPLEX *resultf1, TSIL_COMPLEX *resultf4)
{
  TSIL_COMPLEX Myyzz0, Uz000, S0yz, Tbar0yz, Ty0z, Tz0y, Byzp, Byz, Ay, Az;

  TSIL_Manalytic(y,y,z,z,0,x, &Myyzz0);
  TSIL_Sanalytic(0,y,z,x,qq, &S0yz); 
  TSIL_Tanalytic(z,0,y,x,qq, &Tz0y); 
  TSIL_Tbaranalytic(0,y,z,x,qq, &Tbar0yz);
  Byz = TSIL_B (y,z,x,qq);
  Byzp = TSIL_Bp (z,y,x,qq);
  Ay = TSIL_A (y,qq); 
  Az = TSIL_A (z,qq); 

  if (y < TSIL_TOL)
  {
    TSIL_Uanalytic(z,0,0,0,x,qq, &Uz000);

    if (TSIL_FABS((x-z)/(x+z)) > TSIL_TOL)
    {
      *resultf1 = -(x-z)*(x-z)*Myyzz0 - 5.L*(x-z)*Uz000/2.L
                + S0yz/2.L + x*Tz0y + x*Byz*Byz
                +((x/z - 2.L)*Az + 3.5L*z -1.5L*x)*Byz
                -3.L*Az*Az/z + 5.5L*Az - 13.L*x/8.L - 2.L*z;
    }
    else
    {
      *resultf1 = x*Tz0y + S0yz/2.L + Az * (5.5L - Byz)
         + x * Byz * (2.L + Byz) - 3.L * Az * Az/z - 29.L*z/8.L;
    }

    *resultf4 = 0.L;
  }
  else
  {
    TSIL_Tanalytic(y,0,z,x,qq, &Ty0z); 

    *resultf1 = Byz*Byz*x - Ay*Ay/(2.L*y) + Tz0y*(x + y) 
                + Az*(2.L + 3.L*Byzp*(x + y - z)/2.L - 2.L*y/z) 
                + Az*Byz*(-0.5L + x/z + y/z) 
                - 6.L*Byz*(x + y - z) + 2.L*Tbar0yz*(x + y - z) 
                + Ty0z*(x + 3.L*y - z)/2.L 
                + Ay*Byz*(x + 3.L*y - z)/(2.L*y) + 3.L*Ay*Az/z 
                - 3.L*Az*Az/(2.L*z) + (-33.L*x - 12.L*y + 20.L*z)/8.L 
                + Myyzz0*(-x*x + y*y + 2*x*z - z*z) 
                + Ay*(0.5L + 3.L*Byzp*(x - y + z)) 
                + 2.L*Byzp*((y - z)*(y - z) - x*(y + z)) + S0yz/2.L;

    *resultf4 = -6.L - 14.L*Byz + Byz*Byz + 4.L*Tbar0yz + 2.L*Ty0z 
                + 2.L*Tz0y + 5.L*Ay*Byz/y + Az*(3.L*Byzp - 2.L/z) 
                + 2.L*Az*Byz/z + 3.L*Ay*Az/(y*z) - 2.L*Byzp*(x - y + z) 
                + Ay*(-3.L/y - 3.L*Byzp*(-x + y + z)/y)
                + 2.L*Myyzz0*(-x + y + z); 
  }

  return 0;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.61), (4.63), (4.67), (4.69), (4.71), (4.72) of hep-ph/0509115 */
int f2f3f5f6 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq, 
              TSIL_COMPLEX *resultf2, TSIL_COMPLEX *resultf3,
              TSIL_COMPLEX *resultf5, TSIL_COMPLEX *resultf6)
{
  TSIL_DATA result;
  TSIL_COMPLEX M0yxzy, M0zxyz, Ux0yy, Ux0zz, Uyzxy, Uzyxz, Uz000;
  TSIL_COMPLEX S0yz, Sxyy, Sxzz, Ty0z, Tz0y, Tbar0yz;
  TSIL_COMPLEX Byz, Byzp, Ax, Ay, Az, Ixyz, Ipxyz;

  TSIL_Sanalytic(0,y,z,x,qq, &S0yz); 
  TSIL_Tanalytic(z,0,y,x,qq, &Tz0y); 
  TSIL_Tanalytic(y,0,z,x,qq, &Ty0z); 
  TSIL_Tbaranalytic(0,y,z,x,qq, &Tbar0yz);
  Byz = TSIL_B (y,z,x,qq);
  Byzp = TSIL_Bp (z,y,x,qq);
  Ax = TSIL_A (x,qq); 
  Ay = TSIL_A (y,qq); 
  Az = TSIL_A (z,qq); 
  Ixyz = TSIL_I2 (x,y,z,qq);
  Ipxyz = TSIL_I2p (x,y,z,qq);

  TSIL_SetParameters (&result, 0,y,x,z,y,qq);
  TSIL_Evaluate (&result, x);
  M0yxzy = TSIL_GetFunction(&result, "M");
  Ux0yy = TSIL_GetFunction(&result, "Uzxyv");
  Uyzxy = TSIL_GetFunction(&result, "Uyuzv");
  Sxyy = TSIL_GetFunction(&result, "Svyz");

  TSIL_SetParameters (&result, 0,z,x,y,z,qq);
  TSIL_Evaluate (&result, x);
  M0zxyz = TSIL_GetFunction(&result, "M");
  Ux0zz = TSIL_GetFunction(&result, "Uzxyv");
  Uzyxz = TSIL_GetFunction(&result, "Uyuzv");
  Sxzz = TSIL_GetFunction(&result, "Svyz");

  if (y < TSIL_TOL)
  {
    TSIL_Uanalytic(z,0,0,0,x,qq, &Uz000);
    if (TSIL_FABS((x-z)/(x+z)) > TSIL_TOL)
    {
      *resultf2 = Ax - 13.L*Az/4.L + Sxzz/2.L + 3.L*x/16.L + 3.L*Az*Az/(2.L*z) 
                + 2*z + Ux0zz*(-x + z) + Az*(1.L - x/(2.L*z))*Byz 
                + (x - 5.L*z)*Byz/4.L + (x-z)*(x-z)*M0yxzy 
                + (-x*x + z*z)*M0zxyz - Sxyy/2.L + 5.L*S0yz/4.L 
                + (-x/2.L + z)*Tz0y - x*Uyzxy + (x - z)*Ux0yy 
                + (-x + z)*Uz000/4.L;
  
      *resultf3 = -6.L*Ax - Sxzz/2.L + 4.L*Ax*Az/x + Ux0zz*(x - z) 
                - 1.5L*Az*Az/z + Az*(-9.L/4.L - z/(2.L*x)) 
                + (87.L*x - 2.L*z + 8.L*z*z/x)/16.L 
                + Az*(-1.L + x/(2.L*z))*Byz + (-7.L*x + 3.L*z)*Byz/4.L 
                + Ax*(1.5L + 2.5L*z/x)*Byz + (4.5L - z/(2.L*x))*Ixyz 
                + (x-z)*(x-z)*M0yxzy + (x*x - z*z)*M0zxyz - Sxyy/2.L 
                + (-5.L/4.L + z/x)*S0yz + (x*x - x*z + z*z)*Tz0y/(2.L*x) 
                - x*Uyzxy + (x - z)*Ux0yy - 9.L*(x - z)*Uz000/4.L;
    }
    else 
    {
      *resultf2 = -9.L*Az/4.L + Sxzz/2.L + 35.L*x/16.L + 3.L*Az*Az/(2.L*z)
         + 0.5L*Az*Byz - z*Byz + 3.L*S0yz/4.L + 0.5L*z*Tz0y - x*Uyzxy;

      *resultf3 = -8.75L*Ax - Sxzz/2.L + 2.5L*Ax*Az/x + 93.L*x/16.L
         + 3.5L*Az*Byz -x*Byz + 4.L*Ixyz- 0.75L*Sxyy + x*Tz0y/2.L - x*Uyzxy;
    }

    *resultf5 = 0.L;
    *resultf6 = 0.L; 
  }
  else
  {
    *resultf2 = Ax + 1.25L*S0yz - Sxyy/2.L + Sxzz/2.L + Uyzxy*(-x - y) 
                - Ay*Ay/(4.L*y) - Uzyxz*y + M0zxyz*(-x*x + (y-z)*(y-z)) 
                + Az*(-1.5L - 0.75L*Byzp*(x + y - z) - y/z) 
                + Ux0yy*(x + y - z) + M0yxzy*(x + y - z)*(x + y - z) 
                + Ty0z*(x + 3.L*y - z)/4.L + Ay*Byz*(x + 3.L*y - z)/(4.L*y) 
                + 1.5L*Ay*Az/z + 0.75L*Az*Az/z 
                + Az*Byz*(-2.L*x - 2.L*y + z)/(4.L*z) 
                + Ux0zz*(-x - y + z) + Byzp*(y - z)*(-x + y + z) 
                + Tz0y*(-x - y + 2.L*z)/2.L 
                + Ay*(0.25L + 3.L*Byzp*(x - y + z)/2.L) 
                + (-x + 20.L*(y + z))/16.L;

    *resultf3 = -Sxyy/2.L - Sxzz/2.L - Ax*Ay/x + Uyzxy*(-x - y) 
                - Ay*Ay/(4.L*y) + Uzyxz*y + M0zxyz*(x*x - (y - z)*(y - z)) 
                + Ax*Az*(1.L/x - 3.L/z) 
                - Ax*Byz*(3.L*x + 5.L*y - 5.L*z)/(2.L*x) 
                - 6.L*Byz*(x + y - z) + 3.L*Ipxyz*(x + y - z) 
                + 2.L*Tbar0yz*(x + y - z) + Ux0yy*(x + y - z) 
                + Ux0zz*(x + y - z) + M0yxzy*(x + y - z)*(x + y - z) 
                + Ixyz*(3.L*x + y - z)/(2.L*x) 
                + Ay*Byz*(x + 3.L*y - z)/(4.L*y) 
                + 1.5L*Ay*Az/z - 3.L*Az*Az/(4.L*z) 
                + Az*Byz*(2.L*x + 2.L*y - z)/(4.L*z) 
                + 3.L*Ax*Byzp*(-x + y + z) 
                + Tz0y*(x*x + x*y - x*z - y*z + z*z)/(2.L*x) 
                + (3.L*x*x - 8.L*y*y + 6.L*x*(y - 5.L*z) + 8.L*z*z)/(16.L*x) 
                + S0yz*(-5.L/4.L + (-y + z)/x) 
                + Ay*((x + 2.L*y - 2.L*z)/(4.L*x) 
                + 3.L*Byzp*(x - y + z)/2.L) 
                + Byzp*(-(y-z)*(y-z) + x*(y + z)) 
                + (Ty0z*(x*x + 2.L*y*(-y + z) - x*(3.L*y + z)))/(4.L*x) 
                + Az*((3.L*Byzp*(x + y - z))/4.L 
                + ((y - z)*z + 2.L*x*(y + z))/(2.L*x*z));

    *resultf5 = 1.L - Byz + Ty0z - Tz0y + 2.L*Ux0yy - 2.L*Ux0zz 
                - 2.L*Uyzxy - Uzyxz + 2.5L*Ay*Byz/y 
                - Az*(1.5L*Byzp + 1.L/z) + 2.L*M0yxzy*(x + y - z) 
                - Az*Byz/z + 1.5L*Ay*Az/(y*z) - 2.L*M0zxyz*(x - y + z) 
                + Byzp*(-x + y + 3.L*z) 
                + Ay*(-1.5L/y - 1.5L*Byzp*(-x + y + z)/y);

    *resultf6 = -11.L*Byz + 6.L*Ipxyz + 4.L*Tbar0yz + 2.L*Ux0yy 
                + 2.L*Ux0zz - 2.L*Uyzxy + Uzyxz - Ax*Byz/x + 2.L*Ixyz/x 
                - 4.L*S0yz/x + 2.5L*Ay*Byz/y + Ty0z*(1.L - 2.L*y/x) 
                + Az*(1.5L*Byzp + 2.L/x + 1.L/z) 
                + 2.L*M0yxzy*(x + y - z) + Az*Byz/z 
                - 3.L*Ax*Az/(x*z) + 1.5L*Ay*Az/(y*z) + Byzp*(x - y + z) 
                + 2.L*M0zxyz*(x - y + z) + Tz0y*(1.L - 2.L*z/x) 
                + Ax*(-3.L/x - 3.L*Byzp*(x - y + z)/x) 
                - (x + 4.L*(y + z))/(2.L*x) + Ay*(2.L/x - 1.5L/y 
                - 1.5L*Byzp*(-x + y + z)/y );
  }
  
  return 0;
}

/* ------------------------------------------------------------------- */

/* eq. (4.85) of hep-ph/0509115 */
TSIL_REAL F1 (TSIL_REAL x, TSIL_REAL qq)
{
  TSIL_REAL lnbarx, result;
  
  lnbarx = TSIL_LOG(x/qq);
  result = x * (41.L/4.L + 10.L * PI2 - 27.L * lnbarx 
           +18.L * lnbarx * lnbarx + 24.L * Zeta3 - 16.L * PI2 * ln2);
  return result;
}

/* ------------------------------------------------------------------- */

/* eq. (4.87) of hep-ph/0509115 */
TSIL_REAL F2 (TSIL_REAL x, TSIL_REAL qq)
{
  TSIL_REAL lnbarx, result;
  
  lnbarx = TSIL_LOG(x/qq);
  result = x * (1093.L/12.L - 8.L * PI2/3.L - (179.L/3.L) * lnbarx 
           +11.L * lnbarx * lnbarx - 12.L * Zeta3 + 8.L * PI2 * ln2);
  return result;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.89), (4.93) of hep-ph/0509115 */
TSIL_REAL F3 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL qq)
{
  TSIL_REAL lnbarx, lnxoy, result;

  lnbarx = TSIL_LOG(x/qq);

  if (TSIL_FABS(y/x) < TSIL_TOL) {
    result = x * ((26.L/3.L) * lnbarx - 2.L * lnbarx * lnbarx 
                  -(37.L/3.L) - 4.L * PI2/3.L);
    return result;
  }

  lnxoy = TSIL_LOG(x/y);
  result = -(37.L/3.L) * x - 12.L * y + (26.L/3.L) * x * lnbarx
           -2.L * x * lnbarx * lnbarx + 4.L * y * lnxoy
           -2.L * lnxoy * lnxoy * y*y/x 
           +8.L * (x+y) * f(TSIL_SQRT(y/x))
           -4.L * (x + y*y/x) * (PI2/6.L + TSIL_Dilog(1.L - y/x));

  return result;
}

/* ------------------------------------------------------------------- */

/* eqs. (4.91), (4.96) of hep-ph/0509115 */
TSIL_REAL F4 (TSIL_REAL x, TSIL_REAL y, TSIL_REAL qq)
{
  TSIL_REAL lnbarx, lnxoy, result;

  lnbarx = TSIL_LOG(x/qq);

  if (TSIL_FABS(y/x) < TSIL_TOL) {
    result = x * ((19.L/3.L) * lnbarx - lnbarx * lnbarx 
                  -(125.L/12.L) - 2.L * PI2/3.L);
    return result;
  }

  lnxoy = TSIL_LOG(x/y);
  result = -(125.L/12.L) * x + 14.L * y + (19.L/3.L) * x * lnbarx
           -x * lnbarx * lnbarx - 6.L * y * lnxoy + lnxoy * lnxoy * y*y/x 
           +8.L * (x-y) * f(TSIL_SQRT(y/x))
           +2.L * (y*y/x - x) * (PI2/6.L + TSIL_Dilog(1.L - y/x));

  return result;
}

/* ------------------------------------------------------------------- */

/* eq. (4.92) of hep-ph/0509115 */
TSIL_REAL f (TSIL_REAL r)
{
  TSIL_REAL temp, result;

  temp = (1.L - r)/(1.L + r);
  result = r * (TSIL_Dilog (temp) - TSIL_Dilog (-temp) + PI2/4.L);
  return result;
}

