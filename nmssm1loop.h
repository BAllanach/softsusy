/* extern "C" int gettadS(double g,double gp, double ll, double kk, double ht, double hb, double htau, double v1, double v2, double xx, double Ak, double Al, double At, double Ab, double Atau, double Q, double tadS); */

/* extern "C" int getPiSS(double g,double gp, double ll, double kk, double ht, double hb, double htau, double v1, double v2, double xx, double Ak, double Al, double At, double Ab, double Atau, double p, double Q, double piSS); */

/* extern "C" int getPiPP(double *g,double *gp, double *ll, double *kk, double *ht, double *hb, double *htau, double *v1, double *v2, double *xx, double *Ak, double *Al, double *At, double *Ab, double *Atau, double *p, double *Q, DoubleMatrix piPP); */

/* extern "C" int getpipp_(double *g,double *gp, double *ll, double *kk, double *ht, double *hb, double *htau, double *v1, double *v2, double *xx, double *Ak, double *Al, double *At, double *Ab, double *Atau, double *p, double *Q, double piPP[][4][4]); */


/* extern "C" int getpipp_(double *g,double *gp, double *ll, double *kk, double *ht, double *hb, double *htau, double *v1, double *v2, double *xx, double *Ak, double *Al, double *At, double *Ab, double *Atau, double *p, double *Q, double *pipp11, double *pipp12, double *pipp13, double *pipp22, double *pipp23, double *pipp33); */


 extern "C" int treemasses_(double *g,double *gp,double *ll,double *kk,double *ht,double *hb,double *htau,double *v1,double *v2,double *xx,double *M1,double *M2, double *Ak,double *Al,double *At,double *Ab,double *Atau,double *mQ3,double *mtr,double *mbr,double *mQ,double *mur,double *mdr,double *mL3,double *mtaur,double *mL,double *mer,double *Q,bool *errmass);

 extern "C" int getpiss_(double *g,double *gp, double *ll, double *kk, double *ht, double *hb, double *htau, double *v1, double *v2, double *xx, double *Ak, double *Al, double *At, double *Ab, double *Atau, double *p, double *Q, double piSS[][3][3]); 


extern "C" int getpipp_(double *g,double *gp, double *ll, double *kk, double *ht, double *hb, double *htau, double *v1, double *v2, double *xx, double *Ak, double *Al, double *At, double *Ab, double *Atau, double *p, double *Q, double piPP11[][3][3]);


extern "C" int getpizz_(double *g, double *gp, double *ht, double *hb, double *htau, double *v1, double *v2, double *p, double *Q, double *piZZ);

extern "C" int getpiww_(double *g, double *gp, double *ht, double *hb, double *htau, double *v1, double *v2, double *p, double *Q, double *piWW);

/* extern "C" int treemasses(double g, double gp, double ll, double kk, double ht, double hb, double htau, double v1, double v2, double xx, double M1, double M2,  double Ak, double Al, double At, double Ab, double Atau, double mQ3, double mtr, double mbr, double mQ, double mur, double mdr, double mL3, double mtaur, double mL, double mer,  double Q, bool errmass); */

/* extern "C" int  tree_charginos( double g, double ll, double v1, double v2, double M2, double xx, double xmc, double u, double v); */

/* extern "C" int tree_neutralinos( double g, double gp, double ll, double kk, double v1, double v2, double xx, double M1, double M2, double xmn, double Z); */


/* extern "C" int tree_higgses( double g, double gp, double ll, double kk, double v1, double v2, double xx, double Ak, double Al,  double mss, double maa, double mhc, double RS, double RP, double RC, bool errhiggs); */



/* extern "C" int tree_sfermions( double g, double gp, double ll, double ht, double hb, double htau, double v1, double v2, double xx, double At, double Ab, double Atau,  double mQ3, double mtr, double mbr, double mQ, double mur, double mdr, double mL3, double mtaur, double mL, double mer, double mstop, double msbot,  double mstau, double Rt, double Rb, double Rtau, double msup, double msdown, double msel, double msnutau, double msnue, bool errsfer); */



/* extern "C" int diagsfe( double n, double hf, double g, double gp, double ll, double v1, double v2, double xx, double Af, double mL, double mR, double mass, double R,bool error); */

/* extern "C" int coupl_s_sf( double g, double gp, double ll, double ht, double hb, double htau, double v1, double v2, double xx, double At, double Ab, double Atau,  double Rt, double Rb, double Rtau, double lsstt, double lssbb, double lsstata, double lssntnt, double lssuu, double lssdd, double lssee, double lssnn, double lstt, double lsbb, double lstata, double lsntnt, double lsuu, double lsdd, double lsee, double lsnn); */

/* extern "C" int coupl_s_hh( double g, double gp, double ll, double kk, double v1, double v2, double xx, double Al, double Ak, double RS, double RP, double lsshhl, double lssaa, double lshh, double lsaa, double lsscc, double lscc); */


/* extern "C" int coupl_s_ino( double g, double gp, double ll, double kk, double NN, double UU, double VV, double lsnene, double lschch); */


/*    extern "C" int coupl_p_sf( double g, double gp, double ll, double ht, double hb, double htau, double v1, double v2, double xx, double At, double Ab, double Atau, double  Rt, double Rb, double Rtau, double lpptt, double lppbb, double lpptata, double lppntnt, double lppuu, double lppdd, double lppee, double lppnn, double lptt, double lpbb, double lptata, double lpntnt, double lpuu, double lpdd, double lpee, double lpnn); */


/* extern "C" int coupl_p_hh( double g, double gp, double ll, double kk, double v1, double v2, double xx, double Al, double Ak, double RS, double RP,  double lpphh, double lppaa, double lpah, double lppcc, double lpcc); */

/* extern "C" int coupl_p_ino( double g, double gp, double ll, double kk, double NN, double UU, double VV, double lpnene, double lpchch); */



/*    extern "C" int coupl_Z_ino( double g, double gp, double NN, double UU, double VV, double lznene, double azchch, double bzchch); */



/* extern "C" int coupl_W_ino( double g, double NN, double UU, double VV, double awnech, double bwnech); */




/*    extern "C" int jacobi(double a, int n,int np, double d, double v, int nrot); */


