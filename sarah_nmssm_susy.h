
#include "rge.h"

class NMSSMSusyPars : public RGE {
public:
   NMSSMSusyPars();
   NMSSMSusyPars(double scale_, double loops_, double thresholds_, DoubleMatrix Yu_, DoubleMatrix Yd_, DoubleMatrix Ye_, double Lambdax_,
   double Kappa_, double g1_, double g2_, double g3_, double vd_, double vu_,
   double vS_
);
   virtual ~NMSSMSusyPars() {}
   virtual DoubleVector beta() const;
   virtual const DoubleVector display() const;
   virtual void set(const DoubleVector&);

   NMSSMSusyPars calcBeta() const;

   void set_Yu(const DoubleMatrix& Yu_) { Yu = Yu_; }
   void set_Yd(const DoubleMatrix& Yd_) { Yd = Yd_; }
   void set_Ye(const DoubleMatrix& Ye_) { Ye = Ye_; }
   void set_Lambdax(double Lambdax_) { Lambdax = Lambdax_; }
   void set_Kappa(double Kappa_) { Kappa = Kappa_; }
   void set_g1(double g1_) { g1 = g1_; }
   void set_g2(double g2_) { g2 = g2_; }
   void set_g3(double g3_) { g3 = g3_; }
   void set_vd(double vd_) { vd = vd_; }
   void set_vu(double vu_) { vu = vu_; }
   void set_vS(double vS_) { vS = vS_; }

   const DoubleMatrix& get_Yu() const { return Yu; }
   const DoubleMatrix& get_Yd() const { return Yd; }
   const DoubleMatrix& get_Ye() const { return Ye; }
   double get_Lambdax() const { return Lambdax; }
   double get_Kappa() const { return Kappa; }
   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_vd() const { return vd; }
   double get_vu() const { return vu; }
   double get_vS() const { return vS; }

   DoubleMatrix get_SqSq() const;
   DoubleMatrix get_SlSl() const;
   double get_SHdSHd() const;
   double get_SHuSHu() const;
   DoubleMatrix get_SdRSdR() const;
   DoubleMatrix get_SuRSuR() const;
   DoubleMatrix get_SeRSeR() const;
   double get_SsRSsR() const;


protected:
   DoubleMatrix Yu;
   DoubleMatrix Yd;
   DoubleMatrix Ye;
   double Lambdax;
   double Kappa;
   double g1;
   double g2;
   double g3;
   double vd;
   double vu;
   double vS;


private:
   static const int numberOfParameters = 35;
};
