
#include "rge.h"
#include "sarah_nmssm_susy.h"

class NMSSMSoftPars : public NMSSMSusyPars {
public:
   NMSSMSoftPars();
   NMSSMSoftPars(const NMSSMSusyPars& , DoubleMatrix TYu_, DoubleMatrix TYd_, DoubleMatrix TYe_, double TLambdax_,
   double TKappa_, DoubleMatrix mq2_, DoubleMatrix ml2_, double mHd2_, double
   mHu2_, DoubleMatrix md2_, DoubleMatrix mu2_, DoubleMatrix me2_, double ms2_,
   double MassB_, double MassWB_, double MassG_
);
   virtual ~NMSSMSoftPars() {}
   virtual DoubleVector beta() const;
   virtual const DoubleVector display() const;
   virtual void set(const DoubleVector&);

   NMSSMSoftPars calcBeta() const;

   void set_TYu(const DoubleMatrix& TYu_) { TYu = TYu_; }
   void set_TYd(const DoubleMatrix& TYd_) { TYd = TYd_; }
   void set_TYe(const DoubleMatrix& TYe_) { TYe = TYe_; }
   void set_TLambdax(double TLambdax_) { TLambdax = TLambdax_; }
   void set_TKappa(double TKappa_) { TKappa = TKappa_; }
   void set_mq2(const DoubleMatrix& mq2_) { mq2 = mq2_; }
   void set_ml2(const DoubleMatrix& ml2_) { ml2 = ml2_; }
   void set_mHd2(double mHd2_) { mHd2 = mHd2_; }
   void set_mHu2(double mHu2_) { mHu2 = mHu2_; }
   void set_md2(const DoubleMatrix& md2_) { md2 = md2_; }
   void set_mu2(const DoubleMatrix& mu2_) { mu2 = mu2_; }
   void set_me2(const DoubleMatrix& me2_) { me2 = me2_; }
   void set_ms2(double ms2_) { ms2 = ms2_; }
   void set_MassB(double MassB_) { MassB = MassB_; }
   void set_MassWB(double MassWB_) { MassWB = MassWB_; }
   void set_MassG(double MassG_) { MassG = MassG_; }

   const DoubleMatrix& get_TYu() const { return TYu; }
   const DoubleMatrix& get_TYd() const { return TYd; }
   const DoubleMatrix& get_TYe() const { return TYe; }
   double get_TLambdax() const { return TLambdax; }
   double get_TKappa() const { return TKappa; }
   const DoubleMatrix& get_mq2() const { return mq2; }
   const DoubleMatrix& get_ml2() const { return ml2; }
   double get_mHd2() const { return mHd2; }
   double get_mHu2() const { return mHu2; }
   const DoubleMatrix& get_md2() const { return md2; }
   const DoubleMatrix& get_mu2() const { return mu2; }
   const DoubleMatrix& get_me2() const { return me2; }
   double get_ms2() const { return ms2; }
   double get_MassB() const { return MassB; }
   double get_MassWB() const { return MassWB; }
   double get_MassG() const { return MassG; }


protected:
   DoubleMatrix TYu;
   DoubleMatrix TYd;
   DoubleMatrix TYe;
   double TLambdax;
   double TKappa;
   DoubleMatrix mq2;
   DoubleMatrix ml2;
   double mHd2;
   double mHu2;
   DoubleMatrix md2;
   DoubleMatrix mu2;
   DoubleMatrix me2;
   double ms2;
   double MassB;
   double MassWB;
   double MassG;


private:
   static const int numberOfParameters = 115;
};
