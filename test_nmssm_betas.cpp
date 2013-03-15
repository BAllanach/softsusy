
#include "sarah_nmssm_softpars.h"
#include "nmssmsoftsusy.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

const double max_dev = 1.0e-12;

template <typename T>
bool is_zero(T a)
{
   return std::fabs(a) < std::numeric_limits<T>::epsilon();
}

template <typename T>
bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   return std::fabs(a - b) < prec;
}

bool is_equal(const DoubleMatrix& a, const DoubleMatrix& b, double max_dev)
{
   if (a.displayRows() != b.displayRows()) {
      std::cout << "matrices have different number of rows\n";
      return false;
   }
   if (a.displayCols() != b.displayCols()) {
      std::cout << "matrices have different number of columns\n";
      return false;
   }
   for (int i = 1; i <= a.displayRows(); ++i) {
      for (int l = 1; l <= a.displayCols(); ++l) {
         if (!is_equal(a(i,l), b(i,l), max_dev))
            return false;
      }
   }
   return true;
}

template <typename T>
void check_equality(T a, T b, const std::string& testMsg, T max_dev)
{
   if (!is_equal(a, b, max_dev))
      std::cout << "test failed: " << testMsg << ": "
                << a << " != " << b
                << " (diff.: " << (a-b) << ", rel. diff.: "
                << 100. * (a-b)/a << "%)\n";
}

void check_equality(int a, int b, const std::string& testMsg, double)
{
   if (a != b)
      std::cout << "test failed: " << testMsg << ": "
                << a << " != " << b << ")\n";
}

void check_equality(const DoubleMatrix& a, const DoubleMatrix& b,
                    const std::string& testMsg, double max_dev)
{
   if (!is_equal(a, b, max_dev)) {

      DoubleMatrix diff(a-b);
      int i,l;
      double max_diff = diff.max(i,l);
      for (int i = 1; i <= a.displayRows(); ++i)
         for (int l = 1; l <= a.displayCols(); ++l) {
            if (!is_zero(a(i,l)))
               diff(i,l) = std::fabs(diff(i,l) / a(i,l));
            else
               diff(i,l) = 0.0;
         }
      double max_rel_diff = 100. * diff.max(i,l);

      std::cout << "test failed: " << testMsg
                << " (max. diff.: " << max_diff
                << ", max. rel. diff.: " << max_rel_diff << "% in element "
                << i << "," << l << ")\n";
   }
}

void check_condition(bool condition, const std::string& testMsg)
{
   if (!condition)
      std::cout << "test failed: " << testMsg << "\n";
}

template <typename T>
void check_greater_than(T a, T b, const std::string& testMsg)
{
   if (!(a > b))
      std::cout << "test failed: " << testMsg << ": " << a << " > "
                << b << "\n";
}

#define TEST_EQUALITY(a, b) check_equality(a, b, #a " == " #b, max_dev)
#define TEST_CLOSE(a, b, dev) check_equality(a, b, #a " == " #b, dev)
#define TEST_GREATER(a, b) check_greater_than(a, b, #a " > " #b)
#define TEST(condition) check_condition(condition, #condition);

void compare_anomalous_dimensions(const SoftParsNmssm& a, const NMSSMSoftPars& b)
{
  DoubleMatrix gEE(3,3),gLL(3,3),gQQ(3,3),gDD(3,3),gUU(3,3);
  double gH1H1 = 0.0, gH2H2 = 0.0, gSS = 0.0;
  DoubleVector dg(1,3);
  nmsBrevity brevity;
  a.anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg, gH1H1, gH2H2, gSS, brevity);

  TEST_EQUALITY(a.displayLoops(), b.displayLoops());
  TEST_EQUALITY(a.displayMu(), b.displayMu());
  TEST_EQUALITY(a.displayThresholds(), b.displayThresholds());

  TEST_EQUALITY(gEE, b.get_SeRSeR());
  TEST_EQUALITY(gLL, b.get_SlSl());
  TEST_EQUALITY(gQQ, b.get_SqSq());
  TEST_EQUALITY(gUU, b.get_SuRSuR());
  TEST_EQUALITY(gDD, b.get_SdRSdR());
  TEST_EQUALITY(gH1H1, b.get_SHdSHd());
  TEST_EQUALITY(gH2H2, b.get_SHuSHu());
  TEST_EQUALITY(gSS  , b.get_SsRSsR());
}

void test_parameter_equality(const SoftParsNmssm& a, const NMSSMSoftPars& b)
{
   TEST_EQUALITY(a.displayLoops(), b.displayLoops());
   TEST_EQUALITY(a.displayMu(), b.displayMu());
   TEST_EQUALITY(a.displayThresholds(), b.displayThresholds());

   TEST_EQUALITY(a.displayGaugeCoupling(1), b.get_g1());
   TEST_EQUALITY(a.displayGaugeCoupling(2), b.get_g2());
   TEST_EQUALITY(a.displayGaugeCoupling(3), b.get_g3());

   TEST_EQUALITY(a.displayYukawaMatrix(YU), b.get_Yu());
   TEST_EQUALITY(a.displayYukawaMatrix(YD), b.get_Yd());
   TEST_EQUALITY(a.displayYukawaMatrix(YE), b.get_Ye());

   TEST_EQUALITY(a.displayLambda(), b.get_Lambdax());
   TEST_EQUALITY(a.displayKappa(), b.get_Kappa());

   TEST_EQUALITY(a.displayGaugino(1), b.get_MassB());
   TEST_EQUALITY(a.displayGaugino(2), b.get_MassWB());
   TEST_EQUALITY(a.displayGaugino(3), b.get_MassG());

   TEST_EQUALITY(a.displayMh1Squared(), b.get_mHd2());
   TEST_EQUALITY(a.displayMh2Squared(), b.get_mHu2());
   TEST_EQUALITY(a.displaySoftMassSquared(mQl), b.get_mq2());
   TEST_EQUALITY(a.displaySoftMassSquared(mUr), b.get_mu2());
   TEST_EQUALITY(a.displaySoftMassSquared(mDr), b.get_md2());
   TEST_EQUALITY(a.displaySoftMassSquared(mLl), b.get_ml2());
   TEST_EQUALITY(a.displaySoftMassSquared(mEr), b.get_me2());

   TEST_EQUALITY(a.displayTrilinear(UA), b.get_TYu());
   TEST_EQUALITY(a.displayTrilinear(DA), b.get_TYd());
   TEST_EQUALITY(a.displayTrilinear(EA), b.get_TYe());

   const double tanBeta = b.get_vu() / b.get_vd();
   const double vev = sqrt(sqr(b.get_vu()) + sqr(b.get_vd()));
   TEST_EQUALITY(a.displayTanb(), tanBeta);
   TEST_EQUALITY(a.displayHvev(), vev);
   TEST_EQUALITY(a.displaySvev(), b.get_vS());
}

void test_beta_function_equality(const SoftParsNmssm& a, const NMSSMSoftPars& b)
{
   SoftParsNmssm beta_a(a.beta2());
   NMSSMSoftPars beta_b(b.calcBeta());

   TEST_EQUALITY(beta_a.displayLoops(), beta_b.displayLoops());
   TEST_EQUALITY(beta_a.displayMu(), beta_b.displayMu());
   TEST_EQUALITY(beta_a.displayThresholds(), beta_b.displayThresholds());

   TEST_EQUALITY(beta_a.displayGaugeCoupling(1), beta_b.get_g1());
   TEST_EQUALITY(beta_a.displayGaugeCoupling(2), beta_b.get_g2());
   TEST_EQUALITY(beta_a.displayGaugeCoupling(3), beta_b.get_g3());

   TEST_EQUALITY(beta_a.displayYukawaMatrix(YU), beta_b.get_Yu());
   TEST_EQUALITY(beta_a.displayYukawaMatrix(YD), beta_b.get_Yd());
   TEST_EQUALITY(beta_a.displayYukawaMatrix(YE), beta_b.get_Ye());

   TEST_EQUALITY(beta_a.displayLambda(), beta_b.get_Lambdax());
   TEST_EQUALITY(beta_a.displayKappa(), beta_b.get_Kappa());

   TEST_EQUALITY(beta_a.displayGaugino(1), beta_b.get_MassB());
   TEST_EQUALITY(beta_a.displayGaugino(2), beta_b.get_MassWB());
   TEST_EQUALITY(beta_a.displayGaugino(3), beta_b.get_MassG());

   TEST_EQUALITY(beta_a.displayMh1Squared(), beta_b.get_mHd2());
   TEST_EQUALITY(beta_a.displayMh2Squared(), beta_b.get_mHu2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mQl), beta_b.get_mq2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mUr), beta_b.get_mu2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mDr), beta_b.get_md2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mLl), beta_b.get_ml2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mEr), beta_b.get_me2());

   TEST_EQUALITY(beta_a.displayTrilinear(UA), beta_b.get_TYu());
   TEST_EQUALITY(beta_a.displayTrilinear(DA), beta_b.get_TYd());
   TEST_EQUALITY(beta_a.displayTrilinear(EA), beta_b.get_TYe());

   const double vu = b.get_vu();
   const double vd = b.get_vd();
   const double tanBeta = vu / vd;
   const double beta_tanBeta = tanBeta * (beta_b.get_vu()/vu - beta_b.get_vd() / vd);
   const double vev = sqrt(sqr(vu) + sqr(vd));
   const double beta_vev = (vu * beta_b.get_vu() + vd * beta_b.get_vd()) / vev;

   TEST_EQUALITY(a.displayTanb(), tanBeta);
   TEST_EQUALITY(beta_a.displayTanb(), beta_tanBeta);
   TEST_EQUALITY(a.displayHvev(), vev);
   TEST_EQUALITY(beta_a.displayHvev(), beta_vev);
   TEST_EQUALITY(beta_a.displaySvev(), beta_b.get_vS());
}

void compare_rges(int loopLevel)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * PI * alpha1);
   const double gY = g1*sqrt(0.6);
   const double g2 = sqrt(4 * PI * alpha2);
   const double g3 = sqrt(4 * PI * ALPHASMZ);
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double vs = 1500.0;
   const double lambda = 0.3;
   const double kappa = 0.4;
   DoubleMatrix Yu(3,3), Yd(3,3), Ye(3,3);
   Yu(3,3) = 165.0   * root2 / (vev * sinBeta);
   Yd(3,3) = 2.9     * root2 / (vev * cosBeta);
   Ye(3,3) = 1.77699 * root2 / (vev * cosBeta);
   DoubleMatrix ID(3, 3), mm0(3, 3);
   for (int i=1; i<=3; i++) ID(i, i) = 1.0;
   mm0 = ID * sqr(m0);

   NMSSMSoftPars sarah;
   sarah.setMu(91);
   sarah.setLoops(loopLevel);
   sarah.set_g1(g1);
   sarah.set_g2(g2);
   sarah.set_g3(g3);
   sarah.set_Yu(Yu);
   sarah.set_Yd(Yd);
   sarah.set_Ye(Ye);
   sarah.set_Lambdax(lambda);
   sarah.set_Kappa(kappa);
   sarah.set_MassB(M12);
   sarah.set_MassG(M12);
   sarah.set_MassWB(M12);
   sarah.set_mq2(mm0);
   sarah.set_ml2(mm0);
   sarah.set_md2(mm0);
   sarah.set_mu2(mm0);
   sarah.set_me2(mm0);
   sarah.set_mHd2(sqr(m0));
   sarah.set_mHu2(sqr(m0));
   sarah.set_TYu(a0 * Yu);
   sarah.set_TYd(a0 * Yd);
   sarah.set_TYe(a0 * Ye);
   sarah.set_vu(vu);
   sarah.set_vd(vd);
   sarah.set_vS(vs);

   SoftParsNmssm softSusy;
   softSusy.setMu(91);
   softSusy.setLoops(loopLevel);
   softSusy.setGaugeCoupling(1, g1);
   softSusy.setGaugeCoupling(2, g2);
   softSusy.setGaugeCoupling(3, g3);
   softSusy.setYukawaMatrix(YU, Yu);
   softSusy.setYukawaMatrix(YD, Yd);
   softSusy.setYukawaMatrix(YE, Ye);
   softSusy.setLambda(lambda);
   softSusy.setKappa(kappa);
   softSusy.setGauginoMass(1, M12);
   softSusy.setGauginoMass(2, M12);
   softSusy.setGauginoMass(3, M12);
   softSusy.setSoftMassMatrix(mQl, mm0);
   softSusy.setSoftMassMatrix(mUr, mm0);
   softSusy.setSoftMassMatrix(mDr, mm0);
   softSusy.setSoftMassMatrix(mLl, mm0);
   softSusy.setSoftMassMatrix(mEr, mm0);
   softSusy.setMh1Squared(sqr(m0));
   softSusy.setMh2Squared(sqr(m0));
   softSusy.setTrilinearMatrix(UA, a0 * Yu);
   softSusy.setTrilinearMatrix(DA, a0 * Yd);
   softSusy.setTrilinearMatrix(EA, a0 * Ye);
   softSusy.setHvev(vev);
   softSusy.setTanb(tanBeta);
   softSusy.setSvev(vs);
   // the following parameters are not peresent in the SARAH classes
   softSusy.setZeta(0.0);
   softSusy.setMu_s(0.0);

   std::cout << "comparing parameters ... ";
   test_parameter_equality(softSusy, sarah);
   std::cout << "done\n";
   std::cout << "comparing beta functions ... ";
   test_beta_function_equality(softSusy, sarah);
   std::cout << "done\n";
   std::cout << "comparing anomalous dimensions ... ";
   compare_anomalous_dimensions(softSusy, sarah);
   std::cout << "done\n";
}

int main()
{

   std::cout << "====================\n";
   std::cout << "compare 1-loop level\n";
   std::cout << "====================\n";
   compare_rges(1);

   std::cout << "====================\n";
   std::cout << "compare 2-loop level\n";
   std::cout << "Note: SoftSusy misses two-loop g^4 terms for the vev\n";
   std::cout << "Note: SARAH misses delta Zhat two-loop terms for the vevs\n";
   std::cout << "====================\n";
   compare_rges(2);

   return 0;
}
