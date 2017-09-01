#ifndef Himalaya_interface_HPP
#define Himalaya_interface_HPP

#include "version.hpp"
#include <complex>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
#include <algorithm>

namespace himalaya {

typedef Eigen::Matrix<double,2,1> V2;
typedef Eigen::Matrix<double,2,2> RM22;
typedef Eigen::Matrix<double,3,3> RM33;

/**
 * 	The Himalaya interface struct
 */
struct Parameters {
   // DR-bar parameters
   double scale{};				/**< renormalization scale */
   double mu{};					/**< mu parameter */
   double g3{};					/**< gauge coupling g3 SU(3) */
   double vd{};					/**< VEV of down Higgs */
   double vu{};					/**< VEV of up Higgs */
   RM33 mq2{RM33::Zero()};			/**< soft-breaking squared left-handed squark mass parameters */
   RM33 md2{RM33::Zero()};			/**< soft-breaking squared right-handed down-squark mass parameters */
   RM33 mu2{RM33::Zero()};			/**< soft-breaking squared right-handed up-squark mass parameters */
   double At{};					/**< trilinear stop-Higgs coupling */
   double Ab{};					/**< trilinear sbottom-Higgs coupling */

   // DR-bar masses
   double MG{};					/**< gluino */
   double MW{};					/**< W */
   double MZ{};					/**< Z */
   double Mt{};					/**< top-quark */
   double Mb{};					/**< down-quark */
   double MA{};					/**< CP-odd Higgs */
   V2 MSt{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};	/**< stops */
   V2 MSb{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};	/**< sbottoms */

   // DR-bar mixing angles
   double s2t = std::numeric_limits<double>::quiet_NaN();		/**< sine of 2 times the stop mixing angle */
   double s2b = std::numeric_limits<double>::quiet_NaN();		/**< sine of 2 times the sbottom mixing angle */
   
   /**
    * 	Checks if the stop/sbottom masses and mixing angles are provided. Otherwise calculate them.
    * 	Checks if the stop/sbottom masses are ordered in the right way. If these masses are wrongly ordered
    * 	the right ordering will be introduced.
    * 	Checks if the stops/sbottom masses are degenerated and introduce a small shift to the 1st stop/sbottom mass in this case.
    * 	@param verbose a bool which suppresses the information of the calculation if set to flase
    */
   void validate(const bool verbose){
      // check if stop/sbottom masses and/or mixing angles are nan. If so, calculate these quantities.
      if(std::isnan(MSt(0)) || std::isnan(MSt(1)) || std::isnan(s2t)){
	 double beta = atan(vu / vd);
	 double Xt = Mt * (At - mu * 1 / tan(beta));
	 double sw2 = 1 - MW * MW / MZ / MZ;
	 Eigen::Matrix2d stopMatrix;
	 stopMatrix << mq2(2, 2) + Mt * Mt + (1/2. - 2/3. * sw2) * MZ * MZ * cos(2 * beta), Xt,
	    Xt, mu2(2, 2) + Mt * Mt + 2 / 3. * sw2 * MZ * MZ * cos(2 * beta);
	 // solve eigenvalues and sort them
	 Eigen::EigenSolver<Eigen::Matrix2d> es(stopMatrix);
	 std::vector<double> sortedEigenvalues = {sqrt(std::real(es.eigenvalues()(0))), sqrt(std::real(es.eigenvalues()(1)))};
         std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
	 // set stop masses
	 MSt << sortedEigenvalues.at(0), sortedEigenvalues.at(1);
	 // extraxt mixing angle
	 double delta1 = std::abs(acos(std::real(es.eigenvectors().col(0)(0))) + asin(std::real(es.eigenvectors().col(0)(1))));
	 double delta2 = std::abs(acos(std::real(es.eigenvectors().col(1)(0))) + asin(std::real(es.eigenvectors().col(1)(1))));
	 
	 double theta;
	 if(delta1 < delta2){
	    theta = acos(std::real(es.eigenvectors().col(0)(0)));
	 }
	 else{
	    theta = -acos(std::real(es.eigenvectors().col(1)(0)));
	 }
	 s2t = sin(2 * theta);
	 if(verbose){
	    std::cout << "\033[1;34m Info:\033[0m Stop masses or mixing angle not provided. Calculated values:\n" <<
	       "\tstop masses: " << MSt(0) << " GeV, " << MSt(1) << " GeV,\n" << 
	       "\tmixing angle: " << theta << ".\n";
	 }
      }
      if(std::isnan(MSb(0)) || std::isnan(MSb(1)) || std::isnan(s2b)){
	 double beta = atan(vu / vd);
	 double Xb = Mb * (Ab - mu * tan(beta));
	 double sw2 = 1 - MW * MW / MZ / MZ;
	 Eigen::Matrix2d sbottomMatrix;
	 sbottomMatrix << mq2(2, 2) + Mb * Mb - (1/2. - 1/3. * sw2) * MZ * MZ * cos(2 * beta), Xb,
	    Xb, md2(2, 2) + Mb * Mb - 1/3. * sw2 * MZ * MZ * cos(2 * beta);
	 // solve eigenvalues and sort them
	 Eigen::EigenSolver<Eigen::Matrix2d> es(sbottomMatrix);
	 std::vector<double> sortedEigenvalues = {sqrt(std::real(es.eigenvalues()(0))), sqrt(std::real(es.eigenvalues()(1)))};
         std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
	 // set sbottom masses
	 MSb << sortedEigenvalues.at(0), sortedEigenvalues.at(1);
	 // extract mixing angle
	 double delta1 = std::abs(acos(std::real(es.eigenvectors().col(0)(0))) + asin(std::real(es.eigenvectors().col(0)(1))));
	 double delta2 = std::abs(acos(std::real(es.eigenvectors().col(1)(0))) + asin(std::real(es.eigenvectors().col(1)(1))));
	 
	 double theta;
	 if(delta1 < delta2){
	    theta = acos(std::real(es.eigenvectors().col(0)(0)));
	 }
	 else{
	    theta = -acos(std::real(es.eigenvectors().col(1)(0)));
	 }
	 s2b = sin(2 * theta);
	 if(verbose){
	    std::cout << "\033[1;34m Info:\033[0m Sbottom masses or mixing angle not provided. Calculated values:\n" <<
	       "\tsbottom masses: " << MSb(0) << " GeV, " << MSb(1) << " GeV,\n" << 
	       "\tmixing angle: " << theta << ".\n";
	 }
      }
      // check the ordering of the stop/sbottom quarks
      if (MSt(0) > MSt(1)) {
	 std::swap(MSt(0), MSt(1));
	 s2t *= -1;
      }

      if (MSb(0) > MSb(1)) {
	 std::swap(MSb(0), MSb(1));
	 s2b *= -1;
      }

      // check if the stop/sbottom masses are degenerated. If this is the case one could get spurious poles
      // in Pietro's code. To avoid this numerical issue we shift the stop/bottom 1 mass by a relative (but small)
      // value.
      if(std::abs(MSt(0) - MSt(1)) < 1.0E-5){
	MSt(0) = MSt(1) / (1. + 1.0E-5);
      }

      if(std::abs(MSb(0) - MSb(1)) < 1.0E-5){
	MSb(0) = MSb(0) / (1. + 1.0E-5);
      }
   };
};

}	//	himalaya

#endif	//	Himalaya_interface_HPP
