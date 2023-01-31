#define Pi M_PI

#include "HierarchyCalculator.hpp"
#include "H3.hpp"
#include "H32q2g.hpp"
#include "H3q22g.hpp"
#include "H4.hpp"
#include "H5.hpp"
#include "H5g1.hpp"
#include "H6.hpp"
#include "H6b.hpp"
#include "H6b2qg2.hpp"
#include "H6bq22g.hpp"
#include "H6bq2g2.hpp"
#include "H6g2.hpp"
#include "H9.hpp"
#include "H9q2.hpp"
#include "Utils.hpp"
#include <stdexcept>
#include <iomanip>
#include <type_traits>

extern "C" void dszhiggs_(double *t, double *mg, double *T1, double *T2, double *st, double *ct, double *q, double *mu, double *tanb,
      double *v2, double *gs, int *OS, double *S11, double *S22, double *S12);

static bool isInfoPrinted; /**< If this bool is true, than no info will be printed in further runs */

/**
 * 	Define static variables
 */
namespace himalaya {
   // Hierarchy keys
   const int HierarchyCalculator::h3 		= 0;	/**< The key to hierarchy h3 */
   const int HierarchyCalculator::h32q2g 	= 1;	/**< The key to hierarchy h32q2g */
   const int HierarchyCalculator::h3q22g 	= 2;	/**< The key to hierarchy h3q22g */
   const int HierarchyCalculator::h4 		= 3;	/**< The key to hierarchy h4 */
   const int HierarchyCalculator::h5 		= 4;	/**< The key to hierarchy h5 */
   const int HierarchyCalculator::h5g1 		= 5;	/**< The key to hierarchy h5g1 */
   const int HierarchyCalculator::h6 		= 6;	/**< The key to hierarchy h6 */
   const int HierarchyCalculator::h6b 		= 7;	/**< The key to hierarchy h6b */
   const int HierarchyCalculator::h6b2qg2 	= 8;	/**< The key to hierarchy h6b2qg2 */
   const int HierarchyCalculator::h6bq22g	= 9;	/**< The key to hierarchy h6bq22g */
   const int HierarchyCalculator::h6bq2g2	= 10;	/**< The key to hierarchy h6bq2g2 */
   const int HierarchyCalculator::h6g2		= 11;	/**< The key to hierarchy h6g2 */
   const int HierarchyCalculator::h9 		= 12;	/**< The key to hierarchy h9 */
   const int HierarchyCalculator::h9q2		= 13;	/**< The key to hierarchy h9q2 */
   // expansion depth flags
   const unsigned int HierarchyCalculator::xx				= 14;	/**< This flag can truncate the two loop expansion at the three loop expansion depth */
   const unsigned int HierarchyCalculator::xxMst			= 15;	/**< This flag can truncate the expansion depth of the stop/sbottom masses by one order */
   const unsigned int HierarchyCalculator::xxDmglst1			= 16;	/**< This flag can truncate the expansion depth of the difference of stop/sbottom 1 mass and the gluino mass by one order*/
   const unsigned int HierarchyCalculator::xxDmsqst1			= 17;	/**< This flag can truncate the expansion depth of the difference of the stop/sbottom 1 mass and the average squark mass by one order*/
   const unsigned int HierarchyCalculator::xxDmst12			= 18; 	/**< This flag can truncate the expansion depth of the difference of the stop/sbottom masses by one order*/
   const unsigned int HierarchyCalculator::xxAt				= 19;	/**< This flag can truncate the expansion depth of At/Ab by one order*/
   const unsigned int HierarchyCalculator::xxlmMsusy			= 20;	/**< This flag can truncate the expansion depth of log(Msusy) by one order*/
   const unsigned int HierarchyCalculator::xxMsq			= 21;	/**< This flag can truncate the expansion depth of the average squark mass by one order*/
   const unsigned int HierarchyCalculator::xxMsusy			= 22;	/**< This flag can truncate the expansion depth of the average SUSY mass by one order*/
   const unsigned int HierarchyCalculator::xxDmglst2			= 23;	/**< This flag can truncate the expansion depth of the difference of the stop/sbottom 2 mass and the gluino mass by one order*/
   const unsigned int HierarchyCalculator::xxDmsqst2			= 24;	/**< This flag can truncate the expansion depth of the difference of the average squark mass and the stop/sbottom 2 mass by one order*/
   const unsigned int HierarchyCalculator::xxMgl			= 25;	/**< This flag can truncate the expansion depth of the gluino mass by one order*/
   
}

/**
 * 	Constructor 
 * 	@param p a HimalayaInterface struct
 * 	@param verbose a bool which suppresses the information of the calculation if set to flase
 */
himalaya::HierarchyCalculator::HierarchyCalculator(const Parameters& p, const bool verbose){
   if(!isInfoPrinted && verbose){
      printInfo();
      isInfoPrinted = true;
   }
   this -> p = p;
   this -> p.validate(verbose);

   // imaginary unit
   const std::complex<double> I(0., 1.);

   // Riemann-Zeta
   z2 = pow2(Pi)/6.;

   // init common variables
   init();
}

/**
 * 	Initializes all common variables.
 */
void himalaya::HierarchyCalculator::init(){
   // fill flag list
   flagMap.clear();
   for(unsigned int i = xx; i <= xxMgl; i++){
      flagMap.insert(std::pair<unsigned int, unsigned int> (i, 1));
   }
   // beta
   const double beta = atan(p.vu / p.vd);

   //sw2
   const double sw2 = 1 - pow2(p.MW / p.MZ);

   // Al4p
   Al4p = pow2(p.g3 / (4 * Pi));

   // MGl
   Mgl = p.MG;

   // Msq, checked
   Msq = (2 * sqrt(p.mq2(0, 0)) + sqrt(p.mu2(0, 0)) + sqrt(p.md2(0, 0))	// sup and sdown
      + 2 * sqrt(p.mq2(1, 1)) + sqrt(p.mu2(1, 1)) + sqrt(p.md2(1, 1))	// scharm and sstrange
      // sbottom
      + sqrt(p.mq2(2, 2) + pow2(p.Mb) - (1 / 2. - 1 / 3. * sw2) * pow2(p.MZ) * cos(2 * beta))
      + sqrt(p.md2(2, 2) + pow2(p.Mb) - 1 / 3. * sw2 * pow2(p.MZ) * cos(2 * beta))) / 10.;
   // lmMsq, checked
   lmMsq = log(pow2(p.scale / Msq));

   // lmMgl, checked
   lmMgl = log(pow2(p.scale / Mgl));

   // prefactor, GF = 1/(sqrt(2) * (vu^2 + vd^2)) (here, GF is calculated in the DRbar scheme, checked)
   prefac = (3. / (sqrt(2) * (pow2(p.vu) + pow2(p.vd)) * sqrt(2) * pow2(Pi) * pow2(sin(beta))));
}

/**
 * 	Calculates the 3-loop mass matrix and other information of the hierarchy selection process.
 * 	@param isAlphab a bool which determines if the returned object is proportinal to alpha_b.
 * 	@return A HierarchyObject which holds all information of the calculation.
 */
himalaya::HierarchyObject himalaya::HierarchyCalculator::calculateDMh3L(bool isAlphab, const int mdrFlag){
   HierarchyObject ho (isAlphab);
   
   // set mdrFlag
   ho.setMDRFlag(mdrFlag);
   
   // compare hierarchies and get the best fitting hierarchy
   compareHierarchies(ho);
   
   // calculate the DR to MDR shift with the obtained hierarchy
   if(mdrFlag == 1){
      ho.setDRToMDRShift(calcDRbarToMDRbarShift(ho, true, true));
   }
   else{
      ho.setDRToMDRShift(calcDRbarToMDRbarShift(ho, false, false));
   }
   
   // calculate the 3-loop Higgs mass matrix for the obtained hierarhy
   ho.setDMh(3, calculateHierarchy(ho, 0, 0, 1));
   
   // set the alpha_x contributions
   ho.setDMh(1, getMt41L(ho, mdrFlag, mdrFlag));
   
   // set the alpha_x*alpha_s contributions
   ho.setDMh(2, getMt42L(ho, mdrFlag, mdrFlag));
   
   // estimate the uncertainty of the expansion at 3-loop level
   ho.setExpUncertainty(3, getExpansionUncertainty(ho,
						   ho.getDMh(0) + ho.getDMh(1) + ho.getDMh(2), 0, 0, 1));
   
   // set the uncertainty of the expansion at 1-loop level to 0 by convention, if the user needs this value getExpansionUncertainty should be called
   ho.setExpUncertainty(1, 0.);
   
   return ho;
}

/**
 * 	Compares deviation of all hierarchies with the exact two-loop result and returns the hierarchy which minimizes the error.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@return A integer which is identified with the suitable hierarchy.
 */
int himalaya::HierarchyCalculator::compareHierarchies(himalaya::HierarchyObject& ho){
   // set flags to truncate the expansion
   flagMap.at(xx) = 0;
   flagMap.at(xxMst) = 0;
   double error = -1.;
   int suitableHierarchy = -1;
   // sine of 2 times beta
   const double s2b = sin(2*atan(p.vu/p.vd));
   const double tbeta = p.vu/p.vd;
   // tree level Higgs mass matrix
   Eigen::Matrix2d treelvl;
   treelvl (0,0) = s2b/2.*(pow2(p.MZ) / tbeta + pow2(p.MA) * tbeta);
   treelvl (1,0) = s2b/2.*(-pow2(p.MZ) - pow2(p.MA));
   treelvl (0,1) = treelvl (1,0);
   treelvl (1,1) = s2b/2.*(pow2(p.MZ) * tbeta + pow2(p.MA) / tbeta);
   
   ho.setDMh(0, treelvl);

   // compare the exact higgs mass at 2-loop level with the expanded expressions to find a suitable hierarchy
   for(int hierarchy = h3; hierarchy <= h9q2; hierarchy ++){
      // first, check if the hierarchy is suitable to the mass spectrum
      ho.setSuitableHierarchy(hierarchy);
      if(isHierarchySuitable(ho)){
	 // calculate the exact 1-loop result (only alpha_t/b)
	 Eigen::Matrix2d Mt41L = getMt41L(ho, ho.getMDRFlag(), 0);
	 
	 // call the routine of Pietro Slavich to get the alpha_s alpha_t/b corrections with the MDRbar masses
	 Eigen::Matrix2d Mt42L = getMt42L(ho, ho.getMDRFlag(), 0);
	 
	 // Note: spurious poles are handled by the validate method
	 // of the Himalaya_Interface struct
	 
	 //DEPRECATED calc 1-loop shift for DRbar -> MDRbar
	 //calc difference of Mt41L or Mt41L in the MDRbar scheme directly
	 //it seems that in H3m the sign of the function getShift is wrong (Mt4LDRbar - Mt4LMDRbar)????
	 //to be consistent in the MDRbar-scheme we should subtract Mt4LDRbar, shouldn't we?
	 //Eigen::Matrix2d shift = getShift(hierarchyMap.at(hierarchy), isAlphab);
	 
	 //calculate the exact Higgs mass at 2-loop (only up to alpha_s alpha_t/b)
	 Eigen::EigenSolver<Eigen::Matrix2d> es2L (treelvl + Mt41L + Mt42L);
	 double Mh2l = sortEigenvalues(es2L).at(0);

	 // calculate the expanded 2-loop expression with the specific hierarchy
	 Eigen::EigenSolver<Eigen::Matrix2d> esExpanded (treelvl + Mt41L + calculateHierarchy(ho, 0, 1, 0));
	 
	 // calculate the higgs mass in the given mass hierarchy and compare the result to estimate the error
	 double Mh2LExpanded = sortEigenvalues(esExpanded).at(0);

	 // estimate the error
	 double twoLoopError = fabs((Mh2l - Mh2LExpanded));

	 // estimate the uncertainty of the expansion
	 double expUncertainty = getExpansionUncertainty(ho, treelvl + Mt41L, 0, 1, 0);
	 
	 // add these errors to include the error of the expansion in the comparison
	 double currError = sqrt(pow2(twoLoopError) + pow2(expUncertainty));

	 // if the error is negative, it is the first iteration and there is no hierarchy which fits better
	 if(error < 0){
	    error = currError;
	    suitableHierarchy = hierarchy;
	    ho.setAbsDiff2L(twoLoopError);
	    ho.setRelDiff2L(twoLoopError/Mh2l);
	    ho.setExpUncertainty(2, expUncertainty);
	 }
	 // compare the current error with the last error and choose the hierarchy which fits best (lowest error)
	 else if(currError < error){
	    error = currError;
	    suitableHierarchy = hierarchy;
	    ho.setAbsDiff2L(twoLoopError);
	    ho.setRelDiff2L(twoLoopError/Mh2l);
	    ho.setExpUncertainty(2, expUncertainty);
	 }
      }
   }
   ho.setSuitableHierarchy(suitableHierarchy);
   // reset the flags
   flagMap.at(xx) = 1;
   flagMap.at(xxMst) = 1;
   return suitableHierarchy;
}

//TODO: if one is interested in the expansion at one- and two-loop choose a unified choice for the MDR scheme
/**
 * 	Calculates the hierarchy contributions for a specific hierarchy at a specific loop order.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param oneLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded one-loop results to the returned value, respectivley.
 * 	@param twoLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded two-loop results to the returned value, respectivley.
 * 	@param threeLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded three-loop results to the returned value, respectivley.
 * 	@throws runtime_error Throws a runtime_error if the tree-level is requested in terms of hierarchies.
 * 	@return The loop corrected Higgs mass matrix which contains the expanded corrections at the given order.
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::calculateHierarchy(himalaya::HierarchyObject& ho, const int oneLoopFlagIn,
								  const int twoLoopFlagIn, const int threeLoopFlagIn) {
   // get the hierarchy
   const int hierarchy = ho.getSuitableHierarchy();

   // the hierarchy files containing 1-, 2- and 3-loop terms (alpha_s^0 alpha_t/b, alpha_s alpha_t/b, alpha_s^2 alpha_t/b)
   double sigS1Full = 0., sigS2Full = 0., sigS12Full = 0.;

   // common variables
   double At, Mt, s2t, Mst1 = 0., Mst2 = 0.;
   if (!ho.getIsAlphab()) {
      At = p.At;
      Mt = p.Mt;
      s2t = p.s2t;
   }
   else {
      At = p.Ab;
      Mt = p.Mb;
      s2t = p.s2b;
   }

   const double beta = atan(p.vu / p.vd);
   const double lmMt = log(pow2(p.scale / Mt));

   // this loop is needed to calculate the suitable mass shift order by order
   for(int currentLoopOrder = 1; currentLoopOrder <= 3; currentLoopOrder ++){
      bool runThisOrder;
      double curSig1 = 0., curSig2 = 0., curSig12 = 0.;
      int oneLoopFlag = 0, twoLoopFlag = 0, threeLoopFlag = 0;
      switch (currentLoopOrder){
	 case 1:
	    oneLoopFlag = 1;
	    runThisOrder = oneLoopFlag == oneLoopFlagIn;
	 break;
	 case 2:
	    twoLoopFlag = 1;
	    runThisOrder = twoLoopFlag == twoLoopFlagIn;
	 break;
	 case 3:
	    threeLoopFlag = 1;
	    runThisOrder = threeLoopFlag == threeLoopFlagIn;
	 break;
      }
      if(runThisOrder){
	 // set the Msx masses according to MDRFlag
	 if(oneLoopFlag == 1){
	    Mst1 = shiftMst1ToMDR(ho, 0, 0);
	    Mst2 = shiftMst2ToMDR(ho, 0, 0);
	 }
	 else if(twoLoopFlag == 1){
	    Mst1 = shiftMst1ToMDR(ho, ho.getMDRFlag(), 0);
	    Mst2 = shiftMst2ToMDR(ho, ho.getMDRFlag(), 0);
	 }
	 else if(threeLoopFlag == 1){
	    Mst1 = shiftMst1ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
	    Mst2 = shiftMst2ToMDR(ho, ho.getMDRFlag(), ho.getMDRFlag());
	 }
	 else{
	    throw std::runtime_error("There are no tree-level hierarchies included!");
	 }
	 // select the suitable hierarchy for the specific hierarchy and set variables
	 switch(getCorrectHierarchy(hierarchy)){
	    case h3:{
	       double Dmglst1 = Mgl - Mst1;
	       double Dmsqst1 = pow2(Msq) - pow2(Mst1);
	       double Dmst12 = pow2(Mst1) - pow2(Mst2);
	       double lmMst1 = log(pow2(p.scale / Mst1));
	       switch(hierarchy){
		  case h3:{
		     H3 hierarchy3(flagMap, Al4p, beta,
			Dmglst1, Dmst12, Dmsqst1, lmMt, lmMst1,
			Mgl, Mt, Mst1, Mst2, Msq, p.mu,
			s2t, 
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy3.getS1();
		     curSig2 = hierarchy3.getS2();
		     curSig12 = hierarchy3.getS12();
		  }
		  break;
		  case h32q2g:{
		     H32q2g hierarchy32q2g(flagMap, Al4p, beta,
			Dmglst1, Dmst12, Dmsqst1, lmMt, lmMst1,
			Mt, Mst1, Mst2, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy32q2g.getS1();
		     curSig2 = hierarchy32q2g.getS2();
		     curSig12 = hierarchy32q2g.getS12();
		  }
		  break;
		  case h3q22g:{
		     H3q22g hierarchy3q22g(flagMap, Al4p, beta,
			Dmglst1, Dmst12, Dmsqst1, lmMt, lmMst1,
			Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy3q22g.getS1();
		     curSig2 = hierarchy3q22g.getS2();
		     curSig12 = hierarchy3q22g.getS12();
		  }
		  break;
	       }
	    }
	    break;
	    case h4:{
	       double Msusy = (Mst1 + Mst2 + Mgl) / 3.;
	       double lmMsusy = log(pow2(p.scale / Msusy));
	       H4 hierarchy4(flagMap, Al4p, At, beta,
		  lmMt, lmMsq, lmMsusy, Mt, Msusy, Msq,
		  ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
	       curSig1 = hierarchy4.getS1();
	       curSig2 = hierarchy4.getS2();
	       curSig12 = hierarchy4.getS12();
	    }
	    break;
	    case h5:{
	       double Dmglst1 = Mgl - Mst1;
	       double lmMst1 = log(pow2(p.scale / Mst1));
	       double lmMst2 = log(pow2(p.scale / Mst2));
	       switch(hierarchy){
		  case h5:{
		     H5 hierarchy5(flagMap, Al4p, beta, Dmglst1,
			lmMt, lmMst1, lmMst2, lmMsq, Mt, Mst1,
			Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy5.getS1();
		     curSig2 = hierarchy5.getS2();
		     curSig12 = hierarchy5.getS12();
		  }
		  break;
		  case h5g1:{
		     H5g1 hierarchy5g1(flagMap, Al4p, beta, Dmglst1,
			lmMt, lmMst1, lmMst2, lmMsq,
			Mgl, Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy5g1.getS1();
		     curSig2 = hierarchy5g1.getS2();
		     curSig12 = hierarchy5g1.getS12();
		  }
		  break;
	       }
	    }
	    break;
	    case h6:{
	       double Dmglst2 = Mgl - Mst2;
	       double lmMst1 = log(pow2(p.scale / Mst1));
	       double lmMst2 = log(pow2(p.scale / Mst2));
	       switch(hierarchy){
		  case h6:{
		     H6 hierarchy6(flagMap, Al4p, beta, Dmglst2,
			lmMt, lmMst1, lmMst2, lmMsq,
			Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6.getS1();
		     curSig2 = hierarchy6.getS2();
		     curSig12 = hierarchy6.getS12();
		  }
		  break;
		  case h6g2:{
		     H6g2 hierarchy6g2(flagMap, Al4p, beta, Dmglst2,
			lmMt, lmMst1, lmMst2, lmMsq,
			Mgl, Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6g2.getS1();
		     curSig2 = hierarchy6g2.getS2();
		     curSig12 = hierarchy6g2.getS12();
		  }
		  break;
	       }
	    }
	    break;
	    case h6b:{
	       double Dmglst2 = Mgl - Mst2;
	       double Dmsqst2 = Msq - Mst2;
	       double lmMst1 = log(pow2(p.scale / Mst1));
	       double lmMst2 = log(pow2(p.scale / Mst2));
	       switch(hierarchy){
		  case h6b:{
		     H6b hierarchy6b(flagMap, Al4p, beta, Dmglst2,
			Dmsqst2, lmMt, lmMst1, lmMst2,
		       Mt, Mst1, Mst2, p.mu,
		       s2t,
		       ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6b.getS1();
		     curSig2 = hierarchy6b.getS2();
		     curSig12 = hierarchy6b.getS12();
		  }
		  break;
		  case h6b2qg2:{
		     H6b2qg2 hierarchy6b2qg2(flagMap, Al4p, beta, Dmglst2,
			Dmsqst2, lmMt, lmMst1, lmMst2,
			Mgl, Mt, Mst1, Mst2, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6b2qg2.getS1();
		     curSig2 = hierarchy6b2qg2.getS2();
		     curSig12 = hierarchy6b2qg2.getS12();
		  }
		  break;
		  case h6bq22g:{
		     H6bq22g hierarchy6bq22g(flagMap, Al4p, beta, Dmglst2,
			Dmsqst2, lmMt, lmMst1, lmMst2,
			Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6bq22g.getS1();
		     curSig2 = hierarchy6bq22g.getS2();
		     curSig12 = hierarchy6bq22g.getS12();
		  }
		  break;
		  case h6bq2g2:{
		     H6bq2g2 hierarchy6bq2g2(flagMap, Al4p, beta, Dmglst2,
			Dmsqst2, lmMt, lmMst1, lmMst2,
			Mgl, Mt, Mst1,Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy6bq2g2.getS1();
		     curSig2 = hierarchy6bq2g2.getS2();
		     curSig12 = hierarchy6bq2g2.getS12();
		  }
		  break;
	       }
	    }
	    break;
	    case h9:{
	       double lmMst1 = log(pow2(p.scale / Mst1));
	       double Dmst12 = pow2(Mst1) - pow2(Mst2);
	       double Dmsqst1 = pow2(Msq) - pow2(Mst1);
	       switch(hierarchy){
		  case h9:{
		     H9 hierarchy9(flagMap, Al4p, beta, Dmst12, Dmsqst1,
			lmMt, lmMgl, lmMst1,
			Mgl, Mt, Mst1, Mst2, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy9.getS1();
		     curSig2 = hierarchy9.getS2();
		     curSig12 = hierarchy9.getS12();
		  }
		  break;
		  case h9q2:{
		     H9q2 hierarchy9q2(flagMap, Al4p, beta, Dmst12, Dmsqst1,
			lmMt, lmMgl, lmMst1,
			Mgl, Mt, Mst1, Mst2, Msq, p.mu,
			s2t,
			ho.getMDRFlag(), oneLoopFlag, twoLoopFlag, threeLoopFlag);
		     curSig1 = hierarchy9q2.getS1();
		     curSig2 = hierarchy9q2.getS2();
		     curSig12 = hierarchy9q2.getS12();
		  }
		  break;
	       }
	    }
	    break;
	 }
      }
      sigS1Full += curSig1;
      sigS2Full += curSig2;
      sigS12Full += curSig12;
   }
   // add the MDR masses to the hierarchy object only if a 3-loop calculation has to be done, otherwise let the user decide
   if(oneLoopFlagIn == 0 && twoLoopFlagIn == 0 && threeLoopFlagIn == 1){
      Eigen::Matrix<double,2,1> mdrMasses;
      mdrMasses(0) = Mst1;
      mdrMasses(1) = Mst2;
      ho.setMDRMasses(mdrMasses);
   }
   Eigen::Matrix2d higgsMassMatrix;
   higgsMassMatrix(0, 0) = prefac * sigS1Full;
   higgsMassMatrix(0, 1) = prefac * sigS12Full;
   higgsMassMatrix(1, 0) = higgsMassMatrix(0, 1);
   higgsMassMatrix(1, 1) = prefac * sigS2Full;
   return higgsMassMatrix;
}

/**
 * 	Checks if a hierarchy is suitable to the given mass spectrum.
 * 	@param ho a HierarchyObject with constant isAlphab and a hierarchy candidate.
 * 	@returns A bool if the hierarchy candidate is suitable to the given mass spectrum.
 */
bool himalaya::HierarchyCalculator::isHierarchySuitable(const himalaya::HierarchyObject& ho){
   double Mst1, Mst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
   }
   switch (ho.getSuitableHierarchy()){
      case h3:
	 return Mgl > Mst2;
      case h32q2g:
	 return (Mst2 >= Msq) && (Mst2 > Mgl);
      case h3q22g:
	 return (Msq > Mst2) && (Mst2 > Mgl);
      case h4:
	 return (Mst1 < Msq) && (Mst1 >= Mgl);
      case h5:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mgl - Mst1) < std::abs(Mgl - Mst2)) && (Mst2 < Msq) && (Mst1 >= Mgl);
      case h5g1:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mgl - Mst1) < std::abs(Mgl - Mst2)) && (Mst2 < Msq) && (Mgl > Mst1);
      case h6:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 < Msq) && (Mst2 >= Mgl);
      case h6g2:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 < Msq) && (Mgl > Mst2);
      case h6b:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 >= Msq) && (Mst2 >= Mgl);
      case h6b2qg2:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Mst2 >= Msq) && (Mgl > Mst2);
      case h6bq22g:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Msq > Mst2) && (Mst2 >= Mgl);
      case h6bq2g2:
	 return (Mst2 - Mst1 > 0.1*Mst1) && ((Mst2 - Mgl) < std::abs(Mgl - Mst1)) && (Msq > Mst2) && (Mgl > Mst2);
      case h9:
	 return (Mst2 >= Msq) && ((Mst2 - Mst1) < (Mst1 - Mgl));
      case h9q2:
	 return (Msq > Mst2) && ((Mst1 - Mst1) < (Mst1 - Mgl));
   }
   return false;
}

/**
 * 	Shifts Msx1 according to the hierarchy to the MDR scheme.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
 * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
 * 	@return A double which is the MDR sx_1 mass.
 */
double himalaya::HierarchyCalculator::shiftMst1ToMDR(const himalaya::HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag) {
   double Mst1mod = 0., Mst1, Mst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
   }
   double lmMst2 = log(pow2(p.scale) / pow2(Mst2));
   double Dmglst2 = Mgl - Mst2;
   double mdr2mst1ka = (-8. * twoLoopFlag * pow2(Al4p) * (10 * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2) + pow2(Mst2) * (-1 + 2 * lmMst2 + 2 * z2))) / (3. * pow2(Mst1));
   switch (getCorrectHierarchy(ho.getSuitableHierarchy())) {
   case h3:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case h4:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case h5:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   case h6:
      Mst1mod = (144 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) * pow4(Msq) + 27 * (1 + mdr2mst1ka) * pow4(Msq) * pow2(Mst1) +
         twoLoopFlag * pow2(Al4p) * Mgl * (-5 * (67 + 84 * lmMgl - 84 * lmMsq) * pow5(Mgl) - 40 * (43 + 30 * lmMgl - 30 * lmMsq) * pow3(Mgl) * pow2(Msq) +
            288 * Dmglst2 * pow4(Msq) * (1 - 2 * z2) + 12 * Mgl * pow4(Msq) * (79 + 144 * pow2(lmMgl) - 150 * lmMsq +
               90 * pow2(lmMsq) - 90 * lmMgl * (-3 + 2 * lmMsq) + 208 * z2))) / (27. * pow4(Msq) * pow2(Mst1));
      break;
   case h6b:
      Mst1mod = (48 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) + 9 * (1 + mdr2mst1ka) * pow2(Mst1) +
         8 * twoLoopFlag * pow2(Al4p) * (-135 * pow2(Msq) + 12 * Dmglst2 * Mgl * (1 - 22 * z2) +
            pow2(Mgl) * (77 + 135 * lmMgl + 72 * pow2(lmMgl) - 75 * lmMsq -
               90 * lmMgl * lmMsq + 45 * pow2(lmMsq) + 104 * z2))) / (9. * pow2(Mst1));
      break;
   case h9:
      Mst1mod = (1 + mdr2mst1ka);
      break;
   }
   return Mst1 * sqrt(Mst1mod);
}

/**
 * 	Shifts Msx2 according to the hierarchy to the MDR scheme.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
 * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
 * 	@return A double which is the MDR sx_2 mass.
 */
double himalaya::HierarchyCalculator::shiftMst2ToMDR(const himalaya::HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag) {
   double Mst2mod = 0., Mst2;
   if(!ho.getIsAlphab()){
      Mst2 = p.MSt(1);
   }
   else{
      Mst2 = p.MSb(1);
   }
   double Dmglst2 = Mgl - Mst2;
   double mdr2mst2ka = (-80. * twoLoopFlag * pow2(Al4p) * pow2(Msq) * (-1 + 2 * lmMsq + 2 * z2)) / (3. * pow2(Mst2));
   switch (getCorrectHierarchy(ho.getSuitableHierarchy())) {
   case h3:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case h4:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case h5:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   case h6:
      Mst2mod = (144 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) * pow4(Msq) + 27 * (1 + mdr2mst2ka) * pow4(Msq) * pow2(Mst2) +
         twoLoopFlag * pow2(Al4p) * Mgl * (-5 * (67 + 84 * lmMgl - 84 * lmMsq) * pow5(Mgl) - 40 * (43 + 30 * lmMgl - 30 * lmMsq) * pow3(Mgl) * pow2(Msq) +
            288 * Dmglst2 * pow4(Msq) * (1 - 2 * z2) + 12 * Mgl * pow4(Msq) * (79 + 144 * pow2(lmMgl) - 150 * lmMsq +
               90 * pow2(lmMsq) - 90 * lmMgl * (-3 + 2 * lmMsq) + 208 * z2))) / (27. * pow4(Msq) * pow2(Mst2));
      break;
   case h6b:
      Mst2mod = (48 * oneLoopFlag * Al4p * (1 + lmMgl) * pow2(Mgl) + 9 * (1 + mdr2mst2ka) * pow2(Mst2) +
         8 * twoLoopFlag * pow2(Al4p) * (-135 * pow2(Msq) + 12 * Dmglst2 * Mgl * (1 - 22 * z2) +
            pow2(Mgl) * (77 + 135 * lmMgl + 72 * pow2(lmMgl) - 75 * lmMsq -
               90 * lmMgl * lmMsq + 45 * pow2(lmMsq) + 104 * z2))) / (9. * pow2(Mst2));
      break;
   case h9:
      Mst2mod = (1 + mdr2mst2ka);
      break;
   }
   return Mst2 * sqrt(Mst2mod);
}


/**
 * 	Sorts the eigenvalues of a 2x2 matrix.
 * 	@param es the EigenSolver object corresponding to the matrix whose eigenvalues should be sorted.
 * 	@return A sorted vector with the lowest eigenvalue at position 0.
 */
std::vector<double> himalaya::HierarchyCalculator::sortEigenvalues(const Eigen::EigenSolver<Eigen::Matrix2d> es){
  std::vector<double> sortedEigenvalues = {sqrt(std::real(es.eigenvalues()(0))), sqrt(std::real(es.eigenvalues()(1)))};
  std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
  return sortedEigenvalues;
}

/**
 * 	Calculates the loop corrected Higgs mass matrix at the order O(alpha_x). Here, x can be t or b.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
 * 	@param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
 * 	@return The loop corrected Higgs mass matrix at the order O(alpha_x).
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getMt41L(const himalaya::HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop){
   Eigen::Matrix2d Mt41L;
   double GF = 1/(sqrt(2) * (pow2(p.vu) + pow2(p.vd)));
   double Mst1;
   double Mst2;
   double Mt;
   double s2t;
   const double beta = atan(p.vu/p.vd);
   Mst1 = shiftMst1ToMDR(ho, shiftOneLoop, shiftTwoLoop);
   Mst2 = shiftMst2ToMDR(ho, shiftOneLoop, shiftTwoLoop);
   if(!ho.getIsAlphab()){
      s2t = p.s2t;
      Mt = p.Mt;
   }
   else{
      s2t = p.s2b;
      Mt = p.Mb;
   }
   Mt41L (0,0) = (-3*GF*pow2(Mt)*pow2(p.mu)*pow2(1/sin(beta))*
      (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
        pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
        pow2(Mst2)*log(Mst2))*pow2(s2t))/
    (4.*sqrt(2)*(pow2(Mst1) - pow2(Mst2))*pow2(Pi));
   Mt41L (0,1) = (3*GF*pow2(1/sin(beta))*
      (-(pow3(Mt)*p.mu*(log(Mst1) - log(Mst2))*s2t)/2. + 
        (pow2(Mt)*pow2(p.mu)*1/tan(beta)*
           (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
             pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
             pow2(Mst2)*log(Mst2))*pow2(s2t))/
         (4.*(pow2(Mst1) - pow2(Mst2))) + 
        (Mt*p.mu*(-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
             pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
             pow2(Mst2)*log(Mst2))*pow3(s2t))/8.))/
    (sqrt(2)*pow2(Pi));
   Mt41L (1,0) = Mt41L(0,1);
   Mt41L (1,1) =  (3*GF*pow2(1/sin(beta))*
      (pow4(Mt)*(log(Mst1) + log(Mst2) - 2*log(Mt)) + 
        pow3(Mt)*p.mu*1/tan(beta)*(log(Mst1) - log(Mst2))*s2t + 
        (pow2(Mt)*pow2(1/sin(beta))*
           (pow2(Mst1)*pow2(p.mu)*pow2(cos(beta)) - 
             pow2(Mst2)*pow2(p.mu)*pow2(cos(beta)) - 
             pow2(Mst1)*pow2(p.mu)*pow2(cos(beta))*log(Mst1) - 
             pow2(Mst2)*pow2(p.mu)*pow2(cos(beta))*log(Mst1) + 
             pow2(Mst1)*pow2(p.mu)*pow2(cos(beta))*log(Mst2) + 
             pow2(Mst2)*pow2(p.mu)*pow2(cos(beta))*log(Mst2) + 
             2*pow4(Mst1)*log(Mst1)*pow2(sin(beta)) - 
             4*pow2(Mst1)*pow2(Mst2)*log(Mst1)*pow2(sin(beta)) + 
             2*pow4(Mst2)*log(Mst1)*pow2(sin(beta)) - 
             2*pow4(Mst1)*log(Mst2)*pow2(sin(beta)) + 
             4*pow2(Mst1)*pow2(Mst2)*log(Mst2)*pow2(sin(beta)) - 
             2*pow4(Mst2)*log(Mst2)*pow2(sin(beta)))*pow2(s2t))/
         (4.*(pow2(Mst1) - pow2(Mst2))) - 
        (Mt*p.mu*1/tan(beta)*(-pow2(Mst1) + pow2(Mst2) + 
             pow2(Mst1)*log(Mst1) + pow2(Mst2)*log(Mst1) - 
             pow2(Mst1)*log(Mst2) - pow2(Mst2)*log(Mst2))*
           pow3(s2t))/4. - 
        ((pow2(Mst1) - pow2(Mst2))*
           (-pow2(Mst1) + pow2(Mst2) + pow2(Mst1)*log(Mst1) + 
             pow2(Mst2)*log(Mst1) - pow2(Mst1)*log(Mst2) - 
             pow2(Mst2)*log(Mst2))*pow4(s2t))/16.))/
    (sqrt(2)*pow2(Pi));
    return Mt41L;
}

/**
 * 	@deprecated
 * 	Shifts the 1-loop terms to the MDRbar scheme.
 * 	@param ho a HierarchyObject with constant isAlphab and a hierarchy candidate.
 * 	@return A matrix which corresponds to the difference of the DR and MDR scheme.
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getShift(const himalaya::HierarchyObject& ho){
   Eigen::Matrix2d shift;
   double GF = 1/(sqrt(2) * (pow2(p.vu) + pow2(p.vd)));
   double Mst1;
   double Mst2;
   double Mt;
   double s2t;
   const double beta = atan(p.vu/p.vd);
   double deltamst1, deltamst2;
   if(!ho.getIsAlphab()){
      Mst1 = p.MSt(0);
      Mst2 = p.MSt(1);
      s2t = p.s2t;
      Mt = p.Mt;
      deltamst1 = -(Mst1 - shiftMst1ToMDR(ho, 1, 0));
      deltamst2 = -(Mst2 - shiftMst2ToMDR(ho, 1, 0));
   }
   else{
      Mst1 = p.MSb(0);
      Mst2 = p.MSb(1);
      s2t = p.s2b;
      Mt = p.Mb;
      deltamst1 = -(Mst1 - shiftMst1ToMDR(ho, 1, 0));
      deltamst2 = -(Mst2 - shiftMst2ToMDR(ho, 1, 0));
   }
   shift(0, 0) = (3 * GF * (deltamst2 * Mst1 - deltamst1 * Mst2) * pow2(Mt) * pow2(p.mu) * pow2(1 / Pi) *
      pow2(1/sin(beta)) * pow2(1 / (pow2(Mst1) - pow2(Mst2))) * pow2(s2t) *
      (4 * (log(Mst1) - log(Mst2)) * pow2(Mst1) * pow2(Mst2) - pow4(Mst1) +
	 pow4(Mst2))) / (4. * sqrt(2) * Mst1 * Mst2);
   shift(0, 1) = (-3 * GF * Mt * p.mu * pow2(1 / Pi) * pow2(1/sin(beta)) *
      pow2(1 / (pow2(Mst1) - pow2(Mst2))) * s2t *
      (-(pow2(pow2(Mst1) - pow2(Mst2)) *
	    (4 * (-(deltamst2 * Mst1) + deltamst1 * Mst2) * pow2(Mt) +
	       (-2 * Mst1 * Mst2 * (deltamst1 * Mst1 + deltamst2 * Mst2) *
		  (log(Mst1) - log(Mst2)) +
		  (deltamst2 * Mst1 + deltamst1 * Mst2) * pow2(Mst1) -
		  (deltamst2 * Mst1 + deltamst1 * Mst2) * pow2(Mst2)) * pow2(s2t))) + 2 * (deltamst2 * Mst1 - deltamst1 * Mst2) * Mt * p.mu * 1/tan(beta) *
	 (4 * (log(Mst1) - log(Mst2)) * pow2(Mst1) * pow2(Mst2) - pow4(Mst1) +
	    pow4(Mst2)) * s2t)) / (8. * sqrt(2) * Mst1 * Mst2);
   shift(1, 0) = shift(0, 1);
   shift(1, 1) = (3 * GF * pow2(1 / Pi) * pow2(1/sin(beta)) *
      ((Mt * p.mu * 1/tan(beta) * (-(deltamst1 * Mst1) + deltamst2 * Mst2 +
	       2 * deltamst1 * Mst1 * log(Mst1) + 2 * deltamst2 * Mst2 * log(Mst1) -
	       2 * deltamst1 * Mst1 * log(Mst2) - 2 * deltamst2 * Mst2 * log(Mst2) -
	       (deltamst2 * pow2(Mst1)) / Mst2 + (deltamst1 * pow2(Mst2)) / Mst1) *
	    pow3(s2t)) / 4. + (pow2(Mt) * pow2(1 / (pow2(Mst1) - pow2(Mst2))) *
	    pow2(s2t) * (2 * deltamst2 * pow7(Mst1) -
	       2 * deltamst1 * pow6(Mst1) * Mst2 -
	       2 * deltamst2 * Mst1 * pow6(Mst2) + 2 * deltamst1 * pow7(Mst2) -
	       4 * deltamst1 * pow6(Mst1) * Mst2 * log(Mst1) +
	       4 * deltamst2 * Mst1 * pow6(Mst2) * log(Mst1) +
	       4 * deltamst1 * pow6(Mst1) * Mst2 * log(Mst2) -
	       4 * deltamst2 * Mst1 * pow6(Mst2) * log(Mst2) -
	       deltamst2 * pow5(Mst1) * pow2(p.mu) * pow2(1/tan(beta)) -
	       deltamst1 * pow5(Mst2) * pow2(p.mu) * pow2(1/tan(beta)) +
	       2 * deltamst2 * pow2(Mst2) *
	       (pow5(Mst1) * (-3 + 2 * log(Mst1) - 2 * log(Mst2)) +
		  2 * (log(Mst1) - log(Mst2)) * pow2(p.mu) * pow2(1/tan(beta)) *
		  pow3(Mst1)) - 2 * deltamst1 * pow2(Mst1) *
	       (pow5(Mst2) * (3 + 2 * log(Mst1) - 2 * log(Mst2)) +
		  2 * (log(Mst1) - log(Mst2)) * pow2(p.mu) * pow2(1/tan(beta)) *
		  pow3(Mst2)) + deltamst1 * Mst2 * pow2(p.mu) * pow2(1/tan(beta)) *
	       pow4(Mst1) + 6 * deltamst1 * pow3(Mst2) * pow4(Mst1) +
	       8 * deltamst1 * log(Mst1) * pow3(Mst2) * pow4(Mst1) -
	       8 * deltamst1 * log(Mst2) * pow3(Mst2) * pow4(Mst1) +
	       deltamst2 * Mst1 * pow2(p.mu) * pow2(1/tan(beta)) * pow4(Mst2) +
	       6 * deltamst2 * pow3(Mst1) * pow4(Mst2) -
	       8 * deltamst2 * log(Mst1) * pow3(Mst1) * pow4(Mst2) +
	       8 * deltamst2 * log(Mst2) * pow3(Mst1) * pow4(Mst2))) / (4. * Mst1 * Mst2) -
	 ((deltamst2 * Mst1 + deltamst1 * Mst2) * pow4(Mt)) / (Mst1 * Mst2) -
	 ((deltamst2 * pow5(Mst1) + deltamst1 * pow5(Mst2) -
	    4 * deltamst2 * pow2(Mst2) * pow3(Mst1) -
	    4 * deltamst1 * pow2(Mst1) * pow3(Mst2) +
	    3 * deltamst1 * Mst2 * pow4(Mst1) -
	    4 * deltamst1 * Mst2 * log(Mst1) * pow4(Mst1) +
	    4 * deltamst1 * Mst2 * log(Mst2) * pow4(Mst1) +
	    3 * deltamst2 * Mst1 * pow4(Mst2) +
	    4 * deltamst2 * Mst1 * log(Mst1) * pow4(Mst2) -
	    4 * deltamst2 * Mst1 * log(Mst2) * pow4(Mst2)) * pow4(s2t)) /
	 (16. * Mst1 * Mst2) + (-(deltamst1 / Mst1) + deltamst2 / Mst2) * p.mu *
	 1/tan(beta) * pow3(Mt) * s2t)) / sqrt(2);
   return shift;
}


/**
 * 	Calculates the loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s). Here, x can be t or b.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
 * 	@param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
 * 	@return The loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s).
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::getMt42L(const himalaya::HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop){
   Eigen::Matrix2d Mt42L;
   double S11, S12, S22;
   double Mt2;
   double MG = p.MG;
   double Mst12;
   double Mst22;
   double st;
   double ct;
   Mst12 = pow2(shiftMst1ToMDR(ho, shiftOneLoop, shiftTwoLoop));
   Mst22 = pow2(shiftMst2ToMDR(ho, shiftOneLoop, shiftTwoLoop));
   if(!ho.getIsAlphab()){
      const double theta = asin(p.s2t)/2.;
      Mt2 = pow2(p.Mt);
      st = sin(theta);
      ct = cos(theta);
   }
   else{
      const double theta = asin(p.s2b)/2.;
      Mt2 = pow2(p.Mb);
      st = sin(theta);
      ct = cos(theta);
   }
   double scale2 = pow2(p.scale);
   // note the sign difference in mu
   double mu = - p.mu;
   double tanb = p.vu/p.vd;
   double v2 = pow2(p.vu) + pow2(p.vd);
   double gs = p.g3;
   int os = 0;
   dszhiggs_(&Mt2, &MG, &Mst12, &Mst22, &st, &ct, &scale2, &mu, &tanb, &v2, &gs, &os, &S11, &S22, &S12);
   Mt42L(0, 0) = S11;
   Mt42L(1, 0) = S12;
   Mt42L(0, 1) = S12;
   Mt42L(1, 1) = S22;
   return Mt42L;
}

/**
 * 	Calculates the contribution to the order (alpha_x) and (alpha_s alpha_x) as the difference
 * 	of the Higgs mass matrices of the MDR and DR scheme. Here, x can be t or b.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param shiftOneLoop a bool to shift the terms at one-loop level.
 * 	@param shiftTwoLoop a bool to shift the terms at two-loop level.
 * 	@return The loop corrected Higgs mass matrix difference of the MDR and DR scheme at the given order.
 */
Eigen::Matrix2d himalaya::HierarchyCalculator::calcDRbarToMDRbarShift(const himalaya::HierarchyObject& ho, const bool shiftOneLoop, const bool shiftTwoLoop){
   if(shiftOneLoop && shiftTwoLoop){
      return getMt41L(ho, 1, 1) + getMt42L(ho, 1, 1) - getMt41L(ho, 0, 0) - getMt42L(ho, 0, 0);
   }
   else if(shiftOneLoop){
      return getMt41L(ho, 1, 1) - getMt41L(ho, 0, 0);
   }
   else if(shiftTwoLoop){
      return getMt42L(ho, 1, 1) - getMt42L(ho, 0, 0);
   }
   else{
      return Eigen::Matrix2d::Zero();
   }
}


/**
 * 	Estimates the uncertainty of the expansion at a given order.
 * 	@param ho a HierarchyObject with constant isAlphab.
 * 	@param massMatrix the CP-even Higgs mass matrix without the corrections whose uncertainty should be estimated.
 * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the one-loop expansion terms.
 * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the two-loop expansion terms.
 * 	@param threeLoopFlag an integer flag which is 0 or 1 in order to estimte the uncertainty of the three-loop expansion terms.
 * 	@return A double which is the estimated uncertainty.
 */
double himalaya::HierarchyCalculator::getExpansionUncertainty(himalaya::HierarchyObject& ho, const Eigen::Matrix2d& massMatrix, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag, const unsigned int threeLoopFlag){
   double Mh;
   double Mhcut;
   std::vector<double> errors;
   // reset flags
   flagMap.at(xxMst) = 1;
   Eigen::EigenSolver<Eigen::Matrix2d> es;
   switch (getCorrectHierarchy(ho.getSuitableHierarchy())) {
   case h3:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      // truncate the expansion at all variables with one order lower than the expansion depth and evaluate the expansion uncertainty 
      flagMap.at(xxDmglst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmglst1) = 1;
      flagMap.at(xxDmsqst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmsqst1) = 1;
      flagMap.at(xxDmst12) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmst12) = 1;
      break;
   case h4:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxAt) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxAt) = 1;
      flagMap.at(xxlmMsusy) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxlmMsusy) = 1;
      flagMap.at(xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMsq) = 1;
      flagMap.at(xxMsusy) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMsusy) = 1;
      break;
   case h5:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxDmglst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmglst1) = 1;
      flagMap.at(xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMsq) = 1;
      break;
   case h6:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxDmglst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmglst2) = 1;
      flagMap.at(xxMsq) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMsq) = 1;
      break;
   case h6b:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxDmglst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmglst2) = 1;
      flagMap.at(xxDmsqst2) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmsqst2) = 1;
      break;
   case h9:
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mh = sortEigenvalues(es).at(0);
      flagMap.at(xxDmsqst1) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmsqst1) = 1;
      flagMap.at(xxDmst12) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxDmst12) = 1;
      flagMap.at(xxMgl) = 0;
      es.compute(massMatrix + calculateHierarchy(ho, oneLoopFlag, twoLoopFlag, threeLoopFlag), false);
      Mhcut = sortEigenvalues(es).at(0);
      errors.push_back(fabs(Mh - Mhcut));
      flagMap.at(xxMgl) = 1;
      break;
   }
   // evalue the sqrt of the squared errors
   double squaredErrorSum = 0.;
   for(auto const& error: errors){
      squaredErrorSum = squaredErrorSum + pow2(error);
   }
   // set the expansion depth for the next comparison
   flagMap.at(xxMst) = 0;
   flagMap.at(xx) = 0;
   return sqrt(squaredErrorSum);
}

/*
 * 	input parameters for the check
 */
himalaya::Parameters checkTermsXt33(){

   himalaya::Parameters pars;

   pars.scale = 1973.75;
   pars.mu = 1999.82;
   pars.g3 =  1.02907;
   pars.vd = 49.5751;
   pars.vu = 236.115;
   pars.mq2 <<  4.00428e+06 , 0, 0,
               0, 4.00428e+06, 0,
               0, 0, 3.99786e+06;
   pars.md2 << 4.00361e+06, 0, 0,
               0, 4.00361e+06, 0,
               0, 0, 4.00346e+06;
   pars.mu2 << 4.00363e+06 , 0, 0,
               0, 4.00363e+06, 0,
               0, 0, 3.99067e+06;
   pars.Ab = 9996.81;
   pars.At = 6992.34;

   pars.MA = 1992.14;
   pars.MG = 2000.96;
   pars.MW = 76.7777;
   pars.MZ = 88.4219;
   pars.Mt = 147.295;
   pars.Mb = 2.23149;
   pars.MSt << 1745.3 , 2232.1;
   pars.MSb <<  2000.14, 2001.09;
   pars.s2t = -0.999995;
   pars.s2b =-0.550527;

   return pars;
}

/**
 * 	Performs a sanity check of the implemented expansion terms by comparing them to their numerical value at a given parameter point.
 */
void himalaya::HierarchyCalculator::checkTerms(){
   p = checkTermsXt33();
   init();
   himalaya::HierarchyObject ho (false);
   ho.setMDRFlag(1);
   for(int i = h3; i <= h9q2; i++){
      ho.setSuitableHierarchy(i);
      Eigen::Matrix2d oloMat = calculateHierarchy(ho, 1,0,0);
      Eigen::Matrix2d twloMat = calculateHierarchy(ho, 0,1,0);
      Eigen::Matrix2d thloMat = calculateHierarchy(ho, 0,0,1);
      bool ck1LPassed, ck2LPassed, ck3LPassed;
      switch(i){
	 case h3:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.3521101999062) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-13.4248798129061) < 1e-06 && twloMat(1,0) - 10.91323060388626 < 1e-06 && twloMat(1,1) - 1477.154147584478 < 1e-06;
	    ck3LPassed = thloMat(0,0) - 1.163582002655875 < 1e-06 && thloMat(1,0) - 9.897241079348351 < 1e-06 && thloMat(1,1) - 369.9741236956309 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h32q2g:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.3521101999062) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-13.66052379180129) < 1e-06 && twloMat(1,0) - 11.26755617866339 < 1e-06 && twloMat(1,1) - 1477.465656153518 < 1e-06;
	    ck3LPassed = thloMat(0,0) - 1.113051431370291 < 1e-06 && thloMat(1,0) - 9.903809573970422 < 1e-06 && thloMat(1,1) - 369.7408109643386 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h3q22g:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.3521101999062) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-13.66052379180129) < 1e-06 && twloMat(1,0) - 11.26755617866339 < 1e-06 && twloMat(1,1) - 1477.465656153518 < 1e-06;
	    ck3LPassed = thloMat(0,0) - 1.058450932536496 < 1e-06 && thloMat(1,0) - 10.0141272838662 < 1e-06 && thloMat(1,1) - 370.3301180635573 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h4:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 0 < 1e-06 && oloMat(1,0) - 0 < 1e-06 && oloMat(1,1) - 6685.123085628641 < 1e-06;
	    ck2LPassed = twloMat(0,0) - 0 < 1e-06 && twloMat(1,0) - 1183.325484493686 < 1e-06 && twloMat(1,1) - 1458.970501474495 < 1e-06;
	    ck3LPassed = thloMat(0,0) - 162.1379208650191 < 1e-06 && thloMat(1,0) - 326.0219627343553 < 1e-06 && thloMat(1,1) - 431.6926278454841 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h5:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 15921.69462848581 < 1e-06 && oloMat(1,0) - (-388569.2043081555) < 1e-06 && oloMat(1,1) - 7874.401574063407 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-86.77887344841422) < 1e-06 && twloMat(1,0) - (-20625.63783863484) < 1e-06 && twloMat(1,1) - (-42446.62009872038) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 2442.115080578889 < 1e-06 && thloMat(1,0) - (-3859.942907446577) < 1e-06 && thloMat(1,1) - 60593.055768119 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h5g1:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 15921.69462848581 < 1e-06 && oloMat(1,0) - (-388569.2043081556) < 1e-06 && oloMat(1,1) - 7874.401574063407 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-114.6037388932203) < 1e-06 && twloMat(1,0) - (-20341.84471909946) < 1e-06 && twloMat(1,1) - (-42843.48046642416) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 2415.507513838155 < 1e-06 && thloMat(1,0) - (-3766.750163753644) < 1e-06 && thloMat(1,1) - 59380.34497121828 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h6:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702315 < 1e-06 && oloMat(1,0) - (-184.7601614832763) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1078.578574572312) < 1e-06 && twloMat(1,0) - 7096.529601647042 < 1e-06 && twloMat(1,1) - (-1927.791631086123) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 245.4412216221288 < 1e-06 && thloMat(1,0) - 573.1296253278389 < 1e-06 && thloMat(1,1) - 8448.4582538127 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h6b:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702311 < 1e-06 && oloMat(1,0) - (-184.7601614832763) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1078.578574572312) < 1e-06 && twloMat(1,0) - 7096.52960164704 < 1e-06 && twloMat(1,1) - (-1900.197036824461) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 283.0253770519464 < 1e-06 && thloMat(1,0) - 566.2182257407396 < 1e-06 && thloMat(1,1) - 10093.33785879814 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h6b2qg2:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702311 < 1e-06 && oloMat(1,0) - (-184.7601614832759) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1089.201418061661) < 1e-06 && twloMat(1,0) - 7145.267026465748 < 1e-06 && twloMat(1,1) - (-2077.345120153528) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 285.3154791763894 < 1e-06 && thloMat(1,0) - 544.3654284413091 < 1e-06 && thloMat(1,1) - 10336.22756889787 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h6bq22g:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702315 < 1e-06 && oloMat(1,0) - (-184.7601614832763) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1078.578574572311) < 1e-06 && twloMat(1,0) - 7096.529601647042 < 1e-06 && twloMat(1,1) - (-1900.197036824461) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 283.0220052455883 < 1e-06 && thloMat(1,0) - 566.2190953470737 < 1e-06 && thloMat(1,1) - 10093.33986048966 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h6bq2g2:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702315 < 1e-06 && oloMat(1,0) - (-184.7601614832759) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1089.201418061661) < 1e-06 && twloMat(1,0) - 7145.267026465748 < 1e-06 && twloMat(1,1) - (-2077.345120153528) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 285.3120881213721 < 1e-06 && thloMat(1,0) - 544.3662758149513 < 1e-06 && thloMat(1,1) - 10336.23012077387 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h6g2:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - 9272.477351702315 < 1e-06 && oloMat(1,0) - (-184.7601614832761) < 1e-06 && oloMat(1,1) - 7581.278122072418 < 1e-06;
	    ck2LPassed = twloMat(0,0) - (-1089.201418061661) < 1e-06 && twloMat(1,0) - 7145.267026465748 < 1e-06 && twloMat(1,1) - (-2112.642999123034) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 246.0217489966267 < 1e-06 && thloMat(1,0) - 557.451210096066 < 1e-06 && thloMat(1,1) - 8628.076480526881 < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h9:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.352110199906) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - 420.2050380976995 < 1e-06 && twloMat(1,0) - (-554.6021924866435) < 1e-06 && twloMat(1,1) - (-797.8089039452509) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 132.8584579769461 < 1e-06 && thloMat(1,0) - (-171.9326869339159) < 1e-06 && thloMat(1,1) - (-800.8408283898472) < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
	 case h9q2:
	    ck1LPassed = false;
	    ck2LPassed = false;
	    ck3LPassed = false;
	    ck1LPassed = oloMat(0,0) - (-1033.437882123761) < 1e-06 && oloMat(1,0) - (-394.3521101999065) < 1e-06 && oloMat(1,1) - 17633.47392819223 < 1e-06;
	    ck2LPassed = twloMat(0,0) - 420.2050380976993 < 1e-06 && twloMat(1,0) - (-554.6021924866436) < 1e-06 && twloMat(1,1) - (-797.8089039452487) < 1e-06;
	    ck3LPassed = thloMat(0,0) - 132.6358855624267 < 1e-06 && thloMat(1,0) - (-171.4711818838455) < 1e-06 && thloMat(1,1) - (-800.9569014303727) < 1e-06;
	    std::cout << "Hierarchy " << i << " passed checks 1L: " << tf(ck1LPassed) << " 2L: " << tf(ck2LPassed) << " 3L: "<< tf(ck3LPassed) << "." << "\n";
	 break;
      }
   }
}

/**
 * 	A function which maps a boolean to a string.
 * 	@param tf a boolean.
 * 	@return A string which is 'true' if tf is true or 'false' if tf is false.
 */
std::string himalaya::HierarchyCalculator::tf(const bool tf){
   return tf ? "true" : "false";
}

/**
 * 	Maps a hierarchy to it's mother hierarchy.
 * 	@param hierarchy the key to a hierarchy.
 * 	@throws runtime_error Throws a runtime_error if the given hierarchy is not included.
 * 	@returns The key of the mother hierarchy.
 */
int himalaya::HierarchyCalculator::getCorrectHierarchy(const int hierarchy){
   if(hierarchy < 0 || hierarchy > 13){
      throw std::runtime_error("\033[1;31m Error: Hierarchy " + std::to_string(hierarchy) + " not included!\033[0m");
   }
   return hierarchyMap.at(hierarchy);
}

/**
 * 	Prints out some information about Himalaya.
 */
void himalaya::HierarchyCalculator::printInfo(){
   std::cout << "....................................................." << "\n";
   std::cout << "Himalaya " << Himalaya_VERSION_MAJOR << "." << Himalaya_VERSION_MINOR << "." << Himalaya_VERSION_RELEASE << "\tѧѦ ѧ \n";
   std::cout << "Uses code by: P. Slavich et al. (2-loop at*as)." << "\n";
   std::cout << "Uses the 3-loop at*as^2 contributions of Kant et al." << "\n";
   std::cout << "....................................................." << "\n";
}
