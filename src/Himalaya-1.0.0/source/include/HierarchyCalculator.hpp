#ifndef HierarchyCalculator_HPP
#define HierarchyCalculator_HPP

#include "Himalaya_interface.hpp"
#include "HierarchyObject.hpp"
#include "version.hpp"
#include <map>
#include <vector>

namespace himalaya{
   /**
    * The HierarchyCalculatur class 
    */
   class HierarchyCalculator{
   public:
      /**
       * 	Constructor 
       * 	@param p a HimalayaInterface struct
       * 	@param verbose a bool which suppresses the information of the calculation if set to flase
       */
      HierarchyCalculator(const Parameters& p, const bool verbose = true);
      /**
       * 	Calculates the 3-loop mass matrix and other information of the hierarchy selection process.
       * 	@param isAlphab a bool which determines if the returned object is proportinal to alpha_b.
       * 	@return A HierarchyObject which holds all information of the calculation.
       */
      HierarchyObject calculateDMh3L(bool isAlphab, const int mdrFlag = 0);
      /**
       * 	Compares deviation of all hierarchies with the exact two-loop result and returns the hierarchy which minimizes the error.
       * 	@param ho a HierarchyObject with constant isAlphab.
       * 	@param mdrFlag (0) to use the DR-scheme (1) to use the MDR-scheme.
       * 	@return A integer which is identified with the suitable hierarchy.
       */
      int compareHierarchies(HierarchyObject& ho);
      /**
       * 	Calculates the hierarchy contributions for a specific hierarchy at a specific loop order.
       * 	@param ho a HierarchyObject with constant isAlphab.
       * 	@param oneLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded one-loop results to the returned value, respectivley.
       * 	@param twoLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded two-loop results to the returned value, respectivley.
       * 	@param threeLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded three-loop results to the returned value, respectivley.
       * 	@throws runtime_error Throws a runtime_error if the tree-level is requested in terms of hierarchies.
       * 	@return The loop corrected Higgs mass matrix which contains the expanded corrections at the given order.
       */
      Eigen::Matrix2d calculateHierarchy(himalaya::HierarchyObject& ho, const int oneLoopFlagIn, const int twoLoopFlagIn, const int threeLoopFlagIn);
      /**
       * 	Calculates the contribution to the order (alpha_x) and (alpha_s alpha_x) as the difference
       *	of the Higgs mass matrices of the MDR and DR scheme. Here, x can be t or b.
       * 	@param ho a HierarchyObject with constant isAlphab.
       * 	@param shiftOneLoop a bool to shift the terms at one-loop level.
       * 	@param shiftTwoLoop a bool to shift the terms at two-loop level.
       * 	@return The loop corrected Higgs mass matrix difference of the MDR and DR scheme at the given order.
       */
      Eigen::Matrix2d calcDRbarToMDRbarShift(const HierarchyObject& ho, const bool shiftOneLoop, const bool shiftTwoLoop);
      /**
       * 	Calculates the loop corrected Higgs mass matrix at the order O(alpha_x). Here, x can be t or b.
       * 	@param ho a HierarchyObject with constant isAlphab.
       * 	@param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
       * 	@param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
       * 	@return The loop corrected Higgs mass matrix at the order O(alpha_x).
       */
      Eigen::Matrix2d getMt41L(const HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop);
      /**
       * 	Calculates the loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s). Here, x can be t or b.
       * 	@param ho a HierarchyObject with constant isAlphab.
       * 	@param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
       * 	@param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
       * 	@return The loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s).
       */
      Eigen::Matrix2d getMt42L(const HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop);
      /**
       *	Shifts Msx1 according to the hierarchy to the MDR scheme.
       * 	@param ho a HierarchyObject with constant isAlphab.
       * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
       * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
       * 	@return A double which is the MDR sx_1 mass.
       */
      double shiftMst1ToMDR(const HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag);
      /**
       *	Shifts Mst2 according to the hierarchy to the MDR scheme.
       * 	@param ho a HierarchyObject with constant isAlphab.
       * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
       * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
       * 	@return A double which is the MDR stop_2 mass.
       */
      double shiftMst2ToMDR(const HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag);
      /**
       * 	Estimates the uncertainty of the expansion at a given order.
       * 	@param ho a HierarchyObject with constant isAlphab.
       * 	@param massMatrix the CP-even Higgs mass matrix without the corrections whose uncertainty should be estimated.
       * 	@param oneLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the one-loop expansion terms.
       * 	@param twoLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the two-loop expansion terms.
       * 	@param threeLoopFlag an integer flag which is 0 or 1 in order to estimte the uncertainty of the three-loop expansion terms.
       * 	@return A double which is the estimated uncertainty.
       */
      double getExpansionUncertainty(himalaya::HierarchyObject& ho, const Eigen::Matrix2d& massMatrix, const unsigned int oneLoopFlag, 
				     const unsigned int twoLoopFlag, const unsigned int threeLoopFlag);
      /**
       * 	Performs a sanity check of the implemented expansion terms by comparing them to their numerical value at a given parameter point.
       */
      void checkTerms();
      // expansion depth flags
      const static unsigned int xx;		/**< This flag can truncate the two loop expansion at the three loop expansion depth */
      const static unsigned int xxMst;		/**< This flag can truncate the expansion depth of the stop/sbottom masses by one order */
      const static unsigned int xxDmglst1;	/**< This flag can truncate the expansion depth of the difference of stop/sbottom 1 mass and the gluino mass by one order*/
      const static unsigned int xxDmsqst1;	/**< This flag can truncate the expansion depth of the difference of the stop/sbottom 1 mass and the average squark mass by one order*/
      const static unsigned int xxDmst12; 	/**< This flag can truncate the expansion depth of the difference of the stop/sbottom masses by one order*/
      const static unsigned int xxAt;		/**< This flag can truncate the expansion depth of At/Ab by one order*/
      const static unsigned int xxlmMsusy;	/**< This flag can truncate the expansion depth of log(Msusy) by one order*/
      const static unsigned int xxMsq;		/**< This flag can truncate the expansion depth of the average squark mass by one order*/
      const static unsigned int xxMsusy;	/**< This flag can truncate the expansion depth of the average SUSY mass by one order*/
      const static unsigned int xxDmglst2;	/**< This flag can truncate the expansion depth of the difference of the stop/sbottom 2 mass and the gluino mass by one order*/
      const static unsigned int xxDmsqst2;	/**< This flag can truncate the expansion depth of the difference of the average squark mass and the stop/sbottom 2 mass by one order*/
      const static unsigned int xxMgl;		/**< This flag can truncate the expansion depth of the gluino mass by one order*/
   private:
      Parameters p;	/** The HimalayaInterface struct. */
      double Al4p,lmMgl, lmMsq, Mgl, Msq, prefac, z2;/** alpha_s/(4*Pi), log(pow2(p.scale / Mgl)), log(pow2(p.scale / Msq)), Gluino mass, mean Squark mass, prefactor of the Higgs mass matrix, Zeta[2] */
      /**
       * 	Initializes all common variables.
       */
      void init();
      /**
       * 	Checks if a hierarchy is suitable to the given mass spectrum.
       * 	@param ho a HierarchyObject with constant isAlphab and a hierarchy candidate.
       * 	@returns A bool if the hierarchy candidate is suitable to the given mass spectrum.
       */
      bool isHierarchySuitable(const HierarchyObject& ho);
      /**
       * 	Sorts the eigenvalues of a 2x2 matrix.
       * 	@param es the EigenSolver object corresponding to the matrix whose eigenvalues should be sorted.
       * 	@return A sorted vector with the lowest eigenvalue at position 0.
       */
      std::vector<double> sortEigenvalues(const Eigen::EigenSolver<Eigen::Matrix2d> es);
      /**
       * 	@deprecated
       * 	Shifts the 1-loop terms to the MDRbar scheme.
       * 	@param ho a HierarchyObject with constant isAlphab and a hierarchy candidate.
       * 	@return A matrix which corresponds to the difference of the DR and MDR scheme.
       */
      Eigen::Matrix2d getShift(const HierarchyObject& ho);
      std::map<unsigned int, unsigned int> flagMap;	/** A map which holds all hierarchy key value pairs. */
      /**
       * 	A function which maps a boolean to a string.
       * 	@param tf a boolean.
       * 	@return A string which is 'true' if tf is true or 'false' if tf is false.
       */
      std::string tf(const bool tf);
      /**
       * 	Maps a hierarchy to it's mother hierarchy.
       * 	@param hierarchy the key to a hierarchy.
       * 	@throws runtime_error Throws a runtime_error if the given hierarchy is not included.
       * 	@returns The key of the mother hierarchy.
       */
      int getCorrectHierarchy(const int hierarchy);
      /**
       * 	Prints out some information about Himalaya.
       */
      void printInfo();
      //hierarchy keys TODO: use an enum instead?
      static const int h3;	/**< The key to hierarchy h3 */
      static const int h32q2g;	/**< The key to hierarchy h32q2g */
      static const int h3q22g;	/**< The key to hierarchy h3q22g */
      static const int h4;	/**< The key to hierarchy h4 */
      static const int h5;	/**< The key to hierarchy h5 */
      static const int h5g1;	/**< The key to hierarchy h5g1 */
      static const int h6;	/**< The key to hierarchy h6 */
      static const int h6b;	/**< The key to hierarchy h6b */
      static const int h6b2qg2;	/**< The key to hierarchy h6b2qg2 */
      static const int h6bq22g;	/**< The key to hierarchy h6bq22g */
      static const int h6bq2g2;	/**< The key to hierarchy h6bq2g2 */
      static const int h6g2;	/**< The key to hierarchy h6g2 */
      static const int h9;	/**< The key to hierarchy h9 */
      static const int h9q2;	/**< The key to hierarchy h9q2 */
      const std::map<int, int> hierarchyMap = {{h3, h3}, {h32q2g, h3}, {h3q22g, h3}, {h4, h4}, {h5, h5}, {h5g1, h5},
	 {h6, h6}, {h6g2, h6}, {h6b, h6b}, {h6b2qg2, h6b}, {h6bq22g, h6b}, {h6bq2g2, h6b}, {h9, h9}, {h9q2, h9}}; /** The hierarchy map which maps all hierarchies to their mother hierarchies */
  };
}	// himalaya
#endif	// HierarchyCalculator_HPP
