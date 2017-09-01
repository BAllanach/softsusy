#ifndef H4_HPP
#define H4_HPP

#include <map>

namespace himalaya{
   
   class H4{
   public:
      /**
       * 	Constructor
       * 	@param flagMap the flagMap for the truncation of expansion variables
       * 	@param Al4p a double alpha_s/4/Pi
       * 	@param At a double tri-linear breaking term
       * 	@param beta a double which is the mixing angle beta
       * 	@param lmMt a double log((<renormalization scale> / Mt)^2)
       * 	@param lmMsq a double log((<renormalization scale> / Msq)^2)
       * 	@param lmMst1 a double log((<renormalization scale> / Mst1)^2)
       * 	@param Mt a double top/bottom quark mass
       * 	@param Msusy a double (Mst1 + Mst2 + Mgl) / 3.
       * 	@param Msq a double the average squark mass w/o the top squark
       * 	@param mdrFlag an int 0 for DR and 1 for MDR scheme
       * 	@param oneLoopFlag an int flag to consider the one-loop expansion terms
       * 	@param twoLoopFlag an int flag to consider the two-loop expansion terms
       * 	@param threeLoopFlag an int flag to consider the three-loop expansion terms
       */
      H4(std::map<unsigned int, unsigned int> flagMap, double Al4p, double At, double beta,
		 double lmMt, double lmMsq, double lmMsusy, double Mt, double Msusy, double Msq,
		 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag);
      /**
       * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
       */
      double getS1();
      /**
       * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
       */
      double getS2();
      /**
       * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
       */
      double getS12();
   private:
      double s1 = 0., s2 = 0., s12 = 0.;	/**< The Higgs mass matrix elements s1 = (1, 1), s2 = (2, 2), s12 = (1, 2) = (2, 1) */
   };
}	// himalaya
#endif	// H4_HPP