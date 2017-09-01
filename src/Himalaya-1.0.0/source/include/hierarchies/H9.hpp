#ifndef H9_HPP
#define H9_HPP

#include <map>

namespace himalaya{
   
   class H9{
   public:
      /**
       * 	Constructor
       * 	@param flagMap the flagMap for the truncation of expansion variables
       * 	@param Al4p a double alpha_s/4/Pi
       * 	@param beta a double which is the mixing angle beta
       * 	@param Dmst12 a double Mst1^2 - Mst2^2
       * 	@param Dmsqst1  a double Msq1^2 - Mst1^2
       * 	@param lmMt a double log((<renormalization scale> / Mt)^2)
       * 	@param lmMgl a double log((<renormalization scale> / Mgl)^2)
       * 	@param lmMst1 a double log((<renormalization scale> / Mst1)^2)
       * 	@param Mgl a double gluino mass
       * 	@param Mt a double top/bottom quark mass
       * 	@param Mst1 a double stop 1 mass
       * 	@param Mst2 a double stop 2 mass
       * 	@param MuSUSY a double mu parameter
       * 	@param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
       * 	@param mdrFlag an int 0 for DR and 1 for MDR scheme
       * 	@param oneLoopFlag an int flag to consider the one-loop expansion terms
       * 	@param twoLoopFlag an int flag to consider the two-loop expansion terms
       * 	@param threeLoopFlag an int flag to consider the three-loop expansion terms
       */
      H9(std::map<unsigned int, unsigned int> flagMap, double Al4p, double beta, double Dmst12, double Dmsqst1,
		 double lmMt, double lmMgl, double lmMst1,
		 double Mgl, double Mt,  double Mst1, double Mst2, double MuSUSY,
		 double s2t,
		 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag);
      /**
       * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9'
       */
      double getS1();
      /**
       * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9'
       */
      double getS2();
      /**
       * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H9'
       */
      double getS12();
   private:
      double s1 = 0., s2 = 0., s12 = 0.;	/**< The Higgs mass matrix elements s1 = (1, 1), s2 = (2, 2), s12 = (1, 2) = (2, 1) */
   };
}	// himalaya
#endif	// H9_HPP