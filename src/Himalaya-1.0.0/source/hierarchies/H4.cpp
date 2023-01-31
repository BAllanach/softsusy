#define Pi M_PI

#include "H4.hpp"
#include "HierarchyCalculator.hpp"
#include "Utils.hpp"
#include <type_traits>
#include <math.h>

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
himalaya::H4::H4(std::map<unsigned int, unsigned int> flagMap, double Al4p, double At, double beta,
		 double lmMt, double lmMsq, double lmMsusy, double Mt, double Msusy, double Msq,
		 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for cos(beta) and sin(beta)
   double Cbeta = cos(beta);
   double Sbeta = sin(beta);
   // zeta functions
   double z2 = pow2(Pi)/6.;
   double z3 = 1.202056903159594;
   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   int shiftst1 = mdrFlag;
   int shiftst2 = mdrFlag;
   int shiftst3 = mdrFlag;
   // expansion flags
   int xAt = flagMap.at(HierarchyCalculator::xxAt);
   int xMsq = flagMap.at(HierarchyCalculator::xxMsq);
   int xlmMsusy = flagMap.at(HierarchyCalculator::xxlmMsusy);
   int xMsusy = flagMap.at(HierarchyCalculator::xxMsusy);
   s1 = 
   #include "../hierarchies/h4/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h4/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h4/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
 */
double himalaya::H4::getS1(){
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
 */
double himalaya::H4::getS2(){
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H4'
 */
double himalaya::H4::getS12(){
   return s12;
}


