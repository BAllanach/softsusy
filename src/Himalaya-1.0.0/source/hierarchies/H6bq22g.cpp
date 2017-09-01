#define Pi M_PI

#include "H6bq22g.hpp"
#include "HierarchyCalculator.hpp"
#include "Utils.hpp"
#include <type_traits>
#include <math.h>

/**
 * 	Constructor
 * 	@param flagMap the flagMap for the truncation of expansion variables
 * 	@param Al4p a double alpha_s/4/Pi
 * 	@param beta a double which is the mixing angle beta
 * 	@param Dmglst2 a double Mgl - Mst2
 * 	@param Dmsqst2 a double Msq - Mst2
 * 	@param lmMt a double log((<renormalization scale> / Mt)^2)
 * 	@param lmMst1 a double log((<renormalization scale> / Mst1)^2)
 * 	@param lmMst2 a double log((<renormalization scale> / Mst2)^2)
 * 	@param Mt a double top/bottom quark mass
 * 	@param Mst1 a double stop 1 mass
 * 	@param Mst2 a double stop 2 mass
 * 	@param Msq a double average squark mass w/o the stop quark
 * 	@param MuSUSY a double mu parameter
 * 	@param s2t a double 2 times the sine of the stop/sbottom quark mixing angle
 * 	@param mdrFlag an int 0 for DR and 1 for MDR scheme
 * 	@param oneLoopFlag an int flag to consider the one-loop expansion terms
 * 	@param twoLoopFlag an int flag to consider the two-loop expansion terms
 * 	@param threeLoopFlag an int flag to consider the three-loop expansion terms
 */
himalaya::H6bq22g::H6bq22g(std::map<unsigned int, unsigned int> flagMap, double Al4p, double beta, double Dmglst2,
		 double Dmsqst2, double lmMt, double lmMst1, double lmMst2,
		 double Mt, double Mst1, double Mst2, double Msq, double MuSUSY,
		 double s2t,
		 int mdrFlag, int oneLoopFlag, int twoLoopFlag, int threeLoopFlag){
   // abbrev for tan(beta) and sin(beta)
   double Tbeta = tan(beta);
   double Sbeta = sin(beta);
   // zeta functions
   double z2 = pow2(Pi)/6.;
   double z3 = 1.202056903159594;
   double z4 = pow4(Pi)/90.;
   // poly logs
   double pl412 = 0.51747906167389934317668576113647; // PolyLog[4,1/2]
   std::complex<double> pl2expPi3 (0.27415567780803773941206919444, 1.014941606409653625021202554275); // PolyLog[2, Exp[I Pi / 3]]
   std::complex<double> pl3expPi6sqrt3 (0.51928806536375962552715984277228, - 0.33358157526196370641686908633664); // PolyLog[3, Exp[- I Pi / 6] / Sqrt[3]]

   // polylog functions, checked
   double B4 = (-4 * z2 * pow2(log(2)) + 2 / 3.* pow4(log(2)) - 13 / 2. * z4 + 16. * pl412);
   double D3 = 6 * z3 - 15 / 4. * z4 - 6. * pow2(std::imag(pl2expPi3));
   double DN = 6 * z3 - 4 * z2 * pow2(log(2)) + 2 / 3. * pow4(log(2)) - 21 / 2. * z4 + 16. * pl412;
   double OepS2 = - 763 / 32. - (9 * Pi * sqrt(3) * pow2(log(3))) / 16. - (35 * pow3(Pi) * sqrt(3)) / 48.
      + 195 / 16. * z2 - 15 / 4. * z3 + 57 / 16. * z4 + 45 * sqrt(3) / 2. * std::imag(pl2expPi3)
      - 27 * sqrt(3) * std::imag(pl3expPi6sqrt3);
   double S2 = 4 * std::imag(pl2expPi3) / (9. * sqrt(3));
   double T1ep = - 45 / 2. - (Pi * sqrt(3) * pow2(log(3))) / 8. - (35 * pow3(Pi) * sqrt(3)) / 216. - 9 / 2. * z2 + z3 
      + 6. * sqrt(3) * std::imag(pl2expPi3) - 6. * sqrt(3) * std::imag(pl3expPi6sqrt3);
   // mdr flags, indicates if one wants to shift the dr stop mass to the mdr stop mass
   int shiftst1 = mdrFlag;
   int shiftst3 = mdrFlag;
   int xDR2DRMOD = mdrFlag;
   // expansion flags
   int xDmglst2 = flagMap.at(HierarchyCalculator::xxDmglst2);
   int xDmsqst2 = flagMap.at(HierarchyCalculator::xxDmsqst2);
   int xMst = flagMap.at(HierarchyCalculator::xxMst);
   s1 = 
   #include "../hierarchies/h6bq22g/sigS1Full.inc"
   ;
   s2 = 
   #include "../hierarchies/h6bq22g/sigS2Full.inc"
   ;
   s12 = 
   #include "../hierarchies/h6bq22g/sigS12Full.inc"
   ;
}

/**
 * 	@return The diagonal (1, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq22g'
 */
double himalaya::H6bq22g::getS1(){
   return s1;
}

/**
 * 	@return The diagonal (2, 2) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq22g'
 */
double himalaya::H6bq22g::getS2(){
   return s2;
}

/**
 * 	@return The off-diagonal (1, 2) = (2, 1) matrix element of the Higgs mass matrix as a double for the hierarchy 'H6bq22g'
 */
double himalaya::H6bq22g::getS12(){
   return s12;
}


