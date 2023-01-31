#include "HierarchyObject.hpp"

/**
 * 	A constructor.
 * 	@param isAlphab the boolean which determines wether the members are proportinal to alpha_b or alpha_t.
 */
himalaya::HierarchyObject::HierarchyObject(bool isAlphab){
   this -> isAlphab = isAlphab;
}

/**
 * 	@return The value of isAlphab.
 */
bool himalaya::HierarchyObject::getIsAlphab() const{
   return isAlphab;
}

/**
 * 	Sets the suitable hierarchy
 * 	@param hierarchy the integer key of the hierarchy.
 */
void himalaya::HierarchyObject::setSuitableHierarchy(int hierarchy){
   this -> hierarchy = hierarchy;
}

/**
 * 	@return The key to the suitable hierarchy
 */
int himalaya::HierarchyObject::getSuitableHierarchy() const{
   return hierarchy;
}

/**
 * 	Sets the absolute difference of the Higgs masses at two-loop level
 * 	@param absDiff2L the absolute difference of the Higgs masses as a double.
 */
void himalaya::HierarchyObject::setAbsDiff2L(double absDiff2L){
   this -> absDiff2L = absDiff2L;
}

/**
 * 	@return The absolute difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
 */
double himalaya::HierarchyObject::getAbsDiff2L() const{
   return absDiff2L;
}

/**
 * 	Sets the relative difference ot the Higgs masses at two-loop level
 * 	@param relDiff2L the relative difference of the Higgs masses as a double.
 */
void himalaya::HierarchyObject::setRelDiff2L(double relDiff2L){
   this -> relDiff2L = relDiff2L;
}

/**
 * 	@return The relative difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
 */
double himalaya::HierarchyObject::getRelDiff2L() const{
   return relDiff2L;
}


/**
 * 	Sets the uncertainty of the expansion at a given loop level.
 * 	@param loops the integer value of the corresponding loops. Can be 1, 2 or 3.
 * 	@param uncertainty the expansion untertainty at the given loop order as a double.
 */
void himalaya::HierarchyObject::setExpUncertainty(int loops, double uncertainty){
   if(loops > 0 && loops <=3){
      expUncertainties.insert(std::pair<int, double> (loops, uncertainty));
   }
   else {
      throw std::runtime_error("Expansion uncertainty for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

/**
 * 	@param loops an integer which can be 1, 2 or 3.
 * 	@return A double which is the expansion uncertainty for the given loop order.
 */
double himalaya::HierarchyObject::getExpUncertainty(int loops) const{
   if(loops > 0 && loops <=3){
      return expUncertainties.at(loops);
   }
   else {
      throw std::runtime_error("Expansion uncertainty for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

/**
 * 	Sets the delta of the CP-even Higgs mass matrix
 * 	@param loops the integer value of the corresponding loops. Can be 0, 1, 2 or 3. 0 corresponds to the tree-level.
 * 	@param dMh the delta of the mass matrix.
 */
void himalaya::HierarchyObject::setDMh(int loops, const Eigen::Matrix2d& dMh){
   if(loops >= 0 && loops <=3){
      dMhMap.insert(std::pair<int, Eigen::Matrix2d> (loops, dMh));
   }
   else {
      throw std::runtime_error("Higgs mass matrix for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

/**
 * 	@param loops an integer which can be 0, 1, 2, 3. Here, 0 corresponds to the tree-level matrix.
 * 	@return The CP-even Higgs mass matrix at the given loop order.
 */
Eigen::Matrix2d himalaya::HierarchyObject::getDMh(int loops) const{
   if(loops >= 0 && loops <=3){
      return dMhMap.at(loops);
   }
   else {
      throw std::runtime_error("Higgs mass matrix for " + std::to_string(loops) + " loop(s) is not available.");
   }
}

/**
 * 	Sets the DR -> MDR shift
 * 	@param dMh2L the DR -> MDR shiftet matrix of the form M(MDR) - M(DR).
 */
void himalaya::HierarchyObject::setDRToMDRShift(const Eigen::Matrix2d& mdrShift){
   this -> mdrShift = mdrShift;
}

/**
 * 	@return The matrix M(MDR) - M(DR) at the order O(alpha_x + alpha_x*alpha_s)
 */
Eigen::Matrix2d himalaya::HierarchyObject::getDRToMDRShift() const{
   return mdrShift;
}

/**
 * 	Sets the MDR masses
 * 	@param mdrMasses a vector containting the MDR masses with the lightest particle at position 0.
 */
void himalaya::HierarchyObject::setMDRMasses(Eigen::Matrix<double,2,1>& mdrMasses){
   this -> mdrMasses = sortVector(mdrMasses);
}

/**
 * 	@return A vector of the MDR stop/sbottom masses. The 0th entry corresponds to the lighter particle.
 */
Eigen::Matrix<double,2,1> himalaya::HierarchyObject::getMDRMasses() const{
   return mdrMasses;
}

/**
 * 	Sets the mdrFlag to calculate the corretions in the DR (0) or MDR (1) scheme
 * 	@param mdrFlag an int. (0) for DR- and (1) for MDR-scheme
 * 	@throws runtime_exception if the flag is neither 0 or 1 an exception is thrown.
 */
void himalaya::HierarchyObject::setMDRFlag(int mdrFlag){
   if(mdrFlag != 0 && mdrFlag != 1)
      throw std::runtime_error("The MDR-flag has to be 0 (DR-scheme) or 1 (MDR-scheme). Input: " + std::to_string(mdrFlag) + ".");
   this -> mdrFlag = mdrFlag;
}

/**
 * 	@return the MDRFlag integer.
 */
int himalaya::HierarchyObject::getMDRFlag() const{
   return mdrFlag;
}



/**
 * 	Sorts a vector.
 * 	@param vector The vector which should be sorted.
 *	@return Returns a vector the lightest entry at position 0.
 */
Eigen::Matrix<double,2,1> himalaya::HierarchyObject::sortVector(Eigen::Matrix<double,2,1>& vector){
   // checks if all variables are ordered in the right way
   if (vector(0) > vector(1)) {
      std::swap(vector(0), vector(1));
   }
   return vector;
}
