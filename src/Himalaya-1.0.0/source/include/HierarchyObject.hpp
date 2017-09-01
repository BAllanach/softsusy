#ifndef HierarchyObject_HPP
#define HierarchyObject_HPP

#include "version.hpp"
#include <Eigen/Eigenvalues>
#include <vector>
#include <map>

namespace himalaya{
   /**
    * 	The HierarchyObject class.
    */
   class HierarchyObject{
   public:
      /**
       * 	A constructor.
       * 	@param isAlphab the boolean which determines wether the members are proportinal to alpha_b or alpha_t.
       */
      HierarchyObject(bool isAlphab);
      /**
       * 	@return the value of isAlphab.
       */
      bool getIsAlphab() const;
      /**
       * 	@return The key to the suitable hierarchy
       */
      int getSuitableHierarchy() const;
      /**
       * 	@return The absolute difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
       */
      double getAbsDiff2L() const;
      /**
       * 	@return The relative difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
       */
      double getRelDiff2L() const;
      /**
       * 	@param loops an integer which can be 1, 2 or 3.
       * 	@return A double which is the expansion uncertainty for the given loop order.
       */
      double getExpUncertainty(int loops) const;
      /**
       * 	@param loops an integer which can be 0, 1, 2, 3. Here, 0 corresponds to the tree-level matrix.
       * 	@return The CP-even Higgs mass matrix at the given loop order.
       */
      Eigen::Matrix2d getDMh(int loops) const;
      /**
       * 	@return The matrix M(MDR) - M(DR) at the order O(alpha_x + alpha_x*alpha_s)
       */
      Eigen::Matrix2d getDRToMDRShift() const;
      /**
       * 	@return A vector of the MDR stop/sbottom masses. The 0th entry corresponds to the lighter particle.
       */
      Eigen::Matrix<double, 2, 1> getMDRMasses() const;
      /**
       * 	@return the MDRFlag integer.
       */
      int getMDRFlag() const;
      /**
       * 	Sets the suitable hierarchy
       * 	@param hierarchy the integer key of the hierarchy.
       */
      void setSuitableHierarchy(int hierarchy);
      /**
       * 	Sets the absolute difference of the Higgs masses at two-loop level
       * 	@param absDiff2L the absolute difference of the Higgs masses as a double.
       */
      void setAbsDiff2L(double absDiff2L);
      /**
       * 	Sets the relative difference ot the Higgs masses at two-loop level
       * 	@param relDiff2L the relative difference of the Higgs masses as a double.
       */
      void setRelDiff2L(double relDiff2L);
      /**
       * 	Sets the uncertainty of the expansion at a given loop level.
       * 	@param loops the integer value of the corresponding loops. Can be 1, 2 or 3.
       * 	@param uncertainty the expansion untertainty at the given loop order as a double.
       */
      void setExpUncertainty(int loops, double uncertainty);
      /**
       * 	Sets the DR -> MDR shift
       * 	@param dMh2L the DR -> MDR shiftet matrix of the form M(MDR) - M(DR).
       */
      void setDRToMDRShift(const Eigen::Matrix2d& dMh2L);
      /**
       * 	Sets the MDR masses
       * 	@param mdrMasses a vector containting the MDR masses with the lightest particle at position 0.
       */
      void setMDRMasses(Eigen::Matrix<double, 2, 1>& mdrMasses);
      /**
       * 	Sets the delta of the CP-even Higgs mass matrix
       * 	@param loops the integer value of the corresponding loops. Can be 0, 1, 2 or 3. 0 corresponds to the tree-level.
       * 	@param dMh the delta of the mass matrix.
       */
      void setDMh(int loops, const Eigen::Matrix2d& dMh);
      /**
       * 	Sets the mdrFlag to calculate the corretions in the DR (0) or MDR (1) scheme
       * 	@param mdrFlag an int. (0) for DR- and (1) for MDR-scheme.
       * 	@throws runtime_exception if the flag is neither 0 or 1 an exception is thrown.
       */
      void setMDRFlag(int mdrFlag);
   private:
      bool isAlphab;									/**< the bool isAlphab */
      int hierarchy;									/**< the suitable hierarchy */
      int mdrFlag = 0;									/**< the MDR-scheme flag */
      double absDiff2L;									/**< the absolute difference of the two loop Higgs masses */
      double relDiff2L;									/**< the relative difference of the two loop Higgs masses */
      std::map<int, double> expUncertainties;						/**< the map which holds the expansion uncertainties, the keys are the loop order: 1, 2, 3 */
      std::map<int, Eigen::Matrix2d> dMhMap;						/**< the map which holds all mass matrices at the given loop order */
      Eigen::Matrix2d mdrShift;								/**< the mass matrix of the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s */
      Eigen::Matrix<double, 2, 1> mdrMasses;						/**< the 'vector' which holds the MDR masses */
      /**
       * 	Sorts a vector.
       * 	@param vector The vector which should be sorted.
       *	@return Returns a vector the lightest entry at position 0.
       */
      Eigen::Matrix<double, 2, 1> sortVector(Eigen::Matrix<double, 2, 1>& vector);
   };
}	// himalaya
#endif	// HierarchyObject_HPP
