
/** \file def.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief switches (options) and parameters such as default fermion masses,
   required accuracy etc 
*/

#ifndef DEF_H
#define DEF_H

#include <cmath>
namespace softsusy{
  const char SOFTSUSY_VERSION[] = "3.3.8";

  /// uncomment if you want checking of vector/matrices bounds: slows code
  /// down. It also now checks over/underflows in matrix multiplication etc
  ///  #define ARRAY_BOUNDS_CHECKING 

  /// Make true if you want to include the 2-loop RGE corrections to scalar mass
  /// squared parameters and trilinear terms: they slow it down by a factor of
  /// 3. Note that gaugino and Higgs mass parameters are evolved to 2-loops by
  /// default anyway.
  extern bool INCLUDE_2_LOOP_SCALAR_CORRECTIONS;
  /// Set to number of loops to use for calculation of Higgs mass 
  /// (currently up to 2, the default)
  extern int numHiggsMassLoops;
  /// Set to number of loops to use for REWSB condition up to the default of 2
  extern int numRewsbLoops;
  
  const double EPSTOL = 1.0e-11; ///< underflow accuracy
  const double PI = atan(1.0) * 4.0; ///< or 3.141592653589793 longhand;
  const double root2 = sqrt(2.0);

  extern double GMU;
  extern double MZ; 
  
  /// LEPEWWG central value 14/06/06. Is just used for intialisation etc
  const double MW = 80.404; 
  /// particle data book 2004 central value. Is just used for intialisation etc
  const double MZCENT = 91.1876;
  /// variable for level of output and amount of quark: 0-3, higher numbers
  /// giving more diagnostics. Set by user in file "massIn"
  extern int PRINTOUT;
  /// quark mixing flag: set by user in file "massIn":
  /// 0=no quark mixing, 1=in up sector, 2=in down sector, -1=3rd family
  /// approximation (all at MZ)
  extern int MIXING; 
  /// overall accuracy required
  extern double TOLERANCE;
  /// SUSY breaking scale - if set by user
  extern double QEWSB;
  /// For RPV: do you want to fix the RPV SUSY couplings at MSUSY?
  extern bool susyRpvBCatMSUSY;
  /// Again, for RPV: do you want old school (SLHA-1 like) output?
  extern bool forceSlha1;
  /// For flavour violation: records which SCKM trilinears have been set in
  /// SLHA2 
  extern bool slha2setTrilinear[];
  /// Flag which can be set to switch on producing spectrum output even for
  /// theoretically ruled out regions of parameter space
  extern bool printRuledOutSpectra;
  /// Tries really hard to get convergence when the going gets tough, at the
  /// cost of time for those points (takes up to 43 sec on my computer)
  extern bool tryToConvergeHard;
  /// If =0 (default), sets tachyonic mA=0, otherwise resets mA=sqrt(|mA|^2)
  extern bool mAFlag;
  /// If true, it prints out special developer debugging
  extern bool printDEBUG; 
}
#endif
