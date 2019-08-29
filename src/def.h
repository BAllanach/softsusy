
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
  /// uncomment if you want to use looptools (slower) rather than SOFTSUSY's
  /// defined loop functions: good for testing
  /// #define USE_LOOPTOOLS

  /// uncomment if you want checking of vector/matrices bounds: slows code
  /// down. It also now checks over/underflows in matrix multiplication etc
  ///  #define ARRAY_BOUNDS_CHECKING 

  /// Set to true (default) if you want to include 2-loop SM contributions to
  /// electroweak gauge coupling threshold effects
  //  extern bool twoLEW;
  
  extern bool NMSSMTools;
  extern bool SoftHiggsOut;
  /// Set to number of loops to use for calculation of Higgs mass 
  /// (currently up to 3 with Himalaya; the default is 2)
  extern int numHiggsMassLoops;
  /// Set to number of loops to use for REWSB condition up to the default of 2
  extern int numRewsbLoops;
  
  const double EPSTOL = 1.0e-11; ///< underflow accuracy
  const double PI = atan(1.0) * 4.0; ///< or 3.141592653589793 longhand;
  const double root2 = sqrt(2.0);
  ///< used to flag diabolical problems
  const double numberOfTheBeast = 6.66e66; 
  
  extern double GMU; ///< Fermi constant
  extern double MZ;  ///< Z boson mass
  /// LEPEWWG central value 14/06/06. Is just used for intialisation etc
  const double MW = 80.404; 
  /// particle data book 2004 central value. Is just used for intialisation etc
  const double MZCENT = 91.1876;
  /// variable for level of output and amount of quark: 0-3, higher numbers
  /// giving more diagnostics. Set by user in file "massIn"
  extern int PRINTOUT;
  /// overall accuracy required
  extern double TOLERANCE;
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
  /// If =0 (default), sets tachyonic mA=0, otherwise resets mA=sqrt(|mA|^2)
  extern bool mAFlag;
  /// A handy global variable for random number generator
  extern long idummySave;
  /// minimum branching ratio that will get printed out
  extern double minBR;
  /// If false, don't calculate the three-body decays
  extern bool threeBodyDecays;
  /// If true, output partial widths of decays in SLHA comments
  extern bool outputPartialWidths;
  /// If true, calculate decays
  extern bool calcDecays;
  
  /// Includes the evaluation of leading two-loop thresholds corrections
  /// to the strong coupling constant and to the third family of fermion masses 
  extern bool USE_TWO_LOOP_GAUGE_YUKAWA; 
  /// just implements decoupling procedure "consistently" for
  /// the case of b-quark mass. It requires the external momentum to be zero. 
  /// However, the difference between the p^2 = 0 and p^2 = mb^2 cases
  /// is expected to be of O((mb/MSUSY)^2), which we can formally neglect.
  extern bool MB_DECOUPLING;

  enum { ENABLE_TWO_LOOP_MT_AS  = 0x1,    
	 ENABLE_TWO_LOOP_AS_AS_YUK = 0x2,
	 ENABLE_TWO_LOOP_MB_AS  = 0x4, 
	 ENABLE_TWO_LOOP_MB_YUK = 0x8,    
	 ENABLE_TWO_LOOP_MTAU_YUK = 0x10    
  };    

  /// Various two-loop thresholds, eg 2-loop QCD corrections to m_gluino
  extern bool USE_TWO_LOOP_SPARTICLE_MASS;
  /// If = 0, no expansion for gluino. If = 1, expand around gluino and squark
  /// pole masses. If = 2, expand only around gluino pole mass rather than the
  /// tree-level mass and iterate. 1 and 2 are a little slower. 
  extern int expandAroundGluinoPole;
  /// Quick start bool to trap svev setting rather than xiS setting
  extern bool setsVevEWSB;
}

#endif
