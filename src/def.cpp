
/** \file def.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "def.h"

/// global variable declaration
namespace softsusy {
  /// no verbose debug output
  int PRINTOUT = 0;
  /// fractional accuracy required - may need to make smaller for accurate
  double TOLERANCE = 1.0e-4;
  /// decay constant of muon
  double GMU = 1.16637e-5; 
  //If true then the EWSB conditions will output soft Higgs masses
  //Will be inconsistent with constrained models
  //but can be useful for non-universal Higgs cases 
  bool SoftHiggsOut = false;
  /// If true we give the output needed for nmssmTools
  /// otherwise normal nmssm softsusy output
  bool NMSSMTools = false;
  /// number of loops used to calculate Higgs mass and tadpoles. They should be
  /// identical for a consistent calculation
  int numHiggsMassLoops = 2, numRewsbLoops = 2;
  /// global pole mass of MZ in GeV - MZCENT is defined in def.h
  double MZ = MZCENT;
  /// default is to fix RPV parameters at the GUT scale
  bool susyRpvBCatMSUSY = false;
  /// default: use SLHA1 conventions for RPV output
  bool forceSlha1 = false;
  /// default is to have no trilinears set by SLHA2 conventions
  bool slha2setTrilinear[] = {false, false, false, false, false, false, false, 
			      false, false, false, false, false, false, false, 
			      false, false, false, false, false, false, false, 
			      false, false, false, false, false, false};
  /// default is to *not* print out theoretically excluded spectra
  bool printRuledOutSpectra = false;
  /// default is to set tree-level tachyonic A masses to 0 in loops
  bool mAFlag = false;
  /// random number seed
  long idummySave = -22;
  /// Default: use SOFTSUSY conventions for masses of sparticles in loops, ie
  /// tree-level masses computed with the 2-loop Higgs potential
  bool sphenoMassConv = false;
  /// Some parameters used in the computation
  double sw2 = 1.0 - MW * MW / (MZ * MZ), gnuL = 0.5, 
    guL = 0.5 - 2.0 * sw2 / 3.0,
    gdL = -0.5 + sw2 / 3.0, geL = -0.5 + sw2, guR = 2.0 * sw2 / 3.0,
    gdR = -sw2 / 3.0, geR = -sw2, yuL = 1.0 / 3.0, yuR = -4.0 / 3.0,
    ydL = 1.0 / 3.0, ydR = 2.0 / 3.0, yeL = -1.0, yeR = 2.0, ynuL = -1.0;
  /// Default: print out branching ratios bigger than \f$10^{-6}\f$
  double minBR = 1.0e-6;
  /// Default: calculate three-body decays
  bool threeBodyDecays = true;
  /// Default: don't output partial widths in decays
  bool outputPartialWidths = false;
  /// Default: don't calculate decays
  bool calcDecays = false;
  
  /// Includes the evaluation of leading two-loop thresholds corrections
  /// to the strong coupling constant and to the third family of fermion masses 
  bool USE_TWO_LOOP_GAUGE_YUKAWA = false;
  /// just implements decoupling procedure "consistently" for
  /// the case of b-quark mass. It requires the external momentum to be zero. 
  /// However, the difference between the p^2 = 0 and p^2 = mb^2 cases
  /// is expected to be of O((mb/MSUSY)^2), which we can formally neglect.
  bool MB_DECOUPLING = false;

  /// Various two-loop thresholds, eg 2-loop QCD corrections to m_gluino
  bool USE_TWO_LOOP_SPARTICLE_MASS = false;
  /// Default: expand around gluino and squark pole masses
  int expandAroundGluinoPole = 3;
  /// Do you want to calculate svev from EWSB (and input xiS) as opposed to
  /// the converse?
  bool setsVevEWSB = false;
  /// Thresholds at 2-loop SM to electroweak gauge couplings at QED x QCD /
  /// MSSM matching scale. If the QEDxQCD matching scale is not set to MZ,
  /// these two-loop corrections will be ineffective
  //  bool twoLEW = true; 
}
/// end of global variable declaration
