
/** \file def.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "def.h"

/// global variable declaration
/// no quark mixing (dominant third family approx), and no verbose output
namespace softsusy {
  int MIXING = 1, PRINTOUT = 0;
  /// fractional accuracy required - may need to make smaller for accurate
  /// neutrino masses
  double TOLERANCE = 1.0e-4;
  /// decay constant of muon
  double GMU = 1.16637e-5; 
  /// there are two possible conventions: if QEWSB > MZ, its value is assumed
  /// in GeV and used as a constant MSUSY. Otherwise, it MULTIPLIES the usual 
  /// MSUSY value, of root(mstop1 mstop2)
  double QEWSB = 1.0; 
  /// Do we include 2-loop RGEs of *all* scalar masses and A-terms, or only the
  /// scalar mass Higgs parameters? (Other quantities all 2-loop anyway): the
  /// default in SOFTSUSY 3 is to include all 2-loop terms
  bool INCLUDE_2_LOOP_SCALAR_CORRECTIONS = true;
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
  /// default is to do a quick calculation instead of trying hard to get
  /// convergence 
  bool tryToConvergeHard = false;
  /// default is to set tree-level tachyonic A masses to 0 in loops
  bool mAFlag = false;
  /// default is that developer-level debugging printout flag is OFF
  bool printDEBUG = false;
  /// random number seed
  long idummySave = -22;
  /// Default is to not include 2-loop MSSM corrections in degenerate
  /// gluino/squark mass limit 
  bool includeTwoLoopMssmCorrectionsToMt = false;
  /// Default: use SOFTSUSY conventions for masses of sparticles in loops, ie
  /// tree-level masses computed with the 2-loop Higgs potential
  bool sphenoMassConv = false;
}
/// end of global variable declaration
