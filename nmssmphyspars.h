
/** \file physpars.h
   - Project:     SOFTSUSY 
   - File:        physpars.h
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Flags for potential problems in sProblem structure, and structure
   for containing physical MSSM parameters such as pole masses in sPhysical. 
*/

#ifndef NMSSMPHYSPARSH
#define NMSSMPHYSPARSH

#include <iostream>
#include "linalg.h"
using namespace softsusy;
using namespace std;

typedef enum {none=0, selectron, smuon, stau, sup, scharm, stop, sdown, sstrange, sbottom, h0, A0, hpm, snue, snumu, snutau} tachyonType;

const static char* tachyonNames[16] = {
  "none", "selectron", "smuon", "stau", "sup", "scharm", "stop", "sdown", 
  "sstrange", "sbottom", "h0", "A0", "hpm", "snue", "snumu", "snutau"};

/// Masses of the physical particles. 
struct nmsPhysical {
  DoubleVector mh0, mA0;
  double mHpm;
  //PA: new nmssm masses (mH0, mA0 no longer correspond to phyical states but leaving them in at least for just for sanity checks)
   double mH1, mH2, mA1, mA2;
  /// sneutrino masses
  DoubleVector msnu;
  /// chargino/neutralino masses: ordered
   DoubleVector mch, mneut, mhiggs;
  /// Gluino mass
  double mGluino;
  /// neutralino mixing
   DoubleMatrix mixNeut, mixHiggs;
  /// chargino and third family mixing angles
  double thetaL, thetaR, thetat, thetab, thetatau;
  /// sparticle masses in order (i=L/R, family)
  DoubleMatrix mu, md, me;
   /// Higgs mixing angle (alpha)
   double thetaH;
   //PA: new mixing angle.
   double thetaHP;
  
  /// DRbar tadpoles evaluated at MSusy: these tadpoles include the 2-loop
  /// pieces
  double t1OV1Ms, t2OV2Ms;

  /// These DRbar tadpoles evaluated at MSusy only include 1-loop pieces. One
  /// must keep them in order to get the Higgs mass calculations correct to
  /// two loops because of the way that they are organised.
  double t1OV1Ms1loop, t2OV2Ms1loop;
  
  /// Feynman rules
  //  DoubleMatrix aChi0ChicW, bChi0ChicW;
  
  nmsPhysical(); ///< Constructor: initialises with zeroes
  nmsPhysical(const nmsPhysical &); ///< Constructor copies another object
  /// Displays the object in situations where a const is required
  const nmsPhysical & displaynmsPhysical() const { return *this; };
  /// Sets whole contents to those of another object s
  void setnmsPhysical(const nmsPhysical &s) { *this = s; };
  /// Sets whole contents to those of another object s
  const nmsPhysical & operator = (const nmsPhysical & s);
  
  /// Displays contents in a C-style convention *a (starts with index zero)
  void display(double *a) const;
}; 

/// Formatted printout
std::ostream & operator <<(std::ostream &, const nmsPhysical &); 

/// Various boolean values to flag any problems with the parameter point
struct NMsProblem {
  bool mgutOutOfBounds;
  bool badConvergence; ///< Nowhere near a decent solution
  bool irqfp; ///< Infra-red quasi fixed point breached
  bool noRhoConvergence; ///< Couldn't calculate electroweak rho parameter
  bool noConvergence; ///< Iteration did not converge: not always serious
  tachyonType tachyon; ///< Tachyonic point
  bool muSqWrongSign; ///< mu^2 came out with wrong sign; no REWSB
  bool m3sq; ///< m3sq came out with wrong sign; no REWSB
  bool higgsUfb; ///< Higgs potential inconsistent with a good minimum
  bool nonperturbative; ///< Running went non-perturbative
  bool noMuConvergence; ///< mu couldn't be calculated
  /// Higgs mass is potentially inaccurate and cant be trusted
  bool inaccurateHiggsMass; 
  bool problemThrown;
   bool szero;
   /// returns true if there's any problem 
  bool test() const 
  {return (mgutOutOfBounds || irqfp || noConvergence || tachyon || 
	   muSqWrongSign || higgsUfb || nonperturbative || noRhoConvergence || 
	   noMuConvergence || m3sq || badConvergence || inaccurateHiggsMass ||
	   problemThrown);}; 
  /// Only returns true if there's a serious problem
  bool testSeriouNMsProblem() const 
  {return (irqfp || tachyon || muSqWrongSign || higgsUfb || nonperturbative 
	   || noRhoConvergence || noMuConvergence || m3sq || badConvergence ||
	   mgutOutOfBounds || problemThrown);}; 

  inline NMsProblem(); ///< constructor full of false values
  /// Constructor that sets flags equal to those of s
  inline NMsProblem(const NMsProblem & s);
  /// Sets flags equal to those of s
  const NMsProblem & operator = (const NMsProblem &);
  
};
/// Formatted output, but won't print unflagged problems
std::ostream & operator <<(std::ostream &st, const NMsProblem & p);
/// Formatted input of physical parameters
std::istream & operator >>(std::istream & left, nmsPhysical &s);

/// DRbar values of masses and mixings in MSSM
struct NMdrBarPars: public nmsPhysical { 
  double mz, mw;       /// Running electroweak gauge boson masses
  double mt, mb, mtau; /// Running top, bottom and tau mass
  double ht, hb, htau; /// Running diagonal top, bottom, tau Yukawa couplings
  double ut, ub, utau; /// Running diagonal top, bottom, tau trlinear couplings
  /// BPMZ convention mixing matrices for neutralinos and charginos
  ComplexMatrix nBpmz, uBpmz, vBpmz; 
  /// positive definite masses for neutralinos and charginos
  DoubleVector mnBpmz, mchBpmz; 

  inline NMdrBarPars(); ///< Initialises with zero values
  inline NMdrBarPars(const NMdrBarPars &); ///< Initialises with another object
  /// Sets contents equal to those of another object
  const NMdrBarPars & operator = (const NMdrBarPars &s);

  /// Returns mixing matrix o and neutralino masses mn in the MPZ convention
  /// (hep-ph/9606211), n is 4 by 4 and mneut is 1->4, ie
  /// Calculates mnBpmz, nBpmz
  void mpzNeutralinos(); 
  /// Calculates uBpmz, vBpmz, mchBpmz, ie
  /// Returns mixing matrices u,v and neutralino masses mneut in the MPZ
  /// convention (hep-ph/9606211),  u+v are (2,2) and mch is 1->2.
  void mpzCharginos();
};

std::ostream & operator <<(std::ostream &, const NMdrBarPars &); 

// ---------------------- inline class members ------------------------
inline NMdrBarPars::NMdrBarPars()
  : nmsPhysical(), mz(0.0), mw(0.0), 
    mt(0.0), mb(0.0), mtau(0.0), ht(0.0), hb(0.0), htau(0.0), 
    ut(0.0), ub(0.0), utau(0.0),
    nBpmz(4, 4), uBpmz(2, 2), 
    vBpmz(2, 2), mnBpmz(4), mchBpmz(2)
{}

inline NMdrBarPars::NMdrBarPars(const NMdrBarPars &s)
  : nmsPhysical(s.displaynmsPhysical()), mz(s.mz), mw(s.mw),
    mt(s.mt), mb(s.mb), mtau(s.mtau), ht(s.ht), hb(s.hb), htau(s.htau),
    ut(s.ut), ub(s.ub), utau(s.utau),
    nBpmz(s.nBpmz), uBpmz(s.uBpmz), vBpmz(s.vBpmz), mnBpmz(s.mnBpmz), 
    mchBpmz(s.mchBpmz)
{}

inline nmsPhysical::nmsPhysical()
  : mh0(2), mA0(1), mHpm(0.), mH1(0.), mH2(0.), mA1(0.), mA2(0.),
    msnu(3), mch(2), mneut(5), mhiggs(3), mGluino(0.0),
    mixNeut(5, 5), mixHiggs(3,3), 
    thetaL(0.0), thetaR(0.0), thetat(0.0), thetab(0.0),
    thetatau(0.0), mu(2, 3), md(2, 3), me(2, 3), thetaH(0.0), 
    t1OV1Ms(0.0), t2OV2Ms(0.0), t1OV1Ms1loop(0.), t2OV2Ms1loop(0.)
{}

inline nmsPhysical::nmsPhysical(const nmsPhysical & s)
   : mh0(s.mh0), mA0(s.mA0), mHpm(s.mHpm), 
    mH1(s.mH1), mH2(s.mH2), mA1(s.mA1), mA2(s.mA2), 
    msnu(s.msnu), mch(s.mch), 
    mneut(s.mneut), mhiggs(s.mhiggs), mGluino(s.mGluino), 
    mixNeut(s.mixNeut), mixHiggs(s.mixHiggs), thetaL(s.thetaL),
    thetaR(s.thetaR), thetat(s.thetat), thetab(s.thetab),
    thetatau(s.thetatau), mu(s.mu), md(s.md), me(s.me), thetaH(s.thetaH),
    t1OV1Ms(s.t1OV1Ms), t2OV2Ms(s.t2OV2Ms),    
    t1OV1Ms1loop(s.t1OV1Ms1loop), t2OV2Ms1loop(s.t2OV2Ms1loop)
{}

inline NMsProblem::NMsProblem()
  : mgutOutOfBounds(false), badConvergence(false), 
    irqfp(false), noRhoConvergence(false), noConvergence(false),
    tachyon(none), muSqWrongSign(false), m3sq(false), higgsUfb(false), 
    nonperturbative(false), noMuConvergence(false),     
    inaccurateHiggsMass(false), problemThrown(false)
{}

inline NMsProblem::NMsProblem(const NMsProblem & s)
  : mgutOutOfBounds(s.mgutOutOfBounds), badConvergence(s.badConvergence), 
    irqfp(s.irqfp), noRhoConvergence(s.noRhoConvergence), 
    noConvergence(s.noConvergence),
    tachyon(s.tachyon), muSqWrongSign(s.muSqWrongSign), m3sq(s.m3sq),
    higgsUfb(s.higgsUfb), nonperturbative(s.nonperturbative), 
    noMuConvergence(s.noMuConvergence), 
    inaccurateHiggsMass(s.inaccurateHiggsMass), problemThrown(s.problemThrown)
{}

#endif
