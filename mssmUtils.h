
#ifndef MSSM_UTILS_H
#define MSSM_UTILS_H

#include <iosfwd>
#include <string>

class sBrevity;
class MssmSusy;

template<class SoftPars> class Softsusy;
template<class Susy, class Brevity> class SoftPars;

typedef SoftPars<MssmSusy, sBrevity> SoftParsMssm;
typedef Softsusy<SoftParsMssm> MssmSoftsusy;

/// Pietro's fit to the totality of LEP2 data: gives whether the lightest
/// Higgs is allowed or not - depends only upon the mass and beta-alpha.
/// true means allowed, false means ruled out. 
bool testLEPHiggs(const MssmSoftsusy & r, double error = 3.0);

/// Given mu paramer, tau Yukawa htau, family number examined, finds height of
/// potential for temp at |H_2|=h2.
double ufb3fn(double mu, double htau, double h2, int family, const MssmSoftsusy
	      & temp);

/// For UFB-3direction, returns scale at which one-loop corrections are
/// smallest. IO parameters: inminTol is fractional accuracy with which
/// minimum is found, eR is value of RH selectron field, h2 is value of H2
/// field, Lisq=|L_i|^2 slepton VEV value, mx=high BC-scale
double getQhat(double inminTol,double eR, double h2, double Lisq, double mx,
	       MssmSoftsusy & temp);

/// Calculates fractional difference in Drbar masses between in and out
double sumTol(const MssmSoftsusy & in, const MssmSoftsusy & out, int numTries);

/// Prints out the identity of the LSP to standard output. 
/// temp, j are defined from lsp function 
std::string recogLsp(int temp, int j);

/// Formatted output
std::ostream & operator<<(std::ostream&, const MssmSoftsusy&); 

#endif
