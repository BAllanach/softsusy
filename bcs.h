
#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

class DoubleVector;
class sBrevity;
class MssmSusy;

template<class SoftPars> class Softsusy;
template<class Susy, class Brevity> class SoftPars;

typedef SoftPars<MssmSusy, sBrevity> SoftParsMssm;
typedef Softsusy<SoftParsMssm> MssmSoftsusy;

/// Allows user to specify a boundary condition where ALL SUSY breaking
/// parameters are specified in inputParameters
void generalBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Sets all soft parameters in m except for mh1sq or mh2sq: it is intended
/// for the case where mu and M_A^0(pole) is specified
void generalBcs2(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// non-universal mSUGRA boundary conditions including mH1^2 and mH2^2
void extendedSugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// User supplied routine. Inputs m at the unification scale, and uses
/// inputParameters vector to output m with high energy soft boundary
/// conditions. 
void sugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Non-universal higgs mass conditions. Paramaters are, in order: m0,m12,mH,a0
void nuhmI(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Non-universal higgs mass conditions. Paramaters are, in order:
/// m0,m12,mH1,mH2,a0
void nuhmII(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Adds 2-loop AMSB boundary conditions onto m
void amsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Large Volume string compactification boundary conditions
void lvsBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// One-loop GMSB boundary conditions
void gmsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// For the user to define....
void userDefinedBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);

void nonUniGauginos(MssmSoftsusy & m, const DoubleVector & inputParameters);

void splitGmsb(MssmSoftsusy & m, const DoubleVector & inputParameters);

#endif
