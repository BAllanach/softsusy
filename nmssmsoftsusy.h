/** \file nmssmsoftsusy.cpp
    Project: NMSSMSOFTSUSY 
    Author: Ben Allanach, Peter Athron, Lewis Tunstall, Alexander Voigt 
    Manual: TBW
    Webpage:  https://github.com/Expander/softsusy.git 
*/


#ifndef NMSSMSOFTSUSY_H
#define NMSSMSOFTSUSY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <def.h>
#include <utils.h>
#include <numerics.h>
#include <physpars.h> 
#include <lowe.h>
#include <nmssmsoftpars.h>
#include <softsusy.h>
#include "bcs.h"
#include "mssmUtils.h"
//#include <nmssm2loop.h>
using namespace softsusy;
using namespace std;
/* class NmssmSoftsusy;  */
/* std::istream & operator >>(std::istream &left, NmssmSoftsusy &s);\ */

/// Contains all supersymmetric NMSSM parameters, incorporating R_p NMSSM
class NmssmSoftsusy: public Softsusy<SoftParsNmssm> {
private:
  
   double tSOVSMs;  ///< New Nmssm DRbar tadpole(MSusy): incl 2 loops
   double tSOVSMs1loop; ///<New Nmssm DRbar tadpole(MSusy): excl 2 loops

public:
//  void (*boundaryCondition)(NmssmSoftsusy &, const DoubleVector &);
  /// Default constructor fills object with zeroes
  NmssmSoftsusy();
  /// Constructor sets SUSY parameters only from another object
  NmssmSoftsusy(const NmssmSusy &);
  /// Constructor copies another object
  NmssmSoftsusy(const NmssmSoftsusy &);
  /// Sets all parameters from s, sp, mu is the mu superpotential parameter, l
  /// is the number of loops used for RG evolution, t is the thresholds
  /// accuracy parameter, mg is the gravitino mass, hv is the Higgs VEV
  /// parameter.
  NmssmSoftsusy(const SoftParsNmssm & s, const sPhysical & sp, double mu, int l, int t, double hv);
  /// Set all data in the object equal to another
  const NmssmSoftsusy & operator=(const NmssmSoftsusy & s);
  
 
  /// Organises tree-level calculation of all sparticle masses and mixings
  virtual void calcDrBarPars();
  
  /// Returns tree-level up squark mass matrix in "mass" for NMSSM.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mtrun=DR bar top
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
void treeUpSquark(DoubleMatrix & mass, double mtrun, 
                    double pizztMS, double sinthDRbarMS, 
                    int family);

/// Returns tree-level down squark mass matrix in "mass" for NMSSM.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mbrun=DR bar bottom
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
  void treeDownSquark(DoubleMatrix & mass, double mbrun, double pizztMS, 
		double sinthDRbarMS, int family);
  /// Returns tree-level down squark mass matrix in "mass" for NMSSM.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mTrun=DR bar tau
  /// mass, family=generation of slepton, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w
  void treeChargedSlepton(DoubleMatrix & mass, double mTrun, double pizztMS, 
		double sinthDRbarMS, int family);
  /// LCT: new routine to allocate NMSSM chargino masses, 
  //Returns tree-level chargino mass matrix in the NMSSM 
  void calcDrBarCharginos(DoubleMatrix & mass, double beta, double mw);
  
/// LCT: new routine for NMSSM neutralino masses, 
  //Returns tree-level Neutralino mass matrix in the NMSSM 
  void calcDrBarNeutralinos(DoubleMatrix & mass, double beta, double mz, double mw, double sinthDRbar);

  void calcDrBarHiggs(double beta, double mz2, double mw2, double sinthDRbar, drBarPars & eg);
 
 //PA: A print method used in development.  I find it useful and easier to read than couting the normal display function or calling printlong etc.    
  void printall();

};
inline NmssmSoftsusy::NmssmSoftsusy()
		     : Softsusy(), tSOVSMs(0.0), tSOVSMs1loop(0.0)  {}


inline NmssmSoftsusy::NmssmSoftsusy(const NmssmSoftsusy & s)
		     : Softsusy<SoftParsNmssm>(s),
 tSOVSMs1loop(s.tSOVSMs1loop ), tSOVSMs(s.tSOVSMs) {
    setPars(121);   
}


inline NmssmSoftsusy::NmssmSoftsusy(const NmssmSusy &s)
  : Softsusy<SoftParsNmssm>(s),tSOVSMs(0.0), tSOVSMs1loop(0.0)  {
      setPars(121);
}


inline NmssmSoftsusy::NmssmSoftsusy
(const SoftParsNmssm & s, const sPhysical & sp, double mu, int l, int t,
 double hv): Softsusy<SoftParsNmssm>(s, sp, mu, l, t, hv),tSOVSMs(0.0), tSOVSMs1loop(0.0)  {
 setPars(121);
  
}

#endif
