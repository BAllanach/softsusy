/** \file mssm_twoloop_mt.h
    - Project:     SOFTSUSY 
    - Author:      Alex Voigt
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
    
    \brief Two loop corrections to the top mass.
    This file has been generated at Fri 31 Mar 2017 14:09:29
    with the script "tquark_to_cpp.m". This file is part of FlexibleSUSY.
    
*/


#ifndef MSSM_TWO_LOOP_SQCD_MT_H
#define MSSM_TWO_LOOP_SQCD_MT_H

#include "numerics.h"
#include "def.h"
#include <cmath>
#include <limits>
 
namespace softsusy {

  struct mtParameters {
    mtParameters() = default;
    mtParameters(double g3_, double mt_, double mg_, double mst1_,
               double mst2_, double msusy_, double xt_, double Q_)
      : g3(g3_), mt(mt_), mg(mg_), mst1(mst1_)
      , mst2(mst2_), msusy(msusy_), xt(xt_), Q(Q_)
    {}
    
    double g3{};    ///< MSSM strong gauge coupling DR-bar
    double mt{};    ///< MSSM top mass DR-bar
    double mg{};    ///< MSSM gluino mass DR-bar
    double mst1{};  ///< MSSM light stop mass DR-bar
    double mst2{};  ///< MSSM heavy stop mass DR-bar
    double msusy{}; ///< MSSM SUSY particle mass scale
    double xt{};    ///< MSSM stop mixing parameter DR-bar
    double Q{};     ///< renormalization scale
  };
  
  /// 2-loop QCD contributions to Delta Mt over mt [hep-ph/0507139]
  double dMt_over_mt_2loop_qcd(const mtParameters&);
  /// 2-loop SUSY contributions to Delta Mt over mt [hep-ph/0507139]
  double dMt_over_mt_2loop_susy(const mtParameters&);
  /// 2-loop full SQCD contributions to Delta Mt over mt [hep-ph/0507139].
  /// Here, msusy should be mu(1, 2)
  double dMt_over_mt_2loop(double g3, double mt, double mg, double mst1,
			   double mst2, double msusy, double thetat,
			   double q);
  
} ///< namespace softsusy

#endif
