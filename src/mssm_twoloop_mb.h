/** \file mssm_twoloop_mb.h
    - Project:     SOFTSUSY 
    - Author:      Alex Voigt
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
    
    \brief Two loop corrections to the bottom mass.
    This file has been generated at Fri 7 Apr 2017 21:00:53
    with the script "bquark_to_cpp.m". It is part of FlexibleSUSY.    
*/

#ifndef MSSM_TWO_LOOP_SQCD_MB_H
#define MSSM_TWO_LOOP_SQCD_MB_H

#include "numerics.h"
#include "def.h"
#include <array>
#include <cmath>
#include <limits>

namespace softsusy {

  struct Parameters {
    Parameters() = default;
    Parameters(double g3_, double mt_, double mb_, double mg_,
               double mst1_, double mst2_,
               double msb1_, double msb2_,
               double msusy_,
               double xt_, double xb_, double Q_)
      : g3(g3_), mt(mt_), mb(mb_), mg(mg_)
      , mst1(mst1_), mst2(mst2_)
      , msb1(msb1_), msb2(msb2_)
      , msusy(msusy_)
      , xt(xt_), xb(xb_), Q(Q_)
    {}
    
    double g3{};    ///< MSSM strong gauge coupling DR-bar
    double mt{};    ///< MSSM top mass DR-bar
    double mb{};    ///< SM   bottom mass DR-bar
    double mg{};    ///< MSSM gluino mass DR-bar
    double mst1{};  ///< MSSM light stop mass DR-bar
    double mst2{};  ///< MSSM heavy stop mass DR-bar
    double msb1{};  ///< MSSM light sbottom mass DR-bar
    double msb2{};  ///< MSSM heavy sbottom mass DR-bar
    double msusy{}; ///< MSSM light squark masses DR-bar
    double xt{};    ///< MSSM sbottom mixing parameter DR-bar
    double xb{};    ///< MSSM stop mixing parameter DR-bar
    double Q{};     ///< renormalization scale
  };
  
  /// 2-loop full SQCD contributions to mb [arXiv:0707.0650]
  double delta_mb_2loop(const Parameters&);
  /// 2-loop full SQCD contributions to mb [arXiv:0707.0650]
  double delta_mb_2loop(double g3, double mt, double mb, double mg,
			double mst1, double mst2, double msb1, double msb2,
			double msusy, double thetat, double thetab, double q);
  
} ///< namespace softsusy

#endif
