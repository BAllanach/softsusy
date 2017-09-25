/** \file mssm_twoloop_as.h
    - Project:     SOFTSUSY 
    - Author:      Ben Allanach and Alex Voigt
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
    
    \brief Two loop corrections to the top mass.
    This file has been generated at Fri 14 Jul 2017 22:12:37
    with the script "as2_to_cpp.m".
*/


#ifndef MSSM_TWO_LOOP_AS_H
#define MSSM_TWO_LOOP_AS_H

#include <iosfwd>

namespace flexiblesusy {
namespace mssm_twoloop_as {

using Real = long double;
  
struct Parameters {
    Real g3{};    ///< MSSM strong gauge coupling DR-bar
    Real yt{};    ///< MSSM top Yukawa coupling DR-bar
    Real yb{};    ///< MSSM bottom Yukawa coupling DR-bar
    Real mt{};    ///< MSSM top mass DR-bar
    Real mb{};    ///< MSSM bottom mass DR-bar
    Real mg{};    ///< MSSM gluino mass DR-bar
    Real mst1{};  ///< MSSM light stop mass DR-bar
    Real mst2{};  ///< MSSM heavy stop mass DR-bar
    Real msb1{};  ///< MSSM light sbottom mass DR-bar
    Real msb2{};  ///< MSSM heavy sbottom mass DR-bar
    Real msd1{};  ///< MSSM light sdown mass DR-bar
    Real msd2{};  ///< MSSM heavy sdown mass DR-bar
    Real xt{};    ///< MSSM stop mixing parameter DR-bar
    Real xb{};    ///< MSSM sbottom mixing parameter DR-bar
    Real mw{};    ///< MSSM W boson mass DR-bar
    Real mz{};    ///< MSSM Z boson mass DR-bar
    Real mh{};    ///< MSSM light CP-even Higgs mass DR-bar
    Real mH{};    ///< MSSM heavy CP-even Higgs mass DR-bar
    Real mC{};    ///< MSSM charged Higgs mass DR-bar
    Real mA{};    ///< MSSM CP-odd Higgs mass DR-bar
    Real mu{};    ///< MSSM mu superpotential parameter DR-bar
    Real tb{};    ///< MSSM tan(beta) DR-bar
    Real Q{};     ///< renormalization scale
};

/// 2-loop O(alpha_s^2) contributions to Delta alpha_s [hep-ph/0509048,arXiv:0810.5101]
Real delta_alpha_s_2loop_as_as(const Parameters&);

/// 2-loop O(alpha_t*alpha_s) contributions to Delta alpha_s [arXiv:1009.5455]
Real delta_alpha_s_2loop_at_as(const Parameters&);

/// 2-loop O(alpha_b*alpha_s) contributions to Delta alpha_s [arXiv:1009.5455]
Real delta_alpha_s_2loop_ab_as(const Parameters&);

std::ostream& operator<<(std::ostream&, const Parameters&);

} // namespace mssm_twoloop_as
} // namespace flexiblesusy

#endif
