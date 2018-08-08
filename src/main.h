
/** \file main.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief a main C++ program to calculate Higgs masses as a function of tan
   beta 
*/

/** \mainpage Detailed SOFTSUSY Documentation 

    \section install Installation or downloads
    For installation instructions or a download, please go to the 
    <a href="http://projects.hepforge.org/softsusy/">
    SOFTSUSY Homepage</a>

    \section manual Official manuals
    If you use SOFTSUSY to write a paper, please cite <br>
    [Allanach:2001kg] B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145,

    which is the SOFTSUSY manual for the R-parity conserving MSSM. 
    If you use the decay calculations, please cite [Allanach:2001kg] and <br>
    [Allanach:2017hcf] B.C. Allanach and T. Cridge, Comput. Phys. Commun. 220 (2017) 417, arxiv:1703.09717

    If you calculate in the NMSSM, please cite [Allanach:2001kg] and <br>
    [Allanach:2013kza] B.C. Allanach, P. Athron, L. Tunstall, A. Voigt and A. Williams, Comput. Phys. Comm. 185 (2014) 2322, arXiv:1311.7659.

    If you use the R-parity violating aspects, please cite [Allanach:2001kg] and <br>    
    [Allanach:2009bv] B.C. Allanach and M.A. Bernhardt, Comput. Phys. Commun. 181 (2010) 232, arXiv:0903.1805.
    
    If you use it to calculate neutrino masses and mixings, please cite
    [Allanach:2001kg], [Allanach:2009bv] and <br>
    [Allanach:2011de] B.C. Allanach, M. Hanussek and C.H. Kom, Comput. Phys. Commun. 183 (2012) 785, arXiv:1109.3735.
    
    If you use the three-loop RGEs or two-loop threshold corrections to
    gauge/Yukawa couplings, please cite [Allanach:2001kg] and <br>
    [Allanach:2014nba] B.C. Allanach, A. Bednyakov and R. Ruiz de Austri, Comput. Phys. Commun. (2015) 192, arXiv:1407.6130.
    
    If you use the two-loop SUSY QCD corrections to squark and gluino pole masses, please cite [Allanach:2001kg] and<br>
    [Allanach:2016rxd] B.C. Allanach, Stephen P. Martin, David G. Robertson
    and Roberto Ruiz de Austri, Comput.Phys.Commun. 219 (2017) 339, 
    arXiv:1601.06657.
    
    \section documentation Documentation
    These web-pages contain the documentation of the latest SOFTSUSY code.
    There are class diagrams and cross-referenced links a la doxygen to help 
    you navigate.

    \section updates Official Updates
    Updates will be posted on the    
    <a href="https://softsusy.hepforge.org/">
    SOFTSUSY Homepage</a>, and the manuals in the doc subdirectory
    will also be updated and released with the distribution.
 */

#include <iostream>
#include <string>
#include "mycomplex.h"
#include "linalg.h"
#include "lowe.h"
#include "rge.h"
#include "softsusy.h"
#include "softpars.h"
#include "susy.h"
#include "utils.h"
#include "numerics.h"
using namespace softsusy;

namespace NR {
  extern int nn;
  extern DoubleVector fvec;
}

  extern void (*nrfuncv)(int n, DoubleVector v, DoubleVector & f);

