
/** \file mainDecay.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief a main C++ program to calculate Higgs masses as a function of tan
   beta 
*/

#include <iostream>
#include "mycomplex.h"
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "rge.h"
#include "softsusy.h"
#include "softpars.h"
#include "susy.h"
#include "utils.h"
#include "numerics.h"
#include "nmssmsoftsusy.h"
#include "decays.h"
using namespace softsusy;

namespace NR {
  extern int nn;
  extern DoubleVector fvec;
}

  extern void (*nrfuncv)(int n, DoubleVector v, DoubleVector & f);

