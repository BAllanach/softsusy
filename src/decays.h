/** \file decays.h
   - Project:     SOFTSUSY 
   - Author:      Tom Cridge, Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Code calculates decay modes and prints out an SLHA format file with 
          them in. For R-parity conserving NMSSM/MSSM. See arXiv:1703.09717
*/

#ifndef DECAYS_H
#define DECAYS_H

#include "nmssmsoftsusy.h"
#include "softsusy.h"
#include "physpars.h"
#include "lowe.h"
#include "def.h"
#include "softpars.h"
#include "softsusy.h"
#include "flavoursoft.h"
#include "susy.h"
#include "particle.h"
#include "twoBodyDecays.h"
#include "threeBodyDecays.h"
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include <complex>

using namespace std;

/// Calculate Decays does all decay table calculations and outputs
int calculateDecays(ostream & out, MssmSoftsusy * r,
		    const NmssmSoftsusy & nmssm, bool nmssmIsIt);
  
#endif
