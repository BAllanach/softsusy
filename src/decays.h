/** \file decays.h
   - Project:     SOFTSUSY 
   - Author:      Tom Cridge, Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Code calculates decay modes and prints out an SLHA format file with 
          them in 
*/

#ifndef DECAYS_H
#define DECAYS_H

#include "nmssmsoftsusy.h"
#include "softsusy.h"
//#include "Particle.h"
#include "physpars.h"
#include "lowe.h"
#include "softpars.h"
#include "softsusy.h"
#include "flavoursoft.h"
#include "susy.h"
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <stdlib.h>
#include <vector>


using namespace std;

/// Particle class definition
class Particle {
 public:
  string name;
  double mass;
  double PDG;
  double No_of_Decays;
  double No_1to2_Decays;
  double No_1to3_Decays;
  double No_grav_Decays;
  double total_width;
  double two_width;
  double three_width;
  vector<vector<double> > Array_Decays;
  vector<string> Array_Comments;
};

void calculateDecays(MssmSoftsusy * r);
void calculateDecaysNmssm(const NmssmSoftsusy & r);

//static MssmSoftsusy *useInDGauss;

double fdgauss(double x);

#endif
