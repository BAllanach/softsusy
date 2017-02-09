/** \file particle.h
   - Project:     SOFTSUSY 
   - Author:      Tom Cridge, Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Definitions of particles and container for decay widths.
*/

#ifndef PARTICLE_H
#define PARTICLE_H

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
  double No_NMSSM_Decays;
  double total_width;
  double two_width;
  double three_width;
  vector<vector<double> > Array_Decays;
  vector<string> Array_Comments;
};

#endif
