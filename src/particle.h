/** \file particle.h
   - Project:     SOFTSUSY 
   - Author:      Tom Cridge, Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Definitions of particles and container for decay widths.
*/

#ifndef PARTICLE_H
#define PARTICLE_H

/// Particle class definition for decays
class Particle {
 public:
  string name;                             ///< name of the parent
  double mass;                             ///< mass of the parent
  double PDG;                              ///< PDG code of parent
  double No_of_Decays;                     ///< How many decay modes
  double No_1to2_Decays;                   ///< How many 2-body decay modes
  double No_1to3_Decays;                   ///< How many 3-body decay modes
  double No_grav_Decays;                   ///< How many decays to gravitinos
  double No_NMSSM_Decays;                  ///< How many NMSSM specific modes
  double total_width;                
  double two_width;                        ///< Total width of 2-body modes
  double three_width;                      ///< Total width of 3-body modes
  vector<vector<double> > Array_Decays;    ///< Partial width for each decay
  vector<string> Array_Comments;           ///< Comments for SLHA file
};

#endif
