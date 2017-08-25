/** \file particle.h
   - Project:     SOFTSUSY 
   - Author:      Tom Cridge, Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Definitions of particles and container for decay widths.
*/

#ifndef PARTICLE_H
#define PARTICLE_H

/// Define Particle PDG codes - for SLHA output of decay tables
const int PDGdown = 1, PDGup = 2, PDGstrange = 3, PDGcharm = 4, PDGbottom = 5,
  PDGtop = 6;
const int PDGelectron = 11, PDGnuelectron = 12, PDGmuon = 13, PDGnumuon = 14,
  PDGtau = 15, PDGnutau = 16, PDGgluon = 21, PDGphoton = 22, PDGZboson = 23,
  PDGWplus = 24, PDGh0 = 25, PDGH0 = 35, PDGA0 = 36, PDGHplus = 37;
const int PDGsdownL = 1000001, PDGsupL = 1000002, PDGsstrangeL = 1000003,
  PDGscharmL = 1000004, PDGsbottom1 = 1000005, PDGstop1 = 1000006;
const int PDGselectronL = 1000011, PDGnuselectronL = 1000012,
  PDGsmuonL = 1000013, PDGnusmuonL = 1000014, PDGstau1 = 1000015,
  PDGnustauL = 1000016;
const int PDGgluino = 1000021, PDGneutralino1 = 1000022,
  PDGneutralino2 = 1000023, PDGchargino1 = 1000024, PDGneutralino3 = 1000025,
  PDGneutralino4 = 1000035, PDGchargino2 = 1000037;
const int PDGsdownR = 2000001, PDGsupR = 2000002, PDGsstrangeR = 2000003,
  PDGscharmR = 2000004, PDGsbottom2 = 2000005, PDGstop2 = 2000006;
const int PDGselectronR = 2000011, PDGsmuonR = 2000013, PDGstau2 = 2000015;
/// const int PDGnuselectronR = 2000012, PDGnusmuonR = 2000014, PDGnustauR =
/// 2000016 - for use later with right handed neutrinos
const int PDGgravitino = 1000039;
/// useful meson codes
const int PDGpi0 = 111, PDGpiPlus = 211;

///PDG codes for extra NMSSM particles:
const int PDGA2 = 46, PDGH3 = 45, PDGneutralino5 = 1000045;

typedef enum {gluino = 0, sdownL, sdownR, supL, supR, strangeL, strangeR,
	      scharmL, scharmR, sbottom1, sbottom2, stop1, stop2, selectronL,
	      selectronR, smuonL, smuonR, snue, snumu, stau1, stau2, snutau,
	      chargino1, chargino2, neutralino1, neutralino2, neutralino3,
	      neutralino4, h0, H0, A0, Hplus, A2, H3, neutralino5} particleType;

/// Particle class definition for decays
class Particle {
 public:
  string name;                             ///< name of the parent
  double mass;                             ///< mass of the parent
  int PDG;                              ///< PDG code of parent
  int No_of_Decays;                     ///< How many decay modes
  int No_1to2_Decays;                   ///< How many 2-body decay modes
  int No_1to3_Decays;                   ///< How many 3-body decay modes
  int No_1to4_Decays;                   ///< How many 4-body decay modes
  int No_grav_Decays;                   ///< How many decays to gravitinos
  int No_NMSSM_Decays;                  ///< How many NMSSM specific modes
  double total_width;                
  double two_width;                        ///< Total width of 2-body modes
  double three_width;                      ///< Total width of 3-body modes
  double four_width;                      ///< Total width of 4-body modes  
  vector<vector<double> > Array_Decays;    ///< Partial width for each decay
  vector<string> Array_Comments;           ///< Comments for SLHA file
};

#endif
