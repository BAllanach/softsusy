
/** \file softpoint.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach, Markus Bernhardt 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: main calling program: command line interface
   \brief Main program: reads input in SLHA or command line format
*/

#include <iostream>
#include <sstream>
#include <cstring>
#include <src/mycomplex.h>
#include <src/def.h>
#include <src/linalg.h>
#include <src/lowe.h>
#include <src/rge.h>
#include <src/softsusy.h>
#include <src/flavoursoft.h>
#include <src/nmssmsoftsusy.h>
#include <src/softpars.h>
#include <src/physpars.h>
#include <src/susy.h>
#include <src/utils.h>
#include <src/numerics.h>
#include <src/twoloophiggs.h>
#include <src/dilogwrap.h>
#include <src/rpvneut.h>
using namespace softsusy;

/// string routine for options
bool starts_with(const std::string& str,
		 const std::string& prefix) {
  return !str.compare(0, prefix.size(), prefix);
}

double get_value(const std::string& str, const std::string& prefix) {
   return atof(str.substr(prefix.size()).c_str());
}

int get_valuei(const std::string& str, const std::string& prefix) {
   return atoi(str.substr(prefix.size()).c_str());
}

bool contains_only_whitespace(const std::string& str) {
  return str.find_first_not_of(" \t") == string::npos;
}

namespace softsusy {
   /// Requested by CMS
   void splitGmsb(MssmSoftsusy & m, const DoubleVector & inputParameters);
}

/// Does the user require gauge unification or not -- gaugeUnification changed
/// to be correct value
inline double mgutCheck(char * a, bool & gaugeUnification, 
			bool & ewsbBCscale) { 
  gaugeUnification = false; ewsbBCscale = false;
  if (!strcmp(a, "--mgut=?") || !strcmp(a,"--mgut=unified")) {
    gaugeUnification = true; 
    return 2.0e16;
  }
  if (!strcmp(a, "--mgut=msusy")) {
    ewsbBCscale = true;
    return 1.0e3;
  }
  else return get_value(a, "--mgut");
}

/// Incorrect input: gives advice on how to supply it
void errorCall();
