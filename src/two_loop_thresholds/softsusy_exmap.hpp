#include "softsusy.h"

#ifndef _SOFTSUSY_EXMAP_H
#define _SOFTSUSY_EXMAP_H
//#include <ginac/ginac.h>
//#include "mssmparam.hpp"
//#include "softsusy.h"
//#include "utils.h"


//class MssmSoftsusy;
//using namespace softsusy;

//#include <softsusy.h>

namespace SoftSusy_helpers_ {

using namespace GiNaC;

 /// Return GiNac exmap with MSSM (Yukawa Model) parameters
  GiNaC::exmap drBarPars_exmap(const MssmSoftsusy & mssy_);

} // namespace SoftSusy_helpers_
#endif

// vi:ts=2:sw=2

