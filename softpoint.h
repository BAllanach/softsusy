
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
#include <mycomplex.h>
#include <def.h>
#include <linalg.h>
#include <lowe.h>
#include <rge.h>
#include <softsusy.h>
#include <flavoursoft.h>
#include <nmssmsoftsusy.h>
#include <softpars.h>
#include <physpars.h>
#include <susy.h>
#include <utils.h>
#include <numerics.h>
#include <twoloophiggs.h>
#include <dilogwrap.h>
#include <rpvneut.h>
#include <cassert>
using namespace softsusy;

/// Requested by CMS
void splitGmsb(MssmSoftsusy & m, const DoubleVector & inputParameters);

/// Does the user require gauge unification or not -- gaugeUnification changed
/// to be correct value
inline double mgutCheck(char * a, bool & gaugeUnification, 
			bool & ewsbBCscale) { 
  gaugeUnification = false; ewsbBCscale = false;
  if (!strcmp(a, "?") || !strcmp(a,"unified")) {
    gaugeUnification = true; 
    return 2.0e16;
  }
  if (!strcmp(a, "msusy")) {
    ewsbBCscale = true;
    return 1.0e3;
  }
  else return atof(a);
}

/// Incorrect input: gives advice on how to supply it
void errorCall();

/// class for NMSSM input parameters
class NMSSM_input {
public:
   enum NMSSM_parameters {
      tanBeta,               // MINPAR entry 3
      mHd2,                  // EXTPAR entry 21
      mHu2,                  // EXTPAR entry 22
      mu,                    // EXTPAR entry 23
      BmuOverCosBetaSinBeta, // m3^2/(cos(Beta)sin(Beta)) EXTPAR entry 24
      lambda,                // EXTPAR entry 61
      kappa,                 // EXTPAR entry 62
      Alambda,               // EXTPAR entry 63
      Akappa,                // EXTPAR entry 64
      lambdaS,               // lambda * <S> EXTPAR entry 65
      xiF,                   // EXTPAR entry 66
      xiS,                   // EXTPAR entry 67
      muPrime,               // EXTPAR entry 68
      mPrimeS2,              // EXTPAR entry 69
      mS2,                   // EXTPAR entry 70
      NUMBER_OF_NMSSM_INPUT_PARAMETERS
   };

   NMSSM_input()
      : parameter()    // sets all parameters to zero
      , has_been_set() // sets all values to zero (false)
   {}

   /// set parameter to given value
   void set(NMSSM_parameters par, double value) {
      assert(par < NUMBER_OF_NMSSM_INPUT_PARAMETERS);
      parameter[par] = value;
      has_been_set[par] = true;
   }

   /// get value of parameter
   double get(NMSSM_parameters par) const {
      assert(par < NUMBER_OF_NMSSM_INPUT_PARAMETERS);
      return parameter[par];
   }

   /// returns vector with supersymmetric NMSSM parameters
   DoubleVector get_nmpars() const {
      DoubleVector nmpars(5);
      nmpars(1) = get(NMSSM_input::lambda);
      nmpars(2) = get(NMSSM_input::kappa);
      if (!close(get(NMSSM_input::lambda),0.,EPSTOL)) {
         nmpars(3) = get(NMSSM_input::lambdaS) / get(NMSSM_input::lambda);
      } else {
         cout << "# Error: you set lambda * <S> = "
              << get(NMSSM_input::lambdaS) << ", but lambda is zero. "
            "Please set lambda (EXTPAR entry 61) to a non-zero value.\n";
      }
      nmpars(4) = get(NMSSM_input::xiF);
      nmpars(5) = get(NMSSM_input::muPrime);
      return nmpars;
   };

   /// returns the number of set NMSSM parameters
   unsigned get_number_of_set_parameters() const {
      unsigned num = 0;
      for (unsigned i = 0; i < NUMBER_OF_NMSSM_INPUT_PARAMETERS; i++)
         num += is_set(static_cast<NMSSM_parameters>(i));
      return num;
   }

   /// returns true if parameter was set, false otherwise
   bool is_set(NMSSM_parameters par) const { return has_been_set[par]; }

   /// returns true if input parameter set defines a Z3 symmetric NMSSM
   bool is_Z3_symmetric() const {
      return (!has_been_set[mu]      || close(parameter[mu]      , 0., EPSTOL))
         && (!has_been_set[BmuOverCosBetaSinBeta]
             || close(parameter[BmuOverCosBetaSinBeta], 0., EPSTOL))
         && (!has_been_set[muPrime]  || close(parameter[muPrime] , 0., EPSTOL))
         && (!has_been_set[mPrimeS2] || close(parameter[mPrimeS2], 0., EPSTOL))
         && (!has_been_set[xiF]      || close(parameter[xiF]     , 0., EPSTOL));
   }

   bool check_setup() {
       // check if the user has given enough input parameters
       if (get_number_of_set_parameters() == 6) {
          if (!check_6_parameter_input())
             return false;
       } else if (get_number_of_set_parameters() != 12) {
          cout << "# Error: " << get_number_of_set_parameters() << " NMSSM"
             " parameters given: " << (*this) << "  Please select either"
             " 6 or 12 parameters.\n";
          return false;
       }

       // check if the EWSB output parameter set is supported
       if (!ewsb_output_parameters_are_supported())
          return false;

       return true;
   }

   /// checks the 6 input parameters and sets the non-standard parameters to zero
   bool check_6_parameter_input() {
      if (get_number_of_set_parameters() != 6)
         return false;

      // check that the 6 set parameters are the ones from Eq. (60)
      // arxiv.org/abs/0801.0045
      const bool are_reduced_parameters = !is_set(mu) &&
         !is_set(BmuOverCosBetaSinBeta) && !is_set(muPrime) &&
         !is_set(mPrimeS2) && !is_set(xiF) && !is_set(xiS);

      if (!are_reduced_parameters) {
         cout << "# Error: 6 parameter set, but not all are \"standard"
            " parameters\".  Please select 6 parameters from {"
            "tanBeta, mHd2, mHu2, lambda, kappa, Alambda, Akappa, lambdaS,"
            " mS2}" << endl;
         return false;
      }

      set(mu, 0.);
      set(BmuOverCosBetaSinBeta, 0.);
      set(muPrime, 0.);
      set(mPrimeS2, 0.);
      set(xiF, 0.);
      set(xiS, 0.);

      return true;
   }

   /// checks if the unset parameters can be used as EWSB output
   bool ewsb_output_parameters_are_supported() {
      if (get_number_of_set_parameters() != 12)
         return false;

      bool supported = false;

      // check supported cases
      if (is_Z3_symmetric()) {
         if (!has_been_set[lambdaS] && !has_been_set[kappa] &&
             !has_been_set[mS2])
            supported = true;
      } else {
         if (!has_been_set[mu] && !has_been_set[BmuOverCosBetaSinBeta] &&
             !has_been_set[xiS])
            supported = true;
      }

      if (!supported) {
         cout << "# Error: the specified set of EWSB output parameters is"
            " is currently not supported: ";
         for (unsigned i = 0; i < NUMBER_OF_NMSSM_INPUT_PARAMETERS; i++) {
            if (!has_been_set[static_cast<NMSSM_parameters>(i)])
               cout << parameter_names[i] << ", ";
         }
         cout << '\n';
      }

      return supported;
   }

   friend std::ostream& operator<<(std::ostream& lhs, const NMSSM_input& rhs) {
      for (unsigned i = 0; i < NUMBER_OF_NMSSM_INPUT_PARAMETERS; i++) {
         if (rhs.is_set(static_cast<NMSSM_parameters>(i))) {
            lhs << rhs.parameter_names[i] << " = "
                << rhs.get(static_cast<NMSSM_parameters>(i)) << ", ";
         }
      }
      return lhs;
   }

private:
   double parameter[NUMBER_OF_NMSSM_INPUT_PARAMETERS];  ///< NMSSM parameters
   bool has_been_set[NUMBER_OF_NMSSM_INPUT_PARAMETERS];
   static char const * const parameter_names[NUMBER_OF_NMSSM_INPUT_PARAMETERS];
};
