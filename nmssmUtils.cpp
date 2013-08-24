
#include "nmssmUtils.h"
#include "def.h"
#include "utils.h"

#include <cassert>
#include <sstream>

char const * const NMSSM_input::parameter_names[NUMBER_OF_NMSSM_INPUT_PARAMETERS] = {
   "tanBeta", "mHd2", "mHu2", "mu", "BmuOverCosBetaSinBeta", "lambda",
   "kappa", "Alambda", "Akappa", "lambda*S", "xiF", "xiS", "muPrime",
   "mPrimeS2", "mS2"
};

NMSSM_input::NMSSM_input()
  : parameter()    // sets all parameters to zero
  , has_been_set() // sets all values to zero (false)
{}

void NMSSM_input::set(NMSSM_parameters par, double value) {
   assert(par < NUMBER_OF_NMSSM_INPUT_PARAMETERS);
   parameter[par] = value;
   has_been_set[par] = true;
}

double NMSSM_input::get(NMSSM_parameters par) const {
   assert(par < NUMBER_OF_NMSSM_INPUT_PARAMETERS);
   return parameter[par];
}

DoubleVector NMSSM_input::get_nmpars() const {
   DoubleVector nmpars(5);
   nmpars(1) = get(NMSSM_input::lambda);
   nmpars(2) = get(NMSSM_input::kappa);
   if (is_set(NMSSM_input::lambdaS)) {
      if (!close(get(NMSSM_input::lambda),0.,EPSTOL)) {
         nmpars(3) = get(NMSSM_input::lambdaS) / get(NMSSM_input::lambda);
      } else {
         std::string msg =
            "# Error: you set lambda * <S> to a non-zero value"
            ", but lambda is zero.  "
            "Please set lambda (EXTPAR entry 61) to a non-zero value.";
         throw msg;
      }
   }
   nmpars(4) = get(NMSSM_input::xiF);
   nmpars(5) = get(NMSSM_input::muPrime);
   return nmpars;
};

unsigned NMSSM_input::get_number_of_set_parameters() const {
   unsigned num = 0;
   for (unsigned i = 0; i < NUMBER_OF_NMSSM_INPUT_PARAMETERS; i++)
      num += is_set(static_cast<NMSSM_parameters>(i));
   return num;
}

bool NMSSM_input::is_set(NMSSM_parameters par) const {
   assert(par < NUMBER_OF_NMSSM_INPUT_PARAMETERS);
   return has_been_set[par];
}

bool NMSSM_input::is_Z3_symmetric() const {
   return (!is_set(mu)      || close(parameter[mu]      , 0., EPSTOL))
      && (!is_set(BmuOverCosBetaSinBeta)
          || close(parameter[BmuOverCosBetaSinBeta], 0., EPSTOL))
      && (!is_set(muPrime)  || close(parameter[muPrime] , 0., EPSTOL))
      && (!is_set(mPrimeS2) || close(parameter[mPrimeS2], 0., EPSTOL))
      && (!is_set(xiF)      || close(parameter[xiF]     , 0., EPSTOL));
}

void NMSSM_input::check_setup() {
   const unsigned number_of_set_parameter = get_number_of_set_parameters();

   if (number_of_set_parameter == 6) {
      check_6_input_parameters();
   } else if (number_of_set_parameter != 12) {
      std::ostringstream msg;
      msg << "# Error: " << get_number_of_set_parameters() << " NMSSM"
         " parameters given: " << (*this) << "  Please select either"
         " 6 or 12 parameters.";
      throw msg.str();
   }

   check_ewsb_output_parameters();
}

void NMSSM_input::check_6_input_parameters() {
   if (get_number_of_set_parameters() != 6) {
      std::string msg =
         "# Internal error: check_6_input_parameters() called but"
         " not 6 paramters set";
      throw msg;
   }

   // check that the 6 set parameters are the ones from Eq. (60)
   // arxiv.org/abs/0801.0045
   const bool are_reduced_parameters = !is_set(mu) &&
      !is_set(BmuOverCosBetaSinBeta) && !is_set(muPrime) &&
      !is_set(mPrimeS2) && !is_set(xiF) && !is_set(xiS);

   if (!are_reduced_parameters) {
      std::string msg =
         "# Error: 6 parameter set, but not all are \"standard"
         " parameters\".  Please select 6 parameters from {"
         "tanBeta, mHd2, mHu2, lambda, kappa, Alambda, Akappa, lambdaS,"
         " mS2}";
      throw msg;
   }

   set(mu, 0.);
   set(BmuOverCosBetaSinBeta, 0.);
   set(muPrime, 0.);
   set(mPrimeS2, 0.);
   set(xiF, 0.);
   set(xiS, 0.);
}

void NMSSM_input::check_ewsb_output_parameters() const {
   if (get_number_of_set_parameters() != 12) {
      std::string msg =
         "# Internal error: check_ewsb_output_parameters()"
         " called but not 12 paramters set";
      throw msg;
   }

   bool supported = false;

   // check supported cases
   const bool Z3_symmetric = is_Z3_symmetric();
   if (Z3_symmetric) {
      if (!is_set(lambdaS) && !is_set(kappa) && !is_set(mS2))
         supported = true;
   } else {
      if (!is_set(mu) && !is_set(BmuOverCosBetaSinBeta) && !is_set(xiS))
         supported = true;
   }

   if (!supported) {
      std::ostringstream msg;
      msg << "# Error: the specified set of EWSB output parameters is"
         " currently not supported for a Z3 "
          << (Z3_symmetric ? "symmetric" : "violating") << " NMSSM: ";
      for (unsigned i = 0; i < NUMBER_OF_NMSSM_INPUT_PARAMETERS; i++) {
         if (!is_set(static_cast<NMSSM_parameters>(i)))
            msg << parameter_names[i] << ", ";
      }
      msg << '\n';
      throw msg.str();
   }
}

std::ostream& operator<<(std::ostream& lhs, const NMSSM_input& rhs) {
   for (unsigned i = 0; i < NMSSM_input::NUMBER_OF_NMSSM_INPUT_PARAMETERS; i++) {
      if (rhs.is_set(static_cast<NMSSM_input::NMSSM_parameters>(i))) {
         lhs << rhs.parameter_names[i] << " = "
             << rhs.get(static_cast<NMSSM_input::NMSSM_parameters>(i)) << ", ";
      }
   }
   return lhs;
}
