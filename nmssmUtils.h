
#ifndef NMSSM_UTILS_H
#define NMSSM_UTILS_H

#include <iosfwd>
#include <string>

#include "linalg.h"

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

   NMSSM_input();
   ~NMSSM_input() {}

   /// set parameter to given value
   void set(NMSSM_parameters, double);
   /// get value of parameter
   double get(NMSSM_parameters) const;
   /// returns vector with supersymmetric NMSSM parameters
   DoubleVector get_nmpars() const;
   /// returns the number of set NMSSM parameters
   unsigned get_number_of_set_parameters() const;
   /// returns true if parameter was set, false otherwise
   bool is_set(NMSSM_parameters par) const;
   /// returns true if input parameter set defines a Z3 symmetric NMSSM
   bool is_Z3_symmetric() const;
   /// checks the NMSSM parameter setup, throws if not SLHA2 confrom
   void check_setup();

   friend std::ostream& operator<<(std::ostream&, const NMSSM_input&);

private:
   double parameter[NUMBER_OF_NMSSM_INPUT_PARAMETERS];  ///< NMSSM parameters
   bool has_been_set[NUMBER_OF_NMSSM_INPUT_PARAMETERS];
   static char const * const parameter_names[NUMBER_OF_NMSSM_INPUT_PARAMETERS];

   /// checks the 6 input parameters and sets the non-standard parameters to zero
   void check_6_input_parameters();
   /// checks if the unset parameters can be used as EWSB output
   void check_ewsb_output_parameters() const;
};

#endif
