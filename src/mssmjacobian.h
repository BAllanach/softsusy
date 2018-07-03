/** \file mssmjacobian.h

    \brief Provides routines for calculating Jacobian fine-tuning.
*/

#ifndef MSSMJACOBIAN_H
#define MSSMJACOBIAN_H

#include "linalg.h"

#include <utility>

namespace softsusy {

  class MssmSoftsusy;

  class MssmJacobian {
  public:

    explicit MssmJacobian(bool doTop = false);
    ~MssmJacobian();

    /// \brief Displays whether the top Yukawa is included in the tuning.
    /// \return true if the top Yukawa is included
    bool displayIncludeTopFlag() const { return includeTop; }

    /// \brief Sets whether to include top Yukawa derivatives
    /// \param[in] flag if true, include top Yukawa
    void setIncludeTopFlag(bool flag) { includeTop = flag; }

    /// \brief Displays whether running or pole masses are being used.
    /// \return true if running masses are used instead of pole masses
    bool displayUseRunningMassesFlag() const { return useRunningMasses; }

    /// \brief Sets whether to use running or pole masses in derivatives.
    /// \param[in] flag if true, use running masses instead of pole masses
    void setUseRunningMassesFlag(bool flag) { useRunningMasses = flag; }

    /// \brief Displays whether SUGRA style trilinears are assumed
    /// \return true if SUGRA style trilinears are being assumed
    bool displayUseSugraTrilinearsFlag() const { return useSugraTrilinears; }

    /// \brief Sets whether to assume SUGRA style trilinears
    /// \param[in] flag if true, assume SUGRA style trilinears
    void setUseSugraTrilinearsFlag(bool flag) { useSugraTrilinears = flag; }

    /// \brief Resets the error status to false.
    void resetError() { hasError = false; }

    /// \brief Displays whether an error occurred during the calculation.
    /// \return true if an error occurred
    bool displayError() const { return hasError; }

    /// \brief Displays the Jacobian for transforming from low- to high-scale
    /// parameters.
    /// \return the Jacobian matrix of derivatives for the transformation
    DoubleMatrix displayInverseRGFlowJacobian() const { return invJacRGFlow; }

    /// \brief Displays the error in the Jacobian for transforming from low-
    /// to high-scale parameters.
    /// \return the absolute error on the Jacobian matrix
    DoubleMatrix displayInverseRGFlowJacobianErrors() const {
      return invJacRGFlowErrors;
    }

    /// \brief Displays the Jacobian for transforming observables to
    /// low-scale parameters.
    /// \return the Jacobian matrix of derivatives for the transformation
    DoubleMatrix displayInverseEWSBJacobian() const { return invJacEWSB; }

    /// \brief Displays the error in the Jacobian for transforming observables
    /// to low-scale parameters.
    /// \return the absolute error on the Jacobian matrix
    DoubleMatrix displayInverseEWSBJacobianErrors() const {
      return invJacEWSBErrors;
    }

    /// \brief Calculates the Jacobian transforming observables to parameters.
    /// \param[in] m the model to calculate the Jacobian for
    /// \return the value of the Jacobian
    double calcFTInverseJacobian(const MssmSoftsusy& m);

    /// \brief Calculates the Jacobian transforming observables to parameters.
    /// \param[in] m the model to calculate the Jacobian for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the Jacobian
    double calcFTInverseJacobian(const MssmSoftsusy& m, double mx);

    /// \brief Calculates the fine-tuning using the Jacobian measure.
    /// \param[in] m the model to calculate the fine-tuning for
    /// \return the value of the fine-tuning. It returns NaN if it had trouble
    /// calculating the fine-tuning
    double calcDeltaJ(const MssmSoftsusy& m);

    /// \brief Calculates the fine-tuning using the Jacobian measure.
    /// \param[in] m the model to calculate the fine-tuning for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the fine-tuning
    double calcDeltaJ(const MssmSoftsusy& m, double mx);

  private:
    enum Parameters { Mzsq, Tanb, Mtsq, SMu, M3Sq, Mh1Sq, Mh2Sq, Yt };

    DoubleMatrix invJacRGFlow;
    DoubleMatrix invJacEWSB;

    DoubleMatrix invJacRGFlowErrors;
    DoubleMatrix invJacEWSBErrors;

    bool includeTop;
    bool useRunningMasses;
    bool useSugraTrilinears;
    bool hasError;

    struct EWSBPars {
      MssmSoftsusy* model;
      Parameters independent;
      Parameters dependent;
      DoubleVector outputs;
      bool useRunningMasses;
      bool useSugraTrilinears;

      EWSBPars()
        : model(0), independent(Mzsq), dependent(SMu),
          outputs(2), useRunningMasses(false),
          useSugraTrilinears(false)
        {}
    };

    struct RGFlowPars {
      MssmSoftsusy* model;
      Parameters independent;
      Parameters dependent;
      double toScale;
      bool useSugraTrilinears;
    };

    static double calcMz(MssmSoftsusy& model, bool getRunningMass = false);
    static double calcMt(MssmSoftsusy& model, bool getRunningMass = false);

    static double calcRunningParameter(double x, void* parameters);
    std::pair<double,double> calcRGDerivative(MssmSoftsusy& model,
                                              Parameters dep, Parameters indep,
                                              double toScale);
    double calcRGFlowJacobian(MssmSoftsusy& model, double startScale,
                              double endScale);
    static double calcEWSBOutput(double x, void* parameters);
    std::pair<double,double> calcEWSBOutputDerivative(MssmSoftsusy& model,
                                                      Parameters dep,
                                                      Parameters indep);
    double calcInverseEWSBJacobian(MssmSoftsusy& model);
  };

} /// namespace softsusy

#endif
