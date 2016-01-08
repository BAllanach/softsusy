/* -*- mode: c++; c-basic-offset: 2; -*- */
/** \file nmssmjacobian.h

    \brief Provides routines for calculating Jacobian fine-tuning.
 */

#ifndef NMSSMJACOBIAN_H
#define NMSSMJACOBIAN_H

#include "linalg.h"

namespace softsusy {

  class NmssmSoftsusy;

  /// \brief Class for calculating Jacobian fine-tuning measure.
  class NmssmJacobian {
  public:

    explicit NmssmJacobian(bool doTop = false);
    virtual ~NmssmJacobian();

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

    /// \brief Resets the error status to false.
    void resetError() { hasError = false; }

    /// \brief Displays whether an error occurred during the calculation.
    /// \return true if an error occurred
    bool displayError() const { return hasError; }

    /// \brief Displays the Jacobian for transforming from high- to low-scale
    /// parameters.
    /// \return the Jacobian matrix of derivatives for the transformation
    DoubleMatrix displayRGFlowJacobian() const { return jacRGFlow; }

    /// \brief Displays the Jacobian for transforming low-scale parameters
    /// to observables.
    /// \return the Jacobian matrix of derivatives for the transformation
    DoubleMatrix displayEWSBJacobian() const { return jacEWSB; }

    /// \brief Displays the Jacobian for transforming from low- to high-scale
    /// parameters.
    /// \return the Jacobian matrix of derivatives for the transformation
    DoubleMatrix displayInverseRGFlowJacobian() const { return invJacRGFlow; }

    /// \brief Displays the Jacobian for transforming observables to
    /// low-scale parameters.
    /// \return the Jacobian matrix of derivatives for the transformation
    DoubleMatrix displayInverseEWSBJacobian() const { return invJacEWSB; }

    /// \brief Calculates the Jacobian transforming parameters to observables.
    /// \param[in] model the model to calculate the Jacobian for
    /// \return the value of the Jacobian
    double calcFTJacobian(NmssmSoftsusy& model);

    /// \brief Calculates the Jacobian transforming parameters to observables.
    /// \param[in] model the model to calculate the Jacobian for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the Jacobian
    double calcFTJacobian(NmssmSoftsusy& model, double mx);

    /// \brief Calculates the Jacobian transforming observables to parameters.
    /// \param[in] model the model to calculate the Jacobian for
    /// \return the value of the Jacobian
    double calcFTInverseJacobian(NmssmSoftsusy& model);

    /// \brief Calculates the Jacobian transforming observables to parameters.
    /// \param[in] model the model to calculate the Jacobian for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the Jacobian
    double calcFTInverseJacobian(NmssmSoftsusy& model, double mx);

    /// \brief Calculates the fine-tuning using the Jacobian measure.
    /// \param[in] model the model to calculate the fine-tuning for
    /// \return the value of the fine-tuning
    double calcDeltaJ(NmssmSoftsusy& model);

    /// \brief Calculates the fine-tuning using the Jacobian measure.
    /// \param[in] model the model to calculate the fine-tuning for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the fine-tuning
    double calcDeltaJ(NmssmSoftsusy& model, double mx);

    static double calcMz(NmssmSoftsusy* m, bool getRunningMass = false);
    static double calcMt(NmssmSoftsusy* m, bool getRunningMass = false);

  private:
    enum Parameters { Mzsq, Tanb, Svev, Mtsq, Lambda, Kappa,
                      SMu, M3Sq, XiS, Mh1Sq, Mh2Sq, MsSq, Yt };

    DoubleMatrix jacRGFlow;
    DoubleMatrix jacEWSB;
    DoubleMatrix invJacRGFlow;
    DoubleMatrix invJacEWSB;
    bool includeTop;
    bool useRunningMasses;
    bool hasError;

    struct EWSBPars {
      NmssmSoftsusy* model;
      Parameters independent;
      Parameters dependent;
      DoubleVector outputs;
      bool useRunningMasses;

      EWSBPars()
        : model(0), independent(Mzsq), dependent(Lambda),
          outputs(3), useRunningMasses(false)
        {}
    };

    struct RGFlowPars {
      NmssmSoftsusy* model;
      Parameters independent;
      Parameters dependent;
      double toScale;
    };

    static double calcRunningParameter(double x, void* parameters);
    double calcRGDerivative(NmssmSoftsusy& model, Parameters dep,
                            Parameters indep, double toScale);
    double calcRGFlowJacobian(NmssmSoftsusy& model, double startScale,
                              double endScale);
    static int ewsbOutputErrors(const DoubleVector & guess, void* parameters,
                                DoubleVector & errors);
    static void fixEWSBOutputs(EWSBPars* pars, int & err);
    static double calcEWSBParameter(double x, void* parameters);
    static double calcEWSBOutput(double x, void* parameters);
    double calcEWSBOutputDerivative(NmssmSoftsusy& model, Parameters dep,
                                    Parameters indep);
    double calcEWSBParameterDerivative(NmssmSoftsusy& model, Parameters dep,
                                       Parameters indep);
    double calcEWSBJacobian(NmssmSoftsusy& model);
    double calcInverseEWSBJacobian(NmssmSoftsusy& model);
  };

} /// namespace softsusy

#endif
