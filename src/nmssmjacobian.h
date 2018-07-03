/** \file nmssmjacobian.h

    \brief Provides routines for calculating Jacobian fine-tuning.
 */

#ifndef NMSSMJACOBIAN_H
#define NMSSMJACOBIAN_H

#include "linalg.h"

#include <utility>

namespace softsusy {

  class NmssmSoftsusy;

  /// \brief Class for calculating Jacobian fine-tuning measure.
  ///
  /// This class provides methods for calculating the Jacobian
  /// fine-tuning measure as defined in arXiv:1312.4150,
  /// \f$\Delta_J = | \partial O_j / \partial p_i |\f$, where
  /// \f$\{O_j\}\f$ refers to a set of observables and
  /// \f$\{p_i\}\f$ to a set of model parameters, and \f$ \partial
  /// O_j / \partial p_i\f$ is the Jacobian for the transformation
  /// from the observables to the parameters.
  ///
  /// The sets \f$\{O_j\}\f$ and \f$\{p_i\}\f$ used depend on
  /// the global flags softsusy::Z3 and softsusy::SoftHiggsOut.
  /// If softsusy::SoftHiggsOut is true, the parameters used are
  /// the soft Higgs masses at an input scale \c mx, \f$\{ m_{H_1,0}^2,
  /// m_{H_2,0}^2, m_{S_0}^2\} \f$, and the observables are
  /// \f$\{M_Z^2, \tan\beta, s\}\f$, irrespective of the value of
  /// softsusy::Z3.  If softsusy::SoftHiggsOut is false and
  /// softsusy::Z3 is true, the parameters are taken to be
  /// \f$\{\lambda_0, \kappa_0, m_{S_0}^2\}\f$ at \c mx.  The
  /// observables in this case are \f$\{M_Z^2, \tan\beta,
  /// \lambda\}\f$.  If both flags are false, the parameters
  /// used are \f$\{\mu_0, m_{3_0}^2, \xi_{S_0}\}\f$ and the set of
  /// observables is \f$\{M_Z^2, \tan\beta, s\}\f$.
  ///
  /// The top Yukawa \f$y_t\f$ can be included in the set of parameters,
  /// and \f$M_t^2\f$ in the set of observables, by calling
  /// setIncludeTopFlag() with the desired flag.  Whether \f$M_Z^2\f$
  /// and \f$M_t^2\f$ are taken to be the pole or running masses may be
  /// set by calling setUseRunningMassesFlag() with the appropriate flag.
  /// By default the pole masses are used.
  ///
  /// The treatment of the soft SUSY breaking trilinears can be modified
  /// using the setUseSugraTrilinearsFlag routine.  If set to true, the
  /// independent parameters are taken to be the capital \f$A\f$ trilinears,
  /// e.g. \f$A_\lambda = a_\lambda / \lambda\f$.  In this case, it is
  /// the quantity \f$A\f$ that is held fixed in partial derivatives, rather
  /// than the quantity \f$a\f$.  If this flag is false, the trilinears
  /// \f$a\f$ are treated as independent of the corresponding coupling, and
  /// are held fixed in the partial derivatives.  By default this flag is
  /// false, i.e. the \f$a\f$ couplings are taken to be the independent
  /// parameters.
  class NmssmJacobian {
  public:

    explicit NmssmJacobian(bool doTop = false);
    ~NmssmJacobian();

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

    /// \brief Resets the warning status to false
    void resetWarnings() {
      failedVevIterationWarning = false;
      inaccurateHiggsMassWarning = false;
    }

    /// \brief Displays whether a warning was flagged during the calculation.
    /// \return true if a warning was flagged
    bool displayWarning() const {
      return failedVevIterationWarning || inaccurateHiggsMassWarning; }

    /// \brief Displays whether the VEV iteration failed
    /// \return true if the VEV iteration failed
    bool displayFailedVevIterationWarning() const {
      return failedVevIterationWarning;
    }

    /// \brief Displays whether an inaccurate Higgs mass was encountered
    /// \return true if an inaccurate Higgs mass was flagged
    bool displayInaccurateHiggsMassWarning() const {
      return inaccurateHiggsMassWarning;
    }

    /// \brief Returns the worst achieved accuracy for the Higgs masses
    /// \return the worst case accuracy for the Higgs masses
    double displayWorstCaseHiggsAccuracy() const {
      return worstCaseHiggsAccuracy;
    }

    /// \brief Displays the Jacobian for transforming from high- to low-scale
    /// parameters.
    /// \return the Jacobian matrix of derivatives for the transformation
    DoubleMatrix displayRGFlowJacobian() const { return jacRGFlow; }

    /// \brief Displays the error in the Jacobian for transforming from high-
    /// to low-scale parameters.
    /// \return the absolute error on the Jacobian matrix
    DoubleMatrix displayRGFlowJacobianErrors() const { return jacRGFlowErrors; }

    /// \brief Displays the Jacobian for transforming low-scale parameters
    /// to observables.
    /// \return the Jacobian matrix of derivatives for the transformation
    DoubleMatrix displayEWSBJacobian() const { return jacEWSB; }

    /// \brief Displays the error in the Jacobian for transforming low-scale
    /// parameters to observables.
    /// \return the absolute error on the Jacobian matrix
    DoubleMatrix displayEWSBJacobianErrors() const { return jacEWSBErrors; }

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

    /// \brief Calculates the Jacobian transforming parameters to observables.
    /// \param[in] m the model to calculate the Jacobian for
    /// \return the value of the Jacobian
    double calcFTJacobian(const NmssmSoftsusy& m);

    /// \brief Calculates the Jacobian transforming parameters to observables.
    /// \param[in] m the model to calculate the Jacobian for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the Jacobian
    double calcFTJacobian(const NmssmSoftsusy& m, double mx);

    /// \brief Calculates the Jacobian transforming observables to parameters.
    /// \param[in] m the model to calculate the Jacobian for
    /// \return the value of the Jacobian
    double calcFTInverseJacobian(const NmssmSoftsusy& m);

    /// \brief Calculates the Jacobian transforming observables to parameters.
    /// \param[in] m the model to calculate the Jacobian for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the Jacobian
    double calcFTInverseJacobian(const NmssmSoftsusy& m, double mx);

    /// \brief Calculates the Jacobian transforming observables to parameters.
    /// \param[in] m the model to calculate the Jacobian for
    /// \return the value of the Jacobian
    double calcFTInverseJacobianLogs(const NmssmSoftsusy& m);

    /// \brief Calculates the Jacobian transforming observables to parameters.
    /// \param[in] m the model to calculate the Jacobian for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the Jacobian
    double calcFTInverseJacobianLogs(const NmssmSoftsusy& m, double mx);

    /// \brief Calculates the fine-tuning using the Jacobian measure.
    /// \param[in] m the model to calculate the fine-tuning for
    /// \return the value of the fine-tuning. It returns NaN if it had trouble
    /// calculating the fine-tuning
    double calcDeltaJ(const NmssmSoftsusy& m);

    /// \brief Calculates the fine-tuning using the Jacobian measure.
    /// \param[in] m the model to calculate the fine-tuning for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the fine-tuning
    double calcDeltaJ(const NmssmSoftsusy& m, double mx);

    /// \brief Calculates the fine-tuning using logarithmic derivatives.
    /// \param[in] m the model to calculate the fine-tuning for
    /// \return the value of the fine-tuning
    double calcDeltaJLogs(const NmssmSoftsusy& m);

    /// \brief Calculates the fine-tuning using logarithmic derivatives.
    /// \param[in] m the model to calculate the fine-tuning for
    /// \param[in] mx the scale at which the input parameters are defined
    /// \return the value of the fine-tuning
    double calcDeltaJLogs(const NmssmSoftsusy& m, double mx);

  private:
    enum Parameters { Mzsq, Tanb, Svev, Mtsq, Lambda, Kappa,
                      SMu, M3Sq, XiS, Mh1Sq, Mh2Sq, MsSq, Yt };

    DoubleMatrix jacRGFlow;
    DoubleMatrix jacEWSB;
    DoubleMatrix invJacRGFlow;
    DoubleMatrix invJacEWSB;

    DoubleMatrix jacRGFlowErrors;
    DoubleMatrix jacEWSBErrors;
    DoubleMatrix invJacRGFlowErrors;
    DoubleMatrix invJacEWSBErrors;

    bool includeTop;
    bool useRunningMasses;
    bool useSugraTrilinears;
    bool hasError;
    bool failedVevIterationWarning;
    bool inaccurateHiggsMassWarning;

    double worstCaseHiggsAccuracy;
    
    struct EWSBPars {
      NmssmSoftsusy* model;
      Parameters independent;
      Parameters dependent;
      DoubleVector outputs;
      bool useRunningMasses;
      bool useSugraTrilinears;
      bool failedVevIteration;
      bool inaccurateHiggsMass;
      double higgsMassAccuracy;

      EWSBPars()
        : model(0), independent(Mzsq), dependent(Lambda),
          outputs(3), useRunningMasses(false),
          useSugraTrilinears(false), failedVevIteration(false),
          inaccurateHiggsMass(false), higgsMassAccuracy(0.)
        {}
    };

    struct RGFlowPars {
      NmssmSoftsusy* model;
      Parameters independent;
      Parameters dependent;
      double toScale;
      bool useSugraTrilinears;
    };

    static double calcMz(NmssmSoftsusy& model, bool getRunningMass = false);
    static double calcMt(NmssmSoftsusy& model, bool getRunningMass = false);

    static double calcRunningParameter(double x, void* parameters);
    static double calcLogRunningParameter(double x, void* parameters);
    std::pair<double,double> calcRGDerivative(NmssmSoftsusy& model,
                                              Parameters dep, Parameters indep,
                                              double toScale) const;
    std::pair<double,double> calcRGLogDerivative(NmssmSoftsusy& model,
                                                 Parameters dep, Parameters indep,
                                                 double toScale) const;
    double calcRGFlowJacobian(NmssmSoftsusy& model, double startScale,
                              double endScale);
    double calcRGFlowJacobianLogs(NmssmSoftsusy& model, double startScale,
                                  double endScale);
    static int ewsbOutputErrors(const DoubleVector & guess, void* parameters,
                                DoubleVector & errors);
    static void fixEWSBOutputs(EWSBPars* pars, int & err);
    static double calcEWSBParameter(double x, void* parameters);
    static double calcEWSBOutput(double x, void* parameters);
    static double calcLogEWSBOutput(double x, void* parameters);
    std::pair<double,double> calcEWSBOutputDerivative(NmssmSoftsusy& model,
                                                      Parameters dep,
                                                      Parameters indep);
    std::pair<double,double> calcEWSBOutputLogDerivative(NmssmSoftsusy& model,
                                                         Parameters dep,
                                                         Parameters indep);
    std::pair<double,double> calcEWSBParameterDerivative(NmssmSoftsusy& model,
                                                         Parameters dep,
                                                         Parameters indep);
    double calcEWSBJacobian(NmssmSoftsusy& model);
    double calcInverseEWSBJacobian(NmssmSoftsusy& model);
    double calcInverseEWSBJacobianLogs(NmssmSoftsusy& model);

    bool checkDerivativeErrors(DoubleMatrix& derivs, const DoubleMatrix& errors,
                               const std::vector<double>& paramValues) const;
  };

} /// namespace softsusy

#endif
