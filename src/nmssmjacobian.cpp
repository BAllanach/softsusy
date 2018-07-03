/** \file nmssmjacobian.cpp

    \brief Implementation of routines for calculating Jacobian fine-tuning.
 */

#include "config.h"
#include "def.h"
#include "nmssmjacobian.h"
#include "nmssmsoftsusy.h"
#include "numerics.h"
#include "utils.h"

#include <iostream>
#include <limits>
#include <vector>

#ifdef ENABLE_GSL
#include <gsl/gsl_deriv.h>
#endif

namespace softsusy {

  NmssmJacobian::NmssmJacobian(bool doTop)
     : jacRGFlow(3,3), jacEWSB(3,3)
     , invJacRGFlow(3,3), invJacEWSB(3,3)
     , jacRGFlowErrors(3,3), jacEWSBErrors(3,3)
     , invJacRGFlowErrors(3,3), invJacEWSBErrors(3,3)
     , includeTop(doTop), useRunningMasses(false)
     , useSugraTrilinears(false), hasError(false)
     , failedVevIterationWarning(false)
     , inaccurateHiggsMassWarning(false)
     , worstCaseHiggsAccuracy(0.) {
  }

  NmssmJacobian::~NmssmJacobian() {}

  /// Calculates the fine-tuning for the model, using the Jacobian
  /// measure.  The parameters are taken to be defined at the input
  /// scale \c mx obtained from a call to \c m.displayMxBC(),
  /// while the observables are calculated at the SUSY scale
  /// returned by \c m.calcMs().  The Jacobian matrix is
  /// calculated by calling calcFTInverseJacobian() internally, and
  /// can be accessed after the calculation.
  /// \see calcFTInverseJacobian(const NmssmSoftsusy&)
  double NmssmJacobian::calcDeltaJ(const NmssmSoftsusy& m) {
    return calcDeltaJ(m, m.displayMxBC());
  }

  /// Calculates the fine-tuning as in calcDeltaJ(const NmssmSoftsusy& m),
  /// but with \c mx set to the given value instead of that returned by
  /// \c m.displayMxBC().
  double NmssmJacobian::calcDeltaJ(const NmssmSoftsusy& m, double mx) {
    NmssmSoftsusy model(m);

    model.runto(model.calcMs());

    model.calcDrBarPars();
    const double mtrun = model.displayDrBarPars().mt;
    const double sinthDRbar = model.calcSinthdrbar();
    model.doTadpoles(mtrun, sinthDRbar);

    const double determinant = calcFTInverseJacobian(model, mx);

    const double mz2 = sqr(calcMz(model, useRunningMasses));
    const double tb = model.displayTanb();

    double denominator = mz2 * tb;
    if ( m.displayZ3() ) {
      denominator *= model.displayLambda();
    } else {
      denominator *= model.displaySvev();
    }
    if (includeTop) {
      denominator *= sqr(calcMt(model, useRunningMasses));
    }

    model.runto(mx);

    double numerator;
    if (SoftHiggsOut) {
      numerator = model.displayMh1Squared() * model.displayMh2Squared()
        * model.displayMsSquared();
    } else if ( m.displayZ3() ) {
      numerator = model.displayLambda() * model.displayKappa()
        * model.displayMsSquared();
    } else {
      numerator = model.displaySusyMu() * model.displayM3Squared()
        * model.displayXiS();
    }
    if (includeTop) numerator *= model.displayYukawaElement(YU, 3, 3);

    double ans = fabs(numerator * determinant / denominator);
    if (displayError()) ans = asin(2.0);
    
    return ans;
  }

  /// Calculates the fine-tuning for the model, using the Jacobian
  /// measure.  The Jacobian matrix is calculated using derivatives
  /// computed with respect to the logarithms of the parameters and
  /// observables.  The parameters are taken to be defined at the
  /// input scale \c mx obtained from a call to \c m.displayMxBC(),
  /// while the observables are calculated at the SUSY scale
  /// returned by \c m.calcMs().   Internally, the Jacobian matrix
  /// is calculated by calling calcFTInverseJacobianLogs(), and
  /// can be accessed after the calculation.
  /// \see calcFTInverseJacobianLogs(const NmssmSoftsusy&)
  double NmssmJacobian::calcDeltaJLogs(const NmssmSoftsusy& m) {
    return calcDeltaJLogs(m, m.displayMxBC());
  }

  /// Calculates the fine-tuning as in calcDeltaJLogs(const NmssmSoftsusy& m),
  /// but with \c mx set to the given value instead of that returned by
  /// \c m.displayMxBC().
  double NmssmJacobian::calcDeltaJLogs(const NmssmSoftsusy& m, double mx) {
    NmssmSoftsusy model(m);

    model.runto(model.calcMs());

    model.calcDrBarPars();
    const double mtrun = model.displayDrBarPars().mt;
    const double sinthDRbar = model.calcSinthdrbar();
    model.doTadpoles(mtrun, sinthDRbar);

    // derivatives due to RG flow are not calculated logarithmically
    double prefactor;
    if (SoftHiggsOut) {
      prefactor = 1. / (model.displayMh1Squared() * model.displayMh2Squared()
                         * model.displayMsSquared());
      if (includeTop) prefactor /= model.displayYukawaElement(YU, 3, 3);
      model.runto(mx);
      prefactor *= model.displayMh1Squared() * model.displayMh2Squared()
        * model.displayMsSquared();
      if (includeTop) prefactor *= model.displayYukawaElement(YU, 3, 3);
    } else if ( m.displayZ3() ) {
      prefactor = 1. / (model.displayLambda() * model.displayKappa()
                        * model.displayMsSquared());
      if (includeTop) prefactor /= model.displayYukawaElement(YU, 3, 3);
      model.runto(mx);
      prefactor *= model.displayLambda() * model.displayKappa()
        * model.displayMsSquared();
      if (includeTop) prefactor *= model.displayYukawaElement(YU, 3, 3);
    } else {
      prefactor = 1. / (model.displaySusyMu() * model.displayM3Squared()
                        * model.displayXiS());
      if (includeTop) prefactor /= model.displayYukawaElement(YU, 3, 3);
      model.runto(mx);
      prefactor *= model.displaySusyMu() * model.displayM3Squared()
         * model.displayXiS();
      if (includeTop) prefactor *= model.displayYukawaElement(YU, 3, 3);
    }

    const double determinant = calcFTInverseJacobianLogs(model, mx);

    return fabs(prefactor * determinant);
  }

  /// Calculates the Jacobian of the form \f$ J^{-1} = |
  /// \partial O / \partial p | \f$.  The transformation is done
  /// in two stages.  In the first, the observables at the
  /// SUSY scale, obtained from \c m.calcMs(), are traded for
  /// the parameters at this scale using the EWSB conditions.
  /// The Jacobian matrix for this transformation can be
  /// accessed afterwards using displayInverseEWSBJacobian().
  /// Then, the parameters at this scale are transformed to
  /// parameters at the scale \c mx, obtained from calling
  /// \c m.displayMxBC(), using the RGEs.  The Jacobian matrix for
  /// this transformation may be accessed by calling
  /// displayInverseRGFlowJacobian().
  double NmssmJacobian::calcFTInverseJacobian(const NmssmSoftsusy& m) {
    return calcFTInverseJacobian(m, m.displayMxBC());
  }

  /// Calculates the Jacobian as in
  /// calcFTInverseJacobian(const NmssmSoftsusy& m), but with \c mx set to
  /// the given value instead of that returned by \c m.displayMxBC().
  double NmssmJacobian::calcFTInverseJacobian(const NmssmSoftsusy& m, double mx) {
    NmssmSoftsusy model(m);

    model.runto(model.calcMs());

    const double scale = model.displayMu();
    const double rgDet = calcRGFlowJacobian(model, mx, scale);
    const double ewsbDet = calcInverseEWSBJacobian(model);

    return ewsbDet * rgDet;
  }

  /// Calculates the Jacobian of the form \f$ J^{-1} = |
  /// \partial \ln O / \partial \ln p | \f$.  The transformation is
  /// done in two stages.  In the first, the observables at the
  /// SUSY scale, obtained from \c m.calcMs(), are traded for
  /// the parameters at this scale using the EWSB conditions.
  /// The Jacobian matrix for this transformation can be
  /// accessed afterwards using displayInverseEWSBJacobian().
  /// Then, the parameters at this scale are transformed to
  /// parameters at the scale \c mx, obtained from calling
  /// \c m.displayMxBC(), using the RGEs.  The Jacobian matrix for
  /// this transformation may be accessed by calling
  /// displayInverseRGFlowJacobian().
  double NmssmJacobian::calcFTInverseJacobianLogs(const NmssmSoftsusy& m) {
    return calcFTInverseJacobianLogs(m, m.displayMxBC());
  }

  /// Calculates the Jacobian as in
  /// calcFTInverseJacobianLogs(const NmssmSoftsusy& m), but with \c mx set to
  /// the given value instead of that returned by \c m.displayMxBC().
  double NmssmJacobian::calcFTInverseJacobianLogs(const NmssmSoftsusy& m, double mx) {
    NmssmSoftsusy model(m);

    model.runto(model.calcMs());

    const double rgDet = calcRGFlowJacobian(model, mx, model.displayMu());
    const double ewsbDet = calcInverseEWSBJacobianLogs(model);

    return ewsbDet * rgDet;
  }

  /// Calculates the Jacobian of the form \f$ J = |
  /// \partial p / \partial O | \f$.  The transformation,
  /// which is the inverse of that calculated by
  /// calcFTInverseJacobian(), is done in two stages.  In the
  /// first, the parameters at the input scale \c mx, obtained
  /// from a call to \c m.displayMxBC(), are transformed to parameters
  /// \f$\{q\}\f$ at the SUSY scale, obtained from \c m.calcMs(),
  /// using the RGEs.  The Jacobian matrix for this
  /// transformation may be accessed afterwards by calling
  /// displayRGFlowJacobian().  Then the parameters
  /// \f$\{q\}\f$ are traded for the observables at the same
  /// scale using the EWSB conditions.  The Jacobian matrix for
  /// this second transformation can be accessed using
  /// displayEWSBJacobian().
  double NmssmJacobian::calcFTJacobian(const NmssmSoftsusy& m) {
    return calcFTJacobian(m, m.displayMxBC());
  }

  /// Calculates the Jacobian as in calcFTJacobian(const NmssmSoftsusy& m),
  /// but with \c mx set to the given value instead of that returned by
  /// \c m.displayMxBC().
  double NmssmJacobian::calcFTJacobian(const NmssmSoftsusy& m, double mx) {
    NmssmSoftsusy model(m);

    model.runto(model.calcMs());

    const double scale = model.displayMu();
    const double rgDet = calcRGFlowJacobian(model, scale, mx);
    const double ewsbDet = calcEWSBJacobian(model);

    return ewsbDet * rgDet;
  }

  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double NmssmJacobian::calcMz(NmssmSoftsusy& model, bool getRunningMass) {

    double mz = model.displayMzRun();
    if (!getRunningMass) {
      const double scale = model.displayMu();
      const double pizzt = model.piZZT(mz, scale);
      mz = sqrt(fabs(sqr(mz) - pizzt));
    }

    return mz;
  }

  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double NmssmJacobian::calcMt(NmssmSoftsusy& model, bool getRunningMass) {

    double mt = model.displayDrBarPars().mt;
    if (!getRunningMass) {
      const double scale = model.displayMu();

      // workaround hard-coded external momentum dependence
      QedQcd savedData(model.displayDataSet());
      QedQcd tempData;
      tempData.setPoleMt(mt);
      model.setData(tempData);

      const double oneLoopQCD = model.calcRunMtQCD();
      const double stopGluino = model.calcRunMtStopGluino();
      const double higgs = model.calcRunMtHiggs();
      const double neutralinos = model.calcRunMtNeutralinos();
      const double charginos = model.calcRunMtCharginos();

      double resigmat = mt * (stopGluino + higgs + neutralinos
                              + charginos + oneLoopQCD) / (16.0 * sqr(PI));

      const double g3Sq = sqr(model.displayGaugeCoupling(3));
      const double logMtSqOverQSq = 2.0 * log(mt / scale);
      const double twoLoopQCD = sqr(g3Sq) * (0.005191204615668296
        - 0.0032883224409535764 * logMtSqOverQSq + 0.0008822328500119351
        * sqr(logMtSqOverQSq));
      resigmat -= mt * twoLoopQCD;

      mt -= resigmat;

      model.setData(savedData);
    }

    return mt;
  }

  double NmssmJacobian::calcRunningParameter(double x, void* parameters) {
    RGFlowPars* pars = static_cast<RGFlowPars*>(parameters);

    NmssmSoftsusy* tempModel = pars->model;
    Parameters independent = pars->independent;
    Parameters dependent = pars->dependent;
    const double endScale = pars->toScale;

    const DoubleVector savedPars(tempModel->display());
    const double startScale = tempModel->displayMu();

    ostringstream msg;
    if (PRINTOUT > 1) msg << "# Qstart = " << startScale << ", ";

    switch (independent) {
    case Lambda: {
      if (pars->useSugraTrilinears) {
        const double Alambda = tempModel->displaySoftAlambda();
        tempModel->setTrialambda(x * Alambda);
      }
      tempModel->setLambda(x);
      if (PRINTOUT > 1) msg << "lambda = " << x << ", ";
      break;
    }
    case Kappa: {
      if (pars->useSugraTrilinears) {
        const double Akappa = tempModel->displaySoftAkappa();
        tempModel->setTriakappa(x * Akappa);
      }
      tempModel->setKappa(x);
      if (PRINTOUT > 1) msg << "kappa = " << x << ", ";
      break;
    }
    case SMu: {
      tempModel->setSusyMu(x);
      if (PRINTOUT > 1) msg << "mu = " << x << ", ";
      break;
    }
    case M3Sq: {
      tempModel->setM3Squared(x);
      if (PRINTOUT > 1) msg << "m3sq = " << x << ", ";
      break;
    }
    case XiS: {
      tempModel->setXiS(x);
      if (PRINTOUT > 1) msg << "xiS = " << x << ", ";
      break;
    }
    case Mh1Sq: {
      tempModel->setMh1Squared(x);
      if (PRINTOUT > 1) msg << "mH1Sq = " << x << ", ";
      break;
    }
    case Mh2Sq: {
      tempModel->setMh2Squared(x);
      if (PRINTOUT > 1) msg << "mH2Sq = " << x << ", ";
      break;
    }
    case MsSq: {
      tempModel->setMsSquared(x);
      if (PRINTOUT > 1) msg << "mSsq = " << x << ", ";
      break;
    }
    case Yt: {
      if (pars->useSugraTrilinears) {
        const double At = tempModel->displaySoftA(UA, 3, 3);
        tempModel->setTrilinearElement(UA, 3, 3, x * At);
      }
      tempModel->setYukawaElement(YU, 3, 3, x);
      if (PRINTOUT > 1) msg << "ht = " << x << ", ";
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcRunningParameter called with incorrect"
         << " independent parameter " << independent << '\n';
      throw ii.str();
    }
    }

    tempModel->runto(endScale);

    if (PRINTOUT > 1) msg << "Qend = " << endScale << ", ";

    double output;
    switch (dependent) {
    case Lambda: {
      output = tempModel->displayLambda();
      if (PRINTOUT > 1) msg << "lambda = " << output << '\n';
      break;
    }
    case Kappa: {
      output = tempModel->displayKappa();
      if (PRINTOUT > 1) msg << "kappa = " << output << '\n';
      break;
    }
    case SMu: {
      output = tempModel->displaySusyMu();
      if (PRINTOUT > 1) msg << "mu = " << output << '\n';
      break;
    }
    case M3Sq: {
      output = tempModel->displayM3Squared();
      if (PRINTOUT > 1) msg << "m3sq = " << output << '\n';
      break;
    }
    case XiS: {
      output = tempModel->displayXiS();
      if (PRINTOUT > 1) msg << "xiS = " << output << '\n';
      break;
    }
    case Mh1Sq: {
      output = tempModel->displayMh1Squared();
      if (PRINTOUT > 1) msg << "mH1Sq = " << output << '\n';
      break;
    }
    case Mh2Sq: {
      output = tempModel->displayMh2Squared();
      if (PRINTOUT > 1) msg << "mH2Sq = " << output << '\n';
      break;
    }
    case MsSq: {
      output = tempModel->displayMsSquared();
      if (PRINTOUT > 1) msg << "mSsq = " << output << '\n';
      break;
    }
    case Yt: {
      output = tempModel->displayYukawaElement(YU, 3, 3);
      if (PRINTOUT > 1) msg << "ht = " << output << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcRunningParameter called with incorrect"
         << " dependent parameter " << dependent << '\n';
      throw ii.str();
    }
    }

    if (PRINTOUT > 1) cout << msg.str();

    tempModel->setMu(startScale);
    tempModel->set(savedPars);

    return output;
  }

  std::pair<double,double> NmssmJacobian::calcRGDerivative(
    NmssmSoftsusy& model, Parameters dep, Parameters indep, double toScale) const {

    double x = 0.;
    double h = 0.01;

    switch (indep) {
    case Lambda: {
      x = model.displayLambda();
      h = 0.0005 * x;
      break;
    }
    case Kappa: {
      x = model.displayKappa();
      h = 0.0005 * x;
      break;
    }
    case SMu: {
      x = model.displaySusyMu();
      h = 0.01 * x;
      break;
    }
    case M3Sq: {
      x = model.displayM3Squared();
      h = 0.01 * x;
      break;
    }
    case XiS: {
      x = model.displayXiS();
      h = 0.01 * x;
      break;
    }
    case Mh1Sq: {
      x = model.displayMh1Squared();
      h = 0.01 * x;
      break;
    }
    case Mh2Sq: {
      x = model.displayMh2Squared();
      h = 0.01 * x;
      break;
    }
    case MsSq: {
      x = model.displayMsSquared();
      h = 0.01 * x;
      break;
    }
    case Yt: {
      x = model.displayYukawaElement(YU, 3, 3);
      h = 0.0005 * x;
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcRGDerivative called with incorrect"
         << " independent parameter " << indep << '\n';
      throw ii.str();
    }
    }

#ifdef ENABLE_GSL
    h = sqrt(std::numeric_limits<double>::epsilon()) * maximum(fabs(x), 1.0);
#endif

    volatile const double temp = x + h;
    h = temp - x;

    RGFlowPars pars;
    pars.model = &model;
    pars.independent = indep;
    pars.dependent = dep;
    pars.toScale = toScale;
    pars.useSugraTrilinears = useSugraTrilinears;

    double derivative = 0.;
    double err = 0.;

    bool has_error = false;

#ifdef ENABLE_GSL
    gsl_function func;
    func.function = &calcRunningParameter;
    func.params = &pars;

    gsl_deriv_central(&func, x, h, &derivative, &err);

    if (fabs(x) > 1.0e-10 &&
        (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
      has_error = fabs(err / derivative) > 1.0;
    }

    if (has_error) {
      // attempt to use a forward or backward difference in
      // case we are at a boundary at parameter space
      gsl_deriv_forward(&func, x, h, &derivative, &err);

      if (fabs(x) > 1.0e-10 &&
          (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
        has_error = fabs(err / derivative) > 1.0;
      }

      if (has_error) {
        gsl_deriv_backward(&func, x, h, &derivative, &err);

        if (fabs(x) > 1.0e-10 &&
            (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
          has_error = fabs(err / derivative) > 1.0;
        }
      }
    }
#else
    if (fabs(x) > 1.0e-10) {
      derivative = calcDerivative(calcRunningParameter, x, h, &err, &pars);

      if (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10) {
        has_error = fabs(err / derivative) > 1.0;
      }
    }
#endif

    if (has_error) {
      derivative = -numberOfTheBeast;
    }

    if (PRINTOUT > 1) {
      const double scale = model.displayMu();

      ostringstream msg;
      msg << "# derivative[d";

      switch (dep) {
      case Lambda: msg << "lambda(Q=" << toScale << ")/d"; break;
      case Kappa: msg << "kappa(Q=" << toScale << ")/d"; break;
      case SMu: msg << "mu(Q=" << toScale << ")/d"; break;
      case M3Sq: msg << "m3sq(Q=" << toScale << ")/d"; break;
      case XiS: msg << "xiS(Q=" << toScale << ")/d"; break;
      case Mh1Sq: msg << "mH1Sq(Q=" << toScale << ")/d"; break;
      case Mh2Sq: msg << "mH2Sq(Q=" << toScale << ")/d"; break;
      case MsSq: msg << "mSsq(Q=" << toScale << ")/d"; break;
      case Yt: msg << "ht(Q=" << toScale << ")/d"; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcRGDerivative called with incorrect"
           << " dependent parameter " << dep << '\n';
        throw ii.str();
      }
      }

      switch (indep) {
      case Lambda: msg << "lambda(Q=" << scale << ")] = "; break;
      case Kappa: msg << "kappa(Q=" << scale << ")] = "; break;
      case SMu: msg << "mu(Q=" << scale << ")] = "; break;
      case M3Sq: msg << "m3sq(Q=" << scale << ")] = "; break;
      case XiS: msg << "xiS(Q=" << scale << ")] = "; break;
      case Mh1Sq: msg << "mH1Sq(Q=" << scale << ")] = "; break;
      case Mh2Sq: msg << "mH2Sq(Q=" << scale << ")] = "; break;
      case MsSq: msg << "mSsq(Q=" << scale << ")] = "; break;
      case Yt: msg << "ht(Q=" << scale << ")] = "; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcRGDerivative called with incorrect"
           << " independent parameter " << indep << '\n';
        throw ii.str();
      }
      }

      msg << derivative << ", error = " << err << '\n';
      cout << msg.str();
    }

    return std::pair<double,double>(derivative, err);
  }

  double NmssmJacobian::calcLogRunningParameter(double x, void* parameters) {
    RGFlowPars* pars = static_cast<RGFlowPars*>(parameters);

    NmssmSoftsusy* tempModel = pars->model;
    Parameters independent = pars->independent;
    Parameters dependent = pars->dependent;
    const double endScale = pars->toScale;

    const DoubleVector savedPars(tempModel->display());
    const double startScale = tempModel->displayMu();

    ostringstream msg;
    if (PRINTOUT > 1) msg << "# Qstart = " << startScale << ", ";

    const double value = std::exp(x);
    switch (independent) {
    case Lambda: {
      if (pars->useSugraTrilinears) {
        const double Alambda = tempModel->displaySoftAlambda();
        tempModel->setTrialambda(value * Alambda);
      }
      tempModel->setLambda(value);
      if (PRINTOUT > 1) msg << "lambda = " << value << ", ";
      break;
    }
    case Kappa: {
      if (pars->useSugraTrilinears) {
        const double Akappa = tempModel->displaySoftAkappa();
        tempModel->setTriakappa(value * Akappa);
      }
      tempModel->setKappa(value);
      if (PRINTOUT > 1) msg << "kappa = " << value << ", ";
      break;
    }
    case SMu: {
      tempModel->setSusyMu(value);
      if (PRINTOUT > 1) msg << "mu = " << value << ", ";
      break;
    }
    case M3Sq: {
      tempModel->setM3Squared(value);
      if (PRINTOUT > 1) msg << "m3sq = " << value << ", ";
      break;
    }
    case XiS: {
      tempModel->setXiS(value);
      if (PRINTOUT > 1) msg << "xiS = " << value << ", ";
      break;
    }
    case Mh1Sq: {
      tempModel->setMh1Squared(value);
      if (PRINTOUT > 1) msg << "mH1Sq = " << value << ", ";
      break;
    }
    case Mh2Sq: {
      tempModel->setMh2Squared(value);
      if (PRINTOUT > 1) msg << "mH2Sq = " << value << ", ";
      break;
    }
    case MsSq: {
      tempModel->setMsSquared(value);
      if (PRINTOUT > 1) msg << "mSsq = " << value << ", ";
      break;
    }
    case Yt: {
      if (pars->useSugraTrilinears) {
        const double At = tempModel->displaySoftA(UA, 3, 3);
        tempModel->setTrilinearElement(UA, 3, 3, value * At);
      }
      tempModel->setYukawaElement(YU, 3, 3, value);
      if (PRINTOUT > 1) msg << "ht = " << value << ", ";
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcLogRunningParameter called with incorrect"
         << " independent parameter " << independent << '\n';
      throw ii.str();
    }
    }

    tempModel->runto(endScale);

    if (PRINTOUT > 1) msg << "Qend = " << endScale << ", ";

    double output;
    switch (dependent) {
    case Lambda: {
      output = tempModel->displayLambda();
      if (PRINTOUT > 1) msg << "lambda = " << output << '\n';
      break;
    }
    case Kappa: {
      output = tempModel->displayKappa();
      if (PRINTOUT > 1) msg << "kappa = " << output << '\n';
      break;
    }
    case SMu: {
      output = tempModel->displaySusyMu();
      if (PRINTOUT > 1) msg << "mu = " << output << '\n';
      break;
    }
    case M3Sq: {
      output = tempModel->displayM3Squared();
      if (PRINTOUT > 1) msg << "m3sq = " << output << '\n';
      break;
    }
    case XiS: {
      output = tempModel->displayXiS();
      if (PRINTOUT > 1) msg << "xiS = " << output << '\n';
      break;
    }
    case Mh1Sq: {
      output = tempModel->displayMh1Squared();
      if (PRINTOUT > 1) msg << "mH1Sq = " << output << '\n';
      break;
    }
    case Mh2Sq: {
      output = tempModel->displayMh2Squared();
      if (PRINTOUT > 1) msg << "mH2Sq = " << output << '\n';
      break;
    }
    case MsSq: {
      output = tempModel->displayMsSquared();
      if (PRINTOUT > 1) msg << "mSsq = " << output << '\n';
      break;
    }
    case Yt: {
      output = tempModel->displayYukawaElement(YU, 3, 3);
      if (PRINTOUT > 1) msg << "ht = " << output << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcLogRunningParameter called with incorrect"
         << " dependent parameter " << dependent << '\n';
      throw ii.str();
    }
    }

    if (PRINTOUT > 1) cout << msg.str();

    tempModel->setMu(startScale);
    tempModel->set(savedPars);

    return std::log(output);
  }

  std::pair<double,double> NmssmJacobian::calcRGLogDerivative(
    NmssmSoftsusy& model, Parameters dep, Parameters indep, double toScale) const {

    double x = 0.;
    double h = 0.01;

    switch (indep) {
    case Lambda: {
      x = std::log(model.displayLambda());
      h = 0.0005 * x;
      break;
    }
    case Kappa: {
      x = std::log(model.displayKappa());
      h = 0.0005 * x;
      break;
    }
    case SMu: {
      x = std::log(model.displaySusyMu());
      h = 0.01 * x;
      break;
    }
    case M3Sq: {
      x = std::log(model.displayM3Squared());
      h = 0.01 * x;
      break;
    }
    case XiS: {
      x = std::log(model.displayXiS());
      h = 0.01 * x;
      break;
    }
    case Mh1Sq: {
      x = std::log(model.displayMh1Squared());
      h = 0.01 * x;
      break;
    }
    case Mh2Sq: {
      x = std::log(model.displayMh2Squared());
      h = 0.01 * x;
      break;
    }
    case MsSq: {
      x = std::log(model.displayMsSquared());
      h = 0.01 * x;
      break;
    }
    case Yt: {
      x = std::log(model.displayYukawaElement(YU, 3, 3));
      h = 0.0005 * x;
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcRGLogDerivative called with incorrect"
         << " independent parameter " << indep << '\n';
      throw ii.str();
    }
    }

#ifdef ENABLE_GSL
    h = sqrt(std::numeric_limits<double>::epsilon()) * maximum(fabs(x), 1.0);
#endif

    volatile const double temp = x + h;
    h = temp - x;

    RGFlowPars pars;
    pars.model = &model;
    pars.independent = indep;
    pars.dependent = dep;
    pars.toScale = toScale;
    pars.useSugraTrilinears = useSugraTrilinears;

    double derivative = 0.;
    double err = 0.;

    bool has_error = false;

#ifdef ENABLE_GSL
    gsl_function func;
    func.function = &calcLogRunningParameter;
    func.params = &pars;

    gsl_deriv_central(&func, x, h, &derivative, &err);

    if (fabs(x) > 1.0e-10 &&
        (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
      has_error = fabs(err / derivative) > 1.0;
    }

    if (has_error) {
      // attempt to use a forward or backward difference in
      // case we are at a boundary in parameter space
      gsl_deriv_forward(&func, x, h, &derivative, &err);

      if (fabs(x) > 1.0e-10 &&
          (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
        has_error = fabs(err / derivative) > 1.0;
      }

      if (has_error) {
        gsl_deriv_backward(&func, x, h, &derivative, &err);

        if (fabs(x) > 1.0e-10 &&
            (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
          has_error = fabs(err / derivative) > 1.0;
        }
      }
    }
#else
    if (fabs(x) > 1.0e-10) {
      derivative = calcDerivative(calcLogRunningParameter, x, h, &err, &pars);

      if (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10) {
        has_error = fabs(err / derivative) > 1.0;
      }
    }
#endif

    if (has_error) {
      derivative = -numberOfTheBeast;
    }

    if (PRINTOUT > 1) {
      const double scale = model.displayMu();

      ostringstream msg;
      msg << "# derivative[d";

      switch (dep) {
      case Lambda: msg << "log(lambda(Q=" << toScale << "))/d"; break;
      case Kappa: msg << "log(kappa(Q=" << toScale << "))/d"; break;
      case SMu: msg << "log(mu(Q=" << toScale << "))/d"; break;
      case M3Sq: msg << "log(m3sq(Q=" << toScale << "))/d"; break;
      case XiS: msg << "log(xiS(Q=" << toScale << "))/d"; break;
      case Mh1Sq: msg << "log(mH1Sq(Q=" << toScale << "))/d"; break;
      case Mh2Sq: msg << "log(mH2Sq(Q=" << toScale << "))/d"; break;
      case MsSq: msg << "log(mSsq(Q=" << toScale << "))/d"; break;
      case Yt: msg << "log(ht(Q=" << toScale << "))/d"; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcRGLogDerivative called with incorrect"
           << " dependent parameter " << dep << '\n';
        throw ii.str();
      }
      }

      switch (indep) {
      case Lambda: msg << "log(lambda(Q=" << scale << "))] = "; break;
      case Kappa: msg << "log(kappa(Q=" << scale << "))] = "; break;
      case SMu: msg << "log(mu(Q=" << scale << "))] = "; break;
      case M3Sq: msg << "log(m3sq(Q=" << scale << "))] = "; break;
      case XiS: msg << "log(xiS(Q=" << scale << "))] = "; break;
      case Mh1Sq: msg << "log(mH1Sq(Q=" << scale << "))] = "; break;
      case Mh2Sq: msg << "log(mH2Sq(Q=" << scale << "))] = "; break;
      case MsSq: msg << "log(mSsq(Q=" << scale << "))] = "; break;
      case Yt: msg << "log(ht(Q=" << scale << "))] = "; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcRGLogDerivative called with incorrect"
           << " independent parameter " << indep << '\n';
        throw ii.str();
      }
      }

      msg << derivative << ", error = " << err << '\n';
      cout << msg.str();
    }

    return std::pair<double,double>(derivative, err);
  }

  double NmssmJacobian::calcRGFlowJacobian(NmssmSoftsusy& model,
                                           double startScale, double endScale) {
    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();

    model.runto(startScale);

    const int numPars = includeTop ? 4 : 3;

    DoubleMatrix jac(numPars, numPars);
    DoubleMatrix jacErrors(numPars, numPars);

    std::vector<Parameters> indepPars;
    std::vector<double> paramValues;

    if (SoftHiggsOut) {
      indepPars.push_back(Mh1Sq);
      indepPars.push_back(Mh2Sq);
      indepPars.push_back(MsSq);

      paramValues.push_back(model.displayMh1Squared());
      paramValues.push_back(model.displayMh2Squared());
      paramValues.push_back(model.displayMsSquared());
    } else if ( model.displayZ3() ) {
      indepPars.push_back(Lambda);
      indepPars.push_back(Kappa);
      indepPars.push_back(MsSq);

      paramValues.push_back(model.displayLambda());
      paramValues.push_back(model.displayKappa());
      paramValues.push_back(model.displayMsSquared());
    } else {
      indepPars.push_back(SMu);
      indepPars.push_back(M3Sq);
      indepPars.push_back(XiS);

      paramValues.push_back(model.displaySusyMu());
      paramValues.push_back(model.displayM3Squared());
      paramValues.push_back(model.displayXiS());
    }
    if (includeTop) {
      indepPars.push_back(Yt);
      paramValues.push_back(model.displayYukawaElement(YU, 3, 3));
    }

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      std::pair<double,double> result;
      if (SoftHiggsOut) {
        result = calcRGDerivative(model, Mh1Sq, indepPars[i], endScale);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcRGDerivative(model, Mh2Sq, indepPars[i], endScale);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcRGDerivative(model, MsSq, indepPars[i], endScale);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      } else if ( model.displayZ3() ) {
        result = calcRGDerivative(model, Lambda, indepPars[i], endScale);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcRGDerivative(model, Kappa, indepPars[i], endScale);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcRGDerivative(model, MsSq, indepPars[i], endScale);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      } else {
        result = calcRGDerivative(model, SMu, indepPars[i], endScale);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcRGDerivative(model, M3Sq, indepPars[i], endScale);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcRGDerivative(model, XiS, indepPars[i], endScale);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      }
      if (includeTop) {
        result = calcRGDerivative(model, Yt, indepPars[i], endScale);
        jac(i + 1, 4) = result.first;
        jacErrors(i + 1, 4) = result.second;
      }
    }

    // check for errors
    hasError = checkDerivativeErrors(jac, jacErrors, paramValues);

    // save calculated matrix
    // convention: inverse refers to case where transformation
    // is from high-scale to low-scale parameters
    if (startScale == endScale) {
      if (invJacRGFlow.displayRows() != numPars
          || invJacRGFlow.displayCols() != numPars) {
        invJacRGFlow.resize(numPars, numPars);
      }
      invJacRGFlow = jac;

      if (invJacRGFlowErrors.displayRows() != numPars
          || invJacRGFlowErrors.displayCols() != numPars) {
        invJacRGFlowErrors.resize(numPars, numPars);
      }
      invJacRGFlowErrors = jacErrors;

      if (jacRGFlow.displayRows() != numPars
          || jacRGFlow.displayCols() != numPars) {
        jacRGFlow.resize(numPars, numPars);
      }
      jacRGFlow = jac;

      if (jacRGFlowErrors.displayRows() != numPars
          || jacRGFlowErrors.displayCols() != numPars) {
        jacRGFlowErrors.resize(numPars, numPars);
      }
      jacRGFlowErrors = jacErrors;
    } else if (startScale > endScale) {
      if (invJacRGFlow.displayRows() != numPars
          || invJacRGFlow.displayCols() != numPars) {
        invJacRGFlow.resize(numPars, numPars);
      }
      invJacRGFlow = jac;

      if (invJacRGFlowErrors.displayRows() != numPars
          || invJacRGFlowErrors.displayCols() != numPars) {
        invJacRGFlowErrors.resize(numPars, numPars);
      }
      invJacRGFlowErrors = jacErrors;
    } else {
      if (jacRGFlow.displayRows() != numPars
          || jacRGFlow.displayCols() != numPars) {
        jacRGFlow.resize(numPars, numPars);
      }
      jacRGFlow = jac;

      if (jacRGFlowErrors.displayRows() != numPars
          || jacRGFlowErrors.displayCols() != numPars) {
        jacRGFlowErrors.resize(numPars, numPars);
      }
      jacRGFlowErrors = jacErrors;
    }

    model.setMu(scale);
    model.set(savedPars);

    return jac.determinant();
  }

  // @todo handle negative parameters properly
  double NmssmJacobian::calcRGFlowJacobianLogs(
    NmssmSoftsusy& model, double startScale, double endScale) {
    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();

    model.runto(startScale);

    const int numPars = includeTop ? 4 : 3;

    DoubleMatrix jac(numPars, numPars);
    DoubleMatrix jacErrors(numPars, numPars);

    std::vector<Parameters> indepPars;
    std::vector<double> paramValues;

    if (SoftHiggsOut) {
      indepPars.push_back(Mh1Sq);
      indepPars.push_back(Mh2Sq);
      indepPars.push_back(MsSq);

      paramValues.push_back(std::log(model.displayMh1Squared()));
      paramValues.push_back(std::log(model.displayMh2Squared()));
      paramValues.push_back(std::log(model.displayMsSquared()));
    } else if ( model.displayZ3() ) {
      indepPars.push_back(Lambda);
      indepPars.push_back(Kappa);
      indepPars.push_back(MsSq);

      paramValues.push_back(std::log(model.displayLambda()));
      paramValues.push_back(std::log(model.displayKappa()));
      paramValues.push_back(std::log(model.displayMsSquared()));
    } else {
      indepPars.push_back(SMu);
      indepPars.push_back(M3Sq);
      indepPars.push_back(XiS);

      paramValues.push_back(std::log(model.displaySusyMu()));
      paramValues.push_back(std::log(model.displayM3Squared()));
      paramValues.push_back(std::log(model.displayXiS()));
    }
    if (includeTop) {
      indepPars.push_back(Yt);
      paramValues.push_back(std::log(model.displayYukawaElement(YU, 3, 3)));
    }

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      std::pair<double,double> result;
      if (SoftHiggsOut) {
        result = calcRGLogDerivative(model, Mh1Sq, indepPars[i], endScale);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcRGLogDerivative(model, Mh2Sq, indepPars[i], endScale);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcRGLogDerivative(model, MsSq, indepPars[i], endScale);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      } else if ( model.displayZ3() ) {
        result = calcRGLogDerivative(model, Lambda, indepPars[i], endScale);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcRGLogDerivative(model, Kappa, indepPars[i], endScale);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcRGLogDerivative(model, MsSq, indepPars[i], endScale);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      } else {
        result = calcRGLogDerivative(model, SMu, indepPars[i], endScale);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcRGLogDerivative(model, M3Sq, indepPars[i], endScale);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcRGLogDerivative(model, XiS, indepPars[i], endScale);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      }
      if (includeTop) {
        result = calcRGLogDerivative(model, Yt, indepPars[i], endScale);
        jac(i + 1, 4) = result.first;
        jacErrors(i + 1, 4) = result.second;
      }
    }

    // check for errors
    hasError = checkDerivativeErrors(jac, jacErrors, paramValues);

    // save calculated matrix
    // convention: inverse refers to case where transformation
    // is from high-scale to low-scale parameters
    if (startScale == endScale) {
      if (invJacRGFlow.displayRows() != numPars
          || invJacRGFlow.displayCols() != numPars) {
        invJacRGFlow.resize(numPars, numPars);
      }
      invJacRGFlow = jac;

      if (invJacRGFlowErrors.displayRows() != numPars
          || invJacRGFlowErrors.displayCols() != numPars) {
        invJacRGFlowErrors.resize(numPars, numPars);
      }
      invJacRGFlowErrors = jacErrors;

      if (jacRGFlow.displayRows() != numPars
          || jacRGFlow.displayCols() != numPars) {
        jacRGFlow.resize(numPars, numPars);
      }
      jacRGFlow = jac;

      if (jacRGFlowErrors.displayRows() != numPars
          || jacRGFlowErrors.displayCols() != numPars) {
        jacRGFlowErrors.resize(numPars, numPars);
      }
      jacRGFlowErrors = jacErrors;
    } else if (startScale > endScale) {
      if (invJacRGFlow.displayRows() != numPars
          || invJacRGFlow.displayCols() != numPars) {
        invJacRGFlow.resize(numPars, numPars);
      }
      invJacRGFlow = jac;

      if (invJacRGFlowErrors.displayRows() != numPars
          || invJacRGFlowErrors.displayCols() != numPars) {
        invJacRGFlowErrors.resize(numPars, numPars);
      }
      invJacRGFlowErrors = jacErrors;
    } else {
      if (jacRGFlow.displayRows() != numPars
          || jacRGFlow.displayCols() != numPars) {
        jacRGFlow.resize(numPars, numPars);
      }
      jacRGFlow = jac;

      if (jacRGFlowErrors.displayRows() != numPars
          || jacRGFlowErrors.displayCols() != numPars) {
        jacRGFlowErrors.resize(numPars, numPars);
      }
      jacRGFlowErrors = jacErrors;
    }

    model.setMu(scale);
    model.set(savedPars);

    return jac.determinant();
  }

  double NmssmJacobian::calcEWSBOutput(double x, void* parameters) {
    EWSBPars* pars = static_cast<EWSBPars*>(parameters);

    NmssmSoftsusy* tempModel = pars->model;
    Parameters independent = pars->independent;
    Parameters dependent = pars->dependent;

    const DoubleVector savedPars(tempModel->display());
    const drBarPars savedDrBarPars(tempModel->displayDrBarPars());
    // reset problems to avoid subsequent iterations skipping
    // calculating tadpoles
    const sProblem savedProblems(tempModel->displayProblem());
    tempModel->setProblem(sProblem());

    ostringstream msg;
    if (PRINTOUT > 1) msg << "# ";

    switch (independent) {
    case Lambda: {
      if (pars->useSugraTrilinears) {
        const double Alambda = tempModel->displaySoftAlambda();
        tempModel->setTrialambda(x * Alambda);
      }
      tempModel->setLambda(x);
      if (PRINTOUT > 1) msg << "lambda = " << x << ", ";
      break;
    }
    case Kappa: {
      if (pars->useSugraTrilinears) {
        const double Akappa = tempModel->displaySoftAkappa();
        tempModel->setTriakappa(x * Akappa);
      }
      tempModel->setKappa(x);
      if (PRINTOUT > 1) msg << "kappa = " << x << ", ";
      break;
    }
    case SMu: {
      tempModel->setSusyMu(x);
      if (PRINTOUT > 1) msg << "mu = " << x << ", ";
      break;
    }
    case M3Sq: {
      tempModel->setM3Squared(x);
      if (PRINTOUT > 1) msg << "m3sq = " << x << ", ";
      break;
    }
    case XiS: {
      tempModel->setXiS(x);
      if (PRINTOUT > 1) msg << "xiS = " << x << ", ";
      break;
    }
    case Mh1Sq: {
      tempModel->setMh1Squared(x);
      if (PRINTOUT > 1) msg << "mH1Sq = " << x << ", ";
      break;
    }
    case Mh2Sq: {
      tempModel->setMh2Squared(x);
      if (PRINTOUT > 1) msg << "mH2Sq = " << x << ", ";
      break;
    }
    case MsSq: {
      tempModel->setMsSquared(x);
      if (PRINTOUT > 1) msg << "mSsq = " << x << ", ";
      break;
    }
    case Yt: {
      if (pars->useSugraTrilinears) {
        const double At = tempModel->displaySoftA(UA, 3, 3);
        tempModel->setTrilinearElement(UA, 3, 3, x * At);
      }
      tempModel->setYukawaElement(YU, 3, 3, x);
      if (PRINTOUT > 1) msg << "ht = " << x << ", ";
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBOutput called with incorrect"
         << " independent parameter " << independent << '\n';
      throw ii.str();
    }
    }

    tempModel->calcDrBarPars();
    tempModel->doTadpoles(tempModel->displayDrBarPars().mt,
                          tempModel->calcSinthdrbar());

    DoubleVector vevs(3);
    vevs(1) = tempModel->displayHvev();
    vevs(2) = tempModel->displayTanb();
    vevs(3) = tempModel->displaySvev();

    int error = 0;
    tempModel->predVevs(vevs, error);

    if (error != 0) {
      pars->failedVevIteration = true;
      if (PRINTOUT > 0) {
        cout << "# Warning: could not solve for VEVs\n";
      }
    }

    tempModel->setHvev(vevs(1));
    tempModel->setTanb(vevs(2));
    tempModel->setSvev(vevs(3));

    tempModel->calcDrBarPars();
    const double mtrun = tempModel->displayDrBarPars().mt;
    const double sinthDRbar = tempModel->calcSinthdrbar();
    tempModel->doTadpoles(mtrun, sinthDRbar);

    if (tempModel->displayProblem().inaccurateHiggsMass) {
      pars->inaccurateHiggsMass = true;
      if (fabs(tempModel->displayDrBarHiggsAccuracy())
          > fabs(pars->higgsMassAccuracy))
        pars->higgsMassAccuracy = tempModel->displayDrBarHiggsAccuracy();
    }

    double output;
    switch (dependent) {
    case Mzsq: {
      output = sqr(calcMz(*tempModel, pars->useRunningMasses));
      if (PRINTOUT > 1) msg << "MZ = " << sqrt(output) << '\n';
      break;
    }
    case Tanb: {
      output = tempModel->displayTanb();
      if (PRINTOUT > 1) msg << "tanb = " << output << '\n';
      break;
    }
    case Svev: {
      output = tempModel->displaySvev();
      if (PRINTOUT > 1) msg << "Svev = " << output << '\n';
      break;
    }
    case Mtsq: {
      output = sqr(calcMt(*tempModel, pars->useRunningMasses));
      if (PRINTOUT > 1) msg << "MT = " << sqrt(output) << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBOutput called with incorrect"
         << " dependent parameter " << dependent << '\n';
      throw ii.str();
    }
    }

    if (PRINTOUT > 1) cout << msg.str();

    tempModel->set(savedPars);
    tempModel->setDrBarPars(savedDrBarPars);
    tempModel->setProblem(savedProblems);

    return output;
  }

  double NmssmJacobian::calcLogEWSBOutput(double x, void* parameters) {
    EWSBPars* pars = static_cast<EWSBPars*>(parameters);

    NmssmSoftsusy* tempModel = pars->model;
    Parameters independent = pars->independent;
    Parameters dependent = pars->dependent;

    const DoubleVector savedPars(tempModel->display());
    const drBarPars savedDrBarPars(tempModel->displayDrBarPars());
    // reset problems to avoid subsequent iterations skipping
    // calculating tadpoles
    const sProblem savedProblems(tempModel->displayProblem());
    tempModel->setProblem(sProblem());

    ostringstream msg;
    if (PRINTOUT > 1) msg << "# ";

    const double value = std::exp(x);
    switch (independent) {
    case Lambda: {
      if (pars->useSugraTrilinears) {
        const double Alambda = tempModel->displaySoftAlambda();
        tempModel->setTrialambda(value * Alambda);
      }
      tempModel->setLambda(value);
      if (PRINTOUT > 1) msg << "lambda = " << value << ", ";
      break;
    }
    case Kappa: {
      if (pars->useSugraTrilinears) {
        const double Akappa = tempModel->displaySoftAkappa();
        tempModel->setTriakappa(value * Akappa);
      }
      tempModel->setKappa(value);
      if (PRINTOUT > 1) msg << "kappa = " << value << ", ";
      break;
    }
    case SMu: {
      tempModel->setSusyMu(value);
      if (PRINTOUT > 1) msg << "mu = " << value << ", ";
      break;
    }
    case M3Sq: {
      tempModel->setM3Squared(value);
      if (PRINTOUT > 1) msg << "m3sq = " << value << ", ";
      break;
    }
    case XiS: {
      tempModel->setXiS(value);
      if (PRINTOUT > 1) msg << "xiS = " << value << ", ";
      break;
    }
    case Mh1Sq: {
      tempModel->setMh1Squared(value);
      if (PRINTOUT > 1) msg << "mH1Sq = " << value << ", ";
      break;
    }
    case Mh2Sq: {
      tempModel->setMh2Squared(value);
      if (PRINTOUT > 1) msg << "mH2Sq = " << value << ", ";
      break;
    }
    case MsSq: {
      tempModel->setMsSquared(value);
      if (PRINTOUT > 1) msg << "mSsq = " << value << ", ";
      break;
    }
    case Yt: {
      if (pars->useSugraTrilinears) {
        const double At = tempModel->displaySoftA(UA, 3, 3);
        tempModel->setTrilinearElement(UA, 3, 3, value * At);
      }
      tempModel->setYukawaElement(YU, 3, 3, value);
      if (PRINTOUT > 1) msg << "ht = " << value << ", ";
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBOutput called with incorrect"
         << " independent parameter " << independent << '\n';
      throw ii.str();
    }
    }

    tempModel->calcDrBarPars();
    tempModel->doTadpoles(tempModel->displayDrBarPars().mt,
                          tempModel->calcSinthdrbar());

    DoubleVector vevs(3);
    vevs(1) = tempModel->displayHvev();
    vevs(2) = tempModel->displayTanb();
    vevs(3) = tempModel->displaySvev();

    int error = 0;
    tempModel->predVevs(vevs, error);

    if (error != 0) {
      pars->failedVevIteration = true;
      if (PRINTOUT > 0) {
        cout << "# Warning: could not solve for VEVs\n";
      }
    }

    tempModel->setHvev(vevs(1));
    tempModel->setTanb(vevs(2));
    tempModel->setSvev(vevs(3));

    tempModel->calcDrBarPars();
    const double mtrun = tempModel->displayDrBarPars().mt;
    const double sinthDRbar = tempModel->calcSinthdrbar();
    tempModel->doTadpoles(mtrun, sinthDRbar);

    if (tempModel->displayProblem().inaccurateHiggsMass) {
      pars->inaccurateHiggsMass = true;
      if (fabs(tempModel->displayDrBarHiggsAccuracy())
          > fabs(pars->higgsMassAccuracy))
        pars->higgsMassAccuracy = tempModel->displayDrBarHiggsAccuracy();
    }

    double output;
    switch (dependent) {
    case Mzsq: {
      output = sqr(calcMz(*tempModel, pars->useRunningMasses));
      if (PRINTOUT > 1) msg << "MZ = " << sqrt(output) << '\n';
      break;
    }
    case Tanb: {
      output = tempModel->displayTanb();
      if (PRINTOUT > 1) msg << "tanb = " << output << '\n';
      break;
    }
    case Svev: {
      output = tempModel->displaySvev();
      if (PRINTOUT > 1) msg << "Svev = " << output << '\n';
      break;
    }
    case Mtsq: {
      output = sqr(calcMt(*tempModel, pars->useRunningMasses));
      if (PRINTOUT > 1) msg << "MT = " << sqrt(output) << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcLogEWSBOutput called with incorrect"
         << " dependent parameter " << dependent << '\n';
      throw ii.str();
    }
    }

    if (PRINTOUT > 1) cout << msg.str();

    tempModel->set(savedPars);
    tempModel->setDrBarPars(savedDrBarPars);
    tempModel->setProblem(savedProblems);

    return std::log(output);
  }

  int NmssmJacobian::ewsbOutputErrors(const DoubleVector & guess,
                                      void* parameters,
                                      DoubleVector & errors) {
    EWSBPars* pars = static_cast<EWSBPars*>(parameters);

    NmssmSoftsusy* m = pars->model;
    const DoubleVector outputs = pars->outputs;
    const int numOutputs = outputs.displayEnd();

    if (SoftHiggsOut) {
      m->setMh1Squared(guess(1));
      m->setMh2Squared(guess(2));
      m->setMsSquared(guess(3));
    } else if ( m->displayZ3() ) {
      if (pars->useSugraTrilinears) {
        const double Alambda = m->displaySoftAlambda();
        const double Akappa = m->displaySoftAkappa();
        m->setTrialambda(guess(1) * Alambda);
        m->setTriakappa(guess(2) * Akappa);
      }
      m->setLambda(guess(1));
      m->setKappa(guess(2));
      m->setMsSquared(guess(3));
    } else {
      m->setSusyMu(guess(1));
      m->setM3Squared(guess(2));
      m->setXiS(guess(3));
    }
    if (numOutputs > 3) {
      if (pars->useSugraTrilinears) {
        const double At = m->displaySoftA(UA, 3, 3);
        m->setTrilinearElement(UA, 3, 3, guess(4) * At);
      }
      m->setYukawaElement(YU, 3, 3, guess(4));
    }

    DoubleVector vevs(3);
    vevs(1) = m->displayHvev();
    vevs(2) = m->displayTanb();
    vevs(3) = m->displaySvev();

    int error = 0;
    m->predVevs(vevs, error);

    if (error != 0) {
      pars->failedVevIteration = true;
      if (PRINTOUT > 0) {
        cout << "# Warning: could not solve for VEVs\n";
      }
    }

    m->setHvev(vevs(1));
    m->setTanb(vevs(2));
    m->setSvev(vevs(3));

    m->calcDrBarPars();

    if (m->displayProblem().inaccurateHiggsMass) {
      pars->inaccurateHiggsMass = true;
      if (fabs(m->displayDrBarHiggsAccuracy()) > fabs(pars->higgsMassAccuracy))
        pars->higgsMassAccuracy = m->displayDrBarHiggsAccuracy();
    }

    if (m->displayZ3() && !SoftHiggsOut) {
      errors(1) = 1.0 - (sqr(calcMz(*m, pars->useRunningMasses)) / outputs(1));
      errors(2) = 1.0 - (m->displayTanb() / outputs(2));
      errors(3) = 1.0 - (m->displayLambda() / outputs(3));
    } else {
      errors(1) = 1.0 - (sqr(calcMz(*m, pars->useRunningMasses)) / outputs(1));
      errors(2) = 1.0 - (m->displayTanb() / outputs(2));
      errors(3) = 1.0 - (m->displaySvev() / outputs(3));
    }

    error = error && testNan(errors(1)) && testNan(errors(2))
      && testNan(errors(3));

    if (numOutputs > 3) {
      errors(4) = 1.0 - (sqr(calcMt(*m, pars->useRunningMasses)) / outputs(4));
      error = error && testNan(errors(4));
    }

    return error;
  }

  void NmssmJacobian::fixEWSBOutputs(EWSBPars* pars, int & err) {

    NmssmSoftsusy* m = pars->model;
    const int numOutputs = pars->outputs.displayEnd();
    // initial guess
    DoubleVector guess(numOutputs);
    if (SoftHiggsOut) {
      guess(1) = m->displayMh1Squared();
      guess(2) = m->displayMh2Squared();
      guess(3) = m->displayMsSquared();
    } else if ( m->displayZ3() ) {
      guess(1) = m->displayLambda();
      guess(2) = m->displayKappa();
      guess(3) = m->displayMsSquared();
    } else {
      guess(1) = m->displaySusyMu();
      guess(2) = m->displayM3Squared();
      guess(3) = m->displayXiS();
    }
    if (numOutputs > 3) {
      guess(4) = m->displayYukawaElement(YU, 3, 3);
    }

    const double oldTol = TOLERANCE;
    TOLERANCE *= 1.0e-3;

    bool error = newt(guess, ewsbOutputErrors, pars);
    err = error ? 1 : 0;

    TOLERANCE = oldTol;
  }

  double NmssmJacobian::calcEWSBParameter(double x, void* parameters) {
    EWSBPars* pars = static_cast<EWSBPars*>(parameters);

    NmssmSoftsusy* tempModel = pars->model;
    Parameters independent = pars->independent;
    Parameters dependent = pars->dependent;

    const DoubleVector savedPars(tempModel->display());
    const double startScale = tempModel->displayMu();
    const drBarPars savedDrBarPars(tempModel->displayDrBarPars());
    const DoubleVector savedOutputs(pars->outputs);
    // reset problems to avoid subsequent iterations skipping
    // calculating tadpoles
    const sProblem savedProblems(tempModel->displayProblem());
    tempModel->setProblem(sProblem());

    ostringstream msg;
    if (PRINTOUT > 1) msg << "# ";

    switch (independent) {
    case Lambda: {
      pars->outputs(3) = x;
      if (PRINTOUT > 1) msg << "lambda = " << x << ", ";
      break;
    }
    case Mzsq: {
      pars->outputs(1) = x;
      if (PRINTOUT > 1) msg << "MZ = " << sqrt(x) << ", ";
      break;
    }
    case Tanb: {
      pars->outputs(2) = x;
      if (PRINTOUT > 1) msg << "tanb = " << x << ", ";
      break;
    }
    case Svev: {
      pars->outputs(3) = x;
      if (PRINTOUT > 1) msg << "Svev = " << x << ", ";
      break;
    }
    case Mtsq: {
      pars->outputs(4) = x;
      if (PRINTOUT > 1) msg << "MT = " << sqrt(x) << ", ";
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBParameter called with incorrect"
         << " independent parameter " << independent << '\n';
      throw ii.str();
    }
    }

    int error = 0;
    fixEWSBOutputs(pars, error);

    if (error != 0) {
      pars->failedVevIteration = true;
      if (PRINTOUT > 0) {
        cout << "# Warning: could not solve for VEVs\n";
      }
    }

    double output;
    switch (dependent) {
    case Lambda: {
      output = tempModel->displayLambda();
      if (PRINTOUT > 1) msg << "lambda = " << output << '\n';
      break;
    }
    case Kappa: {
      output = tempModel->displayKappa();
      if (PRINTOUT > 1) msg << "kappa = " << output << '\n';
      break;
    }
    case SMu: {
      output = tempModel->displaySusyMu();
      if (PRINTOUT > 1) msg << "mu = " << output << '\n';
      break;
    }
    case M3Sq: {
      output = tempModel->displayM3Squared();
      if (PRINTOUT > 1) msg << "m3sq = " << output << '\n';
      break;
    }
    case XiS: {
      output = tempModel->displayXiS();
      if (PRINTOUT > 1) msg << "xiS = " << output << '\n';
      break;
    }
    case Mh1Sq: {
      output = tempModel->displayMh1Squared();
      if (PRINTOUT > 1) msg << "mH1Sq = " << output << '\n';
      break;
    }
    case Mh2Sq: {
      output = tempModel->displayMh2Squared();
      if (PRINTOUT > 1) msg << "mH2Sq = " << output << '\n';
      break;
    }
    case MsSq: {
      output = tempModel->displayMsSquared();
      if (PRINTOUT > 1) msg << "mSsq = " << output << '\n';
      break;
    }
    case Yt: {
      output = tempModel->displayYukawaElement(YU, 3, 3);
      if (PRINTOUT > 1) msg << "ht = " << output << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBParameter called with incorrect"
         << " dependent parameter " << dependent << '\n';
      throw ii.str();
    }
    }

    if (PRINTOUT > 1) cout << msg.str();

    tempModel->setMu(startScale);
    tempModel->set(savedPars);
    tempModel->setDrBarPars(savedDrBarPars);
    tempModel->setProblem(savedProblems);
    pars->outputs = savedOutputs;

    return output;
  }

  std::pair<double,double> NmssmJacobian::calcEWSBOutputDerivative(
    NmssmSoftsusy& model, Parameters dep, Parameters indep) {

    double x = 0.;
    double h = 0.01;

    switch (indep) {
    case Lambda: {
      x = model.displayLambda();
      h = 0.0005 * x;
      break;
    }
    case Kappa: {
      x = model.displayKappa();
      h = 0.0005 * x;
      break;
    }
    case SMu: {
      x = model.displaySusyMu();
      h = 0.01 * x;
      break;
    }
    case M3Sq: {
      x = model.displayM3Squared();
      h = 0.01 * x;
      break;
    }
    case XiS: {
      x = model.displayXiS();
      h = 0.01 * x;
      break;
    }
    case Mh1Sq: {
      x = model.displayMh1Squared();
      h = 0.01 * x;
      break;
    }
    case Mh2Sq: {
      x = model.displayMh2Squared();
      h = 0.01 * x;
      break;
    }
    case MsSq: {
      x = model.displayMsSquared();
      h = 0.01 * x;
      break;
    }
    case Yt: {
      x = model.displayYukawaElement(YU, 3, 3);
      h = 0.0005 * x;
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBOutputDerivative called with incorrect"
         << " independent parameter " << indep << '\n';
      throw ii.str();
    }
    }

#ifdef ENABLE_GSL
    h = sqrt(std::numeric_limits<double>::epsilon()) * maximum(fabs(x), 1.0);
#endif

    volatile const double temp = x + h;
    h = temp - x;

    EWSBPars pars;
    pars.model = &model;
    pars.independent = indep;
    pars.dependent = dep;
    pars.useRunningMasses = useRunningMasses;
    pars.useSugraTrilinears = useSugraTrilinears;

    double derivative = 0.;
    double err = 0.;

    bool has_error = false;

#ifdef ENABLE_GSL
    gsl_function func;
    func.function = &calcEWSBOutput;
    func.params = &pars;

    gsl_deriv_central(&func, x, h, &derivative, &err);

    if (fabs(x) > 1.0e-10 &&
        (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
      has_error = fabs(err / derivative) > 1.0;
    }

    if (has_error) {
      // attempt to use a forward or backward difference in
      // case we are at a boundary at parameter space
      gsl_deriv_forward(&func, x, h, &derivative, &err);

      if (fabs(x) > 1.0e-10 &&
          (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
        has_error = fabs(err / derivative) > 1.0;
      }

      if (has_error) {
        gsl_deriv_backward(&func, x, h, &derivative, &err);

        if (fabs(x) > 1.0e-10 &&
            (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
          has_error = fabs(err / derivative) > 1.0;
        }
      }
    }
#else
    if (fabs(x) > 1.0e-10) {
      derivative = calcDerivative(calcEWSBOutput, x, h, &err, &pars);

      if (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10) {
        has_error = fabs(err / derivative) > 1.0;
      }
    }
#endif

    if (has_error) {
      derivative = -numberOfTheBeast;
    }

    failedVevIterationWarning = pars.failedVevIteration;
    inaccurateHiggsMassWarning = pars.inaccurateHiggsMass;
    if (fabs(pars.higgsMassAccuracy) > fabs(worstCaseHiggsAccuracy))
      worstCaseHiggsAccuracy = pars.higgsMassAccuracy;

    if (PRINTOUT > 1) {
      ostringstream msg;
      msg << "# derivative[d";

      switch (dep) {
      case Mzsq: msg << "MZsq/d"; break;
      case Tanb: msg << "tanb/d"; break;
      case Svev: msg << "Svev/d"; break;
      case Mtsq: msg << "MTsq/d"; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcEWSBOutputDerivative called with incorrect"
           << " dependent parameter " << dep << '\n';
        throw ii.str();
      }
      }

      switch (indep) {
      case Lambda: msg << "lambda] = "; break;
      case Kappa: msg << "kappa] = "; break;
      case SMu: msg << "mu] = "; break;
      case M3Sq: msg << "m3sq] = "; break;
      case XiS: msg << "xiS] = "; break;
      case Mh1Sq: msg << "mH1Sq] = "; break;
      case Mh2Sq: msg << "mH2Sq] = "; break;
      case MsSq: msg << "mSsq] = "; break;
      case Yt: msg << "ht] = "; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcEWSBOutputDerivative called with incorrect"
           << " independent parameter " << indep << '\n';
        throw ii.str();
      }
      }
      msg << derivative << ", error = " << err << '\n';
      cout << msg.str();
    }

    return std::pair<double,double>(derivative, err);
  }

  std::pair<double,double> NmssmJacobian::calcEWSBOutputLogDerivative(
    NmssmSoftsusy& model, Parameters dep, Parameters indep) {

    double x = 0.;

    switch (indep) {
    case Lambda: {
      x = std::log(model.displayLambda());
      break;
    }
    case Kappa: {
      x = std::log(model.displayKappa());
      break;
    }
    case SMu: {
      x = std::log(model.displaySusyMu());
      break;
    }
    case M3Sq: {
      x = std::log(model.displayM3Squared());
      break;
    }
    case XiS: {
      x = std::log(model.displayXiS());
      break;
    }
    case Mh1Sq: {
      x = std::log(model.displayMh1Squared());
      break;
    }
    case Mh2Sq: {
      x = std::log(model.displayMh2Squared());
      break;
    }
    case MsSq: {
      x = std::log(model.displayMsSquared());
      break;
    }
    case Yt: {
      x = std::log(model.displayYukawaElement(YU, 3, 3));
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBOutputLogDerivative called with incorrect"
         << " independent parameter " << indep << '\n';
      throw ii.str();
    }
    }

    double h = sqrt(std::numeric_limits<double>::epsilon()) * maximum(fabs(x), 1.0);

    volatile const double temp = x + h;
    h = temp - x;

    EWSBPars pars;
    pars.model = &model;
    pars.independent = indep;
    pars.dependent = dep;
    pars.useRunningMasses = useRunningMasses;
    pars.useSugraTrilinears = useSugraTrilinears;

    double derivative = 0.;
    double err = 0.;

    bool has_error = false;

#ifdef ENABLE_GSL
    gsl_function func;
    func.function = &calcLogEWSBOutput;
    func.params = &pars;

    gsl_deriv_central(&func, x, h, &derivative, &err);

    if (fabs(x) > 1.0e-10 &&
        (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
      has_error = fabs(err / derivative) > 1.0;
    }

    if (has_error) {
      // attempt to use a forward or backward difference in
      // case we are at a boundary at parameter space
      gsl_deriv_forward(&func, x, h, &derivative, &err);

      if (fabs(x) > 1.0e-10 &&
          (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
        has_error = fabs(err / derivative) > 1.0;
      }

      if (has_error) {
        gsl_deriv_backward(&func, x, h, &derivative, &err);

        if (fabs(x) > 1.0e-10 &&
            (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
          has_error = fabs(err / derivative) > 1.0;
        }
      }
    }
#else
    if (fabs(x) > 1.0e-10) {
      derivative = calcDerivative(calcLogEWSBOutput, x, h, &err, &pars);

      if (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10) {
        has_error = fabs(err / derivative) > 1.0;
      }
    }
#endif

    if (has_error) {
      derivative = -numberOfTheBeast;
    }

    failedVevIterationWarning = pars.failedVevIteration;
    inaccurateHiggsMassWarning = pars.inaccurateHiggsMass;
    if (fabs(pars.higgsMassAccuracy) > fabs(worstCaseHiggsAccuracy))
      worstCaseHiggsAccuracy = pars.higgsMassAccuracy;

    if (PRINTOUT > 1) {
      ostringstream msg;
      msg << "# derivative[d";

      switch (dep) {
      case Mzsq: msg << "log(MZsq)/d"; break;
      case Tanb: msg << "log(tanb)/d"; break;
      case Svev: msg << "log(Svev)/d"; break;
      case Mtsq: msg << "log(MTsq)/d"; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcEWSBOutputLogDerivative called with incorrect"
           << " dependent parameter " << dep << '\n';
        throw ii.str();
      }
      }

      switch (indep) {
      case Lambda: msg << "log(lambda)] = "; break;
      case Kappa: msg << "log(kappa)] = "; break;
      case SMu: msg << "log(mu)] = "; break;
      case M3Sq: msg << "log(m3sq)] = "; break;
      case XiS: msg << "log(xiS)] = "; break;
      case Mh1Sq: msg << "log(mH1Sq)] = "; break;
      case Mh2Sq: msg << "log(mH2Sq)] = "; break;
      case MsSq: msg << "log(mSsq)] = "; break;
      case Yt: msg << "log(ht)] = "; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcEWSBOutputLogDerivative called with incorrect"
           << " independent parameter " << indep << '\n';
        throw ii.str();
      }
      }
      msg << derivative << ", error = " << err << '\n';
      cout << msg.str();
    }

    return std::pair<double,double>(derivative, err);
  }

  std::pair<double,double> NmssmJacobian::calcEWSBParameterDerivative(
    NmssmSoftsusy& model, Parameters dep, Parameters indep) {

    const int numOutputs = includeTop ? 4 : 3;
    DoubleVector outputs(numOutputs);
    outputs(1) = sqr(calcMz(model, useRunningMasses));
    outputs(2) = model.displayTanb();
    outputs(3) = (model.displayZ3() && !SoftHiggsOut) ? model.displayLambda() :
      model.displaySvev();
    if (numOutputs > 3) outputs(4) = sqr(calcMt(model, useRunningMasses));

    double x = 0.;
    double h = 0.01;

    switch (indep) {
    case Lambda: {
      x = outputs(3);
      h = 0.0005 * x;
      break;
    }
    case Mzsq: {
      x = outputs(1);
      h = 0.01 * x;
      break;
    }
    case Tanb: {
      x = outputs(2);
      h = 0.001 * x;
      break;
    }
    case Svev: {
      x = outputs(3);
      h = 0.01 * x;
      break;
    }
    case Mtsq: {
      x = outputs(4);
      h = 0.01 * x;
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBParameterDerivative called with incorrect"
         << " independent parameter " << indep << '\n';
      throw ii.str();
    }
    }

#ifdef ENABLE_GSL
    h = sqrt(std::numeric_limits<double>::epsilon()) * maximum(fabs(x), 1.0);
#endif

    volatile const double temp = x + h;
    h = temp - x;

    EWSBPars pars;
    pars.model = &model;
    pars.independent = indep;
    pars.dependent = dep;
    pars.outputs = outputs;
    pars.useRunningMasses = useRunningMasses;
    pars.useSugraTrilinears = useSugraTrilinears;

    double derivative = 0.;
    double err = 0.;

    bool has_error = false;

#ifdef ENABLE_GSL
    gsl_function func;
    func.function = &calcEWSBParameter;
    func.params = &pars;

    gsl_deriv_central(&func, x, h, &derivative, &err);

    if (fabs(x) > 1.0e-10 &&
        (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
      has_error = fabs(err / derivative) > 1.0;
    }

    if (has_error) {
      // attempt to use a forward or backward difference in
      // case we are at a boundary at parameter space
      gsl_deriv_forward(&func, x, h, &derivative, &err);

      if (fabs(x) > 1.0e-10 &&
          (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
        has_error = fabs(err / derivative) > 1.0;
      }

      if (has_error) {
        gsl_deriv_backward(&func, x, h, &derivative, &err);

        if (fabs(x) > 1.0e-10 &&
            (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
          has_error = fabs(err / derivative) > 1.0;
        }
      }
    }
#else
    if (fabs(x) > 1.0e-10) {
      derivative = calcDerivative(calcEWSBParameter, x, h, &err, &pars);

      if (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10) {
        has_error = fabs(err / derivative) > 1.0;
      }
    }
#endif

    if (has_error) {
      derivative = -numberOfTheBeast;
    }

    failedVevIterationWarning = pars.failedVevIteration;
    inaccurateHiggsMassWarning = pars.inaccurateHiggsMass;
    if (fabs(pars.higgsMassAccuracy) > fabs(worstCaseHiggsAccuracy))
      worstCaseHiggsAccuracy = pars.higgsMassAccuracy;

    if (PRINTOUT > 1) {
      ostringstream msg;
      msg << "# derivative[d";

      switch (dep) {
      case Lambda: msg << "lambda/d"; break;
      case Kappa: msg << "kappa/d"; break;
      case SMu: msg << "mu/d"; break;
      case M3Sq: msg << "m3sq/d"; break;
      case XiS: msg << "xiS/d"; break;
      case Mh1Sq: msg << "mH1Sq/d"; break;
      case Mh2Sq: msg << "mH2Sq/d"; break;
      case MsSq: msg << "mSsq/d"; break;
      case Yt: msg << "ht/d"; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcEWSBParameterDerivative called with incorrect"
           << " dependent parameter " << dep << '\n';
        throw ii.str();
      }
      }

      switch (indep) {
      case Lambda: msg << "lambda] = "; break;
      case Mzsq: msg << "MZsq] = "; break;
      case Tanb: msg << "tanb] = "; break;
      case Svev: msg << "Svev] = "; break;
      case Mtsq: msg << "MTsq] = "; break;
      default: {
        ostringstream ii;
        ii << "NmssmJacobian:calcEWSBParameterDerivative called with incorrect"
           << " independent parameter " << indep << '\n';
        throw ii.str();
      }
      }

      msg << derivative << ", error = " << err << '\n';
      cout << msg.str();
    }

    return std::pair<double,double>(derivative, err);
  }

  double NmssmJacobian::calcEWSBJacobian(NmssmSoftsusy& model) {

    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();
    const drBarPars savedDrBarPars(model.displayDrBarPars());

    const int numPars = includeTop ? 4 : 3;

    // reset warnings on each call
    failedVevIterationWarning = false;
    inaccurateHiggsMassWarning = false;
    worstCaseHiggsAccuracy = 0.;

    DoubleMatrix jac(numPars, numPars);
    DoubleMatrix jacErrors(numPars, numPars);

    std::vector<Parameters> indepPars;
    std::vector<double> paramValues;

    if (model.displayZ3() && !SoftHiggsOut) {
      indepPars.push_back(Mzsq);
      indepPars.push_back(Tanb);
      indepPars.push_back(Lambda);

      paramValues.push_back(sqr(calcMz(model, useRunningMasses)));
      paramValues.push_back(model.displayTanb());
      paramValues.push_back(model.displayLambda());
    } else {
      indepPars.push_back(Mzsq);
      indepPars.push_back(Tanb);
      indepPars.push_back(Svev);

      paramValues.push_back(sqr(calcMz(model, useRunningMasses)));
      paramValues.push_back(model.displayTanb());
      paramValues.push_back(model.displaySvev());
    }
    if (includeTop) {
      indepPars.push_back(Mtsq);
      paramValues.push_back(sqr(calcMt(model, useRunningMasses)));
    }

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      std::pair<double,double> result;
      if (SoftHiggsOut) {
        result = calcEWSBParameterDerivative(model, Mh1Sq, indepPars[i]);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcEWSBParameterDerivative(model, Mh2Sq, indepPars[i]);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcEWSBParameterDerivative(model, MsSq, indepPars[i]);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      } else if ( model.displayZ3() ) {
        if (indepPars[i] == Lambda) {
          jac(i + 1, 1) = 1.;
          jacErrors(i + 1, 1) = 0.;
          jac(i + 1, 2) = 0.;
          jacErrors(i + 1, 2) = 0.;
          jac(i + 1, 3) = 0.;
          jacErrors(i + 1, 3) = 0.;
        } else {
          jac(i + 1, 1) = 0.;
          jacErrors(i + 1, 1) = 0.;

          result = calcEWSBParameterDerivative(model, Kappa, indepPars[i]);
          jac(i + 1, 2) = result.first;
          jacErrors(i + 1, 2) = result.second;

          result = calcEWSBParameterDerivative(model, MsSq, indepPars[i]);
          jac(i + 1, 3) = result.first;
          jacErrors(i + 1, 3) = result.second;
        }
      } else {
        result = calcEWSBParameterDerivative(model, SMu, indepPars[i]);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcEWSBParameterDerivative(model, M3Sq, indepPars[i]);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcEWSBParameterDerivative(model, XiS, indepPars[i]);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      }
      if (includeTop) {
        if ((model.displayZ3() && !SoftHiggsOut) && indepPars[i] == Lambda) {
          jac(i + 1, 4) = 0.;
          jacErrors(i + 1, 4) = 0.;
        } else {
          result = calcEWSBParameterDerivative(model, Yt, indepPars[i]);
          jac(i + 1, 4) = result.first;
          jacErrors(i + 1, 4) = result.second;
        }
      }
    }

    // check for errors
    hasError = checkDerivativeErrors(jac, jacErrors, paramValues);

    // save calculated matrix
    if (jacEWSB.displayRows() != numPars
        || jacEWSB.displayCols() != numPars) {
      jacEWSB.resize(numPars, numPars);
    }
    jacEWSB = jac;

    if (jacEWSBErrors.displayRows() != numPars
        || jacEWSBErrors.displayCols() != numPars) {
      jacEWSBErrors.resize(numPars, numPars);
    }
    jacEWSBErrors = jacErrors;

    model.setMu(scale);
    model.set(savedPars);
    model.setDrBarPars(savedDrBarPars);

    return jac.determinant();
  }

  double NmssmJacobian::calcInverseEWSBJacobian(NmssmSoftsusy& model) {

    const int numPars = includeTop ? 4 : 3;

    // reset warnings on each call
    failedVevIterationWarning = false;
    inaccurateHiggsMassWarning = false;
    worstCaseHiggsAccuracy = 0.;

    DoubleMatrix jac(numPars, numPars);
    DoubleMatrix jacErrors(numPars, numPars);

    std::vector<Parameters> indepPars;
    std::vector<double> paramValues;

    if (SoftHiggsOut) {
      indepPars.push_back(Mh1Sq);
      indepPars.push_back(Mh2Sq);
      indepPars.push_back(MsSq);

      paramValues.push_back(model.displayMh1Squared());
      paramValues.push_back(model.displayMh2Squared());
      paramValues.push_back(model.displayMsSquared());
    } else if ( model.displayZ3() ) {
      indepPars.push_back(Lambda);
      indepPars.push_back(Kappa);
      indepPars.push_back(MsSq);

      paramValues.push_back(model.displayLambda());
      paramValues.push_back(model.displayKappa());
      paramValues.push_back(model.displayMsSquared());
    } else {
      indepPars.push_back(SMu);
      indepPars.push_back(M3Sq);
      indepPars.push_back(XiS);

      paramValues.push_back(model.displaySusyMu());
      paramValues.push_back(model.displayM3Squared());
      paramValues.push_back(model.displayXiS());
    }
    if (includeTop) {
      indepPars.push_back(Yt);
      paramValues.push_back(model.displayYukawaElement(YU, 3, 3));
    }

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      std::pair<double,double> result;
      if (model.displayZ3() && !SoftHiggsOut) {
        if (indepPars[i] == Lambda) {
          jac(i + 1, 1) = 0.;
          jacErrors(i + 1, 1) = 0.;
          jac(i + 1, 2) = 0.;
          jacErrors(i + 1, 2) = 0.;
          jac(i + 1, 3) = 1.;
          jacErrors(i + 1, 3) = 0.;
        } else {
          result = calcEWSBOutputDerivative(model, Mzsq, indepPars[i]);
          jac(i + 1, 1) = result.first;
          jacErrors(i + 1, 1) = result.second;

          result = calcEWSBOutputDerivative(model, Tanb, indepPars[i]);
          jac(i + 1, 2) = result.first;
          jacErrors(i + 1, 2) = result.second;

          jac(i + 1, 3) = 0.;
          jacErrors(i + 1, 3) = 0.;
        }
      } else {
        result = calcEWSBOutputDerivative(model, Mzsq, indepPars[i]);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcEWSBOutputDerivative(model, Tanb, indepPars[i]);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcEWSBOutputDerivative(model, Svev, indepPars[i]);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      }
      if (includeTop) {
        if ((model.displayZ3() && !SoftHiggsOut) && indepPars[i] == Lambda) {
          jac(i + 1, 4) = 0.;
          jacErrors(i + 1, 4) = 0.;
        } else {
          result = calcEWSBOutputDerivative(model, Mtsq, indepPars[i]);
          jac(i + 1, 4) = result.first;
          jacErrors(i + 1, 4) = result.second;
        }
      }
    }

    // check for errors
    hasError = checkDerivativeErrors(jac, jacErrors, paramValues);

    // save calculated matrix
    if (invJacEWSB.displayRows() != numPars
        || invJacEWSB.displayCols() != numPars) {
      invJacEWSB.resize(numPars, numPars);
    }
    invJacEWSB = jac;

    if (invJacEWSBErrors.displayRows() != numPars
        || invJacEWSBErrors.displayCols() != numPars) {
      invJacEWSBErrors.resize(numPars, numPars);
    }
    invJacEWSBErrors = jacErrors;

    return jac.determinant();
  }

  // @todo handle negative parameters properly
  double NmssmJacobian::calcInverseEWSBJacobianLogs(NmssmSoftsusy& model) {

    const int numPars = includeTop ? 4 : 3;

    // reset warnings on each call
    failedVevIterationWarning = false;
    inaccurateHiggsMassWarning = false;
    worstCaseHiggsAccuracy = 0.;

    DoubleMatrix jac(numPars, numPars);
    DoubleMatrix jacErrors(numPars, numPars);

    std::vector<Parameters> indepPars;
    std::vector<double> paramValues;

    if (SoftHiggsOut) {
      indepPars.push_back(Mh1Sq);
      indepPars.push_back(Mh2Sq);
      indepPars.push_back(MsSq);

      paramValues.push_back(std::log(model.displayMh1Squared()));
      paramValues.push_back(std::log(model.displayMh2Squared()));
      paramValues.push_back(std::log(model.displayMsSquared()));
    } else if ( model.displayZ3() ) {
      indepPars.push_back(Lambda);
      indepPars.push_back(Kappa);
      indepPars.push_back(MsSq);

      paramValues.push_back(std::log(model.displayLambda()));
      paramValues.push_back(std::log(model.displayKappa()));
      paramValues.push_back(std::log(model.displayMsSquared()));
    } else {
      indepPars.push_back(SMu);
      indepPars.push_back(M3Sq);
      indepPars.push_back(XiS);

      paramValues.push_back(std::log(model.displaySusyMu()));
      paramValues.push_back(std::log(model.displayM3Squared()));
      paramValues.push_back(std::log(model.displayXiS()));
    }
    if (includeTop) {
      indepPars.push_back(Yt);
      paramValues.push_back(std::log(model.displayYukawaElement(YU, 3, 3)));
    }

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      std::pair<double,double> result;
      if (model.displayZ3() && !SoftHiggsOut) {
        if (indepPars[i] == Lambda) {
          jac(i + 1, 1) = 0.;
          jacErrors(i + 1, 1) = 0.;
          jac(i + 1, 2) = 0.;
          jacErrors(i + 1, 2) = 0.;
          jac(i + 1, 3) = 1.;
          jacErrors(i + 1, 3) = 0.;
        } else {
          result = calcEWSBOutputLogDerivative(model, Mzsq, indepPars[i]);
          jac(i + 1, 1) = result.first;
          jacErrors(i + 1, 1) = result.second;

          result = calcEWSBOutputLogDerivative(model, Tanb, indepPars[i]);
          jac(i + 1, 2) = result.first;
          jacErrors(i + 1, 2) = result.second;

          jac(i + 1, 3) = 0.;
          jacErrors(i + 1, 3) = 0.;
        }
      } else {
        result = calcEWSBOutputLogDerivative(model, Mzsq, indepPars[i]);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcEWSBOutputLogDerivative(model, Tanb, indepPars[i]);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;

        result = calcEWSBOutputLogDerivative(model, Svev, indepPars[i]);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      }
      if (includeTop) {
        if ((model.displayZ3() && !SoftHiggsOut) && indepPars[i] == Lambda) {
          jac(i + 1, 4) = 0.;
          jacErrors(i + 1, 4) = 0.;
        } else {
          result = calcEWSBOutputLogDerivative(model, Mtsq, indepPars[i]);
          jac(i + 1, 4) = result.first;
          jacErrors(i + 1, 4) = result.second;
        }
      }
    }

    // check for errors
    hasError = checkDerivativeErrors(jac, jacErrors, paramValues);

    // save calculated matrix
    if (invJacEWSB.displayRows() != numPars
        || invJacEWSB.displayCols() != numPars) {
      invJacEWSB.resize(numPars, numPars);
    }
    invJacEWSB = jac;

    if (invJacEWSBErrors.displayRows() != numPars
        || invJacEWSBErrors.displayCols() != numPars) {
      invJacEWSBErrors.resize(numPars, numPars);
    }
    invJacEWSBErrors = jacErrors;

    return jac.determinant();
  }

  bool NmssmJacobian::checkDerivativeErrors(
     DoubleMatrix& derivs, const DoubleMatrix& errors,
     const std::vector<double>& paramValues) const {
    bool error = false;
    const int numDependent = derivs.displayCols();
    for (int i = 0, numPars = paramValues.size(); i < numPars; ++i) {
      if (fabs(paramValues[i]) <= 1.0e-10)
        continue;
      for (int j = 0; j < numDependent; ++j) {
        if (fabs(derivs(i + 1, j + 1)) > 1.0e-10
            || fabs(errors(i + 1, j + 1)) > 1.0e-10) {
          if (fabs(errors(i + 1, j + 1) / derivs(i + 1, j + 1)) > 1.0) {
            error = true;
            derivs(i + 1, j + 1) = -numberOfTheBeast;
          }
        }
      }
    }
    return error;
  }

} /// namespace softsusy
