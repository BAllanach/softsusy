/** \file mssmjacobian.cpp

    \brief Implementation of routines for calculating Jacobian fine-tuning.
*/

#include "config.h"
#include "def.h"
#include "mssmjacobian.h"
#include "softsusy.h"
#include "numerics.h"
#include "utils.h"

#include <iostream>
#include <limits>
#include <vector>

#ifdef ENABLE_GSL
#include <gsl/gsl_deriv.h>
#endif

namespace softsusy {

  MssmJacobian::MssmJacobian(bool doTop)
     : invJacRGFlow(3,3), invJacEWSB(3,3)
     , invJacRGFlowErrors(3,3), invJacEWSBErrors(3,3)
     , includeTop(doTop), useRunningMasses(false)
     , useSugraTrilinears(false), hasError(false) {
  }

  MssmJacobian::~MssmJacobian() {}

  /// Calculates the fine-tuning for the model, using the Jacobian
  /// measure.  The parameters are taken to be defined at the input
  /// scale \c mx obtained from a call to \c m.displayMxBC(),
  /// while the observables are calculated at the SUSY scale
  /// returned by \c m.calcMs().  The Jacobian matrix is
  /// calculated by calling calcFTInverseJacobian() internally, and
  /// can be accessed after the calculation.
  /// \see calcFTInverseJacobian(const MssmSoftsusy&)
  double MssmJacobian::calcDeltaJ(const MssmSoftsusy& m) {
    return calcDeltaJ(m, m.displayMxBC());
  }

  /// Calculates the fine-tuning as in calcDeltaJ(const MssmSoftsusy& m),
  /// but with \c mx set to the given value instead of that returned by
  /// \c m.displayMxBC().
  double MssmJacobian::calcDeltaJ(const MssmSoftsusy& m, double mx) {
    MssmSoftsusy model(m);

    model.runto(model.calcMs());

    model.calcDrBarPars();
    const double mtrun = model.displayDrBarPars().mt;
    const double sinthDRbar = model.calcSinthdrbar();
    model.doTadpoles(mtrun, sinthDRbar);

    const double determinant = calcFTInverseJacobian(model, mx);

    const double mz2 = sqr(calcMz(model, useRunningMasses));
    const double tb = model.displayTanb();

    double denominator = mz2 * tb;
    if (includeTop) {
      denominator *= sqr(calcMt(model, useRunningMasses));
    }

    model.runto(mx);

    double numerator;
    if (model.displayAltEwsb()) {
      numerator = model.displayMh1Squared() * model.displayMh2Squared();
    } else {
      numerator = model.displaySusyMu() * model.displayM3Squared();
    }
    if (includeTop) numerator *= model.displayYukawaElement(YU, 3, 3);

    double ans = fabs(numerator * determinant / denominator);
    if (displayError()) ans = asin(2.0);
    
    return ans;
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
  double MssmJacobian::calcFTInverseJacobian(const MssmSoftsusy& m) {
    return calcFTInverseJacobian(m, m.displayMxBC());
  }

  /// Calculates the Jacobian as in
  /// calcFTInverseJacobian(const MssmSoftsusy& m), but with \c mx set to
  /// the given value instead of that returned by \c m.displayMxBC().
  double MssmJacobian::calcFTInverseJacobian(const MssmSoftsusy& m, double mx) {
    MssmSoftsusy model(m);

    model.runto(model.calcMs());

    const double rgDet = calcRGFlowJacobian(model, mx, model.displayMu());
    const double ewsbDet = calcInverseEWSBJacobian(model);

    return ewsbDet * rgDet;
  }

  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double MssmJacobian::calcMz(MssmSoftsusy& model, bool getRunningMass) {

    double mz = model.displayMzRun();
    if (!getRunningMass) {
      const double scale = model.displayMu();
      const double pizzt = model.piZZT(mz, scale);
      mz = sqrt(fabs(sqr(mz) - pizzt));
    }

    return mz;
  }

  // @todo make sure this is correct for MSSM
  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double MssmJacobian::calcMt(MssmSoftsusy& model, bool getRunningMass) {

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

  double MssmJacobian::calcRunningParameter(double x, void* parameters) {
    RGFlowPars* pars = static_cast<RGFlowPars*>(parameters);

    MssmSoftsusy* tempModel = pars->model;
    Parameters independent = pars->independent;
    Parameters dependent = pars->dependent;
    const double endScale = pars->toScale;

    const DoubleVector savedPars(tempModel->display());
    const double startScale = tempModel->displayMu();

    ostringstream msg;
    if (PRINTOUT > 1) msg << "# Qstart = " << startScale << ", ";

    switch (independent) {
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
      ii << "MssmJacobian:calcRunningParameter called with incorrect"
         << " independent parameter " << independent << '\n';
      throw ii.str();
    }
    }

    tempModel->runto(endScale);

    if (PRINTOUT > 1) msg << "Qend = " << endScale << ", ";

    double output;
    switch (dependent) {
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
    case Yt: {
      output = tempModel->displayYukawaElement(YU, 3, 3);
      if (PRINTOUT > 1) msg << "ht = " << output << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "MssmJacobian:calcRunningParameter called with incorrect"
         << " dependent parameter " << dependent << '\n';
      throw ii.str();
    }
    }

    if (PRINTOUT > 1) cout << msg.str();

    tempModel->setMu(startScale);
    tempModel->set(savedPars);

    return output;
  }

  std::pair<double,double> MssmJacobian::calcRGDerivative(
    MssmSoftsusy& model, Parameters dep, Parameters indep, double toScale) {

    double x = 0.;
    double h = 0.01;

    switch (indep) {
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
    case Yt: {
      x = model.displayYukawaElement(YU, 3, 3);
      h = 0.0005 * x;
      break;
    }
    default: {
      ostringstream ii;
      ii << "MssmJacobian:calcRGDerivative called with incorrect"
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
      has_error = false;
      gsl_deriv_forward(&func, x, h, &derivative, &err);

      if (fabs(x) > 1.0e-10 &&
          (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
        has_error = fabs(err / derivative) > 1.0;
      }

      if (has_error) {
        has_error = false;
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
      hasError = true;
    }

    if (PRINTOUT > 1) {
      const double scale = model.displayMu();

      ostringstream msg;
      msg << "# derivative[d";

      switch (dep) {
      case SMu: msg << "mu(Q=" << toScale << ")/d"; break;
      case M3Sq: msg << "m3sq(Q=" << toScale << ")/d"; break;
      case Mh1Sq: msg << "mH1Sq(Q=" << toScale << ")/d"; break;
      case Mh2Sq: msg << "mH2Sq(Q=" << toScale << ")/d"; break;
      case Yt: msg << "ht(Q=" << toScale << ")/d"; break;
      default: {
        ostringstream ii;
        ii << "MssmJacobian:calcRGDerivative called with incorrect"
           << " dependent parameter " << dep << '\n';
        throw ii.str();
      }
      }

      switch (indep) {
      case SMu: msg << "mu(Q=" << scale << ")] = "; break;
      case M3Sq: msg << "m3sq(Q=" << scale << ")] = "; break;
      case Mh1Sq: msg << "mH1Sq(Q=" << scale << ")] = "; break;
      case Mh2Sq: msg << "mH2Sq(Q=" << scale << ")] = "; break;
      case Yt: msg << "ht(Q=" << scale << ")] = "; break;
      default: {
        ostringstream ii;
        ii << "MssmJacobian:calcRGDerivative called with incorrect"
           << " independent parameter " << indep << '\n';
        throw ii.str();
      }
      }

      msg << derivative << ", error = " << err << '\n';
      cout << msg.str();
    }

    return std::pair<double,double>(derivative, err);
  }

  double MssmJacobian::calcRGFlowJacobian(MssmSoftsusy& model,
                                           double startScale, double endScale) {
    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();

    model.runto(startScale);

    const int numPars = includeTop ? 3 : 2;

    DoubleMatrix jac(numPars, numPars);
    DoubleMatrix jacErrors(numPars, numPars);

    std::vector<Parameters> indepPars;

    if (model.displayAltEwsb()) {
      indepPars.push_back(Mh1Sq);
      indepPars.push_back(Mh2Sq);
    } else {
      indepPars.push_back(SMu);
      indepPars.push_back(M3Sq);
    }
    if (includeTop) indepPars.push_back(Yt);

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      std::pair<double,double> result;
      if (model.displayAltEwsb()) {
        result = calcRGDerivative(model, Mh1Sq, indepPars[i], endScale);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcRGDerivative(model, Mh2Sq, indepPars[i], endScale);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;
      } else {
        result = calcRGDerivative(model, SMu, indepPars[i], endScale);
        jac(i + 1, 1) = result.first;
        jacErrors(i + 1, 1) = result.second;

        result = calcRGDerivative(model, M3Sq, indepPars[i], endScale);
        jac(i + 1, 2) = result.first;
        jacErrors(i + 1, 2) = result.second;
      }
      if (includeTop) {
        result = calcRGDerivative(model, Yt, indepPars[i], endScale);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      }
    }

    // save calculated matrix
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

    model.setMu(scale);
    model.set(savedPars);

    return jac.determinant();
  }

  double MssmJacobian::calcEWSBOutput(double x, void* parameters) {
    EWSBPars* pars = static_cast<EWSBPars*>(parameters);

    MssmSoftsusy* tempModel = pars->model;
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
      ii << "MssmJacobian:calcEWSBOutput called with incorrect"
         << " independent parameter " << independent << '\n';
      throw ii.str();
    }
    }

    tempModel->calcDrBarPars();
    tempModel->doTadpoles(tempModel->displayDrBarPars().mt,
                          tempModel->calcSinthdrbar());

    DoubleVector vevs(2);
    vevs(1) = tempModel->displayHvev();
    vevs(2) = tempModel->displayTanb();

    int error = 0;
    tempModel->predVevs(vevs, error);

    if (error != 0) {
      if (PRINTOUT > 0) {
        cout << "# Warning: could not solve for VEVs\n";
      }
    }

    tempModel->setHvev(vevs(1));
    tempModel->setTanb(vevs(2));

    tempModel->calcDrBarPars();
    const double mtrun = tempModel->displayDrBarPars().mt;
    const double sinthDRbar = tempModel->calcSinthdrbar();
    tempModel->doTadpoles(mtrun, sinthDRbar);

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
    case Mtsq: {
      output = sqr(calcMt(*tempModel, pars->useRunningMasses));
      if (PRINTOUT > 1) msg << "MT = " << sqrt(output) << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "MssmJacobian:calcEWSBOutput called with incorrect"
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

  std::pair<double,double> MssmJacobian::calcEWSBOutputDerivative(
    MssmSoftsusy& model, Parameters dep, Parameters indep) {

    double x = 0.;
    double h = 0.01;

    switch (indep) {
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
    case Yt: {
      x = model.displayYukawaElement(YU, 3, 3);
      h = 0.0005 * x;
      break;
    }
    default: {
      ostringstream ii;
      ii << "MssmJacobian:calcEWSBOutputDerivative called with incorrect"
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
      has_error = false;
      gsl_deriv_forward(&func, x, h, &derivative, &err);

      if (fabs(x) > 1.0e-10 &&
          (fabs(derivative) > 1.0e-10 || fabs(err) > 1.0e-10)) {
        has_error = fabs(err / derivative) > 1.0;
      }

      if (has_error) {
        has_error = false;
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
      hasError = true;
    }

    if (PRINTOUT > 1) {
      ostringstream msg;
      msg << "# derivative[d";

      switch (dep) {
      case Mzsq: msg << "MZsq/d"; break;
      case Tanb: msg << "tanb/d"; break;
      case Mtsq: msg << "MTsq/d"; break;
      default: {
        ostringstream ii;
        ii << "MssmJacobian:calcEWSBOutputDerivative called with incorrect"
           << " dependent parameter " << dep << '\n';
        throw ii.str();
      }
      }

      switch (indep) {
      case SMu: msg << "mu] = "; break;
      case M3Sq: msg << "m3sq] = "; break;
      case Mh1Sq: msg << "mH1Sq] = "; break;
      case Mh2Sq: msg << "mH2Sq] = "; break;
      case Yt: msg << "ht] = "; break;
      default: {
        ostringstream ii;
        ii << "MssmJacobian:calcEWSBOutputDerivative called with incorrect"
           << " independent parameter " << indep << '\n';
        throw ii.str();
      }
      }
      msg << derivative << ", error = " << err << '\n';
      cout << msg.str();
    }

    return std::pair<double,double>(derivative, err);
  }

  double MssmJacobian::calcInverseEWSBJacobian(MssmSoftsusy& model) {

    const int numPars = includeTop ? 3 : 2;

    DoubleMatrix jac(numPars, numPars);
    DoubleMatrix jacErrors(numPars, numPars);

    std::vector<Parameters> indepPars;

    if (model.displayAltEwsb()) {
      indepPars.push_back(Mh1Sq);
      indepPars.push_back(Mh2Sq);
    } else {
      indepPars.push_back(SMu);
      indepPars.push_back(M3Sq);
    }
    if (includeTop) indepPars.push_back(Yt);

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      std::pair<double,double> result;
      result = calcEWSBOutputDerivative(model, Mzsq, indepPars[i]);
      jac(i + 1, 1) = result.first;
      jacErrors(i + 1, 1) = result.second;

      result = calcEWSBOutputDerivative(model, Tanb, indepPars[i]);
      jac(i + 1, 2) = result.first;
      jacErrors(i + 1, 2) = result.second;

      if (includeTop) {
        result = calcEWSBOutputDerivative(model, Mtsq, indepPars[i]);
        jac(i + 1, 3) = result.first;
        jacErrors(i + 1, 3) = result.second;
      }
    }

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

} /// namespace softsusy
