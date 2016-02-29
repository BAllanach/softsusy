/** \file nmssmjacobian.cpp

    \brief Implementation of routines for calculating Jacobian fine-tuning.
 */

#include "def.h"
#include "nmssmjacobian.h"
#include "nmssmsoftsusy.h"
#include "numerics.h"
#include "utils.h"

#include <vector>

namespace softsusy {

  NmssmJacobian::NmssmJacobian(bool doTop)
     : jacRGFlow(3,3), jacEWSB(3,3)
     , invJacRGFlow(3,3), invJacEWSB(3,3)
     , includeTop(doTop), useRunningMasses(false), hasError(false) {
  }

  NmssmJacobian::~NmssmJacobian() {}

  /// Calculates the fine-tuning for the model, using the Jacobian
  /// measure.  The parameters are taken to be defined at the input
  /// scale \c mx obtained from a call to \c model.displayMxBC(),
  /// while the observables are calculated at the SUSY scale
  /// returned by \c model.calcMs().  The Jacobian matrix is
  /// calculated by calling calcFTInverseJacobian() internally, and
  /// can be accessed after the calculation.
  /// \see calcFTInverseJacobian(NmssmSoftsusy&)
  double NmssmJacobian::calcDeltaJ(NmssmSoftsusy& model) {
    return calcDeltaJ(model, model.displayMxBC());
  }

  /// Calculates the fine-tuning as in calcDeltaJ(NmssmSoftsusy& model),
  /// but with \c mx set to the given value instead of that returned by
  /// \c model.displayMxBC().
  double NmssmJacobian::calcDeltaJ(NmssmSoftsusy& model, double mx) {

    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();
    const sPhysical savedPhys(model.displayPhys());

    model.runto(model.calcMs());

    const double determinant = calcFTInverseJacobian(model, mx);

    const double mz2 = sqr(calcMz(model, useRunningMasses));
    const double tb = model.displayTanb();

    double denominator = mz2 * tb;
    if (Z3) {
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
    } else if (Z3) {
      numerator = model.displayLambda() * model.displayKappa()
        * model.displayMsSquared();
    } else {
      numerator = model.displaySusyMu() * model.displayM3Squared()
        * model.displayXiS();
    }
    if (includeTop) numerator *= model.displayYukawaElement(YU, 3, 3);

    model.setMu(scale);
    model.set(savedPars);
    model.setPhys(savedPhys);

    return fabs(numerator * determinant / denominator);
  }

  /// Calculates the Jacobian of the form \f$ J^{-1} = |
  /// \partial O / \partial p | \f$.  The transformation is done
  /// in two stages.  In the first, the observables at the
  /// SUSY scale, obtained from \c model.calcMs(), are traded for
  /// the parameters at this scale using the EWSB conditions.
  /// The Jacobian matrix for this transformation can be
  /// accessed afterwards using displayInverseEWSBJacobian().
  /// Then, the parameters at this scale are transformed to
  /// parameters at the scale \c mx, obtained from calling
  /// \c model.displayMxBC(), using the RGEs.  The Jacobian matrix for
  /// this transformation may be accessed by calling
  /// displayInverseRGFlowJacobian().
  double NmssmJacobian::calcFTInverseJacobian(NmssmSoftsusy& model) {
    return calcFTInverseJacobian(model, model.displayMxBC());
  }

  /// Calculates the Jacobian as in calcFTInverseJacobian(NmssmSoftsusy& model),
  /// but with \c mx set to the given value instead of that returned by
  /// \c model.displayMxBC().
  double NmssmJacobian::calcFTInverseJacobian(NmssmSoftsusy& model, double mx) {

    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();
    const sPhysical savedPhys(model.displayPhys());

    model.runto(model.calcMs());

    const double rgDet = calcRGFlowJacobian(model, mx, model.displayMu());
    const double ewsbDet = calcInverseEWSBJacobian(model);

    model.setMu(scale);
    model.set(savedPars);
    model.setPhys(savedPhys);

    return ewsbDet * rgDet;
  }

  /// Calculates the Jacobian of the form \f$ J = |
  /// \partial p / \partial O | \f$.  The transformation,
  /// which is the inverse of that calculated by
  /// calcFTInverseJacobian(), is done in two stages.  In the
  /// first, the parameters at the input scale \c mx, obtained
  /// from a call to \c model.displayMxBC(), are transformed to parameters
  /// \f$\{q\}\f$ at the SUSY scale, obtained from \c model.calcMs(),
  /// using the RGEs.  The Jacobian matrix for this
  /// transformation may be accessed afterwards by calling
  /// displayRGFlowJacobian().  Then the parameters
  /// \f$\{q\}\f$ are traded for the observables at the same
  /// scale using the EWSB conditions.  The Jacobian matrix for
  /// this second transformation can be accessed using
  /// displayEWSBJacobian().
  double NmssmJacobian::calcFTJacobian(NmssmSoftsusy& model) {
    return calcFTJacobian(model, model.displayMxBC());
  }

  /// Calculates the Jacobian as in calcFTJacobian(NmssmSoftsusy& model),
  /// but with \c mx set to the given value instead of that returned by
  /// \c model.displayMxBC().
  double NmssmJacobian::calcFTJacobian(NmssmSoftsusy& model, double mx) {

    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();
    const sPhysical savedPhys(model.displayPhys());

    model.runto(model.calcMs());

    const double rgDet = calcRGFlowJacobian(model, model.displayMu(), mx);
    const double ewsbDet = calcEWSBJacobian(model);

    model.setMu(scale);
    model.set(savedPars);
    model.setPhys(savedPhys);

    return ewsbDet * rgDet;
  }

  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double NmssmJacobian::calcMz(NmssmSoftsusy& model, bool getRunningMass) {

    double mz = model.displayMzRun();
    if (!getRunningMass) {
      const double scale = model.displayMu();

      // @note currently this is using the given pole top mass, not the
      // running mass, check impact
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

    if (PRINTOUT > 1) cout << "# Qstart= " << startScale << ' ';

    switch (independent) {
    case Lambda: {
      tempModel->setLambda(x);
      if (PRINTOUT > 1) cout << "lambda= " << x << ' ';
      break;
    }
    case Kappa: {
      tempModel->setKappa(x);
      if (PRINTOUT > 1) cout << "kappa= " << x << ' ';
      break;
    }
    case SMu: {
      tempModel->setSusyMu(x);
      if (PRINTOUT > 1) cout << "mu= " << x << ' ';
      break;
    }
    case M3Sq: {
      tempModel->setM3Squared(x);
      if (PRINTOUT > 1) cout << "m3sq= " << x << ' ';
      break;
    }
    case XiS: {
      tempModel->setXiS(x);
      if (PRINTOUT > 1) cout << "xiS= " << x << ' ';
      break;
    }
    case Mh1Sq: {
      tempModel->setMh1Squared(x);
      if (PRINTOUT > 1) cout << "mH1Sq= " << x << ' ';
      break;
    }
    case Mh2Sq: {
      tempModel->setMh2Squared(x);
      if (PRINTOUT > 1) cout << "mH2Sq= " << x << ' ';
      break;
    }
    case MsSq: {
      tempModel->setMsSquared(x);
      if (PRINTOUT > 1) cout << "mSsq= " << x << ' ';
      break;
    }
    case Yt: {
      tempModel->setYukawaElement(YU, 3, 3, x);
      if (PRINTOUT > 1) cout << "ht= " << x << ' ';
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

    if (PRINTOUT > 1) cout << "Qend= " << endScale << ' ';

    double output;
    switch (dependent) {
    case Lambda: {
      output = tempModel->displayLambda();
      if (PRINTOUT > 1) cout << "lambda= " << output << '\n';
      break;
    }
    case Kappa: {
      output = tempModel->displayKappa();
      if (PRINTOUT > 1) cout << "kappa= " << output << '\n';
      break;
    }
    case SMu: {
      output = tempModel->displaySusyMu();
      if (PRINTOUT > 1) cout << "mu= " << output << '\n';
      break;
    }
    case M3Sq: {
      output = tempModel->displayM3Squared();
      if (PRINTOUT > 1) cout << "m3sq= " << output << '\n';
      break;
    }
    case XiS: {
      output = tempModel->displayXiS();
      if (PRINTOUT > 1) cout << "xiS= " << output << '\n';
      break;
    }
    case Mh1Sq: {
      output = tempModel->displayMh1Squared();
      if (PRINTOUT > 1) cout << "mH1Sq= " << output << '\n';
      break;
    }
    case Mh2Sq: {
      output = tempModel->displayMh2Squared();
      if (PRINTOUT > 1) cout << "mH2Sq= " << output << '\n';
      break;
    }
    case MsSq: {
      output = tempModel->displayMsSquared();
      if (PRINTOUT > 1) cout << "mSsq= " << output << '\n';
      break;
    }
    case Yt: {
      output = tempModel->displayYukawaElement(YU, 3, 3);
      if (PRINTOUT > 1) cout << "ht= " << output << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcRunningParameter called with incorrect"
         << " dependent parameter " << dependent << '\n';
      throw ii.str();
    }
    }

    tempModel->setMu(startScale);
    tempModel->set(savedPars);

    return output;
  }

  double NmssmJacobian::calcRGDerivative(NmssmSoftsusy& model,
                                         Parameters dep, Parameters indep,
                                         double toScale) {

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

    volatile const double temp = x + h;
    h = temp - x;

    RGFlowPars pars;
    pars.model = &model;
    pars.independent = indep;
    pars.dependent = dep;
    pars.toScale = toScale;

    double derivative = 0.;
    double err = 0.;
    if (fabs(x) > 1.0e-10) {
      derivative = calcDerivative(calcRunningParameter, x, h, &err, &pars);
    }

    if (PRINTOUT > 1)
      cout << "derivative=" << derivative << " error=" << err << '\n';

    const bool has_error
      = fabs(x) > 1.0e-10 && fabs(err / derivative) > 1.0;

    if (has_error) {
      derivative = -numberOfTheBeast;
      hasError = true;
    }

    return derivative;
  }

  double NmssmJacobian::calcRGFlowJacobian(NmssmSoftsusy& model,
                                           double startScale, double endScale) {
    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();

    model.runto(startScale);

    const int numPars = includeTop ? 4 : 3;

    DoubleMatrix jac(numPars, numPars);

    vector<Parameters> indepPars;

    if (SoftHiggsOut) {
      indepPars.push_back(Mh1Sq);
      indepPars.push_back(Mh2Sq);
      indepPars.push_back(MsSq);
    } else if (Z3) {
      indepPars.push_back(Lambda);
      indepPars.push_back(Kappa);
      indepPars.push_back(MsSq);
    } else {
      indepPars.push_back(SMu);
      indepPars.push_back(M3Sq);
      indepPars.push_back(XiS);
    }
    if (includeTop) indepPars.push_back(Yt);

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      if (SoftHiggsOut) {
        jac(i + 1, 1) = calcRGDerivative(model, Mh1Sq, indepPars[i], endScale);
        jac(i + 1, 2) = calcRGDerivative(model, Mh2Sq, indepPars[i], endScale);
        jac(i + 1, 3) = calcRGDerivative(model, MsSq, indepPars[i], endScale);
      } else if (Z3) {
        jac(i + 1, 1) = calcRGDerivative(model, Lambda, indepPars[i], endScale);
        jac(i + 1, 2) = calcRGDerivative(model, Kappa, indepPars[i], endScale);
        jac(i + 1, 3) = calcRGDerivative(model, MsSq, indepPars[i], endScale);
      } else {
        jac(i + 1, 1) = calcRGDerivative(model, SMu, indepPars[i], endScale);
        jac(i + 1, 2) = calcRGDerivative(model, M3Sq, indepPars[i], endScale);
        jac(i + 1, 3) = calcRGDerivative(model, XiS, indepPars[i], endScale);
      }
      if (includeTop)
        jac(i + 1, 4) = calcRGDerivative(model, Yt, indepPars[i], endScale);
    }

    // save calculated matrix
    // convention: inverse refers to case where transformation
    // is from high-scale to low-scale parameters
    if (startScale == endScale) {
      if (invJacRGFlow.displayRows() != numPars
          || invJacRGFlow.displayCols() != numPars) {
        invJacRGFlow.resize(numPars, numPars);
      }
      invJacRGFlow = jac;

      if (jacRGFlow.displayRows() != numPars
          || jacRGFlow.displayCols() != numPars) {
        jacRGFlow.resize(numPars, numPars);
      }
      jacRGFlow = jac;
    } else if (startScale > endScale) {
      if (invJacRGFlow.displayRows() != numPars
          || invJacRGFlow.displayCols() != numPars) {
        invJacRGFlow.resize(numPars, numPars);
      }
      invJacRGFlow = jac;
    } else {
      if (jacRGFlow.displayRows() != numPars
          || jacRGFlow.displayCols() != numPars) {
        jacRGFlow.resize(numPars, numPars);
      }
      jacRGFlow = jac;
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
    const double startScale = tempModel->displayMu();
    const drBarPars savedDrBarPars(tempModel->displayDrBarPars());

    if (PRINTOUT > 1) cout << '#';

    switch (independent) {
    case Lambda: {
      tempModel->setLambda(x);
      if (PRINTOUT > 1) cout << "lambda= " << x << ' ';
      break;
    }
    case Kappa: {
      tempModel->setKappa(x);
      if (PRINTOUT > 1) cout << "kappa= " << x << ' ';
      break;
    }
    case SMu: {
      tempModel->setSusyMu(x);
      if (PRINTOUT > 1) cout << "mu= " << x << ' ';
      break;
    }
    case M3Sq: {
      tempModel->setM3Squared(x);
      if (PRINTOUT > 1) cout << "m3sq= " << x << ' ';
      break;
    }
    case XiS: {
      tempModel->setXiS(x);
      if (PRINTOUT > 1) cout << "xiS= " << x << ' ';
      break;
    }
    case Mh1Sq: {
      tempModel->setMh1Squared(x);
      if (PRINTOUT > 1) cout << "mH1Sq= " << x << ' ';
      break;
    }
    case Mh2Sq: {
      tempModel->setMh2Squared(x);
      if (PRINTOUT > 1) cout << "mH2Sq= " << x << ' ';
      break;
    }
    case MsSq: {
      tempModel->setMsSquared(x);
      if (PRINTOUT > 1) cout << "mSsq= " << x << ' ';
      break;
    }
    case Yt: {
      tempModel->setYukawaElement(YU, 3, 3, x);
      if (PRINTOUT > 1) cout << "ht= " << x << ' ';
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBOutput called with incorrect"
         << " independent parameter " << independent << '\n';
      throw ii.str();
    }
    }

    DoubleVector vevs(3);
    vevs(1) = tempModel->displayHvev();
    vevs(2) = tempModel->displayTanb();
    vevs(3) = tempModel->displaySvev();

    int error = 0;
    tempModel->iterateVevs(vevs, error);

    if (error != 0) {
      if (PRINTOUT > 0) {
        cout << "Warning: could not solve for VEVs\n";
      }
    }

    tempModel->setHvev(vevs(1));
    tempModel->setTanb(vevs(2));
    tempModel->setSvev(vevs(3));

    tempModel->calcDrBarPars();
    const double mtrun = tempModel->displayDrBarPars().mt;
    const double sinthDRbar = tempModel->calcSinthdrbar();
    tempModel->doTadpoles(mtrun, sinthDRbar);

    double output;
    switch (dependent) {
    case Mzsq: {
      output = sqr(calcMz(*tempModel, pars->useRunningMasses));
      if (PRINTOUT > 1) cout << "MZ=" << sqrt(output) << '\n';
      break;
    }
    case Tanb: {
      output = tempModel->displayTanb();
      if (PRINTOUT > 1) cout << "tanb=" << output << '\n';
      break;
    }
    case Svev: {
      output = tempModel->displaySvev();
      if (PRINTOUT > 1) cout << "Svev=" << output << '\n';
      break;
    }
    case Mtsq: {
      output = sqr(calcMt(*tempModel, pars->useRunningMasses));
      if (PRINTOUT > 1) cout << "MT=" << sqrt(output) << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBOutput called with incorrect"
         << " dependent parameter " << dependent << '\n';
      throw ii.str();
    }
    }

    tempModel->setMu(startScale);
    tempModel->set(savedPars);
    tempModel->setDrBarPars(savedDrBarPars);

    return output;
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
    } else if (Z3) {
      m->setLambda(guess(1));
      m->setKappa(guess(2));
      m->setMsSquared(guess(3));
    } else {
      m->setSusyMu(guess(1));
      m->setM3Squared(guess(2));
      m->setXiS(guess(3));
    }
    if (numOutputs > 3) {
      m->setYukawaElement(YU, 3, 3, guess(4));
    }

    DoubleVector vevs(3);
    vevs(1) = m->displayHvev();
    vevs(2) = m->displayTanb();
    vevs(3) = m->displaySvev();

    int error = 0;
    m->iterateVevs(vevs, error);

    if (error != 0) {
      if (PRINTOUT > 0) {
        cout << "Warning: could not solve for VEVs\n";
      }
    }

    m->setHvev(vevs(1));
    m->setTanb(vevs(2));
    m->setSvev(vevs(3));

    m->calcDrBarPars();

    if (Z3 && !SoftHiggsOut) {
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
    } else if (Z3) {
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

    if (PRINTOUT > 1) cout << '#';

    switch (independent) {
    case Lambda: {
      pars->outputs(3) = x;
      if (PRINTOUT > 1) cout << "lambda= " << x << ' ';
      break;
    }
    case Mzsq: {
      pars->outputs(1) = x;
      if (PRINTOUT > 1) cout << "MZ= " << sqrt(x) << ' ';
      break;
    }
    case Tanb: {
      pars->outputs(2) = x;
      if (PRINTOUT > 1) cout << "tanb= " << x << ' ';
      break;
    }
    case Svev: {
      pars->outputs(3) = x;
      if (PRINTOUT > 1) cout << "Svev= " << x << ' ';
      break;
    }
    case Mtsq: {
      pars->outputs(4) = x;
      if (PRINTOUT > 1) cout << "MT= " << sqrt(x) << ' ';
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
      if (PRINTOUT > 0) {
        cout << "Warning: could not solve for VEVs\n";
      }
    }

    double output;
    switch (dependent) {
    case Lambda: {
      output = tempModel->displayLambda();
      if (PRINTOUT > 1) cout << "lambda= " << output << '\n';
      break;
    }
    case Kappa: {
      output = tempModel->displayKappa();
      if (PRINTOUT > 1) cout << "kappa= " << output << '\n';
      break;
    }
    case SMu: {
      output = tempModel->displaySusyMu();
      if (PRINTOUT > 1) cout << "mu= " << output << '\n';
      break;
    }
    case M3Sq: {
      output = tempModel->displayM3Squared();
      if (PRINTOUT > 1) cout << "m3sq= " << output << '\n';
      break;
    }
    case XiS: {
      output = tempModel->displayXiS();
      if (PRINTOUT > 1) cout << "xiS= " << output << '\n';
      break;
    }
    case Mh1Sq: {
      output = tempModel->displayMh1Squared();
      if (PRINTOUT > 1) cout << "mH1Sq= " << output << '\n';
      break;
    }
    case Mh2Sq: {
      output = tempModel->displayMh2Squared();
      if (PRINTOUT > 1) cout << "mH2Sq= " << output << '\n';
      break;
    }
    case MsSq: {
      output = tempModel->displayMsSquared();
      if (PRINTOUT > 1) cout << "mSsq= " << output << '\n';
      break;
    }
    case Yt: {
      output = tempModel->displayYukawaElement(YU, 3, 3);
      if (PRINTOUT > 1) cout << "ht= " << output << '\n';
      break;
    }
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcEWSBParameter called with incorrect"
         << " dependent parameter " << dependent << '\n';
      throw ii.str();
    }
    }

    tempModel->setMu(startScale);
    tempModel->set(savedPars);
    tempModel->setDrBarPars(savedDrBarPars);
    pars->outputs = savedOutputs;

    return output;
  }

  double NmssmJacobian::calcEWSBOutputDerivative(
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

    volatile const double temp = x + h;
    h = temp - x;

    EWSBPars pars;
    pars.model = &model;
    pars.independent = indep;
    pars.dependent = dep;
    pars.useRunningMasses = useRunningMasses;

    double derivative = 0.;
    double err = 0.;
    if (fabs(x) > 1.0e-10) {
      derivative = calcDerivative(calcEWSBOutput, x, h, &err, &pars);
    }

    if (PRINTOUT > 1)
      cout << "derivative=" << derivative << " error=" << err << '\n';

    const bool has_error
      = fabs(x) > 1.0e-10 && fabs(err / derivative) > 1.0;

    if (has_error) {
      derivative = -numberOfTheBeast;
      hasError = true;
    }

    return derivative;
  }

  double NmssmJacobian::calcEWSBParameterDerivative(
    NmssmSoftsusy& model, Parameters dep, Parameters indep) {

    double x = 0.;
    double h = 0.01;

    const int numOutputs = includeTop ? 4 : 3;
    DoubleVector outputs(numOutputs);
    outputs(1) = sqr(calcMz(model, useRunningMasses));
    outputs(2) = model.displayTanb();
    outputs(3) = (Z3 && !SoftHiggsOut) ? model.displayLambda() :
      model.displaySvev();
    if (numOutputs > 3) outputs(4) = sqr(calcMt(model, useRunningMasses));

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

    volatile const double temp = x + h;
    h = temp - x;

    EWSBPars pars;
    pars.model = &model;
    pars.independent = indep;
    pars.dependent = dep;
    pars.outputs = outputs;
    pars.useRunningMasses = useRunningMasses;

    double derivative = 0.;
    double err = 0.;
    if (fabs(x) > 1.0e-10) {
      derivative = calcDerivative(calcEWSBParameter, x, h, &err, &pars);
    }

    if (PRINTOUT > 1)
      cout << "derivative=" << derivative << " error=" << err << '\n';

    const bool has_error
      = fabs(x) > 1.0e-10 && fabs(err / derivative) > 1.0;

    if (has_error) {
      derivative = -numberOfTheBeast;
      hasError = true;
    }

    return derivative;
  }

  double NmssmJacobian::calcEWSBJacobian(NmssmSoftsusy& model) {

    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();
    const drBarPars savedDrBarPars(model.displayDrBarPars());

    const int numPars = includeTop ? 4 : 3;

    DoubleMatrix jac(numPars, numPars);

    vector<Parameters> indepPars;

    if (Z3 && !SoftHiggsOut) {
      indepPars.push_back(Mzsq);
      indepPars.push_back(Tanb);
      indepPars.push_back(Lambda);
    } else {
      indepPars.push_back(Mzsq);
      indepPars.push_back(Tanb);
      indepPars.push_back(Svev);
    }
    if (includeTop) indepPars.push_back(Mtsq);

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      if (SoftHiggsOut) {
        jac(i + 1, 1) = calcEWSBParameterDerivative(model, Mh1Sq, indepPars[i]);
        jac(i + 1, 2) = calcEWSBParameterDerivative(model, Mh2Sq, indepPars[i]);
        jac(i + 1, 3) = calcEWSBParameterDerivative(model, MsSq, indepPars[i]);
      } else if (Z3) {
        if (indepPars[i] == Lambda) {
          jac(i + 1, 1) = 1.;
          jac(i + 1, 2) = 0.;
          jac(i + 1, 3) = 0.;
        } else {
          jac(i + 1, 1) = 0.;
          jac(i + 1, 2) = calcEWSBParameterDerivative(model, Kappa, indepPars[i]);
          jac(i + 1, 3) = calcEWSBParameterDerivative(model, MsSq, indepPars[i]);
        }
      } else {
        jac(i + 1, 1) = calcEWSBParameterDerivative(model, SMu, indepPars[i]);
        jac(i + 1, 2) = calcEWSBParameterDerivative(model, M3Sq, indepPars[i]);
        jac(i + 1, 3) = calcEWSBParameterDerivative(model, XiS, indepPars[i]);
      }
      if (includeTop) {
        if ((Z3 && !SoftHiggsOut) && indepPars[i] == Lambda) {
          jac(i + 1, 4) = 0.;
        } else {
          jac(i + 1, 4) = calcEWSBParameterDerivative(model, Yt, indepPars[i]);
        }
      }
    }

    // save calculated matrix
    if (jacEWSB.displayRows() != numPars
        || jacEWSB.displayCols() != numPars) {
      jacEWSB.resize(numPars, numPars);
    }
    jacEWSB = jac;

    model.setMu(scale);
    model.set(savedPars);
    model.setDrBarPars(savedDrBarPars);

    return jac.determinant();
  }

  double NmssmJacobian::calcInverseEWSBJacobian(NmssmSoftsusy& model) {

    const DoubleVector savedPars(model.display());
    const double scale = model.displayMu();
    const drBarPars savedDrBarPars(model.displayDrBarPars());

    const int numPars = includeTop ? 4 : 3;

    DoubleMatrix jac(numPars, numPars);

    vector<Parameters> indepPars;

    if (SoftHiggsOut) {
      indepPars.push_back(Mh1Sq);
      indepPars.push_back(Mh2Sq);
      indepPars.push_back(MsSq);
    } else if (Z3) {
      indepPars.push_back(Lambda);
      indepPars.push_back(Kappa);
      indepPars.push_back(MsSq);
    } else {
      indepPars.push_back(SMu);
      indepPars.push_back(M3Sq);
      indepPars.push_back(XiS);
    }
    if (includeTop) indepPars.push_back(Yt);

    for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
      if (Z3 && !SoftHiggsOut) {
        if (indepPars[i] == Lambda) {
          jac(i + 1, 1) = 0.;
          jac(i + 1, 2) = 0.;
          jac(i + 1, 3) = 1.;
        } else {
          jac(i + 1, 1) = calcEWSBOutputDerivative(model, Mzsq, indepPars[i]);
          jac(i + 1, 2) = calcEWSBOutputDerivative(model, Tanb, indepPars[i]);
          jac(i + 1, 3) = 0.;
        }
      } else {
        jac(i + 1, 1) = calcEWSBOutputDerivative(model, Mzsq, indepPars[i]);
        jac(i + 1, 2) = calcEWSBOutputDerivative(model, Tanb, indepPars[i]);
        jac(i + 1, 3) = calcEWSBOutputDerivative(model, Svev, indepPars[i]);
      }
      if (includeTop) {
        if ((Z3 && !SoftHiggsOut) && indepPars[i] == Lambda) {
          jac(i + 1, 4) = 0.;
        } else {
          jac(i + 1, 4) = calcEWSBOutputDerivative(model, Mtsq, indepPars[i]);
        }
      }
    }

    // save calculated matrix
    if (invJacEWSB.displayRows() != numPars
        || invJacEWSB.displayCols() != numPars) {
      invJacEWSB.resize(numPars, numPars);
    }
    invJacEWSB = jac;

    model.setMu(scale);
    model.set(savedPars);
    model.setDrBarPars(savedDrBarPars);

    return jac.determinant();
  }

} /// namespace softsusy
