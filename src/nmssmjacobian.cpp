// -*- mode: c++; c-basic-offset: 2; -*-
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

  NmssmJacobian::NmssmJacobian(NmssmSoftsusy* m)
     : model(m), jacRGFlow(3,3), jacEWSB(3,3)
     , invJacRGFlow(3,3), invJacEWSB(3,3), useRunningMasses(false) {
  }

  NmssmJacobian::~NmssmJacobian() {}

  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double NmssmJacobian::calcMz(NmssmSoftsusy* m, bool getRunningMass) {
    double mzpole = 0.;

    if (m) {
      const double mzrun = m->displayMzRun();
      if (getRunningMass) {
        mzpole = mzrun;
      } else {
        const double scale = m->displayMu();

        // @note currently this is using the given pole top mass, not the
        // running mass, check impact
        const double pizzt = m->piZZT(mzrun, scale);

        mzpole = sqrt(fabs(sqr(mzrun) - pizzt));
      }
    }

    return mzpole;
  }

  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double NmssmJacobian::calcMt(NmssmSoftsusy* m, bool getRunningMass) {
    double mtpole = 0.;

    if (m) {
      const double mtrun = m->displayDrBarPars().mt;
      if (getRunningMass) {
        mtpole = mtrun;
      } else {
        const double scale = m->displayMu();

        // workaround hard-coded external momentum dependence
        QedQcd savedData(m->displayDataSet());
        QedQcd tempData;
        tempData.setPoleMt(mtrun);
        m->setData(tempData);

        const double stopGluino = m->calcRunMtStopGluino();
        const double higgs = m->calcRunMtHiggs();
        const double neutralinos = m->calcRunMtNeutralinos();
        const double charginos = m->calcRunMtCharginos();

        double resigmat = mtrun * (stopGluino + higgs + neutralinos
                                   + charginos) / (16.0 * sqr(PI));

        const double g3Sq = sqr(m->displayGaugeCoupling(3));
        const double logMtSqOverQSq = 2.0 * log(mtrun / scale);
        const double oneLoopQCD = 4.0 * g3Sq * (5.0 - 3.0 * logMtSqOverQSq)
          / (3.0 * 16.0 * sqr(PI));
        const double twoLoopQCD = sqr(g3Sq) * (0.005191204615668296
          - 0.0032883224409535764 * logMtSqOverQSq + 0.0008822328500119351 *
          sqr(logMtSqOverQSq));
        resigmat -= mtrun * (oneLoopQCD + twoLoopQCD);

        mtpole = mtrun - resigmat;

        m->setData(savedData);
      }
    }

    return mtpole;
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

  double NmssmJacobian::calcRGDerivative(Parameters dep, Parameters indep,
                                         double toScale) {
    double derivative = 0.;

    if (model) {
      double x = 0.;
      double h = 0.01;

      switch (indep) {
      case Lambda: {
        x = model->displayLambda();
        h = 0.0005 * x;
        break;
      }
      case Kappa: {
        x = model->displayKappa();
        h = 0.0005 * x;
        break;
      }
      case SMu: {
        x = model->displaySusyMu();
        h = 0.01 * x;
        break;
      }
      case M3Sq: {
        x = model->displayM3Squared();
        h = 0.01 * x;
        break;
      }
      case XiS: {
        x = model->displayXiS();
        h = 0.01 * x;
        break;
      }
      case Mh1Sq: {
        x = model->displayMh1Squared();
        h = 0.01 * x;
        break;
      }
      case Mh2Sq: {
        x = model->displayMh2Squared();
        h = 0.01 * x;
        break;
      }
      case MsSq: {
        x = model->displayMsSquared();
        h = 0.01 * x;
        break;
      }
      case Yt: {
        x = model->displayYukawaElement(YU, 3, 3);
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
      pars.model = model;
      pars.independent = indep;
      pars.dependent = dep;
      pars.toScale = toScale;

      double err = 0.;
      if (fabs(x) > 1.0e-10) {
        derivative = calcDerivative(calcRunningParameter, x, h, &err, &pars);
      }

      if (PRINTOUT > 1)
        cout << "derivative=" << derivative << " error=" << err << '\n';

      const bool has_error
        = fabs(x) > 1.0e-10 && fabs(err / derivative) > 1.0;

      if (has_error) derivative = -numberOfTheBeast;
    }

    return derivative;
  }

  double NmssmJacobian::calcRGFlowJacobian(double startScale, double endScale,
                                           bool doTop) {
    double rgDet = 0.;

    if (model) {
      const DoubleVector savedPars(model->display());
      const double scale = model->displayMu();

      model->runto(startScale);

      const int numPars = doTop ? 4 : 3;

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
      if (doTop) indepPars.push_back(Yt);

      for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
        if (SoftHiggsOut) {
          jac(i + 1, 1) = calcRGDerivative(Mh1Sq, indepPars[i], endScale);
          jac(i + 1, 2) = calcRGDerivative(Mh2Sq, indepPars[i], endScale);
          jac(i + 1, 3) = calcRGDerivative(MsSq, indepPars[i], endScale);
        } else if (Z3) {
          jac(i + 1, 1) = calcRGDerivative(Lambda, indepPars[i], endScale);
          jac(i + 1, 2) = calcRGDerivative(Kappa, indepPars[i], endScale);
          jac(i + 1, 3) = calcRGDerivative(MsSq, indepPars[i], endScale);
        } else {
          jac(i + 1, 1) = calcRGDerivative(SMu, indepPars[i], endScale);
          jac(i + 1, 2) = calcRGDerivative(M3Sq, indepPars[i], endScale);
          jac(i + 1, 3) = calcRGDerivative(XiS, indepPars[i], endScale);
        }
        if (doTop)
          jac(i + 1, 4) = calcRGDerivative(Yt, indepPars[i], endScale);
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

      rgDet = jac.determinant();

      model->setMu(scale);
      model->set(savedPars);
    }

    return rgDet;
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

    // @todo better error handling
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
      output = sqr(calcMz(tempModel, pars->useRunningMasses));
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
      output = sqr(calcMt(tempModel, pars->useRunningMasses));
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

    m->setHvev(vevs(1));
    m->setTanb(vevs(2));
    m->setSvev(vevs(3));

    m->calcDrBarPars();

    if (Z3 && !SoftHiggsOut) {
      errors(1) = 1.0 - (sqr(calcMz(m, pars->useRunningMasses)) / outputs(1));
      errors(2) = 1.0 - (m->displayTanb() / outputs(2));
      errors(3) = 1.0 - (m->displayLambda() / outputs(3));
    } else {
      errors(1) = 1.0 - (sqr(calcMz(m, pars->useRunningMasses)) / outputs(1));
      errors(2) = 1.0 - (m->displayTanb() / outputs(2));
      errors(3) = 1.0 - (m->displaySvev() / outputs(3));
    }

    error = error && testNan(errors(1)) && testNan(errors(2))
      && testNan(errors(3));

    if (numOutputs > 3) {
      errors(4) = 1.0 - (sqr(calcMt(m, pars->useRunningMasses)) / outputs(4));
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

    bool error = newt(guess, ewsbOutputErrors, pars);
    err = error ? 1 : 0;
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

    // @todo better error handling
    if (error != 0) {
      if (PRINTOUT > 0) {
        cout << "Warning: could not set EWSB outputs\n";
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

  double NmssmJacobian::calcEWSBOutputDerivative(Parameters dep, Parameters indep) {
    double derivative = 0.;

    if (model) {
      double x = 0.;
      double h = 0.01;

      switch (indep) {
      case Lambda: {
        x = model->displayLambda();
        h = 0.0005 * x;
        break;
      }
      case Kappa: {
        x = model->displayKappa();
        h = 0.0005 * x;
        break;
      }
      case SMu: {
        x = model->displaySusyMu();
        h = 0.01 * x;
        break;
      }
      case M3Sq: {
        x = model->displayM3Squared();
        h = 0.01 * x;
        break;
      }
      case XiS: {
        x = model->displayXiS();
        h = 0.01 * x;
        break;
      }
      case Mh1Sq: {
        x = model->displayMh1Squared();
        h = 0.01 * x;
        break;
      }
      case Mh2Sq: {
        x = model->displayMh2Squared();
        h = 0.01 * x;
        break;
      }
      case MsSq: {
        x = model->displayMsSquared();
        h = 0.01 * x;
        break;
      }
      case Yt: {
        x = model->displayYukawaElement(YU, 3, 3);
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
      pars.model = model;
      pars.independent = indep;
      pars.dependent = dep;
      pars.useRunningMasses = useRunningMasses;

      double err = 0.;
      if (fabs(x) > 1.0e-10) {
        derivative = calcDerivative(calcEWSBOutput, x, h, &err, &pars);
      }

      if (PRINTOUT > 1)
        cout << "derivative=" << derivative << " error=" << err << '\n';

      const bool has_error
        = fabs(x) > 1.0e-10 && fabs(err / derivative) > 1.0;

      if (has_error) derivative = -numberOfTheBeast;
    }

    return derivative;
  }

  double NmssmJacobian::calcEWSBParameterDerivative(Parameters dep, Parameters indep, bool doTop) {
    double derivative = 0.;

    if (model) {
      double x = 0.;
      double h = 0.01;

      const int numOutputs = doTop ? 4 : 3;
      DoubleVector outputs(numOutputs);
      outputs(1) = sqr(calcMz(model, useRunningMasses));
      outputs(2) = model->displayTanb();
      outputs(3) = (Z3 && !SoftHiggsOut) ? model->displayLambda() :
        model->displaySvev();
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
      pars.model = model;
      pars.independent = indep;
      pars.dependent = dep;
      pars.outputs = outputs;
      pars.useRunningMasses = useRunningMasses;

      double err = 0.;
      if (fabs(x) > 1.0e-10) {
        derivative = calcDerivative(calcEWSBParameter, x, h, &err, &pars);
      }

      if (PRINTOUT > 1)
        cout << "derivative=" << derivative << " error=" << err << '\n';

      const bool has_error
        = fabs(x) > 1.0e-10 && fabs(err / derivative) > 1.0;

      if (has_error) derivative = -numberOfTheBeast;
    }

    return derivative;
  }

  double NmssmJacobian::calcEWSBJacobian(bool doTop) {
    double ewsbDet = 0.;

    if (model) {
      const DoubleVector savedPars(model->display());
      const double scale = model->displayMu();
      const drBarPars savedDrBarPars(model->displayDrBarPars());

      const int numPars = doTop ? 4 : 3;

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
      if (doTop) indepPars.push_back(Mtsq);

      for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
        if (SoftHiggsOut) {
          jac(i + 1, 1) = calcEWSBParameterDerivative(Mh1Sq, indepPars[i], doTop);
          jac(i + 1, 2) = calcEWSBParameterDerivative(Mh2Sq, indepPars[i], doTop);
          jac(i + 1, 3) = calcEWSBParameterDerivative(MsSq, indepPars[i], doTop);
        } else if (Z3) {
          if (indepPars[i] == Lambda) {
            jac(i + 1, 1) = 1.;
            jac(i + 1, 2) = 0.;
            jac(i + 1, 3) = 0.;
          } else {
            jac(i + 1, 1) = 0.;
            jac(i + 1, 2) = calcEWSBParameterDerivative(Kappa, indepPars[i], doTop);
            jac(i + 1, 3) = calcEWSBParameterDerivative(MsSq, indepPars[i], doTop);
          }
        } else {
          jac(i + 1, 1) = calcEWSBParameterDerivative(SMu, indepPars[i], doTop);
          jac(i + 1, 2) = calcEWSBParameterDerivative(M3Sq, indepPars[i], doTop);
          jac(i + 1, 3) = calcEWSBParameterDerivative(XiS, indepPars[i], doTop);
        }
        if (doTop) {
          if ((Z3 && !SoftHiggsOut) && indepPars[i] == Lambda) {
            jac(i + 1, 4) = 0.;
          } else {
            jac(i + 1, 4) = calcEWSBParameterDerivative(Yt, indepPars[i], doTop);
          }
        }
      }

      // save calculated matrix
      if (jacEWSB.displayRows() != numPars
          || jacEWSB.displayCols() != numPars) {
        jacEWSB.resize(numPars, numPars);
      }
      jacEWSB = jac;

      ewsbDet = jac.determinant();

      model->setMu(scale);
      model->set(savedPars);
      model->setDrBarPars(savedDrBarPars);
    }

    return ewsbDet;
  }

  double NmssmJacobian::calcInverseEWSBJacobian(bool doTop) {
    double ewsbDet = 0.;

    if (model) {
      const DoubleVector savedPars(model->display());
      const double scale = model->displayMu();
      const drBarPars savedDrBarPars(model->displayDrBarPars());

      const int numPars = doTop ? 4 : 3;

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
      if (doTop) indepPars.push_back(Yt);

      for (int i = 0, numIndep = indepPars.size(); i < numIndep; ++i) {
        if (Z3 && !SoftHiggsOut) {
          if (indepPars[i] == Lambda) {
            jac(i + 1, 1) = 0.;
            jac(i + 1, 2) = 0.;
            jac(i + 1, 3) = 1.;
          } else {
            jac(i + 1, 1) = calcEWSBOutputDerivative(Mzsq, indepPars[i]);
            jac(i + 1, 2) = calcEWSBOutputDerivative(Tanb, indepPars[i]);
            jac(i + 1, 3) = 0.;
          }
        } else {
          jac(i + 1, 1) = calcEWSBOutputDerivative(Mzsq, indepPars[i]);
          jac(i + 1, 2) = calcEWSBOutputDerivative(Tanb, indepPars[i]);
          jac(i + 1, 3) = calcEWSBOutputDerivative(Svev, indepPars[i]);
        }
        if (doTop) {
          if ((Z3 && !SoftHiggsOut) && indepPars[i] == Lambda) {
            jac(i + 1, 4) = 0.;
          } else {
            jac(i + 1, 4) = calcEWSBOutputDerivative(Mtsq, indepPars[i]);
          }
        }
      }

      // save calculated matrix
      if (invJacEWSB.displayRows() != numPars
          || invJacEWSB.displayCols() != numPars) {
        invJacEWSB.resize(numPars, numPars);
      }
      invJacEWSB = jac;

      ewsbDet = jac.determinant();

      model->setMu(scale);
      model->set(savedPars);
      model->setDrBarPars(savedDrBarPars);
    }

    return ewsbDet;
  }

  /// The Jacobian calculated is of the form \f$ J^{-1} = |
  /// \partial O / \partial p | \f$.  The transformation is done
  /// in two stages.  In the first, the observables at the
  /// current scale are traded for the parameters at this
  /// scale using the EWSB conditions.  The Jacobian matrix for
  /// this transformation can be accessed afterwards using
  /// displayInverseEWSBJacobian().  Then, the parameters
  /// at this scale are transformed to parameters at the
  /// scale \c mx using the RGEs.  The Jacobian matrix for this
  /// transformation may be accessed by calling
  /// displayInverseRGFlowJacobian().
  ///
  /// The sets of observables \f$\{O\}\f$ and parameters
  /// \f$\{p\}\f$ are selected using the global flags softsusy::Z3
  /// and softsusy::SoftHiggsOut.  If softsusy::SoftHiggsOut is true,
  /// the parameters used are the soft Higgs masses at the scale
  /// \c mx, \f$\{ m_{H_1,0}^2, m_{H_2,0}^2, m_{S_0}^2\} \f$, and
  /// the observables are \f$\{M_Z^2, \tan\beta, s\}\f$, irrespective
  /// of the value of softsusy::Z3.  If softsusy::SoftHiggsOut is false
  /// and softsusy::Z3 is true, the parameters are taken to be
  /// \f$\{\lambda_0, \kappa_0, m_{S_0}^2\}\f$ at \c mx.  The
  /// observables in this case are \f$\{M_Z^2, \tan\beta,
  /// \lambda\}\f$.  If both flags are false, the parameters
  /// used are \f$\{\mu_0, m_{3_0}^2, \xi_{S_0}\}\f$ and the set of observables
  /// is \f$\{M_Z^2, \tan\beta, s\}\f$.
  ///
  /// If \c doTop is true, the top Yukawa \f$y_t\f$ will also be included
  /// in the set of parameters, and \f$M_t^2\f$ in the set of observables.
  /// Whether \f$M_Z^2\f$ and \f$M_t^2\f$ are taken to be the pole or
  /// running masses may be set by calling setUseRunningMassesFlag()
  /// with the appropriate flag.  By default the pole masses are used.
  double NmssmJacobian::calcFTInverseJacobian(double mx, bool doTop) {
    double determinant = 0.;

    if (model) {
      const DoubleVector savedPars(model->display());
      const double scale = model->displayMu();
      const sPhysical savedPhys(model->displayPhys());

      const double rgDet = calcRGFlowJacobian(mx, scale, doTop);
      const double ewsbDet = calcInverseEWSBJacobian(doTop);

      determinant = ewsbDet * rgDet;

      model->setMu(scale);
      model->set(savedPars);
      model->setPhys(savedPhys);
    }

    return determinant;
  }

  /// The Jacobian calculated is of the form \f$ J = |
  /// \partial p / \partial O | \f$.  The transformation,
  /// which is the inverse of that calculated by
  /// calcFTInverseJacobian(), is done in two stages.  In the
  /// first, the parameters at the scale \c mx are transformed
  /// to parameters \f$\{q\}\f$ at the current scale using the
  /// RGEs.  The Jacobian matrix for this
  /// transformation may be accessed afterwards by calling
  /// displayRGFlowJacobian().  Then the parameters
  /// \f$\{q\}\f$ are traded for the observables at the same
  /// scale using the EWSB conditions.  The Jacobian matrix for
  /// this second transformation can be accessed using
  /// displayEWSBJacobian().
  ///
  /// The selection of the sets of observables and parameters
  /// is based on the global flags softsusy::Z3 and
  /// softsusy::SoftHiggsOut, and is the same as is used in
  /// calcFTInverseJacobian().  The effect of \c doTop is
  /// the same as for calcFTInverseJacobian() as well.
  double NmssmJacobian::calcFTJacobian(double mx, bool doTop) {
    double determinant = 0.;

    if (model) {
      const DoubleVector savedPars(model->display());
      const double scale = model->displayMu();
      const sPhysical savedPhys(model->displayPhys());

      const double rgDet = calcRGFlowJacobian(scale, mx, doTop);
      const double ewsbDet = calcEWSBJacobian(doTop);

      determinant = ewsbDet * rgDet;

      model->setMu(scale);
      model->set(savedPars);
      model->setPhys(savedPhys);
    }

    return determinant;
  }

  double NmssmJacobian::calcDeltaJ(double mx, bool doTop) {

    double tuning = 0.;

    if (model) {
      const DoubleVector savedPars(model->display());
      const double scale = model->displayMu();
      const sPhysical savedPhys(model->displayPhys());

      const double determinant = calcFTInverseJacobian(mx, doTop);

      const double mz2 = sqr(calcMz(model, useRunningMasses));
      const double tb = model->displayTanb();
      const double mt2 = sqr(calcMt(model, useRunningMasses));

      double denominator = mz2 * tb;
      if (Z3) {
        denominator *= model->displayLambda();
      } else {
        denominator *= model->displaySvev();
      }
      if (doTop) denominator *= mt2;

      model->runto(mx);

      double numerator;
      if (SoftHiggsOut) {
        numerator = model->displayMh1Squared() * model->displayMh2Squared()
          * model->displayMsSquared();
      } else if (Z3) {
        numerator = model->displayLambda() * model->displayKappa()
          * model->displayMsSquared();
      } else {
        numerator = model->displaySusyMu() * model->displayM3Squared()
          * model->displayXiS();
      }
      if (doTop) numerator *= model->displayYukawaElement(YU, 3, 3);

      model->setMu(scale);
      model->set(savedPars);
      model->setPhys(savedPhys);

      tuning = fabs(numerator * determinant / denominator);
    }

    return tuning;
  }

} /// namespace softsusy
