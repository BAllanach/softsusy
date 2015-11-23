// -*- mode: c++; c-basic-offset: 2; -*-
/** \file nmssmjacobian.cpp
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
     , invJacRGFlow(3,3), invJacEWSB(3,3) {
  }

  NmssmJacobian::~NmssmJacobian() {}

  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double NmssmJacobian::calcMzPole(NmssmSoftsusy* m) {
    double mzpole = 0.;

    if (m) {
      const double scale = m->displayMu();
      const double mzrun = m->displayMzRun();
      // @note currently this is using the given pole top mass, not the running
      // mass, check impact
      const double pizzt = m->piZZT(mzrun, scale);

      mzpole = sqrt(fabs(sqr(mzrun) - pizzt));
    }

    return mzpole;
  }

  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double NmssmJacobian::calcMtPole(NmssmSoftsusy* m) {
    double mtpole = 0.;

    if (m) {
      const double mtrun = m->displayDrBarPars().mt;
      const double scale = m->displayMu();

      // workaround hard-coded external momentum dependence
      QedQcd savedData(m->displayDataSet());
      QedQcd tempData;
      tempData.setPoleMt(mtrun);
      m->setData(tempData);

      const double qcd = m->calcRunMtQCD();
      const double stopGluino = m->calcRunMtStopGluino();
      const double higgs = m->calcRunMtHiggs();
      const double neutralinos = m->calcRunMtNeutralinos();
      const double charginos = m->calcRunMtCharginos();

      double resigmat = mtrun * (qcd + stopGluino + higgs + neutralinos
                                 + charginos) / (16.0 * sqr(PI));

      const double logMtSqOverQSq = 2.0 * log(mtrun / scale);
      const double twoLoopQCD = sqr(sqr(m->displayGaugeCoupling(3))) *
        (-0.538314 + 0.181534*logMtSqOverQSq - 0.0379954*sqr(logMtSqOverQSq));
      resigmat += mtrun * twoLoopQCD / (16.0 * sqr(PI));

      mtpole = mtrun - resigmat;

      m->setData(savedData);
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
      if (startScale > endScale) {
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
      output = sqr(calcMzPole(tempModel));
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
      output = sqr(calcMtPole(tempModel));
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

  double NmssmJacobian::calcEWSBDerivative(Parameters dep, Parameters indep) {
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
        ii << "NmssmJacobian:calcEWSBDerivative called with incorrect"
           << " independent parameter " << indep << '\n';
        throw ii.str();
      }
      }

      EWSBPars pars;
      pars.model = model;
      pars.independent = indep;
      pars.dependent = dep;

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

  double NmssmJacobian::calcEWSBJacobian(bool doTop) {
    double ewsbDet = 0.;

    if (model) {

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
            jac(i + 1, 1) = calcEWSBDerivative(Mzsq, indepPars[i]);
            jac(i + 1, 2) = calcEWSBDerivative(Tanb, indepPars[i]);
            jac(i + 1, 3) = 0.;
          }
        } else {
          jac(i + 1, 1) = calcEWSBDerivative(Mzsq, indepPars[i]);
          jac(i + 1, 2) = calcEWSBDerivative(Tanb, indepPars[i]);
          jac(i + 1, 3) = calcEWSBDerivative(Svev, indepPars[i]);
        }
        if (doTop) {
          if ((Z3 && !SoftHiggsOut) && indepPars[i] == Lambda) {
            jac(i + 1, 4) = 0.;
          } else {
            jac(i + 1, 4) = calcEWSBDerivative(Mtsq, indepPars[i]);
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

  double NmssmJacobian::calcDeltaJ(double mx, bool doTop) {

    double tuning = 0.;

    if (model) {
      const DoubleVector savedPars(model->display());
      const double scale = model->displayMu();
      const sPhysical savedPhys(model->displayPhys());

      const double determinant = calcFTInverseJacobian(mx, doTop);

      const double mz2 = sqr(calcMzPole(model));
      const double tb = model->displayTanb();
      const double mt2 = sqr(calcMtPole(model));

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
