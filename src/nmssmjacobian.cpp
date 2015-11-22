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
  double NmssmJacobian::calcMzPole() const {
    double mzpole = 0.;

    if (model) {
      const double scale = model->displayMu();
      const double mzrun = model->displayMzRun();
      // @note currently this is using the given pole top mass, not the running
      // mass, check impact
      const double pizzt = model->piZZT(mzrun, scale);

      mzpole = sqrt(fabs(sqr(mzrun) - pizzt));
    }

    return mzpole;
  }

  // @todo make sure this is consistent (e.g. is scale dependence right?)
  double NmssmJacobian::calcMtPole() {
    double mtpole = 0.;

    if (model) {
      const double mtrun = model->displayDrBarPars().mt;
      const double scale = model->displayMu();

      // workaround hard-coded external momentum dependence
      QedQcd savedData(model->displayDataSet());
      QedQcd tempData;
      tempData.setPoleMt(mtrun);
      model->setData(tempData);

      const double qcd = model->calcRunMtQCD();
      const double stopGluino = model->calcRunMtStopGluino();
      const double higgs = model->calcRunMtHiggs();
      const double neutralinos = model->calcRunMtNeutralinos();
      const double charginos = model->calcRunMtCharginos();

      double resigmat = mtrun * (qcd + stopGluino + higgs + neutralinos
                                 + charginos) / (16.0 * sqr(PI));

      const double logMtSqOverQSq = 2.0 * log(mtrun / scale);
      const double twoLoopQCD = sqr(sqr(model->displayGaugeCoupling(3))) * 
        (-0.538314 + 0.181534*logMtSqOverQSq - 0.0379954*sqr(logMtSqOverQSq));
      resigmat += mtrun * twoLoopQCD / (16.0 * sqr(PI));

      mtpole = mtrun - resigmat;

      model->setData(savedData);
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

    if (PRINTOUT > 1) cout << '#';

    switch (independent) {
    case Lambda: tempModel->setLambda(x); break;
    case Kappa: tempModel->setKappa(x); break;
    case SMu: tempModel->setSusyMu(x); break;
    case M3Sq: tempModel->setM3Squared(x); break;
    case XiS: tempModel->setXiS(x); break;
    case Mh1Sq: tempModel->setMh1Squared(x); break;
    case Mh2Sq: tempModel->setMh2Squared(x); break;
    case MsSq: tempModel->setMsSquared(x); break;
    case Yt: tempModel->setYukawaElement(YU, 3, 3, x); break;
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcRunningParameter called with incorrect"
         << " independent parameter " << independent << '\n';
      throw ii.str();
    }
    }

    tempModel->runto(endScale);

    double output;
    switch (dependent) {
    case Lambda: output = tempModel->displayLambda(); break;
    case Kappa: output = tempModel->displayKappa(); break;
    case SMu: output = tempModel->displaySusyMu(); break;
    case M3Sq: output = tempModel->displayM3Squared(); break;
    case XiS: output = tempModel->displayXiS(); break;
    case Mh1Sq: output = tempModel->displayMh1Squared(); break;
    case Mh2Sq: output = tempModel->displayMh2Squared(); break;
    case MsSq: output = tempModel->displayMsSquared(); break;
    case Yt: output = tempModel->displayYukawaElement(YU, 3, 3); break;
    default: {
      ostringstream ii;
      ii << "NmssmJacobian:calcRunningParameter called with incorrect"
         << " dependent parameter " << independent << '\n';
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
        ii << "NmssmJacobian:calcRunningParameter called with incorrect"
           << " dependent parameter " << indep << '\n';
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
      const sPhysical savedPhys(model->displayPhys());

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
      model->setPhys(savedPhys);
    }

    return rgDet;
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
      const sPhysical savedPhys(model->displayPhys());

      const int numPars = doTop ? 4 : 3;

      DoubleVector pars(numPars);

      if (SoftHiggsOut) {
        pars(1) = model->displayMh1Squared();
        pars(2) = model->displayMh2Squared();
        pars(3) = model->displayMsSquared();
      } else if (Z3) {
        pars(1) = model->displayLambda();
        pars(2) = model->displayKappa();
        pars(3) = model->displayMsSquared();
      } else {
        pars(1) = model->displaySusyMu();
        pars(2) = model->displayM3Squared();
        pars(3) = model->displayXiS();
      }
      if (doTop) {
        pars(4) = model->displayYukawaElement(YU, 3, 3);
      }

      DoubleMatrix jac(numPars, numPars);



      // save calculated matrix
      if (invJacEWSB.displayRows() != numPars
          || invJacEWSB.displayCols() != numPars) {
        invJacEWSB.resize(numPars, numPars);
      }
      invJacEWSB = jac;

      ewsbDet = jac.determinant();

      model->setMu(scale);
      model->set(savedPars);
      model->setPhys(savedPhys);
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

      const double mz2 = sqr(calcMzPole());
      const double tb = model->displayTanb();
      const double mt2 = sqr(calcMtPole());

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
