/* -*- mode: c++; c-basic-offset: 2; -*- */
/** \file nmssmjacobian.h

    \brief provides routines for calculating Jacobian fine tuning
 */

#ifndef NMSSMJACOBIAN_H
#define NMSSMJACOBIAN_H

#include "linalg.h"

namespace softsusy {

  class NmssmSoftsusy;

  class NmssmJacobian {
  public:

    explicit NmssmJacobian(NmssmSoftsusy* m);
    virtual ~NmssmJacobian();

    DoubleMatrix displayRGFlowJacobian() const { return jacRGFlow; }
    DoubleMatrix displayEWSBJacobian() const { return jacEWSB; }
    DoubleMatrix displayInverseRGFlowJacobian() const { return invJacRGFlow; }
    DoubleMatrix displayInverseEWSBJacobian() const { return invJacEWSB; }

    double calcFTJacobian(double mx, bool doTop = false);
    double calcFTInverseJacobian(double mx, bool doTop = false);
    double calcDeltaJ(double mx, bool doTop = false);

  private:
    // used for derivatives
    enum Parameters { Mzsq, Tanb, Svev, Mtsq, Lambda, Kappa,
                      SMu, M3Sq, XiS, Mh1Sq, Mh2Sq, MsSq, Yt };

    // @note have model as member (pointer/reference?) or not?
    NmssmSoftsusy* model;
    DoubleMatrix jacRGFlow;
    DoubleMatrix jacEWSB;
    DoubleMatrix invJacRGFlow;
    DoubleMatrix invJacEWSB;

    // helper functions getting pole Z and top mass
    double calcMzPole() const;
    double calcMtPole();

    struct RGFlowPars {
      NmssmSoftsusy* model;
      Parameters independent;
      Parameters dependent;
      double toScale;
    };

    static double calcRunningParameter(double x, void* parameters);
    double calcRGDerivative(Parameters dep, Parameters indep, double toScale);
    double calcRGFlowJacobian(double startScale, double endScale, bool doTop);
    double calcEWSBJacobian(bool doTop);
    double calcInverseEWSBJacobian(bool doTop);
  };

} /// namespace softsusy

#endif
