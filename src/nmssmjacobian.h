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
    static double calcMzPole(NmssmSoftsusy*);
    static double calcMtPole(NmssmSoftsusy*);

    struct EWSBPars {
      NmssmSoftsusy* model;
      Parameters independent;
      Parameters dependent;
      DoubleVector outputs;

      EWSBPars()
        : model(0), independent(Mzsq), dependent(Lambda), outputs(3)
        {}
    };

    struct RGFlowPars {
      NmssmSoftsusy* model;
      Parameters independent;
      Parameters dependent;
      double toScale;
    };

    static double calcRunningParameter(double x, void* parameters);
    double calcRGDerivative(Parameters dep, Parameters indep, double toScale);
    double calcRGFlowJacobian(double startScale, double endScale, bool doTop);
    static int ewsbOutputErrors(const DoubleVector & guess, void* parameters,
                                DoubleVector & errors);
    static void fixEWSBOutputs(EWSBPars* pars, int & err);
    static double calcEWSBParameter(double x, void* parameters);
    static double calcEWSBOutput(double x, void* parameters);
    double calcEWSBOutputDerivative(Parameters dep, Parameters indep);
    double calcEWSBParameterDerivative(Parameters dep, Parameters indep, bool doTop);
    double calcEWSBJacobian(bool doTop);
    double calcInverseEWSBJacobian(bool doTop);
  };

} /// namespace softsusy

#endif
