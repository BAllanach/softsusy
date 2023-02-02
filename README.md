# SOFTSUSY - supersymmetric/Higgs spectrum and decays

by: B C Allanach, P Athron, A Bednyakov, M Bernhardt, T Cridge, D Grellscheid,
                  M Hanussek, C H Kom, S Martin, D Robertson, R Ruiz de Austri,
		  P Slavich, L Tunstall, A Voigt and A G Williams 

## Summary

This program provides a SUSY spectrum in the NMSSM, or the MSSM including flavour violation and with or without R-parity consistent with input Standard Model fermion mass/mixings and electroweak/strong coupling data. The R-parity violating mode can calculate neutrino masses and mixings to 1 loop. SOFTSUSY
can be used in conjunction with other programs for [many different particle
physics calculations](https://arxiv.org/abs/0805.2088). SOFTSUSY now has a mode with 3 loop RGEs and some 2-loop threshold correction and 2-loop SUSY QCD corrections to gluino and squark pole masses. SOFTSUSY *now computes decay branching ratios for the MSSM and NMSSM*. It also ships with (and links to) Himalaya-1.0 for three-loop corrections to mh0. Legacy code can be found on the old [softsusy website](https://softsusy.hepforge.org/).

## Quick Installation and Run Test

See [INSTALL.md](INSTALL.md) for quick installation and run-test instructions.

For other ultra-basic instructions, see the [introduction video](https://www.youtube.com/watch?v=avRPn9uUKJI&ab_channel=BenAllanach). Otherwise, see a quick [tutorial](https://softsusy.hepforge.org/softsusyTutorial.pdf) given at BUSSTEPP 2012.

Note that the executables are actually wrapper scripts, the "true" executables lie in the directory `.libs/`.

## References

If you use SOFTSUSY to write a paper, *please cite* (see [MCnet guidelines](https://www.montecarlonet.org/publications_guidelines/)) - collected in [soft.bib](soft.bib)

> [1] [B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145](https://arxiv.org/abs/hep-ph/0104145)

which is the SOFTSUSY manual for the R-parity conserving MSSM.
If you use the *decay* calculations, please cite [1] and

> [2] [B.C. Allanach and T. Cridge, Comput. Phys. Comm. 220 (2017) 417, arxiv:1703.09717](https://arxiv.org/abs/1703.09717)

If you
calculate in the NMSSM, please cite [1] and 

> [3] [B.C. Allanach, P. Athron, L. Tunstall, A. Voigt and A. Williams,
Comput. Phys. Comm. 185 (2014) 2322,
arXiv:1311.7659](https://arxiv.org/abs/1311.7659)

If you use the R-parity violating aspects, please cite [1] and

> [4] [B.C. Allanach and M.A. Bernhardt, Comput. Phys. Commun. 181 (2010) 232,
arXiv:0903.1805](https://arxiv.org/abs/0903.1805)

If you use it to calculate neutrino masses and mixings, please cite [1], [4] and

> [5] [B.C. Allanach, M. Hanussek and C.H. Kom, Comput. Phys. Commun. 183 (2012)
785, arXiv:1109.3735](https://arxiv.org/abs/1109.3735)

If you use the three-loop RGEs or two-loop threshold corrections, please cite
[1] and 

> [6] [B.C. Allanach, A. Bednyakov and R. Ruiz de Autri,
Comput. Phys. Commun. 189 (2015) 192, 
arXiv:1407.6130](https://arxiv.org/abs/1407.6130)

If you use the two-loop SUSY QCD corrections to squark and gluino pole masses,
please cite [1] and 

> [7] [B.C. Allanach, Stephen P. Martin, David G. Robertson and Roberto Ruiz de
Austri, Comput. Phys. Commun. 05 (2017) 006, arXiv:1601.06657](https://arxiv.org/abs/1601.06657)

## Particle Decays

An example point including the calculation of sparticle decays, neglecting modes with a branching ratio of less than 1.0e-5, and outputting the partial widths in the comments:
```bash
./softpoint.x gmsb --n5=2 --mMess=1.0e6 --LAMBDA=5.0e5 --tanBeta=10 --sgnMu=1 --decays --minBR=1.0e-5 --outputPartialWidths
```
For queries regarding decay calculations please contact Tom Cridge

## Three Loop Corrections to the Lightest CP Even Higgs Mass

For Himalaya-1.0 three-loop corrections to mh0, you must first install the package [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page). Then do:
```bash
./configure CPPFLAGS="-I/usr/include/eigen3" --enable-two-loop-gauge-yukawa-compilation --enable-himalaya; make
```
After this, setting the `SLHA SOFTSUSY BLock` parameter 7 (number of Higgs mass loops) to 3 will include the corrections. If you use this option, you should cite [1], [6] and
* [Robert V. Harlander, Jonas Klappert, Alexander Voigt, arXiv:1708.05720](https://arxiv.org/abs/1708.05720)
* [P Kant, R Harlander, L Mihaila, M Steinhauser, JHEP 1008 (2010) 104, arXiv:1005.5709](https://arxiv.org/abs/1005.5709)


## SOFTSUSY-specific input for SUSY Les Houches Accord input files

```bash
Block SOFTSUSY           # SOFTSUSY specific inputs
  1   <TOLERANCE>        # desired fractional accuracy in output
  2   <MIXING>           # quark mixing option
  3   <PRINTOUT>         # gives additional verbose output during calculation
  4   <QEWSB>            # change electroweak symmetry breaking scale
  5   <INCLUDE_2_LOOP_SCALAR_CORRECTIONS>    # Full 2-loop running in RGEs
  6   <PRECISION>        # number of significant figures in SLHA output
  7   <numHiggsLoops>    # number of loops in REWSB/mh calculation
  8   <susyRpvBCatMSUSY> # Switch MSUSY-scale RPV boundary conditions ON
  9   <invertedOutput>   # RPV neutrino output uses normal hierarchy (=0.0) or inverted (=1.0)
 10   <forceSlha1>       # if =1, tries to force output into SLHA *1* format
 11   <m32>              # sets gravitino mass to m32
 12   <printSpectra>     # Prints spectrum even when point is theoretically excluded if=1
 13   <mAFlag>           # If=0 (default), sets tachyonic mA=0, otherwise mA=sqrt(|mA|^2)
 15   <NMSSMTools>       # If=1, enables NMSSMTools compatible SLHA2 output
 16   <MICROMEGAS>       # Micromegas options for NMSSMTools use: 1=RD, 2=DD, 3=ID, 4=both
 17   <NMSDECAY>         # If=1, flags for sparticle decays to be calculated via NMSDECAY
 18   <SoftHiggsOut>     # If=1, then the EWSB conditions output soft Higgs masses in NMSSM
 19   <threeLoopRGEs>    # If=1, then 3-loop MSSM RGEs included (default of 0 to disable)
 20   <gyThresholds>     # If>0, switch on gauge/Yukawa two-loop thresholds (see manual [6] for details). 
                           If=31, they all are switched on (default 0 to disable).
 22    <2-loop squark/gluino>   # Include 2-loop terms in gluino/squark masses (default of 0 to disable)
 23    <expandAroundGluinoPole> # sets expandAroundGluinoPole parameter (default 3)
24   <minBR>                  # If decay BR is below this number, don't output that mode
 25   <threeBodyDecays>        # If set to 0, don't calculate 3-body decays (1=default)
 26   <outputPartialWidths>    # If set to 1, output partial widths (0=default)
 27   <MQEDxQCD>    # Set scale at which QEDxQCD is matched to MSSM (mt=default)```

## High orders mode

If the 2-loop SUSY QCD corrections to squark and gluino masses are required, do

```bash
./configure --enable-two-loop-sparticle-mass-compilation
make 
```
An example point using the higher order terms can be run with, for example,
```bash
./softpoint.x sugra --tol=1.0e-4 --m0=1000 --m12=1000 --a0=0 --tanBeta=10 --sgnMu=1 --two-loop-sparticle-masses --two-loop-sparticle-mass-method=1
```

## High accuracy mode

If the high accuracy mode with 3-loop RGEs and some 2-loop threshold corrections is required, do
```bash
./configure --enable-full-susy-threshold-compilation --enable-three-loop-rge-compilation
make 
```

An example point using the high accuracy mode can be run with, for example,

```bash
./softpoint.x sugra --tol=1.0e-5 --m0=7240 --m12=800 --a0=-6000 --tanBeta=50 --sgnMu=1 --mt=173.2 --alpha_s=0.1187 --mbmb=4.18 --two-loop-susy-thresholds --three-loop-rges
```

See [6] for more details.

## Comparisons with other SUSY spectrum generators

There are detailed comparisons between SOFTSUSY and other publicly available
codes in 
* Uncertainties in the Lightest CP Even Higgs Boson Mass Prediction in the Minimal Supersymmetric Standard Model: Fixed Order Versus Effective Field Theory Prediction, [B.C. Allanach and A. Voigt, Eur.Phys.J. C78 (2018) no.7, arxiv:1804.09410](https://arxiv.org/abs/1804.09410)
* The Calculation of Sparticle and Higgs Decays in the Minimal and Next-to-Minimal Supersymmetric Standard Models: SOFTSUSY4.0, [B.C. Allanach and T. Cridge, Comput. Phys. Comm. 220 (2017) 417, arxiv:1703.09717](https://arxiv.org/abs/1703.09717)
* Precise Determination of the Neutral Higgs Boson Masses in the MSSM, [B.C. Allanach, A. Djouadi, J.L. Kneur, W. Porod, P. Slavich, JHEP 0409 (2004) 044, hep-ph/0406166](https://arxiv.org/abs/hep-ph/0406166) 
* Theoretical uncertainties in sparticle mass predictions from computational tools, [B.C. Allanach, S. Kraml, W. Porod, JHEP 03 (2003) 045, hep-ph/0302102](https://arxiv.org/abs/hep-ph/0302102)  

and comparisons with NMSSM generators in
* [B.C. Allanach and T. Cridge, Comput. Phys. Comm. 220 (2017) 417, arxiv:1703.09717](https://arxiv.org/abs/1703.09717)
* Higgs mass predictions of public NMSSM spectrum generators, [Staub et al, Comp. Phys. Comm. 202 (2016) 113, arXiv:1507.05093](https://arxiv.org/abs/1507.05093)
* Improved predictions for intermediate and heavy Supersymmetry in the MSSM and beyond [Staub and Porod, Eur. Phys. J. C (2017) 77, arXiv:1703.03267](https://arxiv.org/abs/1703.03267)


## Executable files: after installation

* `softpoint.x`: command-line interface. GMSB, AMSB, mSUGRA and general boundary conditions possible, icluding SLHA. Main program: `src/softpoint.cpp`
* `softsusy.x`: example C++ test program - calculates spectrum of SPS1a mSUGRA point with varying tan beta. Main program: `src/main.cpp`
* `softsusy-nmssm.x`: example NMSSM test program - loops over tan beta. Main program: `src/main-nmssm.cpp`
* `rpvsoftsusy.x`: example C++ test program - calculates spectrum of SPS1a mSUGRA point with varying lambda'_{331}(M_GUT). Main program: `src/rpvmain.cpp`
* `rpvneut.x`: example neutrino mass calculating R-parity violating test program. Main program `src/rpvNeut.cpp`

## Documentation

Full code documentation can be obtained from: [softsusy](http://softsusy.hepforge.org/). Manuals are found in the `doc/` subdirectory: `rpcManual.pdf` is the main one, and the other `.pdf` files describe various extensions to this.

## Input and information files
* `README.md` contains these instructions 
* `inOutFiles/outputTest` is the output from the test program
* `inOutFiles/slha2Input` is an alternative input file in the SUSY Les Houches Accord 2 format for SPS1a' 
* `inOutFiles/slha2Output` is the result of running with the above input file and includes flavour violation, for inclusion into codes like SusyBsg1.3 which include flavour corrections     
* `inOutFiles/lesHouchesInput` is an alternative input file in the SUSY Les Houches Accord format
* `inOutFiles/rpvHouchesInput` is an alternative input file in the SUSY Les Houches Accord format for R-parity violation
* `inOutFiles/rpvHouchesOutput` is the output from the R-parity violating test program rpvmain.cpp

## Files included in this distribution

Source files are to be found in the `src/` subdirectory. The `inOutFiles/` directory contains input and output files. `doc/` contains the manuals (see above).

## Licence

    SOFTSUSY Copyright (C) 2007 B.C. Allanach

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    See [licenses](http://www.gnu.org/licenses/)  

