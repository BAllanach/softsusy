# SOFTSUSY - supersymmetric/Higgs spectrum and decays

by: [B C Allanach](https://www.damtp.cam.ac.uk/user/bca20/), P Athron, A Bednyakov, M Bernhardt, T Cridge, D Grellscheid,
                  M Hanussek, C H Kom, S Martin, D Robertson, R Ruiz de Austri,
		  P Slavich, L Tunstall, A Voigt and A G Williams 

## Summary

This program provides a SUSY spectrum in the NMSSM, or the MSSM including flavour violation and with or without R-parity consistent with input Standard Model fermion mass/mixings and electroweak/strong coupling data. The R-parity violating mode can calculate neutrino masses and mixings to 1 loop. SOFTSUSY
can be used in conjunction with other programs for [many different particle
physics calculations](https://arxiv.org/abs/0805.2088). SOFTSUSY now has a mode with 3 loop RGEs and some 2-loop threshold correction and 2-loop SUSY QCD corrections to gluino and squark pole masses. SOFTSUSY *now computes decay branching ratios for the MSSM and NMSSM*. It also ships with (and links to) Himalaya-1.0 for three-loop corrections to mh0. 

## Download
* <a id="raw-url" href="https://raw.githubusercontent.com/BAllanach/softsusy/master/tags/softsusy-4.1.22.tar.gz">SOFTSUSY-4.1.22</a> (07/08/25): bug-fixed neutralino_i -> chargino_i three-body decays: added missing CP-conjugate decay modes. Fixed naming of antiparticles and inos in comments of decay file to be more universal. Thanks to K Zhang for reporting the bug and finding the code location of it.
* <a id="raw-url" href="https://raw.githubusercontent.com/BAllanach/softsusy/master/tags/softsusy-4.1.20.tar.gz">SOFTSUSY-4.1.20</a> (23/8/23): Fixed SLHA output which was giving the wrong RPV inputs. SOFTSUSY-specific SLHA input options given automatically printed out now in SLHA output. Fixed dm^2(atm) in neutrino output for RPV mode. Also, corrected RPV pMSSM SLHA output, adding rpvPmssmInput to repository.
* <a id="raw-url" href="https://raw.githubusercontent.com/BAllanach/softsusy/master/tags/softsusy-4.1.19.tar.gz">SOFTSUSY-4.1.19</a> (16/8/23): improved 2/3 body neutralino/chargino decays. Previously, there were occasionally numerical problems yielding junk for a (very heavy) RH sneutrino contribution. The RH sneutrino parts have now been explicitly decoupled. Thanks to S Kraml and T Pascal for reporting the bug.
* <a id="raw-url" href="https://raw.githubusercontent.com/BAllanach/softsusy/master/tags/softsusy-4.1.18.tar.gz">SOFTSUSY-4.1.18</a> (25/7/23): improved 2/3 body neutralino decays. Occasionally, SOFTSUSY would find neutralino2 to be stable when in fact it really had 3 body decays. This occured when the decay of the running mass of the Z would allow 2 body decays but not the pole mass. Thanks to S Kraml and T Pascal for reporting the bug.
* <a id="raw-url" href="https://raw.githubusercontent.com/BAllanach/softsusy/master/tags/softsusy-4.1.17.tar.gz">SOFTSUSY-4.1.17</a> (2/2/23): improved documentation for github

Previous releases (4.1.13 and before) can be obtained from the [old website](https://softsusy.hepforge.org/)

## Quick Installation and Run Test

See [INSTALL.md](INSTALL.md) for quick installation and run-test instructions.

For other ultra-basic instructions, see the [introduction video](https://www.youtube.com/watch?v=avRPn9uUKJI&ab_channel=BenAllanach). Otherwise, see a quick [tutorial](https://softsusy.hepforge.org/softsusyTutorial.pdf) given at BUSSTEPP 2012.

Note that the executables are actually wrapper scripts, the "true" executables lie in the directory `.libs/`.

## References/manuals

If you use SOFTSUSY to write a paper, *please cite* (see [MCnet guidelines](https://www.montecarlonet.org/publications_guidelines/)) - collected in [soft.bib](soft.bib) the manuals (up-to-date manuals included with distribution in `doc` subdirectory):

> 1. [`rpcManual.pdf`](doc/rpcManual.pdf) is the *main* one for the R-parity conserving MSSM: the base for all others: [B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145](https://arxiv.org/abs/hep-ph/0104145) 
> 2. [`decays.pdf`](doc/rpvManual.pdf) details calculations of sparticle and Higgs decays in the MSSM/NMSSM: [B.C. Allanach and T. Cridge, Comput. Phys. Comm. 220 (2017) 417, arxiv:1703.09717](https://arxiv.org/abs/1703.09717) 
> 3. [`nmssmManual.pdf`](doc/nmssmManual.pdf) describes the NMSSM implementation: [B.C. Allanach, P. Athron, L. Tunstall, A. Voigt and A. Williams, Comput. Phys. Comm. 185 (2014) 2322, arXiv:1311.7659](https://arxiv.org/abs/1311.7659) 
> 4. [`rpvManual.pdf`](doc/rpvManual.pdf) for R-parity violating generalisation: [B.C. Allanach and M.A. Bernhardt, Comput. Phys. Commun. 181 (2010) 232,
arXiv:0903.1805](https://arxiv.org/abs/0903.1805)
> 5. [`neutManual.pdf`](doc/neutManual.pdf) on the one-loop calculation of neutrino masses and lepton mixing in the R-parity violating MSSM: [B.C. Allanach, M. Hanussek and C.H. Kom, Comput. Phys. Commun. 183 (2012) 785, arXiv:1109.3735](https://arxiv.org/abs/1109.3735)
> 6. [`threeLoop.pdf`](doc/threeLoop.pdf) describes the inclusion of three-loop MSSM RGEs and two-loop threshold corrections: [B.C. Allanach, A. Bednyakov and R. Ruiz de Autri, Comput. Phys. Commun. 189 (2015) 192, arXiv:1407.6130](https://arxiv.org/abs/1407.6130)
> 7. [`ho.pdf`](doc/ho.pdf) describes the inclusion of two-loop SUSYQCD corrections to squark and gluino pole masses: [B.C. Allanach, Stephen P. Martin, David G. Robertson and Roberto Ruiz de Austri, Comput. Phys. Commun. 05 (2017) 006, arXiv:1601.06657](https://arxiv.org/abs/1601.06657)

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
* [Robert V. Harlander, Jonas Klappert, Alexander Voigt, Eur. Phys. J. C77 (2017) 814, arXiv:1708.05720](https://arxiv.org/abs/1708.05720)
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
 20   <gyThresholds>     # If>0, switch on gauge/Yukawa two-loop thresholds (see manual [6] for details). If=31, they all are switched on (default 0 to disable).
 22   <2-loop squark/gluino>   # Include 2-loop terms in gluino/squark masses (default of 0 to disable)
 23   <expandAroundGluinoPole> # sets expandAroundGluinoPole parameter (default 3)
 24   <minBR>                  # If decay BR is below this number, don't output that mode
 25   <threeBodyDecays>        # If set to 0, don't calculate 3-body decays (1=default)
 26   <outputPartialWidths>    # If set to 1, output partial widths (0=default)
 27   <MQEDxQCD>    # Set scale at which QEDxQCD is matched to MSSM (mt=default)
```

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

## Executable files: after installation

* `softsusy.x`: example C++ test program - calculates spectrum of SPS1a mSUGRA point with varying tan beta. Main program: [`src/main.cpp`](src/main.cpp). Output: [`inOutFiles/outputTest`](inOutFiles/outputTest)
* `softsusy-nmssm.x`: example NMSSM test program - loops over tan beta. Main program: [`src/main-nmssm.cpp`](src/main-nmssm.cpp). Output: [`inOutFiles/outputTest-nmssm`](inOutFiles/outputTest-nmssm)
* `rpvsoftsusy.x`: example C++ test program - calculates spectrum of SPS1a mSUGRA point with varying lambda'_{331}(M_GUT). Main program: [`src/rpvmain.cpp`](src/rpvmain.cpp). Output [`inOutFiles/rpvOutputTest`](inOutFiles/rpvOutputtest)
* `rpvneut.x`: example neutrino mass calculating R-parity violating test program. Main program [`src/rpvNeut.cpp`](src/rpvNeut.cpp). Output [`inOutFiles/neutOutputTest`](inOutFiles/neutOuptutTest)
* `softpoint.x`: command-line interface. GMSB, AMSB, mSUGRA and general boundary conditions possible, icluding SLHA. Main program: [`src/softpoint.cpp`](src/softpoint.cpp). See [INSTALL.md](INSTALL.md) for examples of reading in SLHA files and producing output.

### Input and information files
* `README.md` contains these instructions 
* [`inOutFiles/lesHouchesInput`](inOutFileslesHouchesInput) is an alternative input file in the [SUSY Les Houches Accord](https://arxiv.org/abs/hep-ph/0311123) (SLHA) format
* [`inOutFiles/nmssmSLHAnoZ3Input`](inOutFiles/nmssmSLHAnoZ3Input) is an [SLHA2](https://arxiv.org/abs/0801.0045) NMSSM input file *without* assuming Z3 symmetry
* [`inOutFiles/nmssmSLHAZ3Input`](inOutFiles/nmssmSLHAZ3Input) is an SLHA NMSSM input file with the Z3 assumption
* [`inOutFiles/slha2Input`](inOutFiles/slha2Input) is an alternative input file in the SUSY Les Houches Accord 2 format for SPS1a' 
* [`inOutFiles/rpvHouchesInput`](inOutFiles/rpvHouchesInput) is an alternative input file in the SUSY Les Houches Accord format for R-parity violation


## Files included in this distribution

Source files are to be found in the `src/` subdirectory. The `inOutFiles/` directory contains input and output files. `doc/` contains the manuals (see above).

## Comparisons with other SUSY spectrum generators

There are detailed comparisons between SOFTSUSY and other publicly available programs in 
* Uncertainties in the Lightest CP Even Higgs Boson Mass Prediction in the Minimal Supersymmetric Standard Model: Fixed Order Versus Effective Field Theory Prediction, [B.C. Allanach and A. Voigt, Eur.Phys.J. C78 (2018) no.7, arxiv:1804.09410](https://arxiv.org/abs/1804.09410)
* The Calculation of Sparticle and Higgs Decays in the Minimal and Next-to-Minimal Supersymmetric Standard Models: SOFTSUSY4.0, [B.C. Allanach and T. Cridge, Comput. Phys. Comm. 220 (2017) 417, arxiv:1703.09717](https://arxiv.org/abs/1703.09717)
* Precise Determination of the Neutral Higgs Boson Masses in the MSSM, [B.C. Allanach, A. Djouadi, J.L. Kneur, W. Porod, P. Slavich, JHEP 0409 (2004) 044, hep-ph/0406166](https://arxiv.org/abs/hep-ph/0406166) 
* Theoretical uncertainties in sparticle mass predictions from computational tools, [B.C. Allanach, S. Kraml, W. Porod, JHEP 03 (2003) 045, hep-ph/0302102](https://arxiv.org/abs/hep-ph/0302102)  

and comparisons with NMSSM generators in
* [B.C. Allanach and T. Cridge, Comput. Phys. Comm. 220 (2017) 417, arxiv:1703.09717](https://arxiv.org/abs/1703.09717)
* Higgs mass predictions of public NMSSM spectrum generators, [Staub et al, Comp. Phys. Comm. 202 (2016) 113, arXiv:1507.05093](https://arxiv.org/abs/1507.05093)
* Improved predictions for intermediate and heavy Supersymmetry in the MSSM and beyond [Staub and Porod, Eur. Phys. J. C (2017) 77, arXiv:1703.03267](https://arxiv.org/abs/1703.03267)


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

