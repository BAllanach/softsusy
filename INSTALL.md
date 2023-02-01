# Installation and quick run test

The following releases contain a test program (`main.cpp`) and the SOFTSUSY library (`libsoft.a,` link with `-L.libs -lsoft`). In linux, just unpack the files with (eg for `softsusy-4.0`)
```bash
gunzip softsusy-4.0.tar.gz
tar -xvf softsusy-4.0.tar 
cd softsusy-4.0
```

Then, for simplest installation, to compile the code:
```bash
./configure
make programs
```

To run SOFTSUSY, you should need only standard `C++` and `fortran` libraries.

There are four C++ test programs, which can be run by the commands 
```bash
./softsusy.x
./rpvsoftsusy.x 
./rpvneut.x
./softsusy-nmssm.x
```
The output from these commands can be checked against `outputTest`,
`rpvOutputTest`, `neutOutputTest` and `outputTest-nmssm`.

You can run the SUSY Les Houches Accord input provided by running the commands
```bash
./softpoint.x leshouches < inOutFiles/lesHouchesInput > inOutFiles/lesHouchesOutput
./softpoint.x leshouches < inOutFiles/nmssmSLHAnoZ3Input > inOutFiles/nmssmSLHAnoZ3Output
./softpoint.x leshouches < inOutFiles/nmssmSLHAZ3Input > inOutFiles/nmssmSLHAZ3Output
./softpoint.x leshouches < inOutFiles/rpvHouchesInput > inOutFiles/rpvHouchesOutput
./softpoint.x leshouches < inOutFiles/slha2Input > inOutFiles/slha2Output
```
You may check the output of these commands against the output files
in directory `inOutFiles/`.

All of the output files mentioned above are produced by the `Makefile` automatically.
*SOFTSUSY executables use no input or output files except for standard input or standard output.*
