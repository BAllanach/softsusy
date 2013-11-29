SUBDIRS = src doc

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

OUTPUTFILES = inOutFiles/slha2Output inOutFiles/lesHouchesOutput \
inOutFiles/nmssmSLHAZ3Output \
inOutFiles/nmssmSLHAnoZ3Output
TESTOUTPUTFILES = .test/slha2Output .test/rpvHouchesOutput \
.test/lesHouchesOutput .test/nmssmSLHAZ3Output .test/nmssmSLHAnoZ3Output
INPUTFILES = inOutfiles/slha2Input inOutFiles/lesHouchesInput \
inOutFiles/rpvHouchesInput inOutFiles/nmssmSLHAZ3Input \
inOutFiles/nmssmSLHAnoZ3Input

test: softpoint.x 
	git checkout $(OUTPUTFILES)
	cp $(OUTPUTFILES) .test; rm $(OUTPUTFILES); make $(OUTPUTFILES); \
	ndiff --relative-error 1.0e-3 inOutFiles/slha2Output .test/slha2Output; \
	ndiff --relative-error 1.0e-3 inOutFiles/lesHouchesOutput .test/lesHouchesOutput; \
	ndiff --relative-error 1.0e-3 inOutFiles/rpvHouchesOutput .test/rpvHouchesOutput; \
	ndiff --relative-error 1.0e-3 inOutFiles/nmssmSLHAZ3Output .test/nmssmSLHAZ3Output; \
	ndiff --relative-error 1.0e-3 inOutFiles/nmssmSLHAnoZ3Output .test/nmssmSLHAnoZ3Output

inOutFiles/%Output : inOutFiles/%Input softpoint.x 
	./softpoint.x leshouches < $< > $@

inOutFiles/outputTest: softsusy.x
	./softsusy.x > inOutFiles/outputTest

inOutFiles/rpvOutputTest: rpvsoftsusy.x
	./rpvsoftsusy.x > inOutFiles/rpvOutputTest

inOutFiles/neutOutputTest: rpvneut.x
	./rpvneut.x > inOutFiles/neutOutputTest

inOutFiles/gaugeUn.out: gaugeUn.x
	./gaugeUn.x > inOutFiles/gaugeUn.out 
