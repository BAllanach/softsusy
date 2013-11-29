SUBDIRS = src

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = softsusy.x softsusy-nmssm.x softpoint.x rpvsoftsusy.x rpvneut.x gaugeUn.x doc/nmssmManual.pdf doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf

softsusy_x_SOURCES = src/main.cpp
softsusy_x_LDADD = src/libsoft.la  

softsusy_nmssm_x_SOURCES = src/main-nmssm.cpp
softsusy_nmssm_x_LDADD = src/libsoft.la

softpoint_x_SOURCES = src/softpoint.cpp 
softpoint_x_LDADD = src/libsoft.la

rpvsoftsusy_x_SOURCES = src/rpvmain.cpp
rpvsoftsusy_x_LDADD = src/libsoft.la 

rpvneut_x_SOURCES = src/rpvneutmain.cpp
rpvneut_x_LDADD = src/libsoft.la 

gaugeUn_x_SOURCES = src/gaugeUn.cpp
gaugeUn_x_LDADD = src/libsoft.la 

DOCFILES=doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf \
doc/nmssmManual.pdf
EXTRA_DIST=$(DOCFILES)

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

doc: doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf \
	doc/nmssmManual.pdf

doc/rpcManual.pdf: doc/rpcManual.tex
	cd doc; \
	latex rpcManual.tex; latex rpcManual.tex; \
	dvips -f rpcManual > rpcManual.ps; 	\
	ps2pdf rpcManual.ps; cd ..

doc/rpvManual.pdf: doc/rpvManual.tex
	cd doc; \
	latex rpvManual.tex; latex rpvManual.tex; \
	dvips -f rpvManual > rpvManual.ps; 	\
	ps2pdf rpvManual.ps; cd ..

doc/neutManual.pdf: doc/neutManual.tex
	cd doc; \
	latex neutManual.tex; latex neutManual.tex; \
	dvips -f neutManual > neutManual.ps; 	\
	ps2pdf neutManual.ps; cd ..

doc/nmssmManual.pdf: doc/nmssmManual.tex
	cd doc; \
	latex nmssmManual.tex; latex nmssmManual.tex; \
	dvips -f nmssmManual > nmssmManual.ps; 	\
	ps2pdf nmssmManual.ps; cd ..
