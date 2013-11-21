bin_PROGRAMS = softsusy.x softpoint.x rpvsoftsusy.x rpvneut.x gaugeUn.x Devtest.x TestNM.x doc/nmssmManual.pdf doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf

lib_LTLIBRARIES = libsoft.la 

softsusy_x_SOURCES = main.cpp
softsusy_x_LDADD = libsoft.la  

softpoint_x_SOURCES = softpoint.cpp 
softpoint_x_LDADD = libsoft.la

rpvsoftsusy_x_SOURCES = rpvmain.cpp
rpvsoftsusy_x_LDADD = libsoft.la 

rpvneut_x_SOURCES = rpvneutmain.cpp
rpvneut_x_LDADD = libsoft.la 

gaugeUn_x_SOURCES = gaugeUn.cpp
gaugeUn_x_LDADD = libsoft.la 

Devtest_x_SOURCES = mainDev.cpp
Devtest_x_LDADD = libsoft.la

TestNM_x_SOURCES = mainTestNM.cpp
TestNM_x_LDADD = libsoft.la 

pkginclude_HEADERS = \
def.h main.h numerics.h \
softpars.h susy.h linalg.h \
physpars.h softpoint.h twoloophiggs.h \
nmssm2loop.h dilogwrap.h \
lowe.h mycomplex.h rge.h rgeroutines.h softsusy.h \
utils.h mssmUtils.h rpvsoft.h rpvsusypars.h \
tensor.h rpvmain.h flavoursoft.h \
xpr-base.h xpr-matrix.h xpr-vector.h rpvneut.h \
nmssmsusy.h nmssmsoftpars.h  \
nmssmsoftsusy.h nmssmUtils.h

libsoft_la_SOURCES = def.cpp linalg.cpp lowe.cpp numerics.cpp  rpvneut.cpp \
physpars.cpp rge.cpp rgeroutines.cpp nmssm2loop.f dilogwrap.f \
 nmssmsoftsusy.cpp nmssmsusy.cpp nmssmsoftpars.cpp nmssmUtils.cpp \
susy.cpp utils.cpp mssmUtils.cpp \
tensor.cpp rpvsusypars.cpp rpvsoft.cpp flavoursoft.cpp twoloophiggs.f 
libsoft_la_LDFLAGS = $(FLIBS)

DOCFILES=doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf \
doc/nmssmManual.pdf
EXTRA_DIST=$(DOCFILES)

OUTPUTFILES = slha2Output lesHouchesOutput rpvHouchesOutput \
 nmssmSLHAZ3Output nmssmSLHAnoZ3Output
TESTOUTPUTFILES = .test/slha2Output .test/rpvHouchesOutput .test/lesHouchesOutput
INPUTFILES = slha2Input lesHouchesInput rpvHouchesInput

test: softpoint.x 
	git checkout $(OUTPUTFILES)
	cp $(OUTPUTFILES) .test; rm $(OUTPUTFILES); make $(OUTPUTFILES); \
	ndiff --relative-error 1.0e-3 slha2Output .test/slha2Output; \
	ndiff --relative-error 1.0e-3 lesHouchesOutput .test/lesHouchesOutput; \
	ndiff --relative-error 1.0e-3 rpvHouchesOutput .test/rpvHouchesOutput; \
	ndiff --relative-error 1.0e-3 nmssmSLHAZ3Output .test/nmssmSLHAZ3Output; \
	ndiff --relative-error 1.0e-3 nmssmSLHAnoZ3Output .test/nmssmSLHAnoZ3Output

%Output : %Input softpoint.x 
	./softpoint.x leshouches < $< > $@

outputTest: softsusy.x
	./softsusy.x > outputTest

rpvOutputTest: rpvsoftsusy.x
	./rpvsoftsusy.x > rpvOutputTest

neutrinoOutputTest: rpvneut.x
	./rpvneut.x > neutrinoOutputTest

gaugeUn.out: gaugeUn.x
	./gaugeUn.x > gaugeUn.out 


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
