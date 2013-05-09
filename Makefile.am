bin_PROGRAMS = softsusy.x softpoint.x rpvsoftsusy.x rpvneut.x gaugeUn.x doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf runnmssm.x nmssmsoftsusy.x .test/test_nmssm_betas.x .test/test_matrix_resize.x

lib_LTLIBRARIES = libsoft.la libnmsoft.la .test/libsarah.la

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

nmssmsoftsusy_x_SOURCES = main-nmssm.cpp
nmssmsoftsusy_x_LDADD = libnmsoft.la

runnmssm_x_SOURCES = runnmssm.cpp
runnmssm_x_LDADD = libnmsoft.la

_test_test_nmssm_betas_x_SOURCES = .test/test_nmssm_betas.cpp
_test_test_nmssm_betas_x_LDADD = .test/libsarah.la

_test_test_matrix_resize_x_SOURCES = .test/test_matrix_resize.cpp
_test_test_matrix_resize_x_LDADD = libsoft.la

pkginclude_HEADERS = \
def.h main.h numerics.h \
softpars.h susy.h linalg.h \
physpars.h softpoint.h twoloophiggs.h \
lowe.h mycomplex.h rge.h rgeroutines.h softsusy.h \
utils.h rpvsoft.h rpvsusypars.h \
tensor.h rpvmain.h flavoursoft.h \
xpr-base.h xpr-matrix.h xpr-vector.h rpvneut.h \
nmssmsusy.h nmssmsoftpars.h nmssmphyspars.h nmssmsoftsusy.h

libsoft_la_SOURCES = def.cpp linalg.cpp lowe.cpp numerics.cpp  rpvneut.cpp \
physpars.cpp rge.cpp rgeroutines.cpp susy.cpp utils.cpp bcs.cpp mssmUtils.cpp \
tensor.cpp rpvsusypars.cpp rpvsoft.cpp flavoursoft.cpp twoloophiggs.f
libsoft_la_LDFLAGS = $(FLIBS)

libnmsoft_la_SOURCES = def.cpp linalg.cpp lowe.cpp numerics.cpp \
nmssmphyspars.cpp rge.cpp nmssmsoftsusy.cpp nmssmsusy.cpp nmssmsoftpars.cpp \
 utils.cpp twoloophiggs.f tensor.cpp susy.cpp
libnmsoft_la_LDFLAGS = $(FLIBS)

_test_libsarah_la_SOURCES = def.cpp linalg.cpp lowe.cpp numerics.cpp \
nmssmphyspars.cpp rge.cpp nmssmsoftsusy.cpp nmssmsusy.cpp nmssmsoftpars.cpp \
utils.cpp twoloophiggs.f tensor.cpp susy.cpp \
.test/sarah_nmssm_susy.cpp .test/sarah_nmssm_softpars.cpp \
.test/mathematica_wrappers.cpp
_test_libsarah_la_LDFLAGS = $(FLIBS)

DOCFILES=doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf
EXTRA_DIST=$(DOCFILES)

OUTPUTFILES = slha2Output lesHouchesOutput rpvHouchesOutput 
TESTOUTPUTFILES = .test/slha2Output .test/rpvHouchesOutput .test/lesHouchesOutput
INPUTFILES = slha2Input lesHouchesInput rpvHouchesInput

test: softpoint.x 
	git checkout $(OUTPUTFILES)
	cp $(OUTPUTFILES) .test; rm $(OUTPUTFILES); make $(OUTPUTFILES); \
	ndiff --relative-error 1.0e-3 slha2Output .test/slha2Output; \
	ndiff --relative-error 1.0e-3 lesHouchesOutput .test/lesHouchesOutput; \
	ndiff --relative-error 1.0e-3 rpvHouchesOutput .test/rpvHouchesOutput

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

doc: doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf

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


