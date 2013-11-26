bin_PROGRAMS = softsusy.x softsusy-nmssm.x softpoint.x rpvsoftsusy.x rpvneut.x 

lib_LTLIBRARIES = libsoft.la 

softsusy_x_SOURCES = main.cpp
softsusy_x_LDADD = libsoft.la  

softsusy_nmssm_x_SOURCES = main-nmssm.cpp
softsusy_nmssm_x_LDADD = libsoft.la

softpoint_x_SOURCES = softpoint.cpp 
softpoint_x_LDADD = libsoft.la

rpvsoftsusy_x_SOURCES = rpvmain.cpp
rpvsoftsusy_x_LDADD = libsoft.la 

rpvneut_x_SOURCES = rpvneutmain.cpp
rpvneut_x_LDADD = libsoft.la 

pkginclude_HEADERS = \
def.h main.h numerics.h \
softpars.h susy.h linalg.h \
physpars.h softpoint.h twoloophiggs.h \
nmssm2loop.h dilogwrap.h \
lowe.h mycomplex.h rge.h softsusy.h \
utils.h mssmUtils.h rpvsoft.h rpvsusypars.h \
tensor.h rpvmain.h flavoursoft.h \
xpr-base.h xpr-matrix.h xpr-vector.h rpvneut.h \
nmssmsusy.h nmssmsoftpars.h  \
nmssmsoftsusy.h nmssmUtils.h

libsoft_la_SOURCES = def.cpp linalg.cpp lowe.cpp numerics.cpp  rpvneut.cpp \
physpars.cpp rge.cpp nmssm2loop.f dilogwrap.f \
 nmssmsoftsusy.cpp nmssmsusy.cpp nmssmsoftpars.cpp nmssmUtils.cpp \
susy.cpp utils.cpp mssmUtils.cpp \
tensor.cpp rpvsusypars.cpp rpvsoft.cpp flavoursoft.cpp twoloophiggs.f 
libsoft_la_LDFLAGS = $(FLIBS)

DOCFILES=doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf \
doc/nmssmManual.pdf
EXTRA_DIST=$(DOCFILES)

OUTPUTFILES = slha2Output lesHouchesOutput rpvHouchesOutput \
 nmssmSLHAZ3Output nmssmSLHAnoZ3Output neutOutputTest \
benchmarks/cmssm10.1.1.spec benchmarks/cmssm10.1.2.spec \
benchmarks/cmssm10.1.3.spec benchmarks/cmssm10.1.4.spec \
benchmarks/cmssm10.1.5.spec \
benchmarks/cmssm10.2.1.spec benchmarks/cmssm10.2.2.spec \
benchmarks/cmssm10.2.3.spec benchmarks/cmssm10.2.4.spec \
benchmarks/cmssm10.2.5.spec \
benchmarks/cmssm10.3.1.spec benchmarks/cmssm10.3.2.spec \
benchmarks/cmssm10.3.3.spec benchmarks/cmssm10.3.4.spec \
benchmarks/cmssm10.3.5.spec \
benchmarks/cmssm10.4.1.spec benchmarks/cmssm10.4.2.spec \
benchmarks/cmssm10.4.3.spec benchmarks/cmssm10.4.4.spec \
benchmarks/cmssm10.4.5.spec \
benchmarks/cmssm40.1.1.spec benchmarks/cmssm40.1.2.spec \
benchmarks/cmssm40.1.3.spec benchmarks/cmssm40.1.4.spec \
benchmarks/cmssm40.1.5.spec \
benchmarks/cmssm40.2.1.spec benchmarks/cmssm40.2.2.spec \
benchmarks/cmssm40.2.3.spec benchmarks/cmssm40.2.4.spec \
benchmarks/cmssm40.2.5.spec \
benchmarks/cmssm40.3.1.spec benchmarks/cmssm40.3.2.spec \
benchmarks/cmssm40.3.3.spec benchmarks/cmssm40.3.4.spec \
benchmarks/cmssm40.3.5.spec \
benchmarks/mamsb1.1.spec benchmarks/mamsb1.2.spec benchmarks/mamsb1.3.spec \
benchmarks/mgmsb1.1.spec benchmarks/mgmsb1.2.spec benchmarks/mgmsb1.3.spec \
benchmarks/mgmsb1.4.spec benchmarks/mgmsb1.5.spec benchmarks/mgmsb1.6.spec \
benchmarks/mgmsb2.1.1.spec benchmarks/mgmsb2.1.2.spec benchmarks/mgmsb2.1.3.spec \
benchmarks/mgmsb2.1.4.spec benchmarks/mgmsb2.1.5.spec benchmarks/mgmsb2.1.6.spec \
benchmarks/mgmsb2.2.1.spec benchmarks/mgmsb2.2.2.spec benchmarks/mgmsb2.2.3.spec \
benchmarks/mgmsb2.2.4.spec benchmarks/mgmsb2.2.5.spec benchmarks/mgmsb2.2.6.spec \
benchmarks/pmssm1.1.spec benchmarks/pmssm1.2.spec benchmarks/pmssm1.3.spec \
benchmarks/pmssm1.4.spec benchmarks/pmssm1.5.spec benchmarks/pmssm1.6.spec \
benchmarks/pmssm2.1.spec benchmarks/pmssm2.2.spec benchmarks/pmssm2.3.spec \
benchmarks/pmssm2.4.spec benchmarks/pmssm2.5.spec benchmarks/pmssm2.6.spec \
benchmarks/rpv1.1.spec benchmarks/rpv1.2.spec benchmarks/rpv1.3.spec \
benchmarks/rpv1.4.spec benchmarks/rpv1.5.spec \
benchmarks/rpv2.1.spec benchmarks/rpv2.2.spec benchmarks/rpv2.3.spec \
benchmarks/rpv2.4.spec benchmarks/rpv2.5.spec \
benchmarks/rpv3.1.spec benchmarks/rpv3.2.spec benchmarks/rpv3.3.spec \
benchmarks/rpv3.4.spec benchmarks/rpv3.5.spec 

INPUTFILES = slha2Input lesHouchesInput rpvHouchesInput \
benchmarks/cmssm10.1.1.in benchmarks/cmssm10.1.2.in \
benchmarks/cmssm10.1.3.in benchmarks/cmssm10.1.4.in \
benchmarks/cmssm10.1.5.in \
benchmarks/cmssm10.2.1.in benchmarks/cmssm10.2.2.in \
benchmarks/cmssm10.2.3.in benchmarks/cmssm10.2.4.in \
benchmarks/cmssm10.2.5.in \
benchmarks/cmssm10.3.1.in benchmarks/cmssm10.3.2.in \
benchmarks/cmssm10.3.3.in benchmarks/cmssm10.3.4.in \
benchmarks/cmssm10.3.5.in \
benchmarks/cmssm10.4.1.in benchmarks/cmssm10.4.2.in \
benchmarks/cmssm10.4.3.in benchmarks/cmssm10.4.4.in \
benchmarks/cmssm10.4.5.in \
benchmarks/cmssm40.1.1.in benchmarks/cmssm40.1.2.in \
benchmarks/cmssm40.1.3.in benchmarks/cmssm40.1.4.in \
benchmarks/cmssm40.1.5.in \
benchmarks/cmssm40.2.1.in benchmarks/cmssm40.2.2.in \
benchmarks/cmssm40.2.3.in benchmarks/cmssm40.2.4.in \
benchmarks/cmssm40.2.5.in \
benchmarks/cmssm40.3.1.in benchmarks/cmssm40.3.2.in \
benchmarks/cmssm40.3.3.in benchmarks/cmssm40.3.4.in \
benchmarks/cmssm40.3.5.in \
benchmarks/mamsb1.1.in benchmarks/mamsb1.2.in benchmarks/mamsb1.3.in \
benchmarks/mgmsb1.1.in benchmarks/mgmsb1.2.in benchmarks/mgmsb1.3.in \
benchmarks/mgmsb1.4.in benchmarks/mgmsb1.5.in benchmarks/mgmsb1.6.in \
benchmarks/mgmsb2.1.1.in benchmarks/mgmsb2.1.2.in benchmarks/mgmsb2.1.3.in \
benchmarks/mgmsb2.1.4.in benchmarks/mgmsb2.1.5.in benchmarks/mgmsb2.1.6.in \
benchmarks/mgmsb2.2.1.in benchmarks/mgmsb2.2.2.in benchmarks/mgmsb2.2.3.in \
benchmarks/mgmsb2.2.4.in benchmarks/mgmsb2.2.5.in benchmarks/mgmsb2.2.6.in \
benchmarks/pmssm1.1.in benchmarks/pmssm1.2.in benchmarks/pmssm1.3.in \
benchmarks/pmssm1.4.in benchmarks/pmssm1.5.in benchmarks/pmssm1.6.in \
benchmarks/pmssm2.1.in benchmarks/pmssm2.2.in benchmarks/pmssm2.3.in \
benchmarks/pmssm2.4.in benchmarks/pmssm2.5.in benchmarks/pmssm2.6.in \
benchmarks/rpv1.1.in benchmarks/rpv1.2.in benchmarks/rpv1.3.in \
benchmarks/rpv1.4.in benchmarks/rpv1.5.in \
benchmarks/rpv2.1.in benchmarks/rpv2.2.in benchmarks/rpv2.3.in \
benchmarks/rpv2.4.in benchmarks/rpv2.5.in \
benchmarks/rpv3.1.in benchmarks/rpv3.2.in benchmarks/rpv3.3.in \
benchmarks/rpv3.4.in benchmarks/rpv3.5.in 

dist_pkgdata_DATA = $(INPUTFILES) $(OUTPUTFILES) CHANGES GUIDELINES COPYING $(DOCFILES)
CLEANFILES=$(OUTPUTFILES) *.x

%Output : %Input softpoint.x 
	./softpoint.x leshouches < $< > $@

outputTest: softsusy.x
	./softsusy.x > outputTest

rpvOutputTest: rpvsoftsusy.x
	./rpvsoftsusy.x > rpvOutputTest

neutOutputTest: rpvneut.x
	./rpvneut.x > neutOutputTest

benchmarks/%.spec: benchmarks/%.in softpoint.x
	$(install_sh) -d benchmarks
	./softpoint.x leshouches < $< > $@

updateWeb: dist 
	scp softsusy-*.tar.gz allanach@login.hepforge.org:~/softsusy/downloads; 
	scp *{Input,Output}* allanach@login.hepforge.org:~/softsusy/public_html/
	scp {def.h,outputTest,README} allanach@login.hepforge.org:~/softsusy/pub
	scp Makefile.alt allanach@login.hepforge.org:~/softsusy/public_html/Make
	scp -r doc allanach@login.hepforge.org:~/softsusy/public_html/; \
	doxygen softsusy.cfg &&\
	scp -r html allanach@login.hepforge.org:~/softsusy/public_html; \
	rm -fr html *.tar.gz
