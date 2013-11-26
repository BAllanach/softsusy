AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = softsusy.x softsusy-nmssm.x softpoint.x rpvsoftsusy.x rpvneut.x 

lib_LTLIBRARIES = libsoft.la 

softsusy_x_SOURCES = src/main.cpp
softsusy_x_LDADD = libsoft.la  

softsusy_nmssm_x_SOURCES = src/main-nmssm.cpp
softsusy_nmssm_x_LDADD = libsoft.la

softpoint_x_SOURCES = src/softpoint.cpp 
softpoint_x_LDADD = libsoft.la

rpvsoftsusy_x_SOURCES = src/rpvmain.cpp
rpvsoftsusy_x_LDADD = libsoft.la 

rpvneut_x_SOURCES = src/rpvneutmain.cpp
rpvneut_x_LDADD = libsoft.la 

pkginclude_HEADERS = \
src/def.h src/main.h src/numerics.h \
src/softpars.h src/susy.h src/linalg.h \
src/physpars.h src/softpoint.h src/twoloophiggs.h \
src/nmssm2loop.h src/dilogwrap.h \
src/lowe.h src/mycomplex.h src/rge.h src/softsusy.h \
src/utils.h src/mssmUtils.h src/rpvsoft.h src/rpvsusypars.h \
src/tensor.h src/rpvmain.h src/flavoursoft.h \
src/xpr-base.h src/xpr-matrix.h src/xpr-vector.h src/rpvneut.h \
src/nmssmsusy.h src/nmssmsoftpars.h  \
src/nmssmsoftsusy.h src/nmssmUtils.h

libsoft_la_SOURCES = src/def.cpp src/linalg.cpp src/lowe.cpp src/numerics.cpp  src/rpvneut.cpp \
src/physpars.cpp src/rge.cpp src/nmssm2loop.f src/dilogwrap.f \
 src/nmssmsoftsusy.cpp src/nmssmsusy.cpp src/nmssmsoftpars.cpp src/nmssmUtils.cpp \
src/susy.cpp src/utils.cpp src/mssmUtils.cpp \
src/tensor.cpp src/rpvsusypars.cpp src/rpvsoft.cpp src/flavoursoft.cpp src/twoloophiggs.f 
libsoft_la_LDFLAGS = $(FLIBS)

DOCFILES=doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf \
doc/nmssmManual.pdf
EXTRA_DIST=$(DOCFILES)

OUTPUTFILES = inOutFiles/slha2Output inOutFiles/lesHouchesOutput inOutFiles/rpvHouchesOutput \
 inOutFiles/nmssmSLHAZ3Output inOutFiles/nmssmSLHAnoZ3Output inOutFiles/neutOutputTest \
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

INPUTFILES = inOutFiles/slha2Input inOutFiles/lesHouchesInput inOutFiles/rpvHouchesInput \
inOutFiles/nmssmSLHAZ3Input inOutFiles/nmssmSLHAnoZ3Input \
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

inOutFiles/%Output : inOutFiles/%Input softpoint.x 
	./softpoint.x leshouches < $< > $@

inOutFiles/outputTest: softsusy.x
	./softsusy.x > inOutFiles/outputTest

inOutFiles/rpvOutputTest: rpvsoftsusy.x
	./rpvsoftsusy.x > inOutFiles/rpvOutputTest

inOutFiles/neutOutputTest: rpvneut.x
	./rpvneut.x > inOutFiles/neutOutputTest

benchmarks/%.spec: benchmarks/%.in softpoint.x
	$(install_sh) -d benchmarks
	./softpoint.x leshouches < $< > $@

updateWeb: dist 
	scp softsusy-*.tar.gz allanach@login.hepforge.org:~/softsusy/downloads; 
	scp inOutFiles/*{Input,Output}* allanach@login.hepforge.org:~/softsusy/public_html/code
	scp {def.h,inOutFiles/outputTest,README} allanach@login.hepforge.org:~/softsusy/public_html/code
	scp Makefile.alt allanach@login.hepforge.org:~/softsusy/public_html/Makefile
	scp -r doc allanach@login.hepforge.org:~/softsusy/public_html/; \
	doxygen softsusy.cfg &&\
	scp -r html allanach@login.hepforge.org:~/softsusy/public_html; \
	scp -r CHANGES $cedar:~/softsusy/public_html/CHANGES.html.3.4.0; \
	rm -fr html *.tar.gz
