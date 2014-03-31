lib_LTLIBRARIES = 
bin_PROGRAMS =
CLEANFILES =
BUILT_SOURCES=
EXTRA_DIST=
pkginclude_HEADERS=
dist_pkgdata_DATA = 
#SUBDIRS = src doc
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} 

AUTOMAKE_OPTIONS = subdir-objects

include src/softsusy.amk
include doc/doc.amk
include benchmarks/benchmarks.amk


TESTOUTPUT = inOutFiles/slha2Output inOutFiles/lesHouchesOutput inOutFiles/rpvHouchesOutput \
 inOutFiles/nmssmSLHAZ3Output inOutFiles/nmssmSLHAnoZ3Output

OUTPUTFILES = $(TESTOUTPUT) \
 inOutFiles/neutOutputTest inOutFiles/rpvOutputTest \
 inOutFiles/outputTest inOutFiles/outputTest-nmssm 

INPUTFILES = inOutFiles/slha2Input inOutFiles/lesHouchesInput inOutFiles/rpvHouchesInput \
inOutFiles/nmssmSLHAZ3Input inOutFiles/nmssmSLHAnoZ3Input


dist_pkgdata_DATA += softsusy_nmssmtools.x.in setup_nmssmtools.sh Makefile.nmssmtools $(INPUTFILES) $(OUTPUTFILES) CHANGES GUIDELINES COPYING $(DOCFILES)

CLEANFILES += $(OUTPUTFILES) *.x

inOutFiles/%Output : inOutFiles/%Input softpoint.x 
	$(install_sh) -d inOutFiles
	./softpoint.x leshouches < $< > $@

inOutFiles/outputTest: softsusy.x
	$(install_sh) -d inOutFiles
	./softsusy.x > inOutFiles/outputTest

inOutFiles/outputTest-nmssm: softsusy-nmssm.x
	$(install_sh) -d inOutFiles
	./softsusy-nmssm.x > inOutFiles/outputTest-nmssm

inOutFiles/rpvOutputTest: rpvsoftsusy.x
	$(install_sh) -d inOutFiles
	./rpvsoftsusy.x > inOutFiles/rpvOutputTest

inOutFiles/neutOutputTest: rpvneut.x
	$(install_sh) -d inOutFiles
	./rpvneut.x > inOutFiles/neutOutputTest

.PHONY: doc test updateWeb

test: softpoint.x
	git checkout $(TESTOUTPUT)
	cp $(TESTOUTPUT) .test; rm $(TESTOUTPUT); make $(TESTOUTPUT); \
	ndiff --relative-error 1.0e-3 inOutFiles/slha2Output .test/slha2Output; \
	ndiff --relative-error 1.0e-3 inOutFiles/lesHouchesOutput .test/lesHouchesOutput; \
	ndiff --relative-error 1.0e-3 inOutFiles/rpvHouchesOutput .test/rpvHouchesOutput; \
	ndiff --relative-error 1.0e-3 inOutFiles/nmssmSLHAZ3Output .test/nmssmSLHAZ3Output; \
	ndiff --relative-error 1.0e-3 inOutFiles/nmssmSLHAnoZ3Output .test/nmssmSLHAnoZ3Output

#doc:
#	@echo Making $@ in $(top_builddir)/doc
#	cd $(top_builddir)/doc && make $@

updateWeb: dist 
	scp softsusy-${VERSION}.tar.gz allanach@login.hepforge.org:~/softsusy/downloads; 
	scp inOutFiles/*{Input,Output}* allanach@login.hepforge.org:~/softsusy/public_html/code
	scp {src/def.h,inOutFiles/outputTest,README} allanach@login.hepforge.org:~/softsusy/public_html/code
	scp Makefile.alt allanach@login.hepforge.org:~/softsusy/public_html/Makefile
	scp -r doc allanach@login.hepforge.org:~/softsusy/public_html/; \
	doxygen softsusy.cfg &&\
	scp -r html allanach@login.hepforge.org:~/softsusy/public_html; \
	scp -r CHANGES allanach@login.hepforge.org:~/softsusy/public_html/CHANGES.${VERSION}.html; \
	rm -fr html *.tar.gz

install-data-local:
	rm -f softsusy.x softsusy-nmssm.x softpoint.x rpvsoftsusy.x rpvneut.x *.la
