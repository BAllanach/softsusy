SUBDIRS = src doc

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

OUTPUTFILES = inOutFiles/slha2Output inOutFiles/lesHouchesOutput inOutFiles/rpvHouchesOutput \
 inOutFiles/nmssmSLHAZ3Output inOutFiles/nmssmSLHAnoZ3Output 

INPUTFILES = inOutFiles/slha2Input inOutFiles/lesHouchesInput inOutFiles/rpvHouchesInput \
inOutFiles/nmssmSLHAZ3Input inOutFiles/nmssmSLHAnoZ3Input 

dist_pkgdata_DATA = softsusy_nmssmtools.x.in setup_nmssmtools.sh Makefile.nmssmtools $(INPUTFILES) $(OUTPUTFILES) CHANGES GUIDELINES COPYING $(DOCFILES)
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

.PHONY: doc test updateWeb

test: softpoint.x
	git checkout $(OUTPUTFILES)
	cp $(OUTPUTFILES) .test; rm $(OUTPUTFILES); make $(OUTPUTFILES); \
	ndiff --relative-error 1.0e-3 inOutFiles/slha2Output .test/slha2Output; \
	ndiff --relative-error 1.0e-3 inOutFiles/lesHouchesOutput .test/lesHouchesOutput; \
	ndiff --relative-error 1.0e-3 inOutFiles/rpvHouchesOutput .test/rpvHouchesOutput; \
	ndiff --relative-error 1.0e-3 inOutFiles/nmssmSLHAZ3Output .test/nmssmSLHAZ3Output; \
	ndiff --relative-error 1.0e-3 inOutFiles/nmssmSLHAnoZ3Output .test/nmssmSLHAnoZ3Output

doc:
	@echo Making $@ in $(top_builddir)/doc
	cd $(top_builddir)/doc && make $@

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
