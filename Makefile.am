# Adding comment
lib_LTLIBRARIES = 
bin_PROGRAMS = 
CLEANFILES =
BUILT_SOURCES=
EXTRA_DIST=
pkginclude_HEADERS=
dist_pkgdata_DATA = 
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

AUTOMAKE_OPTIONS = subdir-objects


include src/softsusy.amk
DOCFILES = doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf \
	   doc/nmssmManual.pdf doc/threeLoop.pdf doc/ho.pdf doc/decays.pdf
include doc/doc.amk
include inOutFiles/inOut.amk

dist_pkgdata_DATA += CHANGES GUIDELINES COPYING $(DOCFILES)

CLEANFILES += $(OUTPUTFILES) *.x

.PHONY: doc test updateWeb install

programs: softpoint.x rpvsoftsusy.x softsusy.x softsusy-nmssm.x rpvneut.x decay.x

updateWeb: dist
	scp -P 222 softsusy-${VERSION}.tar.gz allanach@login.hepforge.org:~/softsusy/downloads; 
	scp -P 222 inOutFiles/*{Input,Output}* allanach@login.hepforge.org:~/softsusy/public_html/code
	scp -P 222 {src/def.h,inOutFiles/outputTest,README} allanach@login.hepforge.org:~/softsusy/public_html/code
	scp -P 222 -r doc/*.pdf allanach@login.hepforge.org:~/softsusy/public_html/doc; \
	doxygen softsusy.cfg &&\
	scp -P 222 -r html allanach@login.hepforge.org:~/softsusy/public_html; \
	scp -P 222 -r CHANGES allanach@login.hepforge.org:~/softsusy/public_html/CHANGES.${VERSION}.html; \
	rm -fr html *.tar.gz

install-data-local:
	rm -f higher.x softsusy.x softsusy-nmssm.x softpoint.x rpvsoftsusy.x rpvneut.x *.la

# For distributions: don't re-make the docs
dist-hook:
	sed -i -e 's/^include doc/#include doc/g' $(distdir)/Makefile.am
	rm $(distdir)/doc/doc.amk
	cd $(distdir); autoreconf -vif; cd ..
