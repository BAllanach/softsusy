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
#Comment the next line for a distribution - don't want doc sources
#include doc/doc.amk
include inOutFiles/inOut.amk

DOCFILES = doc/rpcManual.pdf doc/rpvManual.pdf doc/neutManual.pdf \
	   doc/nmssmManual.pdf doc/threeLoop.pdf doc/ho.pdf doc/decays.pdf

dist_pkgdata_DATA += softsusy_nmssmtools.x.in setup_nmssmtools.sh \
NMSSMTools_4.2.1.tgz Makefile.nmssmtools CHANGES GUIDELINES COPYING $(DOCFILES)

CLEANFILES += $(OUTPUTFILES) *.x

.PHONY: doc test updateWeb install

programs: softpoint.x rpvsoftsusy.x softsusy.x softsusy-nmssm.x rpvneut.x decay.x

updateWeb: dist 
	scp softsusy-${VERSION}.tar.gz allanach@login.hepforge.org:~/softsusy/downloads; 
	scp inOutFiles/*{Input,Output}* allanach@login.hepforge.org:~/softsusy/public_html/code
	scp {src/def.h,inOutFiles/outputTest,README} allanach@login.hepforge.org:~/softsusy/public_html/code
	scp -r doc/*.pdf allanach@login.hepforge.org:~/softsusy/public_html/doc; \
	doxygen softsusy.cfg &&\
	scp -r html allanach@login.hepforge.org:~/softsusy/public_html; \
	scp -r CHANGES allanach@login.hepforge.org:~/softsusy/public_html/CHANGES.${VERSION}.html; \
	rm -fr html *.tar.gz

install-data-local:
	rm -f higher.x softsusy.x softsusy-nmssm.x softpoint.x rpvsoftsusy.x rpvneut.x *.la

