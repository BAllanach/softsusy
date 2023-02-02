# NB REMINDER: Making an official tag:
# Make sure CHANGES.md is up to date.
# Make sure version number in configure.ac is up-to-date
# Then, comment include doc/doc.amk line
# below. Then,
# git push
# git tag -a v4.X.Y -m "Official version v4.X.Y"; git push origin --tags
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
# Comment this line before making an official tag
include doc/doc.amk
include inOutFiles/inOut.amk

dist_pkgdata_DATA += CHANGES.md $(DOCFILES)

CLEANFILES += $(OUTPUTFILES) *.x

.PHONY: doc test updateWeb install

nodoc: $(OUTPUTFILES) programs

programs: softpoint.x rpvsoftsusy.x softsusy.x softsusy-nmssm.x rpvneut.x 

updateWeb: dist
	scp -P 222 softsusy-${VERSION}.tar.gz allanach@login.hepforge.org:~/softsusy/downloads; 
	scp -P 222 inOutFiles/*{Input,Output}* allanach@login.hepforge.org:~/softsusy/public_html/code
	scp -P 222 {src/def.h,inOutFiles/outputTest,README.md} allanach@login.hepforge.org:~/softsusy/public_html/code
	scp -P 222 -r doc/*.pdf allanach@login.hepforge.org:~/softsusy/public_html/doc; \
	doxygen softsusy.cfg &&\
	scp -P 222 -r html allanach@login.hepforge.org:~/softsusy/public_html; \
	scp -P 222 -r CHANGES allanach@login.hepforge.org:~/softsusy/public_html/CHANGES.${VERSION}.shtml; \
	rm -fr html *.tar.gz

install-data-local:
	rm -f softsusy.x softsusy-nmssm.x softpoint.x rpvsoftsusy.x rpvneut.x *.la

# For distributions: don't re-make the docs
dist-hook:
	sed -i -e 's/^include doc/#include doc/g' $(distdir)/Makefile.am
	rm $(distdir)/doc/doc.amk
	cd $(distdir); autoreconf -vif; cd ..
