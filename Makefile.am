# NB REMINDER: Making an official tag:
# Make sure CHANGES.md is up to date.
# Make sure version number in configure.ac is up-to-date
# make -j8 dist
# mv softsusy-4.X.Y.tar.gz tags/
# git status
# git add tags/softsusy-4.X.Y.tar.gz
# git commit
# git push
# git tag -a v4.X.Y -m "Official version v4.X.Y"
# git push origin --tags
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

# NB new version should be ready: remember to update README.md
updateWeb: dist
	mv softsusy-${VERSION}.tar.gz tags
	git add tags/softsusy-${VERSION}.tar.gz
	git commit -a -m "Added dist to tags/ directory"
	git tag -a v${VERSION} -m "Official version v${VERSION}"
	git commit -a -m "Official version v${VERSION}"
	git push origin --tags
	git push

install-data-local:
	rm -f softsusy.x softsusy-nmssm.x softpoint.x rpvsoftsusy.x rpvneut.x *.la

# For distributions: don't re-make the docs
dist-hook:
	sed -i -e 's/^include doc/#include doc/g' $(distdir)/Makefile.am
	rm $(distdir)/doc/doc.amk
	cd $(distdir); autoreconf -vif; cd ..
