ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

lib_LTLIBRARIES = libsoft.la 

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

libsoft_la_LDFLAGS = @FLIBS@
