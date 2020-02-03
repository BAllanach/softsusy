dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([softsusy], [4.1.9], [B.C.Allanach@damtp.cam.ac.uk], [softsusy], [http://softsusy.hepforge.org/])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_HEADERS([Config/config.h])
AC_CONFIG_MACRO_DIR([m4])

dnl Put default fortran compiler for this variable
dnl F77=gfortran
dnl === disable debug symbols by default =====
if test "x$CXXFLAGS" = "x"; then
    CXXFLAGS="-O2 -std=c++14"
fi
if test "x$CFLAGS" = "x"; then
    CFLAGS="-O2"
fi
if test "x$FFLAGS" = "x"; then
    FFLAGS="-O2"
fi
dnl ==========================================

AC_LANG(C++)

AM_INIT_AUTOMAKE([foreign])
dnl Checks for programs.
AC_PROG_CXX

AM_PROG_AR

if test x"$CXX" = x; then
 AC_MSG_ERROR([could not find C++ compiler])
fi

AC_PROG_CC

if test x"$CC" = x; then
 AC_MSG_ERROR([could not find CC compiler])
fi

#check whether the C compiler is an Apple one
$CC -v &>ver
a=`grep -m 1 -o 'apple' ver`
if test "$a" = "apple"; then
  CFLAGS="-O0 -Wall"
fi

AC_PROG_INSTALL
AC_PROG_MAKE_SET

## needed for twoloophiggs
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
##

## If you want shared libraries, configure with --enable-shared
AC_PROG_LIBTOOL

## Switch on two-loop threshold corrections to third family Yukawa couplings 
## and strong gauge coupling 
AC_ARG_ENABLE([two-loop-gauge-yukawa-compilation], 
               AC_HELP_STRING([--enable-two-loop-gauge-yukawa-compilation],[Compile two-loop SUSY thresholds to the third generation fermionic Yukawa couplings and the strong gauge coupling]),
             [AC_DEFINE([COMPILE_TWO_LOOP_GAUGE_YUKAWA],[1],[Compile leading two-loop SUSY threshold corrections to fermionic third generation Yukawa couplings])
              case "${enableval}" in
              yes) two_loop_gauge_yukawa_compilation=true ;;
              no)  two_loop_gauge_yukawa_compilation=false ;;
              *) AC_MSG_ERROR([bad value ${enableval} for --enable-two_loop_gauge_yukawa_compilation]) ;;
              esac],[two_loop_gauge_yukawa_compilation=false])
AM_CONDITIONAL(COMPILE_TWO_LOOP_GAUGE_YUKAWA, [test x$two_loop_gauge_yukawa_compilation = xtrue])

# Use Himalaya for 3L contributions to Mh
AC_ARG_ENABLE([himalaya], 
              AC_HELP_STRING([--enable-himalaya],[Include 3-loop Higgs mass contributions from Himalaya]),
             [AC_DEFINE([ENABLE_HIMALAYA],[1],[Include 3-loop Higgs mass contributions from Himalaya])
              case "${enableval}" in
              yes) enable_himalaya=yes ;;
              no)  enable_himalaya=no ;;
              *) AC_MSG_ERROR([bad value ${enableval} for --enable-himalaya]) ;;
              esac],[enable_himalaya=no])
AM_CONDITIONAL(ENABLE_HIMALAYA, [test x$enable_himalaya = xyes])

AS_IF([test "x$enable_himalaya" = "xyes"],
      [AC_LANG_PUSH([C++])
       AC_CHECK_HEADERS([Eigen/Dense], [],
                        [AC_MSG_ERROR([Eigen headers not found])]
                       )
       AC_LANG_POP([C++])]
     )

# Use looptools instead of SOFTSUSY's built in 1-loop functions
AC_ARG_ENABLE([looptools], 
               AC_HELP_STRING([--enable-looptools],[Looptools evaluation of Pasarino-Veltman functions]),
             [AC_DEFINE([LOOPTOOLS],[1],[Looptools evaluation of Pasarino-Veltman functions])
              case "${enableval}" in
              yes) looptools_compilation=true ;;
              no)  looptools_compilation=false ;;
              *) AC_MSG_ERROR([bad value ${enableval} for --enable-looptools]) ;;
              esac],[looptools_compilation=false])
AM_CONDITIONAL(LOOPTOOLS, [test x$looptools_compilation = xtrue])

# AC_ARG_ENABLE([debug], 
#               AC_HELP_STRING([--enable-debug],[enables debug printing in two-loop threshold corrections]), 
# 	      AC_DEFINE([DEBUG],[1],[Threshold correction debugging output]),[])

AC_ARG_ENABLE([gsl],
              AC_HELP_STRING([--enable-gsl],[enables linking against GSL (default: no)]),
              [],[enable_gsl=no])

AS_IF([test "x$enable_gsl" = "xyes"],
   [AC_CHECK_LIB([m],[cos],[],
                [AC_MSG_ERROR([libm not found, but required for enabling GSL])])
   AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],
                [AC_MSG_ERROR([libgslcblas not found, but required for enabling GSL])])
   AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],
                [AC_MSG_ERROR([libgsl not found, but required for enabling GSL])])
   AC_CHECK_HEADERS([gsl/gsl_multiroots.h gsl/gsl_deriv.h],[],
                   [AC_MSG_ERROR([GSL must be installed to use GSL routines])], [])
   AC_DEFINE([ENABLE_GSL],[1],[Enable GSL routines])]
)

AM_CONDITIONAL([ENABLE_GSL], [test "x$enable_gsl" = "xyes"])

#Autogen
#AG_PATH_AUTOOPTS([5.17], [], [AC_MSG_ERROR([No suitable version of AutoGen could be found])])

# TSIL for 2-loop sparticle mass functions. But at the moment we are including 
# it, so no worries.
#if test x"$two_loop_sparticle_mass_compilation" = x"yes"; then
# PKG_CHECK_MODULES([TSIL], tsil >= 1.2.1)
#
# LIBS="$LIBS $TSIL_LIBS"
# CPPFLAGS="$CPPFLAGS TSIL_CFLAGS"
#fi
# For TSIL
# TSIL_LIB="-L/usr/lib -ltsil" TSIL_CFLAGS="-I/usr/include" ./configure --enable-two-loop-sparticle-mass-compilation

AC_SUBST(FLIBS)

AC_CONFIG_FILES([Makefile])

AC_OUTPUT

# To run with looptools, run: ./configure FFLAGS="-DQUAD=1" --enable-looptools CPPFLAGS=-I/home/bca20/code/LoopTools-2.13/install/include LDFLAGS=-L/home/bca20/code/LoopTools-2.13/build

# For GINAC/CLN
# GINAC_LIBS="-L/usr/lib -lginac" GINAC_CFLAGS="-I/usr/include" CLN_LIBS="-L/usr/lib -lcln" CLN_CFLAGS="-I/usr/include" ./configure --enable-two-loop-gauge-yukawa-compilation --enable-three-loop-rge-compilation

