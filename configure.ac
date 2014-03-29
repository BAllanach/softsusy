dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([softsusy],[3.5.0],[softsusy@projects.hepforge.org],[softsusy])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_HEADERS([Config/config.h])
AC_CONFIG_MACRO_DIR([m4])

dnl Put default fortran compiler for this variable
dnl if test "x$CXXFLAGS" = "x"; then
dnl dnl     F77=gfortran
dnl fi
dnl === disable debug symbols by default =====
if test "x$CXXFLAGS" = "x"; then
    CXXFLAGS="-O2 -Wall"
fi
if test "x$FFLAGS" = "x"; then
    FFLAGS="-O2 -Wall"
fi
dnl ==========================================

AC_LANG(C++)

AM_INIT_AUTOMAKE([1.9 foreign])
dnl Checks for programs.
AC_PROG_CXX

if test x"$CXX" = x; then
 AC_MSG_ERROR([could not find C++ compiler])
fi

AC_PROG_INSTALL
AC_PROG_MAKE_SET
#AC_PROG_LN_S

## needed for twoloophiggs
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
##

## If you want shared libraries, configure with --enable-shared
#AC_DISABLE_STATIC
AC_PROG_LIBTOOL

#AM_CONDITIONAL(NEED_APPLE_FIXES,
#		[test "xx${host/darwin/foundit}xx" != "xx${host}xx"])

AC_ARG_ENABLE([full_three_loop], 
               AC_HELP_STRING([--enable-full_three_loop],[Full three-loop RGEs with Yukawa matrices (with mixing)]), 
 	      AC_DEFINE([FULL_THREE_LOOP],[1],[Full 3-loop RGE with mixing]),[])

AC_ARG_ENABLE([full_susy_threshold], 
               AC_HELP_STRING([--enable-full_susy_threshold],[Full two-loop SUSY threshold corr. to fermions third generation]),
             [AC_DEFINE([FULL_SUSY_THRESHOLD],[1],[Full two-loop SUSY threshold corr. to fermions third generation]),
              case "${enableval}" in
              yes) full_susy_threshold=true ;;
              no)  full_susy_threshold=false ;;
              *) AC_MSG_ERROR([bad value ${enableval} for --enable-full_susy_threshold]) ;;
              esac],[full_susy_threshold=false])

AM_CONDITIONAL(FULL_SUSY_THRESHOLD, [test x$full_susy_threshold = xtrue])


AC_ARG_ENABLE([debug], 
               AC_HELP_STRING([--enable-debug],[enables debug printing]), 
 	      AC_DEFINE([DEBUG],[1],[Debugging output]),[])

#Autogen
#AG_PATH_AUTOOPTS([5.17], [], [AC_MSG_ERROR([No suitable version of AutoGen could be found])])

# GiNaC/cln for two-loop decoupling corrections
if test x"$enable_full_susy_threshold" = x"yes"; then
 PKG_CHECK_MODULES([CLN], cln >= 1.3.1)  
 PKG_CHECK_MODULES([GINAC], ginac >= 1.3.5)  

 LIBS="$LIBS $CLN_LIBS $GINAC_LIBS"
 CPPFLAGS="$CPPFLAGS $CLN_CFLAGS $GINAC_CFLAGS"
fi

AC_SUBST(FLIBS)

AC_CONFIG_FILES([Makefile])

AC_OUTPUT

# To run with looptools, run: ./configure CPPFLAGS=-I/home/bca20/code/looptools/LoopTools/install/include LDFLAGS=-L/home/bca20/code/looptools/LoopTools/install/lib/ F77=gfortran and add -looptools to Makefile.am libs
