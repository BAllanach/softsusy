dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([softsusy],[3.3.7],[softsusy@projects.hepforge.org],[softsusy])
AC_CONFIG_SRCDIR([softpoint.h])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_HEADERS([Config/config.h])

dnl === disable debug symbols by default =====
if test "x$CXXFLAGS" = "x"; then
   CXXFLAGS=-O2
fi
if test "x$FFLAGS" = "x"; then
   FFLAGS=-O2 
fi
dnl ==========================================

AC_LANG(C++)

AM_INIT_AUTOMAKE([1.9 foreign])
dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
#AC_PROG_LN_S

## needed for twoloophiggs
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
##

## If you want shared libraries, configure with --enable-shared
AC_DISABLE_STATIC
AC_PROG_LIBTOOL

#AM_CONDITIONAL(NEED_APPLE_FIXES,
#		[test "xx${host/darwin/foundit}xx" != "xx${host}xx"])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
