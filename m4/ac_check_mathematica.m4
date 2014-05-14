dnl @synopsis AC_CHECK_MATHEMATICA
dnl 
dnl This macro tests if Mathematica executable
dnl is installed. If math is installed it sets 
dnl $MATHEMATICA to the right value
dnl
dnl
AC_DEFUN([AC_CHECK_MATHEMATICA], [
 AC_PATH_PROG(MATHEMATICA, math, no)
 if test "x$MATHEMATICA" = "xno"; then
   ifelse([$3], , :,[$3])
 else
   mathematica_min_version=ifelse([$1], ,6.0,[$1])
   AC_MSG_CHECKING([for Mathematica version >= $mathematica_min_version])
   mathematica_installed_version=`echo -e \\$VersionNumber | $MATHEMATICA -noprompt | tr -d "\"\n"`
   math_found=`expr $mathematica_installed_version \> $mathematica_min_version`
   if test "x$math_found" == "x1"; then
	mathlink_base_directory=`echo -e \\$InstallationDirectory | $MATHEMATICA -noprompt | tr -d "\"\n"`
	mathlink_dev_kit=$mathlink_base_directory/SystemFiles/Links/MathLink/DeveloperKit
	mathlink_system_id=`echo -e \\$SystemID | $MATHEMATICA -noprompt | tr -d "\"\n"`
	mathlink_compiler_additions=$mathlink_dev_kit/$mathlink_system_id/CompilerAdditions
   	AC_MSG_RESULT([yes])
	AC_PATH_PROGS(MPREP,mprep,no,path= $mathlink_compiler_additions)
	AC_PATH_PROGS(MCC,mcc,no,path= $mathlink_compiler_additions)

	MATHLINK_CFLAGS=-I$mathlink_compiler_additions
	case $mathlink_system_id in
	Linux)
		mathlink_lib=-lML32i3
		mathlink_extra_libs="-lm -lpthread -lrt -lstdc++"	
		;;	
	Linux-x86-64) 
		mathlink_lib=-lML64i3
		mathlink_extra_libs="-lm -lpthread -lrt -lstdc++"	
		;;	
	MacOSX)
		mathlink_lib=-lMLi3
		mathlink_extra_libs="-lstdc++ -framework Fondation"	
		;;
	*)
		echo "Mathlink $mathlink_system_id not implemented yet"
	esac
	MATHLINK_LDFLAGS="-L$mathlink_compiler_additions $mathlink_lib $mathlink_extra_libs"
	ifelse([$2], , :,[$2])
   else
        AC_MSG_RESULT([no]) 
	ifelse([$3], , :,[$3])
   fi
 fi
AC_SUBST(MATHEMATICA)
AC_SUBST(MCC)
AC_SUBST(MPREP)
AC_SUBST(MATHLINK_CFLAGS)
AC_SUBST(MATHLINK_LDFLAGS)
]) 
