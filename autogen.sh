#! /bin/sh

DIE=0
PROJECT="Softsusy"
CONFIG_DIR="Config"

(aclocal --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have aclocal installed to compile $PROJECT."
	echo "Download the appropriate package for your distribution,"
	echo "or get the source tarball at ftp://ftp.gnu.org/pub/gnu/"
	DIE=1
}

(autoconf --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have autoconf installed to compile $PROJECT."
	echo "Download the appropriate package for your distribution,"
	echo "or get the source tarball at ftp://ftp.gnu.org/pub/gnu/"
	DIE=1
}

(autoheader --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have autoheader installed to compile $PROJECT."
	echo "Download the appropriate package for your distribution,"
	echo "or get the source tarball at ftp://ftp.gnu.org/pub/gnu/"
	DIE=1
}

(automake --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have automake installed to compile $PROJECT."
	echo "Get ftp://ftp.gnu.org/pub/gnu/automake-1.3.tar.gz"
	echo "(or a newer version if it is available)"
	DIE=1
}

if test "$DIE" -eq 1; then
	exit 1
fi


rm -f config.cache
aclocal
if [ ! -d "$CONFIG_DIR" ]
then
	echo "Creating Config/"
	mkdir "$CONFIG_DIR"
fi
autoheader
# create file Config/ltmain.sh
case `uname` in Darwin*) glibtoolize -c -f ;;
  *) libtoolize -c -f ;; esac
automake --add-missing --copy
autoconf

exit
