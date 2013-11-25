#! /bin/sh

DIE=0
PROJECT="Softsusy"
CONFIG_DIR="Config"

exists_in_path() {
    # This function will try to locate an executable [$1] in $PATH.
    #
    # The result of the search is stored in cmd, which should be
    # immediately copied, since the variables value will be
    # overwritten at next invocation of this function.

    # Assert that we got enough arguments
    if test $# -ne 1 ; then
        echo "exists_in_path: Exactly one argument required"
        return 1
    fi

    cmd=$(command -v -- "$1")
    case "$cmd" in
	/*) return 0 ;;
	alias\ *) return 1 ;; # alias
	*) return 1 ;; # built-in or function
    esac
}

call_libtoolize() {
    # This function calls libtoolize.
    libtoolize_cmd="libtoolize"

    # on Darwin systems prefer glibtoolize over libtoolize
    case `uname` in
        Darwin*)
            exists_in_path glibtoolize
            if [ -n "$cmd" ] ; then
                libtoolize_cmd="glibtoolize"
            else
                libtoolize_cmd="libtoolize"
            fi
            ;;
    esac

    exists_in_path "${libtoolize_cmd}"
    if [ -n "$cmd" ] ; then
        ${libtoolize_cmd} -c -f
    else
        echo "Error: ${libtoolize_cmd} not found."
        echo "You must have libtool installed to run this script."
        exit 1
    fi
}

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
call_libtoolize
automake --add-missing --copy
autoconf

exit
