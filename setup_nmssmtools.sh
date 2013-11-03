#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

# absolute path to this script
ABSBASEDIR=$(cd $BASEDIR; pwd)

enable_compile_nmssmtools=no
nmssmtools_dir=""

usage() {
cat <<EOF
Usage: ./`basename $0` path-to-nmssmtools-dir [options]
Options:
  --compile|-c       compile NMSSMTools
  --nmssmtools-dir=  Path to NMSSMTools
  --help|-h          Print this message and exit
EOF
}

check_nmssmtools() {
    if ! test -d "${nmssmtools_dir}"; then
        echo "Error: Cannot find NMSSMTools directory \"${nmssmtools_dir}\""
        exit 1
    fi
    if ! test -r "${nmssmtools_dir}/main/Makefile"; then
        echo "Error: ${nmssmtools_dir}/main/Makefile does not exist."
        exit 1
    fi
    if ! test -r "${nmssmtools_dir}/main/nmspec.f"; then
        echo "Error: ${nmssmtools_dir}/main/nmspec.f does not exist."
        exit 1
    fi
}

copy_file() {
    if test $# -ne 2 ; then
        echo "Internal error: copy_file not called with two arguments"
    fi
    printf "Backup file $2 to $2~ ..."
    if cp "$2" "$2~"; then
        echo " done"
    else
        echo " failed"
        exit 1
    fi
    printf "Copying $1 to $2 ..."
    if cp "$1" "$2"; then
        echo " done"
    else
        echo " failed"
        exit 1
    fi
}

copy_files() {
    copy_file ${BASEDIR}/nmProcessSpec.f ${nmssmtools_dir}/main/nmProcessSpec.f
    copy_file ${BASEDIR}/Makefile.nmssmtools ${nmssmtools_dir}/main/Makefile
}

compile_nmssmtools() {
    printf "Compiling NMSSMTools ..."
    if cd ${nmssmtools_dir} && make init > make.nmssmtools 2>&1 && make >> make.nmssmtools 2>&1; then
        echo " done"
    else
        echo " error"
        exit 1
    fi
}

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case "$1" in
            --compile|-c)            enable_compile_nmssmtools=yes ;;
            --nmssmtools-dir=*)      nmssmtools_dir="$optarg" ;;
            --help|-h)               usage; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

check_nmssmtools
copy_files
if test "x${enable_compile_nmssmtools}" = "xyes"; then
    compile_nmssmtools
else
    echo ""
    echo "Next steps: recompile NMSSMTools"
    echo "  $ cd ${nmssmtools_dir}"
    echo "  $ make init"
    echo "  $ make"
fi
