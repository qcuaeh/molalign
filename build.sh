#!/bin/bash -e
shopt -s nullglob
unalias -a

toarray() {
   IFS=\  read -r -a "$1" <<< "${!1}"
}

clean_build() {
   if test -d build; then
      pushd build >/dev/null
      for file in *.f90 *.mod *.o; do
         rm "$file"
      done
      popd >/dev/null
   fi
   unset obj_files
   unset f2py_files
}

compile() {
   pic=false
   debug=false
   toarray std_flags
   toarray pic_flags
   toarray optim_flags
   toarray debug_flags
   options=$(getopt -o '' -al pic,debug -- "$@") || exit
   eval set -- "$options"
   while true; do
      case "$1" in
      --pic) pic=true; shift ;;
      --debug) debug=true; shift ;;
      --) shift; break ;;
      *) exit
      esac
   done
   srcdir=$topdir/$1
   if test ! -d "$srcdir"; then
      echo Error: $srcdir does not exist
      exit 1
   fi
   pushd "$buildir" >/dev/null
   flags=("${std_flags[@]}")
   if $pic; then
      flags+=("${pic_flags[@]}")
   fi
   if $debug; then
      flags+=("${debug_flags[@]}")
   else
      flags+=("${optim_flags[@]}")
   fi
   while IFS= read -r srcfile; do
      prefix=${srcfile%.*}
      objfile=$prefix.o
      if ! test -e "$objfile" \
      || ! test -e "$srcfile" \
      || ! diff -q "$srcfile" "$srcdir/$srcfile" >/dev/null
      then
         echo Compiling $prefix...
         if test "$srcdir" != "$PWD"; then
            cp -f "$srcdir/$srcfile" "$srcfile"
         fi
         "$F90" "${flags[@]}" -c "$srcfile" -o "$objfile"
      fi
      obj_files+=("$objfile")
   done < <(grep -v ^# "$srcdir/f90_files")
   popd >/dev/null
   if test -f "$srcdir/f2py_files"; then
      while IFS= read -r f2pyfile; do
         f2py_files+=("$f2pyfile")
      done < <(grep -v ^# "$srcdir/f2py_files")
   fi
}

makeprog() {
   name=$1
   executable=$buildir/$1
   if test -z "$name"; then
      echo Error: name is empty
      exit 1
   fi
   echo Linking program...
   pushd "$buildir" > /dev/null
   "$F90" -o "$name" "${obj_files[@]}" -llapack
   popd > /dev/null
   echo Done
}

makelib() {
   name=$1
   if test -z "$name"; then
      echo Error: name is empty
      exit 1
   fi
   echo Linking dynamic library...
   pushd "$buildir" >/dev/null
   "$F90" -shared -o "$name.so" "${obj_files[@]}" -llapack
   popd >/dev/null
   echo Done
}

makepyext() {
   name=$1
   if test -z "$name"; then
      echo Error: name is empty
      exit 1
   fi
   if ! type "$F2PY" &>/dev/null; then
      echo Error: F2PY executable not found 
      exit 1
   fi
   pushd "$buildir" >/dev/null
   echo Linking extension module...
   echo "$f2py_f2cmap" > "$buildir/.f2py_f2cmap"
   "$F2PY" -h "$name.pyf" -m "$name" --overwrite-signature "${f2py_files[@]}" --quiet
   "$F2PY" -c "$name.pyf" --fcompiler=gnu95 --link-lapack "${obj_files[@]}" --quiet
   popd >/dev/null
   echo Done
}

runtests() {
   testdir=$topdir/$1
   shift
   if test -z "$executable"; then
      echo Error: executable is not set
      exit 1
   fi
   for file in "$testdir"/*.out; do
      name=$(basename "$file")
      echo -n Running test ${name%.out}.xyz...
      if diff -bB "$file" <("$executable" "$testdir/${name%.out}.xyz" "$@") > /dev/null; then
          echo \ pass
      else
          echo \ failed
      fi
   done
   for file in aligned_*.xyz; do
      rm "$file"
   done
   echo Done
}

if test ! -e ./build.env; then
   echo Error: build.env does not exist
   exit 1
elif test ! -f ./build.env; then
   echo Error: build.env does exist but is not a file
   exit 1
fi

topdir=$(dirname "$(readlink -e "$0")")
buildir=$topdir/build

if test ! -e "$buildir"; then
   mkdir "$buildir"
elif test ! -d "$buildir"; then
   echo Error: $buildir does exist but is not a directory
   exit 1
fi

# Set environment
while IFS= read -r line; do
   var=${line%%=*}
   value=${line#*=}
   declare -- "$var"="$value"
done < <(grep -v ^# ./build.env)

# Build program
compile molalignlib
compile molalign
makeprog molalign
clean_build

# Run tests
runtests tests/r05 -test -rec 10 -sort -fast -tol 0.17
#runtests tests/r10 -test -rec 10 -sort -fast -tol 0.35
#runtests tests/r20 -test -rec 10 -sort -fast -tol 0.69

# Build dynamic library
#compile -pic molalignlib
#makelib molalignlib
#clean_build

# Build python extension module
#compile -pic molalignlib
#makepyext molalignlibext
#clean_build
