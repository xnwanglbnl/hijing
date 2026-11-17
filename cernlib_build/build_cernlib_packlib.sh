#!/bin/sh
# Build CERNLIB 2006 packlib + kernlib from github.com/apc-llc/cernlib
# Usage:  sh build_cernlib_packlib.sh
# Result: ./cernlib/2006/lib/libpacklib.a and libkernlib.a
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
for t in gfortran gcc make perl git imake makedepend; do
  command -v $t >/dev/null 2>&1 || { echo "missing tool: $t (imake/makedepend come from xutils-dev or imake package)"; exit 1; }
done
if [ ! -d cernlib ]; then
  git clone --depth 1 https://github.com/apc-llc/cernlib.git
fi
cd cernlib
chmod +x build64 gcc g77 gmake

# gfortran >= 10 flags
grep -q "fallow-invalid-boz" g77 || \
  sed -i 's|/usr/bin/gfortran -std=legacy|/usr/bin/gfortran -std=legacy -fallow-invalid-boz -fallow-argument-mismatch|' g77

# glibc >= 2.32: sys_errlist -> strerror
f=2006/src/packlib/cspack/sysreq/serror.c
grep -q "sys_nerr" $f && {
  sed -i 's|^extern int      sys_nerr;.*|#include <string.h>|; s|^extern char \*sys_errlist\[\];.*||' $f
  sed -i 's|if ((n>0) \&\& (n<sys_nerr)) {|if (n>0) {|; s|return(sys_errlist\[n\]);|return(strerror(n));|' $f
}
f=2006/src/packlib/cspack/sysreq/socket.c
grep -q "sys_errlist\[errno\]" $f && {
  sed -i 's|^extern int      sys_nerr;.*|#include <string.h>|; s|^extern char \*sys_errlist\[\];.*||' $f
  sed -i 's|return(sys_errlist\[errno\]);|return(strerror(errno));|' $f
}
f=2006/src/packlib/kernlib/kernbit/z268/systems.c
grep -q "sys_errlist\[errno\]" $f && {
  sed -i 's|^extern char     \*sys_errlist\[\] ;.*|#include <string.h>|; s|sys_errlist\[errno\]|strerror(errno)|g' $f
}

# no Motif: make the default MotifDependantMakeVar depend on HAVE_MOTIF (unset below)
f=2006/src/config/Imake.rules
grep -q "NOMOTIF_PATCH" $f || awk '
/^variable=value/ && prev ~ /MotifDependantMakeVar/ {
  print "ifneq ($(strip $(HAVE_MOTIF)),)\t\t\t\t\t\t@@\\"
  print "    variable=value\t\t\t\t\t\t\t@@\\"
  print "endif\t\t\t\t\t\t\t\t\t@@\\"
  print "/* NOMOTIF_PATCH */"
  prev=$0; next }
{ print; prev=$0 }' $f > $f.new && mv $f.new $f
unset HAVE_MOTIF

# binutils >= 2.35: drop the 'l' modifier of ar
cat > ar <<'EOS'
#!/usr/bin/perl
my(@a) = @ARGV;
if (@a && $a[0] !~ /^-/) { $a[0] =~ s/l//g; }
exit(system("/usr/bin/ar", @a) >> 8);
EOS
chmod +x ar

# build packlib only
export CERN_LEVEL=2006
export CERN=$(pwd)
export CERN_ROOT=$CERN/$CERN_LEVEL
export CVSCOSRC=$CERN_ROOT/src
export PATH=$CERN_ROOT/bin:$CERN_ROOT/src/scripts:$CERN:$PATH
cd $CERN_ROOT
mkdir -p build bin lib build/log
cd $CERN_ROOT/build
[ -f packlib/Makefile ] || make tree -j8 TOP=$CERN_ROOT/build > log/tree 2>&1 || {
  echo "imake tree failed, see $CERN_ROOT/build/log/tree"; exit 1; }
make bin/kuipc TOP=$CERN_ROOT/build > log/kuipc 2>&1 || {
  echo "kuipc build failed, see $CERN_ROOT/build/log/kuipc"; exit 1; }
cd $CERN_ROOT/build/packlib
C_INCLUDE_PATH=/usr/include/freetype2 make -j8 install.lib TOP=$CERN_ROOT/build > ../log/packlib 2>&1 || {
  echo "packlib build failed, see $CERN_ROOT/build/log/packlib"; exit 1; }
ls -la $CERN_ROOT/lib/libpacklib.a $CERN_ROOT/lib/libkernlib.a
echo "CERNLIB packlib/kernlib built: $CERN_ROOT/lib"
