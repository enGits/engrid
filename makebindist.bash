#!/bin/bash
set -e -x -u

# Check if all parameters are present
# If no, exit
if [ $# -ne 3 ]
then
        echo "usage :"
        echo "`basename $0` VERSION SRCDIR ARCHITECTURE"
        exit 0
fi

VERSION=$1
SRCDIR=$(readlink -f $2)
#SRCDIR=$2
ARCH=$3

ORIG_WD=$(pwd)

#define variables
tarname="enGrid_linux"$ARCH"bit_"$VERSION".tar"
gzname="enGrid_linux"$ARCH"bit_"$VERSION".tar.gz"
dirname="enGrid_linux"$ARCH"bit_"$VERSION

#TMP=$(mktemp -d)
TMP=/tmp/$dirname

if [ -e $TMP ]
then
        echo "$TMP already exists."
        echo "rm -rfv $TMP ?(y/n/q)"
	read ans
	case $ans in
	  y|Y|yes) rm -rfv $TMP;;
	  q) exit 0;;
	  *) echo "proceeding without removing";;
	esac
fi

mkdir -p $TMP

#get all dependencies
$SRCDIR/engrid -distbin
mv -v enGrid_bin.tar.gz $TMP

#clean up
cd $SRCDIR
qmake
make clean
cd netgen_svn
qmake
make clean

#add scripts + dependencies
cp -rv $SRCDIR/setup $TMP
cp -rv $SRCDIR/start_engrid $TMP

#change back to SRCDIR + add source files
cd $SRCDIR
cp -rv $SRCDIR/*.h $TMP
cp -rv $SRCDIR/*.cpp $TMP
cp -rv $SRCDIR/*.cxx $TMP
cp -rv $SRCDIR/*.ui $TMP
cp -rv $SRCDIR/licence.txt $TMP
cp -rv  $SRCDIR/resources $TMP
cp -rv  $SRCDIR/engrid.pro $TMP
cp -rv  $SRCDIR/engrid.qrc $TMP
cp -rv  $SRCDIR/math/*.h $TMP
cp -rv  $SRCDIR/netgen_svn/netgen-mesher/netgen $TMP
cp -rv $SRCDIR/netgen_svn/ng.pro $TMP

cd /tmp
tar -czvf $gzname $dirname
mv -v $gzname $ORIG_WD

#exit 0

#extract .tar in ./tmp
#mkdir tmp
#cd tmp
#tar xf ../$tarname

#rename the extracted dir
#mv $SRCDIR $dirname
#tar czf $gzname $dirname
#mv $gzname ..
#cd ..
#rm -rf tmp
#rm $tarname
