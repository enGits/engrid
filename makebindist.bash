#!/bin/bash
set -e -x -u

# Check if all parameters are present
# If no, exit
if [ $# -ne 3 ]
then
        echo "usage :"
        echo "`basename $0` VERSION SRCDIR ARCHITECTURE"
	echo "You must run this script from SRCDIR/.."
        exit 0
fi

VERSION=$1
#SRCDIR=$(readlink -f $2)
SRCDIR=$2
ARCH=$3

ORIG_WD=$(pwd)

#define variables
tarname="enGrid_linux"$ARCH"bit_"$VERSION".tar"
gzname="enGrid_linux"$ARCH"bit_"$VERSION".tar.gz"
dirname="enGrid_linux"$ARCH"bit_"$VERSION

saferemove()
{
	TARGET=$(readlink -f $1)
	if [ -e $TARGET ]
	then
	        echo "$TARGET already exists."
	        echo "rm -v $TARGET ?(y/n/q)"
		read ans
		case $ans in
		  y|Y|yes) rm -v $TARGET;;
		  q) exit 0;;
		  *) echo "proceeding without removing";;
		esac
	fi
}

saferemove_recursive()
{
	TARGET=$(readlink -f $1)
	if [ -e $TARGET ]
	then
	        echo "$TARGET already exists."
	        echo "rm -rfv $TARGET ?(y/n/q)"
		read ans
		case $ans in
		  y|Y|yes) rm -rfv $TARGET;;
		  q) exit 0;;
		  *) echo "proceeding without removing";;
		esac
	fi
}

saferemove_recursive tmp
saferemove $tarname

#TMP=$(mktemp -d)
#TMP=/tmp/toto

#if [ -e $TMP ]
#then
#        echo "$TMP already exists."
#        echo "rm -rfv $TMP ?(y/n/q)"
#	read ans
#	case $ans in
#	  y|Y|yes) rm -rfv $TMP;;
#	  q) exit 0;;
#	  *) echo "proceeding without removing";;
#	esac
#fi

#mkdir -p $TMP

#get all dependencies
cd $SRCDIR
./engrid -distbin
#mv -v enGrid_bin.tar.gz $TMP

#clean up
#cd $SRCDIR
qmake
make clean
cd netgen_svn
qmake
make clean
cd ../..

#add scripts + dependencies
pwd
tar -f $tarname -r $SRCDIR/setup
#tar -f $tarname -r $SRCDIR/start_engrid
tar -f $tarname -r $SRCDIR/enGrid_bin.tar.gz
tar -f $tarname -r $SRCDIR/README

#change back to SRCDIR + add source files
tar -f $tarname -r $SRCDIR/*.h
tar -f $tarname -r $SRCDIR/*.cpp
tar -f $tarname -r $SRCDIR/*.cxx
tar -f $tarname -r $SRCDIR/*.ui
tar -f $tarname -r $SRCDIR/licence.txt
tar -f $tarname -r $SRCDIR/resources
tar -f $tarname -r $SRCDIR/engrid.pro
tar -f $tarname -r $SRCDIR/engrid.qrc
tar -f $tarname -r $SRCDIR/math/*.h
tar -f $tarname -r $SRCDIR/netgen_svn/netgen-mesher/netgen
tar -f $tarname -r $SRCDIR/netgen_svn/ng.pro

#cd /tmp
#tar -czvf $gzname $TMP
#mv -v $gzname $ORIG_WD

#exit 0

#extract .tar in ./tmp
pwd
mkdir -p tmp
cd tmp
tar xf ../$tarname
pwd

#rename the extracted dir
mv $SRCDIR $dirname
tar czf $gzname $dirname
mv $gzname ..
cd ..
saferemove_recursive tmp
saferemove $tarname
