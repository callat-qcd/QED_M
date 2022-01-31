#!/bin/bash

LIB=noncompact_u1
SRCDIR=$PWD
BUILDDIR=$PWD
SRCDIR2=$SRCDIR/src
BUILDDIR2=$BUILDDIR/objs
FFTWDIR=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-5/fftw-3.3.8-azzdjlzx2j6dpqvzdir2nwvxypohyfq4/lib

CC=gcc
CXX=g++
AR=ar

LDFLAGS=""
CXXFLAGS="-O2"
ARFLAGS="ruv"

INCLUDE_FLAGS="-I$SRCDIR/include -I$FFTWDIR/include"

if [ -d $BUILDDIR2 ]; then
  rm -r $BUILDDIR2
fi

mkdir $BUILDDIR2
cd $BUILDDIR2

$CXX $CXXFLAGS $LDFLAGS $INCLUDE_FLAGS -c $SRCDIR2/*.C
$AR $ARFLAGS $BUILDDIR/$LIB.a $BUILDDIR2/*.o
ranlib  $BUILDDIR/$LIB.a
