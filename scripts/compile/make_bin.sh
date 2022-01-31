#!/bin/bash

LIB=noncompact_u1
SRCDIR=../../
BUILDDIR=../../
FFTWDIR=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-5/fftw-3.3.8-azzdjlzx2j6dpqvzdir2nwvxypohyfq4/lib

CXX=g++

#LDFLAGS="-L$FFTWDIR/lib -lm  -Wl,-framework -Wl,Accelerate -lfftw3"# mac
LDFLAGS="-L$FFTWDIR/lib -lm -lfftw3"
CXXFLAGS="-O3"
ARFLAGS=

INCLUDE_FLAGS="-I$SRCDIR/include -I$FFTWDIR/include"

MAIN=$1
$CXX $CXXFLAGS $MAIN $BUILDDIR/$LIB.a $LDFLAGS $INCLUDE_FLAGS
