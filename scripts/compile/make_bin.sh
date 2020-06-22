#!/bin/bash

LIB=noncompact_u1
SRCDIR=../../
BUILDDIR=../../
FFTWDIR=$HOME/Documents/libraries/fftw-3.3.3

CXX=g++

LDFLAGS="-L$FFTWDIR/lib -lm  -Wl,-framework -Wl,Accelerate -lfftw3"
CXXFLAGS="-O2"
ARFLAGS=

INCLUDE_FLAGS="-I$SRCDIR/include -I$FFTWDIR/include"

MAIN=$1
$CXX $CXXFLAGS $MAIN $BUILDDIR/$LIB.a $LDFLAGS $INCLUDE_FLAGS
