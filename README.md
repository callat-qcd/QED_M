# noncompact_u1
C/C++ code to generate electro-quenched U1 fields for LQCD.  
Original author: Michael G. Endres.  
Porting to new repo: Ben Hoerz and Andre Walker-Loud

- [Installation](#installation)
- [Debugging](#debugging)
- [Acknowledgements](#acknowledgements)


## Installation

### 1) To create the library edit the FFTWLIB flag in the files:  
`scripts/compile/make_lib.sh`  
and  
`scripts/compile/make_bin.sh`  
to point to where your FFTW library is located.

Note that these files (along with scripts/compile/make_clean.sh) are my lazy alternative to using a Makefile; they will be used to compile the library and test/production programs.

### 2) Next run:  
`./scripts/compile/make_lib.sh`

This should create the library file called noncompact_u1.a (as well as an obj directory containing some *.o files which were created during the comilation)

### 3) To compile one of the test programs, ENTER the test program directory, e.g.,  
`cd tests/gauge_invariance/`  
and then run the script:  
`../../scripts/compile/make_bin.sh main.C`  
This will create your new binary called a.out

### 4 To run the test program, type:

./a.out

### 5) Production code is compiled and run similarly to test code


## Debugging
====DEBUGGING AND MEMORY USE====

A useful tool for detecting memory leaks/use and profiling is valgrind (no recompilation of your code is required!):

http://valgrind.org/

Some of the tools included are:

### 1) To profile memory use, run:  
`valgrind --tool=massif ./a.out`  

### 2) To detect memory leaks, access of unallocated memory locations, etc.,, run:  
`valgrind --leak-check=yes ./a.out`

### 3) Profiling can be done with  
`valgrind --tool=callgrind ./a.out`  
`callgrind_annotate [options] callgrind.out.<pid>`



## Acknowledgements
If you use this code, please cite

Massive photons: an infrared regularization scheme for lattice QCD+QED  
M. G. Endres, A. Shindler, B. C. Tiburzi, A. Walker-Loud  
[Phys.Rev.Lett. 117 (2016) 7, 072002](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.117.072002) [[arXiv:1507.08916](https://arxiv.org/abs/1507.08916)]
