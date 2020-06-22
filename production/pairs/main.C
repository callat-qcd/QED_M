#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <ctime>

using namespace std;
#include "lattice.h"
#include "ranlxd.h"

#define BETA 1.30  // beta = 1/g^2
#define MASS2 0.0 // photon mass-squared

#define NT 8
#define NX 4
#define NY 4
#define NZ 4

#define NCONF 1 // Number of uncorrelated configurations desired
#define BINARY_EXPORT 1
#define HUMAN_READABLE_EXPORT 1
#define MEASURE 1

#define RNG_SEED 1997

// This code will generate field configurations, write them to a file, then perform a random gauge transformation on the field and write them to a file a second time.
// This code is useful for checking the U(1) gauge-invariance of observables written in chroma_51.
int main()
{

  cout << setprecision(15);
  if (RNG_SEED<0) {
    rng.Init( time(0) );
  } else {
    rng.Init( RNG_SEED );
  }

  int sites[4] = { NT, NX, NY, NZ };

  Lattice lattice(sites);

  clock_t clock_base = clock();
  double clock_elapsed;

  string basename("cfgs/");

  for (int i=0; i<NCONF; ++i) {

    clock_elapsed = 0.0; 
    clock_elapsed -= clock();

    if (MASS2>0.0) {
      cout << "Massive photon update with MASS2 = " << MASS2 << endl;
      lattice.MassiveUpdate(BETA, MASS2);
    } else {
      cout << "Coulomb gauge-fixed update" << endl;
      lattice.CoulombGaugeUpdate(BETA);
    }


    if (BINARY_EXPORT) {
      ostringstream filename;
      filename << basename;
      filename << "nt" << NT;
      filename << "nx" << NX;
      filename << "ny" << NY;
      filename << "nz" << NZ;
      filename << "beta" << BETA;
      filename << "msq" << MASS2;
      filename << "-" << i;
      filename << "-A.bin";
      lattice.Export( filename.str() );
    }

    if (HUMAN_READABLE_EXPORT) {
      ostringstream filename;
      filename << basename;
      filename << "nt" << NT;
      filename << "nx" << NX;
      filename << "ny" << NY;
      filename << "nz" << NZ;
      filename << "beta" << BETA;
      filename << "msq" << MASS2;
      filename << "-" << i;
      filename << "-A.dat";
      int ordering[4] = {1,2,3,0};
      lattice.Export( filename.str() , ordering);
    }


    lattice.RandomGauge();

    if (BINARY_EXPORT) {
      ostringstream filename;
      filename << basename;
      filename << "nt" << NT;
      filename << "nx" << NX;
      filename << "ny" << NY;
      filename << "nz" << NZ;
      filename << "beta" << BETA;
      filename << "msq" << MASS2;
      filename << "-" << i;
      filename << "-B.bin";
      lattice.Export( filename.str() );
    }

    if (HUMAN_READABLE_EXPORT) {
      ostringstream filename;
      filename << basename;
      filename << "nt" << NT;
      filename << "nx" << NX;
      filename << "ny" << NY;
      filename << "nz" << NZ;
      filename << "beta" << BETA;
      filename << "msq" << MASS2;
      filename << "-" << i;
      filename << "-B.dat";
      int ordering[4] = {1,2,3,0};
      lattice.Export( filename.str() , ordering);
    }

    clock_elapsed += clock();
    clock_elapsed /= CLOCKS_PER_SEC;
    cout << "Elapsed time / configuration : " << clock_elapsed << "s" << endl;

  }

  exit(EXIT_SUCCESS);

}
