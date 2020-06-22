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
#include "block.h"

#define WORKING_DIR "./cfgs/" // The last slash is required!!

#define BETA 1.0  // beta = 1/e^2 = 137/(4pi) 
                  // **Note that we decided to keep beta=1 in photon field generation
                  // and specify its value instead in propagator generation

#define XI 1.0    // Anisotropy factor

#define MASS2 2.3 // photon mass-squared

#define NT 16
#define NX 16
#define NY 16
#define NZ 16

#define NCONF 10 // Number of uncorrelated configurations desired
#define BINARY_EXPORT 1 // Write lattice to a binary file
#define HUMAN_READABLE_EXPORT 1 // Write lattice to a text file
#define MEASURE 1 // Perform various observables

#define RNG_SEED 1997 // Random seed: positive (fixed seed) negative (time based)

void PrintParameters(); // This just outputs the above parameters so they can be stored in
                        // a log file of some sort
string Basename(); // This just generates a standard basename for config files
                   // with all parameters specified in the name

int main()
{

  PrintParameters();

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

  string path(WORKING_DIR);
  string basename( Basename() );

  for (int i=0; i<NCONF; ++i) {

    clock_elapsed = 0.0; 
    clock_elapsed -= clock();

    cout << "Configuration " << i << " :" << endl;

    if (MASS2>0.0) {
      cout << "Massive photon update..." << endl;
      lattice.MassiveUpdate(BETA, XI, MASS2);
    } else {
      cout << "Coulomb gauge-fixed update..." << endl;
      lattice.CoulombGaugeUpdate(BETA, XI);
    }

    if (BINARY_EXPORT) {
      ostringstream filename;
      filename << path << basename;
      filename << "-" << i;
      filename << ".bin";
      lattice.Export( filename.str() );
    }

    if (HUMAN_READABLE_EXPORT) {
      ostringstream filename;
      filename << path  << basename;
      filename << "-" << i;
      filename << ".dat";
      int ordering[4] = {1,2,3,0}; // This makes the ordering of the written field identical to Chroma
      lattice.Export( filename.str() , ordering);
    }

    if (MEASURE) {
      if (i==0) {
        lattice.PrintObservables(1); // Be verbose the first time
      } else {
        lattice.PrintObservables(0); // Then be silent
      }
    }

    clock_elapsed += clock();
    clock_elapsed /= CLOCKS_PER_SEC;
    cout << "Elapsed time / configuration : " << clock_elapsed << "s" << endl;

  }

  exit(EXIT_SUCCESS);

}

void PrintParameters()
{
  cout << "NT = " <<  NT << endl;
  cout << "NX = " <<  NX << endl;
  cout << "NY = " <<  NY << endl;
  cout << "NZ = " <<  NZ << endl;
  cout << "BETA = " <<  BETA << endl;
  cout << "XI = " <<  XI << endl;
  cout << "MASS2 = " <<  MASS2  << endl;
}

string Basename()
{
  ostringstream basename;
  basename << "nt" << NT;
  basename << "nx" << NX;
  basename << "ny" << NY;
  basename << "nz" << NZ;
  basename << "beta" << BETA;
  basename << "xi" << XI;
  basename << "msq" << MASS2;
  return basename.str();
}
