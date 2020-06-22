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

#define BETA 2.1  // beta = 1/e^2 = 137/(4pi)
#define XI 1.7
//#define MASS2 1.8 // photon mass-squared
#define MASS2 0.0 // photon mass-squared

#define NT 32
#define NX 16
#define NY 8
#define NZ 4

#define NCONF 1000000 // Number of uncorrelated configurations desired
#define BINARY_EXPORT 0 // Write lattice to a binary file
#define HUMAN_READABLE_EXPORT 0 // Write lattice to a text file
#define MEASURE 1 // Perform various observables

#define RNG_SEED 1997 // Random seed: positive (fixed seed) negative (time based)
#define VERBOSE 0 // When on, an explaination of the observables is given
#define BLOCK 0 // Block gauge-field 

void PrintParameters();
string Basename();

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
  int half_sites[4];
  for (int mu=0; mu<4; ++mu) { half_sites[mu] = sites[mu]/2; }

  Lattice lattice(sites);
  Lattice* half_lattice;
  if (BLOCK) half_lattice = new Lattice(half_sites);

  clock_t clock_base = clock();
  double clock_elapsed;

  string path("cfgs/");
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

    if (0) { lattice.RandomGauge(); }

    if (BLOCK) {
      cout << "Blocking..." << endl;
      BlockPAC_CS(lattice,*half_lattice);
    }

    if (BINARY_EXPORT) {
      ostringstream filename;
      filename << path << basename;
      filename << "-" << i;
      filename << ".bin";
      lattice.Export( filename.str() );
    }

    if (BLOCK)
    if (BINARY_EXPORT) {
      ostringstream filename;
      filename << path << "half-" <<  basename;
      filename << "-" << i;
      filename << ".bin";
      half_lattice->Export( filename.str() );
    }

    if (HUMAN_READABLE_EXPORT) {
      ostringstream filename;
      filename << path  << basename;
      filename << "-" << i;
      filename << ".dat";
      int ordering[4] = {1,2,3,0};
      lattice.Export( filename.str() , ordering);
    }

    if (BLOCK)
    if (HUMAN_READABLE_EXPORT) {
      ostringstream filename;
      filename << path  << "half-" <<  basename;
      filename << "-" << i;
      filename << ".dat";
      int ordering[4] = {1,2,3,0};
      half_lattice->Export( filename.str() , ordering);
    }

    if (0) {
      lattice.LorenzGauge();
      lattice.CoulombGauge(0);
    }

    if (MEASURE) {
      if (i==0) {
        lattice.PrintObservables(1);
        if (BLOCK) half_lattice->PrintObservables(1);
      } else {
        lattice.PrintObservables(0);
        if (BLOCK) half_lattice->PrintObservables(0);
      }
    }

//  if (0) {
//    int r[4] = {0,0,0,0};
//    for (r[3]=0; r[3]<sites[3]; ++r[3]) {
//      complex<double> poly = lattice.MeanPolyakovLoop(0,1,r);
//      cout << real(poly) << " ";
//    }
//    cout << endl;
//  }

//  if (0) {
//    int r[4] = {0,0,0,0};
//    complex<double> poly;
//    for (r[3]=0; r[3]<sites[3]; ++r[3])
//    for (r[2]=0; r[2]<sites[2]; ++r[2])
//    for (r[1]=0; r[1]<sites[1]; ++r[1])
//    for (r[0]=0; r[0]<sites[0]; ++r[0])
//    for (int mu=0; mu<4; ++mu)
//    for (int nu=0; nu<4; ++nu) {
//      poly = lattice.MeanPolyakovLoop(mu,nu,r);
//      cout << real(poly) << " ";
//    }
//    cout << endl;
//  }

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
