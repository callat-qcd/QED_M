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

int main()
{

  cout << setprecision(15);
  rng.Init( 407 );

  int sites[4] = { 4, 4, 4, 4 };
  int volume = 1;
  for (int mu=0; mu<4; ++mu) { volume *= sites[mu]; }
  Lattice lattice(sites);

  double beta = 0.73;
  double xi = 0.73;
  double mass2 = 0.73;

  lattice.MassiveUpdate(beta,mass2); 
  lattice.MassiveUpdate(beta,xi,mass2); 

  double *a = new double[4*volume];
  for (int idx=0; idx<4*volume; ++idx) { a[idx] = lattice.GaugeField(idx); }

  cout << "Before gauge transformation:" << endl;
  cout << "---------------------------:" << endl;
  lattice.PrintObservables(1);

  cout << "Performing multiple gauge fixings followed by Coulomb gauge fixing..." << endl;
  cout << endl;

  for (int n=0; n<100; ++n) {
    lattice.RandomGauge();
    lattice.LorenzGauge();
    lattice.CoulombGauge(3);
    lattice.CoulombGauge(2);
    lattice.CoulombGauge(1);
  }
  lattice.CoulombGauge(0);


  for (int idx=0; idx<4*volume; ++idx) { a[idx] -= lattice.GaugeField(idx); }

  cout << "After gauge transformation:" << endl;
  cout << "--------------------------:" << endl;
  lattice.PrintObservables(1);

  double diff = 0.0;
  for (int idx=0; idx<4*volume; ++idx) { diff += a[idx]*a[idx]; }
  cout << "Difference beween before and after: " << diff;
  cout << endl ;
  cout << endl ;

  delete []a;

  exit(EXIT_SUCCESS);

}
