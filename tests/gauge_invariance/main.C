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

  int sites[4] = { 8, 8, 8, 8 };
  Lattice lattice(sites);

  cout << "Before gauge transformation:" << endl;
  cout << "---------------------------:" << endl;
  lattice.PrintObservables(1);

  cout << "Performing random gauge transformation..." << endl;
  cout << endl;
  lattice.RandomGauge();

  cout << "After gauge transformation:" << endl;
  cout << "--------------------------:" << endl;
  lattice.PrintObservables(1);


  exit(EXIT_SUCCESS);

}
