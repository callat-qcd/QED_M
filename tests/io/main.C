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
//  rng.Init( time(0) );

  int sites[4] = { 32, 4, 8, 16 };
  Lattice lattice0(sites);
  Lattice lattice1(sites);

  lattice0.MassiveUpdate(2.2,2.3);
  lattice1.MassiveUpdate(2.2,2.3);

  cout << "Lattice 0 : " << lattice0.MeanPlaquetteSq(1,1) << endl;
  cout << "Lattice 1 : " << lattice1.MeanPlaquetteSq(1,1) << endl;

  cout << "Exporting Lattice to file..." << endl;
  lattice0.Export("config0.bin");
  lattice1.Export("config1.bin");

  int ordering[4] = {1,2,3,0}; // Fastest (right) to slowest (left)
  //int ordering[4] = {0,1,2,3}; // Fastest (right) to slowest (left)
  lattice0.Export("config0.dat",ordering);
  lattice1.Export("config1.dat",ordering);

  cout << "Importing Lattices from file..." << endl;
  lattice1.Import("config0.bin");
  lattice0.Import("config1.bin");

  cout << "Lattice 0 : " << lattice0.MeanPlaquetteSq(1,1) << endl;
  cout << "Lattice 1 : " << lattice1.MeanPlaquetteSq(1,1) << endl;

  exit(EXIT_SUCCESS);

}
