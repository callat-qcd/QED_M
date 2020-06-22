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
  rng.Init( time(0) );

  int sites[4] = { 32, 16, 8, 4 };
  Lattice lattice(sites);

  cout << "0 : " << lattice.Sites(0) << endl;
  cout << "1 : " << lattice.Sites(1) << endl;
  cout << "2 : " << lattice.Sites(2) << endl;
  cout << "3 : " << lattice.Sites(3) << endl;

  cout << "The following should give an error message:" << endl;
  lattice.Sites(4);

  exit(EXIT_SUCCESS);

}
