#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <ctime>

using namespace std;
#include "block.h"
#include "lattice.h"
#include "ranlxd.h"

#define BETA 2.3

int main()
{

  cout << setprecision(15);
  rng.Init( 971 );

  int big_sites[4] = { 32, 32, 32, 32 };
  int small_sites[4] = { 16, 16, 16, 16 };

  Lattice big_lattice(big_sites);
  Lattice small_lattice(small_sites);

  big_lattice.CoulombGaugeUpdate(BETA);
  cout << big_lattice.MeanPlaquetteSq(1,1) << endl;

  for (int level=0; level<4; ++level) {
    BlockPAC_CS(big_lattice,small_lattice,level);
    cout << small_lattice.MeanPlaquetteSq(1,1) << endl;
  }


  exit(EXIT_SUCCESS);

}
