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
  int medium_sites[4] = { 16, 16, 16, 16 };
  int small_sites[4] = { 8, 8, 8, 8 };
  int tiny_sites[4] = { 4, 4, 4, 4 };
  int miniscule_sites[4] = { 2, 2, 2, 2 };

  Lattice big_lattice(big_sites);
  Lattice medium_lattice(medium_sites);
  Lattice small_lattice(small_sites);
  Lattice tiny_lattice(tiny_sites);
  Lattice miniscule_lattice(miniscule_sites);

  int two_volume = 2*2*2*2;
  double big_beta = BETA;
  double medium_beta = big_beta*two_volume;
  double small_beta = medium_beta*two_volume;
  double tiny_beta = small_beta*two_volume;
  double miniscule_beta = tiny_beta*two_volume;

  big_lattice.CoulombGaugeUpdate(big_beta);
  medium_lattice.CoulombGaugeUpdate(medium_beta);
  small_lattice.CoulombGaugeUpdate(small_beta);
  tiny_lattice.CoulombGaugeUpdate(tiny_beta);
  miniscule_lattice.CoulombGaugeUpdate(miniscule_beta);

  cout << "Before blocking:" << endl;
  cout << "---------------------------:" << endl;

  cout << big_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << medium_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << small_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << tiny_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << miniscule_lattice.MeanPlaquetteSq(1,1) << endl;

  cout << "Before gauge transformation:" << endl;
  cout << "---------------------------:" << endl;

  BlockPAC_CS(big_lattice,medium_lattice);
  BlockPAC_CS(medium_lattice,small_lattice);
  BlockPAC_CS(small_lattice,tiny_lattice);
  BlockPAC_CS(tiny_lattice,miniscule_lattice);

  cout << big_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << medium_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << small_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << tiny_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << miniscule_lattice.MeanPlaquetteSq(1,1) << endl;

  cout << "Performing random gauge transformation..." << endl;
//  cout << endl;
  big_lattice.RandomGauge();

  cout << "After gauge transformation:" << endl;
  cout << "--------------------------:" << endl;

  BlockPAC_CS(big_lattice,medium_lattice);
  BlockPAC_CS(medium_lattice,small_lattice);
  BlockPAC_CS(small_lattice,tiny_lattice);
  BlockPAC_CS(tiny_lattice,miniscule_lattice);

  cout << big_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << medium_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << small_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << tiny_lattice.MeanPlaquetteSq(1,1) << endl;
  cout << miniscule_lattice.MeanPlaquetteSq(1,1) << endl;

  exit(EXIT_SUCCESS);

}
