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

  if (0) {
    int s[4];
    int idx=0;
    for ( s[0]=0; s[0]<sites[0]; ++s[0] )
    for ( s[1]=0; s[1]<sites[1]; ++s[1] )
    for ( s[2]=0; s[2]<sites[2]; ++s[2] )
    for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
      cout << s[0] << " ";
      cout << s[1] << " ";
      cout << s[2] << " ";
      cout << s[3] << " ";
      cout << lattice.SiteIndex(s) << " ";
      cout << endl;
      idx++;
    }
  }

  if (1) {
    int s[4];
    int idx=0;
    for ( s[0]=0; s[0]<sites[0]; ++s[0] )
    for ( s[1]=0; s[1]<sites[1]; ++s[1] )
    for ( s[2]=0; s[2]<sites[2]; ++s[2] )
    for ( s[3]=0; s[3]<sites[3]; ++s[3] )
    for (int mu=0; mu<4; ++mu) {
      cout << s[0] << " ";
      cout << s[1] << " ";
      cout << s[2] << " ";
      cout << s[3] << " ";
      cout << mu   << " ";
      cout << lattice.GaugeIndex(mu,s) << " ";
      cout << lattice.GaugeField(mu,s) << " ";
      cout << lattice.GaugeField(idx)  << " ";
      cout << endl;
      idx++;
    }
  }

  exit(EXIT_SUCCESS);

}
