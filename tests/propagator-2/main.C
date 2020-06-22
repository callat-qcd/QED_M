#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <ctime>

using namespace std;
#include "propagator.h"

int main()
{

  cout << setprecision(15);

  int sites[4] = { 4, 6, 8, 10 };
  int volume = 1;
  for (int mu=0; mu<4; ++mu) { volume *= sites[mu]; }

  complex<double> inner[4][4];
  for (int mu=0; mu<4; ++mu)
  for (int nu=0; nu<4; ++nu) {
    inner[mu][nu] = 0.0;
  }

  int s[4];
  complex<double> a, b, c;
  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) { 

    Propagator propagator(s,sites);

    for (int mu=0; mu<4; ++mu)
    for (int nu=0; nu<4; ++nu) {
      for (int rho=0; rho<4; ++rho) {
        a = propagator.SimilarityTransformationMatrix(rho, mu);
        b = propagator.SimilarityTransformationMatrix(rho, nu);
        inner[mu][nu] +=conj( a ) * b;
      }
    }

  }

  for (int mu=0; mu<4; ++mu)
  for (int nu=0; nu<4; ++nu) {
    inner[mu][nu] /= volume;
    if ( abs( inner[mu][nu] ) < 1.0e-15 ) { inner[mu][nu] = 0.0; }
    cout << mu << " " << nu << " " << inner[mu][nu] << endl;
  }

  exit(EXIT_SUCCESS);

}
