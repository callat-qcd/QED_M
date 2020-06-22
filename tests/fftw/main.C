#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <ctime>

using namespace std;
#include "ranlxd.h"
#include "fftw3.h"

int main()
{

  cout << setprecision(15);
  rng.Init( time(0) );

  int sites[4] = { 8, 8, 8, 8 };
  int vol = 1;

  for (int mu=0; mu<4; ++mu) { vol *= sites[mu]; }
  cout << "Volume : " << vol << endl;

  fftw_complex *in = new fftw_complex [vol];
  fftw_complex *out = new fftw_complex [vol];
  fftw_complex *in_again = new fftw_complex [vol];

  fftw_plan dft_forward = fftw_plan_dft(4, sites, in, out, FFTW_FORWARD, FFTW_PATIENT);
  fftw_plan dft_backward = fftw_plan_dft(4, sites, out, in_again, FFTW_BACKWARD, FFTW_PATIENT);

  for (int i=0; i<vol; ++i) {
    in[i][0] = rng.Normal(0.0, 1.0);
    in[i][1] = rng.Normal(0.0, 1.0);
    out[i][0] = 0.0;
    out[i][1] = 0.0;
    in_again[i][0] = 0.0;
    in_again[i][1] = 0.0;
  }

  fftw_execute(dft_forward);
  fftw_execute(dft_backward);

  double err = 0.0;
  double diff = 0.0;
  for (int i=0; i<vol; ++i) {
    diff = in_again[i][0]/vol - in[i][0];
    err += diff*diff;
    diff = in_again[i][1]/vol - in[i][1];
    err += diff*diff;
  }
  cout << "|initial - final|^2 : " << err << endl;

  fftw_destroy_plan(dft_forward);
  fftw_destroy_plan(dft_backward);
  fftw_cleanup();


  delete[] in;
  delete[] out;
  delete[] in_again;

  exit(EXIT_SUCCESS);

}
