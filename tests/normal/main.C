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
  rng.Init( 407 );

  complex<double> z;
  double x = 0.0;
  double y = 0.0;
  double xx = 0.0;
  double yy = 0.0;
  double xy = 0.0;

  int configs = 100000000;
  cout << "Generating " << configs << " complex random numbers z = x+iy..." << endl;

  for (int i=0; i<configs; ++i) {
    z = rng.Normal(1.0);
    x  += z.real();
    y  += z.imag();
    xx += z.real() * z.real();
    yy += z.imag() * z.imag();
    xy += z.real() * z.imag();
  }

  x  /= configs;  
  y  /= configs;  
  xx /= configs;  
  yy /= configs;  
  xy /= configs;  

  cout << "x = " << x << " " << endl;
  cout << "y = " << y << " " << endl;
  cout << "<x>^2 - <x>^2 = " << xx-x*x << " " << endl;
  cout << "<y^2> - <y>^2 = " << yy-y*y << " " << endl;
  cout << "<xy> - <x><y> = " << xy-x*y << " " << endl;

  exit(EXIT_SUCCESS);

}
