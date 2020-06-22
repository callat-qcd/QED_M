#include <iostream>
#include <stdlib.h>
#include <math.h>
using namespace std;
#include "propagator.h"

Propagator::Propagator(int* n, int* sites)
: PI(3.141592653589793),
  I(0,1)
{

  // Complex momentum vectors
  p_sq = 0.0;
  p3vec_sq = 0.0;
  for (int mu=0; mu<4; ++mu) {
    double q = 2.0*PI*n[mu]/sites[mu];
    p[mu]  = 2.0 * exp( I*q/2.0 ) * sin(q/2.0);
    p_sq += 4.0*sin(q/2)*sin(q/2);
    if (mu>0) { p3vec_sq += 4.0*sin(q/2)*sin(q/2); }
  }

  for (int mu=0; mu<4; ++mu)
  for (int nu=0; nu<4; ++nu) {
    s[mu][nu] = 0.0;
  }

  if ( n[1]==0&&n[2]==0&&n[3]==0 ) {
    for (int mu=0; mu<4; ++mu) {
      s[mu][mu] = 1.0;
    }
  } else {

    // Column 0
    for (int mu=0; mu<4; ++mu) { s[mu][0] = p[mu]; }


    // Column 1
    s[0][1] = p[0]*conj(p[0]) - p_sq;
    for (int j=1; j<4; ++j) { s[j][1] = p[j]*conj(p[0]); }

    // Column 2
    s[0][2]  = 0.0;
    if ( n[2]==0 && n[3]==0 ) {
      s[1][2]  =  conj( p[2] );
      s[2][2]  = -conj( p[1] );
      s[3][2]  = 0.0;
    } else {
      s[1][2]  = 0.0;
      s[2][2]  =  conj( p[3] );
      s[3][2]  = -conj( p[2] );
    }

    // Column 3
    s[0][3]  =  0.0;
    s[1][3]  =  p[3]*s[2][2] - p[2]*s[3][2];
    s[2][3]  = -p[3]*s[1][2] + p[1]*s[3][2];
    s[3][3]  =  p[2]*s[1][2] - p[1]*s[2][2];
    for (int mu=0; mu<4; ++mu) {
      s[mu][3] = conj( s[mu][3] );
    }

    // Normalize column vectors
    double norm;
    for (int nu=0; nu<4; ++nu) {
        norm = 0.0;
      for (int mu=0; mu<4; ++mu) {
        norm += real( conj( s[mu][nu] ) * s[mu][nu] );
      }
      norm = sqrt(norm);
      for (int mu=0; mu<4; ++mu) {
        s[mu][nu] /= norm;
      }
    }

  }

}

Propagator::~Propagator()
{
}

double Propagator::ThreeMomentumSquared() { return p3vec_sq; };

double Propagator::FourMomentumSquared() { return p_sq; };

double Propagator::FourMomentumSquared(double xi) { return xi*p_sq + (1.0/xi - xi)*p3vec_sq;  };

complex<double>& Propagator::SimilarityTransformationMatrix(int mu, int nu) {return s[mu][nu]; }

void Propagator::PrintMomentum()
{
  for (int mu=0; mu<4; ++mu) {
    cout << mu << " : " << p[mu] << endl;
  }
}

void Propagator::PrintSimilarityTransformationMatrix()
{
  for (int mu=0; mu<4; ++mu)
  for (int nu=0; nu<4; ++nu) {
    cout << mu << "," << nu << " : " << s[mu][nu] << endl;
  }
}


complex<double> Propagator::Momentum(int mu)
{
  return p[mu];
}
