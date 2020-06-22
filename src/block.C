#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;
#include "block.h"

void BlockPAC_CS(Lattice & big, Lattice & small)
{
  // do all four levels of averaging in the p-cell orthogonal to mu
  BlockPAC_CS(big, small, 99);
}

void BlockPAC_CS(Lattice & big, Lattice & small, int level)
{

  // This blocking routine follows the method of
  // Aoki, et al., Phys. Rev. D 86, 034507 (2012)
  // and given explicitly in Eq. 8
  // Averaging is done over links within the p-cell orthogonal to mu
  // where level specifies the maximal p included in the average.

  int big_sites[4];
  int small_sites[4];

  // Confirm that the lattice dimensions are double
  for (int mu=0; mu<4;++mu) {
    big_sites[mu] = big.Sites(mu);
    small_sites[mu] = small.Sites(mu);
    if ( big_sites[mu] - 2*small_sites[mu] != 0 ) {
      cout << "Big lattice must have dimensions twice that of small lattice." << endl;
      exit(EXIT_FAILURE);
    }
  }

  // Now for the fun part
  int b[4];
  int s[4];
  double a;
  double a0;
  double a1;
  double a2;
  double a3;

  for ( s[0]=0; s[0]<small_sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<small_sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<small_sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<small_sites[3]; ++s[3] ) {
    for (int mu=0; mu<4; ++mu) { b[mu] = 2*s[mu]; }
    for (int mu=0; mu<4; ++mu) {

      // 0-cell
      a0 = 0.0;

      a0 += big.GaugeField(mu,b);
      b[mu]+=1;
      a0 += big.GaugeField(mu,b);
      b[mu]-=1;

      // 1-cell
      for (int nu=0; nu<4; ++nu)
      if (nu!=mu) {

        a1 = 0.0;

        a1 += big.GaugeField(nu,b);
        b[nu]+=1;
        a1 += big.GaugeField(mu,b);
        b[mu]+=1;
        a1 += big.GaugeField(mu,b);
        b[mu]+=1;
        b[nu]-=1;
        a1 -= big.GaugeField(nu,b);
        b[mu]-=2;

        b[nu]-=1;
        a1 -= big.GaugeField(nu,b);
        a1 += big.GaugeField(mu,b);
        b[mu]+=1;
        a1 += big.GaugeField(mu,b);
        b[mu]+=1;
        a1 += big.GaugeField(nu,b);
        b[nu]+=1;
        b[mu]-=2;

      }

      // 2-cell
      for (int nu=0; nu<4; ++nu)
      for (int rho=0; rho<4; ++rho)
      if ( (nu!=mu)&&(rho!=mu)&&(nu!=rho) ) {

        a2 = 0.0;

        a2 += big.GaugeField(rho,b);
        b[rho]+=1;
        a2 += big.GaugeField(nu,b);
        b[nu]+=1;
        a2 += big.GaugeField(mu,b);
        b[mu]+=1;
        a2 += big.GaugeField(mu,b);
        b[mu]+=1;
        b[nu]-=1;
        a2 -= big.GaugeField(nu,b);
        b[rho]-=1;
        a2 -= big.GaugeField(rho,b);
        b[mu]-=2;

        a2 += big.GaugeField(rho,b);
        b[rho]+=1;
        b[nu]-=1;
        a2 -= big.GaugeField(nu,b);
        a2 += big.GaugeField(mu,b);
        b[mu]+=1;
        a2 += big.GaugeField(mu,b);
        b[mu]+=1;
        a2 += big.GaugeField(nu,b);
        b[nu]+=1;
        b[rho]-=1;
        a2 -= big.GaugeField(rho,b);
        b[mu]-=2;

        b[rho]-=1;
        a2 -= big.GaugeField(rho,b);
        a2 += big.GaugeField(nu,b);
        b[nu]+=1;
        a2 += big.GaugeField(mu,b);
        b[mu]+=1;
        a2 += big.GaugeField(mu,b);
        b[mu]+=1;
        b[nu]-=1;
        a2 -= big.GaugeField(nu,b);
        a2 += big.GaugeField(rho,b);
        b[rho]+=1;
        b[mu]-=2;

        b[rho]-=1;
        a2 -= big.GaugeField(rho,b);
        b[nu]-=1;
        a2 -= big.GaugeField(nu,b);
        a2 += big.GaugeField(mu,b);
        b[mu]+=1;
        a2 += big.GaugeField(mu,b);
        b[mu]+=1;
        a2 += big.GaugeField(nu,b);
        b[nu]+=1;
        a2 += big.GaugeField(rho,b);
        b[rho]+=1;
        b[mu]-=2;

      }

      // 3-cell
      for (int nu=0; nu<4; ++nu)
      for (int rho=0; rho<4; ++rho)
      for (int sigma=0; sigma<4; ++sigma)
      if ( (nu!=mu)&&(rho!=mu)&&(sigma!=mu)&&(rho!=nu)&&(sigma!=nu)&&(sigma!=rho) ) {

        a3 = 0.0;

        a3 += big.GaugeField(sigma,b);
        b[sigma]+=1;
        a3 += big.GaugeField(rho,b);
        b[rho]+=1;
        a3 += big.GaugeField(nu,b);
        b[nu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        b[nu]-=1;
        a3 -= big.GaugeField(nu,b);
        b[rho]-=1;
        a3 -= big.GaugeField(rho,b);
        b[sigma]-=1;
        a3 -= big.GaugeField(sigma,b);
        b[mu]-=2;

        a3 += big.GaugeField(sigma,b);
        b[sigma]+=1;
        a3 += big.GaugeField(rho,b);
        b[rho]+=1;
        b[nu]-=1;
        a3 -= big.GaugeField(nu,b);
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(nu,b);
        b[nu]+=1;
        b[rho]-=1;
        a3 -= big.GaugeField(rho,b);
        b[sigma]-=1;
        a3 -= big.GaugeField(sigma,b);
        b[mu]-=2;

        a3 += big.GaugeField(sigma,b);
        b[sigma]+=1;
        b[rho]-=1;
        a3 -= big.GaugeField(rho,b);
        a3 += big.GaugeField(nu,b);
        b[nu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        b[nu]-=1;
        a3 -= big.GaugeField(nu,b);
        a3 += big.GaugeField(rho,b);
        b[rho]+=1;
        b[sigma]-=1;
        a3 -= big.GaugeField(sigma,b);
        b[mu]-=2;

        a3 += big.GaugeField(sigma,b);
        b[sigma]+=1;
        b[rho]-=1;
        a3 -= big.GaugeField(rho,b);
        b[nu]-=1;
        a3 -= big.GaugeField(nu,b);
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(nu,b);
        b[nu]+=1;
        a3 += big.GaugeField(rho,b);
        b[rho]+=1;
        b[sigma]-=1;
        a3 -= big.GaugeField(sigma,b);
        b[mu]-=2;

        b[sigma]-=1;
        a3 -= big.GaugeField(sigma,b);
        a3 += big.GaugeField(rho,b);
        b[rho]+=1;
        a3 += big.GaugeField(nu,b);
        b[nu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        b[nu]-=1;
        a3 -= big.GaugeField(nu,b);
        b[rho]-=1;
        a3 -= big.GaugeField(rho,b);
        a3 += big.GaugeField(sigma,b);
        b[sigma]+=1;
        b[mu]-=2;

        b[sigma]-=1;
        a3 -= big.GaugeField(sigma,b);
        a3 += big.GaugeField(rho,b);
        b[rho]+=1;
        b[nu]-=1;
        a3 -= big.GaugeField(nu,b);
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(nu,b);
        b[nu]+=1;
        b[rho]-=1;
        a3 -= big.GaugeField(rho,b);
        a3 += big.GaugeField(sigma,b);
        b[sigma]+=1;
        b[mu]-=2;

        b[sigma]-=1;
        a3 -= big.GaugeField(sigma,b);
        b[rho]-=1;
        a3 -= big.GaugeField(rho,b);
        a3 += big.GaugeField(nu,b);
        b[nu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        b[nu]-=1;
        a3 -= big.GaugeField(nu,b);
        a3 += big.GaugeField(rho,b);
        b[rho]+=1;
        a3 += big.GaugeField(sigma,b);
        b[sigma]+=1;
        b[mu]-=2;

        b[sigma]-=1;
        a3 -= big.GaugeField(sigma,b);
        b[rho]-=1;
        a3 -= big.GaugeField(rho,b);
        b[nu]-=1;
        a3 -= big.GaugeField(nu,b);
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(mu,b);
        b[mu]+=1;
        a3 += big.GaugeField(nu,b);
        b[nu]+=1;
        a3 += big.GaugeField(rho,b);
        b[rho]+=1;
        a3 += big.GaugeField(sigma,b);
        b[sigma]+=1;
        b[mu]-=2;

      }

      // The average should be thought of as an average over links
      // extending in the mu-direction. These links can be identifed with the
      // sites (1), links (6), faces (12) and cubes (8) of the 3-cube orthogonal
      // to mu. In the above sums, these appear with degeneracies (1,1,2,6).
      // Thus the appropriate average is given by:
//      a = a0 + a1 + a2/2.0 + a3/6.0;
//      a /= 27.0; 

      double p = 1.0;
      a = a0/1.0;
      if (level>0) { a += a1; p+=6.0; }
      if (level>1) { a += a2/2.0; p+=12.0; }
      if (level>2) { a += a3/6.0; p+=8.0; }
      a /= p;

      small.GaugeField(mu,s) = a;
 
    }
  }
}

