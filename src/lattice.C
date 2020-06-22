#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;
#include "lattice.h"
#include "ranlxd.h"

Lattice::Lattice(int* my_sites)
{

  // sites = { NT, NX, NY, NZ }
  volume = 1;
  for (int mu=0; mu<4; ++mu) {
     sites[mu] = my_sites[mu];
     volume *= sites[mu];
  }

  gauge_field_size = 4*volume;
  gauge_field = new double [gauge_field_size];

  propagator = new Propagator* [volume];
  int s[4];
  int idx;
  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    idx = SiteIndex(s);
    propagator[idx] = new Propagator(s, sites); // Create an instance for every momentum vector
  }

  // You will have to read the FFTW documentation about the following
  // I had to  perform a "strided" FFT due to the layout of the gauge field
  // ( the way thisngs are set up, mu is the fastest running index)
  in = (fftw_complex *) fftw_malloc( gauge_field_size * sizeof(fftw_complex) );
  out = (fftw_complex *) fftw_malloc( gauge_field_size * sizeof(fftw_complex) );

  int rank = 4;
  int howmany = 4;
  int stride = 4;
  int dist = 1;
  cout << "Creating fftw plan (this could take a while)........" << flush;
  dft_forward  = fftw_plan_many_dft(rank, sites, howmany, in, NULL, stride, dist, out, NULL, stride, dist, FFTW_FORWARD, FFTW_PATIENT);
  dft_backward = fftw_plan_many_dft(rank, sites, howmany, out, NULL, stride, dist, in, NULL, stride, dist, FFTW_BACKWARD, FFTW_PATIENT);
  cout << "DONE." << endl;

  InitializeFields();

}

Lattice::~Lattice()
{
  delete[] gauge_field;

  for (int idx=0; idx<volume; ++idx) {
    delete propagator[idx];
  }
  delete[] propagator;

  fftw_free(in);
  fftw_free(out);
  fftw_destroy_plan(dft_forward); 
  fftw_destroy_plan(dft_backward); 
  fftw_cleanup();
}


void Lattice::InitializeFields()
{
  // It doesn't really matter what the initial field is. For debuging perposes, I chose
  // a Normal distribution. this gets wiped out once you update the lattice.
  for (int idx=0; idx<gauge_field_size; ++idx) { gauge_field[idx] = rng.Normal(0,1); }
}

void Lattice::RandomGauge()
{
  double* alpha = new double [volume];
  int s[4];
  int idx0, idx1;

  for (int idx=0; idx<volume; ++idx) { alpha[idx] = rng.Normal(0,1);  }

  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    idx0 = SiteIndex(s);

    for (int mu=0; mu<4; ++mu) {
      s[mu] += 1;
      idx1 = SiteIndex(s);
      s[mu] -= 1;
      GaugeField(mu,s) += alpha[idx1] - alpha[idx0];
    }

  }

  delete[] alpha;

}

void Lattice::CoulombGauge(int nu)
{

  int s[4];
  int idx;
  complex<double> p[4];
  complex<double> a, b;
  double p_sq;
  double eps = 1e-123;

  for (idx=0; idx<gauge_field_size; ++idx) {
    in[idx][0] = gauge_field[idx];
    in[idx][1] = 0.0;
  }

  fftw_execute(dft_forward);  

  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    idx = SiteIndex(s);
    a = 0.0;
    p_sq = 0.0;
    for (int mu=0; mu<4; ++mu) {
      p[mu] = propagator[idx]->Momentum(mu);
      if (mu!=nu) {
        a += conj( p[mu] ) * complex<double>(out[mu+4*idx][0], out[mu+4*idx][1] );
        p_sq += real( p[mu]*conj(p[mu]) );
      }
    }
    for (int mu=0; mu<4; ++mu) {
      b = a*p[mu];
      b /= p_sq+eps;
      b /= volume;
      out[mu+4*idx][0] = real(b); 
      out[mu+4*idx][1] = imag(b); 
    }
  }

  fftw_execute(dft_backward);  

  for (idx=0; idx<gauge_field_size; ++idx) {
    gauge_field[idx] -= in[idx][0];
  }


}




void Lattice::LorenzGauge() { CoulombGauge(99); }

double Lattice::MeanDivSq(int nu)
{
  double div;
  double mean_div_sq = 0.0;
  int s[4];
  int idx;

  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    idx = SiteIndex(s);
    div = 0.0;
    for (int mu=0; mu<4; ++mu) {
      if (mu!=nu) {
        div  += GaugeField(mu,s);
        s[mu] -= 1;
        div -= GaugeField(mu,s);
        s[mu] += 1;
      }
    }
    mean_div_sq += div*div;
  }
  return mean_div_sq / volume;
 
}



double Lattice::MeanDivSq()
{
  return MeanDivSq(99); 
}

double Lattice::MeanCosPlaquette(int mu, int nu, int m, int n)
{

  double cos_plaquette = 0.0;
  double plaq;
  int s[4]; 

  for ( s[0]=0; s[0]<sites[0]; ++s[0] ) 
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    plaq = 0.0;
    for (int i=0; i<m; ++i) { plaq += GaugeField(mu,s); s[mu]+=1; }
    for (int i=0; i<n; ++i) { plaq += GaugeField(nu,s); s[nu]+=1; }
    for (int i=0; i<m; ++i) { s[mu]-=1; plaq -= GaugeField(mu,s); }
    for (int i=0; i<n; ++i) { s[nu]-=1; plaq -= GaugeField(nu,s); }
    cos_plaquette += cos(plaq);
  }
  cos_plaquette /= volume;
  return cos_plaquette;
}

double Lattice::MeanCosPlaquette(int m, int n)
{

  double cos_plaquette = 0.0;

  for (int mu=0; mu<4; ++mu) 
  for (int nu=0; nu<4; ++nu) {
    cos_plaquette += 1.0 - MeanCosPlaquette(mu, nu, m, n);
  }
  cos_plaquette /= -12.0;
  cos_plaquette += 1.0;
  return cos_plaquette;
}



double Lattice::MeanPlaquetteSq(int mu, int nu, int m, int n)
{

  double plaquette_sq = 0.0;
  double plaq;
  int s[4]; 

  for ( s[0]=0; s[0]<sites[0]; ++s[0] ) 
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    plaq = 0.0;
    for (int i=0; i<m; ++i) { plaq += GaugeField(mu,s); s[mu]+=1; }
    for (int i=0; i<n; ++i) { plaq += GaugeField(nu,s); s[nu]+=1; }
    for (int i=0; i<m; ++i) { s[mu]-=1; plaq -= GaugeField(mu,s); }
    for (int i=0; i<n; ++i) { s[nu]-=1; plaq -= GaugeField(nu,s); }
    plaquette_sq += plaq*plaq;

  }
  plaquette_sq /= volume;
  return plaquette_sq ;
}

double Lattice::MeanPlaquetteSq(int m, int n)
{

  double plaquette_sq = 0.0;

  for (int mu=0; mu<4; ++mu) 
  for (int nu=0; nu<4; ++nu) {
    plaquette_sq += MeanPlaquetteSq(mu, nu, m, n);
  }
  plaquette_sq /= 12.0;
  return plaquette_sq ;
}

double Lattice::MeanLink()
{
  double mean_link = 0.0;
  for (int idx=0; idx<gauge_field_size; ++idx) { mean_link += gauge_field[idx]; }
  return mean_link / gauge_field_size; 
}

double Lattice::MeanLink(int mu)
{
  int s[4];
  double mean_link = 0.0;
  for ( s[0]=0; s[0]<sites[0]; ++s[0] ) 
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    mean_link += GaugeField(mu, s);
  }
  return mean_link / volume; 
}

double Lattice::MeanLinkSq()
{
  double mean_link_sq = 0.0;
  for (int idx=0; idx<gauge_field_size; ++idx) { mean_link_sq += gauge_field[idx]*gauge_field[idx]; }
  return mean_link_sq / gauge_field_size; 
}


void Lattice::MassiveUpdate(double beta, double mass2)
{
  MassiveUpdate(beta, 1.0, mass2);
}

void Lattice::MassiveUpdate(double beta, double xi, double mass2)
{

  int s[4];
  int idx;

  complex<double> a[4];
  complex<double> b[4];
  double p_sq, p3vec_sq, pXi_sq;

  double sigma[4];

  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {

    idx = SiteIndex(s);

    p_sq = propagator[idx]->FourMomentumSquared(); // |p|^2
    p3vec_sq = propagator[idx]->ThreeMomentumSquared(); // |\vec p|^2
    pXi_sq = propagator[idx]->FourMomentumSquared(xi);   // |p_\xi|^2 = xi |p_0|^2 + 1/xi |\vec p|^2

    sigma[0] = 1.0/sqrt( mass2 );
    sigma[1] = 1.0/sqrt( mass2 +  xi*p_sq  );
    for (int mu=2; mu<4; ++mu) {
      sigma[mu] = 1.0/sqrt( mass2 + pXi_sq  );
    }
    for (int mu=0; mu<4; ++mu) {
      a[mu] = rng.Normal( sigma[mu]/sqrt(beta) );
    }
    for (int mu=0; mu<4; ++mu) {
      b[mu] = 0.0;
      for (int nu=0; nu<4; ++nu) {
        b[mu] += propagator[idx]->SimilarityTransformationMatrix(mu,nu)*a[nu];
      }
      b[mu] /= sqrt(volume);
    }

    for (int mu=0; mu<4; ++mu) {
      out[mu+4*idx][0] = real( b[mu] );
      out[mu+4*idx][1] = imag( b[mu] );
    }

  }

  fftw_execute(dft_backward);

  for (idx=0; idx<gauge_field_size; ++idx) {
    gauge_field[idx] = in[idx][0];
  }

}


void Lattice::MassiveRxiGaugeUpdate(double beta, double xi, double mass2)
{

  if (mass2<=0) {
    cout << "mass2 must be greater than zero!" << endl;
    exit(EXIT_FAILURE);
  }


  int s[4];
  int idx;

  complex<double> a[4];
  complex<double> b[4];
  double p_sq;
  double eps = 1e-123;

  double sigma[4];

  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {

    idx = SiteIndex(s);

    p_sq = propagator[idx]->FourMomentumSquared(); // |p|^2

    sigma[0] = sqrt(2.0*xi+eps) / sqrt( (2*xi+eps)*mass2 + p_sq  );
    for (int mu=1; mu<4; ++mu) {
      sigma[mu] = 1.0/sqrt( mass2 + p_sq );
    }
    for (int mu=0; mu<4; ++mu) {
      a[mu] = rng.Normal( sigma[mu]/sqrt(beta) );
    }
    for (int mu=0; mu<4; ++mu) {
      b[mu] = 0.0;
      for (int nu=0; nu<4; ++nu) {
        b[mu] += propagator[idx]->SimilarityTransformationMatrix(mu,nu)*a[nu];
      }
      b[mu] /= sqrt(volume);
    }

    for (int mu=0; mu<4; ++mu) {
      out[mu+4*idx][0] = real( b[mu] );
      out[mu+4*idx][1] = imag( b[mu] );
    }

  }

  fftw_execute(dft_backward);

  for (idx=0; idx<gauge_field_size; ++idx) {
    gauge_field[idx] = in[idx][0];
  }

}

int Lattice::SiteIndex(int* s)
{
  int t[4];
  for (int nu=0; nu<4; ++nu) { t[nu] = ( s[nu]+sites[nu] ) % sites[nu]; } // Impose PBCs
  int idx = t[3]+sites[3]*( t[2] + sites[2]*( t[1]+ sites[1]*t[0] ) );
  return idx;
}

int Lattice::GaugeIndex(int mu, int* s) { return mu+4*SiteIndex(s); }

double& Lattice::GaugeField(int mu, int* s)
{
  int idx = GaugeIndex(mu, s);
  return gauge_field[idx];
}

double& Lattice::GaugeField(int idx) { return gauge_field[idx]; }

int Lattice::Sites(int mu)
{
  if (mu<0||mu>3) {
    cout << "Direction must be 0, 1, 2 or 3" << endl;
    exit(EXIT_FAILURE);
  }
  return sites[mu];
}

void Lattice::Import(string filename)
{

  ifstream file;
  file.open( filename.c_str(), ios::in | ios::binary );
  if (!file) {
    cerr << "Can't open file " << filename << endl;
    exit(EXIT_FAILURE);
  }
  file.read( (char*)gauge_field, gauge_field_size*sizeof(double) );
  file.close();
}


void Lattice::Export(string filename)
{

  ofstream file;
  file.open( filename.c_str(), ios::out | ios::binary );
  if (!file) {
    cerr << "Can't open file " << filename << endl;
    exit(EXIT_FAILURE);
  }
  file.write( (char*)gauge_field, gauge_field_size*sizeof(double) );
  file.close();
}

void Lattice::Export(string filename, int* ordering)
{

  int check[4];
  for (int mu=0; mu<4; ++mu) { check[mu] = 0; }

  for (int mu=0; mu<4; ++mu)
  for (int nu=0; nu<4; ++nu) {
    if ( ordering[mu] == nu ) {
      check[nu] += 1;
    }
  }

  for (int mu=0; mu<4; ++mu) {
    if (check[mu] != 1 ) {
      cerr << "Ordering is out of bounds or not one-to-one." << endl;
      exit(EXIT_FAILURE);
    }
  }

  ofstream file;
  file.open( filename.c_str(), ios::out);
  if (!file) {
    cerr << "Can't open file " << filename << endl;
    exit(EXIT_FAILURE);
  }
  file << setprecision(15);
  int s[4];

  int i = ordering[0];
  int j = ordering[1];
  int k = ordering[2];
  int l = ordering[3];

  for ( s[i]=0; s[i]<sites[i]; ++s[i] )
  for ( s[j]=0; s[j]<sites[j]; ++s[j] )
  for ( s[k]=0; s[k]<sites[k]; ++s[k] )
  for ( s[l]=0; s[l]<sites[l]; ++s[l] )
  for ( int mu=0; mu<4; ++mu ) {
    file << s[0] << " ";
    file << s[1] << " ";
    file << s[2] << " ";
    file << s[3] << " ";
    file << ordering[mu] << " ";
    file << GaugeField( ordering[mu], s);
    file << endl;
  }

  file.close();

}


complex<double> Lattice::PolyakovLoop(int mu, int* s)
{

  complex<double> I(0.0,1.0);
  double loop;
  int idx;
  int t[4];
  for (int nu=0; nu<4; ++nu) { t[nu] = s[nu]; }

  loop = 0.0;
  for ( t[mu]=0; t[mu]<sites[mu]; ++t[mu] ) {
    idx = SiteIndex(t);
    loop += GaugeField(mu,t);
  }

  return exp(I*loop);

}

complex<double> Lattice::MeanPolyakovLoop(int mu)
{

  complex<double> mean(0.0);
  int s[4];

  s[0] = 0;
  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    mean += PolyakovLoop(mu,s);
  }
  mean /= volume;
  return mean;

}


complex<double> Lattice::MeanPolyakovLoop(int mu, int nu, int* r)
{

  complex<double> mean(0.0);
  int s[4];
  int t[4];

  s[0] = 0;
  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    for (int rho=0; rho<4; ++rho) { t[rho] = ( s[rho]+r[rho]+sites[rho] ) % sites[rho]; } 
    mean += PolyakovLoop(mu,s) * conj( PolyakovLoop(nu,t) );
  }
  mean /= volume;
  return mean;
  
}

void Lattice::CoulombGaugeUpdate(double beta)
{
  CoulombGaugeUpdate(beta, 1.0);
}


void Lattice::CoulombGaugeUpdate(double beta, double xi)
{

  int s[4];
  int idx;

  double p3vec_sq, pXi_sq, norm_sq;
  complex<double> a[4];
  complex<double> b[4];
  complex<double> p[4];

  for ( s[0]=0; s[0]<sites[0]; ++s[0] )
  for ( s[1]=0; s[1]<sites[1]; ++s[1] )
  for ( s[2]=0; s[2]<sites[2]; ++s[2] )
  for ( s[3]=0; s[3]<sites[3]; ++s[3] ) {
    idx = SiteIndex(s);
    for (int mu=0; mu<4; ++mu) {
      p[mu] = propagator[idx]->Momentum(mu);
      a[mu] = 0.0;
    }

    p3vec_sq = propagator[idx]->ThreeMomentumSquared(); // |\vec p|^2
    pXi_sq = propagator[idx]->FourMomentumSquared(xi);   // |p_\xi|^2 = xi |p_0|^2 + 1/xi |\vec p|^2

    // case : mu==0 
//    if ( s[0]!=0 ) {
      if ( (s[1]!=0)||(s[2]!=0)||(s[3]!=0) ) {
        a[0] = rng.Normal( 1.0/sqrt(beta*xi*p3vec_sq) );
     }
//    }

    // case : mu!=0 
    if ( s[3]!=0 ) {
      norm_sq = real( p[1]*conj(p[1]) + p[2]*conj(p[2]) );
      if ( (s[1]==0) && (s[2]==0) ) {
        a[1] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
        a[2] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
      } else {
        b[1] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
        b[2] = rng.Normal( 1.0/sqrt( beta*pXi_sq*( 1.0+norm_sq/real(p[3]*conj(p[3]) )) ) );
        if ( (s[1]!=0) && (s[2]!=0) ) {
          a[1] = p[1]*conj(p[2])*( -b[1]/abs(p[1]) + b[2]/abs(p[2]));
          a[1] /= sqrt(norm_sq);
          a[2] = abs(p[1])*b[1] + abs(p[2])*b[2];
          a[2] /= sqrt(norm_sq);
        }
        if ( (s[1]!=0) && (s[2]==0) ) {
          a[1] = p[1]*b[2]/abs(p[1]);
          a[2] = b[1];
        }
        if ( (s[1]==0) && (s[2]!=0) ) {
          a[1] = -conj(p[2])*b[1]/abs(p[2]);
          a[2] = b[2];
        }
        a[3] = -( conj(p[1])*a[1]+conj(p[2])*a[2] );
        a[3] /= conj(p[3]);
      }
    } else {
      if ( (s[1]!=0)&&(s[2]!=0) ) {
        a[1] = rng.Normal( 1.0/sqrt(beta*pXi_sq*(1.0+real(p[1]*conj(p[1]))/real(p[2]*conj(p[2])))) );
        a[2] = -conj(p[1])*a[1]/conj(p[2]);
        a[3] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
      }
      if ( (s[1]!=0)&&(s[2]==0) ) {
        a[2] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
        a[3] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
      }
      if ( (s[1]==0)&&(s[2]!=0) ) {
        a[1] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
        a[3] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
      }
      if ( (s[1]==0)&&(s[2]==0)&&(s[0]!=0) ) {
        a[1] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
        a[2] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
        a[3] = rng.Normal( 1.0/sqrt(beta*pXi_sq) );
      }
    }
    for (int mu=0; mu<4; ++mu) {
      a[mu] /= sqrt(volume);
      out[mu+4*idx][0] = real(a[mu]);
      out[mu+4*idx][1] = imag(a[mu]);
    }
  }

  fftw_execute(dft_backward);

  for (idx=0; idx<gauge_field_size; ++idx) {
    gauge_field[idx] = in[idx][0];
  }

}



void Lattice::PrintObservables(int verbose)
{

  if (verbose) {
     cout << endl;
     cout << "O = Gauge invariant quantity" << endl;
     cout << "X = Gauge variant quantity" << endl;
  }

  if (verbose) cout << endl << "(O) Link <A> : ";
  cout << MeanLink() << " ";
  if (verbose) cout << endl << "(X) Link Squared <AA> : ";
  cout << MeanLinkSq() << " ";
  if (verbose) cout << endl << "(X) Four-divergence Squared <(d.A)^2> : ";
  cout << MeanDivSq() << " ";
  if (verbose) cout << endl << "(X) Three-divergence Squared (0) <(d.A)^2> : ";
  cout << MeanDivSq(0) << " ";
  if (verbose) cout << endl << "(X) Three-divergence Squared (1) <(d.A)^2> : ";
  cout << MeanDivSq(1) << " ";
  if (verbose) cout << endl << "(X) Three-divergence Squared (2) <(d.A)^2> : ";
  cout << MeanDivSq(2) << " ";
  if (verbose) cout << endl << "(X) Three-divergence Squared (3) <(d.A)^2> : ";
  cout << MeanDivSq(3) << " ";

  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <FF> : ";
  cout << MeanPlaquetteSq(1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <F12F12> : ";
  cout << MeanPlaquetteSq(1,2,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <F13F13> : ";
  cout << MeanPlaquetteSq(1,3,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <F23F23> : ";
  cout << MeanPlaquetteSq(2,3,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <F01F01> : ";
  cout << MeanPlaquetteSq(0,1,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <F02F02> : ";
  cout << MeanPlaquetteSq(0,2,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <F03F03> : ";
  cout << MeanPlaquetteSq(0,3,1,1) << " ";

  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <cos(F)> : ";
  cout << MeanCosPlaquette(1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <cos(F12)> : ";
  cout << MeanCosPlaquette(1,2,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <cos(F13)> : ";
  cout << MeanCosPlaquette(1,3,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <cos(F23)> : ";
  cout << MeanCosPlaquette(2,3,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <cos(F01)> : ";
  cout << MeanCosPlaquette(0,1,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <cos(F02)> : ";
  cout << MeanCosPlaquette(0,2,1,1) << " ";
  if (verbose) cout << endl << "(O) 1x1 Plaquette Squared <cos(F03)> : ";
  cout << MeanCosPlaquette(0,3,1,1) << " ";

  complex<double> polyakov;
  polyakov = MeanPolyakovLoop(0);
  if (verbose) cout << endl << "(O) Polyakov loop re( <P_0> ) : ";
  cout << real(polyakov) << " ";
  if (verbose) cout << endl << "(O) Polyakov loop im( <P_0> ): ";
  cout << imag(polyakov) << " ";
  polyakov = MeanPolyakovLoop(1);
  if (verbose) cout << endl << "(O) Polyakov loop re( <P_1> ) : ";
  cout << real(polyakov) << " ";
  if (verbose) cout << endl << "(O) Polyakov loop im( <P_1> ): ";
  cout << imag(polyakov) << " ";
  polyakov = MeanPolyakovLoop(2);
  if (verbose) cout << endl << "(O) Polyakov loop re( <P_2> ) : ";
  cout << real(polyakov) << " ";
  if (verbose) cout << endl << "(O) Polyakov loop im( <P_2> ): ";
  cout << imag(polyakov) << " ";
  polyakov = MeanPolyakovLoop(3);
  if (verbose) cout << endl << "(O) Polyakov loop re( <P_3> ) : ";
  cout << real(polyakov) << " ";
  if (verbose) cout << endl << "(O) Polyakov loop im( <P_3> ): ";
  cout << imag(polyakov) << " ";

  if (verbose) cout << endl << "(O) Link <A_0> : ";
  cout << MeanLink(0) << " ";
  if (verbose) cout << endl << "(O) Link <A_1> : ";
  cout << MeanLink(1) << " ";
  if (verbose) cout << endl << "(O) Link <A_2> : ";
  cout << MeanLink(2) << " ";
  if (verbose) cout << endl << "(O) Link <A_3> : ";
  cout << MeanLink(3) << " ";

//  int r[4] = {0,1,2,3};
//  for (int mu=0;mu<4;++mu)
//  for (int nu=0;nu<4;++nu) {
//    cout << "(O) Polyakov loop 2-pt; (mu,nu)=(" << mu << "," <<nu << ") and  r=("<< r[0] << "," << r[1] << "," << r[2] << ","<<r[3] << ") : " << MeanPolyakovLoop(mu,nu,r) << endl;
//  }

  cout << endl;
}
