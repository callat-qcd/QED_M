#include <stdlib.h>
#include <complex>
#include "fftw3.h"
#include "propagator.h"

#ifndef INCLUDED_LATTICE
#define INCLUDED_LATTICE

// This class contains a memory allocation for the noncompact U(1) gauge field.
// It also contains class functions which operate on the field, such as updating
// the field, and performing gauge transformations on the field.
// In addition to these, the class contains functions ofr performing measurements of
// observables, and configuration I/O
class Lattice
{
  private:
    int sites[4]; // Dimensions of the lattice ordered as (T,X,Y,Z)

    int volume; // Product of the dimensions V = T*X*Y*Z

    int gauge_field_size; // This is (4 components * V sites)
    double* gauge_field; // Pointer to the gauge field configuration

    Propagator** propagator; // Pointer to pointers to similarity transformations used to diagonalize the
                             // gauge field propagator in momentum space, there is an instance of this
                             // class created for each momentum vector (pt,px,py,pz)

    fftw_complex* in; // These are the in and out vectors for used to perform FFTs
    fftw_complex* out;
    fftw_plan dft_forward; // The FFT is most effient when a "plan" is created, specific to the problem
    fftw_plan dft_backward; // it is created once when the Lattice class is intantiated, and then used 
                            // repeatedly during FFTs afterwards

    Lattice& operator=(const Lattice&); // This prevents assignment of Lattice objects 
                                        // since doing the assignments properly requires care
                                        // when memory is allocated within a constructor; I'm lazy
    Lattice(const Lattice&); // Prevents copy construction

  public:
    explicit Lattice(int*); // Instantiation of a lattice object requires the dimensions (T,X,Y,Z)
                            // as an argument
    ~Lattice();

    void InitializeFields();

    void RandomGauge(); // Performs a random gauge transformation on the field
    void LorenzGauge(); // Puts the field in Lorenz Gauge
    void CoulombGauge(int); // Puts the field in a Coulomb-like gauge, where the preferred direction is
                            // specified by an integer T=0, X=1, Y=2, Z=3

    double MeanDivSq(); // Space-time averaged mean four-divergence squared
    double MeanDivSq(int); // Space-time averaged mean three-divergence squared, with preferred direction
    double MeanLink(); // Space-time average of 1/4 \sum_\mu A_\mu
    double MeanLink(int); // Space-time average of A_\mu
    double MeanLinkSq(); // Space-time average of 1/4 \sum_\mu A_\mu^2

    double MeanPlaquetteSq(int,int,int,int); // Space-time averaged P_{\mu,\nu}^2 with n (m) steps in direction mu (nu)
    double MeanPlaquetteSq(int,int); // Space-time averaged 1/6 \sum_{\mu>\nu} P_{\mu,\nu}^2 with n (m) steps in direction mu (nu)

    double MeanCosPlaquette(int,int,int,int); // Space-time averaged cos( P_{\mu,\nu} ) with n (m) steps in direction mu (nu)
    double MeanCosPlaquette(int,int); // Space-time averaged 1/6 \sum_{\mu\nu} cos( P_{\mu,\nu} ) with n (m) steps in direction mu (nu)

    complex<double> PolyakovLoop(int,int*); // Polyakov loop at a specified site
    complex<double> MeanPolyakovLoop(int); // Space-time averaged polyakov loop
    complex<double> MeanPolyakovLoop(int,int,int*); // Space-time averaged polyakov loops in the mu and nu directions, separated by a distance

    void PrintObservables(int);

    void MassiveUpdate(double, double); // beta, mass, xi=1
    void MassiveUpdate(double, double, double); // beta, xi, mass, where xi is an anisotropy factor
    void CoulombGaugeUpdate(double); // beta, xi=1
    void CoulombGaugeUpdate(double, double); // beta, xi
    void MassiveRxiGaugeUpdate(double, double, double); // beta, xi, mass !! Here, xi is the R-xi xi!! Not anisotropy.

    int Sites(int); // Returns the number of sites in a particular direction
    int SiteIndex(int*); // Gives the collective index associated with the site (t,x,y,z)
    int GaugeIndex(int, int*); // Gives the collective index associated with direction mu and site (t,x,y,z)

    double& GaugeField(int, int*); // Returns the gauge field given mu and (t,x,y,z)
    double& GaugeField(int); // Returns the gauge field given the Gauge index

    void Export(string); // I/O to a binary file specifed by the string
    void Import(string); // Import from a binary file specifed by the string
    void Export(string, int*); // Write a human readable file, the last argument is used to specify
                               // ordering of the output

};
#endif
