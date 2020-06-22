#include <stdlib.h>
#include <complex>

#ifndef INCLUDED_PROPAGATOR
#define INCLUDED_PROPAGATOR

// To understand what this class does, I'd advise reading my note about diagonalization
// of the photon propagator in momentum space in the anisotropic case
// This class generates and holds the appropriate similarity transformaion, given a momentum vector
// and dimensions of the lattice. An intance must be created for every momentum vector 
class Propagator
{
  private:
    
    complex<double> p[4];     //  p_\mu = 2 e^{i p/2} sin(p/2)
    complex<double> s[4][4]; // similarity transformation matrix that diagonalizes the propagator

    double p_sq; // four-momentum squared
    double p3vec_sq; // three-momentum squared

    double PI;
    complex<double> I;

    Propagator& operator=(const Propagator&); // Prevents assignment
    Propagator(const Propagator&); // Prevents copy construction

  public:
    explicit Propagator(int*, int*); // arguments are momentum p = (pt,px,py,pz) and sites = (T,X,Y,Z)
    ~Propagator();


    void PrintMomentum(); // for debugging purposes
    void PrintSimilarityTransformationMatrix(); // for debugging purposes

    double ThreeMomentumSquared(); // three-momentum squared 
    double FourMomentumSquared(); // four-momentum squared
    double FourMomentumSquared(double); // argument is the anisotropy factor xi
                                        // returns p^2 = xi p0^2 + 1/xi (p1^2+p2^2+p3^2)

    complex<double> Momentum(int); // returns a complex momentum associated with a site index
    complex<double>& SimilarityTransformationMatrix(int, int); // Gives the matrix elements of the Matrix


};
#endif
