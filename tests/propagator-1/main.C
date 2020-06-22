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

  int n[4] = { 1, 4, 4, 2 };
  int sites[4] = { 8, 8, 8, 8 };

  Propagator propagator(n,sites);

  cout << "Momentum Vector :" << endl;
  propagator.PrintMomentum();

  cout << "Momentum Vector (obtained from class function) :" << endl;
  for (int mu=0; mu<4; ++mu) {
    cout << mu << " : " << propagator.Momentum(mu) << endl;
  }

  cout << "Similarity Transformation Matrix :" << endl;
  propagator.PrintSimilarityTransformationMatrix();

  cout << "Momentum Squared :" << endl;
  cout << "p_sq : " << propagator.FourMomentumSquared() << endl;
  cout << "pXi_sq (xi=1.1) : " << propagator.FourMomentumSquared(1.1) << endl;
  cout << "p3vec_sq : " << propagator.ThreeMomentumSquared() << endl;

  exit(EXIT_SUCCESS);

}
