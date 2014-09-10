#include <iostream>
#include "taylor.hpp"

// This is a very small example of taking derivatives
// using libtaylor. 

using namespace std;

// Define a function we want to differentiate: f(x,y)
template<class T>
T f(const T &x, const T &y)
{
  return sin(log(7*x)+exp(y))+9;
}

int main(void)
{
  // Compute a directional derivative of f(x,y) in the direction
  // (1,2), at point (x,y) = (3,4). This is equivalent to computing the
  // taylor expansion of f(3+1*eps, 4+2*eps) in the variable eps.
  const int Ndeg = 5; // Order of expansion.
  const int Nvar = 1; // Only one differentiation variable in this example
  taylor<double,Nvar,Ndeg> eps(0,0); // Set up seed variable.
  taylor<double,Nvar,Ndeg> fexpansion = f(3+eps,4+2*eps);
  // Now fexpansion contains Taylor coefficients. If we want
  // derivatives we have to multiply by the appropriate factorials
  fexpansion.deriv_facs();
  cout << "Directional derivative: " << fexpansion << endl;
  return 0;
}
