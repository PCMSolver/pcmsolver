#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#ifndef pi
  #define pi 3.1415926535897932385
#endif

// Constants regarding solvent properties
const double epsilon = 78.39;   ///< dielectric constant of the solvent
const double kappaIS = 0.0;       ///< ion screening

extern unsigned int minLevel;

const unsigned int delta = 10; ///< constant for memory allocation

const double eps = 1e-6; ///< accuracy for point equality

const double		op = -1; ///< order of the operator

const double		a = 1.25; ///< compression constant,  a > 1
const double		b = 0.001; ///< compression constant, 0 < b < 1

// quadrature
const double		scalingFactor = 0.7071;///< size of relative outer radius
const unsigned int minQuadratureLevel =2; ///< minimal quadrature level
const int g_max = 15; ///< the maximum number of Cubatures calculated

const double q = 1; ///<

#endif
