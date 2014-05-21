#ifndef RANDWERTE_HPP
#define RANDWERTE_HPP

class Vector3;

typedef struct {
	Vector3 *Chi;     ///< the value of the interpolation
	Vector3 *n_Chi;   ///< the normal derivative of the interpolation
	double *det_dChi; ///< ???
  unsigned int noP; ///< number of points
} Randwerte;

#endif
