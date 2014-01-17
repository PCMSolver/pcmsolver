#ifndef SOLVERDATA_HPP
#define SOLVERDATA_HPP

#include "Config.hpp"

class GreensFunction;

/*! @struct solverData
 *  @brief Contains all data defined from user input in the solver section.
 *  @var solverData::gfInside 
 *  The Green's function inside the cavity.
 *  @var solverData::gfOutside
 *  The Green's function outside the cavity.
 *  @var solverData::correction
 *  The correction factor to be use in a CPCM calculation.
 *  @var solverData::integralEquation
 *  The type of integral equation to solve, relevant only for wavelet solvers.
 */

struct solverData
{
	GreensFunction * gfInside;
	GreensFunction * gfOutside;
	double correction;
	int integralEquation;
	solverData(GreensFunction * _gfInside, GreensFunction * _gfOutside, double _correction = 0.0,  int _integralEquation = 1) :
		gfInside(_gfInside), gfOutside(_gfOutside), correction(_correction), integralEquation(_integralEquation) {}
};

#endif // SOLVERDATA_HPP
