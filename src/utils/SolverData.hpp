#ifndef SOLVERDATA_HPP
#define SOLVERDATA_HPP

#include <vector>

#include "Config.hpp"

class GreensFunction;

/*! @struct solverData
 *  @brief Contains all data defined from used input in the cavity section.
 *  @var solverData::spheres 
 *  Contains the list of generating spheres.
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
