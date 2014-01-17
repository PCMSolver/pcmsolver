#ifndef GREENDATA_HPP
#define GREENDATA_HPP

#include <vector>

#include "Config.hpp"

/*! @struct greenData
 *  @brief Contains all data defined from user input in the green section.
 *  @var greenData::how
 *  The way the derivatives of the Green's function are evaluated.
 *  @var greenData::epsilon
 *  The permittivity.
 *  @var greenData::kappa
 *  Inverse of the Debye length.
 *  @var greenData::epsilonReal
 *  Real part of the permittivity of a metal sphere.
 *  @var greenData::epsilonImaginary
 *  Imaginary part of the permittivity of a metal sphere.
 *  @var greenData::spherePosition
 *  Coordinates of the metal sphere center.
 *  @var greenData::sphereRadius
 *  Radius of the the metal sphere.
 */

struct greenData
{
	int how;
	double epsilon;
	double kappa;
	double epsilonReal;
	double epsilonImaginary;
	std::vector<double> spherePosition;
	double sphereRadius;

	greenData(int _how, double _epsilon = 1.0, double _kappa = 0.0, double _epsReal = 0.0, double _epsImaginary = 0.0, 
		const std::vector<double> & _sphere = std::vector<double>(), double _sphRadius = 0.0) :
		how(_how), epsilon(_epsilon), kappa(_kappa), epsilonReal(_epsReal), epsilonImaginary(_epsImaginary), 
		spherePosition(_sphere), sphereRadius(_sphRadius) {}
};

#endif // GREENDATA_HPP
