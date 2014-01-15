#ifndef CAVITYDATA_HPP
#define CAVITYDATA_HPP

#include <vector>

#include "Config.hpp"

#include "Sphere.hpp"

/*! @struct cavityData
 *  @brief Contains all data defined from used input in the cavity section.
 *  @var cavityData::spheres 
 *  Contains the list of generating spheres.
 */

struct cavityData
{
	std::vector<Sphere> spheres;
	double area;
	double probeRadius;
	double minDistance;
	int derOrder;
	bool addSpheres;
	int patchLevel;
	double coarsity;
	cavityData(const std::vector<Sphere> & _spheres, double _area, double _probeRadius = 0.0, 
		   double _minDistance = 0.1, int _derOrder = 4, bool _addSpheres = false, 
		   int _patchLevel = 2, double _coarsity = 0.5) :
	spheres(_spheres), area(_area), probeRadius(_probeRadius), 
	minDistance(_minDistance), derOrder(_derOrder), addSpheres(_addSpheres), 
	patchLevel(_patchLevel), coarsity(_coarsity) {}
};

#endif // CAVITYDATA_HPP
