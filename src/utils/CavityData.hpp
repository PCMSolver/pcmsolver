#ifndef CAVITYDATA_HPP
#define CAVITYDATA_HPP

#include <vector>

#include "Config.hpp"

#include "Sphere.hpp"

/*! @struct cavityData
 *  @brief Contains all data defined from user input in the cavity section.
 *  @var cavityData::spheres 
 *  Contains the list of generating spheres.
 *  @var cavityData::area
 *  The average tesserae area. Relevant for GePolCavity.
 *  @var cavityData::probeRadius
 *  The radius of the spherical probe representing the solvent.
 *  @var cavityData::minDistance
 *  The minimal distance between two sampling 
 *  points on different spheres. Relevant for TsLessCavity.
 *  @var cavityData::derOrder
 *  The maximum derivative order to be used in the definition
 *  of the smoothing function. Relevant for TsLessCavity.
 *  @var cavityData::addSpheres
 *  Triggers the addition of spheres not centered on atoms.
 *  Relevant for GePolCavity.
 *  @var cavityData::patchLevel
 *  Relevant for WaveletCavity.
 *  @var cavityData::coarsity
 *  Relevant for WaveletCavity.
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
