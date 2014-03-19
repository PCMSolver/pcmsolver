#ifndef CAVITYDATA_HPP
#define CAVITYDATA_HPP

#include <string>
#include <vector>

#include "Config.hpp"

#include "Sphere.hpp"
#include "Symmetry.hpp"

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
 *  @var cavityData::minimalRadius
 *  Triggers the addition of spheres not centered on atoms.
 *  Relevant for GePolCavity.
 *  @var cavityData::patchLevel
 *  Relevant for WaveletCavity.
 *  @var cavityData::coarsity
 *  Relevant for WaveletCavity.
 *  @var cavityData::filename
 *  Name of the file containing the cavity
 *  specification for a restart.
 *  @var cavityData::symmetry
 *  Integer specifying the point group.
 */

struct cavityData
{
	std::vector<Sphere> spheres;
	double area;
	double probeRadius;
	double minDistance;
	int derOrder;
	double minimalRadius;
	int patchLevel;
	double coarsity;
	std::string filename;
	Symmetry symmetry;
	cavityData(const std::vector<Sphere> & _spheres, double _area, double _probeRadius, 
		   double _minDistance, int _derOrder, double _minRadius, 
		   int _patchLevel, double _coarsity, const std::string & _fname,
		   const Symmetry & _symmetry) :
	spheres(_spheres), area(_area), probeRadius(_probeRadius), 
	minDistance(_minDistance), derOrder(_derOrder), minimalRadius(_minRadius), 
	patchLevel(_patchLevel), coarsity(_coarsity), filename(_fname), symmetry(_symmetry) {}
};

#endif // CAVITYDATA_HPP
