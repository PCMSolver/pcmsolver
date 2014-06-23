#define BOOST_TEST_MODULE NumericalQuadrature

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>
#include <boost/function.hpp>

#include "GePolCavity.hpp"
#include "MathUtils.hpp"
#include "PhysicalConstants.hpp"
#include "Sphere.hpp"
#include "Symmetry.hpp"

typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> integrand;

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_sphere tests numerical quadrature function on a sphere of radius 1.0
 */
BOOST_AUTO_TEST_CASE(sphere)
{
    struct f {
	double operator()(const Eigen::Vector3d & s, const Eigen::Vector3d & p) { return 1.0; }
    };

    Eigen::Vector3d origin(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(origin,  1.0);
    spheres.push_back(sph1);
    double area = 10.0;
    // C1
    Symmetry pGroup = buildGroup(0, 0, 0, 0);
    GePolCavity cavity(spheres, area, 0.0, 100.0, pGroup);
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());
    
    integrand F = f();
    
    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator(F, cavity.elements(i));
    }

    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), cavity.elementArea(i), 1.0e-12);
    }
}

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_molecule tests numerical quadrature function on H3+ cavity
 */
BOOST_AUTO_TEST_CASE(molecule)
{
    struct f {
	double operator()(const Eigen::Vector3d & s, const Eigen::Vector3d & p) { return 1.0; }
    };

    Eigen::Vector3d H1( 0.735000, 0.000000, -1.333333);
    Eigen::Vector3d H2(-0.735000, 0.000000, -1.333333);
    Eigen::Vector3d H3( 0.000000, 0.000000,  2.666667);
    std::vector<Sphere> spheres;
    double radiusH = (1.20 * 1.20) / convertBohrToAngstrom;
    Sphere sph2(H1, radiusH);
    Sphere sph3(H2, radiusH);
    Sphere sph4(H3, radiusH);
    spheres.push_back(sph2);
    spheres.push_back(sph3);
    spheres.push_back(sph4);
    double area = 0.2 / convertBohr2ToAngstrom2;
    double probeRadius = 1.385 / convertBohrToAngstrom;
    double minRadius = 0.2 / convertBohrToAngstrom;
    // C1
    Symmetry pGroup = buildGroup(0, 0, 0, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, pGroup);
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());
    
    integrand F = f();
    
    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator(F, cavity.elements(i));
    }

    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), cavity.elementArea(i), 1.0e-12);
    }
}
