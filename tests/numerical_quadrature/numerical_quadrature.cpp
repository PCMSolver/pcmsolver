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

#include "cnpyPimpl.hpp"
#include "Element.hpp"
#include "GePolCavity.hpp"
#include "MathUtils.hpp"
#include "PhysicalConstants.hpp"
#include "Sphere.hpp"
#include "Symmetry.hpp"

typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
singleLayerIntegrand;
typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &, 
		        const Eigen::Vector3d &)> doubleLayerIntegrand;

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_sphere tests numerical quadrature on a sphere for integrand = 1.0
 */
BOOST_AUTO_TEST_CASE(sphere)
{
    struct f {
	double operator()(const Eigen::Vector3d & s, const Eigen::Vector3d & p) { return 1.0; }
    };

    Eigen::Vector3d origin(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(origin, 1.55);
    spheres.push_back(sph1);
    double area = 0.4;
    // C1
    Symmetry pGroup = buildGroup(0, 0, 0, 0);
    GePolCavity cavity(spheres, area, 0.0, 100.0, pGroup);
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());
    
    singleLayerIntegrand F = f();
    
    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator<32, 16>(F, cavity.elements(i));
    }

    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), cavity.elementArea(i), 1.0e-12);
    }
}

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_sphere_1r tests numerical quadrature on a sphere for integrand 1.0/r
 */
BOOST_AUTO_TEST_CASE(sphere_1r)
{
    Eigen::Vector3d origin(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(origin, 1.55);
    spheres.push_back(sph1);
    double area = 0.4;
    // C1
    Symmetry pGroup = buildGroup(0, 0, 0, 0);
    GePolCavity cavity(spheres, area, 0.0, 100.0, pGroup);

    struct f {
	double operator()(const Eigen::Vector3d & s, const Eigen::Vector3d & p) { double radius = 1.55; return (1.0/radius); }
    };

    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());
    
    singleLayerIntegrand F = f();
    
    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator<32, 16>(F, cavity.elements(i));
    }

    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), (cavity.elementArea(i)/sph1.radius()), 1.0e-12);
    }
}

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_molecule tests numerical quadrature function on H2 cavity
 */
BOOST_AUTO_TEST_CASE(molecule)
{
    struct f {
	double operator()(const Eigen::Vector3d & s, const Eigen::Vector3d & p) { return 1.0; }
    };

    Eigen::Vector3d H1( 0.735000, 0.000000, 0.000000);
    Eigen::Vector3d H2(-0.735000, 0.000000, 0.000000);
    std::vector<Sphere> spheres;
    double radiusH = 1.20;
    Sphere sph1(H1, radiusH);
    Sphere sph2(H2, radiusH);
    spheres.push_back(sph1);
    spheres.push_back(sph2);
    double area = 0.2;
    double probeRadius = 1.385;
    double minRadius = 0.2;
    // C1
    Symmetry pGroup = buildGroup(0, 0, 0, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, pGroup);
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());
    
    singleLayerIntegrand F = f();
    
    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator<64, 16>(F, cavity.elements(i));
	double diff = results(i) - cavity.elementArea(i); 
	if (diff > 1.0e-12) {
	    BOOST_TEST_MESSAGE("Test versus area for H2 molecule");
            BOOST_TEST_MESSAGE("Tessera n. " << i+1);
	    BOOST_TEST_MESSAGE("diff = " << results(i) - cavity.elementArea(i)); 
	}
    }
    /* 
    // In case you need to update the reference files...
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("molecule.npy", results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_ref = cnpy::npy_load("molecule.npy");
    int dim = raw_ref.shape[0];
    Eigen::VectorXd reference = Eigen::VectorXd::Zero(dim);
    reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_ref.data);
    
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), reference(i), 1.0e-12);
    }
}

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_molecule_1r tests numerical quadrature function on H2 cavity
 */
BOOST_AUTO_TEST_CASE(molecule_1r)
{
    struct f {
	double operator()(const Eigen::Vector3d & s, const Eigen::Vector3d & p) { double radius = 1.55; return (1.0/radius); }
    };

    Eigen::Vector3d H1( 0.735000, 0.000000, 0.000000);
    Eigen::Vector3d H2(-0.735000, 0.000000, 0.000000);
    std::vector<Sphere> spheres;
    double radiusH = 1.20;
    Sphere sph1(H1, radiusH);
    Sphere sph2(H2, radiusH);
    spheres.push_back(sph1);
    spheres.push_back(sph2);
    double area = 0.2;
    double probeRadius = 1.385;
    double minRadius = 0.2;
    // C1
    Symmetry pGroup = buildGroup(0, 0, 0, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, pGroup);
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());
    
    singleLayerIntegrand F = f();
    
    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator<64, 16>(F, cavity.elements(i));
	double diff = results(i) - (cavity.elementArea(i)/sph1.radius()); 
	if (diff > 1.0e-12) {
	    BOOST_TEST_MESSAGE("Test versus area divided by radius for H2 molecule");
            BOOST_TEST_MESSAGE("Tessera n. " << i+1);
	    BOOST_TEST_MESSAGE("diff = " << results(i) - cavity.elementArea(i)); 
	}
    }
    /*
    // In case you need to update the reference files...
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("molecule_1r.npy", results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_ref = cnpy::npy_load("molecule_1r.npy");
    int dim = raw_ref.shape[0];
    Eigen::VectorXd reference = Eigen::VectorXd::Zero(dim);
    reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_ref.data);
    
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), reference(i), 1.0e-12);
    }
}
