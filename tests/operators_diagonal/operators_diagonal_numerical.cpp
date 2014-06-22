#define BOOST_TEST_MODULE NumericalIntegrator 

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "AnisotropicLiquid.hpp"
#include "Element.hpp"
#include "GePolCavity.hpp"
#include "NumericalIntegrator.hpp"
#include "PhysicalConstants.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

struct NumericalIntegratorTest {
    double radius;
    Eigen::Vector3d epsilon;
    Eigen::Vector3d euler;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    GePolCavity cavity;
    NumericalIntegrator * diag;
    NumericalIntegratorTest() { SetUp(); }
    void SetUp() {
/*    	epsilon << 2.0, 80.0, 15.0;
    	euler << 6.0, 40.0, 15.0;*/
	epsilon << 80.0, 80.0, 80.0;
	euler << 0.0, 0.0, 0.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        
	Eigen::Vector3d origin(0.0, 0.0, 0.0);
        std::vector<Sphere> spheres;
	radius = 1.0;//1.44 / convertBohrToAngstrom;
        Sphere sph1(origin,  radius);
        spheres.push_back(sph1);
        double area = 10.0;
        // C1
        Symmetry pGroup = buildGroup(0, 0, 0, 0);
        cavity = GePolCavity(spheres, area, 0.0, 100.0, pGroup);
	
	diag = new NumericalIntegrator();
    }
};

/*! \class NumericalIntegrator 
 *  \test \b NumericalIntegratorTest_vacuum tests the numerical evaluation of the vacuum diagonal elements of S and D
 */
/*BOOST_FIXTURE_TEST_CASE(vacuum, NumericalIntegratorTest)
{
    Vacuum<double> gf(diag);

    double S, S_sum = 0.0;
    for (int i = 0; i < cavity.size(); ++i) {
    	S = gf.diagonalS(cavity.elements(i)); 
        S_sum += S;
    }
    std::cout << "S_sum = " << S_sum << std::endl;
    double ref_S = 4*M_PI*radius; 
    std::cout << "ref_S " << ref_S << std::endl;
    BOOST_REQUIRE_CLOSE(ref_S, S, 1.0e-12);
}*/

/*! \class NumericalIntegrator 
 *  \test \b NumericalIntegratorTest_uniformdielectric tests the numerical evaluation of the uniform dielectric diagonal elements of S and D
 */
/*BOOST_FIXTURE_TEST_CASE(uniformdielectric, NumericalIntegratorTest)
{
    UniformDielectric<double> gf(80.0, diag);

    double S; 
    for (int i = 0; i < cavity.size(); ++i) {
    	S = gf.diagonalS(cavity.elements(0)); 
     	std::cout << "S_{" << i << ", " << i << "} = " << S << std::endl;
    }
    double ref_S = 0.0; 
    BOOST_REQUIRE_CLOSE(ref_S, S, 1.0e-12);
}*/

/*! \class NumericalIntegrator 
 *  \test \b NumericalIntegratorTest_anisotropic tests the numerical evaluation of the anisotropic liquid diagonal elements of S and D
 */
BOOST_FIXTURE_TEST_CASE(anisotropic, NumericalIntegratorTest)
{
    AnisotropicLiquid<double> gf(epsilon, euler, diag);

    double S, accum; 
    accum = 0.0;
    for (int i = 0; i < cavity.size(); ++i) {
	std::cout << "                                        Tessera n. " << i+1 << std::endl;
    	S = gf.diagonalS(cavity.elements(i)); 
    }
    double ref_S = 0.0; 
    BOOST_REQUIRE_CLOSE(ref_S, S, 1.0e-12);
}
