#define BOOST_TEST_MODULE WaveletCavityNH3

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include <Eigen/Dense>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "PWCSolver.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "WaveletCavity.hpp"

struct WaveletCavityNH3Test {
    WaveletCavity cavity;
    WaveletCavityNH3Test() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d N( -0.000000000,   -0.104038047,    0.000000000);
        Eigen::Vector3d H1(-0.901584415,    0.481847022,   -1.561590016);
        Eigen::Vector3d H2(-0.901584415,    0.481847022,    1.561590016);
        Eigen::Vector3d H3( 1.803168833,    0.481847022,    0.000000000);
        std::vector<Sphere> spheres;
        Sphere sph1(N,  2.929075493);
        Sphere sph2(H1, 2.267671349);
        Sphere sph3(H2, 2.267671349);
        Sphere sph4(H3, 2.267671349);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        spheres.push_back(sph4);
        double probeRadius = 1.385; // Probe Radius for water
        int patchLevel = 2;
        double coarsity = 0.5;
        cavity = WaveletCavity(spheres, probeRadius, patchLevel, coarsity);
        cavity.readCavity("molec_dyadic.dat");

	CollocationIntegrator * diag = new CollocationIntegrator();
        double permittivity = 78.39;
        Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
        UniformDielectric<AD_directional> * gfOutside = new
        UniformDielectric<AD_directional>(permittivity, diag);
        int firstKind = 0;
        PWCSolver solver(gfInside, gfOutside, firstKind);
        solver.buildSystemMatrix(cavity);
        cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_(), false);
    }
};

/*! \class WaveletCavity
 *  \test \b WaveletCavityNH3Test_size tests Wavelet cavity size for ammonia
 */
BOOST_FIXTURE_TEST_CASE(size, WaveletCavityNH3Test)
{
    int size = 4288;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class WaveletCavity
 *  \test \b WaveletCavityNH3Test_area tests Wavelet cavity surface area for ammonia
 */
BOOST_FIXTURE_TEST_CASE(area, WaveletCavityNH3Test)
{
    double area = 146.41490284471513;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class WaveletCavity
 *  \test \b WaveletCavityNH3Test_volume tests Wavelet cavity volume for ammonia
 */
BOOST_FIXTURE_TEST_CASE(volume, WaveletCavityNH3Test)
{
    double volume = 153.3346491517182;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}
