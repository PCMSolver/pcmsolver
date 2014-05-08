#define BOOST_TEST_MODULE WaveletCavity

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "PWCSolver.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "WaveletCavity.hpp"

struct WaveletCavityTest {
    WaveletCavity cavity;
    WaveletCavityTest() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d N(0.0, 0.0, 0.0);
        std::vector<Sphere> spheres;
        Sphere sph1(N, 1.0);
        spheres.push_back(sph1);
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
 *  \test \b WaveletCavityTest_size tests Wavelet cavity size for a point charge
 */
BOOST_FIXTURE_TEST_CASE(size, WaveletCavityTest)
{
    int size = 4864;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class WaveletCavity
 *  \test \b WaveletCavityTest_area tests Wavelet cavity surface area for a point charge
 */
BOOST_FIXTURE_TEST_CASE(area, WaveletCavityTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class WaveletCavity
 *  \test \b WaveletCavityTest_volume tests Wavelet cavity volume for a point charge
 */
BOOST_FIXTURE_TEST_CASE(volume, WaveletCavityTest)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}
