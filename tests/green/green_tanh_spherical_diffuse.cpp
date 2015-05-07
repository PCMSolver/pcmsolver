#define BOOST_TEST_MODULE GreensFunctionTanhSphericalDiffuse

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include "Config.hpp"

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "TanhSphericalDiffuse.hpp"

struct TanhSphericalDiffuseTest {
    double eps1, eps2, sphereRadius, width;
    double inside_reference, outside_reference;
    double deriv_inside_reference, deriv_outside_reference;
    Eigen::Vector3d sphereCenter;
    Eigen::Vector3d source1, probe1, sourceNormal1, probeNormal1;
    Eigen::Vector3d source2, probe2, sourceNormal2, probeNormal2;
    TanhSphericalDiffuseTest() { SetUp(); }
    void SetUp() {
        // High dielectric constant inside
        eps1 = 80.0;
        // Low dielectric constant outside
        eps2 = 2.0;
	    sphereCenter << 0.0, 0.0, 0.0;
        sphereRadius = 100.0;
	    width = 5.0;
	    // Evaluation inside the sphere
        source1 << 1.0, 0.0, 0.0;
        sourceNormal1 = source1; // + Eigen::Vector3d::Random();
        sourceNormal1.normalize();
        probe1 << 2.0, 0.0, 0.0;
        probeNormal1 = probe1; // + Eigen::Vector3d::Random();
        probeNormal1.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        inside_reference = 0.012377848015483936;
        deriv_inside_reference = -1.0601898887311156;
	    // Evaluation outside the sphere
        source2 << 150.0, 150.0, 150.0;
        sourceNormal2 = source2; // + Eigen::Vector3d::Random();
        sourceNormal2.normalize();
        probe2 << 151.0, 150.0, 150.0;
        probeNormal2 = probe2; // + Eigen::Vector3d::Random();
        probeNormal2.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        outside_reference = 0.50008829802731714;
        deriv_outside_reference = -0.5799074235842071;
    }
};

BOOST_FIXTURE_TEST_SUITE(TanhSphericalDiffuse1, TanhSphericalDiffuseTest)

/*! \class TanhSphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_inside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(inside, TanhSphericalDiffuseTest)
{
    TanhSphericalDiffuse gf(eps1, eps2, width, sphereRadius);
    double value = inside_reference;
    double gf_value = gf.function(source1, probe1);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, value, 1.0e-08);

    double deriv = deriv_inside_reference;
    double gf_deriv = gf.derivative(probeNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("deriv    = " << std::setprecision(std::numeric_limits<long double>::digits10) << deriv);
    BOOST_TEST_MESSAGE("gf_deriv = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_deriv);
    BOOST_REQUIRE_CLOSE(gf_deriv, deriv, 1.0e-06);

    /*
    double derSource = resultInside(2);
    double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("derSource    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-06);
    */
}
/*! \class TanhSphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_outside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(outside, TanhSphericalDiffuseTest)
{
    TanhSphericalDiffuse gf(eps1, eps2, width, sphereRadius);
    double value = outside_reference;
    double gf_value = gf.function(source2, probe2);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, value, 1.0e-08);

    double deriv = deriv_outside_reference;
    double gf_deriv = gf.derivative(probeNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("deriv    = " << std::setprecision(std::numeric_limits<long double>::digits10) << deriv);
    BOOST_TEST_MESSAGE("gf_deriv = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_deriv);
    BOOST_REQUIRE_CLOSE(gf_deriv, deriv, 1.0e-06);

    /*
    double derSource = resultOutside(2);
    double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("derSource    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-06);
    */
}
BOOST_AUTO_TEST_SUITE_END()
