#define BOOST_TEST_MODULE GreensFunctionTanhSphericalDiffuse

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include "Config.hpp"

#include <Eigen/Dense>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "SphericalDiffuse.hpp"
#include "OneLayerTanh.hpp"

struct TanhSphericalDiffuseTest {
    double eps1, eps2, sphereRadius, width;
    double inside_reference, outside_reference;
    double der_probe_inside_reference, der_source_inside_reference;
    double der_probe_outside_reference, der_source_outside_reference;
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
        der_probe_inside_reference = -0.0125248562596845525;
        der_source_inside_reference = 0.0125498847608398328;
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
        der_probe_outside_reference = -0.289954122292268579;
        der_source_outside_reference = 0.288674001099886723;
    }
};

BOOST_FIXTURE_TEST_SUITE(TanhSphericalDiffuse, TanhSphericalDiffuseTest)
/*! \class SphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_inside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(inside, TanhSphericalDiffuseTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(eps1, eps2, width, sphereRadius, sphereCenter);
    double value = inside_reference;
    double gf_value = gf.kernelS(source1, probe1);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, value, 1.0e-08);

    double der_probe = der_probe_inside_reference;
    double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("der_probe   = " << std::setprecision(std::numeric_limits<long double>::digits10) << der_probe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_REQUIRE_CLOSE(gf_derProbe, der_probe, 1.0e-06);

    double der_source = der_source_inside_reference;
    double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("der_source   = " << std::setprecision(std::numeric_limits<long double>::digits10) << der_source);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(gf_derSource, der_source, 1.0e-06);
}
/*! \class SphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_outside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(outside, TanhSphericalDiffuseTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(eps1, eps2, width, sphereRadius, sphereCenter);
    double value = outside_reference;
    double gf_value = gf.kernelS(source2, probe2);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, value, 1.0e-08);

    double der_probe = der_probe_outside_reference;
    double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("der_probe   = " << std::setprecision(std::numeric_limits<long double>::digits10) << der_probe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_REQUIRE_CLOSE(gf_derProbe, der_probe, 1.0e-06);

    double der_source = der_source_outside_reference;
    double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("der_source   = " << std::setprecision(std::numeric_limits<long double>::digits10) << der_source);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(gf_derSource, der_source, 1.0e-06);
}

BOOST_AUTO_TEST_SUITE_END()

struct TanhSphericalDiffuseShiftedTest {
    double eps1, eps2, sphereRadius, width;
    double inside_reference, outside_reference;
    double der_probe_inside_reference, der_source_inside_reference;
    double der_probe_outside_reference, der_source_outside_reference;
    Eigen::Vector3d sphereCenter;
    Eigen::Vector3d source1, probe1, sourceNormal1, probeNormal1;
    Eigen::Vector3d source2, probe2, sourceNormal2, probeNormal2;
    TanhSphericalDiffuseShiftedTest() { SetUp(); }
    void SetUp() {
        // High dielectric constant inside
        eps1 = 80.0;
        // Low dielectric constant outside
        eps2 = 2.0;
	    sphereCenter << 0.0, 0.0, 300.0;
        sphereRadius = 100.0;
	    width = 5.0;
	    // Evaluation inside the sphere
        source1 << 0.0, 0.0, 301.0;
        sourceNormal1 = source1; // + Eigen::Vector3d::Random();
        sourceNormal1.normalize();
        probe1 << 0.0, 0.0, 302.0;
        probeNormal1 = probe1; // + Eigen::Vector3d::Random();
        probeNormal1.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        inside_reference = 0.012377848015483936;
        der_probe_inside_reference = -0.0125248562565360294;
        der_source_inside_reference = 0.0125498847577173306;
	    // Evaluation outside the sphere
        source2 << 150.0, 150.0, 150.0;
        sourceNormal2 = source2; // + Eigen::Vector3d::Random();
        sourceNormal2.normalize();
        probe2 << 150.0, 150.0, 151.0;
        probeNormal2 = probe2; // + Eigen::Vector3d::Random();
        probeNormal2.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        outside_reference = 0.5000891977645715;
        der_probe_outside_reference = -0.289954254408808509;
        der_source_outside_reference = 0.288674994271542751;
    }
};

BOOST_FIXTURE_TEST_SUITE(TanhSphericalDiffuseShifted1, TanhSphericalDiffuseShiftedTest)

/*! \class SphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_inside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(inside, TanhSphericalDiffuseShiftedTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(eps1, eps2, width, sphereRadius, sphereCenter);
    double value = inside_reference;
    double gf_value = gf.kernelS(source1, probe1);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, value, 1.0e-08);

    double der_probe = der_probe_inside_reference;
    double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("der_probe   = " << std::setprecision(std::numeric_limits<long double>::digits10) << der_probe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_REQUIRE_CLOSE(gf_derProbe, der_probe, 1.0e-06);

    double der_source = der_source_inside_reference;
    double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("der_source   = " << std::setprecision(std::numeric_limits<long double>::digits10) << der_source);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(gf_derSource, der_source, 1.0e-06);
}
/*! \class SphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_outside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(outside, TanhSphericalDiffuseShiftedTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(eps1, eps2, width, sphereRadius, sphereCenter);
    double value = outside_reference;
    double gf_value = gf.kernelS(source2, probe2);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, value, 1.0e-08);

    double der_probe = der_probe_outside_reference;
    double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("der_probe   = " << std::setprecision(std::numeric_limits<long double>::digits10) << der_probe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_REQUIRE_CLOSE(gf_derProbe, der_probe, 1.0e-06);

    double der_source = der_source_outside_reference;
    double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("der_source   = " << std::setprecision(std::numeric_limits<long double>::digits10) << der_source);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(gf_derSource, der_source, 1.0e-06);
}

BOOST_AUTO_TEST_SUITE_END()
