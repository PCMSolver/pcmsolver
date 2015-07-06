#define BOOST_TEST_MODULE GreensFunctionTanhSphericalDiffuse

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include "Config.hpp"

#include <Eigen/Core>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "SphericalDiffuse.hpp"
#include "OneLayerTanh.hpp"

struct TanhSphericalDiffuseTest {
    int maxL;
    double eps1, eps2, sphereRadius, width;
    double inside_reference, outside_reference;
    double der_probe_inside_reference, der_source_inside_reference;
    double der_probe_outside_reference, der_source_outside_reference;
    Eigen::Vector3d sphereCenter;
    Eigen::Vector3d source1, probe1, sourceNormal1, probeNormal1;
    Eigen::Vector3d source2, probe2, sourceNormal2, probeNormal2;
    TanhSphericalDiffuseTest() { SetUp(); }
    void SetUp() {
        maxL = 3;
        // High dielectric constant inside
        eps1 = 80.0;
        // Low dielectric constant outside
        eps2 = 2.0;
	    sphereCenter << 0.0, 0.0, 0.0;
        sphereRadius = 100.0;
	    width = 5.0;
	    // Evaluation inside the sphere
        source1 << 1.0, 0.0, 0.0;
        sourceNormal1 = source1;
        sourceNormal1.normalize();
        probe1 << 2.0, 0.0, 0.0;
        probeNormal1 = probe1;
        probeNormal1.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        inside_reference = 0.012507311388168523;
        der_probe_inside_reference = -0.012506305835640469;
	    der_source_inside_reference = 0.012498621932118328;
	    // Evaluation outside the sphere
        source2 << 150.0, 150.0, 150.0;
        sourceNormal2 = source2;
        sourceNormal2.normalize();
        probe2 << 151.0, 150.0, 150.0;
        probeNormal2 = probe2;
        probeNormal2.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        outside_reference = 0.50004567416494572;
        der_probe_outside_reference = -0.29005549287308696;
        der_source_outside_reference = 0.28879776627577236;
    }
};

BOOST_FIXTURE_TEST_SUITE(TanhSphericalDiffuse, TanhSphericalDiffuseTest)
/*! \class SphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_inside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(inside, TanhSphericalDiffuseTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
    double ref_value = inside_reference;
    double gf_value  = gf.kernelS(source1, probe1);
    BOOST_TEST_MESSAGE("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_value);
    BOOST_TEST_MESSAGE("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, ref_value, 1.0e-08);

    double ref_derProbe = der_probe_inside_reference;
    double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("ref_derProbe   = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_derProbe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_REQUIRE_CLOSE(gf_derProbe, ref_derProbe, 1.0e-06);

    double ref_derSource = der_source_inside_reference;
    double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("ref_derSource   = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(gf_derSource, ref_derSource, 1.0e-06);
}
/*! \class SphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_outside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(outside, TanhSphericalDiffuseTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
    double ref_value = outside_reference;
    double gf_value = gf.kernelS(source2, probe2);
    BOOST_TEST_MESSAGE("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_value);
    BOOST_TEST_MESSAGE("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, ref_value, 1.0e-08);

    double ref_derProbe = der_probe_outside_reference;
    double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("ref_derProbe   = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_derProbe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_REQUIRE_CLOSE(gf_derProbe, ref_derProbe, 1.0e-06);

    double ref_derSource = der_source_outside_reference;
    double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("ref_derSource   = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(gf_derSource, ref_derSource, 1.0e-06);
}

BOOST_AUTO_TEST_SUITE_END()

struct TanhSphericalDiffuseShiftedTest {
    int maxL;
    double eps1, eps2, sphereRadius, width;
    double inside_reference, outside_reference;
    double der_probe_inside_reference, der_source_inside_reference;
    double der_probe_outside_reference, der_source_outside_reference;
    Eigen::Vector3d sphereCenter;
    Eigen::Vector3d source1, probe1, sourceNormal1, probeNormal1;
    Eigen::Vector3d source2, probe2, sourceNormal2, probeNormal2;
    TanhSphericalDiffuseShiftedTest() { SetUp(); }
    void SetUp() {
        maxL = 3;
        // High dielectric constant inside
        eps1 = 80.0;
        // Low dielectric constant outside
        eps2 = 2.0;
	    sphereCenter << 25.0, 0.0, 0.0;
        sphereRadius = 100.0;
	    width = 5.0;
	    // Evaluation inside the sphere
        source1 << 0.0, 0.0, 1.0;
        sourceNormal1 = source1;
        sourceNormal1.normalize();
        probe1 << 0.0, 0.0, 2.0;
        probeNormal1 = probe1;
        probeNormal1.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        inside_reference = 0.012520055180345364;
        der_probe_inside_reference = -0.01250488850888104;
	    der_source_inside_reference = 0.012504949286531314;
	    // Evaluation outside the sphere
        source2 << 150.0, 150.0, 150.0;
        sourceNormal2 = source2;
        sourceNormal2.normalize();
        probe2 << 151.0, 150.0, 150.0;
        probeNormal2 = probe2;
        probeNormal2.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        outside_reference = 0.5000329900631173;
        der_probe_outside_reference = -0.28999970912518824;
        der_source_outside_reference = 0.28865758495810745;
    }
};

BOOST_FIXTURE_TEST_SUITE(TanhSphericalDiffuseShifted1, TanhSphericalDiffuseShiftedTest)

/*! \class SphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_inside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(inside, TanhSphericalDiffuseShiftedTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
    double ref_value = inside_reference;
    double gf_value = gf.kernelS(source1, probe1);
    BOOST_TEST_MESSAGE("ref_value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, ref_value, 1.0e-08);

    double ref_derProbe = der_probe_inside_reference;
    double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("ref_derProbe   = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_derProbe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_REQUIRE_CLOSE(gf_derProbe, ref_derProbe, 1.0e-06);

    double ref_derSource = der_source_inside_reference;
    double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("ref_derSource   = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(gf_derSource, ref_derSource, 1.0e-06);
}
/*! \class SphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_outside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(outside, TanhSphericalDiffuseShiftedTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
    double ref_value = outside_reference;
    double gf_value = gf.kernelS(source2, probe2);
    BOOST_TEST_MESSAGE("ref_value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_REQUIRE_CLOSE(gf_value, ref_value, 1.0e-08);

    double ref_derProbe = der_probe_outside_reference;
    double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("ref_derProbe   = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_derProbe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_REQUIRE_CLOSE(gf_derProbe, ref_derProbe, 1.0e-06);

    double ref_derSource = der_source_outside_reference;
    double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("ref_derSource   = " << std::setprecision(std::numeric_limits<long double>::digits10) << ref_derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_REQUIRE_CLOSE(gf_derSource, ref_derSource, 1.0e-06);
}

BOOST_AUTO_TEST_SUITE_END()
