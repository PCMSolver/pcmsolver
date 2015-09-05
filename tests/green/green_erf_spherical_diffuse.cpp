#define BOOST_TEST_MODULE GreensFunctionErfSphericalDiffuse

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
#include "OneLayerErf.hpp"

struct ErfSphericalDiffuseTest {
    int maxL;
    double eps1, eps2, sphereRadius, width;
    double inside_reference, outside_reference;
    double der_probe_inside_reference, der_source_inside_reference;
    double der_probe_outside_reference, der_source_outside_reference;
    Eigen::Vector3d sphereCenter;
    Eigen::Vector3d source1, probe1, sourceNormal1, probeNormal1;
    Eigen::Vector3d source2, probe2, sourceNormal2, probeNormal2;
    ErfSphericalDiffuseTest() { SetUp(); }
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
        inside_reference = 0.012507311377769814;
        der_probe_inside_reference = -0.012506340818360315;
	    der_source_inside_reference = 0.012498621813619368;
	    // Evaluation outside the sphere
        source2 << 150.0, 150.0, 150.0;
        sourceNormal2 = source2;
        sourceNormal2.normalize();
        probe2 << 151.0, 150.0, 150.0;
        probeNormal2 = probe2;
        probeNormal2.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        outside_reference = 0.49991229576650942;
        der_probe_outside_reference = -0.28997553534915177;
        der_source_outside_reference = 0.28871292628573908;
    }
};

BOOST_FIXTURE_TEST_SUITE(ErfSphericalDiffuse, ErfSphericalDiffuseTest)
/*! \class SphericalDiffuse
 *  \test \b ErfSphericalDiffuseTest_inside tests the evaluation of the ErfSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(inside, ErfSphericalDiffuseTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerErf> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
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
 *  \test \b ErfSphericalDiffuseTest_outside tests the evaluation of the ErfSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(outside, ErfSphericalDiffuseTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerErf> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
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

struct ErfSphericalDiffuseShiftedTest {
    int maxL;
    double eps1, eps2, sphereRadius, width;
    double inside_reference, outside_reference;
    double der_probe_inside_reference, der_source_inside_reference;
    double der_probe_outside_reference, der_source_outside_reference;
    Eigen::Vector3d sphereCenter;
    Eigen::Vector3d source1, probe1, sourceNormal1, probeNormal1;
    Eigen::Vector3d source2, probe2, sourceNormal2, probeNormal2;
    ErfSphericalDiffuseShiftedTest() { SetUp(); }
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
        source1 << 1.0, 0.0, 0.0;
        sourceNormal1 = source1;
        sourceNormal1.normalize();
        probe1 << 2.0, 0.0, 0.0;
        probeNormal1 = probe1;
        probeNormal1.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        inside_reference = 0.012523344896520634;
        der_probe_inside_reference = -0.012502439280500169;
        der_source_inside_reference = 0.012489903372303254;
	    // Evaluation outside the sphere
        source2 << 150.0, 150.0, 150.0;
        sourceNormal2 = source2;
        sourceNormal2.normalize();
        probe2 << 151.0, 150.0, 150.0;
        probeNormal2 = probe2;
        probeNormal2.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        outside_reference = 0.49989736527661349;
        der_probe_outside_reference = -0.28995345687399254;
        der_source_outside_reference = 0.28869402146192158;
    }
};

BOOST_FIXTURE_TEST_SUITE(ErfSphericalDiffuseShifted1, ErfSphericalDiffuseShiftedTest)

/*! \class SphericalDiffuse
 *  \test \b ErfSphericalDiffuseTest_inside tests the evaluation of the ErfSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(inside, ErfSphericalDiffuseShiftedTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerErf> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
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
 *  \test \b ErfSphericalDiffuseTest_outside tests the evaluation of the ErfSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(outside, ErfSphericalDiffuseShiftedTest)
{
    SphericalDiffuse<CollocationIntegrator, OneLayerErf> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
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
