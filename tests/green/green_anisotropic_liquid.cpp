#define BOOST_TEST_MODULE GreensFunctionAnisotropicLiquid

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "AnalyticEvaluate.hpp"
#include "DerivativeTypes.hpp"
#include "AnisotropicLiquid.hpp"

struct AnisotropicLiquidTest {
    Eigen::Vector3d epsilon;
    Eigen::Vector3d euler;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    AnisotropicLiquidTest() { SetUp(); }
    void SetUp() {
    	epsilon << 2.0, 80.0, 15.0;
    	euler << 6.0, 40.0, 15.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticAnisotropicLiquid(epsilon, euler, sourceNormal, source, probeNormal, probe);
    }
};

BOOST_FIXTURE_TEST_SUITE(Anisotropic, AnisotropicLiquidTest)

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidTest_numerical tests the numerical evaluation of the AnisotropicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(numerical, AnisotropicLiquidTest)
{
    AnisotropicLiquid<double> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-05);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-05);
}

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidTest_directional_AD tests the automatic evaluation (directional derivative only)
 *  of the AnisotropicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(directional_AD, AnisotropicLiquidTest)
{
    AnisotropicLiquid<AD_directional> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);
}

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidTest_gradient_AD tests the automatic evaluation (full gradient)
 *  of the AnisotropicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(gradient_AD, AnisotropicLiquidTest)
{
    AnisotropicLiquid<AD_gradient> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);
}

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidTest_hessian_AD tests the automatic evaluation (full hessian)
 *  of the AnisotropicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(hessian_AD, AnisotropicLiquidTest)
{
    AnisotropicLiquid<AD_hessian> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);

    /*	double hessian = result(4);
    	double gf_hessian = gf.hessian(sourceNormal, source, probeNormal, probe);
    	BOOST_REQUIRE_CLOSE(hessian, gf_hessian, 1.0e-12);*/
}

BOOST_AUTO_TEST_SUITE_END()

// Define a uniform dielectric as an anisotropic dielectric and test the evaluation of the
// Green's function and derivatives against the analytic result for the uniform dielectric
struct AnisotropicLiquidUniformTest {
    Eigen::Vector3d epsilon;
    Eigen::Vector3d euler;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    AnisotropicLiquidUniformTest() { SetUp(); }
    void SetUp() {
    	epsilon << 80.0, 80.0, 80.0;
    	euler << 0.0, 0.0, 0.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
	result = analyticUniformDielectric(epsilon(0), sourceNormal, source, probeNormal, probe);
    }
};

BOOST_FIXTURE_TEST_SUITE(AnisotropicVsUniform, AnisotropicLiquidUniformTest)

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidUniformTest_numerical tests the numerical evaluation of the AnisotropicLiquid Green's function against analytical result for a uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(numerical, AnisotropicLiquidUniformTest)
{
    AnisotropicLiquid<double> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-05);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-05);
}

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidUniformTest_directional_AD tests the automatic evaluation (directional derivative only)
 *  of the AnisotropicLiquid Green's function against analytical result for a uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(directional_AD, AnisotropicLiquidUniformTest)
{
    AnisotropicLiquid<AD_directional> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);
}

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidUniformTest_gradient_AD tests the automatic evaluation (full gradient)
 *  of the AnisotropicLiquid Green's function against analytical result for a uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(gradient_AD, AnisotropicLiquidUniformTest)
{
    AnisotropicLiquid<AD_gradient> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);
}

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidUniformTest_hessian_AD tests the automatic evaluation (full hessian)
 *  of the AnisotropicLiquid Green's function against analytical result for a uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(hessian_AD, AnisotropicLiquidUniformTest)
{
    AnisotropicLiquid<AD_hessian> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);

    /*	double hessian = result(4);
    	double gf_hessian = gf.hessian(sourceNormal, source, probeNormal, probe);
    	BOOST_REQUIRE_CLOSE(hessian, gf_hessian, 1.0e-12);*/
}

BOOST_AUTO_TEST_SUITE_END()
