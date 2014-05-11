#define BOOST_TEST_MODULE GreensFunctionAnisotropicLiquid

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "AnisotropicLiquid.hpp"

struct AnisotropicLiquidTest {
    Eigen::Array4d analyticEvaluate(double eps, double k,
                                    const Eigen::Vector3d & spNormal,
                                    const Eigen::Vector3d & sp,
                                    const Eigen::Vector3d & ppNormal, const Eigen::Vector3d & pp) {
        Eigen::Array4d result = Eigen::Array4d::Zero();
        double distance = (sp - pp).norm();
        double distance_3 = std::pow(distance, 3);
        double distance_5 = std::pow(distance, 5);

        // Value of the function
        result(0) = std::exp(- k * distance) / (eps * distance);
        // Value of the directional derivative wrt probe
        result(1) = (sp - pp).dot(ppNormal) * (1 + k * distance ) * std::exp(
                        - k * distance) / (eps * distance_3);
        // Directional derivative wrt source
        result(2) = - (sp - pp).dot(spNormal) * (1 + k * distance ) * std::exp(
                        - k * distance) / (eps * distance_3);
        // Value of the Hessian
        result(3) = spNormal.dot(ppNormal) * (1 + k * distance) * std::exp(
                        - k * distance) / (eps * distance_3)
                    - std::pow(k, 2) * (sp - pp).dot(spNormal) * (sp - pp).dot(
                        ppNormal) * std::exp(- k * distance) / (eps * distance_3)
                    - 3 * (sp - pp).dot(spNormal) * (sp - pp).dot(
                        ppNormal) * (1 + k * distance) * std::exp(- k * distance) /
                    (eps * distance_5);

        return result;
    }
    Eigen::Vector3d epsilon;
    Eigen::Vector3d euler;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    AnisotropicLiquidTest() { SetUp(); }
    void SetUp() {
    	epsilon << 10.0, 5.0, 3.0;
    	euler << 40.0, 50.0, 70.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticEvaluate(epsilon(0), euler(0), sourceNormal, source, probeNormal, probe);
    }
};

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidTest_numerical tests the numerical evaluation of the AnisotropicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(numerical, AnisotropicLiquidTest)
{
    Eigen::Array4d result = analyticEvaluate(epsilon(0), euler(0), sourceNormal, source,
                            probeNormal,
                            probe);

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
    Eigen::Array4d result = analyticEvaluate(epsilon(0), euler(0), sourceNormal, source,
                            probeNormal,
                            probe);

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
    Eigen::Array4d result = analyticEvaluate(epsilon(0), euler(0), sourceNormal, source,
                            probeNormal,
                            probe);

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
    Eigen::Array4d result = analyticEvaluate(epsilon(0), euler(0), sourceNormal, source,
                            probeNormal,
                            probe);

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
