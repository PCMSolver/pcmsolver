#define BOOST_TEST_MODULE GreensFunctionVacuum

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "Vacuum.hpp"

struct VacuumTest {
    Eigen::Array4d analyticEvaluate(const Eigen::Vector3d & spNormal,
                                    const Eigen::Vector3d & sp,
                                    const Eigen::Vector3d & ppNormal, const Eigen::Vector3d & pp) {
        Eigen::Array4d result = Eigen::Array4d::Zero();
        double distance = (sp - pp).norm();
        double distance_3 = std::pow(distance, 3);
        double distance_5 = std::pow(distance, 5);

        // Value of the function
        result(0) = 1.0 / distance;
        // Value of the directional derivative wrt probe
        result(1) = (sp - pp).dot(ppNormal) / distance_3 ;
        // Directional derivative wrt source
        result(2) = - (sp - pp).dot(spNormal) / distance_3;
        // Value of the Hessian
        result(3) = spNormal.dot(ppNormal) / distance_3 - 3 * ((
                        sp - pp).dot(spNormal))*((sp - pp).dot(
                                    ppNormal)) / distance_5;

        return result;
    }
    VacuumTest() { SetUp(); }
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    void SetUp() {
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticEvaluate(sourceNormal, source, probeNormal, probe);
    }
};

/*! \class Vacuum
 *  \test \b VacuumTest_numerical tests the numerical evaluation of the Vacuum Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(numerical, VacuumTest)
{
    Vacuum<double> gf;
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-06);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-06);
}

/*! \class Vacuum
 *  \test \b VacuumTest_directional_AD tests the automatic evaluation (directional derivative only)
 *  of the Vacuum Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(directional_AD, VacuumTest)
{
    Vacuum<AD_directional> gf;
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

/*! \class Vacuum
 *  \test \b VacuumTest_gradient_AD tests the automatic evaluation (full gradient)
 *  of the Vacuum Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(gradient_AD, VacuumTest)
{
    Vacuum<AD_gradient> gf;
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

/*! \class Vacuum
 *  \test \b VacuumTest_hessian_AD tests the automatic evaluation (full hessian)
 *  of the Vacuum Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(hessian_AD, VacuumTest)
{
    Vacuum<AD_hessian> gf;
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
