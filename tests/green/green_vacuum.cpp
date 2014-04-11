#include <cmath>
#include <iostream>

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "Vacuum.hpp"

// Disable obnoxious warnings from Google Test headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#include "gtest/gtest.h"
#pragma GCC diagnostic pop
#endif

class VacuumTest : public ::testing::Test
{
public:
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
protected:
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    virtual void SetUp() {
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticEvaluate(sourceNormal, source, probeNormal, probe);
    }
};

TEST_F(VacuumTest, numerical)
{
    Vacuum<double> gf;
    double value = result(0);
    double gf_value = gf.function(source, probe);
    EXPECT_DOUBLE_EQ(value, gf_value);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    EXPECT_NEAR(derProbe, gf_derProbe, 1.0e-09);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    EXPECT_NEAR(derSource, gf_derSource, 1.0e-09);
}

TEST_F(VacuumTest, AD_directional)
{
    Vacuum<AD_directional> gf;
    double value = result(0);
    double gf_value = gf.function(source, probe);
    EXPECT_DOUBLE_EQ(value, gf_value);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    EXPECT_DOUBLE_EQ(derProbe, gf_derProbe);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    EXPECT_DOUBLE_EQ(derSource, gf_derSource);
}

TEST_F(VacuumTest, AD_gradient)
{
    Vacuum<AD_gradient> gf;
    double value = result(0);
    double gf_value = gf.function(source, probe);
    EXPECT_DOUBLE_EQ(value, gf_value);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    EXPECT_DOUBLE_EQ(derProbe, gf_derProbe);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    EXPECT_DOUBLE_EQ(derSource, gf_derSource);
}

TEST_F(VacuumTest, AD_hessian)
{
    Vacuum<AD_hessian> gf;
    double value = result(0);
    double gf_value = gf.function(source, probe);
    EXPECT_DOUBLE_EQ(value, gf_value);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    EXPECT_DOUBLE_EQ(derProbe, gf_derProbe);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    EXPECT_DOUBLE_EQ(derSource, gf_derSource);

    /*	double hessian = result(4);
    	double gf_hessian = gf.hessian(sourceNormal, source, probeNormal, probe);
    	EXPECT_DOUBLE_EQ(hessian, gf_hessian);
    	EXPECT_NEAR(hessian, gf_hessian, 1.0e-13);*/
}
