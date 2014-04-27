#include <cmath>
#include <iostream>

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "UniformDielectric.hpp"

#include "gtestPimpl.hpp"

class UniformDielectricTest : public ::testing::Test
{
public:
    Eigen::Array4d analyticEvaluate(double eps, const Eigen::Vector3d & spNormal,
                                    const Eigen::Vector3d & sp,
                                    const Eigen::Vector3d & ppNormal, const Eigen::Vector3d & pp) {
        Eigen::Array4d result = Eigen::Array4d::Zero();
        double distance = (sp - pp).norm();
        double distance_3 = std::pow(distance, 3);
        double distance_5 = std::pow(distance, 5);

        // Value of the function
        result(0) = 1.0 / (eps * distance);
        // Value of the directional derivative wrt probe
        result(1) = (sp - pp).dot(ppNormal) / (eps * distance_3);
        // Directional derivative wrt source
        result(2) = - (sp - pp).dot(spNormal) / (eps * distance_3);
        // Value of the Hessian
        result(3) = spNormal.dot(ppNormal) / (eps * distance_3) - 3 * ((
                        sp - pp).dot(spNormal))*((sp - pp).dot(
                                    ppNormal)) / (eps * distance_5);

        return result;
    }
protected:
    double epsilon;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    virtual void SetUp() {
        epsilon = 60.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticEvaluate(epsilon, sourceNormal, source, probeNormal, probe);
    }
};

TEST_F(UniformDielectricTest, numerical)
{
    UniformDielectric<double> gf(epsilon);
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

TEST_F(UniformDielectricTest, AD_directional)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, sourceNormal, source, probeNormal,
                            probe);

    UniformDielectric<AD_directional> gf(epsilon);
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

TEST_F(UniformDielectricTest, AD_gradient)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, sourceNormal, source, probeNormal,
                            probe);

    UniformDielectric<AD_gradient> gf(epsilon);
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

TEST_F(UniformDielectricTest, AD_hessian)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, sourceNormal, source, probeNormal,
                            probe);

    UniformDielectric<AD_hessian> gf(epsilon);
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
