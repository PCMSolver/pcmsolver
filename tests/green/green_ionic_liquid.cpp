#include <cmath>
#include <iostream>

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "IonicLiquid.hpp"

#include "gtestPimpl.hpp"

class IonicLiquidTest : public ::testing::Test
{
public:
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
protected:
    double epsilon;
    double kappa;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    virtual void SetUp() {
        epsilon = 60.0;
        kappa = 5.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticEvaluate(epsilon, kappa, sourceNormal, source, probeNormal, probe);
    }
};

/*! \class IonicLiquid
 *  \test \b IonicLiquidTest_numerical tests the numerical evaluation of the IonicLiquid Green's function against analytical result
 */
TEST_F(IonicLiquidTest, numerical)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, kappa, sourceNormal, source,
                            probeNormal,
                            probe);

    IonicLiquid<double> gf(epsilon, kappa);
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

/*! \class IonicLiquid
 *  \test \b IonicLiquidTest_AD_directional tests the automatic evaluation (directional derivative only) 
 *  of the IonicLiquid Green's function against analytical result
 */
TEST_F(IonicLiquidTest, AD_directional)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, kappa, sourceNormal, source,
                            probeNormal,
                            probe);

    IonicLiquid<AD_directional> gf(epsilon, kappa);
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

/*! \class IonicLiquid
 *  \test \b IonicLiquidTest_AD_gradient tests the automatic evaluation (full gradient) 
 *  of the IonicLiquid Green's function against analytical result
 */
TEST_F(IonicLiquidTest, AD_gradient)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, kappa, sourceNormal, source,
                            probeNormal,
                            probe);

    IonicLiquid<AD_gradient> gf(epsilon, kappa);
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

/*! \class IonicLiquid
 *  \test \b IonicLiquidTest_AD_hessian tests the automatic evaluation (full hessian) 
 *  of the IonicLiquid Green's function against analytical result
 */
TEST_F(IonicLiquidTest, AD_hessian)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, kappa, sourceNormal, source,
                            probeNormal,
                            probe);

    IonicLiquid<AD_hessian> gf(epsilon, kappa);
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
