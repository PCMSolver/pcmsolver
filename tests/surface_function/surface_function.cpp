#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "SurfaceFunction.hpp"

#include "gtestPimpl.hpp"

class SurfaceFunctionTest : public ::testing::Test
{
public:
    int nPoints;
    Eigen::VectorXd values1;
    Eigen::VectorXd values2;
protected:
    SurfaceFunction func1;
    SurfaceFunction func2;
    virtual void SetUp() {
        nPoints = 1000;
        values1.resize(nPoints);
        values2.resize(nPoints);
        values1 = Eigen::VectorXd::Random(nPoints);
        values2 = Eigen::VectorXd::Random(nPoints);
        double * values1_ptr = values1.data();
        double * values2_ptr = values2.data();
        func1 = SurfaceFunction("TestFunction1", nPoints, values1_ptr);
        func2 = SurfaceFunction("TestFunction2", nPoints, values2_ptr);
    }
};

/*! \class SurfaceFunction
 *  \test \b SurfaceFunctionTest_addition tests addition of two SurfaceFunction
 */
TEST_F(SurfaceFunctionTest, addition)
{
    SurfaceFunction addition = func1 + func2;
    EXPECT_EQ("TestFunction1+TestFunction2", addition.getName());
    EXPECT_EQ(nPoints, addition.getNPoints());
    Eigen::VectorXd result(nPoints);
    result = values1 + values2;
    for (int i = 0; i < nPoints; ++i) {
        EXPECT_DOUBLE_EQ(result(i), addition.getValue(i));
    }
}

/*! \class SurfaceFunction
 *  \test \b SurfaceFunctionTest_subtraction tests subtraction of two SurfaceFunction
 */
TEST_F(SurfaceFunctionTest, subtraction)
{
    SurfaceFunction subtraction = func1 - func2;
    EXPECT_EQ("TestFunction1-TestFunction2", subtraction.getName());
    EXPECT_EQ(nPoints, subtraction.getNPoints());
    Eigen::VectorXd result(nPoints);
    result = values1 - values2;
    for (int i = 0; i < nPoints; ++i) {
        EXPECT_DOUBLE_EQ(result(i), subtraction.getValue(i));
    }
}

/*! \class SurfaceFunction
 *  \test \b SurfaceFunctionTest_multiply_by_scalar tests multiplication of a SurfaceFunction by a scalar
 */
TEST_F(SurfaceFunctionTest, multiply_by_scalar)
{
    SurfaceFunction scaled1 = 2.5 * func1;
    func2 *= 0.5;
    EXPECT_EQ("2.5*TestFunction1", scaled1.getName());
    EXPECT_EQ("0.5*TestFunction2", func2.getName());
    EXPECT_EQ(nPoints, scaled1.getNPoints());
    EXPECT_EQ(nPoints, func2.getNPoints());
    Eigen::VectorXd result1(nPoints), result2(nPoints);
    result1 = 2.5 * values1;
    result2 = 0.5 * values2;
    for (int i = 0; i < nPoints; ++i) {
        EXPECT_DOUBLE_EQ(result1(i), scaled1.getValue(i));
        EXPECT_DOUBLE_EQ(result2(i), func2.getValue(i));
    }
}

/*! \class SurfaceFunction
 *  \test \b SurfaceFunctionTest_multiply tests dot product of two SurfaceFunction
 */
TEST_F(SurfaceFunctionTest, multiply)
{
    double product = func1 * func2;
    double expected_product = values1.dot(values2);
    EXPECT_DOUBLE_EQ(expected_product, product);
}
