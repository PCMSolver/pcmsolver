#include <Eigen/Dense>

#include "SurfaceFunction.h"
#include "gtest/gtest.h"

class SurfaceFunctionTest : public ::testing::Test
{
	protected:
		SurfaceFunction func1;
		SurfaceFunction func2;
		virtual void SetUp()
		{
			int nPoints = 10;
			Eigen::VectorXd values1(nPoints);
			Eigen::VectorXd values2(nPoints);
			values1 = Eigen::VectorXd::Constant(10,1,1.0);
			values2 = Eigen::VectorXd::Constant(10,1,4.0);
			double * values1_ptr = values1.data();
			double * values2_ptr = values2.data();
			func1 = SurfaceFunction("TestFunction1", nPoints, values1_ptr);
			func2 = SurfaceFunction("TestFunction2", nPoints, values2_ptr);
		}
};

TEST_F(SurfaceFunctionTest, addition)
{
	SurfaceFunction addition = func1 + func2;
	EXPECT_EQ("TestFunction1+TestFunction2", addition.getName());
	EXPECT_EQ(10, addition.getNPoints());
	Eigen::VectorXd result(10);
	result = Eigen::VectorXd::Constant(10, 1, 5.0);
	for (int i = 0; i < 10; ++i)
	{
		EXPECT_EQ(result(i), addition.getValue(i));
	}
}

TEST_F(SurfaceFunctionTest, subtraction)
{
	SurfaceFunction subtraction = func1 - func2;
	EXPECT_EQ("TestFunction1-TestFunction2", subtraction.getName());
	EXPECT_EQ(10, subtraction.getNPoints());
	Eigen::VectorXd result(10);
	result = Eigen::VectorXd::Constant(10, 1, -3.0);
	for (int i = 0; i < 10; ++i)
	{
		EXPECT_EQ(result(i), subtraction.getValue(i));
	}
}

TEST_F(SurfaceFunctionTest, multiply_by_scalar)
{
	SurfaceFunction scaled1 = 2.5 * func1;
	func2 *= 0.5;
	EXPECT_EQ("2.500000*TestFunction1", scaled1.getName());
	EXPECT_EQ("0.500000*TestFunction2", func2.getName());
	EXPECT_EQ(10, scaled1.getNPoints());
	EXPECT_EQ(10, func2.getNPoints());
	Eigen::VectorXd result1(10), result2(10);
	result1 = Eigen::VectorXd::Constant(10, 1, 2.5);
	result2 = Eigen::VectorXd::Constant(10, 1, 2.0);
	for (int i = 0; i < 10; ++i)
	{
		EXPECT_EQ(result1(i), scaled1.getValue(i));
		EXPECT_EQ(result2(i), func2.getValue(i));
	}
}

TEST_F(SurfaceFunctionTest, multiply)
{
	double product = func1 * func2;
	Eigen::VectorXd result1(10), result2(10);
	result1 = Eigen::VectorXd::Constant(10, 1, 1.0);
	result2 = Eigen::VectorXd::Constant(10, 1, 4.0);
	double expected_product = result1.dot(result2);
	EXPECT_EQ(expected_product, product);
}
