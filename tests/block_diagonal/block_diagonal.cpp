#include <iostream>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

#include "MathUtils.hpp"

// Disable obnoxious warnings from Google Test headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include "gtest/gtest.h"
#pragma GCC diagnostic pop
#endif

class BlockDiagonalMatrixTest : public ::testing::Test
{
    public:
        int nPoints;
	Eigen::MatrixXd block1;
	Eigen::MatrixXd block2;
	Eigen::MatrixXd block3;
    protected:
	BlockDiagonalMatrix<double> bd;
	virtual void SetUp()
	{
            nPoints = 5;
	    block1 = Eigen::MatrixXd::Random(nPoints, nPoints);
	    block2 = Eigen::MatrixXd::Random(nPoints, nPoints);
	    block3 = Eigen::MatrixXd::Random(nPoints, nPoints);
	    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(3*nPoints, 3*nPoints);
	    matrix.block(0, 0, nPoints, nPoints) = block1;
	    matrix.block(nPoints, nPoints, nPoints, nPoints) = block2;
	    matrix.block(2*nPoints, 2*nPoints, nPoints, nPoints) = block3;
	    std::cout << matrix << std::endl;
	    bd = BlockDiagonalMatrix<double>(matrix, 3, nPoints);
	}
};

TEST_F(BlockDiagonalMatrixTest, printout)
{
	std::cout << bd << std::endl;
}

/*
TEST_F(SurfaceFunctionTest, addition)
{
	SurfaceFunction addition = func1 + func2;
	EXPECT_EQ("TestFunction1+TestFunction2", addition.getName());
	EXPECT_EQ(nPoints, addition.getNPoints());
	Eigen::VectorXd result(nPoints);
    result = values1 + values2;
	for (int i = 0; i < nPoints; ++i)
	{
		EXPECT_EQ(result(i), addition.getValue(i));
	}
}

TEST_F(SurfaceFunctionTest, subtraction)
{
	SurfaceFunction subtraction = func1 - func2;
	EXPECT_EQ("TestFunction1-TestFunction2", subtraction.getName());
	EXPECT_EQ(nPoints, subtraction.getNPoints());
	Eigen::VectorXd result(nPoints);
    result = values1 - values2;
	for (int i = 0; i < nPoints; ++i)
	{
		EXPECT_EQ(result(i), subtraction.getValue(i));
	}
}

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
	for (int i = 0; i < nPoints; ++i)
	{
		EXPECT_EQ(result1(i), scaled1.getValue(i));
		EXPECT_EQ(result2(i), func2.getValue(i));
	}
}

TEST_F(SurfaceFunctionTest, multiply)
{
	double product = func1 * func2;
	double expected_product = values1.dot(values2);
	EXPECT_EQ(expected_product, product);
}*/
