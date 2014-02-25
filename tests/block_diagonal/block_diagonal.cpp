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
	int nBlocks;
        int nPoints;
	int fullDim;
	Eigen::MatrixXd block1;
	Eigen::MatrixXd block2;
	Eigen::MatrixXd block3;
    protected:
	BlockDiagonalMatrix<double, 3, 15> bd1;
	BlockDiagonalMatrix<double, 3, 15> bd2;
	virtual void SetUp()
	{
	    nBlocks = 3;
	    nPoints = 15;
	    fullDim = nBlocks * nPoints;
	    block1 = Eigen::MatrixXd::Random(nPoints, nPoints);
	    block2 = Eigen::MatrixXd::Random(nPoints, nPoints);
	    block3 = Eigen::MatrixXd::Random(nPoints, nPoints);
	    Eigen::MatrixXd matrix1 = Eigen::MatrixXd::Zero(3*nPoints, 3*nPoints);
	    matrix1.block(0, 0, nPoints, nPoints) = block1;
	    matrix1.block(nPoints, nPoints, nPoints, nPoints) = block2;
	    matrix1.block(2*nPoints, 2*nPoints, nPoints, nPoints) = block3;
	    bd1 = BlockDiagonalMatrix<double, 3, 15>(matrix1);
	    Eigen::MatrixXd matrix2 = Eigen::MatrixXd::Zero(3*nPoints, 3*nPoints);
	    matrix2.block(0, 0, nPoints, nPoints) = block3;
	    matrix2.block(nPoints, nPoints, nPoints, nPoints) = block2;
	    matrix2.block(2*nPoints, 2*nPoints, nPoints, nPoints) = block1;
	    bd2 = BlockDiagonalMatrix<double, 3, 15>(matrix2);
	}
};

TEST_F(BlockDiagonalMatrixTest, addition)
{
	const int nBlocks = 3;
        const int nPoints = 15;
	BlockDiagonalMatrix<double, nBlocks, nPoints> addition = bd1 + bd2;
	EXPECT_EQ(fullDim, addition.fullDim());
	Eigen::MatrixXd result1(nPoints, nPoints);
	result1 = block1 + block3;
	Eigen::MatrixXd result2(nPoints, nPoints);
	result2 = block2 + block2;
	Eigen::MatrixXd result3(nPoints, nPoints);
	result3 = block3 + block1;
	Eigen::MatrixXd tmp1(nPoints, nPoints);
	tmp1 = addition.block(0);
	Eigen::MatrixXd tmp2(nPoints, nPoints);
	tmp2 = addition.block(1);
	Eigen::MatrixXd tmp3(nPoints, nPoints);
	tmp3 = addition.block(2);
	for (int i = 0; i < nPoints; ++i)
	{
		for (int j = 0; j < nPoints; ++j)
		{
			EXPECT_EQ(result1(i, j), tmp1(i, j));
			EXPECT_EQ(result2(i, j), tmp2(i, j));
			EXPECT_EQ(result3(i, j), tmp3(i, j));
		}
	}
}

TEST_F(BlockDiagonalMatrixTest, subtraction)
{
	const int nBlocks = 3;
        const int nPoints = 15;
	BlockDiagonalMatrix<double, nBlocks, nPoints> subtraction = bd1 - bd2;
	EXPECT_EQ(fullDim, subtraction.fullDim());
	Eigen::MatrixXd result1(nPoints, nPoints);
	result1 = block1 - block3;
	Eigen::MatrixXd result2(nPoints, nPoints);
	result2 = block2 - block2;
	Eigen::MatrixXd result3(nPoints, nPoints);
	result3 = block3 - block1;
	Eigen::MatrixXd tmp1(nPoints, nPoints);
	tmp1 = subtraction.block(0);
	Eigen::MatrixXd tmp2(nPoints, nPoints);
	tmp2 = subtraction.block(1);
	Eigen::MatrixXd tmp3(nPoints, nPoints);
	tmp3 = subtraction.block(2);
	for (int i = 0; i < nPoints; ++i)
	{
		for (int j = 0; j < nPoints; ++j)
		{
			EXPECT_EQ(result1(i, j), tmp1(i, j));
			EXPECT_EQ(result2(i, j), tmp2(i, j));
			EXPECT_EQ(result3(i, j), tmp3(i, j));
		}
	}
}
/*
TEST_F(BlockDiagonalMatrixTest, multiply_by_scalar)
{
	BlockDiagonalMatrix scaled1 = 2.5 * func1;
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

TEST_F(BlockDiagonalMatrixTest, multiply)
{
	double product = func1 * func2;
	double expected_product = values1.dot(values2);
	EXPECT_EQ(expected_product, product);
}*/
