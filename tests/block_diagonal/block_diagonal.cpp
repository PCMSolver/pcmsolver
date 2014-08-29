/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *     
 *     This file is part of PCMSolver.
 *     
 *     PCMSolver is free software: you can redistribute it and/or modify       
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *     
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *     
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *     
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#define BOOST_TEST_MODULE BlockDiagonalMatrix

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include "BlockDiagonalMatrix.hpp"

struct BlockDiagonalMatrixTest {
    int nBlocks;
    int nPoints;
    int fullDim;
    Eigen::MatrixXd block1;
    Eigen::MatrixXd block2;
    Eigen::MatrixXd block3;
    BlockDiagonalMatrixTest() { SetUp(); }
    BlockDiagonalMatrix<double, 3, 15> bd1;
    BlockDiagonalMatrix<double, 3, 15> bd2;
    void SetUp() {
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

/*! \class BlockDiagonalMatrix
 *  \test \b BlockDiagonalMatrixTest_addition tests addition of two compatible block diagonal matrices
 */
BOOST_FIXTURE_TEST_CASE(addition, BlockDiagonalMatrixTest)
{
    const int nBlocks = 3;
    const int nPoints = 15;
    BlockDiagonalMatrix<double, nBlocks, nPoints> addition = bd1 + bd2;
    BOOST_REQUIRE_EQUAL(fullDim, addition.fullDim());
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
    for (int i = 0; i < nPoints; ++i) {
        for (int j = 0; j < nPoints; ++j) {
            BOOST_REQUIRE_CLOSE(result1(i, j), tmp1(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result2(i, j), tmp2(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result3(i, j), tmp3(i, j), 1.0e-12);
        }
    }
}

/*! \class BlockDiagonalMatrix
 *  \test \b BlockDiagonalMatrixTest_subtraction tests subtraction of two compatible block diagonal matrices
 */
BOOST_FIXTURE_TEST_CASE(subtraction, BlockDiagonalMatrixTest)
{
    const int nBlocks = 3;
    const int nPoints = 15;
    BlockDiagonalMatrix<double, nBlocks, nPoints> subtraction = bd1 - bd2;
    BOOST_REQUIRE_EQUAL(fullDim, subtraction.fullDim());
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
    for (int i = 0; i < nPoints; ++i) {
        for (int j = 0; j < nPoints; ++j) {
            BOOST_REQUIRE_CLOSE(result1(i, j), tmp1(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result2(i, j), tmp2(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result3(i, j), tmp3(i, j), 1.0e-12);
        }
    }
}

/*! \class BlockDiagonalMatrix
 *  \test \b BlockDiagonalMatrixTest_multiply_by_scalar1 tests multiplication of a block diagonal matrix by a scalar
 */
BOOST_FIXTURE_TEST_CASE(multiply_by_scalar1, BlockDiagonalMatrixTest)
{
    const int nBlocks = 3;
    const int nPoints = 15;
    BlockDiagonalMatrix<double, nBlocks, nPoints> scaled1 = 2.5 * bd1;
    BOOST_REQUIRE_EQUAL(fullDim, scaled1.fullDim());
    Eigen::MatrixXd result1(nPoints, nPoints);
    result1 = 2.5 * block1;
    Eigen::MatrixXd result2(nPoints, nPoints);
    result2 = 2.5 * block2;
    Eigen::MatrixXd result3(nPoints, nPoints);
    result3 = 2.5 * block3;
    Eigen::MatrixXd tmp1(nPoints, nPoints);
    tmp1 = scaled1.block(0);
    Eigen::MatrixXd tmp2(nPoints, nPoints);
    tmp2 = scaled1.block(1);
    Eigen::MatrixXd tmp3(nPoints, nPoints);
    tmp3 = scaled1.block(2);
    for (int i = 0; i < nPoints; ++i) {
        for (int j = 0; j < nPoints; ++j) {
            BOOST_REQUIRE_CLOSE(result1(i, j), tmp1(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result2(i, j), tmp2(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result3(i, j), tmp3(i, j), 1.0e-12);
        }
    }
}

/*! \class BlockDiagonalMatrix
 *  \test \b BlockDiagonalMatrixTest_multiply_by_scalar2 tests multiplication of a block diagonal matrix by a scalar
 */
BOOST_FIXTURE_TEST_CASE(multiply_by_scalar2, BlockDiagonalMatrixTest)
{
    int nPoints = 15;
    bd2 *= 0.5;
    BOOST_REQUIRE_EQUAL(fullDim, bd2.fullDim());
    Eigen::MatrixXd result1(nPoints, nPoints);
    result1 = 0.5 * block3;
    Eigen::MatrixXd result2(nPoints, nPoints);
    result2 = 0.5 * block2;
    Eigen::MatrixXd result3(nPoints, nPoints);
    result3 = 0.5 * block1;
    Eigen::MatrixXd tmp1(nPoints, nPoints);
    tmp1 = bd2.block(0);
    Eigen::MatrixXd tmp2(nPoints, nPoints);
    tmp2 = bd2.block(1);
    Eigen::MatrixXd tmp3(nPoints, nPoints);
    tmp3 = bd2.block(2);
    for (int i = 0; i < nPoints; ++i) {
        for (int j = 0; j < nPoints; ++j) {
            BOOST_REQUIRE_CLOSE(result1(i, j), tmp1(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result2(i, j), tmp2(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result3(i, j), tmp3(i, j), 1.0e-12);
        }
    }
}

/*! \class BlockDiagonalMatrix
 *  \test \b BlockDiagonalMatrixTest_multiply1 tests matrix multiplication of two compatible block diagonal matrices
 */
BOOST_FIXTURE_TEST_CASE(multiply1, BlockDiagonalMatrixTest)
{
    const int nBlocks = 3;
    const int nPoints = 15;
    BlockDiagonalMatrix<double, nBlocks, nPoints> prod = bd1 * bd2;
    Eigen::MatrixXd result1(nPoints, nPoints);
    result1 = block1 * block3;
    Eigen::MatrixXd result2(nPoints, nPoints);
    result2 = block2 * block2;
    Eigen::MatrixXd result3(nPoints, nPoints);
    result3 = block3 * block1;
    Eigen::MatrixXd tmp1(nPoints, nPoints);
    tmp1 = prod.block(0);
    Eigen::MatrixXd tmp2(nPoints, nPoints);
    tmp2 = prod.block(1);
    Eigen::MatrixXd tmp3(nPoints, nPoints);
    tmp3 = prod.block(2);
    for (int i = 0; i < nPoints; ++i) {
        for (int j = 0; j < nPoints; ++j) {
            BOOST_REQUIRE_CLOSE(result1(i, j), tmp1(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result2(i, j), tmp2(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result3(i, j), tmp3(i, j), 1.0e-12);
        }
    }
}

/*! \class BlockDiagonalMatrix
 *  \test \b BlockDiagonalMatrixTest_multiply2 tests matrix multiplication of two compatible block diagonal matrices
 */
BOOST_FIXTURE_TEST_CASE(multiply2, BlockDiagonalMatrixTest)
{
    bd1 *= bd2;
    Eigen::MatrixXd result1(nPoints, nPoints);
    result1 = block1 * block3;
    Eigen::MatrixXd result2(nPoints, nPoints);
    result2 = block2 * block2;
    Eigen::MatrixXd result3(nPoints, nPoints);
    result3 = block3 * block1;
    Eigen::MatrixXd tmp1(nPoints, nPoints);
    tmp1 = bd1.block(0);
    Eigen::MatrixXd tmp2(nPoints, nPoints);
    tmp2 = bd1.block(1);
    Eigen::MatrixXd tmp3(nPoints, nPoints);
    tmp3 = bd1.block(2);
    for (int i = 0; i < nPoints; ++i) {
        for (int j = 0; j < nPoints; ++j) {
            BOOST_REQUIRE_CLOSE(result1(i, j), tmp1(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result2(i, j), tmp2(i, j), 1.0e-12);
            BOOST_REQUIRE_CLOSE(result3(i, j), tmp3(i, j), 1.0e-12);
        }
    }
}
