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

#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

#include <algorithm>
#include <bitset>
#include <cmath>
#include <functional>
#include <iterator>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

/*! \fn inline int parity(std::bitset<nBits> bitrep)
 *  \param[in] bitrep a bitset
 *  \tparam nBits lenght of the input bitset
 *
 *  Calculate the parity of the bitset as defined by:
 *     bitrep[0] XOR bitrep[1] XOR ... XOR bitrep[nBits-1]
 */
template <size_t nBits>
inline int parity(std::bitset<nBits> bitrep)
{
    int parity = 0;
    for (size_t i = 0; i < bitrep.size(); ++i) {
        parity ^= bitrep[i];
    }
    return parity;
}

/*! \fn inline double parity(unsigned int i)
 *  \param[in] i an integer, usually an index for an irrep or a symmetry operation
 *
 * Returns parity of input integer.
 * The parity is defined as the result of using XOR on the bitrep
 * of the given integer. For example:
 *   2 -> 010 -> 0^1^0 = 1 -> -1.0
 *   6 -> 110 -> 1^1^0 = 0 ->  1.0
 *
 * It can also be interpreted as the action of a given operation on the
 * Cartesian axes:
 *      zyx         Parity
 *   0  000    E      1.0
 *   1  001   Oyz    -1.0
 *   2  010   Oxz    -1.0
 *   3  011   C2z     1.0
 *   4  100   Oxy    -1.0
 *   5  101   C2y     1.0
 *   6  110   C2x     1.0
 *   7  111    i     -1.0
 */
inline double parity(unsigned int i)
{
    // Use a ternary if construct. If the bitset is odd return -1.0 Return +1.0 otherwise.
    return (parity(std::bitset<3>(i)) ? -1.0 : 1.0);
}

/*! \fn inline bool isZero(double value, double threshold)
 *  \param[in] value     the value to be checked
 *  \param[in] threshold the threshold
 *
 *  Returns true if value is less or equal to threshold
 */
inline bool isZero(double value, double threshold)
{
    return (std::abs(value) <= threshold);
}

/*! \fn inline bool numericalZero(double value)
 *  \param[in] value the value to be checked
 *
 *  Returns true if value is less than 1.0e-14
 */
inline bool numericalZero(double value)
{
    return (isZero(value, 1.0e-14));
}

/*! \fn inline void symmetryBlocking(Eigen::MatrixXd & matrix, int cavitySize, int ntsirr, int nr_irrep)
 *  \param[out] matrix the matrix to be block-diagonalized
 *  \param[in]  cavitySize the size of the cavity (size of the matrix)
 *  \param[in]  ntsirr     the size of the irreducible portion of the cavity (size of the blocks)
 *  \param[in]  nr_irrep   the number of irreducible representations (number of blocks)
 */
inline void symmetryBlocking(Eigen::MatrixXd & matrix, int cavitySize, int ntsirr,
                             int nr_irrep)
{
    // This function implements the simmetry-blocking of the PCM
    // matrix due to point group symmetry as reported in:
    // L. Frediani, R. Cammi, C. S. Pomelli, J. Tomasi and K. Ruud, J. Comput.Chem. 25, 375 (2003)
    // u is the character table for the group (t in the paper)
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nr_irrep, nr_irrep);
    for (int i = 0; i < nr_irrep; ++i) {
        for (int j = 0; j < nr_irrep; ++j) {
            u(i, j) = parity(i&j);
        }
    }
    // Naming of indices:
    //     a, b, c, d   run over the total size of the cavity (N)
    //     i, j, k, l   run over the number of irreps (n)
    //     p, q, r, s   run over the irreducible size of the cavity (N/n)
    // Instead of forming U (T in the paper) and then perform the dense
    // multiplication, we multiply block-by-block using just the u matrix.
    //      matrix = U * matrix * Ut; U * Ut = Ut * U = id
    // First half-transformation, i.e. first_half = matrix * Ut
    Eigen::MatrixXd first_half = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    for (int i = 0; i < nr_irrep; ++i) {
        int ioff = i * ntsirr;
        for (int k = 0; k < nr_irrep; ++k) {
            int koff = k * ntsirr;
            for (int j = 0; j < nr_irrep; ++j) {
                int joff = j * ntsirr;
                double ujk = u(j, k) / nr_irrep;
                for (int p = 0; p < ntsirr; ++p) {
                    int a = ioff + p;
                    for (int q = 0; q < ntsirr; ++q) {
                        int b = joff + q;
                        int c = koff + q;
                        first_half(a, c) += matrix(a, b) * ujk;
                    }
                }
            }
        }
    }
    // Second half-transformation, i.e. matrix = U * first_half
    matrix.setZero(cavitySize, cavitySize);
    for (int i = 0; i < nr_irrep; ++i) {
        int ioff = i * ntsirr;
        for (int k = 0; k < nr_irrep; ++k) {
            int koff = k * ntsirr;
            for (int j = 0; j < nr_irrep; ++j) {
                int joff = j * ntsirr;
                double uij = u(i, j);
                for (int p = 0; p < ntsirr; ++p) {
                    int a = ioff + p;
                    int b = joff + p;
                    for (int q = 0; q < ntsirr; ++q) {
                        int c = koff + q;
                        matrix(a, c) += uij * first_half(b, c);
                    }
                }
            }
        }
    }
    // Traverse the matrix and discard numerical zeros
    for (int a = 0; a < cavitySize; ++a) {
        for (int b = 0; b < cavitySize; ++b) {
            if (numericalZero(matrix(a, b))) {
                matrix(a, b) = 0.0;
            }
        }
    }
}

/*! \fn inline void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix, const Eigen::MatrixXd & fullMatrix, int nrBlocks, int dimBlock)
 *  \param[out] blockedMatrix the result of packing fullMatrix
 *  \param[in]  fullMatrix the matrix to be packed
 *  \param[in]  dimBlock the dimension of the square blocks
 *  \param[in]  nrBlocks the number of square blocks
 */
inline void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix,
                            const Eigen::MatrixXd & fullMatrix, int dimBlock, int nrBlocks)
{
    // This function packs the square block diagonal fullMatrix with nrBlocks of dimension dimBlock
    // into a std::vector<Eigen::MatrixXd> containing nrBlocks square matrices of dimenion dimBlock.
    int j = 0;
    for (int i = 0; i < nrBlocks; ++i) {
        blockedMatrix.push_back(fullMatrix.block(j, j, dimBlock, dimBlock));
        j += dimBlock;
    }
}

/*! \fn inline void hermitivitize(Eigen::MatrixBase<Derived> & matrix_)
 *  \param[out] matrix_ the matrix to be hermitivitized
 *  \tparam     Derived the numeric type of matrix_ elements
 *
 *  Given matrix_ returns 0.5 * (matrix_ + matrix_^dagger)
 */
template <typename Derived>
inline void hermitivitize(Eigen::MatrixBase<Derived> & matrix_)
{
    // We need to use adjoint().eval() to avoid aliasing issues, see:
    // http://eigen.tuxfamily.org/dox/group__TopicAliasing.html
    // The adjoint is evaluated explicitly into an intermediate.
    matrix_ = 0.5 * (matrix_ + matrix_.adjoint().eval());
}

/*! \brief Returns an Eigen matrix of type T, with dimensions _rows*_columns.
 *  \param _rows the number of rows.
 *  \param _columns the number of columns.
 *  \param _inData the raw data buffer.
 *  \tparam T the data type of the matrix to be returned.
 *
 *  Warning! This function assumes that the raw buffer is in column-major order
 *  as in Fortran.
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getFromRawBuffer(
    size_t _rows, size_t _columns, void * _inData)
{
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _outData;
    _outData.resize(_rows, _columns);
    T * data = reinterpret_cast<T*>(_inData);
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > mapped_data(data,
            _rows, _columns);
    _outData = mapped_data;
    return _outData;
}

/*! \brief Return value of function defined on grid at an arbitrary point
 *  \param[in] point     where the function has to be evaluated
 *  \param[in] grid      holds points on grid where function is known
 *  \param[in] function  holds known function values
 *
 *  This function finds the nearest values for the given point
 *  and performs a linear interpolation.
 *  \warning This function assumes that grid has already been sorted!
 */
inline double linearInterpolation(const double point, const std::vector<double> & grid,
        const std::vector<double> & function)
{
    // Find nearest points on grid to the arbitrary point given
    size_t index = std::distance(grid.begin(), std::lower_bound(grid.begin(), grid.end(), point));

    // Parameters for the interpolating line
    double f_idx1 = function[index+1], f_idx = function[index];
    double g_idx1 = grid[index+1], g_idx = grid[index];
    double m = (f_idx1 - f_idx) / (g_idx1 - g_idx);

    return (m * (point - g_idx) + f_idx);
}

/*! \typedef DifferentiableFunction
 *  \brief sort of a function pointer to a function of a pair of vectors that can be numerically differentiated
 */
typedef std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> DifferentiableFunction;

/*! \brief Calculate directional derivative using a three-point stencil
 *  \param[in] func function to be differentiated
 *  \param[in] arg_1 first point, the directional derivative is calculated with respect to this point
 *  \param[in] arg_2 second point
 *  \param[in] direction the direction in which the directional derivative is to be evaluated
 *  \param[in] step finite difference value for the stencil
 */
inline double threePointStencil(const DifferentiableFunction & func,
        const Eigen::Vector3d & arg_1, const Eigen::Vector3d & arg_2,
        const Eigen::Vector3d & direction, double step)
{
    // f(x-h)
    Eigen::Vector3d delta_m1 = arg_1 - direction * step / direction.norm();
    // f(x+h)
    Eigen::Vector3d delta_1  = arg_1 + direction * step / direction.norm();

    Eigen::Vector2d stencil;
    stencil << -0.5, 0.5;
    Eigen::Vector2d function_values;
    function_values << func(delta_m1, arg_2), func(delta_1,  arg_2);

    return (function_values.dot(stencil) / step);
}

/*! \brief Calculate directional derivative using a five-point stencil
 *  \param[in] func function to be differentiated
 *  \param[in] arg_1 first point, the directional derivative is calculated with respect to this point
 *  \param[in] arg_2 second point
 *  \param[in] direction the direction in which the directional derivative is to be evaluated
 *  \param[in] step finite difference value for the stencil
 */
inline double fivePointStencil(const DifferentiableFunction & func,
        const Eigen::Vector3d & arg_1, const Eigen::Vector3d & arg_2,
        const Eigen::Vector3d & direction, double step)
{
    // f(x-2h)
    Eigen::Vector3d delta_m2 = arg_1 - 2.0 * direction * step / direction.norm();
    // f(x-h)
    Eigen::Vector3d delta_m1 = arg_1 - direction * step / direction.norm();
    // f(x+h)
    Eigen::Vector3d delta_1  = arg_1 + direction * step / direction.norm();
    // f(x+2h)
    Eigen::Vector3d delta_2  = arg_1 + 2.0 * direction * step / direction.norm();

    Eigen::Vector4d stencil;
    stencil << 1.0/12.0, -2.0/3.0, 2.0/3.0, -1.0/12.0;
    Eigen::Vector4d function_values;
    function_values << func(delta_m2, arg_2), func(delta_m1, arg_2),
                       func(delta_1,  arg_2), func(delta_2,  arg_2);

    return (function_values.dot(stencil) / step);
}

/*! \brief Calculate directional derivative using a seven-point stencil
 *  \param[in] func function to be differentiated
 *  \param[in] arg_1 first point, the directional derivative is calculated with respect to this point
 *  \param[in] arg_2 second point
 *  \param[in] direction the direction in which the directional derivative is to be evaluated
 *  \param[in] step finite difference value for the stencil
 */
inline double sevenPointStencil(const DifferentiableFunction & func,
        const Eigen::Vector3d & arg_1, const Eigen::Vector3d & arg_2,
        const Eigen::Vector3d & direction, double step)
{
    // f(x-3h)
    Eigen::Vector3d delta_m3 = arg_1 - 3.0 * direction * step / direction.norm();
    // f(x-2h)
    Eigen::Vector3d delta_m2 = arg_1 - 2.0 * direction * step / direction.norm();
    // f(x-h)
    Eigen::Vector3d delta_m1 = arg_1 - direction * step / direction.norm();
    // f(x+h)
    Eigen::Vector3d delta_1  = arg_1 + direction * step / direction.norm();
    // f(x+2h)
    Eigen::Vector3d delta_2  = arg_1 + 2.0 * direction * step / direction.norm();
    // f(x+3h)
    Eigen::Vector3d delta_3  = arg_1 + 3.0 * direction * step / direction.norm();

    Eigen::Matrix<double, 6, 1> stencil;
    stencil << -1.0/60.0, 3.0/20.0, -3.0/4.0, 3.0/4.0, -3.0/20.0, 1.0/60.0;
    Eigen::Matrix<double, 6, 1> function_values;
    function_values << func(delta_m3, arg_2), func(delta_m2, arg_2), func(delta_m1, arg_2),
                       func(delta_1,  arg_2), func(delta_2,  arg_2), func(delta_3, arg_2);

    return (function_values.dot(stencil) / step);
}

#endif // MATHUTILS_HPP
