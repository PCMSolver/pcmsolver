/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#pragma once

#include <algorithm>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "SplineFunction.hpp"
#include "Symmetry.hpp"
#include "cnpy.hpp"

/*! \file MathUtils.hpp */

namespace pcm {
namespace utils {
/*! \fn inline int parity(std::bitset<nBits> bitrep)
 *  \param[in] bitrep a bitset
 *  \tparam nBits lenght of the input bitset
 *
 *  Calculate the parity of the bitset as defined by:
 *     bitrep[0] XOR bitrep[1] XOR ... XOR bitrep[nBits-1]
 */
template <size_t nBits> inline int parity(std::bitset<nBits> bitrep) {
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
inline double parity(unsigned int i) {
  // Use a ternary if construct. If the bitset is odd return -1.0 Return +1.0
  // otherwise.
  return (parity(std::bitset<3>(i)) ? -1.0 : 1.0);
}

/*! \fn inline bool isZero(double value, double threshold)
 *  \param[in] value     the value to be checked
 *  \param[in] threshold the threshold
 *
 *  Returns true if value is less or equal to threshold
 */
inline bool isZero(double value, double threshold) {
  return (std::abs(value) <= threshold);
}

/*! \fn inline bool numericalZero(double value)
 *  \param[in] value the value to be checked
 *
 *  Returns true if value is less than 1.0e-14
 */
inline bool numericalZero(double value) { return (isZero(value, 1.0e-14)); }

/*! \fn inline void symmetryBlocking(Eigen::MatrixXd & matrix, int cavitySize, int
 * ntsirr, int nr_irrep)
 *  \param[out] matrix the matrix to be block-diagonalized
 *  \param[in]  cavitySize the size of the cavity (size of the matrix)
 *  \param[in]  ntsirr     the size of the irreducible portion of the cavity (size of
 * the blocks)
 *  \param[in]  nr_irrep   the number of irreducible representations (number of
 * blocks)
 */
inline void symmetryBlocking(Eigen::MatrixXd & matrix,
                             PCMSolverIndex cavitySize,
                             PCMSolverIndex ntsirr,
                             int nr_irrep) {
  // This function implements the simmetry-blocking of the PCM
  // matrix due to point group symmetry as reported in:
  // L. Frediani, R. Cammi, C. S. Pomelli, J. Tomasi and K. Ruud, J. Comput.Chem. 25,
  // 375 (2003)
  // u is the character table for the group (t in the paper)
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nr_irrep, nr_irrep);
  for (int i = 0; i < nr_irrep; ++i) {
    for (int j = 0; j < nr_irrep; ++j) {
      u(i, j) = parity(i & j);
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
  for (PCMSolverIndex a = 0; a < cavitySize; ++a) {
    for (PCMSolverIndex b = 0; b < cavitySize; ++b) {
      if (numericalZero(matrix(a, b))) {
        matrix(a, b) = 0.0;
      }
    }
  }
}

/*! \fn inline void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix,
 * const Eigen::MatrixXd & fullMatrix, int nrBlocks, int dimBlock)
 *  \param[out] blockedMatrix the result of packing fullMatrix
 *  \param[in]  fullMatrix the matrix to be packed
 *  \param[in]  dimBlock the dimension of the square blocks
 *  \param[in]  nrBlocks the number of square blocks
 */
inline void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix,
                            const Eigen::MatrixXd & fullMatrix,
                            int dimBlock,
                            int nrBlocks) {
  // This function packs the square block diagonal fullMatrix with nrBlocks of
  // dimension dimBlock
  // into a std::vector<Eigen::MatrixXd> containing nrBlocks square matrices of
  // dimenion dimBlock.
  int j = 0;
  for (int i = 0; i < nrBlocks; ++i) {
    blockedMatrix.push_back(fullMatrix.block(j, j, dimBlock, dimBlock));
    j += dimBlock;
  }
}

/*! \fn inline void hermitivitize(Eigen::MatrixBase<Derived> & obj_)
 *  \param[out] obj_ the Eigen object to be hermitivitized
 *  \tparam     Derived the numeric type of obj_ elements
 *
 *  Given obj_ returns 0.5 * (obj_ + obj_^dagger)
 *  \note We check if a matrix or vector was given, since in the latter
 *  case we only want the complex conjugation operation to happen.
 */
template <typename Derived>
inline void hermitivitize(Eigen::MatrixBase<Derived> & obj_) {
  // We need to use adjoint().eval() to avoid aliasing issues, see:
  // http://eigen.tuxfamily.org/dox/group__TopicAliasing.html
  // The adjoint is evaluated explicitly into an intermediate.
  obj_ = 0.5 * (obj_ + obj_.adjoint().eval());
}

/*! \fn inline void eulerRotation(Eigen::Matrix3d & R_, const Eigen::Vector3d &
 *eulerAngles_)
 *  \brief Build rotation matrix between two reference frames given the Euler angles.
 *  \param[out] R_ the rotation matrix
 *  \param[in]  eulerAngles_ the Euler angles, in degrees, describing the rotation
 *
 *  We assume the convention \f$ R = Z_3 X_2 Z_1 \f$ for the ordering of the
 *extrinsic
 *  elemental rotations (see http://en.wikipedia.org/wiki/Euler_angles)
 *  The Euler angles are given in the order \f$ \phi, \theta, \psi \f$.
 *  If we write \f$ c_i, s_i \,\, i = 1, 3 \f$ for their cosines and sines the
 *rotation
 *  matrix will be:
 *  \f[
 *  	R = \begin{pmatrix}
 *  	      c_1c_3 - s_1c_2s_3 & -s_1c_3 - c_1c_2s_3 &  s_2s_3 \\
 *  	      c_1s_3 + s_1c_2c_3 & -s_1s_3 + c_1c_2c_3 & -s_2c_3 \\
 *  	      s_1s_2             & c_1s_2              &  c_2
 *  	    \end{pmatrix}
 *  \f]
 *  Eigen's geometry module is used to calculate the rotation matrix
 */
inline void eulerRotation(Eigen::Matrix3d & R_,
                          const Eigen::Vector3d & eulerAngles_) {
  double to_radians = M_PI / 180.0;
  double phi = eulerAngles_(0) * to_radians;
  double theta = eulerAngles_(1) * to_radians;
  double psi = eulerAngles_(2) * to_radians;
  R_ = Eigen::AngleAxis<double>(psi, Eigen::Vector3d::UnitZ()) *
       Eigen::AngleAxis<double>(theta, Eigen::Vector3d::UnitX()) *
       Eigen::AngleAxis<double>(phi, Eigen::Vector3d::UnitZ());
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
inline double linearInterpolation(const double point,
                                  const std::vector<double> & grid,
                                  const std::vector<double> & function) {
  // Find nearest points on grid to the arbitrary point given
  size_t index = std::distance(grid.begin(),
                               std::lower_bound(grid.begin(), grid.end(), point)) -
                 1;

  // Parameters for the interpolating line
  double y_1 = function[index], y_0 = function[index - 1];
  double x_1 = grid[index], x_0 = grid[index - 1];
  double m = (y_1 - y_0) / (x_1 - x_0);

  return (m * (point - x_0) + y_0);
}

/*! \brief Return value of function defined on grid at an arbitrary point
 *  \param[in] point     where the function has to be evaluated
 *  \param[in] grid      holds points on grid where function is known
 *  \param[in] function  holds known function values
 *
 *  This function finds the nearest values for the given point
 *  and performs a cubic spline interpolation.
 *  \warning This function assumes that grid has already been sorted!
 */
inline double splineInterpolation(const double point,
                                  const std::vector<double> & grid,
                                  const std::vector<double> & function) {
  // Find nearest points on grid to the arbitrary point given
  size_t index = std::distance(grid.begin(),
                               std::lower_bound(grid.begin(), grid.end(), point)) -
                 1;

  // Parameters for the interpolating spline
  Eigen::Vector3d x, y;
  x << grid[index - 1], grid[index], grid[index + 1];
  y << function[index - 1], function[index], function[index + 1];
  SplineFunction s(x, y);

  return s(point);
}

/*! \brief Prints Eigen object (matrix or vector) to file
 *  \param[in] matrix Eigen object
 *  \param[in] fname  name of the file
 *  \tparam Derived template parameters of the MatrixBase object
 *
 *  \note This is for debugging only, the format is in fact rather ugly.
 *      Row index     Column index      Matrix entry
 *         0               0              0.0000
 */
template <typename Derived>
inline void print_eigen_matrix(const Eigen::MatrixBase<Derived> & matrix,
                               const std::string & fname) {
  std::ofstream fout;
  fout.open(fname.c_str());
  fout << " Row index " << '\t' << " Column index " << '\t' << " Matrix entry "
       << std::endl;
  int rows = matrix.rows();
  int cols = matrix.cols();
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      fout << i << '\t' << j << '\t' << matrix(i, j) << std::endl;
    }
  }
  fout.close();
}
} // namespace utils
} // namespace pcm

namespace cnpy {
/*! Custom overloads for cnpy load and save functions */
namespace custom {
/*! \brief Save Eigen object to NumPy array file
 *  \param fname name of the NumPy array file
 *  \param obj Eigen object to be saved, either a matrix or a vector
 *  \tparam Scalar the data type of the matrix to be returned. Default is double
 *  \tparam Rows number of rows in the Eigen object. Default is dynamic
 e  \tparam Cols number of columns in the Eigen object. Default is dynamic
 */
template <typename Scalar, int Rows, int Cols>
inline void npy_save(const std::string & fname,
                     const Eigen::Matrix<Scalar, Rows, Cols> & obj) {
  unsigned int rows = static_cast<unsigned int>(obj.rows());
  unsigned int cols = static_cast<unsigned int>(obj.cols());
  const unsigned int shape[] = {rows, cols};
  bool colMajor = obj.IsRowMajor ? false : true;
  cnpy::npy_save(fname, obj.data(), shape, 2, "w", colMajor);
}

/*! \brief Save Eigen object to a compressed NumPy file
 *  \param fname name of the compressed NumPy file
 *  \param name tag for the given object in the compressed NumPy file
 *  \param obj  Eigen object to be saved, either a matrix or a vector
 *  \param overwrite if file exists, overwrite. Appends by default.
 *  \tparam Scalar the data type of the matrix to be returned. Default is double
 *  \tparam Rows number of rows in the Eigen object. Default is dynamic
 *  \tparam Cols number of columns in the Eigen object. Default is dynamic
 */
template <typename Scalar, int Rows, int Cols>
inline void npz_save(const std::string & fname,
                     const std::string & name,
                     const Eigen::Matrix<Scalar, Rows, Cols> & obj,
                     bool overwrite = false) {
  unsigned int rows = static_cast<unsigned int>(obj.rows());
  unsigned int cols = static_cast<unsigned int>(obj.cols());
  const unsigned int shape[] = {rows, cols};
  bool colMajor = obj.IsRowMajor ? false : true;
  std::string mode = overwrite ? "w" : "a";
  cnpy::npz_save(fname, name, obj.data(), shape, 2, mode, colMajor);
}

/*! \brief Load NpyArray object into Eigen object
 *  \param npy_array the NpyArray object
 *  \tparam Scalar the data type of the matrix to be returned. Default is double
 *  \return An Eigen object (matrix or vector) with the data
 *
 *  \todo Extend to read in also data in row-major (C) storage order
 *  \warning We check that the rank of the object read is not more than 2
 *  Eigen cannot handle general tensors.
 */
template <typename Scalar>
inline Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> npy_to_eigen(
    const NpyArray & npy_array) {
  if (npy_array.shape.size() > 2)
    PCMSOLVER_ERROR("Only vectors and matrices can be read into Eigen objects.");
  return Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> >(
      reinterpret_cast<Scalar *>(npy_array.data),
      npy_array.shape[0],
      npy_array.shape[1]);
}

/*! \brief Load NumPy array file into Eigen object
 *  \param fname name of the NumPy array file
 *  \tparam Scalar the data type of the matrix to be returned. Default is double
 *  \return An Eigen object (matrix or vector) with the data
 *
 *  \todo Extend to read in also data in row-major (C) storage order
 */
template <typename Scalar>
inline Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> npy_load(
    const std::string & fname) {
  return npy_to_eigen<Scalar>(cnpy::npy_load(fname));
}
} // namespace custom
} // namespace cnpy
