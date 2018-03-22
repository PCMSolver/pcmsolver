/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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

#include <cmath>
#include <fstream>
#include <iomanip>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "cnpy.hpp"

/*! \file MathUtils.hpp */

namespace pcm {
namespace utils {
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

/*! \fn inline int sign(T val)
 *  \param[in] val value whose sign should be determined
 *  \tparam    T of the parameter val
 *
 *  This function implements the signum function and returns the sign of the passed
 * value: -1, 0 or 1
 */
template <typename T> inline int sign(T val) { return (T(0) < val) - (val < T(0)); }

/*! \brief Enforce Hermitian symmetry on a matrix
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
 *    R = \begin{pmatrix}
 *          c_1c_3 - s_1c_2s_3 & -s_1c_3 - c_1c_2s_3 &  s_2s_3 \\
 *          c_1s_3 + s_1c_2c_3 & -s_1s_3 + c_1c_2c_3 & -s_2c_3 \\
 *          s_1s_2             & c_1s_2              &  c_2
 *        \end{pmatrix}
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
