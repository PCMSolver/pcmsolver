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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>
#include <boost/function.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/map.hpp>

#include "Element.hpp"
#include "QuadratureRules.hpp"
#include "Sphere.hpp"
#include "Symmetry.hpp"

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
            u(i, j) = Symmetry::parity(i&j);
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

/*! \fn inline void eulerRotation(Eigen::Matrix3d & R_, const Eigen::Vector3d & eulerAngles_)
 *  \brief Build rotation matrix between two reference frames given the Euler angles.
 *  \param[out] R_ the rotation matrix
 *  \param[in]  eulerAngles_ the Euler angles, in degrees, describing the rotation
 *
 *  We assume the convention \f$ R = Z_3 X_2 Z_1 \f$ for the ordering of the extrinsic
 *  elemental rotations (see http://en.wikipedia.org/wiki/Euler_angles)
 *  The Euler angles are given in the order \f$ \phi, \theta, \psi \f$.
 *  If we write \f$ c_i, s_i \,\, i = 1, 3 \f$ for their cosines and sines the rotation
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
inline void eulerRotation(Eigen::Matrix3d & R_, const Eigen::Vector3d & eulerAngles_)
{
    double to_radians = M_PI / 180.0;
    double phi   = eulerAngles_(0) * to_radians;
    double theta = eulerAngles_(1) * to_radians;
    double psi   = eulerAngles_(2) * to_radians;
    R_ = Eigen::AngleAxis<double>(psi,   Eigen::Vector3d::UnitZ())
         * Eigen::AngleAxis<double>(theta, Eigen::Vector3d::UnitX())
         * Eigen::AngleAxis<double>(phi,   Eigen::Vector3d::UnitZ());
}

typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
singleLayerIntegrand;

template <int PhiPoints_, int ThetaPoints_>
inline double integrator(const singleLayerIntegrand & F, const Element & e)
{
    double result = 0.0;

    // Get the quadrature rules for azimuthal and polar integrations
    namespace mpl = boost::mpl;
    typedef typename mpl::at<rules_map, mpl::int_<PhiPoints_> >::type PhiPolicy;
    typedef typename mpl::at<rules_map, mpl::int_<ThetaPoints_> >::type ThetaPolicy;
    QuadratureRule<PhiPolicy> phiRule;
    QuadratureRule<ThetaPolicy> thetaRule;
    int upper_phi = PhiPoints_ / 2; // Upper limit for loop on phi points
    int upper_theta = ThetaPoints_ / 2; // Upper limit for loop on theta points

    // Extract relevant data from Element
    int nVertices = e.nVertices();
    Eigen::Vector3d normal = e.normal();
    Sphere sph = e.sphere();
    Eigen::Matrix3Xd vertices = e.vertices();
    Eigen::Matrix3Xd arcs = e.arcs();

    // Calculation of the tangent and the bitangent (binormal) vectors
    // Tangent, Bitangent and Normal form a local reference frame:
    // T <-> x; B <-> y; N <-> z
    Eigen::Vector3d tangent, bitangent;
    tangent_and_bitangent(normal, tangent, bitangent);

    std::vector<double> theta(nVertices), phi(nVertices), phinumb(nVertices+1);
    std::vector<int> numb(nVertices+1);
    // Clean-up heap crap
    std::fill_n(theta.begin(),   nVertices,   0.0);
    std::fill_n(phi.begin(),     nVertices,   0.0);
    std::fill_n(numb.begin(),    nVertices+1, 0);
    std::fill_n(phinumb.begin(), nVertices+1, 0.0);
    // Populate arrays and redefine tangent and bitangent
    e.spherical_polygon(tangent, bitangent, theta, phi, phinumb, numb);

    // Actual integration occurs here
    for (int i = 0; i < nVertices; ++i) { // Loop on edges
        double phiLower = phinumb[i];
        double phiUpper = phinumb[i+1];
        double phiA= (phiUpper - phiLower) / 2.0;
        double phiB = (phiUpper + phiLower) / 2.0;
        double thetaLower = theta[numb[i]];
        double thetaUpper = theta[numb[i+1]];
        double thetaMax = 0.0;
        Eigen::Vector3d oc = (arcs.col(i) - sph.center()) / sph.radius();
        double oc_norm = oc.norm();
        double oc_norm2 = std::pow(oc_norm, 2);
        for (int j = 0; j < upper_phi; ++j) { // Loop on Gaussian points: phi integration
            for (int k = 0; k <= 1; ++k) {
                double ph = (2*k - 1) * phiA * phiRule.abscissa(j) + phiB;
                double cos_ph = std::cos(ph);
                double sin_ph = std::sin(ph);
		// We need to calculate the upper bound for the integration on theta, which depends on phi
                if (oc_norm2 < 1.0e-07) { // This means that the edge is centered on the same sphere the tessera belongs to
                    double cotg_thmax = (std::sin(ph-phiLower) / std::tan(thetaUpper) + std::sin(phiUpper-ph) / std::tan(
                                             thetaLower)) / std::sin(phiUpper - phiLower);
                    thetaMax = std::atan(1.0 / cotg_thmax);
                } else {
                    Eigen::Vector3d scratch;
                    scratch << tangent.dot(oc), bitangent.dot(oc), normal.dot(oc);
                    double aa = std::pow(tangent.dot(oc)*cos_ph + bitangent.dot(oc)*sin_ph,
                                         2) + std::pow(normal.dot(oc), 2);
                    double bb = -normal.dot(oc) * oc_norm2;
                    double cc = std::pow(oc_norm2,
                                         2) - std::pow(tangent.dot(oc)*cos_ph + bitangent.dot(oc)*sin_ph, 2);
                    double ds = std::pow(bb, 2) - aa*cc;
                    if (ds < 0.0) ds = 0.0;
                    double cs = (-bb + std::sqrt(ds)) / aa;
                    if (cs > 1.0) cs = 1.0;
                    if (cs < -1.0) cs = 1.0;
                    thetaMax = std::acos(cs);
                }
                double scratch = 0.0;
                if (!(thetaMax < 1.0e-08)) {
                double thetaA = thetaMax / 2.0;
                    for (int l = 0; l < upper_theta; ++l) { // Loop on Gaussian points: theta integration
                        for (int m = 0; m <= 1; ++m) {
                            double th = (2*m - 1) * thetaA * thetaRule.abscissa(l) + thetaA;
                            double cos_th = std::cos(th);
                            double sin_th = std::sin(th);
                            Eigen::Vector3d point;
                            point(0) = tangent(0) * sin_th * cos_ph
                                       + bitangent(0) * sin_th * sin_ph
                                       + normal(0) * (cos_th - 1.0);
                            point(1) = tangent(1) * sin_th * cos_ph
                                       + bitangent(1) * sin_th * sin_ph
                                       + normal(1) * (cos_th - 1.0);
                            point(2) = tangent(2) * sin_th * cos_ph
                                       + bitangent(2) * sin_th * sin_ph
                                       + normal(2) * (cos_th - 1.0);
                            double value = F(point,
                                             Eigen::Vector3d::Zero()); // Evaluate integrand at Gaussian point
                            scratch += std::pow(sph.radius(), 2) * value * sin_th * thetaA * thetaRule.weight(l);
                        }
                    }
                    result += scratch * phiA * phiRule.weight(j);
                }
            }
        }
    }
    return result;
}

typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &, 
		        const Eigen::Vector3d &)> doubleLayerIntegrand;

template <int PhiPoints_, int ThetaPoints_>
inline double integrator(const doubleLayerIntegrand & F, const Element & e)
{
    double result = 0.0;

    // Get the quadrature rules for azimuthal and polar integrations
    namespace mpl = boost::mpl;
    typedef typename mpl::at<rules_map, mpl::int_<PhiPoints_> >::type PhiPolicy;
    typedef typename mpl::at<rules_map, mpl::int_<ThetaPoints_> >::type ThetaPolicy;
    QuadratureRule<PhiPolicy> phiRule;
    QuadratureRule<ThetaPolicy> thetaRule;
    int upper_phi = PhiPoints_ / 2; // Upper limit for loop on phi points
    int upper_theta = ThetaPoints_ / 2; // Upper limit for loop on theta points

    // Extract relevant data from Element
    int nVertices = e.nVertices();
    Eigen::Vector3d normal = e.normal();
    Sphere sph = e.sphere();
    Eigen::Matrix3Xd vertices = e.vertices();
    Eigen::Matrix3Xd arcs = e.arcs();

    // Calculation of the tangent and the bitangent (binormal) vectors
    // Tangent, Bitangent and Normal form a local reference frame:
    // T <-> x; B <-> y; N <-> z
    Eigen::Vector3d tangent, bitangent;
    tangent_and_bitangent(normal, tangent, bitangent);

    std::vector<double> theta(nVertices), phi(nVertices), phinumb(nVertices+1);
    std::vector<int> numb(nVertices+1);
    // Clean-up heap crap
    std::fill_n(theta.begin(),   nVertices,   0.0);
    std::fill_n(phi.begin(),     nVertices,   0.0);
    std::fill_n(numb.begin(),    nVertices+1, 0);
    std::fill_n(phinumb.begin(), nVertices+1, 0.0);
    // Populate arrays and redefine tangent and bitangent
    e.spherical_polygon(tangent, bitangent, theta, phi, phinumb, numb);

    // Actual integration occurs here
    for (int i = 0; i < nVertices; ++i) { // Loop on edges
        double phiLower = phinumb[i]; // Lower vertex of edge
        double phiUpper = phinumb[i+1]; // Upper vertex of edge
        double phiA = (phiUpper - phiLower) / 2.0;
        double phiB = (phiUpper + phiLower) / 2.0;
        double thetaLower = theta[numb[i]];
        double thetaUpper = theta[numb[i+1]];
        double thetaMax = 0.0;
        Eigen::Vector3d oc = (arcs.col(i) - sph.center()) / sph.radius();
        double oc_norm = oc.norm();
        double oc_norm2 = std::pow(oc_norm, 2);
        for (int j = 0; j < upper_phi; ++j) { // Loop on Gaussian points: phi integration
            for (int k = 0; k <= 1; ++k) {
                double ph = (2*k - 1) * phiA * phiRule.abscissa(j) + phiB;
                double cos_phi = std::cos(ph);
                double sin_phi = std::sin(ph);
                if (oc_norm2 < 1.0e-07) { // This should check if oc_norm2 is zero
                    double cotg_thmax = (std::sin(ph-phiLower) / std::tan(thetaUpper) + std::sin(phiUpper-ph) / std::tan(
                                             thetaLower)) / std::sin(phiUpper - phiLower);
                    thetaMax = std::atan(1.0 / cotg_thmax);
                } else {
                    Eigen::Vector3d scratch;
                    scratch << tangent.dot(oc), bitangent.dot(oc), normal.dot(oc);
                    double aa = std::pow(tangent.dot(oc)*cos_phi + bitangent.dot(oc)*sin_phi,
                                         2) + std::pow(normal.dot(oc), 2);
                    double bb = -normal.dot(oc) * oc_norm2;
                    double cc = std::pow(oc_norm2,
                                         2) - std::pow(tangent.dot(oc)*cos_phi + bitangent.dot(oc)*sin_phi, 2);
                    double ds = std::pow(bb, 2) - aa*cc;
                    if (ds < 0.0) ds = 0.0;
                    double cs = (-bb + std::sqrt(ds)) / aa;
                    if (cs > 1.0) cs = 1.0;
                    if (cs < -1.0) cs = 1.0;
                    thetaMax = std::acos(cs);
                }
                double thetaA = thetaMax / 2.0;
                double scratch = 0.0;
                if (!(thetaMax < 1.0e-08)) {
                    for (int l = 0; l < upper_theta; ++l) { // Loop on Gaussian points: theta integration
                        for (int m = 0; m <= 1; ++m) {
                            double th = (2*m - 1) * thetaA  * thetaRule.abscissa(l) + thetaA;
                            double cos_theta = std::cos(th);
                            double sin_theta = std::sin(th);
                            Eigen::Vector3d point;
                            point(0) = tangent(0) * sin_theta * cos_phi
                                       + bitangent(0) * sin_theta * sin_phi
                                       + normal(0) * (cos_theta - 1.0);
                            point(1) = tangent(1) * sin_theta * cos_phi
                                       + bitangent(1) * sin_theta * sin_phi
                                       + normal(1) * (cos_theta - 1.0);
                            point(2) = tangent(2) * sin_theta * cos_phi
                                       + bitangent(2) * sin_theta * sin_phi
                                       + normal(2) * (cos_theta - 1.0);
                            double value = F(e.normal(),
                                             Eigen::Vector3d::Zero(),
					     point); // Evaluate integrand at Gaussian point
                            scratch += std::pow(sph.radius(), 2) * value * sin_theta * thetaA * thetaRule.weight(l);
                        }
                    }
                    result += scratch * phiA * phiRule.weight(j);
                }
            }
        }
    }
    return result;
}

#endif // MATHUTILS_HPP
