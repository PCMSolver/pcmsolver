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

#include "Numerical.hpp"

#include <cmath>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include <boost/mpl/at.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/map.hpp>

#include "BIOperatorData.hpp"
#include "cavity/Element.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/QuadratureRules.hpp"

namespace pcm {
using cavity::Element;
using cavity::detail::tangent_and_bitangent;

namespace bi_operators {
Eigen::MatrixXd Numerical::computeS_impl(const std::vector<Element> & elems,
                                         const IGreensFunction & gf) const {
  PCMSolverIndex cavitySize = elems.size();
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    Element source = elems[i];
    S(i, i) = integrateS<32, 16>(gf.exportKernelS(), source);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      Element probe = elems[j];
      if (i != j)
        S(i, j) = gf.kernelS(source.center(), probe.center());
    }
  }
  return S;
}

Eigen::MatrixXd Numerical::computeD_impl(const std::vector<Element> & elems,
                                         const IGreensFunction & gf) const {
  PCMSolverIndex cavitySize = elems.size();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    Element source = elems[i];
    D(i, i) = integrateD<32, 16>(gf.exportKernelD(), source);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      Element probe = elems[j];
      if (i != j)
        D(i, j) =
            gf.kernelD(probe.normal().normalized(), source.center(), probe.center());
    }
  }
  return D;
}

IBoundaryIntegralOperator * createNumerical(const BIOperatorData & /* data */) {
  return new Numerical();
}

template <int PhiPoints, int ThetaPoints>
double integrateS(const KernelS & F, const Element & e) {
  double result = 0.0;

  // Get the quadrature rules for azimuthal and polar integrations
  namespace mpl = boost::mpl;
  typedef typename mpl::at<rules_map, mpl::int_<PhiPoints> >::type PhiPolicy;
  typedef typename mpl::at<rules_map, mpl::int_<ThetaPoints> >::type ThetaPolicy;
  QuadratureRule<PhiPolicy> phiRule;
  QuadratureRule<ThetaPolicy> thetaRule;
  int upper_phi = PhiPoints / 2;     // Upper limit for loop on phi points
  int upper_theta = ThetaPoints / 2; // Upper limit for loop on theta points

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

  std::vector<double> theta(nVertices), phi(nVertices), phinumb(nVertices + 1);
  std::vector<int> numb(nVertices + 1);
  // Clean-up heap crap
  std::fill_n(theta.begin(), nVertices, 0.0);
  std::fill_n(phi.begin(), nVertices, 0.0);
  std::fill_n(numb.begin(), nVertices + 1, 0);
  std::fill_n(phinumb.begin(), nVertices + 1, 0.0);
  // Populate arrays and redefine tangent and bitangent
  e.spherical_polygon(tangent, bitangent, theta, phi, phinumb, numb);

  // Actual integration occurs here
  for (int i = 0; i < nVertices; ++i) { // Loop on edges
    double phiLower = phinumb[i];
    double phiUpper = phinumb[i + 1];
    double phiA = (phiUpper - phiLower) / 2.0;
    double phiB = (phiUpper + phiLower) / 2.0;
    double thetaLower = theta[numb[i]];
    double thetaUpper = theta[numb[i + 1]];
    double thetaMax = 0.0;
    Eigen::Vector3d oc = (arcs.col(i) - sph.center) / sph.radius;
    double oc_norm = oc.norm();
    double oc_norm2 = std::pow(oc_norm, 2);
    for (int j = 0; j < upper_phi; ++j) { // Loop on Gaussian points: phi integration
      for (int k = 0; k <= 1; ++k) {
        double ph = (2 * k - 1) * phiA * phiRule.abscissa(j) + phiB;
        double cos_ph = std::cos(ph);
        double sin_ph = std::sin(ph);
        // We need to calculate the upper bound for the integration on theta, which
        // depends on phi
        if (oc_norm2 < 1.0e-07) { // This means that the edge is centered on the same
                                  // sphere the tessera belongs to
          double cotg_thmax = (std::sin(ph - phiLower) / std::tan(thetaUpper) +
                               std::sin(phiUpper - ph) / std::tan(thetaLower)) /
                              std::sin(phiUpper - phiLower);
          thetaMax = std::atan(1.0 / cotg_thmax);
        } else {
          Eigen::Vector3d scratch;
          scratch << tangent.dot(oc), bitangent.dot(oc), normal.dot(oc);
          double aa =
              std::pow(tangent.dot(oc) * cos_ph + bitangent.dot(oc) * sin_ph, 2) +
              std::pow(normal.dot(oc), 2);
          double bb = -normal.dot(oc) * oc_norm2;
          double cc =
              std::pow(oc_norm2, 2) -
              std::pow(tangent.dot(oc) * cos_ph + bitangent.dot(oc) * sin_ph, 2);
          double ds = std::pow(bb, 2) - aa * cc;
          if (ds < 0.0)
            ds = 0.0;
          double cs = (-bb + std::sqrt(ds)) / aa;
          if (cs > 1.0)
            cs = 1.0;
          if (cs < -1.0)
            cs = 1.0;
          thetaMax = std::acos(cs);
        }
        double scratch = 0.0;
        if (!(thetaMax < 1.0e-08)) {
          double thetaA = thetaMax / 2.0;
          for (int l = 0; l < upper_theta;
               ++l) { // Loop on Gaussian points: theta integration
            for (int m = 0; m <= 1; ++m) {
              double th = (2 * m - 1) * thetaA * thetaRule.abscissa(l) + thetaA;
              double cos_th = std::cos(th);
              double sin_th = std::sin(th);
              Eigen::Vector3d point;
              point(0) = tangent(0) * sin_th * cos_ph +
                         bitangent(0) * sin_th * sin_ph + normal(0) * (cos_th - 1.0);
              point(1) = tangent(1) * sin_th * cos_ph +
                         bitangent(1) * sin_th * sin_ph + normal(1) * (cos_th - 1.0);
              point(2) = tangent(2) * sin_th * cos_ph +
                         bitangent(2) * sin_th * sin_ph + normal(2) * (cos_th - 1.0);
              double value =
                  F(point,
                    Eigen::Vector3d::Zero()); // Evaluate integrand at Gaussian point
              scratch += std::pow(sph.radius, 2) * value * sin_th * thetaA *
                         thetaRule.weight(l);
            }
          }
          result += scratch * phiA * phiRule.weight(j);
        }
      }
    }
  }
  return result;
}

template <int PhiPoints, int ThetaPoints>
double integrateD(const KernelD & F, const Element & e) {
  double result = 0.0;

  // Get the quadrature rules for azimuthal and polar integrations
  namespace mpl = boost::mpl;
  typedef typename mpl::at<rules_map, mpl::int_<PhiPoints> >::type PhiPolicy;
  typedef typename mpl::at<rules_map, mpl::int_<ThetaPoints> >::type ThetaPolicy;
  QuadratureRule<PhiPolicy> phiRule;
  QuadratureRule<ThetaPolicy> thetaRule;
  int upper_phi = PhiPoints / 2;     // Upper limit for loop on phi points
  int upper_theta = ThetaPoints / 2; // Upper limit for loop on theta points

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

  std::vector<double> theta(nVertices), phi(nVertices), phinumb(nVertices + 1);
  std::vector<int> numb(nVertices + 1);
  // Clean-up heap crap
  std::fill_n(theta.begin(), nVertices, 0.0);
  std::fill_n(phi.begin(), nVertices, 0.0);
  std::fill_n(numb.begin(), nVertices + 1, 0);
  std::fill_n(phinumb.begin(), nVertices + 1, 0.0);
  // Populate arrays and redefine tangent and bitangent
  e.spherical_polygon(tangent, bitangent, theta, phi, phinumb, numb);

  // Actual integration occurs here
  for (int i = 0; i < nVertices; ++i) { // Loop on edges
    double phiLower = phinumb[i];       // Lower vertex of edge
    double phiUpper = phinumb[i + 1];   // Upper vertex of edge
    double phiA = (phiUpper - phiLower) / 2.0;
    double phiB = (phiUpper + phiLower) / 2.0;
    double thetaLower = theta[numb[i]];
    double thetaUpper = theta[numb[i + 1]];
    double thetaMax = 0.0;
    Eigen::Vector3d oc = (arcs.col(i) - sph.center) / sph.radius;
    double oc_norm = oc.norm();
    double oc_norm2 = std::pow(oc_norm, 2);
    for (int j = 0; j < upper_phi; ++j) { // Loop on Gaussian points: phi integration
      for (int k = 0; k <= 1; ++k) {
        double ph = (2 * k - 1) * phiA * phiRule.abscissa(j) + phiB;
        double cos_phi = std::cos(ph);
        double sin_phi = std::sin(ph);
        if (oc_norm2 < 1.0e-07) { // This should check if oc_norm2 is zero
          double cotg_thmax = (std::sin(ph - phiLower) / std::tan(thetaUpper) +
                               std::sin(phiUpper - ph) / std::tan(thetaLower)) /
                              std::sin(phiUpper - phiLower);
          thetaMax = std::atan(1.0 / cotg_thmax);
        } else {
          Eigen::Vector3d scratch;
          scratch << tangent.dot(oc), bitangent.dot(oc), normal.dot(oc);
          double aa =
              std::pow(tangent.dot(oc) * cos_phi + bitangent.dot(oc) * sin_phi, 2) +
              std::pow(normal.dot(oc), 2);
          double bb = -normal.dot(oc) * oc_norm2;
          double cc =
              std::pow(oc_norm2, 2) -
              std::pow(tangent.dot(oc) * cos_phi + bitangent.dot(oc) * sin_phi, 2);
          double ds = std::pow(bb, 2) - aa * cc;
          if (ds < 0.0)
            ds = 0.0;
          double cs = (-bb + std::sqrt(ds)) / aa;
          if (cs > 1.0)
            cs = 1.0;
          if (cs < -1.0)
            cs = 1.0;
          thetaMax = std::acos(cs);
        }
        double thetaA = thetaMax / 2.0;
        double scratch = 0.0;
        if (!(thetaMax < 1.0e-08)) {
          for (int l = 0; l < upper_theta;
               ++l) { // Loop on Gaussian points: theta integration
            for (int m = 0; m <= 1; ++m) {
              double th = (2 * m - 1) * thetaA * thetaRule.abscissa(l) + thetaA;
              double cos_theta = std::cos(th);
              double sin_theta = std::sin(th);
              Eigen::Vector3d point;
              point(0) = tangent(0) * sin_theta * cos_phi +
                         bitangent(0) * sin_theta * sin_phi +
                         normal(0) * (cos_theta - 1.0);
              point(1) = tangent(1) * sin_theta * cos_phi +
                         bitangent(1) * sin_theta * sin_phi +
                         normal(1) * (cos_theta - 1.0);
              point(2) = tangent(2) * sin_theta * cos_phi +
                         bitangent(2) * sin_theta * sin_phi +
                         normal(2) * (cos_theta - 1.0);
              double value = F(e.normal(),
                               Eigen::Vector3d::Zero(),
                               point); // Evaluate integrand at Gaussian point
              scratch += std::pow(sph.radius, 2) * value * sin_theta * thetaA *
                         thetaRule.weight(l);
            }
          }
          result += scratch * phiA * phiRule.weight(j);
        }
      }
    }
  }
  return result;
}

template double integrateS<32, 16>(const KernelS & F, const Element & e);
template double integrateD<32, 16>(const KernelD & F, const Element & e);

template double integrateS<64, 16>(const KernelS & F, const Element & e);
template double integrateD<64, 16>(const KernelD & F, const Element & e);
} // namespace bi_operators
} // namespace pcm
