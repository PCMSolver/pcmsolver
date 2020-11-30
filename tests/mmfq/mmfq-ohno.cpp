/**
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

#include "catch.hpp"

#include <iostream>

#include <Eigen/Core>

#include "mmfq/FQOhno.hpp"
#include "utils/MMFQ.hpp"
#include "utils/MathUtils.hpp"

using namespace pcm;
using mmfq::FQOhno;
using utils::MMFQ;

TEST_CASE("Test MMFQ for a pair of water fragments and the Ohno kernel",
          "[mmfq][ohno][mmfq-ohno]") {
  MMFQ fragments;
  fragments.nFragments = 2;
  fragments.nSitesPerFragment = 3;
  int dim = fragments.nFragments * fragments.nSitesPerFragment;
  fragments.chi = Eigen::VectorXd::Zero(dim);
  fragments.chi(0) = 0.189194;
  fragments.chi(1) = 0.012767;
  fragments.chi(2) = 0.012767;
  fragments.chi(3) = 0.189194;
  fragments.chi(4) = 0.012767;
  fragments.chi(5) = 0.012767;

  fragments.eta = Eigen::VectorXd::Zero(dim);
  fragments.eta(0) = 0.623700;
  fragments.eta(1) = 0.637512;
  fragments.eta(2) = 0.637512;
  fragments.eta(3) = 0.623700;
  fragments.eta(4) = 0.637512;
  fragments.eta(5) = 0.637512;

  fragments.sites = Eigen::Matrix3Xd::Zero(3, dim);
  Eigen::Matrix3Xd tmp = Eigen::Matrix3Xd::Zero(3, dim);
  tmp << 1.9602622018, 1.2008654437, 1.7053244255, 1.5747753414, 1.9507811044,
      2.3218071081, 0.2129635729, -0.4083959416, 0.8826907262, 0.2954197242,
      0.4993935756, 0.0746990811, -1.0911273996, -1.0103265807, -1.7421569411,
      1.6434710664, 0.7568211662, 2.2184343370;
  fragments.sites = tmp * angstromToBohr();
  FQOhno ohno(fragments);
  Eigen::VectorXd potential = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd charge = ohno.computeCharge(potential);

  Eigen::VectorXd reference = cnpy::custom::npy_load<double>("FQ-2xH2O.npy");
  for (int i = 0; i < dim; ++i) {
    REQUIRE(reference(i) == Approx(charge(i)));
  }

  double ref_energy = -0.12133410792682;
  double energy = 0.5 * (charge.dot(fragments.chi));
  REQUIRE(ref_energy == Approx(energy));
  CAPTURE(ref_energy - energy);
}
