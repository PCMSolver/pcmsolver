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

#include "TsLessCavity.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "CavityData.hpp"
#include "utils/Sphere.hpp"
#include "utils/Symmetry.hpp"

namespace pcm {
namespace cavity {
TsLessCavity::TsLessCavity(const Molecule & molec,
                           double a,
                           double pr,
                           double minR,
                           double minD,
                           int der)
    : ICavity(molec),
      averageArea_(a),
      probeRadius_(pr),
      minimalRadius_(minR),
      minDistance_(minD),
      derOrder_(der) {
  TIMER_ON("TsLessCavity::build from Molecule");
  build(10000, 200, 25000);
  TIMER_OFF("TsLessCavity::build from Molecule");
}

TsLessCavity::TsLessCavity(const Sphere & sph,
                           double a,
                           double pr,
                           double minR,
                           double minD,
                           int der)
    : ICavity(sph),
      averageArea_(a),
      probeRadius_(pr),
      minimalRadius_(minR),
      minDistance_(minD),
      derOrder_(der) {
  TIMER_ON("TsLessCavity::build from a single sphere");
  build(10000, 200, 25000);
  TIMER_OFF("TsLessCavity::build from a single sphere");
}

TsLessCavity::TsLessCavity(const std::vector<Sphere> & sph,
                           double a,
                           double pr,
                           double minR,
                           double minD,
                           int der)
    : ICavity(sph),
      averageArea_(a),
      probeRadius_(pr),
      minimalRadius_(minR),
      minDistance_(minD),
      derOrder_(der) {
  TIMER_ON("TsLessCavity::build from list of spheres");
  build(10000, 200, 25000);
  TIMER_OFF("TsLessCavity::build from list of spheres");
}

void TsLessCavity::build(int maxts, int maxsph, int maxvert) {
  // This is a wrapper for the generatecavity_cpp_ function defined in the Fortran
  // code TsLess. Here we allocate the necessary arrays to be passed to TsLess, in
  // particular we allow for the insertion of additional spheres as in the most
  // general formulation of the GePol algorithm.

  int lwork = maxts * maxsph;
  double * xtscor = new double[maxts];
  double * ytscor = new double[maxts];
  double * ztscor = new double[maxts];
  double * ar = new double[maxts];
  double * xsphcor = new double[maxts];
  double * ysphcor = new double[maxts];
  double * zsphcor = new double[maxts];
  double * rsph = new double[maxts];
  double * work = new double[lwork];

  // Clean-up possible heap-crap
  std::fill_n(xtscor, maxts, 0.0);
  std::fill_n(ytscor, maxts, 0.0);
  std::fill_n(ztscor, maxts, 0.0);
  std::fill_n(ar, maxts, 0.0);
  std::fill_n(xsphcor, maxts, 0.0);
  std::fill_n(ysphcor, maxts, 0.0);
  std::fill_n(zsphcor, maxts, 0.0);
  std::fill_n(rsph, maxts, 0.0);
  std::fill_n(work, lwork, 0.0);

  int nts = 0;
  int ntsirr = 0;

  // If there's an overflow in the number of spheres TsLess will die.
  // The maximum number of spheres in PEDRA is set to 200 (primitive+additional)
  // so the integer here declared is just to have enough space C++ side to pass
  // everything back.
  int maxAddedSpheres = 200;

  // maximum number of spheres we allow the algorithm to add to our original set.
  // If this number is exceeded, then the algorithm crashes (should look into
  // this...) After the cavity is generated we will update ALL the class data
  // members, both related to spheres and finite elements so that the cavity is fully
  // formed.

  Eigen::VectorXd xv = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
  Eigen::VectorXd yv = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
  Eigen::VectorXd zv = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
  Eigen::VectorXd radii_scratch =
      Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres); // Not to be confused with
                                                          // the data member
                                                          // inherited from Cavity!!!

  for (int i = 0; i < nSpheres_; ++i) {
    for (int j = 0; j < 3; ++j) {
      xv(i) = sphereCenter_(0, i);
      yv(i) = sphereCenter_(1, i);
      zv(i) = sphereCenter_(2, i);
    }
    radii_scratch(i) = sphereRadius_(i);
  }

  double * xe = xv.data();
  double * ye = yv.data();
  double * ze = zv.data();

  double * rin = radii_scratch.data();
  double * mass = new double[molecule_.nAtoms()];
  for (size_t i = 0; i < molecule_.nAtoms(); ++i) {
    mass[i] = molecule_.masses(i);
  }

  addedSpheres = 0;
  // Number of generators and generators of the point group
  int nr_gen = molecule_.pointGroup().nrGenerators();
  int gen1 = molecule_.pointGroup().generators(0);
  int gen2 = molecule_.pointGroup().generators(1);
  int gen3 = molecule_.pointGroup().generators(2);

  int weightFunction = 1;
  // Go TsLess, Go!
  TIMER_ON("TsLessCavity::tsless_driver");
  tsless_driver(&maxts,
                &maxsph,
                &maxvert,
                &nSpheres_,
                &nts,
                &ntsirr,
                &addedSpheres,
                xtscor,
                ytscor,
                ztscor,
                ar,
                xsphcor,
                ysphcor,
                zsphcor,
                rsph,
                xe,
                ye,
                ze,
                rin,
                mass,
                &nr_gen,
                &gen1,
                &gen2,
                &gen3,
                &averageArea_,
                &minDistance_,
                &derOrder_,
                &weightFunction,
                &probeRadius_,
                work);
  TIMER_OFF("TsLessCavity::tsless_driver");

  // The "intensive" part of updating the spheres related class data members will be
  // of course executed iff addedSpheres != 0
  if (addedSpheres != 0) {
    // Save the number of original spheres
    int orig = nSpheres_;
    // Update the nSpheres
    nSpheres_ += addedSpheres;
    // Resize sphereRadius and sphereCenter...
    sphereRadius_.resize(nSpheres_);
    sphereCenter_.resize(Eigen::NoChange, nSpheres_);
    // Transfer radii and centers.
    // Eigen has no push_back function, so we need to traverse all the spheres...
    sphereRadius_ = radii_scratch.head(nSpheres_);
    for (int i = 0; i < nSpheres_; ++i) {
      sphereCenter_(0, i) = xv(i);
      sphereCenter_(1, i) = yv(i);
      sphereCenter_(2, i) = zv(i);
    }
    // Now grow the vector<Sphere> containing the list of spheres
    for (int i = orig; i < nSpheres_; ++i) {
      spheres_.push_back(Sphere(sphereCenter_.col(i), sphereRadius_(i)));
    }
  }

  nElements_ = static_cast<int>(nts);
  nIrrElements_ = static_cast<int>(ntsirr);
  elementCenter_.resize(Eigen::NoChange, nElements_);
  elementSphereCenter_.resize(Eigen::NoChange, nElements_);
  elementNormal_.resize(Eigen::NoChange, nElements_);
  elementArea_.resize(nElements_);
  elementRadius_.resize(nElements_);
  for (int i = 0; i < nElements_; ++i) {
    elementCenter_(0, i) = xtscor[i];
    elementCenter_(1, i) = ytscor[i];
    elementCenter_(2, i) = ztscor[i];
    elementArea_(i) = ar[i];
    elementSphereCenter_(0, i) = xsphcor[i];
    elementSphereCenter_(1, i) = ysphcor[i];
    elementSphereCenter_(2, i) = zsphcor[i];
    elementRadius_(i) = rsph[i];
  }

  elementNormal_ = elementCenter_ - elementSphereCenter_;
  for (int i = 0; i < nElements_; ++i) {
    elementNormal_.col(i) /= elementNormal_.col(i).norm();
  }

  // Fill elements_ vector
  for (int i = 0; i < nElements_; ++i) {
    bool irr = false;
    int nv = 1; // TsLess does not generate spherical polygons!!
    // TsLess puts the irreducible tesserae first (? Check with Cris!)
    if (i < nIrrElements_)
      irr = true;
    Sphere sph(elementSphereCenter_.col(i), elementRadius_(i));
    Eigen::Matrix3Xd vertices, arcs;
    vertices.resize(Eigen::NoChange, nv);
    arcs.resize(Eigen::NoChange, nv);
    // FIXME index of the sphere the element belongs to
    elements_.push_back(Element(nv,
                                0,
                                elementArea_(i),
                                elementCenter_.col(i),
                                elementNormal_.col(i),
                                irr,
                                sph,
                                vertices,
                                arcs));
  }

  delete[] xtscor;
  delete[] ytscor;
  delete[] ztscor;
  delete[] ar;
  delete[] xsphcor;
  delete[] ysphcor;
  delete[] zsphcor;
  delete[] rsph;
  delete[] work;
  delete[] mass;

  built = true;
}

std::ostream & TsLessCavity::printCavity(std::ostream & os) {
  os << "Cavity type: TsLess" << std::endl;
  os << "Average point weight = " << averageArea_ << " AU^2" << std::endl;
  os << "Minimal distance between sampling points = " << minDistance_ << " AU"
     << std::endl;
  os << "Switch function is of class C^" << derOrder_ << std::endl;
  os << "Addition of extra spheres enabled" << std::endl;
  os << "Probe radius = " << probeRadius_ << " AU" << std::endl;
  os << "Number of spheres = " << nSpheres_
     << " [initial = " << nSpheres_ - addedSpheres << "; added = " << addedSpheres
     << "]" << std::endl;
  os << "Number of finite elements = " << nElements_;
  os << "Number of irreducible finite elements = " << nIrrElements_ << std::endl;
  os << "============ Spheres list (in Angstrom)" << std::endl;
  os << " Sphere   on   Radius   Alpha       X            Y            Z     \n";
  os << "-------- ---- -------- ------- -----------  -----------  -----------\n";
  // Print original set of spheres
  int original = nSpheres_ - addedSpheres;
  Eigen::IOFormat CleanFmt(6, Eigen::DontAlignCols, "     ", "\n", "", "");
  for (int i = 0; i < original; ++i) {
    os << std::setw(4) << i + 1;
    os << "      " << molecule_.atoms()[i].symbol << "    ";
    os << std::fixed << std::setprecision(4) << molecule_.atoms()[i].radius;
    os << std::fixed << std::setprecision(2) << "   "
       << molecule_.atoms()[i].radiusScaling << "     ";
    os << (molecule_.geometry().col(i).transpose() * bohrToAngstrom())
              .format(CleanFmt);
    os << std::endl;
  }
  // Print added spheres
  for (int j = 0; j < addedSpheres; ++j) {
    int idx = original + j;
    os << std::setw(4) << idx + 1;
    os << "      Du    ";
    os << std::fixed << std::setprecision(4)
       << sphereRadius_(idx) * bohrToAngstrom();
    os << std::fixed << std::setprecision(2) << "   1.00";
    os << (sphereCenter_.col(idx).transpose() * bohrToAngstrom()).format(CleanFmt);
    os << std::endl;
  }
  return os;
}

ICavity * createTsLessCavity(const CavityData & data) {
  return new TsLessCavity(data.molecule,
                          data.area,
                          data.probeRadius,
                          data.minimalRadius,
                          data.minDistance,
                          data.derOrder);
}
} // namespace cavity
} // namespace pcm
