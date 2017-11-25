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

#include "ICavity.hpp"

#include <iostream>

#include "Config.hpp"

#include <Eigen/Core>

#include "utils/MathUtils.hpp"
#include "utils/Sphere.hpp"
#include "utils/Symmetry.hpp"
#include "utils/cnpy.hpp"

namespace pcm {
ICavity::ICavity() : nElements_(0), built(false) {}

ICavity::ICavity(const Sphere & sph) : built(false) {
  spheres_.push_back(sph);
  molecule_ = Molecule(spheres_);
  nSpheres_ = spheres_.size();
  transfer_spheres(spheres_, sphereCenter_, sphereRadius_);
}

ICavity::ICavity(const std::vector<Sphere> & sph) : spheres_(sph), built(false) {
  molecule_ = Molecule(spheres_);
  nSpheres_ = spheres_.size();
  transfer_spheres(spheres_, sphereCenter_, sphereRadius_);
}

ICavity::ICavity(const Molecule & molec)
    : spheres_(molec.spheres()), molecule_(molec), built(false) {
  nSpheres_ = spheres_.size();
  transfer_spheres(spheres_, sphereCenter_, sphereRadius_);
}

void ICavity::saveCavity(const std::string & fname) {
  // Write everything in a single .npz binary file
  unsigned int dim = static_cast<unsigned int>(nElements_);
  // Write the number of elements, it will be used to check sanity of the save/load
  // operations.
  const unsigned int shape[] = {1};
  cnpy::npz_save(fname, "elements", &dim, shape, 1, "w", false);
  // Write weights
  cnpy::custom::npz_save(fname, "weights", elementArea_);
  // Write element sphere center
  cnpy::custom::npz_save(fname, "elSphCenter", elementSphereCenter_);
  // Write element radius
  cnpy::custom::npz_save(fname, "elRadius", elementRadius_);
  // Write centers
  cnpy::custom::npz_save(fname, "centers", elementCenter_);
  // Write normals
  cnpy::custom::npz_save(fname, "normals", elementNormal_);
  for (PCMSolverIndex i = 0; i < nElements_; ++i) {
    // Write vertices
    cnpy::custom::npz_save(
        fname, "vertices_" + pcm::to_string(i), elements_[i].vertices());
    // Write arcs
    cnpy::custom::npz_save(fname, "arcs_" + pcm::to_string(i), elements_[i].arcs());
  }
}

void ICavity::loadCavity(const std::string & fname) {
  // We initialize molecule_ to a dummy Molecule
  molecule_ = Molecule();
  // Load the .npz binary file and then traverse it to get the data needed to rebuild
  // the cavity.
  cnpy::npz_t loaded_cavity = cnpy::npz_load(fname);
  // 0. Get the number of elements
  cnpy::NpyArray raw_ele = loaded_cavity["elements"];
  int * ne = reinterpret_cast<int *>(raw_ele.data);
  nElements_ = *ne;
  // Set the size of the irreducible portion of the cavity
  // it will be the same as the total size, since a restarted cavity is always C1
  nIrrElements_ = nElements_;

  // 1. Get the weights
  cnpy::NpyArray raw_weights = loaded_cavity["weights"];
  if (raw_weights.shape[0] != nElements_)
    PCMSOLVER_ERROR("elementArea_: incoherent dimensions read in");
  elementArea_ = cnpy::custom::npy_to_eigen<double>(raw_weights);
  // 2. Get the element sphere center
  if (loaded_cavity.find("elSphCenter") == loaded_cavity.end()) {
    // Element sphere center was not found on file, fill it up with zeros
    elementSphereCenter_ = Eigen::Matrix3Xd::Zero(3, nElements_);
  } else {
    cnpy::NpyArray raw_elSphCenter = loaded_cavity["elSphCenter"];
    if (raw_elSphCenter.shape[1] != nElements_)
      PCMSOLVER_ERROR("elementSphereCenter_: incoherent dimensions read in");
    elementSphereCenter_ = cnpy::custom::npy_to_eigen<double>(raw_elSphCenter);
  }
  // 3. Get the element radius
  cnpy::NpyArray raw_elRadius = loaded_cavity["elRadius"];
  if (raw_elRadius.shape[0] != nElements_)
    PCMSOLVER_ERROR("elementRadius_: incoherent dimensions read in");
  elementRadius_ = cnpy::custom::npy_to_eigen<double>(raw_elRadius);
  // 4. Get the centers
  cnpy::NpyArray raw_centers = loaded_cavity["centers"];
  if (raw_centers.shape[1] != nElements_)
    PCMSOLVER_ERROR("elementCenter_: incoherent dimensions read in");
  elementCenter_ = cnpy::custom::npy_to_eigen<double>(raw_centers);
  // 5. Get the normal vectors
  cnpy::NpyArray raw_normals = loaded_cavity["normals"];
  if (raw_normals.shape[1] != nElements_)
    PCMSOLVER_ERROR("elementNormal_: incoherent dimensions read in");
  elementNormal_ = cnpy::custom::npy_to_eigen<double>(raw_normals);

  bool has_arcs = loaded_cavity.find("arcs_0") == loaded_cavity.end() ? false : true;
  bool has_vertices =
      loaded_cavity.find("vertices_0") == loaded_cavity.end() ? false : true;
  // Reconstruct the elements_ vector
  for (PCMSolverIndex i = 0; i < nElements_; ++i) {
    bool irr = false;
    // PEDRA puts the irreducible tesserae first
    if (i < nIrrElements_)
      irr = true;
    Sphere sph(elementSphereCenter_.col(i), elementRadius_(i));
    Eigen::Matrix3Xd vertices, arcs;
    // 6. Get vertices and arcs
    if (has_vertices) {
      cnpy::NpyArray raw_vertices = loaded_cavity["vertices_" + pcm::to_string(i)];
      vertices = cnpy::custom::npy_to_eigen<double>(raw_vertices);
    } else {
      // Vertices were not found on file, fill them up with zeros
      vertices = Eigen::Matrix3Xd::Zero(3, 3);
    }
    if (has_arcs) {
      cnpy::NpyArray raw_arcs = loaded_cavity["arcs_" + pcm::to_string(i)];
      arcs = cnpy::custom::npy_to_eigen<double>(raw_arcs);
    } else {
      // Arcs were not found on file, fill them up with zeros
      arcs = Eigen::Matrix3Xd::Zero(3, 3);
    }
    if (arcs.cols() != vertices.cols())
      PCMSOLVER_ERROR("Inconsistent number of vertices read from file for element " +
                      pcm::to_string(i));
    int nv = vertices.cols();
    // Populate vertices and arcs
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
}
} // namespace pcm
