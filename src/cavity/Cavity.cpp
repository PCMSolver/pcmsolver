/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#include "Cavity.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>

#include "Config.hpp"

#include <Eigen/Core>
#include "utils/cnpy.hpp"

#include "utils/MathUtils.hpp"
#include "utils/Symmetry.hpp"

void Cavity::saveCavity(const std::string & fname)
{
  /*
     std::ofstream weights("weights.txt", std::ios_base::out);
  // First line in the weights file is the number of elements.
  // This is for a sanity-check on the save/load operations.
  weights << std::setprecision(std::numeric_limits<double>::digits10) << nElements_ << std::endl;
  weights << std::setprecision(std::numeric_limits<double>::digits10) << elementArea_ << std::endl;
  std::ofstream elRadius("element_radius.txt", std::ios_base::out);
  elRadius << std::setprecision(std::numeric_limits<double>::digits10) << elementRadius_ << std::endl;
  std::ofstream centers("centers.txt", std::ios_base::out);
  centers << std::setprecision(std::numeric_limits<double>::digits10) << elementCenter_ << std::endl;
  std::ofstream normals("normals.txt", std::ios_base::out);
  normals << std::setprecision(std::numeric_limits<double>::digits10) << elementNormal_ << std::endl;
  */

  // Write everything in a single .npz binary file
  unsigned int dim = static_cast<unsigned int>(nElements_);
  // Write the number of elements, it will be used to check sanity of the save/load operations.
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
  // Write vertices TODO!!!!
  //const unsigned int vertices_shape[] = {dim, 10, 3};
  // Write arcs     TODO!!!!
  //const unsigned int arcs_shape[] = {dim, 10, 3};
}

void Cavity::loadCavity(const std::string & fname)
{
  // We initialize molecule_ to a dummy Molecule
  molecule_ = Molecule();
  // Load the .npz binary file and then traverse it to get the data needed to rebuild the cavity.
  cnpy::npz_t loaded_cavity = cnpy::npz_load(fname);
  // 0. Get the number of elements
  cnpy::NpyArray raw_ele = loaded_cavity["elements"];
  int * ne = reinterpret_cast<int*>(raw_ele.data);
  nElements_ = *ne;
  // Set the size of the irreducible portion of the cavity
  // it will be the same as the total size, since a restarted cavity is always C1
  nIrrElements_ = nElements_;

  // 1. Get the weights
  cnpy::NpyArray raw_weights = loaded_cavity["weights"];
  if (raw_weights.shape[0] != nElements_) PCMSOLVER_ERROR("elementArea_: incoherent dimensions read in", BOOST_CURRENT_FUNCTION);
  elementArea_ = cnpy::custom::npy_to_eigen(raw_weights);
  // 2. Get the element sphere center
  cnpy::NpyArray raw_elSphCenter = loaded_cavity["elSphCenter"];
  if (raw_elSphCenter.shape[1] != nElements_) PCMSOLVER_ERROR("elementSphereCenter_: incoherent dimensions read in", BOOST_CURRENT_FUNCTION);
  elementSphereCenter_ = cnpy::custom::npy_to_eigen(raw_elSphCenter);
  // 3. Get the element radius
  cnpy::NpyArray raw_elRadius = loaded_cavity["elRadius"];
  if (raw_elRadius.shape[0] != nElements_) PCMSOLVER_ERROR("elementRadius_: incoherent dimensions read in", BOOST_CURRENT_FUNCTION);
  elementRadius_ = cnpy::custom::npy_to_eigen(raw_elRadius);
  // 4. Get the centers
  cnpy::NpyArray raw_centers = loaded_cavity["centers"];
  if (raw_centers.shape[1] != nElements_) PCMSOLVER_ERROR("elementCenter_: incoherent dimensions read in", BOOST_CURRENT_FUNCTION);
  elementCenter_ = cnpy::custom::npy_to_eigen(raw_centers);
  // 5. Get the normal vectors
  cnpy::NpyArray raw_normals = loaded_cavity["normals"];
  if (raw_normals.shape[1] != nElements_) PCMSOLVER_ERROR("elementNormal_: incoherent dimensions read in", BOOST_CURRENT_FUNCTION);
  elementNormal_ = cnpy::custom::npy_to_eigen(raw_normals);

  // Reconstruct the elements_ vector
  for (size_t i = 0; i < nElements_; ++i) {
    bool irr = false;
    // PEDRA puts the irreducible tesserae first
    if (i < nIrrElements_) irr = true;
    Sphere sph(elementSphereCenter_.col(i), elementRadius_(i));
    int nv = 3; // BOGUS!!!
    Eigen::Matrix3Xd vertices, arcs;
    vertices.resize(Eigen::NoChange, nv); // BOGUS!!!
    arcs.resize(Eigen::NoChange, nv); // BOGUS!!
    // Populate vertices and arcs
    elements_.push_back(Element(nv, 0,
          elementArea_(i),
          elementCenter_.col(i),
          elementNormal_.col(i),
          irr, sph,
          vertices, arcs));
  }
}
