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

#include "Cavity.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <stdexcept>

#include "Config.hpp"


#include <Eigen/Dense>
#include "cnpyPimpl.hpp"

#include "Symmetry.hpp"

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
    const unsigned int weights_shape[] = {dim};
    cnpy::npz_save(fname, "weights", elementArea_.data(), weights_shape, 1, "a", true);
    // Write element radius
    const unsigned int elRadius_shape[] = {dim};
    cnpy::npz_save(fname, "elRadius", elementRadius_.data(), elRadius_shape, 1, "a",
                   true);
    // Write centers
    const unsigned int centers_shape[] = {3, dim};
    cnpy::npz_save(fname, "centers", elementCenter_.data(), centers_shape, 2, "a", true);
    // Write normals
    const unsigned int normals_shape[] = {3, dim};
    cnpy::npz_save(fname, "normals", elementNormal_.data(), normals_shape, 2, "a", true);
}

void Cavity::loadCavity(const std::string & fname)
{
    // If the cavity has been loaded from file, the point group is C1
    pointGroup_ = buildGroup(0, 0, 0, 0);
    // Load the .npz binary file and then traverse it to get the data needed to rebuild the cavity.
    cnpy::npz_t loaded_cavity = cnpy::npz_load(fname);
    // 0. Get the number of elements
    cnpy::NpyArray raw_ele = loaded_cavity["elements"];
    int * ne = reinterpret_cast<int*>(raw_ele.data);
    nElements_ = *ne;
    // Set the size of the irreducible portion of the cavity
    // it will be the same as the total size, since a restarted cavity is always C1
    nIrrElements_ = (nElements_ / pointGroup_.nrIrrep());

    // 1. Get the weights
    cnpy::NpyArray raw_weights = loaded_cavity["weights"];
    int dim = raw_weights.shape[0];
    if (dim != nElements_) {
        throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of weights vector!");
    } else {
        elementArea_ = cnpy::getFromRawBuffer<double>(dim, 1, raw_weights.data);
    }

    // 2. Get the element radius
    cnpy::NpyArray raw_elRadius = loaded_cavity["elRadius"];
    dim = raw_elRadius.shape[0];
    if (dim != nElements_) {
        throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of element radius vector!");
    } else {
        elementRadius_ = cnpy::getFromRawBuffer<double>(dim, 1, raw_elRadius.data);
    }

    // 3. Get the centers
    cnpy::NpyArray raw_centers = loaded_cavity["centers"];
    dim = raw_centers.shape[1];
    if (dim != nElements_) {
        throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of centers matrix!");
    } else {
        elementCenter_ = cnpy::getFromRawBuffer<double>(3, dim, raw_centers.data);
    }

    // 4. Get the normal vectors
    cnpy::NpyArray raw_normals = loaded_cavity["normals"];
    dim = raw_normals.shape[1];
    if (dim != nElements_) {
        throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of normals matrix!");
    } else {
        elementNormal_ = cnpy::getFromRawBuffer<double>(3, dim, raw_normals.data);
    }
}
