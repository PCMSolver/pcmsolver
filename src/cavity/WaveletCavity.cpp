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

#include "WaveletCavity.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

#include <boost/foreach.hpp>

#include "Vector2.hpp"
#include "Vector3.hpp"
#include "Interpolation.hpp"
#include "Cubature.hpp"
#include "GaussSquare.hpp"

void WaveletCavity::writeInput(std::string &fileName)
{
    std::ofstream output;
    output.open(fileName.c_str(), std::fstream::out);

    output << nSpheres_ << std::endl;
    output.setf(std::ios_base::showpoint);
    output.precision(12);
    for(int i=0; i < nSpheres_; i++) {
        output << sphereCenter_(0,i) << " ";
        output << sphereCenter_(1,i) << " ";
        output << sphereCenter_(2,i) << " ";
        output << sphereRadius_(i) << std::endl;
    }
    output.close();
    uploadedDyadic_ = false;
}

extern "C"
{
    int waveletCavityDrv_(double pr, double c, int pl, const char * infile);
}

void WaveletCavity::makeCavity()
{
    int check = 0;
    std::string infile = "cavity.inp";
    writeInput(infile);
#if defined (WAVELET_DEVELOPMENT)
    check = waveletCavityDrv_(probeRadius_, coarsity_, patchLevel_, infile.c_str());
#else
    check = 1;
#endif
    if (check != 0) {
        throw std::runtime_error("Problem with the wavelet cavity inside makeCavity method");
    }
}


void WaveletCavity::readCavity(const std::string & filename)
{
    size_t i, j, k;
    double x, y, z;

    std::ifstream file;
    file.open(filename);
    if (file.is_open()) {
        file >> nLevels_ >> nPatches_;           
                                                
        int nNodes = (1 << nLevels_) + 1;
                                                
        nPoints_ = nPatches_ * nNodes * nNodes;

                                                
        for (size_t l = 0; l < nPoints_; ++l) {
            file >> i >> j >> k >> x >> y >> z;
            Eigen::Vector3i index(i, k, j);
            nodeIndex_.push_back(index);
            Eigen::Vector3d point(x, y, z);
            nodePoint_.push_back(point);
        }
                                                
        file.close();
        uploadedDyadic_ = true;
    } else {
        throw std::runtime_error("Dyadic file not opened properly");
    }
}

void WaveletCavity::scaleCavity(const double scalingFactor)
{
    BOOST_FOREACH(Eigen::Vector3d & point, nodePoint_) {
        point *= scalingFactor;
    }
}

void WaveletCavity::uploadPoints(int quadLevel, Interpolation *interp)
{
    if (not uploadedDyadic_) {
        throw std::runtime_error("Dyadic file must be uploaded first.");
    }
    Vector2 s, t;
    Vector3 point;
    Vector3 norm;
    int n = 1 << nLevels_;
    double h = 1.0 / n;
    Cubature *Q;
    initGaussSquare(&Q, quadLevel + 1);

    nElements_ = nPatches_ * n * n * Q[quadLevel].noP;
    // We set nIrrElements_ to nElements_ since there is no symmetry!
    nIrrElements_ = nElements_;

    elementCenter_.resize(Eigen::NoChange, nElements_);
    elementNormal_.resize(Eigen::NoChange, nElements_);
    elementArea_.resize(nElements_);
    // The following is to guarantee that writing cavity to .npz file works
    elementRadius_.setZero(nElements_);

    int j = 0;
    for (size_t i1 = 0; i1 < nPatches_; ++i1) {
        for (int i2 = 0; i2 < n; ++i2) {
            s.y = h * i2;
            for (int i3=0; i3 < n; ++i3) {
                s.x = h * i3;
                for (size_t k = 0; k < Q[quadLevel].noP; k++) {
                    t = vector2Add(s,vector2SMul(h,Q[quadLevel].xi[k]));
                    point = interp->Chi(t,i1);
                    norm = interp->n_Chi(t,i1);
                    Eigen::Vector3d center(point.x, point.y, point.z);
                    Eigen::Vector3d normal(norm.x,  norm.y,  norm.z);
                    normal.normalize();
                    double area = h * h * Q[quadLevel].weight[k] * vector3Norm(norm);
                    elementCenter_.col(j) = center.transpose();
                    elementNormal_.col(j) = normal.transpose();
                    elementArea_(j) = area;
                    ++j;
                }
            }
        }
    }
    freeGaussSquare(&Q,quadLevel+1);
    built = true;
}

std::ostream & WaveletCavity::printCavity(std::ostream & os)
{
    // Calculation of quadrature level for printout
    os << "Cavity type: Wavelet" << std::endl;
    os << "Probe Radius =  " << probeRadius_ << std::endl;
    os << "Coarsity     =  " << coarsity_ << std::endl;
    os << "Patch Level  =  " << patchLevel_ << std::endl;
    os << "Number of spheres = " << nSpheres_ << std::endl;
    os << "Number of patches = " << nPatches_ << std::endl;
    os << "Number of levels  = " << nLevels_  << std::endl;
    os << "Number of potential points = " << nElements_;
    /*
    os << std::endl; 
    os << "-------------- Potential points printout " << std::endl;
    for(int i = 0; i < nElements_; i++) {
    	os << std::endl;
    	os << "Point " << i+1 << " \n";
    	os << elementCenter_(0,i) << " ";
    	os << elementCenter_(1,i) << " ";
    	os << elementCenter_(2,i) << " ";
    	os << elementArea_(i) << " ";
    }
    */
    /*
    for(int i = 0; i < nSpheres_; i++) {
    	os << endl;
    	os << i+1 << " ";
    	os << sphereCenter_(0,i) << " ";
    	os << sphereCenter_(1,i) << " ";
    	os << sphereCenter_(2,i) << " ";
    	os << sphereRadius_(i) << " ";
    }
    */
    /*	
    if (uploadedDyadic_) {
        os << "Printing nodes" << std::endl;
    	for(int i = 0; i < nPoints_; i++) {
    	    os << std::endl;
    	    os << i+1 << " ";
    	    os << nodeIndex_[i].transpose() << " " << nodePoint_[i].transpose() << " ";
    	}
    }
    */
    return os;
}
