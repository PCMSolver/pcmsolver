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

extern "C"
{
//#include "WEM.h"
//#include "read_points.h"
#include "vector2.h"
#include "vector3.h"
#include "interpolate.h"
#include "interpolate_pwl.h"
//#include "topology.h"
//#include "kern.h"
//#include "compression.h"
//#include "postproc.h"
//#include "WEMRHS.h"
//#include "WEMPCG.h"
//#include "WEMPGMRES.h"
//#include "dwt.h"
#include "cubature.h"
#include "gauss_square.h"
//#include "constants.h"
}

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
    size_t i, j;
    double x, y, z;

    std::ifstream file;
    file.open(filename.c_str());
    file >> nLevels_ >> nPatches_;

    int nNodes = (1 << nLevels_) + 1;

    nPoints_ = nPatches_ * nNodes * nNodes;

    for (size_t k = 0; k < nPoints_; ++k) {
        file >> i >> j >> k >> x >> y >> z;
        Eigen::Vector3i index(i, j, k);
        nodeIndex_.push_back(index);
        Eigen::Vector3d point(x, y, z);
        nodePoint_.push_back(point);
    }

    file.close();
    uploadedDyadic_ = true;
}

void WaveletCavity::uploadPoints(int quadLevel, vector3 **** T_, bool isPWL)
{
    if(isPWL) {
        uploadPointsPWL(quadLevel, T_);
    } else {
        uploadPointsPWC(quadLevel, T_);
    }
}

void WaveletCavity::uploadPointsPWC(int quadLevel, vector3 **** T_)
{
    if (not uploadedDyadic_) {
        throw std::runtime_error("Dyadic file must be uploaded first.");
    }
    vector2 s, t;
    vector3 point;
    vector3 norm;
    int n = 1 << nLevels_;
    double h = 1.0 / n;
    cubature *Q;
    init_Gauss_Square(&Q, quadLevel + 1);

    nElements_ = nPatches_ * n * n * Q[quadLevel].nop;

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
                for (size_t k = 0; k < Q[quadLevel].nop; ++k) {
                    t = vector2_add(s,vector2_Smul(h,Q[quadLevel].xi[k]));
                    point = Chi(t,T_[i1], nLevels_);
                    norm = n_Chi(t,T_[i1], nLevels_);
                    Eigen::Vector3d center(point.x, point.y, point.z);
                    Eigen::Vector3d normal(norm.x,  norm.y,  norm.z);
                    normal.normalize();
                    double area = h * h * Q[quadLevel].w[k] * vector3_norm(n_Chi(t, T_[i1], nLevels_));
                    elementCenter_.col(j) = center.transpose();
                    elementNormal_.col(j) = normal.transpose();
                    elementArea_(j) = area;
                    ++j;
                }
            }
        }
    }
    free_Gauss_Square(&Q,quadLevel+1);
    built = true;
}

void WaveletCavity::uploadPointsPWL(int quadLevel, vector3 **** T_)
{
    if (not uploadedDyadic_) {
        throw std::runtime_error("Dyadic file must be uploaded first.");
    }
    vector2 s, t;
    vector3 point;
    vector3 norm;
    int n = 1 << nLevels_;
    double h = 1.0 / n;
    cubature *Q;
    init_Gauss_Square(&Q, quadLevel + 1);

    nElements_ = nPatches_ * n * n * Q[quadLevel].nop;

    elementCenter_.resize(Eigen::NoChange, nElements_);
    elementNormal_.resize(Eigen::NoChange, nElements_);
    elementArea_.resize(nElements_);

    int j = 0;
    for (size_t i1 = 0; i1 < nPatches_; ++i1) {
        for (int i2 = 0; i2 < n; ++i2) {
            s.y = h * i2;
            for (int i3=0; i3 < n; ++i3) {
                s.x = h * i3;
                for (size_t k = 0; k < Q[quadLevel].nop; k++) {
                    t = vector2_add(s,vector2_Smul(h,Q[quadLevel].xi[k]));
                    point = Chi_pwl(t,T_[i1], nLevels_);
                    norm = n_Chi_pwl(t,T_[i1], nLevels_);
                    Eigen::Vector3d center(point.x, point.y, point.z);
                    Eigen::Vector3d normal(norm.x,  norm.y,  norm.z);
                    normal.normalize();
                    double area = h * h * Q[quadLevel].w[k] * vector3_norm(n_Chi_pwl(t, T_[i1],
                                  nLevels_));
                    elementCenter_.col(j) = center.transpose();
                    elementNormal_.col(j) = normal.transpose();
                    elementArea_(j) = area;
                    ++j;
                }
            }
        }
    }
    free_Gauss_Square(&Q,quadLevel+1);
    built = true;
}

std::ostream & WaveletCavity::printCavity(std::ostream & os)
{
    os << "Cavity type: Wavelet" << std::endl;
    os << "Probe Radius =  " << probeRadius_ << std::endl;
    os << "Coarsity     =  " << coarsity_ << std::endl;
    os << "Patch Level  =  " << patchLevel_ << std::endl;
    os << "Number of spheres = " << nSpheres_ << std::endl;
    os << "Number of finite elements = " << nElements_;
    /*for(int i = 0; i < nElements_; i++)
    {
    	os << std::endl;
    	os << i+1 << " ";
    	os << elementCenter_(0,i) << " ";
    	os << elementCenter_(1,i) << " ";
    	os << elementCenter_(2,i) << " ";
    	os << elementArea_(i) << " ";
           }
       	for(int i = 0; i < nSpheres_; i++)
    {
    	os << endl;
    	os << i+1 << " ";
    	os << sphereCenter_(0,i) << " ";
    	os << sphereCenter_(1,i) << " ";
    	os << sphereCenter_(2,i) << " ";
    	os << sphereRadius_(i) << " ";
       	}*/
    /*	if (uploadedDyadic_)
    	{
    		os << "Printing nodes" << std::endl;
    		for(int i = 0; i < nPoints_; i++)
    		{
    			os << std::endl;
    			os << i+1 << " ";
    			os << nodeIndex_[i].transpose() << " " << nodePoint_[i].transpose() << " ";
    		}
    	}*/
    return os;
}

