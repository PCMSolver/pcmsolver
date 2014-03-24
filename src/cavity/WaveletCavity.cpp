#include "WaveletCavity.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

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
    uploadedDyadic = false;
}

extern "C"
{
    int waveletCavityDrv_(double probeRadius, double coarsity, int patchLevel,
                          const char * infile);
}

void WaveletCavity::makeCavity()
{
    int check = 0;
    std::string infile = "cavity.inp";
    writeInput(infile);
#if defined (WAVELET_DEVELOPMENT)
    check = waveletCavityDrv_(probeRadius, coarsity, patchLevel, infile.c_str());
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
    file >> nLevels >> nPatches;

    int nNodes = (1 << nLevels) + 1;

    nPoints = nPatches * nNodes * nNodes;

    for (size_t k = 0; k < nPoints; ++k) {
        file >> i >> j >> k >> x >> y >> z;
        Eigen::Vector3i index(i, j, k);
        nodeIndex.push_back(index);
        Eigen::Vector3d point(x, y, z);
        nodePoint.push_back(point);
    }

    file.close();
    uploadedDyadic = true;
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
    if (not uploadedDyadic) {
        throw std::runtime_error("Dyadic file must be uploaded first.");
    }
    vector2 s, t;
    vector3 point;
    vector3 norm;
    int n = 1 << nLevels;
    double h = 1.0 / n;
    cubature *Q;
    init_Gauss_Square(&Q, quadLevel + 1);

    nElements_ = nPatches * n * n * Q[quadLevel].nop;

    elementCenter_.resize(Eigen::NoChange, nElements_);
    elementNormal_.resize(Eigen::NoChange, nElements_);
    elementArea_.resize(nElements_);
    // The following is to guarantee that writing cavity to .npz file works
    elementRadius_.setZero(nElements_);

    int j = 0;
    for (size_t i1 = 0; i1 < nPatches; ++i1) {
        for (int i2 = 0; i2 < n; ++i2) {
            s.y = h * i2;
            for (int i3=0; i3 < n; ++i3) {
                s.x = h * i3;
                for (size_t k = 0; k < Q[quadLevel].nop; ++k) {
                    t = vector2_add(s,vector2_Smul(h,Q[quadLevel].xi[k]));
                    point = Chi(t,T_[i1], nLevels);
                    norm = n_Chi(t,T_[i1], nLevels);
                    Eigen::Vector3d center(point.x, point.y, point.z);
                    Eigen::Vector3d normal(norm.x,  norm.y,  norm.z);
                    normal.normalize();
                    double area = h * h * Q[quadLevel].w[k] * vector3_norm(n_Chi(t, T_[i1], nLevels));
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
    if (not uploadedDyadic) {
        throw std::runtime_error("Dyadic file must be uploaded first.");
    }
    vector2 s, t;
    vector3 point;
    vector3 norm;
    int n = 1 << nLevels;
    double h = 1.0 / n;
    cubature *Q;
    init_Gauss_Square(&Q, quadLevel + 1);

    nElements_ = nPatches * n * n * Q[quadLevel].nop;

    elementCenter_.resize(Eigen::NoChange, nElements_);
    elementNormal_.resize(Eigen::NoChange, nElements_);
    elementArea_.resize(nElements_);

    int j = 0;
    for (size_t i1 = 0; i1 < nPatches; ++i1) {
        for (int i2 = 0; i2 < n; ++i2) {
            s.y = h * i2;
            for (int i3=0; i3 < n; ++i3) {
                s.x = h * i3;
                for (size_t k = 0; k < Q[quadLevel].nop; k++) {
                    t = vector2_add(s,vector2_Smul(h,Q[quadLevel].xi[k]));
                    point = Chi_pwl(t,T_[i1], nLevels);
                    norm = n_Chi_pwl(t,T_[i1], nLevels);
                    Eigen::Vector3d center(point.x, point.y, point.z);
                    Eigen::Vector3d normal(norm.x,  norm.y,  norm.z);
                    normal.normalize();
                    double area = h * h * Q[quadLevel].w[k] * vector3_norm(n_Chi_pwl(t, T_[i1],
                                  nLevels));
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
    os << "Probe Radius =  " << probeRadius << std::endl;
    os << "Coarsity =      " << coarsity << std::endl;
    os << "Patch Level =   " << patchLevel << std::endl;
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
    /*	if (uploadedDyadic)
    	{
    		os << "Printing nodes" << std::endl;
    		for(int i = 0; i < nPoints; i++)
    		{
    			os << std::endl;
    			os << i+1 << " ";
    			os << nodeIndex[i].transpose() << " " << nodePoint[i].transpose() << " ";
    		}
    	}*/
    return os;
}

