#include "WaveletCavity.hpp"

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

        output << nSpheres << std::endl;
	output.setf(std::ios_base::showpoint);
	output.precision(12);
        for(int i=0; i < nSpheres; i++) 
	{
		output << sphereCenter(0,i) << " ";
		output << sphereCenter(1,i) << " ";
		output << sphereCenter(2,i) << " ";
		output << sphereRadius(i) << std::endl;
        }
        output.close();
	uploadedDyadic = false;
}

extern "C" 
{
	int waveletCavityDrv_(double probeRadius, double coarsity, int patchLevel, const char * infile);
}

void WaveletCavity::makeCavity() 
{
//	int dummy = 0; 
	int check = 0;
	std::string infile = "cavity.inp";
	writeInput(infile);
	check = waveletCavityDrv_(probeRadius, coarsity, patchLevel, infile.c_str());
	if (check != 0) 
	{
		throw std::runtime_error("Problem with the wavelet cavity inside makeCavity method");
	}
}


void WaveletCavity::readCavity(const std::string & filename) 
{
	int i, j, k;
	double x, y, z;

	std::ifstream file;
	file.open(filename.c_str());
	file >> nLevels >> nPatches;

	int nNodes = (1 << nLevels) + 1;

	nPoints = nPatches * nNodes * nNodes;
	
	for (int k = 0; k < nPoints; ++k) 
	{
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
	if(isPWL) 
	{
		uploadPointsPWL(quadLevel, T_);
	} 
	else 
	{
		uploadPointsPWC(quadLevel, T_);
	}
}

void WaveletCavity::uploadPointsPWC(int quadLevel, vector3 **** T_) 
{
	if (not uploadedDyadic) 
	{
		throw std::runtime_error("Dyadic file must be uploaded first.");
	}
	vector2 s, t;
	vector3 point;
	vector3 norm;
	int n = 1 << nLevels;
	double h = 1.0 / n;
	cubature *Q;
	init_Gauss_Square(&Q, quadLevel + 1);

	nElements = nPatches * n * n * Q[quadLevel].nop;

	elementCenter.resize(Eigen::NoChange, nElements);
	elementNormal.resize(Eigen::NoChange, nElements);
	elementArea.resize(nElements);

	int j = 0;
	for (int i1 = 0; i1 < nPatches; ++i1)
	{
		for (int i2 = 0; i2 < n; ++i2)
		{
			s.y = h * i2;
			for (int i3=0; i3 < n; ++i3)
			{
				s.x = h * i3;
				for (int k = 0; k < Q[quadLevel].nop; ++k)
				{
					t = vector2_add(s,vector2_Smul(h,Q[quadLevel].xi[k]));
					point = Chi(t,T_[i1], nLevels);
					norm = n_Chi(t,T_[i1], nLevels);
					Eigen::Vector3d center(point.x, point.y, point.z);	 
					Eigen::Vector3d normal(norm.x,  norm.y,  norm.z);	 
					normal.normalize();
					double area = h * h * Q[quadLevel].w[k] * vector3_norm(n_Chi(t, T_[i1], nLevels));
					elementCenter.col(j) = center.transpose();
					elementNormal.col(j) = normal.transpose();
					elementArea(j) = area;
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
	if (not uploadedDyadic) 
	{
		throw std::runtime_error("Dyadic file must be uploaded first.");
	}
	vector2 s, t;
	vector3 point;
	vector3 norm;
	int n = 1 << nLevels;
	double h = 1.0 / n;
	cubature *Q;
	init_Gauss_Square(&Q, quadLevel + 1);

	nElements = nPatches * n * n * Q[quadLevel].nop;

	elementCenter.resize(Eigen::NoChange, nElements);
	elementNormal.resize(Eigen::NoChange, nElements);
	elementArea.resize(nElements);

	int j = 0;
	for (int i1 = 0; i1 < nPatches; ++i1)
	{
		for (int i2 = 0; i2 < n; ++i2)
		{
			s.y = h * i2;
			for (int i3=0; i3 < n; ++i3)
			{
				s.x = h * i3;
				for (int k = 0; k < Q[quadLevel].nop; k++)
				{
					t = vector2_add(s,vector2_Smul(h,Q[quadLevel].xi[k]));
					point = Chi_pwl(t,T_[i1], nLevels);
					norm = n_Chi_pwl(t,T_[i1], nLevels);
					Eigen::Vector3d center(point.x, point.y, point.z);	 
					Eigen::Vector3d normal(norm.x,  norm.y,  norm.z);	 
					normal.normalize();
					double area = h * h * Q[quadLevel].w[k] * vector3_norm(n_Chi_pwl(t, T_[i1], nLevels));
					elementCenter.col(j) = center.transpose();
					elementNormal.col(j) = normal.transpose();
					elementArea(j) = area;
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
	os << "Number of spheres = " << nSpheres << std::endl;
        os << "Number of finite elements = " << nElements << std::endl;
        /*for(int i = 0; i < nElements; i++) 
	{
		os << std::endl;
		os << i+1 << " ";
		os << elementCenter(0,i) << " ";
		os << elementCenter(1,i) << " ";
		os << elementCenter(2,i) << " ";
		os << elementArea(i) << " ";
        }
    	for(int i = 0; i < nSpheres; i++) 
	{
		os << endl;
		os << i+1 << " ";
		os << sphereCenter(0,i) << " ";
		os << sphereCenter(1,i) << " ";
		os << sphereCenter(2,i) << " ";
		os << sphereRadius(i) << " ";
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

