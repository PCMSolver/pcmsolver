#ifndef WAVELETCAVITY
#define WAVELETCAVITY

#include <iostream>
#include <string>

#include "Config.h"

#include <Eigen/Dense>

#include "vector3.h"
#include "Getkw.h"
#include "Cavity.h"

class WaveletCavity : public Cavity {
 public:
    WaveletCavity(){}
    //    WaveletCavity(string &filename);
    WaveletCavity(const Getkw & input, const string path = "Cavity");
    WaveletCavity(const Section & cavity);
    WaveletCavity(int _patchLevel, const std::vector<Sphere> & _spheres, double _coarsity, double _probeRadius) : 
     patchLevel(_patchLevel), coarsity(_coarsity), probeRadius(_probeRadius) 
           {
		// Initialize data members inherited from CavityOfSpheres...
		spheres = _spheres;
                nSpheres = spheres.size();
		sphereCenter.resize(Eigen::NoChange, nSpheres);
		sphereRadius.resize(nSpheres);
		for (int i = 0; i < nSpheres; ++i) {
			sphereCenter.col(i) = spheres[i].getSphereCenter();
			sphereRadius(i) = spheres[i].getSphereRadius();
		}
		uploadedDyadic = false;
		// ...and build the cavity!
		makeCavity();
           }
    ~WaveletCavity(){};
    void makeCavity();
    void readCavity(const string & filename);
    void uploadPoints(int quadLevel, vector3 **** T_, bool isPWL);
//    VectorXd & getTessRadius(){return tessRadius;};
//    VectorXd & getSphereRadius(){return sphereRadius;};
//    int getNSpheres(){return nSpheres;};
//    Matrix3Xd & getSphereCenter(){return sphereCenter;};
//    Matrix3Xd & getTessSphereCenter(){return tessSphereCenter;};
//    double getTessRadius(int i){return tessRadius(i);};
    unsigned int getNPatches(){return nPatches;}
    unsigned int getNLevels(){return nLevels;}
    unsigned int getNPoints(){return nPoints;}
    Vector3d getNodePoint(int i){return nodePoint[i];}
    Vector3i getNodeIndex(int i){return nodeIndex[i];}
    friend std::ostream& operator<<(std::ostream &o, const WaveletCavity &c);
    void compFakePotential();
 private:
    void uploadPointsPWC(int quadLevel, vector3 **** T_);
    void uploadPointsPWL(int quadLevel, vector3 **** T_);
    std::vector<Vector3d> nodePoint;
    std::vector<Vector3i> nodeIndex;
    unsigned int nPatches;
    unsigned int nLevels;
    unsigned int nPoints;
    bool uploadedDyadic;
    void writeInput(string &fileName);
 //   int nSpheres;
 //   Matrix3Xd sphereCenter;
 //   Matrix3Xd tessSphereCenter;
 //   VectorXd sphereRadius;
 //   VectorXd tessRadius;
    int patchLevel;
    double probeRadius;
    double coarsity;
};

#endif
