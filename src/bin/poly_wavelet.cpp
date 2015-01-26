
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "PWCSolver.hpp"
#include "PWLSolver.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "WaveletCavity.hpp"
#include "PhysicalConstants.hpp"

#include "LoggerInterface.hpp"

void read_data(const std::string &filename, Compression *comp,  std::vector<double> *charge_, std::vector<Eigen::Vector3d> *atoms_){

    unsigned int nAtoms;
    double charge;
    double x,y,z;
    double a,dp,b, eps;

    std::ifstream file;
    file.open(filename);

    if(file.is_open()) {
        LOG(">>> DATA_FILE");
        file >> a >> dp >> b >> eps;
        LOG(a, " ", dp, " ", b, " ", eps);
    	comp->aPrioriA = a;
    	comp->aPrioridPrime = dp;
    	comp->aPosterioriB = b;
        file >> nAtoms;
        for(unsigned int i = 0; i < nAtoms; ++i){
            file >> charge >> x >> y >> z;
            LOG(charge, " ", x," ", y, " ", z);
            Eigen::Vector3d position(x,y,z);
            atoms_->push_back(position);
            charge_->push_back(charge);
        }
        LOG("<<< DATA_FILE");
    } else {
        throw std::runtime_error("Data file could not be opened");
    }
}



int main(int argc, char* argv[]) {
    //read_sphere();
    Compression comp;
    std::vector<double> charge_;
    std::vector<Eigen::Vector3d> atoms_;
    read_data(argv[2], &comp, &charge_, &atoms_);
    
    LOG(">>> GEOMETRY_FILE");
    LOG(argv[1]);
    LOG("<<< GEOMETRY_FILE");
    // Set up cavity, read it from Maharavo's file 
    WaveletCavity cavity(argv[1]);
    //cavity.scaleCavity(1./convertBohrToAngstrom);

    double permittivity = 78.39;
    double totalASC = 0;
    for(unsigned int i = 0; i < charge_.size(); ++i){
        totalASC -= charge_[i];
    }
    totalASC *= (permittivity-1)/permittivity; 

    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>();
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity);
    int firstKind = 0;
#ifdef DEBUG2
    FILE* debugFile = fopen("debug.out","w");
    fclose(debugFile);
#endif
    LOG(">>> PCMSOLVER");
    PWCSolver solver(gfInside, gfOutside, comp, firstKind);
    solver.buildSystemMatrix(cavity);
    cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_());

    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        fake_mep(i) = 0.0;
        Eigen::Vector3d center = cavity.elementCenter(i);
        for( unsigned int j = 0; j < atoms_.size(); ++j){
            fake_mep(i) += charge_[j]/(center - atoms_[j]).norm();
        }
    }
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalFakeASC = fake_asc.sum();

    double energy = 0.5 * (fake_mep.dot(fake_asc));
    energy *= 0.5291772;
   
    std::stringstream solverString;
    solverString << solver; 
    LOG(solverString.str());
    LOG(cavity);
    LOG("totalASC     = ", totalASC);
    LOG("totalFakeASC = ", totalFakeASC);
    LOG("Delta        = ", totalASC - totalFakeASC);
    LOG("Energy       = ", energy);
    LOG("<<< PCMSOLVER");
    LOG(">>> TIMING");
    LOG_TIME;
    LOG("<<< TIMING");
    LOG("# vim: foldmarker=>>>,<<< foldlevel=0 foldmethod=marker");
}
