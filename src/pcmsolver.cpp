#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "TraditionalPCMOperator.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "MetalSphere.h"
#include "GreensFunctionSum.h"

#include "PCMSolver.h"

int main(){
    VectorXd potential, charges;
    std::string fname("cavity.out");
    
    Vector3d p1(0.0, 10.1, 0.0);
    Vector3d p2(0.0, 10.2, 0.0);
    Vector3d ps(0.0,  0.0, 0.0);
    
    UniformDielectric water(78.39);
    UniformDielectric cyclohexane(2.0);
    Vacuum vacuum;
    
    PCMSolver<Vacuum, UniformDielectric> waterSolver(vacuum, water);
    GreensFunction &water2 = waterSolver.getGreenOutside();
    double green = water2.evalf(p1,p2);
    cout << " green " << green << endl;
    waterSolver.readCavity(fname);
    waterSolver.buildPCMMatrix();
}

//
//  atomicinput.open("inp_atoms", fstream::in);
//  if (atomicinput.eof()){
//    cout << "Unexpected end of file." ;
//    exit(1);
//  }
//  if (atomicinput.bad() != 0 ){
//    cout << "Type mismatch or file corrupted." ;
//    exit(1);
//  }
//  atomicinput >> numberOfAtoms;
//  for(int i = 0; i < numberOfAtoms; i++) {
//    if (atomicinput.eof()){
//      cout << "Unexpected end of file." ;
//      exit(1);
//    }
//    if (atomicinput.bad() != 0 ){
//      cout << "Type mismatch or file corrupted." ;
//      exit(1);
//    }
//
//    atomicinput >> Rwater(0);
//    atomicinput >> Rwater(1);
//    atomicinput >> Rwater(2);
//    atomicinput >> atomicCharges(i);
//
//    atoms.push_back(Rwater);
//
//  }
//
//  atomicinput.close();
//  for(unsigned int i = 0; i < points.size(); i++){
//    for(int j = 0; j < numberOfAtoms; j++){
//      potential[i] = potential[i] + atomicCharges[j]/(points[i] - atoms[j]).norm();
//    }
//  }
//
//  op.solveForPotential(potential, charges);
//
//  
//  for(unsigned int i = 0; i < points.size(); i++){
//    output << points[i] ;
//    output << " " << charges[i] << endl;
//  }
//  
//  /*
//  for(unsigned int i = 0; i < points.size(); i++){
//    cout << points[i];
//    cout << charges[i] << endl;
//  }
//  */
//
//  double charge2 = 0.0;
//  for(unsigned int i = 0; i < points.size(); i++){
//    totalCharge = totalCharge + charges[i];
//    charge2 += charges[i]*areas[i];
//  }
//
//  /*
//  for(unsigned int i = 0; i < points.size(); i++){
//    cout << "nr chg pot: " << i << " " << charges[i] << " " << potential[i] << endl;
//  }
//  */
//  
//  cout << "Villes value: " << 9.748 << " Calculated total charge: " << totalCharge << endl;
//  cout << "Calculated total charge*areas: " << charge2 << endl;
//
//  output.close();
//
//
//}
//
