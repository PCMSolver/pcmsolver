#include <string>
#include <iostream>
#include <fstream>

using namespace std;


#include "TraditionalPCMOperator.h"
#include "GreensFunction.h"
#include "UniformDielectric.h"
#include "MetalSphere.h"
#include "GreensFunctionSum.h"

#include "PCMSolver.h"

int main(){
  std::vector<Vector3d> points, normals, atoms;
  std::vector<double> areas;
  VectorXd potential, charges;
  TraditionalPCMOperator op(78.39);
  //  GreensFunctionInterface iface;
  std::string fname("output_cav");
  //  VectorXd *potential, *charges;



  // potential testing variables
  Vector3d Rwater;
  int numberOfAtoms;
  double totalCharge;
  ifstream atomicinput;
  Vector3d atomicCharges;
  ofstream output;
  //  output.open("output_solver", fstream::out);

  //  op.readCavity(fname);
  //  op.getCavityDefs(points, normals, areas);
  //  op.test();
  //  op.constructSystemMatrix();

  //  potential = new VectorXd(points.size(),1);
  //  charges = new VectorXd(points.size(),1);

    // here I nastily hard code coordinates for water, just for a test case. To be edited in future. KM

  //  potential.setZero(points.size());
  //  charges.setZero(points.size());

  totalCharge = 0;

  double p1[3] = {0.0, 10.1, 0.0};
  double p2[3] = {0.0, 10.2, 0.0};
  double ps[3] = {0.0, 0.0, 0.0};

  UniformDielectric water(78.39);
  MetalSphere metal(78.39, -1000000.0, 0.0, ps, 10.0);
  GreensFunctionSum waterMetal(water, metal);

  PCMSolver<UniformDielectric> waterSolver(water);
  PCMSolver<MetalSphere> metalSolver(metal);
  PCMSolver<GreensFunctionSum> waterMetalSolver(waterMetal);
  
  GreensFunction &water2 = waterSolver.getGreensFunction();
  GreensFunction &metal2 = metalSolver.getGreensFunction();
  GreensFunction &waterMetal2 = waterMetalSolver.getGreensFunction();

  double green = water2.evalf(p1,p2);
  cout << " green " << green << endl;
  double green2 = metal2.evalf(p1,p2);
  cout << " green2 " << green2 << endl;
  double green3 = waterMetal2.evalf(p1,p2);
  cout << " green3 " << green3 << endl;

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
