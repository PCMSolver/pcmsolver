/*! \file TraditionalPCMOperator.cpp
\brief Class implementation for the traditional pcm solver
*/

#include "TraditionalPCMOperator.h"

/*
  #include <iostream>
  #include <fstream>
  #include <ostream>
  #include <sstream> 
*/

/*
Traditional PCM Operator program written by Krzysztof Mozgawa and Ville Weijo, 2010
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/LU>


using namespace std;

bool TraditionalPCMOperator::readCavity(string &filename){
  
  /*
    Vector3d v;
    v(0) = x;
    v(1) = y;
    v(2) = z;

    points.push_back(v);
  */

  Vector3d v;
  Vector3d sph;
  double area;
  double radius;
  int TesNum;
  ifstream input;
  input.open(filename.c_str(), fstream::in);
  if (input.eof()){
    cout << "Unexpected end of file." ;
    exit(1);
  }
  if (input.bad() != 0 ){
    cout << "Type mismatch or file corrupted." ;
    exit(1);
  }
  input >> TesNum;
  for(int i = 0; i < TesNum; i++) {
    if (input.eof()){
      cout << "Unexpected end of file." ;
      exit(1);
    }
    if (input.bad() != 0 ){
      cout << "Type mismatch or file corrupted." ;
      exit(1);
    }
    input >> v(0);
    input >> v(1);
    input >> v(2);
    input >> area;
    input >> sph(0);
    input >> sph(1);
    input >> sph(2);
    input >> radius;
    points_.push_back(v);
    areas_.push_back(area);
    correspondingSpheres_.push_back(sph);
    sphereRadii_.push_back(radius*1.889726127);
    }
  /*
    double area;
    area = something;
    areas.push_back(area);
  */
  input.close();
  return false;
}

void TraditionalPCMOperator::getCavityDefs(std::vector<Vector3d> &points, 
                                           std::vector<Vector3d> &normals,
					   std::vector<double> &areas){

  Vector3d work;
  areas = areas_;
  points = points_;
  for (unsigned int i=0; i< points_.size(); i++){
    work=points_[i] - correspondingSpheres_[i];
    work.normalize();
    normals.push_back(work);
  }
  normals_ = normals;
}

//void TraditionalPCMOperator::constructSystemMatrix(GreensFunctionInterface &iface){
void TraditionalPCMOperator::constructSystemMatrix(){

  /* THE OLD VERSION STARTS HERE

  MatrixXd S((int)points_.size(), (int)points_.size());
  MatrixXd D((int)points_.size(), (int)points_.size());
  MatrixXd Ainv((int)points_.size(), (int)points_.size());
  MatrixXd IEFPCM((int)points_.size(), (int)points_.size());  

  
  Vector3d loadingVector;

  S.setZero();
  D.setZero();
  Ainv.setZero();
  IEFPCM.setZero();
  
  for(unsigned int i = 0; i < points_.size(); i++){
    for (unsigned int j=0; j < points_.size(); j++){
	     if (i != j) {
	       S(i,j) = 1.0 / (points_[i] - points_[j]).norm();
	     }
	     else {
	       S(i,j) = 1.07 * sqrt((4*M_PI)/areas_[i]);
	     }
    }
  }   

  for(unsigned int i = 0; i < points_.size(); i++){
    for (unsigned int j=0; j < points_.size(); j++){
	     if (i != j) {
	       D(i,j) = (points_[i] - points_[j]).dot(normals_[j]);
	       D(i,j) = D(i,j) / pow((points_[i] - points_[j]).norm(),3);
	     }
	     else {
	       D(i,j) = 1.07 * sqrt(4*M_PI*areas_[i]);
	       D(i,j) = D(i,j) / (2.0 * sphereRadii_[i]);
	     }
    }
  } 
   
  for(unsigned int i = 0; i < points_.size(); i++){
    Ainv(i,i) = 1/areas_[i];
      }
 

  Ainv = 2*M_PI*Ainv;
  IEFPCM = Ainv*((PCM_E+1)/(PCM_E-1));
  IEFPCM = IEFPCM - D;
  IEFPCM = IEFPCM * S; 

  IEFPCM = IEFPCM.inverse();

  IEFPCM = IEFPCM * (Ainv - D);


  systemMatrix_ = IEFPCM;
  

  THE OLD VERSION ENDS HERE */


  MatrixXd S((int)points_.size(), (int)points_.size());
  MatrixXd D((int)points_.size(), (int)points_.size());
  MatrixXd IEFPCM((int)points_.size(), (int)points_.size());
  MatrixXd Ainv((int)points_.size(), (int)points_.size());

  S.setZero();
  D.setZero();
  Ainv.setZero();
  IEFPCM.setZero();

  for(unsigned int i = 0; i < points_.size(); i++){
    Ainv(i,i) = areas_[i];
  }

  for(unsigned int i = 0; i < points_.size(); i++){
    for (unsigned int j=0; j < points_.size(); j++){
      if (i != j) {
	S(i,j) = areas_[j] / (4*M_PI*(points_[i] - points_[j]).norm());
      }
      else {
	S(i,j) = 1.07 * sqrt(areas_[i]/(4*M_PI));
      }
    }
  }

  for(unsigned int i = 0; i < points_.size(); i++){
    for (unsigned int j=0; j < points_.size(); j++){
      if (i != j) {
	D(i,j) = areas_[i]*areas_[j]*(points_[i] - points_[j]).dot(normals_[j]);;
	D(i,j) = D(i,j) / (4*M_PI*pow((points_[i] - points_[j]).norm(),3));
      }
      else {
	D(i,j) = -1.07 * areas_[i] * sqrt(areas_[i]/(4*M_PI));
	D(i,j) = D(i,j) / (2.0 * sphereRadii_[i]);
      }
    }
  }

  /*

  for(unsigned int i = 0; i < points_.size(); i++){
    for (unsigned int j= 0; j < points_.size(); j++){
      output << i+1 << " " << j+1 << " " << D(i,j) << " " 
	     << S(i,j) << " " << points_[j] << " "
	     << areas_[i] << " " << areas_[j] << " "
	     << (points_[i] - points_[j]).norm() 
	     <<	" " << sphereRadii_[i] << endl;
    }
  }
  */
  /*
  for(unsigned int i = points_.size(); i > 0; i--){
    output << "Column: " << i << " ";
    for (unsigned int j= 0; j < points_.size(); j++){
      output << S(i-1,j) << " " ;
    }
    output << endl;
  }
  */

  
  /*
  for(unsigned int i = 0; i < points_.size(); i++){
    for(unsigned int k = 0; k < points_.size(); k++){
      for (unsigned int j=0; j < points_.size(); j++){
	if (i != j) {
	  IEFPCM(i,k) = IEFPCM(i,k) - 4 * M_PI * D(i,j) * S(j,k);
	}
	else {
	  IEFPCM(i,k) = IEFPCM(i,k) + 4 * M_PI * ( Ainv(i,i) * (((PCM_E+1)/(PCM_E-1))/2) - D(i,j))*S(j,k);
	}
      }
    }
  }

  */

  IEFPCM = 4.0*M_PI*(Ainv*(((PCM_E+1.0)/(PCM_E-1.0))/2.0) - D)*S;
  /*
  for(unsigned int i = points_.size(); i > 0; i--){
    cout << "Column: " << i << " ";
    for (unsigned int j= 0; j < points_.size(); j++){
      cout << j << " " << IEFPCM(i-1,j);
    }
    cout << endl;
  }

  */
  IEFPCM = IEFPCM.inverse()*(Ainv/2.0 - D);

  systemMatrix_ = IEFPCM;

  /*  
  for(unsigned int i = 0; i < points.size(); i++){
    for (unsigned int j=0; j < points.size(); j++){
      cout << "next value:" << IEFPCM(i,j);
    }
    cout << endl;
  }
  IEFPCMtest = IEFPCMtest*IEFPCM;

  for(unsigned int i = 0; i < points.size(); i++){
    for (unsigned int j=0; j < points.size(); j++){
      cout << "next value:" << IEFPCMtest(i,j);
    }
    cout << endl;
  }
  */

  
}

void TraditionalPCMOperator::solveForPotential(VectorXd &potential, VectorXd &charges){

  //void TraditionalPCMOperator::solveForPotential(){


  //  VectorXd Q((int)points_.size(), 1);


  charges = systemMatrix_ * potential;



  //  Q = potential * systemMatrix_;

  //  charges = Q;



}

void TraditionalPCMOperator::printInfo(std::ostream &out){

  out << "Used epsilon value: "<< PCM_E << endl << "Number of tesserae " 
      << points_.size() << endl; 

}
/*
BEMInterface *TraditionalPCMOperator::newBEMInterface(){

  return new TraditionalPCMOperator(PCM_E);

}
*/
/*	     		     

void TraditionalPCMOperator::test(){
  cout << "Number of points: " << points_.size() << endl;
  for(unsigned int i = 0; i < points_.size(); i++)
    cout << points_[i].transpose() << endl;
}
*/
