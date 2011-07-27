#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "Cavity.h"

/*

Methods for basic cavity class
written by Krzysztof Mozgawa, 2011

*/

void Cavity::writeOutput(string &filename){
    ofstream output;
    output.open(filename.c_str(), fstream::out);
    output << nTess << endl;
    for(int i=0; i < nTess; i++) {
		output << tessCenter(i,0) << " ";
		output << tessCenter(i,1) << " ";
		output << tessCenter(i,2) << " ";
		output << tessArea(i) << " ";
    }
    output.close();
}

ostream & operator<<(ostream &os, const Cavity &cavity) {
	os << "Molecular cavity" << endl;
	os << "Nr. of tesserae: " << cavity.nTess;
    for(int i = 0; i < cavity.nTess; i++) {
		os << endl;
		os << i+1 << " ";
		os << cavity.tessCenter(i,0) << " ";
		os << cavity.tessCenter(i,1) << " ";
		os << cavity.tessCenter(i,2) << " ";
		os << cavity.tessArea(i);
    }
	return os;
}
