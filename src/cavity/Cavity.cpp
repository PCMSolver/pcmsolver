#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Cavity.h"

/*

Methods for basic cavity class
written by Krzysztof Mozgawa, 2011

*/


void Cavity::writeOutput(string &filename){
    ofstream output;
    output.open(filename.c_str(), fstream::out);
    output << nElements << endl;
    for(int i=0; i < nElements; i++) {
		output << elementCenter(0,i) << " ";
		output << elementCenter(1,i) << " ";
		output << elementCenter(2,i) << " ";
		output << elementArea(i) << " ";
    }
    output.close();
}

ostream & operator<<(ostream & os, Cavity & cavity) {
	return cavity.printObject(os);
}

ostream & Cavity::printObject(ostream & os) {
	os << "Molecular cavity" << endl;
	os << "Nr. of tesserae: " << nElements;
    for(int i = 0; i < nElements; i++) {
		os << endl;
		os << i+1 << " ";
		os << elementCenter(0,i) << " ";
		os << elementCenter(1,i) << " ";
		os << elementCenter(2,i) << " ";
		os << elementArea(i) << " ";
    }
	return os;
}

void Cavity::setMode(const string & type) {
	if (type == "Atoms") {
		setMode(Atoms);
	} else if (type == "Implicit") {
		setMode(Implicit);
	} else if (type == "Explicit") {
		setMode(Explicit);
	} else {
		exit(-1);
	}
}

void Cavity::setMode(int type) {
	switch (type) {
	case Atoms :
		mode = Atoms;
		break;
	case Implicit :
		mode = Implicit;
		break;
	case Explicit :
		mode = Explicit;
		break;
	default :
		exit(-1);
	}
}
