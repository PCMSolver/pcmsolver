#include <iostream>
#include <fstream> 
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "CavityOfSpheres.h"

void CavityOfSpheres::setMode(const string & type) {
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

void CavityOfSpheres::setMode(int type) {
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
