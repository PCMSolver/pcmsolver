#include <iostream>
#include <fstream> 
#include <string>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Atom.h"
#include "Sphere.h"
#include "Cavity.h"
#include "CavityOfSpheres.h"

void GePolCavity::setMode(const string & type) {
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

void GePolCavity::setMode(int type) {
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
