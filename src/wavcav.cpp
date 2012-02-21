#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "WaveletCavity.h"

int main(int argc, char** argv){

	const char *infile = 0;
	if (argc == 1) {
		infile = "STDIN";
	} else if (argc == 2) {
		infile = argv[1];
	} else {
		cout << "Invalid nr. of arguments" << endl;
		exit(1);
	}
	
	Getkw Input = Getkw(infile, false, true);
	const Section & WaveletCavitySection = Input.getSect("Cavity<wavelet>");
    WaveletCavity wavcav(WaveletCavitySection);
	wavcav.makeCavity();
}
