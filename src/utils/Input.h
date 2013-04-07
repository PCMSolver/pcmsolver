#ifndef INPUT_H
#define INPUT_H

#include "Getkw.h"
#include "Sphere.h"

/**
 * The Input class is a wrapper class for the Getkw Library C++ bindings. 
 * It is to be used as the unique point of access to user-provided input:
 *  input ---> parsed input (Python script) ---> Input object (contains all the input data)
 * Definition of input parameters is to be done in the Python script and in this class.
 *
 * At runtime a unique Input object in created, this object contains all input parameters
 * processed as to be directly usable in the module.
 * 
 * They must be specified as private data members with public accessor methods (get-ters).
 * In general, no mutator methods (set-ters) should be needed, exceptions to this rule
 * should be carefully considered.
 */

class Input 
{
	private:
		std::string type;
		int patchLevel;
		double coarsity;
		double area;
		bool scaling;
		bool addSpheres;
		std::string mode;
		std::vector<int> atoms;
		std::vector<double> radii;
		std::vector<Sphere> spheres;
		std::string solvent;
		std::string solverType;
		std::string equationType;
		double probeRadius;
		std::string greenType;
		std::string derivativeType;
		double epsilon;
		double epsilonReal;
		double epsilonImaginary;
		std::vector<double> spherePosition;
		double sphereRadius;
		
	public:
		Input();
	        Input(const char *parsedInputFile);
		/// Copy constructor
		Input(const Input &other);
        	~Input(){};
		// Accessor methods
		// Cavity section input
		std::string getType();
		int getPatchLevel();
		double getCoarsity();
		double getArea();
		bool getScaling(); 
		bool getAddSpheres();
		std::string getMode(); // vector<Sphere> should be setup here. Could this work in Implicit mode?
		std::vector<int> getAtoms();
		std::vector<double> getRadii();
		std::vector<Sphere> getSpheres();
		// Medium section input
		std::string getSolvent();
		std::string getSolverType();
		std::string getEquationType();
		double getProbeRadius();
		std::string getGreenType();
		std::string getDerivativeType();
		double getEpsilon(); 
		double getEpsilonReal();
		double getEpsilonImaginary();
		std::vector<double> getSpherePosition();
		double getSphereRadius();
		/// Operators
		/// Assignment operator.
		Input & operator=(const Input &other);
		/// operator<<
                friend std::ostream & operator<<(std::ostream &os, const Input &input);
                /// @}
};

#endif // INPUT_H
