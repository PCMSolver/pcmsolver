#ifndef INPUT_H
#define INPUT_H

#include "Getkw.h"
#include "Sphere.h"

/**
 * The Input class is a wrapper class for the Getkw Library C++ bindings. 
 * It is implemented as a Singleton: only one instance of it is needed 
 * and therefore admitted (see Alexandrescu "Modern C++ Design").
 * 
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
		Input();
		Input(const char *parsedInputFile);
		/// Copy constructor is made private
		Input(const Input &other);
		Input& operator=(const Input &other);
        	~Input(){};
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
		double correction;
		double probeRadius;
		std::string greenType;
		std::string derivativeType;
		double epsilon;
		double epsilonReal;
		double epsilonImaginary;
		std::vector<double> spherePosition;
		double sphereRadius;
		
	public:
		static Input& CreateInput(const char *parsedInputFile) 
		{
			static Input obj(parsedInputFile);
			return obj;
		}
		// Accessor methods
		// Cavity section input
		std::string getType(){ return this->type; };
		int getPatchLevel(){ return this->patchLevel; };
		double getCoarsity(){ return this->coarsity; };
		double getArea(){ return this->area; };
		bool getScaling(){ return this->scaling; }; 
		bool getAddSpheres(){ return this->addSpheres; };
		std::string getMode(){ return this->mode; }; // vector<Sphere> should be setup here. Could this work in Implicit mode?
		std::vector<int> getAtoms(){ return this->atoms; };
		std::vector<double> getRadii(){ return this->radii; };
		std::vector<Sphere> getSpheres(){ return this->spheres; };
		// Medium section input
		std::string getSolvent(){ return this->solvent; };
		std::string getSolverType(){ return this->solverType; };
		std::string getEquationType(){ return this->equationType; };
		double getCorrection(){ return this->correction; };
		double getProbeRadius(){ return this->probeRadius; };
		std::string getGreenType(){ return this->greenType; };
		std::string getDerivativeType(){ return this->derivativeType; };
		double getEpsilon(){ return this->epsilon; }; 
		double getEpsilonReal(){ return this->epsilonReal; };
		double getEpsilonImaginary(){ return this->epsilonImaginary; };
		std::vector<double> getSpherePosition(){ return this->spherePosition; };
		double getSphereRadius(){ return this->sphereRadius; };
		/// Operators
		/// operator<<
                friend std::ostream & operator<<(std::ostream &os, const Input &input);
                /// @}
};

#endif // INPUT_H
