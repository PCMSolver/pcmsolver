#ifndef INPUT_H
#define INPUT_H

#include "Getkw.h"
#include "Sphere.h"
#include "Solvent.h"

/*!
 * \brief A wrapper class for the Getkw Library C++ bindings.
 * \author Roberto Di Remigio
 * \date 2013
 *
 * Implemented as a Singleton: only one instance of it is needed and therefore admitted.
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
		//Input(const char *parsedInputFile);
		/// Copy constructor is made private
		Input(const Input &other);
		Input& operator=(const Input &other);
        	~Input(){}
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
		Solvent solvent;
		std::string solverType;
		std::string equationType;
		double correction;
		double probeRadius;
		std::string greenInsideType;
		std::string greenOutsideType;
		std::string derivativeInsideType;
		std::string derivativeOutsideType;
		double epsilonInside;
		double epsilonOutside;
		double epsilonReal;
		double epsilonImaginary;
		std::vector<double> spherePosition;
		double sphereRadius;
		
	public:
		//static Input& TheInput(const char *parsedInputFile) 
		static Input& TheInput() 
		{
			static Input obj;
			return obj;
		}
		// Accessor methods
		// Cavity section input
		std::string getCavityType(){ return type; }
		int getPatchLevel(){ return patchLevel; }
		double getCoarsity(){ return coarsity; }
		double getArea(){ return area; }
		bool getScaling(){ return scaling; } 
		bool getAddSpheres(){ return addSpheres; }
		std::string getMode(){ return mode; }
		std::vector<int> getAtoms(){ return atoms; }
		std::vector<double> getRadii(){ return radii; }
		std::vector<Sphere> getSpheres(){ return spheres; }
		void setSpheres(const std::vector<Sphere> & _spheres){ spheres = _spheres; }
		// Medium section input
		Solvent getSolvent(){ return solvent; }
		std::string getSolverType(){ return solverType; }
		std::string getEquationType(){ return equationType; }
		double getCorrection(){ return correction; }
		double getProbeRadius(){ return probeRadius; }
		std::string getGreenInsideType(){ return greenInsideType; }
		std::string getGreenOutsideType(){ return greenOutsideType; }
		std::string getDerivativeInsideType(){ return derivativeInsideType; }
		std::string getDerivativeOutsideType(){ return derivativeOutsideType; }
		double getEpsilonInside(){ return epsilonInside; } 
		double getEpsilonOutside(){ return epsilonOutside; }
		double getEpsilonReal(){ return epsilonReal; }
		double getEpsilonImaginary(){ return epsilonImaginary; }
		std::vector<double> getSpherePosition(){ return spherePosition; }
		double getSphereRadius(){ return sphereRadius; }
		/// Operators
		/// operator<<
                friend std::ostream & operator<<(std::ostream &os, const Input &input);
                /// @}
};

#endif // INPUT_H
