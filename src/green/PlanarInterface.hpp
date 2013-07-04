#ifndef PLANARINTERFACE_HPP
#define PLANARINTERFACE_HPP

#include "Config.hpp"

#include "GreensFunction.hpp"

class PlanarInterface : public GreensFunction
{
	public:
		PlanarInterface(double eps1, double eps2, double pos, double width);
		~PlanarInterface(){};
    		double evalf(double* p1, double* p2);
	
	private:
		double eps1;
		double eps2;
		double pos;
		double width;
		bool computed;
};

#endif // PLANARINTERFACE_HPP
