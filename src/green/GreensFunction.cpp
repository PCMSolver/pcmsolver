#include <map>
#include <stdexcept>

#include "GreensFunction.h"

Eigen::Array4d GreensFunction::evaluate(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const
{
        enum                           	 
        {
		NumericalDirectional, 
        	AnalyticDirectional, 
        	AutomaticDirectional, 
        	AutomaticGradient,
        	AutomaticHessian
        };
	std::map<std::string, int> StringToInt;
	StringToInt.insert(std::map<std::string, int>::value_type("Numerical", 0));
	StringToInt.insert(std::map<std::string, int>::value_type("Analytic", 1));
	StringToInt.insert(std::map<std::string, int>::value_type("Derivative", 2));
	StringToInt.insert(std::map<std::string, int>::value_type("Gradient", 3));
	StringToInt.insert(std::map<std::string, int>::value_type("Hessian", 4));

	Eigen::Array4d result = Eigen::Array4d::Zero();

	switch(StringToInt[how])                                                               			
	{                                                                          			
		case NumericalDirectional:                                         			
			result = numericalDirectional(sourceNormal_, source_, probeNormal_, probe_);
			break;
		case AnalyticDirectional:                                          			
			result = analyticDirectional(sourceNormal_, source_, probeNormal_, probe_); 			
			break;
		case AutomaticDirectional:                                         			
			result = automaticDirectional(sourceNormal_, source_, probeNormal_, probe_);			
			break;
		case AutomaticGradient:                                            			
			result = automaticGradient(sourceNormal_, source_, probeNormal_, probe_);   			
			break;
		case AutomaticHessian:                                             			
			result = automaticHessian(sourceNormal_, source_, probeNormal_, probe_);     			
			break;
		default:
			throw std::runtime_error("In GreensFunction.h an unknown Green's function evaluation strategy occurred.");
	}                                                                          			
	return result;                                                             			
}                  

void GreensFunction::compOffDiagonal(const Eigen::Matrix3Xd & elementCenter_, const Eigen::Matrix3Xd & elementNormal_, 
		                     Eigen::MatrixXd & S_, Eigen::MatrixXd & D_) const
{
	int size = S_.rows(); // which is, of course, equal to S_.cols()

	for (int i = 0; i < size; ++i)
	{
		Eigen::Vector3d source = elementCenter_.col(i);
		Eigen::Vector3d sourceNormal = elementNormal_.col(i);
		sourceNormal.normalize();
		for (int j = 0; j < size; ++j)
		{
			Eigen::Vector3d probe = elementCenter_.col(j);
			Eigen::Vector3d probeNormal = elementNormal_.col(j);
			probeNormal.normalize();
			if (i != j)
			{
				Eigen::Array4d tmp = evaluate(sourceNormal, source, probeNormal, probe);
				S_(i, j) = tmp(0);
				D_(i, j) = tmp(1);
			}
		}
	}
}
