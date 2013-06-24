#ifndef ANALYTICALGRADIENT_H
#define ANALYTICALGRADIENT_H

#include "Evaluate.h"

/*! \file AnalyticalGradient.h
 *  \class AnalyticalGradient
 *  \brief Implements an analytical gradient. 
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function and of its gradient 
 *  using analytic differentiation formulae.
 */

class AnalyticalGradient : public Evaluate
{
	public:
		virtual void evaluate()
		{
		}
};

#endif // ANALYTICALGRADIENT_H
