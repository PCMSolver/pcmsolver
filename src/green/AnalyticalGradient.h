#ifndef ANALYTICALGRADIENT_H
#define ANALYTICALGRADIENT_H

#include "Evaluate.h"

/*! \file AnalyticalGradient.h
 *  \class AnalyticalGradient
 *  \brief Implements a directional derivative using the full analytical gradient. 
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function and of its directinal
 *  derivative using analytic differentiation formulae.
 *  The full gradient is calculated and contracted with the desired direction:
 *  @f[ \mathbf{\nabla}_{\mathbf{u}} f(\mathbf{x}) = \mathbf{u}^{t}  \mathbf{\nabla} f(\mathbf{x}) @f]
 */

class AnalyticalGradient : public Evaluate
{
	public:
		virtual void evaluate()
		{
		}
};

#endif // ANALYTICALGRADIENT_H
