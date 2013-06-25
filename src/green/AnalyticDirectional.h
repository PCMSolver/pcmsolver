#ifndef ANALYTICDIRECTIONAL_H
#define ANALYTICDIRECTIONAL_H

#include "Evaluate.h"

/*! \file AnalyticDirectional.h
 *  \class AnalyticDirectional
 *  \brief Implements an analytic directional derivative. 
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function and of its directional
 *  derivative using analytic differentiation formulae:
 *  \f$ \mathbf{\nabla}_{\mathbf{u}} f(\mathbf{x}) = \mathbf{u}^{t}  \mathbf{\nabla} f(\mathbf{x}) \f$
 */

class AnalyticDirectional : public Evaluate
{
	public:
		virtual void evaluate()
		{
		}
};

#endif // ANALYTICDIRECTIONAL_H
