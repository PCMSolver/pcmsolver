#ifndef NUMERICALDIRECTIONAL_H
#define NUMERICALDIRECTIONAL_H

#include <Evaluate.h>

/*! \file NumericalDirectional.h
 *  \class NumericalDirectional
 *  \brief Implements a numerical directional derivative.
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function and of its directional
 *  derivative using numerical differentiation formulae.
 */

class NumericalDerivative : public Evaluate
{
	public:
		virtual void evaluate() 
		{
		}
};

#endif // NUMERICALDIRECTIONAL_H
