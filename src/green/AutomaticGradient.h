#ifndef AUTOMATICGRADIENT_H
#define AUTOMATICGRADIENT_H

#include "taylor.hpp"

#include "Evaluate.h"

/*! \file AutomaticGradient.h
 *  \class AutomaticGradient
 *  \brief Implements an automatic gradient. 
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function and of its gradient 
 *  using automatic differentiation.
 *  The automatic differentiation engine is Ulf Ekström's libtaylor library.
 */

class AutomaticGradient : public Evaluate
{
	public:
		virtual void evaluate()
		{
		}
};

#endif // AUTOMATICGRADIENT_H
