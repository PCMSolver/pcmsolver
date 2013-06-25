#ifndef AUTOMATICGRADIENT_H
#define AUTOMATICGRADIENT_H

#include "taylor.hpp"

#include "Evaluate.h"

/*! \file AutomaticGradient.h
 *  \class AutomaticGradient
 *  \brief Implements a directional derivative using the full automatic gradient.
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function and of its directinal
 *  derivative using automatic differentiation.
 *  The full gradient is calculated and contracted with the desired direction:
 *  @\f[ \mathbf{\nabla}_{\mathbf{u}} f(\mathbf{x}) = \mathbf{u}^{t} \mathbf{\nabla} f(\mathbf{x}) @\f]
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
