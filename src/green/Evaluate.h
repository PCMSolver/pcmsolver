#ifndef EVALUATE_H
#define EVALUATE_H

/*! \file Evaluate.h
 *  \class Evaluate
 *  \brief Encapsulated behaviour for Green's functions evaluation.
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This Abstract Base Class provides the basic facility to encapsulate the
 *  evaluation of a Green's function and its derivatives. We are applying the
 *  Strategy Pattern.
 */

class Evaluate
{
	public:
		virtual void evaluate() = 0;
};

#endif // EVALUATE_H
