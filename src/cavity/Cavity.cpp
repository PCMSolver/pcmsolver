#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Dense>

#include "Cavity.hpp"

/*

Methods for basic cavity class
written by Krzysztof Mozgawa, 2011

*/


void Cavity::writeOutput(std::string &filename)
{
	std::ofstream output;
	output.open(filename.c_str(), std::fstream::out);
    	output << nElements << std::endl;
    	for(int i=0; i < nElements; i++) 
	{
		output << elementCenter(0,i) << " ";
		output << elementCenter(1,i) << " ";
		output << elementCenter(2,i) << " ";
		output << elementArea(i) << " ";
   	}
   	output.close();
}

std::ostream & operator<<(std::ostream & os, Cavity & cavity) 
{
	return cavity.printObject(os);
}

std::ostream & Cavity::printObject(std::ostream & os) 
{
	os << "Molecular cavity" << std::endl;
	os << "Nr. of tesserae: " << nElements;
        for(int i = 0; i < nElements; i++) 
	{
		os << std::endl;
		os << i+1 << " ";
		os << elementCenter(0,i) << " ";
		os << elementCenter(1,i) << " ";
		os << elementCenter(2,i) << " ";
		os << elementArea(i) << " ";
        }
	return os;
}

void Cavity::setMode(const std::string & type) 
{
	if (type == "Atoms") 
	{
		setMode(Atoms);
	} 
	else if (type == "Implicit") 
	{
		setMode(Implicit);
	} 
	else if (type == "Explicit") 
	{
		setMode(Explicit);
	} 
	else 
	{
		exit(-1);
	}
}

void Cavity::setMode(int type) 
{
	switch (type) 
	{
	case Atoms :
		mode = Atoms;
		break;
	case Implicit :
		mode = Implicit;
		break;
	case Explicit :
		mode = Explicit;
		break;
	default :
		exit(-1);
	}
}
