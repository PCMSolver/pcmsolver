#include "RestartCavity.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Config.hpp"


std::ostream & RestartCavity::printCavity(std::ostream & os)
{
    os << "Cavity type: Restart" << std::endl;
    os << "Number of finite elements = " << nElements_;
    return os;
}
