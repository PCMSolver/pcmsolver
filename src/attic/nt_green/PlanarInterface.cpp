#include "PlanarInterface.hpp"

#include <cmath>
#include <stdexcept>

double PlanarInterface::evalf(double* p1, double* p2)
{
    if (computed) {
        throw std::runtime_error("Method not implemented yet");
    } else {
        throw std::runtime_error("Compute green´s function before using it!");
    }
}
