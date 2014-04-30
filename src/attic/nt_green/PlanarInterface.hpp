#ifndef PLANARINTERFACE_HPP
#define PLANARINTERFACE_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GreensFunction.hpp"

class PlanarInterface : public GreensFunction
{
public:
    PlanarInterface(double eps1_, double eps2_, double pos_, double width_)
        : eps1(eps1_), eps2(eps2_), pos(pos_), width(width_), computed(false) {}
    ~PlanarInterface() {};
    double evalf(double* p1, double* p2);
    friend std::ostream & operator<<(std::ostream & os, PlanarInterface & green) {
        return green.printGreensFunction(os);
    }
private:
    double eps1;
    double eps2;
    double pos;
    double width;
    bool computed;
    virtual std::ostream & printGreensFunction(std::ostream & os);
};

#endif // PLANARINTERFACE_HPP
