#ifndef SOLVERDATA_HPP
#define SOLVERDATA_HPP

#include "Config.hpp"


class IGreensFunction;

/*! @struct solverData
 *  @brief Contains all data defined from user input in the solver section.
 *  @var solverData::gfInside
 *  The Green's function inside the cavity.
 *  @var solverData::gfOutside
 *  The Green's function outside the cavity.
 *  @var solverData::correction
 *  The correction factor to be use in a CPCM calculation.
 *  @var solverData::integralEquation
 *  The type of integral equation to solve, relevant only for wavelet solvers.
 *  @var solverData::hermitivitize
 *  Triggers hermitivitization of the PCM matrix obtained by collocation.
 */

struct solverData {
    IGreensFunction * gfInside;
    IGreensFunction * gfOutside;
    double correction;
    int integralEquation;
    bool hermitivitize;
    solverData(IGreensFunction * _gfInside, IGreensFunction * _gfOutside,
               double _correction = 0.0,  int _integralEquation = 1, bool _symm = true) :
        gfInside(_gfInside), gfOutside(_gfOutside), correction(_correction),
        integralEquation(_integralEquation), hermitivitize(_symm) {}
};

#endif // SOLVERDATA_HPP
