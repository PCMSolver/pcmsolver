#ifndef PWLSOLVER_HPP
#define PWLSOLVER_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

extern "C"
{
#include "vector3.h"
#include "sparse2.h"
#include "intvector_pwl.h"
#include "basis_pwl.h"
}

#include "IGreensFunction.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"
#include "WEMSolver.hpp"

/*! \file PWLSolver.hpp
 *  \class PWLSolver
 *  \brief Wavelet solver, piecewise linear.
 *  \author Luca Frediani
 *  \date 2012
 */

class PWLSolver : public WEMSolver
{
private:
    virtual void initPointers();
    virtual std::ostream & printSolver(std::ostream & os);
public:
    PWLSolver(IGreensFunction * gfInside_,
              IGreensFunction * gfOutside_) : WEMSolver(gfInside_, gfOutside_, FirstKind) {
        initPointers();
    }
    PWLSolver(IGreensFunction * gfInside_, IGreensFunction * gfOutside_,
              int integralEquation_) : WEMSolver(gfInside_, gfOutside_, integralEquation_) {
        initPointers();
    }
    //PWLSolver(Section solver);
    virtual ~PWLSolver();
    friend std::ostream & operator<<(std::ostream & os, PWLSolver & solver) {
        return solver.printSolver(os);
    }
private:
    virtual void initInterpolation();
    virtual void constructWavelets();
    virtual void constructSi();
    virtual void constructSe();
    virtual void solveFirstKind(const Eigen::VectorXd & potential,
                                Eigen::VectorXd & charge);
    virtual void solveSecondKind(const Eigen::VectorXd & potential,
                                 Eigen::VectorXd & charge);
    virtual void solveFull(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
    element_pwl *elementTree; //*E_; Hierarchical element list
    wavelet_pwl *waveletList; //*W_; List of wavelets
};

namespace
{
    PCMSolver * createPWLSolver(const solverData & _data)
    {
        return new PWLSolver(_data.gfInside, _data.gfOutside, _data.integralEquation);
    }
    const std::string PWLSOLVER("Linear");
    const bool registeredPWLSolver = SolverFactory::TheSolverFactory().registerSolver(
                                         PWLSOLVER, createPWLSolver);
}

#endif // PWLSOLVER_HPP
