#ifndef PWCSOLVER_HPP
#define PWCSOLVER_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

extern "C"
{
#include "vector3.h"
#include "sparse2.h"
#include "intvector.h"
#include "basis.h"
}

#include "IGreensFunction.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"
#include "WEMSolver.hpp"

/*! \file PWCSolver.hpp
 *  \class PWCSolver
 *  \brief Wavelet solver, piecewise constant.
 *  \author Luca Frediani
 *  \date 2012
 */

class PWCSolver : public WEMSolver
{
private:
    virtual void initPointers();
    virtual std::ostream & printSolver(std::ostream & os);
public:
    PWCSolver(IGreensFunction * gfInside_,
              IGreensFunction * gfOutside_) : WEMSolver(gfInside_, gfOutside_, Full) {
        initPointers();
    }
    PWCSolver(IGreensFunction * gfInside_, IGreensFunction * gfOutside_,
              int integralEquation_) : WEMSolver(gfInside_, gfOutside_, integralEquation_) {
        initPointers();
    }
    //PWCSolver(const Section & solver);
    virtual ~PWCSolver();
    friend std::ostream & operator<<(std::ostream & os, PWCSolver & solver) {
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
    element *elementTree; //*E_; Hierarchical element list
    wavelet *waveletList; //*W_; List of wavelets
};

namespace
{
    PCMSolver * createPWCSolver(const solverData & _data)
    {
        return new PWCSolver(_data.gfInside, _data.gfOutside, _data.integralEquation);
    }
    const std::string PWCSOLVER("Wavelet");
    const bool registeredPWCSolver = SolverFactory::TheSolverFactory().registerSolver(
                                         PWCSOLVER, createPWCSolver);
}

#endif // PWCSOLVER_HPP
