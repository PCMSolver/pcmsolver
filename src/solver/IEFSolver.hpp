#ifndef IEFSOLVER_HPP
#define IEFSOLVER_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"


class Cavity;
class GreensFunction;

#include "PCMSolver.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"

/*!
 * \file IEFSolver.hpp
 * \class IEFSolver
 * \brief Traditional solver.
 * \author Luca Frediani
 * \date 2011
 */

class IEFSolver : public PCMSolver
{
private:
    bool builtIsotropicMatrix;
    bool builtAnisotropicMatrix;
    bool hermitivitize_;
    Eigen::MatrixXd fullPCMMatrix;
    std::vector<Eigen::MatrixXd> blockPCMMatrix;
    /*! \brief Builds PCM matrix for an anisotropic environment
     *  \param[in] cavity the cavity to be used.
     */
    void buildAnisotropicMatrix(const Cavity & cavity);
    /*! \brief Builds PCM matrix for an isotropic environment
     *  \param[in] cavity the cavity to be used.
     */
    void buildIsotropicMatrix(const Cavity & cavity);
    virtual std::ostream & printSolver(std::ostream & os);
public:
    IEFSolver() {}
    IEFSolver(GreensFunction * gfInside, GreensFunction * gfOutside, bool symm)
        : PCMSolver(gfInside, gfOutside), builtIsotropicMatrix(false),
          builtAnisotropicMatrix(false), hermitivitize_(symm) {}
    virtual ~IEFSolver() {}
    virtual void buildSystemMatrix(const Cavity & cavity);
    virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge,
                            int irrep = 0);
    friend std::ostream & operator<<(std::ostream & os, IEFSolver & solver) {
        return solver.printSolver(os);
    }
};

namespace
{
    PCMSolver * createIEFSolver(const solverData & _data)
    {
        return new IEFSolver(_data.gfInside, _data.gfOutside, _data.hermitivitize);
    }
    const std::string IEFSOLVER("IEFPCM");
    const bool registeredIEFSolver = SolverFactory::TheSolverFactory().registerSolver(
                                         IEFSOLVER, createIEFSolver);
}

#endif // IEFSOLVER_HPP
