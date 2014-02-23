#include "IEFSolver.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

#include "Cavity.hpp"
#include "GreensFunction.hpp"

void matrixSymmetryBlocking(Eigen::MatrixXd & matrix, Cavity & cav)
{
	// This function implements the simmetry-blocking of the PCM
	// matrix due to point group symmetry as reported in:
	// L. Frediani, R. Cammi, C. S. Pomelli, J. Tomasi and K. Ruud, J. Comput.Chem. 25, 375 (2003)
	int cavitySize = cav.size();
	int nr_irrep = cav.pointGroup().nontrivialOps() + 1;
	// u is the character table for the group (t in the paper)
	Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nr_irrep, nr_irrep);
	for (int i = 0; i < nr_irrep; ++i)
	{
		for (int j = 0; j < nr_irrep; ++j)
		{
			u(i, j) = Symmetry::parity(i&j);
		}
	}
	// Naming of indices:
	//     a, b, c, d   run over the total size of the cavity (N)
	//     i, j, k, l   run over the number of irreps (n)
	//     p, q, r, s   run over the irreducible size of the cavity (N/n)
	// Instead of forming U (T in the paper) and then perform the dense
	// multiplication, we multiply block-by-block using just the u matrix.
	//      matrix = U * matrix * Ut; U * Ut = Ut * U = id
	// First half-transformation, i.e. first_half = matrix * Ut
	Eigen::MatrixXd first_half = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
	int ntsirr = cav.irreducible_size();
	for (int i = 0; i < nr_irrep; ++i)
	{
		int ioff = i * ntsirr;
		for (int k = 0; k < nr_irrep; ++k)
		{
			int koff = k * ntsirr;
			for (int j = 0; j < nr_irrep; ++j)
			{
				int joff = j * ntsirr;
				double ujk = u(j, k) / nr_irrep;
				for (int p = 0; p < ntsirr; ++p)
				{
					int a = ioff + p;
					for (int q = 0; q < ntsirr; ++q)
					{
						int b = joff + q;
						int c = koff + q;
						first_half(a, c) += matrix(a, b) * ujk;
					}
				}
			}
		}
	}
	// Second half-transformation, i.e. matrix = U * first_half
	matrix.setZero(cavitySize, cavitySize);
	for (int i = 0; i < nr_irrep; ++i)
	{
		int ioff = i * ntsirr;
		for (int k = 0; k < nr_irrep; ++k) 
		{
			int koff = k * ntsirr;
			for (int j = 0; j < nr_irrep; ++j)
			{
				int joff = j * ntsirr;
				double uij = u(i, j);
				for (int p = 0; p < ntsirr; ++p)
				{
					int a = ioff + p;
					int b = joff + p;
					for (int q = 0; q < ntsirr; ++q)
					{
						int c = koff + q;
						matrix(a, c) += uij * first_half(b, c);
					}
				}
			}
		}
	}
	// Traverse the matrix and discard numerical zeros
	for (int a = 0; a < cavitySize; ++a)
	{
		for (int b = 0; b < cavitySize; ++b)
		{
			if ( std::abs(matrix(a, b)) < 1.0e-14 )
			{
				matrix(a, b) = 0.0;	        				
			}
		}
	}
}

void IEFSolver::buildSystemMatrix(Cavity & cavity) 
{
	if (greenInside->isUniform() && greenOutside->isUniform())
	{
		buildIsotropicMatrix(cavity);
	}
	else
	{
		buildAnisotropicMatrix(cavity);
	}
}

void IEFSolver::buildAnisotropicMatrix(Cavity & cav)
{
	int cavitySize = cav.size();
    	Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
	Eigen::MatrixXd SE = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd DE = Eigen::MatrixXd::Zero(cavitySize, cavitySize);

    	// This is the very core of PCMSolver
    	greenInside->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SI, DI);
    	greenInside->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SI, DI);
    	greenOutside->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SE, DE);
    	greenOutside->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SE, DE);
/*    for(int i = 0; i < cavitySize; i++)
    {
	    Eigen::Vector3d p1 = cav.getElementCenter(i); 
		double area = cav.getElementArea(i);
		double radius = cav.getElementRadius(i);
		SI(i,i) =  greenInside->compDiagonalElementS(area); 
		SE(i,i) = greenOutside->compDiagonalElementS(area); 
		DI(i,i) =  greenInside->compDiagonalElementD(area, radius); 
		DE(i,i) = greenOutside->compDiagonalElementD(area, radius); 
		for (int j = 0; j < cavitySize; j++){
			Eigen::Vector3d p2 = cav.getElementCenter(j);
			Eigen::Vector3d n2 = cav.getElementNormal(j);
			n2.normalize();
			if (i != j) {
				SI(i,j) = greenInside->evalf(p1, p2);
				SE(i,j) = greenOutside->evalf(p1, p2);
				DI(i,j) = greenInside->evald(n2, p1, p2);
				DE(i,j) = greenOutside->evald(n2, p1, p2);
			}
		}
    }*/
    	Eigen::MatrixXd a = cav.getElementArea().asDiagonal();// Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd aInv = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	aInv = 2 * M_PI * a.inverse();
    //a.setZero();
    //aInv.setZero();
/*    for (int i = 0; i < cavitySize; i++) {
		a(i,i) = cav.getElementArea(i);
		aInv(i,i) = 2 * M_PI / cav.getElementArea(i);
    }*/
    	PCMMatrix = ((aInv - DE) * a * SI + SE * a * (aInv + DI.transpose().eval()));
    	PCMMatrix = PCMMatrix.inverse();
    	PCMMatrix *= ((aInv - DE) - SE * SI.inverse() * (aInv - DI));
    	PCMMatrix = PCMMatrix * a;
    	builtAnisotropicMatrix = true;
    	builtIsotropicMatrix = false;
}


void IEFSolver::buildIsotropicMatrix(Cavity & cav)
{
    	std::ofstream matrixOut("PCM_matrix");
	bool symmetry = cav.pointGroup().groupInteger();
	double epsilon = greenOutside->getDielectricConstant();
        int cavitySize = cav.size();
    	Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);

    	// This is the very core of PCMSolver
    	greenInside->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SI, DI);
    	greenInside->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SI, DI);
	if (symmetry)
	{
		matrixSymmetryBlocking(DI, cav);
		matrixSymmetryBlocking(SI, cav);
	}

    	Eigen::MatrixXd a = cav.getElementArea().asDiagonal();
    	Eigen::MatrixXd aInv = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	aInv = a.inverse();

	// Tq = -Rv -> q = -(Tinv * R)v
	// T = (2 * M_PI * fact * aInv - DI) * a * SI; R = (2 * M_PI * aInv - DI)
	// PCMMatrix = Tinv * R * a
	// 1. Form T
	double fact = (epsilon + 1.0)/(epsilon - 1.0);
    	PCMMatrix = (2 * M_PI * fact * aInv - DI) * a * SI;
	// 2. Invert T
    	PCMMatrix = PCMMatrix.inverse();
	// 3. Multiply Tinv and R
    	PCMMatrix *= (2 * M_PI * aInv - DI);
	// 4. Multiply by a
    	PCMMatrix = PCMMatrix * a;
	// 5. Symmetrize K := (K + K+)/2
    	Eigen::MatrixXd PCMAdjoint(cavitySize, cavitySize); 
    	PCMAdjoint = PCMMatrix.adjoint().eval(); // See Eigen doc for the reason of this
    	PCMMatrix = 0.5 * (PCMMatrix + PCMAdjoint);
	// Diagonalize and print out some useful info.
	// This should be done only when debugging... 
    	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(PCMMatrix);
    	if (solver.info() != Eigen::Success) 
		abort();
    	matrixOut << "PCM matrix printout" << std::endl;
 	matrixOut << "Number of Tesserae: " << cavitySize << std::endl;
    	matrixOut << "Largest Eigenvalue: " << solver.eigenvalues()[cavitySize-1] << std::endl;
    	matrixOut << "Lowest Eigenvalue: " << solver.eigenvalues()[0] << std::endl;
    	matrixOut << "Average of Eigenvalues: " << (solver.eigenvalues().sum() / cavitySize)<< std::endl;
    	matrixOut << "List of Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    	matrixOut.close();

	builtIsotropicMatrix = true;
	builtAnisotropicMatrix = false;
}

void IEFSolver::compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) 
{
	if (builtIsotropicMatrix or builtAnisotropicMatrix) 
	{
		charge = - PCMMatrix * potential;
	} 
	else 
	{
		throw std::runtime_error("PCM matrix not initialized!");
	}
}

std::ostream & IEFSolver::printSolver(std::ostream & os) 
{
	std::string type;
	if (builtAnisotropicMatrix) 
	{
		type = "IEFPCM, anisotropic";
	} 
	else 
	{
		type = "IEFPCM, isotropic";
	}
	os << "Solver Type: " << type;
	return os;
}

