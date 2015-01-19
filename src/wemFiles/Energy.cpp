/**
 * @file Energy.cpp
 *
 * @brief calculates the potential energy
 */

#include <stdio.h>
#include "Vector2.hpp"
#include "Cubature.hpp"
#include "GaussSquare.hpp"
#include "GenericAnsatzFunction.hpp"
#include "Constants.hpp"

/*
#include "molecule.h"
double energy(double *u, GenericAnsatzFunction *af) {
	unsigned int	n = 1 << af->nLevels;	///< n*n elements per patch
	unsigned int	zi;                         ///< row index: zi = i1*(n*n)+i2*n+i3
	Cubature        *Q;                       ///< integration formulas for a 2D domain
	Vector2		s;                              ///< left lower corner of element zi
	Vector2		t;                              ///< evaluation points of gauss quadrature
	double		U;                              ///<  evaluation of density in quadrature point
	double		E;                              ///< energy
	double 		h = 1./n;                       ///< step
	double		C = 0;                          ///< surface meassure
	double		D = 0;                          ///< total charge

	// initialization
	initGaussSquare(&Q,af->quadratureLevel_+1);	// quadrature formulas 2D
	
  // calculating the error?
	E = zi = 0;
	for (unsigned int i1=0; i1<af->nPatches; i1++) {
		s.y = 0;
		for (unsigned int i2=0; i2<n; i2++) {
			s.x = 0;
      // rowwise numbering of patch zi = (i1, i2, i3)
			for (unsigned int i3=0; i3<n; i3++) {
				for (unsigned int l=0; l<Q[af->quadratureLevel_].noP; l++){
					t = vector2Add(s,vector2SMul(h,Q[af->quadratureLevel_].xi[l]));
					U = af->calculateUEnergy(u, Q[af->quadratureLevel_].xi[l], zi);
					E += Q[af->quadratureLevel_].weight[l]*U*potmol(af->interCoeff->Chi(t,i1));
					// total charge calculated using integration of charge density
          C += Q[af->quadratureLevel_].weight[l]*U;
				}
				s.x += h;
				zi++;
			}
			s.y += h;
		}
	}
	E = -0.5*h*E;	// correct scaling

	// output
	printf("Computing the energy:            %g\n",E*0.5291772);

	for (unsigned int l=0; l<k; l++) D += alpha[l];	// total charge
	// C = total charge * (epsilon-1)/epsilon
	printf("Testing total charge (= 0):      %g\n",h*C-D*(epsilon-1)/epsilon);
	freeGaussSquare(&Q,af->quadratureLevel_+1);
	return(E);
}
*/

double energy_ext(double *u, double *potential, GenericAnsatzFunction *af)
{
    unsigned int n = 1 << af->nLevels;    /* n*n Patches pro Parametergebiet             */
    unsigned int i1, i2, i3;    /* Laufindizes fuer Ansatzfunktion             */
    unsigned int zi = 0;        /* Zeilenindex hieraus: zi = i1*(n*n)+i2*n+i3  */
    Cubature *Q;                /* Kubatur-Formeln                             */
    unsigned int l;             /* Laufindex fuer Quadratur                    */
    Vector2 s;                  /* Linker, unterer Eckpunkt des Patches zi     */
    Vector2 t;                  /* Auswertepunkte der Gaussquadratur           */
    double E = 0;               /* Energie                                     */
    double h = 1. / n;          /* Schrittweite                                */
    unsigned int index;
    af->quadratureLevel_ = 1;

/* Initialisierung */
    initGaussSquare(&Q, af->quadratureLevel_ + 1);       /* Kubatur-Formeln */

/* Berechnung des Fehlers */
    for(unsigned int i = af->nPatches*(n*n-1)/3; i <af->elementTree.totalSizeElementList; ++i){
      for (l = 0; l < Q[af->quadratureLevel_].noP; l++) {
        t.x = h*(af->elementTree.element[i].index_s+Q[af->quadratureLevel_].xi[l].x);
        t.y = h*(af->elementTree.element[i].index_t+Q[af->quadratureLevel_].xi[l].y);
        index = (af->elementTree.element[i].patch*n*n + af->elementTree.element[i].index_t*n+af->elementTree.element[i].index_s);
        E += Q[af->quadratureLevel_].weight[l] * af->calculateUEnergy(u, Q[af->quadratureLevel_].xi[l],index)
          * potential[Q[af->quadratureLevel_].noP*index+l];
        }
      }
      E = -0.5 * h * E;           /* correct scaling */
      E *= 0.5291772;

/* Datenausgabe */
    printf("PWC Computed energy:            %f10\n", E);
	  freeGaussSquare(&Q,af->quadratureLevel_+1);
    return (E);
}

double charge_ext(double *u, double *charge, GenericAnsatzFunction *af)
{
    unsigned int n = 1 << af->nLevels;    /* n*n Patches pro Parametergebiet             */
    unsigned int i1, i2, i3;    /* Laufindizes fuer Ansatzfunktion             */
    unsigned int zi = 0;        /* Zeilenindex hieraus: zi = i1*(n*n)+i2*n+i3  */
    Cubature *Q;                /* Kubatur-Formeln                             */
    unsigned int l;             /* Laufindex fuer Quadratur                    */
    Vector2 s;                  /* Linker, unterer Eckpunkt des Patches zi     */
    Vector2 t;                  /* Auswertepunkte der Gaussquadratur           */
    double h = 1. / n;          /* Schrittweite                                */
    double C = 0;               /* surface meassure                            */
    unsigned int index;
    af->quadratureLevel_ = 1;

/* Initialisierung */
    initGaussSquare(&Q, af->quadratureLevel_ + 1);       /* Kubatur-Formeln */

/* Berechnung des Fehlers */
	zi = 0;
	for (unsigned int i1=0; i1<af->nPatches; i1++) {
		s.y = 0;
		for (unsigned int i2=0; i2<n; i2++) {
			s.x = 0;
      // rowwise numbering of patch zi = (i1, i2, i3)
			for (unsigned int i3=0; i3<n; i3++) {
				for (unsigned int l=0; l<Q[af->quadratureLevel_].noP; l++){
					t = vector2Add(s,vector2SMul(h,Q[af->quadratureLevel_].xi[l]));
          charge[zi*Q[af->quadratureLevel_].noP+l] = Q[af->quadratureLevel_].weight[l] * af->calculateUEnergy(u, Q[af->quadratureLevel_].xi[l],zi) * h;
          C += charge[zi*Q[af->quadratureLevel_].noP+l];
#ifdef DEBUG
          //printf("%d %g\n",zi*Q[af->quadratureLevel_].noP+l,charge[zi*Q[af->quadratureLevel_].noP+l]);
#endif
				}
				s.x += h;
				zi++;
			}
			s.y += h;
		}
	}

/* Datenausgabe */
//    printf("PWC Computed charge:            %f10\n", C);
	  freeGaussSquare(&Q,af->quadratureLevel_+1);
    return (C);
}
