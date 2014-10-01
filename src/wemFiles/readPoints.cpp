/**
 * @file readPoints.hpp
 *
 * @brief reads the points from file
 */
#ifdef UNITTEST_READPOINTS_CPP
  #include "RDTSCTest.hpp"
  #include <iostream>
#endif

#include "readPoints.hpp"
#include "Vector3.hpp"
#include <stdio.h>
#include <stdlib.h>

void readPoints(const char * fileName, Vector3 ****ppPointsIn, unsigned int *pnoPatch, unsigned int *pnLevels) {
	unsigned int	n = 0;      ///< n*n elemente per patch on level pnLevels
	unsigned int	k;          ///< index
	unsigned int	j1, j2, j3;	///< point indeces
	float		x, y, z;          ///< point coordinates
	FILE		*file;            ///< input file
	
	// open file with geometric data
	file = fopen(fileName,"r");
	if (file == NULL) {
		printf("ERROR: Datei kann nicht geoeffnet werden! \n");
	} else {
		// allocate memory
		fscanf(file,"%d\n%d\n",pnLevels,pnoPatch);

		n = 1<<(*pnLevels);
		/// @note allocate memory in one chunk, seems faster most of the times
		(*ppPointsIn) = (Vector3***) malloc( ((*pnoPatch)*sizeof(Vector3**)) + ((*pnoPatch)*(n+1)*sizeof(Vector3*)) + ((*pnoPatch)*(n+1)*(n+1)*sizeof(Vector3)) );
		
		for(j1 = 0; j1 < *pnoPatch; ++j1){
			(*ppPointsIn)[j1] = (Vector3 **) ((*ppPointsIn) + (*pnoPatch)) + j1 * (n+1);
			for(j2 = 0; j2 < n+1; ++j2){
				(*ppPointsIn)[j1][j2] = (Vector3 *) ((*ppPointsIn) + (*pnoPatch) +(*pnoPatch)*(n+1)) + j1*(n+1)*(n+1) + j2*(n+1);
			}
		}
		
    // read data
		for (k=0; k<*pnoPatch*(n+1)*(n+1); k++) {
			fscanf(file,"%d %d %d %g %g %g\n",&j1,&j3,&j2,&x,&y,&z);
			// point corresponds to (patch,y-index,x-index)
			(*ppPointsIn)[j1][j2][j3] = Vector3(x,y,z);
		}
		fclose(file);
	}
	return;
}

void allocPoints(Vector3 ****ppPointsIn, unsigned int noPatch, unsigned int nLevels) {
    unsigned int n = 1<<(nLevels);
		/// @note allocate memory in one chunk, seems faster most of the times
		(*ppPointsIn) = (Vector3***) malloc( ((noPatch)*sizeof(Vector3**)) + ((noPatch)*(n+1)*sizeof(Vector3*)) + ((noPatch)*(n+1)*(n+1)*sizeof(Vector3)) );
		
		for(unsigned int j1 = 0; j1 < noPatch; ++j1){
			(*ppPointsIn)[j1] = (Vector3 **) ((*ppPointsIn) + (noPatch)) + j1 * (n+1);
			for(unsigned int j2 = 0; j2 < n+1; ++j2){
				(*ppPointsIn)[j1][j2] = (Vector3 *) ((*ppPointsIn) + (noPatch) +(noPatch)*(n+1)) + j1*(n+1)*(n+1) + j2*(n+1);
			}
		}
	return;
}


void freePoints(Vector3 ****ppPointsIn) {
	free(*ppPointsIn);
	return;
}

#ifdef UNITTEST_READPOINTS_CPP
int main() {
	Vector3 ***pPointsIn;       // node points
	unsigned int noPatch;       // number of patches
	unsigned int nLevels;	// 2^M*2^M elements per patch
	
	unsigned long long startTime = 0;
	unsigned long long functionTime = 0;
	startTime = startRDTSC();

	readPoints(&pPointsIn,&noPatch,&nLevels);

	functionTime = stopRDTSCP() - startTime;
	std::cout << functionTime << std::endl;
}
#endif

