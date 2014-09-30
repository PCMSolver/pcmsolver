/**
 * @file readPoints.hpp
 *
 * @brief reads the points from file
 */

 #ifndef READPOINTS_HPP
 #define READPOINTS_HPP
 
class Vector3;

void readPoints(const char * fileName, Vector3 ****ppPointsIn, unsigned int *pnoPatch, unsigned int *pnLevels);
void allocPoints(Vector3 ****ppPointsIn, unsigned int noPatch, unsigned int nLevels);
void freePoints(Vector3 ****ppPointsIn);

#endif
