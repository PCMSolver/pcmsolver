#ifndef BOUNDING_BOX_HPP
#define BOUNDING_BOX_HPP
/**
 * @file BoundingBox.hpp
 *
 * @brief the structure declarations of the bounding boxes, for the constant and
 * linear case
 */
typedef struct{
  double rx, ry, rz; ///< the radius of the bounding box in x, y, and z respectively
  double mx, my, mz; ///< coordinates of the midpoint
} BoundingBox;

typedef struct{
  unsigned int patch;
  double rx, ry; ///< the radius if the bounding box in the reference domain in x and y respectively
  double mx, my; ///< the midpoint in the reference domain
} BoundingBoxSquare;

#endif

