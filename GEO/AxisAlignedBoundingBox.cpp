/*
 * \file AxisAlignedBoundingBox.cpp
 *
 *  Created on: April 22, 2010
 *      Author: TF
 */

#include <limits>
#include <cstddef>
#include "AxisAlignedBoundingBox.h"

namespace GEOLIB {

AABB::AABB ()
{
	for (std::size_t k(0); k<3; k++) {
		_min[k] = std::numeric_limits<double>::max();
		_max[k] = std::numeric_limits<double>::min();
	}
}

void AABB::update (double x, double y, double z)
{
	if (x < _min[0]) _min[0] = x;
	if (_max[0] < x) _max[0] = x;
	if (y < _min[1]) _min[1] = y;
	if (_max[1] < y) _max[1] = y;
	if (z < _min[2]) _min[2] = z;
	if (_max[2] < z) _max[2] = z;
}

bool AABB::containsPoint (const double *pnt) const
{
	return containsPoint (pnt[0], pnt[1], pnt[2]);
}

bool AABB::containsPoint (double x, double y, double z) const
{
	if (_min[0] <= x && x <= _max[0]) {
		if (_min[1] <= y && y <= _max[1]) {
			if (_min[2] <= z && z <= _max[2]) {
				return true;
			} else return false;
		} else return false;
	} else return false;
}

}
