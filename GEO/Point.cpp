/*
 * Point.cpp
 *
 *  Created on: Jun 22, 2010
 *      Author: TF
 */


#include <cmath>
#include <limits>

#include "Point.h"

bool operator<= (const GEOLIB::Point& p0, const GEOLIB::Point& p1)
{
	double tol (sqrt (std::numeric_limits<double>::min()));

	if (fabs (p0[0]-p1[0]) > tol * fabs(p0[0])) {
		if (p0[0] < p1[0]) return true;
		else return false;
	} else {
		// assume p0[0] == p1[0]
		if (fabs (p0[1]-p1[1]) > tol * fabs(p0[1])) {
			if (p0[1] < p1[1]) return true;
			else return false;
		} else {
			// assume p0[1] == p1[1] and p0[0] == p1[0]
			if (p0[2] < p1[2]) return true;
			else return false;
		}
	}
}
