/***************************************************************************
 *  Point.h
 *  Created on: Jan 12, 2010
 *      Author: TF
**************************************************************************/

#ifndef POINT_H_
#define POINT_H_

#include "TemplatePoint.h"

namespace GEOLIB {

/**
 * \ingroup GEOLIB
 */

typedef TemplatePoint<double> Point;
}

/**
 * lexicographic comparison of points
 */
bool operator<= (const GEOLIB::Point& p0, const GEOLIB::Point& p1);

#endif /* POINT_H_ */
