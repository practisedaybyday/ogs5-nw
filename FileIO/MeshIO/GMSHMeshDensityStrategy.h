/*
 * GMSHMeshDensityStrategy.h
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
 */

#ifndef GMSHMESHDENSITYSTRATEGY_H_
#define GMSHMESHDENSITYSTRATEGY_H_

#include <ostream>
#include <vector>

// GEOLIB
#include "Point.h"

namespace FileIO
{
/**
 * virtual base class GMSHMeshDensityStrategy
 */
class GMSHMeshDensityStrategy
{
public:
	virtual void init(std::vector<GEOLIB::Point*> const&) = 0;
	virtual void getMeshDensityAtPoint(GEOLIB::Point const*const, std::ostream&) const = 0;
};

} // end namespace


#endif /* GMSHMESHDENSITYSTRATEGY_H_ */
