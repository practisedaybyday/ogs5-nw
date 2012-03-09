/*
 * GMSHNoMeshDensity.h
 *
 *  Created on: Mar 5, 2012
 *      Author: fischeth
 */

#ifndef GMSHNOMESHDENSITY_H_
#define GMSHNOMESHDENSITY_H_

#include "GMSHMeshDensityStrategy.h"

namespace FileIO {

class GMSHNoMeshDensity: public FileIO::GMSHMeshDensityStrategy {
public:
	GMSHNoMeshDensity() {};
	void init(std::vector<GEOLIB::Point*> const& vec)
	{
		// to avoid a warning here:
		(void)(vec);
	}

	void getMeshDensityAtPoint(GEOLIB::Point const*const pnt, std::ostream& out) const
	{
		// to avoid a warning here:
		(void)(pnt);
		(void)(out);
	}
};

}

#endif /* GMSHNOMESHDENSITY_H_ */
