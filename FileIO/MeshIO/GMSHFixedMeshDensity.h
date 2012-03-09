/*
 * GMSHFixedMeshDensity.h
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
 */

#ifndef GMSHFIXEDMESHDENSITY_H_
#define GMSHFIXEDMESHDENSITY_H_

#include "GMSHMeshDensityStrategy.h"

namespace FileIO {

class GMSHFixedMeshDensity : public GMSHMeshDensityStrategy
{
public:
	GMSHFixedMeshDensity(double mesh_density);
	void init(std::vector<GEOLIB::Point*> const& vec);
	void getMeshDensityAtPoint(GEOLIB::Point const*const, std::ostream&) const;

private:
	double _mesh_density;
};

}

#endif /* GMSHFIXEDMESHDENSITY_H_ */
