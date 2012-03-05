/*
 * GMSHFixedMeshDensityStrategy.h
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
 */

#ifndef GMSHFIXEDMESHDENSITYSTRATEGY_H_
#define GMSHFIXEDMESHDENSITYSTRATEGY_H_

#include "GMSHMeshDensityStrategy.h"

namespace FileIO {

class GMSHFixedMeshDensityStrategy : public GMSHMeshDensityStrategy {
public:
	GMSHFixedMeshDensityStrategy(double mesh_density);
private:
	double _mesh_density;
};

}

#endif /* GMSHFIXEDMESHDENSITYSTRATEGY_H_ */
