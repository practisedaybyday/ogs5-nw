/*
 * GMSHAdaptiveMeshDensity.h
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
 */

#ifndef GMSHADAPTIVEMESHDENSITY_H_
#define GMSHADAPTIVEMESHDENSITY_H_

// FileIO
#include "GMSHMeshDensityStrategy.h"

// GEOLIB
#include "Point.h"
#include "QuadTree.h"

namespace FileIO {

class GMSHAdaptiveMeshDensity: public GMSHMeshDensityStrategy {
public:
	GMSHAdaptiveMeshDensity(double pnt_density, double station_density, size_t max_pnts_per_leaf);
	virtual ~GMSHAdaptiveMeshDensity();
	void addPoints(std::vector<GEOLIB::Point*> const& pnts);

private:
	double _pnt_density;
	double _station_density;
	size_t _max_pnts_per_leaf;
	GEOLIB::QuadTree<GEOLIB::Point> *_quad_tree;
};

} // end namespace FileIO

#endif /* GMSHADAPTIVEMESHDENSITY_H_ */
