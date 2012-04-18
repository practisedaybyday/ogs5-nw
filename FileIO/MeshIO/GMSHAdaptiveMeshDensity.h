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

namespace GEOLIB {
class Polygon;
}

namespace FileIO {

class GMSHAdaptiveMeshDensity: public GMSHMeshDensityStrategy {
public:
	GMSHAdaptiveMeshDensity(double pnt_density, double station_density, size_t max_pnts_per_leaf);
	virtual ~GMSHAdaptiveMeshDensity();
	void init(std::vector<GEOLIB::Point const*> const& pnts);
	double getMeshDensityAtPoint(GEOLIB::Point const*const pnt) const;
	void addPoints(std::vector<GEOLIB::Point const*> const& pnts);
	double getMeshDensityAtStation(GEOLIB::Point const*const) const;
	void writeSteinerPoints (GEOLIB::Polygon const*const bounding_polygon,
					size_t &pnt_idx_offset, size_t sfc_idx_offset, std::ostream & out) const;

private:
	double _pnt_density;
	double _station_density;
	size_t _max_pnts_per_leaf;
	GEOLIB::QuadTree<GEOLIB::Point> *_quad_tree;
};

} // end namespace FileIO

#endif /* GMSHADAPTIVEMESHDENSITY_H_ */
